import h5py
import numpy as np
import os
import argparse

# Constants
BOX_SIZE = 100.0  # Set box size value for forestWeight adjustment

# Input to output mapping
forest_halos_mapping = {
    "TreeFirstDescendant": "descendantIndex",
    "TreeFirstHaloInFOFgroup": "hostIndex",
    "TreeIndex": "nodeIndex",
    "SubhaloMass": "nodeMass",
    "SubhaloLen": "particleCount",
    "SubhaloPos": "position",
    "SubhaloSpin": "spin",
    "SubhaloVel": "velocity",
    "SubhaloHalfmassRad": "halfMassRadius"
}

output_fields_1d = [
    "descendantIndex", "hostIndex", "nodeIndex",
    "nodeMass", "particleCount", "halfMassRadius"
]

# Placeholder for forest index mappings
forest_index_mapping = {
    "Length": "numberOfNodes",
    "StartOffset": "firstNode",
    "TreeID": "forestIndex"
}


def add_attributes_to_group(group, attributes):
    for key, value in attributes.items():
        group.attrs[key] = value

def read_parameters_from_input(input_filename):
    with h5py.File(input_filename, 'r') as input_file:
        params_group = input_file['/Parameters']
        return dict(params_group.attrs)

def convert(input_filename, output_filename):
    attributes = read_parameters_from_input(input_filename)

    with h5py.File(output_filename, 'w') as f:
        f.attrs['formatVersion'] = 2

        # Create required groups
        cosmology = f.create_group('cosmology')
        group_finder = f.create_group('groupFinder')
        simulation = f.create_group('simulation')
        f.create_group('forestHalos')
        f.create_group('forests')
        f.create_group('particles')
        f.create_group('provenance')
        f.create_group('forestIndex')
        f.create_group('units')

        # Track attributes that have been used
        used_attributes = set()

        # Populate cosmology group
        cosmology_mapping = {
            'HubbleParam': 'HubbleParam',
            'Omega0': 'OmegaMatter',
            'OmegaLambda': 'OmegaLambda',
            'OmegaBaryon': 'OmegaBaryon'
        }
        cosmology_attrs = {}
        for input_key, output_key in cosmology_mapping.items():
            if input_key in attributes:
                cosmology_attrs[output_key] = attributes[input_key]
                used_attributes.add(input_key)
        # Add hardcoded values
        cosmology_attrs.update({
            'powerSpectrumIndex': 1.0,
            'sigma_8': 0.9,
            'transferFunction': "CAMB"
        })
        add_attributes_to_group(cosmology, cosmology_attrs)

        # Populate groupFinder group
        group_finder.attrs['COMMENT'] = "Group finder parameters."
        group_finder_attrs = {
            'code': "SUBFIND",
            'linkingLength': 0.2,
            'minimumParticleNumber': 20
        }
        add_attributes_to_group(group_finder, group_finder_attrs)

        # Populate simulation group
        simulation.attrs['COMMENT'] = "Simulation parameters."
        simulation_attrs = {'code': 'Gadget-4'}  # hardcoded

        for attr_name, attr_value in attributes.items():
            if attr_name in used_attributes:
                continue  # skip already used attributes
            if attr_name == 'BoxSize':
                simulation_attrs['boxSize'] = attr_value  # special case
            else:
                simulation_attrs[attr_name] = attr_value  # everything else

        add_attributes_to_group(simulation, simulation_attrs)

def read_input_file(input_filename):
    """
    Reads an input HDF5 file, returns the required datasets and header information.
    """
    with h5py.File(input_filename, 'r') as f:
        header = f['/Header']
        nhalos = header.attrs['Nhalos_Total']
        ntrees = header.attrs['Ntrees_Total']
        num_files = header.attrs['NumFiles']

        # Read TreeTable group (mapping input to output)
        tree_data = None
        if '/TreeTable' in f:
            tree_table = f['/TreeTable']
            tree_data = {}
            for input_field, output_field in forest_index_mapping.items():
                if input_field in tree_table:
                    tree_data[output_field] = tree_table[input_field][:]
        
        # Read TreeHalos group (start of forest-related datasets)
        tree_halos = f['/TreeHalos']
        halo_data  = {}
        for input_field, output_field in forest_halos_mapping.items():
            if input_field in tree_halos:
                halo_data[output_field] = tree_halos[input_field][:]
        
        return halo_data, tree_data, nhalos, ntrees, num_files

def process_forest_index(output_file, all_forest_data):
    """
    Process and append data to forestIndex group in the output file.
    """
    with h5py.File(output_file, 'a') as f:
        if 'forestIndex' not in f:
            forest_group = f.create_group('forestIndex')
        else:
            forest_group = f['forestIndex']

        for data in all_forest_data:
            print(data.keys())

        # Concatenate all forest data across multiple input files
        first_node_all      = np.concatenate([data['firstNode'] for data in all_forest_data])
        number_of_nodes_all = np.concatenate([data['numberOfNodes'] for data in all_forest_data])
        forest_index_all    = np.concatenate([data['forestIndex'] for data in all_forest_data])
        forest_weight_all   = np.concatenate([np.ones_like(data['forestIndex'], dtype=np.float64) / BOX_SIZE for data in all_forest_data])

        forest_group.create_dataset('firstNode',     data=first_node_all)
        forest_group.create_dataset('numberOfNodes', data=number_of_nodes_all)
        forest_group.create_dataset('forestIndex',   data=forest_index_all)
        forest_group.create_dataset('forestWeight',  data=forest_weight_all)

def process_forest_halos(output_file, all_tree_data):
    """
    Process and append data to forestHalos group in the output file.
    """
    with h5py.File(output_file, 'a') as f:
        if 'forestHalos' not in f:
            forest_group = f.create_group('forestHalos')
        else:
            forest_group = f['forestHalos']

        # Concatenate all tree data across multiple input files
        for output_field in output_fields_1d:
            data_all = np.concatenate([data[output_field] for data in all_tree_data])
            forest_group.create_dataset(output_field, data=data_all)

def generate_input_filenames(base_filename):
    """
    Generates the list of input filenames based on the base filename.
    The number of files is read from the header of the first file.
    """
    with h5py.File(base_filename, 'r') as f:
        # Read the number of files from the header
        num_files = f['/Header'].attrs['NumFiles']

    # Extract base name and extension from the file path
    dirname, basename = os.path.split(base_filename)
    base_name_idx, ext = os.path.splitext(basename)
    base_name, idx = os.path.splitext(base_name_idx)
    base_name = base_name.strip('.')
    ext       = ext.strip('.')

    # Generate the list of input filenames
    input_filenames = [os.path.join(dirname, f"{base_name}.{i}.{ext}") for i in range(num_files)]
    return input_filenames

def read_and_process_files(input_base_filename, output_filename):
    """
    Reads and processes multiple input files, then writes the combined output.
    """
    # Generate the full list of input filenames based on the base input filename
    print('Reading metadata from {:s}'.format(input_base_filename))

    input_filenames = generate_input_filenames(input_base_filename)

    all_tree_data = []
    all_forest_data = []

    # Read all input files
    for input_filename in input_filenames:
        print('Reading {:s}'.format(input_filename))
        tree_data, forest_data, nhalos, ntrees, num_files = read_input_file(input_filename)
        if tree_data is not None:
            all_tree_data.append(tree_data)
        if forest_data is not None:
            all_forest_data.append(forest_data)

    # Process and append forestIndex data
    process_forest_index(output_filename, all_forest_data)
    
    # Process and append forestHalos data
    process_forest_halos(output_filename, all_tree_data)

if __name__ == '__main__':
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process and combine HDF5 input files into one output file.')
    parser.add_argument('input_base_filename', type=str, help='Base name of the input files (e.g., input.0.hdf5)')
    parser.add_argument('output_filename', type=str, help='Output filename to write combined data to')
    args = parser.parse_args()

    # Read and process files
    read_and_process_files(args.input_base_filename, args.output_filename)
