import numpy as np
import h5py
import argparse

def add_attributes_to_group(group, attributes):
    """
    Adds attributes to a given HDF5 group from a dictionary.

    Parameters:
    group (h5py.Group): The HDF5 group to add attributes to.
    attributes (dict): Dictionary of attribute names and values.
    """
    for key, value in attributes.items():
        group.attrs[key] = value

def read_parameters_from_input(input_filename):
    """
    Reads attributes from the input file under the Parameters group.

    Parameters:
    input_filename (str): The name of the input HDF5 file.

    Returns:
    dict: Dictionary of attributes from the Parameters group.
    """
    with h5py.File(input_filename, 'r') as input_file:
        params_group = input_file['/Parameters']
        attributes = {}
        
        # Read all attributes in the /Parameters group
        for attr_name, attr_value in params_group.attrs.items():
            attributes[attr_name] = attr_value
        
        return attributes

def convert(input_filename, output_filename):
    """
    Creates a new HDF5 file and populates the cosmology, groupFinder, and simulation
    groups using attributes read from the input file.

    Parameters:
    input_filename (str): The input HDF5 file to read attributes from.
    output_filename (str): The output HDF5 file to create and populate.
    """
    attributes = read_parameters_from_input(input_filename)

    with h5py.File(output_filename, 'w') as f:
        # Set root-level attribute
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

        # Cosmology group attributes mapping from input file
        cosmology_attrs = {
            'HubbleParam': attributes.get('HubbleParam', None),
            'OmegaMatter': attributes.get('Omega0', None),
            'OmegaLambda': attributes.get('OmegaLambda', None),
            'OmegaBaryon': attributes.get('OmegaBaryon', None),
            'powerSpectrumIndex': 1.0,  # Assuming fixed value
            'sigma_8': 0.9,              # Assuming fixed value
            'transferFunction': "CAMB"   # Assuming fixed value
        }
        add_attributes_to_group(cosmology, cosmology_attrs)

        # GroupFinder group attributes mapping from input file
        group_finder.attrs['COMMENT'] = "Group finder parameters."
        group_finder_attrs = {
            'code': "SUBFIND",  # Assuming fixed value
            'linkingLength': 0.2,  # Assuming fixed value
            'minimumParticleNumber': 20  # Assuming fixed value
        }
        add_attributes_to_group(group_finder, group_finder_attrs)

        # Simulation group attributes
        simulation.attrs['COMMENT'] = "Simulation parameters."
        
        # Any attribute that is not in cosmology or groupFinder is stored in simulation
        simulation_attrs = {}
        for attr_name, attr_value in attributes.items():
            if attr_name == 'BoxSize':
                # Special treatment for 'BoxSize' to store it as 'boxSize'
                simulation_attrs['boxSize'] = attr_value
            elif attr_name == 'code':
                # Special treatment for 'code' to store it as 'Gadget-4'
                simulation_attrs['code'] = "Gadget-4"
            elif attr_name not in cosmology_attrs and attr_name not in group_finder_attrs:
                simulation_attrs[attr_name] = attr_value
        
        add_attributes_to_group(simulation, simulation_attrs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create an HDF5 file from an input file with specific group attributes.")
    parser.add_argument("input_filename", help="The input HDF5 file to read attributes from.")
    parser.add_argument("output_filename", help="The output HDF5 file to create.")
    args = parser.parse_args()

    convert(args.input_filename, args.output_filename)

