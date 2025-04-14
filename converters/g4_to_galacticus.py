import numpy as np
import h5py
import argparse

import h5py
import argparse

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create an HDF5 file from an input file with specific group attributes.")
    parser.add_argument("input_filename", help="The input HDF5 file to read attributes from.")
    parser.add_argument("output_filename", help="The output HDF5 file to create.")
    args = parser.parse_args()

    convert(args.input_filename, args.output_filename)
