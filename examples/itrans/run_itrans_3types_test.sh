#!/bin/bash

# All paths relative to executable
rm -v ./output_tree*.hdf5
./pchtrees --ntrees 10 --mphalo 1e12 --nlev 128 --zmax 30.0  --params ./itrans_tabulated_quick_test.toml
./pchtrees --ntrees 10 --mphalo 1e12 --nlev 128 --zmax 30.0  --params ./itrans_tabulated_quick_test_h1.toml
./pchtrees --ntrees 10 --mphalo 1e12 --nlev 128 --zmax 30.0  --params ./itrans_eh_quick_test.toml
./pchtrees --ntrees 10 --mphalo 1e12 --nlev 128 --zmax 30.0  --params ./itrans_eh_quick_test_h1.toml
./pchtrees --ntrees 10 --mphalo 1e12 --nlev 128 --zmax 30.0  --params ./itrans_bbks_quick_test.toml
./pchtrees --ntrees 10 --mphalo 1e12 --nlev 128 --zmax 30.0  --params ./itrans_bbks_quick_test_h1.toml
