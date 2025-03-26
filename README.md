# PCH Tree Code

This repository contains a branch of the Extended Press-Schechter (EPS) merger tree generation code described by Parkinson, Cole & Helly (2008), MNRAS, 383, 557. 

The reference version of this code is available from Shaun Cole:
https://astro.dur.ac.uk/~cole/merger_trees/

The version in this repository adds the following features:

- Specification of parametrs via a parameter file and command-line arguments
- Output of trees to a file.

The algorithm is exactly the same, although minor edits have been made to parts of the original code for readbility, fixing compiler warnings etc.

For the original README, see [`./docs/README`](https://github.com/nthu-ga/pchtrees/blob/main/docs/README_original).

## Acknowledgement

Work that makes use of this code should reference Parkinson, Cole and Helly 2008 MNRAS, 383, 557 and also acknowledge the "GALFORM Team" for making the code available.

The code in this repository is made available under the same conditions as the original PCH code. 

Note that the code in this repository maskes use of the TinyTOML parser by Thomas Marks under the MIT License (see LICENSE). Link: https://github.com/archermarx/TinyTOML

## Building

**The code requires gcc (gfortran) 9.4.0 or later.**

After cloning this repository, the code can be built with `make all`, which will produce the `./pchtrees` executable. 

Before building, you may want to edit the makefile to set the following parameters:

* `BUILD_TYPE = OPT` builds optimized code;
* `HDF5_DIR` should be set to the root directory of an HDF5 (fortran) library installation (with the libraries under `/lib` below this path.

If `HDF5_DIR` is not set, the code will be built wihtout support for HDF5 output.

## Running

`./pchtrees` alone prints a usage message with information about the command line options.

`./pchtrees --defaults` dumps the default parameters to `STDOUT`. You can copy this output to create a fiducial parameter file, which you can customize.




