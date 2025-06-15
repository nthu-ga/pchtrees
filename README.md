# PCH Tree Code

This repository contains a branch of the Extended Press-Schechter (EPS) merger
tree generation code described by [Parkinson, Cole & Helly (2008), MNRAS, 383,
557](https://ui.adsabs.harvard.edu/abs/2008MNRAS.383..557P/abstract). 

The reference version of this code is available from Shaun Cole:
https://astro.dur.ac.uk/~cole/merger_trees/

The version in this repository adds the following features:

- Specification of parametrs via a parameter file and command-line arguments
- Output of trees to a file.

The algorithm is exactly the same, although minor edits have been made to parts
of the original code for readbility, fixing compiler warnings etc.

For the original README, see
[`./docs/README`](https://github.com/nthu-ga/pchtrees/blob/main/docs/README_original).

## Acknowledgement

> [!IMPORTANT]  
> Work that makes use of this code should reference Parkinson, Cole and Helly
2008 MNRAS, 383, 557 and also acknowledge the "GALFORM Team" for making the
code available.

The code in this repository is made available under the same conditions as the
original PCH code. There is, strictly speaking, no need to acknowledge this repository (although it is appreciated). Please note that I am not an author of the orginal code in any way, shape or form. The scientific algorithms in this repository are the work of H. Parkinson, S. Cole and J. Helly.

The code in this repository makes use of a modified version of the TinyTOML parser by
Thomas Marks under the MIT License (see LICENSE). Link:
https://github.com/archermarx/TinyTOML

## Building

> [!IMPORTANT]  
> **The code requires gcc (gfortran) 9.4.0 or later.**

After cloning this repository, the code can be built with `make all`, which
will produce the `./pchtrees` executable. 

Before building, you may want to edit the makefile to set the following
parameters:

* `BUILD_TYPE = OPT` builds optimized code;
* `HDF5_DIR` should be set to the root directory of an HDF5 (fortran) library
  installation (with the libraries under `/lib` below this path.

If `HDF5_DIR` is not set, the code will be built wihtout support for HDF5 output.

## Running

`./pchtrees` : Prints a usage message with information about the command line
options.

`./pchtrees --defaults` : Dumps the default parameters to `STDOUT`. You can
copy this output to create a fiducial parameter file, which you can customize.

Example: `./pchtrees --defaults > my_params.toml`

`./pchtrees parameter_file_path ntrees mphalo ahalo zmax [options]`

Options (positional or by keyword):
* path   (--path  ) : path to parameter file in TOML format
* ntrees (--ntrees) : integer number of trees to generate (1)
* mphalo (--mphalo) : target mass of tree root notes (1e12 Msol)
* ahalo  (--ahalo)  : Expansion factor at root of tree (1.0)
* zmax   (--zmax)   : highest redshift in tree (4.0)

Options (keyword only)
* --mmax (value) : if given (with a value), sample ntreees trees between mphalo and mmax
* --loguniform : if given with mmax, sample uniformly in log10 mass rather than mass

## Mass sampling 

If the command line option "--mmax" is given, the code generates ntrees trees
with random uniform mass in the range mphalo < mass < mmax.

If the option "--loguniform" is also give, random sample in log10 mass instead. 
