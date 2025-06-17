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
original PCH code. There is, strictly speaking, no need to acknowledge this
repository (although it is appreciated). Please note that I (the author of this
wrapper code)  am *not* an author of the orginal code. All the 'scientific'
algorithms in this repository are the work of H. Parkinson, S. Cole and J.
Helly.

The code in this repository makes use of a modified version of the TinyTOML
parser by Thomas Marks under the MIT License (see LICENSE). Link:
https://github.com/archermarx/TinyTOML

## Building

> [!IMPORTANT]  
> **The code requires gcc (gfortran) 9.4.0 or later.**

After cloning this repository, the code can be built with `make all`, which
will produce the `./pchtrees` executable. 

Before building, you may want to edit the makefile to set the following
parameters:

* `BUILD_TYPE = OPT` builds optimized code (almost certainly this is what you want, unless you are developing the code);
* `HDF5_DIR` should be set to the root directory of an HDF5 (fortran) library
  installation (with the libraries under `/lib` below this path.

If `HDF5_DIR` is not set, the code will be built wihtout support for HDF5 output. Many of the newer features will not work in that case.

## Running

`./pchtrees` : Prints a usage message with information about the command line
options.

`./pchtrees --defaults` : Dumps the default parameters to `STDOUT`. You can
copy this output to create a fiducial parameter file, which you can customize.

Example: `./pchtrees --defaults > my_params.toml`

`./pchtrees parameter_file_path ntrees mphalo ahalo zmax [options]`

Options (positional or by keyword):
* `path   (--path  )` : path to parameter file in TOML format
* `ntrees (--ntrees)` : integer number of trees to generate (1)
* `mphalo (--mphalo)` : target mass of tree root notes (1e12 Msol)
* `ahalo  (--ahalo)`  : Expansion factor at root of tree (1.0)
* `zmax   (--zmax)`   : highest redshift in tree (4.0)

Options (keyword only)
* `--verbose` : print more output (for testing/debugging)
* `--mmax (value)` : if given (with a value), sample ntreees trees between mphalo and mmax (see "Mass sampling")
* `--loguniform` : if given with mmax, sample uniformly in log10 mass rather than mass (see "Mass sampling")
* `--nlev (value)`: if given (with a value, override the number of tree levels in the parameter file
* `--no-output-trees`: if given, do not write merger tree output
* `--process-first-order-progenitors`: see "Progenitor information".

## Mass sampling 

If the command line option `--mmax` is given, the code generates `ntrees` trees
with masses drawn randomly from uniform distribution in the range `mphalo <
mass < mmax`.

If the option `--loguniform` is also given, the random sampling will be uniform
in `log10 mass` instead. 

## Progenitor information

The command line option `--process-first-order-progentiors` writes additional
output. At present, this information concerns information on the "first
order" progenitors of each tree -- those that are accreted directly on to the
main branch.

This is likely most useful in combinations with `--no-output-trees`, to
generate only the first order progenitor information without saving a great
deal of additional data on each tree. If the trees are saved, information
about the progenitors can be retrieved by postprocesing. However, there may
also be cases where it is more efficient to have this code summarise the properties of first-order mergers directly.

The behavior of this option is controlled by a section `[pfop]` in
the paramter file which must be present if this command line option is used.

```
[pfop]
 file_path  = './output_pfop'
 mass_limit = 1e10
```

`file_path`: Location of the PFOP output.
`mass_limit`: Only write progenitors above this mass.

The PFOP output file has the following structure:

`/Header`: _Values stored as group attributes_

- `LastSnapShotNr` : Number of tree output levels
- `Nhalos_ThisFile` : Number of halos (in this case, progenitors / length of `Progenitors` datasets
- `Nhalos_Total` : Number of halos (in this case, progenitors) over all files
- `Ntrees_ThisFile` : Number of trees / length of `TreeTable` datasets
- `Ntrees_Total` : Number of trees over all files
- `NumFiles` : Number of files

`/Mainbranch` 

- `/Mainbranch/MainbranchMass`:  _The mass of the main branch at each output level (NTREES, NLEV)_

`/OutputTimes`: _Each dataset has one row per output level_
	
- `/OutputTimes/DeltaCrit`: _The critical density at each output time_
- `/OutputTimes/ExpansionFactor`: _The expansion factor at each output time_
- `/OutputTimes/Redshift`: _The redshift at each output time_

`/Progenitors`: _Each dataset has one row per merger event, concatenating events in all trees_
	
- `/Progenitors/HostMass`: _The main branch mass immediately before the merger_
- `/Progenitors/MergedMass`: _The main branch mass immediately after the merger_
- `/Progenitors/MergedZred`: _The redshift immediately after the merger_
- `/Progenitors/ProgenitorMass`:  _The progenitor mass immediately before the merger_
- `/Progenitors/ProgenitorZred`: _The redshift immediately before the merger_
- `/Progenitors/TreeID`: _The index in the TreeTable associated with this event_

`/TreeTable`: _Each dataset has one row per tree_
	
- `/TreeTable/NFirstOrderProg`: _The number of first order progenitors in this tree_
- `/TreeTable/RootMass`: _The mass of the main branch at the final output level_
- `/TreeTable/StartOffset`: _The first offset in the `Progenitors` table for the tree_
- `/TreeTable/TreeID`: _The index of this tree (simply the row number)_

