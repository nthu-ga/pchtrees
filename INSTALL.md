
### Quick start

- Edit HDF5 path in makefile
- make (see below)
- run `./pchtrees` with no argumetns for usage
- run `./pchtrees --defaults > parameters.toml` to generate an example parameter file.

Edit the example parameter file as needed and run with

`./pchtress

### Building

The code requires a Fortran 2003 compiler.

Compilers known to work:     gfortran 9.4.0
Compilers known not to work: gfortran 8.3.0

### Third party libraries included

- TinyTOML (Thomas A. Marks) https://github.com/archermarx/TinyTOML
