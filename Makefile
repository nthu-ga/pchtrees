# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Commands
RM:=rm -f

# Compiler options inclusion
# Compile and Compiler flags
FC     = gfortran
FFLAGS = -O0 -fbacktrace -Wall -Wextra -Wpedantic -fimplicit-none

HDF5_DIR   = /cluster/software/hdf5/1.10.5/gcc--8.3.0/serial
HDF5_FLAGS = -I$(HDF5_DIR)/include
HDF5_LIBS  = -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5

# General rules
.SUFFIXES:
.SUFFIXES: .exe .o .F90 .f90 .mod .f

.F90.o:
	$(FC) $(FFLAGS) $(HDF5_FLAGS) -c $< -o $*.o 
.f90.o:
	$(FC) $(FFLAGS) $(HDF5_FLAGS) -c $< -o $*.o 
.f.o:
	$(FC) $(FFLAGS) $(HDF5_FLAGS) -c $< -o $*.o 

# Object code list
OBJECTS := defined_types.o kind_numbers.o ran3.o indexxx.o interp.o locate.o hyperbolic.o
OBJECTS += deltcrit.o halo_mass_function.o sigmacdm_spline.o transfer_function.o  
OBJECTS += memory_modules.o memory.o run_statistics.o spline.o
OBJECTS += tree_routines.o modified_merger_tree.o split_PCH.o  make_tree.o unresolved_mass.o
OBJECTS += parameter_file.o runtime_parameters.o power_spectrum_parameters.o cosmological_parameters.o num_pars.o time_parameters.o
OBJECTS += commandline.o io.o 
OBJECTS += tinytoml.o 

# Declare all public targets
.PHONY: all clean
all:	pchtrees
clean:
	$(RM) *.mod *.smod $(OBJECTS) pchtrees

#Dependencies
indexxx.o: num_pars.o
memory.o: memory_modules.o
parameter_file.o: tinytoml.o
sigmacdm_spline.o: num_pars.o cosmological_parameters.o power_spectrum_parameters.o parameter_file.o
deltcrit.o: num_pars.o cosmological_parameters.o parameter_file.o
defined_types.o: kind_numbers.o
tree_routines.o: defined_types.o
split_PCH.o: time_parameters.o
trees.o: defined_types.o memory_modules.o tree_routines.o modified_merger_tree.o cosmological_parameters.o runtime_parameters.o

#Rule for making the executable
pchtrees: $(OBJECTS) trees.o
	$(FC) $(FFLAGS) -o pchtrees $^ $(HDF5_LIBS)

