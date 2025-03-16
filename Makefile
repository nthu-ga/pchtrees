# Compiler options inclusion
# Compile and Compiler flags
FC     = gfortran
FFLAGS = -O2 -Wall -Wextra -Wpedantic -fimplicit-none

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
TREE_OBJS := defined_types.o kind_numbers.o halo_mass_function.o run_statistics.o 
TREE_OBJS += memory_modules.o tree_routines.o modified_merger_tree.o deltcrit.o memory.o sigmacdm_spline.o interp.o locate.o hyperbolic.o
TREE_OBJS += split_PCH.o ran3.o spline.o make_tree.o indexxx.o transfer_function.o unresolved_mass.o
TREE_OBJS += parameter_file.o runtime_parameters.o power_spectrum_parameters.o cosmological_parameters.o num_pars.o time_parameters.o
TREE_OBJS += commandline.o io.o 
TREE_OBJS += tinytoml.o 

all:	trees.exe

clean:
	\rm -f ./*.o *.o *.mod tree.exe

#Dependencies
memory.o: memory_modules.o
parameter_file.o: tinytoml.o
sigmacdm_spline.o: num_pars.o cosmological_parameters.o power_spectrum_parameters.o parameter_file.o
deltcrit.o: num_pars.o cosmological_parameters.o parameter_file.o
defined_types.o: kind_numbers.o
tree_routines.o: defined_types.o
split_PCH.o: time_parameters.o
trees.o: defined_types.o memory_modules.o tree_routines.o modified_merger_tree.o cosmological_parameters.o runtime_parameters.o

#Rule for making the executable
trees.exe: $(TREE_OBJS) trees.o
	$(FC) $(FFLAGS) -o pchtrees $^ $(HDF5_LIBS)

