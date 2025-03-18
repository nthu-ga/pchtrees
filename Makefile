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

# List of all source files
SRC_DIR := ./src
SRCS := $(notdir $(wildcard $(SRC_DIR)/*.F90))

# List of all object files
BUILD_DIR := ./build
OBJECTS := $(addprefix $(BUILD_DIR)/, $(addsuffix .o, $(basename $(SRCS))))

MOD_DIR := ./build

$(OBJECTS): $(BUILD_DIR)/%.o : $(SRC_DIR)/%.F90
	$(FC) -c $(FFLAGS) -J $(MOD_DIR) -I$(HDF5_DIR)/include -o $@ $<

# Object code list
#OBJECTS := defined_types.o kind_numbers.o ran3.o indexxx.o interp.o locate.o hyperbolic.o
#OBJECTS += deltcrit.o halo_mass_function.o sigmacdm_spline.o transfer_function.o  
#OBJECTS += memory_modules.o memory.o run_statistics.o spline.o
#OBJECTS += tree_routines.o modified_merger_tree.o split_PCH.o  make_tree.o unresolved_mass.o
#OBJECTS += parameter_file.o runtime_parameters.o power_spectrum_parameters.o cosmological_parameters.o num_pars.o time_parameters.o
#OBJECTS += commandline.o io.o 
#OBJECTS += tinytoml.o 

# Declare all public targets
.PHONY: all clean debug
all:	pchtrees
clean:
	$(RM) $(addprefix $(BUILD_DIR)/, *.mod *.smod) $(OBJECTS) pchtrees

debug:
	@echo "SRC_DIR = $(SRC_DIR)"
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJECTS)"
	#@echo "MODS = $(MODS)"
	#@echo "MOD_OBJS = $(MOD_OBJS)"
	#@echo "PROGRAM = $(PROGRAM)"
	#@echo "PRG_OBJ = $(PRG_OBJ)"

#Dependencies
${BUILD_DIR}/defined_types.o: $(addprefix $(BUILD_DIR)/,kind_numbers.o)
${BUILD_DIR}/memory_modules.o: $(addprefix $(BUILD_DIR)/,defined_types.o)
${BUILD_DIR}/io.o: $(addprefix $(BUILD_DIR)/, memory.o memory_modules.o tree_routines.o cosmological_parameters.o runtime_parameters.o time_parameters.o deltcrit.o)
${BUILD_DIR}/indexxx.o: $(addprefix $(BUILD_DIR)/, num_pars.o)
${BUILD_DIR}/memory.o: $(addprefix $(BUILD_DIR)/, memory_modules.o)
${BUILD_DIR}/parameter_file.o: $(addprefix $(BUILD_DIR)/, tinytoml.o)
${BUILD_DIR}/sigmacdm_spline.o: $(addprefix $(BUILD_DIR)/, num_pars.o cosmological_parameters.o power_spectrum_parameters.o parameter_file.o)
${BUILD_DIR}/deltcrit.o: $(addprefix $(BUILD_DIR)/, num_pars.o cosmological_parameters.o parameter_file.o)
${BUILD_DIR}/tree_routines.o: $(addprefix $(BUILD_DIR)/, defined_types.o)
${BUILD_DIR}/split_PCH.o: $(addprefix $(BUILD_DIR)/, time_parameters.o)
${BUILD_DIR}/trees.o: $(addprefix $(BUILD_DIR)/, defined_types.o memory_modules.o tree_routines.o modified_merger_tree.o cosmological_parameters.o runtime_parameters.o)

#Rule for making the executable
pchtrees: $(OBJECTS) 
	$(FC) $(FFLAGS) -o pchtrees $^ $(HDF5_LIBS)

