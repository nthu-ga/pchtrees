# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Commands
RM := rm -f

# Build options
BUILD_TYPE := OPT
#BUILD_TYPE := OPT_PROFILE
#BUILD_TYPE:= DEVELOP
#BUILD_TYPE:= DEVELOP_FIXES
#BUILD_TYPE := DEBUG

# Compiler
FC = gfortran
COMPILER   := $(strip $(COMPILER))

# Comment following line to disable HDF5
HDF5_DIR := /cluster/software/hdf5/1.10.5/gcc--9.4.0/serial

# Report build type
$(info BUILD_TYPE = '$(BUILD_TYPE)')

# Watch out for trailing whitespace in variables set by the end user...
ifeq ($(strip $(BUILD_TYPE)), DEBUG)
    FPP_FLAGS += -DDEBUG
endif

ifeq ($(strip $(BUILD_TYPE)), DEBUG)
    FC_FLAGS := -O0 -g -fbacktrace -Wall -Wextra -Wpedantic -fimplicit-none -fbounds-check -fcheck=all
endif

ifeq ($(strip $(BUILD_TYPE)), DEVELOP)
    FC_FLAGS := -O0 -g -fbacktrace -Wno-maybe-uninitialized -Wall -Wextra -Wpedantic -fimplicit-none  -fbounds-check
endif

ifeq ($(strip $(BUILD_TYPE)), DEVELOP_FIXES)
    FPP_FLAGS := -DSTRICT_REAL_EQ
    FC_FLAGS := -O0 -g -fbacktrace -Wno-maybe-uninitialized -Wall -Wextra -Wpedantic -fimplicit-none  -fbounds-check
endif


ifeq ($(strip $(BUILD_TYPE)), OPT)
    FC_FLAGS := -O3 -fimplicit-none
endif

ifeq ($(strip $(BUILD_TYPE)), OPT_PROFILE)
    FC_FLAGS := -O3 -fimplicit-none -pg --coverage
endif


HDF5_DIR := $(strip $(HDF5_DIR))
ifdef HDF5_DIR
     FPP_FLAGS += -DWITH_HDF5
     HDF5_FLAGS := -I$(HDF5_DIR)/include
     HDF5_LIBS  := -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5
     HDF5_INCL_DIR := '-I$(HDF5_DIR)/include'
endif

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

# Declare all public targets
.PHONY: all clean debug prepare

all: prepare pchtrees
	@echo ''
ifdef HDF5_DIR
	@echo 'Built with HDF5 support'
else
	@echo 'Built without HDF5 support'
endif

prepare:
	@mkdir -p $(BUILD_DIR)

clean:
	$(RM) $(addprefix $(BUILD_DIR)/, *.mod *.smod) $(OBJECTS) pchtrees

debug:
	@echo "BUILD_TYPE = $(BUILD_TYPE)"
	@echo "FC_FLAGS = $(FC_FLAGS)"
	@echo "FPP_FLAGS = $(FPP_FLAGS)"
	@echo "HDF5_DIR = $(HDF5_DIR)"
	@echo "SRC_DIR = $(SRC_DIR)"
	@echo
	@echo "SRCS = $(SRCS)"
	@echo
	@echo "OBJS = $(OBJECTS)"

$(OBJECTS): $(BUILD_DIR)/%.o : $(SRC_DIR)/%.F90
	$(FC) -c $(FC_FLAGS) $(FPP_FLAGS) -J $(MOD_DIR) $(HDF5_INCL_DIR) -o $@ $<

# Dependencies
${BUILD_DIR}/defined_types.o: $(addprefix $(BUILD_DIR)/,kind_numbers.o)
${BUILD_DIR}/memory_modules.o: $(addprefix $(BUILD_DIR)/,defined_types.o)
${BUILD_DIR}/io.o: $(addprefix $(BUILD_DIR)/, memory.o memory_modules.o tree_routines.o cosmological_parameters.o runtime_parameters.o time_parameters.o deltcrit.o power_spectrum_parameters.o)
${BUILD_DIR}/indexxx.o: $(addprefix $(BUILD_DIR)/, num_pars.o)
${BUILD_DIR}/memory.o: $(addprefix $(BUILD_DIR)/, memory_modules.o)
${BUILD_DIR}/interp.o: $(addprefix $(BUILD_DIR)/, real_comparison.o)
${BUILD_DIR}/transfer_function.o: $(addprefix $(BUILD_DIR)/, real_comparison.o)
${BUILD_DIR}/parameter_file.o: $(addprefix $(BUILD_DIR)/, tinytoml.o)
${BUILD_DIR}/sigmacdm_spline.o: $(addprefix $(BUILD_DIR)/, num_pars.o cosmological_parameters.o power_spectrum_parameters.o parameter_file.o)
${BUILD_DIR}/deltcrit.o: $(addprefix $(BUILD_DIR)/, num_pars.o cosmological_parameters.o parameter_file.o real_comparison.o)
${BUILD_DIR}/tree_routines.o: $(addprefix $(BUILD_DIR)/, defined_types.o)
${BUILD_DIR}/split_PCH.o: $(addprefix $(BUILD_DIR)/, time_parameters.o run_statistics.o real_comparison.o)
${BUILD_DIR}/make_tree.o: $(addprefix $(BUILD_DIR)/, run_statistics.o real_comparison.o)
${BUILD_DIR}/trees.o: $(addprefix $(BUILD_DIR)/, defined_types.o memory_modules.o tree_routines.o modified_merger_tree.o cosmological_parameters.o runtime_parameters.o)

# Rule for making the executable
pchtrees: $(OBJECTS)
	$(FC) $(FC_FLAGS) $(FPP_FLAGS) -o pchtrees $^ $(HDF5_LIBS)
