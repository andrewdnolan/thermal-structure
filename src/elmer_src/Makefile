#Default compiler
FC=gfortran
#Elmer compiler
EC=elmerf90
# folder to store compiled files
BIN_DIR = ../../bin
# folder with fortran source
SRC_DIR = .
# fortran flags
FFLAGS=-fcheck=all

# List of elmer code to be built
SIMPLE_OBJS := $(BIN_DIR)/SlipCoefficent \
			   $(BIN_DIR)/DummySolver \
			   $(BIN_DIR)/FrictionLoads 
NC_OBJS     := $(BIN_DIR)/NetcdfUGRIDOutputSolver
LINKED_OBJS := $(BIN_DIR)/SurfaceBoundary \
			   $(BIN_DIR)/mass_balance \
			   $(BIN_DIR)/Thermodynamics 
			   
# Modules to be built, and linked to elmer source
MOD_SRCS := $(SRC_DIR)/surfacetemperature_mod.f90 \
						$(SRC_DIR)/PositiveDegreeDays.f90 \

MOD_OBJS := $(patsubst $(SRC_DIR)/%.f90, $(BIN_DIR)/%.o,  $(MOD_SRCS))

# Concat all the objs to one list
OBJS := $(SIMPLE_OBJS) $(NC_OBJS) $(LINKED_OBJS)

# global rule for making
all: $(OBJS) $(MOD_OBJS)

# Flags for targets that need special flags
# ref: https://stackoverflow.com/a/9355185/10221482
$(NC_OBJS)    : FFLAGS := `nf-config --fflags --flibs`
$(LINKED_OBJS): FFLAGS :=  -I$(BIN_DIR)

# Create module files from Fortran source
$(MOD_OBJS): $(BIN_DIR)/%.o : $(SRC_DIR)/%.f90
	$(EC) -c $^ -o $@ -J $(BIN_DIR)

# Create object files from Fortran source
$(OBJS): $(BIN_DIR)/% : $(SRC_DIR)/%.f90
	$(EC) -o $@ $^ $(FFLAGS) -J $(BIN_DIR)

# Define all module interdependencies
special_funcs.mod      := $(BIN_DIR)/special_funcs.o
fitpack_interface.mod  := $(BIN_DIR)/fitpack_interface.o \
													$(BIN_DIR)/fpbspl.o            \
													$(BIN_DIR)/splev.o
surfacetemperature.mod := $(BIN_DIR)/surfacetemperature_mod.o
PositiveDegreeDays.mod := $(BIN_DIR)/PositiveDegreeDays.o\
 													$(BIN_DIR)/special_funcs.o

# Rules for creating executables
$(BIN_DIR)/mass_balance            : $(fitpack_interface.mod)
$(BIN_DIR)/SurfaceBoundary         : $(surfacetemperature.mod) $(PositiveDegreeDays.mod)
