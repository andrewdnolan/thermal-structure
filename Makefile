# folder to store compiled files
BIN_DIR = bin
# folder with fortran source
SRC_DIR = src/elmer_UDF
# Find all the source files
SRC  := $(wildcard src/elmer_UDF/*.f90)
# compiled executables used in .sif files
EXEC := $(SRC:src/elmer_UDF/%.f90=$(BIN_DIR)/%)
# compiled fitpack source code
F77_OBJ := $(BIN_DIR)/splev.o $(BIN_DIR)/fpbspl.o
# fortran flags
FFLAGS=-fcheck=all

all: $(EXEC) elmer2nc $(BIN_DIR)/mass_balance $(BIN_DIR)/SurfaceBoundary

# compile the net balance boundary conditions functions, while linking to the
# necessary interface and fitpack source code objects
$(BIN_DIR)/mass_balance: $(SRC_DIR)/mass_balance.f90 $(BIN_DIR)/fitpack_interface.o $(F77_OBJ)
	elmerf90 $^ -o $@ -I$(BIN_DIR)

# compile the fitpack source code and f90 interface
$(BIN_DIR)/fitpack_interface.o: $(SRC_DIR)/mass_balance.f90
	$(MAKE) -C include/fitpack

# compile surface boundary condtions and link surf temp module
$(BIN_DIR)/SurfaceBoundary: $(SRC_DIR)/SurfaceBoundary.f90 $(BIN_DIR)/SurfaceTemperature_mod.o
	elmerf90 $^ -o $@ -I$(BIN_DIR)

# compile surface temperature module neeed for surface boundary conditions
$(BIN_DIR)/SurfaceTemperature_mod.o: $(SRC_DIR)/SurfaceTemperature_mod.f90
	elmerf90 -c $^ -o $@ -J $(BIN_DIR)

# compile NetCDF output solver from fgillet
$(BIN_DIR)/NetcdfUGRIDOutputSolver: $(SRC_DIR)/NetcdfUGRIDOutputSolver.f90
	elmerf90 $(FFLAGS) $^ -o $@ `nf-config --fflags --flibs`

# compile the *.F90 files with the `elmerf90` alias
$(BIN_DIR)/%: $(SRC_DIR)/%.f90
	@if [ $@ = "bin/mass_balance" ] || [ $@ = "bin/SurfaceBoundary" ] || [ $@ = "bin/SurfaceTemperature_mod" ] || [ $@ = "bin/NetcdfUGRIDOutputSolver" ]; then\
		continue; \
	 else \
		elmerf90 $^ -o $@ ; \
	 fi

# make the elmer2nc .result parser using thr makefile in it's source folder
elmer2nc: $(wildcard src/elmer2nc/*.f90)
	$(MAKE) -C src/elmer2nc

clean:
	ls bin/* | grep -v "\." | xargs rm && \
	rm bin/*.o bin/*.mod
