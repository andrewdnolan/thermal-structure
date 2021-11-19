# folder to store compiled files
BIN_DIR = bin
# folder with fortran source
SRC_DIR = src
# Find all the source files
SRC  := $(wildcard src/*.f90)
# compiled executables used in .sif files
EXEC := $(SRC:src/%.f90=$(BIN_DIR)/%)
# compiled fitpack source code
F77_OBJ := $(BIN_DIR)/splev.o $(BIN_DIR)/fpbspl.o

all: $(EXEC) elmer2nc $(BIN_DIR)/mass_balance

# compile the net balance boundary conditions functions, while linking to the
# necessary interface and fitpack source code objects
$(BIN_DIR)/mass_balance: $(SRC_DIR)/mass_balance.f90 $(BIN_DIR)/fitpack_interface.o $(F77_OBJ)
	elmerf90 $^ -o $@ -I$(BIN_DIR)

# compile the fitpack source code and f90 interface
$(BIN_DIR)/fitpack_interface.o: $(SRC_DIR)/mass_balance.f90
	$(MAKE) -C include/fitpack

# compile the *.F90 files with the `elmerf90` alias
$(BIN_DIR)/%: $(SRC_DIR)/%.f90
	@if [ $@ = "bin/mass_balance" ]; then\
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
