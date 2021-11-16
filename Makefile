# folder to store compiled files
BIN_DIR = bin
# folder with fortran source
SRC_DIR = src
# Find all the source files
SRC  := $(wildcard src/*.f90)
# compiled executables used in .sif files
EXEC := $(SRC:src/%.f90=$(BIN_DIR)/%)

all: fitpack $(EXEC) elmer2nc

# compile the fitpack source code and f90 interface
fitpack: $(wildcard src/elmer2nc/*.f90)
	$(MAKE) -C include/fitpack

# compile the *.F90 files with the `elmerf90` alias
$(BIN_DIR)/%: $(SRC_DIR)/%.f90
	@if [ $@ = "bin/mass_balance" ]; then\
		elmerf90 $^ -o $@ $(BIN_DIR)/fitpack_interface.o -I$(BIN_DIR); \
	 else \
		elmerf90 $^ -o $@ ; \
	 fi

# make the elmer2nc .result parser using thr makefile in it's source folder
elmer2nc: $(wildcard src/elmer2nc/*.f90)
	$(MAKE) -C src/elmer2nc

clean:
	ls bin/* | grep -v "\." | xargs rm && \
	rm bin/*.o bin/*.mod
