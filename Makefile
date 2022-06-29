# Global rule for all the various makefiles
all: fitpack_interface numerical_recipes elmer_src elmer2nc

# compile the fitpack source code and f90 interface
fitpack_interface: $(wildcard include/fitpack/*.f90)
	$(MAKE) -C include/fitpack

# compile the numerical recipes source code
numerical_recipes: $(wildcard include/numerical_recipes/*.f90)
	$(MAKE) -C include/numerical_recipes

# compile the elmer source code and link all the modules
elmer_src: $(wildcard src/elmer_src/*.f90)
	$(MAKE) -C src/elmer_src

# make the elmer2nc .result parser using thr makefile in it's source folder
elmer2nc: $(wildcard src/elmer2nc/*.f90)
	$(MAKE) -C src/elmer2nc

clean:
	ls bin/* | grep -v "\." | xargs rm && \
	rm bin/*.o bin/*.mod
