#####################################
## libraries and compilers
## update to your installation
####################################
# NETCDF INC and LIB dirs
NETCDF_INC=/usr/local/include
NETCDF_LIB=/usr/local/lib

# fortrangis INC and LIB dir
PROJ_INC=/usr/local/include
PROJ_LIB=/usr/local/lib

# fortran compiler
F90=gfortran
# fortran flag
FFLAGS=-fcheck=all

# source directory
SRC=../elmersrc

#######################################
LIBS= NetcdfUGRIDOutputSolver USFIcyMask

NetcdfUGRIDOutputSolver: NetcdfUGRIDOutputSolver.F90
	elmerf90 $(FFLAGS) $^ -o $@ `nf-config --fflags --flibs`

clean:
	rm -f $(LIBS)
