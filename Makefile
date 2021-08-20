
#https://stackoverflow.com/questions/488333/gfortran-how-to-control-output-directory-for-mod-files
FC=gfortran
FCFLAGS=-O3 -Wall -Wextra -Wconversion -pedantic
EXE=bin/elmer2nc
OBJS=bin/main.o bin/nodes_mod.o bin/results_mod.o bin/variable_mod.o bin/utils_mod.o


all: $(EXE)

bin/nodes_mod.o: bin/utils_mod.o
bin/variable_mod.o: bin/utils_mod.o
bin/results_mod.o: bin/variable_mod.o bin/utils_mod.o
bin/main.o: src/elmer2nc/main.f90 bin/nodes_mod.o bin/results_mod.o bin/variable_mod.o bin/utils_mod.o

$(EXE): $(OBJS)
	$(FC) $(FCFLAGS) $(OBJS) -o $@ `nf-config --fflags --flibs`

bin/%.o: src/elmer2nc/%.f90
	gfortran -Wall -Wextra -Wconversion -pedantic -c $< -o $@ -J bin/ `nf-config --fflags --flibs`

clean:
	rm $(EXE) bin/*.mod bin/*.o
