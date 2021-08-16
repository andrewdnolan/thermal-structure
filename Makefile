
#https://stackoverflow.com/questions/488333/gfortran-how-to-control-output-directory-for-mod-files
FC=gfortran
FCFLAGS=-O3 -Wall
EXE=bin/elmer2nc
OBJS=bin/main.o bin/nodes_mod.o bin/results_mod.o bin/variable_mod.o


all: $(EXE)

bin/results_mod.o: bin/variable_mod.o
bin/main.o: src/elmer2nc/main.f90 bin/nodes_mod.o bin/results_mod.o bin/variable_mod.o

$(EXE): $(OBJS)
	$(FC) $(FCFLAGS) $(OBJS) -o $@

bin/%.o: src/elmer2nc/%.f90
	gfortran -c $< -o $@ -J bin/

clean:
	rm $(EXE) *.mod bin/*.o
