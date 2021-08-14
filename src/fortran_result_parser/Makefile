
#https://stackoverflow.com/questions/488333/gfortran-how-to-control-output-directory-for-mod-files
FC=gfortran
FCFLAGS=-O3 -Wall
EXE=parse_result
OBJS=bin/main.o bin/nodes_mod.o bin/results_mod.o bin/variable_mod.o


all: $(EXE)

bin/results_mod.o: bin/variable_mod.o
bin/main.o: src/main.f90 bin/nodes_mod.o bin/results_mod.o bin/variable_mod.o

$(EXE): $(OBJS)
	$(FC) $(FCFLAGS) $(OBJS) -o $@

bin/%.o: src/%.f90
	gfortran -c $< -o $@

clean:
	rm parse_result *.mod bin/*.o
