FC=gfortran
FCFLAGS=-O3 -Wall -Wextra -Wconversion -pedantic
BIN_DIR=../../bin
EXE=elmer2nc
OBJS_local := main.o nodes_mod.o results_mod.o variable_mod.o utils_mod.o
OBJS       := $(foreach obj,$(OBJS_local),$(BIN_DIR)/$(obj))


all: $(BIN_DIR)/$(EXE)

$(BIN_DIR)/nodes_mod.o:     $(BIN_DIR)/utils_mod.o
$(BIN_DIR)/variable_mod.o:  $(BIN_DIR)/utils_mod.o
$(BIN_DIR)/results_mod.o:   $(BIN_DIR)/variable_mod.o $(BIN_DIR)/utils_mod.o
$(BIN_DIR)/main.o: main.f90 $(BIN_DIR)/nodes_mod.o $(BIN_DIR)/results_mod.o $(BIN_DIR)/variable_mod.o $(BIN_DIR)/utils_mod.o

$(BIN_DIR)/$(EXE): $(OBJS)
	$(FC) $(FCFLAGS) $(OBJS) -o $@ `nf-config --fflags --flibs`

$(BIN_DIR)/%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@ -J $(BIN_DIR) `nf-config --fflags --flibs`

clean:
	rm $(EXE) bin/*.mod bin/*.o
