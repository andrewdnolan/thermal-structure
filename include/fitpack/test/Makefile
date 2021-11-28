FC=gfortran
FCFLAGS=-O3 -Wall -Wextra -Wconversion -pedantic
BIN_DIR=bin

#gfortran -c splev.f -o bin/splev.o
#gfortran -c fpbspl.f -o bin/fpbspl.o
#gfortran -c fitpack_interface.f90 -o bin/fitpack_interface.o -J bin/ bin/*.o
#gfortran test.f90 -o test bin/fitpack_interface.o bin/fpbspl.o bin/splev.o -I./bin/

OBJS_local := fpbspl.o splev.o
OBJS       := $(foreach obj,$(OBJS_local),$(BIN_DIR)/$(obj))
MOD         = $(BIN_DIR)/fitpack_interface.mod

test: $(BIN_DIR)/test_fitpack $(OBJS) $(MOD)

$(BIN_DIR)/test_fitpack: test_fitpack.f90 $(BIN_DIR)/fitpack_interface.mod
	gfortran test_fitpack.f90 -o $(BIN_DIR)/test_fitpack -J $(BIN_DIR) $(OBJS)

# Complile and link interface module with f77 code
$(BIN_DIR)/fitpack_interface.mod: ../fitpack_interface.f90 $(OBJS)
	gfortran -c ../fitpack_interface.f90 -o $(BIN_DIR)/fitpack_interface.o \
				   -J $(BIN_DIR) -shared $(OBJS)

# Compile and assmble, but do not link the fitpack source code (f77 code)
$(BIN_DIR)/%.o: ../%.f
	gfortran -c $< -o $@

clean:
	ls bin/* | grep -v "\." | xargs rm && \
	rm bin/*.o bin/*.mod
