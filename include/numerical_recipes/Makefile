#Directories
BIN_DIR=../../bin
SRC_DIR=.

# List of all source files
SRCS := $(SRC_DIR)/special_funcs.f90
OBJS := $(patsubst $(SRC_DIR)/%.f90, $(BIN_DIR)/%.o,  $(SRCS))

# Complile and link module
$(OBJS) : $(SRCS)
	elmerf90 -c $< -o $@ -J $(BIN_DIR)
