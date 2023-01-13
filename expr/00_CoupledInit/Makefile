# set the default shell to bash
SHELL=/bin/bash

# directory on the scratch file system where jobs will be run
destination:="/home/anolan/scratch/thermal-structure/expr/00_CoupledInit"

local_src   := $(wildcard initialize.*)
scratch_src := $(local_src:%=$(destination)/%)

# list of entire directories which need to be copied 
local_dirs := run/ sifs/ params/
scratch_dirs := $(local_dirs:%=$(destination)/%)

# list of directories which need their older structure copied
local_skl := logs/ result/
scratch_skl := $(local_skl:%=$(destination)/%)

all : $(scratch_dirs) $(scratch_skl)

# folder structures to copy to scratch
$(scratch_skl)  : $(destination)/%/ : %/
	# add filter to exclude everything that isn't a folder
	rsync -a --filter="-! */" $^ $@
.PHONY: $(scrathc_skl)


# Whole directories to copy to scratch
$(scratch_dirs) : $(destination)/%/ : %/ 
	rsync -rupE $^ $@ 
# run rsync every time, even if the folder already exists on scratch
.PHONY: $(scrathc_dirs)

$(destination)/initialize.% : initialize.%
	@echo cp -p $^ $@
