#!/usr/bin/env bash

# Example https://github.com/NOAA-EMC/global-workflow/tree/develop/modulefiles

# Load Elmer and neccesary prereqs
module load StdEnv/2020
module load gcc/9.3.0
module load openmpi/4.0.3
module load elmerfem/9.0

# for netcdf postprocessing in fortran
module load netcdf/4.7.4
module load netcdf-fortran/4.6.0

# for json parsing from the command line
module load perl

# Load python 3
module load python/3.11.5
module load scipy-stack/2023b
source $HOME/python_envs/py4elmer/bin/activate
