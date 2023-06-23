#!/bin/bash
#SBATCH --job-name=netcdf_io_gridding
#SBATCH --time=00:20:00           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000M
#SBATCH --mail-type=ALL                      # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu # email to spend updates too
#SBATCH --output=netcdf_io_gridding_%A.out              # standard output
#SBATCH --error=netcdf_io_gridding_%A.err               # standard error


# numbers of cores each job in the array will have
export NUM_WORKERS=16
# use a single thread per cpu core
export THREADS_PER_WORKER=1

source ../../config/modulefile.cc.cedar

# load the source functions
source ./periodic_surge.sh


KEY='crmpt12'

# set up the dask cluter for the instance of the job array
create_dask_cluster

run_name="crmpt12_dx_50_TT_0--3ka_MB_-0.37_OFF_Tma_-8.5_B_1.778e-04_SP_2_QP_58"

# run the post processing commands
post_proccess