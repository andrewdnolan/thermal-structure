#!/bin/bash
#SBATCH --job-name=dask_gridding
#SBATCH --time=01:00:00           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8000M
#SBATCH --mail-type=ALL                      # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu # email to spend updates too
#SBATCH --output=dask_%A_%a.out              # standard output
#SBATCH --error=dask_%A_%a.err               # standard error

# numbers of cores each job in the array will have
export NUM_WORKERS=32
# use a single thread per cpu core
export THREADS_PER_WORKER=1

source ../../config/modulefile.cc.cedar

# load the source functions
source ./periodic_surge.sh

KEY='crmpt12'

# # hard code betas, b/c I only want to process a subset of the files
# betas=(4.019e-04 3.831e-04 3.652e-04 3.481e-04 3.318e-04)

# set up the dask cluter for the instance of the job array
create_dask_cluster


base="crmpt12_dx_50_TT_15000.0_MB_-0.37_OFF_Tma_-8.5_B_1.000e-03_SP_2_QP_58.nc"
file=$(find ./result/${KEY}/nc/ -name "${base}") 

# get the base filename ,with no path info 
fn="${file##*/}"
# strip the file extension, to get the runname 
run_name="${fn%%.nc}"
# run the post processing commands
post_proccess


# for beta in ${betas[@]}; do
#     # post process the surge files
#     for file in $(find ./result/${KEY}/nc/ -name "crmpt12_dx_50_TT_6000.0*B_${beta}*.nc"); do 
#         # get the base filename ,with no path info 
#         fn="${file##*/}"
#         # strip the file extension, to get the runname 
#         run_name="${fn%%.nc}"
#         # run the post processing commands
#         post_proccess
#     done 
# done
