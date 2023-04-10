#!/bin/bash
#SBATCH --job-name=dask_gridding
#SBATCH --time=03:00:00           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000M
#SBATCH --mail-type=ALL                      # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu # email to spend updates too
#SBATCH --output=dask_%A_%a.out              # standard output
#SBATCH --error=dask_%A_%a.err               # standard error

# numbers of cores each job in the array will have
export NUM_WORKERS=16
# use a single thread per cpu core
export THREADS_PER_WORKER=1

source ../../config/modulefile.cc.cedar

# load the source functions
source ./periodic_surge.sh

KEY='crmpt12'

# hard code betas, b/c I only want to process a subset of the files
betas=(1.000e-03 7.499e-04 5.623e-04 4.217e-04
       3.162e-04 2.371e-04 1.778e-04)

# set up the dask cluter for the instance of the job array
create_dask_cluster

for beta in ${betas[@]}; do
  # post process the surge files
  for file in $(find ./result/${KEY}/nc/ -name "*B_${beta}*.nc"); do 
      # get the base filename ,with no path info 
      fn="${file##*/}"
      # strip the file extension, to get the runname 
      run_name="${fn%%.nc}"
      # run the post processing commands
      post_proccess
  done 
done
