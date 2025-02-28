#!/bin/bash
#SBATCH --job-name=dask_gridding
#SBATCH --time=00:40:00           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000M
#SBATCH --mail-type=END                      # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu # email to spend updates too
#SBATCH --output=dask_%A_%a.out              # standard output
#SBATCH --error=dask_%A_%a.err               # standard error

# use a single thread per cpu core
export THREADS_PER_WORKER=8

source ../../config/modulefile.cc.cedar

# load the source functions
source ./surge2steady.sh

KEY='crmpt12'

# make the thinned dir for each task in the job array
mkdir "${SLURM_TMPDIR}/thinned"

# check to make sure a thinned dir exists to dump the tar files
if [ ! -d "result/${KEY}/thinned/" ]; then
  mkdir  "result/${KEY}/thinned/"
fi

# post process the surge files
for file in $(find ./result/${KEY}/nc/ -name "*.nc" -ctime -1); do 
    # get the base filename ,with no path info 
    fn="${file##*/}"
    # strip the file extension, to get the runname 
    run_name="${fn%%.nc}"

    # run the post processing commands
    post_proccess
done 

