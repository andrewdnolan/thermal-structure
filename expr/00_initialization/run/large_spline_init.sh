#!/bin/bash
#SBATCH --array=1-74%25                  # 74 jobs that run 25 at a time
#SBATCH --job-name=large_spline_init           # base job name for the array
#SBATCH --mem-per-cpu=2250M                     # maximum 2250MMB per job
#SBATCH --time=10:00:00                      # maximum walltime per job
#SBATCH --nodes=1                                  # Only one node is needed
#SBATCH --ntasks=1                                 # These are serial jobs
#SBATCH --mail-type=ALL                            # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu       # email to spend updates too
#SBATCH --output=logs/large_spline_init_%A_%a.out  # standard output
#SBATCH --error=logs/large_spline_init_%A_%a.err   # standard error
# in the previous two lines %A" is replaced by jobID and "%a" with the array index

#Load cedar module file
source ../../config/modulefile.cc.cedar

# Get the command to create run specific .sif file
CREATE=$( sed -n "${SLURM_ARRAY_TASK_ID}p" ./run/large.in )

# create the .sif file, filename is 'returend' (kinda) by the function
SIF=$(./prepare2submit $CREATE)

# Get the command to convert from .result to NetCDF
CONVERT=$( sed -n "${SLURM_ARRAY_TASK_ID}p" ./run/large.out )

# Start the timer
start=$(date +%s.%N)

# Execute the .sif file
ElmerSolver $SIF

# End the timer
end=$(date +%s.%N)

# Execution time of the solver
runtime=$(awk -v start=$start -v end=$end 'BEGIN {print end - start}')

# Log the execution time
../../src/utils/elmer_log.sh $CREATE --ET $runtime

# Do post-processing file conversion
$CONVERT

# remove the sif
rm -f $SIF