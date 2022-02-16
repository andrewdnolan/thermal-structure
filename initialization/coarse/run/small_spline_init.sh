#!/bin/bash
#SBATCH --array=1-55%10                  # 55 jobs that run 10 at a time
#SBATCH --job-name=small_spline_init           # base job name for the array
#SBATCH --mem-per-cpu=300MB                     # maximum 300MBMB per job
#SBATCH --time=1:45:00                      # maximum walltime per job
#SBATCH --nodes=1                                  # Only one node is needed
#SBATCH --ntasks=1                                 # These are serial jobs
#SBATCH --mail-type=ALL                            # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu       # email to spend updates too
#SBATCH --output=logs/small_spline_init_%A_%a.out  # standard output
#SBATCH --error=logs/small_spline_init_%A_%a.err   # standard error
# in the previous two lines %A" is replaced by jobID and "%a" with the array index

#Load cedar module file
source ../../config/modulefile.cc.cedar

# Get the command to create run specific .sif file
CREATE=$( sed -n "${SLURM_ARRAY_TASK_ID}p" <in_fp> )

# strip .sif file name from the creation command
SIF=$(awk '{split($0, array, " "); print $NF}' <<< "$CREATE")

# create the .sif file
./prepare2submit $CREATE

# Get the command to convert from .result to NetCDF
CONVERT=$( sed -n "${SLURM_ARRAY_TASK_ID}p" <in_fp> )

# Execute the .sif file
ElmerSolver $SIF

# Do post-processing file conversion
$CONVERT

# remove the sif
rm -f $SIF
