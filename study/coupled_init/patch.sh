#!/bin/bash
#SBATCH --array=1-29                               # 442 jobs that run
#SBATCH --job-name=fix_gridding                    # base job name for the array
#SBATCH --mem-per-cpu=2250M                        # maximum 2250MMB per job
#SBATCH --time=0:20:00                             # maximum walltime per job
#SBATCH --nodes=1                                  # Only one node is needed
#SBATCH --ntasks=1                                 # These are serial jobs
#SBATCH --mail-type=ALL                            # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu       # email to spend updates too
#SBATCH --output=logs/fix_gridding_%a.out          # standard output
#SBATCH --error=logs/fix_gridding_%a.err           # standard error
# in the previous two lines %A" is replaced by jobID and "%a" with the array index

#Load cedar module file
source ../../config/modulefile

# Get the file to fix
target=$( sed -n "${SLURM_ARRAY_TASK_ID}p" bad_files )
# Get the original filename
source="${target/_gridded.nc/.nc}"
# get rid of the problem target files
rm $target
# re-grid the NetCDF file written by the NetcdfUGRIDOutputSolver
python3 grid_data.py $source
