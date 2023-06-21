#!/bin/bash
#SBATCH --job-name=netcdf_io_error                 # base job name for the array
#SBATCH --mem-per-cpu=4000M                        # maximum 200M per job
#SBATCH --time=4-4:00:00                                # maximum walltime per job
#SBATCH --nodes=1                                  # Only one node is needed
#SBATCH --ntasks=1                                 # These are serial jobs
#SBATCH --mail-type=ALL                            # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu       # email to spend updates too
#SBATCH --output=logs/crmpt12/netcdf_io_error_%A.out # standard output
#SBATCH --error=logs/crmpt12/netcdf_io_error_%A.err  # standard error
# in the previous two lines %A" is replaced by jobID and "%a" with the array index


# load the neccessary packages
source ../../config/modulefile.cc.cedar

./periodic_surge.py -k "crmpt12" -SP 2 -QP 58 -beta 1.778e-04 -TT 3000 -T0 0

# Exit status of the code:
STATUS=$?

# save the job statuses to file with unique filename per submission
# echo $SLURM_ARRAY_TASK_ID $STATUS >> ./run/crmpt12/crmpt12_$SLURM_ARRAY_JOB_ID.status
