#!/bin/bash
#SBATCH --array=1-9%9
#SBATCH --job-name=crmpt12_%A                      # base job name for the array
#SBATCH --mem-per-cpu=4000M                        # maximum 200M per job
#SBATCH --time=48:00:00                                # maximum walltime per job
#SBATCH --nodes=1                                  # Only one node is needed
#SBATCH --ntasks=1                                 # These are serial jobs
#SBATCH --mail-type=END                            # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu       # email to spend updates too
#SBATCH --output=logs/crmpt12/crmpt12_%A_%a.out   # standard output
#SBATCH --error=logs/crmpt12/crmpt12_%A_%a.err    # standard error
# in the previous two lines %A" is replaced by jobID and "%a" with the array index

echo "Starting task $SLURM_ARRAY_TASK_ID"

# load the neccessary packages
source ../../config/modulefile.cc.cedar

CMD=$(  sed -n "${SLURM_ARRAY_TASK_ID}p" ./run/crmpt12/crmpt12.commands)

# excute the ith model run
eval "$CMD"

# Exit status of the code:
STATUS=$?

# save the job statuses to file with unique filename per submission
echo $SLURM_ARRAY_TASK_ID $STATUS >> ./run/crmpt12/crmpt12_$SLURM_ARRAY_JOB_ID.status
