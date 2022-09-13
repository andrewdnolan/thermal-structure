#!/bin/bash
#SBATCH --array=1-82%40                           # 82 jobs that run
#SBATCH --job-name=glc1-a_coupled_init             # base job name for the array
#SBATCH --mem-per-cpu=3500M                        # maximum 2250MMB per job
#SBATCH --time=12:00:00                             # maximum walltime per job
#SBATCH --nodes=1                                  # Only one node is needed
#SBATCH --ntasks=1                                 # These are serial jobs
#SBATCH --mail-type=ALL                            # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu       # email to spend updates too
#SBATCH --output=logs/glc1-a/coupled_init_%A_%a.out  # standard output
#SBATCH --error=logs/glc1-a/coupled_init_%A_%a.err   # standard error
# in the previous two lines %A" is replaced by jobID and "%a" with the array index

#Load cedar module file
source ../../config/modulefile.cc.cedar

# load the souce functions
source initialize.sh

# parse the parameters from the json files
parse_json "params/glc1-a.json"

# Add one to SLURM_ARRAY_TASK_ID to slip header
count=$((SLURM_ARRAY_TASK_ID + 1))
# get the mean annual airtemp and mass balance offset
T_ma=$(sed -n "${count}p" "./run/glc1-a.incomplete" | cut -d$'	' -f 1)
offset=$(sed -n "${count}p" "./run/glc1-a.incomplete" | cut -d$'	' -f 2)

# at the end of the array, remove the *.incomplete file
if [[ $SLURM_ARRAY_TASK_ID == 82 ]]; then
    rm "./run/glc1-a.incomplete"
fi

# run the initialization sequence for one air temp and mass balance
# combination.
./initialize.py -dx 50 --key "glc1-a" -t_f 2000 -off $offset -T_ma $T_ma
