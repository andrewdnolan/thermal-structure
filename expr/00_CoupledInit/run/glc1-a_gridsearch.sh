#!/bin/bash
#SBATCH --array=1-441%40                           # 441 jobs that run
#SBATCH --job-name=glc1-a_coupled_init             # base job name for the array
#SBATCH --mem-per-cpu=2250M                        # maximum 2250MMB per job
#SBATCH --time=6:30:00                             # maximum walltime per job
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

# get ith offset and T_ma corresponing to SLURM_ARRAY_TASK_ID.
# Nested for loops each job in the job-array is junky.
# The cartesian product of the two arrays would be a more elegant solution,
# but I was having trouble gettig that to work.
# Refs:
#   - https://unix.stackexchange.com/questions/97814/array-cartesian-product-in-bash
#   - https://stackoverflow.com/questions/23363003/how-to-produce-cartesian-product-in-bash
#   - https://rosettacode.org/wiki/Cartesian_product_of_two_or_more_lists#UNIX_Shell

count=1
# loop over the mass balance offsets
for offset in $(seq -w $MB_0 $MB_s $MB_f); do
  # loop over mean annual air temps
  for T_ma in $(seq -w $T_ma_0 $T_ma_s $T_ma_f); do
    # check if counter equals SLURM_ARRAY_TASK_ID
    if [[ $count -eq $SLURM_ARRAY_TASK_ID ]]; then
        break 2
    fi
    count=$((count+1))
  done
done

# run the initialization sequence for one air temp and mass balance
# combination.
./initialize.py -dx 50 --key "glc1-a" -t_f 2000 -off $offset -T_ma $T_ma
