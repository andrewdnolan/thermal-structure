#!/usr/bin/env python3

import os
import sys
import json
import argparse
import numpy as np

def N_points(j,f):
    # find the number of grid points for a given field (f).
    # Need to add one because it's inclusive of the
    # final point in the grid search
    x = j[f]['range']
    return int((x[2] - x[0]) / x[1]) + 1

def make_script(KEY):
    """Make the SLURM job array submission script.
       All the parameter information is contained in the params/{KEY}.json files

       KEY (str) ---> glacier identifier
    """

    #################################################
    # Get all the parameter values from the json file
    #################################################
    # open the json file
    with open(f'params/{KEY}.json') as file:
        j = json.load(file)

    # integration length [a]
    t_f = j['TT']
    # horizontal mesh resolution [m]
    dx = j['dx']
    # time step [a]
    dt = j['dt']
    # (W)all (T)ime HH:MM:SS
    WT  = j['walltime']
    # Memory request in megabytes (M)
    MEM = j['memory']
    # number of temperature gridpoints
    N_t = N_points(j, 'T_ma')
    # number of mass balance gridpoints
    N_b = N_points(j, "MB")
    # Total number of gridpoints
    N_P = N_t*N_b
    # (J)ob (S)tride
    J_s = j['stride']


    #################################################
    # dump the parameter values in the template
    #################################################
    script = f"""#!/bin/bash
#SBATCH --array=1-{N_P}%{J_s}                           # {N_P} jobs that run
#SBATCH --job-name={KEY}_coupled_init             # base job name for the array
#SBATCH --mem-per-cpu={MEM}                        # maximum 2250MMB per job
#SBATCH --time={WT}                             # maximum walltime per job
#SBATCH --nodes=1                                  # Only one node is needed
#SBATCH --ntasks=1                                 # These are serial jobs
#SBATCH --mail-type=ALL                            # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu       # email to spend updates too
#SBATCH --output=logs/{KEY}/coupled_init_%A_%a.out  # standard output
#SBATCH --error=logs/{KEY}/coupled_init_%A_%a.err   # standard error
# in the previous two lines %A" is replaced by jobID and "%a" with the array index

#Load cedar module file
source ../../config/modulefile.cc.cedar

# load the souce functions
source initialize.sh

# parse the parameters from the json files
parse_json "params/{KEY}.json"

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
./initialize.py -dx {dx} --key "{KEY}" -t_f {t_f} -off $offset -T_ma $T_ma
"""


    ##################################
    # Write the bash submission script
    ##################################
    with open(f'run/{KEY}_gridsearch.sh', 'w') as f:
        f.write(script)


def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('key',
                        help='Glacier identifier to create submission scipt for.')

    args, _ = parser.parse_known_args(argv)

    key = args.key

    if key.endswith('.json'):
        # if the json was passed instead of the key, then just get the file name
        file_name = os.path.basename(key)
        # .split is more cautious then .strip incase key resembels json
        key = file_name.split('.json')[0]


    #with the key, make the submission script
    make_script(key)

if __name__ == '__main__':
    main(sys.argv[1:])
