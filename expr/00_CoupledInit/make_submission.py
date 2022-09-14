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

def make_resubmission_script(KEY, WT=None, MEM=None):
    """ After parsing through which model runs failed this function makes a
        script to resubmit them.
    """

    # make sure the "./run/${KEY}.incomplete" file exists
    if not os.path.exists(f"run/{KEY}.incomplete"):
        raise FileNotFoundError(f"No './run/{KEY}.incomplete' file exists")

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
    # number of time integration points
    NT = int(t_f / dt)
    # Total number of gridpoints
    N_P = sum(1 for _ in open(f"run/{KEY}.incomplete")) - 1
    # (J)ob (S)tride
    J_s = j['stride']

    if WT == None:
        # (W)all (T)ime HH:MM:SS
        WT  = j['walltime']

        print(f'''WARNING: Using walltime from "params/${KEY}.json".
                           this has previsouly cause runs to fail. Have you updated "params/${KEY}.json"?? ''')
    if MEM == None:
        # Memory request in megabytes (M)
        MEM = j['memory']

        print(f'''WARNING: Using memory from "params/{KEY}.json".
                           this has previsouly cause runs to fail. Have you updated "params/${KEY}.json"?? ''')

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

# Add one to SLURM_ARRAY_TASK_ID to slip header
count=$((SLURM_ARRAY_TASK_ID + 1))
# get the mean annual airtemp and mass balance offset
T_ma=$(sed -n "${{count}}p" "./run/{KEY}.incomplete" | cut -d$'\\t' -f 1)
offset=$(sed -n "${{count}}p" "./run/{KEY}.incomplete" | cut -d$'\\t' -f 2)

# remove the existing *.nc/*.result files since they will cause writing errors
# diagnostic
rm "result/{KEY}/*/{KEY}_dx_{dx}_MB_${{offset}}_OFF_Tma_${{T_ma}}_diag.nc"
rm "result/{KEY}/mesh_dx{dx}/{KEY}_dx_{dx}_MB_${{offset}}_OFF_Tma_${{T_ma}}_diag.result"
# prognostic
rm "result/{KEY}/*/{KEY}_dx_{dx}_NT_{NT}_dt_{dt}_MB_${{offset}}_OFF_Tma_${{T_ma}}_prog.nc"
rm "result/{KEY}/mesh_dx{dx}/{KEY}_dx_{dx}_NT_{NT}_dt_{dt}_MB_${{offset}}_OFF_Tma_${{T_ma}}_prog.result"



# at the end of the array, remove the *.incomplete file
if [[ $SLURM_ARRAY_TASK_ID == {N_P} ]]; then
    rm "./run/{KEY}.incomplete"
fi

# run the initialization sequence for one air temp and mass balance
# combination.
./initialize.py -dx {dx} --key "{KEY}" -t_f {t_f} -off $offset -T_ma $T_ma
"""


    ##################################
    # Write the bash submission script
    ##################################
    with open(f'run/{KEY}_resubmit.sh', 'w') as f:
        f.write(script)


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


def parse_resubmission(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('-key', help='Glacier identifier to create submission scipt for')
    parser.add_argument('-mem', '--memory', help='Amount of Memory (MB) to request e.g. "3000M"')
    parser.add_argument('-WT', '--walltime', help='walltime to request from SLURM e.g. "12:00:00"')

    args, _ = parser.parse_known_args(argv)

    key = args.key
    mem = args.memory
    WT  = args.walltime

    if key.endswith('.incomplete'):
        # if the json was passed instead of the key, then just get the file name
        file_name = os.path.basename(key)
        # .split is more cautious then .strip incase key resembels json
        key = file_name.split('.incomplete')[0]

    #with the key, make the resubmission script
    make_resubmission_script(key, MEM=mem, WT=WT)


def parse_gridsearch(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('-key', help='Glacier identifier to create submission scipt for')

    args, _ = parser.parse_known_args(argv)

    key = args.key

    if key.endswith('.json'):
        # if the json was passed instead of the key, then just get the file name
        file_name = os.path.basename(key)
        # .split is more cautious then .strip incase key resembels json
        key = file_name.split('.json')[0]


    #with the key, make the submission script
    make_script(key)


def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('run_type', help='Type of submission script to make', nargs='?', choices=('gridsearch', 'resubmission'))
    # TO DO: fix help string
    # https://docs.python.org/3.2/library/argparse.html#sub-commands
    args, _ = parser.parse_known_args(argv)

    run_type = args.run_type

    if run_type == "gridsearch":
        parse_gridsearch(argv)
    elif run_type == "resubmission":
        parse_resubmission(argv)
    else:
        raise NotImplementedError('Unaccepted run_type')

if __name__ == '__main__':
    main(sys.argv[1:])
