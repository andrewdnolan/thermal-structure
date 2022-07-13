#!/usr/bin/env bash
#SBATCH --job-name=subset                          # base job name for the array
#SBATCH --mem-per-cpu=20000M                       # maximum 2250MMB per job
#SBATCH --time=2:00:00                             # maximum walltime per job
#SBATCH --nodes=1                                  # Only one node is needed
#SBATCH --ntasks=1                                 # These are serial jobs
#SBATCH --mail-type=ALL                            # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu       # email to spend updates too
#SBATCH --output=logs/post_proc.out                # standard output
#SBATCH --error=logs/post_proc.err                 # standard error
# in the previous two lines %A" is replaced by jobID and "%a" with the array index

#Load cedar module file
source ../../config/modulefile.cc.cedar

python3 subset_test.py
