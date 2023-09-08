#!/bin/bash
#SBATCH --job-name=Adrien_Edits                      # base job name for the array
#SBATCH --mem-per-cpu=4000M                        # maximum 200M per job
#SBATCH --time=15-4:00:00                          # maximum walltime per job
#SBATCH --nodes=1                                  # Only one node is needed
#SBATCH --ntasks=1                                 # These are serial jobs
#SBATCH --mail-type=ALL                            # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu       # email to spend updates too
#SBATCH --output=logs/crmpt12/Adrien_Edits.out   # standard output
#SBATCH --error=logs/crmpt12/Adrien_Edits.err    # standard error
# in the previous two lines %A" is replaced by jobID and "%a" with the array index

# load the neccessary packages
source ../../config/modulefile.cc.cedar

./periodic_surge.py -k "crmpt12" -SP 2 -QP 28 -beta 1.000e-04 -TT 9000 -T0 0 
