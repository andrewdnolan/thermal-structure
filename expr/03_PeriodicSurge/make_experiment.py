#!/usr/bin/env python3

import numpy as np
from itertools import product

# glacier identifer
key  = 'crmpt12'
# simulation length [yr]
TT   = 3000
# walltime in dd-hh:mm:ss
WT   = '6-12:00:00'
# memory allocation
MEM  = '4000M'
# fixed 2-year surge period [a]
SP = 2

# ##################################################
# # Reference Experiment:
# ##################################################
# slip coefficients to test (3 order of magnitude)
betas  = np.logspace(-3, -4, 9)
# length of the surge cycle [a]
cycles = [15, 30, 60]

# ##################################################
# # Bisect to see periodicity emerge Experiment:
# ##################################################
# # slip coefficients to test (3 order of magnitude)
# betas  = np.geomspace(0.0004217,0.00031623, 7)[1:-1]
# # length of the surge cycle [a]
# cycles = [30]

# base command used to execute a single simulation
cmd = "./periodic_surge.py -k \"{KEY}\" -SP {SP} -QP {QP} -beta {beta:1.3e} -TT {TT}"

# open the commands file used by the job array
with open(f'./run/{key}/{key}.commands', 'w') as f:

    # itterate over betas and surge-cycle periods [yr]
    for i, (beta, CP) in enumerate(product(betas, cycles)):

        # calculate quiescent period [yr]
        QP = CP-SP

        # dump the command in the text file
        f.write(cmd.format(KEY=key, SP=SP, QP=QP, TT=TT, beta=beta))
        f.write('\n')

# account for pythons zero indexing
i += 1

# read the submission template
with open(f'./run/submission.template', 'r') as f:
    template = f.read()

# fill the submission template with the set parameters and write to disk
with open(f'./run/{key}/{key}_submit.sh', 'w') as f:
    f.write(template.format(KEY=key, M=i, S=i, run_name=key, WT=WT, MEM=MEM))
