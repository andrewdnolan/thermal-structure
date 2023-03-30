#!/usr/bin/env python3

import numpy as np
# from itertools import product

# glacier identifer
key  = 'crmpt12'
# simulation length [yr]
TT   = 3000
# slip coeffiecent  [???]
beta = 0.001
# walltime in hh:mm:ss
WT   = '75:00:00'
# memory allocation
MEM  = '4G'

# base command used to execute a single simulation
cmd = "./surge2steady.py -k \"{KEY}\" -SP 2 -QP 3000 -beta {beta}"

# open the commands file used by the job array
with open(f'./run/{key}/{key}.commands', 'w') as f:

    # itterate various beta values
    for i, beta in enumerate(np.logspace(-6, -8, 9)**0.5):

        # dump the command in the text file
        f.write(cmd.format(KEY=key, beta=beta))
        f.write('\n')
# account for pythons zero indexing
i += 1

# read the submission template
with open(f'./run/submission.template', 'r') as f:
    template = f.read()

# fill the submission template with the set parameters and write to disk
with open(f'./run/{key}/{key}_submit.sh', 'w') as f:
    f.write(template.format(KEY=key, M=i, S=i, run_name=key, WT=WT, MEM=MEM))
