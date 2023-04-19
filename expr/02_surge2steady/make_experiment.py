#!/usr/bin/env python3

import numpy as np
# from itertools import product

# glacier identifer
key  = 'crmpt12'
# walltime in hh:mm:ss
WT   = '48:00:00'
# memory allocation
MEM  = '4000M'

# base command used to execute a single simulation
cmd = "./surge2steady.py -k \"{KEY}\" -SP 2 -QP 4000 -beta {beta:1.3e}"\
      "  -cycle2 -RESTART crmpt12_dx_50_NT_40_dt_0.05_MB_-0.37_OFF_Tma_-8.5_B_{beta:1.3e}_pseudo_NT_2000_recovery.result"

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
