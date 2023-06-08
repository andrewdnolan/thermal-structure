#!/usr/bin/env python3

import numpy as np
from itertools import product

# glacier identifer
key  = 'crmpt12'
# simulation length [yr]
TT   = 3000
# walltime in dd-hh:mm:ss
WT   = '2-00:00:00'
# memory allocation
MEM  = '4000M'

sim_name= "V_versus_FT"

offsets = np.linspace(-0.65, -0.51, 15)

# base command used to execute a single simulation
cmd = "./initialize.py -dx 50 -k \"{KEY}\" -t_f {TT} -dt 0.1 -Dynamic_int 10 -off {offset:1.2f} -T_ma -8.5 "

# open the commands file used by the job array
with open(f'./run/{key}/{sim_name}.commands', 'w') as f:

    # itterate over betas and surge-cycle periods [yr]
    for i, offset in enumerate(offsets):

        # dump the command in the text file
        f.write(cmd.format(KEY=key, TT=TT, offset=offset))
        f.write('\n')

# account for pythons zero indexing
i += 1

# read the submission template
with open(f'./run/submission.template', 'r') as f:
    template = f.read()

# fill the submission template with the set parameters and write to disk
with open(f'./run/{key}/{sim_name}_submit.sh', 'w') as f:
    f.write(template.format(KEY=key, M=i, S=i, search_type=sim_name, run_name=key, WT=WT, MEM=MEM))