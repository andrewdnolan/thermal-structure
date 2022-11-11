#!/usr/bin/env python3

from itertools import product

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
cmd = "./periodic_surge.py -k \"{KEY}\" -SP {SP} -QP {QP} -beta {beta} -TT {TT}"

# open the commands file used by the job array
with open(f'./run/{key}/{key}.commands', 'w') as f:

    # itterate over active phase and surge-cycle periods [yr]
    for i, (SP, CP) in enumerate(product([2,5], [15, 30, 60])):

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
