#!/usr/bin/env python3

import subprocess

# for i, key in enumerate(['glc1-a', 'lilk-a', 'twds-a']):
#
#     bfp   = './input_data/{}_bed.dat'.format(key)
#     cmd   = "cd  ../.. &&"
#     cmd  += " ./scripts/make_mesh.sh 100"
#     cmd  += " ./study/MB_tune/result/{} {}".format(key,bfp)
#
#     print(cmd)
#     result = subprocess.run(
#                             [cmd],
#                             shell=True,
#                             capture_output=True,
#                             text=True
#                             )
#     #print("stdout:", result.stdout)
#     print("stderr:", result.stderr)

for dx in [100]:
    key   = 'crmpt12'
    bfp   = './input_data/{}_bed.dat'.format(key)
    cmd   = "cd  .. &&"
    cmd  += " ./scripts/make_mesh.sh {}".format(dx)
    cmd  += " ./initialization/coarse/{}/mesh_dx{} {}".format(key,dx,bfp)

    print(cmd)
    result = subprocess.run(
                            [cmd],
                            shell=True,
                            capture_output=True,
                            text=True
                            )
    #print("stdout:", result.stdout)
    print("stderr:", result.stderr)
