# Solver Timing 

This directory contains a small, self contained, simulation to profile the `CPU` time used by each of the solvers in an initialization simulation.
At all places possible symbolic links to the corresponding files/folders from the `00_CoupledInit` directory are used to minimize the amount of duplicate code. 

We use the "reference glacier" initialization to profile the time used, and there we need to run
```bash
./initialize.py -dx 50 --key "crmpt12"            \\
                -t_f 3000 -dt 0.1 -Dynamic_int 10 \\
                -off -0.37 -T_ma -8.5 
```