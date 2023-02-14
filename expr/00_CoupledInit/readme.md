# Coupled Initialization

Coupled model runs to initialize our synthetic geometries and solve for steady-
state thermal structures.

### Single initialization

The basis of this repository is the `initialize.py` script, which is a wrapper for the `bash` functions in the `initialize.sh` file.
An example usage of script is:
```bash
./initialize.py -dx 50 --key "crmpt12" -t_f 3000 -dt 0.1 -Dynamic_int 10 -off -0.50 -T_ma -9.0
```
where:
  - `-dx 50`
  - `--key "crmpt12"`
  - `t_f 3000`
  - `-dt 0.1`
  - `-Dynamic_int 10`
  - `-off -0.50`
  - `-T_ma -9.0`

This is useful for running one off simulations and doing some preliminary inquiry into parameter space.

### Gridsearch

Assuming the `params/${KEY}.json` file is update with the adequate parameter values
```bash
./make_submission.py gridsearch -key ${KEY}
```

Then on `cedar` run
```bash
sbatch ./run/${key}_gridsearch.sh
```
to conduct the gridsearch over air temperature and mass balance.

### Resubmission

It's inevitable that some of the grid search will fail.
There are various reason for some runs to fail, primarily running out of walltime and under allocation of memory.


Load the initialization functions
```bash
source initialize.sh
```
so they can be run from command line.
Then run,
```bash
find_incomplete '${KEY}'
```
which will make a file `run/${KEY}.incomplete` file with all the air temperature and mass balance parameter combination which did not finish.

Then create the resubmission
```bash
./make_submission.py resubmit -key ${KEY} -mem "" -wt ""
```


If (or when) a subset of model runs fail, there are a series of helper functions to determine what subset of the original runs did not execute or did not finish.
