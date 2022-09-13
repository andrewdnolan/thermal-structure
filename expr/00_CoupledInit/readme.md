# Coupled Initialization

Coupled model runs to initialize our synthetic geometries and solve for steady-
state thermal structures.

__Before you start__:   

The `*.py` scripts should already be executable, but you should double check.
If not make them executable with the command:
```bash
chmod +x ./initialize.py
chmod +x ./make_submission.py
```

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
