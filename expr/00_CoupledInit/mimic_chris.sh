#!/bin/bash
#SBATCH --job-name=DaskSingleNode
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000M
#SBATCH --output=mimic_chris_%A.out   # standard output
#SBATCH --error=mimic_chris_%A.err    # standard error


module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 python/3.8.10 scipy-stack

# Create virtual env on fast local disk of node
virtualenv $SLURM_TMPDIR/venv
source $SLURM_TMPDIR/venv/bin/activate
pip install --no-index --upgrade pip
pip install --no-index dask distributed xarray

export SCHEDULER_FILE=${SLURM_JOB_ID}-scheduler.json
dask-scheduler --host 127.0.0.1 --no-dashboard --scheduler-file $SCHEDULER_FILE &
sleep 15
dask-worker --scheduler-file $SCHEDULER_FILE --no-dashboard --no-nanny --nworkers 1 --nthreads 1 &
dask-worker --scheduler-file $SCHEDULER_FILE --no-dashboard --no-nanny --nworkers 1 --nthreads 1 &
dask-worker --scheduler-file $SCHEDULER_FILE --no-dashboard --no-nanny --nworkers 1 --nthreads 1 &
dask-worker --scheduler-file $SCHEDULER_FILE --no-dashboard --no-nanny --nworkers 1 --nthreads 1 &

# Bogus input/output
touch in
touch out

PYTHONPATH=$PYTHONPATH:../../src/thermal python3 "../../src/thermal/grid_data.py "./result/crmpt12/nc/crmpt12_dx_50_NT_30000_dt_0.1_MB_-0.50_OFF_Tma_-8.5_prog.nc" \
                                    -out_fn "./result/crmpt12/gridded/crmpt12_dx_50_NT_30000_dt_0.1_MB_-0.50_OFF_Tma_-8.5_prog.zarr"
