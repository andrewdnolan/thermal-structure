#!/bin/bash
#SBATCH --job-name=DaskSingleNode
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G

export NUM_WORKERS=12
export THREADS_PER_WORKER=1
export MEM_PER_WORKER=2GB

# This file (created by the scheduler) communicates Dask network
# information between the scheduler, workers and the client program
export SCHEDULER_FILE=${SLURM_JOB_ID}-scheduler.json

# load the computing environment
source ../../config/modulefile.cc.cedar

# Start a scheduler and wait for it to spin up
dask-scheduler --host 127.0.0.1 \
               --no-dashboard \
               --scheduler-file $SCHEDULER_FILE &
sleep 15

# Start a bunch of workers connected to the scheduler and wait
for worker in `seq $NUM_WORKERS`
do
dask-worker tcp://192.168.1.27:8786 \
             --no-dashboard \
             --no-nanny \
             --nworkers 1 \
             --nthreads $THREADS_PER_WORKER \
             --memory-limit $MEM_PER_WORKER &
done

time python3 ../../src/thermal/grid_data.py "./result/crmpt12/nc/crmpt12_dx_50_NT_30000_dt_0.1_MB_-0.50_OFF_Tma_-8.5_prog.nc" \
                                    -out_fn "./result/crmpt12/gridded/crmpt12_dx_50_NT_30000_dt_0.1_MB_-0.50_OFF_Tma_-8.5_prog.zarr"
