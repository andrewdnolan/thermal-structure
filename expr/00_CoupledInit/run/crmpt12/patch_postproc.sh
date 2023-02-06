#!/bin/bash
#SBATCH --job-name=DaskSingleNode
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000M
#SBATCH --output=16_cores_2GB_JID_%A.out   # standard output
#SBATCH --error=16_cores_2GB_JID_%A.err    # standard error

export NUM_WORKERS=16
export THREADS_PER_WORKER=1


parse_params()
{
  IFS=" " read -r offset T_ma \
    <<< $(sed -n "${1}p" ./run/crmpt12/gridsearch.commands | cut -d " " -f 13,15)
}

source ../../config/modulefile.cc.cedar

export SCHEDULER_FILE=${SLURM_JOB_ID}-scheduler.json
dask-scheduler --host 127.0.0.1 --no-dashboard --scheduler-file $SCHEDULER_FILE &
sleep 15

for worker in $(seq $NUM_WORKERS); do
dask worker --scheduler-file $SCHEDULER_FILE \
            --no-dashboard \
            --no-nanny \
            --nworkers 1 \
            --nthreads 1 &
done
sleep 15


CHUNK_SIZE=48

# # start the dask client for each job in the array

# for SLURM_ID in $(seq 0 5); do 
#   start=$((SLURM_ID * CHUNK_SIZE))
#   stop=$((start + CHUNK_SIZE))
#   for j in $(seq $start $stop); do 
#     if [[ $j -le 286 ]]; then 
#       echo $j
#     fi
#   done 
# done


mkdir "${SLURM_TMPDIR}/thinned"

if [ -d "result/crmpt12/thinned/" ]; then
  mkdir  "result/crmpt12/thinned/"
fi

# What need to get done for each file: 
# 
#   1. rsync file to $SLURM_TMPDIR, i.e. on the local (fast) SSD
#   2. grid_data.py file, writing zarr to $SLURM_TMPDIR
#   3. downsample.py file, again writing to $SLURM_TMPDIR
#   4. tar whole zarr file
#   5. tar downsampled zarr file
#   6. mv whole zarr file to scratch
#   7. mv downsampled zarr file to 
#   8. delete any local copies in  $SLURM_TMPDIR to make room for next run. 

# get the j-th air temp and mass balance values
parse_params $CHUNK_SIZE

# from the parameters fill out the run_name
run_name="crmpt12_dx_50_NT_30000_dt_0.1_MB_${offset}_OFF_Tma_${T_ma}_prog"

# create parameter json (dictionary) for gridding Zarr files
param_dict="{\"T_ma\"   : ${T_ma},
             \"offset\" : ${offset}}"

# copy the source file from scratch to local (compute node's) SSD
time rsync -ah "result/crmpt12/nc/${run_name}.nc" "${SLURM_TMPDIR}"

# grid the NetCDF file written by the NetcdfUGRIDOutputSolver, 
# convert from NetCDF to Zarr file format
time grid_data.py -i "${SLURM_TMPDIR}/${run_name}.nc" \
             -o "${SLURM_TMPDIR}/${run_name}.zarr" \
             -p "${param_dict}"

# run the subsampling script, write a years worth of data every 10 years
time downsample.py -i "${SLURM_TMPDIR}/${run_name}.zarr" \
              -o  "${SLURM_TMPDIR}/thinned/${run_name}.zarr" \
              --value --years_worth 10

# tar the full zarr file, and write the tar to scratch
time tar -cf "result/crmpt12/gridded/${run_name}.zarr.tar" -C "${SLURM_TMPDIR}" "${run_name}.zarr"

# tar the thinned zarr file, and write the tar to scratch
time tar -cf "result/crmpt12/thinned/${run_name}.zarr.tar" -C "${SLURM_TMPDIR}/thinned" "${run_name}.zarr"

# delete files from SSD to make room for next files
rm "${SLURM_TMPDIR}/${run_name}.nc"
rm "${SLURM_TMPDIR}/${run_name}.zarr"
rm "${SLURM_TMPDIR}/thinned/${run_name}.zarr"

# echo $run_name
