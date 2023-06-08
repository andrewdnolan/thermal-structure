#!/bin/bash
#SBATCH --job-name=dask_gridding
#SBATCH --time=02:00:00           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000M
#SBATCH --mail-type=ALL                      # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu # email to spend updates too
#SBATCH --output=dask_%A_%a.out              # standard output
#SBATCH --error=dask_%A_%a.err               # standard error


# numbers of cores each job in the array will have
export NUM_WORKERS=16
# use a single thread per cpu core
export THREADS_PER_WORKER=1

post_proccess()
{
  T_ma=-8.5
  offset=$1

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
  rm -r "${SLURM_TMPDIR}/${run_name}.zarr"
  rm -r "${SLURM_TMPDIR}/thinned/${run_name}.zarr"
}


# make the thinned dir for each task in the job array
mkdir "${SLURM_TMPDIR}/thinned"

for off in $(seq -w -0.65 0.01 -0.51) -0.31 -0.36 -0.27; do 
    # run the post processing commands
    post_proccess $off
done