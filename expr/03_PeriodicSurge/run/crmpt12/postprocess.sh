#!/bin/bash
#SBATCH --array=0-9
#SBATCH --job-name=dask_gridding
#SBATCH --time=01:30:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4096M
#SBATCH --output=dask_%A.out
#SBATCH --error=dask_%A.err

export KEY='crmpt12'

# number of files each member of the job array will process
CHUNK_SIZE=11

# because mostly numeric dominated (e.g. numpy) use a large number of threads
export THREADS_PER_WORKER=16

source ../../config/modulefile.cc.cedar

parse_params()
{   
  IFS=" " read -r QP beta NT T_0<<< $(sed -n "${1}p" "run/${KEY}/${KEY}.commands" | cut -d " " -f 7,9,11,13)
}

post_proccess()
{
  parse_params $1

  # get the run name based on parsed parameters for the i-th run
  # NOTE: default value which aren't varied are hard coded
  run_name="${KEY}_dx_50_TT_${NT}.0_MB_-0.36_OFF_Tma_-8.5_B_${beta}_SP_2_QP_${QP}"

  # based on the start time (T_0) and time intergration length (NT)
  # caluculate the final time in kiloyears
  T_f=$(awk -v T_0=$T_0 -v NT=$NT 'BEGIN {print (T_0 + NT)/1e3}')
  # convert the start time from years to kiloyears
  T_0=$(awk -v T_0=$T_0 'BEGIN {print T_0/1e3}')

  # rename the file based on the restart point and integration length
  new_name="${KEY}_dx_50_TT_${T_0}--${T_f}ka_MB_-0.36_OFF_Tma_-8.5_B_${beta}_SP_2_QP_${QP}"

  # rename restart files
  mv "result/${KEY}/mesh_dx50/${run_name}.result"\
     "result/${KEY}/mesh_dx50/${new_name}.result"

  # rename the "raw" netcdf files
  mv "result/${KEY}/nc/${run_name}.nc"\
     "result/${KEY}/nc/${new_name}.nc"

  # overwrite the runname with the updated name,
  # we postprocess after renaming so the files packed into the zarr archieves have 
  # the correct filenames when unpacked
  run_name="${new_name}"

  # copy the source file from scratch to local (compute node's) SSD
  time rsync -ah "result/${KEY}/nc/${run_name}.nc" "${SLURM_TMPDIR}"

  # grid the NetCDF file written by the NetcdfUGRIDOutputSolver, 
  # convert from NetCDF to Zarr file format
  time grid_data.py -i "${SLURM_TMPDIR}/${run_name}.nc" \
                    -o "${SLURM_TMPDIR}/${run_name}.zarr"  #\
                    # -p "${param_dict}"

  # tar the full zarr file, and write the tar to scratch
  time tar -cf "result/${KEY}/gridded/${run_name}.zarr.tar" -C "${SLURM_TMPDIR}" "${run_name}.zarr"

  # delete files from SSD to make room for next files
  rm "${SLURM_TMPDIR}/${run_name}.nc"
  rm -r "${SLURM_TMPDIR}/${run_name}.zarr"

}

# get the number of simulations to process
N=$(wc -l < "run/${KEY}/${KEY}.commands")

# loop of the corresponding input files
start=$((SLURM_ARRAY_TASK_ID * CHUNK_SIZE))
stop=$((start + CHUNK_SIZE))
for j in $(seq $start $stop); do 
  # add one to j since sed is indexed at 1
  j=$((j + 1))
  # only ran 286 runs, so make sure we don't go over 
  if [[ $j -le $N ]]; then 
    # run the post processing commands
    post_proccess $j
  fi
done
