#!/usr/bin/env bash

parse_reference()
{
  #-----------------------------------------------------------------------------
  # Parse reference model parameter from the json file
  # perl code from: https://stackoverflow.com/a/27127626/10221482
  #-----------------------------------------------------------------------------
  # $1 (var) --->   string 
  #-----------------------------------------------------------------------------
  export var=$1
  # parse MB grid search start
  reference=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{params}->{$ENV{var}}->{reference}
    ' params/crmpt12.json )
}


parse_params()
{
  IFS=" " read -r C_firn f_dd w_en w_aq IC \
    <<< $(sed -n "${1}p" ./run/crmpt12/gridsearch.commands | cut -d " " -f 17,19,21,23,25)
}

create_dask_cluster()
{
  export SCHEDULER_FILE="${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}-scheduler.json"
  dask scheduler --host 127.0.0.1 --no-dashboard --scheduler-file $SCHEDULER_FILE &
  sleep 15

  for worker in $(seq $NUM_WORKERS); do
  dask worker --scheduler-file $SCHEDULER_FILE \
              --no-dashboard \
              --no-nanny \
              --nworkers 1 \
              --nthreads 1 &
  done
  sleep 15

}

post_proccess()
{
  # get the j-th air temp and mass balance values
  parse_params $1

  # based on the set parameter values, create a unique run name 
  run_name=$(make_runname)

  # create parameter json (dictionary) for gridding Zarr files
  param_dict="{\"param\"   : ${param},
              \"${param}\" : ${value}}"

  # copy the source file from scratch to local (compute node's) SSD
  time rsync -ah "result/${KEY}/nc/${run_name}.nc" "${SLURM_TMPDIR}"

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
  time tar -cf "result/${KEY}/gridded/${run_name}.zarr.tar" -C "${SLURM_TMPDIR}" "${run_name}.zarr"

  # tar the thinned zarr file, and write the tar to scratch
  time tar -cf "result/${KEY}/thinned/${run_name}.zarr.tar" -C "${SLURM_TMPDIR}/thinned" "${run_name}.zarr"

  # delete files from SSD to make room for next files
  rm "${SLURM_TMPDIR}/${run_name}.nc"
  rm -r "${SLURM_TMPDIR}/${run_name}.zarr"
  rm -r "${SLURM_TMPDIR}/thinned/${run_name}.zarr"
}

calc_enthalpy()
{ #-----------------------------------------------------------------------------
  # Convert from temperature [C] to enthalpy [J kg-1]
  #
  # Variables:
  # ----------
  # $1  ---> temperature [C]
  #-----------------------------------------------------------------------------

  # simple python script to convert from temp to enthalpy
  script="from thermal.utils import calc_enthalpy; print(calc_enthalpy(273.15 + float(${1})))"

  # run script and print result
  python -c "${script}"
}

prognostic_run()
{
  #-----------------------------------------------------------------------------
  # Run a prognostic simulation. This function:
  #         1) upates the template sif file,
  #         2) runs the prognostic Elmer/Ice simulation,
  #         3) adds the .sif file and params as attributes to the NetCDF
  #
  # Variables:
  # ---------
  #  KEY       ---> glacier identifer
  #  dx        ---> mesh resolution. ${KEY}/mesh_dx${dx}/mesh.* should exists
  #  T_ma      ---> mean annual air temp. [C] @ $z_ref from "params/ref_params.sif" file
  #  offset    ---> Mass balance anomoly [m i.e. yr-1]
  #  SS_itters ---> Number of S.S. itterations for diagnostic simulation
  #  run_name  ---> unique simulation identifer
  #
  #  dt        ---> timestep             [yr]
  #  NT        ---> Number of timesteps  [-]
  #
  #-----------------------------------------------------------------------------

  # Update the .SIF FILE with the model run specifc params
  sed "s#<C_firn>#"$C_firn"#g;
       s#<w_en>#"$w_en"#g;
       s#<w_aq>#"$w_aq"#g;
       s#<f_dd>#"$f_dd"#g;
       s#<enth_IC>#"$enth_IC"#g;
       s#<DX>#"$dx"#g;
       s#<dt>#"$dt"#g;
       s#<NT>#"$NT"#g;
       s#<KEY>#"$KEY"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<run_name>#"$run_name"#g;
       s#<Dynamic_int>#"$Dynamic_int"#g
       s#<SS_itters>#"$SS_itters"#g;" "./sifs/prognostic.sif" > "./sifs/${run_name}.sif"

  # Run the model
  ElmerSolver "./sifs/${run_name}.sif" | tee $log_file

  # add the sif as a global attribute to the netcdf file
  python3 ../../src/thermal/add_attr.py -f "./sifs/${run_name}.sif" \
                                        -a "sif"                \
                                           "result/${KEY}/nc/${run_name}.nc"
  
  # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

make_runname()
{
  # NOTE, these are NOT the evauated variables but their string "names"
  for var in  C_firn f_dd w_en w_aq IC; do 

    # get the reference value for current parameter
    parse_reference $var
    # store current parameter value in dummy varible
    eval test_value='$'$var
    
    # compare the current parameter value to the reference
    if awk "BEGIN {exit !($reference != $test_value)}"; then

      # if we need to do simultaneous parameter variation, this should be a list 
      # that we append varibale names and values to 
      param=$var
      value=$test_value
    fi
  done

  echo "${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_1aTST_${param}_${value}"

  #note: param and value are also exported
}

log_runtime()
{

  local OUT_fp="result/${KEY}/${KEY}.sensitivity.time_profile"

  if [ ! -f "$OUT_fp" ]; then
      echo "#run_name runtime" |
      awk -v OFS='\t' '{print $1 "\t" $2}' >> \
                      $OUT_fp
  fi

  echo "${run_name} ${runtime}" |
  awk -v OFS='\t'  '{print $1 "\t" $2}' >> \
                  $OUT_fp
}


test_sensitivity()
{ # Variables that need to be set:
  #    Dynamic_int
  #    SS_itters
  #    KEY
  #    dx
  #    offset
  #    T_ma
  #    dt
  #    t_f
  #-----------------------------------------------------------------------------
  # NOTE: these variables are most likely set wihtin the python wrapper. 
  #-----------------------------------------------------------------------------


  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  #  Prognostic parameteric sensitivity tests
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  # Number of time interval based on dt
  NT=$(awk -v dt=$dt -v t_f=$t_f 'BEGIN {OFMT = "%.0f"; print (t_f/dt)}')

  # convert IC from temp [C] to enthalpy [J kg-1]
  enth_IC=$(calc_enthalpy $IC)

  # based on the set parameter values, create a unique run name 
  run_name=$(make_runname)

  # Start the timer
  start=$(date +%s.%N)

  # fill in the template sif file, and run the model
  prognostic_run

  # End the timer
  end=$(date +%s.%N)

  # Execution time of the solver
  runtime=$(awk -v start=$start -v end=$end 'BEGIN {print end - start}')

  # log the runtime info 
  log_runtime $run_name $runtime
}
