#!/usr/bin/env bash


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TO DO:
#    [ ] set global flag for wether to write log files. No need on westgrid becuase
#        of the *.out and *.err files
#
#    [ ] set optional flag for executing the Elmer commands from the docker container
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# source ../../src/utils/elmer_log.sh
source ../../src/utils/initialization.sh

parse_json()
{
  #-----------------------------------------------------------------------------
  # Parse model parameters from input .json file needed for initialization
  # perl code from: https://stackoverflow.com/a/27127626/10221482
  #-----------------------------------------------------------------------------
  # $1 (json_fp) --->   file path to json parameter file
  #-----------------------------------------------------------------------------

  # parse MB grid search start
  MB_0=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{MB}->{range}->[0]
    ' $1 )
  # parse MB grid search stride
  MB_s=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{MB}->{range}->[1]
    ' $1 )
  # parse MB grid search end
  MB_f=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{MB}->{range}->[2]
    ' $1 )
  # parse T_ma grid search start
  T_ma_0=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{T_ma}->{range}->[0]
    ' $1 )
  # parse T_ma grid search stride
  T_ma_s=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{T_ma}->{range}->[1]
    ' $1 )
  # parse T_ma grid search end
  T_ma_f=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{T_ma}->{range}->[2]
    ' $1 )
  # parse horizontal gird spacing
  dx=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{dx}
    ' $1 )
  # parse timestep size
  dt=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{dt}
    ' $1 )
  # number of time integration steps
  t_f=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{TT}
    ' $1 )
  # parse MB curve fit type
  FIT=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{fit}
    ' $1 )
  # parse MB curve fit type
  k=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{k}
    ' $1 )
}

diagnostic_run()
{
  # Update the .SIF FILE with the model run specifc params
  sed "s#<DX>#"$dx"#g;
       s#<FIT>#"$FIT"#g;
       s#<KEY>#"$KEY"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<run_name>#"$run_name"#g;
       s#<SS_itters>#"$SS_itters"#g;" "./sifs/diagnostic.sif" > "./sifs/${run_name}.sif"

  # filepath to log file
  log_file="logs/${KEY}/${run_name}.log"

  # Run the model

  # TO DO: if log, else
  ElmerSolver "./sifs/${run_name}.sif" #| tee $log_file

  # add the params as a global attribute to the netcdf file
  python3 ../../src/thermal/add_attr.py -f "params/ref_params.sif" \
                                        -a "params"                \
                                           "result/${KEY}/nc/${run_name}.nc"

  # add the sif as a global attribute to the netcdf file
  python3 ../../src/thermal/add_attr.py -f "./sifs/${run_name}.sif" \
                                        -a "sif"                \
                                           "result/${KEY}/nc/${run_name}.nc"

  # # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

prognostic_run()
{
  # Update the .SIF FILE with the model run specifc params
  sed "s#<DX>#"$dx"#g;
       s#<dt>#"$dt"#g;
       s#<NT>#"$NT"#g;
       s#<KEY>#"$KEY"#g;
       s#<FIT>#"$FIT"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<RESTART>#"$RESTART"#g
       s#<run_name>#"$run_name"#g;
       s#<Dynamic_int>#"$Dynamic_int"#g
       s#<SS_itters>#"$SS_itters"#g;" "./sifs/prognostic.sif" > "./sifs/${run_name}.sif"

  # filepath to log file
  log_file="logs/${KEY}/${run_name}.log"

  # Run the model
  ElmerSolver "./sifs/${run_name}.sif" #| tee $log_file

  # add the params as a global attribute to the netcdf file
  python3 ../../src/thermal/add_attr.py -f "params/ref_params.sif" \
                                        -a "params"                \
                                           "result/${KEY}/nc/${run_name}.nc"
  # add the sif as a global attribute to the netcdf file
  python3 ../../src/thermal/add_attr.py -f "./sifs/${run_name}.sif" \
                                        -a "sif"                \
                                           "result/${KEY}/nc/${run_name}.nc"

  # # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

log_runtime()
{

  local OUT_fp="result/${KEY}/${KEY}.coupled_init.time_profile"

  if [ ! -f "$OUT_fp" ]; then
      echo "#dx dt NT offset T_ma runtime" |
      awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' >> \
                      $OUT_fp
  fi

  echo "${dx} ${dt} ${NT} ${offset} ${T_ma} ${runtime}" |
  awk -v OFS='\t'  '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' >> \
                  $OUT_fp
}



full_initialization(){
  # Variables that need to be set:
  #    diag_SS_itters
  #    FIT
  #    KEY
  #    dx
  #    offset
  #    T_ma
  #    dt
  #    t_f
  #    diag_SS_itters

  # parameter dictionary for griding NetCDFs
  param_dict="{\"T_ma\"   : ${T_ma},
               \"offset\" : ${offset}}"

  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # steady-state run
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  # let's allow for X number of steady state itterations
  SS_itters=$diag_SS_itters

  # make the run name based on model params
  run_name="${KEY}_dx_${dx}_MB_${offset}_OFF_Tma_${T_ma}_diag"

  # run the model for a given offset
  diagnostic_run $dx $FIT $KEY $offset $run_name $SS_itters

  # grid the NetCDF file written by the NetcdfUGRIDOutputSolver
  python3 ../../src/thermal/grid_data.py "result/${KEY}/nc/${run_name}.nc"      \
                                 -out_fn "result/${KEY}/gridded/${run_name}.nc" \
                                 -params "${param_dict}"

  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # transient run
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # Number of time interval based on dt
  NT=$(awk -v dt=$dt -v t_f=$t_f 'BEGIN {OFMT = "%.0f"; print (t_f/dt)}')
  # Execute interval for dynamics solvers,
  Dynamic_int=$(awk -v dt=$dt 'BEGIN {OFMT = "%.0f"; print (1.0/dt)}')

  # limit to 10 S.S. itters for transient runs
  SS_itters=$prog_SS_itters
  # diagnostic run is now restart variable
  RESTART="${run_name}.result"
  # prognostic run name
  run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_prog"

  # Start the timer
  start=$(date +%s.%N)

  # run the transient model with diagnostic solution as restart fiedl
  prognostic_run $dx $FIT $KEY $offset $run_name $SS_itters $restart $NT $dt

  # grid the NetCDF file written by the NetcdfUGRIDOutputSolver
  python3 ../../src/thermal/grid_data.py "result/${KEY}/nc/${run_name}.nc"      \
                                 -out_fn "result/${KEY}/gridded/${run_name}.nc" \
                                 -params "${param_dict}"
  # End the timer
  end=$(date +%s.%N)

  # Execution time of the solver
  runtime=$(awk -v start=$start -v end=$end 'BEGIN {print end - start}')

  # record the prognostic runtimes for future info
  log_runtime $dx $dt $NT $offset $T_ma $runtime
}

find_final_timestep()
{ #-----------------------------------------------------------------------------
  # Get the final timestep from a given NetCDF file.
  # Needed to determine which model runs ran out of walltime or memory alloc.
  #-----------------------------------------------------------------------------
  # $1 (NetCDF_fp) --->   file path to NetCDF file
  #-----------------------------------------------------------------------------

  final_timestep=$(ncks -d 't',-1 -v 't' $1 | tac | grep -m 1 "t = " | tr -s ' ' | cut -d " " -f 4)
  # get the final timestep from the input file |
  # reverse cat the file, is needed because time dimension format look similar to timestep info
  # get the final timestep |
  # squeeze any duplicate blank spaces |
  # get the actual time info by splitting
}

log_incomplete()
{ #-----------------------------------------------------------------------------
  # Log the current offset and T_ma as an incomplete run.
  #-----------------------------------------------------------------------------
  # $1 (incomp_file) --->   file path to ./run/${KEY}.incomplete file
  #-----------------------------------------------------------------------------

  # if the first time add comment at header
  if [[ ! -f $1 ]]; then
    echo "#T_ma offset" |
    awk -v OFS='\t' '{print $1 "\t" $2}' >> $1
  fi

  # print the current value to the passed file
  echo "${T_ma} ${offset}" |
  awk -v OFS='\t' '{print $1 "\t" $2}' >> $1
}

find_incomplete()
{ #-----------------------------------------------------------------------------
  # Create a list of incomplete model runs, which can be resubmitted
  #
  # Results are added to the "./run/${KEY}.incomplete" file
  #-----------------------------------------------------------------------------
  # $1 (KEY) --->   Glacier identifer
  #-----------------------------------------------------------------------------
  KEY=$1

  # parse the parameters from the json files
  parse_json "params/${KEY}.json"

  # set the incomplete files filename
  incomp_file="run/${KEY}.incomplete"

  #WARNING: hardcoding timstep
  dt=1.0

  # loop over the mass balance offsets
  for offset in $(seq -w $MB_0 $MB_s $MB_f); do
    # loop over mean annual air temps
    for T_ma in $(seq -w $T_ma_0 $T_ma_s $T_ma_f); do

      # Number of time interval based on dt
      NT=$(awk -v dt=$dt -v t_f=$t_f 'BEGIN {OFMT = "%.0f"; print (t_f/dt)}')
      # get the transient runname
      run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_prog"

      # gridded transient results don't exist
      if [[ ! -f "result/${KEY}/gridded/${run_name}.nc" ]]; then
        # record current offset and T_ma, which needs to be re-run
        log_incomplete ${incomp_file}
      else
        # since the file exists, check to see if the run finished
        find_final_timestep "result/${KEY}/gridded/${run_name}.nc"

        if [[ $final_timestep != $t_f ]]; then
          # record current offset and T_ma, which needs to be re-run
          log_incomplete ${incomp_file}
	      fi
      fi
    done
  done
}
