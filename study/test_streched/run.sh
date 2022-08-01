#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run.sh:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set +x

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
  # parse horizontal gird spacing
  dx=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{dx}
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
       s#<Cfirn>#"$Cfirn"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<run_name>#"$run_name"#g;
       s#<SS_itters>#"$SS_itters"#g;
       s#<limit_type>#"$limit_type"#g;" "./sifs/diagnostic.sif" > "./sifs/${run_name}.sif"

  # filepath to log file
  log_file="${KEY}/logs/${run_name}.log"

  # Run the model
  ElmerSolver "./sifs/${run_name}.sif" | tee $log_file

  # number of steady state itteration required
  n_itter=$(tac "./${KEY}/mesh_dx${dx}/${run_name}.result" | \
             grep -m1 "Time" | tr -s ' ' | cut -d " " -f 2)

  # Convert result files into NetCDFs
  ../../src/elmer2nc/elmer2nc.sh -r "./${KEY}/mesh_dx${dx}/${run_name}.result" \
                                -m "./${KEY}/mesh_dx${dx}/" \
                                -t $n_itter                 \
                                -o "./${KEY}/nc/"

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
       s#<Cfirn>#"$Cfirn"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<RESTART>#"$RESTART"#g
       s#<run_name>#"$run_name"#g;
       s#<SS_itters>#"$SS_itters"#g;
       s#<limit_type>#"$limit_type"#g;" "./sifs/prognostic.sif" > "./sifs/${run_name}.sif"

  # filepath to log file
  log_file="${KEY}/logs/${run_name}.log"

  # Run the model
  ElmerSolver "./sifs/${run_name}.sif" | tee $log_file

  # # Convert result files into NetCDFs
  # ../../src/elmer2nc/elmer2nc.sh -r "./${KEY}/mesh_dx${dx}/${run_name}.result" \
  #                               -m "./${KEY}/mesh_dx${dx}/" \
  #                               -t $NT               \
  #                               -o "./${KEY}/nc/"

  # # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

# parse the parameters from the json files
parse_json "params/glc1-a.json"

# let's allow for X number of steady state itterations
SS_itters=25
# set the glacier key
KEY='glc1-a'
#Mean annual air temperate
T_ma=-9.00
offset=-1.225
limit_type='rho' # 'surf' or 'rho'


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# steady-state run
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# make the run name based on model params
run_name="${KEY}_dx_${dx}_MB_${offset}_OFF_Tma_${T_ma}_diag"

# run the model for a given offset
diagnostic_run $dx $KEY $offset $run_name $SS_itters

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# transient run
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
dt=1.0
NT=250
# limit to 10 S.S. itters for transient runs
SS_itters=10
# diagnostic run is now restart variable
RESTART="${run_name}.result"
# prognostic run name
run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_prog"

# run the transient model with diagnostic solution as restart fiedl
prognostic_run $dx $KEY $offset $run_name $SS_itters $restart $NT $dt
