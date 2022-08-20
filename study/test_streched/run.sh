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
       s#<heat_source>#"$heat_source"#g;" "./sifs/diagnostic.sif" > "./sifs/${run_name}.sif"


  docker exec elmerenv /bin/sh -c "cd shared_directory/Thesis/thermal-structure/study/test_streched;
                                   ElmerSolver ./sifs/${run_name}.sif "

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
       s#<SS_itters>#"$SS_itters"#g;
       s#<Dynamic_int>#"$Dynamic_int"#g
       s#<heat_source>#"$heat_source"#g;" "./sifs/prognostic.sif" > "./sifs/${run_name}.sif"

  docker exec elmerenv /bin/sh -c "cd shared_directory/Thesis/thermal-structure/study/test_streched;
                                  ElmerSolver ./sifs/${run_name}.sif "
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
T_ma=-7.00
offset=-1.225
# 500 year prognostic runs
T_f=1000

# for heat_source in "volumetric" "mass";do
for heat_source in "volumetric";do

  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # steady-state run
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  # make the run name based on model params
  run_name="${KEY}_dx_${dx}_${heat_source}_diag"

  # run the model for a given offset
  diagnostic_run $dx $KEY $offset $run_name $SS_itters

  # parameter dictionary for griding NetCDFs
  param_dict="{\"T_ma\"      : ${T_ma},
               \"Delta_MB\"  : ${offset},
               \"heat_source\": \"${heat_source}\"}"

  # grid the NetCDF file written by the NetcdfUGRIDOutputSolver
  python3 ../../src/thermal/grid_data.py "glc1-a/nc/${run_name}.nc" \
                                         -params "${param_dict}"

  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # transient run
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  # test multiple timesteps [a]
  # for dt in 1.0 0.5 0.25 0.1; do
  for dt in 0.25; do
    # Number of time interval based on dt
    NT=$(awk -v dt=$dt -v T_f=$T_f 'BEGIN {OFMT = "%.0f"; print (T_f/dt)}')
    # Execute interval for dynamics solvers,
    Dynamic_int=$(awk -v dt=$dt 'BEGIN {OFMT = "%.0f"; print (1.0/dt)}')

    # limit to 10 S.S. itters for transient runs
    SS_itters=10
    # diagnostic run is now restart variable
    RESTART="${run_name}.result"
    # prognostic run name
    # run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_${heat_source}_prog"
    run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_vol_prog"


    # run the transient model with diagnostic solution as restart fiedl
    prognostic_run $dx $KEY $offset $run_name $SS_itters $restart $NT $dt

    # parameter dictionary for griding NetCDFs
    param_dict="{\"T_ma\"      : ${T_ma},
                 \"Delta_MB\"  : ${offset},
                 \"dt\"        : ${dt},
                 \"heat_source\": \"${heat_source}\"}"

    # grid the NetCDF file written by the NetcdfUGRIDOutputSolver
    python3 ../../src/thermal/grid_data.py "glc1-a/nc/${run_name}.nc" \
                                           -params "${param_dict}"
  done
done
