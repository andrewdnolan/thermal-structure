#!/usr/bin/env bash

diagnostic_run()
{
  # Update the .SIF FILE with the model run specifc params
  sed "s#<DX>#$dx#g;
       s#<FIT>#$FIT#g;
       s#<KEY>#$KEY#g;
       s#<T_mean>#$T_ma#g;
       s#<offset>#$offset#g;
       s#<run_name>#$run_name#g;
       s#<w_limiter>#$w_limiter#g;
       s#<SS_itters>#$SS_itters#g;" "./sifs/diagnostic.sif" > "./sifs/${run_name}.sif"

  # Run the model
  docker exec elmerenv /bin/sh -c "cd shared_directory/Thesis/thermal-structure/study/max_w_prescription;
                                   ElmerSolver ./sifs/${run_name}.sif "

  # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

prognostic_run()
{
  # Update the .SIF FILE with the model run specifc params
  sed "s#<DX>#$dx#g;
       s#<dt>#$dt#g;
       s#<NT>#$NT#g;
       s#<KEY>#$KEY#g;
       s#<FIT>#$FIT#g;
       s#<T_mean>#$T_ma#g;
       s#<offset>#$offset#g;
       s#<RESTART>#$RESTART#g
       s#<run_name>#$run_name#g;
       s#<w_limiter>#$w_limiter#g;
       s#<SS_itters>#$SS_itters#g;" "./sifs/prognostic.sif" > "./sifs/${run_name}.sif"

  # Run the model
  docker exec elmerenv /bin/sh -c "cd shared_directory/Thesis/thermal-structure/study/max_w_prescription;
                                   ElmerSolver ./sifs/${run_name}.sif "

  # # Remove the sif file
  rm "./sifs/${run_name}.sif"

}
source ../../src/utils/initialization.sh

# set the glacier key
KEY='glc1-a'

# parse the parameters from the json files
parse_json "params/${KEY}.json"

T_ma=-9.0
offset=-1.225

# different methods for prescribing max water content
for wmax_limit in "Density" "Depth"; do

  # parameter dictionary for griding NetCDFs
  param_dict="{\"T_ma\"      : ${T_ma},
               \"Delta_MB\"  : ${offset},
               \"wmax_limit\": \"${wmax_limit}\"}"

  if [[ $wmax_limit == "Density" ]]; then
    limiter="
    Enthalpy_h Upper Limit = Variable Densi, Phase Change Enthalpy \n
     \t\t\t\t\t\t\t\t\t\t\t\t\t
     real Procedure \"../../bin/Thermodynamics\" \"Limit_Enthalpy_rho\""
  elif [[ $wmax_limit == "Depth" ]]; then
    limiter="
    Enthalpy_h Upper Limit = Variable Depth, Phase Change Enthalpy \n
    \t\t\t\t\t\t\t\t\t\t\t\t\t
    real Procedure \"../../bin/Thermodynamics\" \"Limit_Enthalpy_surf\""
  fi

  # echo the variable so it's all one line and works with sed
  # ref: https://stackoverflow.com/a/42356201/10221482
  w_limiter=$(echo $limiter)

  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # steady-state run
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # let's allow for 25 number of steady state itterations
  SS_itters=25

  # make the run name based on model params
  run_name="${KEY}_dx_${dx}_MB_${offset}_OFF_Tma_${T_ma}_wlimit_${wmax_limit}_diag"

  # run the model for a given offset
  diagnostic_run $dx $KEY $offset $run_name $SS_itters
  # grid the NetCDF file written by the NetcdfUGRIDOutputSolver
  python3 ../../src/thermal/grid_data.py "glc1-a/nc/${run_name}.nc" -params "${param_dict}"

  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # transient run
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  dt=1.0
  NT=1000
  TT=1000
  # limit to 10 S.S. itters for transient runs
  SS_itters=10
  # diagnostic run is now restart variable
  RESTART="${run_name}.result"
  # prognostic run name
  run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_wlimit_${wmax_limit}_prog"

  # run the transient model with diagnostic solution as restart fiedl
  prognostic_run $dx $KEY $offset $run_name $SS_itters $restart $NT $dt
  # grid the NetCDF file written by the NetcdfUGRIDOutputSolver
  python3 ../../src/thermal/grid_data.py "glc1-a/nc/${run_name}.nc" -params "${param_dict}"

done
