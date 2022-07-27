#!/bin/bash


source ../../src/utils/initialization.sh

pseudo_run()
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
       s#<SS_itters>#"$SS_itters"#g;" "./sifs/pseudo_surge.sif" > "./sifs/${run_name}.sif"

  # filepath to log file
  log_file="${KEY}/logs/${run_name}.log"

  # Run the model
  ElmerSolver "./sifs/${run_name}.sif" | tee $log_file

  # # Remove the sif file
  # rm "./sifs/${run_name}.sif"

}


# parse the parameters from the json files
parse_json "params/glc1-a.json"

offset=-1.225
T_ma=-9.00
KEY='glc1-a'

# limit to 10 S.S. itters for transient runs
SS_itters=10

dt=0.1
NT=100
TT=10

# diagnostic run is now restart variable
RESTART="glc1-a_dx_50_NT_2000_dt_1.0_MB_-1.225_OFF_Tma_-9.00_prog.result"

# prognostic run name
run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_pseudo"

# run the transient model with diagnostic solution as restart fiedl
pseudo_run $dx $KEY $offset $run_name $SS_itters $restart $NT $dt

# grid the NetCDF file written by the NetcdfUGRIDOutputSolver
python3 ../couple_init/grid_data.py "${KEY}/nc/${run_name}.nc"
