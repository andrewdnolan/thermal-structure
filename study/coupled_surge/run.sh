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
       s#<beta>#"$beta"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<RESTART>#"$RESTART"#g
       s#<run_name>#"$run_name"#g;
       s#<SS_itters>#"$SS_itters"#g;" "./sifs/pseudo_surge.sif" > "./sifs/${run_name}.sif"

  # filepath to log file
  log_file="${KEY}/logs/${run_name}.log"

  # Run the model
  docker exec elmerenv /bin/sh -c "cd shared_directory/Thesis/thermal-structure/study/coupled_surge;
                                   ElmerSolver ./sifs/${run_name}.sif "
  # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

# parse the parameters from the json files
parse_json "params/crmpt12.json"

offset=-0.41
T_ma=-8.5
KEY='crmpt12'

# limit to 10 S.S. itters for transient runs
SS_itters=10

dt=0.1
NT=200
TT=20
beta=0.003

# diagnostic run is now restart variable
RESTART="crmpt12_dx_50_NT_2000_dt_1.0_MB_-0.41_OFF_Tma_-8.5_prog.result"

# # prognostic run name
run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_B_${beta}_pseudo"
# # run_name="test_noseasonality"

# run the transient model with diagnostic solution as restart fiedl
pseudo_run $dx $KEY $offset $run_name $SS_itters $restart $NT $dt

# grid the NetCDF file written by the NetcdfUGRIDOutputSolver
python3 ../../src/thermal/grid_data.py "result/${KEY}/nc/${run_name}.nc"      \
                               -out_fn "result/${KEY}/gridded/${run_name}.nc"

dt=1.0
NT=100
TT=100

# recovery run name
RESTART="${run_name}.result"
run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_B_${beta}_recovery"

beta=1.0
# run the transient model with diagnostic solution as restart fiedl
pseudo_run $dx $KEY $offset $run_name $SS_itters $restart $NT $dt

# grid the NetCDF file written by the NetcdfUGRIDOutputSolver
python3 ../../src/thermal/grid_data.py "result/${KEY}/nc/${run_name}.nc"      \
                               -out_fn "result/${KEY}/gridded/${run_name}.nc"
