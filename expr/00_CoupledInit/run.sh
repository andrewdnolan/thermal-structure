#!/usr/bin/env bash

# source ../../src/utils/elmer_log.sh
source ../../src/utils/initialization.sh


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
  ElmerSolver "./sifs/${run_name}.sif" #| tee $log_file

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

  # # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

log_runtime()
{

  local OUT_fp="result/${KEY}/${KEY}.coupled_spinup.time_profile"

  if [ ! -f "$OUT_fp" ]; then
      echo "#dx dt NT offset runtime" |
      awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' >> \
                      $OUT_fp
  fi

  echo "${dx} ${dt} ${NT} ${offset} ${runtime}" |
  awk -v OFS='\t'  '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' >> \
                  $OUT_fp
}

# parse the parameters from the json files
parse_json "params/glc1-a.json"

# # make a list of mass balance offsets
# offsets=($(seq -w $MB_0 $MB_s $MB_f))
# # bash arrays zero indexed
# N=$((SLURM_ARRAY_TASK_ID-1))
# # get ith offset corresponing to SLURM_ARRAY_TASK_ID
# offset=${offsets[$N]}

# get ith offset and T_ma corresponing to SLURM_ARRAY_TASK_ID, nested for loops
# each script is junky. The cartesian product of the two arrays would be a more
# elegant solution, but I was having trouble gettig that to work.
# Refs:
#   - https://unix.stackexchange.com/questions/97814/array-cartesian-product-in-bash
#   - https://stackoverflow.com/questions/23363003/how-to-produce-cartesian-product-in-bash
#   - https://rosettacode.org/wiki/Cartesian_product_of_two_or_more_lists#UNIX_Shell

# hard code mb final because i accidently wrote over data
MB_f=-1.275

# count=1
# # loop over the mass balance offsets
# for offset in $(seq -w $MB_0 $MB_s $MB_f); do
#   # loop over mean annual air temps
#   for T_ma in $(seq -w -9.00 0.1 -7.00); do
#     # check if counter equals SLURM_ARRAY_TASK_ID
#     if [[ $count -eq $SLURM_ARRAY_TASK_ID ]]; then
#         break 2
#     fi
#     count=$((count+1))
#   done
# done

T_ma=-7.0
offset=-1.5
# set the glacier key
KEY='glc1-a'
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# steady-state run
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# let's allow for X number of steady state itterations
SS_itters=25

# make the run name based on model params
run_name="${KEY}_dx_${dx}_MB_${offset}_OFF_Tma_${T_ma}_diag_44intpoints"

# run the model for a given offset
diagnostic_run $dx $KEY $offset $run_name $SS_itters


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# transient run
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
dt=1.0
T_f=2000
# Number of time interval based on dt
NT=$(awk -v dt=$dt -v T_f=$T_f 'BEGIN {OFMT = "%.0f"; print (T_f/dt)}')
# Execute interval for dynamics solvers,
Dynamic_int=$(awk -v dt=$dt 'BEGIN {OFMT = "%.0f"; print (1.0/dt)}')

# limit to 10 S.S. itters for transient runs
SS_itters=10
# diagnostic run is now restart variable
RESTART="${run_name}.result"
# prognostic run name
run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_prog"

# Start the timer
start=$(date +%s.%N)

# run the transient model with diagnostic solution as restart fiedl
prognostic_run $dx $KEY $offset $run_name $SS_itters $restart $NT $dt

# grid the NetCDF file written by the NetcdfUGRIDOutputSolver
# python3 grid_data.py "${KEY}/nc/${run_name}.nc"

# End the timer
end=$(date +%s.%N)

# Execution time of the solver
runtime=$(awk -v start=$start -v end=$end 'BEGIN {print end - start}')

# record the prognostic runtimes for future info
log_runtime $dx $dt $NT $offset $runtime
