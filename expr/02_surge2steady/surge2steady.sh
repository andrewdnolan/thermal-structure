#!/usr/bin/env bash

evenly_divisible()
{ # helper function to determine if two numbers are evenly divisible
  #               1 : true  | 0 : false
  echo $(awk -v a=$1 -v b=$2 'BEGIN{print((a / b % 1) == 0)}')
}

log_runtime()
{

  local OUT_fp="result/${KEY}/${KEY}.surge2steady.time_profile"

  if [ ! -f "$OUT_fp" ]; then
      echo "#dx dt NT offset T_ma beta runtime" |
      awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' >> \
                      $OUT_fp
  fi

  echo "${dx} ${dt} ${NT} ${offset} ${T_ma} ${beta} ${runtime}" |
  awk -v OFS='\t'  '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' >> \
                  $OUT_fp
}

prognostic_run()
{ 
  #-----------------------------------------------------------------------------
  # Run a prognostic simulation. This function can either:
  #        (1) run an $((NT/dt)) pseduo surge with the prescribed $beta value
  #                              or
  #        (2) run an $((NT/dt)) recovery (i.e $beta==1) from a pseudo surge 
  #
  # Variables:
  # ---------
  #  KEY       ---> glacier identifer
  #  dx        ---> mesh resolution. ${KEY}/mesh_dx${dx}/mesh.* should exists
  #  T_ma      ---> mean annual air temp. [C] @ $z_ref from "params/ref_params.sif" file
  #  offset    ---> Mass balance anomoly [m i.e. yr-1]
  #  SS_itters ---> Number of S.S. itterations for diagnostic simulation
  #  run_name  ---> unique simulation identifer
  #  beta      ---> Basal slip coefficent [??]
  #  dt        ---> timestep              [yr]
  #  NT        ---> Number of timesteps   [-]
  #  RESTART   ---> filename of the simulation to use as the intial condition
  #
  #-----------------------------------------------------------------------------

  # Update the .SIF FILE with the model run specifc params
  sed "s#<DX>#"$dx"#g;
       s#<dt>#"$dt"#g;
       s#<NT>#"$NT"#g;
       s#<KEY>#"$KEY"#g;
       s#<beta>#"$beta"#g;
       s#<z_lim>#"$z_lim"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<RESTART>#"$RESTART"#g
       s#<run_name>#"$run_name"#g;
       s#<Dynamic_int>#"$Dynamic_int"#g
       s#<SS_itters>#"$SS_itters"#g;" "./sifs/pseudo_surge.sif" > "./sifs/${run_name}.sif"

  # filepath to log file
  log_file="${KEY}/logs/${run_name}.log"

  # Run the model
  ElmerSolver ./sifs/${run_name}.sif #| tee $log_file

  # Remove the sif file
  rm "./sifs/${run_name}.sif"

}



# RESTART! 
# cycle 2, flag variable which should be true is passed

# FOR SURGE
# beta --> (passed over cli)
# S_TT --> surge length (passed over cli)
# S_dt --> surge timestep, mechanics and thermodynamics are executed at the same timestep, reasonable deafult would be nice
# Dynamic_int=1

# FOR RECOVERY
#
#

surge2steady()
{ 
  # Run two simultions: (1) the surge, (2) the recovery

  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  #             Pseudo Surge (Part 1)                   #
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  # Number of time interval based on dt
  NT=$(awk -v dt=$ST_dt -v t_f=$SP 'BEGIN {OFMT = "%.0f"; print (t_f/dt)}')
  # Dynamic interval, how often to execute thermal vesus stokes solvers
  Dynamic_int=$(awk -v T_dt=$ST_dt -v D_dt=$SD_dt 'BEGIN {OFMT = "%.0f"; print (D_dt/T_dt)}')

  # set the timestep equal to the thermodynamic timestep, assumed to lower the stokes
  dt=$ST_dt

  # prognostic run name
  if [[ $cycle2 = True ]]; then 
    run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_B_${beta}_pseudo_C2"
  else 
    run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_B_${beta}_pseudo"
  fi

  # Start the timer
  start=$(date +%s.%N)

  # run the transient model with diagnostic solution as restart fiedl
  prognostic_run $dx $KEY $offset $run_name $SS_itters $restart $NT $dt

  # End the timer
  end=$(date +%s.%N)

  # Execution time of the solver
  runtime=$(awk -v start=$start -v end=$end 'BEGIN {print end - start}')

  # record the prognostic runtimes for future info
  log_runtime $dx $dt $NT $offset $T_ma $beta $runtime


  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  #            Recovery from surge (Part 2)             #
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  # Number of time interval based on dt
  NT=$(awk -v dt=$QT_dt -v t_f=$QP 'BEGIN {OFMT = "%.0f"; print (t_f/dt)}')
  # Dynamic interval, how often to execute thermal vesus stokes solvers
  Dynamic_int=$(awk -v T_dt=$QT_dt -v D_dt=$QD_dt 'BEGIN {OFMT = "%.0f"; print (D_dt/T_dt)}')

  # set the timestep equal to the thermodynamic timestep, assumed to lower the stokes
  dt=$QT_dt

  # recovery run name
  RESTART="${run_name}.result"
  run_name="${run_name}_NT_${TT}_recovery"

  beta=1.0
  
  # Start the timer
  start=$(date +%s.%N)
  # run the transient model with diagnostic solution as restart fiedssl
  prognostic_run $dx $KEY $offset $run_name $SS_itters $restart $NT $dt

  # End the timer
  end=$(date +%s.%N)

  # Execution time of the solver
  runtime=$(awk -v start=$start -v end=$end 'BEGIN {print end - start}')

  # record the prognostic runtimes for future info
  log_runtime $dx $dt $NT $offset $T_ma $beta $runtime

}