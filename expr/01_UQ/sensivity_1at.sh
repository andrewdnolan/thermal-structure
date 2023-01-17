#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TO DO:
#    [ ] set global flag for wether to write log files. No need on westgrid becuase
#        of the *.out and *.err files
#
#    [ ] set optional flag for executing the Elmer commands from the docker container
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
diagnostic_run()
{
  #-----------------------------------------------------------------------------
  # Run a diagnostic simulation. This function:
  #         1) upates the template sif file,
  #         2) runs the diagnostic Elmer/Ice simulation,
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
  #-----------------------------------------------------------------------------

  # Update the .SIF FILE with the model run specifc params
  sed "s#<C_firn>#"$C_firn"#g;
       s#<w_en>#"$w_en"#g;
       s#<w_aq>#"$w_aq"#g;
       s#<h_aq>#"$h_aq"#g;
       s#<f_dd>#"$f_dd"#g;
       s#<r_frac>#"$r_frac"#g;
       s#<Q_geo>#"$Q_geo"#g;
       s#<DX>#"$dx"#g;
       s#<FIT>#"$FIT"#g;
       s#<KEY>#"$KEY"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<run_name>#"$run_name"#g;
       s#<RESTART>#"$restart"#g;
       s#<SS_itters>#"$SS_itters"#g;" "./sifs/diagnostic.sif" > "./sifs/${run_name}.sif"

  # filepath to log file
  log_file="logs/${KEY}/${run_name}.log"

  # Run the model

  # TO DO: if log, else
  ElmerSolver "./sifs/${run_name}.sif" | tee $log_file

  # # add the sif as a global attribute to the netcdf file
  # python3 ../../src/thermal/add_attr.py -f "./sifs/${run_name}.sif" \
  #                                       -a "sif"                \
  #                                          "result/${KEY}/nc/${run_name}.nc"
  #
  # # Remove the sif file
  rm "./sifs/${run_name}.sif"

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
       s#<h_aq>#"$h_aq"#g;
       s#<f_dd>#"$f_dd"#g;
       s#<r_frac>#"$r_frac"#g;
       s#<Q_geo>#"$Q_geo"#g;
       s#<DX>#"$dx"#g;
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
  ElmerSolver "./sifs/${run_name}.sif" | tee $log_file

  # # add the params as a global attribute to the netcdf file
  # python3 ../../src/thermal/add_attr.py -f "params/ref_params.sif" \
  #                                       -a "params"                \
  #                                          "result/${KEY}/nc/${run_name}.nc"
  # # add the sif as a global attribute to the netcdf file
  # python3 ../../src/thermal/add_attr.py -f "./sifs/${run_name}.sif" \
  #                                       -a "sif"                \
  #                                          "result/${KEY}/nc/${run_name}.nc"
  #
  # # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

# # parse the i-th line of the parameters file and unpack the variables
# IFS="," read -r id H_aq f_dd r C_firn w_en w_aq Q_geo  \
#                 <<< "$(sed -n "${i}p" seed_123456_d7_n16.samples)"



KEY='crmpt12'
dx=50
T_ma=-8.5
offset=-0.35
SS_itters=25
Dynamic_int=10
dt=0.1
NT=30000

# parse the i-th line of the parameters file and unpack the variables
IFS="," read -r id h_aq f_dd r_frac C_firn w_en w_aq Q_geo  \
                <<< "1,3.0,0.006,0.3,0.05,0.01,0.10,0.055"


for C_firn in 0.15; do

  run_name="crmpt12-diag_C_${C_firn}_25itters_Tma_${T_ma}"

  SS_itters=25

  diagnostic_run

  # diagnostic run is now restart variable
  RESTART="${run_name}.result"

  # run_name="w_aq_depth_haq${h_aq}_fdd_${f_dd}_prog_seasonal_rho"

  run_name="crmpt12_prog_C_${C_firn}_3kya_Tma_${T_ma}"

  SS_itters=10
  prognostic_run

done
