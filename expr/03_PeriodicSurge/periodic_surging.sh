#!/usr/bin/env bash

evenly_divisible()
{ # helper function to determine if two numbers are evenly divisible
  #               1 : true  | 0 : false
  echo $(awk -v a=$1 -v b=$2 'BEGIN{print((a / b % 1) == 0)}')
}

cycle_period()
{ # perform float addition
  echo $(awk -v a=$1 -v b=$2 'BEGIN{printf("%.3f", a + b)}')
}

number_of_cycles()
{ # Calculate the number of surge cycles for input params:
  ##############################################################################
  # $1 --> simulation length  [a]
  # $2 --> surge cycle period [a]
  ##############################################################################

  # test if simulation length is evenly divisible by the surge cycle period
  local remain=$(evenly_divisible $1 $2)

  # if not evenly divisible throw error
  if [ $remain -ne 1 ]; then
    echo
    echo "ERROR simulation length not evenly divisible by surge cycle"
    echo
    exit 1
  fi

  # divide the simulation length by the surge cycle period
  echo $(awk -v a=$1 -v b=$2 'BEGIN{printf("%.0f", a / b)}')
}

fill_timestep_sizes()
{ # fill the timestep_size array
  ##############################################################################
  # $1 --> number of surge cycles
  # $2 --> surging   timestep [a]
  # $3 --> quiescent timestep [a]
  ##############################################################################

  # timestep array lenghts
  local M=$((2*$1))

  # initialize the empty array
  local timestep_sizes=()

  for ((i = 1 ; i <= M ; i+=2)); do
    timestep_sizes+=($2)
    timestep_sizes+=($3)
  done

  # return the timestep sizes array
  echo "${timestep_sizes[*]}"
}

fill_timestep_intervals()
{ # fill the timestep_intervals array
  ##############################################################################
  # $1 --> number of surge cycles
  # $2 --> surging   timestep [a]
  # $3 --> quiescent timestep [a]
  # $4 --> surging   period   [a]
  # $5 --> quiescent period   [a]
  ##############################################################################

  # test if simulation length is evenly divisible by the surge cycle period
  local surge_remain=$(evenly_divisible $4 $2)
  local quisc_remain=$(evenly_divisible $5 $3)

  # if not evenly divisible throw error
  if [ $surge_remain -ne 1 ]; then
    echo
    echo "ERROR: surge period not evenly divisible by surging timestep"
    echo
    exit 1
  elif [ $quisc_remain -ne 1 ]; then
    echo
    echo "ERROR: quiescent period not evenly divisible by quiescent timestep"
    echo
    exit 1
  fi

  # timestep array lenghts
  local M=$((2*$1))

  # initialize the empty array
  local timestep_intervals=()

  for ((i = 1 ; i <= M ; i+=2)); do
    timestep_intervals+=($(awk -v a=$4 -v b=$2 'BEGIN{printf("%.0f", a / b)}'))
    timestep_intervals+=($(awk -v a=$5 -v b=$3 'BEGIN{printf("%.0f", a / b)}'))
  done

  # return the timestep sizes array
  echo "${timestep_intervals[*]}"
}

fill_dynamicexec_intervals()
{ # fill the timestep_intervals array
  ##############################################################################
  # $1 --> number of surge cycles
  # $2 --> dynamic surging   timestep [a]
  # $3 --> dynamic quiescent timestep [a]
  # $4 --> thermal surging   timestep [a]
  # $5 --> thermal quiescent timestep [a]
  ##############################################################################

  # test if simulation length is evenly divisible by the surge cycle period
  local surge_remain=$(evenly_divisible $2 $4)
  local quisc_remain=$(evenly_divisible $3 $5)

  # if not evenly divisible throw error
  if [ $surge_remain -ne 1 ]; then
    echo
    echo "ERROR (surging)  : dynamic dt not evenly divisible by quiescent dt"
    echo
    exit 1
  elif [ $quisc_remain -ne 1 ]; then
    echo
    echo "ERROR (quiescent): dynamic dt not evenly divisible by quiescent dt"
    echo
    exit 1
  fi

  # timestep array lenghts
  local M=$((2*$1))

  # initialize the empty array
  local dynamicexec_intervals=()

  for ((i = 1 ; i <= M ; i+=2)); do
    dynamicexec_intervals+=($(awk -v a=$2 -v b=$4 'BEGIN{printf("%.0f", a / b)}'))
    dynamicexec_intervals+=($(awk -v a=$3 -v b=$5 'BEGIN{printf("%.0f", a / b)}'))
  done

  # return the timestep sizes array
  echo "${dynamicexec_intervals[*]}"
}

calc_number_of_timesteps()
{ # sum the timestep array to calculate the number of timesteps
  sum=0
  for i in ${1[@]}; do
    let sum+=$i
  done

  echo $sum
}

parse_json()
{ # Parse model parameters from .json config file

  ##############################################################################
  # $1 (json_fp) --->   file path to json parameter file
  ##############################################################################

  # parse reference MB value
  MB=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{MB}
    ' $1 )
  # parse reference T_ma value
  T_ma=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{T_ma}
    ' $1 )
  # surging dynamic timestep
  SD_dt=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{SD_dt}
    ' $1 )
  # quiescent dynamic timestep
  QD_dt=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{QD_dt}
    ' $1 )
  # surging thermal timestep
  ST_dt=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{ST_dt}
    ' $1 )
  # quiescent thermal timestep
  QT_dt=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{QT_dt}
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
  TT=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{TT}
    ' $1 )
}

log_runtime()
{ # Log the simulation runtime as function of pertinent params
  #-----------------------------------------------------------------------------
  #
  # Variables:
  # ---------
  #  KEY     ---> glacier identifer
  #  dx      ---> mesh resolution. ${KEY}/mesh_dx${dx}/mesh.* should exists
  #  T_ma    ---> air temp. [C] @ $z_ref from "params/ref_params.sif" file
  #  offset  ---> Mass balance anomoly [m i.e. yr-1]
  #  ST_dt   ---> surging thermal timestep   [yr]
  #  SD_dt   ---> surging dynamic timestep   [yr]
  #  QT_dt   ---> quies.  thermal timestep   [yr]
  #  QD_dt   ---> quies.  dynamic timestep   [yr]
  #  SP      ---> active phase length        [yr]
  #  QP      ---> quiescent period length    [yr]
  #  TT      ---> Total length of simulation [yr]
  #  beta    ---> slip coeffiecent value
  #  runtime ---> Total simulation length    [sec]
  #-----------------------------------------------------------------------------
  local OUT_fp="result/${KEY}/${KEY}.periodic_surge.time_profile"

  if [ ! -f "$OUT_fp" ]; then
      echo "#dx offset T_ma ST_dt SD_dt QT_dt QD_dt SP QP TT beta runtime" |
      awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' >> \
      $OUT_fp
  fi

  echo "${dx} ${offset} ${T_ma} ${ST_dt} ${SD_dt} ${QT_dt} ${QD_dt} ${S_P} ${Q_P} ${TT} ${beta} ${runtime}" |
  awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' >> \
  $OUT_fp
}

periodic_run()
{ # Run a periodic surging simulation
  #-----------------------------------------------------------------------------
  # This function:
  #         1) upates the template sif file,
  #         2) runs the prognostic Elmer/Ice simulation,
  #         3) adds the .sif file and params as attributes to the NetCDF
  #
  # Variables:
  # ---------
  #  KEY          ---> glacier identifer
  #  dx           ---> mesh resolution. ${KEY}/mesh_dx${dx}/mesh.* should exists
  #  T_ma         ---> air temp. [C] @ $z_ref from "params/ref_params.sif" file
  #  offset       ---> Mass balance anomoly [m i.e. yr-1]
  #  SS_itters    ---> Number of S.S. itterations for diagnostic simulation
  #  run_name     ---> unique simulation identifer
  #  M            ---> length of timestep related arrays
  #  dt_arr       ---> timestep array of length $M            [yr]
  #  NT_arr       ---> timstep intervals array of length $M   [-]
  #  dyn_exec_arr ---> interval to exec stokes solver array of length $M [-]
  #-----------------------------------------------------------------------------

  # Update the .SIF FILE with the model run specifc params
  sed "s#<M>#"$M"#g;
       s#<DX>#"$dx"#g;
       s#<KEY>#"$KEY"#g;
       s#<beta>#"$beta"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<dt_arr>#""${dt_arr[*]}""#g;
       s#<NT_arr>#""${NT_arr[*]}""#g;
       s#<RESTART>#"$RESTART"#g
       s#<run_name>#"$run_name"#g;
       s#<SS_itters>#"$SS_itters"#g;
       s#<dyn_exec_int>#""${dyn_exec_arr[*]}""#g;" "./sifs/periodic_surge.sif" > "./sifs/${run_name}.sif"

  # filepath to log file
  log_file="${KEY}/logs/${run_name}.log"

  ElmerSolver ./sifs/${run_name}.sif

  # add the params as a global attribute to the netcdf file
  python3 ../../src/thermal/add_attr.py -f "params/ref_params.sif" \
                                        -a "params"                \
                                           "result/${KEY}/nc/${run_name}.nc"
  # add the sif as a global attribute to the netcdf file
  python3 ../../src/thermal/add_attr.py -f "./sifs/${run_name}.sif" \
                                        -a "sif"                \
                                           "result/${KEY}/nc/${run_name}.nc"

  # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

periodic_simulation()
{
  # parse the parameters from the json files
  parse_json "params/${KEY}.json"

  # check for any overwrites of the default parameters
  

}

# parse the parameters from the json files
parse_json "params/crmpt12.json"

offset=-0.41
T_ma=-8.5
KEY='crmpt12'
# diagnostic run is now restart variable
RESTART="crmpt12_dx_50_NT_100_dt_0.05_MB_-0.41_OFF_Tma_-8.5_B_0.001_pseudo_dt_1.0_NT_2000_recovery.result"
# limit to 10 S.S. itters for transient runs
SS_itters=10

beta=0.001

# run the model for 2000 years
TT=2000

# Surge intervals
S_P=2.00                       # (s)urge  (p)eriod
Q_P=18.0                       # (q)uies. (p)eriod

# calculate surge cycle period; from the set surge/quiescent interval
C_P=$(cycle_period $S_P $Q_P)  # (c)ycle  (p)eriod

# number of surge cycles
NC=$(number_of_cycles $TT $C_P)

# timestep array lenghts
M=$((2*NC))
# get the timestep sizes array
dt_arr=($(fill_timestep_sizes $NC $ST_dt $QT_dt))
# get the timestep intervals array
NT_arr=($(fill_timestep_intervals $NC $ST_dt $QT_dt $S_P $Q_P))
# get the array of intervals to execute the dynamic solver at
dyn_exec_arr=($(fill_dynamicexec_intervals $NC $SD_dt $QD_dt $ST_dt $QT_dt))

# prognostic run name
run_name="${KEY}_dx_${dx}_TT_${TT}_MB_${offset}_OFF_Tma_${T_ma}_B_${beta}_SP_${S_P}_QP_${Q_P}"

# create file form template
periodic_run
