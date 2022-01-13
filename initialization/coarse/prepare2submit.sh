#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prepare2submit.sh
#   - create neccesary files and job submission scripts for uncoupled init
#     of our input geometries.
#
#   Note:
#     - Files and submission scipts are group by geometry size classification
#       since runtime and memory usage should be similar amougnst geometries of
#       the same size. This flawed, but straight forward approach means the mem.
#       and runtime allocation in each `submit` script is dictated by the largest
#       geometry in each group.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_RUN_name(){
  #-----------------------------------------------------------------------------
  # Make the `RUN` named based on the required parameters passed
  # ----------------------------------------------------------------------------
  # $1 (KEY)      --->   glacier id key
  # $2 (TT)       --->   length of the simulation                [y]
  # $3 (dt)       --->   time step                               [y]
  # $4 (dx)       --->   mesh resolution                         [m]
  # $5 (OFFSET)   --->   offset to mass balance curve            [m i.e.q yr^-1]
  # $6 (FIT)      --->   fit type to Young et al 2020. MB data
  #-----------------------------------------------------------------------------

  local RUN="$1_$2a_dt_$3_dx_$4_MB_$5_OFF_$6"
  echo $RUN
}

make_sif(){
  #-----------------------------------------------------------------------------
  # Pass parameters to be replaced in sif template. Function calll results in an
  # `.sif` file that can be run by `ElmerSolver`.
  #
  # This function can be called internally, but the `check_args` function has
  # been written so that it can be called directly from the command line.
  #
  # ----------------------------------------------------------------------------
  # $1 (dx)       --->   mesh resolution                         [m]
  # $2 (dt)       --->   time step                               [y]
  # $3 (TT)       --->   length of the simulation                [y]
  # $4 (KEY)      --->   glacier id key
  # $5 (FIT)      --->   fit type to Young et al 2020. MB data
  # $6 (OFFSET)   --->   offset to mass balance curve            [m i.e.q yr^-1]
  # $7 (SIF)      --->   path to template sif file
  #-----------------------------------------------------------------------------

  # File paths to input topography data
  local Zb_fp="../../input_data/$4_bed.dat"
  local Zs_fp="../../input_data/$4_surf.dat"

  # Given the passed parametes make the `RUN` name variable
  local RUN=$( make_RUN_name $4 $3 $2 $1 $6 $5 )
  # Calculate total number of time steps (NT) from dt and TT
  local NT=$( awk -v TT=$3 -v dt=$2 "BEGIN { print TT/dt }" )

  # Update the .SIF FILE with the model run specifc params
  sed "s#<DX>#"$1"#g;
       s#<dt>#"$2"#g;
       s#<NT>#"$NT"#g;
       s#<RUN>#"$RUN"#g;
       s#<KEY>#"$4"#g;
       s#<FIT>#"$5"#g;
       s#<Zs_fp>#"$Zs_fp"#g;
       s#<Zb_fp>#"$Zb_fp"#g;
       s#<OFFSET>#"$6"#g" "$7" > "./sifs/${RUN}.sif"
}

# Declare arrays for the individual size classifications
declare -a  small=("crpmt12" "crpmt18-a" "crpmt18-b" "glac1-a" "glac1-b")
declare -a medium=("lilk-a" "lilk-b" "klun-a" "klun-b" "sprg")
declare -a  large=("twds-a" "twds-b" "fish" "klut-a" "klut-b")
# Declare a nested arrays containing the sub-arrays
declare -a groups=("small" "medium" "large")

for group in "${groups[@]}"; do
  len="$group[@]"
  for glac in "${!len}"; do
    echo $glac
  done
done

#make_sif 50 1 1000 lilk-a cubic_spline 0.0 "./sifs/simple_spinup.sif"
