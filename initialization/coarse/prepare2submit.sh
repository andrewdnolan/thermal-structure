#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prepare2submit.sh
#   - create neccesary files and job submission scripts for uncoupled init
#     of our input geometries.
#
#   Note:
#     - This script does not work correctly on OSX!!!!
#       Must be run from linux machine!!!!
#       getopt (which support long opts) only works on linux
#
#     - Files and submission scipts are group by geometry size classification
#       since runtime and memory usage should be similar amougnst geometries of
#       the same size. This flawed, but straight forward approach means the mem.
#       and runtime allocation in each `submit` script is dictated by the largest
#       geometry in each group.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parse_args()
{ #https://gist.github.com/cosimo/3760587
  #https://www.aplawrence.com/Unix/getopts.html
  #https://stackoverflow.com/a/7948533/10221482
  #https://www.bahmanm.com/2015/01/command-line-options-parse-with-getopt.html

  OPTS=`getopt -o fh: --long force,help, -n 'parse-options' -- "$@"`

  if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

  # echo "$OPTS"
  eval set -- "$OPTS"

  # Set deaultsh
  #VERBOSE=false
  HELP=false
  FORCE=false
  #DRY_RUN=false
  STACK_SIZE=0

  while true; do
    case "$1" in
      -f | --force ) FORCE=true; shift ;;
      -h | --help )   HELP=true; shift ;;
      -n | --dry-run ) DRY_RUN=true; shift ;;
      -s | --stack-size ) STACK_SIZE="$2"; shift 2;;
      -- ) shift; break ;;
      * ) break ;;
    esac
  done

  echo VERBOSE=$VERBOSE
  echo HELP=$HELP
  echo DRY_RUN=$DRY_RUN
  echo STACK_SIZE=$STACK_SIZE

}
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
  # parse time step size
  dt=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{dt}
    ' $1 )
  # parse MB curve fit type
  TT=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{TT}
    ' $1 )
  # parse MB curve fit type
  FIT=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{fit}
    ' $1 )
}

make_RUN_name()
{
  #-----------------------------------------------------------------------------
  # Make the `RUN` named based on the required parameters passed
  # ----------------------------------------------------------------------------
  # $KEY       --->   glacier id key
  # $TT        --->   length of the simulation                [y]
  # $dt        --->   time step                               [y]
  # $dx        --->   mesh resolution                         [m]
  # $OFFSET    --->   offset to mass balance curve            [m i.e.q yr^-1]
  # $FIT       --->   fit type to Young et al 2020. MB data
  #-----------------------------------------------------------------------------

  local RUN="${KEY}_${TT}a_dt_${dt}_dx_${dx}_MB_${OFFSET}_OFF_${FIT}"
  echo $RUN
}

make_sif()
{
  #-----------------------------------------------------------------------------
  # Pass parameters to be replaced in sif template. Function calll results in an
  # `.sif` file that can be run by `ElmerSolver`.
  #
  # This function can be called internally, but the `check_args` function has
  # been written so that it can be called directly from the command line.
  #
  # ----------------------------------------------------------------------------
  # $dx       --->   mesh resolution                         [m]
  # $dt       --->   time step                               [y]
  # $TT       --->   length of the simulation                [y]
  # $KEY      --->   glacier id key
  # $FIT      --->   fit type to Young et al 2020. MB data
  # $OFFSET   --->   offset to mass balance curve            [m i.e.q yr^-1]
  # $SIF      --->   path to template sif file
  #-----------------------------------------------------------------------------

  # File paths to input topography data
  local Zb_fp="../../input_data/${KEY}_bed.dat"
  local Zs_fp="../../input_data/${KEY}_surf.dat"

  # Given the passed parametes make the `RUN` name variable
  local RUN=$( make_RUN_name $KEY $TT $dt $dx $OFFSET $FIT )
  # Calculate total number of time steps (NT) from dt and TT
  local NT=$( awk -v TT=$3 -v dt=$2 "BEGIN { print TT/dt }" )

  # Update the .SIF FILE with the model run specifc params
  sed "s#<DX>#"$dx"#g;
       s#<dt>#"$dt"#g;
       s#<NT>#"$NT"#g;
       s#<RUN>#"$RUN"#g;
       s#<KEY>#"$KEY"#g;
       s#<FIT>#"$FIT"#g;
       s#<Zs_fp>#"$Zs_fp"#g;
       s#<Zb_fp>#"$Zb_fp"#g;
       s#<OFFSET>#"$OFFSET"#g" "$7" > "./sifs/${RUN}.sif"
}

check_inout_files()
{
  for ext in "in" "out"; do

    # check if file already exists
    if [ -f "./run/${group}.${ext}" ] && [ $FORCE = false ]; then
      echo
      echo "********************************************************************"
      echo -e "Input file: \"./run${group}.${ext}\" already exists.\n"\
      "  Either manually remove the \"./run/${group}.${ext}\" file or\n"\
      "  run this script with the \" -FORCE \" flag to overwrite files."
      echo "********************************************************************"
      echo
      exit 1

    # if file exists but force flag passed, delete it
  elif [ -f "$./run/{group}.${ext}" ] && [ $FORCE = true ]; then
      # remove the existing file
      rm -f "./run/${group}.${ext}"
    fi

  done
}

make_input_file()
{
  #-----------------------------------------------------------------------------
  # $MB_0  --->   start  of MB gridsearch                 [m i.e.q yr^-1]
  # $MB_f  --->   end    of MB gridsearch                 [m i.e.q yr^-1]
  # $MB_s  --->   stride of MB gridsearch                 [m i.e.q yr^-1]
  # $dx    --->   mesh resolution                         [m]
  # $TT    --->   length of the simulation                [y]
  # $dt    --->   time step                               [y]
  # $KEY   --->   glacier id key
  # $FIT   --->   fit type to Young et al 2020. MB data
  #-----------------------------------------------------------------------------

  for OFFSET in $(seq -w $MB_0 $MB_s $MB_f);do
    echo "make_sif $dx $dt $TT $KEY $FIT $OFFSET "./sifs/simple_spinup.sif" " \
         >> "./run/${group}.in"
  done

}

make_output_file()
{
  #-----------------------------------------------------------------------------
  # $MB_0  --->   start  of MB gridsearch                 [m i.e.q yr^-1]
  # $MB_f  --->   end    of MB gridsearch                 [m i.e.q yr^-1]
  # $MB_s  --->   stride of MB gridsearch                 [m i.e.q yr^-1]
  # $KEY   --->   glacier id key
  # $TT    --->   length of the simulation                [y]
  # $dt    --->   time step                               [y]
  # $dx    --->   mesh resolution                         [m]
  # $FIT   --->   fit type to Young et al 2020. MB data
  #-----------------------------------------------------------------------------

  for OFFSET in $(seq -w $MB_0 $MB_s $MB_f);do

    # Given the passed parametes and the current offset make `RUN` name variable
    local RUN=$( make_RUN_name $KEY $TT $dt $dx $OFFSET $FIT )

    echo  ../../src/elmer2nc/elmer2nc.sh  -r "./result/${KEY}/mesh_dx${dx}/${RUN}.result" \
                                          -m "./result/${KEY}/mesh_dx${dx}/" \
                                          -t $NT \
                                          -o "./result/${KEY}/nc/" \
                                          >> "./run/${group}.out"
  done
}

# make_submit_scipt(){
#
# }
################################################################################
FORCE=false

# Declare arrays for the individual size classifications
declare -a  small=("crmpt12" "crmpt18-a" "crmpt18-b" "glc1-a" "glc1-b")
declare -a medium=("lilk-a" "lilk-b" "klun-a" "klun-b" "sprg")
declare -a  large=("twds-a" "twds-b" "fish" "klut-a" "klut-b")
# Declare a nested arrays containing the sub-arrays
declare -a groups=("small" "medium" "large")

# for group in "${groups[@]}"; do
#   # get elements of the group
#   len="$group[@]"
#
#   # check whether in/out files exist
#   check_inout_files $group $FORCE
#
#   # loop over elements of the group
#   for KEY in "${!len}"; do
#
#     # extract parameter values from the json file
#     parse_json "./params/${KEY}.json"
#
#     # Add the `.sif` creation command to the input .txt file
#     make_input_file $MB_0 $MB_s $MB_f $dx $TT $dt $KEY
#
#     # Add the  NetCDF creation command to the output .txt file
#     make_output_file $MB_0 $MB_s $MB_f $KEY $TT $dt $dx $FIT
#
#   done
# done

# For formating to work correctly, this variable has to be printed within quotes
template=$(cat << EOF
#!/bin/bash
#SBATCH --array=1-<NJ>%<JS>                  # <NJ> jobs that run <JS> at a time
#SBATCH --job-name=<JOB_NAME>           # base job name for the array
#SBATCH --mem-per-cpu=<MEM>                     # maximum <MEM>MB per job
#SBATCH --time=<RUN_TIME>                      # maximum walltime per job
#SBATCH --nodes=1                                  # Only one node is needed
#SBATCH --ntasks=1                                 # These are serial jobs
#SBATCH --mail-type=ALL                            # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu       # email to spend updates too
#SBATCH --output=logs/<JOB_NAME>_%A_%a.out  # standard output
#SBATCH --error=logs/<JOB_NAME>_%A_%a.err   # standard error
# in the previous two lines %A" is replaced by jobID and "%a" with the array index

#Load cedar module file
source ../../config/modulefile.cc.cedar

# Get the command to create run specific .sif file
CREATE=\$( sed -n "\${SLURM_ARRAY_TASK_ID}p" <in_fp> )

# strip .sif file name from the creation command
SIF=\$(awk '{split(\$0, array, " "); print \$NF}' <<< "\$CREATE")

# Get the command to convert from .result to NetCDF
CONVERT=\$( sed -n "\${SLURM_ARRAY_TASK_ID}p" <in_fp> )

# Execute the .sif file
ElmerSolver \$SIF > logs/\${SIF##*/}.log
EOF
)

echo "$template" > "./run/test.sh"

sed "s#<NJ>#10#g" "./run/test.sh"

# CREATE=$( sed -n "1p" "./run/small.in" )
#
# # strip .sif file name from the creation command
# SIF=$(awk '{split($0, array, " "); print $NF}' <<< "$CREATE")
#
# echo $SIF
# # echo "$template"
