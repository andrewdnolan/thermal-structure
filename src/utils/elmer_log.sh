#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# elmer_log.sh
#   - create entry in spingup time profile
#
#   Note:
#     - This script does not work correctly on OSX!!!!
#       Must be run from linux machine!!!!
#       getopt (which support long opts) only works on linux
#
#     - this take the same arguments (and than some) as the initialization
#       'prepare2submit' script. This is for the sake of convience to
#       unnecassary variale parsing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

check_vars_defined()
{
  local not_set=false
  for var in $2; do
    if [ -z ${!var} ]; then
      echo
      echo "Error:"
      echo "    ${var} is NOT set"
      local not_set=true
    fi
  done

  if [ $not_set = true ]; then
    echo
    echo "********************************************************************"
    echo "function: ${1} can not be run, not all necessary vars are set."
    echo "********************************************************************"
    echo
    exit 1
  fi
}

parse_make_sif_args()
{ #https://gist.github.com/cosimo/3760587
  #https://www.aplawrence.com/Unix/getopts.html
  #https://stackoverflow.com/a/7948533/10221482
  #https://www.bahmanm.com/2015/01/command-line-options-parse-with-getopt.html

  OPTS=`getopt -o hx:t:T:k:f:p:o:s:e: --long help,dx:,dt:,tt:,key:,fit:,ord:,off:,sif:,ET: -n 'parse-options' -- "$@"`

  if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

  usage="
  $(basename "$0") [-h] [--dx dx --dt dt --tt tt --key key --fit fit  --off off --sif sif ]
  --
  log the runtime info for an exection of ElmerSolver based on input params

  NOTE: this take the same arguments (and than some) as the initialization
        'prepare2submit' script. This is for the sake of convience to
        unnecassary variale parsing

  where:
      -h | --help ) show this help text
      -x | --dx   ) mesh resolution                         [m]
      -t | --dt   ) time step                               [y]
      -T | --tt   ) length of the simulation                [y]
      -k | --key  ) glacier id key
      -f | --fit  ) fit type to Young et al 2020. MB data
      -p | --ord  ) Order of spline used to prescribe MB
      -o | --off  ) offset to mass balance curve            [m i.e.q yr^-1]
      -s | --sif  ) path to template sif file
      -e | --ET   ) execution time                          [s]
      "
  # set default arg for those that need one
  local HELP=false

  eval set -- "$OPTS"
  while true; do
    case "$1" in
      -h | --help ) echo "$usage";exit 0;;
      -x | --dx   ) dx="$2";     shift 2;;
      -t | --dt   ) dt="$2";     shift 2;;
      -T | --tt   ) TT="$2";     shift 2;;
      -k | --key  ) KEY="$2";    shift 2;;
      -f | --fit  ) FIT="$2";    shift 2;;
      -p | --ord  ) k="$2";      shift 2;;
      -o | --off  ) OFFSET="$2"; shift 2;;
      -s | --sif  ) SIF="$2";    shift 2;;
      -e | --ET   ) runtime="$2";shift 2;;
      -- ) shift; break ;;
      * ) break ;;
    esac
  done

  SIF=$(sed -e 's/^"//' -e 's/"$//' <<<"$SIF")
  OFFSET=$(sed -e 's/^"//' -e 's/"$//' <<<"$OFFSET")

  check_vars_defined "parse_make_sif_args" "dx dt TT KEY FIT OFFSET SIF k runtime"
}


log_runtime()
{

  check_vars_defined "make_sif" "dx dt TT KEY FIT OFFSET SIF k runtime"

  OUT_fp="result/${KEY}/${KEY}.spinup.time_profile"

  if [ ! -f "$OUT_fp" ]; then
      echo "#dx dt NT OFFSET runtime" |
      awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' >> \
                      $OUT_fp
  fi

  # Calculate total number of time steps (NT) from dt and TT
  local NT=$( awk -v TT=$3 -v dt=$2 "BEGIN { print TT/dt }" )

  echo "${dx} ${dt} ${NT} ${OFFSET} ${runtime}" |
  awk -v OFS='\t'  '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' >> \
                  $OUT_fp
}


main()
{
  parse_make_sif_args "$@"

  # If all the check are passed, finally actually call the function
  log_runtime $dx $dt $TT $KEY $FIT $OFFSET $SIF
}

#run the whole thing.
main "$@"
