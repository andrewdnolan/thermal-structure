#!/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# elmer2nc:
#   wrapper to parse, pre-process, and pass arguments to the fortran impelmentation
#   of `elmer2nc`
# https://www.golinuxcloud.com/bash-getopts/
# https://gitlab.awi.de/sicopolis/sicopolis/-/blob/master/get_input_files.sh
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function usage()
{
  echo -e "\n ./$(basename $0) -h --> shows usage"
  exit 1
}


function parse_time()
{ # Parse lines (with line numbers) from .result file. Piped to temporary file
  # which is then proccessed by the fortran program
  grep -n "Time" $src_fn | awk -F: '{print $1 " " $3}' > results/mesh/timesteps.dat
}


# list of arguments expected in the input
optstring=":h"

while getopts ${optstring} arg; do
  case ${arg} in
    h)
      echo "showing usage"
      usage
      ;;
    :)
      echo "$0: Must supply an argument to -$OPTARG/" >&2
      exit 1
      ;;
    ?)
      echo "Invalid option: -${OPTARG}."
      usage
      ;;
  esac
done
