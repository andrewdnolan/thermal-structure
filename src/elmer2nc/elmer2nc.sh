#!/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# elmer2nc:
#   wrapper to parse, pre-process, and pass arguments to the fortran impelmentation
#   of `elmer2nc`
# https://www.golinuxcloud.com/bash-getopts/
# https://gitlab.awi.de/sicopolis/sicopolis/-/blob/master/get_input_files.sh
#
# TODO: this is a much better way to parse, lets us have mutli cahracter cmd line
# options:
#   - https://github.com/amanzi/amanzi/blob/master/bootstrap.sh#L603
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function usage()
{
  echo -e "\n Usage ./$(basename $0) -r <result_fp> [further options...]\n"\
  "     Convert Elmer's .result file to a NetCDF using Fortran code found in the src dir.\n"\
  "     By deafault outfile file will have same name as input but with \`.nc\` file extension\n\n"\
  "     [-r <fp>]  => file path to .result file we  want to converts\n"\
  "     [-m <dir>] => mesh directory. If not supplied infered from <result_fp>\n"\
  "     [-t <NT>]  => Number of timesteps within <result_fp> {integer}\n"\
  "     [-o <dir>] => output directory to write .nc file. If not supplied will write to -m <dir> \n"\
  "     [-s      ] => "

  exit 1
}


function parse_time()
{ # Parse lines (with line numbers) from .result file. Piped to temporary file
  # which is then proccessed by the fortran program
  grep -n "Time" $src_fn | awk -F: '{print $1 " " $3}' > results/mesh/timesteps.dat
}

function check_args()
{
  # list of arguments expected in the input
  optstring=":r:m:t:o:h"

  while getopts ${optstring} arg; do
    case ${arg} in
      r) result_fn=$OPTARG ;;
      m) mesh_db=$OPTARG   ;;
      t) NT=$OPTARG;;
      o) out_fp=$OPTARG ;;
      :) echo "$0: Must supply an argument to -$OPTARG/" >&2; exit 1 ;;
    h|?) usage; exit 1     ;;
    esac
  done

  # Make sure a .result file was passed, and that it exists
  if [[ ! "$result_fn" ]]; then
    echo -e "\n Error: need to supply -r <result_fp> \n"; exit 1
  else
    if [[ ! -f "${result_fn}" ]]; then
      echo "\n Error: -r <result_fp> supplied does not exist \n"; exit 1
    fi
  fi

  # if mesh_db is not passed, then infer the mesh_db from the input .result file path
  # this assumes that the .result file is in the same directory as the mesh.* files
  # if this is not the case then this will throw and error
  if [[ ! "$mesh_db" ]]; then
    #stip filename from filepath, and use as our mesh_db
    mesh_db=${result_fn%/*}/
  fi

  # test to make sure $mesh_db contains neccessary mesh.* files
  if [[ -z "$(find $mesh_db -name "mesh.*")" ]]; then
    echo "\n No mesh.* file were found in mesh_db: ${mesh_db} \n"; exit 1
  fi

  # make sure NT was passed (for now)
  if [[ ! "$NT" ]]; then
    echo "\n Error: need to supply -t <NT> \n"; exit 1
  fi


  # check if out_fp was specified, if so make sure it exists. If not set to mesh_db
  if [[ ! "$out_fp" ]]; then
    out_fp="${result_fn%.*}.nc"
  else
    if [[ (! -d "$out_fp") && (! -z "$(echo ${out_fp} | grep -o ".nc")") ]]; then
      echo -e "\n Error: you supplied an output filename not path. Remove the file name from "\
              " -o <dir> and just specify path. Filename is will be the same as the input \n"
      exit 1;
    elif [[ ! -d "$out_fp" ]]; then
      echo "\n Error: supplied -o <dir> does not exist. \`mkdir\` it yourself if you need to \n"
      exit 1;
    else
      fn=${result_fn##*/}         # get the base input filename
      out_fp=${out_fp}${fn%.*}.nc # concat out_fp dir with new out filename
    fi

  fi
}

################################################################################
check_args $*

./bin/elmer2nc $result_fn $mesh_db $out_fp $NT
################################################################################
