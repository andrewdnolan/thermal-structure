#!/bin/bash

################################################################################
#
################################################################################

usage () {
  echo -e "\n Usage ./$(basename $0) dx out_dir \n"\
	  "[ dx ]     => horizontal mesh resolution in meters \n"\
	  "[out_dir]  => name of folder to write mesh too \n"
    "[bedfile]  => filepath (from top dir) to .dat file to be meshed \n"

  exit 1
}

update_grd () {
  echo "start=${start}"
  echo "end=${end}"
  echo "Nx=${Nx}"
  echo "grdfile=${grdfile}"
  sed "s/^Subcell Limits 1 =.*/Subcell Limits 1 = "${start}" "${end}"/g;
       s/^Element Divisions 1 =.*/Element Divisions 1 = "${Nx}"/g" "${grdfile}" > "tmp.grd"
}

main () {
  dx=$1      # horizontal grid cell resolution [m]
  out_dir=$2 # directory name to build mesh in
  bedfile=$3 # relative filepath to .dat file to mesh

  # Get the start and end of the bedfile
  start=$(awk 'NR==1 {print $1}' $bedfile)
    end=$(awk 'END {print $1}'   $bedfile)

  # Given the desired dx, back calculate the numer of x gridcells
  Nx=$(awk -v L=$end -v dx=$dx 'BEGIN {OFMT = "%.0f"; print (L/dx)}')

  # update .grd file with parsed parameters
  update_grd $start $end $Nx $grdfile

  # Actually the mesh.*
  ElmerGrid 1 2 tmp.grd -out $out_dir -autoclean

  # Remove tmp.grd file to reduce clutter
  rm tmp.grd
}

#-----------------------------
# Edit these paths as need be:
#-----------------------------
grdfile="./src/mesh.grd"                    # path to grd file
bedfile="./input_data/lilk-a_bed.dat"       # path to bed file

if [ $# -eq 0 ]; then
  usage
fi

# Run the main loop
main $1 $2 $3
