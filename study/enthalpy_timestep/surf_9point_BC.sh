#!/usr/bin/env bash

T=50  #[a]

# "discrete" 
for case in "mean"; do


  # Place appropriate variables in
  sed "s#<case>#"$case"#g" "./src/surf_9point_BC.sif" > "surf_9point_BC_${case}.sif"

  # Run the model
  ElmerSolver "surf_9point_BC_${case}.sif"

  # Convert from .result to .nc
  ../../src/elmer2nc/elmer2nc.sh  -r "./mesh/surf_9point_BC_${case}.result" \
                                  -m "./mesh/" \
                                  -t 450 \
                                  -o "./nc/"

  rm "surf_9point_BC_${case}.sif"
done
