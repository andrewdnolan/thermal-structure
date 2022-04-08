#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run.sh:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set +x

#-------------------------------------------------------------------------------
# Run the recovery from the surge
#-------------------------------------------------------------------------------

dx=100                           # time step size
KEY=klun-b
RESTART="klun-b_2000a_dt_0.5_dx_100_MB_-0.20_OFF_spline_k2.result"


# Update the .SIF FILE with the model run specifc params
sed "s#<DX>#"$dx"#g;
     s#<KEY>#"$KEY"#g;
     s#<RESTART>#"$RESTART"#g" "./sifs/diagnostic.sif" > "./sifs/SS.sif"

# Run the model
ElmerSolver "./sifs/SS.sif"

# Convert result files into NetCDFs
../../src/elmer2nc/elmer2nc.sh -r "./${KEY}/mesh_dx50/SS_cup_test.result" \
                               -m "./${KEY}/mesh_dx50/" \
                               -t 1 \
                               -o "./${KEY}/nc/"

# Remove the sif file
rm "./sifs/SS.sif"
