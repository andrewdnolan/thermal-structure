#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run.sh:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set +x

#-------------------------------------------------------------------------------
# Run the recovery from the surge
#-------------------------------------------------------------------------------

dt=0.11111111                           # time step size
NT=4500                                 # number of time step
TT=$((NT*dt))                           # total time of simulation

# Update the .SIF FILE with the model run specifc params
sed "s#<dt>#"$dt"#g;
     s#<NT>#"$NT"#g "./sifs/SS_cupped_netbalance.sif" > "./sifs/SS.sif"

# Run the model
ElmerSolver "./sifs/SS.sif"

# Convert result files into NetCDFs
../../src/elmer2nc/elmer2nc.sh -r "./crmpt18-b/mesh_dx50/SS_cup_test.result" \
                               -m "./crmpt18-b/mesh_dx50/" \
                               -t $NT \
                               -o "./crmpt18-b/nc/"

# Remove the sif file
rm "./sifs/SS.sif"
