#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample_run.sh:
#   - Mass balance grid-search to find steady state positions for 3
#     test glaciers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set +x

dt=1                                    # time step size
NT=1000                                 # number of time step
TT=$((NT*dt))                           # total time of simulation

SIF='./simple_spinup.sif'               # template SIF file
BED='farinotti_corrected'               # Mesh DB for the given bed config
KEY='glc1-a'


for OFFSET in $(seq -w -2.0 0.1 -0.1);do
  # Model RUN identifier
  RUN="${KEY}_${TT}a_dt_${dt}_dx_100_MB_${OFFSET}_OFF"

  # File paths to input data
  Zb_fp="../../input_data/${KEY}_bed.dat"
  Zs_fp="../../input_data/${KEY}_surf.dat"

  # Update the .SIF FILE with the model run specifc params
  sed "s#<NT>#"$NT"#g;
       s#<RUN>#"$RUN"#g;
       s#<KEY>#"$KEY"#g;
       s#<Zs_fp>#"$Zs_fp"#g;
       s#<Zb_fp>#"$Zb_fp"#g;
       s#<OFFSET>#"$OFFSET"#g" "$SIF" > "./sifs/${RUN}.sif"

  ElmerSolver "./sifs/${RUN}.sif"

 # Convert result files into NetCDFs
 ../../src/elmer2nc/elmer2nc.sh  -r "./result/${KEY}/${RUN}.result" \
                                 -m "./result/${KEY}/" \
                                 -t $NT \
                                 -o "./result/nc/"

  # Remove the sif file
  rm "./sifs/${RUN}.sif"
done
