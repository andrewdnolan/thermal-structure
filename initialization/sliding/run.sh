#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample_run.sh:
#   - Mass balance grid-search to find steady state positions for 3
#     test glaciers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set +x

KEYS=('crmpt12' 'crmpt18-a' 'glc1-a' 'lilk-a'
      'klun-a' 'sprg' 'fish' 'klut-a' 'twds-a')
RES=(50 50 50 100 100 100 200 200 200)
#-------------------------------------------------------------------------------
# Numerical parameters
#-------------------------------------------------------------------------------
BETA=1.0                                # slip Coefficient, i.e. no sliding
#-------------------------------------------------------------------------------
# input data parameters
#-------------------------------------------------------------------------------
SIF='./sifs/diagnostic_sliding.sif'          # template SIF file
#KEY='sprg'                           # glacier key for input data

for i in "${!KEYS[@]}"; do
  # i-th resolution and glacier key
  dx="${RES[i]}"
  KEY="${KEYS[i]}"

  # Model RUN identifier
  RUN="${KEY}_dx_${dx}_Beta_${BETA}"

  # File paths to input data
  Zb_fp="../../input_data/${KEY}_bed.dat"
  Zs_fp="../../input_data/${KEY}_surf.dat"

  # Update the .SIF FILE with the model run specifc params
  sed "s#<DX>#"$dx"#g;
       s#<RUN>#"$RUN"#g;
       s#<KEY>#"$KEY"#g;
       s#<BETA>#"$BETA"#g;
       s#<Zs_fp>#"$Zs_fp"#g;
       s#<Zb_fp>#"$Zb_fp"#g" "$SIF" > "./sifs/${RUN}.sif"

 # Run the model
 ElmerSolver "./sifs/${RUN}.sif" | tee "result/${KEY}/logs/${RUN}.log"


 # Convert result files into NetCDFs
 ../../src/elmer2nc/elmer2nc.sh  -r "./result/${KEY}/mesh_dx${dx}/${RUN}.result" \
                                 -m "./result/${KEY}/mesh_dx${dx}/" \
                                 -t 2 \
                                 -o "./result/${KEY}/nc/"

  # Remove the sif file
  rm "./sifs/${RUN}.sif"
done
