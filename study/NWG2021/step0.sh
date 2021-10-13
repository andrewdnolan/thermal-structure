#!/usr/bin/env bash

for accum in $(seq 2000 50 2450); do

  if [[ $accum == 2100 ]]; then
    continue
  fi

  # replace variable with correct value in template sif
  sed "s#<A_mean>#${accum}#g;" \
      "sifs/step0.freesurfacerelax.sif" > \
      "relax_${accum}.sif"

  # execute the sif file
  ElmerSolver "relax_${accum}.sif"

  # Convert result files into NetCDFs
  ../../src/elmer2nc/elmer2nc.sh  -r "./mesh/relax_${accum}.result" \
                                  -m "./mesh/" \
                                  -t 1000 \
                                  -o "./nc/"

  rm "relax_${accum}.sif"
done
