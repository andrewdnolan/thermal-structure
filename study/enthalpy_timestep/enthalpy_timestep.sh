#!/usr/bin/env bash

T=200  #[a]

for dt_m in 1 2 6 12; do
  # Total number of timesteps
  NT=$(((T*12)/dt_m))
  dt_a=$(awk -v dt=$dt_m 'BEGIN {print dt/12.0}')
  NT_a=$((12/dt_m))

  # Place appropriate variables in
  sed "s#<NT>#"$NT"#g;
       s#<dt_m>#"$dt_m"#g
       s#<H_dt>#"$dt_a"#g;
       s#<NT_a>#"$NT_a"#g" "enthalpy_timestep.sif" > "20a_dt_${dt_m}m.sif"

  # Run the model
  ElmerSolver "20a_dt_${dt_m}m.sif"

  # Convert from .result to .nc
  ../../src/elmer2nc/elmer2nc.sh  -r "./mesh/20a_dt_${dt_m}m.result" \
                                  -m "./mesh/" \
                                  -t "${NT}" \
                                  -o "./nc/"

  rm "20a_dt_${dt_m}m.sif"
done
