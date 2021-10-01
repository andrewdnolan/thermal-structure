#!/usr/bin/env bash

T=50  #[a]

# echo "dt_m runtime" |
# awk -v OFS='\t' '{print $1 "\t" $2}' > out.enthalpy_timestep

for dt_m in 1 2 6 12; do
  # Total number of timesteps
  NT=$(((T*12)/dt_m))
  dt_a=$(awk -v dt=$dt_m 'BEGIN {print dt/12.0}')
  NT_a=$((12/dt_m))

  # Place appropriate variables in
  sed "s#<T>#"${T}"#g;
       s#<NT>#"$NT"#g;
       s#<dt_m>#"$dt_m"#g
       s#<H_dt>#"$dt_a"#g;
       s#<NT_a>#"$NT_a"#g" "enthalpy_timestep.sif" > "${T}a_dt_${dt_m}m.sif"

  # # Start the timer
  # start=$(date +%s.%N)
  #
  # # Run the model
  # ElmerSolver "${T}a_dt_${dt_m}m.sif"
  #
  # # End the timer
  # end=$(date +%s.%N)
  #
  # # Execution time of the solver
  # runtime=$(awk -v start=$start -v end=$end 'BEGIN {print end - start}')
  #
  # echo "${dt_m} ${runtime}" |
  # awk -v OFS='\t' '{print $1 "\t" $2}' >> out.enthalpy_timestep

  # # Convert from .result to .nc
  # ../../src/elmer2nc/elmer2nc.sh  -r "./mesh/${T}a_dt_${dt_m}m.result" \
  #                                 -m "./mesh/" \
  #                                 -t "${NT}" \
  #                                 -o "./nc/"

  # rm "${T}a_dt_${dt_m}m.sif"
done
