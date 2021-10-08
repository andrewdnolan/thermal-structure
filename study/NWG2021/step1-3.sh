#!/usr/bin/env bash

src_sif="./sifs/step1-3.heatsources.sif"

run_simulation () {
  ElmerSolver "step${step}.sif" 2>>"test.out" | tee "logs/step${step}.log"

  n_itter=$(grep "Steady state iteration"  "logs/step${step}.log" | wc -l)

  echo -e "Step finsihed in: \t ${n_itter} itterations" >> "test.out" 2>&1

  # Convert result files into NetCDFs
  ../../src/elmer2nc/elmer2nc.sh  -r "./mesh/heat_sources_step${step}.result" \
                                  -m "./mesh/" \
                                  -t 2 \
                                  -o "./nc/" >> "test.out" 2>&1
  rm "step${step}.sif"
}

for step in 1 2; do

  if [[ $step == 1 ]]; then
    echo -e "------------------------------------------------------------------------\n" \
            "Step 1: Dirichlet B.C. for surface, no heatflux at bed\n"\
            "------------------------------------------------------------------------\n" \
            > "test.out" 2>&1

    sed "s#<step>#step1#g;" $src_sif > "step${step}.sif"

  elif [[ $step == 2 ]]; then
    echo -e "------------------------------------------------------------------------\n" \
            "Step 2: Dirichlet B.C. for surface and geothermal heat flux at bed \n"\
            "------------------------------------------------------------------------\n" \
            >> "test.out" 2>&1

    sed "s#<step>#step2#g;
         s#Enthalpy Heat Flux BC =.*#Enthalpy Heat Flux BC = Logical True#g;"\
           $src_sif > "step${step}.sif"
  fi


  run_simulation
done
