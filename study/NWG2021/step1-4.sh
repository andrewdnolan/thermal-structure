#!/usr/bin/env bash

src_sif="./sifs/step1-3.heatsources.sif"

run_simulation () {
  ElmerSolver "step${step}.sif" 2>>"step1-4.log" | tee "logs/step${step}.log"

  n_itter=$(grep "Steady state iteration"  "logs/step${step}.log" | wc -l)

  echo -e "Step finsihed in: \t ${n_itter} itterations" >> "step1-4.log" 2>&1

  # Convert result files into NetCDFs
  ../../src/elmer2nc/elmer2nc.sh  -r "./mesh/heat_sources_step${step}.result" \
                                  -m "./mesh/" \
                                  -t 2 \
                                  -o "./nc/" >> "step1-4.log" 2>&1
  rm "step${step}.sif"
}

for step in 4; do
#for step in 1 2 3 4; do

  if [[ $step == 1 ]]; then
    echo -e "------------------------------------------------------------------------\n" \
            "Step 1: Dirichlet B.C. for surface, no heatflux at bed\n"\
            "------------------------------------------------------------------------\n" \
            > "step1-4.log" 2>&1

    sed "s#<step>#step${step}#g;" $src_sif > "step${step}.sif"

  elif [[ $step == 2 ]]; then
    echo -e "------------------------------------------------------------------------\n" \
            "Step 2: Dirichlet B.C. for surface and geothermal heat flux at bed \n"\
            "------------------------------------------------------------------------\n" \
            >> "step1-4.log" 2>&1

    sed "s#<step>#step${step}#g;
         s#Enthalpy Heat Flux BC =.*#Enthalpy Heat Flux BC = Logical True#g;"\
           $src_sif > "step${step}.sif"

  elif [[ $step == 3 ]]; then
    echo -e "------------------------------------------------------------------------\n" \
            "Step 3: Dirichlet B.C. for surface and geothermal heat flux at bed \n"\
            "        and Deformational heating. \n"\
            "------------------------------------------------------------------------\n" \
            >> "step1-4.log" 2>&1

    sed "s#<step>#step${step}#g;
         s#Enthalpy Heat Flux BC =.*#Enthalpy Heat Flux BC = Logical True#g;
         s#Exec Solver .* = String \"Never\" ! <DeformationalHeat>#Exec Solver = String \"Always\"#g;"\
         $src_sif > "step${step}.sif"

  fi

  if [[ $step == 4 ]]; then
    ElmerSolver "sifs/step4.heatsources.sif" 2>>"step1-4.log" | tee "logs/step${step}.log"

    #n_itter=$(grep "Steady state iteration"  "logs/step${step}.log" | wc -l)

    echo -e "Step ran for: \t 90 timesteps" >> "step1-4.log" 2>&1

    # Convert result files into NetCDFs
    ../../src/elmer2nc/elmer2nc.sh  -r "./mesh/heat_sources_step${step}.result" \
                                    -m "./mesh/" \
                                    -t 90 \
                                    -o "./nc/" >> "step1-4.log" 2>&1
  else
    run_simulation
  fi
done
