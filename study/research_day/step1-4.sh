#!/usr/bin/env bash

src_sif="./sifs/steps1-4.sif"

run_simulation () {
  ElmerSolver "step${step}.sif" 2>>"./logs/step1-4_T_${temp}.log" | tee "./logs/step${step}_T_${temp}.log"

  n_itter=$(grep "Steady state iteration"  "logs/step${step}_T_${temp}.log" | wc -l)

  echo -e "Step finsihed in: \t ${n_itter} itterations" >> "./logs/step1-4_T_${temp}.log" 2>&1

  # Convert result files into NetCDFs
  ../../src/elmer2nc/elmer2nc.sh  -r "./glc1-a/mesh_dx50/heat_sources_step${step}_T_${temp}.result" \
                                  -m "./glc1-a/mesh_dx50/" \
                                  -t $n_itter \
                                  -o "./glc1-a/nc/" >> "./logs/step1-4.log" 2>&1
  rm "step${step}.sif"
}
for temp in -9.02 -7.02; do
  for step in 1 2 3 4; do
    if [[ $step == 1 ]]; then
      echo -e "------------------------------------------------------------------------\n" \
              "Step 1: Dirichlet B.C. for surface, no heatflux at bed\n"\
              "------------------------------------------------------------------------\n" \
              > "./logs/step1-4.log" 2>&1

      sed "s#<temp>#${temp}#g;
           s#<step>#step${step}#g;" $src_sif > "step${step}.sif"

    elif [[ $step == 2 ]]; then
      echo -e "------------------------------------------------------------------------\n" \
              "Step 2: Dirichlet B.C. for surface and geothermal heat flux at bed \n"\
              "------------------------------------------------------------------------\n" \
              >> "./logs/step1-4.log" 2>&1

      sed "s#<temp>#${temp}#g;
           s#<step>#step${step}#g;
           s#Enthalpy Heat Flux BC =.*#Enthalpy Heat Flux BC = Logical True#g;"\
             $src_sif > "step${step}.sif"

    elif [[ $step == 3 ]]; then
      echo -e "------------------------------------------------------------------------\n" \
              "Step 3: Dirichlet B.C. for surface and geothermal heat flux at bed \n"\
              "        and Deformational heating. \n"\
              "------------------------------------------------------------------------\n" \
              >> "./logs/step1-4.log" 2>&1

      sed "s#<temp>#${temp}#g;
           s#<step>#step${step}#g;
           s#Enthalpy Heat Flux BC =.*#Enthalpy Heat Flux BC = Logical True#g;
           s#Exec Solver .* = String \"Never\" ! <DeformationalHeat>#Exec Solver = String \"Always\"#g;"\
           $src_sif > "step${step}.sif"

    elif [[ $step == 4 ]]; then
      echo -e "------------------------------------------------------------------------\n" \
              "Step 4: Dirichlet B.C. for surface and geothermal heat flux at bed \n"\
              "        Deformational heating, and latnet heat from meltwater refreezing.\n"\
              "------------------------------------------------------------------------\n" \
              >> "./logs/step1-4.log" 2>&1

      sed "s#<temp>#${temp}#g;
           s#<step>#step${step}#g;
           s#Enthalpy Heat Flux BC =.*#Enthalpy Heat Flux BC = Logical True#g;
           s#Exec Solver .* = String \"Never\" ! <DeformationalHeat>#Exec Solver = String \"Always\"#g;
           s#! <Latent_Heat_Source># Heat Source = Variable Q_lat, Densi \n\t\t\tReal lua \"(tx[0] / tx[1])\"#g;"\
           $src_sif > "step${step}.sif"
    fi

    # Actually run the model
    # run_simulation
  done
done
