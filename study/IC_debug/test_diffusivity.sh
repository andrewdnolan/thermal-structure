#!/usr/bin/env bash

NT=50

# "constant"
for diff in "variable_fixed_diff"; do

  echo -e "------------------------------------------------------------------------\n" \
          "Running test case for ${diff} temperate diffusivity\n"\
          "------------------------------------------------------------------------\n" \
          > "test.out" 2>&1


  if   [[ "$diff" == "constant" ]]; then

    sed "s#<run>#"$diff"#g;
         s#<diffusion>#real \#9.92e-4*spy#g" "diffusivity.sif" > "${diff}.sif"

  elif [[ "$diff" == "variable" ]] || [[ "$diff" == "variable_fixed_diff" ]]; then
    string="Variable Temperature \n \t\t\t Real procedure \"../../bin/thermodynamics\" \"Diffusivity\" "

    sed "s#<run>#"$diff"#g;
         s#<diffusion>#$string#g" "diffusivity.sif" > "${diff}.sif"
  fi


  #  Run the model
  ElmerSolver "${diff}.sif" 2>>"test.out" | tee "logs/${diff}.log"


  # Convert from .result to .nc
  ../../src/elmer2nc/elmer2nc.sh  -r "./mesh/${diff}.result" \
                                  -m "./mesh/" \
                                  -t "${NT}" \
                                  -o "./nc/" >> "test.out" 2>&1

  rm "${diff}.sif"
done
