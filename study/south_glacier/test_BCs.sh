#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test_BCs.sh:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set +x

source ../../src/utils/initialization.sh

diagnostic_run()
{
  # Update the .SIF FILE with the model run specifc params
  sed "s#<DX>#"$dx"#g;
       s#<FIT>#"$FIT"#g;
       s#<KEY>#"$KEY"#g;
       s#<ELA>#"$ELA"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<run_name>#"$run_name"#g;
       s#<SS_itters>#"$SS_itters"#g;
       s#<heat_source>#"$heat_source"#g;
       s#<smart_heater>#"$smart_heater"#g;
       s#<enthalpy_SBC>#$enthalpy_SBC#g;
       s#<Source_at_Surface>#"$Source_at_Surface"#g;
       s#<temperate_headwall>#$temp_headwall#g;" "./sifs/diagnostic.sif" > "./sifs/${run_name}.sif"


  docker exec elmerenv /bin/sh -c "cd shared_directory/Thesis/thermal-structure/study/south_glacier;
                                   ElmerSolver ./sifs/${run_name}.sif "

  # # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

prognostic_run()
{
  # Update the .SIF FILE with the model run specifc params
  sed "s#<DX>#"$dx"#g;
       s#<dt>#"$dt"#g;
       s#<NT>#"$NT"#g;
       s#<KEY>#"$KEY"#g;
       s#<FIT>#"$FIT"#g;
       s#<ELA>#"$ELA"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<RESTART>#"$RESTART"#g
       s#<run_name>#"$run_name"#g;
       s#<SS_itters>#"$SS_itters"#g;
       s#<Dynamic_int>#"$Dynamic_int"#g
       s#<heat_source>#"$heat_source"#g;
       s#<smart_heater>#"$smart_heater"#g;
       s#<enthalpy_SBC>#$enthalpy_SBC#g;
       s#<Source_at_Surface>#"$Source_at_Surface"#g;
       s#<temperate_headwall>#$temp_headwall#g;" "./sifs/prognostic.sif" > "./sifs/${run_name}.sif"

  docker exec elmerenv /bin/sh -c "cd shared_directory/Thesis/thermal-structure/study/south_glacier;
                                  ElmerSolver ./sifs/${run_name}.sif "
  # # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

get_enthalpy_SBC()
{

  # echo the variable so it's all one line and works with sed
  # ref: https://stackoverflow.com/a/42356201/10221482
  if [[ $1 == "DSBC" ]]; then
    tmp="\t Enthalpy_h = Equals \"Surface_Enthalpy\""
    enthalpy_SBC=$(echo $tmp)
  else
    tmp="
    \t External Enthalpy_h = Equals \"Surface_Enthalpy\" \n
    \t Heat Transfer Coefficient  = Variable \"Enthalpy_h\", \"Phase Change Enthalpy\", \"Enthalpy Heat Diffusivity\", \"Enthalpy Water Diffusivity\" \n
    \t\t\t
    Real LUA \"IfThenElse((tx[0]<tx[1]), tx[2]/0.1, tx[3]/0.1)\""
    enthalpy_SBC=$(echo $tmp)
  fi

}

get_heatsource_node()
{
  if [[ $1 == "1bfs" ]]; then
    Source_at_Surface="False"
  else
    Source_at_Surface="True"
  fi
}

get_temperate_headwall()
{
  if [[ $1 == "True" ]]; then
    # echo "well shit"
    temp_headwall="\t Enthalpy_h = Real 134404.00"
  else
    temp_headwall=""
  fi
}
get_run_name(){

  if [[ $temperate_headwall == "True" ]];then
    head="TempHead"
  else
    head="ColdHead"
  fi

  if [[ $smart_heater == "True" ]];then
    heater="_SH"
  else
    heater=""
  fi


  # provide a prefix to which parameter values are attached too
  run_name="${1}_${SBC}_vhs_${node}_${head}${heater}"
}

# parse the parameters from the json files
parse_json "params/glc1-a.json"
# let's allow for X number of steady state itterations
SS_itters=25
# set the glacier key
KEY='glc1-fetal2011'
T_ma=-7.0
offset=-1.5
ELA=0
dt=1.0
T_f=250
# Number of time interval based on dt
NT=$(awk -v dt=$dt -v T_f=$T_f 'BEGIN {OFMT = "%.0f"; print (T_f/dt)}')
# Execute interval for dynamics solvers,
Dynamic_int=$(awk -v dt=$dt 'BEGIN {OFMT = "%.0f"; print (1.0/dt)}')

# test various (s)urface (b)oundary (c)onditions
for SBC in "RSBC" "DSBC"; do

  # call fucntion which return var $enthalpy_SBC
  get_enthalpy_SBC $SBC

  # test different vertical nodes to prescribe latent heat source
  for node in "1bfs" "atfs"; do
    # (1) (b)elow (f)ree (s)urface
    # (at)        (f)ree (s)urface
    #-----------------------------

    # call fucntion which returns var $Source_at_Surface
    get_heatsource_node $node

    # test if forcing the headwall to temperate has any effect
    for temperate_headwall in True False; do

      # call fucntion which returns var $temperate_headwall
      get_temperate_headwall $temperate_headwall

      # try using the elmer smart heater keyword
      for smart_heater in True False; do
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        # diagnostic run
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        # diagnostic model run_name
        get_run_name "${KEY}_dx_${dx}_MB_${offset}_Tma_${T_ma}"
        # run_name="test_diag"
        # run the model for a parameter combination
        diagnostic_run $dx $KEY $offset $run_name $SS_itters

        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        # prognostic run
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        # diagnostic run is now restart variable
        RESTART="${run_name}.result"
        # don't need as many SS iterations for transient runs
        SS_itters=25
        # prognostic model run name
        get_run_name "${KEY}_dx_${dx}_Tf_${T_f}_dt_${dt}_MB_${offset}_Tma_${T_ma}"

        # run_name="test_prog"

        # run the model for a parameter combination
        prognostic_run $dx $KEY $offset $run_name $SS_itters $restart $NT $dt

      done
    done
  done
done
