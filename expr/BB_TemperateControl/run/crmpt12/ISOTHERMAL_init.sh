#!/usr/bin/env bash

########################################
# Isothermal initialization gridsearch #
########################################

dx=50         # horizontal mesh resolution [m]
dt=1.0        # annual timesteps [a]
NT=3000       # 3ka initialization
t_f=$NT        # b/c dt=1.0, NT==t_f
key="crmpt12" # glacier identifier
T_ma=0.0      # Dummy value for the air temperature

for offset in $(seq -w -0.5 0.01 -0.35); do

    # create parameter json (dictionary) for gridding NetCDF files
    param_dict="{\"offset\" : ${offset}}"

    # run_name from the coupled experiments
    run_name="${key}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_prog"
    # file name for the gridded results, which denotes the isotheraml temperature
    gridded_fn="${key}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_ISOTHERM_prog"

    # run the isothermal simulation for a given mass balance offset
    ./initialize.py -dx $dx --key $key \
                    -t_f $t_f -dt $dt -Dynamic_int 1 \
                    -off $offset  -T_ma $T_ma

    # grid and write the data as netcdfs
    # isothermal results have problem with .zarr file format for some reason
    grid_data.py -i "result/crmpt12/nc/${run_name}.nc" \
                 -o "result/crmpt12/gridded/${gridded_fn}.nc" \
                 -p "${param_dict}"
done 