#!/usr/bin/env bash

run(){

    wrk_dir=shared_directory/Thesis/thermal-structure/expr/00_CoupledInit/

    docker exec elmerenv /bin/sh -c "cd ${wrk_dir}; ElmerSolver ./sifs/${run_name}.sif"

    python3 ../../src/thermal/add_attr.py -f "./sifs/${run_name}.sif" \
                                          -a "sif" \
                                            "result/${KEY}/nc/${run_name}.nc"

    grid_data.py  -i "result/${KEY}/${run_name}.nc" \
                  -o "result/${KEY}/${run_name}.zarr"
}

KEY='crmpt12'

run_names=(crmpt12_dx_50_NT_30000_dt_0.1_MB_-0.36_OFF_Tma_-8.5_prog_fdd_2.0_NoStrainHeating
           crmpt12_dx_50_NT_30000_dt_0.1_MB_-0.36_OFF_Tma_-8.5_prog_fdd_2.0_WithStrainHeating
           crmpt12_dx_50_NT_30000_dt_0.1_MB_-0.50_OFF_Tma_-9.0_prog_NoStrainHeating
           )

for run_name in ${run_names[*]};do

    run

done

