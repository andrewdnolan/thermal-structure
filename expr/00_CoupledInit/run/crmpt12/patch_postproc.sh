

parse_params()
{
  IFS=" " read -r offset T_ma \
    <<< $(sed -n "${1}p" ./run/crmpt12/gridsearch.commands | cut -d " " -f 13,15)
}

mkdir -p "sub_from_0.5_to_3000_by_10_YEARS/nc/"
mkdir -p "result/crmpt12/sub_from_0.5_to_3000_by_10_YEARS/zarr/"

parse_params 10


# parameter dictionary for griding NetCDFs
param_dict="{\"T_ma\"   : ${T_ma},
             \"offset\" : ${offset}}"

# from the parameters fill out the run_name
run_name="crmpt12_dx_50_NT_30000_dt_0.1_MB_${offset}_OFF_Tma_${T_ma}_prog"

# run the subsampling script
downsample -i "./result/crmpt12/nc/${run_name}.nc" \
           -o "result/crmpt12/sub_from_0.5_to_3000_by_10_YEARS/nc/${run_name}.nc" \ 
           --start 4 --stop -1 --stride 100

# grid the NetCDF file written by the NetcdfUGRIDOutputSolver
python3 ../../src/thermal/grid_data.py "result/crmpt12/sub_from_0.5_to_3000_by_10_YEARS/nc/${run_name}.nc" \
                               -out_fn "result/crmpt12/sub_from_0.5_to_3000_by_10_YEARS/zarr/${run_name}.zarr" \
                               -params "${param_dict}"
echo $run_name
