# Periodic Surges

## Pseudo surge simulation 
The basis of the this experiment folder is the `periodic_surge.py` script, which is a wrapper for the functions in the `periodic_surge.sh` file. An example usage of the script would be: 
```bash
./periodic_surge.py -k "crmpt12" -SP 2 -QP 28 -beta 1.000e-03 -TT 6000
```

## Split simulations 
For computational and/or numerical stability reasons it's often useful to split the 9ka simulations into
more managable 3ka chunks. The `./periodic_surge.py` script supports restart files, but the `$run_name` variable
is unaware of the restarts / time splitting. So, it's neccessary to rename the resulting files to contain
more descriptive information about the corresponding time span run. So for example, if you've run the first
3ka of a 9ka simulation, you should run: 
```bash
for file in $(find result/crmpt12/ -name "crmpt12_dx_50_TT_3000.0_MB*"); do  
    mv $file ${file/3000.0/0--3ka}
done
```

## Postprocessing 
For analysis and plotting, it's convenient to extract the timeseries of interest from the full `zarr` files. 
To do so, run: 
```bash 
python3 run/crmpt12/postprocess_timeseries.py \
    --gridded_dir "/Volumes/thermal/Thesis/thermal-structure/expr/03_PeriodicSurge/result/crmpt12/gridded/"
```
