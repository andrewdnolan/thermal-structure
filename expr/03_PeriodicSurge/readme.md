# Periodic Surges

## Pseudo surge simulation 
The basis of the this experiment folder is the `periodic_surge.py` script, which is a wrapper for the functions in the `periodic_surge.sh` file. An example usage of the script would be: 
```bash
./periodic_surge.py -k "crmpt12" -SP 2 -QP 28 -beta 1.000e-03 -TT 6000
```

## Postprocessing 
For analysis and plotting, it's convenient to extract the timeseries of interest from the full `zarr` files. 
To do so, run: 
```bash 
python3 run/crmpt12/postprocess_timeseries.py \
    --gridded_dir "/Volumes/thermal/Thesis/thermal-structure/expr/03_PeriodicSurge/result/crmpt12/gridded/"
```
