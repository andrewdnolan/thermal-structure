# Uncertainty Quantification    

Parametric sensitivity tests.  

## One-off Test's

```bash
./sensitivity.py --key "crmpt12" \
                 -t_f 3000 \
                 -dx 50\ 
                 -dt 0.1 
                 -Dynamic_int 10 
                 -off -0.37 \
                 -T_ma -8.5 -C_firn 0.0500 -f_dd 0.006 -w_en 1.0 -w_aq 5.0 -IC 0.0
```

## "Production" Runs 

When the goal is to fully sample the parameter space (of a single parameter as time), there are helper scripts to automatically generate `SLURM` job submission scripts. 
Update the parameter dictionaries found in the `./params/` folder. 

To generate the model run commands and SLURM job array submission script, run: 
```bash
./make_jobscript.py --key 'crmpt12'
```

On a SLURM cluster (e.g. `cedar`):
```bash
sbatch ./run/crmpt12/parametric_sensitivity_submit.sh
```

### Bug in making `runname` for reference parameter runs 

A bug in the name `make_runname` function in `sensitivity.sh` script results
in all reference parameter runs having the runname
`crmpt12_dx_50_NT_30000_dt_0.1_1aTST__`. So as a stopgap solution run the following
```bash
source sensitivity.sh

for var in C_firn f_dd w_en w_aq IC; do 

  # get the reference value for current parameter
  parse_reference $var
  runname="crmpt12_dx_50_NT_30000_dt_0.1_1aTST_${var}_${reference}"
    
  for dir in "gridded" "thinned"; do
    pushd result/crmpt12/${dir}
    echo cp crmpt12_dx_50_NT_30000_dt_0.1_1aTST__.zarr.tar ${runname}.zarr.tar
    popd
  done
done
```
