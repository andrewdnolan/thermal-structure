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