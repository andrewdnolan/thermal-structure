# Uncertainty Quantification    

Parametric sensitivity tests.  

To generate the model run commands and SLURM job array submission script, run: 
```bash
./make_jobscript.py --key 'crmpt12
```

On a SLURM cluster (e.g. `cedar`):
```bash
sbatch ./run/crmpt12/parametric_sensitivity_submit.sh
```