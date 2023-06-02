# Numerical Experiments  

## Experiment Descriptions: 
### `00_CoupledInit`
  - Steady state initialization of the synthetic geometries. Contains workflow to conduct a grid search across surface (i.e., air temperature and mass balance) forcing. 

### `01_UQ` 
  - Prognostic one at a time parametric sensitivity tests. Contains workflow to conduct grid search across defined parametric ranges, using a reference glacier identified through the previous experiment (`00_CoupledInit`). 

### `02_surge2steady`
  - Prognostic simulation where the reference glacier is subject to a pseudo surge, of varying magnitudes, and the allowed to fully recover from the perturbation. 

### `03_PeriodicSurge`
  - 

### `AA_SolverTiming`
  - 

### `BB_TemperateControl`
  - 

# Workflow on Compute Canada Computing Resources 

