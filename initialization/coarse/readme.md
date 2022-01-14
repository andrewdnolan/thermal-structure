# Initialization

Here follows run scripts, results, and figures from uncoupled initialization experiments.

## `prepare2submit.sh`:
  - This script create necessary files to submit `SLURM` job array on cedar.
  - Geometries and submission scripts are group by size since runtime and memory
    usage should be similar within the groups.

__TO DO__:
  - [x] Folder structure to create hold the 10+ `.sif` files generated per glacier.
  - [x] Can we set this up so that the `.sif` files are only created as needed:
      1. To not clutter up the repo
      2. I think HPC platforms slow down when writing ASCII files

  - [x] function to create the `Inputs.txt` file

  - [x] function to create the `Outputs.txt` file:
    - [x] call `elmer2nc.sh` to conver form `.result` to `NetCDF`

  - Actual submission script:
    - [ ] count number of jobs per group in existing nested `for` loops
    - [ ] create group specific run name
    - [ ] Need to pass three runtimes and three memory request one for each size class
    - [ ] check that paths work if submitted from  `run` directory 
    - [ ] Make the scripts executable?


  - [ ] Can we plot, but only after all the runs are completed for a single glacier?
    - maybe this needs to be it's own script which waits for the job arrays to complete?
