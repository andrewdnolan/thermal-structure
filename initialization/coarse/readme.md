# Initialization

Here follows run scripts, results, and figures from uncoupled initialization experiments.

```
.
├── figs
  ├──${glac_key}/
├── params               # folder w/ .json param files for each glacier
  ├──${glac_key}.json    # json dictionary with param values
├── result
  ├──${glac_key}  
    ├──nc/
    ├──mesh_dx${DX}/
└── sifs
    ├──simple_spinup.sif #
├── make_meshes.py       # 
├── makefile
├── prepare2submit       #
├── readme.md
├── run                  #
├── sample_run.sh        #

```

## `prepare2submit.sh`:

- This script create necessary files to submit `SLURM` job array on cedar.
- Geometries and submission scripts are group by size since runtime and memory
  usage should be similar within the groups.

__Example usage__:
```
./prepare2submit --force          \ # overwrite .in/.out files if already exist
                 --mem "2250M"    \ # mem  request per job in array
                 --time "9:30:00" \ # time request per job in array
                 --group "medium" \ # glacier size class
                 --stride 20        # num jobs in the SLURM array to run at once
```


<!-- __TO DO__:
  - [x] Folder structure to create hold the 10+ `.sif` files generated per glacier.
  - [x] Can we set this up so that the `.sif` files are only created as needed:
      1. To not clutter up the repo
      2. I think HPC platforms slow down when writing ASCII files

  - [x] function to create the `Inputs.txt` file

  - [x] function to create the `Outputs.txt` file:
    - [x] call `elmer2nc.sh` to conver form `.result` to `NetCDF`

  - Actual submission script:
    - [x] count number of jobs per group in existing nested `for` loops
    - [x] create group specific run name
    - [x] parse `MEM` and `RUN_TIME` from command line
    - [x] wrap whole process into function to be run per size group
    - [x] parse group from command line and appropriate control flow.
    - [ ] check that paths work if submitted from  `run` directory
    - [x] Make the scripts executable?

  - Make sense to run script once per group.
    - [ ] `makefile` with rule for each group size

  - [ ] Can we plot, but only after all the runs are completed for a single glacier?
    - maybe this needs to be it's own script which waits for the job arrays to complete? -->
