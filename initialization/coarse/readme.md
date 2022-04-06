# Initialization

This folder contains run scripts, and the base directory structure for storing results and figures from uncoupled initialization experiments.
Here we seek to find the appropriate mass balance offset to the KMR curve which produces a glacier of the same size as the initial condition (i.e. observed).  

---

The basic workflow for how these results are as follows:

First, execute
```bash
python3 make_meshes.py
```
in an environment with `Elmer` executables available (i.e. within the docker container).
This script will create the meshes and appropriate directory structure for each glacier, based on the `dx` parameter in the respective glaciers `.json` file (in the `params` folder).


Then on either server side or locally `make` runscripts by executing:
```bash
make
```
which creates the `.in`/`.out` files and submission scripts for each size class, by calling the `prepare2submit` script. To change the runtime or memory allocation either change the `makefile` or run the `prepare2submit` script directly:
```bash
./prepare2submit --force            \ # overwrite .in/.out files if already exist
                 --mem    "2250M"   \ # mem  request per job in array
                 --time   "9:30:00" \ # time request per job in array
                 --group  "medium"  \ # glacier size class
                 --stride 20          # num jobs in the SLURM array to run at once
```
If either of the previous two commands were run locally, then the changes need to transferred onto cedar (i.e. commit locally and pull on cedar).


Once you're on cedar navigate back to this directory. To run the initialization
experiment for one size class (in this case medium) execute:
```bash
# must be executed on cedar
sbatch run/medium_init_spline.sh
```
which queues our job submission scripts. Edit appropriate (`prepare2submit`) to include your email to receive information about when jobs begin and end.

To get results back on your local computer, run:
```bash
fp2src="/home/user/scratch/thermal-structure/initialization/coarse/result"
rsync -avP user@cedar.computecanada.ca:${fp2src}/${glac_key}/nc/*.nc  result/${glac_key}/nc/
```
where `${glac_key}` should be replaced with whatever glacier you'd like to pull results for. (Easy enough to write for loop over glacier keys).

Finally, if you'd like to visualize the spin-up results run:
```
python3 ../../src/plotting/plot_spinup.py -p "params/${glac_key}.json"
```
where `${glac_key}` is replaced with the glacier which you want to plot. This script will produce both relative volume and final free surface plots, named based on the model run params in the `fig/${glac_key}` directory.

<!-- ## `prepare2submit.sh`:

- This script create necessary files to submit `SLURM` job array on Cedar.
- Submission scripts are group by size since runtime and memory
  usage should be similar within the groups.

__Example usage__:
```bash
./prepare2submit --force            \ # overwrite .in/.out files if already exist
                 --mem    "2250M"   \ # mem  request per job in array
                 --time   "9:30:00" \ # time request per job in array
                 --group  "medium"  \ # glacier size class
                 --stride 20          # num jobs in the SLURM array to run at once
```

## Plotting
For quick visualization of the of the model results use:
```bash
python3 ../../src/plotting/plot_spinup.py -p "params/${glac_key}.json"
```
where `${glac_key}` is replaced with the glacier which you want to plot. -->

<!-- ```
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

``` -->


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
