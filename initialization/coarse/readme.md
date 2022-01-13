# Initialization

Here follows run scripts, results, and figures from uncoupled initialization experiments.

## `prepare2submit.sh`:
  - This script create necessary files to submit `SLURM` job array on cedar.
  - Geometries and submission scripts are group by size since runtime and memory
    usage should be similar within the groups.

  __TO DO__:
    - [ ] Folder structure to create hold the 10+ `.sif` files generated per glacier.
      - For the past project these were created all at once before the job was submitted.
        Can we set this up so that the `.sif` files are only created as needed:
          1. To not clutter up the repo
          2. I think HPC platforms slow down when writing ASCII files

    - [ ] function to create the `Inputs.txt` file
      - [ ] At least need to include file path to the generated `.sif` file.
       - Maybe this should instead be a call to `bash` function to create this file.
    - [ ] function to create the `Outputs.txt` file:
      - [ ] call `elmer2nc.sh` to conver form `.result` to `NetCDF`
      - [ ] Can we plot, but only after all the runs are completed for a single glacier?
