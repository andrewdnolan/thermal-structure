# `thermal-structure`

This repository contains code used to investigate how glacier surges alter glacier thermal structure. 
In pursuit of our research objectives, we've employed a thermomechanically coupled numerical ice-flow model, `Elmer/Ice` in a 2D ($x-z$) configuration.
We conduct a suite of experiments ranging from initialization through periodic perturbations to the  basal slip coefficient, mean to mimic observed surges. 
This repository contains code used pre-process the necessary input data, execute and post-process the model runs, and do some plotting and analysis. 


https://user-images.githubusercontent.com/32367657/227734840-fa76a693-07a7-4c8d-a49b-55b8ac3761b4.mp4


### Project layout: 

The repository structure is as follows: 

```
├── Makefile    # Global makefile used to compile code throughout the repository
├── README.md  
├── bin         # Folder to hold the compiled code (mainly Elmer user functions)
├── config      # Repository wide configuration files 
├── expr        # Various experiments we conducted using `Elmer/Ice`
├── include     # External code used in our analysis
├── input_data  # Surface mass balance and topography data need to run `Elmer` simulations
├── notebooks   # Collection of jupyter notebooks doing data pre-processing and analysis
├── src         # source code folder, containing `Elmer` boundary condition code and more
└── study
```

### Getting started
To get started first clone the repository by running, 
```bash
$ git clone git@github.com:andrewdnolan/thermal-structure.git
```
Then, navigate into the top directory (i.e. `cd thermal-strucure`) and compile/install the code by running 
```bash
$ make
```
which will execute the global Makefile that compiles external `FORTRAN` code from the `include` folder, compiles the `Elmer` user functions (`USF`s) from the `src/elmer_src` folder, and builds the legacy `NetCDF` post-processing `FORTRAN` program contained in the `src/elmer2nc` folder.

Finally, you'll also need to run
```
pip install --editable src/thermal/
```
which will create a locally editable copy of the `python` post-processing library `thermal`. `thermal` is built around `xarray` and `dask` to enable efficient, out of memory, and parallel postprocessing of the terabytes of data produced in the various experiments. 

### Installing `Elmer/Ice`  

If working from a Linux machine, follow the [compilation instructions](https://elmerice.elmerfem.org/wiki/doku.php?id=compilation:compilationcmake) from the `Elmer/Ice` documentation. Compiling `Elmer/Ice` on non-Linux machines (e.g. `OSX` and Windows) is notoriously challenging. To circumvent this problem, I have written a [`Docker` container](https://hub.docker.com/r/andrewdnolan/elmerice). Instructions on how to install the `Docker` can be found in the associated `readme`. All within code/files in this repository assume the Elmer executables (e.g. `ElmerSolver`, `ElmerGird`) are available. 


### Note on reproducibility  

Given the computational cost of thermomechanically coupled numerical modeling, most of the model execution has been done on high performance computing (HPC) resources (e.g. Compute Canada's `cedar` cluster). Therefore, many of the files various `expr` directories are not meant to be run on a local laptop/workstation. That being said, command line scripts to run one-off simulations, which still could take < 24 hours, are available in each of the `expr` directories. 

<!-- __To Do__:
  - Set up additional solvers, so that at each times step we have a record of the
    amount of heat contributed by each source term in the governing equation.
    - Diffusive Flux:
      - http://www.nic.funet.fi/pub/sci/physics/elmer/doc/ElmerModelsManual.pdf#page=229
      - How do I pass a variable diffusivity to the solver??
        - I can easily write another solver, but is there a way with the existing
          elmer variables?

  - Write a solver to calculate the peclet number as field variable for each timestep.
    - Should also write the "Brinkman Number" (see [Meyer and Minchew, 2018](https://www-sciencedirect-com.proxy.lib.sfu.ca/science/article/pii/S0012821X18303790?via%3Dihub#se0080) for example of it being used.)


Notes from NWG:
  - from GEF during card ride:
    - we need to quantify how (and if) the changes in surge vigor during periodic surges are results of a less temperate area along the bed, or difference in driving stress -->
