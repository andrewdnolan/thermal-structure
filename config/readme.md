# `config`

In this folder there are a some useful configuration files, for managing the various software packages required.  

---

### `thermal.yml`
This is the configuration file needed for creating the `conda` environment with all the necessary `python3` packages for post processing model results. To install, assuming you have `conda` or `miniconda` installed, run:
```bash
conda env create -f thermal.yml
```

### `modulefile.cc.cedar`
This file loads all the necessary `modules` on `cedar` needed for model execution. To use, run the command:
```bash
source modulefile.cc.cedar
```
All the automatically generated submission scripts, in the `expr/` folder, will load this file as their first step.


### `build_py4elmer.sh`
Will build the `python` virtual env needed for post processing on `cedar`. Unfortunately, `conda` env aren't available so we need a seperate package management approach for `cedar`. This isn't a vert elegant solution, and whenever the `thermal.yml` file is update the appropriate packages need to be added to this script. To run,
```bash
bash ./build_py4elmer.sh
```
