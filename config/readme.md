# `config`

In this folder there are a some useful configuration files, for managing the various software packages required.  

---

### `dev-environment.txt`
This is the list of python packages needed for running all workflows.
The txt file can be used to create a `conda` environment with all the necessary
`python3` packages for post processing model results.
To install, assuming you have `miniforge` installed, run:
```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -y -n thermal python=3.12.4
conda install -n thermal --file dev-environment.txt

conda activate thermal

cd ../src/thermal/
pip install --editable .
```

### `modulefile.cc.cedar`
This file loads all the necessary `modules` on `cedar` needed for model execution. To use, run the command:
```bash
source modulefile.cc.cedar
```
All the automatically generated submission scripts, in the `expr/` folder, will load this file as their first step.


### `build_py4elmer.sh`
Will build the `python` virtual env needed for post processing on `cedar`. 
Unfortunately, `conda` env aren't available on `cedar` but this build script
uses the same `dev-environment.txt` file so we can guarantee consistent packages
(but not necessarily versions) between `cedar` and our `conda` env. To install, run:
```bash
bash ./build_py4elmer.sh
```

### `link2scratch
Will create the symbolic links and build out the repo structure on the `scratch` filesystem.
This is a workaround to have the same repo structure on `scratch` and `project`
or any other location, without needing to maintain the `.git` repo on `scratch`. 
```
make -f config/link2scratch
```
