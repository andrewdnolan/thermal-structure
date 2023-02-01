#!/usr/bin/env bash

# script for building python environment on cedar

# Load python 3
module load python/3.10.2
# Lazy scipy stack
module load scipy-stack

if [[ -d "$HOME/python_envs/py4elmer" ]]; then
  rm -r "$HOME/python_envs/py4elmer"
fi

# create virtualenv
virtualenv --no-download $HOME/python_envs/py4elmer
# activate env
source $HOME/python_envs/py4elmer/bin/activate
# update pip
pip install --no-index --upgrade pip
# install the various packages needed ontop of scipy stack
pip install --no-index xarray jupyterlab dask_jobqueue bokeh seaborn zarr dask==2023.1.0 distributed==2023.1.0

# install the local pyton module (thermal)
pip install --editable ../src/thermal/

# Steps three and four form https://docs.alliancecan.ca/wiki/Advanced_Jupyter_configuration
echo -e '#!/bin/bash\nunset XDG_RUNTIME_DIR\njupyter lab --ip $(hostname -f) --no-browser' \
        > $VIRTUAL_ENV/bin/jupyterlab.sh
chmod u+x $VIRTUAL_ENV/bin/jupyterlab.sh

deactivate
