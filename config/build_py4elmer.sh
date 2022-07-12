#!/usr/bin/env bash

# script for building python environment on cedar

# Load python 3
module load python/3.8.10
# Lazy scipy stack
module load scipy-stack

if [[ -d "$HOME/py4elmer_new" ]]; then
  rm -r "$HOME/py4elmer_new"
fi

# create virtualenv
virtualenv --no-download $HOME/py4elmer_new
# activate env
source $HOME/py4elmer_new/bin/activate
# update pip
pip install --no-index --upgrade pip
# install the various packages needed ontop of scipy stack
pip install xarray jupyterlab dask dask_jobqueue seaborn --no-index

deactivate
