#!/usr/bin/env bash

# script for building python environment on cedar

# Load python 3
module load python/3.8.10
# Lazy scipy stack
module load scipy-stack

if [[ -d "$HOME/py4elmer" ]]; then
  rm -r "$HOME/py4elmer"
fi

# create virtualenv
virtualenv --no-download $HOME/py4elmer
# activate env
source $HOME/py4elmer/bin/activate
# update pip
pip install --no-index --upgrade pip
# install the various packages needed ontop of scipy stack
pip install --no-index xarray jupyterlab dask dask_jobqueue bokeh seaborn

# Steps three and four form https://docs.alliancecan.ca/wiki/Advanced_Jupyter_configuration
echo -e '#!/bin/bash\nunset XDG_RUNTIME_DIR\njupyter lab --ip $(hostname -f) --no-browser' \
        > $VIRTUAL_ENV/bin/jupyterlab.sh
chmod u+x $VIRTUAL_ENV/bin/jupyterlab.sh

deactivate
