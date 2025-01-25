#!/usr/bin/env bash

# script for building python environment on cedar

# Load python 3
module load python/3.12.4

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
pip install --no-index -r dev-environment.txt

# install the local pyton module (thermal)
pip install --editable $HOME/project/thermal-structure/src/thermal/ --no-deps

# Steps three and four form https://docs.alliancecan.ca/wiki/Advanced_Jupyter_configuration
echo -e '#!/bin/bash\nunset XDG_RUNTIME_DIR\njupyter lab --ip $(hostname -f) --no-browser' \
        > $VIRTUAL_ENV/bin/jupyterlab.sh
chmod u+x $VIRTUAL_ENV/bin/jupyterlab.sh

deactivate
