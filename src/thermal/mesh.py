#!/usr/bin/env python3

import re
import os
import subprocess
import numpy as np

def update_template(L, Nx, Ny, grd_fp='../mesh.grd'):
    """Update a templare ElmerGird file with approraite mesh parameters.

    Inputs:
        grd_dp    (str)   --> file path to reference ElmerGrid (.grd) file
        L         (float) --> Length of the domain (meters)
        Nx        (int)   --> Number of horizontal mesh nodes
        Ny        (int)   --> Number of vertical mesh nodes
    Outputs:
        ElmerGrid (str)   --> ElmerGrid with approraite parameters values
    """

    # Open the template file, and read into memory as a single string
    with open(grd_fp) as f:
        template = f.read()

    # Replace the length of the domain
    ElmerGrid = re.sub('Subcell Limits 1 =.*',
                       'Subcell Limits 1 = {:e} {:e}'.format(0, L),
                       template)

    # Replace the length of the domain
    ElmerGrid = re.sub('Element Divisions 1 =.*',
                       'Element Divisions 1 = {:d}'.format(Nx),
                       ElmerGrid)

    # Replace the length of the domain
    ElmerGrid = re.sub('Element Divisions 2 =.*',
                       'Element Divisions 2 = {:d}'.format(Ny),
                       ElmerGrid)

    return ElmerGrid


# make mesh function

key = 'sprg'
dx  = 200
out_fp = f"initialization/coarse/result/{key}/mesh_dx{dx}"

# Make sure we are in top directory of repo
#os.chdir('../..')

if not os.path.exists(out_fp):
    os.makedirs(os.path.join('../..', out_fp))


# find the x coordinate of the approraite bed file
x = np.loadtxt(f'../../input_data/{key}_bed.dat')[:,0]
L  = x.max() - x.min()
Nx = int(round(L/dx))

# create the glacier specific mesh file
ElmerGird = update_template(L, Nx, 15)


# write the ElmerGrid string to a temporary file
with open("./../../tmp.grd", 'w+') as f:
    f.write(ElmerGird)

cmd   = "cd ../.. &&"
cmd  += "ElmerGrid 1 2 tmp.grd -out {} -autoclean".format(out_fp)

result = subprocess.run(
                        [cmd],
                        shell=True,
                        capture_output=True,
                        text=True
                        )

if "ElmerGrid: command not found" in result.stderr:
    raise OSError('Elmer/Ice executables are not available. Make sure you are '\
                  'in the correct computing environment')
