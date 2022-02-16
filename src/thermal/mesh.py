#!/usr/bin/env python3

import re
import os
import glob
import subprocess
import numpy as np

# Deal with relative path problems when using a module
thermal_fp = os.path.dirname(__file__)
git_top    = os.path.join(thermal_fp, '../..')

def update_template(L, Nx, Ny, grd_fp=os.path.join(thermal_fp, '../mesh.grd')):
    """Update a templare ElmerGird file with approraite mesh parameters.

    Inputs:
        grd_dp    (str)   --> file path to reference ElmerGrid (.grd) file
        L         (float) --> Length of the domain (meters)
        Nx        (int)   --> Number of horizontal mesh nodes
        Ny        (int)   --> Number of vertical mesh nodes
    Outputs:
        ElmerGrid (str)   --> ElmerGrid with approraite parameters values
    """

    # Make sure we can find the template file
    if not os.path.exists(grd_fp):
        raise FileNotFoundError(f'\n Trouble finding .grd file: {grd_fp} \n')

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



def make_mesh(key, dx, out_fp, force=False):
    """Make the mesh directory. Currently only suports first-order
       quadrilateral elments (Elmer 404 nodes)

    Inputs:
        key    (str)  --> glacier identifier (see input_data repo)
        dx     (int)  --> horzontal mesh resolution
        out_fp (str)  --> file path to directory to place the mesh.* files
        force  (bool) --> whether to overwrite mesh if it already exists

    return:
        result (...)  --> instance of subprocess.CompletedProcess class
                          has attributes stderr and stdout which record the
                          stderr and stdout respectively from the mesh generation
                          call. If you mesh isn't made correctly check these attrs
    """

    # Really only need to check about one directory up.....
    # ElmerGrid can create the right most directory (i.e. .../mesh_dx{dx})
    git_top_fp = '/'.join(out_fp.split('/')[:-1])
    if not os.path.exists(os.path.join(git_top, git_top_fp)):
        os.makedirs(os.path.join(git_top, git_top_fp))

    # check if the mesh files alrady exists.  no need to create the mesh if it
    # already exists To overwrite exiting mesh set force == True
    if ( os.path.exists(os.path.join(git_top, out_fp)) ) and \
       ( len(glob.glob( os.path.join(git_top, out_fp, 'mesh.*'))) == 4 ) and \
       ( not force ):
       return


    # find the x coordinate of the approraite bed file
    x = np.loadtxt(os.path.join(git_top,f'input_data/topography/{key}_bed.dat'))[:,0]
    L  = x.max() - x.min()
    Nx = int(round(L/dx))

    # create the glacier specific mesh file
    ElmerGird = update_template(L, Nx, 15)

    # write the ElmerGrid string to a temporary file
    with open(os.path.join(git_top,"tmp.grd"), 'w+') as f:
        f.write(ElmerGird)

    pwd = os.getcwd()
    os.chdir(git_top)

    cmd   = "ElmerGrid 1 2 tmp.grd -out {} -autoclean".format(out_fp)
    cmd  += " && rm tmp.grd"
    result = subprocess.run([cmd],
                            shell=True,
                            capture_output=True,
                            text=True
                            )
    os.chdir(pwd)

    if "ElmerGrid: command not found" in result.stderr:
        raise OSError('Elmer/Ice executables are not available. Make sure you are '\
                      'in the correct computing environment')

    return result
