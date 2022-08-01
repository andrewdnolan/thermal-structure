#!/usr/bin/env python3

import os
import sys
sys.path.append('../../src')
import thermal.mesh as meshing

# input glaciers and desired mesh resolutions for coarse grid searches
glaciers = { # small
             'glc1-a'    : 50,  # Δx [m]
        }


# filepath relative to top of git repo
rel_fp = "./study/test_streched"

for key in glaciers:
    Δx     = glaciers[key]
    out_fp = os.path.join(rel_fp, key, "mesh_dx{}".format(Δx))
    nc_fp  = f'{key}/nc/'

    #-----
    # Mesh
    #-----
    # Make the mesh
    stat = meshing.make_mesh(key , Δx, out_fp)

    #-------
    # NetCDF
    #-------
    # Make the NetCDF folder
    if not os.path.exists(nc_fp):
        os.mkdir(nc_fp)

    # Add .gitkeep file to NetCDF folder so dir struc is preserved on GitHub
    if not os.path.exists(os.path.join(nc_fp, '.gitkeep')):
        with open(os.path.join(nc_fp, '.gitkeep'), 'w') as f:
            pass

    #--------
    # figures
    #--------
    # Make figures dir
    if not os.path.exists(f'figs/{key}/'):
        os.mkdir(f'figs/{key}/')

    # Add .gitkeep file to figures folder so dir struc is preserved on GitHub
    if not os.path.exists(os.path.join(f'figs/{key}/', '.gitkeep')):
        with open(os.path.join(f'figs/{key}/', '.gitkeep'), 'w') as f:
            pass
