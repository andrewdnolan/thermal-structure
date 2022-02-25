#!/usr/bin/env python3

import os
import sys
sys.path.append('../../src')
import thermal.mesh as meshing

# input glaciers and desired mesh resolutions for coarse grid searches
glaciers = { # small
             'crmpt12'   : 50,  # Δx [m]
             'crmpt18-a' : 50,  # Δx [m]
             'crmpt18-b' : 50,  # Δx [m]
             'glc1-a'    : 50,  # Δx [m]
             'glc1-b'    : 50,  # Δx [m]
             # medium
             'lilk-a'    : 100, # Δx [m]
             'lilk-b'    : 100, # Δx [m]
             'klun-a'    : 100, # Δx [m]
             'klun-b'    : 100, # Δx [m]
             'sprg'      : 100, # Δx [m]
             # large
             # 'fish'      : 500, # Δx [m]
             # 'klut-b'    : 500, # Δx [m]
             # 'twds-b'    : 250, # Δx [m]
            }


# filepath relative to top of git repo
rel_fp = "./initialization/coarse/result"

for key in glaciers.keys():
    Δx     = glaciers[key]
    out_fp = os.path.join(rel_fp, key, "mesh_dx{}".format(Δx))
    nc_fp  = f'result/{key}/nc/'

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
