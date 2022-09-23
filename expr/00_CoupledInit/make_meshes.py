#!/usr/bin/env python3

import os
import sys
sys.path.append('../../src')
import thermal.mesh as meshing

# input glaciers and desired mesh resolutions for coarse grid searches
glaciers = { # small
             'crmpt12'   : 50,  # Δx [m]
             # 'crmpt18-a' : 50,  # Δx [m]
             # 'crmpt18-b' : 50,  # Δx [m]
             'glc1-a'    : 50,  # Δx [m]
             # 'glc1-b'    : 50,  # Δx [m]
             # # medium
             # 'lilk-a'    : 100, # Δx [m]
             # 'lilk-b'    : 100, # Δx [m]
             # 'klun-a'    : 100, # Δx [m]
             # 'klun-b'    : 100, # Δx [m]
             'sprg'      : 100, # Δx [m]
             # # large
             # 'fish'      : 200, # Δx [m]
             # 'klut-a'    : 200, # Δx [m]
             # 'klut-b'    : 200, # Δx [m]
             'twds-b'    : 200, # Δx [m]
            }


# filepath relative to top of git repo
rel_fp = "./expr/00_CoupledInit/result"

for key in glaciers.keys():
    Δx     = glaciers[key]
    out_fp = os.path.join(rel_fp, key, "mesh_dx{}".format(Δx))
    nc_fp  = f'result/{key}/nc/'
    gridded_fp = f'result/{key}/gridded/'

    #-----
    # Mesh
    #-----
    # Make the mesh
    stat = meshing.make_mesh(key , Δx, out_fp)

    #-------
    # NetCDF
    #-------
    for nc_folder in [nc_fp, gridded_fp]:
        # Make the NetCDF folder
        if not os.path.exists(nc_folder):
            os.mkdir(nc_folder)

        # Add .gitkeep file to NetCDF folder so dir struc is preserved on GitHub
        if not os.path.exists(os.path.join(nc_folder, '.gitkeep')):
            with open(os.path.join(nc_folder, '.gitkeep'), 'w') as f:
                pass

    #--------
    # figures
    #--------
    # Make figures dir
    if not os.path.exists(f'figs/'):
        os.mkdir('figs/')
    if not os.path.exists(f'figs/{key}/'):
        os.mkdir(f'figs/{key}/')


    # Add .gitkeep file to figures folder so dir struc is preserved on GitHub
    if not os.path.exists(os.path.join(f'figs/{key}/', '.gitkeep')):
        with open(os.path.join(f'figs/{key}/', '.gitkeep'), 'w') as f:
            pass


    #-----
    # logs
    #-----
    # Make logs dir
    if not os.path.exists(f'logs/'):
        os.mkdir('logs/')
    if not os.path.exists(f'logs/{key}/'):
        os.mkdir(f'logs/{key}/')


    # Add .gitkeep file to logs folder so dir struc is preserved on GitHub
    if not os.path.exists(os.path.join(f'logs/{key}/', '.gitkeep')):
        with open(os.path.join(f'logs/{key}/', '.gitkeep'), 'w') as f:
            pass

    #-----
    # run
    #-----
    # Make run dir
    if not os.path.exists(f'run/'):
        os.mkdir('run/')

    # Add .gitkeep file to run folder so dir struc is preserved on GitHub
    if not os.path.exists(os.path.join(f'run/', '.gitkeep')):
        with open(os.path.join(f'run/', '.gitkeep'), 'w') as f:
            pass
