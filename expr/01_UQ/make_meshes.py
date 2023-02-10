#!/usr/bin/env python3

import os
import sys
# sys.path.append('../../src')
# import thermal.mesh as meshing

# all the sensitivty test are done with crmpt12
glacier = 'crmpt12'

# instead of glaciers, we have a dictionary of different sensitivty test
experiments = {
                'crmpt12' : [50,  15], # Δx, Nz [m]
              }

# filepath relative to top of git repo
rel_fp = "./expr/01_UQ/result"

#
for exp in experiments:
    Δx, Nz = experiments[exp]
    out_fp = os.path.join(rel_fp, exp, "mesh_dx{}".format(Δx, Nz))
    nc_fp  = f'result/{exp}/nc/'
    gridded_fp = f'result/{exp}/gridded/'

    # -----
    # Mesh
    # -----
    if type(Nz) is list:
        for nz  in Nz:
            # append to the mesh name
            out_fp += f"_nz{nz}"
            # make the mesh
            stat = meshing.make_mesh(glacier , Δx, out_fp, Nz=nz)

    else:
        stat = meshing.make_mesh(glacier , Δx, out_fp, Nz=Nz)

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
    if not os.path.exists(f'figs/{exp}/'):
        os.mkdir(f'figs/{exp}/')


    # Add .gitkeep file to figures folder so dir struc is preserved on GitHub
    if not os.path.exists(os.path.join(f'figs/{exp}/', '.gitkeep')):
        with open(os.path.join(f'figs/{exp}/', '.gitkeep'), 'w') as f:
            pass


    #-----
    # logs
    #-----
    # Make logs dir
    if not os.path.exists(f'logs/'):
        os.mkdir('logs/')
    if not os.path.exists(f'logs/{exp}/'):
        os.mkdir(f'logs/{exp}/')


    # Add .gitkeep file to logs folder so dir struc is preserved on GitHub
    if not os.path.exists(os.path.join(f'logs/{exp}/', '.gitkeep')):
        with open(os.path.join(f'logs/{exp}/', '.gitkeep'), 'w') as f:
            pass

    #-----
    # run
    #-----
    # Make run dir
    if not os.path.exists(f'run/'):
        os.mkdir('run/')
    if not os.path.exists(f'run/{exp}/'):
        os.mkdir(f'run/{exp}/')

    # Add .gitkeep file to run folder so dir struc is preserved on GitHub
    if not os.path.exists(os.path.join(f'run/{exp}/', '.gitkeep')):
        with open(os.path.join(f'run/{exp}/', '.gitkeep'), 'w') as f:
            pass
