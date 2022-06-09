#!/usr/bin/env python3

import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

sys.path.append('../../src/thermal')
sys.path.append('../../src/plotting')
import animate_Zs
import open as test_open


# Set global matplotlib style parameters
plt.rcParams.update({'text.usetex': True,
                     'animation.html': 'jshtml',
                     'figure.facecolor': 'w',
                     'savefig.bbox':'tight'})



def calc_perc_temp(src):
    # number of valid nodes above minimum ice-thickness
    # NN = src.height.where(src.height < 10).count(dim=('nMesh_node'))#.count(dim=('coord_1', 'coord_2'))
    NN = src.height.where(src.height < 10).count(dim=('coord_1', 'coord_2'))

    # mask for determining where nodes are temperate
    mask = src['enthalpy_h'] > src['phase change enthalpy']
    # number of nodes where ethalpy exceeds phase chnage enthalpy
    NT = src.enthalpy_h.where(mask).count(dim=('coord_1', 'coord_2'))
    return NT/NN * 100


Delta_mb = -1.000
dt       = 1.0
src_fp = "glc1-a/nc/glc1-a_dx_50_NT_100_dt_1.0_MB_-1.000_OFF_prog_local_NetcdfOutPutSolve.nc"

ds = test_open.dataset(src_fp)


perc_temp = calc_perc_temp(ds)


plt.plot(perc_temp.t, perc_temp)
plt.show()


fig, anim = animate_Zs.animate_2D_field(ds,
                                        field=["vel_m", "temperature"],
                                        title="$\Delta \dot b$ = {:.3f}, $\Delta t$ = {:.1f}".format(Delta_mb, dt),
                                        cbar_label=['Velocity [m a$^{-1}$]', 'Temperature [C]'],
                                        vmin=[0, -12.5],
                                        vmax=[25, 0.0],
                                        cmap=['viridis', 'plasma'],
                                        stride=1,
                                        interval=100)

plt.close()

anim.save(f'./figs/Delta_b{Delta_mb}_dt{dt}_anim_netcdfout.mp4')
