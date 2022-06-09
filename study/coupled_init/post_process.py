#!/usr/bin/env python3

import sys
import numpy as np
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.append('../../src/thermal')
sys.path.append('../../src/plotting')
import animate_Zs
import open as test_open


def make_colorbar(array):
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    #---------------------------------------------------------------------------
    # For Seting up the colorbar:
    #    - http://csc.ucdavis.edu/~jmahoney/matplotlib-tips.html
    #---------------------------------------------------------------------------
    cmap = cm.plasma
    norm = mcolors.Normalize(vmin=np.min(array),
                             vmax=np.max(array))

    s_map = cm.ScalarMappable(norm=norm, cmap=cmap)
    s_map.set_array(np.linspace(np.min(array), np.max(array), len(array)+1))

    # If color parameters is a linspace, we can set boundaries in this way
    halfdist = (array[1] - array[0]) / 2.0
    bounds   = np.linspace(np.min(array)  - halfdist,
                           np.max(array)  + halfdist,
                           len(array) + 1)

    return cmap, norm, s_map, bounds


# Set global matplotlib style parameters
plt.rcParams.update({'text.usetex': True,
                     'animation.html': 'jshtml',
                     'figure.facecolor': 'w',
                     'savefig.bbox':'tight'})

Delta_mb = -1.0
diag_fp = "glc1-a/nc/glc1-a_dx_50_MB_-1.000_OFF_diag.nc"
src_fp = "glc1-a/nc/glc1-a_dx_50_NT_20000_dt_0.1_MB_-1.000_OFF_prog_NetcdfOutPutSolve.nc"


with xr.open_dataset(diag_fp) as diag:
    # diag["height"] = xr.where(diag.height <= 10, 0, diag.height)
    diag["Z"]     = diag.zbed + diag.Z * diag.height
    diag['vel_m'] = np.sqrt(diag['velocity 1']**2 + diag['velocity 2']**2)

fig, ax = plt.subplots(1,1)

idx = 125
T = diag['temperature'].isel(t=-1, coord_1=idx)

dTdz = np.gradient(diag.temperature.isel(t=-1, coord_1=idx),
                   diag.Z.isel(t=-1, coord_1=idx),
                   edge_order=2)

# ax.plot(dTdz, ds.depth.sel(t=t).isel(coord_1=100))
ax.plot(dTdz, diag.coord_2 / diag.coord_2.size)

# ax.invert_yaxis()
plt.show()

# ds = test_open.dataset(src_fp)
#
# fig, ax = plt.subplots(1,1)
#
# for t in ds.t[5:500:10]:
#
#     T = ds['enthalpy_h'].sel(t=t).isel(coord_1=100)
#
#     dTdz = np.gradient(ds.enthalpy_h.sel(t=t).isel(coord_1=100),
#                        ds.Z.sel(t=t).isel(coord_1=100),
#                        edge_order=2)
#
#     # ax.plot(dTdz, ds.depth.sel(t=t).isel(coord_1=100))
#     ax.plot(dTdz, ds.coord_2 / ds.coord_2.size)
#
# # ax.invert_yaxis()
# plt.show()
#
#
# cmap, norm, s_map, bounds = make_colorbar(ds.t[5::10])
#
#
# fig, ax = plt.subplots(1,1)
#
# for t in ds.t[5::10]:
#
#     color = cmap(norm(t))
#
#     T = ds['enthalpy_h'].sel(t=t).isel(coord_2=-1)
#
#     # dTdz = np.gradient(ds.enthalpy_h.sel(t=t).isel(coord_1=100),
#     #                    ds.Z.sel(t=t).isel(coord_1=100),
#     #                    edge_order=2)
#
#     # ax.plot(dTdz, ds.depth.sel(t=t).isel(coord_1=100))
#     ax.plot(ds.X.isel(coord_2=-1), T, color=color)
#
# cbar = fig.colorbar(s_map,
#                     spacing='proportional',
#                     ticks=np.linspace(ds.t[5::10].values.min(),
#                                       ds.t[5::10].values.max(),
#                                       10),
#                     ax=ax,
#                     boundaries=bounds,
#                     drawedges=True,
#                     format='%2.{}f'.format(3),
#                     aspect=35)
# # ax.invert_yaxis()
# plt.show()
#
# test = [np.gradient(ds.temperature.sel(coord_1=i).isel(t=0).values, ds.Z.sel(coord_1=i).isel(t=0), edge_order=2) for i in ds.coord_1]


# fig, anim = animate_Zs.animate_2D_field(ds.isel(t=slice(0,210)),
#                                         field=["vel_m", "temperature"],
#                                         title="$\Delta \dot b$ = {:.3f}".format(Delta_mb),
#                                         cbar_label=['Velocity [m a$^{-1}$]', 'Temperature [C]'],
#                                         vmin=[0, -12.5],
#                                         vmax=[25, 0.0],
#                                         cmap=['viridis', 'plasma'],
#                                         stride=1,
#                                         interval=300)
#
# plt.close()
#
# anim.save(f'./figs/Delta_b{Delta_mb}_anim_netcdfout.mp4')


# with xr.open_dataset(src_fp) as src:
#     src["height"] = xr.where(src.height <= 10, 0, src.height)
#     src["Z"]      = src.zbed + src.Z * src.height
#     src['vel_m']  = np.sqrt(src['velocity 1']**2 + src['velocity 2']**2)
#
#
# Vol = src.height.isel(coord_2=-1).integrate("coord_1") /\
#       src.height.isel(coord_2=-1).isel(t=0).integrate("coord_1")
#
#
# # plt.plot(src.coord_1,
# #          src['surface_enthalpy'].isel(t=1).isel(coord_2=-1), 'x-')
#
#
# plt.pcolormesh(src.t,
#          src.coord_1/1e3,
#          src['mass balance'].isel(coord_2=0).T,#.differentiate("t").T,
#          shading='auto')
#
# plt.colorbar()
# plt.show()
