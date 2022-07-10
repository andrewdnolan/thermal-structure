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

def animation_2D_field(ds, anim_fn, title):
    fig, anim = animate_Zs.animate_2D_field(ds,
                                            field=["vel_m", "temperature"],
                                            title=title,
                                            cbar_label=['Velocity [m a$^{-1}$]',
                                                        'Temperature [C]'],
                                            vmin=[0, -12.5],
                                            vmax=[25, 0.0],
                                            cmap=['viridis', 'plasma'],
                                            stride=1,
                                            interval=100)

    plt.close()

    anim.save(anim_fn)

limit_type='surf' # 'surf' or 'rho'


for Delta_mb in ['-1', '-1.5']:
    for T_ma in ['-9.02', '-7.02']:
    # if Delta_mb == "-1": continue

        src_fp = f"glc1-a/nc/glc1-a_dx_50_NT_250_dt_1.0_MB_{Delta_mb}_OFF_Tma_{T_ma}_limit_rho_Cfirn_0.05_prog.nc"

        ds = test_open.dataset(src_fp)

        out_fn = f'./figs/glc1-a_MB_{Delta_mb}_OFF_250a_Tma_{T_ma}.mp4'
        title  = "$\Delta \dot b$ = {:.3f}, $T_m$ = {:.2f}".format(float(Delta_mb), float(T_ma))

        animation_2D_field(ds, out_fn, title)

        print('mp4 file succesfully written to:')
        print(f'\t {out_fn}')

    # fig, ax = plt.subplots()
    #
    # for t in ds.t:
    #
    #     dTdz = np.gradient(ds.enthalpy_h.sel(t=t).isel(coord_1=100) / 1e3,
    #                        ds.Z.sel(t=t).isel(coord_1=100),
    #                        edge_order=2)
    #
    #     ax.plot(dTdz, ds.coord_2 / ds.coord_2.size)
    #
    #
    # plt.show()


# # def animate_vertical_gradient():
#
# from matplotlib import animation
#
#
# i  = 100              # x-coordinate index to plot
# NT = ds.t.size        # Number of timestep
# NZ = ds.coord_2.size  # Number of verical nodes
#
# # create empty array for storing enthalpy gradient
# dEdz = np.zeros((NT,NZ))
#
# # calculate the vertical enthalpy gradient for all timesteps
# for j, t in enumerate(ds.t):
#     dEdz[j,:] = np.gradient(ds.enthalpy_h.sel(t=t).isel(coord_1=i) / 1e3,
#                             ds.Z.sel(t=t).isel(coord_1=i),
#                             edge_order=2)
#
# fig, ax = plt.subplots(figsize=(3,6), dpi=300, constrained_layout=True)
#
#
# line1, = ax.plot([], [], lw=2, color='tab:blue')
# line   = [line1]
#
# ax.set_xlim(dEdz.min(),dEdz.max())
# ax.set_ylim(0 , 1)
#
# ax.set_ylabel(r'$\zeta = z/H [-]$')
# ax.set_xlabel(r'$\frac{\partial E}{\partial z} [\rm{kJ} \, \rm{kg}^{-1} \, \rm{m}^{-1}]$')
#
# # Function to be animated
# def animate(i):
#     # plot the free surface
#     line[0].set_data(dEdz[i,:],
#                      ds.coord_2 / ds.coord_2.size)
#     string = r'$\frac{\partial E}{\partial z}('
#     line[0].set_label(string+'t={:.1f})$'.format(float(ds.t.isel(t=i))))
#
#     ax.legend(loc=3)
#     return line
#
# anim = animation.FuncAnimation(fig, animate,
#                                frames=np.arange(0, NT),
#                                interval=50,
#                                blit=True)
#
#
# anim.save(f'./figs/test.mp4')
#
#
# fig, ax = plt.subplots()
#
# for t in ds.t:
#
#     dTdz = np.gradient(ds.enthalpy_h.sel(t=t).isel(coord_1=100) / 1e3,
#                        ds.Z.sel(t=t).isel(coord_1=100),
#                        edge_order=2)
#
#     ax.plot(dTdz, ds.coord_2 / ds.coord_2.size)
#
#     out_fn = f'./figs/glc1-a_MB_{Delta_mb}_OFF_100a_hightempres.mp4'
#     title  = "$\Delta \dot b$ = {:.3f}".format(float(Delta_mb))
#
#     animation_2D_field(ds, out_fn, title)
#
#     print('mp4 file succesfully written to:')
#     print(f'\t {out_fn}')

# fig, anim = animate_Zs.animate_2D_field(ds,
#                                         field=["vel_m", "temperature"],
#                                         title="$\Delta \dot b$ = {:.3f}, $\Delta t$ = {:.1f}".format(Delta_mb, dt),
#                                         cbar_label=['Velocity [m a$^{-1}$]', 'Temperature [C]'],
#                                         vmin=[0, -12.5],
#                                         vmax=[25, 0.0],
#                                         cmap=['viridis', 'plasma'],
#                                         stride=1,
#                                         interval=100)
#
# plt.close()
#
# anim.save(f'./figs/Delta_b{Delta_mb}_dt{dt}_anim_netcdfout.mp4')
