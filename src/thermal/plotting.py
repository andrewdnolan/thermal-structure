#!/usr/bin/env python3

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.animation import FuncAnimation

plt.rcParams['animation.html'] = 'jshtml'

class EnthalpyNormalizer(colors.Normalize):
    # https://matplotlib.org/stable/tutorials/colors/colormapnorms.html#twoslopenorm-different-mapping-on-either-side-of-a-center
    def __init__(self, T_min=-8.0, W_max=1.0, vcenter=0.0, ncenter=2/3, clip=False):
        # value which is the midpoint
        self.vcenter = vcenter
        # normalized location of the midpoint, normally 0.5
        self.ncenter = ncenter
        super().__init__(T_min, W_max, clip)

    def __call__(self, value, clip=None):
        # Note also that we must extrapolate beyond vmin/vmax
        x, y = [self.vmin, self.vcenter, self.vmax], [0, self.ncenter, 1.]
        return np.ma.masked_array(np.interp(value, x, y,
                                            left=-np.inf, right=np.inf))
    def inverse(self, value):
        y, x = [self.vmin, self.vcenter, self.vmax], [0, self.ncenter, 1]
        return np.interp(value, x, y, left=-np.inf, right=np.inf)

def custom_RdBu(norm):
    # https://stackoverflow.com/a/45787677/10221482
    cold_colors = plt.cm.Reds_r( np.linspace(0.0, 1.0,128))
    warm_colors = plt.cm.Blues(np.linspace(0, 1.,128))
    poly_colors = list(zip(np.linspace(0.00, norm.ncenter, 128), cold_colors)) +\
                  list(zip(np.linspace(norm.ncenter, 1.00, 128), warm_colors))
    poly_cmap   = colors.LinearSegmentedColormap.from_list('mycmap', poly_colors)

    return poly_cmap


def enthalpy_pcolormesh(src, i, axes=None, T_min=-8.0, W_max=1.0, contour_CTS=False):
    """
    i is time index (ith timestep)
    """

    if axes is None:
        fig, ax = plt.subplots(1,1, figsize=(6,3), constrained_layout=True)
    else:
        ax = axes
    # create custom norm and colormap for plotting the enthalpy field
    norm = EnthalpyNormalizer(T_min=T_min, W_max=W_max)
    cmap = custom_RdBu(norm)

    # based on phase change enthaply index appropriate variable
    T_and_w = xr.where(src.enthalpy_h >= src['phase change enthalpy'],
                       src['water content']*100,
                       src['temperature'] )

    # plot the temperature/water content field
    im = ax.pcolormesh(src.X[:,::-1], src.Z.isel(t=i), T_and_w.isel(t=i),
                       shading='gouraud', norm=norm, cmap=cmap)

    # contour the CTS if the flag is passed
    if contour_CTS:
        ax.contour(src.X[:,::-1],
                   src.Z.isel(t=i),
                   (src['water content']*100).isel(t=i),
                   levels=0.0, colors='k', linewidths=0.75)


    # ax.set_xlim(None, 5.5e3)
    # ax.set_ylim(1900, None)
    #
    # # set aspect ratio of plot, i.e. vertical exheration
    # ratio = 1/2
    # xleft, xright = ax.get_xlim()
    # ybottom, ytop = ax.get_ylim()
    # ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)

    if axes is None:
        # FIXME: These need to be set programtically..........
        # cb = plt.colorbar(im,shrink=ratio*1.2)
        cb = plt.colorbar(im)
        cb.set_ticks(np.concatenate((np.linspace(-8, 0, 9),
                                     np.linspace(0.2, 1.0, 5))))

    # # FIXME: These need to be set programtically..........
    # ax.set_xlim(0, 500)
    # ax.set_ylim(2600, 3000)
    # ax.set_xlim(None, 5.5e3)
    # ax.set_ylim(1900, None)
    if axes is None:
        return fig, ax, cb
    else:
        return im


def animate_enthalpy(src, i_0=0, i_f=None, i_s=10, interval=100):

    # create the initial frame
    fig, ax, cb = enthalpy_pcolormesh(src, 0)

    # Add time annotation
    label = ax.text(0.9, 0.9, r"$t=0$",  transform=ax.transAxes,
                       ha='center', va='center')

    def update(i):
        # wipe the axis
        ax.clear()
        # update the axis with the current timestep
        enthalpy_pcolormesh(src, i, axes=ax)
        # update time annotation
        label = ax.text(0.9, 0.9, "t={}".format(src.t[i].values), transform=ax.transAxes,
                           ha='center', va='center')

    # if final index is passed, plot the entire time acis
    if i_f is None:
        i_f = src.t.size

    ani = FuncAnimation(fig, update, frames=np.arange(i_0, i_f, i_s), interval=interval)

    return ani
