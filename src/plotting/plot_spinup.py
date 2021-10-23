
#!/usr/bin/env python3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting script
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import sys
import glob
import argparse
import numpy as np
import xarray as xr
import pandas as pd
import scipy.linalg as LA
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib import animation, rc
import matplotlib.colors as mcolors


# Set some matplotlib parameters
plt.rcParams['text.usetex']    = False
plt.rcParams['animation.html'] = 'jshtml'

def make_colorbar(mf_dataset):
    #---------------------------------------------------------------------------
    # For Seting up the colorbar:
    #    - http://csc.ucdavis.edu/~jmahoney/matplotlib-tips.html
    #---------------------------------------------------------------------------
    cmap = cm.plasma
    norm = mcolors.Normalize(vmin=np.min(mf_dataset.Delta_MB),
                             vmax=np.max(mf_dataset.Delta_MB))

    s_map = cm.ScalarMappable(norm=norm, cmap=cmap)
    s_map.set_array(mf_dataset.Delta_MB)

    # If color parameters is a linspace, we can set boundaries in this way
    halfdist = (mf_dataset.Delta_MB[1] - mf_dataset.Delta_MB[0]) / 2.0
    bounds   = np.linspace(mf_dataset.Delta_MB[0]   - halfdist,
                           mf_dataset.Delta_MB[-1]  + halfdist,
                           len(mf_dataset.Delta_MB) + 1)

    return cmap, norm, s_map, bounds

def plot_volume(mf_dataset, precision=3, title=''):

    # Make a volume xarray
    Vol = mf_dataset.height.isel(coord_2=-1).integrate("coord_1") /\
          mf_dataset.height.isel(coord_2=-1).isel(t=0).integrate("coord_1")

    # Make a colormap and all the associated var names
    cmap, norm, s_map, bounds = make_colorbar(mf_dataset)


    fig, ax = plt.subplots(figsize=(7, 5))

    for delta_mb in Vol.Delta_MB:
        color = cmap(norm(delta_mb))
        ax.plot(Vol.t[1:], Vol.sel(Delta_MB=delta_mb)[1:], color=color)

    ax.axhline(1.0,c='k',ls=':',lw=1, alpha=0.5)

    cbar = fig.colorbar(s_map,
                    spacing='proportional',
                    ticks=mf_dataset.Delta_MB,
                    ax=ax,
                    boundaries=bounds,
                    drawedges=True,
                    format='%2.{}f'.format(precision))


    ax.set_title(title)

    # annotate the figures axes
    ax.set_ylabel('Relative Volume per Unit Width (km$^2$)')
    ax.set_xlabel('Time (a)')
    # annotate the colorbar axes
    cbar.set_label('$\Delta \dot b$ (m a$^{-1}$)', rotation=270, labelpad=20)
    cbar.ax.tick_params(labelsize=7)

    fig.tight_layout()

    return fig, ax

def plot_final_z_s(mf_dataset, precision=3, title=''):
    # Make a colormap and all the associated var names
    cmap, norm, s_map, bounds = make_colorbar(mf_dataset)

    fig, ax = plt.subplots(figsize=(9, 5))

    for delta_mb in mf_dataset.Delta_MB:
        color = cmap(norm(delta_mb))
        ax.plot(mf_dataset.coord_1/1000.,
                mf_dataset.isel(t=-1, coord_2=-1).z_s.sel(Delta_MB=delta_mb),
                color=color)

    ax.plot(mf_dataset.coord_1/1000.,
            mf_dataset.isel(t=0,Delta_MB=0,coord_2=-1).zbed,
            color='k', label=r'$z_{\rm b}$')

    ax.plot(mf_dataset.coord_1/1000.,
            mf_dataset.isel(t=1,Delta_MB=0,coord_2=-1).z_s,
            color='k', ls=':', lw=0.5, alpha = 0.5, label=r'$z_{\rm s}(t=0)$')

    ax.fill_between(mf_dataset.coord_1/1000.,
                    mf_dataset.isel(t=0,Delta_MB=0,coord_2=-1).zbed,
                    color='gray', alpha=0.5)

    cbar = fig.colorbar(s_map,
                        spacing='proportional',
                        ticks=mf_dataset.Delta_MB,
                        ax=ax,
                        boundaries=bounds,
                        drawedges=True,
                        format='%2.2f')

    ax.legend(loc=2)

    ax.set_title(title)

    ax.set_xlabel('Length (km)')
    ax.set_ylabel('m a.s.l.')
    ax.set_xlim(0,np.max(mf_dataset.coord_1)/1000.)
    ax.set_ylim( mf_dataset.z_s.min() - \
                 (mf_dataset.Z.max() - mf_dataset.Z.min()) / 20, None)

    cbar.set_label('$\Delta \dot b$ (m a$^{-1}$)', rotation=270, labelpad=20)
    cbar.ax.tick_params(labelsize=6)

    fig.tight_layout()

    return fig, ax

def main(argv):

    #---------------------------------------------------------------------------
    # Specify command line arguments
    #---------------------------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("src_path", metavar="path", type=str,
                        help = "Path to .nc files to be plotted"\
                               "enclose in quotes, accepts * as wildcard for directories or filenames")
    parser.add_argument('-mb','--mb_range', nargs='+',
                        help = "mimics 'seq' unix commands where:"\
                               " first value is start"\
                               " middle values is stride"\
                               " last value is stop")
    parser.add_argument('-T','--title', type=str,
                        help = "string for the title of the plot")
    parser.add_argument('-V','--plot_volume', action='store_true',
                        help = "volume convergence plots after mass balance grid search")
    parser.add_argument('-Z_s','--plot_Z_s',  action='store_true',
                        help = "final z_s after mass balance grid search")
    parser.add_argument('-out_fn','--output_filename', type=str,
                        help = "full path to the output figure")

    args, _ = parser.parse_known_args(argv)

    volume_plot = args.plot_volume
    z_s_plot    = args.plot_Z_s
    out_fn      = args.output_filename
    title       = args.title
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    # Load and concatenate the .nc files
    #---------------------------------------------------------------------------

    # Glob the file paths and return list of files
    files = sorted(glob.glob(args.src_path))

    # Raise error if glob didn't work
    if not files:
        raise OSError('value passed for "src_path" is invalid')

    # Create array of mass balance values used in spin-up
    MB, dx = np.linspace(float(args.mb_range[0]),
                         float(args.mb_range[2]),
                         len(files),
                         retstep=True)

    # Check the the MB stride is the same as the stride that was specified
    if not np.isclose(dx, float(args.mb_range[1])):
        raise OSError('MB stride passed does not match that of files found')

    # Make an empty list to store the read in .nc files
    xarrays = []

    # Iterate over each .nc file and read in with xarray
    for file in files:
        with xr.open_dataset(file) as src:
                # correct for minimum ice thickness
                src["depth"] = xr.where(src.depth <= 10, 0, src.depth)
                # apply sigma coordinate transform for vertical coordinate
                src["z_s"]     = src.zbed + src.Z * src.height
                # Calculate the magnitude of the velocity vectors
                src['vel_m'] = np.sqrt(src['velocity 1']**2 + src['velocity 2']**2)

        xarrays.append(src)

    # Concatenate the .nc files via their mass balance offset
    mf_dataset = xr.concat(xarrays,
                           pd.Index(data = MB, name='Delta_MB'))
    #---------------------------------------------------------------------------

    out_fn = args.output_filename

    if volume_plot:
        fig, _  = plot_volume( mf_dataset,
                           precision=len(args.mb_range[1])-2,
                           title=title)
    if z_s_plot:
        fig, _  = plot_final_z_s( mf_dataset,
                              precision=len(args.mb_range[1])-2,
                              title=title)

    # Write the plot to a file
    fig.savefig(out_fn, dpi=400, bbox_inches='tight', facecolor='w')

    # Don't know if I actually need this
    plt.close()

if __name__ == '__main__':
    main(sys.argv[1:])
