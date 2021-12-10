
#!/usr/bin/env python3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting script
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import re
import os
import sys
import glob
import argparse
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def make_colorbar(mf_dataset):
    #---------------------------------------------------------------------------
    # For Seting up the colorbar:
    #    - http://csc.ucdavis.edu/~jmahoney/matplotlib-tips.html
    #---------------------------------------------------------------------------
    cmap = cm.plasma
    norm = mcolors.Normalize(vmin=np.min(mf_dataset.Delta_MB.values),
                             vmax=np.max(mf_dataset.Delta_MB.values))

    s_map = cm.ScalarMappable(norm=norm, cmap=cmap)
    s_map.set_array(np.linspace(mf_dataset.Delta_MB.values.min(),
                                mf_dataset.Delta_MB.values.max(),
                                mf_dataset.Delta_MB.size))

    # If color parameters is a linspace, we can set boundaries in this way
    halfdist = (mf_dataset.Delta_MB[1] - mf_dataset.Delta_MB[0]) / 2.0
    bounds   = np.linspace(mf_dataset.Delta_MB.values.min()  - halfdist,
                           mf_dataset.Delta_MB.values.max()  + halfdist,
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
        color = cmap(1-norm(delta_mb))
        ax.plot(Vol.t[1:], Vol.sel(Delta_MB=delta_mb)[1:], color=color)

    ax.axhline(1.0,c='k',ls=':',lw=1, alpha=0.5)

    cbar = fig.colorbar(s_map,
                    spacing='proportional',
                    ticks=np.linspace(mf_dataset.Delta_MB.values.min(),
                                      mf_dataset.Delta_MB.values.max(),
                                      mf_dataset.Delta_MB.size),
                    ax=ax,
                    boundaries=bounds,
                    drawedges=True,
                    format='%2.{}f'.format(precision))


    ax.set_title(title)

    # annotate the figures axes
    ax.set_ylabel('Relative Volume per Unit Width')
    ax.set_xlabel('Time (yr)')
    # annotate the colorbar axes
    cbar.set_label('$\Delta \dot b$ (m i.e.q. yr$^{-1}$)', rotation=270, labelpad=20)
    #cbar.ax.tick_params(labelsize=7)

    fig.tight_layout()

    return fig, ax

def plot_final_z_s(mf_dataset, precision=3, title=''):
    # Make a colormap and all the associated var names
    cmap, norm, s_map, bounds = make_colorbar(mf_dataset)

    fig, ax = plt.subplots(figsize=(10, 3))

    for delta_mb in mf_dataset.Delta_MB:
        color = cmap(1-norm(delta_mb))
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
                    y2=np.minimum(mf_dataset.isel(t=0,Delta_MB=0,coord_2=-1).zbed.min(), 0.0),
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
    ax.set_ylabel('Elevation  [m a.s.l.]')
    ax.set_xlim(0,np.max(mf_dataset.coord_1)/1000.)
    ax.set_ylim( mf_dataset.z_s.min() - \
                 (mf_dataset.Z.max() - mf_dataset.Z.min()) / 20, None)

    cbar.set_label('$\Delta \dot b$ (m i.e.q. yr$^{-1}$)', rotation=270, labelpad=20)
    #cbar.ax.tick_params(labelsize=6)

    fig.tight_layout()

    ratio = 1/3
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
    return fig, ax

def main(argv):

    #---------------------------------------------------------------------------
    # Specify command line arguments
    #---------------------------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("src_path", metavar="path", type=str, nargs='+',
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
    parser.add_argument('-TeX','--use_LaTeX', action='store_true',
                        help = "use LaTeX for text rendeing ")

    args, _ = parser.parse_known_args(argv)

    volume_plot = args.plot_volume
    z_s_plot    = args.plot_Z_s
    out_fn      = args.output_filename
    title       = args.title
    #---------------------------------------------------------------------------

    if args.use_LaTeX:
        plt.rcParams['text.usetex'] = True

    #---------------------------------------------------------------------------
    # Load and concatenate the .nc files
    #---------------------------------------------------------------------------

    if type(args.src_path) == list:
        if "*" in args.src_path[0]:
            # Glob the file paths and return list of files
            files = sorted(glob.glob(args.src_path[0]))
    else:
        files = args.src_path

    # Raise error if glob didn't work
    if not files:
        raise OSError('value passed for "src_path" is invalid')

    # We need to do some more filtering cause there will be lots of misc files
    # in th nc directories which might be matched by a overly general glob
    stride_length = len(args.mb_range[1].split(".")[1])

    # loop over the globed files
    for file in files[:]:
        match = []
        # loop over mb offset values
        for MB in np.arange(float(args.mb_range[0]),
                            float(args.mb_range[2]) + float(args.mb_range[1]),
                            float(args.mb_range[1])):

            # if MB offset is not in any of the globed filenames will return flase
            match.append(str(np.round(MB, stride_length)) in file)

            # Deal wih -0.0 error
            if np.abs(np.round(MB, stride_length)) == 0.0:
                match.append(str(np.round(np.abs(MB), stride_length)) in file)
            else:
                match.append(str(np.round(MB, stride_length)) in file)

        # if none of the MB offsets are contained in the filename drop it
        if not any(match):
            files.remove(file)

    # Sorting isn't guaranted to work corretly, so if pattern is matched
    # do special sorting ensure it's done correctly
    regex = re.search('MB_-*\d*\d\.\d*\d_OFF', files[0])
    if regex:
        files.sort(key = lambda x: float(x.split('MB_')[-1].split('_OFF')[0]),
                   reverse = True)


    # Create array of mass balance values used in spin-up
    MB, dx = np.linspace(float(args.mb_range[0]),
                         float(args.mb_range[2]),
                         len(files),
                         retstep=True)


    # Check the the MB stride is the same as the stride that was specified
    if not np.isclose(dx, float(args.mb_range[1])):

        print(dx)
        raise OSError('MB stride passed does not match that of files found')

    # Make an empty list to store the read in .nc files
    xarrays = []

    # Iterate over each .nc file and read in with xarray
    for file in files:
        with xr.open_dataset(file) as src:
                # correct for minimum ice thickness
                src["height"] = xr.where(src.height <= 10, 0, src.height)
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
