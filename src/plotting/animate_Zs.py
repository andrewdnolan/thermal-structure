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
import matplotlib.pyplot as plt
from matplotlib import animation, rc

# Set some matplotlib parameters
plt.rcParams['text.usetex']    = False
plt.rcParams['animation.html'] = 'jshtml'

def animate_Zs(result, title=''):
    fig, ax = plt.subplots(figsize=(9,3), constrained_layout=True)

    # Set the x and y limits, which do not change throughout animation
    ax.set_xlim(0,np.max(result.coord_1)/1000.)
    ax.set_ylim(result.z_s.min(),
                result.z_s.isel(t=0).max() + \
                (result.z_s.isel(t=0).max() - result.z_s.isel(t=0).min())/ 10 )

    # Set axes labels, which do not change throughout animation
    ax.set_ylabel('Elevation [m a.s.l.]')
    ax.set_xlabel('Distance along flowline [km]')

    # Plot the bed which is constant throughout animation
    ax.plot(result.coord_1/1000.,
            result.isel(t=0,coord_2=-1).zbed,
            color='k', label=r'$z_{\rm b}$')
    # Fill between bed and bottom of plot
    ax.fill_between(result.coord_1/1000.,
                    result.isel(t=0,coord_2=-1).zbed,
                    y2=np.minimum(result.isel(t=0,coord_2=-1).zbed.min(), 0.0),
                    color='gray', alpha=0.5)

    ax.legend(loc=2)
    ax.set_title(title)

    line1, = ax.plot([], [], lw=2, color='tab:blue', label='$z_s(t=0.0)$',)
    line   = [line1]

    # Function to be animated
    def animate(i):
        # plot the free surface
        line[0].set_data(result.coord_1/1000.,
                         result.isel(t=i,coord_2=-1).z_s)
        line[0].set_label('$z_s(t={{{:.1f}}})$'.format(result.t.isel(t=i).values))

        ax.fill_between(result.coord_1/1000.,
                        result.isel(t=0,coord_2=-1).zbed,
                        y2=np.minimum(result.isel(t=0,coord_2=-1).zbed.min(), 0.0),
                        color='gray', alpha=0.5)
        ax.legend(loc=2)
        return line

    NT = result.t.size
    anim = animation.FuncAnimation(fig, animate,
                                   frames=np.arange(0, NT, 10),
                                   interval=50,
                                   blit=True)

    return anim
def main(argv):
    #---------------------------------------------------------------------------
    # Specify command line arguments
    #---------------------------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("src_path", metavar="path", type=str,
                        help = "Path to .nc file to be plotted"\
                               "enclose in quotes")
    parser.add_argument('-T','--title', type=str,
                        help = "string for the title of the plot")
    parser.add_argument('-out_fn','--output_filename', type=str,
                        help = "full path to the output figure")

    args, _ = parser.parse_known_args(argv)

    out_fn      = args.output_filename
    title       = args.title
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    # Load the .nc files
    #---------------------------------------------------------------------------


    # Raise error if glob didn't work
    if not os.path.exists(args.src_path):
        raise OSError('\n value passed for "src_path" is invalid \n')

    with xr.open_dataset(args.src_path) as src:
        # correct for minimum ice thickness
        src["height"] = xr.where(src.height <= 10, 0, src.height)
        # apply sigma coordinate transform for vertical coordinate
        src["z_s"]     = src.zbed + src.Z * src.height
        # Calculate the magnitude of the velocity vectors
        src['vel_m'] = np.sqrt(src['velocity 1']**2 + src['velocity 2']**2)

    fig = animate_Zs(src, title=title)

    # Write the plot to a file
    fig.save(out_fn, dpi=400, fps=10)

    # Don't know if I actually need this
    plt.close()

if __name__ == '__main__':
    main(sys.argv[1:])
