#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import animation
from .plotting import enthalpy_pcolormesh, get_axis_limits

class AnimateEnthalpy: 

    def __init__(self, src, interval=150, frames=np.arange(0,2000,10)):

        # Initialize the figure and axes objects
        self.create_figure()
        # Make the input dataset available to the underlying functions
        self.src = src
        # # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update,
                                           frames=frames,
                                           interval=interval,
                                           init_func=self.setup_plot,
                                           blit=False, repeat=False,
                                           cache_frame_data=False)
    def create_figure(self):
        self.fig, self.ax = plt.subplots(figsize=(6,3), constrained_layout=True)

    def setup_plot(self):
        """Initial plot."""
        self.enthalpy_h =enthalpy_pcolormesh(self.src, 0, axes=self.ax)
        
        # create the colorbar
        self.H_cb = self.fig.colorbar(self.enthalpy_h, ax=self.ax, extend='both')

        # special colorbar formatting to deal with broken (enthalpy) axis
        self.H_cb.set_ticks(np.concatenate((np.linspace(-8, 0, 5),
                                            np.linspace(0.1, 0.5, 3))))
        self.H_cb.ax.tick_params(labelsize='small')
        self.H_cb.set_label('$\omega$ [\%]  /  $T\'$ [$^\circ$C]', rotation=270, labelpad=25)

        # Set the axis labels
        self.ax.set_ylabel('Elevation [m a.s.l.]')
        self.ax.set_xlabel('Distance [km]')

        # get the axis bounds for the pcolormesh
        xlim, ylim = get_axis_limits(self.src)
        # set the pcolormesh axis limits
        self.ax.set_xlim(*xlim)
        self.ax.set_ylim(*ylim)

        # Add time annotation
        self.time_annot = self.ax.text(0.85, 0.9, "$t={{{:6.0f}}}$".format(float(self.src.t[0])),
                                       transform=self.ax.transAxes, ha='center', va='center')

        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.enthalpy_h,

    def update(self, i):
        """Update the scatter plot."""

        # update the pcolormesh
        self.enthalpy_h.remove()
        self.enthalpy_h = enthalpy_pcolormesh(self.src, i, axes=self.ax)

        # update the time annotation
        self.time_annot.set_text("$t={{{:6.0f}}}$".format(float(self.src.t[i])-1.0))

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.enthalpy_h,