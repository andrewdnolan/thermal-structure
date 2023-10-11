#!/usr/bin/env python3

"""
collection of functions for visulaizing the posterior distibutions
"""

import numpy as np
import arviz as az
import xarray as xr
from scipy import stats
import matplotlib.pyplot as plt


def plot_components(z, RAM):
    pass


def plot_BAMR(z: np.ndarray, ppc: tuple, obs: tuple, out_fp: str, axis_labels=None):
    """Plot net (B)alance, (A)ccumulation, and (M)elt (BAM).
    """

    if axis_labels == None:
        axis_labels = ('Refreezing (m a$^{-1}$)','Accumulation (m a$^{-1}$)','Melt (m a$^{-1}$)', 'Mass Balance (m a$^{-1}$)')

    fig, ax = plt.subplots(4, 1, figsize=(5,7), sharex=True)

    for i, (pp, obs) in enumerate(zip(ppc, obs)):

        ax[i].scatter(z, obs, s=1.0, c='#7570b3')

        ax[i].plot(z, np.mean(pp, axis=(0, 1)), lw=1.0, color='#1b9e77')#,
                   #label=labels[i])

        az.plot_hdi(z, pp, hdi_prob=0.6, ax=ax[i], color='#1b9e77',
                        fill_kwargs={"alpha": 0.5, "linewidth": 0})

        az.plot_hdi(z, pp, hdi_prob=0.95, ax=ax[i], color='#1b9e77',
                        fill_kwargs={"alpha": 0.33, "linewidth": 0})

        ax[i].set_ylabel(axis_labels[i])

        ax[i].legend()

    ax[0].set_ylim(None, 0.5)
    ax[3].set_xlabel('Elevation [m a.s.l.]')
    plt.tight_layout()
    fig.savefig(out_fp, dpi=300)



def plot_hierarchical(obs_df, pred, keys, labels=None):

    if len(pred.shape) == 3:
        axis = (0,1)
    elif len(pred.shape) == 2:
        axis = 0

    # get one and two sigma from posterior prediction
    one_σ  = np.quantile(pred, 0.68, axis=axis)
    two_σ  = np.quantile(pred, 0.95, axis=axis)
    # find the mean of the posterior prediction
    mean   = np.mean(pred, axis=axis)


    # number of glacier "keys"
    n_glac = len(obs_df.key.unique())


    fig, axes = plt.subplots(n_glac, 1, figsize=(5,8), sharex=True, sharey=True, 
                             constrained_layout=True)

    bbox_props  = dict(boxstyle="round4, pad=0.3", fc="white", ec="black", lw=0.5)

    if n_glac == 1:
        axes = [axes]

    for i, key in enumerate(obs_df.key.unique()):

        idxs   = obs_df[(obs_df["key"] == key)].index
        z_obs  = obs_df[(obs_df["key"] == key)].z
        B_obs  = obs_df[(obs_df["key"] == key)].MB

        # get the predicitions, depending on the shape of posterior predicted array
        if len(pred.shape) == 3:
            pp = pred[:,:,idxs]
        elif len(pred.shape) == 2:
            pp = pred[:,idxs]

        B_pred = mean[idxs]
        B_1σ   = one_σ[idxs]
        B_2σ   = two_σ[idxs]

        # scatter plots the observations
        axes[i].scatter(z_obs, B_obs, s=1.0, c='#7570b3')

        axes[i].plot(z_obs, B_pred, lw=1.0, color='#1b9e77')

        az.plot_hdi(z_obs, pp, hdi_prob=0.68, ax=axes[i], color='#1b9e77',
                        fill_kwargs={"alpha": 0.5, "linewidth": 0})

        az.plot_hdi(z_obs, pp, hdi_prob=0.95, ax=axes[i], color='#1b9e77',
                        fill_kwargs={"alpha": 0.33, "linewidth": 0})


        if labels: 
            label = labels[i]
        else: 
            label = key

        axes[i].grid()
        axes[i].text(0.05, 0.85, label, 
                     bbox=bbox_props,
                     ha="left", va="center", 
                     transform=axes[i].transAxes)

        axes[i].set_ylabel("$B$ (m a$^{-1}$)")

    axes[i].set_xlabel("Elevation (m a.s.l.)")

    return fig, axes
