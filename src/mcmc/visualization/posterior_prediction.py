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


def plot_hierarchical(obs_df, pred, keys):

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

    fig, axes = plt.subplots(n_glac, 1, figsize=(4,8), sharex=True)

    if n_glac == 1:
        axes = [axes]
    print(axes)
    for i, key in enumerate(obs_df.key.unique()):

        idxs   = obs_df[(obs_df["key"] == key)].index
        z_obs  = obs_df[(obs_df["key"] == key)].z
        B_obs  = obs_df[(obs_df["key"] == key)].MB

        pp = pred[:,:,idxs]
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

        axes[i].text(0.15, 0.9, key, ha="center", va="center", transform=axes[i].transAxes)

    plt.tight_layout()
