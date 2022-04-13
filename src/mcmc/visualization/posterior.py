#!/usr/bin/env python3

"""
collection of functions for visulaizing the posterior distibutions
"""

import numpy as np
import arviz as az
import xarray as xr
from scipy import stats
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from matplotlib.ticker  import NullFormatter


plt.rcParams.update({'text.usetex': True,
                     'figure.facecolor': 'w',
                     'savefig.bbox': 'tight'})

def dist_bounds(samples, n=30):
    """Find mix, max and (n) bins for an array of samples"""
    min_val, max_val = samples.min(), samples.max()
    bins = np.linspace(min_val, max_val, n)

    return max_val, max_val, bins

def chain_plot():
    # fig 2 from Brinkerhoff et al. 2016
    # https://dbrinkerhoff.files.wordpress.com/2016/10/subglacial_hydro_bayes.pdf
    pass

def pair_plot(posterior, var_keys, out_fp=None):
    # Most of this function is an exact replica of :
    #https://github.com/pism/pism-emulator/blob/main/speedemulator/plot_posterior.py

    # if arviz InferenceData object was passed, get the posterior attribute
    if isinstance(posterior, az.InferenceData):
        posterior = posterior.posterior

    # make sure the posterior passed is a xarray object
    if not isinstance(posterior, xr.Dataset):
        raise TypeError("posterior must be xarray Dataset")

    # collapse all chains and draws along one dimenesion for plotting
    if "draw" and "chain" in posterior.coords:
        posterior = posterior.stack(samples=("chain", "draw"))

    keys = list(var_keys)
    n_params = len(var_keys)

    #-----------------------------------------------------
    # compute correlation coeff for upper diagonal pathces
    #-----------------------------------------------------
    # convert to pandas df, much easier to work with np.corrcoeff
    df = posterior[keys].to_dataframe().reset_index(drop=True)
    # compute unbiased correlation coeff
    C_0  = np.corrcoef((df - df.mean(axis=0)).T)
    # rescale corrcoeff from [-1, 1] --> [0, 1] for cmap plotting
    Cn_0 = (np.sign(C_0) * C_0 ** 2 + 1) / 2.0
    #-----------------------------------------------------

    fig, axs = plt.subplots(nrows=n_params, ncols=n_params, figsize=(6,6))
    fig.subplots_adjust(hspace=0.0, wspace=0.0)

    # plotting loop
    for i in range(n_params):
        for j in range(n_params):

            # below the digonal (scatter plots)
            if i > j:
                axs[i,j].scatter(posterior[keys[j]],
                                 posterior[keys[i]],
                                 s=0.1,
                                 alpha=0.05,
                                 rasterized=True)
                # gaussian kde contours
                min_val, max_val, bins_j = dist_bounds(posterior[keys[j]])
                min_val, max_val, bins_i = dist_bounds(posterior[keys[i]])

                kernel = stats.gaussian_kde(posterior[[keys[j], keys[i]]].to_array())

                # rescale from bin edges to bin centers
                bj = 0.5 * (bins_j[1:] + bins_j[:-1])
                bi = 0.5 * (bins_i[1:] + bins_i[:-1])
                Bj, Bi = np.meshgrid(bj, bi)

                axs[i, j].contour(Bj,
                                  Bi,
                                  kernel(np.vstack((Bj.ravel(),
                                                    Bi.ravel()))).reshape(Bj.shape),
                                  5,
                                  linewidths=0.5,
                                  colors="black",
                )

                # set axis limits
                axs[i, j].set_xlim(posterior[keys[j]].min(),
                                   posterior[keys[j]].max())
                axs[i, j].set_ylim(posterior[keys[i]].min(),
                                   posterior[keys[i]].max())

            # digonal (dist plot)
            elif i == j:
                min_val, max_val, bins = dist_bounds(posterior[keys[i]])
                ith_post_hist, b = np.histogram(posterior[keys[i]], bins, density=True)

                # rescale from bin edges to bin centers
                b = 0.5 * (b[1:] + b[:-1])

                axs[i,j].plot(b, ith_post_hist*0.5, c='k')

            # above the diagonal (correlation patches)
            elif i < j:
                patch_upper = Polygon(
                    np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]]),
                    facecolor=plt.cm.seismic_r(Cn_0[i, j]),
                )

                axs[i, j].add_patch(patch_upper)

                axs[i, j].text(
                    0.5,
                    0.5,
                    "{0:+.2f}".format(C_0[i, j]),
                    horizontalalignment="center",
                    verticalalignment="center",
                    transform=axs[i, j].transAxes,
                    color='k',
                )
    #------------------
    # formatting loops:
    #------------------
    # loop over first column and set y-label
    for i, ax in enumerate(axs[:, 0]):
        ax.set_ylabel(var_keys[keys[i]])

    # loop over last row and set x-label
    for j, ax in enumerate(axs[-1, :]):
        ax.set_xlabel(var_keys[keys[j]])
        ax.tick_params("x", rotation=45)

        # turn off y-axis for all but first column
        if j > 0:
            ax.yaxis.set_minor_formatter(NullFormatter())
            ax.yaxis.set_major_formatter(NullFormatter())
            ax.tick_params(axis="y", which="both", length=0)

    # Loop over all rows (except last) in first column and turn off x-axis
    for ax in axs[:-1, 0].ravel():
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis="x", which="both", length=0)

    # Loop over interior rows and cols, turning off both x and y axis
    for ax in axs[:-1, 1:].ravel():
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis="both", which="both", length=0)

    #Save the figure
    if out_fp:
        print(f"Saving figure to {out_fp}")
        fig.savefig(out_fp, dpi=300)
    else:
        plt.show()

def compare_prior_and_posterior():
    pass
