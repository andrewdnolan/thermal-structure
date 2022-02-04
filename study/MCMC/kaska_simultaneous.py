#!/usr/bin/env python3

import PDD
import theano
import warnings
import numpy as np
import pymc3 as pm
import arviz as az
import xarray as xr
import matplotlib.pyplot as plt

warnings.filterwarnings(action='ignore')
theano.config.compute_test_value = 'warn'
plt.rcParams.update({'figure.facecolor': 'w', 'savefig.bbox':'tight'})

def plot_BAM(z: np.ndarray, ppc: tuple, obs: tuple, out_fp: str, axis_labels=None):
    """Plot net (B)alance, (A)ccumulation, and (M)elt (BAM).
    """

    if axis_labels == None:
        axis_labels = ('Mass Balance (m/yr)','Accumulation (m/yr)','Melt (m/yr)')

    fig, ax = plt.subplots(3, 1, figsize=(5,6), sharex=True)

    for i, (pp, obs) in enumerate(zip(ppc, obs)):

        ax[i].scatter(z, obs, s=1.0, c='#7570b3')

        ax[i].plot(z, np.mean(pp, axis=(0, 1)), lw=1.0, color='#1b9e77')

        az.plot_hdi(z, pp, hdi_prob=0.6, ax=ax[i], color='#1b9e77',
                        fill_kwargs={"alpha": 0.5, "linewidth": 0})

        az.plot_hdi(z, pp, hdi_prob=0.95, ax=ax[i], color='#1b9e77',
                        fill_kwargs={"alpha": 0.33, "linewidth": 0})

        ax[i].set_ylabel(axis_labels[i])

    ax[2].set_xlabel('Elevation [m a.s.l.]')
    plt.tight_layout()
    fig.savefig(out_fp, dpi=300)

def plot_posterior(trace, vars: list, out_fp: str):
    """Wrapper for arviz trace plot
    """
    _ = az.plot_trace(trace, var_names=vars)
    plt.tight_layout()
    plt.savefig(out_fp, dpi=300)

def read_observation():
    """Read and return EMY/KMR results
    """

    # file paths to the NetCDF files
    nc_fp = "../MB_tune/Young_etal_2020_ref_MB.nc"

    with xr.open_dataset(nc_fp) as MB_ref:
        # average over the 25 reference model runs and stack along elevation
        stacked = MB_ref.mean('ref_run').stack(z=('x','y'))

        # Sort the indexes by elevation
        idxs  = stacked.dropna('z').z.values
        elev  = stacked.Elevation.dropna('z').values
        idxs  = idxs[np.argsort(elev)]

        # Calculate the Kaskawulsh ELA
        ELA_idxs = np.argpartition(np.abs(stacked.MB.values), 5)
        z_ELA    = float(stacked.isel(z=ELA_idxs[:5]).Elevation.mean())

        # surface elevation of observed values
        z_obs = stacked.Elevation.sel(z=idxs).values
        # Observed Accumulation [ m i.e. / yr ]
        A_obs = stacked.Accumulation.sel(z=idxs).values
        # Observed melt         [ m i.e. / yr ]
        M_obs = stacked.Melt.sel(z=idxs).values
        # Observed mass balance [m i.e. / yr ]
        B_obs = stacked.MB.sel(z=idxs).values

    return z_obs, A_obs, M_obs, B_obs, z_ELA

def fit_LA(z_obs, A_obs, M_obs, B_obs, z_ELA, draws=4000, tune=2000, cores=1):
    """Simultaneous fitting the linear accumulation model
    """
def fit_PWA(z_obs, A_obs, M_obs, B_obs, z_ELA, draws=4000, tune=2000, cores=1):
    """Simultaneous fitting the Piecewise accumulation model
    """

    # dictionary of PDD model parameters
    const =  dict(T_m = 0.0, T_rs = 1.0, α = 10.5, T_ma = -7.77, ΔTΔz  = 6.5E-3,
                  T_p = 196, ref_z = 2193, T_σ = 8.6, z_ELA = z_ELA)

    # initalize the PDD melt model class
    PDD_PWA = PDD.PDD_PWA(**const)

    # Define Priors
    with pm.Model() as model:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ----> Mass balance Model (physical priors)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Degree Day factor: Braithwaite (2008) / Rounce et al. 2020
        f_s_prior = pm.TruncatedNormal("f_s", mu=4.1, sigma=1.5, lower=0.0)
        # Somewhat base of Aschwanden et al. 2019
        C_prior = pm.TruncatedNormal("C", mu=2.0, sigma=0.75, lower=1, upper=4.0)
        # Somewhat base of Aschwanden et al. 2019
        f_r_prior = pm.TruncatedNormal('f_r', mu=0.5, sigma=0.2, lower=0.0, upper=1)
        #Accumulation function priors
        z_max_prior  = pm.Normal("z_max",  mu = 3000,   sigma = 500)
        P0_prior     = pm.Normal("P0",     mu = 0.5,    sigma = 0.5)
        ΔPΔz_1_prior = pm.Normal("ΔPΔz_1", mu = 1.5e-4, sigma = 1e-4)
        ΔPΔz_2_prior = pm.Normal("ΔPΔz_2", mu = 2.5e-3, sigma = 1e-3)
        P_max_prior  = pm.Normal("P_max",  mu = 4.0,    sigma = 1.0)

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ----> Hyperparameters (likelihood related priors)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        A_sigma = pm.HalfCauchy("A_sigma", 1)
        M_sigma = pm.HalfCauchy("M_sigma", 5)
        B_sigma = pm.HalfCauchy("B_sigma", 2)

    # Define Forward model (wrapped through theano)
    with model:
        # components of mass balance model
        A, R, M = PDD_PWA.tt_forward(z_obs, f_s_prior, C_prior, f_r_prior,
                  z_max_prior, P0_prior, ΔPΔz_1_prior, ΔPΔz_2_prior, P_max_prior)
        # net balance [m i.e. / yr ]
        B = A + R - M

    # Define likelihood (function?)
    with model:
        # Individual likelihood functions for each component
        B_est = pm.Normal("B_est", mu=B, sigma=B_sigma, observed=B_obs)
        A_est = pm.Normal("A_est", mu=A, sigma=A_sigma, observed=A_obs)
        M_est = pm.Normal("M_est", mu=M, sigma=M_sigma, observed=M_obs)

        potential = pm.Potential("obs", B_est.sum()*A_est.sum()*M_est.sum())

    # run inference: Sample
    with model:
        trace = pm.sample(init='adapt_diag',
                          draws=draws,
                          tune=tune,
                          cores=cores,
                          target_accept=0.9,
                          return_inferencedata=True);

    # do posterior predictive inference
    with model:
        ppc = pm.sample_posterior_predictive(trace,
                                             var_names=["B_est", "A_est", "M_est"],
                                             keep_size=True)
        pp_B = ppc["B_est"]
        pp_A = ppc["A_est"]
        pp_M = ppc["M_est"]

    return model, trace, pp_B, pp_A, pp_M

def run_PWA(draws=4000, tune=2000, cores=1):

    # load observations
    z_obs, A_obs, M_obs, B_obs, z_ELA = read_observation()

    # fit the PWA model with MCMC
    model, trace, pp_B, pp_A, pp_M = fit_PWA(z_obs, A_obs, M_obs, B_obs,
                                             z_ELA, draws, tune, cores)

    # plot the traces
    vars   = ['z_max', 'P0', 'ΔPΔz_1', 'ΔPΔz_2', 'P_max', 'f_s', 'C', 'f_r']
    out_fp = './simultaneous/PWA/trace.png'
    plot_posterior(trace, vars=vars, out_fp=out_fp)

    # plot the mass balance results and components
    out_fp = './simultaneous/PWA/BAM_PWA.png'
    plot_BAM(z_obs, (pp_B, pp_A, pp_M), (B_obs, A_obs, M_obs), out_fp=out_fp)

    # save the trace
    az.to_netcdf(trace,'./simultaneous/PWA/trave.nc')


if __name__ == '__main__':
    run_PWA(2000, 2000, 2)
