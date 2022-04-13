#!/usr/bin/env python3

import json
import numpy as np
from data import utils
from model import MCMC, PDD
from visualization import posterior

if __name__ == '__main__':
    # Load the KMR "observational" data
    z_obs, R_obs, A_obs, M_obs, B_obs, z_ELA = utils.load_KMR_DETIM()
    # Downsample the observational data, to efficency
    z_obs, R_obs, A_obs, M_obs, B_obs = z_obs[::10], R_obs[::10], A_obs[::10], M_obs[::10], B_obs[::10]
    # stack the observations into a matrix need for multivariate model
    RAM_obs = np.array([R_obs, A_obs, M_obs]).T


    # dictionary of PDD model parameters
    p = [ 8.29376332e-05, -3.45256005e-02,  6.31076200e+00]
    doy = np.arange(1,366)
    fancy_std = np.polyval(p, doy)[:, np.newaxis]

    const =  dict(T_m   = 0.0,    T_rs  = 1.0,
                  α     = 10.8,   T_ma  = -9.02,
                  ΔTΔz  = 6.5E-3, T_p   = 196,
                  ref_z = 2193,   T_σ   = fancy_std)

    # initialize the PDD melt model class
    PDD_forward = PDD.PositiveDegreeDayModel(**const)

    # Load the dictionary of physical priors
    with open('priors.json') as f:
        priors = json.load(f)

    MCMC_dict = dict(forward=PDD_forward)

    # initialize the MCMC inversion class
    model = MCMC.MultivariateInversion(**MCMC_dict)

    # stack the observations into a matrix
    RAM_obs = np.array([R_obs, A_obs, M_obs]).T

    fit_kwargs = {'init'  : 'advi+adapt_diag',
                  'draws' : 2000,
                  'tune'  : 2000,
                  'cores' : 1,
                  'chains': 10,
                  'random_seed' : 1234567,
                  'target_accept': 0.99,
                  'return_inferencedata': True}

    fit_type = 'Normal_multidim'

    # fit the MCMC model with "C"
    model.fit(z_obs, RAM_obs, priors=priors, fit_type=fit_type, **fit_kwargs)

    var_keys = { "f_s" : "$f_s$",
                  "C"   : "$C$",
                  "f_r" : "$f_r$",
                  "grad_A" : r"$\frac{d P}{d z}$",
                  "A_mean" : r"$\bar{A}$"}

    out_fp = f'/Users/andrewnolan/Desktop/{fit_type}_C.png'
    # plot the
    posterior.pair_plot(model.trace, var_keys=var_keys, out_fp=out_fp)


    # # fit the MCMC model with "f_i"
    # priors.pop("C", None)
    # var_keys.pop("C", None)
    #
    # # replace C with f_i
    # priors["f_i"] = {'type': 'TruncatedNormal',
    #                  'kwargs': {'mu': 8.2, 'sigma': 1.5, 'lower': 4.1}}
    #
    # # fit the MCMC model with "C"
    # model.fit(z_obs, RAM_obs, priors=priors, fit_type=fit_type, **fit_kwargs)
    #
    # out_fp = f'/Users/andrewnolan/Desktop/{fit_type}_f_i.png'
    #
    # var_keys = { "f_s" : "$f_s$",
    #              "f_i" : "$f_i$",
    #              "f_r" : "$f_r$",
    #              "grad_A" : r"$\frac{d P}{d z}$",
    #              "A_mean" : r"$\bar{A}$"}
    # # plot the
    # posterior.pair_plot(model.trace, var_keys=var_keys, out_fp=out_fp)
