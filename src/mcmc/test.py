#!/usr/bin/env python3

import json
import pickle
import numpy as np
import arviz as az
from data import utils
from model import MCMC, PDD
from visualization import posterior

if __name__ == '__main__':
    # Load the KMR "observational" data
    z_obs, R_obs, A_obs, M_obs, B_obs, z_ELA = utils.load_KMR_DETIM()
    # Downsample the observational data, to efficency
    z_obs, R_obs, A_obs, M_obs, B_obs = z_obs[::5], R_obs[::5], A_obs[::5], M_obs[::5], B_obs[::5]
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

    var_keys = { "f_s" : "$f_s$",
                  "C"   : "$C$",
                  "f_r" : "$f_r$",
                  "grad_A" : r"$\frac{d P}{d z}$",
                  "A_mean" : r"$\bar{A}$"}

    # fit the MCMC model with "C"
    model.fit(z_obs, RAM_obs, priors=priors, fit_type=fit_type, **fit_kwargs)


    # out_fp = f'/Users/andrewnolan/Thesis/PDD_inversion/figures/pairplot_{fit_type}_C.pdf'
    # # plot the pair plot
    # posterior.pair_plot(model.trace, var_keys=var_keys, out_fp=out_fp)

    # # save the trace as a netcdf file 
    # out_fp = f'/Users/andrewnolan/Thesis/PDD_inversion/data/{fit_type}_C_trace.nc'
    # az.to_netcdf(model.trace, out_fp)


    # save the model object 
    pkl_fn = f"/Users/andrewnolan/Thesis/PDD_inversion/data/{fit_type}_C.pkl"
    # write the pickle 
    with open(pkl_fn, 'wb') as file:
        pickle.dump(model, file)




    # # fit the MCMC model with "f_i"
    # priors.pop("C", None)
    # var_keys.pop("C", None)
    
    # # replace C with f_i
    # priors["f_i"] = {'type': 'TruncatedNormal',
    #                  'kwargs': {'mu': 8.2, 'sigma': 1.5, 'lower': 4.1}}
    
    # # fit the MCMC model with "C"
    # model.fit(z_obs, RAM_obs, priors=priors, fit_type=fit_type, **fit_kwargs)
    
    # out_fp = f'/Users/andrewnolan/Thesis/PDD_inversion/figures/pairplot_{fit_type}_f_i.pdf'

    # var_keys = { "f_s" : "$f_s$",
    #              "f_i" : "$f_i$",
    #              "f_r" : "$f_r$",
    #              "grad_A" : r"$\frac{d P}{d z}$",
    #              "A_mean" : r"$\bar{A}$"}
    # # plot the
    # posterior.pair_plot(model.trace, var_keys=var_keys, out_fp=out_fp)

    # # save the trace as a netcdf file 
    # out_fp = f'/Users/andrewnolan/Thesis/PDD_inversion/data/{fit_type}_f_i_trace.nc'
    # az.to_netcdf(model.trace, out_fp)


    # # save the model object 
    # pkl_fn = f"/Users/andrewnolan/Thesis/PDD_inversion/data/{fit_type}_f_i.pkl"
    # # write the pickle 
    # with open(pkl_fn, 'wb') as file:
    #     pickle.dump(model, file)