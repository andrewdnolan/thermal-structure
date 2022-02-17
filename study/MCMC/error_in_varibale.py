#!/usr/bin/env python3

import PDD
import theano
import warnings
import numpy as np
import pymc3 as pm
import arviz as az
import xarray as xr
import theano.tensor as tt
import matplotlib.pyplot as plt

warnings.filterwarnings(action='ignore')
theano.config.compute_test_value = 'warn'

class PDD_MCMC:
    def __init__(self, α, T_ma, ΔTΔz, T_p, ref_z, T_m, T_rs, T_σ):
        # These a arguments which all models will need. Any of the model
        # pamaeters which are dependent on the model formulation are passed
        # as kwargs.
        self.α      = α       # anual air temp. amplitude    [K]
        self.T_ma   = T_ma    # Mean annual air temp @ ref_z [K]
        self.ΔTΔz   = ΔTΔz    # air temp lapse rate          [K m^-1]
        self.T_p    = T_p     # DOY of annual temp peak      [DOY]
        self.ref_z  = ref_z   # reference surface elevation  [m a.s.l.]
        self.T_m    = T_m     # Melting temp. threshold      [K]
        self.T_rs   = T_rs    # T_
        self.T_σ    = T_σ

    def _air_temp(self, z):
        """"Evaluate the surface air temperature at a given evlevation for a given
            day of the year

        Inputs:
            z   (float or Nx1 ndarry) ---> Nodal surface elevation [m a.s.l.]

        Outputs:
            T   (floar or Nx365 ndarray) ---> Nodal surface elevation for each
                                              day of the yeat      [C]
        """
        doy  = np.arange(1,366)[:, np.newaxis]
        Temp = self.α*np.cos( 2*np.pi*(doy-self.T_p)/365 ) + \
               self.ΔTΔz*(self.ref_z-z)+self.T_ma + \
               np.random.normal(0, self.T_σ, (365,1))

        return Temp

    def forward(self, z, f_snow, C, f_r, grad_a, A_mean):
        """Theano implementation of the forward model which supports shared
           variables as input.

           This is an intermidiate function which returns a list of the three
           components of the PDD model (i.e. Accumulation, Refreezing and melt).

        Inputs:
            z    (ndarray) --> Nx1 array of elevations         [m a.s.l.]
            f_snow (float) --> degree-day factor for snow      [kg m^-2 yr^-1 K^-1]
            C      (float) --> factor relating f_ice to f_snow [-]
                                1 <= C < ?
            f_r    (float) --> refreezing factor               [-]
            **kwargs
        Outputs:
            [A_snow, R, M_melt]  ([theano.tt, theano.tt, theano.tt]) -->
                List of theano tensors without numeric values for each component
                of the mass balance model [m i.e yr^-1]
        """

        f_ice  = C*f_snow

        # temperature and PDDs calc
        T      = self._air_temp(z)
        PDDs   = tt.switch(tt.gt(T, self.T_m), T, 0.0).sum(axis=0)

        # # accumulation calc
        # A_snow = self.__tt_accumulation(z, grad_a, T)

        A_days = tt.switch(tt.lt(T, self.T_rs), 1/365., 0.0).sum(axis=0)
        A_snow = tt.maximum((A_days*A_mean)*(1+(z-self.ref_z)*grad_a), 0.0)

        # calculate local surface melt assuming f_m = f_snow
        melt_local = PDDs * f_snow

        # calculate refreezing
        R = tt.minimum(f_r*A_snow, melt_local)


        r_s2m = tt.switch(tt.eq(melt_local, 0.0), 1.0, A_snow/melt_local)

        f_m = tt.switch(tt.ge(r_s2m, 1.), f_snow,
                        f_ice - (f_ice - f_snow)*r_s2m)

        # calculate surface melt [kg m^{-2} yr^{-1}] with f_m
        M_melt = f_m*PDDs

        # Return individual components of the mass balance model in [m i.e. / y]
        return R * (1/910), A_snow * (1/910), M_melt * (1/910)

def read_observation():
    """Read and return EMY/KMR results
    """

    # file paths to the NetCDF files
    nc_fp = "../../input_data/mass_balance/Kaskawulsh_NetBalance.nc"

    with xr.open_dataset(nc_fp) as MB_obs:
        # average over the 25 reference model runs and stack along elevation
        stacked = MB_obs.stack(z=('x','y'))

        # Sort the indexes by elevation
        idxs  = stacked.dropna('z').z.values
        elev  = stacked.Z.dropna('z').values
        idxs  = idxs[np.argsort(elev)]

        # Calculate the Kaskawulsh ELA
        ELA_idxs = np.argpartition(np.abs(stacked.B.values), 5)
        z_ELA    = float(stacked.isel(z=ELA_idxs[:5]).Z.mean())

        # surface elevation of observed values
        z_obs = stacked.Z.sel(z=idxs).values
        # Observed Refreezing   [m  i.e. / yr]
        R_obs = stacked.R.sel(z=idxs).values
        # Observed Accumulation [ m i.e. / yr ]
        A_obs = stacked.A.sel(z=idxs).values
        # Observed melt         [ m i.e. / yr ]
        M_obs = stacked.M.sel(z=idxs).values
        # Observed mass balance [m i.e. / yr ]
        B_obs = stacked.B.sel(z=idxs).values

    return z_obs, R_obs, A_obs, M_obs, B_obs

def plot_posterior(trace, vars: list, out_fp: str):
    """Wrapper for arviz trace plot
    """
    _ = az.plot_trace(trace, var_names=vars)
    plt.tight_layout()
    plt.savefig(out_fp, dpi=300)

def plot_BAMR(z: np.ndarray, ppc: tuple, obs: tuple, out_fp: str, axis_labels=None):
    """Plot net (B)alance, (A)ccumulation, and (M)elt (BAM).
    """

    if axis_labels == None:
        axis_labels = ('Refreezing (m/yr)','Accumulation (m/yr)','Melt (m/yr)', 'Mass Balance (m/yr)')

    labels = (r'$X_2(z, \theta_3)$',
              r'$X_1(z, \theta_1, \theta_2)$',
              r'$X_3(z, \theta_4, \theta_5)$',
              r'$Y(z, \theta_1, \theta_2, \theta_3, \theta_4,  \theta_5)$')


    fig, ax = plt.subplots(4, 1, figsize=(5,7), sharex=True)

    for i, (pp, obs) in enumerate(zip(ppc, obs)):

        ax[i].scatter(z, obs, s=1.0, c='#7570b3')

        ax[i].plot(z, np.mean(pp, axis=(0, 1)), lw=1.0, color='#1b9e77',
                   label=labels[i])

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

if __name__ == '__main__':

    fit_type = 'new'
    draws = 2000
    tune = 2000
    cores = 2


    # load observations
    z_obs, R_obs, A_obs, M_obs, B_obs = read_observation()

    # dictionary of PDD model parameters
    p = [ 8.29376332e-05, -3.45256005e-02,  6.31076200e+00]
    doy = np.arange(1,366)
    fancy_std = np.polyval(p, doy)[:, np.newaxis]

    const =  dict(T_m = 0.0, T_rs = 1.0, α = 10.8, T_ma = -9.02, ΔTΔz  = 6.5E-3,
                  T_p = 196, ref_z = 2193, T_σ = fancy_std)

    # initialize the PDD melt model class
    PDD_forward = PDD_MCMC(**const)

    # Define Priors
    with pm.Model() as model:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ----> Mass balance Model (physical priors)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Degree Day factor: Braithwaite (2008) / Rounce et al. 2020
        f_s_prior  = pm.TruncatedNormal("f_s",    mu=4.1,  sigma=1.5,  lower=0.0)
        # Somewhat base of Aschwanden et al. 2019
        C_prior    = pm.TruncatedNormal("C",      mu=2.0,  sigma=0.75,  lower=1)
        # Need a reference for the accumulation grad distibution
        grad_a     = pm.Uniform("grad_a", lower=0.5e-4, upper=5e-4)
        #
        A_m_prior  = pm.Normal("A_mean", mu=1300, sigma=100)
        # Somewhat base of Aschwanden et al. 2019
        f_r_prior  = pm.TruncatedNormal('f_r',    mu=0.5,  sigma=0.2,  lower=0.0, upper=1)

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ----> Hyperparameters (likelihood related priors)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        R_sigma = pm.HalfCauchy("R_sigma", 3)
        A_sigma = pm.HalfCauchy("A_sigma", 3)
        M_sigma = pm.HalfCauchy("M_sigma", 3)
        B_sigma = pm.HalfCauchy("B_sigma", 3)

    # Define Forward model (wrapped through theano)
    with model:
        R, A, M = PDD_forward.forward(z_obs, f_s_prior, C_prior, f_r_prior, grad_a, A_m_prior)
        # net balance [m i.e. / yr ]
        B = A + R - M

    with model:
        # Individual likelihood functions for each component
        R_like = pm.Normal("R_like", mu=R, sigma=R_sigma, observed=R_obs)
        A_like = pm.Normal("A_like", mu=A, sigma=A_sigma, observed=A_obs)
        M_like = pm.Normal("M_like", mu=M, sigma=M_sigma, observed=M_obs)
        # likelihood uses latent x, not observed
        B_like = pm.Normal("B_like", mu=A + R - M, sigma=B_sigma, observed=B_obs)

    
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
                                             var_names=["R_like", "A_like", "M_like", "B_like"],
                                             keep_size=True)
        pp_R = ppc["R_like"]
        pp_A = ppc["A_like"]
        pp_M = ppc["M_like"]
        pp_B = ppc["B_like"]


    # plot the traces
    vars   = ['f_s', 'C', 'f_r', 'grad_a', 'A_mean']
    out_fp = f'./result/error_in_var/trace_LA_{fit_type}.png'

    plot_posterior(trace, vars=vars, out_fp=out_fp)

    # plot the mass balance results and components
    out_fp = f'./result/error_in_var/BAM_LA_{fit_type}.png'


    plot_BAMR(z_obs, (pp_R, pp_A, pp_M, pp_B), (R_obs, A_obs, M_obs, B_obs), out_fp=out_fp)

    # save the trace
    out_fp = f'./result/error_in_var/trace_LA_{fit_type}.nc'
    az.to_netcdf(trace, out_fp)
