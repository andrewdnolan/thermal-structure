#!/usr/bin/env python3

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

SEED = 1234567



def calc_air_temp(self, z):
    """"Evaluate the surface air temperature at a given elev. for a given
        day of the year

    Inputs:
        z   (float or Nx1 ndarry) ---> Nodal surface elevation [m a.s.l.]

    Outputs:
        T   (floar or Nx365 ndarray) ---> Nodal surface elevation for each
                                          day of the yeat      [C]
    """
    # Need to seed the Generator every function call to ensure
    # all MCMC steps have the same temperature forcing
    random = np.random.default_rng(SEED)

    # Use array broadcasting to calc temp as function of elevation and doy
    doy  = np.arange(1,366)[:, np.newaxis]
    Temp = self.α * np.cos( 2*np.pi*(doy-self.T_p)/365 ) + \
           self.ΔTΔz * (self.ref_z-z) + self.T_ma + \
           random.normal(0, self.T_σ, (365,1))

    return Temp



IC_fp = '../../initialization/coarse/result/glc1-b/nc/glc1-b_3000a_dt_1_dx_50_MB_-1.675_OFF_spline_k2.nc'

with xr.open_dataset('./crmpt18-b/nc/steady_state.nc') as src:
    # correct for minimum ice thickness
    src["height"] = xr.where(src.height <= 10, 0, src.height)
    # apply sigma coordinate transform for vertical coordinate
    src["Z"]      = src.zbed + src.Z * src.height
    # Calculate the magnitude of the velocity vectors
    src['vel_m']  = np.sqrt(src['velocity 1']**2 + src['velocity 2']**2)
