#!/usr/bin/env python3

import numpy as np
import xarray as xr
from scipy.optimize import minimize
import matplotlib.pyplot as plt


def net_balance_func(f_snow, f_ice, net_flux):
    """Net balance function to be minimzed
    Inputs:
        f_snow (float) --> degree-day factor for snow   [kg m^{-2} yr^{-1} K^{-1}]
        f_ice  (float) --> degree-day factor for ice    [kg m^{-2} yr^{-1} K^{-1}]

    Returns:
        net_balance (float) -->
    """

    z         = np.loadtxt('../../input_data/lilk-a_surf.dat')[:,1][np.newaxis,:]
    ref_z     = 2000.0  # reference surface elevation  [m a.s.l.]
    α         = 12.0    # anual air temp. amplitude    [K]
    ΔTΔz      = 6.5E-3  # air temp lapse rate          [K m^-1]
    grad_A    = 0.3332  # accum. gradient              [(kg m^{-2} yr^{-1}) m^{-1}]
    T_mean    = 265.15  # Mean annual air temp @ ref_z [K]
    A_mean    = 578.56  # Mean annual accum.  @ ref_z  [kg m^-2 yr^-1]
    temp_peak = 365./2. # DOY of annual temp peak      [doy]
    f_r       = 0.20    # refreezing factor            [-]


    DOY = np.linspace(1,365,365)[:,np.newaxis]
    T          = α*np.cos( 2*np.pi*(DOY-temp_peak)/365 )+ΔTΔz*(ref_z-z)+T_mean
    PDDs       = np.where(T>273.15, T-273.15, 0).sum(axis=0)
    accum_days = np.where(T<274.15, 1/365.,   0).sum(axis=0)

    # calculate snow accumulation
    A_snow=np.maximum((accum_days*A_mean)*(1 + (z-ref_z)*grad_A), 0.0)
    # calculate local surface melt assuming f_m = f_snow
    melt_local = PDDs * f_snow
    # calculate refreezing
    R = np.minimum(f_r*A_snow, melt_local)
    # compute the ratio b/w accumulated snow and total melt assuming f_m = f_snow
    r_s2m = np.where(melt_local==0.0, 1, A_snow / melt_local)
    # Compute nodal specific degree day factor
    f_m = np.where(r_s2m >= 1, f_snow, f_ice - (f_ice - f_snow)*r_s2m)

    # calculate surface melt [kg m^{-2} yr^{-1}] with f_m
    M_melt = f_m*PDDs
    # calculate the mass balance [m yr^{-1}]
    MB = (A_snow + R - M_melt) * (1 / 910.)

    net_balance = np.trapz(MB, z)

    # print('np.trpaz', net_balance)
    # print('np.cumsum', np.cumsum(MB, axis=1))
    return (net_balance - net_flux) / 1e3**2


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in .nc file from dignostic run
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

with xr.open_dataset("./nc/div_flux_lk.nc") as src:
    src["Z"]     = src.zbed + src.Z * src.height
    src['vel_m'] = np.sqrt(src['velocity 1']**2 + src['velocity 2']**2)

    # Calculate the Flux [m^2 a^-1]
    src["q_x"]   = src["velocity 2"].mean('coord_2') * src['height'].isel(coord_2=-1)

    # Calculate the flux Divergence (since 1-d, really only horizontal gradient)
    src["div q_x"] = src['q_x'].differentiate('coord_1')

    # intergate the flux w/ repsect to Z along free surface
    net_flux = np.trapz(src["div q_x"].isel(t=-1),
                        x = src.Z.isel(t=-1, coord_2=-1))

test = []
f_ice_vals  = np.arange(400,  1e3, 20) / 910.
f_snow_vals = np.arange(0,  87.25, 1.0) / 350.
net_bal     = np.zeros((f_ice_vals.size, f_snow_vals.size))

for i, f_ice in enumerate(f_ice_vals):
    for j, f_snow in enumerate(f_snow_vals):

        net_bal[i,j] = net_balance_func(f_snow*350., f_ice*910., net_flux)


fig, ax = plt.subplots(1,1)

im = ax.contourf(f_snow_vals, f_ice_vals, net_bal, cmap='RdBu')

ax.contour(f_snow_vals,f_ice_vals, net_bal, levels=[0], colors='k')
ax.scatter(*np.meshgrid(f_snow_vals, f_ice_vals), marker='x', s=0.1, c='w')

ax.set_ylabel('$f_{\\rm ice}$  [m yr$^{-1}$ K$^{-1}$]')
ax.set_xlabel('$f_{\\rm snow}$ [m yr$^{-1}$ K$^{-1}$]')
fig.colorbar(im, label='$\\dot b_{\\rm i} - \\nabla \\cdot \\vec{q}$ [$km^2$]')

fig.savefig('MB_tune_0.png', dpi=300, bbox_inches='tight')
plt.show()
