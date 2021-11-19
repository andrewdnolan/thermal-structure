import numpy as np
import xarray as xr
import pandas as pd
from scipy import signal
from scipy import linalg as LA
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from IPython.display import display, HTML


def fit_airtemp(ds, plot=True):
    """
    Fit the air temperature data
    """
    DOY   = np.linspace(1,366,366)
    z_ref = ds.stack(z=('x', 'y')).Elevation.mean()

    def func(d, α, temp_peak, T_mean, ΔTΔz=6.5E-3, z=2193.0, ref_z=2193.0):
        ''' Air temperature function'''
        return α*np.cos( 2*np.pi*(d-temp_peak)/365 )+ΔTΔz*(ref_z-z)+T_mean



    # In this case we assume a functional form for the curve fit that was used to generate the synthetic data.
    popt, pcov = optimize.curve_fit(func,
                                    DOY,
                                    ds.temp.groupby("time.dayofyear").mean().values,
                                    p0=[7.5, 365./1.85, -8.0])

    if plot:

        fig, ax = plt.subplots(1,1, figsize=(10,5))

        ax.plot(ds.temp.groupby("time.dayofyear").mean().dayofyear,
                ds.temp.groupby("time.dayofyear").mean().values,
                label='2007-2018 Average')

        ax.plot(DOY,
                func(DOY, *popt),
                label=r'fit: $\alpha$ = {:.1f}, $\hat{{T}}$ = {:1.0f}, $\bar{{T}}$ = {:.2f}'.format(*popt))

        ax.set_xlabel('Day of Year')
        ax.set_ylabel('$\degree$C')

        plt.legend()

        Air_Temp = {'ΔTΔz'      : {'Value': 6.5E-3,  'Tuned': False},
                    'α'         : {'Value': popt[0], 'Tuned': True},
                    'temp_peak' : {'Value': popt[1], 'Tuned': True},
                    'T_mean'    : {'Value': popt[2], 'Tuned': True} }

        display(HTML(pd.DataFrame(Air_Temp).T.to_html()))
    else:
        return popt, pcov

def net_balance_func(z, DOY, f_snow, f_ice, f_r, α, temp_peak, T_mean, A_mean,
                     grad_A, ΔTΔz=6.5E-3, ref_z=2193., T_m=273.15, T_rs=274.15):
    """Net balance function to be minimzed
    Inputs:
        z         (array) --> Nx1 array of elevations      [m a.s.l.]
        DOY       (array) --> 1xN array of Days of Year    [DOY]
        f_snow    (float) --> degree-day factor for snow   [kg m^-2 yr^-1 K^-1]
        f_ice     (float) --> degree-day factor for ice    [kg m^-2 yr^-1 K^-1]
        f_r       (float) --> refreezing factor            [-]
        α         (float) --> anual air temp. amplitude    [K]
        temp_peak  (int)  --> DOY of annual temp peak      [DOY]
        T_mean    (float) --> Mean annual air temp @ ref_z [K]
        A_mean    (float) --> Mean annual accum.   @ ref_z [kg m^-2 yr^-1]
        grad_A    (float) --> alititudinal precip factor   [m^-1]
        ΔTΔz      (float) --> air temp lapse rate          [K m^-1]
        ref_z     (float) --> reference surface elevation  [m a.s.l.]
        T_m       (float) --> Melting temp. threshold      [K]
        T_rs      (float) --> Rain versus snow threshold   [K]
    Returns:
        net_balance (float) -->
    """

    T          = α*np.cos( 2*np.pi*(DOY-temp_peak)/365 )+ΔTΔz*(ref_z-z)+T_mean
    PDDs       = np.where(T>T_m,  T-T_m,    0).sum(axis=0)
    accum_days = np.where(T<T_rs, 1/365.,   0).sum(axis=0)

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

    return MB

class optimizer:
    def __init__(self, popt):
        self.x_true = src["zs accumulation flux 2"].isel(t=-1, coord_2=-1).values
        self.z      = src.Z.isel(t=-1, coord_2=-1).values[np.newaxis, :]
        self.popt   = popt
    def forward(self,x):
        α         = self.popt[0]
        temp_peak = self.popt[1]
        T_mean    = self.popt[2]+273.15
        ΔTΔz      = 6.5E-3
        ref_z     = 2193.
        T_m       = 273.15
        T_rs      = 274.15
        DOY       = np.linspace(1,365,365)[:, np.newaxis]
        MB = net_balance_func(self.z, DOY, x[0],
                              x[1], x[2], α,
                              temp_peak, T_mean,
                              A_mean,x[3], ΔTΔz,
                              ref_z,  T_m, T_rs)
        return MB.flatten()

    def objective(self, x):
        α         = self.popt[0]
        temp_peak = self.popt[1]
        T_mean    = self.popt[2]+273.15
        ΔTΔz      = 6.5E-3
        ref_z     = 2193.
        T_m       = 273.15
        T_rs      = 275.15
        DOY       = np.linspace(1,365,365)[:, np.newaxis]
        MB = net_balance_func(self.z, DOY, x[0],
                              x[1], x[2], α,
                              temp_peak, T_mean,
                              A_mean,x[3], ΔTΔz,
                              ref_z,  T_m, T_rs)

        return LA.norm(self.x_true-MB,2)

# Simple Fitting of air temperature forcing parameters to Katie's Model runs
Young2020 = xr.open_dataset("Young_etal_2020_Delta_T_-0.9_C.nc")

z_ref = Young2020.stack(z=('x', 'y')).Elevation.mean().values

A_mean = Young2020.stack(z=('x', 'y')).Accumulation.values[
          np.argpartition(np.abs(Young2020.stack(z=('x', 'y')).Elevation.values - z_ref), 25)][:25].mean() * 910.0

with xr.open_dataset('../../../initialization/coarse/result/lilk-a/nc/lilk-a_1000a_dt_1_dx_200_MB_00.0_OFF.nc') as src:
    # correct for minimum ice thickness
    src["depth"] = xr.where(src.depth <= 10, 0, src.depth)
    # apply sigma coordinate transform for vertical coordinate
    src["Z"]     = src.zbed + src.Z * src.height
    # Calculate the magnitude of the velocity vectors
    src['vel_m'] = np.sqrt(src['velocity 1']**2 + src['velocity 2']**2)

Vol = src.height.isel(coord_2=-1).integrate("coord_1") /\
      src.height.isel(coord_2=-1).isel(t=0).integrate("coord_1")
