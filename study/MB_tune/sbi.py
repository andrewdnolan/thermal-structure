#!/usr/bin/env python3

import numpy as np
import xarray as xr
from dataclasses import dataclass
import scipy.optimize as optimize

def fit_airtemp(ds, plot=True):
    """
    Fit the air temperature data
    """
    DOY   = np.linspace(1,366, 366)
    z_ref = ds.stack(z=('x', 'y')).Elevation.mean()

    def func(d, α, temp_peak, T_mean, ΔTΔz=6.5E-3, z=2193.0, ref_z=2193.0):
        ''' Air temperature function'''
        return α*np.cos( 2*np.pi*(d-temp_peak)/365 )+ΔTΔz*(ref_z-z)+T_mean



    # In this case we assume a functional form for the curve fit that was used to generate the synthetic data.
    popt, pcov = optimize.curve_fit(func,
                                    DOY,
                                    ds.temp.groupby("time.dayofyear").mean().values,
                                    p0=[7.5, 365./1.85, -8.0])

    parameters = {
            'ΔTΔz'      : 6.5E-3,
            'α'         : popt[0],
            'T_peak'    : int(np.round(popt[1])),
            'T_mean'    : popt[2]
            }

    return parameters

@dataclass
class AirTemp:
    """Class for storing parameter values and evaluating surface air temperature
       forcing.

    Args:
        α      (float): anual air temp. amplitude                     [K]
        T_mean (float): Mean annual air temp @ ref_z                  [K]
        T_peak (int)  : DOY of annual temp peak                       [DOY]
        ΔTΔz   (float): air temp lapse rate         (default:6.5E-3)  [K m^-1]
        ref_z  (float): reference surface elevation (default:2193)    [m a.s.l.]
    """
    α     : float
    T_mean: float
    T_peak: int
    ΔTΔz  : float = 6.5E-3
    ref_z : float = 2193.

    def eval(self, z, DOY):
        """Evaluate the surface air temperature at a given evlevation for a given
           day of the year

        Inputs:
            z   (float or Nx1 ndarry) ---> Nodal surface elevation [m a.s.l.]
            DOY (int   or 1xN ndarry) --->  Days of Year           [DOY]
        """

        T = self.α*np.cos( 2*np.pi*(DOY-self.T_peak)/365 ) + \
            self.ΔTΔz*(self.ref_z-z) + self.T_mean+273.15

        return T

@dataclass
class MassBalance:
    """Class for storing parameter values and running the Gilber et al. 2016
       surface mass balance model.

       By defaults the model is run for one year at a daily timestep, characterisitc
       of most positive degree day models.

    Inputs:
        z    (np.ndarray) --> Nx1 array of elevations      [m a.s.l.]
        Temp    (AirTemp) --> instance of AirTemp class
        f_snow    (float) --> degree-day factor for snow   [kg m^-2 yr^-1 K^-1]
        f_r       (float) --> refreezing factor            [-]
        A_mean    (float) --> Mean annual accum.   @ ref_z [kg m^-2 yr^-1]
        grad_A    (float) --> alititudinal precip factor   [m^-1]
        ref_z     (float) --> reference surface elevation  [m a.s.l.]
        T_m       (float) --> Melting temp. threshold      [K]
        T_rs      (float) --> Rain versus snow threshold   [K]

    """
    Z      : np.ndarray
    Temp   : AirTemp
    f_snow : float = None
    #f_ice  : float = None
    f_r    : float = None
    A_mean : float = None
    grad_A : float = None
    T_m    : float = 273.15
    T_rs   : float = 274.15
    ref_z  : float = 2193.
    arrays : dict  = None

    def __post_init__(self):
        pass
        # assert (len(self.Z.shape)) == 2 & (self.Z.shape[1] ==1 ):
        # print(len(self.Z.shape))

    def eval(self, **kwargs):

        # Get a list of the object attributes (i.e. model parameters)
        params = [a for a in dir(self) if not a.startswith('__')
                                       and not callable(getattr(self, a))
                                       and a not in ['Z', 'Temp', 'arrays']]

        # kwargs passed to the function (i.e. model parameters to be it)
        keys = kwargs.keys()

        # Loop over the attributes (i.e. model parameters) and make sure that
        # they are either set an initialization or are passed to the eval call
        for param in params:
            # All model parameters are converted to private attributes so that
            # parameter used for one evaluation don't affect the next
            if self.__getattribute__(param) != None:
               self.__setattr__("_"+param, self.__getattribute__(param))

            elif (self.__getattribute__(param) == None) & (param in keys):
               self.__setattr__("_"+param, kwargs[param])

            elif (self.__getattribute__(param)==None)&(param not in keys):
                raise AttributeError(
                     f"{param}: Not set at initialization and/or "+\
                      "passed to the tunning procedure"
                      )

        ########################################################################
        # Evaluate the forward model
        ########################################################################

        # Rounce et al. 2020
        self._f_ice = 2.0 * self._f_snow

        DOY = np.linspace(1,365,365,dtype=int)[:, np.newaxis]

        T          = self.Temp.eval(self.Z, DOY)

        PDDs       = np.where(T>self._T_m,  T-self._T_m, 0).sum(axis=0)
        accum_days = np.where(T<self._T_rs,     1/365., 0).sum(axis=0)

        # calculate snow accumulation
        A_snow=np.maximum((accum_days*self._A_mean)*
                          (1+(self.Z-self._ref_z)*self._grad_A), 0.0)

        # calculate local surface melt assuming f_m = f_snow
        melt_local = PDDs * self._f_snow
        # calculate refreezing
        R = np.minimum(self._f_r*A_snow, melt_local)
        # compute the ratio b/w accumulated snow and total melt assuming f_m = f_snow
        with np.errstate(divide='ignore', invalid='ignore'):
            r_s2m = np.where(melt_local==0.0, 1.0, A_snow / melt_local)

        # Compute nodal specific degree day factor
        f_m = np.where(r_s2m >= 1, self._f_snow,
                       self._f_ice - (self._f_ice - self._f_snow)*r_s2m)
        # calculate surface melt [kg m^{-2} yr^{-1}] with f_m
        M_melt = f_m*PDDs
        # calculate the mass balance [m yr^{-1}]
        MB = (A_snow + R - M_melt) * (1 / 910.)

        self.arrays = dict(MB = MB, A = A_snow, M = M_melt, R=R)
        return MB

if __name__ == "__main__":

    # Read in Katies model runs
    Young2020 = xr.open_dataset("study/MB_tune/notebooks/Young_etal_2020_Delta_T_-0.9_C.nc")
    # Fit the air temperature function to Katies results
    AirTemp_params = fit_airtemp(Young2020)
    # instantiate the AirTemp class using the parameters from the curve fitting
    Temp = AirTemp(**AirTemp_params)

    # Read in the lilk model results
    with xr.open_dataset('initialization/coarse/result/lilk-a/nc/lilk-a_1000a_dt_1_dx_200_MB_00.0_OFF.nc') as src:
        # correct for minimum ice thickness
        src["depth"] = xr.where(src.depth <= 10, 0, src.depth)
        # apply sigma coordinate transform for vertical coordinate
        src["Z"]     = src.zbed + src.Z * src.height
        # Calculate the magnitude of the velocity vectors
        src['vel_m'] = np.sqrt(src['velocity 1']**2 + src['velocity 2']**2)

    # Surface Elevation Vector
    Z      = src.Z.isel(t=-1, coord_2=-1).values[np.newaxis, :]
    z_ref  = Young2020.stack(z=('x', 'y')).Elevation.mean().values
    A_mean = Young2020.stack(z=('x', 'y')).Accumulation.values[
              np.argpartition(np.abs(Young2020.stack(z=('x', 'y')).Elevation.values - z_ref), 25)][:25].mean() * 910.0

    with open("dict.pickle","rb") as src:
        params = pickle.load(src)

    Z      = params["Z"]
    z_ref  = params["z_ref"]
    A_mean = params["A_mean"]
    # instantiate the MassBalance class
    model = MassBalance(Z=Z, Temp=Temp, A_mean=A_mean, f_r=0.6, grad_A=3.713e-3)

    def simulation_wrapper(params):
        f_snow  = params[0]
        f_ice   = params[1]
        return model.eval(f_ice=f_ice, f_snow=f_snow)
