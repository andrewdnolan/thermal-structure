#!/usr/bin/env python3

import numpy as np
import xarray as xr
from constants import *

def calc_T_prime(P):
    """[K]
    """
    return T_ptr - β * (P - P_ptr)

def calc_enthalpy(T, ω=0.0):
    """[J kg-1]
    """
    return (CapA / 2.0 * (T**2 - T_ref**2) + CapB * (T - T_ref)) + ω * L_heat

def calc_H_fusion(P):
    """[J kg-1]
    """
    return calc_enthalpy(T=calc_T_prime(P), ω = 0.0)

def calc_watercont(H, P=P_surf):
    """[-]
    """
    # water content [-]
    ω = (H - calc_H_fusion(P)) / L_heat
    # water content only defineed where H is greater than H of fusion
    return np.where(ω < 0, 0.0, ω )

def calc_Temp(H, P):
    """[K]
    """
    # Temperature [K]
    T = (-CapB + (CapB**2 + 2 * CapA * ((CapA / 2) * \
          T_ref**2 + CapB * T_ref + H))**(1/2)) / CapA
    # Temperature relative to pressure melting
    T_prime = calc_T_prime(P)
    # temperature only defineed where H is less than H of fusion
    return np.where(T < T_prime, T_prime, T )

class surface_AirTemp:

    def __init__(self, α = 10.8, dTdz = 6.5E-3, z_ref = 2193.0, T_mean = -9.02,
                 T_peak = 196, T_σ = (8.2938e-05, -3.45256e-02, 6.3108e+00),
                 SEED=12345678
                 ):
        """
        """


        self.α       = α       #
        self.dTdz    = dTdz    #
        self.z_ref   = z_ref
        self.T_mean  = T_mean
        self.T_peak  = T_peak
        self.T_σ     = T_σ
        self.SEED    = SEED


    def __call__(self, z, **kwargs):
        """"Evaluate the surface air temperature at a given elev. for all 365
            days of the year
        Inputs:
            z   (float or Nx1 ndarry) ---> Nodal surface elevation [m a.s.l.]
        Outputs:
            T   (floar or Nx365 ndarray) ---> Nodal surface elevation for each
                                              day of the yeat      [C]
        """

        # allow attributes to be overwritten by passing them to call
        if kwargs:
            for kwarg in kwargs:
                if self.__getattribute__(kwarg) != None:
                    self.__setattr__(kwarg, kwargs[kwarg])


        # Need to seed the Generator every function call to ensure same temp.
        random = np.random.default_rng(self.SEED)

        if hasattr(self.T_σ, '__len__'):
            std = np.polyval(self.T_σ, np.arange(0,365)[:, np.newaxis])
        else:
            std = self.T_σ

        # Use array broadcasting to calc temp as function of elevation and doy
        doy  = np.arange(0,365)[:, np.newaxis]
        Temp = self.α * np.cos( 2*np.pi*(doy-self.T_peak)/365 ) + \
               self.dTdz * (self.z_ref-z) + self.T_mean + \
               random.normal(0, std, (365,1))

        return Temp

class heat_pump:
    """Implentation of Dirichlet surface boundary condition for Enthalpy Eqn.

    """

    def __init__(self, h_aq=3.0, f_dd=0.003, r_frac=0.3, w_max_aq=0.1,
                 w_max_en=0.03, T_melt=0.0, c_firn=0.05):

        self.h_aq     = h_aq         # [m]
        self.f_dd     = f_dd         # [m K-1 d-1]
        self.r_frac   = r_frac       # [-]
        self.w_max_aq = w_max_aq     # [-]
        self.w_max_en = w_max_en     # [-]
        self.T_melt   = T_melt       # [C]
        self.c_firn   = c_firn       # [m-1]

    def __create_xarray(self, z, mb, Surf_E, T_surf, doy_ip1):

        # create a new dataset to populate
        ds = xr.Dataset(
            {
                "Enthalpy"    : (["time", "c1"], Surf_E),
                "Temperature" : (["time", "c1"], T_surf),
                "Water_Cont." : (["time", "c1"], calc_watercont(Surf_E)),
                "Mass_Balance": (["c1"], mb)
            },
            coords={
                "z": (["c1"], z),
                "time": [doy_ip1],
            }
        )

        return ds

    def __call__(self, z, mb, airtemp, doy_i=0, dt=365, airtemp_kwarg={}):

        def clavo_greeve_PDDs(T, σ=0.0):
            import scipy.special as sp

            with np.errstate(divide='ignore', invalid='ignore'):
                T_norm = T / (np.sqrt(2)*σ)

            clavo_greeve = σ / np.sqrt(2*np.pi) * np.exp(-T_norm**2) + T/2*sp.erfc(-T_norm)

            return clavo_greeve

        assert z.size == mb.size

        doy_i    = doy_i
        doy_ip1  = doy_i+dt
        doy  = np.arange(0,365)


        NX = z.size                # number of x gridcells
        DT = doy_ip1 - doy_i       # number of timesteps
        NT = 1

        # PDDs   = np.zeros((NT,NX))    # [C d]
        Q_lat  = np.zeros((NT,NX))    # [J kg-1]
        H_max  = np.zeros((NT,NX))    # [J kg-1]
        T_surf = np.zeros((NT,NX))    # [K]
        Surf_E = np.zeros((NT,NX))    # [J kg-1]
        H_surf = np.zeros((NT,NX))    # [J kg-1]

        σ    = airtemp.T_σ
        T    = airtemp(z, **airtemp_kwarg)                  # [C]
        # PDDs = np.maximum(T - self.T_melt, self.T_melt) * 1 # [C d]
        PDDs = clavo_greeve_PDDs(T, σ=np.polyval(σ, doy)[:,None])

        # Calculte melt in meters of firn equivalent
        melt = self.f_dd * PDDs[doy_i:doy_ip1].sum(axis=0) #* rho_i/rho_s
        # broadcast mass balance array to (NT, NX) matrix
        Q_lat = np.where(np.ones((NT,1))*mb >= 0.0,
                        (1 - self.r_frac) * (rho_w/self.h_aq) * L_heat * melt,
                         0.0) # [J m-3]

        # maxium water content [-] based surface density
        H_max = np.where(np.ones((NT,1))*mb   >= 0.0,
                         calc_H_fusion(P_surf) + self.w_max_aq * L_heat,
                         calc_H_fusion(P_surf) + self.w_max_en * L_heat)

        # surface air temp in kelvin, fixed at 273.15
        T_surf = np.where(T > 0.0, 273.15, T+273.15) # [K]
        # Convert surface air temp to enthalpy
        H_surf = calc_enthalpy(T_surf[doy_i:doy_ip1].mean(axis=0))[None,:]  # [J kg-1]
        # Add the surface heating to the surface enthalpy
        H_surf+= Q_lat / rho_s                         # [J kg-1]
        # Limit the surface enthalpy based on max water content
        Surf_E = np.where(H_surf<H_max, H_surf, H_max)


        ds = self.__create_xarray(z, mb,
                                  Surf_E,
                                  T_surf[doy_i:doy_ip1].mean(axis=0)[None,:],
                                  doy_ip1)

        return ds
