#!/usr/bin/env python3

"""
Constants
"""

# seconds per year                [s a^{-1}]
spy = 365.25 * 24. * 60. * 60.
# ice density                     [kg m^{-3}]
rho_i = 910.0
# snow/surface density            [kg m^{-3}]
rho_s = 350.0
# water  density                  [kg m^{-3}]
rho_w = 1000.0
# Ref. temp. for Enthalpy calc    [K]
T_ref  = 200.0
# Surface atmospheric pressure    [MPa]
P_surf = 0.1013
# Triple point pressure for water [MPa]
P_ptr  = 6.1173e-4
# Temp. @ triple point of water   [K] (Cuffey and Paterson 2010)
T_ptr  = 273.16
# Press. melting slope            [K MPa-1] (Cuffey and Paterson 2010)
β      = 0.098
# Latent heat of fusion           [J kg-1]  (Cuffey and Paterson 2010)
L_heat = 333500.0
#  Cp(T) = A*T + B                [J kg-1 K-2] (Cuffey and Paterson 2010)
CapA   = 7.122
#  Cp(T) = A*T + B                [J kg-1 K-1] (Cuffey and Paterson 2010)
CapB   = 152.5
