#!/usr/bin/env python3

import numpy as np
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

def calc_watercont(H, P):
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
