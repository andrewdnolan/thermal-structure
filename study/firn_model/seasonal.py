#!/usr/bin/env python3

import numpy as np
import xarray as xr
import scipy.linalg as LA
import matplotlib.pyplot as plt

nz  = 400                     # number of vertical gridcells
H   = 20                      # vertical column thickness
k   = 2.1                     # ice conducivity               [J m^-1 K^-1 s^-1]
c   = 2009                    # specific heat capacity of ice [J kg^-1 K^-1]
ρ   = 910                     # ice density                   [kg m^-3]
spd = 60*60*24                # seconds per day
κ   = k / (c*ρ) * spd         # diffusivity
dt  = 1.0

# Annual air temp amplitude
ΔT  = 7
T_d = -13
# vertical grid
z = np.linspace(0, H, nz+1)[:,None]
dz = (H/nz+1)
# temporal grid
doy_v = np.arange(0,365,dt)[None,:]


#
A = np.sqrt(((2*np.pi)/365)/(2*κ))
# numerical solution matrix
T = np.zeros((nz+1, 365))
# analytical solution
t = T_d + ΔT*np.exp(-A*z)*np.sin(((2*np.pi)/365)*doy_v-A*z  )

# set the initial condition
T[:,0] = t[:,0]

# timestepping
for i, doy in enumerate(np.arange(0,364,dt)):

    # fourier number
    α = (κ * dt) / dz

    # Populate matrix
    A = np.diag([α     for __ in range(nz-2)],   k = 1)+\
        np.diag([1-2*α for __ in range(nz-1)], k = 0)+\
        np.diag([α     for __ in range(nz-2)],   k =-1)

    RHS     = T[1:-1,i]
    RHS[0]  = RHS[0]  + α * (T_d+ΔT*np.sin( ((2*np.pi)/365) * (doy+1)))
    RHS[-1] = RHS[-1] + α * (T_d)


    T[0,   i+1] = (T_d+ΔT*np.sin(((2*np.pi)/365)*(doy+1)) )
    T[1:-1,i+1] = np.dot(A, RHS)
    T[-1,  i+1] = T_d


    # # # enforce boundary conditions
    # A[0,0]   = 1.0
    # A[0,1]   = 0.0
    #
    # A[-2,-1] = 0.0
    # A[-1,-1] = 1.0A


fig, ax = plt.subplots()

im = ax.pcolormesh(doy_v, -z, T, shading='nearest', cmap='RdBu')

fig.colorbar(im)
# ax.plot(doy_v, ΔT*np.sin(2*np.pi*doy_v/365))
#
plt.show()
