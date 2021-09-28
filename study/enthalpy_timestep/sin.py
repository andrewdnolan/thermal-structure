#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

z           = 2300.0    # Node Elevation            [m a.s.l.]
α           = 12.       # Annula air temp amp       [K]
ΔTΔz        = 11/2750.  # Air temp lapse rate       [K m^{-1}]
ref_z       = 0.0       # Reference Elevation       [m]
T_mean      = 273.15    # Mean annual temp at ref_z [K]
accum_days  = 0.0

# Define vector to store MB model data
DOY_v  = np.arange(1,366)
T_sin = np.zeros(365)
T_cos = np.zeros(365)

for j, t in enumerate(DOY_v):
    T_sin[j] = α*np.sin(2*np.pi*t/365) + ΔTΔz*(ref_z - z) + T_mean
    T_cos[j] = α*np.cos(2*np.pi*t/365) + ΔTΔz*(ref_z - z) + T_mean

fig, ax = plt.subplots(1,1, sharex=True, sharey=True,
                       constrained_layout=True)

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


ax.plot(DOY_v, T_sin, "k:")

for i, dt_m in enumerate([1.0, 2.0, 3.0, 6.0]):
    dt_a = dt_m / 12.0
    NT_a = 12 / dt_m

    color = colors[i]

    DOY = (np.linspace(0,1, int(NT_a))*365).astype(np.int)
    T_s = α*np.sin(2*np.pi*DOY/365) + ΔTΔz*(ref_z - z) + T_mean


    ax.plot(DOY, T_s, 'x-', c=color, label = 'dt = {} month'.format(int(dt_m)))

ax.legend()
ax.set_xlim(-5,370)
ax.set_xlabel('Day of Year')
ax.set_ylabel('Temperature [K]')

plt.show()
