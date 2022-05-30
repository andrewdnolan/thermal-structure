#!/usr/bin/env python3

import numpy as np
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt

# Set global matplotlib style parameters
plt.rcParams.update({'text.usetex': True,
                     'animation.html': 'jshtml',
                     'figure.facecolor': 'w',
                     'savefig.bbox':'tight'})

drive_fp = Path("/Volumes/thermal/Thesis/thermal-structure/study/coupled_init/")
src_fp   = drive_fp / "glc1-a/nc/glc1-a_dx_50_NT_20000_dt_0.1_MB_-1.500_OFF_prog.nc"


with xr.open_dataset(src_fp) as src:
    src["height"] = xr.where(src.height <= 10, 0, src.height)
    src["Z"]      = src.zbed + src.Z * src.height
    src['vel_m']  = np.sqrt(src['velocity 1']**2 + src['velocity 2']**2)


Vol = src.height.isel(coord_2=-1).integrate("coord_1") /\
      src.height.isel(coord_2=-1).isel(t=0).integrate("coord_1")


# plt.plot(src.coord_1,
#          src['surface_enthalpy'].isel(t=1).isel(coord_2=-1), 'x-')


plt.pcolormesh(src.t,
         src.coord_1/1e3,
         src['mass balance'].isel(coord_2=0).T,#.differentiate("t").T,
         shading='auto')

plt.colorbar()
plt.show()
