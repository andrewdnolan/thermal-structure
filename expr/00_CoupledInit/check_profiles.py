import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

paths = ["/Users/andrewnolan/Thesis/SFU_Thesis/chapters/04_SteadyState/data/",
         "/Users/andrewnolan/Thesis/thermal-structure/expr/00_CoupledInit/"]

fig, ax = plt.subplots(constrained_layout=True)

for path in paths:
    filepath= path + 'crmpt12_dx50_nz15_Z.nc'

    src = xr.open_dataset(filepath).isel(t=-1, coord_2=-1, T_ma=0, offset=-1)
    print(src.offset, src.T_ma)
    ax.plot(src.X, src.Z)

plt.show()
