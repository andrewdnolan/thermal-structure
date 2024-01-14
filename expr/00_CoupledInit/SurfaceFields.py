import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from thermal.plotting import make_colorbar

plt.rcParams['text.usetex'] = True

b_dot = -0.36

dir_path = "/Volumes/thermal/thermal-structure/expr/00_CoupledInit/result/crmpt12/thinned/"
filename = f'crmpt12_dx_50_NT_30000_dt_0.1_MB_{b_dot:1.2f}_OFF_Tma_*_prog.zarr'

file_path = dir_path + filename
# take the summed quantities over the final year of simulation,
# extracted along the free surface
src = xr.open_mfdataset(file_path, engine='zarr').squeeze()
#
src = src.isel(coord_2=-1, t=range(-10,0)).max('t')

fig, ax = plt.subplots(2, 1, figsize=(8,4), sharex=True,
                       constrained_layout=True)

cmap, norm, s_map, bounds = make_colorbar(src.T_ma, cmap='viridis')

for T_ma in src.T_ma:

    color = cmap(norm(T_ma))

    sub = src.sel(T_ma=T_ma)

    ax[0].plot(sub.X[::-1]/1e3, sub.runoff_frac, c=color, lw=1.0)
    ax[1].plot(sub.X[::-1]/1e3, sub.surf_melt, c=color, lw=1.0)


cbar = fig.colorbar(s_map,
                    ax=ax,
                    spacing='proportional',
                    ticks=np.linspace(src.T_ma.min(),
                                      src.T_ma.max(),10),
                                      boundaries=bounds,
                                      drawedges=True,
                                      format='%1.1f')

cbar.set_label('$T_{\\rm ma}$ at $z_{\\rm ref}$  [$^\circ$C]',
                              rotation=270, labelpad=15)
ax[1].set_xlim(0,4.5)

ax[0].set_title(f'$\Delta \dot b$ = {b_dot:1.2f} [m a$^{{-1}}$]')
ax[0].set_ylabel("Runoff Fraction [-]")
ax[1].set_ylabel("Surface Melt [m a$^{-1}$]")
ax[1].set_xlabel("Distance [km]")

plt.savefig(f'SurfaceField_b_dot_{b_dot:1.2f}.png',
            bbox_inches='tight',
            dpi=400)

plt.close()
