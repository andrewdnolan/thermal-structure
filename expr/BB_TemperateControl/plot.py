import xarray as xr
import matplotlib.pyplot as plt
from thermal.plotting import make_colorbar

plt.rcParams['text.usetex'] = False

src = xr.open_mfdataset('result/crmpt12/gridded/*.nc', parallel=True)

fig, ax = plt.subplots(figsize=(6, 3), constrained_layout=True)

cmap, norm, s_map, bounds = make_colorbar(src.offset)

ax.axhline(1.0, ls=':', lw=1.0, c='k')
for off in src.offset:
    color=cmap(norm(off))

    ax.plot(src.t, src.sel(offset=off).relative_volume, c=color, lw=1.0)


cbar = fig.colorbar(s_map,
                spacing='proportional',
                ticks=src.offset,
                ax=ax,
                boundaries=bounds,
                drawedges=True,
                format='%0.{}f'.format(3))


ax.set_ylabel('Relative Volume per Unit Width [-]')
ax.set_xlabel('Time [a]')
# annotate the colorbar axes
cbar.set_label('$\Delta \dot b$ [m a$^{-1}$]', rotation=270, labelpad=20)

plt.show()