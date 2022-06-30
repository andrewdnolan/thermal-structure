
import sys
import dask
import glob
import numpy as np
import xarray as xr

import open as test_open
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sys.path.append('../../src/plotting')
import animate_Zs


def calc_perc_temp(src):
    # number of valid nodes above minimum ice-thickness
    # NN = src.height.where(src.height < 10).count(dim=('nMesh_node'))#.count(dim=('coord_1', 'coord_2'))
    NN = src.height.where(src.height < 10).count(dim=('coord_1', 'coord_2'))

    # mask for determining where nodes are temperate
    mask = src['enthalpy_h'] > src['phase change enthalpy']
    # number of nodes where ethalpy exceeds phase chnage enthalpy
    NT = src.enthalpy_h.where(mask).count(dim=('coord_1', 'coord_2'))
    return NT/NN * 100

def make_colorbar(mf_dataset):
    #---------------------------------------------------------------------------
    # For Seting up the colorbar:
    #    - http://csc.ucdavis.edu/~jmahoney/matplotlib-tips.html
    #---------------------------------------------------------------------------
    cmap = cm.plasma
    norm = mcolors.Normalize(vmin=np.min(mf_dataset.Delta_MB.values),
                             vmax=np.max(mf_dataset.Delta_MB.values))

    s_map = cm.ScalarMappable(norm=norm, cmap=cmap)
    s_map.set_array(np.linspace(mf_dataset.Delta_MB.values.min(),
                                mf_dataset.Delta_MB.values.max(),
                                mf_dataset.Delta_MB.size+1))

    # If color parameters is a linspace, we can set boundaries in this way
    halfdist = (mf_dataset.Delta_MB[1] - mf_dataset.Delta_MB[0]) / 2.0
    bounds   = np.linspace(mf_dataset.Delta_MB.values.min()  - halfdist,
                           mf_dataset.Delta_MB.values.max()  + halfdist,
                           len(mf_dataset.Delta_MB) + 2)

    return cmap, norm, s_map, bounds


test_fp = '../../study/coupled_init/glc1-a/nc/glc1-a_dx_50_NT_20000_dt_0.1_MB_-1.*_OFF_prog_NetcdfOutPutSolve.nc'

file_names = sorted(glob.glob(test_fp))

parse_MB   = lambda x : float(x.split('MB_')[-1].split('_OFF')[0])


def expand_dims(ds, fp):
    offset = float(fp.split('MB_')[-1].split('_OFF')[0])
    ds     = ds.expand_dims("Delta_MB").assign_coords(Delta_MB=('Delta_MB', [offset]))
    return ds

#https://docs.xarray.dev/en/stable/user-guide/io.html#netcdf

# open the files
open_tasks    = [dask.delayed(xr.open_dataset)(f) for f in file_names]
# preprocess according to how the file was generated
preproc_tasks = [dask.delayed(test_open._preprocess)(task) for task in open_tasks]
# add parameter dim for concatenation
expand_tasks  = [dask.delayed(expand_dims)(task, f) for task, f in zip(preproc_tasks, file_names)]

datasets   = dask.compute(expand_tasks)  # get a list of xarray.Datasets


ds = xr.concat(datasets[0],"Delta_MB")

prog_vol = ds.height.isel(coord_2=-1).integrate("coord_1") /\
           ds.height.isel(coord_2=-1).isel(t=0).integrate("coord_1")
prog_perc = calc_perc_temp(ds)

cmap, norm, s_map, bounds = make_colorbar(ds)

fig, ax = plt.subplots(2,1, figsize=(8,5), sharex=True,
                       constrained_layout=True)

for offset in prog_vol.Delta_MB.values:
    color = cmap(norm(offset))

    # Plot Volume
    ax[0].plot(prog_vol.t, prog_vol.sel(Delta_MB=offset), color = color)

    # Plot % temperate
    ax[1].plot(prog_perc.t[:-1], prog_perc.sel(Delta_MB=offset)[:-1], color = color)

    # ax[1].plot(prognostic.coord_1, prognostic.sel(Delta_MB=offset).isel(t=-1, coord_2=-1).Z, color = color)

# ax[1].plot(prognostic.coord_1, prognostic.sel(Delta_MB=offset).isel(t=0, coord_2=-1).Z, color = 'k')

ax[0].axhline(1.0,c='k',ls=':',lw=1, alpha=0.5)

cbar = fig.colorbar(s_map,
                    spacing='proportional',
                    ticks=np.linspace(prog_vol.Delta_MB.values.min(),
                                      prog_vol.Delta_MB.values.max(),
                                      prog_vol.Delta_MB.size),
                    ax=ax,
                    boundaries=bounds,
                    drawedges=True,
                    format='%2.{}f'.format(3),
                    aspect=35)


# annotate the figures axes
ax[0].set_ylabel('Relative Volume per Unit Width')
ax[1].set_ylabel('Percent Temperate [\%]')
ax[1].set_xlabel('Time (yr)')
# annotate the colorbar axes
cbar.set_label('$\Delta \dot b$ (m i.e.q. yr$^{-1}$)', rotation=270, labelpad=20)
#cbar.ax.tick_params(labelsize=7)

plt.show()
plt.close()

# fig.savefig('./figs/Vol_PercTemp_2kya.png', dpi=300)


t_0 = 0
t_f = 250
mb_index = 0
Delta_mb = float(ds.Delta_MB[mb_index].values)

fig, anim = animate_Zs.animate_2D_field(ds.isel(Delta_MB=mb_index, t=slice(t_0,t_f)),
                                        field=["vel_m", "temperature"],
                                        title="$\Delta \dot b$ = {:.3f}".format(Delta_mb),
                                        cbar_label=['Velocity [m a$^{-1}$]', 'Temperature [C]'],
                                        vmin=[0, -12.5],
                                        vmax=[25, 0.0],
                                        cmap=['viridis', 'plasma'],
                                        stride=1,
                                        interval=100)

plt.close()

anim.save(f'Delta_b{Delta_mb}_{t_0}--{t_f}_anim.mp4')
