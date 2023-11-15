import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from thermal.plotting import enthalpy_pcolormesh
plt.rcParams['text.usetex'] = True


def Plot_StrainHeatingComparison(b_dot, T_ma):

    if b_dot == -0.36:
        f_dd = '_fdd_2.0_'
        degree_day_label = '$f_{\\rm dd} = 2.0$ mm d$^{-1}$ K$^{-1}$'
    else:
        f_dd = '_'
        degree_day_label = '$f_{\\rm dd} = 6.0$ mm d$^{-1}$ K$^{-1}$ (Reference)'

    run_names = ['crmpt12_dx_50_NT_30000_dt_0.1_MB_{b_dot:1.2f}_OFF_Tma_{T_ma:1.1f}_prog{f_dd}NoStrainHeating',
                 'crmpt12_dx_50_NT_30000_dt_0.1_MB_{b_dot:1.2f}_OFF_Tma_{T_ma:1.1f}_prog{f_dd}WithStrainHeating']

    labels = ['Without Strain Heating', 'With Strain Heating']

    fig, ax = plt.subplots(2, 1, figsize=(6,4), sharex=True, sharey=True,
                           constrained_layout=True)

    for i, run in enumerate(run_names):
        run_name = run.format(b_dot=b_dot, T_ma=T_ma, f_dd=f_dd)
        fp = f'result/crmpt12/gridded/{run_name}.zarr'

        src = xr.open_zarr(fp).squeeze()

        # plot the final enthalpy field
        im = enthalpy_pcolormesh(src, -1, axes=ax[i])

        # label the subplots
        ax[i].text(0.55, 0.9, s=labels[i] + ',   ' + degree_day_label,
                   va='center', ha='center', transform=ax[i].transAxes)

        ax[i].set_ylabel('Elevation [m a.s.l.]')
    ax[i].set_xlabel('Distance [km]')

    cbar = fig.colorbar(im, ax=ax,)
    # special ticks for enthalpy
    cbar.set_ticks(np.concatenate((np.linspace(-8, 0, 9),
                                   np.linspace(0.1, 0.5, 5))))

    cbar.set_label(r'Water Content ($\%$) / Temp. Relative to PMP ($^\circ$C)',
                   rotation=270, labelpad=12.5)

    fig.savefig(f'StrainHeatingComparison_OFF_{b_dot:1.2f}_Tma_{T_ma:1.1f}.png', dpi=400)

if __name__ == '__main__':

    for param in [(-0.36,-8.5), (-0.50,-9.0)]:
        Plot_StrainHeatingComparison(*param)

