
#############################################
#
# Attempts at quantifing e-folding timescales
#
# NOTE: this was not used, ended up using a
#       simplier approach found in the `notebooks`
#       folder.
#
#############################################
import numpy as np
import xarray as xr
from os import path
from scipy import optimize
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

def amalgamate(base_fp, beta, surge_NT=40):

    # pseudo surge netcdf file
    surge_nc = f'crmpt12_dx_50_NT_{surge_NT}_dt_0.05_MB_-0.37_OFF_Tma_-8.5_B_{beta:1.3e}_pseudo.zarr'
    surge_fp = path.join(base_fp.format('gridded'), surge_nc)
    # # second pseudo surge from recovered state
    surge2_fp = surge_fp.split('.zarr')[0] + '_C2.zarr'

    # surge recovery netcdf file
    recovery_nc = surge_nc.split('.zarr')[0] + '_NT_2000_recovery.zarr'
    recovery_fp = path.join(base_fp.format('gridded'), recovery_nc)
    # # second surge recovery netcdf file
    recovery2_nc = surge_nc.split('.zarr')[0] + '_C2_NT_2000_recovery.zarr'
    recovery2_fp = path.join(base_fp.format('gridded'), recovery2_nc)


    xarrays = []
    for i, fp in enumerate([surge_fp, recovery_fp]): # surge2_fp, recovery2_fp]):

        src = xr.open_zarr(fp)

        src['relative_volume'] = src.relative_volume * src.initial_volume

        if i == 0:
            initial_volume = src.initial_volume.copy(deep=True)

        # drop the initial volume as it's meaningless
        src = src.drop_vars('initial_volume')

        if i == 0:
            xarrays.append(src)
        else:
            # get the final timestep from previous the netcdf file
            t_f = xarrays[i-1].t.isel(t=-1)
            # extend the time coordinate
            src['t'] = src.t +  t_f

            xarrays.append(src)

    ds = xr.concat(xarrays, dim='t')

    # squeeze the x-coordinate array, which is constant in time
    ds['X'] = ds.X.isel(t=0)
    # make the volume realtive
    ds['relative_volume'] = ds.relative_volume / initial_volume
    ds['initial_volume']  = initial_volume

    return ds

def load_data(beta):
    base_fp = '/Volumes/thermal/Thesis/thermal-structure/'
    expr_fp = 'expr/02_surge2steady/result/crmpt12/{}/'
    src_fp  = path.join(base_fp, expr_fp)

    # load the given beta value
    ds = amalgamate(src_fp, beta=beta)
    # down sample to annual timeseries to constant dt=0.1
    ds = ds.interp(t=np.linspace(0.1,  4000.1, 40001), method="linear")

    # only need a subset of the variables for this
    vars = ['relative_volume', 'percent_temperate']

    # subset the variables and take the annual moving average of the data
    ds = ds[vars].rolling(t=10, min_periods=1).mean()

    ds['percent_temperate'] = ds['percent_temperate'] / 100
    return ds


beta = 1e-4

ds = load_data(beta)

fig, ax = plt.subplots(2,2, sharex=True, sharey='row',
                       figsize=(7,4),
                       constrained_layout=True)

for i in range(2):
    ax[0,i].plot(ds.t, ds.relative_volume)
    ax[1,i].plot(ds.t, ds.percent_temperate)


    ax[1,i].set_xlabel('Time [a]')

    for j in range(2):
        ax[i,j].grid()

ax[0,0].set_ylabel("$V'$ [-]")
ax[1,0].set_ylabel('$FT$ [-]')

ax[0,0].set_title(r"$f(t)= {\frac {M}{1+\left({\frac {M-f_{0}}{f_{0}}}\right)e^{-\lambda t}}}$")

ax[0,1].set_title(r"$f(t)= C + A e^{-\lambda t}\cos(\omega t-\varphi )$")


################################
# fit the logarithmic function #
################################

y1_text = [0.925, 0.45]
y2_text = [0.875, 0.35]

for i, var in enumerate(['relative_volume']): #, 'percent_temperate']):

    if var == 'relative_volume':
        peak_idx = ds[var].argmin().values
    elif var == 'percent_temperate':
        peak_idx = ds[var].argmax().values

    # asymptote value
    M = ds[var].isel(t=-1).values
    # value corresponding to peak
    y0 = ds[var].isel(t=peak_idx).values
    # time at the peak
    t_0 = ds.t.isel(t=peak_idx).values
    # end fo the timeseries
    t_f = ds.t.isel(t=-1).values

    # time vector to predict at
    t_vec = np.linspace(t_0, t_f, 1000)

    # # Make our function to be fit
    # def func(t, λ):
    #     ''' Our function to be fit'''
    #     return (y0*M)/(y0+(M-y0)*np.exp(-λ*M*t))

    def func(t, λ):
        ''' Our function to be fit'''
        return M/(1+((M-y0)/y0)*np.exp(-λ*t))

    # only fit from peak to end of timeseries
    subset = ds[var].isel(t=slice(peak_idx, -1))


    # Fit the curve with reasonable inital guesses
    popt, pcov = optimize.curve_fit(func, subset.t-t_0, subset.values, p0=[1/500])


    ax[i, 0].plot(t_vec, func(t_vec-t_0, *popt))

    ax[i, 0].axvline(t_0, c='k', lw=1)

    ax[i, 0].text(2500, y1_text[i], r"$\tau = \frac{1}{\lambda} = $" + f'{1/popt[0]:.1f}' + ' [a]')#, #textcoords='data')
    ax[i, 0].text(2500, y2_text[i], r"$\tau + t_0 = $" + f'{1/popt[0] + t_0:.1f}' + ' [a]')#, #textcoords='data')


    t_2 = t_0 + 500
    y2 = ds[var].sel(t=t_2, method='nearest').values

    T_d = ((t_2 - t_0)) * np.log(2) / np.log(y2/y0)

    print(T_d / np.log(2))
    # print(1/(popt[0]*M))

# ax[i, 0].annotate("$f_0$",
#             xy=(t_0, y0), xycoords='datta',
#             xytext=(600, 0.45), textcoords='data',
#             arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))


# ax[i, 0].annotate("$t_0$",
#             xy=(t_0, 0.35), xycoords='data',
#             xytext=(1250, 0.35), textcoords='data',
#             arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))

# ax[i, 0].annotate("$M$",
#             xy=(t_f, M), xycoords='data',
#             xytext=(3000, 0.5), textcoords='data',
#             arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))



# ax[0, 0].text(2500, 0.925, r"$\tau = \frac{1}{\lambda} = $")#, #textcoords='data')

# ################################
# # fit the damped oscillator    #
# ################################

# for i, var in enumerate(['relative_volume', 'percent_temperate']):

#     if var == 'relative_volume':
#         p0 = [-0.02, 0.003, 0.05, 2, 1.0]
#     elif var == 'percent_temperate':
#         p0 = [-25, 0.003, 0.05, 2, 55.0]

#     # time vector to predict at
#     t_vec = np.linspace(0, 4000, 1000)

#     # Make our function to be fit
#     def func(t, A, λ, ω, ϕ, C):
#         ''' Our function to be fit'''
#         return C + A * np.exp(-λ * t) * np.cos((ω)*t - ϕ)

#     # Fit the curve with reasonable inital guesses
#     popt, pcov = optimize.curve_fit(func, ds.t, ds[var].values, p0=p0)

#     ax[i, 1].plot(t_vec, func(t_vec, *popt))


#     ax[i, 1].text(2500, y1_text[i], r"$\tau = \frac{1}{\lambda} = $" + f'{1/popt[1]:.1f}' + ' [a]')#, #textcoords='data')

# plt.savefig('timescale_funcs.png', dpi=600)
plt.show()
