import numpy as np
import xarray as xr
from os import path
from tqdm.contrib import tenumerate
from thermal.derived_fields import calc_mean_enthalpy, calc_length

def mean_velocity(ds, slice, compute=False): 
    """Compute the average surface velocity [m a-1], in space and time, over a given time slice

    Inputs:
    ds (xr.Dataset) --> Full xr.Dataset which we seek to subset and calculate from
    slice  (slice)  --> `slice` object to be used by `.sel` to subset data
    compute (bool)  --> boolean, whether to load subset into memory. Useful for small time slices

    Outputs:
    vel_mu (xr.DataArray) --> Avg., in x and t, surf. vel. over the given time slice
    """

    if compute: 
        ds = ds.sel(t=slice).compute()
    else: 
        ds = ds.sel(t=slice)

    # compute the glacier length, in order to mask passive nodes
    L = calc_length(ds)
    # flip the x-coord so the glacier flows left to right
    X = ds.X.sortby('coord_1', ascending=False).compute()
    # compute the mean velocity, both space and time, at the active nodes
    vel_mu = ds.vel_m.where(X >= L).mean()
    
    return vel_mu

def amalgamate(base_fp, beta, surge_NT=40):

    # pseudo surge netcdf file
    surge_nc = f'crmpt12_dx_50_NT_{surge_NT}_dt_0.05_MB_-0.37_OFF_Tma_-8.5_B_{beta:1.3e}_pseudo.zarr'
    surge_fp = path.join(base_fp, surge_nc)
    # second pseudo surge from recovered state
    surge2_fp = surge_fp.split('.zarr')[0] + '_C2.zarr'

    # surge recovery netcdf file
    recovery_nc = surge_nc.split('.zarr')[0] + '_NT_2000_recovery.zarr'
    recovery_fp = path.join(base_fp, recovery_nc)
    # second surge recovery netcdf file
    recovery2_nc = surge_nc.split('.zarr')[0] + '_C2_NT_2000_recovery.zarr'
    recovery2_fp = path.join(base_fp, recovery2_nc)

    # empty list to store datasets in
    xarrays = []
    # loop over the four files written as part of the experiment 
    for i, fp in enumerate([surge_fp, recovery_fp, surge2_fp, recovery2_fp]):

        src = xr.open_zarr(fp)
        # redimensionalize the volume to be in m^3
        src['relative_volume'] = src.relative_volume * src.initial_volume
        # store the initial volume, to Nondimensionalize later
        if i == 0:
            initial_volume = src.initial_volume.copy(deep=True)
        # drop the initial volume now that it's redimensionalize
        src = src.drop_vars('initial_volume')  
        # if first file, just add to the list
        if i == 0:
            # add to the list
            xarrays.append(src)
        # otherwise: time is relative to start of simulation, so append as needed
        else:
            # get the final timestep from previous the netcdf file
            t_f = xarrays[i-1].t.isel(t=-1)
            # extend the time coordinate
            src['t'] = src.t +  t_f
            # add to the list
            xarrays.append(src)

    # concatenate the four files along the time dimension
    ds = xr.concat(xarrays, dim='t')
    # squeeze the x-coordinate array, which is constant in time
    ds['X'] = ds.X.isel(t=0)
    # make the volume realtive
    ds['relative_volume'] = ds.relative_volume / initial_volume
    # store the initial relative volume
    ds['initial_volume']  = initial_volume
    # calculate the mean enthalpy in kJ/kg
    ds['mean_enthalpy'] = calc_mean_enthalpy(ds) / 1e3

    return ds

def parse_beta(src): 
    """ Extract beta value from filepath and add as new dimension
    """
    # strip the file extension
    run_name = path.basename(src.encoding['source']).split('.zarr')[0]
    # extract the beta value and turn it into float
    beta = float(run_name.split('B')[-1].split('_')[1])
    # add new dimension to dataset, which we can concat along
    src = src.expand_dims('beta').assign_coords({'beta': ('beta', [beta])})
    
    return src

def main():
    # filepath to external drive
    base_fp = '/Volumes/thermal/Thesis/thermal-structure/'
    # filepath to within repo structure 
    expr_fp = 'expr/02_surge2steady/result/crmpt12/gridded/'
    # joined filepath to data directory
    src_fp  = path.join(base_fp, expr_fp)
    # write file up one directory 
    out_fp  = '/'.join(src_fp.split('/')[:-2])
    # timeseries variables
    vars = ['initial_volume', 'relative_volume', 'percent_temperate', 
            'mean_surge_vel_m', 'mean_quies_vel_m',
            'mean_enthalpy']

    # empty list to store the timeseries datasets
    timeseries = []

    for i, beta in tenumerate(np.logspace(-3, -4, 9), ascii=True):
        # concatenate all the individual files for one beta value
        src = amalgamate(src_fp, beta=beta)

        # average surface velocity in quiescence 
        src['mean_quies_vel_m'] = mean_velocity(src, slice(4000.1, 4001.0))
        # average surface velocity during the surge 
        src['mean_surge_vel_m'] = mean_velocity(src, slice(4002.1, 4004.0))

        # expand the dimension for the ith beta value
        src = parse_beta(src)
        # extract the timeseries variables of interest and load into memory
        src = src[vars].compute()
        # resample to constant sampling of dt=0.1
        src = src.interp(t=np.linspace(0.1,  8000.1, 80001), method="linear")

        # append the preprocessed, filtered, and resampled timeseries to the list
        timeseries.append(src)

    # concatenate timeseries along the beta dimension 
    surge2steady = xr.concat(timeseries, dim='beta')


    # make the output filename 
    out_fn  = 'surge2steady_timeseries.zarr'
    full_fn = path.join(out_fp, out_fn)
    # write the timeseries to disk
    surge2steady.to_zarr(full_fn, mode='a')

    print('\n\n' + '*'*100)
    print(f'Output filename: {out_fn}')
    print(f'File written to: {out_fp}')
    print('*'*100 + '\n\n')

if __name__ == '__main__': 
    main()