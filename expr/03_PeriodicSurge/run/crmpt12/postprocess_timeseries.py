import click
import numpy as np
import xarray as xr
from os import path
from glob import glob
from tqdm.contrib.itertools import product
from thermal.derived_fields import (calc_length, 
                                    calc_volume,
                                    Variable_at_ELA, 
                                    calc_mean_enthalpy)

def calc_mean_velocity(ds, loc='surface'): 
    """Compute the spatial average of the velocity [m a-1], at valid horizontal nodes

    Inputs:
    ds (xr.Dataset) --> Full xr.Dataset which we seek to subset and calculate from
    loc  (str)      --> ['surface'|'bed'] free surface to extract velocity along

    Outputs:
    vel_mu (xr.DataArray) --> Avg. in x vel. along the specified surface
    """
    
    # mapping of valid "locations" to extract velocity and associated vertical indexed
    locations = dict(surface=-1, bed=0)

    # make sure a valid location was specified 
    if loc not in locations:
        raise ValueError(f"`loc` must be one of {list(locations)}" )
    else: 
    # if so, extract the associated index
        idx = locations[loc]

    # compute the glacier length, in order to mask passive nodes
    # returned in km but convert to m to match `coord_1`
    L = calc_length(ds) * 1e3
    # flip the x-coord so the glacier flows left to right
    X = ds.X.sortby('coord_1', ascending=False).compute()
    # compute spatial average of velocity along the specified surface, at the active nodes
    vel_mu = ds.vel_m.where(X <= L).isel(coord_2=idx).mean('coord_1')
    
    return vel_mu


def calc_flux_at_ELA(src): 
    """ Calculate timeseries of flux [m2 a-1] across the ELA 
    """
    # ice thickness at the ELA [m]
    H_ELA = Variable_at_ELA(src, 'height').isel(coord_2=-1)
    # depth averaged velocity [m a-1]
    ubar_ELA = Variable_at_ELA(src, 'vel_m').mean('coord_2')

    # flux at ELA [m2 a-1], Equation (8.1) of Cuffey and Paterson
    return H_ELA*ubar_ELA

def parse_params(fn): 
    """ Extract beta value and surge period from filepath
    """
    # strip the file extension
    run_name = path.basename(fn).split('.zarr')[0]
    # extract the beta value and turn it into float
    beta = float(run_name.split('B')[-1].split('_')[1])
    # extract the surge recurrence, add 2 year surge period to get surge cycle
    SP = float(run_name.split('QP')[-1].split('_')[1]) + 2 
    return np.array([beta, SP])

def expand_dims(src, SP, beta): 
    # add new dimensions to dataset, which we can concat along
    src = src.expand_dims('SP').assign_coords({'SP': ('SP', [SP])})
    src = src.expand_dims('beta').assign_coords({'beta': ('beta', [beta])})  
    return src

@click.command()
@click.option('--gridded_dir', 
              help='Filepath to the the `gridded` folder containing the zarr file to postprocess',
              type=click.Path(exists=True), required=True)
def main(gridded_dir):
    # runname template to search, most params are hard coded but could be pass over cli 
    run_name = "crmpt12_dx_50_TT_*--*ka_MB_-0.35_OFF_Tma_-8.5_B_{beta}_SP_2_QP_{QP}.zarr"

    # source filepath  
    src_fp  =  path.join(gridded_dir, run_name)
    # output filepath, where results will be written up one directory 
    out_fp  = '/'.join(gridded_dir.split('/')[:-3])

    print(out_fp)
    # timeseries variables
    vars = ['initial_volume', 'relative_volume', 'percent_temperate', 'h_ELA',
            'mean_enthalpy',  'flux_across_ELA', 'u_bar_surface', 'u_bar_bed']

    # empty list to store the timeseries datasets
    timeseries = []

    # find all the files in the specified directory
    files = sorted(glob(src_fp.format(beta='*',QP='*')))

    # parse the parameter combinations from the file names
    params = np.vstack(list(map(parse_params,files)))

    betas   = np.unique(params[:,0])
    periods = np.unique(params[:,1])
    
    for beta, SP in product(betas, periods):

        # quiescent period is surge period minus two year active phase
        QP = SP-2

        # check that the beta/QP combination exists
        # specifally, bisected beta values were not run for the 15 & 60 year sugres. 
        found = glob(src_fp.format(beta=f'{beta:.3e}', QP=int(QP)))
        # So, skip the bisected beta values for those surge periods
        if not found: 
            continue

        # open the ith file
        src = xr.open_mfdataset(found, engine='zarr', data_vars='different')

        ##############################################################
        # Calculate diagnostic variables for analysis 
        ##############################################################
        # recompute the volume based on the concatenated timeseires
        volume = calc_volume(src)
        # write the initial volume in m^2
        src['initial_volume']  = volume.isel(t=0)
        # write the time dependent "relative volume"
        src['relative_volume'] = volume / volume.isel(t=0)    

        # calculate the area weighted average enthalpy (kJ kg-3)
        src['mean_enthalpy'] = calc_mean_enthalpy(src) / 1e3
        # calculate the flux [m2 a-1] at the ELA as a function of time
        src['flux_across_ELA'] = calc_flux_at_ELA(src)
        # calculate the ice thickness at the ELA [m]
        src['h_ELA'] = Variable_at_ELA(src, 'height').isel(coord_2=-1)
        # extract spatial average of the surface velocity [m a-1]
        src['u_bar_surface'] = calc_mean_velocity(src, loc='surface')
        # extract spatial average of the basal velocity [m a-1]
        src['u_bar_bed'] = calc_mean_velocity(src, loc='bed')

        # extract the timeseries variables of interest and load into memory
        src = src[vars]
        # resample to constant sampling of dt=0.1
        src = src.interp(t=np.linspace(0.1,  9000.1, 90001), method="linear")

        # expand the dimension for the ith beta value
        src = expand_dims(src, SP, beta)

        # append the preprocessed, filtered, and resampled timeseries to the list
        timeseries.append(src)

    # combine the timeseries along the beta and SP dims 
    periodic_surges = xr.combine_by_coords(timeseries)

    # make the output filename 
    out_fn  = 'PeriodicSurges_Timeseries_9ka.zarr'
    full_fn = path.join(out_fp, out_fn)
    # write the timeseries to disk
    periodic_surges.to_zarr(full_fn, mode='a')

    print('\n\n' + '*'*100)
    print(f'Output filename: {out_fn}')
    print(f'File written to: {out_fp}')
    print('*'*100 + '\n\n')

if __name__ == '__main__': 
    main()
