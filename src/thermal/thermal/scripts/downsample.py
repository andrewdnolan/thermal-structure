#!/usr/bin/env python3

import os
import click
import warnings
import numpy as np
import xarray as xr
from pandas.errors import InvalidIndexError

def find_open_method(fp): 
    """Figure out how to open source file, based on file extension
    """
    # split and check the input file extension
    __, file_ext = os.path.splitext(fp)

    # call approraite function based on file extension
    if file_ext == ".zarr":
        open_method = 'open_zarr'
    elif file_ext == ".nc":
        open_method = 'open_dataset'
    else:
        raise NotImplementedError('Only .nc and .zarr formats supported')

    return open_method

def find_save_method(fp): 
    """Figure out how to save subsampled file, based on file extension
    """
    # split and check the input file extension
    __, file_ext = os.path.splitext(fp)

    # call approraite function based on file extension
    if file_ext == ".zarr":
        save_method = 'to_zarr'
    elif file_ext == ".nc":
        save_method = 'to_netcdf'
    else:
        raise NotImplementedError('Only .nc and .zarr formats supported')

    return save_method

@click.command()
@click.option("-i", "--in_fp",
                help="filepath to NetCDF/Zarr file to subsample",
                type=click.Path(exists=True), required=True)
@click.option("-o", "--out_fp",
                help="filepath to resulting NetCDF/Zarr",
                type=click.Path(), required=True)
@click.option('--index', 'selection_type', flag_value='isel',
                help="use index to subsample (i.e. `isel`)",
                default=True, show_default=True)
@click.option('--value', 'selection_type', flag_value='sel',
                help="use value to subsample (i.e. `sel`)",)
@click.option('--start', help="index start", 
                default=0, show_default=True)
@click.option('--stop',  help="index stop, use -1 for final",
                default=-1, show_default=True)
@click.option('--stride', help="slice stride",)
@click.option('--years_worth', type=click.INT, help="Write a full years worth of data very N years or N timesteps",)


def downsample(in_fp, out_fp, selection_type, start, stop, stride, years_worth):
    """ Downsample the input NetCDF file, along the 'time' (or 't') dimension.
    """

    if years_worth and (selection_type=='isel'): 
        warnings.warn('--years_worth flag only works in "--value" (.sel) mode')
        # force selection type to value (i.e. sel)
        selection_type == 'sel'
    if years_worth and stride: 
        raise KeyError('Use either --years_worth or --stride flags, not both')


    # set the index datatypes based on indexing method:
    if selection_type == 'isel':
        idx_type = int
    else:
        idx_type = float

    # call approraite function based on file extension
    open_method = find_open_method(in_fp)
    save_method = find_save_method(out_fp)

    # open the input file with approiate method
    with getattr(xr, open_method)(in_fp, chunks={'time':'auto'}) as src:
        
        # find the time dimensions name
        if 't' in src:
            var='t'
        elif 'time' in src:
            var='time'
        elif 't' and 'time' in src:
            raise KeyError('Both "t" and "time" present, ambigous')
        else:
            raise KeyError('No valid time dimension found ("t" or "time")')

        if stride:
            # if stride provided, use a slice in the dictionary
            dict = {var : slice(*map(idx_type, [start, stop, stride]))}
        else:
            # if no stride provided, just select the start and stop timesteps
            dict = {var : list(map(idx_type, [start, stop]))}

        # (N)umber of (T)ime(s)teps (P)er (Y)ear
        NTsPY  = src[var].where(np.floor(src[var])==2, drop=True).size


        # years worth flag ONLY works with value indexing
        if years_worth: 

            # if final time value wasn't passed find the final time value in the file
            if stop == -1: 
                # get the final time value, rounded to the nearest whole number
                final = round(float(getattr(src[var], 'isel')({var : -1})), 0)
            else: 
                # if passed use passed value
                final = stop 

            # years we want all the data: 
            #   add the stride to deal with non-inclusive indexing
            years = np.arange(0, final + years_worth, years_worth)
            # create dicionary to find time values
            conditional_dict ={'cond' : np.floor(src[var]).isin(years), 'drop' : True}

            # Get the time values which fall wihtin the years we want data
            times_idx = src[var].where(**conditional_dict)

            # subsampling dictionary selecting discrete timeslices
            dict = {var : times_idx}
        
        try: 
            # after wrangling the dictionary arguments, actually do the selection
            subset = getattr(src, selection_type)(dict)

        except InvalidIndexError: 
            # Sometimes when there's blow up from CFL condition violation, the time 
            # dimension will have NaNs (i.e. really big floats). These break the indexing
            # scheme, due to repeat values (i.e. some really big floats). 
            #
            # So, we catch that error here and filter the time axis to remove duplicates, 
            # which should fix things (hopfully)

            # get the unique time indexes. ref: https://stackoverflow.com/a/51077784/10221482
            _, valid_idxs = np.unique(src[var], return_index=True)
            
            # use isel to get unique time values
            src = src.isel({var : valid_idxs})

            # after wrangling the duplicates, try do the selection again
            subset = getattr(src, selection_type)(dict)

        # need to re-chunk the subsetted data to prevent corrupting the data
        subset = subset.chunk("auto")
        
    # after the subsetting, write the file to disk 
    getattr(subset, save_method)(out_fp, mode="w")
