#!/usr/bin/env python3

import click
import xarray as xr

@click.command()
@click.option("-i", "--in_fp",
                help="filepath to NetCDF file to subsample",
                type=click.Path(exists=True), required=True)
@click.option("-o", "--out_fp",
                help="filepath to resulting NetCDF",
                type=click.Path(), required=True)
@click.option('--index', 'selection_type', flag_value='isel',
                help="use index to subsample (i.e. `isel`)",
                default=True, show_default=True)
@click.option('--value', 'selection_type', flag_value='sel',
                help="use value to subsample (i.e. `sel`)",
)
@click.option('--start', help="index start", required=True)
@click.option('--stop',  help="index stop, use -1 for final", required=True)
@click.option('--stride',help="slice stride",)

def downsample(in_fp, out_fp, selection_type, start, stop, stride):
    """ Downsample the input NetCDF file, along the 'time' (or 't') dimension.
    """
    # set the index datatypes based on indexing method:
    if selection_type == 'isel':
        idx_type = int
    else:
        idx_type = float

    print()
    print('opening input file')
    print()

    # open the input file
    with xr.open_dataset(in_fp) as src:

        print()
        print('input file is open')
        print()

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

        print()
        print('Subsetting open file')
        print()

        # after wrangling the dictionary arguments, actually do the selection
        subset = getattr(src, selection_type)(dict)

        print()
        print('Subsetting of open file DONE!')
        print()


    print()
    print('Writing new file!')
    print()
    # after the subsetting, write the file to disk as NetCDF
    subset.to_netcdf(out_fp)
    print()
    print('New file written!')
    print()
