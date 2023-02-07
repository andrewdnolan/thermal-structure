#!/usr/bin/env python3

import os
import sys
import json
import click
from dask.distributed import Client

# functions from local src code
from thermal.open import dataset as open_dataset

def start_cluster():
    """
    Start the dask cluster for parallel support.  
    
    Built in felxibilty to support both local and HPC side proccessing. 
    https://github.com/ualberta-rcg/python-dask/blob/master/cluster-examples/
    """
    # check if scheduler file set, i.e. on slurm cluster
    if 'SCHEDULER_FILE' in os.environ:
        # get the scheduler file
        scheduler_file = os.environ['SCHEDULER_FILE']
        # start the client
        client = Client(scheduler_file=scheduler_file)
    # otherwise on local workstation: just make use of resources avail.
    else:
        client = Client()

    return client

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
                help="filepath to NetCDF file to subsample",
                type=click.Path(exists=True), required=True)
@click.option("-o", "--out_fp",
                help="filepath to resulting NetCDF",
                type=click.Path(), required=True)
@click.option("-p", "--params",
                help="Dictionary of parameters to add to NetCDF object"\
                     "as dimmension. e.g. \'{\"T_ma\":-9.00, \"Delta_MB\":-1.5}\'",
                type=click.STRING)
def grid_data(in_fp, out_fp, params):
    """ Postprocess in the input NetCDF file written by Elmer.
    """
    # Check is parameter dictionary was passed
    if params:
        # should probably add some better error handling here
        param_dict = json.loads(params)
    else:
        # set to empty dict, loop below will not itterate
        param_dict = dict()

    # start the dask distributed cluster
    client = start_cluster()

    # print the cluster info 
    print('\n', client, '\n')
    # prevent buffered output
    sys.stdout.flush()

    # open and preprocess (i.e grid) the dataset, out of memory
    ds = open_dataset(in_fp, chunks={'time':'auto', 'nMesh_node':-1})
    # Note: need single chunk along the reshaping dimension (nMesh_node)
    # https://github.com/pydata/xarray/discussions/7217#discussioncomment-4150763

    # loop over ensemble parameter values
    for key in param_dict:
        val = param_dict[key]
        # expand the dimmension along parameter values to concatentate along
        ds = ds.expand_dims(key).assign_coords({key: (key, [val])})

    # call approraite function based on the requested file extension
    save_method = find_save_method(out_fp)

    # call approraite function and save file to disk
    getattr(ds, save_method)(out_fp, mode="w")
