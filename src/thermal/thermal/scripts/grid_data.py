#!/usr/bin/env python3

import click
import dask
import json
import os
import sys

from dask.distributed import Client, LocalCluster

from thermal.open import dataset as open_dataset

# custom dask config needed for scipt to execute on compute canada systems
cedar_configs = {
    "distributed.worker.memory.pause": False,
    "distributed.worker.memory.spill": False,
    "distributed.worker.memory.target": False,
    "distributed.worker.memory.terminate": False, 
    "distributed.nanny.pre-spawn-environ.MALLOC_TRIM_THRESHOLD_": 0,
}

def start_cluster():
    """
    Start a LocalCuster for parallel support processing

    LocalCluster will work for both a laptop and SLURM job, assuming it's a 
    a shared memory job (i.e. single node).
    """
    
    timeout = 600

    if "SLURM_JOB_ID" in os.environ:
        tmpdir = os.environ["SLURM_TMPDIR"]
        n_cpus = int(os.environ["SLURM_CPUS_PER_TASK"])

        if "THREADS_PER_WORKER" in os.environ:
            threads_per_worker = int(os.environ["THREADS_PER_WORKER"])
        else:
            threads_per_worker = 4

        n_workers = n_cpus // threads_per_worker

        with dask.config.set(cedar_configs):
            cluster = LocalCluster(n_workers=n_workers,
                                   threads_per_worker=threads_per_worker,
                                   #processes=False,
                                   host="localhost", dashboard_address=":8787",
                                   # worker kwargs
                                   local_directory=tmpdir)

            client = Client(cluster, timeout=timeout)

    else:
        cluster = LocalCluster()
        client = Client(cluster, timeout=timeout)

    # print the cluster info
    print('\n', client, '\n')
    # prevent buffered output
    sys.stdout.flush()

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
        try:
            param_dict = json.loads(params)
        except json.JSONDecodeError:
            print('Parameter dictionary incorrectly formatted, resulting in JSONDecodeError')
            print(f'\n {params} \n')
            raise
    else:
        # set to empty dict, loop below will not itterate
        param_dict = dict()

    # start the dask distributed cluster
    client = start_cluster()

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
