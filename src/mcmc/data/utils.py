#!/usr/bin/env python3

import os
import glob
import json
import numpy as np
import pandas as pd
import xarray as xr

def load_KMR_DETIM():
    """Load KMR (D)istributed (E)nhanced (T)emperature (I)ndex (M)odel results
    """
    # file paths to the NetCDF files
    nc_fp = "/Users/andrewnolan/Thesis/thermal-structure/input_data/mass_balance/Kaskawulsh_NetBalance.nc"

    with xr.open_dataset(nc_fp) as MB_obs:
        # stack along elevation
        stacked = MB_obs.stack(z=('x','y'))

        # Sort the indexes by elevation
        idxs  = stacked.dropna('z').z.values
        elev  = stacked.Z.dropna('z').values
        idxs  = idxs[np.argsort(elev)]

        # Calculate the Kaskawulsh ELA
        ELA_idxs = np.argpartition(np.abs(stacked.B.values), 5)
        z_ELA    = float(stacked.isel(z=ELA_idxs[:5]).Z.mean())

        # surface elevation of observed values
        z_obs = stacked.Z.sel(z=idxs).values
        # Observed Refreezing   [m  i.e. / yr]
        R_obs = stacked.R.sel(z=idxs).values
        # Observed Accumulation [ m i.e. / yr ]
        A_obs = stacked.A.sel(z=idxs).values
        # Observed melt         [ m i.e. / yr ]
        M_obs = stacked.M.sel(z=idxs).values
        # Observed mass balance [m i.e. / yr ]
        B_obs = stacked.B.sel(z=idxs).values

    return z_obs, R_obs, A_obs, M_obs, B_obs, z_ELA

def load_NetBalance_SpinUps(SpinUp_fp, keys):
    """Load the net-balance curve results from uncoupled Elmer/Ice model runs.

    Data is returned as a pandas Dataframe, which is formatted for hierechecal
    inversion procedure.
    """

    if not os.path.exists(SpinUp_fp):
        raise OSError('Invalid path to params folder')

    json_fp = os.path.join(SpinUp_fp,'params','{}.json')
    nc_fp  = os.path.join(SpinUp_fp, 'result', '{}', 'nc')

    # list to store pandas dataframes in
    dfs = []

    for key in keys:
        param_fp = json_fp.format(key)

        if not os.path.exists(SpinUp_fp):
            raise OSError(f'{key} is invalid')

        # load the parameters json dictinary
        with open(param_fp) as f:
            params = json.load(f)


        nc_folder = nc_fp.format(key)
        NetCDFs_fps = sorted(os.listdir(nc_folder))

        # get rest of the params from the dictinary for filtering
        sub_strings = [str(params[key]) for key in ['dx', 'dt', 'TT', 'fit', 'k']]
        # mass balance offset to KMR MB curve such that argmin(ΔV'(Δb))
        sub_strings.append(str(params["MB"]["value"]))

        # cull files based on params from json file, make sure to loop over copy
        # of list for this to work properly
        for file in NetCDFs_fps[:]:
            # cull files based on .json params
            if not all(sub_str in file for sub_str in sub_strings):
                NetCDFs_fps.remove(file)

        # if first attempt at culling didn't work, actually parse the offset from
        # the fp and compare offsets
        if len(NetCDFs_fps) != 1:
            for file in NetCDFs_fps[:]:
                offset = float(file.split('MB_')[-1].split('_OFF')[0])

                if params["MB"]["value"] != offset:
                    NetCDFs_fps.remove(file)

        # at this point rasise an error if things still didn't work
        if len(NetCDFs_fps) != 1:
            raise OSError('Unable to find unqiue NetCDF file')

        # join the realtive path and the found file
        src_fp = os.path.join(nc_folder, NetCDFs_fps[0])

        # read, correct, and extract from the NetCDF file
        with xr.open_dataset(src_fp) as src:
            # correct for minimum ice thickness
            src["height"] = xr.where(src.height <= 10, 0, src.height)
            # apply sigma coordinate transform for vertical coordinate
            src["z_s"]    = src.zbed + src.Z * src.height
            # Exract data of interest from NetCDF file
            z  = src['zs'].isel(t=-1,coord_2=-1).values
            MB = src['zs accumulation flux 2'].isel(t=-1,coord_2=-1).values

        # list of dataframes to be concatenated
        dfs.append(pd.DataFrame({f'z':z, "MB": MB, 'key':[key]*len(z)}))

    # concatenate and factorize the dataframes
    Elmer_Runs = pd.concat(dfs, ignore_index=True)
    Elmer_Runs["key_factor"] =  pd.factorize(Elmer_Runs.key)[0]

    return Elmer_Runs
