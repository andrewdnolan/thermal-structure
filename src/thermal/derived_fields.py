#!/usr/bin/env python3

import numpy as np
import xarray as xr
from scipy import integrate
from elements import calc_element_area, calc_element_mean

def calc_magnitude(a, b):
    """ Calculate vector magnitude.

        Using xt.apply_ufunc to apply a vectorized function for unlabeled arrays
        on xarray objects.
    """
    func = lambda x, y: np.sqrt(x**2 + y**2)

    return xr.apply_ufunc(func, a, b, dask="allowed")

def calc_volume(ds):
    """ Compute the glacier volume per unit width [m^2], by trapeziod integration

        This function work on both xrray and dask objects. Input ds would be a
        dask object, if it was chunked when it was read in
    """
    return xr.apply_ufunc(integrate.trapz,
                          ds.isel(coord_2=-1).height,
                          ds.isel(coord_2=-1).X,
                          input_core_dims=[["coord_1"], ["coord_1"]],
                          kwargs={"axis": -1},
                          dask="parallelized")

def calc_percent_temperate(src, dz_var='height'):
    """ Calculate the percentage temperate using element areas
    """
    # Calculate the element area using structured grid
    elm_area = calc_element_area(src, dz_var)
    # Calculate mean elemental enthalpy
    elm_enth = calc_element_mean(src, 'enthalpy_h')
    # Calculate mean elemental phase change enthalpy (pce)
    elm_PCE  = calc_element_mean(src, 'phase change enthalpy')

    # mask to determine temperate nodes
    mask = elm_enth >= elm_PCE

    # find total glacier area [m2]
    A_totl = elm_area.sum('element')
    # find total temperate area [m2]
    A_temp = xr.where(mask, elm_area, 0.0).sum('element')

    # convert to percentage [%]
    return (A_temp / A_totl) * 100

def calculate_length(src, H_min=10.):
    """Calculate the glacier length [km]

    Inputs:
        src (xr.Dataset) --> Dataset with 'height' varibales. Supports dask
                             delayed objects

        H_min    (float) --> Minimum ice thicknes used to filter passive nodes

    Outputs:
        (xr.DataArray)   --> Length [km] with 't' dimension / coordinate
    """
    # function to find terminus postion as a function of time
    find_Term = lambda x: x.where(~x.isnull(), drop=True).max('coord_1')

    # (T)otal (D)omain (L)ength [km]
    TDL = float(src.X.max()) / 1e3
    # (N)umber of (H)orizontal (N)odess
    NHN = int(src.coord_1.max())

    # ice thickness [m] along free surface
    H = src.height.isel(coord_2=-1)
    # Mask to find passive nodes
    passive = xr.where(H <= (H_min+1.0), src.coord_1, np.nan)

    # Get indexes of the ice free nodes
    ice_free = passive.diff('coord_1')*passive.coord_1.isel(coord_1=slice(0,-1))
    # Get terminus index
    term_idx = (NHN - find_Term(ice_free)).astype(int)
    # Get the glacier lenght as a function of time [km]
    Length = src.X.isel(coord_1=term_idx, coord_2=-1)/1e3

    return Length
