#!/usr/bin/env python3

import warnings
import numpy as np
import xarray as xr
from scipy import integrate
from .elements import calc_element_area, calc_element_mean

def check_filtered(func):
    """ Decorator to check if the ice thickness has been filtered.
        If not, a warning is raised.
    """
    def inner(src, **kwargs):
        # check if "h" greater than zero in all gricells, suggesting "h" was not
        # filtered correctly, OR the glacier has reached the end of the domain
        # NOTE: bottom most gridcell (coord_2=0) is excluded, since h always = 0
        if (src.height.isel(coord_2=slice(1,None)) < 0).all():
            warnings.warn('Make sure fictitious ice thickness has been removed.')
        return func(src,**kwargs)
    return inner

def calc_magnitude(a, b):
    """ Calculate vector magnitude.

        Using xt.apply_ufunc to apply a vectorized function for unlabeled arrays
        on xarray objects.
    """
    func = lambda x, y: np.sqrt(x**2 + y**2)

    return xr.apply_ufunc(func, a, b, dask="allowed")

@check_filtered
def calc_volume(ds):
    """ Compute the glacier volume per unit width [m^2], by trapezoid integration

        This function work on both xarray and dask objects. Input ds would be a
        dask object, if it was chunked when it was read in
    """
    return xr.apply_ufunc(integrate.trapezoid,
                          ds.isel(coord_2=-1).height,
                          ds.isel(coord_2=-1).X,
                          input_core_dims=[["coord_1"], ["coord_1"]],
                          kwargs={"axis": -1},
                          dask="parallelized")

@check_filtered
def calc_percent_temperate(src, dz_var='Z'):
    """ Calculate the percentage temperate using element areas

    Note: we now use the water content to find temperate ice
          instead of "phase change enthalpy" because enthalpy was 
          within floating point roundoff of "phase change enthalpy" 
          and caused spurious oscillations in the percent temperate
    """
    # Calculate the element area using structured grid
    elm_area  = calc_element_area(src, var=dz_var)
    #################################################################
    # Calculate mean elemental enthalpy 
    elm_enth = calc_element_mean(src, 'enthalpy_h')
    # Calculate mean elemental phase change enthalpy (pce)
    elm_PCE  = calc_element_mean(src, 'phase change enthalpy')
    # Old way of calculating temperate mask: 
    mask = elm_enth >= elm_PCE
    ################################################################

    # # Calculate mean elemental water content 
    # elm_omega = calc_element_mean(src, 'water content')

    # # new way of calculating temperate mask: 
    # # i.e. must have some water content 
    # mask = elm_omega >= 1e-4

    # find total glacier area [m2]
    A_totl = elm_area.sum('element')
    # find total temperate area [m2]
    A_temp = xr.where(mask, elm_area, 0.0).sum('element')

    # convert to percentage [%]
    return (A_temp / A_totl) * 100

@check_filtered
def calc_mean_enthalpy(src, dz_var='Z'):
    """ Calculate the weighted (by element area) mean enthalpy [J kg-1]
    """

    # Calculate the elemental area [m2]
    elm_area = calc_element_area(src, var=dz_var)
    # Calculate mean elemental enthalpy [J kg-1]
    elm_enth = calc_element_mean(src, 'enthalpy_h')
    # Calculate the total glacier area [m2]
    tot_area = elm_area.sum('element')

    # calculate the weighted mean [J kg-1]
    enth_bar = (elm_enth * elm_area/tot_area).sum('element')

    # weighting will return 0 if results are "nan", so check and fix
    return xr.where(enth_bar==0, np.nan, enth_bar)

def calc_length(src, H_min=10.):
    """Calculate the glacier length [km]

    Inputs:
        src (xr.Dataset) --> Dataset with 'height' variables. Supports dask
                             delayed objects

        H_min    (float) --> Minimum ice thickness used to filter passive nodes

    Outputs:
        (xr.DataArray)   --> Length [km] with 't' dimension / coordinate
    """
    # function to find terminus position as a function of time
    find_Term = lambda x: x.where(~x.isnull(), drop=True).max('coord_1')

    # (T)otal (D)omain (L)ength [km]
    TDL = float(src.X.max()) / 1e3
    # (N)umber of (H)orizontal (N)odes
    NHN = int(src.coord_1.max())

    # ice thickness [m] along free surface
    H = src.height.isel(coord_2=-1)

    # check if all values are NAN, then return NAN
    if H.isnull().all():
        # copy coordinates of input dataset, but squeeze coordinates since 
        # dataarray should have been reduced along them (i.e. coord_1, coord_2)
        Length = xr.DataArray(np.nan, coords=H.coords).squeeze()
        
    #otherwise, calculate length as expected
    else: 
        # Mask to find passive nodes
        passive = xr.where(H <= (H_min+1.0), src.coord_1, np.nan)
        # Get indexes of the ice free nodes
        ice_free = passive.diff('coord_1')*passive.coord_1.isel(coord_1=slice(0,-1))
        # Get terminus index
        term_idx = (NHN - find_Term(ice_free)).astype(int)
        # Get the glacier length as a function of time [km]
        Length = src.X.isel(coord_1=term_idx, coord_2=-1)/1e3

    return Length

def _get_ELA_indexes(src): 
    # function to find the indexes of the ELA, as a function of time 
    # Generalized enough to support multiple horizontal indexes being at the ELA
        
    # create a two dimensional 'coord_1' array for boolean masking
    _, coord_1 = xr.broadcast(src.t, src.coord_1)
    
    # the difference of the sign will tell us where the mass balance crosses zero
    mask = xr.apply_ufunc(np.diff, 
                          np.sign(src['mass balance']), 
                          kwargs={"axis": -1, "prepend" : -1},
                          dask="allowed")

    # the difference should be 2 at the ELA, 0 all other places
    ELA_idxs = coord_1.where(mask == 2, drop=False)
    
    # returns nan values everywhere except the horizontal coordinates that 
    # coincide w/ the ELA
    return ELA_idxs

def Variable_at_ELA(src, variable:str): 
    """Extract the specified `variable`'s value along at the ELA as function of time

    Inputs:
        src (xr.Dataset) --> Dataset with 'height' variables. Supports dask
                             delayed objects

        variable  (str) --> Valid variable within `src` dataset to extract

    Outputs:
        (xr.DataArray)   --> 
    """
     
    # get the ELA mask
    ELA_mask = _get_ELA_indexes(src)
    # get the desired variable at the non-nan horizontal coordinates (i.e. ELA)
    # in the case where there are multiple horizontal nodes that correspond to 
    # the ELA, take the average. Otherwise will return the value associated with
    # the single coordinate. 
    var_at_ELA = src[variable].where(~ELA_mask.isnull(),drop=True).mean('coord_1')

    return var_at_ELA

def calc_Peclet(src, dim="X"):
    np.diff(np.sign(src['mass balance'].isel(t=0, coord_2=-1)), prepend=-1) == 2
