#!/usr/bin/env python3

import xarray as xr

from elements import calc_element_area, calc_element_mean

def calc_relative_volume(src):

    rel_vol = src.height.isel(coord_2 = -1).integrate('coord_1') /\
              src.height.isel(coord_2 = -1, t = 0).integrate('coord_1')

    return rel_vol

def calc_percent_temperate(src, dz_var='height'):
    """
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
