#!/usr/bin/env python3

"""
Functions for working with element wise data, designed for a strcuted
quadrilateral mesh
"""

# dictionaries for vertical indexing structured mesh
l_i   = dict(coord_2=slice(0, -1),  coord_1=slice(1,None))
l_ip1 = dict(coord_2=slice(1,None), coord_1=slice(1,None))
r_i   = dict(coord_2=slice(0, -1),  coord_1=slice(0, -1))
r_ip1 = dict(coord_2=slice(1,None), coord_1=slice(0, -1))

def calc_element_mean(src, var):
    """ Calculate element wise mean of variable

    Inputs:
        src (xr.Dataset) --> Dataset (or DataArray) to access field from
        var        (str) --> variable to calculate elment mean of

    Outputs:
        (xr.DataArray)   --> element mean stacked along the "element" coordinate
    """
    if var in src.variables:
        field = src.get(var)
    else:
        raise KeyError(f'{var} not found in src')

    # sum of field on left  side of element
    l_s = field.isel(l_ip1) + field.isel(l_i)
    # sum of field on right side of element
    r_s = field.isel(r_ip1) + field.isel(r_i)

    # return element averaged field
    return ((l_s + r_s) / 4).stack(element=('coord_2', 'coord_1'))


def calc_element_area(src, var='height'):
    """Calculate element area of a structured mesh

    Note: We use the approximate area of the quadrilateral by taking the avergae
          height of the left and right sides of the element, then using the width
          approximate the area as though the quadrilateral were a perfect rectangle.
          A more general formula could be used, but results in insignificant
          differences (<1e-6 [m2])

    Inputs:
        src (xr.Dataset) --> Dataset (or DataArray) to access field from
        var        (str) --> variable used to calcute vetical gridcell spacing

    Outputs:
        (xr.DataArray)   --> element area stacked along the "element" coordinate
    """
    if var in src.variables:
        dz_field = src.get(var)
    else:
        raise KeyError(f'{var} not found in src')

    # vertical spacing (dz) on left side of element
    l_dz = dz_field.isel(l_ip1) - dz_field.isel(l_i)
    # vertical spacing (dz) on right side of element
    r_dz = dz_field.isel(r_ip1) - dz_field.isel(r_i)
    # horizontal spacing (dx), average is used since horizontal grid is never updated
    m_dx = (src.X.isel(l_ip1) - src.X.isel(r_ip1)).mean()

    # return the element area
    return (m_dx * (l_dz + r_dz) / 2).stack(element=('coord_2', 'coord_1'))
