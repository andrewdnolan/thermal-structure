#!/usr/bin/env python3
import numpy as np
import xarray as xr
from scipy import optimize
from scipy import interpolate
import matplotlib.pyplot as plt


def poly_fit(ds, deg=1, ret_p=True):
    # Deal with nan's internally
    x = ds.Elevation.values
    y = ds.MB.values

    # find indexes of nans and sort indexes
    mask = np.where(~np.isnan(x))
    idxs = np.argsort(x[mask])

    # Remove the nans and sort the data
    x_ = x[  mask][idxs]
    y_ = y[mask  ][idxs]

    p   = np.polyfit(x_, y_, deg)
    fit = np.polyval(p, x_)

    if ret_p:
        return x_, fit, p
    else:
        return x_, fit

def spline_fit(ds, deg=3, ret_p=True, s=240):
    # Deal with nan's internally
    x = ds.Elevation.values
    y = ds.MB.values

    # find indexes of nans and sort indexes
    mask = np.where(~np.isnan(x))
    idxs = np.argsort(x[mask])

    # Remove the nans and sort the data and downsample data so there are
    # less verically overlappig points, causing the spline fitting to fail
    x_sub = x[~np.isnan(x)][::10]
    y_sub = y[~np.isnan(x)][::10]

    m = len(x_sub)

    idxs = np.argsort(x_sub, kind='heapsort')

    x_ = x_sub[idxs]
    y_ = y_sub[idxs]


    tck  = interpolate.splrep(x_, y_, k=deg, s=s)
    xnew = np.linspace(x_.min(), x_.max(), 1000)
    fit  = interpolate.splev(xnew, tck)

    if ret_p:
        return xnew, fit, tck
    else:
        return xnew, fit


def pert_ela(z, tck, δb_bounds=(0.0, 0.0), n_δb = 250):
    if type(tck) == tuple:
        b = interpolate.splev(z, tck)[:,np.newaxis]
    else:
        b = np.polyval(tck, z)[:,np.newaxis]

    δb    = np.linspace(δb_bounds[0], δb_bounds[1], n_δb)[np.newaxis, :]
    Δb    = b + δb
    ELAs  = z[np.argpartition(np.abs(Δb),2, axis=0)[:2,:]].mean(axis=0)

    return δb, ELAs, b

# ------------------------------------------------------------------------------
# Archived functions
# ------------------------------------------------------------------------------
def latex_float(f):
    """https://stackoverflow.com/a/13490601/10221482"""
    float_str = "{0:.2e}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

def piecewise_linear(x, x0, y0, k1, k2):
    """https://stackoverflow.com/a/29384899/10221482"""
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

def piecewise_fit(ds, deg=1, ret_p=True):
    # Deal with nan's internally
    x = ds.Elevation.values
    y = ds.Accumulation.values

    # find indexes of nans and sort indexes
    mask = np.where(~np.isnan(x))
    idxs = np.argsort(x[mask])

    # Remove the nans and sort the data
    x_ = x[  mask][idxs]
    y_ = y[mask  ][idxs]

    # Fit the model, optimize the slopes and y-intercepts
    # (and therefore the knot point)
    p , e = optimize.curve_fit(piecewise_linear, x_, y_, p0=[ 2300,0,0,0])

    if ret_p:
        return x_, piecewise_linear(x_, *p), p
    else:
        return x_, piecewise_linear(x_, *p)
