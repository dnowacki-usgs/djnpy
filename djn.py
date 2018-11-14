 # -*- coding: utf-8 -*-
from __future__ import print_function, division
import matplotlib.pyplot as plt
import numpy as np
from plotly import offline as py
import plotly.tools as tls
import scipy

def uv2sd(u, v):
    """Convert east, north components to speed and direction"""

    u = np.asarray(u)
    v = np.asarray(v)

    spds = np.full_like(u, np.nan)
    # dirs = np.full_like(u.values, np.nan)

    # print(spds)
    # print(dirs)

    spds = np.sqrt(u**2 + v**2)
    dirs = np.rad2deg(np.arctan(u/v))
    dirs[np.logical_and(u == 0, v > 0)] = 0
    # dirs[np.logical_and(u > 0, v > 0)] = np.rad2deg(np.arctan(u[np.logical_and(u > 0, v > 0)] / v[np.logical_and(u > 0, v > 0)]))
    dirs[np.logical_and(u > 0, v == 0)] = 90
    dirs[np.logical_and(u > 0, v < 0)] = dirs[np.logical_and(u > 0, v < 0)] + 180
    dirs[np.logical_and(u == 0 , v < 0)] = 180
    dirs[np.logical_and(u < 0, v < 0)] = dirs[np.logical_and(u < 0, v < 0)] + 180
    dirs[np.logical_and(u < 0, v == 0)] = 270
    dirs[np.logical_and(u < 0, v > 0)] = dirs[np.logical_and(u < 0, v > 0)] + 360

    return spds, dirs

def sd2uv(s, d):
    """
    Convert speed and direction to u, v components
    """
    s = np.asarray(s)
    d = np.asarray(d)

    u = s * np.sin(d * np.pi / 180)
    v = s * np.cos(d * np.pi / 180)

    return u, v

def boxoff():
    """
    A Matlab-like boxoff to remove top & right border of plots
    """
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

def twinboxoff():
    """
    A Matlab-like boxoff for twinx plots
    """
    plt.gca().spines['top'].set_visible(False)

def thinspines(lw=0.5):
    for axis in ['top','bottom','left','right']:
        plt.gca().spines[axis].set_linewidth(lw)
        plt.gca().tick_params(width=lw)

def find_nearest(array,value):
    """
    Find nearest value in numpy array
    http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    idx = (np.abs(array-value)).argmin()
    return idx

def middles(edges):
    """
    Make middles vector from edges vector
    """
    diffs = np.median(np.diff(edges));
    edgestart = edges[0] + diffs/2.;
    edgeend = edges[-1] - diffs/2.;

    return np.linspace(edgestart, edgeend, len(edges)-1)

def show():
    """
    Easy way to make an mpl plot into a plotly plot.
    Call djn.show() instead of plt.show()
    """
    py.iplot(tls.mpl_to_plotly(plt.gcf()))

def set_fontsize(fig,fontsize):
    """
    For each text object of a figure fig, set the font size to fontsize
    http://stackoverflow.com/questions/7082597/in-matplotlib-how-do-you-change-the-fontsize-of-a-single-figure
    """
    def match(artist):
        return artist.__module__ == "matplotlib.text"

    for textobj in fig.findobj(match=match):
        textobj.set_fontsize(fontsize)

def getcols():
    """
    Get the default color order
    """
    return plt.rcParams['axes.prop_cycle'].by_key()['color']

def nextcolor(n=1):
    """
    Get the next color in the default matplotlib color order
    """
    for x in range(n):
        next(plt.gca()._get_lines.prop_cycler)['color']

def siegel(x, y):
    """
    Compute robust regression using repeated medians, following Siegel (1982)

    Inputs:
    x, y: x and y locations of points

    Outputs:
    slope, intercept: slope and intercept of robust regression line

    ANDREW F. SIEGEL; Robust regression using repeated medians, Biometrika,
    Volume 69, Issue 1, 1 April 1982, Pages 242–244,
    https://doi.org/10.1093/biomet/69.1.242

    Based on 2-clause BSD licensed code by Vlad Niculae, available at
    http://codegists.com/snippet/python/siegelpy_vene_python
    """

    x = np.asarray(x)
    y = np.asarray(y)
    deltax = x[:, np.newaxis] - x
    deltay = y[:, np.newaxis] - y

    slopes = deltay / deltax

    slope = np.median(np.nanmedian(slopes, axis=0))
    intercept = np.median(y - slope * x)

    # # compute residuals for non-parametric prediction intervals
    # if np.shape(x)[0] > 6:
    #     pred = slope * x + intercept
    #
    #     yhat = np.sort(y - pred)
    #
    #     a=1-.6826
    #     eL = (np.shape(x)[0] + 1) * a / 2
    #     eU = (np.shape(x)[0] + 1) * (1 - a/2)
    #
    #     L1 = np.floor(eL).astype(int) - 1
    #     L2 = L1 + 1
    #
    #     U1 = np.ceil(eU).astype(int) - 1
    #     U2 = U1 - 1
    #
    #     rL = yhat[L1] + (eL - L1) * (yhat[L2] - yhat[L1])
    #     rU = yhat[U1] - (eU - U1) * (yhat[U2] - yhat[U1])
    # else:
    #     rL = np.nan
    #     rU = np.nan
    #
    # print(rL, rU)

    return slope, intercept

def princax(w):
    """
     PRINCAX Principal axis, rotation angle, principal ellipse

       [theta,maj,min,wr]=princax(w)

       Input:  w   = complex vector time series (u+i*v)

       Output: theta = angle of maximum variance, math notation (east == 0, north=90)
               maj   = major axis of principal ellipse
               min   = minor axis of principal ellipse
               wr    = rotated time series, where real(wr) is aligned with
                       the major axis.

     For derivation, see Emery and Thompson, "Data Analysis Methods
       in Oceanography", 1998, Pergamon, pages 325-327.  ISBN 0 08 0314341
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Version 1.0 (12/4/1996) Rich Signell (rsignell@usgs.gov)
     Version 1.1 (4/21/1999) Rich Signell (rsignell@usgs.gov)
         fixed bug that sometimes caused the imaginary part
         of the rotated time series to be aligned with major axis.
         Also simplified the code.
     Version 1.2 (3/1/2000) Rich Signell (rsignell@usgs.gov)
         Simplified maj and min axis computations and added reference
         to Emery and Thompson book
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    """

    ind = np.nonzero(np.isfinite(w))
    wr = w.copy()
    w = w[ind].copy()

    # find covariance matrix
    cv = np.cov([np.real(w), np.imag(w)])

    # find direction of maximum variance
    theta = 0.5 * np.arctan2(2. * cv[1, 0], cv[0, 0] - cv[1, 1])

    # find major and minor axis amplitudes
    term1 = cv[0, 0] + cv[1, 1]
    term2 = np.sqrt((cv[0, 0] - cv[1, 1])**2 + 4. * cv[1, 0]**2)
    majo = np.sqrt(0.5 * (term1 + term2))
    mini = np.sqrt(0.5 * (term1 - term2))

    # rotate into principal ellipse orientation
    wr[ind] = w * np.exp(-1j * theta)
    theta = theta * 180 / np.pi

    # convert from math notation to geographic coordinates
    degrees = (450 - theta) % 360

    return theta, majo, mini, wr, degrees

def rot_earth(u, v, degrees):
    """
    Rotate vectors u, v by given number of degrees using earthwise coordinates
    (0=N, 90=E)
    - Positive degrees results in a counterclockwise (CCW) rotation
    - Negative degrees rotates values clockwise
    """

    up = np.cos(np.deg2rad(degrees)) * u - np.sin(np.deg2rad(degrees)) * v
    vp = np.sin(np.deg2rad(degrees)) * u + np.cos(np.deg2rad(degrees)) * v
    return up, vp

def tidalfilt(inmat, fs, cutoff=48.):
    """
    Low-pass filter data using a 5th order Butterworth filter
    """

    # fs in samples per hour
    b, a = scipy.signal.butter(5, (1./cutoff)/(fs/2.))
    return scipy.signal.filtfilt(b, a, inmat)

def get_nan_block_idxs(a, f=np.isnan, mindiff=0, maxdiff=None):
    """
    Modified from https://stackoverflow.com/a/15200385/3657988
    Returns start and stop indexes of (by default) nan blocks
    Can also return non-nan blocks by modifying f (e.g. by using np.isfinite)
    minlength is minimum size of the block of nans (or non-nans)
    """
    nan_mask = f(a)
    start_nans_mask = np.concatenate((np.resize(nan_mask[...,0],a.shape[:-1]+(1,)),
                                 np.logical_and(np.logical_not(nan_mask[...,:-1]), nan_mask[...,1:])
                                 ), axis=a.ndim-1)
    stop_nans_mask = np.concatenate((np.logical_and(nan_mask[...,:-1], np.logical_not(nan_mask[...,1:])),
                                np.resize(nan_mask[...,-1], a.shape[:-1]+(1,))
                                ), axis=a.ndim-1)

    start_idxs = np.where(start_nans_mask)[0]
    stop_idxs = np.where(stop_nans_mask)[0]
    # return stop_idxs[-1] - start_idxs[-1] + 1
    # return start_idxs, stop_idxs
    idxs = np.vstack([start_idxs, stop_idxs]).T

    d = np.ravel(np.diff(idxs))

    if maxdiff:
        return idxs[(d >= mindiff) & (d <= maxdiff), :]
    else:
        return idxs[d >= mindiff, :]

def trim_max_diff(da, diff):
    da[np.ediff1d(da, to_begin=0) > diff] = np.nan

def trim_min_diff(da, diff):
    da[np.ediff1d(da, to_begin=0) < diff] = np.nan
