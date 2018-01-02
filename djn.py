from __future__ import print_function, division
import matplotlib.pyplot as plt
import numpy as np
from plotly import offline as py
import plotly.tools as tls

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
    mid = np.arange(edgestart, edgeend + diffs, diffs)

    return mid

def show():
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

def nextcolor():
    """
    Get the next color in the default matplotlib color order
    """
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

    # compute residuals for non-parametric prediction intervals
    if np.shape(x)[0] > 6:
        pred = slope * x + intercept

        yhat = np.sort(y - pred)

        a=1-.6826
        eL = (np.shape(x)[0] + 1) * a / 2
        eU = (np.shape(x)[0] + 1) * (1 - a/2)

        L1 = np.floor(eL).astype(int) - 1
        L2 = L1 + 1

        U1 = np.ceil(eU).astype(int) - 1
        U2 = U1 - 1

        rL = yhat[L1] + (eL - L1) * (yhat[L2] - yhat[L1])
        rU = yhat[U1] - (eU - U1) * (yhat[U2] - yhat[U1])
    else:
        rL = np.nan
        rU = np.nan

    print(rL, rU)

    return slope, intercept
