# djnpy
Handy python utilities

## airy.py
`import airy`

Functions for linear wave theory, including celerity, energy density, and group velocity

## djn.py
`import djn`

Useful little tidbits:

* `uv2sd()` convert east, north components to speed, direction
* `sd2uv()` convert speed, direction, to east, north components
* `boxoff()` a Matlab-like `box off` function for matplotlib
* `twinboxoff()` a Matlab-like `box off` function for matplotlib with twin y axes
* `thinspines()` make thinner matplotlib plot spines (default linewidth 0.5)
* `find_nearest()` find nearest value to the requested value in a numpy array
* `middles()` given a vector of bin edges of length n, make a vector of the middle values with length n-1.
* `show()` like `plt.show()` but for plotly. Makes it easy to turn a matplotlib plot into a plotly plot.
* `set_fontsize()` set fontsize of all elements of a figure, like `djn.set_fontsize(plt.gcf(), 14)`
* `getcols()` return default color order for matplotlib plots
* `nextcolor()` return the next color in the default matplotlib color order
* `splabel()` plot a subplot letter label like "(a)" in figures
* `siegel()` compute robust regression using repeated medians, following Siegel (1982). From Vlad Niculae.
* `princax()` compute principal axis, rotation axis, and principal ellipse for velocity data. From Rich Signell.
* `rot_earth()` rotate vectors u, v by given number of degrees using earthwise coordinates (0=N, 90=E)
* `tidalfilt()` low-pass filter data using a 5th order Butterworth filter. Useful for tidally filtering data.
* `get_nan_block_idxs()` return start and stop indexes of (by default) nan blocks`
* `xcorr()` cross-correlate data
* `haversine()` calculate the great circle distance in kilometers between two points on the earth (specified in decimal degrees)
* `argmaxn()` return largest N values from array

## matlabtools.py
`import matlabtools`

Functions to interact with `.mat` files created by MATLAB. Some functions based on code originally written by J.Paul Rinehimer.
* `loadmat()` load `.mat` file
* `savemat()` save `.mat` file
* `matlab2datetime()` convert MATLAB serial time to Python datetime, optionally with timezone
* `datetime2matlab()` convert Python datetime to MATLAB serial time

## noaa.py
`import noaa`

* `get_coops_data()` Read data (water level, air pressure, wind, etc) from NOAA COOPS/Tides and Currents sites via the API as a dict
* `get_long_coops_data()` Read data longer than one month and return as an xarray Dataset

## nwis.py
`import nwis`

* `nwis_json()` Obtain USGS NWIS data (like water level, discharge, meteorological data, etc) via JSON and return it as a Pandas DataFrame

## wavenumber.py
`import wavenumber`

A couple of functions to compute wavenumber from linear wave theory
