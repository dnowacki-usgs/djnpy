# djnpy
Handy python utilities

## airy.py
`import airy`

Functions for linear wave theory, including celerity, energy density, and group velocity

## matlabtools.py
`import matlabtools`

Functions to interact with `.mat` files created by MATLAB. Includes `loadmat()` and `savemat()`

## noaa.py
`import noaa`

Read data (generally water level) from NOAA COOPS/Tides and Currents sites as a dict: `get_coops_data()`

## nwis.py
`import nwis`

Obtain USGS NWIS data (like water level, discharge, meteorological data, etc) via JSON and return it as a Pandas DataFrame: `nwis_json()`

## wavenumber.py
`import wavenumber`

A couple of functions to compute wavenumber from linear wave theory
