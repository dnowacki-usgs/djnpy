from __future__ import division, print_function

import datetime
import json

import numpy as np
import pandas as pd
import requests


def nwis_json(
    site="01646500",
    parm="00065",
    start=None,
    end=None,
    period=None,
    freq="iv",
    xarray=False,
):
    """Obtain NWIS data via JSON.

    Parameters
    ----------
    site : str, optional
        NWIS site code. Default '01646500'
    parm : str, optional
        Parameter code. Default '00065'
        e.g. '00065' for water level, '00060' for discharge.
    freq : str, optional
        Data frequency. Valid values are:
        - 'iv' for unit-value data (default)
        - 'dv' for daily average data
    start : str, optional
    end : str, optional
        Start and end dates in ISO-8601 format (e.g. '2016-07-01').
        If specifying start and end, do not specify period.
    period : str, optional
        Duration following ISO-8601 duration format, e.g. 'P1D' for one day
        Default one day.
        If specifying period, do not specify start and end.
    xarray : bool, optional
        If True, return an xarray Dataset instead of pandas DataFrame

    Returns
    -------
    pandas.DataFrame
        A Pandas DataFrame of the returned data with the following columns:
        - 'time': time in UTC
        - 'sitename': site name
        - 'sitecode': site code
        - 'val': value
        - 'unit': unit
        - 'variableName': variable name
        - 'timelocal': local time
        - 'latitude': site latitude
        - 'longitude': site longitude
        - 'srs': latitude and longitude spatial reference system

        If xarray is True, return an xarray Dataset instead

    Warnings
    --------
    This will fail on Windows systems when trying to retrieve data from before
    1970. See https://bugs.python.org/issue36759 for more details.

    Notes
    -----
    More info about the URL format is at https://waterservices.usgs.gov

    dnowacki@usgs.gov 2016-07
    """

    if start is not None and (pd.Timestamp(start) < pd.Timestamp("1970-01-01")):
        import warnings

        warnings.warn(
            "Requesting data from before 1970 on Windows systems will fail.",
            RuntimeWarning,
        )

    url = (
        "http://waterservices.usgs.gov/nwis/"
        + freq
        + "/?format=json&sites="
        + str(site)
        + "&parameterCd="
        + str(parm)
    )

    if period is None and start is None and end is None:
        period = "P1D"

    if period is None:
        url = url + "&startDT=" + start + "&endDt=" + end
    else:
        url = url + "&period=" + period

    try:
        payload = requests.get(url).json()
    except json.JSONDecodeError:
        raise ValueError(
            f"Error decoding JSON. For more details check the following URL in a browser: <{url}>"
        )

    try:
        pvt = payload["value"]["timeSeries"][0]
    except IndexError:
        raise IndexError(
            f"Data may not be available during the specified time period. For more details check the following URL in a browser: <{url}>"
        )

    dss = []
    for n in range(len(pvt["values"])):
        dss.append(get_single_value(pvt, n, freq, xarray))

    # if only a single value, return it not as a list
    # otherwise return as a list
    if len(dss) == 1:
        dss = dss[0]

    return dss


def get_single_value(pvt, n, freq, xarray):
    v = pvt["values"][n]["value"]
    nwis = {}
    nwis["timelocal"] = pd.to_datetime([v[i]["dateTime"] for i in range(len(v))])
    # Convert local time to UTC if unit values
    if freq == "iv":
        nwis["time"] = np.array(
            [x.astimezone(datetime.timezone.utc) for x in nwis["timelocal"]]
        )
    elif freq == "dv":
        nwis["time"] = nwis["timelocal"]  # keep naive date
    nwis["sitename"] = pvt["sourceInfo"]["siteName"]
    nwis["sitecode"] = pvt["sourceInfo"]["siteCode"][0]["value"]
    nwis["latitude"] = pvt["sourceInfo"]["geoLocation"]["geogLocation"]["latitude"]
    nwis["longitude"] = pvt["sourceInfo"]["geoLocation"]["geogLocation"]["longitude"]
    nwis["srs"] = pvt["sourceInfo"]["geoLocation"]["geogLocation"]["srs"]
    nwis["unit"] = pvt["variable"]["unit"]["unitCode"]
    nwis["variableName"] = pvt["variable"]["variableName"]
    nwis["methodDescription"] = pvt["values"][n]["method"][0]["methodDescription"]
    nwis["val"] = np.array([float(v[i]["value"]) for i in range(len(v))])
    nwis["val"][nwis["val"] == pvt["variable"]["noDataValue"]] = np.nan

    # There can be times when certain values for a variable are empty, but others have data.
    # For example, this is sometimes the case with turbidity.
    if not len(nwis["val"]):
        return None

    df = pd.DataFrame(
        nwis,
        columns=[
            "time",
            "sitename",
            "sitecode",
            "val",
            "unit",
            "variableName",
            "methodDescription",
            "timelocal",
            "latitude",
            "longitude",
            "srs",
        ],
    ).set_index("time")

    if xarray:
        ds = df.to_xarray()
        ds["time"] = pd.DatetimeIndex(ds["time"].values)
        ds["val"].attrs["units"] = ds["unit"].values[0]
        ds["val"].attrs["variableName"] = ds["variableName"].values[0]
        ds["val"].attrs["methodDescription"] = ds["methodDescription"].values[0]
        ds = ds.drop_vars(["unit", "variableName", "methodDescription", "timelocal"])
        for k in ["sitename", "latitude", "longitude", "sitecode", "srs"]:
            ds.attrs[k] = ds[k].values[0]
            ds = ds.drop_vars(k)
        return ds
    else:
        return df
