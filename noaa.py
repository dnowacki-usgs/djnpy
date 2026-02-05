import time
import urllib

import numpy as np
import pandas as pd
import pytz
import requests
import xarray as xr
from dateutil import parser
from requests.adapters import HTTPAdapter, Retry


def get_coops_data(
    station,
    start_date,
    end_date,
    product="hourly_height",
    units="metric",
    datum="MLLW",
    time_zone="GMT",
    interval=False,
    bin=False,
):
    """
    units can be 'english' or 'metric'

    start_date and end_date must be formatted like:
    yyyyMMdd, yyyyMMdd HH:mm, MM/dd/yyyy, or MM/dd/yyyy HH:mm

    product options include 'water_level', 'hourly_height', 'predictions'
    from https://tidesandcurrents.noaa.gov/api/
    Option              Description
    water_level         Preliminary or verified water levels, depending on availability.
    air_temperature     Air temperature as measured at the station.
    water_temperature   Water temperature as measured at the station.
    wind                Wind speed, direction, and gusts as measured at the station.
    air_pressure        Barometric pressure as measured at the station.
    air_gap             Air Gap (distance between a bridge and the water's surface) at the station.
    conductivity        The water's conductivity as measured at the station.
    visibility          Visibility from the station's visibility sensor. A measure of atmospheric clarity.
    humidity            Relative humidity as measured at the station.
    salinity            Salinity and specific gravity data for the station.
    hourly_height       Verified hourly height water level data for the station.
    high_low            Verified high/low water level data for the station.
    daily_mean          Verified daily mean water level data for the station.
    monthly_mean        Verified monthly mean water level data for the station.
    one_minute_water_level  One minute water level data for the station.
    predictions         6 minute predictions water level data for the station.
    datums              datums data for the stations.
    currents            Currents data for currents stations.
    """

    url = (
        "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product="
        + product
        + "&application=NOS.COOPS.TAC.WL&begin_date="
        + urllib.parse.quote(str(start_date))
        + "&end_date="
        + urllib.parse.quote(str(end_date))
        + "&datum="
        + datum
        + "&station="
        + str(station)
        + "&time_zone="
        + time_zone
        + "&units="
        + units
        + "&format=json"
    )

    if interval:
        url = url + "&interval=" + interval

    if bin:
        url = f"{url}&bin={bin}"

    if product == "currents_predictions":
        url = f"{url}&vel_type=speed_dir"

    s = requests.Session()

    retries = Retry(total=5, backoff_factor=0.1, status_forcelist=[500, 502, 503, 504])

    s.mount("https://", HTTPAdapter(max_retries=retries))

    payload = s.get(url, timeout=10).json()

    if "error" in payload.keys():
        raise ValueError("Error in returning dataset: " + payload["error"]["message"])

    t = []
    v = []
    wind = {"s": [], "d": [], "dr": [], "g": [], "f": []}
    cp = {"Speed": [], "Bin": [], "Direction": [], "Depth": []}
    cur = {"s": [], "d": [], "b": []}
    if product in [
        "water_level",
        "hourly_height",
        "air_pressure",
        "air_temperature",
        "water_temperature",
    ]:
        d = payload["data"]
    elif product == "monthly_mean":
        d = payload["data"]
        monthly = {x: [] for x in d[0].keys() if x != "month" and x != "year"}
    elif product == "predictions":
        d = payload["predictions"]
    elif product == "wind":
        d = payload["data"]
    elif product == "currents_predictions":
        d = payload["current_predictions"]["cp"]
    elif product == "currents":
        d = payload["data"]
    elif product == "datums":
        datums = {}
        for k in payload["datums"]:
            datums[k["n"]] = float(k["v"])
        return datums

    for n in range(len(d)):
        if product == "currents_predictions":
            t.append(pytz.utc.localize(parser.parse(d[n]["Time"])))
        elif product != "monthly_mean":
            t.append(pytz.utc.localize(parser.parse(d[n]["t"])))
        else:
            t.append(parser.parse(d[n]["year"] + "-" + d[n]["month"] + "-01"))
        if product == "wind":
            for k in ["s", "d", "dr", "g", "f"]:
                if k == "dr" or k == "f":
                    wind[k].append(d[n][k])
                else:
                    try:
                        wind[k].append(float(d[n][k]))
                    except ValueError:
                        wind[k].append(np.nan)
        elif product == "monthly_mean":
            for k in monthly.keys():
                monthly[k].append(float(d[n][k]))
        elif product == "currents_predictions":
            for k in cp:
                cp[k].append(float(d[n][k]))
        elif product == "currents":
            for k in cur:
                cur[k].append(float(d[n][k]))
        else:
            if d[n]["v"] == "":
                v.append(np.nan)
            else:
                v.append(float(d[n]["v"]))

    ds = xr.Dataset()

    n = {}
    n["time"] = np.array(t)
    if product == "wind":
        for k in ["s", "d", "dr", "g", "f"]:
            n[k] = np.array(wind[k])
    elif product == "monthly_mean":
        for k in monthly.keys():
            n[k] = np.array(monthly[k])
    elif product == "currents_predictions":
        for k in cp:
            n[k] = np.array(cp[k])
        if "units" in payload["current_predictions"]:
            ds.attrs["units"] = payload["current_predictions"]["units"]
    elif product == "currents":
        for k in cur:
            n[k] = np.array(cur[k])
    else:
        n["v"] = np.array(v)

    for k in n:
        ds[k] = xr.DataArray(n[k], dims="time")

    if product == "currents":
        # also get metadata for bins
        bins = requests.get(
            f"https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{station}/bins.json"
        ).json()
        ds["depth"] = bins["bins"][bin - 1]["depth"]

    if "metadata" in payload:
        for k in payload["metadata"]:
            ds.attrs[k] = payload["metadata"][k]

    products = [
        "water_level",
        "hourly_height",
        "high_low",
        "daily_mean",
        "monthly_mean",
        "one_minute_water_level",
        "predictions",
    ]
    if product in products:
        ds.attrs["datum"] = datum

    ds["time"] = pd.DatetimeIndex(ds["time"].values)
    ds["time"] = pd.DatetimeIndex(
        ds["time"].values
    )  # don't know why we need to do it twice, but we do in order for it to return as datetim64[ns]

    return ds


def get_long_coops_data(
    station,
    start_date,
    end_date,
    product="hourly_height",
    units="metric",
    datum="MLLW",
    time_zone="GMT",
    interval=False,
    bin=False,
):
    """
    Get NOAA CO-OPS data for longer than 1 month.
    This function makes recursive calls to get_coops_data() for time ranges
    longer than one month.

    units can be 'english' or 'metric'

    start_date and end_date must be formatted like:
    yyyyMMdd, yyyyMMdd HH:mm, MM/dd/yyyy, or MM/dd/yyyy HH:mm

    product options include 'water_level', 'hourly_height', 'predictions'
    from https://tidesandcurrents.noaa.gov/api/
    Option              Description
    water_level         Preliminary or verified water levels, depending on availability.
    air_temperature     Air temperature as measured at the station.
    water_temperature   Water temperature as measured at the station.
    wind                Wind speed, direction, and gusts as measured at the station.
    air_pressure        Barometric pressure as measured at the station.
    air_gap             Air Gap (distance between a bridge and the water's surface) at the station.
    conductivity        The water's conductivity as measured at the station.
    visibility          Visibility from the station's visibility sensor. A measure of atmospheric clarity.
    humidity            Relative humidity as measured at the station.
    salinity            Salinity and specific gravity data for the station.
    hourly_height       Verified hourly height water level data for the station.
    high_low            Verified high/low water level data for the station.
    daily_mean          Verified daily mean water level data for the station.
    monthly_mean        Verified monthly mean water level data for the station.
    one_minute_water_level  One minute water level data for the station.
    predictions         6 minute predictions water level data for the station.
    datums              datums data for the stations.
    currents            Currents data for currents stations.
    """

    # date ranges in 1-month chunks
    dates = pd.date_range(start_date, end_date, freq="MS")
    # need to add the beginning and end since pd.date_range normalizes to start of month
    if pd.Timestamp(start_date) < dates[0]:
        dates = dates.insert(0, pd.Timestamp(start_date))
    if pd.Timestamp(end_date) > dates[-1]:
        dates = dates.append(pd.DatetimeIndex([end_date]))
    data = []
    for n in range(len(dates) - 1):
        print(dates[n], dates[n + 1])
        try:
            data.append(
                get_coops_data(
                    station,
                    dates[n].strftime("%Y%m%d %H:%M"),
                    dates[n + 1].strftime("%Y%m%d %H:%M"),
                    product=product,
                    units=units,
                    datum=datum,
                    time_zone=time_zone,
                    interval=interval,
                    bin=bin,
                )
            )
        except ValueError as e:
            # sometimes a station goes offline
            if "No data was found" in repr(e):
                print("no data in {} - {}; skipping".format(dates[n], dates[n + 1]))
                continue
            else:
                print(e)
                continue

    ds = xr.concat(data, dim="time")
    # deduplicate times if necessary
    _, index = np.unique(ds["time"], return_index=True)
    ds = ds.isel(time=index)
    ds["time"] = pd.DatetimeIndex(ds["time"].values)

    return ds


def get_coops_metadata(station):
    """
    Get NOAA CO-OPS metadata for a single station
    """
    url = f"https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{station}.json"

    payload = requests.get(url).json()

    return payload


def get_isd(site, years):
    """
    Retrieve NOAA Integrated Surface Data (ISD) from NOAA NCEI

    Formatting from https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf
    INTEGRATED SURFACE DATASET (search online for this)
    https://www.ncei.noaa.gov/access/search/data-search/global-hourly
    """

    # ensure list
    years = [years] if isinstance(years, (str, int)) else years

    tmp = []
    for year in years:
        print(year)
        tmp.append(
            pd.read_csv(
                f"https://www.ncei.noaa.gov/data/global-hourly/access/{year}/{site}.csv",
                low_memory=False,
            )
        )
        time.sleep(2)
    df = pd.concat(tmp, ignore_index=True)

    dfnew = pd.DataFrame()
    dfnew["time"] = pd.DatetimeIndex(df["DATE"])
    dfnew[
        ["wind_direction", "wind_directionq", "wtype", "wind_speed", "wind_speedq"]
    ] = df["WND"].str.split(",", expand=True)
    dfnew[["air_temp", "air_tempq"]] = df["TMP"].str.split(",", expand=True)
    dfnew[["dew_point", "dew_pointq"]] = df["DEW"].str.split(",", expand=True)
    dfnew[["lprecip_period", "lprecip_depth", "lprecip_cond", "lprecipq"]] = df[
        "AA1"
    ].str.split(",", expand=True)
    if "AT1" in df:
        dfnew[["dpwo_source", "dpwo_type", "dpwo_abbrev", "dpwoq"]] = df[
            "AT1"
        ].str.split(",", expand=True)
    dfnew[["sea_level_pressure", "sea_level_pressureq"]] = df["SLP"].str.split(
        ",", expand=True
    )
    if "MA1" in df:
        dfnew[["XXXX", "XXXX", "atmospheric_pressure", "atmospheric_pressureq"]] = df[
            "MA1"
        ].str.split(",", expand=True)
        dfnew.drop(columns="XXXX", inplace=True)

    dfnew.set_index("time", inplace=True)
    for v in dfnew:
        try:
            dfnew[v] = pd.to_numeric(dfnew[v])
        except ValueError as e:
            print("could not make numeric", v, e)
            dfnew[v] = dfnew[v].str.strip()
            continue

    for v in ["lprecip_cond"]:
        dfnew.loc[dfnew[v] == 9, v] = np.nan

    for v in ["lprecip_period"]:
        dfnew.loc[dfnew[v] == 99, v] = np.nan

    for v in ["wind_direction"]:
        dfnew.loc[dfnew[v] == 999, v] = np.nan

    for v in ["wind_speed", "air_temp", "dew_point", "lprecip_depth"]:
        dfnew.loc[dfnew[v] == 9999, v] = np.nan
        dfnew[v] = dfnew[v] / 10

    for v in ["sea_level_pressure", "atmospheric_pressure"]:
        if v in dfnew:
            dfnew.loc[dfnew[v] == 99999, v] = np.nan
            dfnew[v] = dfnew[v] / 10

    ds = dfnew.to_xarray()

    units = {
        "wind_speed": "m s-1",
        "air_temp": "degree_C",
        "sea_level_pressure": "millibar",
        "atmospheric_pressure": "millibar",
        "dew_point": "degree_C",
    }
    for k in units:
        if k in ds:
            ds[k].attrs["units"] = units[k]

    for v in [
        "STATION",
        "SOURCE",
        "LATITUDE",
        "LONGITUDE",
        "ELEVATION",
        "NAME",
        "REPORT_TYPE",
        "CALL_SIGN",
        "QUALITY_CONTROL",
    ]:
        if df[v].eq(df[v][0]).all():
            ds.attrs[v.lower()] = df[v][
                0
            ]  # NAME is a reserved netCDF4 attr, so make them all lowercase:  https://github.com/Unidata/netcdf4-python/issues/1020
            print(v, "all equal")
        else:
            print(v, "not all equal")
            # ds[v] = xr.DataArray(df[v], dims="time")
            # df[v] = df[v].str.strip()
            print(f"  unique values for {v} are", df[v].unique())

    ds = ds.drop_duplicates(dim="time")

    return ds
