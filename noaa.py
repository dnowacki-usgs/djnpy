import requests
import numpy as np
import datetime as dt
import pandas as pd
import pytz
import xarray as xr
from dateutil import parser
import urllib


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

    payload = requests.get(url).json()

    if "error" in payload.keys():
        raise ValueError("Error in returning dataset: " + payload["error"]["message"])

    t = []
    v = []
    wind = {"s": [], "d": [], "dr": [], "g": [], "f": []}
    cp = {"Speed": [], "Bin": [], "Direction": [], "Depth": []}
    cur = {"s": [], "d": [], "b": []}
    if (
        product == "water_level"
        or product == "hourly_height"
        or product == "air_pressure"
    ):
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
                    except:
                        wind[k].append(np.nan)
        elif product == "monthly_mean":
            for k in monthly.keys():
                try:
                    monthly[k].append(float(d[n][k]))
                except:
                    monthly[k].append(np.nan)
        elif product == "currents_predictions":
            for k in cp:
                cp[k].append(float(d[n][k]))
        elif product == "currents":
            for k in cur:
                cur[k].append(float(d[n][k]))
        else:
            try:
                v.append(float(d[n]["v"]))
            except:
                v.append(np.nan)

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
