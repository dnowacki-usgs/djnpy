import numpy as np
import pandas as pd
import xarray as xr


def mesowest_to_xarray(data):
    ds = xr.Dataset()
    ds["time"] = pd.to_datetime(data["STATION"][0]["OBSERVATIONS"]["date_time"])
    ds["time"] = pd.DatetimeIndex(ds["time"].values, tz=None)
    for k in data["STATION"][0]["OBSERVATIONS"].keys():
        if k == "date_time":
            continue
        var = k.split("_set_")
        try:
            ds[var[0]] = xr.DataArray(
                np.array(data["STATION"][0]["OBSERVATIONS"][k]).astype(float),
                dims="time",
            )
        except ValueError:  # for variables which can't be converted to float
            ds[var[0]] = xr.DataArray(
                np.array(data["STATION"][0]["OBSERVATIONS"][k]),
                dims="time",
            )
        ds[var[0]].attrs["units"] = data["UNITS"][var[0]]
    for k in [
        "STATUS",
        "MNET_ID",
        "ELEVATION",
        "NAME",
        "STID",
        "ELEV_DEM",
        "LONGITUDE",
        "STATE",
        "LATITUDE",
        "TIMEZONE",
        "ID",
    ]:
        if k == "NAME":
            ds.attrs["STATION_NAME"] = data["STATION"][0][k]
        else:
            ds.attrs[k] = data["STATION"][0][k]

    return ds
