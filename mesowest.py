import xarray as xr
import pandas as pd
import numpy as np


def mesowest_to_xarray(data):
    ds = xr.Dataset()
    ds["time"] = pd.to_datetime(data["STATION"][0]["OBSERVATIONS"]["date_time"])
    ds["time"] = pd.DatetimeIndex(ds["time"].values, tz=None)
    for k, suffix in zip(
        ["pressure", "wind_speed", "wind_direction", "wind_gust"],
        ["_set_1d", "_set_1", "_set_1", "_set_1"],
    ):
        ds[k] = xr.DataArray(
            np.array(data["STATION"][0]["OBSERVATIONS"][k + suffix]).astype(float),
            dims="time",
        )
        ds[k].attrs["units"] = data["UNITS"][k]
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
