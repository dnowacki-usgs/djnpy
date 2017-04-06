import requests
import numpy as np
import datetime as dt
import pytz
from dateutil import parser

def get_coops_data(station, start_date, end_date, product='hourly_height', units='metric'):
    """
    Other options for product include 'water_level', 'hourly_height'
    units can be 'english' or 'metric'
    """
    url = 'http://tidesandcurrents.noaa.gov/api/datagetter?product=' \
    + product \
    + '&application=NOS.COOPS.TAC.WL&begin_date=' \
    + str(start_date) \
    + '&end_date=' \
    + str(end_date) \
    + '&datum=MLLW&station=' \
    + str(station) \
    + '&time_zone=GMT&units=' \
    + units \
    + '&format=json'

    payload = requests.get(url).json()

    if 'error' in payload.keys():
        raise ValueError('Error in returning dataset. Time requested too long?')

    t = []
    v = []
    for n in range(len(payload['data'])):
        t.append(pytz.utc.localize(parser.parse(payload['data'][n]['t'])))
        try:
            v.append(float(payload['data'][n]['v']))
        except:
            v.append(np.nan)
    t = np.array(t)
    v = np.array(v)

    n = {}
    n['dn'] = t
    n['v'] = v

    return n
