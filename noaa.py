import requests

def get_coops_hourly():
    url = 'http://tidesandcurrents.noaa.gov//api/datagetter?product=hourly_height&application=NOS.COOPS.TAC.WL&begin_date=' \
    + str(20160801) \
    + '&end_date=' \
    + str(20170201) \
    + '&datum=MLLW&station=' \
    + str(8739803) \
    + '&time_zone=GMT&units=english&format=json'

    payload = requests.get(url).json()

    return payload

payload = get_coops_hourly()
# %%

print len(payload['data'])
import numpy as np
import datetime as dt
payload['data'][n]['t']

len(payload['data'])


t = []
for n in range(len(payload['data'])):
    t.append(dt.datetime(payload['data'][n]['t']))
