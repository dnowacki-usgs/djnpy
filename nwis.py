
def nwis_json(site='01646500', parm='00065', start=None, end=None, period=None, freq='iv'):
    '''Obtain NWIS data via JSON.

    Args:
        site (str): NWIS site code. Default '01646500'
        parm (str): Parameter code. Default '00065'
            e.g. '00065' for water level, '00060' for discharge.
        freq: Data frequency. Valid values are:
            'iv' for unit-value data
            'dv' for daily average data

        Date ranges (default one day) can be specified using either:
            period (str): Duration following ISO-8601 duration format
                e.g. 'P1D' for one day
            start (str) and stop (str): begin and end dates/times following
                ISO-8601 format

    Returns:
        nwis, a dict of numpy arrays, containing the following fields:
            ['dn']: datetime
            ['val']: values corresponding to datetimes
            ['sitename']: site name
            ['sitecode']: NWIS site code

    More info about the URL format is at http://waterservices.usgs.gov

    dnowacki@usgs.gov 2016-07
    '''

    import requests
    import numpy as np
    from dateutil import parser
    import pytz

    if period is None and start is None and end is None:
        period='P1D'

    if period is None:
        url = 'http://waterservices.usgs.gov/nwis/' + freq + '/?format=json&sites=' + str(site) + '&startDT=' + start + '&endDt=' + end + '&parameterCd=' + str(parm)
    else:
        url = 'http://waterservices.usgs.gov/nwis/' + freq + '/?format=json&sites=' + str(site) + '&period=' + period + '&parameterCd=' + str(parm)
    payload = requests.get(url).json()
    v = payload['value']['timeSeries'][0]['values'][0]['value']
    nwis = {}
    nwis['dnlocal'] = np.array([parser.parse(v[i]['dateTime']) for i in range(len(v))])
    # Convert local time to UTC
    nwis['dn'] = [x.astimezone(pytz.utc) for x in nwis['dnlocal']]
    nwis['sitename'] = payload['value']['timeSeries'][0]['sourceInfo']['siteName']
    nwis['sitecode'] = payload['value']['timeSeries'][0]['sourceInfo']['siteCode'][0]['value']
    nwis['val'] = np.array([float(v[i]['value']) for i in range(len(v))])
    nwis['val'][nwis['val'] == payload['value']['timeSeries'][0]['variable']['noDataValue']] = np.nan

    return nwis
