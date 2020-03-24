from __future__ import division
import unittest
import nwis
import numpy as np

class TestNwis(unittest.TestCase):

    def test_early_nwis_date(self):
        # try to get it to fail on Windows
        Q = nwis.nwis_json(12101500, parm='00060', start='1914-05-01', end='1914-06-01', freq='dv', xarray=True)

        assert Q['val'][0] == 3070


if __name__ == '__main__':
    unittest.main()
