import datetime
import unittest
from .. import matlabtools
import numpy as np


class TestMatlabtools(unittest.TestCase):
    def test_matlab2datetime(self):
        # Output from Matlab R2021b
        # >> datenum('2021-09-15 12:00')
        #
        # ans =
        #
        #    7.3841e+05
        # (full value is 738414.5)

        assert matlabtools.matlab2datetime(738414.5, tz=False) == datetime.datetime(
            2021, 9, 15, 12
        )

    def test_datetime2matlab(self):

        assert (
            matlabtools.datetime2matlab(datetime.datetime(2021, 9, 15, 12)) == 738414.5
        )


if __name__ == "__main__":
    unittest.main()
