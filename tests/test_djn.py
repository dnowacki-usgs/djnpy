from __future__ import division
import unittest
import math
import djn
import numpy as np

class TestDjn(unittest.TestCase):

    def setUp(self):
        self.s2 = math.sqrt(2)

    def test_uv2sd(self):
        u, v = [1, 0, -1, 0], [0, 1, 0, -1]
        s, d = djn.uv2sd(u, v)

        np.testing.assert_allclose(s, [1, 1, 1, 1])
        np.testing.assert_allclose(d, [90, 0, 270, 180])

        u, v = [self.s2/2, -self.s2/2, -self.s2/2, self.s2/2], \
               [self.s2/2, self.s2/2, -self.s2/2, -self.s2/2]
        s, d = djn.uv2sd(u, v)

        np.testing.assert_allclose(s, [1, 1, 1, 1])
        np.testing.assert_allclose(d, [45, 315, 225, 135])

    def test_sd2uv(self):
        s, d = [1, 1, 1, 1], [90, 0, 270, 180]
        u, v = djn.sd2uv(s, d)

        np.testing.assert_allclose(u, [1, 0, -1, 0], atol=1e-15)
        np.testing.assert_allclose(v, [0, 1, 0, -1], atol=1e-15)

        s, d = [1, 1, 1, 1], [45, 315, 225, 135]
        u, v = djn.sd2uv(s, d)

        np.testing.assert_allclose(u, [self.s2/2, -self.s2/2, -self.s2/2, self.s2/2])
        np.testing.assert_allclose(v, [self.s2/2, self.s2/2, -self.s2/2, -self.s2/2])

if __name__ == '__main__':
    unittest.main()
