from __future__ import division
import numpy as np


def wavenumber(sigma, h):
    """
    Compute wavenumber from sigma and h
    k = wavenumber(sigma, h)
    k is the matrix of same size as sigma and h containing the calculated wave numbers

    sigma is the wave frequencies in rad/s
    h is the water depth

    sigma and h must be scalars,vectors or matricies of the same dimensions

    modified from R.Dalrymple's java code

    D. Nowacki rewritten in Python March 2015
    """

    g = 9.81

    a0 = (sigma ** 2 * h) / g
    b1 = 1.0 / np.tanh(a0 ** (3.0 / 4))
    a1 = a0 * (b1 ** (2.0 / 3))
    da1 = 1000.0

    d1 = np.ones(np.shape(h))

    while np.max(d1) == 1:
        d1 = abs(da1 / a1) > 0.00000001
        th = np.tanh(a1)
        ch = np.cosh(a1)
        f1 = a0 - a1 * th
        f2 = -a1 * (1.0 / ch) ** 2 - th
        da1 = -f1 / f2
        a1 = a1 + da1

    return a1 / h


def qkfs(omega, h):
    """
    Modified from Wiberg & Sherwood 2009; only does 3 iterations.
    Returns only k, not kh

    k = qkfs(omega, h)
    """

    g = 9.81
    x = omega ** 2 * h / g
    y = np.sqrt(x) * (x < 1) + x * (x >= 1)

    t = np.tanh(y)
    y = y - ((y * t - x) / (t + y * (1 - t ** 2)))
    t = np.tanh(y)
    y = y - ((y * t - x) / (t + y * (1 - t ** 2)))
    t = np.tanh(y)
    y = y - ((y * t - x) / (t + y * (1 - t ** 2)))

    return y / h
