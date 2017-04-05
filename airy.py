# Implement functions for linear wave theory, migrated from Matlab codes in
# Dropbox/matlab/+airy

def celerity(omega, d):
    import wavenumber as wn
    # c = celerity(omega, d)
    k = wn.qkfs(omega, d)
    c = omega / k
    return c

def energydens(a, rho=1025, g=9.81):
    # E = airy.energydens(rho, g, a)
    # return energy density following Airy wave theory
    # rho: water density (1025 kg/m^3 if argument is empty)
    # g: gravity (9.81 m/s^2 if argument is empty)
    # a: wave amplitude

    return 0.5 * rho * g * a**2

def groupvel(omega, d):
    import wavenumber  as wn
    import numpy as np

    # cg = airy.groupvel(omega, d)
    # return group velocity following Airy wave theory
    # omega: 2*pi/T
    # d: still water depth


    k = wn.qkfs(omega, d)

    # From Paul & Amos 2011 eq 9 (and elsewhere, obviously)
    # also repeated in Lowe et al 2005 eq 5
    return 0.5 * (1 + 2 * k * d / np.sinh(2 * k * d)) * omega / k


import sys
sys.path
