#!/usr/bin/env python

from netCDF4 import Dataset
import scipy.io as spio

# convert .mat file into simple netcdf file
def mat2nc(matnam, ncnam):
    m = spio.loadmat(matnam) 
    rg = Dataset(ncnam, 'w', format='NETCDF4')
    x = rg.createDimension('x', 736)
    y = rg.createDimension('y', 1555)
    xs = rg.createVariable('x', 'f4', ('x',), zlib=True)
    ys = rg.createVariable('y', 'f4', ('y',), zlib=True)

    for key in m.keys():
        if '__' not in key: # ignore keys with two underscores (no real data there)
            print "Creating", key, "in", ncnam
            var = rg.createVariable(key, 'f4', ('y', 'x',), zlib=True)
            var[:] = m[key]

    rg.close()



mat2nc('disbot.mat', 'disbot.nc')
mat2nc('dissip.mat', 'dissip.nc')
mat2nc('diswcap.mat', 'diswcap.nc')
mat2nc('nplant.mat', 'nplant.nc')
