#!/usr/bin/env python

from argparse import ArgumentParser
from netCDF4 import Dataset as CDF
    
# Set up the option parser
parser = ArgumentParser()
parser.description = '''Script adds CF-conforming mapping variable for EPSG:26710.'''
parser.add_argument("FILE", nargs=1)
options = parser.parse_args()
args = options.FILE

nc = CDF(args[0], 'a')

if not "mapping" in nc.variables.keys():
    mapping = nc.createVariable("mapping", 'c')
else:
    mapping = nc.variables("mapping")
    
mapping.longitude_of_central_meridian = -123.0
mapping.false_easting = 500000.0
mapping.false_northing =  0.0
mapping.grid_mapping_name = "transverse_mercator"
mapping.inverse_flattening = 294.978698213898
mapping.latitude_of_projection_origin = 0.0
mapping.scale_factor_at_central_meridian = 0.9996
mapping.semi_major_axis = 6378206.4
mapping.unit = "metre"

for var in nc.variables.keys():
    if (nc.variables[var].ndim >= 2 and  var not in ['lat', 'lon', 'lat_bnds', 'lon_bnds', 'lat_bounds', 'lon_bounds']):
        nc.variables[var].grid_mapping = "mapping"

nc.close()
