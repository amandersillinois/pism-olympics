#!/usr/bin/env python
# Copyright (C) 2017 Andy Aschwanden

import os
try:
    import subprocess32 as sub
except:
    import subprocess as sub
from glob import glob
import numpy as np
import gdal
from nco import Nco
nco = Nco()
from nco import custom as c
import logging
import logging.handlers
from argparse import ArgumentParser

from netCDF4 import Dataset as NC

def permute_nc(variable, output_order=('time', 'station', 'profile', 'z', 'zb')):
    '''
    Permute dimensions of a NetCDF variable to match the output
    storage order.

    Parameters
    ----------
    variable : a netcdf variable
               e.g. thk = nc.variables['thk']
    output_order: dimension tuple (optional)
                  default ordering is ('time', 'z', 'zb', 'y', 'x')

    Returns
    -------
    var_perm : array_like
    '''

    input_dimensions = variable.dimensions

    # filter out irrelevant dimensions
    dimensions = filter(lambda x: x in input_dimensions,
                        output_order)

    # create the mapping
    mapping = map(lambda x: dimensions.index(x),
                  input_dimensions)

    if mapping:
        return np.transpose(variable[:], mapping)
    else:
        return variable[:]  # so that it does not break processing "mapping"

def permute(array,
               input_order=('time', 'station', 'profile', 'z', 'zb'),
               output_order=('station', 'time', 'profile', 'z', 'zb')):
    '''
    Permute dimensions of an array to match the output
    storage order.

    Parameters
    ----------
    array: array_like
    input_order: dimension tuple (optional)
                  default ordering is ('time', 'station', 'profile', 'z', 'zb')'
    output_order: dimension tuple (optional)
                  default ordering is ('time', 'station', 'profile', 'z', 'zb')'

    Returns
    -------
    var_perm : array_like
    '''

    input_dimensions = input_order

    # filter out irrelevant dimensions
    dimensions = filter(lambda x: x in input_dimensions,
                        output_order)

    # create the mapping
    mapping = map(lambda x: dimensions.index(x),
                  input_dimensions)

    return np.transpose(array, mapping)


def compute_normal_speed(ifile):
    '''
    Compute normal-to-profile speeds
    '''
    nc = NC(ifile, 'a')
    ux = permute_nc(nc.variables['uvelsurf'])
    vy = permute_nc(nc.variables['vvelsurf'])
    nx = permute_nc(nc.variables['nx'])
    ny = permute_nc(nc.variables['ny'])
    dims = nc.variables['uvelsurf'].dimensions
    fill_value = nc.variables['uvelsurf']._FillValue
    nc.createVariable('velsurf_normal', 'd', dimensions=dims, fill_value=fill_value)
    vn = ux * nx + vy * ny
    # filter out irrelevant dimensions
    input_dims = filter(lambda x: x in dims,
                        ('station', 'time', 'profile', 'z', 'zb'))
    nc.variables['velsurf_normal'][:] = permute(vn,
                                                input_order=input_dims,
                                                output_order=dims)
    nc.variables['velsurf_normal'].units = 'm year-1'
    nc.close()


# set up the option parser
parser = ArgumentParser()
parser.description = "Postprocessing files."
parser.add_argument("INDIR", nargs=1,
                    help="main directory", default=None)

options = parser.parse_args()
idir = options.INDIR[0]

# create logger
logger = logging.getLogger('prepare_velocity_observations')
logger.setLevel(logging.DEBUG)

# create file handler which logs even debug messages
fh = logging.handlers.RotatingFileHandler('prepare_velocity_observations.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(module)s:%(lineno)d - %(message)s')

# add formatter to ch and fh
ch.setFormatter(formatter)
fh.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
logger.addHandler(fh)


gdal_gtiff_options = gdal.TranslateOptions(format='GTiff', outputSRS='EPSG:26710')

# Process experiments
dir_gtiff = 'processed_gtiff'
dir_nc = 'processed_nc'

for dir_processed in (dir_gtiff, dir_nc):
    if not os.path.isdir(os.path.join(idir, dir_processed)):
        os.mkdir(os.path.join(idir, dir_processed))

pvars = ('thk', 'usurf', 'velsurf_mag')
fill_value = -2e9
v_str = ' '.join('='.join([x, str(fill_value) + ';']) for x in pvars)
ncap2_str = 'where(thk<10) {{ {} }};'.format(v_str)
exp_files = glob(os.path.join(idir, 'state', '*.nc'))
for exp_file in exp_files:
    exp_basename =  os.path.split(exp_file)[-1].split('.nc')[0]
    exp_nc_wd = os.path.join(idir, dir_nc, exp_basename + '.nc')
    exp_gtiff_wd = os.path.join(idir, dir_gtiff, exp_basename + '.tif')
    logger.info('masking variables where ice thickness < 10m')
    nco.ncap2(input='-s "{}" {}'.format(ncap2_str, exp_file), output=exp_nc_wd, overwrite=True)
    opt = [c.Atted(mode="o", att_name="_FillValue", var_name=myvar, value=fill_value) for myvar in pvars]
    nco.ncatted(input=exp_nc_wd, options=opt)
    for mvar in pvars:
        m_exp_nc_wd = 'NETCDF:{}:{}'.format(exp_nc_wd, mvar)
        m_exp_gtiff_wd = os.path.join(idir, dir_gtiff, mvar + '_' + exp_basename + '.tif')
        logger.info('Converting variable {} to GTiff and save as {}'.format(mvar, m_exp_gtiff_wd))
        gdal.Translate(m_exp_gtiff_wd, m_exp_nc_wd, options=gdal_gtiff_options)
