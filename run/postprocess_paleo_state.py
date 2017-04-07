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

# set up the option parser
parser = ArgumentParser()
parser.description = "Postprocessing files."
parser.add_argument("INDIR", nargs=1,
                    help="main directory", default=None)

options = parser.parse_args()
idir = options.INDIR[0]

# create logger
logger = logging.getLogger('postprocess')
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
