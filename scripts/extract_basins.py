#!/usr/bin/env python
# Copyright (C) 2016-17 Andy Aschwanden


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import ocgis
import os
from cdo import Cdo
cdo = Cdo()

import logging
import logging.handlers

# create logger
logger = logging.getLogger('extract_basins')
logger.setLevel(logging.DEBUG)

# create file handler which logs even debug messages
fh = logging.handlers.RotatingFileHandler('extract.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter
info_formatter = logging.Formatter('%(message)s')
debug_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(module)s:%(lineno)d - %(message)s')

# add formatter to ch and fh
ch.setFormatter(info_formatter)
fh.setFormatter(debug_formatter)

# add ch to logger
logger.addHandler(ch)
logger.addHandler(fh)

script_path = os.path.dirname(os.path.realpath(__file__))
default_basin_file = 'GRE_Basins_IMBIE2_v1.3_ext.shp'

def extract_basins():
    '''
    Extract basin using OCGIS
    '''
    logger.info('Extracting basin {}'.format(basin))

    if GEOM is None:
        select_ugid = None
    else:
        select_geom = filter(lambda x: x['properties']['basin'] == basin,
                             ocgis.GeomCabinetIterator(path=SHAPEFILE_PATH))
        ## this argument must always come in as a list
        select_ugid = [select_geom[0]['properties']['UGID']]
    ## parameterize the operations to be performed on the target dataset
    ops = ocgis.OcgOperations(dataset=rd,
                              geom=SHAPEFILE_PATH,
                              aggregate=False,
                              snippet=False,
                              select_ugid=select_ugid,
                              output_format=output_format,
                              prefix=prefix,
                              dir_output=odir)
    ret = ops.execute()

def calculate_time_series():
    '''
    Calculate scalar time series with CDO
    '''
    ifile = os.path.join(odir, prefix, prefix + '.nc')
    scalar_ofile = os.path.join(odir, prefix, '.'.join(['_'.join(['scalar_fldsum', prefix]), 'nc']))
    logger.info('Calculating field sum and saving to \n {}'.format(scalar_ofile))
    cdo.fldsum(input='-fldsum -selvar,{} {}'.format(','.join(mvar for mvar in mvars), ifile), output=scalar_ofile, overwrite=True)
    scalar_sum_ofile = os.path.join(odir, prefix, '.'.join(['_'.join(['cumsum', prefix]), 'nc']))
    logger.info('Calculating cumulative time sum and saving to \n {}'.format(scalar_sum_ofile))
    cdo.chname(','.join(','.join([mvar, mvars_dict[mvar]]) for mvar in mvars), input='-setattribute,{} -timcumsum {}'.format(','.join('@'.join([mvar, 'units=Gt']) for mvar in mvars), scalar_ofile), output=scalar_sum_ofile, overwrite=True)
    runmean_ofile = os.path.join(odir, prefix, '.'.join(['_'.join(['runmean_{}yr'.format(runmean_steps), prefix]), 'nc']))
    logger.info('Calculating running mean and saving to \n {}'.format(runmean_ofile))
    cdo.runmean('{}'.format(runmean_steps), input=scalar_ofile, output=runmean_ofile, overwrite=True)
    abs_anomaly_ofile = os.path.join(odir, prefix, '.'.join(['_'.join(['abs_anomaly_runmean_{}yr'.format(runmean_steps), prefix]), 'nc']))
    logger.info('Calculating anomalies and saving to \n {}'.format(abs_anomaly_ofile))
    # Choose sign such that a positive anomaly means an increase in discharge
    cdo.mulc('-1', input='-sub {} -timmean -seltimestep,1/{} {}'.format(runmean_ofile, runmean_steps, scalar_ofile), output=abs_anomaly_ofile, overwrite=True)
    rel_anomaly_ofile = os.path.join(odir, prefix, '.'.join(['_'.join(['rel_anomaly_runmean_{}yr'.format(runmean_steps), prefix]), 'nc']))
    logger.info('Calculating anomalies and saving to \n {}'.format(rel_anomaly_ofile))
    cdo.div(input=' {} -timmean -seltimestep,1/{} {}'.format(runmean_ofile, runmean_steps, scalar_ofile), output=rel_anomaly_ofile, overwrite=True)


# set up the option parser
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.description = "Extract basins from continental scale files."
parser.add_argument("FILE", nargs=1)
parser.add_argument("--o_dir", dest="odir",
                    help="output directory", default='.')
parser.add_argument("--shape_file", dest="shape_file",
                    help="Path to shape file with basins", default=os.path.join(script_path, default_basin_file))
parser.add_argument("-v", "--variable", dest="VARIABLE",
                    help="Comma-separated list of variables to be extracted. By default, all variables are extracted.", default=None)
parser.add_argument("--no_extraction", dest="no_extraction", action="store_true",
                    help="Don't extract basins", default=False)
parser.add_argument("--no_timeseries", dest="no_timeseries", action="store_true",
                    help="Don't calculate time-series", default=False)

options = parser.parse_args()
no_extraction = options.no_extraction
no_timeseries = options.no_timeseries

URI = options.FILE[0]
SHAPEFILE_PATH = options.shape_file
if options.VARIABLE is not None:
    VARIABLE=options.VARIABLE.split(',')
else:
    VARIABLE=options.VARIABLE

odir = options.odir
if not os.path.isdir(odir):
    os.mkdir(odir)

ocgis.env.OVERWRITE = True

# Output name
savename=URI[0:len(URI)-3] 

## set the output format to convert to
output_format = 'nc'

## we can either subset the data by a geometry from a shapefile, or convert to
## geojson for the entire spatial domain. there are other options here (i.e. a
## bounding box for tiling or a Shapely geometry).
GEOM = SHAPEFILE_PATH

mvars_dict = {'tendency_of_ice_mass': 'ice_mass',
             'tendency_of_ice_mass_due_to_flow': 'flow_cumulative',
             'tendency_of_ice_mass_due_to_conservation_error': 'conservation_error_cumulative',
             'tendency_of_ice_mass_due_to_basal_mass_flux': 'basal_mass_flux_cumulative',
             'tendency_of_ice_mass_due_to_surface_mass_flux': 'surface_mass_flux_cumulative',
             'tendency_of_ice_mass_due_to_discharge': 'discharge_cumulative'}
mvars = mvars_dict.keys()
cvars = ['pism_config']
basins = ('Hoh', 'Quinault', 'Elwha', 'Queets')
runmean_steps = 10

rd = ocgis.RequestDataset(uri=URI, variable=VARIABLE)
for basin in basins:
    prefix = 'b_{basin}_{savename}'.format(basin=basin, savename=savename)
    rd = ocgis.RequestDataset(uri=URI, variable=VARIABLE)
    if not no_extraction:
        extract_basins()
    if not no_timeseries:
        calculate_time_series()
