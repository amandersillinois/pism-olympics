#!/usr/bin/env python
# Copyright (C) 2016-17 Andy Aschwanden


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import ocgis
import os

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

options = parser.parse_args()

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

basins = ('Hoh', 'Quinault', 'Elwha', 'Queets')

rd = ocgis.RequestDataset(uri=URI, variable=VARIABLE)
for basin in basins:
    prefix = 'b_{basin}_{savename}'.format(basin=basin, savename=savename)
    rd = ocgis.RequestDataset(uri=URI, variable=VARIABLE)
    extract_basins()
