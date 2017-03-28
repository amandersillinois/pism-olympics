#!/usr/bin/env python
import numpy as np
from osgeo import gdal, osr
import itertools
try:
    import subprocess32 as sub
except:
    import subprocess as sub
import csv
import operator
import logging
logger = logging.getLogger('LTOP')

from linear_orog_precip import OrographicPrecipitation, ReadRaster, array2raster

np.seterr(divide='ignore', invalid='ignore')


def read_shapefile(filename):
    '''
    Reads lat / lon from a ESRI shape file.

    Paramters
    ----------
    filename: filename of ESRI shape file.

    Returns
    -------
    lat, lon: array_like coordinates

    '''
    import ogr
    import osr
    driver = ogr.GetDriverByName('ESRI Shapefile')
    data_source = driver.Open(filename, 0)
    layer = data_source.GetLayer(0)
    srs = layer.GetSpatialRef()
    if not srs.IsGeographic():
        print('''Spatial Reference System in % s is not latlon. Converting.'''
              % filename)
        # Create spatialReference, EPSG 4326 (lonlat)
        srs_geo = osr.SpatialReference()
        srs_geo.ImportFromEPSG(4326)
    cnt = layer.GetFeatureCount()
    stations = dict()
    if layer.GetGeomType() == 1:
        for pt in range(0, cnt):
            feature = layer.GetFeature(pt)
            geometry = feature.GetGeometryRef()
            # Transform to latlon if needed
            if not srs.IsGeographic():
                geometry.TransformTo(srs_geo)                
            point = geometry.GetPoint()
            data = dict()
            data['point'] = [point[0], point[1]]
            data['p_obs'] = feature.precip
            data['p_obs_units'] = 'm yr-1'
            stations[pt] = data

    return stations



if __name__ == "__main__":
    print('Linear Theory Orographic Precipitation Model by Smith & Barstad (2004)')

    import itertools
    from argparse import ArgumentParser

    import logging
    import logging.handlers

    # create logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    # create file handler which logs even debug messages
    fh = logging.handlers.RotatingFileHandler('ltop.log')
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(module)s:%(lineno)d - %(message)s')
    # add formatter to ch and fh
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)
    logger.addHandler(fh)
    parser = ArgumentParser()
    parser.add_argument('-o', dest='csv_file',
                        help='CSV output file name', default='ltop_olympics_default.csv')
    parser.add_argument('-i', dest='in_file',
                        help='Gdal-readable DEM', default=None)
    parser.add_argument('-g', '--shape_file', dest='shp_file',
                        help='Point shape file with rain gauges', default=None)

    options = parser.parse_args()
    in_file = options.in_file
    shp_file = options.shp_file
    csv_file = options.csv_file
    lat = 0.
    truncate = True
    ounits = 'm yr-1'
    gd = ReadRaster(in_file)
    X = gd.X
    Y = gd.Y
    Orography = gd.RasterArray

    Theta_m = -6.5     # K / km
    rho_Sref = 7.4e-3  # kg m-3
    gamma = -5.8       # K / km

    tau_values = [1750]
    Nm_values = [0.007]
    Hw_values = [2600]
    magnitude_values = [15]
    direction_values = [220 ]
    precip_scale_factor_values = [ 0.09, 0.11, 0.125]
    background_precip_values_0 = [1]  # mm hr-1
    background_precip_values_1 = [0.057, 0.114]   

    if shp_file is not None:
        stations = read_shapefile(shp_file)
        
    combinations = list(itertools.product(tau_values, Nm_values, Hw_values, magnitude_values, direction_values, precip_scale_factor_values, background_precip_values_0, background_precip_values_1))
    n_exp = len(combinations)
    st_data = dict()
    with open(csv_file, 'w') as f:
        csvwriter = csv.writer(f)
        a = ['p_exp_{}'.format(x) for x in range(n_exp)]
        b = ['p_exp_{}_str'.format(x) for x in range(n_exp)]
        exp_list = []
        for k in range(len(a)):
            exp_list.append(a[k])
            exp_list.append(b[k])
        csvwriter.writerow(['lon', 'lat', 'p_obs', 'p_obs_units'] + exp_list
                           + ['p_exp_units'])
        for exp_no, combination in enumerate(combinations):

            tau, Nm, Hw, magnitude, direction, P_scale, P0, P1 = combination
            out_name = '_'.join(['ltop_olymics_precip', ounits.replace(' ', '_'), 'tau', str(tau), 'Nm', str(Nm), 'Hw', str(Hw), 'mag', str(magnitude), 'dir', str(direction), 'ps', str(P_scale), 'p0', str(P0), 'p1', str(P1)])
            print('Running exp {} {}'.format(exp_no, out_name))

            physical_constants = dict()
            physical_constants['tau_c'] = tau  # conversion time [s]
            physical_constants['tau_f'] = tau  # fallout time [s]
            physical_constants['f'] = 2 * 7.2921e-5 * \
                np.sin(lat * np.pi / 180)  # Coriolis force
            physical_constants['Nm'] = Nm   # moist stability frequency [s-1]
            physical_constants['Cw'] = rho_Sref * Theta_m / \
                gamma  # uplift sensitivity factor [kg m-3]
            physical_constants['Hw'] = Hw         # vapor scale height
            # x-component of wind vector [m s-1]
            physical_constants['u'] = -np.sin(direction * 2 * np.pi / 360) * magnitude
            # y-component of wind vector [m s-1]
            physical_constants['v'] = np.cos(direction * 2 * np.pi / 360) * magnitude
            physical_constants['P0'] = P0   # background precip [mm hr-1]
            physical_constants['P1'] = P1   # background precip [mm hr-1]
            physical_constants['P_scale'] = P_scale   # precip scale factor [1]

            OP = OrographicPrecipitation(
                X,
                Y,
                Orography,
                physical_constants,
                truncate=truncate,
                ounits=ounits)
            P = OP.P
            units = OP.P_units

            out_file = out_name + '.tif'

            array2raster(out_file, gd.geoTrans, gd.proj4, units, P)
            precips = []
            for k in stations:
                station = stations[k]
                lon, lat =  station['point'][0:2]
                cmd = ['gdallocationinfo', '-wgs84', '-valonly', out_file, str(lon), str(lat)]
                result = sub.Popen(cmd, stdout=sub.PIPE)
                precip = float(result.stdout.read().rstrip('\n'))
                precips.append(precip)
            st_data[exp_no] = [precips, ';'.join([str(x) for x in combination])]

        for k in stations:
            station = stations[k]
            station['p_exp'] = dict()
            for exp in st_data:
                print k, exp , st_data[exp][0][k], st_data[exp][1]
                station['p_exp'][exp] = [st_data[exp][0][k], st_data[exp][1]]
            lon, lat =  station['point']
            row = [lon, lat, station['p_obs'], station['p_obs_units']] + sum(station['p_exp'].values(), []) + [ounits]
            csvwriter.writerow(row)
        print('Done writing {}'.format(csv_file))
