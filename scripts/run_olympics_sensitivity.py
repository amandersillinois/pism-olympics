#!/usr/bin/env python
import numpy as np
from osgeo import gdal, osr
import itertools
try:
    import subprocess32 as sub
except:
    import subprocess as sub
import csv
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
    stations = []
    if layer.GetGeomType() == 1:
        for pt in range(0, cnt):
            feature = layer.GetFeature(pt)
            geometry = feature.GetGeometryRef()
            # Transform to latlon if needed
            if not srs.IsGeographic():
                geometry.TransformTo(srs_geo)                
            point = geometry.GetPoint()
            stations.append([point[0], point[1]])

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
    parser.add_argument('-i', dest='in_file',
                        help='Gdal-readable DEM', default=None)
    parser.add_argument('-g', '--shape_file', dest='shp_file',
                        help='Point shape file with rain gauges', default=None)
    parser.add_argument(
        '--background_precip',
        dest='P0',
        type=float,
        help='Background precipitation rate [mm hr-1].',
        default=0.)
    parser.add_argument('--precip_scale_factor', dest='P_scale', type=float,
                        help='Precipitation scale factor.', default=1.)
    parser.add_argument('--no_trunc', dest='truncate', action='store_false',
                        help='Do not truncate precipitation.', default=True)
    parser.add_argument('--latitude', dest='lat', type=float,
                        help='Latitude to compute Coriolis term.', default=45.)
    parser.add_argument('--tau_c', dest='tau_c', type=float,
                        help='conversion time [s].', default=1000)
    parser.add_argument('--tau_f', dest='tau_f', type=float,
                        help='fallout time [s].', default=1000)
    parser.add_argument('--moist_stability', dest='Nm', type=float,
                        help='moist stability frequency [s-1].', default=0.005)
    parser.add_argument('--vapor_scale_height', dest='Hw', type=float,
                        help='Water vapor scale height [m].', default=2500)
    parser.add_argument(
        '--wind_direction',
        dest='direction',
        type=float,
        help='Direction from which the wind is coming.',
        default=270)
    parser.add_argument('--wind_magnitude', dest='magnitude', type=float,
                        help='Magnitude of wind velocity [m/s].', default=15)

    options = parser.parse_args()
    in_file = options.in_file
    shp_file = options.shp_file
    direction = options.direction
    lat = options.lat
    magnitude = options.magnitude
    tau_c = options.tau_c
    tau_f = options.tau_f
    truncate = options.truncate
    Nm = options.Nm
    Hw = options.Hw
    P0 = options.P0
    P_scale = options.P_scale
    ounits = 'm yr-1'
    gd = ReadRaster(in_file)
    X = gd.X
    Y = gd.Y
    Orography = gd.RasterArray

    Theta_m = -6.5     # K / km
    rho_Sref = 7.4e-3  # kg m-3
    gamma = -5.8       # K / km

    tau_c_values = [1000]
    tau_f_values = [1000]
    Nm_values = [0.005]
    Hw_values = [3200]
    magnitude_values = [15, 20]
    direction_values = [248]

    if shp_file is not None:
        stations = read_shapefile(shp_file)
        print stations
        
    combinations = list(itertools.product(tau_c_values, tau_f_values, Nm_values, Hw_values, magnitude_values, direction_values))
    with open('ltop_sensitivity.csv', 'w') as f:
        csvwriter = csv.writer(f)
        # csvwriter.writerow(['lon', 'lat', 'tau_c', 'tau_f', 'Nm', 'Hw', 'wind_speed', 'wind_direction', 'precip', 'precip_unit'])
        for combination in combinations:

            tau_c, tau_f, Nm, Hw, magnitude, direction = combination
            out_name = '_'.join(['ltop_olymics_precip', ounits.replace(' ', '_'), 'tauc', str(tau_c), 'tauf', str(tau_f), 'Nm', str(Nm), 'Hw', str(Hw), 'mag', str(magnitude), 'dir', str(direction)])
            print('Running combination {}'.format(out_name))

            physical_constants = dict()
            physical_constants['tau_c'] = tau_c      # conversion time [s]
            physical_constants['tau_f'] = tau_f      # fallout time [s]
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
            for station in stations:
                lon, lat =  station[0:2]
                cmd = ['gdallocationinfo', '-wgs84', '-valonly', out_file, str(lon), str(lat)]
                result = sub.Popen(cmd, stdout=sub.PIPE)
                precip = float(result.stdout.read().rstrip('\n'))
                precips.append(precip)
            row = [x for x in combination] + precips + [ounits]
            print row
            csvwriter.writerow(row)
