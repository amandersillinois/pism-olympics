#!/usr/bin/env python
import numpy as np
from osgeo import gdal, osr
import itertools
import logging
logger = logging.getLogger('LTOP')
from linear_orog_precip import OrographicPrecipitation, ReadRaster, array2raster

np.seterr(divide='ignore', invalid='ignore')



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
    if in_file is not None:
        gd = ReadRaster(in_file)
        X = gd.X
        Y = gd.Y
        Orography = gd.RasterArray
    else:
        # Reproduce Fig 4c in SB2004
        dx = dy = 750.
        x, y = np.arange(-100e3, 200e3, dx), np.arange(-150e3, 150e3, dy)
        h_max = 500.
        x0 = -25e3
        y0 = 0
        sigma_x = sigma_y = 15e3
        X, Y = np.meshgrid(x, y)
        Orography = h_max * \
            np.exp(-(((X - x0)**2 / (2 * sigma_x**2)) + ((Y - y0)**2 / (2 * sigma_y**2))))

    Theta_m = -6.5     # K / km
    rho_Sref = 7.4e-3  # kg m-3
    gamma = -5.8       # K / km

    tau_c_values = [1000]
    tau_f_values = [1000]
    Nm_values = [0.005]
    Hw_values = [3200]
    magnitude_values = [10, 15, 20, 25]
    direction_values = [244, 246, 248, 250]
    
    combinations = list(itertools.product(tau_c_values, tau_f_values, Nm_values, Hw_values, magnitude_values, direction_values))

    for combination in combinations:

        tau_c, tau_f, Nm, Hw, magnitude, direction = combination
        out_name = '_'.join(['ltop_olymics_precip', ounits, 'tauc', str(tau_c), 'tauf', str(tau_f), 'Nm', str(Nm), 'Hw', str(Hw), 'mag', str(magnitude), 'dir', str(direction)])
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
        
        if in_file is not None:
            array2raster(out_file, gd.geoTrans, gd.proj4, units, P)
        else:
            geoTrans = [0., OP.dx, 0., 0., 0., -OP.dy]
            array2raster(out_file, geoTrans, '', units, P)
