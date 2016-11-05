#!/usr/bin/env python
# Copyright (C) 2016 Andy Aschwanden

import itertools
from collections import OrderedDict
import os
try:
    import subprocess32 as sub
except:
    import subprocess as sub
from argparse import ArgumentParser
import sys
sys.path.append('../resources/')
from resources import *


# set up the option parser
parser = ArgumentParser()
parser.description = "Generating scripts for prognostic simulations."
parser.add_argument("-n", '--n_procs', dest="n", type=int,
                    help='''number of cores/processors. default=2.''', default=2)
parser.add_argument("-w", '--wall_time', dest="walltime",
                    help='''walltime. default: 12:00:00.''', default="12:00:00")
parser.add_argument("-q", '--queue', dest="queue", choices=list_queues(),
                    help='''queue. default=t1standard.''', default='normal')
parser.add_argument("--climate", dest="climate",
                    choices=['elev', 'paleo', 'present'],
                    help="Climate", default='paleo')
parser.add_argument("-d", "--domain", dest="domain",
                    choices=['olympics', 'olympics_mtns'],
                    help="sets the modeling domain", default='olympics')
parser.add_argument("--start_year", dest="start", type=float,
                    help="Start year", default=0)
parser.add_argument("--duration", dest="duration", type=float,
                    help="Duration", default=10)
parser.add_argument("--exstep", dest="exstep",
                    help="Spatial time series writing interval", default=1)
parser.add_argument("-f", "--o_format", dest="oformat",
                    choices=['netcdf3', 'netcdf4_parallel', 'pnetcdf'],
                    help="output format", default='netcdf3')
parser.add_argument("-g", "--grid", dest="grid", type=int,
                    choices=accepted_resolutions(),
                    help="horizontal grid resolution", default=1000)
parser.add_argument("-i", "--input_file", dest="input_file",
                    help="Input file to restart from", default=None)
parser.add_argument("--o_dir", dest="odir",
                    help="output directory. Default: current directory", default='test')
parser.add_argument("--o_size", dest="osize",
                    choices=['small', 'medium', 'big', 'big_2d'],
                    help="output size type", default='big_2d')
parser.add_argument("-s", "--system", dest="system",
                    choices=list_systems(),
                    help="computer system to use.", default='debug')
parser.add_argument("--bed_version", dest="bed_version",
                    help="Version of bed DEM.", default='1')
parser.add_argument("--stress_balance", dest="stress_balance",
                    choices=['sia', 'ssa+sia', 'ssa'],
                    help="stress balance solver", default='sia')


options = parser.parse_args()

nn = options.n
odir = options.odir
oformat = options.oformat
osize = options.osize
queue = options.queue
walltime = options.walltime
system = options.system

bed_version = options.bed_version
climate = options.climate
duration = options.duration
grid = options.grid
stress_balance = options.stress_balance

start = options.start
end  = start + options.duration
exstep = options.exstep
domain = options.domain
pism_exec = generate_domain(domain)
input_file = options.input_file
pism_dataname = 'pism_{domain}_{grid}m_v{version}.nc'.format(domain=domain.capitalize(),
                                                             grid=grid,
                                                             version=bed_version)
if not os.path.isdir(odir):
    os.mkdir(odir)

# Configuration File Setup
pism_config = 'olympics_config'
pism_config_nc = '.'.join([pism_config, 'nc'])
pism_config_cdl = os.path.join('../config', '.'.join([pism_config, 'cdl']))
# Anaconda libssl problem on chinook
if system in ('chinook'):
    ncgen = '/usr/bin/ncgen'
else:
    ncgen = 'ncgen'
cmd = [ncgen, '-o',
       pism_config_nc, pism_config_cdl]
sub.call(cmd)


hydrology = 'diffuse'

# ########################################################
# set up model initialization
# ########################################################


# Model Parameters
ssa_n = (3.0)
ssa_e = (1.0)

# Model Parameters for Sensitivity Studay
sia_e_values = [1.0, 3.0]
ppq_values = [0.50]
tefo_values = [0.020]
plastic_phi_values = [20, 30]
combinations = list(itertools.product(sia_e_values, ppq_values, tefo_values, plastic_phi_values))

tsstep = 'yearly'

scripts = []

for n, combination in enumerate(combinations):

    sia_e, ppq, tefo, plastic_phi = combination

    name_options = OrderedDict()
    name_options['sia_e'] = sia_e
    name_options['plastic_phi'] = plastic_phi
    experiment =  '_'.join([climate, '_'.join(['_'.join([k, str(v)]) for k, v in name_options.items()])])

        
    script = 'cc_{}_g{}m_{}.sh'.format(domain.lower(), grid, experiment)
    scripts.append(script)
    
    for filename in (script):
        try:
            os.remove(filename)
        except OSError:
            pass

    batch_header, batch_system = make_batch_header(system, nn, walltime, queue)
            
    with open(script, 'w') as f:

        f.write(batch_header)

        outfile = '{domain}_g{grid}m_{experiment}_{start}_{end}a.nc'.format(domain=domain.lower(),
                                                                           grid=grid,
                                                                           experiment=experiment,
                                                                           start=int(start),
                                                                            end=int(end))

        prefix = generate_prefix_str(pism_exec)

        # Setup General Parameters
        general_params_dict = OrderedDict()
        if input_file is None:
            general_params_dict['i'] = pism_dataname
            general_params_dict['bootstrap'] = ''
        else:
            general_params_dict['i'] = input_file
        general_params_dict['ys'] = start
        general_params_dict['ye'] = end
        general_params_dict['o'] = os.path.join(odir, outfile)
        general_params_dict['o_format'] = oformat
        general_params_dict['o_size'] = osize
        general_params_dict['config_override'] = pism_config_nc
        
        if input_file is None:
            grid_params_dict = generate_grid_description(grid, accepted_resolutions(), domain)
        else:
            grid_params_dict = generate_grid_description(grid, accepted_resolutions(), domain, restart=True)

        # Setup Stress Balance Paramters
        sb_params_dict = OrderedDict()
        sb_params_dict['sia_e'] = sia_e
        sb_params_dict['ssa_e'] = ssa_e
        sb_params_dict['ssa_n'] = ssa_n
        sb_params_dict['pseudo_plastic_q'] = ppq
        sb_params_dict['till_effective_fraction_overburden'] = tefo
        sb_params_dict['plastic_phi'] = plastic_phi
        sb_params_dict['bed_smoother_range'] = 250.

        stress_balance_params_dict = generate_stress_balance(stress_balance, sb_params_dict)

        # Setup Climate Forcing
        climate_file = 'ltop_climate_olympics_{grid}m_kg_m-2_yr-1.nc'.format(grid=grid)
        climate_params_dict = generate_climate(climate, atmosphere_given_file=climate_file, atmosphere_lapse_rate_file=climate_file)
        # Setup Ocean Forcing
        ocean_params_dict = generate_ocean('null')
        # Setup Hydrology Model
        hydro_params_dict = generate_hydrology(hydrology)

        # Setup Scalar and Spatial Time Series Reporting
        exvars = default_spatial_ts_vars()
        spatial_ts_dict = generate_spatial_ts(outfile, exvars, exstep, odir=odir)
        scalar_ts_dict = generate_scalar_ts(outfile, tsstep, odir=odir)

        # Merge All Parameter Dictionaries
        all_params_dict = merge_dicts(general_params_dict, grid_params_dict, stress_balance_params_dict, climate_params_dict, ocean_params_dict, hydro_params_dict, spatial_ts_dict, scalar_ts_dict)
        all_params = ' '.join([' '.join(['-' + k, str(v)]) for k, v in all_params_dict.items()])
        
        cmd = ' '.join([batch_system['mpido'], prefix, all_params, '> {outdir}/job.${batch}  2>&1'.format(outdir=odir,batch=batch_system['job_id'])])

        f.write(cmd)
        f.write('\n')

    
scripts = uniquify_list(scripts)

print '\n'.join([script for script in scripts])
print('\nwritten\n')
