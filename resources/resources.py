"""
resources
=========

Provides:
  - general resources such as grid constructors, calving, hydrology, etc.
    for the Olympic Peninsula

"""

from collections import OrderedDict
import os


def accepted_resolutions():

    return (50, 100, 200, 250, 500, 1000, 2000, 5000)


def generate_prefix_str(pism_exec):
    '''
    Generate prefix string.

    Returns: string
    '''

    try:
        p = os.environ['PISM_PREFIX']  + pism_exec
    except:
        p  = pism_exec
    
    return p


def generate_domain(domain):
    '''
    Generate domain specific options

    Returns: string
    '''
    
    if domain.lower() in ('olympics', 'olympics_mtns'):
        pism_exec = 'pismr'
    else:
        print('Domain {} not recognized, exiting'.format(domain))
        import sys
        sys.exit(0)

    return pism_exec


def default_spatial_ts_vars():
    '''
    Returns a list of commonly-used extra vars
    '''
    
    exvars = ['air_temp_snapshot',
              'beta',
              'bmelt',
              'cell_area',
              'climatic_mass_balance',
              'dHdt',
              'diffusivity',
              'effective_air_temp',
              'effective_precipitation',
              'ice_mass',
              'ice_surface_temp',
              'mask',
              'lat',
              'lat_bnds',
              'lon',
              'lon_bnds',
              'precipitation',
              'saccum',
              'smelt',
              'srunoff',
              'surface_mass_balance_average',
              'taub_mag',
              'tauc',
              'taud_mag',
              'tempicethk_basal',
              'temppabase',
              'tempsurf',
              'thk',
              'topg',
              'usurf',
              'velbase_mag',
              'velsurf_mag']
    
    return exvars



def generate_spatial_ts(outfile, exvars, step, start=None, end=None, split=None, odir=None):
    '''
    Return dict to generate spatial time series

    Returns: OrderedDict
    '''

    # check if list or comma-separated string is given.
    try:
        exvars = ','.join(exvars)
    except:
        pass

    params_dict = OrderedDict()
    if split is True:
        outfile, ext = os.path.splitext(outfile)
        params_dict['extra_split'] = ''
    if odir is None:
        params_dict['extra_file'] = 'ex_' + outfile
    else:
        params_dict['extra_file'] = os.path.join(odir, 'ex_' + outfile)
    params_dict['extra_vars'] = exvars
        
    if step is None:
        step = 'yearly'

    if (start is not None and end is not None):
        times = '{start}:{step}:{end}'.format(start=start, step=step, end=end)
    else:
        times = step
        
    params_dict['extra_times'] = times
        
  
    return params_dict


def generate_scalar_ts(outfile, step, start=None, end=None, odir=None):
    '''
    Return dict to create scalar time series

    Returns: OrderedDict
    '''

    params_dict = OrderedDict()
    if odir is None:
        params_dict['ts_file'] = 'ts_' + outfile
    else:
        params_dict['ts_file'] = os.path.join(odir, 'ts_' + outfile)
    
    if step is None:
        step = 'yearly'

    if (start is not None and end is not None):
        times = '{start}:{step}:{end}'.format(start=start, step=step, end=end)
    else:
        times = step
    params_dict['ts_times'] = times

    return params_dict


def generate_snap_shots(outfile, times, odir=None):
    '''
    Return dict to generate snap shots

    Returns: OrderedDict
    '''

    params_dict = OrderedDict()
    if odir is None:
        params_dict['save_file'] = 'save_' + outfile.split('.nc')[0]
    else:
        params_dict['save_file'] = os.path.join(odir, 'save_' + outfile.split('.nc')[0])

    params_dict['save_times'] = ','.join(str(e) for e in times)
    params_dict['save_split'] = ''
    params_dict['save_force_output_times'] = ''

    return params_dict


def generate_grid_description(grid_resolution, accepted_resolutions, domain, restart=False):
    '''
    Generate grid description dict

    Returns: OrderedDict
    '''

    if domain.lower() in ('olympics'):
        mx_max = 4000
        my_max = 3600
    elif domain.lower() in ('olympics_mtns'):
        mx_max = 2400
        my_max = 2000
    else:
        print('domain {} not recongnized'.format(domain))
    resolution_max = 50
    

    try:
        grid_resolution in accepted_resolutions
        pass
    except:
        print('grid resolution {}m not recognized'.format(grid_resolution))

    grid_div = (grid_resolution / resolution_max)
              
    mx = mx_max / grid_div
    my = my_max / grid_div

    horizontal_grid = OrderedDict()
    horizontal_grid['Mx'] = mx
    horizontal_grid['My'] = my

    if grid_resolution < 200:
        skip_max = 500
        mz = 101
        mzb = 21
    elif (grid_resolution >= 200) and (grid_resolution <= 500):
        skip_max = 250
        mz = 51
        mzb = 11
    else:
        skip_max = 100
        mz = 26
        mzb = 6

    vertical_grid = OrderedDict()
    vertical_grid['Lz'] = 2000
    vertical_grid['Lbz'] = 2000
    vertical_grid['z_spacing'] = 'equal'
    vertical_grid['Mz'] = mz
    vertical_grid['Mbz'] = mzb

    grid_options = {}
    grid_options['skip'] = ''
    grid_options['skip_max'] = skip_max

    grid_dict = merge_dicts(horizontal_grid, vertical_grid, grid_options)

    if restart is True:
        return grid_options
    else:
        return grid_dict


def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.

    Returns: OrderedDict
    '''
    result = OrderedDict()
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def uniquify_list(seq, idfun=None):
    '''
    Remove duplicates from a list, order preserving.
    From http://www.peterbe.com/plog/uniqifiers-benchmark
    '''

    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen:
            continue
        seen[marker] = 1
        result.append(item)
    return result


def generate_stress_balance(stress_balance, additional_params_dict):
    '''
    Generate stress balance params

    Returns: OrderedDict
    '''

    accepted_stress_balances = ('sia', 'ssa+sia')

    if stress_balance not in accepted_stress_balances:
        print('{} not in {}'.format(stress_balance, accepted_stress_balances))
        print('available stress balance solvers are {}'.format(accepted_stress_balances))
        import sys
        sys.exit(0)

    params_dict = OrderedDict()
    params_dict['stress_balance'] = stress_balance
    if stress_balance in ('ssa+sia'):
        params_dict['options_left'] = ''
        # params_dict['ssafd_pc_type'] = 'asm'
        # params_dict['ssafd_sub_pc_type'] = 'jacobi'
        params_dict['cfbc'] = ''
        params_dict['sia_flow_law'] = 'gpbld3'
        params_dict['pseudo_plastic'] = ''
        params_dict['tauc_slippery_grounding_lines'] = ''

    return merge_dicts(additional_params_dict, params_dict)


def generate_hydrology(hydro, **kwargs):
    '''
    Generate hydrology params

    Returns: OrderedDict
    '''
    
    params_dict = OrderedDict()
    if hydro in ('null'):
        params_dict['hydrology'] = 'null'
    elif hydro in ('diffuse'):
        params_dict['hydrology'] = 'null'
        params_dict['hydrology_null_diffuse_till_water'] = ''
    elif hydro in ('routing'):
        params_dict['hydrology'] = 'routing'
    elif hydro in ('distributed'):
        params_dict['hydrology'] = 'distributed'
    else:
        print('hydrology {} not recognized, exiting'.format(hydro))
        import sys
        sys.exit(0)

    return merge_dicts(params_dict, kwargs)


def generate_calving(calving, **kwargs):
    '''
    Generate calving params

    Returns: OrderedDict
    '''

    params_dict = OrderedDict()
    if calving in ('ocean_kill'):
        params_dict['calving'] = calving
    elif calving in ('eigen_calving', 'vonmises_calving'):
        params_dict['calving'] = '{},thickness_calving'.format(calving)
    elif calving in ('hybrid_calving'):
        params_dict['calving'] = 'eigen_calving,vonmises_calving,thickness_calving'
    elif calving in ('float_kill'):
        params_dict['calving'] = calving
        params_dict['float_kill_margin_only'] = ''
    else:
        print('calving {} not recognized, exiting'.format(calving))
        import sys
        sys.exit(0)

    return merge_dicts(params_dict, kwargs)


def generate_climate(climate, **kwargs):
    '''
    Generate climate params

    Returns: OrderedDict
    '''
    
    params_dict = OrderedDict()
    if climate in ('elev'):
        params_dict['surface'] = 'elevation'
        params_dict['ice_surface_temp'] = '2,-15,0,2000'
        params_dict['climatic_mass_balance'] = '-3.,3,0,800,2500'
    elif climate in ('present'):
        params_dict['atmosphere'] = 'yearly_cycle,lapse_rate'
        if 'atmosphere_yearly_cycle_file' not in kwargs:
            params_dict['atmosphere_yearly_cycle_file'] = 'olympics_climate_1000m.nc'
        else:
            params_dict['atmosphere_yearly_cycle_file'] = kwargs['atmosphere_yearly_cycle_file']
        params_dict['temp_lapse_rate'] = 4.5
        if 'atmosphere_lapse_rate_file' not in kwargs:
            params_dict['atmosphere_lapse_rate_file'] = 'olympics_climate_1000m.nc'
        else:
            params_dict['atmosphere_lapse_rate_file'] = kwargs['atmosphere_lapse_rate_file']
        params_dict['surface'] = 'pdd'
    elif climate in ('paleo'):
        params_dict['atmosphere'] = 'yearly_cycle,lapse_rate,delta_T,frac_P'
        if 'atmosphere_yearly_cycle_file' not in kwargs:
            params_dict['atmosphere_yearly_cycle_file'] = 'olympics_climate_1000m.nc'
        else:
            params_dict['atmosphere_yearly_cycle_file'] = kwargs['atmosphere_yearly_cycle_file']
        params_dict['temp_lapse_rate'] = 4.5
        if 'atmosphere_lapse_rate_file' not in kwargs:
            params_dict['atmosphere_lapse_rate_file'] = 'olympics_climate_1000m.nc'
        else:
            params_dict['atmosphere_lapse_rate_file'] = kwargs['atmosphere_lapse_rate_file']
        if 'atmosphere_delta_T_file' not in kwargs:
            params_dict['atmosphere_delta_T_file'] = 'paleo_modifier.nc'
        else:
            params_dict['atmosphere_delta_T_file'] = kwargs['atmosphere_delta_T_file']
        if 'atmosphere_frac_P_file' not in kwargs:
            params_dict['atmosphere_frac_P_file'] = 'paleo_modifier.nc'
        else:
            params_dict['atmosphere_delta_T_file'] = kwargs['atmosphere_frac_P_file']
        params_dict['surface'] = 'pdd'
    else:
        print('climate {} not recognized, exiting'.format(climate))
        import sys
        sys.exit(0)
        
    return params_dict

        
def generate_ocean(ocean, **kwargs):
    '''
    Generate ocean params

    Returns: OrderedDict
    '''

    params_dict = OrderedDict()
    if ocean in ('null'):
        pass
    elif ocean in ('const'):
        params_dict['ocean'] = 'constant'
    else:
        print('ocean {} not recognized, exiting'.format(ocean))
        import sys
        sys.exit(0)

    return merge_dicts(params_dict, kwargs)


def list_systems():

    '''
    Return a list of supported systems.
    '''
    
    list = ['debug',
            'chinook',
            'fish',
            'pacman',
            'pleiades',
            'pleiades_ivy',
            'pleiades_broadwell']
    
    return list


def list_queues():

    '''
    Return a list of supported queues.
    '''
    
    list = ['debug',
            'devel',
            'gpu',
            'gpu_long',
            'normal',
            'long',
            'standard',
            'standard_16',
            't1standard',
            't2standard',
            't1small']
    
    return list


def make_batch_header(system, cores, walltime, queue):
    '''
    Generate header file for different HPC system.

    Returns: String
    '''
    
    systems = {}
    mpido = 'mpiexec -n {cores}'.format(cores=cores)
    systems['debug'] = {'mpido' : mpido,
                        'submit': 'echo',
                        'job_id' : 'PBS_JOBID'}
    mpido = 'mpiexec -n {cores}'.format(cores=cores)
    systems['fish'] = {'mpido': 'aprun -n {cores}'.format(cores=cores),
                       'submit' : 'qsub',
                       'work_dir' : 'PBS_O_WORKDIR',
                       'job_id' : 'PBS_JOBID',
                       'queue' : {
                           'gpu' : 16,
                           'gpu_long' : 16,
                           'standard' : 12 }}
    mpido = 'mpirun -np {cores}'.format(cores=cores)
    systems['pacman'] = {'mpido' : mpido,
                         'submit' : 'qsub',
                         'work_dir' : 'PBS_O_WORKDIR',
                         'job_id' : 'PBS_JOBID',
                         'queue' : {
                             'standard_16' : 16 }}
    mpido = 'mpirun -np {cores} -machinefile ./nodes_$SLURM_JOBID'.format(cores=cores)                         
    systems['chinook'] = {'mpido' : mpido,
                          'submit' : 'sbatch',
                          'work_dir' : 'SLURM_SUBMIT_DIR',
                          'job_id' : 'SLURM_JOBID',
                          'queue' : {
                              't1standard' : 24,
                              't2standard' : 24,
                              't1small' : 24,
                              'debug' : 24}}
    mpido = 'mpiexec.hydra -n {cores}'.format(cores=cores)
    systems['pleiades'] = {'mpido' : mpido,
                           'submit' : 'qsub',
                           'work_dir' : 'PBS_O_WORKDIR',
                           'job_id' : 'PBS_JOBID',
                           'queue' : {
                               'long' : 20,
                               'normal': 20}}
    systems['pleiades_ivy'] = {'mpido' : mpido,
                           'submit' : 'qsub',
                           'work_dir' : 'PBS_O_WORKDIR',
                           'job_id' : 'PBS_JOBID',
                           'queue' : {
                               'long' : 20,
                               'normal': 20}}
    systems['pleiades_broadwell'] = {'mpido' : mpido,
                           'submit' : 'qsub',
                           'work_dir' : 'PBS_O_WORKDIR',
                           'job_id' : 'PBS_JOBID',
                           'queue' : {
                               'long' : 28,
                               'normal': 28}}

    assert system in systems.keys()
    if system not in 'debug':
        assert queue in systems[system]['queue'].keys()
        assert cores > 0

        ppn = systems[system]['queue'][queue]
        nodes = cores / ppn

    if system in ('debug'):

        header = ''
        
    elif system in ('chinook'):
        
        header = """#!/bin/sh
#SBATCH --partition={queue}
#SBATCH --ntasks={cores}
#SBATCH --tasks-per-node={ppn}
#SBATCH --time={walltime}
#SBATCH --mail-user=aaschwanden@alaska.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=pism.%j

module list

cd $SLURM_SUBMIT_DIR

# Generate a list of compute node hostnames reserved for this job,
# this ./nodes file is necessary for slurm to spawn mpi processes
# across multiple compute nodes
srun -l /bin/hostname | sort -n | awk \'{{print $2}}\' > ./nodes_$SLURM_JOBID

""".format(queue=queue, walltime=walltime, nodes=nodes, ppn=ppn, cores=cores)
    elif system in ('pleiades'):
        
        header = """#PBS -S /bin/bash
#PBS -N cfd
#PBS -l walltime={walltime}
#PBS -m e
#PBS -q {queue}
#PBS -lselect={nodes}:ncpus={ppn}:mpiprocs={ppn}:model=ivy
#PBS -j oe

module list

cd $PBS_O_WORKDIR

""".format(queue=queue, walltime=walltime, nodes=nodes, ppn=ppn, cores=cores)
    elif system in ('pleiades_broadwell'):
        
        header = """#PBS -S /bin/bash
#PBS -N cfd
#PBS -l walltime={walltime}
#PBS -m e
#PBS -q {queue}
#PBS -lselect={nodes}:ncpus={ppn}:mpiprocs={ppn}:model=bro
#PBS -j oe

module list

cd $PBS_O_WORKDIR

""".format(queue=queue, walltime=walltime, nodes=nodes, ppn=ppn, cores=cores)
    else:
        header = """#!/bin/bash
#PBS -q {queue}
#PBS -l walltime={walltime}
#PBS -l nodes={nodes}:ppn={ppn}
#PBS -j oe

module list

cd $PBS_O_WORKDIR

""".format(queue=queue, walltime=walltime, nodes=nodes, ppn=ppn, cores=cores)

    return header, systems[system]


def make_batch_post_header(system):

    if system in ('pleiades', 'pleiades_ivy', 'pleiades_broadwell', 'pleiades_haswell'):

        header = """#PBS -S /bin/bash
#PBS -lselect=1:mem=94GB
#PBS -lwalltime=8:00:00
#PBS -q ldan

module list

cd $PBS_O_WORKDIR

"""
    elif system in ('chinook'):
        header = """#!/bin/bash
#PBS -q t2small
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe

module list

cd $PBS_O_WORKDIR

"""
    else:
        header = """#!/bin/bash

"""
    return header
