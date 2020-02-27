"""
resources
=========

Provides:
  - general resources such as grid constructors, calving, hydrology, etc.
    for the Olympics Peninsula

"""

from collections import OrderedDict
import os, math, sys
import os.path


def accepted_resolutions():

    return (50, 100, 200, 250, 500, 1000, 2000, 5000)


def generate_prefix_str(pism_exec):
    """
    Generate prefix string.

    Returns: string
    """

    try:
        p = os.environ["PISM_PREFIX"] + pism_exec
    except:
        p = pism_exec

    return p


def generate_domain(domain):
    """
    Generate domain specific options

    Returns: string
    """

    if domain.lower() in ("olympics", "olympics_mtns"):
        pism_exec = "pismr"
    else:
        print("Domain {} not recognized, exiting".format(domain))
        import sys

        sys.exit(0)

    return pism_exec


def generate_spatial_ts(outfile, exvars, step, start=None, end=None, split=None, odir=None):
    """
    Return dict to generate spatial time series

    Returns: OrderedDict
    """

    # check if list or comma-separated string is given.
    try:
        exvars = ",".join(exvars)
    except:
        pass

    params_dict = OrderedDict()
    if split is True:
        outfile, ext = os.path.splitext(outfile)
        params_dict["extra_split"] = ""
    if odir is None:
        params_dict["extra_file"] = "ex_" + outfile
    else:
        params_dict["extra_file"] = os.path.join(odir, "ex_" + outfile)
    params_dict["extra_vars"] = exvars

    if step is None:
        step = "yearly"

    if start is not None and end is not None:
        times = "{start}:{step}:{end}".format(start=start, step=step, end=end)
    else:
        times = step

    params_dict["extra_times"] = times

    return params_dict


def generate_scalar_ts(outfile, step, start=None, end=None, odir=None):
    """
    Return dict to create scalar time series

    Returns: OrderedDict
    """

    params_dict = OrderedDict()
    if odir is None:
        params_dict["ts_file"] = "ts_" + outfile
    else:
        params_dict["ts_file"] = os.path.join(odir, "ts_" + outfile)

    if step is None:
        step = "yearly"

    if start is not None and end is not None:
        times = "{start}:{step}:{end}".format(start=start, step=step, end=end)
    else:
        times = step
    params_dict["ts_times"] = times

    return params_dict


def generate_snap_shots(outfile, times, odir=None):
    """
    Return dict to generate snap shots

    Returns: OrderedDict
    """

    params_dict = OrderedDict()
    if odir is None:
        params_dict["save_file"] = "save_" + outfile.split(".nc")[0]
    else:
        params_dict["save_file"] = os.path.join(odir, "save_" + outfile.split(".nc")[0])

    params_dict["save_times"] = ",".join(str(e) for e in times)
    params_dict["save_split"] = ""
    params_dict["save_force_output_times"] = ""

    return params_dict


def generate_grid_description(grid_resolution, accepted_resolutions, domain, restart=False):
    """
    Generate grid description dict

    Returns: OrderedDict
    """

    if domain.lower() in ("olympics"):
        mx_max = 4000
        my_max = 3600
    elif domain.lower() in ("olympics_mtns"):
        mx_max = 2400
        my_max = 2000
    else:
        print("domain {} not recongnized".format(domain))
    resolution_max = 50

    try:
        grid_resolution in accepted_resolutions
        pass
    except:
        print("grid resolution {}m not recognized".format(grid_resolution))

    grid_div = grid_resolution / resolution_max

    mx = mx_max / grid_div
    my = my_max / grid_div

    horizontal_grid = OrderedDict()
    horizontal_grid["Mx"] = mx
    horizontal_grid["My"] = my

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
    vertical_grid["Lz"] = 2000
    vertical_grid["Lbz"] = 2000
    vertical_grid["z_spacing"] = "equal"
    vertical_grid["Mz"] = mz
    vertical_grid["Mbz"] = mzb

    grid_options = {}
    grid_options["skip"] = ""
    grid_options["skip_max"] = skip_max

    grid_dict = merge_dicts(horizontal_grid, vertical_grid, grid_options)

    if restart is True:
        return grid_options
    else:
        return grid_dict


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.

    Returns: OrderedDict
    """
    result = OrderedDict()
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def uniquify_list(seq, idfun=None):
    """
    Remove duplicates from a list, order preserving.
    From http://www.peterbe.com/plog/uniqifiers-benchmark
    """

    if idfun is None:

        def idfun(x):
            return x

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
    """
    Generate stress balance params

    Returns: OrderedDict
    """

    accepted_stress_balances = ("sia", "ssa+sia")

    if stress_balance not in accepted_stress_balances:
        print("{} not in {}".format(stress_balance, accepted_stress_balances))
        print("available stress balance solvers are {}".format(accepted_stress_balances))
        import sys

        sys.exit(0)

    params_dict = OrderedDict()
    params_dict["stress_balance"] = stress_balance
    if stress_balance in ("ssa+sia"):
        params_dict["options_left"] = ""
        # params_dict['cfbc'] = ''
        params_dict["sia_flow_law"] = "gpbld"
        params_dict["pseudo_plastic"] = ""

    return merge_dicts(additional_params_dict, params_dict)


def generate_hydrology(hydro, **kwargs):
    """
    Generate hydrology params

    Returns: OrderedDict
    """

    params_dict = OrderedDict()
    if hydro in ("null"):
        params_dict["hydrology"] = "null"
    elif hydro in ("diffuse"):
        params_dict["hydrology"] = "null"
        params_dict["hydrology_null_diffuse_till_water"] = ""
    elif hydro in ("routing"):
        params_dict["hydrology"] = "routing"
    elif hydro in ("distributed"):
        params_dict["hydrology"] = "distributed"
    else:
        print("hydrology {} not recognized, exiting".format(hydro))
        import sys

        sys.exit(0)

    return merge_dicts(params_dict, kwargs)


def generate_calving(calving, **kwargs):
    """
    Generate calving params

    Returns: OrderedDict
    """

    params_dict = OrderedDict()
    if calving in ("ocean_kill"):
        params_dict["calving"] = calving
    elif calving in ("eigen_calving", "vonmises_calving"):
        params_dict["calving"] = "{},thickness_calving".format(calving)
    elif calving in ("hybrid_calving"):
        params_dict["calving"] = "eigen_calving,vonmises_calving,thickness_calving"
    elif calving in ("float_kill"):
        params_dict["calving"] = calving
        params_dict["float_kill_margin_only"] = ""
    else:
        print("calving {} not recognized, exiting".format(calving))
        import sys

        sys.exit(0)

    return merge_dicts(params_dict, kwargs)


def generate_climate(climate, **kwargs):
    """
    Generate climate params

    Returns: OrderedDict
    """

    params_dict = OrderedDict()
    ice_density = 910.0
    if climate in ("elev"):
        params_dict["surface"] = "elevation"
        params_dict["ice_surface_temp"] = "0,0,-100,5000"
        params_dict["climatic_mass_balance"] = "-3.,3,0,800,2500"
    elif climate in ("present"):
        params_dict["atmosphere"] = "yearly_cycle,elevation_change"
        params_dict["surface.pdd.factor_ice"] = 4.59 / ice_density  # Shea et al (2009)
        params_dict["surface.pdd.factor_snow"] = 3.04 / ice_density  # Shea et al (2009)
        params_dict["surface.pdd.refreeze"] = 0
        params_dict["surface"] = "pdd"
    elif climate in ("paleo", "calib", "constant"):
        params_dict["atmosphere"] = "yearly_cycle,elevation_change,delta_T,precip_scaling"
        params_dict["surface.pdd.factor_ice"] = 4.59 / ice_density  # Shea et al (2009)
        params_dict["surface.pdd.factor_snow"] = 3.04 / ice_density  # Shea et al (2009)
        params_dict["surface.pdd.refreeze"] = 0
        params_dict["surface"] = "pdd"
    elif climate in ("constant_orographic"):
        params_dict["atmosphere"] = "yearly_cycle,elevation_change,delta_T,precip_scalin,orographic_precipitation"
        params_dict["atmosphere.orographic_precipitation.background_precip_post"] = 0.057
        params_dict["atmosphere.orographic_precipitation.background_precip_pre"] = 1.0
        params_dict["atmosphere.orographic_precipitation.conversion_time"] = 1750.0
        params_dict["atmosphere.orographic_precipitation.coriolis_latitude"] = 0.0
        params_dict["atmosphere.orographic_precipitation.fallout_time"] = 1750.0
        params_dict["atmosphere.orographic_precipitation.lapse_rate"] = -5.8
        params_dict["atmosphere.orographic_precipitation.moist_adiabatic_lapse_rate"] = -6.5
        params_dict["atmosphere.orographic_precipitation.moist_stability_frequency"] = 0.007
        params_dict["atmosphere.orographic_precipitation.reference_density"] = 7.4e-3
        params_dict["atmosphere.orographic_precipitation.scale_factor"] = 0.125
        params_dict["atmosphere.orographic_precipitation.truncate"] = "true"
        params_dict["atmosphere.orographic_precipitation.water_vapor_scale_height"] = 2600.0
        params_dict["atmosphere.orographic_precipitation.wind_direction"] = 220
        params_dict["atmosphere.orographic_precipitation.wind_speed"] = 15
        params_dict["surface.pdd.factor_ice"] = 4.59 / ice_density  # Shea et al (2009)
        params_dict["surface.pdd.factor_snow"] = 3.04 / ice_density  # Shea et al (2009)
        params_dict["surface.pdd.refreeze"] = 0
        params_dict["surface"] = "pdd"
    else:
        print("climate {} not recognized, exiting".format(climate))
        import sys

        sys.exit(0)

    return merge_dicts(params_dict, kwargs)


# ../../scripts/linear_orog_precip.py -i ../bed_dem/pism_Olympics_${GRID}m_v1.nc --background_precip_pre 1 --background_precip_post 0.057 --precip_scale_factor 0.125 --tau_c 1750 --tau_f 1750 --wind_direction $dir --wind_magnitude 15 --moist_stability 0.007 --vapor_scale_height 2600


def generate_ocean(ocean, **kwargs):
    """
    Generate ocean params

    Returns: OrderedDict
    """

    params_dict = OrderedDict()
    if ocean in ("null"):
        pass
    elif ocean in ("const"):
        params_dict["ocean"] = "constant"
    else:
        print("ocean {} not recognized, exiting".format(ocean))
        import sys

        sys.exit(0)

    return merge_dicts(params_dict, kwargs)


spatial_ts_vars = {}


spatial_ts_vars["basic"] = [
    "beta",
    "dHdt",
    "diffusivity",
    "effective_precipitation",
    "ice_mass",
    "mask",
    "mass_fluxes",
    "sftgif",
    "thk",
    "topg",
    "usurf",
    "velbase_mag",
    "velsurf_mag",
]


def list_systems():

    """
    Return a list of supported systems.
    """

    list = [
        "debug",
        "chinook",
        "electra_broadwell",
        "fish",
        "keeling",
        "pacman",
        "pleiades",
        "pleiades_haswell",
        "pleiades_ivy",
        "pleiades_broadwell",
    ]

    return list


def list_queues():

    """
    Return a list of supported queues.
    """

    list = ["debug", "devel", "gpu", "gpu_long", "normal", "long", "standard", "standard_16", "t2standard", "t2small"]

    return list


# information about systems
systems = {}

systems["debug"] = {"mpido": "mpiexec -n {cores}", "submit": "echo", "job_id": "PBS_JOBID", "queue": {}}

systems["chinook"] = {
    "mpido": "mpirun -np {cores} -machinefile ./nodes_$SLURM_JOBID",
    "submit": "sbatch",
    "work_dir": "SLURM_SUBMIT_DIR",
    "job_id": "SLURM_JOBID",
    "queue": {"t1standard": 24, "t1small": 24, "t2standard": 24, "t2small": 24, "debug": 24, "analysis": 24},
}

systems["pleiades"] = {
    "mpido": "mpiexec -n {cores}",
    "submit": "qsub",
    "work_dir": "PBS_O_WORKDIR",
    "job_id": "PBS_JOBID",
    "queue": {"long": 20, "normal": 20},
}

systems["pleiades_haswell"] = systems["pleiades"].copy()
systems["pleiades_haswell"]["queue"] = {"long": 24, "normal": 24}

systems["pleiades_ivy"] = systems["pleiades"].copy()
systems["pleiades_ivy"]["queue"] = {"long": 20, "normal": 20}

systems["pleiades_sandy"] = systems["pleiades"].copy()
systems["pleiades_sandy"]["queue"] = {"long": 16, "normal": 16}

systems["pleiades_broadwell"] = systems["pleiades"].copy()
systems["pleiades_broadwell"]["queue"] = {"long": 28, "normal": 28}

systems["electra_broadwell"] = systems["pleiades_broadwell"].copy()

systems["electra_skylake"] = systems["pleiades"].copy()
systems["electra_skylake"]["queue"] = {"long": 40, "normal": 40}


# headers for batch jobs
#
# Available keywords:
#
# cores    - number of cores (MPI tasks)
# queue    - queue (partition) name
# nodes    - number of nodes
# ppn      - number of tasks per node
# walltime - wall time limit

systems["debug"]["header"] = ""

systems["chinook"][
    "header"
] = """#!/bin/sh
#SBATCH --partition={queue}
#SBATCH --ntasks={cores}
#SBATCH --tasks-per-node={ppn}
#SBATCH --time={walltime}
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=pism.%j

module list

umask 007

cd $SLURM_SUBMIT_DIR

# Generate a list of compute node hostnames reserved for this job,
# this ./nodes file is necessary for slurm to spawn mpi processes
# across multiple compute nodes
srun -l /bin/hostname | sort -n | awk '{{print $2}}' > ./nodes_$SLURM_JOBID

ulimit -l unlimited
ulimit -s unlimited
ulimit

"""

systems["chinook"][
    "footer"
] = """
# clean up the list of hostnames
rm -rf ./nodes_$SLURM_JOBID
"""

systems["electra_broadwell"][
    "header"
] = """#PBS -S /bin/bash
#PBS -N cfd
#PBS -l walltime={walltime}
#PBS -m e
#PBS -q {queue}
#PBS -lselect={nodes}:ncpus={ppn}:mpiprocs={ppn}:model=bro_ele
#PBS -j oe

module list

cd $PBS_O_WORKDIR

"""

systems["pleiades"][
    "header"
] = """#PBS -S /bin/bash
#PBS -N cfd
#PBS -l walltime={walltime}
#PBS -m e
#PBS -W group_list=s1878
#PBS -q {queue}
#PBS -lselect={nodes}:ncpus={ppn}:mpiprocs={ppn}:model=ivy
#PBS -j oe

module list

cd $PBS_O_WORKDIR

"""

systems["pleiades_broadwell"][
    "header"
] = """#PBS -S /bin/bash
#PBS -N cfd
#PBS -l walltime={walltime}
#PBS -m e
#PBS -W group_list=s1878
#PBS -q {queue}
#PBS -lselect={nodes}:ncpus={ppn}:mpiprocs={ppn}:model=bro
#PBS -j oe

module list

cd $PBS_O_WORKDIR

"""

systems["pleiades_sandy"][
    "header"
] = """#PBS -S /bin/bash
#PBS -N cfd
#PBS -l walltime={walltime}
#PBS -m e
#PBS -W group_list=s1878
#PBS -q {queue}
#PBS -lselect={nodes}:ncpus={ppn}:mpiprocs={ppn}:model=san
#PBS -j oe

module list

cd $PBS_O_WORKDIR

"""

systems["pleiades_haswell"][
    "header"
] = """#PBS -S /bin/bash
#PBS -N cfd
#PBS -l walltime={walltime}
#PBS -m e
#PBS -W group_list=s1878
#PBS -q {queue}
#PBS -lselect={nodes}:ncpus={ppn}:mpiprocs={ppn}:model=has
#PBS -j oe

module list

cd $PBS_O_WORKDIR

"""

systems["pleiades_ivy"][
    "header"
] = """#PBS -S /bin/bash
#PBS -N cfd
#PBS -l walltime={walltime}
#PBS -m e
#PBS -W group_list=s1878
#PBS -q {queue}
#PBS -lselect={nodes}:ncpus={ppn}:mpiprocs={ppn}:model=ivy
#PBS -j oe

module list

cd $PBS_O_WORKDIR

"""

systems["electra_skylake"][
    "header"
] = """#PBS -S /bin/bash
#PBS -N cfd
#PBS -l walltime={walltime}
#PBS -m e
#PBS -W group_list=s1878
#PBS -q {queue}
#PBS -lselect={nodes}:ncpus={ppn}:mpiprocs={ppn}:model=sky_ele
#PBS -j oe

module list

cd $PBS_O_WORKDIR

"""

systems["debug"][
    "header"
] = """

"""

# headers for post-processing jobs

post_headers = {}
post_headers[
    "default"
] = """#!/bin/bash

"""

post_headers[
    "pbs"
] = """#PBS -S /bin/bash
#PBS -l select=1:mem=94GB
#PBS -l walltime=8:00:00
#PBS -q ldan

cd $PBS_O_WORKDIR

"""

post_headers[
    "slurm"
] = """#!/bin/bash
#SBATCH --partition=analysis
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --output=pism.%j
#SBATCH --mem=214G

cd $SLURM_SUBMIT_DIR

ulimit -l unlimited
ulimit -s unlimited
ulimit

"""


def make_batch_header(system_name, n_cores, walltime, queue):
    """
    Generate header file for different HPC system.

    Returns: String
    """

    # get system info; use "debug" if the requested name was not found
    system = systems.get(system_name, systems["debug"]).copy()

    assert n_cores > 0

    if system_name == "debug":
        # when debugging, assume that all we need is one node
        ppn = n_cores
        nodes = 1
    else:
        try:
            ppn = system["queue"][queue]
        except:
            raise ValueError(
                "There is no queue {} on {}. Pick one of {}.".format(queue, system_name, list(system["queue"].keys()))
            )
        # round up when computing the number of nodes needed to run on 'n_cores' cores
        nodes = int(math.ceil(float(n_cores) / ppn))

        if nodes * ppn != n_cores:
            print(
                (
                    "Warning! Running {n_cores} tasks on {nodes} {ppn}-processor nodes, wasting {N} processors!".format(
                        nodes=nodes, ppn=ppn, n_cores=n_cores, N=ppn * nodes - n_cores
                    )
                )
            )

    system["mpido"] = system["mpido"].format(cores=n_cores)
    system["header"] = system["header"].format(queue=queue, walltime=walltime, nodes=nodes, ppn=ppn, cores=n_cores)
    system["header"] += version_header()

    return system["header"], system


def make_batch_post_header(system):

    v = version_header()

    if system in ("electra_broadwell", "pleiades", "pleiades_ivy", "pleiades_broadwell", "pleiades_haswell"):
        return post_headers["pbs"] + v
    elif system in ("chinook"):
        return post_headers["slurm"] + v
    else:
        return post_headers["default"] + v


def make_batch_header_test():
    "print headers of all supported systems and queues (for testing)"
    for s in list(systems.keys()):
        for q in list(systems[s]["queue"].keys()):
            print("# system: {system}, queue: {queue}".format(system=s, queue=q))
            print(make_batch_header(s, 100, "1:00:00", q)[0])


def version():
    """Return the path to the top directory of the Git repository
    containing this script, the URL of the "origin" remote and the version."""
    import inspect
    import shlex
    import subprocess

    def output(command):
        path = os.path.realpath(os.path.dirname(inspect.stack(0)[0][1]))
        return subprocess.check_output(shlex.split(command), cwd=path).strip()

    return (
        output("git rev-parse --show-toplevel"),
        output("git remote get-url origin"),
        output("git describe --always"),
    )


def version_header():
    "Return shell comments containing version info."
    version_info = version()
    return """
# Generated by {script}
# Command: {command}
# Git top level: {path}
# URL: {url}
# Version: {version}

""".format(
        script=os.path.realpath(sys.argv[0]),
        command=" ".join(sys.argv),
        path=version_info[0],
        url=version_info[1],
        version=version_info[2],
    )
