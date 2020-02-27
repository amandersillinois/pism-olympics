#!/usr/bin/env python
# Copyright (C) 2016-20 Andy Aschwanden

import itertools
from collections import OrderedDict
import numpy as np
import os
import sys
import shlex
from os.path import join, abspath, realpath, dirname

try:
    import subprocess32 as sub
except:
    import subprocess as sub

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys


def current_script_directory():
    import inspect

    filename = inspect.stack(0)[0][1]
    return realpath(dirname(filename))


script_directory = current_script_directory()

sys.path.append(join(script_directory, "../resources"))
from resources import *


def map_dict(val, mdict):
    try:
        return mdict[val]
    except:
        return val


# set up the option parser
parser = ArgumentParser()
parser.description = "Generating scripts for prognostic simulations."
parser.add_argument(
    "-n", "--n_procs", dest="n", type=int, help="""number of cores/processors. default=2.""", default=2
)
parser.add_argument("-w", "--wall_time", dest="walltime", help="""walltime. default: 12:00:00.""", default="12:00:00")
parser.add_argument(
    "-q", "--queue", dest="queue", choices=list_queues(), help="""queue. default=t1standard.""", default="normal"
)
parser.add_argument(
    "--climate",
    dest="climate",
    choices=["elev", "paleo", "present", "calib", "constant"],
    help="Climate",
    default="constant",
)
parser.add_argument(
    "-d",
    "--domain",
    dest="domain",
    choices=["olympics", "olympics_mtns"],
    help="sets the modeling domain",
    default="olympics",
)
parser.add_argument("--start_year", dest="start", type=float, help="Start year", default=0)
parser.add_argument("--duration", dest="duration", type=float, help="Duration", default=1000)
parser.add_argument("--exstep", dest="exstep", type=float, help="Spatial time series writing interval", default=100)
parser.add_argument(
    "-f",
    "--o_format",
    dest="oformat",
    choices=["netcdf3", "netcdf4_parallel", "pnetcdf"],
    help="output format",
    default="netcdf3",
)
parser.add_argument(
    "-g",
    "--grid",
    dest="grid",
    type=int,
    choices=accepted_resolutions(),
    help="horizontal grid resolution",
    default=1000,
)
parser.add_argument("-i", "--input_file", dest="input_file", help="Input file to restart from", default=None)
parser.add_argument("--i_dir", dest="input_dir", help="input directory", default=abspath(join(script_directory, "..")))
parser.add_argument("--o_dir", dest="output_dir", help="output directory", default="test_dir")
parser.add_argument(
    "--o_size", dest="osize", choices=["small", "medium", "big", "big_2d"], help="output size type", default="medium"
)
parser.add_argument(
    "-s", "--system", dest="system", choices=list_systems(), help="computer system to use.", default="debug"
)
parser.add_argument("--bed_version", dest="bed_version", help="Version of bed DEM.", default="1")
parser.add_argument(
    "--stress_balance",
    dest="stress_balance",
    choices=["sia", "ssa+sia", "ssa"],
    help="stress balance solver",
    default="ssa+sia",
)
parser.add_argument(
    "-e",
    "--ensemble_file",
    dest="ensemble_file",
    help="File that has all combinations for ensemble study",
    default="../uncertainty_quantification/constant_climate.csv",
)
parser.add_argument(
    "--spatial_ts", dest="spatial_ts", choices=["basic", "standard", "none"], help="output size type", default="basic",
)


options = parser.parse_args()

nn = options.n
input_dir = abspath(options.input_dir)
output_dir = abspath(options.output_dir)
spatial_tmp_dir = abspath(options.output_dir + "_tmp")

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
end = start + options.duration
exstep = options.exstep
spatial_ts = options.spatial_ts

domain = options.domain
pism_exec = generate_domain(domain)
input_file = options.input_file
ensemble_file = options.ensemble_file
pism_dataname = pism_dataname = "$input_dir/data_sets/bed_dem/pism_{domain}_{grid}m_v{version}.nc".format(
    domain=domain.capitalize(), grid=grid, version=bed_version
)


dirs = {"output": "$output_dir", "spatial_tmp": "$spatial_tmp_dir"}
for d in ["performance", "state", "scalar", "spatial", "jobs", "basins"]:
    dirs[d] = "$output_dir/{dir}".format(dir=d)

if spatial_ts == "none":
    del dirs["spatial"]

# use the actual path of the run scripts directory (we need it now and
# not during the simulation)
scripts_dir = join(output_dir, "run_scripts")
if not os.path.isdir(scripts_dir):
    os.makedirs(scripts_dir)

# use the actual path of the run scripts directory (we need it now and
# not during the simulation)
time_dir = join(output_dir, "time_forcing")
if not os.path.isdir(time_dir):
    os.makedirs(time_dir)

# generate the config file *after* creating the output directory
pism_config = "olympics_config"
pism_config_nc = join(output_dir, pism_config + ".nc")

cmd = "ncgen -o {output} {input_dir}/config/{config}.cdl".format(
    output=pism_config_nc, input_dir=input_dir, config=pism_config
)
sub.call(shlex.split(cmd))

# these Bash commands are added to the beginning of the run scrips
run_header = """# stop if a variable is not defined
set -u
# stop on errors
set -e

# path to the config file
config="{config}"
# path to the input directory (input data sets are contained in this directory)
input_dir="{input_dir}"
# output directory
output_dir="{output_dir}"
# temporary directory for spatial files
spatial_tmp_dir="{spatial_tmp_dir}"

# create required output directories
for each in {dirs};
do
  mkdir -p $each
done

""".format(
    input_dir=input_dir,
    output_dir=output_dir,
    spatial_tmp_dir=spatial_tmp_dir,
    config=pism_config_nc,
    dirs=" ".join(list(dirs.values())),
)

hydrology = "diffuse"


# ########################################################
# set up model initialization
# ########################################################


# Model Parameters
ssa_n = 3.0
ssa_e = 1.0
wind_direction = 220

try:
    combinations = np.loadtxt(ensemble_file, delimiter=",", skiprows=1)
except:
    combinations = np.genfromtxt(ensemble_file, dtype=None, delimiter=",", skip_header=1)

tsstep = "yearly"

scripts = []
scripts_post = []

for n, combination in enumerate(combinations):

    (
        run_id,
        precip_scale_factor,
        temperature_offset,
        temperature_lapse_rate,
        sia_e,
        ppq,
        tefo,
        phi_min,
        phi_max,
        topg_min,
        topg_max,
        bed_smoother_range,
    ) = combination

    ttphi = "{},{},{},{}".format(phi_min, phi_max, topg_min, topg_max)

    name_options = {}
    name_options["id"] = "{}".format(int(run_id))

    experiment = "_".join([climate, "_".join(["_".join([k, str(v)]) for k, v in name_options.items()])])

    atmosphere_paleo_file = "paleo_modifier_{}K.nc".format(temperature_offset)

    script = join(scripts_dir, "cc_{}_g{}m_{}.sh".format(domain.lower(), grid, experiment))
    scripts.append(script)
    script_post = join(scripts_dir, "cc_{}_g{}m_{}_post.sh".format(domain.lower(), grid, experiment))
    scripts_post.append(script_post)

    for filename in script:
        try:
            os.remove(filename)
        except OSError:
            pass

    batch_header, batch_system = make_batch_header(system, nn, walltime, queue)

    with open(script, "w") as f:

        f.write(batch_header)
        f.write(run_header)

        outfile = "{domain}_g{grid}m_{experiment}_{start}_{end}a.nc".format(
            domain=domain.lower(), grid=grid, experiment=experiment, start=int(start), end=int(end)
        )

        pism = generate_prefix_str(pism_exec)

        # Setup General Parameters
        general_params_dict = OrderedDict()
        if input_file is None:
            general_params_dict["i"] = pism_dataname
            general_params_dict["bootstrap"] = ""
        else:
            general_params_dict["i"] = input_file
        general_params_dict["ys"] = start
        general_params_dict["ye"] = end
        general_params_dict["o"] = join(dirs["state"], outfile)

        general_params_dict["o_format"] = oformat
        general_params_dict["o_size"] = osize
        general_params_dict["config_override"] = pism_config_nc

        if input_file is None:
            grid_params_dict = generate_grid_description(grid, accepted_resolutions(), domain)
        else:
            grid_params_dict = generate_grid_description(grid, accepted_resolutions(), domain, restart=True)

        # Setup Stress Balance Paramters
        sb_params_dict = OrderedDict()
        sb_params_dict["sia_e"] = sia_e
        sb_params_dict["ssa_e"] = ssa_e
        sb_params_dict["ssa_n"] = ssa_n
        sb_params_dict["pseudo_plastic_q"] = ppq
        sb_params_dict["till_effective_fraction_overburden"] = tefo
        sb_params_dict["topg_to_phi"] = ttphi
        sb_params_dict["ssa_method"] = "fd"
        sb_params_dict["stress_balance.sia.bed_smoother.range"] = bed_smoother_range

        stress_balance_params_dict = generate_stress_balance(stress_balance, sb_params_dict)

        # Setup Climate Forcing
        climate_file = "../data_sets/climate_forcing/ltop_climate_olympics_{grid}m_dir_{dir}_kg_m-2_yr-1.nc".format(
            grid=grid, dir=wind_direction
        )
        climate_params_dict = generate_climate(
            climate,
            **{
                "atmosphere.yearly_cycle.file": climate_file,
                "atmosphere.elevation_change.file": climate_file,
                "atmosphere.elevation_change.temperature_lapse_rate": temperature_lapse_rate,
                "precip_adjustement": "scale",
                "atmosphere.precip_exponential_factor_for_temperature": precip_scale_factor,
                "atmosphere.delta_T.file": atmosphere_paleo_file,
                "atmosphere.precip_scaling.file": atmosphere_paleo_file,
            }
        )
        # Setup Ocean Forcing
        ocean_params_dict = generate_ocean("null")
        # Setup Hydrology Model
        hydro_params_dict = generate_hydrology(hydrology)

        calving_params_dict = generate_calving("float_kill")

        # Setup Scalar and Spatial Time Series Reporting
        scalar_ts_dict = generate_scalar_ts(outfile, tsstep, odir=dirs["scalar"])

        # Merge All Parameter Dictionaries
        all_params_dict = merge_dicts(
            general_params_dict,
            grid_params_dict,
            stress_balance_params_dict,
            climate_params_dict,
            ocean_params_dict,
            hydro_params_dict,
            calving_params_dict,
            scalar_ts_dict,
        )
        all_params = " ".join([" ".join(["-" + k, str(v)]) for k, v in all_params_dict.items()])

        if not spatial_ts == "none":
            exvars = spatial_ts_vars[spatial_ts]
            spatial_ts_dict = generate_spatial_ts(outfile, exvars, exstep, odir=dirs["spatial_tmp"], split=False)

            all_params_dict = merge_dicts(all_params_dict, spatial_ts_dict)

        all_params = " \\\n  ".join(["-{} {}".format(k, v) for k, v in list(all_params_dict.items())])

        if system == "debug":
            redirect = " "
        else:
            redirect = " > {jobs}/job.${job_id} 2>&1"

        template = "{mpido} {pism} {params}" + redirect

        context = merge_dicts(batch_system, dirs, {"pism": pism, "params": all_params})
        cmd = template.format(**context)

        f.write(cmd)
        f.write("\n")
        f.write("\n")
        f.write("{} {}\n".format(batch_system["submit"], script_post))
        f.write("\n")

    # post_header = make_batch_post_header(system)

    # with open(script_post, "w") as f:

    #     f.write(post_header)

    #     extra_file = spatial_ts_dict["extra_file"]
    #     myfiles = " ".join(["{}_{:.3f}.nc".format(extra_file, k) for k in np.arange(start + exstep, end, exstep)])
    #     myoutfile = extra_file + ".nc"
    #     myoutfile = os.path.join(odir, os.path.split(myoutfile)[-1])
    #     cmd = " ".join(["ncrcat -O -6 -h", myfiles, myoutfile, "\n"])
    #     f.write(cmd)
    #     cmd = " ".join(["ncks -O -4", os.path.join(odir, outfile), os.path.join(odir, outfile), "\n"])
    #     f.write(cmd)


scripts = uniquify_list(scripts)
# scripts_post = uniquify_list(scripts_post)
print("\n".join([script for script in scripts]))
print("\nwritten\n")
