#!/bin/bash
set -x -e 

# Wind direction
dir=220

for GRID  in 100 200 250 500 1000; do
    outfile=pism_ltop_climate_olympics_${GRID}m_dir_${dir}_kg_m-2_yr-1.nc
    python ~/pism/sources/examples/python/orographic_precipitation.py \
           -i ../bed_dem/pism_Olympics_${GRID}m_v1.nc \
           -o ${outfile} \
           -atmosphere.orographic_precipitation.background_precip_post  0.057 \
           -atmosphere.orographic_precipitation.background_precip_pre  1.0 \
           -atmosphere.orographic_precipitation.conversion_time  1750.0 \
           -atmosphere.orographic_precipitation.coriolis_latitude  0.0 \
           -atmosphere.orographic_precipitation.fallout_time  1750.0 \
           -atmosphere.orographic_precipitation.lapse_rate  -5.8 \
           -atmosphere.orographic_precipitation.moist_adiabatic_lapse_rate  -6.5 \
           -atmosphere.orographic_precipitation.moist_stability_frequency  0.007 \
           -atmosphere.orographic_precipitation.reference_density  7.4e-3 \
           -atmosphere.orographic_precipitation.scale_factor  0.125 \
           -atmosphere.orographic_precipitation.truncate  "true" \
           -atmosphere.orographic_precipitation.water_vapor_scale_height  2600.0 \
           -atmosphere.orographic_precipitation.wind_direction  ${dir} \
           -atmosphere.orographic_precipitation.wind_speed  15 \               

    ncap2 -O -s 'air_temp_mean_annual[$time,$y,$x]=9.5f; air_temp_mean_july[$time,$y,$x]=15.5f; air_temp_mean_summer[$time,$y,$x]=15.5f; usurf[$time,$y,$x]=0.f;' $outfile $outfile
    ncatted -a _FillValue,precipitation,d,, -a units,precipitation,o,c,"kg m-2 yr-1" -a units,air_temp_mean_july,o,c,"celsius" -a units,air_temp_mean_summer,o,c,"celsius" -a units,air_temp_mean_annual,o,c,"celsius" -a units,usurf,o,c,"m" -a standard_name,usurf,o,c,"surface_altitude" $outfile
done

ncgen -o paleo_modifier.nc paleo_modifier.cdl
