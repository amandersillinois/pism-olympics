#!/bin/bash
set -x -e 

for GRID  in 50 100 200 250 500 1000; do
    outfile=ltop_climate_olympics_${GRID}m_kg_m-2_yr-1.nc
    python ../../scripts/linear_orog_precip.py -i ../bed_dem/pism_Olympics_${GRID}m_v1.nc --background_precip_pre 1 --background_precip_post 0.114 --precip_scale_factor 0.2 --tau_c 2000 --tau_f 2000 --wind_direction 270 --wind_magnitude 8 --moist_stability 0.005 --vapor_scale_height 2600 -o $outfile --o_units 'kg m-2 yr-1' --water_mass_rate
    ncrename -v Band1,precipitation $outfile
    ncap2 -O -s "air_temp_mean_annual[\$y,\$x]=9.5f; air_temp_mean_july[\$y,\$x]=15.5f; usurf[\$y,\$x]=0.f;" $outfile $outfile
    ncatted -a units,precipitation,o,c,"kg m-2 yr-1" -a units,air_temp_mean_july,o,c,"celsius" -a units,air_temp_mean_annual,o,c,"celsius" -a units,usurf,o,c,"m" -a standard_name,usurf,o,c,"surface_altitude" $outfile
done

exit

# Domain size
xmin=330000
ymin=5201000
xmax=530000
ymax=5381000


infile=preciptif.tif

for GRID  in 1000; do
    outfile=olympics_climate_${GRID}m.nc
    gdalwarp -overwrite -of netCDF -co "FORMAT=NC4" -r average -s_srs EPSG:26710 -t_srs EPSG:26710 -te $xmin $ymin $xmax $ymax -tr $GRID $GRID $infile $outfile
    ncrename -v Band1,precipitation $outfile
    ncap2 -O -s "precipitation=precipitation*(1000./910.); air_temp[\$y,\$x]=10.5f; usurf[\$y,\$x]=0.f;" $outfile $outfile
    ncatted -a units,precipitation,o,c,"mm year-1" -a units,air_temp,o,c,"celsius" -a units,usurf,o,c,"m" -a standard_name,usurf,o,c,"surface_altitude" $outfile
done

ncgen -o paleo_modifier.nc paleo_modifier.cdl
