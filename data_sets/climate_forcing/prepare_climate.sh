#!/bin/bash
set -x -e 

#for GRID  in 50 100 200 250 500 1000 2000; do
for GRID  in 1000; do
    for dir in 220; do
        outfile=ltop_climate_olympics_${GRID}m_dir_${dir}_kg_m-2_yr-1.nc
        python ../../scripts/linear_orog_precip.py -i ../bed_dem/pism_Olympics_${GRID}m_v1.nc --background_precip_pre 1 --background_precip_post 0.057 --precip_scale_factor 0.125 --tau_c 1750 --tau_f 1750 --wind_direction $dir --wind_magnitude 15 --moist_stability 0.007 --vapor_scale_height 2600 -o $outfile --o_units 'kg m-2 yr-1' --water_mass_rate
        gdalwarp -dstnodata 0 -overwrite -of netCDF -cutline ../../shp_files/victoria-cutout.shp $outfile tmp_${GRID}m_dir_${dir}_kg_m-2_yr-1.nc
        ncks -O tmp_${GRID}m_dir_${dir}_kg_m-2_yr-1.nc $outfile
        ncrename -v Band1,precipitation $outfile
        ncap2 -O -s 'defdim("time",1);time[$time]=.5;time@long_name="Time";time@units="yr"; air_temp_mean_annual[$time,$y,$x]=9.5f; air_temp_mean_july[$time,$y,$x]=15.5f; usurf[$time,$y,$x]=0.f;' $outfile $outfile
        ncatted -a _FillValue,precipitation,d,, -a units,precipitation,o,c,"kg m-2 yr-1" -a units,air_temp_mean_july,o,c,"celsius" -a units,air_temp_mean_annual,o,c,"celsius" -a units,usurf,o,c,"m" -a standard_name,usurf,o,c,"surface_altitude" $outfile
    done
done

ncgen -o paleo_modifier.nc paleo_modifier.cdl
