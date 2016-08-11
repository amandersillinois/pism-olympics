#!/bin/bash
set -x -e 

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
