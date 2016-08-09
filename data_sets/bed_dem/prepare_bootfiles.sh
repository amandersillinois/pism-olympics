#!/bin/bash
set -x -e 

# Domain size
xmin=330000
ymin=5201000
xmax=530000
ymax=5381000

infile=Olympics_30m.tif

for GRID  in 50 100 200 500 1000 2000 5000; do
    outfile=pism_Olympics_${GRID}m.nc
    gdalwarp -overwrite -of netCDF -co "FORMAT=NC4" -r average -s_srs EPSG:26710 -t_srs EPSG:26710 -te $xmin $ymin $xmax $ymax -tr $GRID $GRID $infile $outfile
    ncrename -v Band1,topg $outfile
    ncatted -a standard_name,topg,o,c,"bedrock_altitude" -a units,topg,o,c,"m" -a proj4,global,o,c,"+init=EPSG:26710" $outfile
done
