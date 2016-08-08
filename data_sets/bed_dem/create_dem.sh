#!/bin/bash

xmin=330000
ymin=5201000
xmax=530000
ymax=5381000

infile=Olympics_30m.tif

for GRID  in 1000; do
    outfile=pism_Olympics_${GRID}m.nc
    gdalwarp -overwrite -of netCDF -co "FORMAT=NC4" -r average -s_srs EPSG:26710 -t_srs EPSG:26710 -te $xmin $ymin $xmax $ymax -tr $GRID $GRID $infile $outfile

done
