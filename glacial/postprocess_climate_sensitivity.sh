#!/bin/bash

odir=$1

pdir=tiffs
mkdir -p $odir/$pdir

idir=interface
mkdir -p $odir/$idir
for file in $odir/state/*0a.nc; do
    for var in thk velsurf_mag; do
        f="$(basename -- $file)"
        f_ne=${f%.nc}
        gdal_translate NETCDF:$file:$var $odir/$pdir/${var}_${f_ne}.tif
    done
done
