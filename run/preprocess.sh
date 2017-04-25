#!/bin/bash

ncap2 -O -s "delta_T=0.4*delta_T;" pism_dT.nc pism_scaled_dT.nc

for T in -7 -6 -5 -4; do
    ncap2 -O -s "delta_T(0)=${T};" paleo_modifier.nc paleo_modifier_${T}K.nc
done

