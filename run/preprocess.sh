#!/bin/bash

for T in -7 -6 -5 -4; do
    ncap2 -O -s "delta_T(0)=${T};" paleo_modifier.nc paleo_modifier_${T}K.nc
done

