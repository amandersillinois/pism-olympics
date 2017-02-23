#!/bin/bash

for T in -5 -4 -3; do
    ncap2 -O -s "delta_T(0)=${T};" paleo_modifier.nc paleo_modifier_${T}K.nc
done

