#!/usr/bin/env python

import numpy as np
import pandas as pd
import itertools

precip_scale_factor = [0.05, 0.07]
temperature_offset = [-7, -4]
phi_min = [20, 30]
sia_e = [1.0, 3.0]
ppq = [0.5]
tefo = [0.020]
phi_max = [45]
topg_min = [5000]
topg_max = [8000]
temperature_lapse_rate = [5.0, 6.5]


combinations = list(
    itertools.product(
        precip_scale_factor,
        temperature_offset,
        temperature_lapse_rate,
        sia_e,
        ppq,
        tefo,
        phi_min,
        phi_max,
        topg_min,
        topg_max,
    )
)

keys = [
    "precip_scale_factor",
    "temperature_offset",
    "temperature_lapse_rate",
    "sia_e",
    "ppq",
    "tefo",
    "phi_min",
    "phi_max",
    "topg_min",
    "topg_max",
]

# Convert to Pandas dataframe, append column headers, output as csv
df = pd.DataFrame(combinations)
df.to_csv("constant_climate.csv", header=keys, index=True, index_label="id")
