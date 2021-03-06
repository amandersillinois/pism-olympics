#!/usr/bin/env python
import numpy as np
import pandas as pa
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('FILE', nargs=1,
                    help='CSV input file name')

options = parser.parse_args()
csv_file = options.FILE[0]

df = pa.read_csv(csv_file)
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())
