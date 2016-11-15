#!/usr/bin/env python
import numpy as np
import pandas as pa

csv_file = "KUIL_1966080200_2014102623.txt"


df = pa.read_csv(csv_file)

month = []
rh = []

for n in range(0, 2726041):
    str = df.validUTC[n]
    temp = int(str[5:7])
    month.append(temp)
    if df.dwpc[n] != 'M' and df.tmpc[n] != 'M':
           temp2 = np.exp((17.625*float(df.dwpc[n]))/(243.04+float(df.dwpc[n])))
           temp3 = np.exp((17.625*float(df.tmpc[n]))/(243.04+float(df.tmpc[n])))
           temp4 = 100*(temp2/temp3)
           rh.append(temp4)
    else:
           rh.append('M')
                  

df['month'] = month
df['rh'] = rh

df = df.replace('M',np.nan)
df = df.convert_objects(convert_numeric = True)



    
jan = df[df.month == 1]
jan_storms = jan[jan.rh >=85]
jan_storms = jan_storms[jan_storms.height_m >= 1000]
jan_storms = jan_storms[jan_storms.height_m <= 2000]
jan_storms = jan_storms[jan_storms.bearing >= 160]
jan_storms = jan_storms[jan_storms.bearing <= 330]


jan_t = jan_storms["tmpc"].mean()
jan_speed = 0.5144*jan_storms["speed_kts"].mean()
jan_bearing = jan_storms["bearing"].mean()

winter_end = df[df.month >= 10]
winter_start = df[df.month <= 3]
join = [winter_end, winter_start]
winter = pa.concat(join)
winter_storms = winter[winter.rh >=85]
winter_storms = winter_storms[winter_storms.height_m >= 1000]
winter_storms = winter_storms[winter_storms.height_m <= 2000]
winter_storms = winter_storms[winter_storms.bearing >= 160]
winter_storms = winter_storms[winter_storms.bearing <= 330]


winter_t = winter_storms["tmpc"].mean()
winter_speed = 0.5144*winter_storms["speed_kts"].mean()
winter_bearing = winter_storms["bearing"].mean()
