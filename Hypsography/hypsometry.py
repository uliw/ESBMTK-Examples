"""
     ESBMTK-EXAMPLES Hypsography.hypsometry

     Calculate the hypsographic curve from latitude longitude
     registered elevation data, and save the resulting data as csv file
     using this format:
     Elevation, CumSum
     -11000,    0
     
     Use 30 minute resulution to test, see
     https://www.pygmt.org/dev/api/generated/pygmt.datasets.load_earth_relief.html
     for other resolutions

     Copyright (C), 2024 Jingwen (Lisa) Zhong

     This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
 
     You should have received a copy of the GNU General Public License
     along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import time
import pygmt as gmt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from esbmtk.sealevel import process_slice

num_jobs = 8  # Number of CPU cores
res = 5  # Grid resolution in minutes
r = "05m"
dz = 100  # Height interval in meters

elevation_minimum = -11000 -dz  # Relative to sea level
elevation_maximum = 9000
elevations = np.arange(elevation_minimum, elevation_maximum, dz, dtype=int)
dx = res / 60  # Convert resolution to degrees
steps = int(round(180 / dx)) + 1  # +1 to ensure there is no empty array

total_slices = num_jobs

# Load the entire grid once
whole_grid = gmt.datasets.load_earth_relief(
    resolution=r, region="-180/180/-90/90"
).to_numpy()

lat = np.linspace(-90, 90, steps, endpoint=True)

slice_ranges = np.linspace(0, len(lat), total_slices + 1, dtype=int)

start = time.time()

results = Parallel(n_jobs=num_jobs)(
    delayed(process_slice)(
        slice_ranges[i],
        slice_ranges[i + 1],
        lat,
        whole_grid,
        dz,
        elevation_minimum,
        elevation_maximum,
        elevations,
        dx,
    )
    for i in range(total_slices)
)

print(f"duration = {time.time() - start:.2f} seconds")

total_count = np.zeros(len(elevations), dtype=float)
for result in results:
    total_count += result

total_count = total_count / np.sum(total_count)
cum = np.cumsum(total_count)

df = pd.DataFrame({"Elevation": elevations+dz, "CumSum": cum})
df.sort_values(by=["Elevation"], inplace=True)

df.to_csv(f"Hypsometric_Curve_{r}_{dz}.csv", index=False, float_format="%.16f")

fig1, ax1 = plt.subplots()
ax1.plot(1 - cum, elevations)
ax1.set_ylabel("Elevation")
ax1.set_xlabel("Area [cumulative %]")
plt.show()
