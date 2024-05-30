# # """ Calculate the hypsographic curve from latitude longitude
# # registered elevation data, and save the resulting data as csv file
# # using this format:

# # Elevation, CumSum
# # -11000,    0

# # Use 30 minute resulution to test, see
# # https://www.pygmt.org/dev/api/generated/pygmt.datasets.load_earth_relief.html
# # for other resolutions
# # """

import time
import pygmt as gmt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing
from math import cos, sin, radians, pi

def earth_radius(lat: float) -> float:
    """Get earth radius as function of latitude"""
    a = 6378137  # equatorial radius
    b = 6356752  # polar radius

    lat = radians(lat)
    r = (
        ((a**2 * cos(lat)) ** 2 + (b**2 * sin(lat)) ** 2)
        / ((a * cos(lat)) ** 2 + (b * sin(lat)) ** 2)
    ) ** 0.5
    return r

def grid_area(lat: float, size: float) -> float:
    """Calculate the area of a rectangular area of size = 1 deg at a given lat-long position"""
    r = earth_radius(lat) / 1000
    dy = (size * r * pi) / 180
    dx = size / 180 * pi * r * cos(radians(lat))
    return abs(dx * dy)

def slice_count(start: int, end: int, weight: np.ndarray, grid: np.ndarray, elevation_minimum: int, elevation_maximum: int, dz: int) -> np.ndarray:
    """Generate elevation count array for each latitudinal slice"""
    sub_grid = grid[start:end, ...]

    count = np.zeros((elevation_maximum - elevation_minimum) // dz, dtype=float)
    elevation = np.arange(elevation_minimum, elevation_maximum, dz, dtype=int)

    for i, e in enumerate(elevation):
        a = np.sum(np.logical_and(sub_grid > e, sub_grid < e + dz), axis=1)
        count[i] = np.sum(a * weight)

    return count

def process_slice(start: int, end: int, lat: np.ndarray, grid: np.ndarray, dz: int, elevation_minimum: int, elevation_maximum: int, dx: float) -> np.ndarray:
    """Process a slice of the grid"""
    lat_slice = lat[start:end]
    weight = np.array([grid_area(lat_val, dx) for lat_val in lat_slice])
    return slice_count(start, end, weight, grid, elevation_minimum, elevation_maximum, dz)

if __name__ == "__main__":
    num_jobs = 8  # Number of CPU cores available
    res = 30  # Grid resolution in minutes
    r = "30m"
    dz = 1  # Height interval in meters

    elevation_minimum = -11000  # Relative to sea level
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
        delayed(process_slice)(slice_ranges[i], slice_ranges[i + 1], lat, whole_grid, dz, elevation_minimum, elevation_maximum, dx) for i in range(total_slices)
    )
    
    print(f"duration = {time.time() - start:.2f} seconds")

    total_count = np.zeros(len(elevations), dtype=float)
    for result in results:
        total_count += result

    total_count = total_count / np.sum(total_count)
    cum = np.cumsum(total_count)

    df = pd.DataFrame({"Elevation": elevations, "CumSum": cum})
    df.sort_values(by=["Elevation"], inplace=True)

    df.to_csv(f"Hypsometric_Curve_{r}.csv", index=False, float_format="%.16f")

    fig1, ax1 = plt.subplots()
    ax1.plot(1 - cum, elevations)
    ax1.set_ylabel("Elevation")
    ax1.set_xlabel("Area [cumulative %]")
    plt.show()
