# """ Calculate the hypsographic curve from latitude longitude
# registered elevation data, and save the resulting data as csv file
# using this format:

# Elevation, CumSum
# -11000,    0

# Use 30 minute resulution to test, see
# https://www.pygmt.org/dev/api/generated/pygmt.datasets.load_earth_relief.html
# for other resolutions
# """

import time
import pygmt as gmt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed


def earth_radius(lat: float, lon: float) -> float:
    """Get earth radius as function of latitude"""
    from math import cos, sin, radians

    a: int = 6378137  # equatorial radius
    b: int = 6356752  # polar radius

    lat = radians(lat)
    lon = radians(lon)
    r: float = (
        ((a**2 * cos(lat)) ** 2 + (b**2 * sin(lat)) ** 2)
        / ((a * cos(lat)) ** 2 + (b * sin(lat)) ** 2)
    ) ** 0.5
    return r


def grid_area(lat: float, lon: float, size: float) -> float:
    """Calculate the area of a rectangular area of size = 1 deg at a given lat-long position"""
    from math import cos, sin, radians, pi

    lat1 = radians(lat + size)
    lat2 = radians(lat - size)
    lon1 = radians(lon + size)
    lon2 = radians(lon - size)

    r = earth_radius(lat, lon) / 1000
    dy = (size * r * pi) / 180
    dx = size / 180 * pi * r * cos(radians(lat))

    a = abs(dx * dy)
    return a

def slice_count(start, end, weight, grid) -> np.ndarray:
    """Generate elevation count array for each latitudinal slice"""

    # Assuming latitudes span from -90 to 90 degrees
    num_latitudes = grid.shape[0]
    latitudes = np.linspace(-90, 90, num_latitudes)

    # Determine the indices for the given latitude range
    start_index = np.searchsorted(latitudes, start, side="left")
    end_index = np.searchsorted(latitudes, end, side="right")

    # Slice the memmap object based on latitude indices
    sub_grid = grid[start_index:end_index, ...]

    # Calculate the cumulative occurances of each height interval
    # initialize two arrays
    count = np.arange(
        elevation_minimum, elevation_maximum, dz, dtype=float
    ) 
    elevation = np.arange(elevation_minimum, elevation_maximum, dz, dtype=int)

    # loop over elevation intervals
    i = 0
    for e in elevation:
        # this will return a latitudinal transect
        a = np.sum(np.logical_and(sub_grid > e, sub_grid < e + dz), axis=1)
        count[i] = np.sum(a * weight)
        i = i + 1

    return count


# Define the function to process each slice
def process_slice(i,lat):
    lat_slice = lat[(lat >= i) & (lat <= i + 10)]
    weight = np.array([grid_area(lat_val, 0, dx) for lat_val in lat_slice])

    # Compute and return the slice counts
    return slice_count(i, i + 10, weight, whole_grid)


if __name__ == "__main__":
    num_jobs = 8 # number of cpu cores. 1 to 8
    res = 5  # grid resolution in minutes
    r = "05m"
    dz = 1  # height interval in meters

    elevation_minimum = -11000  # relative to sea level
    elevation_maximum = 9000
    elevations = np.arange(elevation_minimum, elevation_maximum, dz, dtype=int)
    dx = res / 60  # convert resolution in degrees
    steps = int(round(180 / dx)) + 1  # +1 to ensure there is no empty array

    # Load the entire grid once
    whole_grid = gmt.datasets.load_earth_relief(
        resolution=r, region="-180/180/-90/90"
    ).to_numpy()

    # Generate lat value every dx degrees from south pole to north pole
    lat = np.linspace(
        -90, 90, steps, endpoint=True
    )

    start = time.time()
    results = Parallel(n_jobs=num_jobs)(
        delayed(process_slice)(i,lat) for i in range(-90, 90, 10))
    
    print(f"duration = {time.time() - start:.2f} seconds")

    # Aggregate the results
    total_count = np.zeros(len(elevations), dtype=float)
    for result in results:
        total_count += result

    # Normalize and calculate the cumulative sum
    total_count = total_count / np.sum(total_count)
    cum = np.cumsum(total_count)

    # Save dataframe
    df = pd.DataFrame({"Elevation": elevations, "CumSum": cum})
    df.sort_values(by=["Elevation"])

    # Save the CSV file
    df.to_csv(f"Hypsometric_Curve_{r}.csv", index=False, float_format="%.16f")

    # Check plot
    fig1, ax1 = plt.subplots()
    ax1.plot(1 - cum, elevations)
    ax1.set_ylabel("Elevation")
    ax1.set_xlabel("Area [cumulative %]")
    plt.show()
