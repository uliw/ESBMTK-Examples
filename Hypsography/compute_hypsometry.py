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
import numpy.typing as npt
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing
from math import cos, sin, radians, pi

# declare numpy float64 array type
NDArrayFloat = npt.NDArray[np.float64]

def earth_radius(
    lat: float,
) -> float:
    """Get earth radius as function of latitude"""
    a = 6378137  # equatorial radius
    b = 6356752  # polar radius

    lat = radians(lat)
    r = (
        ((a**2 * cos(lat)) ** 2 + (b**2 * sin(lat)) ** 2)
        / ((a * cos(lat)) ** 2 + (b * sin(lat)) ** 2)
    ) ** 0.5
    return r


def grid_area(
    lat: float,
    size: float,
) -> float:
    """Calculate the area of a rectangular area of size = 1 deg at a given lat-long position"""
    r = earth_radius(lat) / 1000
    dy = (size * r * pi) / 180
    dx = size / 180 * pi * r * cos(radians(lat))
    return abs(dx * dy)


def slice_count(
    start: int,
    end: int,
    weight: NDArrayFloat,
    grid: NDArrayFloat,
    elevation_minimum: int,
    elevation_maximum: int,
    elevations: NDArrayFloat,
    dz: int,
) -> NDArrayFloat:
    """Generate elevation count array for each latitudinal slice"""
    sub_grid = grid[start:end, ...]

    count = np.zeros((elevation_maximum - elevation_minimum) // dz, dtype=float)

    for i, e in enumerate(elevations):
        a = np.sum(np.logical_and(sub_grid > e, sub_grid < e + dz), axis=1)
        count[i] = np.sum(a * weight)

    return count


def process_slice(
    start: int,
    end: int,
    lat: NDArrayFloat,
    grid: NDArrayFloat,
    dz: int,
    elevation_minimum: int,
    elevation_maximum: int,
    elevations: NDArrayFloat,
    dx: float,
) -> NDArrayFloat:
    """Process a slice of the grid"""
    lat_slice = lat[start:end]
    weight = np.array([grid_area(lat_val, dx) for lat_val in lat_slice])
    return slice_count(
        start, end, weight, grid, elevation_minimum, elevation_maximum, elevations, dz
    )