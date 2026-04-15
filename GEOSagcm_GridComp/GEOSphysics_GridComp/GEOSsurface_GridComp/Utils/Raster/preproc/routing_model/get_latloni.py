import sys  
import numpy as np
import os
from netCDF4 import Dataset
import routing_model_constants

#Main purpose: Computes grid-cell index arrays and per-cell areas for 1-m high-res grid.

lat09_file, lon09_file, lat1m_file, lon1m_file = sys.argv[1:5]

lati09_output_file = "temp/lati_1m_M09.txt"
loni09_output_file = "temp/loni_1m_M09.txt"

# Grid dimensions
nlat1m, nlon1m = routing_model_constants.nlat1m, routing_model_constants.nlon1m
nlat09, nlon09 = routing_model_constants.nlat09, routing_model_constants.nlon09

# Read data
lat1m = np.loadtxt(lat1m_file, dtype=float, max_rows=nlat1m)
lon1m = np.loadtxt(lon1m_file, dtype=float, max_rows=nlon1m)
lat09 = np.loadtxt(lat09_file, dtype=float, max_rows=nlat09)
lon09 = np.loadtxt(lon09_file, dtype=float, max_rows=nlon09)

# Define nearest coordinate function
def ind_nearest_coord(coord_array1, coord_array2):
    """
    Find the index of the nearest value in coord_array2 for each value in coord_array1.
    """
    indices = []
    for coord in coord_array1:
        index = np.argmin(np.abs(coord_array2 - coord))
        indices.append(index)
    return np.array(indices)

# Find nearest coordinates
lati09 = ind_nearest_coord(lat1m, lat09)
loni09 = ind_nearest_coord(lon1m, lon09)

# Save indices to files (1-based index)
np.savetxt(lati09_output_file, lati09 + 1, fmt='%d')
np.savetxt(loni09_output_file, loni09 + 1, fmt='%d')
