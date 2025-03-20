import numpy as np
from netCDF4 import Dataset
import os

# Load data
nt = 112573
nlat = 406
nlon = 964
nc = 291284

# Read input data from text files
lat_bot = np.loadtxt("temp/lat_bottom_M36.txt", dtype=float)
lat_up = np.loadtxt("temp/lat_upper_M36.txt", dtype=float)
lon_left = np.loadtxt("temp/lon_left_M36.txt", dtype=float)
lon_right = np.loadtxt("temp/lon_right_M36.txt", dtype=float)

# Calculate the center latitudes and longitudes
latc = (lat_bot + lat_up) / 2.0
lonc = (lon_left + lon_right) / 2.0

# Read latitudes and longitudes for the grid
lat36m = np.loadtxt("input/lat_M36.txt", dtype=float)
lon36m = np.loadtxt("input/lon_M36.txt", dtype=float)

# Find the nearest coordinates
def ind_nearest_coord(coord_array1, coord_array2):
    """
    Find the index of the nearest value in coord_array2 for each value in coord_array1.
    """
    indices = []
    for coord in coord_array1:
        index = np.argmin(np.abs(coord_array2 - coord))
        indices.append(index)
    return np.array(indices)

lati = ind_nearest_coord(latc, lat36m)
loni = ind_nearest_coord(lonc, lon36m)

# Save the indices to files (1-based index)
np.savetxt("temp/lati_tile_M36.txt", lati + 1, fmt='%d')
np.savetxt("temp/loni_tile_M36.txt", loni + 1, fmt='%d')

# Initialize the map_tile array
map_tile = np.full((nlat, nlon), -9999, dtype=int)

# Fill the map_tile with data
for i in range(nt):
    map_tile[lati[i], loni[i]] = i + 1

# Remove the existing file if it exists
if os.path.exists("temp/map_tile_M36.nc"):
    os.remove("temp/map_tile_M36.nc")

# Create a NetCDF file and write the data
with Dataset("temp/map_tile_M36.nc", "w", format="NETCDF4") as fout:
    # Create dimensions
    fout.createDimension("lat", nlat)
    fout.createDimension("lon", nlon)

    # Create variable to store the map_tile data with fill_value set during creation
    map_tile_var = fout.createVariable("data", "i4", ("lat", "lon"), fill_value=-9999)
    map_tile_var[:] = map_tile

# Print a sample of the map_tile data
#print(map_tile[62, 10])  # Corresponds to map_tile(63-1, 11-1) in NCL (1-based to 0-based)

