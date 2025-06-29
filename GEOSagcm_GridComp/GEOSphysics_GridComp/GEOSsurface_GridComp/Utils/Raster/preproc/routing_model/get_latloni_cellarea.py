import sys  
import numpy as np
import os
from netCDF4 import Dataset

#Main purpose: Computes grid-cell index arrays and per-cell areas for 1-m high-res grid.

lat36_file, lon36_file, lat09_file, lon09_file, lat1m_file, lon1m_file = sys.argv[1:7]

lati36_output_file = "temp/lati_1m_M36.txt"
loni36_output_file = "temp/loni_1m_M36.txt"
lati09_output_file = "temp/lati_1m_M09.txt"
loni09_output_file = "temp/loni_1m_M09.txt"
cellarea_output_file = "temp/cellarea.nc"

# Grid dimensions
nlat36, nlon36 = 406, 964
nlat1m, nlon1m = 10800, 21600
nlat09, nlon09 = 1624, 3856

# Read data
lat36 = np.loadtxt(lat36_file, dtype=float, max_rows=nlat36)
lon36 = np.loadtxt(lon36_file, dtype=float, max_rows=nlon36)
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
lati36 = ind_nearest_coord(lat1m, lat36)
loni36 = ind_nearest_coord(lon1m, lon36)
lati09 = ind_nearest_coord(lat1m, lat09)
loni09 = ind_nearest_coord(lon1m, lon09)

# Save indices to files (1-based index)
np.savetxt(lati36_output_file, lati36 + 1, fmt='%d')
np.savetxt(loni36_output_file, loni36 + 1, fmt='%d')
np.savetxt(lati09_output_file, lati09 + 1, fmt='%d')
np.savetxt(loni09_output_file, loni09 + 1, fmt='%d')

# Compute global grid cell area
def area_global_rectilinear_grid(lat, lon, rearth=6371.22):
    """
    Calculate the approximate area of each grid cell on a global rectilinear grid.
    
    Parameters:
    lat : numpy.ndarray
        Array of latitude values (degrees).
    lon : numpy.ndarray
        Array of longitude values (degrees).
    rearth : float
        Earth radius in kilometers. Default is 6371.22 km.
        
    Returns:
    area_grid : numpy.ndarray
        2D array representing the area of each grid cell (km^2).
    """
    # Convert degrees to radians
    rad = np.pi / 180.0
    rr = rearth * rad

    # Longitude spacing (constant across latitudes)
    dlon = rr * (lon[1] - lon[0])  # Assuming uniform spacing in longitude

    # Compute longitude spacing at each latitude (dx)
    dx = dlon * np.cos(lat * rad)
    
    # Handle rounding issues at poles
    dx[0] = 0.0 if lat[0] < -89.9999 else dx[0]
    dx[-1] = 0.0 if lat[-1] > 89.9999 else dx[-1]

    # Latitude spacing (dy), can be variable
    dy = np.zeros_like(lat)
    dy[0] = (lat[1] - lat[0]) * rr
    dy[1:-1] = (lat[2:] - lat[:-2]) * rr / 2.0
    dy[-1] = (lat[-1] - lat[-2]) * rr

    # Area per latitude band
    area_lat = dx * dy

    # Extend latitude areas to all longitudes
    area_grid = np.outer(area_lat, np.ones(len(lon)))

    # Total area of all grid cells
    area_total = np.sum(area_grid)

    # Total surface area of the sphere
    area_sphere = 4.0 * np.pi * (rearth ** 2)

    # Add metadata as a dictionary
    metadata = {
        "long_name": "area of each grid cell",
        "units": "km^2",
        "area_total": area_total,
        "area_lat": area_lat,
        "rearth": rearth,
        "area_sphere": area_sphere,
        "area_ratio": area_total / area_sphere
    }

    return area_grid, metadata

# Calculate and save cell area
area, metadata = area_global_rectilinear_grid(lat1m, lon1m)
area *= 1.e6  # Convert to m²

# Remove existing file and write new cell area to NetCDF
if os.path.exists(cellarea_output_file):
    os.remove(cellarea_output_file)

with Dataset(cellarea_output_file, "w", format="NETCDF4") as fout:
    # Create dimensions
    fout.createDimension("lat", nlat1m)
    fout.createDimension("lon", nlon1m)
    
    # Create variables for lat and lon
    lat_var = fout.createVariable("lat", "f4", ("lat",))
    lon_var = fout.createVariable("lon", "f4", ("lon",))
    lat_var[:] = lat1m
    lon_var[:] = lon1m
    # Assign units attribute to lat and lon
    lat_var.units = "degrees_north"
    lon_var.units = "degrees_east"

    # Create the area variable
    area_var = fout.createVariable("data", "f8", ("lat", "lon"))
    area_var[:] = area
    area_var.units = "m2"
