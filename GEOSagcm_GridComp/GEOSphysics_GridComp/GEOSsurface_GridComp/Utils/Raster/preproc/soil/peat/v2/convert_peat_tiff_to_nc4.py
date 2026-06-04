#!/usr/bin/env python3

"""
Standalone utility: convert the original GPM 2.0 GeoTIFF to NetCDF4.

This script is NOT required by the peat preprocessing workflow.

The production preprocessing script, build_gpa22_peat_rasters_30arcsec.py,
reads peatGPA22WGS_2cl.tif directly, remaps it to the standard 30-arcsec
global raster grid, and writes the final preprocessing input:

    peatGPA22WGS_2cl_real_30arcsec.nc4

This tif-to-nc4 script is kept only as a convenience tool in case someone
wants a NetCDF4 copy of the original GeoTIFF for inspection or debugging.
"""

import rasterio
import numpy as np
import netCDF4 as nc

# Files
tif_file = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/peat/v2/peatGPA22WGS_2cl.tif"
nc4_file = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/peat/v2/peatGPA22WGS_2cl.nc4"

# Read the GeoTIFF and get some info
with rasterio.open(tif_file) as src:
    print("CRS:", src.crs)
    print("Resolution:", src.res)

    data = src.read(1)

    transform = src.transform
    width, height = src.width, src.height
    res_x, res_y = src.res

# Generate explicit lon-lat coordinates
lon = np.linspace(transform.c + res_x / 2, transform.c + res_x * (width - 0.5), width)
lat = np.linspace(transform.f - res_y / 2, transform.f - res_y * (height - 0.5), height)

# Create NetCDF
with nc.Dataset(nc4_file, 'w', format='NETCDF4') as dst:
    # Create dimensions
    dst.createDimension('lon', len(lon))
    dst.createDimension('lat', len(lat))

    # Create variables
    longitude = dst.createVariable('lon', 'f4', ('lon',))
    latitude = dst.createVariable('lat', 'f4', ('lat',))
    peatland_type = dst.createVariable('peatland_type', 'i1', ('lat', 'lon'), zlib=True)

    # Assign data
    longitude[:] = lon
    latitude[:] = lat
    peatland_type[:, :] = data

    # Add attributes
    dst.description = 'Global Peatland Map 2022 (GPM 2.0) - peat dominated (1), peat in soil mosaic (2)'
    longitude.units = 'degrees_east'
    latitude.units = 'degrees_north'
    peatland_type.description = '1=peat dominated, 2=peat in soil mosaic'

