import rasterio
import numpy as np
import netCDF4 as nc

# Files
tif_file = "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/soil/SOIL-DATA/v2/peatGPA22WGS_2cl.tif"
nc4_file = "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/soil/SOIL-DATA/v2/peatGPA22WGS_2cl.nc4"

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

