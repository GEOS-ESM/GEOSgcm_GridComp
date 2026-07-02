#!/usr/bin/env python3
"""
scrip_to_restart_topo.py

Purpose:
    Convert old-style cubed-sphere NetCDF files with a flat 'ncol' dimension
    into new-style cube files with explicit (lat, lon) dimensions.
    Optionally annotate stretched grids (sg001, sg002) with STRETCH_FACTOR
    and target center information.

Key points:
    - Input: NetCDF with dimension 'ncol' (and possibly 'nrdg')
    - Output: NetCDF with dimensions ('lat','lon'[,'unknown_dim1'])
    - Reshapes variables from (ncol) → (lat, lon) or (nrdg,ncol) → (nrdg,lat,lon)
    - Adds lon/lat index coordinate variables (1..IM, 1..6*IM)
    - Sets global attrs for sg001/sg002 grids
"""
#-------------
# Load modules
#-------------
from netCDF4 import Dataset
import numpy
import argparse
import math

def parse_args():
    p = argparse.ArgumentParser(description='convert old style cube to new style cube input')
    p.add_argument('-i','--input',type=str,help='input file',default=None)
    p.add_argument('-o','--output',type=str,help='output file',default=None)
    p.add_argument('-g','--grid',type=str,help='grid identifier (e.g., sg001, sg002)',default=None)
    return vars(p.parse_args())

#------------------
# Opening the file
#------------------
comm_args    = parse_args()
Input_file   = comm_args['input']
Output_file  = comm_args['output']

ncFid = Dataset(Input_file, mode='r')
ncFidOut = Dataset(Output_file, mode='w', format='NETCDF4')

grid_type = comm_args['grid']

if grid_type == 'sg001':
    ncFidOut.STRETCH_FACTOR = 2.5
    ncFidOut.TARGET_LAT = 39.5
    ncFidOut.TARGET_LON = -98.35
elif grid_type == 'sg002':
    ncFidOut.STRETCH_FACTOR = 3.0
    ncFidOut.TARGET_LAT = 39.5
    ncFidOut.TARGET_LON = -98.35

#---------------------
# Extracting variables
#---------------------

ntiles = len(ncFid.dimensions['ncol'])
haveRdg = False
for dim in ncFid.dimensions:
    if dim == 'nrdg':
           haveRdg = True
           rdgSize = len(ncFid.dimensions['nrdg'])


cRes = ntiles/6
cRes = int(math.sqrt(cRes))

Xdim = ncFidOut.createDimension('lon',cRes)
Ydim = ncFidOut.createDimension('lat',cRes*6)

if haveRdg:
   rdgOut = ncFidOut.createDimension('unknown_dim1',rdgSize)

vXdim = ncFidOut.createVariable('lon','f8',('lon'))
vYdim = ncFidOut.createVariable('lat','f8',('lat'))
setattr(ncFidOut.variables['lon'],'units','degrees_east')
setattr(ncFidOut.variables['lat'],'units','degrees_north')
setattr(ncFidOut.variables['lon'],'long_name','Longitude')
setattr(ncFidOut.variables['lat'],'long_name','Latitude')
vXdim[:]=range(1,cRes+1)
vYdim[:]=range(1,(6*cRes)+1)

temp1d = numpy.zeros([6*cRes,cRes])
if haveRdg: 
   temp2d = numpy.zeros([rdgSize,6*cRes,cRes])

exclude = ['lon', 'lat']
for var in ncFid.variables:
    if var in exclude:
        continue

    v = ncFid.variables[var]
    temp = v[:]
    nd = temp.ndim

    # Choose fill value at creation time
    is_angle = var.upper() in ('ANGLX', 'ANGLL')
    fv = -9999.0 if is_angle else 1.0e15

    # Create destination variable with correct fill immediately
    if nd == 2:  # (nrdg, ncol) -> (unknown_dim1, lat, lon)
        tout = ncFidOut.createVariable(
            var, 'f8', ('unknown_dim1', 'lat', 'lon'),
            fill_value=fv
        )
    elif nd == 1:  # (ncol) -> (lat, lon)
        tout = ncFidOut.createVariable(
            var, 'f8', ('lat', 'lon'),
            fill_value=fv
        )
    else:
        # unexpected rank — skip safely
        continue

    # Copy attributes verbatim EXCEPT _FillValue (already set)
    for att in v.ncattrs():
        if att != '_FillValue':
            setattr(tout, att, getattr(v, att))

    # For angle variables, make the metadata sentinel explicit
    if is_angle:
        setattr(tout, 'missing_value', -9999.0)

    # Simple reshape write (no masking/NaN munging in this tool)
    if nd == 2:
        tout[:, :, :] = numpy.reshape(temp, (len(ncFid.dimensions['nrdg']),
                                             int((len(ncFid.dimensions['ncol'])//6)**0.5)*6,
                                             int((len(ncFid.dimensions['ncol'])//6)**0.5)))
    else:
        tout[:, :] = numpy.reshape(temp,
                                   (int((len(ncFid.dimensions['ncol'])//6)**0.5)*6,
                                    int((len(ncFid.dimensions['ncol'])//6)**0.5)))

#-----------------
# Closing the file
#----------------
ncFidOut.close()
ncFid.close()

