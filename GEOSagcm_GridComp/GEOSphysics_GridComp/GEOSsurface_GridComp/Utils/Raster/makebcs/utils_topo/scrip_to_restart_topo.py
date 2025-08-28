#!/usr/bin/env python

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

exclude = ['lon','lat']
for var in ncFid.variables:
    if var in exclude:
        continue

    invar = ncFid.variables[var]
    temp = invar[:]
    dim_size = temp.ndim

    # choose dims and create output var with the input _FillValue if present
    fillv_in = getattr(invar, '_FillValue', 1.0e15)
    if dim_size == 2:
        tout = ncFidOut.createVariable(var, 'f8', ('unknown_dim1','lat','lon'), fill_value=fillv_in)
    elif dim_size == 1:
        tout = ncFidOut.createVariable(var, 'f8', ('lat','lon'), fill_value=fillv_in)
    else:
        # unexpected shape: skip
        continue

    for att in invar.ncattrs():
        if att != "_FillValue":
            setattr(ncFidOut.variables[var], att, getattr(invar, att))

    # reshape into (6*cRes, cRes) [and lead nrdg dim if present]
    if dim_size == 2:
        data = numpy.reshape(temp, [rdgSize, 6*cRes, cRes]).astype('f8')
    else:
        data = numpy.reshape(temp, [6*cRes, cRes]).astype('f8')

    # --- variance -> std dev for SGH/SGH30 ---
    if var.lower() in ('sgh','sgh30'):
        miss_in = getattr(invar, 'missing_value', None)
        # build mask for invalids and non-physical values
        mask = ~numpy.isfinite(data) | (data < 0.0)
        if miss_in is not None:
            mask |= (data == miss_in)
        if fillv_in is not None:
            mask |= (data == fillv_in)

        out = numpy.empty_like(data)
        out[~mask] = numpy.sqrt(data[~mask])
        out[mask]  = fillv_in

        # write converted values
        if dim_size == 2:
            tout[:,:,:] = out
        else:
            tout[:,:] = out
    else:
        # passthrough
        if dim_size == 2:
            tout[:,:,:] = data
        else:
            tout[:,:] = data

#-----------------
# Closing the file
#----------------
ncFidOut.close()
ncFid.close()

