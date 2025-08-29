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
    if var not in exclude:
        temp = ncFid.variables[var][:]
        dim_size =len(temp.shape)
        
        if dim_size == 2:
            tout = ncFidOut.createVariable(var,'f8',('unknown_dim1','lat','lon'),fill_value=1.0e15)
            for att in ncFid.variables[var].ncattrs():
                if att != "_FillValue":
                   setattr(ncFidOut.variables[var],att,getattr(ncFid.variables[var],att))
            temp2d = numpy.reshape(temp,[rdgSize,cRes*6,cRes])
            tout[:,:,:] = temp2d[:,:,:]

        elif dim_size == 1:
            tout = ncFidOut.createVariable(var,'f8',('lat','lon'),fill_value=1.0e15)
            for att in ncFid.variables[var].ncattrs():
                if att != "_FillValue":
                   setattr(ncFidOut.variables[var],att,getattr(ncFid.variables[var],att))
            temp1d = numpy.reshape(temp,[cRes*6,cRes])
            tout[:,:] = temp1d[:,:]
#-----------------
# Closing the file
#----------------
ncFidOut.close()
ncFid.close()

