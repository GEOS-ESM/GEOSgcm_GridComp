#!/usr/bin/env python

#-------------
# Load modules
#-------------
from netCDF4 import Dataset
import numpy
import argparse

def parse_args():
    p = argparse.ArgumentParser(description='convert old style cube to new style cube input')
    p.add_argument('-i','--input',type=str,help='input file',default=None)
    p.add_argument('-e','--example',type=str,help='example file',default=None)
    p.add_argument('-o','--output',type=str,help='output file',default=None)
    p.add_argument('-v','--vars',type=str,help='output only these',default=None,nargs='+')
    return vars(p.parse_args())

#------------------
# Opening the file
#------------------
comm_args    = parse_args()
Input_file   = comm_args['input']
Output_file  = comm_args['output']
Example_file  = comm_args['example']
only_vars = comm_args['vars']

ncFid = Dataset(Input_file, mode='r')
ncFidEx = Dataset(Example_file, mode='r')
ncFidOut = Dataset(Output_file, mode='w', format='NETCDF4')

#---------------------
# Extracting variables
#---------------------

ntiles = len(ncFid.dimensions['ncol'])
haveRdg = False
for dim in ncFid.dimensions:
    if dim == 'nrdg':
           haveRdg = True
           rdgSize = len(ncFid.dimensions['nrdg'])


cRes = len(ncFidEx.dimensions['Xdim'])

Xdim = ncFidOut.createDimension('Xdim',cRes)
Ydim = ncFidOut.createDimension('Ydim',cRes)
nf = ncFidOut.createDimension('nf',6)
ncontact = ncFidOut.createDimension('contact',4)

if haveRdg:
   rdgOut = ncFidOut.createDimension('nrdg',rdgSize)

vXdim = ncFidOut.createVariable('Xdim','f8',('Xdim'))
vYdim = ncFidOut.createVariable('Ydim','f8',('Ydim'))
setattr(ncFidOut.variables['Xdim'],'units','degrees_east')
setattr(ncFidOut.variables['Ydim'],'units','degrees_north')
setattr(ncFidOut.variables['Xdim'],'long_name','Fake Longitude for GrADS Compatibility')
setattr(ncFidOut.variables['Ydim'],'long_name','Fake Latitude for GrADS Compatibility')
vXdim[:]=range(1,cRes+1)
vYdim[:]=range(1,cRes+1)
vnf = ncFidOut.createVariable('nf','i4',('nf'))
vnf[:]=range(1,7)
setattr(ncFidOut.variables['nf'],'long_name','cubed-sphere face')
setattr(ncFidOut.variables['nf'],'axis','e')
setattr(ncFidOut.variables['nf'],'grads_dim','e')

vchar = ncFidOut.createVariable('cubed_sphere','S1')
setattr(ncFidOut.variables['cubed_sphere'],'grid_mapping_name','gnomonic cubed-sphere')
setattr(ncFidOut.variables['cubed_sphere'],'file_format_version','2.90')
setattr(ncFidOut.variables['cubed_sphere'],'additional_vars','contacts,orientation,anchor')

temp1d = numpy.zeros([6,cRes,cRes])
if haveRdg: 
   temp2d = numpy.zeros([rdgSize,6,cRes,cRes])

if only_vars == None:
   only_vars = ncFid.variables

for var in only_vars:
    temp = ncFid.variables[var][:]
    dim_size =len(temp.shape)
    
    if dim_size == 2:
        tout = ncFidOut.createVariable(var,'f8',('nrdg','nf','Ydim','Xdim'),fill_value=1.0e15)
        for att in ncFid.variables[var].ncattrs():
           if att != "_FillValue":
              setattr(ncFidOut.variables[var],att,getattr(ncFid.variables[var],att))
        setattr(ncFidOut.variables[var],'grid_mapping','cubed_sphere')
        setattr(ncFidOut.variables[var],'coordinates','lons lats')
        temp2d = numpy.reshape(temp,[rdgSize,6,cRes,cRes])
        tout[:,::,:] = temp2d[:,:,:,:]

    elif dim_size == 1: 
        tout = ncFidOut.createVariable(var,'f8',('nf','Ydim','Xdim'),fill_value=1.0e15)
        for att in ncFid.variables[var].ncattrs():
           if att != "_FillValue":
              setattr(ncFidOut.variables[var],att,getattr(ncFid.variables[var],att))
        setattr(ncFidOut.variables[var],'grid_mapping','cubed_sphere')
        setattr(ncFidOut.variables[var],'coordinates','lons lats')
        temp1d = numpy.reshape(temp,[6,cRes,cRes])
        tout[:,:,:] = temp1d[:,:,:]

XCdim = ncFidOut.createDimension('XCdim',cRes+1)
YCdim = ncFidOut.createDimension('YCdim',cRes+1)
center_lons = ncFidOut.createVariable('lons','f8',('nf','Ydim','Xdim'))
setattr(ncFidOut.variables['lons'],'long_name','longitude')
setattr(ncFidOut.variables['lons'],'units','degrees_east')
center_lats = ncFidOut.createVariable('lats','f8',('nf','Ydim','Xdim'))
setattr(ncFidOut.variables['lats'],'long_name','latitude')
setattr(ncFidOut.variables['lats'],'units','degrees_north')

center_lons[:,:,:] = ncFidEx.variables['lons'][:,:,:]
center_lats[:,:,:] = ncFidEx.variables['lats'][:,:,:]


corner_lons = ncFidOut.createVariable('corner_lons','f8',('nf','YCdim','XCdim'))
setattr(ncFidOut.variables['corner_lons'],'long_name','longitude')
setattr(ncFidOut.variables['corner_lons'],'units','degrees_east')
corner_lats = ncFidOut.createVariable('corner_lats','f8',('nf','YCdim','XCdim'))
setattr(ncFidOut.variables['corner_lats'],'long_name','latitude')
setattr(ncFidOut.variables['corner_lats'],'units','degrees_north')

corner_lons[:,:,:] = ncFidEx.variables['corner_lons'][:,:,:]
corner_lats[:,:,:] = ncFidEx.variables['corner_lats'][:,:,:]
#-----------------
# Closing the file
#-----------------
ncFidEx.close()
ncFidOut.close()
ncFid.close()

