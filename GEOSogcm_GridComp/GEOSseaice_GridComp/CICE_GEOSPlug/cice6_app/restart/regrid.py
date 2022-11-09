
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import array
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import glob
import struct
import time
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from scipy import interpolate
import getopt
import string
import datetime
import scipy.interpolate as interp
from scipy.spatial import cKDTree
import scipy.optimize as optm
import scipy.stats as scstats
import os.path
import math
from scipy.io import netcdf
import calendar
import datetime
sys.path.append('/home/bzhao/python_utils')
import read_utils
import plot_utils
import math_utils
import os.path

#### usage:
## run for 1987 jan whole month: 

##     python interpo_aice_hice.py 1987 1 1 30

#####
## this script is modified from /discover/nobackup/zli7/geos5/Sub_seasonal/Obs_to_Tri_interp/20160301_test/interpo_aice_hice.py
## major changes include: write out hice and seaice in different files
#####
##### 
## current version of this scripts are only good until year 2014
## for year 2016 has different input files for hice (heff)
## to run for 2016 onwards, gmask and gm need to be used (uncomment line 360 to 366)
#####
def write_state_aice(fname, AICE):

    nx=AICE.shape[1]
    ny=AICE.shape[0]
    nt=1

    ncfile = Dataset(fname,'w')
    ncfile.createDimension('xaxis_1',nx)
    ncfile.createDimension('yaxis_1',ny)
    ncfile.createDimension('Time',nt)

    aice = ncfile.createVariable('AICE',np.dtype('float32').char,('Time','yaxis_1','xaxis_1'))

    aice[:] = AICE
    ncfile.close()

def write_state_chl(fname, HICE, LON=None, LAT=None, Time=None):

    nx=HICE.shape[2]
    ny=HICE.shape[1]
    nt=HICE.shape[0]

    ncfile = Dataset(fname,'w')
    #ncfile.createDimension('xaxis_1',nx)
    #ncfile.createDimension('yaxis_1',ny)
    ncfile.createDimension('i',nx)
    ncfile.createDimension('j',ny)
    ncfile.createDimension('time',None)

    chl = ncfile.createVariable('chlor_a',np.dtype('float32').char,('time','j','i'), 
                             fill_value=missing)

    chl[:] = HICE

    chl.missing_value = missing
    chl.long_name = "Chlorophyll Concentration, OCI Algorithm"
    chl.standard_name = "mass_concentration_chlorophyll_concentration_in_sea_water" 
    chl.units = "mg m^-3"
    chl.valid_min = np.float32(0.001)
    chl.valid_max = 100. 
    chl.reference = '''Hu, C., Lee Z., and Franz, B.A. (2012). Chlorophyll-a algorithms for o
ligotrophic oceans: A novel approach based on three-band reflectance difference, J. Geophys. Res., 117, C01
011, doi:10.1029/2011JC007395.'''
    

    var_i = ncfile.createVariable('i',np.dtype('float32').char,('i'))
    var_i.long_name = "Grid position along first dimension"
    var_i.cartesian_axis = "X"
    var_i[:] = np.array([x+0.5 for x in range(nx)]) 

    var_j = ncfile.createVariable('j',np.dtype('float32').char,('j'))
    var_j.long_name = "Grid position along second dimension"
    var_j.cartesian_axis = "Y"
    var_j[:] = np.array([y+0.5 for y in range(ny)]) 

    #mask = ncfile.createVariable('MASK',np.dtype('float32').char,('Time','yaxis_1','xaxis_1'))
    #mask[:] = MASK
    if LON is not None:
       lon = ncfile.createVariable('lon',np.dtype('float32').char,('j','i'))
       lon.units = 'degrees_E' 
       lon[:] = LON
    if LAT is not None:
       lat = ncfile.createVariable('lat',np.dtype('float32').char,('j','i'))
       lat.units = 'degrees_N' 
       lat[:] = LAT
    if Time is not None:
       t = ncfile.createVariable('time',np.dtype('float32').char,('time'))
       t.units = "days since 0001-01-01 00:00:00" 
       t.calendar = "noleap" 
       t.long_name = "Time" 
       t.cartesian_axis = "T"
       t.modulo = " " 
       t[:] = Time
    ncfile.close()

def nearest_interp_new(lon, lat, z, LON, LAT):
    lon[lon>80.0]=lon[lon>80.0]-360.0
    xs, ys, zs = lon_lat_to_cartesian(lon.flatten(), lat.flatten())
    xt, yt, zt = lon_lat_to_cartesian(LON.flatten(), LAT.flatten())
    tree = cKDTree(zip(xs, ys, zs))
    #find indices of the nearest neighbors in the flattened array
    #d, inds = tree.query(zip(xt, yt, zt), k = 1)
    #get interpolated 2d field
    #zout = LON.copy().flatten()

    d, inds = tree.query(zip(xt, yt, zt), k = 1)
    #rs=tree.query_ball_point(zip(xt, yt, zt), r)
    #zout = np.sum(w * z.flatten()[inds], axis=1) / np.sum(w, axis=1)
    zout = z.flatten()[inds]
    ''' 
    for i,r in enumerate(rs):                                                   
        if r:
            zout[i] = np.mean(z[r])
        else:
            zout[i] = 0.0 
    '''
    zout.shape = LON.shape
    return zout

def lon_lat_to_cartesian(lon, lat, R = 1):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z



def get_dst_grid(): #reads lat lon for tripolar ocean grid 
    ##ncfile=Dataset('/gpfsm/dnb42/projects/p17/gvernier/SAND_BOXES/PLOT_ODAS/DATA/grid_spec_720x410x40.nc', "r")
    #ncfile=Dataset('/discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM6/1440x1080_V3/MAPL_Tripolar.nc',"r")
    ncfile=Dataset('/discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM6/1440x1080_newtopo/MAPL_Tripolar.nc',"r")
    #ncfile=Dataset('/gpfsm/dnb42/projects/p17/bzhao/cp100/scratch/INPUT/grid_spec.nc',"r")
    LON     = ncfile.variables['lon_centers'][:]
    LAT     = ncfile.variables['lat_centers'][:]
    #wet     = ncfile.variables['mask'][:]
    ncfile.close()
    ncfile=Dataset('/discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM6/1440x1080_newtopo/INPUT/ocean_topog.nc',"r")
    wet     = ncfile.variables['wet'][:]
    ncfile.close()

    return LON, LAT, wet

missing=np.float32(-32767.0)


LON, LAT, wet = get_dst_grid()
#fileout_hice = 'GIOMAS_HICE.201601.nc'
fileout_hice = sys.argv[1] 
filein_hice = '/discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM6/1440x1080/INPUT/seawifs-clim-1997-2010.1440x1080.v20180328.nc' 

#radius = float(sys.argv[2])

ncfile=Dataset(filein_hice,'r')
lons=ncfile.variables['lon'][:]
lats=ncfile.variables['lat'][:]
time=ncfile.variables['time'][:]
hice=ncfile.variables['chlor_a'][:]
print hice.shape
ncfile.close()
#print lonvar.long_name, lonvar.units, lonvar.modulo, lonvar.axis

hicenew=np.zeros(hice.shape, dtype='float32')

n = hice.shape[0]
for k in range(n):
   mask = hice[k].mask 
   hice_in = hice[k][~mask]  
   lon_in = lons[~mask]  
   lat_in = lats[~mask]  
#mask[LAT>=0.0] = 1.0
#mask[LAT<0.0]  = 0.0
   print 'interpolating time slice: ', k
   hice_out = nearest_interp_new(lon_in, lat_in, hice_in, LON, LAT)
  
   hice_out[wet < 0.5] = missing 
   hicenew[k] = hice_out[:]
#hice_out[LAT<0.0] = 0.0

write_state_chl(fileout_hice, hicenew, LON=LON, LAT=LAT, Time=time) 




