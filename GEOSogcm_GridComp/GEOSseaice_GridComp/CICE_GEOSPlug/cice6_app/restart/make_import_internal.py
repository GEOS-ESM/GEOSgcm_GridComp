
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
#sys.path.append('/home/bzhao/python_utils')
#import read_utils
#import plot_utils
#import math_utils
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

    nx=HICE.shape[1]
    ny=HICE.shape[0]

    ncfile = Dataset(fname,'w')
    #ncfile.createDimension('xaxis_1',nx)
    #ncfile.createDimension('yaxis_1',ny)
    ncfile.createDimension('nx',nx)
    ncfile.createDimension('ny',ny)

    chl = ncfile.createVariable('tideamp',np.dtype('float32').char,('ny','nx'), 
                             fill_value=missing)

    chl[:] = HICE

    chl.missing_value = missing
    chl.units = "m s-1"
    

    #mask = ncfile.createVariable('MASK',np.dtype('float32').char,('Time','yaxis_1','xaxis_1'))
    #mask[:] = MASK
    if LON is not None:
       lon = ncfile.createVariable('nx',np.dtype('float64').char,('nx'))
       lon.units = 'degrees_east' 
       lon.cartesian_axis = "X" 
       lon[:] = LON
    if LAT is not None:
       lat = ncfile.createVariable('ny',np.dtype('float64').char,('ny'))
       lat.units = 'degrees_north' 
       lat.cartesian_axis = "Y" 
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
    zout = LON.copy().flatten()

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

class saltwatertile:

    def __init__(self, file): 
         
       header = np.genfromtxt(file, dtype='i4', usecols=(0), max_rows=8)
       #print header
       self.atm = 'x'.join([str(x) for x in header[3:5]])
       self.ocn = 'x'.join([str(x) for x in header[6:]])
       self.nx, self.ny = header[6], header[7]
       print(self.atm, self.ocn) 
       tile=np.genfromtxt(file, dtype=[('type','i1'), ('area','f8'), ('lon','f8'),('lat','f8'), ('gi1','i4'),
                           ('gj1','i4'), ('gw1','f8'),
                           ('idum','i4'), ('gi2','i4'), ('gj2','i4'), ('gw2','f8')], skip_header=8)
       n1=0
       n2=0
       for n in range(1, tile.shape[0]+1, 1):
           if tile[n-1][0] == 0:
               n1 = n
               break
       #print n1
       for n in range(n1, tile.shape[0]+1, 1):
           if tile[n-1][0] != 0:
               n2 = n
               break
       #print n2
       icetile=tile[n1-1:n2-1]
       #print icetile.shape
       #print 'hhh: ',icetile[0][2], icetile[-1][2]
       self.size = icetile.shape[0]
       self.gi = icetile['gi2'][:]
       self.gj = icetile['gj2'][:]
       #print 'hhh: ',self.size,self.gi[-1],self.gj[-1]
       #return icetile

def get_grid(fname): #reads lat lon for tripolar ocean grid 
    ncfile  = Dataset(fname, "r")
    LON     = ncfile.variables['lon_centers'][:]
    LAT     = ncfile.variables['lat_centers'][:]
    ULON    = ncfile.variables['lon_corners'][1:, 1:]
    ULAT    = ncfile.variables['lat_corners'][1:, 1:]
    wet     = ncfile.variables['mask'][:]
    #wet     = ncfile.variables['mask'][:]
    ncfile.close()

    return LON, LAT, ULON, ULAT, wet

missing=np.float32(-1.e20)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputfile', default=None, required=True, help='CICE restart file on source grid')
    parser.add_argument('-ig', '--inputgrid', default=None, required=True, help='source grid file')
    parser.add_argument('-im', '--import', default=None, required=True, help='import restart file on tiles')
    parser.add_argument('-in', '--internal', default=None, required=True, help='internal restart file on tiles')
    parser.add_argument('-ti', '--tilefile', default=None, required=True, help='whether use BL99 fixed salinity profile')
    return parser.parse_args()

def main() -> None:

   args = parse_arguments()    

   LON, LAT, _, _, wet = get_grid(args.inputgrid)

   sw = saltwatertile(args.tilefile)  


if __name__=="__main__":
   main()

