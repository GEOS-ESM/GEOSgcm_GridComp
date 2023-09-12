
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
import argparse
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

def write_import(fname, fr, Time=None):

    ntiles = fr.shape[1]
    ncat   = fr.shape[0]

    ncfile = Dataset(fname,'w')
    ncfile.createDimension('tile', ntiles)
    ncfile.createDimension('unknown_dim1', ncat)
    ncfile.createDimension('subtile', ncat)
    ncfile.createDimension('time', 1)

    chl = ncfile.createVariable('FRACICE',np.dtype('float32').char,('unknown_dim1','tile')) 

    chl[:] = fr 

    chl.units = "1"
    chl.long_name = 'ice_covered_fraction_of_tile'
    

    if Time is not None:
       t = ncfile.createVariable('time',np.dtype('float64').char,('time'))
       t.units =  "minutes since " + Time[:4]+'-'+Time[4:6]+'-'+Time[6:]+" 00:00:00" 
       t.begin_date = int(Time) 
       t.begin_time = 0
       t.long_name = "time" 
       t[:] = np.float64(0.0)
    ncfile.close()

def write_internal(fname, ti, si, Time=None):

    ntiles = ti.shape[1]
    ncat   = ti.shape[0]

    ncfile = Dataset(fname,'w')
    ncfile.createDimension('tile', ntiles)
    ncfile.createDimension('unknown_dim1', ncat)
    ncfile.createDimension('subtile', ncat)
    ncfile.createDimension('time', 1)

    chl = ncfile.createVariable('TSKINI',np.dtype('float32').char,('unknown_dim1','tile')) 
    chl[:] = ti 
    chl.units = "K"
    chl.long_name = 'ice_skin_temperature'

    chl = ncfile.createVariable('SSKINI',np.dtype('float32').char,('tile')) 
    chl[:] = si 
    chl.units = "psu"
    chl.long_name = 'ice_skin_salinity'
    

    if Time is not None:
       t = ncfile.createVariable('time',np.dtype('float64').char,('time'))
       t.units =  "minutes since " + Time[:4]+'-'+Time[4:6]+'-'+Time[6:]+" 00:00:00" 
       t.begin_date = int(Time) 
       t.begin_time = 0
       t.long_name = "time" 
       t[:] = np.float64(0.0)
    ncfile.close()


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
    parser.add_argument('-g', '--inputgrid', default=None, required=True, help='source grid file')
    parser.add_argument('-im', '--imp', default=None, required=True, help='import restart file on tiles')
    parser.add_argument('-in', '--internal', default=None, required=True, help='internal restart file on tiles')
    parser.add_argument('-t', '--tilefile', default=None, required=True, help='whether use BL99 fixed salinity profile')
    parser.add_argument('-d', '--date', default=None, required=True, help='date in format yyyymmdd')
    return parser.parse_args()

def main() -> None:

   args = parse_arguments()    

   LON, LAT, _, _, wet = get_grid(args.inputgrid)

   print(wet.shape)

   sw = saltwatertile(args.tilefile)  

   with Dataset(args.inputfile) as src:
       aicen = src['aicen'][:] 
       tsf = src['Tsfcn'][:] 
  

   fr = np.zeros((aicen.shape[0], sw.size), dtype='float32') 
   si = np.zeros(sw.size, dtype='float32') 
   si[:] = 30.0
   tskin = np.zeros((aicen.shape[0], sw.size), dtype='float32') 

   for k in range(sw.size):
      i, j = sw.gi[k]-1, sw.gj[k]-1
      if wet[j,i] > 0.5:
         fr[:,k]    = np.float32(aicen[:,j,i])       
         tskin[:,k] = np.float32(tsf[:,j,i]) + 273.15       
      else:
         tskin[:,k] = np.float32(273.16)      
 
         
       
   print(aicen.shape)
   write_import(args.imp, fr, Time=args.date)
   write_internal(args.internal, tskin, si, Time=args.date)

if __name__=="__main__":
   main()

