
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
import argparse
import pprint
import scipy.interpolate as interp
from scipy.spatial import cKDTree
import scipy.optimize as optm
import scipy.stats as scstats
import os.path
import math
from scipy.io import netcdf
import calendar
import datetime
import os.path


def nearest_interp_new(lon, lat, z, LON, LAT):
    xs, ys, zs = lon_lat_to_cartesian(lon.flatten(), lat.flatten())
    xt, yt, zt = lon_lat_to_cartesian(LON.flatten(), LAT.flatten())
    #print(xs[:10], ys[:10], zs[:10])
    tree = cKDTree(list(zip(xs, ys, zs)))
    #find indices of the nearest neighbors in the flattened array
    #d, inds = tree.query(zip(xt, yt, zt), k = 1)
    #get interpolated 2d field
    #zout = LON.copy().flatten()

    d, inds = tree.query(list(zip(xt, yt, zt)), k = 1)
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

    x = R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z


def get_src_grid(fname): #reads lat lon for tripolar ocean grid 
    ncfile=Dataset(fname, "r")
    LON     = ncfile.variables['lon'][:]
    LAT     = ncfile.variables['lat'][:]
    bat     = ncfile.variables['Bathymetry'][:]
    #wet     = ncfile.variables['mask'][:]
    ncfile.close()
    #wet = np.zeros(bat.shape)
    #wet[bat>0.0] = 1.0 
    return LON, LAT, bat

def get_dst_grid(fname): #reads lat lon for tripolar ocean grid 
    ncfile=Dataset(fname, "r")
    LON     = ncfile.variables['lon_centers'][:]
    LAT     = ncfile.variables['lat_centers'][:]
    wet     = ncfile.variables['mask'][:]
    #wet     = ncfile.variables['mask'][:]
    ncfile.close()

    return LON, LAT, wet

missing=np.float32(-32767.0)

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputfile', default=None, required=True, help='CICE restart file on source grid')
    parser.add_argument('-ig', '--inputgrid', default=None, required=False, help='source grid file')
    parser.add_argument('-o', '--outputfile', default=None, required=True, help='CICE restart file on target grid')
    parser.add_argument('-og', '--outputgrid', default=None, required=True, help='target grid file')
    return parser.parse_args()



def main() -> None:

   args = parse_arguments()    

   #print(args.inputfile)
   #print(args.outputfile)

   LON, LAT, wet = get_dst_grid(args.outputgrid)

   jm, im = LON.shape
   #print(im, jm)
   if args.inputgrid:
       lons, lats, mask = get_src_grid(args.inputgrid)

   with Dataset(args.inputfile) as src, Dataset(args.outputfile, "w") as dst:
     # copy global attributes all at once via dictionary
      dst.setncatts(src.__dict__)
    # copy dimensions
      for name, dimension in src.dimensions.items():
         if name == 'ni':  
            dst.createDimension(name, (im))
         elif name == 'nj':  
            dst.createDimension(name, (jm))
         else:
            dst.createDimension(
                 name, (len(dimension) if not dimension.isunlimited() else None))
         #print(name, (len(dimension) if not dimension.isunlimited() else None)) 
      #print(dst.dimensions) 
    # copy all file data except for the excluded
      for name, variable in src.variables.items():
        if len(variable.dimensions) == 3:
            x = dst.createVariable(name, variable.datatype,  ('ncat', 'nj', 'ni',))
        else:
            x = dst.createVariable(name, variable.datatype,  ('nj', 'ni',))
        #print(variable.dimensions)   
        #dst[name][:] = src[name][:]
        if 'vel' in name or 'stress' in name or 'strocn' in name:
           dst[name][:] = 0.0 
        else:
           msk = mask.mask 
           lon_in = lons[~msk]  
           lat_in = lats[~msk]  
           if len(variable.dimensions) == 3:
              var = src[name][:]
              for k in range(var.shape[0]):
                  h_in = var[k][~msk]  
                  #print('interpolating time slice:', k)
                  hout = nearest_interp_new(lon_in, lat_in, h_in, LON, LAT)
                  hout[wet<0.5] = 0.0
                  dst[name][k,:,:] = hout    
           else:
              var = src[name][:]
              h_in = var[~msk]  
              hout = nearest_interp_new(lon_in, lat_in, h_in, LON, LAT)
              hout[wet<0.5] = 0.0
              dst[name][:] = hout    
        # copy variable attributes all at once via dictionary
        dst[name].setncatts(src[name].__dict__)


if __name__=="__main__":
   main()

