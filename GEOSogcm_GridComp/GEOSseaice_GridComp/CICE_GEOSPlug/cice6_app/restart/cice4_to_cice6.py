
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import array
import matplotlib.cm as cm
import glob
import struct
import time
import sys
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
import calendar
import datetime
import os.path


c0 = np.float64(0.0)
c1 = np.float64(1.0)
c2 = np.float64(2.0)
c3 = np.float64(3.0)
c4 = np.float64(4.0)
puny = np.float64(1.0e-11)
c1000 = np.float64(1000.0)

rhos      = np.float64(330.0)   # density of snow (kg/m^3)
rhoi      = np.float64(917.0)   # density of ice (kg/m^3)
rhosi     = np.float64(940.0)   # average sea ice density
                    # Cox and Weeks, 1982: 919-974 kg/m^2
rhow      = np.float64(1026.0)  # density of seawater (kg/m^3)
rhofresh  = np.float64(1000.0)  # density of fresh water (kg/m^3)

cp_ice    = np.float64(2106.)  # specific heat of fresh ice (J/kg/K)
cp_ocn    = np.float64(4218.)  # specific heat of ocn    (J/kg/K)
                   # freshwater value needed for enthalpy
depressT = np.float64(0.054)
Lsub      = np.float64(2.835e6) # latent heat, sublimation freshwater (J/kg)
Lvap      = np.float64(2.501e6)  # latent heat, vaporization freshwater (J/kg)
Lfresh    = Lsub-Lvap

# liquidus relation - higher temperature region
az1_liq = np.float64(-18.48)
bz1_liq =  np.float64(0.0)

# liquidus relation - lower temperature region
az2_liq = np.float64(-10.3085)
bz2_liq =     np.float64(62.4)

# liquidus break
Tb_liq = np.float64(-7.6362968855167352) # temperature of liquidus break
Sb_liq =  np.float64(123.66702800276086)    # salinity of liquidus break

# basic liquidus relation constants
az1p_liq = az1_liq / c1000
bz1p_liq = bz1_liq / c1000
az2p_liq = az2_liq / c1000
bz2p_liq = bz2_liq / c1000

ki = np.float64(2.3) # fresh ice conductivity (W m-1 K-1)
kb = np.float64(0.5375) # brine conductivity (W m-1 K-1)


# quadratic constants - higher temperature region
AS1_liq = az1p_liq * (rhow * cp_ocn - rhoi * cp_ice)
AC1_liq = rhoi * cp_ice * az1_liq
BS1_liq = (c1 + bz1p_liq) * (rhow * cp_ocn - rhoi * cp_ice)  \
           + rhoi * Lfresh * az1p_liq
BQ1_liq = -az1_liq
BC1_liq = rhoi * cp_ice * bz1_liq - rhoi * Lfresh * az1_liq
CS1_liq = rhoi * Lfresh * (c1 + bz1p_liq)
CQ1_liq = -bz1_liq
CC1_liq = -rhoi * Lfresh * bz1_liq

# quadratic constants - lower temperature region
AS2_liq = az2p_liq * (rhow * cp_ocn - rhoi * cp_ice)
AC2_liq = rhoi * cp_ice * az2_liq
BS2_liq = (c1 + bz2p_liq) * (rhow * cp_ocn - rhoi * cp_ice)  \
           + rhoi * Lfresh * az2p_liq
BQ2_liq = -az2_liq
BC2_liq = rhoi * cp_ice * bz2_liq - rhoi * Lfresh * az2_liq
CS2_liq = rhoi * Lfresh * (c1 + bz2p_liq)
CQ2_liq = -bz2_liq
CC2_liq = -rhoi * Lfresh * bz2_liq

# break enthalpy constants
D_liq = ((c1 + az1p_liq*Tb_liq + bz1p_liq)    \
          / (   az1_liq*Tb_liq + bz1_liq))    \
        * ((cp_ocn*rhow - cp_ice*rhoi)*Tb_liq + Lfresh*rhoi)
E_liq = cp_ice*rhoi*Tb_liq - Lfresh*rhoi

#  just fully melted enthapy constants
F1_liq = (  -c1000 * cp_ocn * rhow) / az1_liq
G1_liq =    -c1000
H1_liq = (-bz1_liq * cp_ocn * rhow) / az1_liq
F2_liq = (  -c1000 * cp_ocn * rhow) / az2_liq
G2_liq =    -c1000
H2_liq = (-bz2_liq * cp_ocn * rhow) / az2_liq


# warmer than fully melted constants
I_liq = c1 / (cp_ocn * rhow)

def icepack_enthalpy_temperature_bl99(zTin, Tmlt):
    '''
    '''

    zQin = np.zeros(zTin.shape, dtype='float64')

    msk = zTin < c0

    # the T dereived from mushy layer scheme could go above Tmlt
    # reset to just below Tmlt if that is the case
    zTin[zTin >= Tmlt] = Tmlt*1.05

    zQin[msk] = -rhoi*(cp_ice*(Tmlt-zTin[msk]) + 
                   Lfresh*(c1 - Tmlt/zTin[msk]) -cp_ocn*Tmlt)

    return zQin  



def icepack_mushy_temperature_mush(zqin, zSin):
    '''
    taken from icepack_mushy_physics.F90 

    '''

    zTin = np.zeros(zqin.shape, dtype='float64')
    # just melted enthalpy
    # S_low = merge(c1, c0, (zSin < Sb_liq))
    S_low = np.zeros(zSin.shape, dtype='float64')
    S_low[:] = c0
    S_low[zSin < Sb_liq] = c1 

    q0 = ((F1_liq * zSin) / (G1_liq + zSin) + H1_liq) * S_low +  \
         ((F2_liq * zSin) / (G2_liq + zSin) + H2_liq) * (c1 - S_low)
    #q_melt = merge(c1, c0, (zqin > q0))
    q_melt = np.zeros(zqin.shape, dtype='float64')
    q_melt[:] = c0
    q_melt[zqin > q0] = c1 
     

    # break enthalpy
    qb = D_liq * zSin + E_liq
    # t_high = merge(c1, c0, (zqin > qb))
    t_high = np.zeros(zqin.shape, dtype='float64')
    t_high[:] = c0
    t_high[zqin > qb] = c1 
    t_low = c1 - t_high

    # quadratic values
    A = (AS1_liq * zSin                 + AC1_liq) * t_high + \
        (AS2_liq * zSin                 + AC2_liq) * t_low

    B = (BS1_liq * zSin + BQ1_liq * zqin + BC1_liq) * t_high + \
        (BS2_liq * zSin + BQ2_liq * zqin + BC2_liq) * t_low

    C = (CS1_liq * zSin + CQ1_liq * zqin + CC1_liq) * t_high + \
        (CS2_liq * zSin + CQ2_liq * zqin + CC2_liq) * t_low

    zTin = (-B + np.sqrt(np.maximum(B**2 - c4 * A * C,puny))) / (c2 * A)

    # change T if all melted
    zTin = q_melt * zqin * I_liq + (c1 - q_melt) * zTin

    return zTin



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
       self.lons = icetile['lon'][:]
       self.lats = icetile['lat'][:]
       #print 'hhh: ',self.size,self.gi[-1],self.gj[-1]
       #return icetile

def get_src_grid(fname): #reads lat lon for tripolar ocean grid 
    ncfile  = Dataset(fname, "r")
    try:
        LON     = ncfile.variables['lon'][:]
        LAT     = ncfile.variables['lat'][:]
    except:
        pass   
    try:
        LON     = ncfile.variables['TLON'][:]
        LAT     = ncfile.variables['TLAT'][:]
    except:
        pass   
    bat     = ncfile.variables['Bathymetry'][:]
    #wet     = ncfile.variables['mask'][:]
    ncfile.close()
    #wet = np.zeros(bat.shape)
    #wet[bat>0.0] = 1.0 
    return LON, LAT, bat

def get_dst_grid(fname): #reads lat lon for tripolar ocean grid 
    ncfile  = Dataset(fname, "r")
    LON     = ncfile.variables['lon_centers'][:]
    LAT     = ncfile.variables['lat_centers'][:]
    ULON    = ncfile.variables['lon_corners'][1:, 1:]
    ULAT    = ncfile.variables['lat_corners'][1:, 1:]
    wet     = ncfile.variables['mask'][:]
    #wet     = ncfile.variables['mask'][:]
    ncfile.close()

    return LON, LAT, ULON, ULAT, wet

missing=np.float32(-32767.0)

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputfile', default=None, required=True, help='GEOS seaice thermo restart file on tile grid')
    parser.add_argument('-ig', '--inputgrid', default=None, required=True, help='source tile file')
    parser.add_argument('-o', '--outputfile', default=None, required=True, help='CICE restart file on target grid')
    parser.add_argument('-ot', '--outputtemplate', default=None, required=True, help='CICE restart file on target output grid servign as a template')
    parser.add_argument('-og', '--outputgrid', default=None, required=True, help='target grid file')
    return parser.parse_args()



def main() -> None:

   args = parse_arguments()    

   #print(args.inputfile)
   #print(args.outputfile)

   LON, LAT, ULON, ULAT, wet = get_dst_grid(args.outputgrid)

   jm, im = LON.shape
   #print(im, jm)
   #lons, lats, mask = get_src_grid(args.inputgrid)
   sw = saltwatertile(args.inputgrid) 

   #zSin = np.zeros(5, dtype='float64')
   #zQin = np.zeros(5, dtype='float64')
   #zTin = icepack_mushy_temperature_mush(zQin, zSin) 
   #print(zTin)   
   saltmax = np.float64(3.2)
   nsal =  np.float64(0.407)
   msal = np.float64(0.573)
   salinz = np.zeros((nilyr), dtype='float64')
   Tmlt = np.zeros((nilyr), dtype='float64')
   for k in range(nilyr):
      zn = np.float64((k+1-0.5)/nilyr)
      salinz[k] = (saltmax/2.0)*(1.0-np.cos(np.pi*np.power(zn, nsal/(msal+zn))))
      Tmlt[k] = -depressT * salinz[k] 
   print(salinz)


   with Dataset(args.inputfile) as src, Dataset(args.outputfile, "w") as dst, \
        Dataset(args.outputtemplate) as tpl:
     # copy global attributes all at once via dictionary
      dst.setncatts(tpl.__dict__)

      ncat = src.dimensions['subtile'][:]
      nilyr = src.dimensions['unknown_dim2'][:]
      nslyr = src.dimensions['unknown_dim1'][:]
    # copy dimensions
      for name, dimension in tpl.dimensions.items():
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
      for name, variable in tpl.variables.items():
        if len(variable.dimensions) == 3:
            x = dst.createVariable(name, variable.datatype,  ('ncat', 'nj', 'ni',))
        else:
            x = dst.createVariable(name, variable.datatype,  ('nj', 'ni',))
        #print(variable.dimensions)   
        print('processing: ', name)   
        #dst[name][:] = src[name][:]
        if 'vel' in name or 'stress' in name or 'strocn' in name or 'iceu' in name:
           dst[name][:] = c0 
        elif 'ulon' == name:
           dst[name][:] = ULON
        elif 'ulat' == name:
           dst[name][:] = ULAT
        elif 'tlon' == name:
           dst[name][:] = LON
        elif 'tlat' == name:
           dst[name][:] = LAT
        elif 'sice' in name:
           k = int(name[4:]) - 1 # layer index
           var = tpl[name][:]
           for i in range(var.shape[0]):
              dst[name][i,:,:] = salinz[k]
              dst[name][i][wet<0.5] = 0.0
        elif 'qice' in name:
           k = int(name[4:]) - 1 # layer index
           msk = mask.mask 
           lon_in = lons[~msk]  
           lat_in = lats[~msk]  
           var = ti[k]
           for i in range(var.shape[0]):
               h_in = var[i][~msk]  
               hout = nearest_interp_new(lon_in, lat_in, h_in, LON, LAT)
               qin = icepack_enthalpy_temperature_bl99(hout, Tmlt[k]) 
               qin[wet<0.5] = 0.0
               dst[name][i,:,:] = qin    
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

