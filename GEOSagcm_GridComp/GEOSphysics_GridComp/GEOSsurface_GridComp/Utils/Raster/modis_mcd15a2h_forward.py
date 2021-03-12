#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 09:34:21 2020

@author: smahanam
"""

from bs4 import BeautifulSoup
import requests
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from datetime import datetime, timedelta
import os
import re, sys
from pyhdf.SD import SD, SDC
import shutil

#----------------------------------------------------------#
#                 BEGIN USER DEFINED VARIABLES             #

"""
MODIS LAI Product Name (M6_NAME):
    MOD15A2H (Terra   ) 20000218 (DOY : 049) - to date   
    MYD15A2H (Auqa    ) 20020704 (DOY : 185) - to date
    MCD15A2H (combined) 20020704 (DOY : 185) - to date
"""

M6_NAME    = 'MCD15A2H'
FORWARD    = True
# Specify output grid resolution if REGRID
IM      = 43200
JM      = 21600
OUTDIR  = '/discover/nobackup/projects/lis/LS_PARAMETERS/MODIS/' + M6_NAME + '.006/'
MASKFILE= '/discover/nobackup/projects/lis/LS_PARAMETERS/MODIS/MCD15A2H.006/land_mask.nc4'
# If you want to process a single MODIS date, define below
#    YYYYDOY - the format is YYYY + '/' + DOY (for e.g. 2016/065)
# If you want to loop through: YYYYDOY='' (an empty string)
YYYYDOY= ""

#                  END USER DEFINED VARIABLES              # 
#----------------------------------------------------------#

# ---- Global Parameters

NC         = 86400
NR         = 43200
DXY        = np.double(360.)/NC
N_MODIS    = 2400  
CWD        =  os.getcwd()
MAPL_UNDEF = 255
MODIS_PATH = "https://ladsweb.modaps.eosdis.nasa.gov/opendap/hyrax/allData/6/" + M6_NAME + '/'

if M6_NAME == 'MOD15A2H':
    firstdate = '20000222'
if M6_NAME == 'MCD15A2H':
    firstdate = '20020708'

M6_DIR ={'MCD15A2H':'MOTA', 'MOD15A2H':'MOLT','MYD15A2H':'MOLA'}
MODIS_DOWN = "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/" + M6_NAME + '/'
IM = NC
JM = NR

DY = np.double(180.) / JM 
DX = np.double(360.) / IM


    
command = 'module load nco/4.8.1'
os.system(command)
os.chdir(CWD)

#----------------------------------------------------------#
#                     MODULE FUNCTIONS                     #                    
#----------------------------------------------------------#

class DriverFunctions (object):

    def create_netcdf (FILE_NAME,VAR_NAMES):
        import datetime
        ncFidOut = Dataset(FILE_NAME,'w',format='NETCDF4')
        LatDim  = ncFidOut.createDimension('lat', JM)
        LonDim  = ncFidOut.createDimension('lon', IM)
        
        ncFidOut.description = "MODIS " + M6_NAME + " @ " + str(DXY*3600) + ' arc-sec'

        ncFidOut.history     = "Created on " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " by Kristi Arsenault (kristi.r.arsenault@nasa.gov)"
        
        lonout  = ncFidOut.createVariable('lon','f8',('lon',))
        latout  = ncFidOut.createVariable('lat','f8',('lat',))

        for l in range (len(VAR_NAMES)):
            varout = ncFidOut.createVariable(VAR_NAMES[l], 'u1', ('lat','lon'), fill_value=MAPL_UNDEF) 
            setattr(ncFidOut.variables[VAR_NAMES[l]],'missing_value' ,np.uint8(MAPL_UNDEF))
            setattr(ncFidOut.variables[VAR_NAMES[l]],'fmissing_value',np.uint8(MAPL_UNDEF))
         
        latout.units  = 'degrees north'
        lonout.units  = 'degrees east'
        
        lonout [:] = np.array([ -180. + DX/2. + DX*i for i in range(NC)],dtype=np.double) * (np.pi /np.double(180.))
        latout [:] = np.array([DY*i + -90. + DY/2. for i in range(NR)],dtype=np.double) * (np.pi /np.double(180.))
        
        ncFidOut.close()  
        del ncFidOut
    
    def nearest_cell (array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    def fill_gaps (data, fill_value=None, ocean=None):
        NX = min (data.shape)
        data = ma.masked_array(data,data==fill_value)
        odata= ma.masked_array(data,data==ocean)

        for zoom in range (1,NX//2):
            for direction in (-1,1): 
                shift = direction * zoom
                if not np.any(data.mask): break
                for axis in (0,1):                
                    a_shifted = np.roll(data ,shift=shift,axis=axis)
                    o_shifted = np.roll(odata,shift=shift,axis=axis)
                    idx=~a_shifted.mask * data.mask*~o_shifted.mask
                    data[idx]=a_shifted[idx]
    
    def regrid_to_coarse (data):
        global NC, NR, IM, JM
        JX = NR // JM
        IX = NC // IM
        temp = data.reshape((data.shape[0] // JX, JX, data.shape[1] // IX, IX))
        return np.nanmean(temp, axis=(1,3))
    
    def get_tag_list(url,label):
        page = requests.get(url)
        soup = BeautifulSoup(page.content, 'html.parser')
        allfiles = soup.find_all(label)
        thislist = [tag.text for tag in allfiles]
        return thislist


class MAP_FROM_LATLON(object):
    
    # ---- stitch all hdf files together and construct global arrays on the sin-grid.
    # ---- employs forward mapping from the lat lon grid using a mask file.        
    
    def __init__ (self, HDF_FILES, DATADIR, OUTFILE, lc_mask, x_index,y_index):
        
        lai_500 = np.full ((NR,NC),np.uint8(MAPL_UNDEF))
        QC_500  = np.full ((NR,NC),np.uint8(MAPL_UNDEF))
        FP_500  = np.full ((NR,NC),np.uint8(MAPL_UNDEF))
        XC_500  = np.full ((NR,NC),np.uint8(MAPL_UNDEF))

        print (OUTFILE)

        for f in range(len(HDF_FILES)):
            FILE_NAME = DATADIR + HDF_FILES[f]
            if os.path.isfile(FILE_NAME):
                H = int(HDF_FILES[f][HDF_FILES[f].find('h')+1:HDF_FILES[f].find('h')+3])
                V = int(HDF_FILES[f][HDF_FILES[f].find('v')+1:HDF_FILES[f].find('v')+3])
                i1 = H * N_MODIS
                j1 = V * N_MODIS
                
                hdf = SD(FILE_NAME, SDC.READ)
                DATAFIELD_NAME = 'Lai_500m'
                data2D = hdf.select(DATAFIELD_NAME)
                data = data2D[:,:].astype(np.int)
                lai_500 [j1: j1 + N_MODIS,i1: i1 + N_MODIS] = data 
                
                hdf = SD(FILE_NAME, SDC.READ)
                DATAFIELD_NAME = 'FparLai_QC'
                data2D = hdf.select(DATAFIELD_NAME)
                data = data2D[:,:].astype(np.int)
                QC_500 [j1: j1 + N_MODIS,i1: i1 + N_MODIS] = data
                
                hdf = SD(FILE_NAME, SDC.READ)
                DATAFIELD_NAME = 'Fpar_500m'
                data2D = hdf.select(DATAFIELD_NAME)
                data = data2D[:,:].astype(np.int)
                FP_500 [j1: j1 + N_MODIS,i1: i1 + N_MODIS] = data 
                
                hdf = SD(FILE_NAME, SDC.READ)
                DATAFIELD_NAME = 'FparExtra_QC'
                data2D = hdf.select(DATAFIELD_NAME)
                data = data2D[:,:].astype(np.int)
                XC_500 [j1: j1 + N_MODIS,i1: i1 + N_MODIS] = data
                
                hdf.end()
                del hdf
                
            else:
                print ('MISSING hdf FILE ')
                print (FILE_NAME)
                shutil.rmtree(CWD + '/download/')
                sys.exit()                   

        DriverFunctions.create_netcdf (OUTFILE,['Lai_500m','FparLai_QC','Fpar_500m','FparExtra_QC'])

        data_high = np.full((NR,NC),np.uint8(MAPL_UNDEF))
        ncFidOut = Dataset(OUTFILE,mode='a')
        data_high.reshape(NC*NR)[lc_mask] = lai_500[y_index,x_index]
        ncFidOut.variables['Lai_500m'  ][:] = data_high

        data_high[:,:] =  np.uint8(MAPL_UNDEF)
        data_high.reshape(NC*NR)[lc_mask] = QC_500[y_index,x_index]
        ncFidOut.variables['FparLai_QC'][:] = data_high
        
        data_high[:,:] =  np.uint8(MAPL_UNDEF)
        data_high.reshape(NC*NR)[lc_mask] = FP_500[y_index,x_index]
        ncFidOut.variables['Fpar_500m'][:] = data_high  
        
        data_high[:,:] =  np.uint8(MAPL_UNDEF)
        data_high.reshape(NC*NR)[lc_mask] = XC_500[y_index,x_index]
        ncFidOut.variables['FparExtra_QC'][:] = data_high            
        
        ncFidOut.close()
        del lai_500
        del QC_500
        del FP_500
        del XC_500
        del ncFidOut
        del data_high

class MAP_FROM_SINEGRID (object):    
    
    # ---- read MODIS LAI 15A2H granule from downloaded .hdf file.
    # ---- employs inverse mapping from sinusoidal grid to lat/lon
    
    def __init__ (self, FILE_NAME):
        
        hdf = SD(FILE_NAME, SDC.READ)
        
        # ---- Query global attributes
        
        sattr = str(hdf.attributes())
        
        xdim = int(re.findall(r"XDim=(\d+)", sattr)[0])
        ydim = int(re.findall(r"YDim=(\d+)", sattr)[0])
        
        earth_radius = np.double( re.findall(r"ProjParams=\((-?\d+\.\d*|\d*\.\d+)", sattr)[0] )
        upleft   = np.double( re.findall(r"UpperLeftPointMtrs=\((-?\d+\.\d*|\d*\.\d+),(-?\d+\.\d*|\d*\.\d+)", sattr)[0] )
        lowright = np.double( re.findall(r"LowerRightMtrs=\((-?\d+\.\d*|\d*\.\d+),(-?\d+\.\d*|\d*\.\d+)", sattr)[0] )            
        # upleft = np.double( re.findall(r"UpperLeftPointMtrs=\((-?\d+\.\d*|\d*\.\d+),-?(\d+\.\d*|\d*\.\d+)", sattr)[0] )
        # lowright = np.double( re.findall(r"LowerRightMtrs=\((-?\d+\.\d*|\d*\.\d+),-?(\d+\.\d*|\d*\.\d+)", sattr)[0] )
        
        # ---- MODIS pixel geolocation X, Y grid
        dxm = (lowright[0] - upleft[0]) / xdim
        dym = (upleft[1] - lowright[1]) / ydim
        
        xm, ym = np.meshgrid(np.arange(1, xdim+1), np.arange(1, ydim+1))
        x_sin = (xm-xm+1)*(upleft[0] + dxm/2.) + (xm-1)*dxm
        y_sin = (ym-ym+1)*(upleft[1] - dym/2.) - (ym-1)*dym
        
        # ---- Create Longitude and Latitude maps from Sinusoidal projection
        lat = (y_sin / earth_radius * 180.0/np.pi)
        lon = (x_sin / (earth_radius * np.cos(y_sin/earth_radius)) * 180.0/np.pi)
        eastern_hemispehere = np.where(lon < -180.)
        western_hemispehere = np.where(lon >  180.)
        lon [eastern_hemispehere] = lon [eastern_hemispehere] + 360.
        lon [western_hemispehere] = lon [western_hemispehere] - 360.
        self.x_index = np.array(np.floor ((lon + 180.)/DXY),dtype=np.int)         
        self.y_index = np.array(np.floor ((lat + 90.)/DXY) ,dtype=np.int)
        #        self.earth_radius = earth_radius
        
        # ---- Read LAI
        
        
        DATAFIELD_NAME = 'Lai_500m'
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:,:].astype(np.int)
        invalid = data == MAPL_UNDEF
        self.lai_500 =  data
        
        DATAFIELD_NAME = 'FparLai_QC'
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[: ,:].astype(np.int)
        invalid = data == MAPL_UNDEF
        self.QC_500 =  data                
        
        DATAFIELD_NAME = 'Fpar_500m'
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[: ,:].astype(np.int)
        invalid = data == MAPL_UNDEF
        self.FP_500 =  data   
        
        DATAFIELD_NAME = 'FparExtra_QC'
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[: ,:].astype(np.int)
        invalid = data == MAPL_UNDEF
        self.XC_500 =  data
        
        hdf.end()    

if FORWARD:

    class GET_ARRAY_INDICES (object):    
    
        # ---- read MODIS LAI 15A2H granule from downloaded .hdf file.
        # ---- employs inverse mapping from sinusoidal grid to lat/lon
    
        def __init__ (self): 

            print ('Initialize array indices')

            #earth_radius = np.double(6371007.181)
            #plane_width  = np.double(1111950.)
            #pixw         = plane_width / N_MODIS
            #Xmin         = -36.*plane_width/2. + pixw/2.
            #Ymax         =  18.*plane_width/2. - pixw/2.
            
            #xlon = np.array([ -180. + DXY/2. + DXY*i for i in range(NC)],dtype=np.double)*(np.pi /np.double(180.))
            #ylat = np.array([DXY*i + -90. + DXY/2. for i in range(NR)],dtype=np.double)*(np.pi /np.double(180.))
            #xlon = np.array ([xlon,]*NR)
            #xlon = earth_radius * xlon * np.cos(ylat.reshape(NR,1))
            #xx_index= np.floor(((xlon - Xmin)/pixw) - 0.5).astype(int)
            #del xlon
            
            #ylat = earth_radius * np.array ([ylat,]*NC).transpose()
            #yy_index= np.floor(((Ymax - ylat)/pixw) - 0.5).astype(int)
            #del ylat   
            
            #lcFile       = Dataset('/discover/nobackup/projects/lis/LS_PARAMETERS/MODIS/MCD15A2H.006/MCD12Q1_2018_LC_Type1.nc4',mode='r')
            #lc           = np.array(lcFile.variables['LC_Type1'][:])
            #self.lc_mask      = np.flatnonzero ((lc >= 1) & (lc <= 20))
            #lcFile.close
            #del lcFile
            #del lc

            #self.x_index = xx_index.reshape(NC*NR)[self.lc_mask]
            #del xx_index
    
            #self.y_index = yy_index.reshape(NC*NR)[self.lc_mask]
            #del yy_index

            #LANDMASK = Dataset(MASKFILE,'w',format='NETCDF4')
            #NT = self.lc_mask.shape[0]
            #print (NT)
            #lnd  = LANDMASK.createDimension('lnd', NT)
            #lcm  = LANDMASK.createVariable('lc_mask', 'u4', ('lnd'))
            #xind = LANDMASK.createVariable('x_index', 'u4', ('lnd'))
            #yind = LANDMASK.createVariable('y_index', 'u4', ('lnd'))
            #lcm[:] = self.lc_mask
            #xind[:]= self.x_index
            #yind[:]= self.y_index
            #LANDMASK.close()  

            LANDMASK = Dataset(MASKFILE)
            self.lc_mask = np.array(LANDMASK.variables['lc_mask'][:])
            self.x_index = np.array(LANDMASK.variables['x_index'][:])
            self.y_index = np.array(LANDMASK.variables['y_index'][:])
            LANDMASK.close()

#----------------------------------------------------------#   
#        BEGIN PROCESSING MODIS 8-day COMPOSITES           #
#----------------------------------------------------------#   

# ---- single MODIS date
if FORWARD:

    FWI = GET_ARRAY_INDICES() 

if YYYYDOY:
    year     = YYYYDOY[0:-3]
    doy      = YYYYDOY[5:] + '/'
    this_doy = int(doy[0:3])
    date1    = datetime(int(year[0:4]),1,1) + timedelta (days=this_doy-1)
    files    = DriverFunctions.get_tag_list(MODIS_PATH + year + doy,"span")[1:-1]
    
    #DATADIR = CWD + '/hdf_download/' + year + '/'
    DATADIR = OUTDIR + date1.strftime('%Y.%m.%d') + '/'
    #os.makedirs(DATADIR)
    #os.chdir(DATADIR)
    #command = 'wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=5 "' + MODIS_DOWN + year + doy +'" --header "Authorization: Bearer D4FAE892-9F8D-11EA-BB77-CA202F439EFB" -P .'
    #os.system(command)
    #os.chdir(CWD)
    #DATADIR = DATADIR + doy

    OUTFILE  = CWD + '/hdf_download/' + '/' + M6_NAME + '.006_LAI_' + year[0:-1] + doy[0:-1] + '.nc4'

    if FORWARD:
 
        thisset = MAP_FROM_LATLON (files, DATADIR, OUTFILE, FWI.lc_mask, FWI.x_index,FWI.y_index)
        del thisset         
    else:

        lai_high = np.full((NR,NC),np.uint8(MAPL_UNDEF))
        qc_high  = np.full((NR,NC),np.uint8(MAPL_UNDEF))
        fp_high = np.full((NR,NC),np.uint8(MAPL_UNDEF))
        xc_high = np.full((NR,NC),np.uint8(MAPL_UNDEF))

        for f in range(len(files)):
            FILE_NAME = DATADIR + files[f]
            print('PROCESSING : ', FILE_NAME)
            thistile = MAP_FROM_SINEGRID (FILE_NAME)
            
            lai_array= thistile.lai_500.reshape (N_MODIS*N_MODIS)
            xin_array= thistile.x_index.reshape (N_MODIS*N_MODIS)
            yin_array= thistile.y_index.reshape (N_MODIS*N_MODIS)
            qc_array = thistile.QC_500.reshape (N_MODIS*N_MODIS)
            xc_array = thistile.XC_500.reshape (N_MODIS*N_MODIS)
            fp_array = thistile.FP_500.reshape (N_MODIS*N_MODIS)   
            lai_mask = np.where((lai_array >= 0) & (lai_array <= 100))
            qc_high [yin_array[lai_mask],xin_array[lai_mask]] = qc_array [lai_mask]
            lai_high[yin_array[lai_mask],xin_array[lai_mask]] = lai_array[lai_mask]
            fp_high [yin_array[lai_mask],xin_array[lai_mask]] = fp_array [lai_mask]
            xc_high [yin_array[lai_mask],xin_array[lai_mask]] = xc_array [lai_mask]

        
        DriverFunctions.create_netcdf(OUTFILE,['Lai_500m','FparLai_QC','Fpar_500m','FparExtra_QC'])
        ncFidOut = Dataset(OUTFILE,mode='a')
        LAIOUT = ncFidOut.variables['Lai_500m'  ]
        QCOUT  = ncFidOut.variables['FparLai_QC']
        FPOUT  = ncFidOut.variables['Fpar_500m' ]
        XCOUT  = ncFidOut.variables['FparExtra_QC' ]       
        LAIOUT[:] = lai_high
        QCOUT [:] = qc_high
        XCOUT [:] = xc_high
        FPOUT [:] = fp_high
        
        ncFidOut.close()

#    command = 'ncatted -O -a scale_factor,Lai_500m,o,f,0.1 ' + OUTFILE
#    os.system(command)
#    command = 'ncatted -O -a scale_factor,Fpar_500m,o,f,0.01 ' + OUTFILE
#    os.system(command)
#    FINALFILE      = OUTDIR + '/' + M6_NAME + '.006_LAI_' + year[0:-1] + doy[0:#-1] + '.nc4'
#    command = 'nccopy -d6 ' + OUTFILE + ' ' + FINALFILE
#    os.system(command)
    
    sys.exit()

#  ---- loop through years
    
years = DriverFunctions.get_tag_list(MODIS_PATH,"a")[1:-5][-1:]
date0 = datetime(int(firstdate[0:4]), int(firstdate[4:6]), int(firstdate[6:8]))
alldates = []
os.chdir(CWD)

for year in years: 

    doys  = DriverFunctions.get_tag_list(MODIS_PATH + year,"a")[1:-5]

    # ---- DOY loop
    for idx, doy in enumerate(doys):

        this_doy = int(doy[0:3])
        date1 = datetime(int(year[0:4]),1,1) + timedelta (days=this_doy-1)
        if this_doy < 361:
            date2 = date1 + timedelta (days=8)
        else:
            date2 = datetime(int(year[0:4])+1,1,5)
        mday  = date1 + (date2 - date1)/2

        OUTFILE      = OUTDIR + '/' + M6_NAME + '.006_LAI_' + year[0:-1] + doy[0:-1] + '.nc4'

        if not os.path.isfile(OUTFILE):
                
            files    = DriverFunctions.get_tag_list(MODIS_PATH + year + doy,"span")[1:-1]
            #DATADIR = OUTDIR + date1.strftime('%Y.%m.%d') + '/'
            os.chdir(CWD)
            DATADIR = CWD + '/download/' + year + '/'
            os.makedirs(DATADIR)
            os.chdir(DATADIR)
            command = 'wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=5 "' + MODIS_DOWN + year + doy +'" --header "Authorization: Bearer D4FAE892-9F8D-11EA-BB77-CA202F439EFB" -P .'
            os.system(command)
            os.chdir(CWD)
            DATADIR = DATADIR + doy
                
            # ---- Stitch MODIS granules

            TMPFILE  = CWD + '/download/' + '/' + M6_NAME + '.006_LAI_' + year[0:-1] + doy[0:-1] + '.nc4'                
            if FORWARD:

                thisset = MAP_FROM_LATLON (files, DATADIR, TMPFILE, FWI.lc_mask, FWI.x_index,FWI.y_index)
                
            else:
                
                lai_high = np.full((NR,NC),np.uint8(MAPL_UNDEF))
                qc_high = np.full((NR,NC),np.uint8(MAPL_UNDEF))
                fp_high = np.full((NR,NC),np.uint8(MAPL_UNDEF))
                xc_high = np.full((NR,NC),np.uint8(MAPL_UNDEF))
                
                for f in range(len(files)):

                    FILE_NAME = DATADIR + files[f]
                    print('PROCESSING : ', FILE_NAME)
                    if os.path.isfile(FILE_NAME):
                        thistile = MAP_FROM_SINEGRID (FILE_NAME)
                        
                        lai_array= thistile.lai_500.reshape (N_MODIS*N_MODIS)
                        xin_array= thistile.x_index.reshape (N_MODIS*N_MODIS)
                        yin_array= thistile.y_index.reshape (N_MODIS*N_MODIS)
                        
                        #if REGRID:
                        lai_mask = np.where((lai_array >= 0.) & (lai_array <= 10.))
                        lai_high [yin_array[lai_mask],xin_array[lai_mask]] = lai_array[lai_mask]
                        #    else:
                        qc_array = thistile.QC_500.reshape (N_MODIS*N_MODIS)
                        xc_array = thistile.XC_500.reshape (N_MODIS*N_MODIS)
                        fp_array = thistile.FP_500.reshape (N_MODIS*N_MODIS)                    
                        lai_mask = np.where((lai_array >= 0) & (lai_array <= 100))
                        qc_high [yin_array[lai_mask],xin_array[lai_mask]] = qc_array [lai_mask]
                        lai_high[yin_array[lai_mask],xin_array[lai_mask]] = lai_array[lai_mask]
                        fp_high [yin_array[lai_mask],xin_array[lai_mask]] = fp_array [lai_mask]
                        xc_high [yin_array[lai_mask],xin_array[lai_mask]] = xc_array [lai_mask]
                    else:
                        print ('MISSING hdf FILE ')
                        shutil.rmtree(CWD + '/download/')
                        sys.exit()
                        
                    DriverFunctions.create_netcdf(TMPFILE,['Lai_500m','FparLai_QC','Fpar_500m','FparExtra_QC'])
                    ncFidOut = Dataset(TMPFILE,mode='a')
                    LAIOUT = ncFidOut.variables['Lai_500m'  ]
                    QCOUT  = ncFidOut.variables['FparLai_QC']
                    FPOUT  = ncFidOut.variables['Fpar_500m' ]
                    XCOUT  = ncFidOut.variables['FparExtra_QC' ]
                    LAIOUT[:] = lai_high
                    QCOUT [:] = qc_high
                    XCOUT [:] = xc_high
                    FPOUT [:] = fp_high                    
                    ncFidOut.close()

            command = 'ncatted -O -a scale_factor,Lai_500m,o,f,0.1 ' + TMPFILE
            os.system(command)
            command = 'ncatted -O -a scale_factor,Fpar_500m,o,f,0.01 ' + TMPFILE
            os.system(command)
            FINALFILE      = OUTDIR + '/' + M6_NAME + '.006_LAI_' + year[0:-1] + doy[0:-1] + '.nc4'
            command = 'nccopy -d6 ' + TMPFILE + ' ' + OUTFILE
            os.system(command)            
            shutil.rmtree(CWD + '/download/')
            #os.remove (TMPFILE)
            #sys.exit()
    alldates.extend(doys)  

