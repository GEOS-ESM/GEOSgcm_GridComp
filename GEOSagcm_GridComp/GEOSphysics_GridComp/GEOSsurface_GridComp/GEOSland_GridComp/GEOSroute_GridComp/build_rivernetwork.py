#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 12 16:42:19 2019

@author: smahanam
"""
###############
# Load modules
###############

from netCDF4 import Dataset, stringtochar
import numpy as np
import csv
import os
import datetime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.basemap import Basemap
from statistics import mode
from scipy.io import FortranFile
import math
import sys
import gmplot
        
######################################
# Read streamflow station information
######################################

file = input("Enter station lat/lon information file name with full path \
              \n \t (Note : The .csv file format is as follows:\
              \n \t One line header is followed by a seperate line\
              \n \t for each station, please see the below example.)\
              \n \
              \n STA_LAT, STA_LON, STA_NAME\
              \n 38.709372, -91.407608, Missouri Hermann\
              \n 43.065813, -98.563029, Missouri Ft. Randall Dam\
              \n \
              \n ")
sta_info = []
with open(file) as fh:
    rd = csv.DictReader(fh, delimiter=',')
    for row in rd:
        sta_info.append(row)
N_STA = len(sta_info)
lats = np.zeros(shape=(N_STA)); lons = np.zeros(shape=(N_STA)); snames = ["" for x in range(N_STA)]
for x in range(0,N_STA) : lats[x], lons[x], snames[x] = sta_info[x].values()

################################
# Create basin information file
################################

OUTDIR = '../river_basin_infor/'

if not os.path.exists(OUTDIR):
        os.mkdir(OUTDIR)

MAX_CAT_BASIN = 12500
ncFidOut = Dataset(OUTDIR + '/Riverflow_Station_Information.nc4', mode='w', clobber='True',format='NETCDF4')
StaDim  = ncFidOut.createDimension('n_sta' , N_STA)
UniDim  = ncFidOut.createDimension('n_catb', MAX_CAT_BASIN)
StrDim  = ncFidOut.createDimension('strl'  , 40)

# Check whether AGCM or GEOSldas
# ------------------------------

til_file = None
if os.path.isfile('LDAS.rc'):
    til_file = os.popen('ls -1 ../output/*/rc_out/*_tilecoord.bin').read().strip()
    trn_file = os.popen("grep BCS_PATH * | cut -d':' -f3").read().strip() + "/clsm/Grid2Catch_TransferData.nc"

    with FortranFile(til_file, 'r') as tf:
        NTILES = tf.read_ints(np.int32)
        tile_id     = np.array(tf.read_ints(np.int32))
        pfaf_domain = np.array(tf.read_ints(np.int32))
        pfaf_domain = np.array(tf.read_ints(np.int32))
        com_lon     = np.array(tf.read_reals(float))
        com_lat     = np.array(tf.read_reals(float))
        min_lon     = np.array(tf.read_reals(float))
        max_lon     = np.array(tf.read_reals(float))
        min_lat     = np.array(tf.read_reals(float))
        max_lat     = np.array(tf.read_reals(float))
        i_indg      = np.array(tf.read_ints(np.int32))
        j_indg      = np.array(tf.read_ints(np.int32))
                                
    cat_index_domain = np.array(np.unique(np.sort(pfaf_domain)))
    IM = i_indg.max() - i_indg.min() + 1
    JM = j_indg.max() - j_indg.min() + 1

    # read Grid2Catch_TransferData.nc
    Grid2Catch    = Dataset(trn_file,mode='r')
    ncats_in_grid = np.array (Grid2Catch.variables['NCats_in_GRID'][:])
    pfaf_index    = np.array (Grid2Catch.variables['Pfaf_Index'][:])
    Grid2Catch.close()
    sys.exit()
else:    
    print ('Enter output AGCM grid resolution 288x181, 576x361 .. etc : ' )
    imjm = input ('IMxJM : ')
    im, jm = imjm.split('x')     
    IM = int(im)
    JM = int(jm)
    
    DX = 360. / IM
    if (JM % 2) == 1:
        DY = 180. / (JM - 1)
    elif(JM % 2) == 0 :
        DY = 180. / JM    
    xdim    = ncFidOut.createDimension('lon' , IM)
    ydim    = ncFidOut.createDimension('lat' , JM)
    
ncFidOut.description = "Pfafstetter catchment indentifications for river basins"
ncFidOut.history     = "Created " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

latsOut = ncFidOut.createVariable('sta_lat','f4',('n_sta',))
#latsOut.long_name, 'Station Latitude'
latsOut [:] = lats
lonsOut = ncFidOut.createVariable('sta_lon','f4',('n_sta',))
#lonsOut.long_name, 'Station Longitude'
lonsOut [:] = lons
NameOut = ncFidOut.createVariable('sta_name','S1',('n_sta','strl'))
#NameOut.long_name, 'Station Name'
#for nrec in range(nrecs):
data = []
data = np.empty((N_STA,),'S'+repr(40))
for n in range(N_STA):
    tstr = snames[n]
    data[n] = tstr[0:len(tstr)]
datac = stringtochar(data)
NameOut[:] = datac
AreaOut = ncFidOut.createVariable('CompBasinArea','f4',('n_sta',))
#AreaOut.long_name, 'Computed Basin Area'
#AreaOut.units,    'km^2'
NcatOut = ncFidOut.createVariable('NCatB','i4',('n_sta',))
#NcatOut.long_name, 'Number of catchments in the basin'
GCIDOut = ncFidOut.createVariable('GlobalID','i4',('n_catb','n_sta'))
#GCIDOut.long_name, 'Catchment Index numbers of catchments within the basin'

if os.path.isfile('LDAS.rc'):
    NcellsBOut = ncFidOut.createVariable('N_cells_basin','i4',('n_sta',))
    SMAPIDOut  = ncFidOut.createVariable('SMAPID','i4',('n_catb','n_sta'))
else:
    maskOut = ncFidOut.createVariable('basin_mask','i4',('lat','lon','n_sta'))
    lonout  = ncFidOut.createVariable('lon','f4',('lon'))
    varr    = np.full (IM, 0.)
    for i in range (IM):
        varr [i] = -180. + DX*i
    lonout [:] = varr
    latout  = ncFidOut.createVariable('lat','f4',('lat'))
    varr    = np.full (JM,0.)
    varr [0]= -90. + DY/4.
    varr [JM-1] = 90. - DY/4.
    for i in range (1, JM -1):
        varr [i] = DY*i + -90. 
    latout [:] = varr 
    
#######################
# Opening PfafCode file
#######################

# (1) Pfafstetter Index raster
SRTMPfaf = Dataset('SRTM_PfafData.nc',mode='r')
glat     = np.array (SRTMPfaf.variables['latitude'][:])
glon     = np.array (SRTMPfaf.variables['longitude'][:])
catid    = np.array (SRTMPfaf.variables['CatchIndex'][:])
DXY = 360./21600. # Pfafstetter index raster resolution
SRTMPfaf.close()

# (2) River Network information
RNetWork =  Dataset ('RiverNetwork_information.nc4',mode='r')
num_catchs       = np.array (RNetWork.variables['NUM_CATCHS'][:])
downstream_lon   = np.array (RNetWork.variables['DownStream_lon'][:])
downstream_lat   = np.array (RNetWork.variables['DownStream_lat'][:])
upstream_lon     = np.array (RNetWork.variables['UpStream_lon'][:])
upstream_lat     = np.array (RNetWork.variables['UpStream_lat'][:])
catchment_index  = np.array (RNetWork.variables['CatchmentIndex'][:])
upstream_index   = np.array (RNetWork.variables['UPSTRIndex'][:])
downstream_index = np.array (RNetWork.variables['DNSTRIndex'][:])
catch_area       = np.array (RNetWork.variables['CATCH_AREA'][:])

RNetWork.close()

#################
# Define classes
#################
    
class RiverBasin(object):

    def find_pfid_at_station (self,sta_lat, sta_lon):
        global catid, DXY, glat, glon
        j_sta = np.asscalar(np.where(np.logical_and(glat >= sta_lat, glat < sta_lat + DXY))[0])
        i_sta = np.asscalar(np.where(np.logical_and(glon >= sta_lon, glon < sta_lon + DXY))[0])
        loc_index = np.asscalar(catid [j_sta, i_sta])
        if loc_index < 1:
            loc_index = np.max(np.array(catid [j_sta-2:j_sta+3, i_sta-2:i_sta+3]))
        return loc_index 
        
    def get_continent_id (self,loc_index):
        if (1 <= loc_index <= 75368):
            cid = 0                         # Asia
        elif (75369 <= loc_index <= 140751):
            cid = 1                         # Africa                        
        elif (140752 <= loc_index <= 189105):
            cid = 2                         # North America
        elif (189106 <= loc_index <= 229074):
            cid = 3                         # Europe
        elif (229075 <= loc_index <= 267083):
            cid = 4                         # South America
        elif (267084 <= loc_index <= 291284):
            cid = 5                         # Australia
        return cid   
    
    def update_array (self, a, n, v):
        if (a.ndim == 1):
            a[n] = v
        elif(a.ndim > 1):
            a[:,n] = v
        return a
    
    def link_upstream_catchments (self, catch_in, DOM_MIN, upst, mouth_coord, cum_area, cat_area, catch_down):
        if (upst[0,catch_in] < 0):
            if (catch_in != catch_down):
                cum_area[catch_down] = cum_area[catch_down] + cat_area [catch_in]
            return mouth_coord,cum_area
        for i in range(30):
            if (upst[i,catch_in] > 0):
                this_catch = upst[i,catch_in] - DOM_MIN
                mouth_coord[:,this_catch] = mouth_coord[:,catch_in]
                if (catch_in != catch_down):
                    cum_area[this_catch] = cum_area[this_catch] + cat_area [this_catch]
                    cum_area[catch_in]   = cum_area[catch_in]   + cat_area [this_catch]
                    cum_area[catch_down] = cum_area[catch_down] + cat_area [this_catch]
                if(upst[0,this_catch] > 0):
                    mouth_coord,cum_area = self.link_upstream_catchments (this_catch, DOM_MIN, upst, mouth_coord,
                                                              cum_area, cat_area, catch_down)
        return mouth_coord,cum_area
      
    def derive_basin_mask (catid_mask):
        global catid, DXY, glat, glon, IM, JM, DX, DY
        basin_mask    = np.full ((JM, IM), 0)
        dx2 = math.ceil(DX / DXY /2.)
        dy2 = math.ceil(DY / DXY /2.)
        for j in range (1, JM - 1):
            sta_lat = DY*j + -90.
            j_sta = (np.where(np.logical_and(glat >= sta_lat, glat < sta_lat + 2*DXY))[0]).min()
            for i in range (1, IM):
                sta_lon = DX*i -180.
                i_sta = (np.where(np.logical_and(glon >= sta_lon, glon < sta_lon + 2*DXY))[0]).min()
                catid_domain = catid_mask [j_sta - dy2: j_sta + dy2+1, i_sta - dx2: i_sta + dx2+1]
                if any(np.reshape(catid_domain,((2*dx2+1)*(2*dy2+1)))):
                    basin_mask [j,i] = 1
        return basin_mask

    def __init__(self, sta_lat, sta_lon):        
        global num_catchs, catchment_index, catch_area, upstream_index 
        
        pfid_at_station      = self.find_pfid_at_station (sta_lat, sta_lon)
        self.continent_id    = self.get_continent_id (pfid_at_station)
        
        ncat_continent       = num_catchs[self.continent_id]
        domain_index         = np.array (catchment_index[self.continent_id,0:ncat_continent])
        catch_area_continent = np.array (catch_area     [self.continent_id,0:ncat_continent])
        upst                 = np.array (upstream_index [self.continent_id,:,0:ncat_continent])
        
        mouth_coord          = np.full ((2,ncat_continent),-9999.)
        cum_area             = np.full (ncat_continent,0.)
        this_catch           = pfid_at_station - domain_index.min()

        DOM_MIN              = domain_index.min()
        mouth_coord          = self.update_array (mouth_coord, this_catch, [sta_lat, sta_lon])
        cum_area             = self.update_array (cum_area, this_catch, catch_area_continent [this_catch])
         
        mouth_coord,cum_area = self.link_upstream_catchments (this_catch, DOM_MIN, upst, mouth_coord, \
                                                              cum_area, catch_area_continent, this_catch)
        
        self.upst_catchs     = domain_index[(np.where(mouth_coord[0,:] != -9999.))]
        self.comp_basin_area = cum_area[this_catch]
        
#######################
# Loop through stations
#######################

basin_areas    = []
catch_in_basin = np.full ((MAX_CAT_BASIN,N_STA),-9999.)
ncatch_in_basin= []
continent_id   = []
basin_masks    = np.full ((JM, IM, N_STA), 0)
    
for n in range(N_STA):
    this_basin = RiverBasin(lats[n], lons[n])        
    basin_areas.append(this_basin.comp_basin_area)
    ncatch_in_basin.append(np.size(this_basin.upst_catchs))
    continent_id.append(this_basin.continent_id)
    catch_in_basin [0:ncatch_in_basin[n], n] = this_basin.upst_catchs    
    if os.path.isfile('LDAS.rc'):

        catid_mask = np.in1d(catid,[catch_in_basin[0:ncatch_in_basin[n], n]]).reshape(catid.shape)        
    else:       
        catid_mask = np.in1d(catid,[catch_in_basin[0:ncatch_in_basin[n], n]]).reshape(catid.shape)
        basin_masks [:,:,n] = RiverBasin.derive_basin_mask (catid_mask)
    del this_basin
    
AreaOut [:] = basin_areas
NcatOut [:] = ncatch_in_basin
GCIDOut [:] = catch_in_basin
maskOut [:] = basin_masks

ncFidOut.close()

##############################################
#                PLOTTING 
##############################################   
            
class BasinMaps(object):  
    def draw_basinmaps_html (name, lats_down,lats_up,lons_down,lons_up):
        global OUTDIR
        gmap = gmplot.GoogleMapPlotter(np.array([lats_down,lats_up]).mean(),np.array([lons_down,lons_up]).mean(), 5) 
        for i in range(lats_down.size):
            gmap.scatter([lats_down[i], lats_up[i]], [lons_down[i], lons_up [i]], '# FF0000', size = 40, marker = False)
            gmap.plot([lats_down[i], lats_up[i]], [lons_down[i]],'cornflowerblue', edge_width = 2.5)
        gmap.draw( OUTDIR + name + ".html" )
    
    def plot_basinmaps (sta_name, boundary,lats_down,lats_up,lons_down,lons_up, mask = None):
        global IM, JM, DX, DY
        m = Basemap(projection = 'cyl', llcrnrlat= boundary[0] - 0.25,urcrnrlat=boundary[2]+ 0.25, llcrnrlon=boundary[1]- 0.25,urcrnrlon=boundary[3] + 0.25, resolution ='c')
        m.drawcoastlines()
        m.fillcontinents(color='white')
        m.drawmapboundary(fill_color='white')
        m.drawstates(color='black')
        m.drawcountries(color='black')
        for i in range(lats_down.size):
            x, y = m([lons_down[i],lons_up[i]], [lats_down[i], lats_up[i]])
            m.plot(x, y, linewidth=1, color = 'cornflowerblue') 
        if mask.any():
            for j in range (1, JM):
                sta_lat = DY*j + -90. 
                for i in range (1, IM):
                    sta_lon = DX*i -180.    
                    if mask [j,i] == 1:
                        xx = sta_lon + np.array([-1,-1, 1, 1,-1]) * DX/2.
                        yy = sta_lat + np.array([-1, 1, 1,-1,-1]) * DY/2.
#                        m.plot(xx, yy, linewidth=1, color = 'r') 
                        m.plot(xx.mean(), yy.mean(), 'ro', markersize=1)
        plt.title (sta_name)
            
# 1) html table with google maps
################################
    
outfile = open(OUTDIR + "index.html", "w")
outfile.write(""""<html>"
<head>
 <title>Maps of River Basins .. </title>
</head>
<body>
<table border="1">""")

outfile.write ('<tr><th>Name</th><th>Map</th></tr>')

with PdfPages(OUTDIR + 'basin_maps.pdf') as pdf:    
    for n in range(N_STA):
        sname_strip = snames[n].replace(' ','')
        domain_index         = np.array (catchment_index[continent_id[n],0:num_catchs[continent_id[n]]])
        DOM_MIN              = domain_index.min()
        upst_catid           = np.array(catch_in_basin [0:ncatch_in_basin[n], n] - DOM_MIN).astype(int)        
        lats_up              = np.array(upstream_lat[continent_id[n],upst_catid])
        lats_down            = np.array(downstream_lat[continent_id[n],upst_catid])
        lons_up              = np.array(upstream_lon[continent_id[n],upst_catid])
        lons_down            = np.array(downstream_lon[continent_id[n],upst_catid])
        boundary             = []
        boundary.append (np.array([lats_up,lats_down]).min())
        boundary.append (np.array([lons_up,lons_down]).min())
        boundary.append (np.array([lats_up,lats_down]).max())
        boundary.append (np.array([lons_up,lons_down]).max())   
        
        # 1) draw google maps
        BasinMaps.draw_basinmaps_html (sname_strip,lats_down,lats_up,lons_down,lons_up)
        outfile.write ("<tr><td>%s</td><td>%s</td></tr>" % (snames[n], '<a href="' + sname_strip + '.html">' + sname_strip + '</a>'))
        
        # 2) plot maps
        BasinMaps.plot_basinmaps(snames[n], boundary,lats_down,lats_up,lons_down,lons_up, mask = basin_masks[:,:,n])      
        pdf.savefig()
        plt.close()

outfile.write( """</table>
</body></html>""")
outfile.close()


 
    
    






    
    
    
    



