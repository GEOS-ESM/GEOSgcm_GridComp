#!/bin/csh

####################################
##  Setting environment variables ##
####################################

setenv gfile `head -1 clsm/mkCatchParam.log | cut -d 'g' -f2 | cut -d '.' -f1`
setenv workdir `pwd`
setenv NC `head -1 clsm/mkCatchParam.log | cut -d'x' -f2 | cut -d'-' -f1`
setenv NR `head -1 clsm/mkCatchParam.log | cut -d'y' -f2 | cut -d'-' -f1`
setenv ease `head -1 clsm/mkCatchParam.log | cut -d 'e' -f2 | cut -d '-' -f1`

##########################
##  Writing README file ##
##########################

### Start dates for various climatological data
###############################################

set GSWP2_DATES = (1201 0101 0201 0301 0401 0501 0601 0701 0801 0901 1001 1101 1201 0101)
set GEOLAND2_DATES="`printf '1229 0108 0119 0129 0208 0216 0226 0308 0319 0329 0408 0418 0428 0508 0519 0529\\n \
            0608 0618 0628 0708 0719 0729 0808 0819 0829 0908 0918 0928 1008 1019 1029 1108\\n \
            1118 1128 1208 1219 1229 0108'`"
set MODIS8_DATES="`printf '1227 0101 0109 0117 0125 0202 0210 0218 0226 0306 0314 0322 0330 0407 0415 0423\\n \
            0501 0509 0517 0525 0602 0610 0618 0626 0704 0712 0720 0728 0805 0813 0821 0829\\n \
            0906 0914 0922 0930 1008 1016 1024 1101 1109 1117 1125 1203 1211 1219 1227 0101'`"
set MODIS16_DATES="`printf '1219 0101 0117 0202 0218 0306 0322 0407 0423 0509 0525 0610 0626 0712 0728 0813\\n \
            0829 0914 0930 1016 1101 1117 1203 1219 0101'`"
set NDVI_DATES="`printf '1216 0101 0116 0201 0216 0301 0316 0401 0416 0501 0516 0601 0616\\n \
            0701 0716 0801 0816 0901 0916 1001 1016 1101 1116 1201 1216 0101'`"


## Gather information about the data set
########################################
set today=`date +%Y-%m-%d`
set myusage=`head -1 clsm/mkCatchParam.log | tail -1`
set f25tag=`head -2 clsm/mkCatchParam.log | tail -1 | cut -d':' -f2`
set NTILES=`head -1 clsm/catchment.def | tail -1`
set NGLOBAL=`head -1 til/${gfile}.til | cut -c1-12`
set mygrid=`echo $myusage | cut -d'g' -f2 | cut -d '-' -f1`
cvs status src/mkCatchParam.F90 > clsm/TagInfo
#echo GMU_OCT10_SM >> clsm/TagInfo
echo `head -7 clsm/TagInfo | tail -1` | cut -d':' -f2 | cut -d'(' -f1 > clsm/TagName
set MYTAG=`head -1 clsm/TagName | tail -1`
/bin/rm clsm/Tag*

# Set Mask/Topo speifics
########################
set MYMASK=`head -6 clsm/mkCatchParam.log | tail -1`
set NPfafs=291284

set toc_rout="`printf '\\n7. GLOBAL RUNOFF ROUTING MODEL DATA ............................................... 7\\n \
 7.1 Data generation and processing chain\\n \
 7.2 Data files and images\\n \
      7.2.1 Pafafstetter catchment connectivity, channel information\\n \
 7.3 References'`"   
      
      set sec3_veg_des="`printf '(a) Deriving Mosaic vegetation classes:\\n \
\\n \
       Global ESA land cover classification at 10-arcsec resolution is available from\\n \
       European Space Agency (see GLOBCOVER 2009). A simple one-to-one mapping scheme\\n \
       was employed to convert ESA types, at each 10-arcsecond pixel, into Mosaic types.\\n \
       The Mosaic types used for catchment surface elements are determined by computing\\n \
       the dominant Mosaic type of all 10-arcsec pixels within the catchment in question.\\n \
\\n \
       (b) Deriving Catchment-carbon (Catchment-CN) classes:\\n \
\\n \
       The Common Land Model version 4 (CLM4: Oleson et al., 2010) utilizes 17 vegetation\\n \
       classes and the version 4.5 (CLM4.5: Oleson et al., 2013) uses 25 vegetation classes (Table 2).\\n \
       Global arrays of fractional coverages of each of the 17 classes for CLM4 \\n \
       on a 1152x768 grid, and 25 classes of for CLM4.5 on a 7200x3600 grid, \\n \
       as used by CLM4/CLM4.5 was obtained from UCAR.\\n \
\\n \
       The Catchment-CN model’s vegetation classification is based on that used by CLM4, with\\n \
       fractional coverage of each type derived from the Global ESA land cover classification.\\n \
       Each ESA type was mapped into one or more of the CLM4 types. For instance, ESA crop \\n \
       types were mapped one-to-one into CLM4’s crop type, whereas ESA mosaic vegetation was\\n \
       fractionally mapped into region-appropriate grass, shrub, and forest types. This fractional\\n \
       mapping was based on the relative fractions found in the original CLM4/CLM4.5 gridded\\n \
       arrays and on latitude (for differentiating certain types, such as Arctic c3 grass).  \\n \
       Bare soil from the ESA land cover classification is mapped into the broadleaf deciduous \\n \
       shrub type, since bare soil is not an allowed type in our implementation.\\n \
\\n \
       For Catchment-CN, the stress deciduous types (crop and temperate shrubs/grass) utilized \\n \
       by CLM4 is replaced by a mix of two sub-types, one that is seasonally deciduous (with a \\n \
       daylight trigger) and one that is not. Crop type has been further classified to 10 \\n \
       different types in CLM4.5, thus they were not sub-divided into further sub-types, however.  \\n \
       Both sub-types are subject to moisture stress triggers\\n \
       but not to temperature (freezing) stress triggers. The removal of the temperature stress \\n \
       trigger eliminated unnatural swings in leaf carbon during brief temperature stress senescence\\n \
       (Koster et al., 2014).  The relative fractions of the two sub-types applied vary linearly \\n \
       with latitude between 32°-42° in both hemispheres, with 100 percent of the stress deciduous type \\n \
       being replaced by the seasonally deciduous sub-type at 42° and 100 percent replaced by the \\n \
       non-seasonally deciduous type at 32°.\\n \
\\n \
       Within each surface element, the two dominant types are identified; the remaining types are\\n \
       ignored. (Note, however, that in the latitude band 32°-42°, the presence of two sub-types \\n \
       for a potentially dominant type implies that up to four sub-types may be followed.)  The \\n \
       vegetation fractions for the two dominant types are scaled so that they sum to one.\\n \
\\n \
       (c) Deriving Nitrogen Deposition, Annual mean 2m Air Temperature. Soil background albedo \\n \
           for the Catchment-CN Model\\n \
\\n \
       Nitrogen deposition data used in CLM4 (Oleson et al., 2010) were obtained from UCAR. The \\n \
       approximately quarter degree CLM4 deposition data (sum of NHx and NOx species) were spatially\\n \
       interpolated to catchment surface elements and saved (with units of ng m-2 s-1) in  the \\n \
       CLM_NDep_SoilAlb_T2m file (plots/CLM_Ndep_T2m.jpg, top panel). Note that the units of \\n \
       the input deposition data in the Catchment-CN model is g m-2 s-1, and thus the values in \\n \
       the file must be scaled by 1.e-9 before using inside the model. \\n \
\\n \
       The Catchment-CN phenology routine computes a growing-degree-day summation to determine when\\n \
       seasonal and stress deciduous plant functional types (PFT) become active. The critical value \\n \
       for this summation is a function of annual mean 2m air temperature. CLM4 uses the previous \\n \
       years annual mean 2m air temperature (K). In the Catchment-CN model, however, we instead use \\n \
       the climatological annual mean 2m air temperature. Two sources of multi-year global \\n \
       hydrometeorological forcing were used to compute climatological annual mean 2m air temperature\\n \
       separately: 1) GMAO MERRA-2, hourly, 0.625ox0.5o resolution, data from the period 1980-2014, \\n \
       and 2) Sheffield et al (2006), 3-hourly, 1-degree data from the period 1948-2012. The \\n \
       climatological means were computed on the native forcing grids (see middle and bottom panels \\n \
       in plots/CLM_Ndep_T2m.jpg), then spatially interpolated to catchment surface elements, and \\n \
       finally saved in the file CLM_NDep_SoilAlb_T2m.\\n \
\\n \
       In areas of low LAI (less than one), surface albedo in the Catchment-CN model is prescribed\\n \
       using MODIS soil background albedo. For LAI=0, the soil background albedo has full weight; \\n \
       at LAI=1, the albedo is fully determined by mapped catchment vegetation type (using subroutine\\n \
       SIBALB). A linear ramp is used to weight between the MODIS soil background albedo and the \\n \
       SIBALB albedo for 0<LAI<1. There are four albedo components: visible direct (VISDR), visible\\n \
       diffuse (VISDF), near-infrared direct (NIRDR), and near-infrared diffuse (NIRDF). Global \\n \
       fields of soil background VISDR, VISDF, NIRDR and NIRDF at 3-arcmin resolution were obtained\\n \
       from Houldcroft et al. (2009) and spatially aggregated to catchment surface elements for the\\n \
       Catchment-CN model (plots/SoilAlb.jpg).'`"  

set sec3_veg_cite="`printf 'GLOBCOVER 2009 : Products Description and Validation Report, (2011), ESA \\n \
           Technical Note. ,(Available http://due.esrin.esa.int/files/GLOBCOVER2009_\\n \
           Validation_Report_2.2.pdf, last accessed July 2015.)\\n \
       Oleson, K. W., and Co-authors, (2010): Technical description of version 4.0 of\\n \
          the Community Land Model (CLM).  NCAR Technical Note NCAR/TN-478+STR., National\\n \
          Center for Atmospheric Research, P. O. Box 3000, Boulder, Colorado, 80307-3000.\\n \
       Oleson, K. W., and Co-authors, (2013): Technical description of version 4.5 of\\n \
          the Community Land Model (CLM).  NCAR Technical Note NCAR/TN-503+STR., National\\n \
          Center for Atmospheric Research, P. O. Box 3000, Boulder, Colorado, 80307-3000.\\n \
       Koster, R.D., G. K. Walker, G. J. Collatz, and P. E. Thornton, 2014: Hydroclimatic\\n \
          Controls on the Means and Variability of Vegetation Phenology and Carbon Uptake.\\n \
          J. Climate, 27, 56325652. doi: http://dx.doi.org/10.1175/JCLI-D-13-00477.1\\n \
       Houldcroft, Caroline J., William M. F. Grey, Mike Barnsley, Christopher M. Taylor,\\n \
          Sietse O. Los, and Peter R. J. North, 2009: New Vegetation Albedo Parameters and\\n \
          Global Fields of Soil Background Albedo Derived from MODIS for Use in a Climate \\n \
          Model. J. Hydrometeor, 10, 183198. , doi: http://dx.doi.org/10.1175/2008JHM1021.1\\n \
       Sheffield, J., G. Goteti, and E. F. Wood, 2006: Development of a 50-yr high-resolution\\n \
          global dataset of meteorological forcings for land surface modeling, J. Climate, 19\\n \
          (13), 3088-3111\\n \
       Simard, M., N. Pinto, J. B. Fisher, and A. Baccini (2011), Mapping forest canopy height\\n \
          globally with spaceborne lidar, J. Geophys. Res., 116, G04021, doi:10.1029/2011JG001708'`"

if( $MYMASK == GEOS5_10arcsec_mask | $MYMASK == GEOS5_10arcsec_mask.nc | $MYMASK == GEOS5_10arcsec_mask_freshwater-lakes.nc ) then
else
   set MYMASK=global.cat_id.catch.ORIG.DL
   set NPfafs=36716
   set toc_rout=
   set sec3_veg_des="`printf 'Global SiB2 land cover classification at 30-arcsec resolution is avalable from\\n \
       USGS Global Land Cover Characteristics Data Base Version 2.0 (see USGS). A\\n \
       global data set of SiB2 classification at 2.5-arcmin resolution was derived by\\n \
       using dominant types in 30-arcsec USGS data. A simple one-to-one mapping scheme\\n \
       was employed to convert SiB2 types in 2.5-arcmin to Mosaic types and to determine\\n \
       Mosaic types on catchment-tiles.'`"
   set sec3_veg_cite="`printf 'USGS Global Land Cover Characteristics Data Base Version 2.0 \\n \
	    http://edc2.usgs.gov/glcc/globdoc2_0.php'`"
endif

# Set AGCM/SMAP specifics
#########################
set WGRID=AGCM
set int_str1="`printf 'by overlaying the atmospheric grid on ${NPfafs} number of hydraulic catchments \\n in ${MYMASK} mask file.'`"
set sec2_til="`printf 'area, longitude, latitude, ig, jg, cell_frac, integer,   & \\n \
                    pfaf_code, pfaf_index, pfaf_frac'`"
set pfafin_des="`printf 'catchment index (1-$NPfafs) after sorting Pfafstetter codes in ascending order'`"  
set pfaf_des="`printf 'Pfafstetter code of the hydrologic catchment'`"
if( $MYMASK == GEOS5_10arcsec_mask  | $MYMASK == GEOS5_10arcsec_mask.nc | $MYMASK == GEOS5_10arcsec_mask_freshwater-lakes.nc ) set pfaf_des=`echo "${pfafin_des}"`
set pfaf_dest=`echo "${pfaf_des}"`
set sec2_til2="`printf ' (9)    area      [x EarthRadius^2 km2]  tile area\\n\
        (10)   pfaf_frac [-]      fraction of the pfafstetter catchment\\n '`" 
set rout_smap

if(`echo $gfile | cut -d '_' -f1` == SMAP | $ease == EASE) then
   set WGRID=SMAP
   set sec2_til="`printf 'pfaf_code, longitude, latitude, ig, jg, cell_frac, pfaf_index'`"
   if( $MYMASK == GEOS5_10arcsec_mask  | $MYMASK == GEOS5_10arcsec_mask.nc | $MYMASK == GEOS5_10arcsec_mask_freshwater-lakes.nc ) then
    set sec2_til="`printf 'pfaf_index, longitude, latitude, ig, jg, cell_frac, pfaf_index, pfaf_code'`" 
    set pfaf_des=`echo "${pfafin_des}"`
set toc_rout="`printf '\\n7. GLOBAL RUNOFF ROUTING MODEL DATA ............................................... 7\\n \
 7.1 Data generation and processing chain\\n \
 7.2 Data files and images\\n \
      7.2.1 Pafafstetter catchment connectivity, channel information\\n \
      7.2.2 Fractional areas to aggregate from SMAP grid cells to Pfafstetter watersheds\\n \
      7.3 References'`"   

set rout_smap="`printf '\\n \
     7.2.2 Fractional areas to aggregate from SMAP grid cells to Pfafstetter watersheds\\n \
\\n \
       The self-describing SMAP-Catch_TransferData.nc file contains the information needed to \\n \
       convert between SMAP grid cells and watersheds. The information includes the number of \\n \
       Pfafstetter watersheds contributing to the SMAP grid cell and the fractional areas of \\n \
       those contributing watersheds.\\n \
\\n \
       file name : SMAP-Catch_TransferData.nc\\n \
\\n \
       where :\\n \
       NCats_in_SMAP    No. of pfaf catchments contributing to the SMAP cell\\n \
       Pfaf_Index       Pfaf indices (1-291,284) of those contributing catchments\\n \
       Pfaf_Area[km2]   Area of the Pfaf Catchment fraction\\n \
\\n'`" 
 
   endif
   set pfaf_dest="`printf 'Pfafstetter code of the hydrologic catchment'`"
   set int_str1="`printf 'using the land-ice-lakes-ocean mask in ${MYMASK} mask file.'`"
   set sec2_til2=
endif

# Set LAI specifics
###################
set mylai=`head -3 clsm/mkCatchParam.log | tail -1`
if($mylai == GSWP2 | $mylai == GSWPH) then
   set MYLAIDATES="${GSWP2_DATES}"
   set sec4_lai="`printf 'Monthly climatologies of GSWP-2 GrnFrac and LAI were computed by averaging over the\\n \
       17-year period. The computed GrnFrac and LAI climatolgical data wwere spatially\\n \
       interpolated on to catchment-tiles to derive monthly GrnFrac and LAY climaotologies\\n \
       on catchment-tile space.'`"
   set sec4_geo_cite=
endif
if($mylai == GEOLAND2) then
   set MYLAIDATES="${GEOLAND2_DATES}"
endif
if($mylai == MODIS | $mylai == MODGEO) then
   set MYLAIDATES="${MODIS8_DATES}"
	set sec4_lai="`printf 'The Second Global Soil Wetness Project (GSWP-2: Dirmeyer and Oki, 2002)\\n \
       provided monthly Leaf Area Index (LAI) and Greenness Fraction (GrnFrac) data\\n \
       on a 1°×1° grid for the period 1982-1998. A monthly climatology of GrnFrac was\\n \
       computed from these data by temporally averaging over the 17-year period \\n \
       (by month) on the 1°×1° grid and then spatially interpolating the averages \\n \
       onto 30-arcsec pixels. The interpolated GrnFrac data were aggregated over \\n \
       the pixels of each land element to derive a monthly GrnFrac climatology for \\n \
       that land element. \\n \
\\n \
       Global, 10-day averaged LAI data on a 40320×20160 grid for the period \\n \
       1999-2011 are available from GEOLAND2 (Baret et al., 2012 and Camacho et el.\\n \
       2013). In addition, 8-day composites of MOD15A2 v005 MODIS LAI data (MODIS,\\n \
       2008) are available at 30-arcsec (43200×21600) for the period 2000-2013. \\n \
       Preprocessing of the two datasets showed that each had potential flaws, \\n \
       with GEOLAND2 showing questionable seasonal cycles in Siberia, and MODIS \\n \
       showing questionable values over the rain forests. We thus decided to \\n \
       produce a merged LAI data product for GEOS5 to avoid these potential \\n \
       deficiencies. \\n \
\\n \
       The first step in generating the merged product was computing a 10-day \\n \
       climatology of GEOLAND LAI at each GEOLAND2 pixel from the 13 years of \\n \
       GEOLAND2 data and then spatially aggregating the pixel-based climatologies \\n \
       to surface elements. The next step was computing the corresponding 8-day \\n \
       MODIS-based LAI climatology for each surface element from the 14-years of \\n \
       MODIS data. To do this, we used MODIS auxiliary data on surface type to \\n \
       fill in LAI values at certain land pixels as follows: barren, rock, and \\n \
       desert surface types were given LAI values of 0.01, urban-built areas were \\n \
       given values of 0.5; and marshland areas were given values obtained from \\n \
       nearest neighbor pixels.  The filled-in MODIS data had to be projected \\n \
       from a sinusoidal grid to a regular latitude/longitude grid prior to the \\n \
       calculation of the pixel-based 8-day climatologies and the subsequent \\n \
       spatial aggregation of the pixel climatologies to surface element climatologies.\\n \
\\n \
       Note that in high latitudes, both MODIS and GEOLAND2 data are not available\\n \
       over large areal swaths. To address this, we constructed, at every time slice \\n \
       (8 days for MODIS and 10 days for GEOLAND2), a 1°×1° global gridded LAI dataset \\n \
       by spatially aggregating the finer resolution LAI climatological data. Missing \\n \
       LAI values in the finer resolution datasets were filled with the value for the \\n \
       nearest neighbor on the 1°×1° global grid. \\n \
\\n \
       To merge the GEOLAND2 and MODIS data into a single product of 8-day LAI \\n \
       climatologies, the 10-day GEOLAND2 climatological seasonal cycle at each element \\n \
       was first time-interpolated to the 8-day cycle used by MODIS.  MODIS LAI data \\n \
       were then selected for all surface elements except those in South America, Africa \\n \
       and Australia, which took their values from the GEOLAND2 LAI dataset. plots/lai.jpg \\n \
       shows global maps of climatological mean monthly of LAI for the merged LAI data.'`"

      set sec4_geo_cite="`printf 'Baret, F., M. Weiss1, R. Lacaze, F. Camacho, H. Makhmara, P. Pacholcyzk and B. \\n \
          Smets (2012): GEOV1: LAI, FAPAR Essential Climate Variables and FCOVER global\\n \
          time series capitalizing over existing products. Part1: Principles of development\\n \
	  and production. Remote Sensing Environment, 137, 299-309 doi:10.1016/j.rse.2012.12.027\\n \
       Camacho, F., J. Cernicharo, R. Lacaze F. Baret, M. Weiss (2013): GEOV1: LAI, FAPAR\\n \
          Essential Climate Variables and FCover global time series capitalizing over \\n \
          existing products. Part 2: Validation and inter-comparison with reference products\\n \
          Remote Sensing of Enviroment, 137, 310329, doi:10.1016/j.rse.2013.02.030\\n \
       MODIS (2008), MOD15A2 v005, NASA EOSDIS Land Processes DAAC, USGS Earth Resources \\n \
          Observation and Science (EROS) Center, Sioux Falls, South Dakota (https://lpdaac.usgs.gov)\\n \
          , accessed 07-28-2015, at http://e4ftl01.cr.usgs.gov/MOLT/MOD15A2.005/. '`"
endif

# Set Albedo specifics
######################
set myalb=`head -4 clsm/mkCatchParam.log | tail -1`
if($myalb == MODIS1) then
   set AlbFileNames=AlbMap.WS.16-day.tile.0.3_0.7.dat/AlbMap.WS.16-day.tile.0.7_5.0.dat
   set MYALBDATES="${MODIS16_DATES}"
set sec5_mod="`printf 'First, 1-arcmin 16-day composites of MODIS diffused visible (VISDF) and diffused\\n \
       infra-red data (NIRDF) from the period 2000-2004 were processed and the 16-day climatology\\n \
       was computed on each 1-arcmin pixel. Then, on each catchment-tile, MODIS albedo\\n \
       climatology was computed by aggregating 1-arcmin climatological data.'`"
endif
if($myalb == MODIS2) then
   set AlbFileNames=AlbMap.WS.8-day.tile.0.3_0.7.dat/AlbMap.WS.8-day.tile.0.7_5.0.dat
   set MYALBDATES="${MODIS8_DATES}"
set sec5_mod="`printf 'To compute the scaling factors, 30-arcsec 8-day composites of MODIS (MCD43GF, 2014 \\n \
       and Gao et al., 2014) diffuse visible (VISDF) and diffuse near-infrared data (NIRDF) \\n \
       from the period 2001-2011 were temporally averaged into an 8-day climatology. These \\n \
       30-arcsec climatological values were then spatially averaged over a given land surface \\n \
       element’s pixels to produce an 8-day climatology for the land element as a whole.'`"
endif

# Set soil specifics
####################
set mysoil=`head -5 clsm/mkCatchParam.log | tail -1`

## now writing
##############

/bin/rm clsm/README

cat << _EOI_ > clsm/intro

=====================================================================================
||                                                                                 || 
||                        Land Boundary Conditions for the                         ||
||                                                                                 ||
||            Goddard Earth Observing System Model  Version 5 (GEOS-5)             ||
||										   ||
                          ${mygrid} Grid
||                                                                                 ||
||                             ----------------------                              ||
||                                                                                 ||
||                             Data File Descriptions			           ||
||										   ||
||										   ||
||										   ||
||         									   ||
||										   ||
|| _______________________________________________________________________________ ||
||										   ||
||                                         Global Modeling and Assimilation Office ||
||                                                                 610.1 NASA/GSFC ||
||										   ||
||										   ||
|| Author        : Sarith Mahanama (sarith.p.mahanama@nasa.gov)			   ||
||										   ||
|| Data Citation : Mahanama, S.P., R.D. Koster, G.K. Walker, L. Takacs, R.H.       ||
||                 Reichle, G. de Lannoy, Q. Liu, B. Zhao, and M. Suarez (2015) :  ||
||		   Land Boundary Conditions for the Goddard Earth Observing System ||
||		   Model Version 5 (GEOS-5) Climate Modeling System -  Recent      ||
||		   Updates and Data File Descriptions. NASA Technical Report Series||
||		   on Global Modeling and Data Assimilation 104606, v39, 51pp.	   ||
||                 URL: http://gmao.gsfc.nasa.gov/pubs/tm/                         ||
||										   ||
|| Date          : ${today}                                                      ||
||										   ||
=====================================================================================

                                  TABLE OF CONTENTS
                                  ~~~~~~~~~~~~~~~~~

1. INTRODUCTION ................................................................... 1

2. TOPOGRAPHY AND SOIL DATA ....................................................... 2
   2.1 Data generation and processing chain
   2.2 Data files and images
	2.2.1 Tile type, location, area and Pfafstetter catchment mapping
	2.2.2 Western, eastern, southern, northern edges and mean elevation of tiles
	2.2.3 Tile topography - statistics of Compound Topographic Index (CTI) 
	2.2.4 Soil Parameters
   2.3 References

3. VEGETATION CLASSIFICATION DATA ................................................. 3
   3.1 Data generation and processing chain
   3.2 Data files and images
	3.2.1 Mosaic vegetation types and fractions
        3.2.2 vegdyn input data (mosaic primary type, canopy height, and roughness) for GEOS5
	3.2.3 CLM/CLM4.5 and CLM/CLM4.5-carbon vegetation types and fractions 
	3.2.4 CLM Nitrogen Deposition, annual mean T2m, soil back ground albedo
	3.2.5 CLM4.5 ABM, PEATF, GDP, HDM, and soil field capacity
	3.2.6 CLM4.5 lightening frequency climatology
   3.3 References

4. VEGETATION DYNAMIC DATA ........................................................ 4
   4.1 Data generation and processing chain
   4.2 Data files and movies
	4.2.1 Greenness Fraction [-] 
	4.2.2 Leaf Area Index (LAI) [m2/m2]
	4.2.3 Normalized Difference Vegetation Index (NDVI) [-]
   4.3 References

5. SURFACE ALBEDO DATA ............................................................ 5
   5.1 Data generation and processing chain
   5.2 Data files and movies
	5.2.1 MODIS Albedo Climatology [Diffused, Visible (0.3_0.7) and Near-Infrared (0.7_5.0)]
	5.2.2 MODIS Scale Parameters [Diffused, Visible (0.3_0.7) and Near-Infrared (0.7_5.0)] 
   5.3 References

6. CATCHMENT LAND SURFACE MODEL PARAMETES ......................................... 6
   6.1 Data generation and processing chain
   6.2 Data files
	6.2.1 Time scale parameters for moisture transfer between surfexec and rzexec
	6.2.2 Time scale parameters for moisture transfer between rootzone and water table
	6.2.3 Baseflow parameters
	6.2.4 Area fractioning parameters
   6.3 References
`echo "$toc_rout"`

APPENDIX I - mkCatchParam tag, input options, and log ............................ A1

=====================================================================================
====================================== PAGE  1 ======================================
=====================================================================================

1. INTRODUCTION

   This directory contains land boundary data files that are needed by the land models in 
   the GEOS-5 climate modeling system.  The catchment-tiles (computational units at the 
   land surface) have been derived for the GEOS5 ${mygrid} grid 
   `echo "$int_str1"`  The ${mygrid} 
   grid comprises of ${NGLOBAL} number of tiles globally,  out of which ${NTILES} are 
   catchment-tiles.

   An ${NC}x${NR} integer array of tile indices is saved in Fortran binary file  
   "../rst/${gfile}.rst" in little-endian format. 
   In this global raster file, rows are read from the South Pole with the western and 
   sourthern edges of the pixel (1,1) lie along the dateline and the South Pole. Tile index
   numbers 1 - ${NTILES} make up the Earth's land surface. Figure 1 ("plots/US-east.jpg") 
   shows the map of catchment-tiles in the Eastern United States. 

   Soil, vegetation, model specific and other parameters have been computed at each 
   catchment-tile separately and vector-spaced data have been saved as described in below 
   sections. Global maps of important fixed parameters and movies of the daily climatology 
   of seasonal variables are available in "plots" directory.  

=====================================================================================
====================================== PAGE  2 ======================================
=====================================================================================

2. TOPOGRAPHY AND SOIL DATA

   2.1 Data generation and processing chain

         Verdin (2013) produced global 1-arcmin raster arrays of Compound Topographic
	 Index (CTI) statistics: mean, standard deviation and skewness. Verdin (2013) 
	 also delineated the Earth’s land surface into 291,254 hydrologic catchments, 
	 and sub-dividing 30 dateline catchments increased the total number of catchments 
	 to ${NPfafs}. Furthermore, Verdin (2013) produced an associated global 1 arc-minute 
	 raster array of Level 12 Pfafstetter codes to map those hydrologic catchments. 
	 The 1 arc-minute data array was also updated to accommodate the changes stemming
	 from the 30 dateline catchments. 

         The following approach was employed to derive statistics of CTI for each of 
	 the ${NPfafs} catchment using the available CTI statistics at 1 arc-minute 
	 resolution. First, the pixels belonging to a given catchment were identified.  
	 For each of these pixels, a 5000-element array of CTI was constructed using a 
	 3-parameter gamma distribution; the 5000-element array was constructed so as 
	 to have the same mean, standard deviation, and skewness as was identified for 
	 the pixel in the Verdin (2013) dataset.  If the Verdin (2013) data did not 
	 provide a skewness value for the pixel, a Gaussian distribution was used to 
	 construct the 5000-element CTI array.

         Second, the 5000 elements from each of the N pixels making up the given catchment 
	 were combined to form a single, 5000 x N-element array of CTIs.  Finally, the 
	 mean, standard deviation, and skewness of the CTI values in this single larger 
	 array were computed and assigned to the catchment.  "plots/cti.jpg" shows the 
	 resulting global maps of catchment-level CTI statistics ([top] mean, [middle] 
	 standard devisation, and [bottom] skewness).

         The user is referred to De Lannoy et al. (2014) for a complete description of 
	 the procedure used to produce our global datasets of soil type.  The basis of 
	 the datasets is a merging of STATSGO2 (NRCS, 2012) data with the Harmonized 
	 World Soil Data (HWSD, 2009) version 1.21 dataset, with special techniques employed
	 to fill data gaps.  The merging process produced global arrays of sand, clay 
	 and organic matter percentages at 30 arc-second resolution for two soil layers: 
	 0-30cm and 30-100cm. From these data, De Lannoy et al. (2014) defined 253 soil 
	 classes, and from these classes they derived soil hydraulic properties for the 
	 surface and profile layers, separately (Table 1).  The representative (and thus 
	 utilized) soil type for a given hydrological catchment was determined through a 
	 somewhat complex procedure developed by De Lannoy et al. (2014); the upshot of 
	 the procedure is that the chosen type actually appears within the catchment and 
	 has, compared to all of the 30 arc-second pixels in the catchment, roughly median 
	 texture percentages.  The soil hydraulic properties for the 253 classes appear 
	 in Table 1.  See "plots/soil_param.jpg" for corresponding global maps of soil 
	 hydraulic properties (clockwise from upper left: b parameter [-], psis [ m H20],
	 hydraulic conductivity [ms-1], soil depth to bedrock [m], wilting point wetness 
	 [-], and porosity [m3/m3]).

	 The Second Global Soil Wetness Project (Dirmeyer, and Oki, 2002) provided a global,
	 1-degree dataset of soil depth to bedrock. The GSWP-2 soil depth data were spatially
	 interpolated onto 30-arcsec pixels. The interpolated soil depth data were averaged
	 across a given catchment land element to determine that catchment’s effective soil 
	 depth.
         
   2.2 Data files and images

       2.2.1 Tile types, location, area and Pfafstetter catchment mapping
	 file name: ../til/${gfile}.til
	 The 8-line header is followed by ${NGLOBAL} number of rows.
	 do n = 1,${NGLOBAL}
	        read (10,*)type,`echo "${sec2_til}"`
	 end do
	 
	where for each tile:
	 (1)    type      [-]      tile type (100-land; 19-lakes; 20-ice)
	 (2)    longitude [degree] longitude at the centroid of the tile
	 (3)    latitude  [degree] latitude at the centroid of the tile
	 (4)    ig        [-]      i-index of the global grid cell where the tile is located
	 (5)    jg        [-]      j-index of the global grid cell where the tile is located
	 (6)    pfaf_code [-]      ${pfaf_dest} 
	 (7)    pfaf_index[-]      ${pfafin_des} 
	 (8)    cell_frac [-]      fraction of the global grid cell    
	`echo "${sec2_til2}"`
       2.2.2 Western, eastern, southern, northern edges and mean elevation of tiles
	 file name: catchment.def
	 read (10,*) NTILES
	 do n = 1, ${NTILES}
		read (10,'(i8,i8,5(2x,f9.4))') tile_index,pfaf_code,   &
		min_lon,max_lon,min_lat,max_lat, mean_elevation (m) 
         end do
	 
	 where for each tile:
	 (1)	tile_index [-] 	number
	 (2)	pfaf_code [-]	${pfaf_des}
	 (3)    min_lon [degree]   Western edge
	 (4)    max_lon [degree]   Eastern edge
	 (5)    min_lat [degree]   Southern edge
	 (6)    max_lat [degree]   Northern edge
	 (7)    mean_elevation [m] area-averaged elevation
		[Figure 2 : "plots/ELEVATION.jpg"]

       2.2.3 Tile topography - statistics of Compound Topographic Index (CTI) 
	 file name: cti_stats.dat
	 read (10,*) NTILES
	 do n = 1, ${NTILES}
		read (10,'(i8,i8,5(1x,f8.4))') tile_index,pfaf_code,   &
		cti_mean, cti_std, cti_min, cti_max, cti_skew
	 enddo

	 where for each tile:
	 (1)	tile_index [-] 	number
	 (2)	pfaf_code [-]	${pfaf_des}
	 (3)    cti_mean  [log(m)] mean CTI of the underlying hydrologic catchment
		[Figure 3 : "plots/cti.jpg" top panel]
	 (4)    cti_std   standard deviation of CTI of the underlying hydrologic catchment
		[Figure 3 : "plots/cti.jpg" middle panel]
         (5)    cti_min   minimum CTI value in the underlying hydrologic catchment
         (6)    cti_max   maximum CTI value in the underlying hydrologic catchment
         (7)    cti_skew  skewness of CTI of the underlying hydrologic catchment
		[Figure 3 : "plots/cti.jpg" bottom panel]

       2.2.4 Soil Parameters
	 file name: soil_param.dat
	 do n = 1, ${NTILES}
_EOI_
if( $mysoil == HWSD ) then
cat << _EOS1_ > clsm/soil
		read (10,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4)')  &
		tile_index,pfaf_code,soil_class_top,soil_class_com,BEE,              &
		PSIS,POROS,COND, WPWET, soildepth, gravel,OrgCarbon_top,             &
		OrgCarbon_rz,sand_top,clay_top,sand_rz,clay_rz,WPWET_top, POROS_top
	 end do

	 where for each tile:
	 (1)	tile_index [-] 	number
	 (2)	pfaf_code [-]	${pfaf_des}
	 (3)	soil_class_top [-] 	soil class for the surface layer (0-30cm)
	 (4)	soil_class_com  [-]	soil class for the root-zone  (0-100cm)
	 (5)	BEE [-] 	b-parameter of the tension curve
		[Figure 4 : "plots/soil_param.jpg" top-left panel]
	 (6)	PSIS [m H2O] 	air entry pressure for the root-zone (matric potential) 
		[Figure 4 : "plots/soil_param.jpg" top-right panel]
	 (7)	POROS [m3/m3]	soil moisture content at saturation for the root-zone
		[Figure 4 : "plots/soil_param.jpg" middle-left panel]
	 (8)	COND [m/s]	saturated hydraulic conductivity at the surface
		[Figure 4 : "plots/soil_param.jpg" middle-right panel]
	 (9)	WPWET [-]	wilting point/porosity for the root-zone
		[Figure 4 : "plots/soil_param.jpg" bottom-left panel]
	 (10)	soildepth [mm]	depth to bedrock
		[Figure 4 : "plots/soil_param.jpg" bottom-right panel]
	 (11)	gravel [vol% ]	percentage gravel in the surface layer (0-30cm) 
	 (12)	OrgCarbon_top [w%] 	percentage organic carbon in the surface layer (0-30cm)
	 (13)	OrgCarbon_rz [w%] 	percentage organic carbon in the surface layer (0-100cm) 
	 (14)	sand_top [w%]	percentage sand in the surface layer (0-30cm)
	 (15)	clay_top [w%]	percentage clay  in the surface layer (0-30cm)
	 (16)	sand_rz [w%]	percentage sand in the root-zone layer (0-100cm)
	 (17)	clay_rz [w%]	percentage clay in the root-zone layer (0-100cm)
	 (18)	WPWET_top [-]	wilting point/porosity for the surface layer (0-30cm)
	 (19)	POROS_top [m3/m3] soil moisture content at saturation in the surface layer (0-30cm)


	         ========================================================================

                  Sand   Clay   Silt    OC      BEE   PSIS    POROS    WPWET   COND
                  [w%]   [w%]   [w%]   [w%]     [-]  [m H2O] [m3/m3]    [-]    [m/s]

		 ------------------------------------------------------------------------

		 43.333 53.333  3.333  0.258  4.9400 -0.2780  0.4667  0.1301  0.21009E-04
		 36.667 56.667  6.667  0.258  5.4062 -0.3970  0.4969  0.1652  0.14197E-04
		 33.333 53.333 13.333  0.258  6.2117 -0.3664  0.4992  0.1890  0.98520E-05
		 26.667 56.667 16.667  0.258  6.0097 -0.5393  0.5246  0.2050  0.88411E-05
		 23.333 53.333 23.333  0.258  6.3306 -0.4990  0.5214  0.2111  0.73367E-05
 		 16.667 56.667 26.667  0.258  5.9115 -0.7382  0.5463  0.2216  0.68430E-05
 		 13.333 53.333 33.333  0.258  6.0823 -0.6737  0.5406  0.2216  0.58087E-05
  		  6.667 56.667 36.667  0.258  5.5423 -0.9984  0.5655  0.2281  0.54459E-05
 		  3.333 53.333 43.333  0.258  5.6618 -0.8953  0.5581  0.2251  0.45596E-05
 		 53.333 43.333  3.333  0.258  5.0224 -0.1525  0.4237  0.1070  0.20921E-04
 		 46.667 46.667  6.667  0.258  5.8096 -0.1940  0.4497  0.1427  0.14067E-04
 		 43.333 43.333 13.333  0.258  6.4236 -0.1983  0.4535  0.1611  0.95996E-05
 		 36.667 46.667 16.667  0.258  6.5953 -0.2593  0.4750  0.1805  0.86122E-05
 		 33.333 43.333 23.333  0.258  6.6794 -0.2666  0.4730  0.1827  0.70088E-05
		 26.667 46.667 26.667  0.258  6.6552 -0.3492  0.4942  0.1982  0.65632E-05
		 23.333 43.333 33.333  0.258  6.5729 -0.3556  0.4895  0.1946  0.54511E-05
		 16.667 46.667 36.667  0.258  6.4347 -0.4649  0.5111  0.2076  0.51508E-05
		 13.333 43.333 43.333  0.258  6.2961 -0.4671  0.5045  0.2011  0.42112E-05
		  6.667 46.667 46.667  0.258  6.0743 -0.6085  0.5266  0.2120  0.39685E-05
 		  3.333 43.333 53.333  0.258  5.9279 -0.6026  0.5184  0.2037  0.31502E-05
 		 63.333 33.333  3.333  0.258  4.4481 -0.1115  0.3920  0.0773  0.21259E-04
  		 56.667 36.667  6.667  0.258  5.4075 -0.1295  0.4129  0.1116  0.13991E-04
 		 53.333 33.333 13.333  0.258  5.7712 -0.1444  0.4181  0.1250  0.94071E-05
 		 46.667 36.667 16.667  0.258  6.2295 -0.1719  0.4348  0.1462  0.83109E-05
 		 43.333 33.333 23.333  0.258  6.0877 -0.1937  0.4340  0.1451  0.66448E-05
 		 36.667 36.667 26.667  0.258  6.3865 -0.2304  0.4508  0.1630  0.61623E-05
 		 33.333 33.333 33.333  0.258  6.0815 -0.2585  0.4472  0.1565  0.50155E-05
 		 26.667 36.667 36.667  0.258  6.2857 -0.3056  0.4646  0.1728  0.47182E-05
		 23.333 33.333 43.333  0.258  5.9221 -0.3402  0.4589  0.1636  0.37718E-05
		 16.667 36.667 46.667  0.258  6.0565 -0.3991  0.4772  0.1787  0.35560E-05
		 13.333 33.333 53.333  0.258  5.6805 -0.4406  0.4698  0.1677  0.27548E-05
		  6.667 36.667 56.667  0.258  5.7569 -0.5122  0.4889  0.1817  0.25869E-05
		  3.333 33.333 63.333  0.258  5.3927 -0.5605  0.4800  0.1696  0.19300E-05
 		 73.333 23.333  3.333  0.258  3.4528 -0.1022  0.3744  0.0451  0.22840E-04
 		 66.667 26.667  6.667  0.258  4.4073 -0.1107  0.3887  0.0753  0.14465E-04
  		 63.333 23.333 13.333  0.258  4.5541 -0.1324  0.3956  0.0841  0.95901E-05
 		 56.667 26.667 16.667  0.258  5.1544 -0.1471  0.4063  0.1056  0.82079E-05
 		 53.333 23.333 23.333  0.258  4.8747 -0.1788  0.4069  0.1018  0.64538E-05
 		 46.667 26.667 26.667  0.258  5.3580 -0.1983  0.4181  0.1209  0.58360E-05
 		 43.333 23.333 33.333  0.258  4.9350 -0.2414  0.4156  0.1124  0.46596E-05
 		 36.667 26.667 36.667  0.258  5.3435 -0.2656  0.4279  0.1303  0.43009E-05
 		 33.333 23.333 43.333  0.258  4.8657 -0.3231  0.4232  0.1193  0.33652E-05
		 26.667 26.667 46.667  0.258  5.2165 -0.3515  0.4367  0.1363  0.31315E-05
		 23.333 23.333 53.333  0.258  4.7237 -0.4270  0.4301  0.1238  0.23695E-05
		 16.667 26.667 56.667  0.258  5.0266 -0.4581  0.4449  0.1400  0.22088E-05
		 13.333 23.333 63.333  0.258  4.5394 -0.5560  0.4367  0.1267  0.16067E-05
		  6.667 26.667 66.667  0.258  4.8006 -0.5873  0.4527  0.1421  0.14941E-05
 		  3.333 23.333 73.333  0.258  4.3306 -0.7120  0.4432  0.1282  0.10418E-05
 		 83.333 13.333  3.333  0.258  2.3786 -0.1138  0.3741  0.0181  0.26718E-04
  		 76.667 16.667  6.667  0.258  3.1689 -0.1151  0.3803  0.0393  0.16060E-04
 		 73.333 13.333 13.333  0.258  3.1864 -0.1454  0.3890  0.0438  0.10490E-04
 		 66.667 16.667 16.667  0.258  3.7666 -0.1528  0.3922  0.0626  0.85712E-05
 		 63.333 13.333 23.333  0.258  3.4606 -0.1969  0.3944  0.0576  0.66263E-05
 		 56.667 16.667 26.667  0.258  3.9731 -0.2079  0.3986  0.0757  0.57553E-05
 		 53.333 13.333 33.333  0.258  3.5507 -0.2698  0.3975  0.0667  0.45087E-05
 		 46.667 16.667 36.667  0.258  4.0147 -0.2837  0.4033  0.0842  0.40221E-05
		 43.333 13.333 43.333  0.258  3.5437 -0.3706  0.3998  0.0731  0.30816E-05
		 36.667 16.667 46.667  0.258  3.9659 -0.3854  0.4073  0.0901  0.27889E-05
		 33.333 13.333 53.333  0.258  3.4781 -0.5071  0.4018  0.0778  0.20623E-05
		 26.667 16.667 56.667  0.258  3.8626 -0.5190  0.4110  0.0943  0.18813E-05
		 23.333 13.333 63.333  0.258  3.3754 -0.6884  0.4037  0.0814  0.13349E-05
 		 16.667 16.667 66.667  0.258  3.7256 -0.6905  0.4146  0.0973  0.12223E-05
 		 13.333 13.333 73.333  0.258  3.2489 -0.9237  0.4058  0.0842  0.82997E-06
  		  6.667 16.667 76.667  0.258  3.5681 -0.9054  0.4183  0.0993  0.76042E-06
 		  3.333 13.333 83.333  0.258  3.1078 -1.2213  0.4081  0.0862  0.49352E-06
 		 93.333  3.333  3.333  0.258  1.5152 -0.1620  0.3985  0.0043  0.34424E-04
 		 86.667  6.667  6.667  0.258  2.0638 -0.1463  0.3919  0.0135  0.19555E-04
 		 83.333  3.333 13.333  0.258  2.0371 -0.1934  0.4056  0.0153  0.12535E-04
 		 76.667  6.667 16.667  0.258  2.4762 -0.1877  0.3966  0.0265  0.97029E-05
 		 73.333  3.333 23.333  0.258  2.2247 -0.2528  0.4036  0.0227  0.73579E-05
		 66.667  6.667 26.667  0.258  2.6386 -0.2530  0.3963  0.0350  0.60754E-05
		 63.333  3.333 33.333  0.258  2.2989 -0.3442  0.3997  0.0282  0.46641E-05
		 56.667  6.667 36.667  0.258  2.6944 -0.3491  0.3946  0.0413  0.39726E-05
		 53.333  3.333 43.333  0.258  2.3134 -0.4809  0.3953  0.0327  0.29788E-05
		 46.667  6.667 46.667  0.258  2.6895 -0.4881  0.3926  0.0463  0.25869E-05
 		 43.333  3.333 53.333  0.258  2.2907 -0.6827  0.3911  0.0368  0.18693E-05
 		 36.667  6.667 56.667  0.258  2.6458 -0.6860  0.3906  0.0506  0.16455E-05
  		 33.333  3.333 63.333  0.258  2.2430 -0.9768  0.3872  0.0407  0.11391E-05
 		 26.667  6.667 66.667  0.258  2.5756 -0.9632  0.3890  0.0544  0.10125E-05
 		 23.333  3.333 73.333  0.258  2.1774 -1.3983  0.3839  0.0444  0.66957E-06
 		 16.667  6.667 76.667  0.258  2.4873 -1.3439  0.3877  0.0578  0.59923E-06
 		 13.333  3.333 83.333  0.258  2.0991 -1.9899  0.3812  0.0482  0.37816E-06
 		  6.667  6.667 86.667  0.258  2.3867 -1.8539  0.3870  0.0609  0.33987E-06
 		  3.333  3.333 93.333  0.258  2.0118 -2.7985  0.3791  0.0519  0.20465E-06
		 43.333 53.333  3.333  0.463  5.1607 -0.2689  0.4659  0.1363  0.21175E-04
		 36.667 56.667  6.667  0.463  5.6703 -0.3785  0.4953  0.1719  0.14263E-04
		 33.333 53.333 13.333  0.463  6.4853 -0.3529  0.4990  0.1956  0.10015E-04
		 26.667 56.667 16.667  0.463  6.2983 -0.5122  0.5237  0.2119  0.89597E-05
		 23.333 53.333 23.333  0.463  6.6034 -0.4786  0.5218  0.2179  0.75191E-05
 		 16.667 56.667 26.667  0.463  6.1879 -0.6985  0.5459  0.2285  0.69925E-05
 		 13.333 53.333 33.333  0.463  6.3353 -0.6437  0.5415  0.2283  0.59999E-05
  		  6.667 56.667 36.667  0.463  5.7908 -0.9414  0.5658  0.2349  0.56092E-05
 		  3.333 53.333 43.333  0.463  5.8848 -0.8520  0.5597  0.2317  0.47449E-05
 		 53.333 43.333  3.333  0.463  5.1941 -0.1513  0.4253  0.1122  0.21476E-04
 		 46.667 46.667  6.667  0.463  6.0259 -0.1910  0.4508  0.1486  0.14393E-04
 		 43.333 43.333 13.333  0.463  6.6429 -0.1958  0.4557  0.1672  0.99339E-05
 		 36.667 46.667 16.667  0.463  6.8389 -0.2540  0.4767  0.1870  0.88832E-05
 		 33.333 43.333 23.333  0.463  6.9061 -0.2618  0.4758  0.1892  0.73089E-05
		 26.667 46.667 26.667  0.463  6.8971 -0.3404  0.4965  0.2048  0.68220E-05
		 23.333 43.333 33.333  0.463  6.7927 -0.3473  0.4930  0.2012  0.57264E-05
		 16.667 46.667 36.667  0.463  6.6620 -0.4509  0.5139  0.2143  0.53932E-05
		 13.333 43.333 43.333  0.463  6.5004 -0.4536  0.5086  0.2077  0.44549E-05
		  6.667 46.667 46.667  0.463  6.2787 -0.5873  0.5300  0.2185  0.41843E-05
 		  3.333 43.333 53.333  0.463  6.1105 -0.5818  0.5231  0.2102  0.33546E-05
 		 63.333 33.333  3.333  0.463  4.5680 -0.1119  0.3957  0.0814  0.22188E-04
  		 56.667 36.667  6.667  0.463  5.5645 -0.1294  0.4161  0.1167  0.14566E-04
 		 53.333 33.333 13.333  0.463  5.9283 -0.1439  0.4224  0.1304  0.98974E-05
 		 46.667 36.667 16.667  0.463  6.4119 -0.1707  0.4388  0.1520  0.87205E-05
 		 43.333 33.333 23.333  0.463  6.2554 -0.1918  0.4391  0.1509  0.70451E-05
 		 36.667 36.667 26.667  0.463  6.5748 -0.2272  0.4554  0.1691  0.65150E-05
 		 33.333 33.333 33.333  0.463  6.2509 -0.2539  0.4529  0.1626  0.53568E-05
 		 26.667 36.667 36.667  0.463  6.4713 -0.2993  0.4698  0.1792  0.50240E-05
		 23.333 33.333 43.333  0.463  6.0880 -0.3314  0.4653  0.1699  0.40567E-05
		 16.667 36.667 46.667  0.463  6.2336 -0.3881  0.4830  0.1852  0.38123E-05
		 13.333 33.333 53.333  0.463  5.8390 -0.4255  0.4768  0.1740  0.29824E-05
		  6.667 36.667 56.667  0.463  5.9210 -0.4942  0.4954  0.1880  0.27912E-05
		  3.333 33.333 63.333  0.463  5.5403 -0.5364  0.4877  0.1758  0.21026E-05
 		 73.333 23.333  3.333  0.463  3.5376 -0.1029  0.3797  0.0482  0.24151E-04
 		 66.667 26.667  6.667  0.463  4.5158 -0.1110  0.3938  0.0795  0.15279E-04
  		 63.333 23.333 13.333  0.463  4.6643 -0.1321  0.4017  0.0885  0.10227E-04
 		 56.667 26.667 16.667  0.463  5.2815 -0.1465  0.4121  0.1105  0.87411E-05
 		 53.333 23.333 23.333  0.463  4.9923 -0.1768  0.4138  0.1068  0.69388E-05
 		 46.667 26.667 26.667  0.463  5.4915 -0.1957  0.4247  0.1262  0.62641E-05
 		 43.333 23.333 33.333  0.463  5.0549 -0.2361  0.4232  0.1176  0.50491E-05
 		 36.667 26.667 36.667  0.463  5.4787 -0.2595  0.4352  0.1358  0.46511E-05
 		 33.333 23.333 43.333  0.463  4.9858 -0.3123  0.4315  0.1246  0.36739E-05
		 26.667 26.667 46.667  0.463  5.3509 -0.3396  0.4446  0.1419  0.34107E-05
		 23.333 23.333 53.333  0.463  4.8427 -0.4075  0.4392  0.1291  0.26054E-05
		 16.667 26.667 56.667  0.463  5.1581 -0.4377  0.4535  0.1457  0.24220E-05
		 13.333 23.333 63.333  0.463  4.6561 -0.5235  0.4465  0.1319  0.17785E-05
		  6.667 26.667 66.667  0.463  4.9275 -0.5543  0.4620  0.1477  0.16488E-05
 		  3.333 23.333 73.333  0.463  4.4440 -0.6610  0.4535  0.1332  0.11606E-05
 		 83.333 13.333  3.333  0.463  2.4488 -0.1150  0.3806  0.0202  0.28480E-04
  		 76.667 16.667  6.667  0.463  3.2528 -0.1157  0.3867  0.0424  0.17136E-04
 		 73.333 13.333 13.333  0.463  3.2749 -0.1454  0.3964  0.0473  0.11284E-04
 		 66.667 16.667 16.667  0.463  3.8621 -0.1520  0.3994  0.0667  0.92264E-05
 		 63.333 13.333 23.333  0.463  3.5518 -0.1943  0.4026  0.0616  0.71924E-05
 		 56.667 16.667 26.667  0.463  4.0706 -0.2043  0.4067  0.0800  0.62486E-05
 		 53.333 13.333 33.333  0.463  3.6401 -0.2624  0.4066  0.0707  0.49370E-05
 		 46.667 16.667 36.667  0.463  4.1114 -0.2750  0.4122  0.0886  0.44033E-05
		 43.333 13.333 43.333  0.463  3.6302 -0.3545  0.4097  0.0770  0.34030E-05
		 36.667 16.667 46.667  0.463  4.0609 -0.3679  0.4170  0.0944  0.30776E-05
		 33.333 13.333 53.333  0.463  3.5618 -0.4762  0.4124  0.0816  0.22961E-05
		 26.667 16.667 56.667  0.463  3.9558 -0.4872  0.4214  0.0985  0.20921E-05
		 23.333 13.333 63.333  0.463  3.4565 -0.6339  0.4152  0.0849  0.14979E-05
 		 16.667 16.667 66.667  0.463  3.8172 -0.6369  0.4258  0.1013  0.13692E-05
 		 13.333 13.333 73.333  0.463  3.3282 -0.8331  0.4180  0.0873  0.93830E-06
  		  6.667 16.667 76.667  0.463  3.6583 -0.8197  0.4302  0.1030  0.85767E-06
 		  3.333 13.333 83.333  0.463  3.1858 -1.0782  0.4210  0.0889  0.56188E-06
 		 93.333  3.333  3.333  0.463  1.5786 -0.1649  0.4056  0.0053  0.36768E-04
 		 86.667  6.667  6.667  0.463  2.1387 -0.1480  0.3992  0.0155  0.20953E-04
 		 83.333  3.333 13.333  0.463  2.1170 -0.1947  0.4137  0.0177  0.13520E-04
 		 76.667  6.667 16.667  0.463  2.5599 -0.1875  0.4049  0.0295  0.10497E-04
 		 73.333  3.333 23.333  0.463  2.3058 -0.2507  0.4127  0.0256  0.80154E-05
		 66.667  6.667 26.667  0.463  2.7213 -0.2488  0.4055  0.0383  0.66362E-05
		 63.333  3.333 33.333  0.463  2.3760 -0.3349  0.4097  0.0311  0.51315E-05
		 56.667  6.667 36.667  0.463  2.7728 -0.3371  0.4047  0.0446  0.43806E-05
		 53.333  3.333 43.333  0.463  2.3845 -0.4577  0.4063  0.0355  0.33096E-05
		 46.667  6.667 46.667  0.463  2.7624 -0.4613  0.4035  0.0494  0.28793E-05
 		 43.333  3.333 53.333  0.463  2.3553 -0.6337  0.4030  0.0392  0.20972E-05
 		 36.667  6.667 56.667  0.463  2.7131 -0.6334  0.4024  0.0533  0.18481E-05
  		 33.333  3.333 63.333  0.463  2.3012 -0.8822  0.4000  0.0426  0.12900E-05
 		 26.667  6.667 66.667  0.463  2.6382 -0.8673  0.4016  0.0565  0.11472E-05
 		 23.333  3.333 73.333  0.463  2.2302 -1.2266  0.3976  0.0457  0.76529E-06
 		 16.667  6.667 76.667  0.463  2.5462 -1.1784  0.4012  0.0593  0.68466E-06
 		 13.333  3.333 83.333  0.463  2.1476 -1.6929  0.3957  0.0486  0.43603E-06
 		  6.667  6.667 86.667  0.463  2.4431 -1.5812  0.4012  0.0617  0.39144E-06
 		  3.333  3.333 93.333  0.463  2.0575 -2.3066  0.3944  0.0514  0.23795E-06
		 43.333 53.333  3.333  1.123  5.7903 -0.2516  0.4632  0.1531  0.18072E-04
		 36.667 56.667  6.667  1.123  6.4387 -0.3397  0.4905  0.1899  0.12036E-04
		 33.333 53.333 13.333  1.123  7.2677 -0.3247  0.4983  0.2136  0.87933E-05
		 26.667 56.667 16.667  1.123  7.1389 -0.4524  0.5208  0.2303  0.77820E-05
		 23.333 53.333 23.333  1.123  7.3822 -0.4334  0.5230  0.2363  0.67840E-05
 		 16.667 56.667 26.667  1.123  6.9906 -0.6076  0.5449  0.2471  0.62441E-05
 		 13.333 53.333 33.333  1.123  7.0517 -0.5739  0.5447  0.2467  0.55559E-05
  		  6.667 56.667 36.667  1.123  6.5053 -0.8072  0.5667  0.2531  0.51436E-05
 		  3.333 53.333 43.333  1.123  6.5043 -0.7485  0.5647  0.2492  0.45042E-05
 		 53.333 43.333  3.333  1.123  5.6646 -0.1522  0.4305  0.1270  0.19503E-04
 		 46.667 46.667  6.667  1.123  6.6262 -0.1874  0.4542  0.1651  0.12927E-04
 		 43.333 43.333 13.333  1.123  7.2483 -0.1933  0.4629  0.1843  0.92645E-05
 		 36.667 46.667 16.667  1.123  7.5172 -0.2449  0.4820  0.2047  0.81934E-05
 		 33.333 43.333 23.333  1.123  7.5330 -0.2539  0.4851  0.2074  0.69916E-05
		 26.667 46.667 26.667  1.123  7.5692 -0.3226  0.5038  0.2232  0.64543E-05
		 23.333 43.333 33.333  1.123  7.3964 -0.3305  0.5042  0.2199  0.56116E-05
		 16.667 46.667 36.667  1.123  7.2859 -0.4200  0.5232  0.2329  0.52275E-05
		 13.333 43.333 43.333  1.123  7.0512 -0.4237  0.5217  0.2263  0.44667E-05
		  6.667 46.667 46.667  1.123  6.8258 -0.5376  0.5412  0.2365  0.41501E-05
 		  3.333 43.333 53.333  1.123  6.5854 -0.5332  0.5382  0.2279  0.34373E-05
 		 63.333 33.333  3.333  1.123  4.8965 -0.1161  0.4076  0.0940  0.21304E-04
  		 56.667 36.667  6.667  1.123  5.9940 -0.1324  0.4267  0.1316  0.13865E-04
 		 53.333 33.333 13.333  1.123  6.3656 -0.1462  0.4365  0.1464  0.97555E-05
 		 46.667 36.667 16.667  1.123  6.9168 -0.1712  0.4514  0.1690  0.85177E-05
 		 43.333 33.333 23.333  1.123  6.7274 -0.1905  0.4553  0.1685  0.71198E-05
 		 36.667 36.667 26.667  1.123  7.0986 -0.2231  0.4701  0.1874  0.65214E-05
 		 33.333 33.333 33.333  1.123  6.7291 -0.2462  0.4712  0.1812  0.55436E-05
 		 26.667 36.667 36.667  1.123  6.9853 -0.2875  0.4866  0.1981  0.51475E-05
		 23.333 33.333 43.333  1.123  6.5529 -0.3135  0.4857  0.1888  0.42936E-05
		 16.667 36.667 46.667  1.123  6.7163 -0.3643  0.5018  0.2041  0.39930E-05
		 13.333 33.333 53.333  1.123  6.2739 -0.3922  0.4992  0.1929  0.32244E-05
		  6.667 36.667 56.667  1.123  6.3538 -0.4531  0.5161  0.2064  0.29850E-05
		  3.333 33.333 63.333  1.123  5.9299 -0.4812  0.5121  0.1938  0.23191E-05
 		 73.333 23.333  3.333  1.123  3.7828 -0.1080  0.3970  0.0583  0.24231E-04
 		 66.667 26.667  6.667  1.123  4.8248 -0.1153  0.4102  0.0924  0.15268E-04
  		 63.333 23.333 13.333  1.123  4.9912 -0.1352  0.4214  0.1030  0.10544E-04
 		 56.667 26.667 16.667  1.123  5.6513 -0.1484  0.4309  0.1262  0.89683E-05
 		 53.333 23.333 23.333  1.123  5.3491 -0.1757  0.4358  0.1229  0.73425E-05
 		 46.667 26.667 26.667  1.123  5.8861 -0.1928  0.4457  0.1434  0.65905E-05
 		 43.333 23.333 33.333  1.123  5.4247 -0.2272  0.4476  0.1348  0.54776E-05
 		 36.667 26.667 36.667  1.123  5.8819 -0.2482  0.4584  0.1538  0.50120E-05
 		 33.333 23.333 43.333  1.123  5.3589 -0.2903  0.4582  0.1423  0.40811E-05
		 26.667 26.667 46.667  1.123  5.7512 -0.3148  0.4701  0.1603  0.37596E-05
		 23.333 23.333 53.333  1.123  5.2115 -0.3651  0.4681  0.1469  0.29598E-05
		 16.667 26.667 56.667  1.123  5.5453 -0.3923  0.4811  0.1641  0.27276E-05
		 13.333 23.333 63.333  1.123  5.0134 -0.4511  0.4776  0.1494  0.20637E-05
		  6.667 26.667 66.667  1.123  5.2914 -0.4798  0.4917  0.1654  0.18947E-05
 		  3.333 23.333 73.333  1.123  4.7819 -0.5469  0.4868  0.1499  0.13737E-05
 		 83.333 13.333  3.333  1.123  2.6654 -0.1221  0.4016  0.0276  0.29368E-04
  		 76.667 16.667  6.667  1.123  3.5097 -0.1212  0.4076  0.0533  0.17714E-04
 		 73.333 13.333 13.333  1.123  3.5567 -0.1498  0.4202  0.0599  0.11980E-04
 		 66.667 16.667 16.667  1.123  4.1634 -0.1544  0.4229  0.0807  0.98100E-05
 		 63.333 13.333 23.333  1.123  3.8516 -0.1930  0.4292  0.0758  0.78556E-05
 		 56.667 16.667 26.667  1.123  4.3876 -0.2002  0.4327  0.0953  0.68268E-05
 		 53.333 13.333 33.333  1.123  3.9448 -0.2497  0.4357  0.0856  0.55423E-05
 		 46.667 16.667 36.667  1.123  4.4339 -0.2587  0.4408  0.1045  0.49380E-05
		 43.333 13.333 43.333  1.123  3.9344 -0.3214  0.4414  0.0921  0.39226E-05
		 36.667 16.667 46.667  1.123  4.3842 -0.3310  0.4480  0.1105  0.35386E-05
		 33.333 13.333 53.333  1.123  3.8632 -0.4098  0.4467  0.0964  0.27145E-05
		 26.667 16.667 56.667  1.123  4.2769 -0.4180  0.4548  0.1144  0.24632E-05
		 23.333 13.333 63.333  1.123  3.7541 -0.5159  0.4519  0.0992  0.18141E-05
 		 16.667 16.667 66.667  1.123  4.1335 -0.5195  0.4615  0.1166  0.16488E-05
 		 13.333 13.333 73.333  1.123  3.6209 -0.6394  0.4571  0.1007  0.11625E-05
  		  6.667 16.667 76.667  1.123  3.9666 -0.6344  0.4682  0.1174  0.10549E-05
 		  3.333 13.333 83.333  1.123  3.4721 -0.7787  0.4625  0.1011  0.71122E-06
 		 93.333  3.333  3.333  1.123  1.7794 -0.1797  0.4287  0.0097  0.38230E-04
 		 86.667  6.667  6.667  1.123  2.3774 -0.1581  0.4228  0.0235  0.21996E-04
 		 83.333  3.333 13.333  1.123  2.3773 -0.2058  0.4400  0.0273  0.14501E-04
 		 76.667  6.667 16.667  1.123  2.8350 -0.1935  0.4315  0.0410  0.11360E-04
 		 73.333  3.333 23.333  1.123  2.5792 -0.2538  0.4421  0.0369  0.88676E-05
		 66.667  6.667 26.667  1.123  3.0039 -0.2455  0.4351  0.0511  0.74004E-05
		 63.333  3.333 33.333  1.123  2.6478 -0.3212  0.4422  0.0431  0.58547E-05
		 56.667  6.667 36.667  1.123  3.0524 -0.3154  0.4372  0.0576  0.50311E-05
		 53.333  3.333 43.333  1.123  2.6484 -0.4121  0.4418  0.0473  0.38924E-05
		 46.667  6.667 46.667  1.123  3.0351 -0.4063  0.4388  0.0622  0.34030E-05
 		 43.333  3.333 53.333  1.123  2.6090 -0.5313  0.4413  0.0504  0.25406E-05
 		 36.667  6.667 56.667  1.123  2.9779 -0.5220  0.4405  0.0654  0.22457E-05
  		 33.333  3.333 63.333  1.123  2.5448 -0.6844  0.4412  0.0526  0.16083E-05
 		 26.667  6.667 66.667  1.123  2.8957 -0.6654  0.4423  0.0676  0.14316E-05
 		 23.333  3.333 73.333  1.123  2.4648 -0.8759  0.4414  0.0543  0.98080E-06
 		 16.667  6.667 76.667  1.123  2.7974 -0.8385  0.4445  0.0691  0.87634E-06
 		 13.333  3.333 83.333  1.123  2.3751 -1.1086  0.4422  0.0555  0.57369E-06
 		  6.667  6.667 86.667  1.123  2.6893 -1.0407  0.4470  0.0699  0.51317E-06
 		  3.333  3.333 93.333  1.123  2.2797 -1.3816  0.4436  0.0563  0.32091E-06
                 PEAT   (N/A)         >=8.72  3.4130 -1.7600  0.8000  0.2162  0.78600E-06
		 ------------------------------------------------------------------------

		 Table 1: Soil Hydraulic Properties for 253 soil classes (adapted from De 
		          Lannoy et al., 2014).
		
   2.3 References
         NRCS Soil Survey Staff, USDA (2012): General Soil Map (STATSGO2) [United States]. 
	    Available online at http://websoilsurvey.nrcs.usda.gov/. Last accessed July 2015.
	 HWSD (2009): Harmonized World Soil Database Version 1.2 Technical Note, Food and 
            Agric. Organ., Rome  (Available at: http://webarchive.iiasa.ac.at/Research/LUC/
            External-World-soil-database/HTML/).
	 De Lannoy, G.J.M., Koster, R.D., Reichle, R.H., Mahanama, S.P., Liu, Q. (2014). 
	    An Updated Treatment of Soil Texture and Associated Hydraulic Properties in a Global 
	    Land Modeling System. Journal of Advances in Modeling Earth Systems, 06, 
	    doi:10.1002/2014MS000330. 
_EOS1_
else
cat << _EOS2_ > clsm/soil	 
 		read (10,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4)')  &
		tile_index,pfaf_code,soil_class_top,              &
		soil_class_com,BEE, PSIS,POROS,COND,              &
		WPWET,soildepth	
	 end do

	 where for each tile:
	 (1)	tile_index [-] 	number
	 (2)	pfaf_code [-]	${pfaf_des}
	 (3)	soil_class_top [-] 	soil class for the surface layer (0-30cm)
	 (4)	soil_class_com  [-]	soil class for the root-zone  (0-100cm)
	 (5)	BEE [-] 	b-parameter of the tension curve
		[Figure 4 : "plots/soil_param.jpg" top-left panel]
	 (6)	PSIS [m H2O] 	air entry pressure for the root-zone (matric potential) 
		[Figure 4 : "plots/soil_param.jpg" top-right panel]
	 (7)	POROS [m3/m3]	soil moisture content at saturation for the root-zone
		[Figure 4 : "plots/soil_param.jpg" middle-left panel]
	 (8)	COND [m/s]	saturated hydraulic conductivity at the surface
		[Figure 4 : "plots/soil_param.jpg" middle-right panel]
	 (9)	WPWET [-]	wilting point/porosity for the root-zone
		[Figure 4 : "plots/soil_param.jpg" bottom-left panel]
	 (10)	soildepth [mm]	depth to bedrock 
		[Figure 4 : "plots/soil_param.jpg" bottom-right panel]

   2.3 References
         Reynolds, C. A., T. J. Jackson, and W. J. Rawls (2000): Estimating soil 
	    water-holding capacities by linking the Food and Agriculture Organization 
            Soil map of the world with global pedon databases and continuous pedotransfer 
            functions, Water Resour. Res., 36(12), 36533662, doi:10.1029/2000WR900130.
_EOS2_
endif
cat << _EOV1_ > clsm/veg1
	 Dirmeyer, P. and Oki, T. (2002): The Second Global Soil Wetness project (GSWP-2) 
	    Science 2 and Implementation Plan. IGPO Publication Series No. 37, 64p.	
         Verdin, K (2013): Final Report - High Resolution Topographic Analysis for GMAO's 
	    Catchment LSM, pp 21, available from Global Modeling and Assimilation Office, 
	    610.1 NASA/GSFC.

=====================================================================================
====================================== PAGE  3 ======================================
=====================================================================================

3. VEGETATION CLASSIFICATION DATA

   3.1 Data generation and processing chain
  
	 The Mosaic model (Koster and Suarez, 1996) utilizes 8 vegetation types. The current 
	 version of CLSM uses 6 types, with the Mosaic bare soil and desert soil types grouped 
	 with Shrubland. The reduction in type number is explained by the fact that a 
	 shrubland type can be made to act like a bare soil or desert type simply by assigning 
	 an appropriately low value of leaf area index, or LAI, to the surface element.  CLSM 
	 uses mean seasonal cycles of LAI from global data products (Section 4.1), thereby 
	 ensuring that bare soil and desert soil behavior will appear in the correct places, 
	 and in any case albedos are forced to match those for diffuse radiation in MODIS 
	 observations (Section 5.1).  The absence of the explicit definition of the two bare 
	 soil types has an insignificant effect on model behavior. 
            
         `echo "${sec3_veg_des}"`

	 Global 30-arcsec canopy height data, used (in some implementations) in the calculation 
	 of surface roughness, were obtained from NASA’s Jet Propulsion Laboratory (Simrad et al.,
	 2011). These heights were spatially aggregated to catchment surface elements 
	 (plots/Canopy_Height_onTiles.jpg).

	 Global 6km aeolian aerodynamic roughness length data (Prigent et al., 2012) were 
	 pre-prcessed at GMAO to the 0.1x0.1-degree grid, and used to compute ASCAT surface roughness 
       	 length (ASCATZ0) on surface catchment elements.

   3.2 Data files and images
       3.2.1 Mosaic vegetation types and fractions
	 file name: mosaic_veg_typs_fracs
	 do n = 1, ${NTILES}
		read (10,*)tile_index,pfaf_code,   &
		primary_veg_type,secondary_veg_type, primary_veg_frac,         &
		secondary_veg_frac, canopy_height, ASCATZ0
	 end do

	 where for each tile:
	 (1)	tile_index [-] 	number
	 (2)	pfaf_code [-]	${pfaf_des}
	 (3)	primary_veg_type [-]	primary vegetation type    
		[Figure 5 : "plots/mosaic_prim.jpg"]
	 (4)	secondary_veg_type [-]	secondary vegetation type
	 (5)	primary_veg_frac [-] 	primary vegetation fraction
	 (6)	secondary_veg_frac [-] 	secondary vegetation fraction
	 (7)    canopy_height [m]
                [Figure 6 : "plots/Canopy_Height_onTiles.jpg"]
	 (8)    ASCATz0 [m] "plots/ascat_Z0.jpg; plots/icarus_Z0.jpg; and plots/merged_Z0.jpg"

	 with Mosaic types: 1 - Broadleaf Evergreen; 2 - Broadleaf Deciduous; 3 - Needleleaf; 
		4 - Grassland; 5 - Broadleaf Shrubs; 6 - Dwarf Trees 

       3.2.2 vegdyn input data (mosaic primary type, canopy height, and roughness) for GEOS5
       . file name: vegdyn.data or ../vegdyn_*.dat
         file format: fortran binaries, little_endian
             read(10) (primary_veg_type(n),n=1,${NTILES})
             read(10) (canopy_height   (n),n=1,${NTILES})
	     read(10) (ASCATz0         (n),n=1,${NTILES})
	     
_EOV1_
if( $MYMASK == GEOS5_10arcsec_mask | $MYMASK == GEOS5_10arcsec_mask.nc | $MYMASK == GEOS5_10arcsec_mask_freshwater-lakes.nc ) then
cat << _EOV2_ > clsm/veg2

       3.2.3 CLM/CLM4.5, CLM/CLM4.5-carbon, CLM4.5 and CLM4.5-carbon vegetation types and fractions 
	 file names: CLM_veg_typs_fracs and  CLM4.5_veg_typs_fracs
	 do n = 1, ${NTILES} 
	 read (10,'(2I8,4I3,4f7.2,2I3,2f7.2)')       &
            tile_index,pfaf_code,                    &
	    CLM-C_pt1,CLM-C_pt2,CLM-C_st1,CLM-C_st2, &
	    CLM-C_pf1,CLM-C_pf2,CLM-C_sf1,CLM-C_sf2, &
	    CLM_pt, CLM_st, CLM_pf, CLM_sf
	 enddo

	 where for each tile:
	 (1) tile_index [-]	number
	 (2) pfaf_code [-]    ${pfaf_des}
	 (3) CLM-C_pt1 [-]      CLM-Carbon primary type 1  
		[Figure 7a : top panel of "plots/CLM-Carbon_PRIM_veg_typs.jpg and plots/CLM4.5-Carbon_PRIM_veg_typs.jpg"]			    
	 (4) CLM-C_pt2 [-]	CLM-Carbon primary type 2 (moisture stressed only) 
		[Figure 7b : bottom panel of "plots/CLM-Carbon_PRIM_veg_typs.jpg and plots/CLM4.5-Carbon_PRIM_veg_typs.jpg"] 
	 (5) CLM-C_st1 [-]	CLM-Carbon secondary type 1 
		[Figure 8a : top panel of "plots/CLM-Carbon_SEC_veg_typs.jpg and plots/CLM4.5-Carbon_SEC_veg_typs.jpg"]			    
	 (6) CLM-C_st2 [-]	CLM-Carbon secondary type 2 (moisture stressed only) 
		[Figure 8b :  bottom panel of "plots/CLM-Carbon_SEC_veg_typs.jpg and plots/CLM4.5-Carbon_SEC_veg_typs.jpg"]
	 (7) CLM-C_pf1 [-]	CLM-Carbon fraction of 1st primary type
	 (8) CLM-C_pf2 [-]	CLM-Carbon fraction of 2nd primary type (moisture stressed only)
	 (9) CLM-C_sf1 [-]	CLM-Carbon fraction of 1st secondary type
	 (10)CLM-C_sf2 [-]	CLM-Carbon fraction of 2nd secondary type (moisture stressed only)
	 (11)CLM_pt [-]         CLM primary type 
		[Figure 9 : "plots/CLM_PRIM_veg_typs.jpg and plots/CLM4.5_PRIM_veg_typs.jpg"]
	 (12)CLM_st [-]         CLM secondary type 
		[Figure 10: "plots/CLM_SEC_veg_typs.jpg and plots/CLM4.5_SEC_veg_typs.jpg"]
	 (13)CLM_pf [-]         CLM fraction of primary type
	 (14)CLM_sf [-]         CLM fraction of secondary type

	Please see below Table 2 for CLM (CLM4.5)  and CLM-Carbon (CLM4.5-Carbon) land cover classification
		
	===================================================================================       
	Land Cover                                         CLM (CLM4.5)	CLM-Carbon  Map
								     (CLM4.5-Carbon)    
								Class	Class       Legend

	-----------------------------------------------------------------------------------
	Bare                                 			  1  	     -	    BARE 
	Needleleaf evergreen temperate tree  			  2	     1      NLEt 
	Needleleaf evergreen boreal tree     			  3 	     2      NLEB 
	Needleleaf deciduous boreal tree     			  4  	     3      NLDB 
	Broadleaf evergreen tropical tree    			  5 	     4      BLET 
	Broadleaf evergreen temperate tree   			  6 	     5      BLEt 
	Broadleaf deciduous tropical tree    			  7 	     6      BLDT 
	Broadleaf deciduous temperate tree   			  8 	     7      BLDt 
	Broadleaf deciduous boreal tree      			  9 	     8      BLDB 
	Broadleaf evergreen temperate shrub  			 10 	     9      BLEtS
	Broadleaf deciduous temperate shrub  			 11 	    10      BLDtS
	Broadleaf deciduous temperate shrub[moisture stress only] -	    11      BLDtSm
	Broadleaf deciduous boreal shrub     			 12 	    12      BLDBS
	Arctic c3 grass                      			 13 	    13      AC3G 
	Cool c3 grass                        			 14 	    14      CC3G 
	Cool c3 grass [moisture stress only] 			  -         15      CC3Gm
	Warm c4 grass                        			 15 	    16      WC4G 
	Warm c4 grass [moisture stress only]   			  -	    17      WC4Gm
	Crop                                 			 16 	    18      CROP (-) 
	Crop [moisture stress only]          			  -	    19      CROPm(-)
        (C3_crop)                                               (16)       (18)     C3CROP                                 
	(C3_irrigated)                                          (17)       (19)     C3IRR
	(Corn)                                                  (18)       (20)     CORN
	(Irrigated corn)                                        (19)       (21)     ICORN     
	(Spring temperate cereal)                               (20)       (22)     STCER
	(Irrigated spring temperate cereal)                     (21)       (23)     ISTCER
	(winter temperate cereal)                               (22)       (24)     WTCER
        (Irrigated winter temperate cereal)                     (23)       (25)     IWTCER
	(Soybean)                                               (24)       (26)     SOYB
	(Irrigated Soybean)                                     (25)       (27)     ISOYB
	Water                                                    17          -
        -----------------------------------------------------------------------------------
	         Table 2: CLM and CLM-Carbon land cover classification description.  
                          CLM-4.5 and CLM-4.5-Carbon types are in brackets.

       3.2.4 Nitrogen Deposition, annual mean 2m Tair, soil back gorund albedo
	 file name: CLM_Ndep_SoilAlb
	 do n = 1, ${NTILES}
               read (10, '(f10.4,4f7.4,2f8.3)')    &
               NDEP,VISDR, VISDF, NIRDR, NIRDF, T2_M, T2_S 
	enddo

	Where for each tile:

        (1) NDEP [ng m-2 s-1]	Nitrogen deposition 
            [Figure 11a: "plots/CLM_Ndep_T2m.jpg" top panel]
        (2) VISDR [-]	Direct visible soil background albedo
            [Figure 12a: "plots/SoilAlb.jpg" top-left panel]
        (3) VISDF [-]	Diffuse visible soil background albedo
            [Figure 12b: "plots/SoilAlb.jpg" top-right panel]
        (4) NIRDR [-]	Direct near-infrared soil background albedo 
            [Figure 12c: "plots/SoilAlb.jpg" bottom-left panel]
        (5) NIRDF [-]	Diffuse near-infrared soil background albedo
            [Figure 12d: "plots/SoilAlb.jpg" bottom-right panel]
        (6) T2_M [K]	Mean annual 2m air temperature from MERRA-2 (averaged over 1980-2014)
            [Figure 11b: "plots/CLM_Ndep_T2m.jpg" middle panel]
        (7) T2_S [K]	Mean annual 2m air temperature from Sheffield et al. (2006) (averaged over 1980-2014)
            [Figure 11c: "plots/CLM_Ndep_T2m.jpg" bottom panel]

       3.2.5 CLM4.5 Peak month for agri fire (ABM), peat fraction (PEATF), gross domestic product (GDP),
             population density (HDM)
         file name: CLM4.5_abm_peatf_gdp_hdm_fc
         do n = 1, ${NTILES}
               read (10,'(2I8, i3, f8.4, f8.2, f10.2, f8.4)') &
               TID, CID, ABM, PEATF, GDP, HDM, FC
         end do

         where for each tile:

         (1) TID tile_index [-] number
         (2) CID pfaf_index [-]
         (3) ABM [-] peak month for agri fire
         (4) PEATF [-] peat fraction
         (5) GDP [1000 1995 USD / capita] gross domestic product
         (6) HDM [individuals per km-2 (based on 2010 census)] population density
         (7) FC [m3/m3] soil field capacity

       3.2.6 CLM4.5 Lightening Frequency Climatology
         file name: lnfm.dat (or ../lnfm_clim*.data)
         file format: fortran binaries, little_endian
              The read statement should be MAPL_ReadForcing compatible. Each data record
              is preceded by a header containing start and end times of the period that
              the data have been averaged for. Corresponding start dates for climatological
              data records are given below:
              `echo "${GSWP2_DATES}"`

              Loop below through until the last data record:
              read (10) Year_Begin,Month_Begin,Day_Begin,Hour_Begin,Minute_Begin,Secs_Begin,
              Year_End,Month_End,Day_End,Hour_End,Minute_End,Secs_End (Float Numbers)
              read(10) (data(n),n=1,${NTILES})

_EOV2_
endif
cat << _EOF0_ > clsm/README1
   3.3 References
	 Koster, R. and M. Suarez, (1996). Energy and Water Balance Calculations in 
            the Mosaic LSM, NASA Tech. Memo. 104606, Vol. 9. 
	 Prigent, C., C. Jimenez, and J. Catherinot, Comparison of satellite microwave
            backscattering (ASCAT) and visible/near-infrared reflectances (PARASOL) for the
            estimation of aeolian aerodynamic roughness length in arid and semi-arid regions,
            Atmos. Meas. Tech., 5, 2703-2712, 2012
         `echo "${sec3_veg_cite}"`
         
=====================================================================================
====================================== PAGE  4 ======================================
=====================================================================================

4. VEGETATION DYNAMIC DATA

   4.1 Data generation and processing chain

         `echo "${sec4_lai}"`

       The AVHRR NDVI3g data [third generation Global Inventory Modeling and Mapping Studies 
       (GIMMS) Normalizeed difference vegetation index (NDVI) data derived from AVHRR images)] 
       are available two times per month for a 35 year period spanning from 1981 to 2015 at
       5-arcmin resolution (Pinzon et al., 2014). A NDVI climatology data set for the same
       temporal resolution was constructed by temporally averaging over the 35-year period
       on the 5×5 arcmin grid and then aggregating over the pixels of each land element to
       derive NDVI climatology for that land element.

   4.2 Data files and movies
	4.2.1 Greenness Fraction [-] - [Movie 1 : "plots/GREEN.mp4"]
	  file name  : green.dat (or ../green_clim*.data)
	  file format: fortran binaries, little_endian
	      The read statement should be MAPL_ReadForcing compatible. Each data record 
              is preceded by a header containing start and end times of the period that 
              the data have been averaged for. Corresponding start dates for climatological
	      data records are given below:
  	      `echo "${GSWP2_DATES}"`

	      Loop below through until the last data record:
	      read(10) Year_Begin,Month_Begin,Day_Begin,Hour_Begin,Minute_Begin,Secs_Begin,
		Year_End,Month_End,Day_End,Hour_End,Minute_End,Secs_End (Float Numbers)
	      read(10) (data(n),n=1,${NTILES})
	    
	4.2.2 Leaf Area Index (LAI) [m2/m2] - [Movie 2:"plots/LAI.mp4"; Figure 13: "plots/lai.jpg"]
 	  file name  : lai.dat (or ../lai_clim*.data)
	  file format: fortran binaries, little_endian
	      The read statement should be MAPL_ReadForcing compatible. Each data record 
              is preceded by a header containing start and end times of the period that 
              the data have been averaged for. Corresponding start dates for climatological
	      data records are given below:
  	      `echo "${MYLAIDATES}"`

	      Loop below through until the last data record:
	      read(10) Year_Begin,Month_Begin,Day_Begin,Hour_Begin,Minute_Begin,Secs_Begin,
		Year_End,Month_End,Day_End,Hour_End,Minute_End,Secs_End (Float Numbers)
	      read(10) (data(n),n=1,${NTILES})

      	4.2.3 Normalized Difference Vegetation Index (NDVI) [-]
	  file name  : ndvi.dat (or ../ndvi_clim*.data)
	  file format: fortran binaries, little_endian
	      The read statement should be MAPL_ReadForcing compatible. Each data record 
              is preceded by a header containing start and end times of the period that 
              the data have been averaged for. Corresponding start dates for climatological
	      data records are given below:
  	      `echo "${NDVI_DATES}"`

	      Loop below through until the last data record:
	      read(10) Year_Begin,Month_Begin,Day_Begin,Hour_Begin,Minute_Begin,Secs_Begin,
		Year_End,Month_End,Day_End,Hour_End,Minute_End,Secs_End (Float Numbers)
	      read(10) (data(n),n=1,${NTILES})

   4.3 References
	 Dirmeyer, P. and Oki, T. (2002): The Second Global Soil Wetness project (GSWP-2) 
	    Science 2 and Implementation Plan. IGPO Publication Series No. 37, 64p.
         Pinzon, J.E.; Tucker, C.J.	A Non-Stationary 19812012 AVHRR NDVI3g Time Series. 
            Remote Sens. 2014, 6, 6929-6960.	    
         `echo "${sec4_geo_cite}"`

=====================================================================================
====================================== PAGE  5 ======================================
=====================================================================================

5. SURFACE ALBEDO DATA

   5.1 Data generation and processing chain

         The Catchment Land Surface Model (CLSM) computes, at each model time step, 
	 the following quantities: (1) visible (0.3-0.7μm) direct albedo (black sky), 
	 (2) near-infrared (0.7-5.0μm) direct albedo, (3) visible (0.3-0.7μm) diffuse 
	 albedo (white sky), and (4) near-infrared (0.7-5.0μm) diffuse albedo.  
	 Initial diurnally-varying values are first computed using an albedo scheme 
	 (Koster and Suarez 1991) based on the two-stream approximation utilized by 
	 SiB (Sellers et al. 1986).  These values are then scaled so that their 8-day 
	 averages agree with a MODIS-based albedo climatology.  
	 
	 `echo "${sec5_mod}"`

         Meanwhile, the SiB-based albedo scheme was run at a daily time step over a 1-year 
	 period using the vegetation types, greenness fractions, and leaf area indices 
	 established for GEOS-5 for a given distribution of land elements, as described 
	 in sections 3 and 4 above.  Averaging the visible diffuse and near-infrared 
	 diffuse albedos generated by the scheme over 8-day periods produced, in effect, 
	 an 8-day ‘climatology’ of this particular scheme’s diffuse albedos.  The ratio 
	 of the MODIS-based 8-day diffuse visible albedo to the SiB-based diffuse visible 
	 albedo at a given surface element serves as the 8-day ‘scaling factor’ for that 
	 element.  During a full simulation, the time-step values of both the visible 
	 direct and visible diffuse albedos computed with the SiB-based scheme are multiplied 
	 by the diffuse-based scale factor (for the given element and given 8-day period) 
	 prior to being applied to the incoming radiation values.  The same approach is 
	 used to compute the scale factors for the near-infrared albedos.
	                     
   5.2 Data files and movies 
       5.2.1 MODIS Albedo Climatology [Diffused, Visible (0.3_0.7) and Near-Infrared (0.7_5.0);
             Note: GEOS5/CLSM does not read these data
         [Movie 3 : "plots/VISDF.mp4" (Diffused visible) and Movie 4 "plots/NIRDF.mp4" 
                   (Diffused infrared)]
	 file names : ${AlbFileNames}
	 file format: fortran binaries, little_endian
	      The read statement should be MAPL_ReadForcing compatible. Each data record 
              is preceded by a header containing start and end times of the period that 
              the data have been averaged for. Corresponding start dates for climatological
	      data records are given below (MMDD):
  	      `echo "${MYALBDATES}"`

	      Loop below through until the last data record:
	      read(10) Year_Begin,Month_Begin,Day_Begin,Hour_Begin,Minute_Begin,Secs_Begin,
		Year_End,Month_End,Day_End,Hour_End,Minute_End,Secs_End (Float Numbers)
	      read(10) (data(n),n=1,${NTILES})

       5.2.2 MODIS Scale Parameters [Diffused, Visible (0.3_0.7) and Near-Infrared (0.7_5.0)]
	 file names : visdf.dat/nirdf.dat (or ../visdf*dat and ../nirdf*dat)
	 file format: fortran binaries, little_endian
	      The read statement should be MAPL_ReadForcing compatible. Each data record 
              is preceded by a header containing start and end times of the period that 
              the data have been averaged for. Corresponding start dates for climatological
	      data records are given below:
  	      `echo "${MYALBDATES}"`

	      Loop below through until the last data record:
	      read(10) Year_Begin,Month_Begin,Day_Begin,Hour_Begin,Minute_Begin,Secs_Begin,
		Year_End,Month_End,Day_End,Hour_End,Minute_End,Secs_End (Float Numbers)
	      read(10) (data(n),n=1,${NTILES})

   5.3 References
         Gao, F., He, T., Wang, Z., Ghimire, B., Shuai, Y., Masek, J., Schaaf, C. and
	    Williams, C. (2014): Multiscale climatological albedo look-up maps derived 
	    from moderate resolution imaging spectrocradiometer BRDF/albedo products. 
	    Journal of Applied Remote Sensing. 083532-1, Vol. 8
         Koster, R., and M. Suarez, (1991): A simplified treatment of SiB's land surface 
	    albedo parameterization, NASA Tech. Memo. 104538.
         MCD43GF (2014): Gap filled product description. SZN-extended MODIS/Terra+Aqua 
	    30 arc-second global gap-filled snow free BRDF parameters product (Available at:
	    ftp://rsftp.eeos.umb.edu/data02/Gapfilled/readme.pdf, last accessed July 2015.)
         Sellers, P. J., Y. Mintz, Y. C. Sud, and  A. Dalcher, (1986): A simple biosphere 
	    model (SiB) for use within general circulation models, J. Atmos. Sci., 43, 505-531. 

=====================================================================================
====================================== PAGE  6 ======================================
=====================================================================================

6. CATCHMENT LAND SURFACE MODEL PARAMETES

   6.1 Data generation and processing chain

         The Catchment LSM utilizes numerous preprocessed parameters, many of which describe
	 ‘fits’ to the results of highly complex calculations.  For efficiency purposes, 
	 these fits are for used in place of the complex calculations themselves during simulations.
	 The parameters derived for each surface element rely in part on the statistics of 
	 compound topographic index in "cti_stats.dat" (Section 2.2.3) and the soil 
	 hydraulic properties in "soil_param.dat" (Section 2.2.4). The user is referred 
	 to Ducharne et al. (2000) for a description of the CLSM parameter generation process
	 and for definitions of the parameters themselves. 

   6.2 Data files ONLY (NO images)

       In the descriptions below, equation and figure numbers refer to those in Ducharne et al. (2000). 

       6.2.1 Time scale parameters for moisture transfer between 
	 surfexec and rzexec
	 file name : tau_param.dat
         do n = 1, ${NTILES}
		read (10,'(i8,i8,4f10.7)')    &
			tile_index,pfaf_code,atau2,btau2,atau5,btau5
	 end do
	 where:
	 (1) atau2	atau2: Equation (17) for a 2cm surface layer [-]
	 (2) btau2	btau2: Equation (17) for a 2cm surface layer [-]
	 (3) atau5	atau2: Equation (17) for a 5cm surface layer [-]
	 (4) btau5	atau2: Equation (17) for a 5cm surface layer [-]

       6.2.2 Time scale parameters for moisture transfer between
	 root zone and water table
	 file name : ts.dat
   	 do n = 1, ${NTILES}
		read (10,'(i8,i8,f5.2,4(2x,e13.7))')tile_index,pfaf_code,gnu,   &
	             tsa1,tsa2,tsb1,tsb2
      	 end do

 	 where:
	 (1) tsa1	atau1: Equation (16) for +ive root zone excess (Figure 6) [-]
	 (2) tsa2	atau1: Equation (16) for -ive root zone excess (Figure 6) [-] 
	 (3) tsb1	btau1: Equation (16) for +ive root zone excess (Figure 6) [-] 
	 (4) tsb2	btau1: Equation (16) for -ive root zone excess (Figure 6) [-] 

       6.2.3 Baseflow parameters
	 file name : bf.dat
	 do n = 1, ${NTILES}
		read (10,'(i8,i8,f5.2,3(2x,e13.7))')tile_index,pfaf_code,gnu,bf1,bf2,bf3
	 end do

	 where:
	 (1) gnu	 	vertical_transmissivity {Greek nu Equation (8)} [m-1]
	 (2) bf1		A: Equation (9)	[kg m-4]
	 (3) bf2		B: Equation (9)	[m] 
	 (4) bf3         XBAR: Equation (8) [log(m)] == cti_mean in Section 2.2.2

       6.2.4 Area fractioning parameters
	 file name : ar.new
	 do n = 1, ${NTILES}
		read (10,'(i8,i8,f5.2,11(2x,e13.7))')tile_index,pfaf_code,gnu,  &
			ars1,ars2,ars3,ara1,ara2,ara3,ara4,arw1,arw2,arw3,arw4
	 end do

	 where:
	 (1) ars1 	A : Equation (12) for Asat [m+2 kg-1]
	 (2) ars2	B : Equation (12) for Asat [m+2 kg-1]
	 (3) ars3	C : Equation (12) for Asat [m+4 kg-2]
	 (4) ara1	A : Equation (14) of segment1 if skewness < 0.25 [m+2 kg-1] | else straight ara1 = ara3		
	 (5) ara2	B : Equation (14) of segment1 if skewness < 0.25 [-]        | else straight ara2 = ara4		
	 (6) ara3	A : Equation (14) of segment2 if skewness < 0.25 [m+2 kg-1]
	 (7) ara4	B : Equation (14) of segment2 if skewness < 0.25 [-]
	 (8) arw1	A : Equation (12) for THETA0 [m+2 kg-1]
	 (9) arw2	B : Equation (12) for THETA0 [m+2 kg-1]
	 (10)arw3	C : Equation (12) for THETA0 [m+4 kg-2]
	 (11)arw4	Y_infinity : Equation (12) for THETA0 [-]
	
   6.3 References
	Ducharne, A., R. D. Koster, M. J. Suarez, M. Stieglitz, and P. Kumar (2000), 
	    A catchment-based approach to modeling land surface processes in a general 
	    circulation model: 2. Parameter estimation and model demonstration, 
	    J. Geophys. Res., 105(D20), 2482324838, doi:10.1029/2000JD900328.
_EOF0_
if( $MYMASK == GEOS5_10arcsec_mask.nc | $MYMASK == GEOS5_10arcsec_mask | $MYMASK == GEOS5_10arcsec_mask_freshwater-lakes.nc ) then
cat << _EOF1_ > clsm/README2

=====================================================================================
====================================== PAGE  7 ======================================
=====================================================================================

7. GLOBAL RUNOFF ROUTING MODEL DATA

   7.1 Data generation and processing chain

         The Pfafstetter codification (Verdin and Verdin, 1999) assigns a unique multi-digit
	 integer to a given hydrologic catchment within a river basin. The multi-digit 
	 integer, or Pfafstetter code, contains information about connectivity of the 
	 catchment with upstream and downstream catchments; considering all of the catchments’ 
	 codes together allows the construction of a catchment network within the basin. 
	 Verdin (2013) provided global raster arrays of global Level 12 Pfafstetter codes 
	 at 1-arcmin resolution along with information on mean elevation.  These data sets 
	 were used to build the global river channel network slated for use with GEOS-5. 

         The steps used to generate the river network are as follows. Each catchment 
	 (referred to in this discussion as CatchX) has up to two upstream catchments 
	 (and maybe more, as discussed below) but only one downstream catchment (or only 
	 one sink, as an ocean or a lake). The Pfafstetter code helps locate the downstream 
	 catchment. Once it is identified, CatchX automatically becomes one of the upstream 
	 catchments for the downstream catchment. During the Pfafstetter codification process, 
	 sometimes islands can be incorrectly linked across the ocean. Sinks also often get 
	 incorrectly linked with catchments geographically far apart. Thus, a land-water mask 
	 at 1 arc-minute resolution was created using the 10 arc-second mask (Section 2) 
	 created to help identify islands and coastal catchments. Islands, coastal catchments 
	 and sinks were given due attention to ensure that they don’t get linked up incorrectly. 
	 The 1 arc-minute elevation data were used to determine the lowest point (outlet) of 
	 each catchment.

         Once the upstream and downstream catchments for CatchX are identified, the next step
	 is to find the latitudinal and longitudinal coordinates of its upstream and downstream
	 confluences. This is complicated by the fact that the codification at 1-arcmin is 
	 imperfect  in a perfect system, a catchment would only have zero or two upstream 
	 catchments, but because of the discretization to 1-arcmin and the associated loss of 
	 higher resolution information, a catchment might end up having (according to the 
	 discretized codification) one or possibly more than two upstream catchments. To 
	 determine the locations of the upstream confluences, we must therefore consider four 
	 possible cases:

         Case 1: CatchX has no upstream catchments.  In this case, the centroid of CatchX is 
	         used as the location of the upstream confluence, for purposes of computing 
		 river length and elevation difference. 
         Case 2: There is only one upstream catchment.  The nearest 1-arcmin CatchX pixel to 
	         the centroid of the upstream catchment is assumed to be the upstream confluence. 
         Case 3: Exactly 2 upstream catchments are present.  The meeting point of CatchX and 
	         the two upstream catchments is determined and assigned to be the location of 
		 the upstream confluence.
         Case 4: There are more than 2 upstream catchments.  For each 1-arcmin pixel within 
	         CatchX, the sum of the distances between the pixel and the centroids of each 
		 upstream catchment is computed.  The pixel with the minimum sum-of-distances 
		 is assumed to be located at the upstream confluence.

        To determine the location of the downstream confluence, only two cases need to be considered:
        Case 1: There is a downstream catchment.  The location of the downstream confluence of 
	        CatchX is taken to be the same as that of the upstream confluence for that downstream catchment.
        Case 2: There is no downstream catchment. A 1-arcmin CatchX pixel located next to a 
	        water pixel is assumed to be the location of the downstream confluence. 

   7.2 Data files
       7.2.1 Pafafstetter catchment connectivity, channel information
	 file name : /discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/
		    SRTM-TopoData/Pfafcatch-routing.dat
	 read (10,*) NPfafs
	 do n = 1, ${NPfafs}
		read (10,'(i8,i15,4(1x,f9.4),1x,e10.3,3(1x,e9.3),I8,6(1x,f9.4))')         &     
                pfaf_index,pfaf_code,min_lon,max_lon,min_lat,max_lat, mean_elevation,     & 
                cat_area,cum_area, length,ElevDiff, dnst_pfaf_index,DN_long, DN_lat,      &
                UP_lon, UP_lat, mouth_lon, mouth_lat
         end do

         pfaf_index     [-] catchment index (1-$NPfafs) after sorting Pfafstetter codes in ascending order
         pfaf_code      [-] Pfafstetter code of the hydrologic catchment
         min_lon   [degree] Western edge of the catchment
         max_lon   [degree] Eastern edge of the catchment
         min_lat   [degree] Southern edge of the catchment
         max_lat   [degree] Northern edge of the catchment
         mean_elevation [m] area-averaged elevation
         cat_area     [km2] catchment area
         cum_area     [km2] cumulative area upstream of the dwonstream confluence
         length        [km] the length between upstream and downstream confluences
         ElevDiff       [m] the elevation difference between upstream and downstream confluences
         DN_pfaf_index  [-] catchment index  of the downstream catchment (Note : -1 depicts a sink or lake/ocean)
         DN_long   [degree] longitude at the downstream confluence
         DN_lat    [degree] latitude at the downstream confluence
         UP_lon    [degree] longitude at the upstream confluence
         UP_lat    [degree] latitude at the upstream confluence
         mouth_lon [degree] longitude at the river mouth
         mouth_lat [degree] latitude at the river mouth

`echo "${rout_smap}"`
  7.3 References
	 Verdin, K.L., and J.P. Verdin (1999). A topographical system for delineation 
	    and codification of the Earths river basins. J. of Hydrology (218), 1-12.
         Verdin, K (2013): Final Report - High Resolution Topographic Analysis for GMAOs 
	    Catchment LSM, pp 21, available from Global Modeling and Assimilation Office, 
	    610.1 NASA/GSFC.
	    
_EOF1_
endif
cat << _EOF2_ > clsm/README3

=====================================================================================
====================================== PAGE A1 ======================================
=====================================================================================

APPENDIX I - mkCatchParam tag, input options, and log

CVS TAG : $MYTAG

                                  mkCatchParam LOG
                                  ~~~~~~~~~~~~~~~~

_EOF2_

cat << _EOF_ > clsm/back

=====================================================================================
================================ END OF  README FILE ================================
=====================================================================================

_EOF_

sed -e "s/============================================================/ /g"    clsm/mkCatchParam.log > clsm/log
cat clsm/intro clsm/soil clsm/veg1 clsm/veg2 clsm/README1 clsm/README2 clsm/README3 clsm/log clsm/back >> clsm/README

/bin/rm clsm/intro clsm/soil clsm/veg1 clsm/veg2 clsm/README1 clsm/README2 clsm/README3 clsm/log clsm/back

#################################################################################
##  Plotting maps of fixed parameters and making movies of seasonal parameters ##
#################################################################################

mkdir -p clsm/plots
/bin/cp src/clsm_plots.pro clsm/plots/.

cd clsm/plots/

module load tool/idl-8.5

idl  <<EOB

.compile clsm_plots

clsm_plots

exit
EOB

cd $workdir

/bin/rm clsm/plots/clsm_plots.pro

module unload tool/idl-8.5


