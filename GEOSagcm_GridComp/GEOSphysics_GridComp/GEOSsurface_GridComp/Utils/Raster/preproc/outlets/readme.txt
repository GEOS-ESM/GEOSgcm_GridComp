v1 11/08/2023, Yujin Zeng

The "preproc/routing" package is used for creating a 30-arcsec raster file with the locations of river outlets to the ocean.

The output from this package is a binary file "Outlet_latlon.43200x21600".  

The river outlets are located in land or landice tiles as defined in the file "SRTM_PfafData.nc". 

The "Outlet_latlon.43200x21600" file is the input for "mk_runofftbl.F90" in the makebcs package, which further adjusts the outlet locations to be consistent with the ocean model resolution and domain ("mk_runofftbl.F90").

If on NCCS/Discover, the package can be run using the script "run_routing_raster.py". Users should source g5_modules before run to get the necessary env. If not on Discover, please contact yujin.zeng@nasa.gov.

The tasks completed by each f90 program are briefly described as follows:

1. get_finalID_msk.f90:
Get downstream catchment and final destination ID for each catchment, and determine whether it directs to ocean or inland lake. 

2. get_outlets_catchindex.f90:
Get a list of sink catchment IDs.

3. get_outlets_land.f90: 
Get sink points on land or in Greenland (from Lauren Andrews) by picking the point (i.e., 15-arcsec grid cell) within each sink catchment that has the largest drainage area per the HydroSHEDS (https://www.hydrosheds.org/) dataset.

4. get_sinkxy_land.f90: 
Convert outlet locations in degree lat/lon to indices on the 30 arc-sec raster grid.

5. get_outlets_land_allcat.f90: 
Assign outlet locations to all upstream catchments to create a 1d list showing the final x and y indexes for each catchment.

6. get_landocean_Greenland_real.f90: 
Insert the Greenland index map into the catchment index map.

7. Pfaf_to_2d_30s_land.f90: 
Transform the 1d list above to the unformatted Fortran binary file "Outlet_latlon.43200x21600" that can be read directly by "mk_runofftbl.F90" of makebcs.

