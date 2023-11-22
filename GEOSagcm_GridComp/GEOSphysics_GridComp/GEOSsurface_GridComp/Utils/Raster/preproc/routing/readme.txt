v1 11/08/2023, Yujin Zeng

The "preproc/routing" package is used for creating a file with the locations of river outlets to the ocean.

The output from this package is a binary file "Outlet_latlon.43200x21600".  

The river outlets are located in land or landice tiles as defined in the raster file "Pfafstetter.rst" from the makebcs package. 

The "Outlet_latlon.43200x21600" file is the input for "mk_runofftbl.F90" in the makebcs package, which further adjusts the outlet locations to be consistent with the ocean model resolution and domain ("mk_runofftbl.F90").

If on NCCS/Discover, the package can be run using the script "run.sh".  If not on Discover, please contact yujin.zeng@nasa.gov.

The function for each f90 code are briefly described as follows:

1. get_outlets_catchindex.f90:
Get sink catchment IDs.

2. get_outlets_land.f90: 
Get sink points on land or in Greenland (from Lauren Andrews) by picking the point (i.e., 15-arcsec grid cell) within each sink catchment that has the largest drainage area per the HydroSHEDS (https://www.hydrosheds.org/) dataset.

3. get_sinkxy_land.f90: 
Convert outlet locations in degree lat/lon to indices on the 30 arc-sec raster grid.

4. get_outlets_land_allcat.f90: 
Assign outlet locations to all upstream catchments to create a 1d list showing the final x and y indexes for each catchment.

5. get_landocean_Greenland_real.f90: 
Insert the Greenland index map into the catchment index map.

6. Pfaf_to_2d_30s_land.f90: 
Transform the 1d list above to a 30 arc-sec 2d map using the map of indices.

7. read_riveroutlet_land.f90: 
Transform the above 2d maps to the unformatted Fortran binary file "Outlet_latlon.43200x21600" that can be read directly by "mk_runofftbl.F90" of makebcs. 

