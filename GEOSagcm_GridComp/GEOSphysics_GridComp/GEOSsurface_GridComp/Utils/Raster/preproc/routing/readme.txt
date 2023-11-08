v1 11/08/2023
The package is used for creating a file for river outlet locations to ocean used in the GEOS coupled simulation. The outlet locations got from this pre-processing are in the land or glacier area defined by the mk_bcs file Pfafstetter.rst. The output is a binary file Outlet_latlon.43200x21600. This file is the input of the mk_runofftbl.F90 in the mk_bcs package. The inland outlet locations will be further moved to ocean domain by the mk_runofftbl.F90.

If on the Discover, one can simply run the script run.sh to create the output file Outlet_latlon.43200x21600. If not on the Discover, please contact yujin.zeng@nasa.gov .

The function for each f90 code are briefly described as follows:
1. get_outlets_catchindex.f90: get sink catchment IDs.
2. get_outlets_land.f90: get the sink points in land and Greenland (from Lauren Andrews) by picking the point (i.e., 15-sec cell) with the largest drainage area from dataset of HydroSHEDS (https://www.hydrosheds.org/) within each sink catchment.
3. get_sinkxy_land.f90: outlets degree to indexes in the 30 arc-sec map.
4. get_outlets_land_allcat.f90: Assign the outlet locations to all upstream catchments to create a 1d list showing the final x and y indexes for each catchment.
5. get_landocean_Greenland_real.f90: Insert the Greenland index map to the catchment index map.
6. Pfaf_to_2d_30s_land.f90: Transform the 1d list above to a 30s 2d map using the index map.
7. read_riveroutlet_land.f90: Transform above the 2d maps to an unformatted file Outlet_latlon.43200x21600 that can be read directly by the mk_runofftbl.F90 of mk_bcs. 

