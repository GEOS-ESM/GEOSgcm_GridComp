v1 12/05/2025, Yujin Zeng

The "preproc/routing_model" package is used for creating input data related to the GEOS routing model that would be used in the make_bcs.

Usage:
  • On NCCS/Discover:
      1. source g5_modules.sh
      2. build GEOSldas code
      3. python3 run_routing_preproc.py under the build/src directory
  • Off Discover: contact yujin.zeng@nasa.gov

The tasks completed by each F90 or Python program are briefly described as follows:

1. get_Pfaf_file.f90  
   Reads the Pfafstetter code dataset and generates  
   files for the connectivity of catchments in the routing network.

2. get_latloni.py  
   Computes grid-cell index arrays for 1-m high-res grid.

3. get_num_sub_catchment.f90
   Parses high-res map of catchment index to get the area and 
   coordinates (of model grid) of each sub-catchments within each main catchment.

4. get_Qr_clmt.f90  
   Reads SMAP L4 runoff data (2016–2023) from a NetCDF file and computes the climatological  
   mean discharge for each catchment.

5. get_river_length.f90  
   Determines main river channel lengths for each catchment by using HydroSHEDS
   data of distance to sink.

6. get_K_model_calik.f90  
    Calculates the K parameter used in the river routing model.

7. get_dam_data.py  
    Processes reservoir (dam) data: reads dam locations and usage information from GRanD database.

8. read_input_TopoCat.f90  
    Reads lake and lake outlets information from Lake-TopoCa database.

9. process_lake_data.py  
    Processes lake data to be used in the river routing model.

10. create_river_input.py
    Combines all inputs created above to the one nc file route_parameters.nc

The explanations for the input files of this package can be found in the input directory.

The outputs from this package (route_parameters.nc located in the created output folder) will be used as input to the make_bcs. 

route_parameters.nc: River parameters on the catchment space. More detail can be seen in the attributes of each variable in the netcdf file.









