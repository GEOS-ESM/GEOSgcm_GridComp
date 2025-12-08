v1 12/05/2025, Yujin Zeng

The "preproc/routing_model" package is used for creating input data related to the GEOS routing model that would be used in the make_bcs.

Usage:
  • On NCCS/Discover:
      1. source g5_modules.sh
      2. python3 run_routing_preproc.py under the build/src directory
  • Off Discover: contact yujin.zeng@nasa.gov

The tasks completed by each F90 or Python program are briefly described as follows:

1. get_Pfaf_file.f90  
   Reads the Pfafstetter code dataset and generates  
   files for the connectivity of catchments in the routing network.

2. get_latloni_cellarea.py  
   Computes grid-cell index arrays and per-cell areas for 1-m high-res grid.

3. get_num_sub_catchment.f90
   Parses high-res map of catchment index to get the area and 
   coordinates (of model grid) of each sub-catchments within each main catchment.

4. get_lonlat_bond.f90
   Extracts the latitude/longitude boundaries of each catchment-tile from  
   catchment definition files.

5. get_lonlati_maptile.py
   Assigns a catchment‐tile index from catchment definition files to each model grid cell.

6. get_isub.f90
   Assigns a catchment‐tile index from maptile files to each sub-catchment.

7. get_area.f90 
   Gets the area for each catchment-tile.

8. get_Qr_clmt.f90  
   Reads SMAP L4 runoff data (2016–2023) from a NetCDF file and computes the climatological  
   mean discharge for each catchment.

9. get_river_length.f90  
   Determines main river channel lengths for each catchment by using HydroSHEDS
   data of distance to sink.

10. get_K_model_calik.f90  
    Calculates the K parameter used in the river routing model.

11. get_dam_data.py  
    Processes reservoir (dam) data: reads dam locations and usage information from GRanD database.

12. read_input_TopoCat.f90  
    Reads lake and lake outlets information from Lake-TopoCa database.

13. process_lake_data.py  
    Processes lake data to be used in the river routing model.

14. create_river_input.py
    Combines all inputs created above to the one nc file river_input.nc

The explanations for the input files of this package can be found in the input directory.

The outputs from this package (cellarea.nc and river_input.nc) will be used as input to the make_bcs. They are listed as follows:

  cellarea.nc: cell area [m^2] of the 1-min grid.

  river_input.nc: River parameters on the catchment space. More detail can be seen in the attributes of each variable in the netcdf file.









