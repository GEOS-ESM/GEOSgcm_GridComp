v1 05/05/2025, Yujin Zeng

The "preproc/routing_model" package is used for creating input data for the GEOS routing model.

Usage:
  • On NCCS/Discover:
      1. source g5_modules.sh
      2. python3 run_routing_preproc.py
  • Off Discover: contact yujin.zeng@nasa.gov

The tasks completed by each F90 or Python program are briefly described as follows:

1. get_Pfaf_file.f90  
   Reads the Pfafstetter code dataset and generates  
   files for the connectivity of catchments in the routing network.

2. get_latloni_cellarea.py  
   Computes grid-cell index arrays and per-cell areas for 1-m high-res grid.

3. get_num_sub_catchment_M09.f90 / get_num_sub_catchment_M36.f90  
   Parses high-res map of catchment index to get the area and 
   coordinates (of model grid) of each sub-catchments within each main catchment.

4. get_lonlat_bond_M09.f90 / get_lonlat_bond_M36.f90  
   Extracts the latitude/longitude boundaries of each catchment-tile from  
   catchment definition files.

5. get_lonlati_maptile_M09.py / get_lonlati_maptile_M36.py  
   Assigns a catchment‐tile index from catchment definition files to each model grid cell.

6. get_isub_M09.f90 / get_isub_M36.f90  
   Assigns a catchment‐tile index from maptile files to each sub-catchment.

7. get_area_M09.f90 / get_area_M36.f90  
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

The explanations for the input files of this package can be found in the input directory.

The outputs from this package, which are used as input to the river routing model, are listed as follows:

downstream_1D_new_noadj.txt  
  Downstream catchment id for each catchment.

Pfaf_area.txt  
  Catchment area (km^2) for each catchment.

upstream_1D.txt  
  Upstream catchment id for each catchment.

Pfaf_upnum.txt  
  Number of upstream catchments for each catchment.

Pfaf_tosink.txt  
  Number of steps to final sink for each catchment.

Pfaf_nsub_M09.txt / Pfaf_nsub_M36.txt  
  Count of sub‐catchments contained within each catchment at M09 and M36 resolutions, respectively.

Pfaf_xsub_M09.txt / Pfaf_xsub_M36.txt  
  X (longitude) coordinates in M09 and M36 grid for each sub-catchemnt within each catchment, respectively.

Pfaf_ysub_M09.txt / Pfaf_ysub_M36.txt  
  Y (longitude) coordinates in M09 and M36 grid for each sub-catchemnt within each catchment, respectively.

Pfaf_asub_M09.txt / Pfaf_asub_M36.txt  
  Area (km^2) of each sub‐catchment within each catchment at M09 and M36 resolutions, respectively.

Pfaf_isub_M09.txt / Pfaf_isub_M36.txt  
  Tile number (in the catchment definition file) of each sub‐catchment within each catchment at M09 and M36 resolutions, respectively.

area_M09_1d.txt / area_M36_1d.txt  
  Area (km^2) of each tile (in the catchment definition file) for M09 and M36 grid, respectively.

Pfaf_qstr.txt  
  Climatological mean runoff (m^3 s-1) for each catchment.

Pfaf_qri.txt  
  Climatological mean discharge (m^3 s-1) for each catchment.

Pfaf_qin.txt  
  Climatological mean inflow (m^3 s-1) from upstream for each catchment.  

Pfaf_lriv_PR.txt  
  Main river length scale (km) for each catchment. 

Pfaf_lstr_PR.txt  
  Mean local stream length scale (km) for each catchment.

Pfaf_Kstr_PR_fac1_0p35_0p45_0p2_n0p2.txt  
  Calculated K parameters for local streams in each catchment.  

Pfaf_Kv_PR_0p35_0p45_0p2_n0p2.txt  
  Calculated K parameters for main rivers in each catchment. 

area_skm_grand.txt  
  Reservoir surface areas (km^2) for the GRanD dams. 

cap_max_grand.txt  
  Maximum storage capacities (10^6 m^3) for the GRanD dams. 

catid_dam_corr_aca_grand5000.txt  
  Catchment IDs for the GRanD dams.

flag_all_res.txt  
  In-use flags for the GRanD dams.

irr_grand.txt  
  Flags for irrigation use for the GRanD dams.

hydroelec_grand.txt  
  Flags for hydroelectric use for the GRanD dams.

watersupply_grand.txt  
  Flags for water-supply use for the GRanD dams.

nav_grand.txt  
  Flags for navigation use for the GRanD dams.

rec_grand.txt  
  Flags for recreational use for the GRanD dams.

fldmainsec_grand.txt  
  Flags for flood‐control use for the GRanD dams.

other_grand.txt  
  Flags for other use for the GRanD dams.

lake_outlet_lakearea.txt  
  Lake surface area (km^2) for each lake represented in the model.

lake_outlet_flag_valid_2097.txt  
  In-use flags for the lakes.

lake_outlet_catid.txt  
  Catchment IDs for the lakes represented in the model.










