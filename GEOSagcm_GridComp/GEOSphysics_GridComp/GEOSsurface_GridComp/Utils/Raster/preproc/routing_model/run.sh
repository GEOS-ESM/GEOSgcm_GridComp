#!/bin/bash
set -e


#---copy some files----
cp input/area_skm_grand.txt output/
cp input/cap_max_grand.txt output/

#---river--------------
echo get_Pfaf_file.f90
./build get_Pfaf_file.f90
./get_Pfaf_file.out

echo get_latloni_cellarea.py
python3 get_latloni_cellarea.py

echo get_num_sub_catchment_M09.f90
./build get_num_sub_catchment_M09.f90
./get_num_sub_catchment_M09.out

echo get_num_sub_catchment_M36.f90
./build get_num_sub_catchment_M36.f90
./get_num_sub_catchment_M36.out

echo get_lonlat_bond_M09.f90
./build get_lonlat_bond_M09.f90
./get_lonlat_bond_M09.out

echo get_lonlat_bond_M36.f90
./build get_lonlat_bond_M36.f90
./get_lonlat_bond_M36.out

echo get_lonlati_maptile_M09.py
python3 get_lonlati_maptile_M09.py
echo get_lonlati_maptile_M36.py
python3 get_lonlati_maptile_M36.py

echo get_isub_M09.f90
./build get_isub_M09.f90 
./get_isub_M09.out

echo get_isub_M36.f90
./build get_isub_M36.f90
./get_isub_M36.out

echo get_area_M09.f90
./build get_area_M09.f90
./get_area_M09.out

echo get_area_M36.f90
./build get_area_M36.f90
./get_area_M36.out

echo get_Qr_clmt.f90
./build get_Qr_clmt.f90
./get_Qr_clmt.out

echo get_river_length.f90
./build get_river_length.f90
./get_river_length.out

echo get_K_model_calik.f90
./build get_K_model_calik.f90
./get_K_model_calik.out

#--------reservoir-----------
echo get_dam_data.py
python3 get_dam_data.py

#--------lake----------------
echo read_input_TopoCat.f90
./build read_input_TopoCat.f90
./read_input_TopoCat.out

echo process_lake_data.py
python3 process_lake_data.py

