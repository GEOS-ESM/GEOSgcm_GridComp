#!/usr/bin/env python3

#source g5_modules before run to get the necessary env

import os
import shutil
import subprocess
from pathlib import Path

def run(cmd):
    """Run a command and exit on failure."""
    print(f">>> {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def main():
    # -----------------------------
    # Define input and output paths
    # -----------------------------

# Path to "bcs_shared" directory:
    file_path  = "/discover/nobackup/yzeng3/data/river_preproc_input/"     # NCCS Discover
    #file_path = "/nobackup/gmao_SIteam/ModelData/bcs_shared/"      # NAS  


    file_pfafrout       = file_path + "/Pfafcatch-routing.dat"
    file_lat1m          = file_path + "/lat_1m.txt"
    file_lon1m          = file_path + "/lon_1m.txt"
    file_lat            = {"M09": file_path + "/lat_M09.txt",
                           "M36": file_path + "/lat_M36.txt"}
    file_lon            = {"M09": file_path + "/lon_M09.txt",
                           "M36": file_path + "/lon_M36.txt"}
    file_catdef         = {"M09": file_path + "/catchment_M09.def",
                           "M36": file_path + "/catchment_M36.def"}
    file_clmtrunf       = file_path + "/SMAPL4_OL7000_runoff_mean_2016_2023.nc"
    file_pfafmap        = file_path + "/SRTM_PfafData.nc"
    file_ldn            = file_path + "/hyd_glo_ldn_15s.nc"
    file_hyelev         = file_path + "/hyd_glo_dem_15s.nc"
    file_vel            = file_path + "/velocity.txt"
    file_dis            = file_path + "/discharge.txt"
    file_usid           = file_path + "/USGSID.txt"
    file_lats           = file_path + "/lat_for_site.txt"
    file_lons           = file_path + "/lon_for_site.txt"
    file_gage_id        = file_path + "/id_gagesii.txt"
    file_gage_acar      = file_path + "/acar_gagesii.txt"

    file_latdam         = file_path + "/lat_dam_grand.txt"
    file_londam         = file_path + "/lon_dam_grand.txt"
    file_acadam         = file_path + "/catch_acar_grand.txt"
    file_damcat_manfix  = file_path + "/catid_dam_manfix.txt"
    file_dam_manflag    = file_path + "/flag_dam_manfix.txt"
    file_dam_use        = file_path + "/main_use_grand.txt"
    file_damflood       = file_path + "/flood_use_grand.txt"
    file_damarea        = file_path + "/area_skm_grand.txt"
    file_damcap         = file_path + "/cap_max_grand.txt"

    file_lake_area      = file_path + "/Lake_area.csv"
    file_lake_id        = file_path + "/Hylak_id_lake.csv"
    file_lake_aca       = file_path + "/acar_lake.csv"
    file_lakeo_lakeid   = file_path + "/Hylak_id_lakeout.csv"
    file_lakeo_id       = file_path + "/Outlet_id_lakeout.csv"
    file_lakeo_lat      = file_path + "/Outlet_lat_lakeout.csv"
    file_lakeo_lon      = file_path + "/Outlet_lon_lakeout.csv"
    file_lake_mantag    = file_path + "/catid_lake_manfix.txt"
    file_lakecat_manfix = file_path + "/catid_lake_multout_manfix.txt"


    lib_path = "/discover/nobackup/yzeng3/lib"
    old_ld = os.environ.get("LD_LIBRARY_PATH", "")
    os.environ["LD_LIBRARY_PATH"] = f"{lib_path}:{old_ld}"
    
    # -----------------------------
    # Ensure output and temp directory exists
    # -----------------------------
    subprocess.run(["mkdir", "-p", "output"], check=True)
    subprocess.run(["mkdir", "-p", "temp"], check=True) 

    # -----------------------------
    # Copy dam area and capacity files
    # -----------------------------
    shutil.copy(file_damarea, "output/")
    shutil.copy(file_damcap,  "output/")

    # -----------------------------
    # River processing section
    # -----------------------------
    # Compile and run Pfafcatch routing generator
    run(["./build", "get_Pfaf_file.f90"])
    run(["./get_Pfaf_file.out", file_pfafrout])

    # Generate latitude/longitude indices and cell areas
    run([
        "python3", "get_latloni_cellarea.py",
        file_lat["M36"], file_lon["M36"],
        file_lat["M09"], file_lon["M09"],        
        file_lat1m, file_lon1m,
    ])

    # Compute number of sub-catchments for M09 and M36 resolutions
    for res in ("M09", "M36"):
        run(["./build", f"get_num_sub_catchment_{res}.f90"])
        run([f"./get_num_sub_catchment_{res}.out", file_pfafmap])

    # Build longitude-latitude boundary files for each resolution
    for res in ("M09", "M36"):
        run(["./build", f"get_lonlat_bond_{res}.f90"])
        run([f"./get_lonlat_bond_{res}.out", file_catdef[res]])

    # Map tile longitude/latitude for M09 and M36
    for res in ("M09", "M36"):    
        run(["python3", f"get_lonlati_maptile_{res}.py", file_lat[res], file_lon[res]])
    # Build and run isub calculators for both resolutions
    for res in ("M09", "M36"):
        run(["./build", f"get_isub_{res}.f90"])
        run([f"./get_isub_{res}.out"])

    # Calculate area of each catchment
    for res in ("M09", "M36"):
        run(["./build", f"get_area_{res}.f90"])
        run([f"./get_area_{res}.out", file_pfafmap])    

    # Compute climatological runoff
    run(["./build", "get_Qr_clmt.f90"])
    run(["./get_Qr_clmt.out", file_clmtrunf])

    # Determine river lengths
    run(["./build", "get_river_length.f90"])
    run([
        "./get_river_length.out",
        file_pfafmap, file_ldn,
        file_hyelev, file_pfafrout
    ])

    # Calibrate K model using velocity and discharge data
    run(["./build", "get_K_model_calik.f90"])
    run([
        "./get_K_model_calik.out",
        file_vel, file_dis, file_usid,
        file_lats, file_lons,
        file_lat1m, file_lon1m,
        file_pfafmap,
        file_gage_id, file_gage_acar
    ])

    # -----------------------------
    # Reservoir (dam) processing
    # -----------------------------
    run([
        "python3", "get_dam_data.py",
        file_latdam, file_londam,
        file_lat1m, file_lon1m,
        file_pfafmap, file_acadam,
        file_damcat_manfix, file_dam_manflag,
        file_dam_use, file_damflood
    ])

    # -----------------------------
    # Lake processing section
    # -----------------------------
    run(["./build", "read_input_TopoCat.f90"])
    run([
        "./read_input_TopoCat.out",
        file_lake_area, file_lake_id,
        file_lake_aca,
        file_lakeo_lakeid, file_lakeo_id,
        file_lakeo_lat, file_lakeo_lon
    ])

    run([
        "python3", "process_lake_data.py",
        file_lat1m, file_lon1m,
        file_pfafmap,
        file_lake_mantag, file_lakecat_manfix
    ])

if __name__ == "__main__":
    main()