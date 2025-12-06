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
    file_path1  = "/discover/nobackup/projects/gmao/bcs_shared/"
    file_path2  = "/discover/nobackup/yzeng3/data/river_preproc_input/" 


    file_pfafrout       = file_path1 + "/make_bcs_inputs/land/topo/v1/SRTM-TopoData/Pfafcatch-routing.dat"
    file_pfafmap        = file_path1 + "/make_bcs_inputs/land/topo/v1/SRTM-TopoData/SRTM_PfafData.nc"
    file_catdef         = file_path1 + "/fvInput/ExtData/esm/tiles/v14/land/EASEv2_M09/clsm/catchment.def"

    file_lat1m          = file_path2 + "/lat_1m.txt"
    file_lon1m          = file_path2 + "/lon_1m.txt"
    file_lat            = file_path2 + "/lat_SMAPL4.txt"
    file_lon            = file_path2 + "/lon_SMAPL4.txt"
    file_clmtrunf       = file_path2 + "/SMAPL4_OL7000_runoff_mean_2016_2023.nc"
    file_ldn            = file_path2 + "/hyd_glo_ldn_15s.nc"
    file_hyelev         = file_path2 + "/hyd_glo_dem_15s.nc"
    file_vel            = file_path2 + "/velocity.txt"
    file_dis            = file_path2 + "/discharge.txt"
    file_usid           = file_path2 + "/USGSID.txt"
    file_lats           = file_path2 + "/lat_for_site.txt"
    file_lons           = file_path2 + "/lon_for_site.txt"
    file_gage_id        = file_path2 + "/id_gagesii.txt"
    file_gage_acar      = file_path2 + "/acar_gagesii.txt"

    file_latdam         = file_path2 + "/lat_dam_grand.txt"
    file_londam         = file_path2 + "/lon_dam_grand.txt"
    file_acadam         = file_path2 + "/catch_acar_grand.txt"
    file_damcat_manfix  = file_path2 + "/catid_dam_manfix.txt"
    file_dam_manflag    = file_path2 + "/flag_dam_manfix.txt"
    file_dam_use        = file_path2 + "/main_use_grand.txt"
    file_damflood       = file_path2 + "/flood_use_grand.txt"
    file_damarea        = file_path2 + "/area_skm_grand.txt"
    file_damcap         = file_path2 + "/cap_max_grand.txt"

    file_lake_area      = file_path2 + "/Lake_area.csv"
    file_lake_id        = file_path2 + "/Hylak_id_lake.csv"
    file_lake_aca       = file_path2 + "/acar_lake.csv"
    file_lakeo_lakeid   = file_path2 + "/Hylak_id_lakeout.csv"
    file_lakeo_id       = file_path2 + "/Outlet_id_lakeout.csv"
    file_lakeo_lat      = file_path2 + "/Outlet_lat_lakeout.csv"
    file_lakeo_lon      = file_path2 + "/Outlet_lon_lakeout.csv"
    file_lake_mantag    = file_path2 + "/catid_lake_manfix.txt"
    file_lakecat_manfix = file_path2 + "/catid_lake_multout_manfix.txt"


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
    shutil.copy(file_damarea, "temp/")
    shutil.copy(file_damcap,  "temp/")

    # -----------------------------
    # River processing section
    # -----------------------------
    # Compile and run Pfafcatch routing generator
    run(["./build", "get_Pfaf_file.f90"])
    run(["./get_Pfaf_file.out", file_pfafrout])

    # Generate latitude/longitude indices and cell areas
    run([
        "python3", "get_latloni_cellarea.py",
        file_lat, file_lon,        
        file_lat1m, file_lon1m,
    ])

    # Compute number of sub-catchments
    run(["./build", f"get_num_sub_catchment.f90"])
    run([f"./get_num_sub_catchment.out", file_pfafmap])

    # Build longitude-latitude boundary files for each resolution
    run(["./build", f"get_lonlat_bond.f90"])
    run([f"./get_lonlat_bond.out", file_catdef])

    # Map tile longitude/latitude 
    run(["python3", f"get_lonlati_maptile.py", file_lat, file_lon])
    # Build and run isub calculators
    run(["./build", f"get_isub.f90"])
    run([f"./get_isub.out"])

    # Calculate area of each catchment
    run(["./build", f"get_area.f90"])
    run([f"./get_area.out", file_pfafmap])    

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

    run([
        "python3", "create_river_input.py",
    ])

    subprocess.run(["rm", "-rf", "temp"], check=True)
    subprocess.run("rm -f *.out *.mod", shell=True, check=True)

if __name__ == "__main__":
    main()