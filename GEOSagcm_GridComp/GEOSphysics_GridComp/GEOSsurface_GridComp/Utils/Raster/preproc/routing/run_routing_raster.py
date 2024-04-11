#!/usr/bin/env python3

#source g5_modules before run to get the necessary env

import os
import subprocess

#input_path = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/"

file_Pfafcatch="/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/topo/v1/SRTM-TopoData/Pfafcatch-routing.dat"
file_SRTMPfaf="/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/topo/v1/SRTM-TopoData/SRTM_PfafData.nc"
file_Drainage="/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/HydroSHEDS_drainage_area.nc"
file_GrnLat="/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/Greenland_outlets_lat.txt"
file_GrnLon="/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/Greenland_outlets_lon.txt"
file_GrnMap="/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/GreenlandID_30s.nc"


# Remove files and directories
os.system("rm -rf outputs >& /dev/null")
os.system("rm -f Outlet_latlon.43200x21600 >& /dev/null")

# Create directories and symbolic links
os.makedirs("outputs", exist_ok=True)


programs_inputs = {
    "get_finalID_msk.x":[file_Pfafcatch],
    "get_outlets_catchindex.x":[],
    "get_outlets_land.x":[file_SRTMPfaf,file_Drainage,file_GrnLat,file_GrnLon],
    "get_sinkxy_land.x":[],
    "get_outlets_land_allcat.x":[],
    "get_landocean_Greenland_real.x":[file_SRTMPfaf,file_GrnMap],
    "Pfaf_to_2d_30s_land.x":[]
}


for program, input_files in program_inputs.items():
    if input_files: 
        for input_file in input_files:
            command = [program, input_file]
            subprocess.run(command)
    else:
        command = [program]
        subprocess.run(command)

#for out_program in out_programs:
#    print(f"running {out_program}")
#    subprocess.run(f"./{out_program}",shell=True)

print("Outlet_latlon.43200x21600 created!")

# Clean up
print("Removing temporary output files ...")
os.system("rm -rf outputs")

