#!/usr/bin/env python3

#source g5_modules before run to get the necessary env

import os
import subprocess

# Input files
# For NAS users, please replace the Discover base path "/discover/nobackup/projects/gmao/bcs_shared/"
# by the NAS path "/nobackup/gmao_SIteam/ModelData/bcs_shared/"
file_Pfafcatch="/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/topo/v1/SRTM-TopoData/Pfafcatch-routing.dat"
file_SRTMPfaf="/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/topo/v1/SRTM-TopoData/SRTM_PfafData.nc"
file_Drainage="/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/HydroSHEDS_drainage_area.nc"
file_GrnLat="/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/Greenland_outlets_lat.txt"
file_GrnLon="/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/Greenland_outlets_lon.txt"
file_GrnMap="/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/GreenlandID_30s.nc"

name_Pfafcatch = os.path.basename(file_Pfafcatch)
name_SRTMPfaf = os.path.basename(file_SRTMPfaf)
name_Drainage = os.path.basename(file_Drainage)
name_GrnLat =os.path.basename(file_GrnLat)
name_GrnLon = os.path.basename(file_GrnLon)
name_GrnMap = os.path.basename(file_GrnMap)

files=[file_Pfafcatch,file_SRTMPfaf,file_Drainage,file_GrnLat,file_GrnLon,file_GrnMap]
# Remove files and directories
os.system("rm -rf inputs >& /dev/null")
os.system("rm -rf outputs >& /dev/null")
os.system("rm -f Outlet_latlon.43200x21600 >& /dev/null")

# Create directories and symbolic links
os.makedirs("inputs", exist_ok=True)
os.makedirs("outputs", exist_ok=True)
for file in files:
    file_name = os.path.basename(file)
    os.symlink(file, os.path.join("inputs", file_name))

program_inputs = {
    "./get_finalID_msk.x": [f"inputs/{name_Pfafcatch}"],
    "./get_outlets_catchindex.x": [],
    "./get_outlets_land.x": [f"inputs/{name_SRTMPfaf}", f"inputs/{name_Drainage}", f"inputs/{name_GrnLat}", f"inputs/{name_GrnLon}"],
    "./get_sinkxy_land.x": [],
    "./get_outlets_land_allcat.x": [],
    "./get_landocean_Greenland_real.x": [f"inputs/{name_SRTMPfaf}", f"inputs/{name_GrnMap}"],
    "./Pfaf_to_2d_30s_land.x": []
}

for program, input_files in program_inputs.items():
    print(f"running {program}")
    command = [program] + input_files
    subprocess.run(command)

print("Outlet_latlon.43200x21600 created!")

# Clean up
print("Removing temporary output files ...")
os.system("rm -rf outputs")
os.system("rm -rf inputs")
