#!/usr/bin/env python3

#source g5_modules before run to get the necessary env

import os
import subprocess

input_path = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing/"

# Remove files and directories
os.system("rm -rf inputs >& /dev/null")
os.system("rm -rf outputs >& /dev/null")
os.system("rm -f Outlet_latlon.43200x21600 >& /dev/null")

# Create directories and symbolic links
os.makedirs("inputs", exist_ok=True)
os.makedirs("outputs", exist_ok=True)
for file in os.listdir(input_path):
    os.symlink(os.path.join(input_path, file), os.path.join("inputs", file))

out_programs = [
    "get_finalID_msk.x",
    "get_outlets_catchindex.x",
    "get_outlets_land.x",
    "get_sinkxy_land.x",
    "get_outlets_land_allcat.x",
    "get_landocean_Greenland_real.x",
    "Pfaf_to_2d_30s_land.x",
]

for out_program in out_programs:
    print(f"running {out_program}")
    subprocess.run(f"./{out_program}",shell=True)

print("Outlet_latlon.43200x21600 created!")

# Clean up
print("Removing temporary input/output files ...")
os.system("rm -rf outputs")
os.system("rm -rf inputs")

