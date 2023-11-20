#!/usr/bin/env python3

import os
import subprocess

input_path = "/discover/nobackup/yzeng3/work/outlets/inputs"
netcdf_path = "/discover/nobackup/yzeng3/apps/netcdf-4.2.1.1"

# Remove files and directories
os.system("rm -rf inputs >& /dev/null")
os.system("rm -rf outputs >& /dev/null")
os.system("rm -f *.mod >& /dev/null")
os.system("rm -f *.out >& /dev/null")
os.system("rm -f Outlet_latlon.43200x21600 >& /dev/null")

# Create directories and symbolic links
os.makedirs("inputs", exist_ok=True)
os.makedirs("outputs", exist_ok=True)
for file in os.listdir(input_path):
    os.symlink(os.path.join(input_path, file), os.path.join("inputs", file))

# Build and run Fortran programs
programs = [
    "get_outlets_catchindex",
    "get_outlets_land",
    "get_sinkxy_land",
    "get_outlets_land_allcat",
    "get_landocean_Greenland_real",
    "Pfaf_to_2d_30s_land",
]

for program in programs:
    print(f"Building {program} ...")
    #subprocess.run(["./build", program])
    subprocess.run(f"ifort constant.f90 {program}.f90 -I{netcdf_path}/include -L{netcdf_path}/lib -lnetcdf -lnetcdff -o {program}.out",shell=True)

out_programs = [
    "get_outlets_catchindex.out",
    "get_outlets_land.out",
    "get_sinkxy_land.out",
    "get_outlets_land_allcat.out",
    "get_landocean_Greenland_real.out",
    "Pfaf_to_2d_30s_land.out",
]

for out_program in out_programs:
    print(f"running {out_program}")
    subprocess.run(f"./{out_program}",shell=True)

print("Outlet_latlon.43200x21600 created!")

# Clean up
print("Removing temporary input/output files ...")
os.system("rm -rf outputs")
os.system("rm -rf inputs")
print("Removing *.out files ...")
os.system("rm -f *.out")
os.system("rm -f *.mod")


