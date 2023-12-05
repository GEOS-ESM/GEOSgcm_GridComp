#!/usr/bin/env python3

import os
import subprocess

input_path = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing"
baselib = "/discover/swdev/gmao_SIteam/Baselibs/ESMA-Baselibs-7.8.1/x86_64-pc-linux-gnu/ifort_2021.6.0-intelmpi_2021.6.0/Linux"

# Remove files and directories
os.system("rm -rf inputs >& /dev/null")
os.system("rm -rf outputs >& /dev/null")
os.system("rm -f *.mod >& /dev/null")
os.system("rm -f *.out >& /dev/null")
os.system("rm -f Outlet_latlon.43200x21600 >& /dev/null")

# Link basedir
base_lib_path = baselib + "/lib"
for file in os.listdir(base_lib_path):
  os.system(f"rm -f {file} >& /dev/null")
  os.symlink(os.path.join(base_lib_path, file), os.path.join(os.getcwd(), file))
base_inc_path = baselib + "/include/netcdf"
for file in os.listdir(base_inc_path):
  os.system(f"rm -f {file} >& /dev/null")
  os.symlink(os.path.join(base_inc_path, file), os.path.join(os.getcwd(), file))

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
    subprocess.run(f"ifort constant.f90 {program}.f90 -lnetcdf -lnetcdff -o {program}.out",shell=True)

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
base_lib_path = os.getenv("BASELIB") + "/lib/"
for file in os.listdir(base_lib_path):
  os.system(f"rm -f {file} >& /dev/null")
base_inc_path = os.getenv("BASELIB") + "/include/netcdf"
for file in os.listdir(base_inc_path):
  os.system(f"rm -f {file} >& /dev/null")

