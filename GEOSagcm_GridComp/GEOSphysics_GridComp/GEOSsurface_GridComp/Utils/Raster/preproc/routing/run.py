#!/usr/bin/env python3

import os
import subprocess

input_path = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/routing"
#netcdf_path = "/usr/local/other/netcdf4/4.1.2/gcc-4.8.5"
install_path = "../../../../../../../../../../../install/bin"

# Remove files and directories
os.system("rm -rf inputs >& /dev/null")
os.system("rm -rf outputs >& /dev/null")
#os.system("rm -f *.mod >& /dev/null")
#os.system("rm -f *.out >& /dev/null")
os.system("rm -f *.x >& /dev/null")
os.system("rm -f Outlet_latlon.43200x21600 >& /dev/null")

# Create directories and symbolic links
os.makedirs("inputs", exist_ok=True)
os.makedirs("outputs", exist_ok=True)
for file in os.listdir(input_path):
    os.symlink(os.path.join(input_path, file), os.path.join("inputs", file))

# Link and run Fortran programs
programs = [
    "get_outlets_catchindex",
    "get_outlets_land",
    "get_sinkxy_land",
    "get_outlets_land_allcat",
    "get_landocean_Greenland_real",
    "Pfaf_to_2d_30s_land",
]

#for program in programs:
#    print(f"Building {program} ...")
    #subprocess.run(["./build", program])
#    subprocess.run(f"gfortran constant.f90 {program}.f90 -I{netcdf_path}/include -L{netcdf_path}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lcurl -lz -lsz -ldl -o {program}.out",shell=True)
current_working_directory = os.getcwd()
for program in programs:
    os.symlink(os.path.join(install_path, program+".x"), os.path.join(current_working_directory, program+".x"))

out_programs = [
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
print("Removing *.out files ...")
#os.system("rm -f *.out")
os.system("rm -f *.x")
#os.system("rm -f *.mod")

