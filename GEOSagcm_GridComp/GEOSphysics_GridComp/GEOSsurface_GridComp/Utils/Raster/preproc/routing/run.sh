#!/bin/bash
set -e

INPUT=/discover/nobackup/yzeng3/work/outlets/inputs

module load comp/intel/2021.3.0

rm -rf inputs >& /dev/null
rm -rf outputs >& /dev/null
rm -f *.mod >& /dev/null
rm -f *.out >& /dev/null
rm -f Outlet_latlon.43200x21600 >& /dev/null 

mkdir -p inputs outputs
ln -s ${INPUT}/* inputs 

echo "Building get_outlets_catchindex.f90 ..."
./build get_outlets_catchindex.f90
echo "Building get_outlets_land.f90 ..."
./build get_outlets_land.f90
echo "Building get_sinkxy_land.f90 ..."
./build get_sinkxy_land.f90
echo "Building get_outlets_land_allcat.f90 ..."
./build get_outlets_land_allcat.f90
echo "Building get_landocean_Greenland_real.f90 ..."
./build get_landocean_Greenland_real.f90
echo "Building Pfaf_to_2d_30s_land.f90 ..."
./build Pfaf_to_2d_30s_land.f90
echo "Building read_riveroutlet_land.f90 ..."
./build read_riveroutlet_land.f90

echo "Getting the outlet locations in land:"
echo "running get_outlets_catchindex.out"
./get_outlets_catchindex.out
echo "running get_outlets_land.out"
./get_outlets_land.out
echo "running get_sinkxy_land.out"
./get_sinkxy_land.out

echo "Finalizing the outlet files for use in the mk_bcs:"
echo "running get_outlets_land_allcat.out"
./get_outlets_land_allcat.out
echo "running get_landocean_Greenland_real.out"
./get_landocean_Greenland_real.out
echo "running Pfaf_to_2d_30s_land.out"
./Pfaf_to_2d_30s_land.out
echo "running read_riveroutlet_land.out"
./read_riveroutlet_land.out
echo "Outlet_latlon.43200x21600 created!"

echo "Removing temporary input/output files ..."
rm -rf outputs
rm -rf inputs
echo "Removing *.out files ..."
rm -f *.out 
rm -f *.mod
