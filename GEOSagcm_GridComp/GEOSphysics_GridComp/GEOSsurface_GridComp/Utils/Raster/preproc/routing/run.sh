#!/bin/bash
set -e

INPUT=/discover/nobackup/yzeng3/work/outlets/inputs

module load comp/intel/2021.3.0
module load ncl

rm -rf inputs >& /dev/null
rm -rf outputs >& /dev/null
rm -f *.mod >& /dev/null
rm -f *.out >& /dev/null
rm -f Outlet_latlon.43200x21600 >& /dev/null 

mkdir -p inputs outputs
ln -s ${INPUT}/* inputs 

echo "Building get_outlets_land.f90 ..."
./build get_outlets_land.f90
echo "Building get_landocean_Greenland_real_TM0072xTM0036.f90 ..."
./build get_landocean_Greenland_real_TM0072xTM0036.f90
echo "Building get_mask_MAPL_2d.f90 ..."
./build get_mask_MAPL_2d.f90
echo "Building get_mask_TM0072xTM0036.f90 ..."
./build get_mask_TM0072xTM0036.f90
echo "Building get_oceanbond_TM0072xTM0036_mask.f90 ..."
./build get_oceanbond_TM0072xTM0036_mask.f90
echo "Building get_oceanbond_points_TM0072xTM0036_mask.f90 ..."
./build get_oceanbond_points_TM0072xTM0036_mask.f90
echo "Building mv_outlets_ocean.f90 ..."
./build mv_outlets_ocean.f90
echo "Building get_outlets_ocean_allcat.f90 ..."
./build get_outlets_ocean_allcat.f90
echo "Building Pfaf_to_2d_30s.f90 ..."
./build Pfaf_to_2d_30s.f90
echo "Building read_riveroutlet.f90 ..."
./build read_riveroutlet.f90

echo "STEP ONE:"
echo "running get_outlets_catchindex.ncl"
ncl get_outlets_catchindex.ncl
echo "running get_outlets_land.out"
./get_outlets_land.out
echo "running get_sinkxy_land.ncl"
ncl get_sinkxy_land.ncl

echo "STEP TWO:"
echo "running get_landocean_Greenland_real_TM0072xTM0036.out"
./get_landocean_Greenland_real_TM0072xTM0036.out
echo "running get_mask_MAPL_1d.ncl"
ncl get_mask_MAPL_1d.ncl
echo "running get_mask_MAPL_2d.out"
./get_mask_MAPL_2d.out
echo "running get_mask_TM0072xTM0036.out"
./get_mask_TM0072xTM0036.out
echo "running get_oceanbond_TM0072xTM0036_mask.out"
./get_oceanbond_TM0072xTM0036_mask.out
echo "running get_oceanbond_points_TM0072xTM0036_mask.out"
./get_oceanbond_points_TM0072xTM0036_mask.out
echo "running mv_outlets_ocean.out"
./mv_outlets_ocean.out
echo "running get_sinkxy_ocean.ncl"
ncl get_sinkxy_ocean.ncl

echo STEP THREE:""
echo "running get_outlets_ocean_allcat.out"
./get_outlets_ocean_allcat.out
echo "running Pfaf_to_2d_30s.out"
./Pfaf_to_2d_30s.out
echo "running read_riveroutlet.out"
./read_riveroutlet.out
echo "Outlet_latlon.43200x21600 created!"

echo "Removing temporary input/output files ..."
rm -rf outputs
rm -rf inputs
echo "Removing *.out files ..."
rm -f *.out 
rm -f *.mod
