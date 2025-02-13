#!/usr/bin/env python3
#
# source install/bin/g5_modules

import os
import glob

BUILT_ON_SLES15 = "@BUILT_ON_SLES15@"
def get_script_head() :

  head =  """#!/bin/csh -x

#SBATCH --output={EXPDIR}/{TMP_DIR}/logs/{GRIDNAME}/{GRIDNAME}.log
#SBATCH --error={EXPDIR}/{TMP_DIR}/logs/{GRIDNAME}/{GRIDNAME}.err
#SBATCH --account={account}
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --job-name={GRIDNAME}.j
"""
  constraint = '#SBATCH --constraint="[mil|cas]"'
  #if 'TRUE' not in BUILT_ON_SLES15:
  #   constraint = "#SBATCH --constraint=sky"

  head = head + constraint + """
echo "-----------------------------" 
echo "make_bcs starts date/time" 
echo `date` 
echo "-----------------------------" 

cd {SCRATCH_DIR}

if( ! -d bin ) then
  /bin/ln -s {bin_dir}
endif

source bin/g5_modules
module load nco
setenv MASKFILE {MASKFILE}
setenv MAKE_BCS_INPUT_DIR {MAKE_BCS_INPUT_DIR}
limit stacksize unlimited

if( ! -d geometry ) then
  mkdir -p geometry land/shared til rst data/MOM5 data/MOM6 clsm/plots
endif
"""
  return head

def get_change_til_file(grid_type):
  script = ""

  if grid_type == "Stretched_CS" or grid_type == "Cubed-Sphere" :

       script = """

cd geometry/{GRIDNAME}/
if( {TRIPOL_OCEAN} == True ) then
  sed -i 's/{OCEAN_VERSION}{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter/PE{imo}x{jmo}-{OCEAN_VERSION}/g' {GRIDNAME}{RS}.til
  sed -i 's/CF{NC}x6C/PE{nc}x{nc6}-CF/g' {GRIDNAME}{RS}.til
  ncatted -a Grid_Name,global,o,c,'PE{nc}x{nc6}-CF' {GRIDNAME}{RS}.nc4
  ncatted -a Grid_ocn_Name,global,o,c,'PE{imo}x{jmo}-{OCEAN_VERSION}' {GRIDNAME}{RS}.nc4
endif
if( {CUBED_SPHERE_OCEAN} == True ) then
  sed -i 's/{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter/OC{nc}x{nc6}-CF/g' {GRIDNAME}{RS}.til
  sed -i 's/CF{NC}x6C{SGNAME}/PE{nc}x{nc6}-CF/g' {GRIDNAME}{RS}.til
  ncatted -a Grid_Name,global,o,c,'PE{nc}x{nc6}-CF' {GRIDNAME}{RS}.nc4
  ncatted -a Grid_ocn_Name,global,o,c,'OC{nc}x{nc6}-CF' {GRIDNAME}{RS}.nc4
endif
if( {LATLON_OCEAN} == True ) then
  sed -i 's/{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter/PE{imo}x{jmo}-{DATENAME}/g' {GRIDNAME}{RS}.til
  sed -i 's/CF{NC}x6C{SGNAME}/PE{nc}x{nc6}-CF/g' {GRIDNAME}{RS}.til
  ncatted -a Grid_Name,global,o,c,'PE{nc}x{nc6}-CF' {GRIDNAME}{RS}.nc4
  ncatted -a Grid_ocn_Name,global,o,c,'PE{imo}x{jmo}-{DATENAME}' {GRIDNAME}{RS}.nc4
endif
cd ../../

"""
  if grid_type == "Lat-Lon" :

     script = """
cd geometry/{GRIDNAME}/
sed -i 's/{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter/PE{imo}x{jmo}-{DATENAME}/g' {GRIDNAME}{RS}.til
sed -i 's/DC{IM}xPC{JM}/PC{im}x{jm}-DC/g' {GRIDNAME}{RS}.til
ncatted -a Grid_Name,global,o,c,'PE{nc}x{nc6}-CF' {GRIDNAME}{RS}.nc4
ncatted -a Grid_ocn_Name,global,o,c,'PE{imo}x{jmo}-{DATENAME}' {GRIDNAME}{RS}.nc4
cd ../../

"""
  return script

def get_script_mv(grid_type):

   mv_template = """

mkdir -p geometry/{GRIDNAME}
/bin/mv {GRIDNAME}.j geometry/{GRIDNAME}/.
/bin/cp til/{GRIDNAME}{RS}.til geometry/{GRIDNAME}/.
/bin/cp til/{GRIDNAME}{RS}.nc4 geometry/{GRIDNAME}/.
if( {TRIPOL_OCEAN} == True ) /bin/cp til/{GRIDNAME}{RS}.TRN geometry/{GRIDNAME}/.

/bin/mv rst til geometry/{GRIDNAME}/.

mkdir -p land/{GRIDNAME}/clsm

/bin/mv clsm/irrig.dat   land/{GRIDNAME}/irrigation_{RC}.dat
/bin/mv clsm/vegdyn.data land/{GRIDNAME}/vegdyn_{RC}.dat
/bin/mv clsm/nirdf.dat   land/{GRIDNAME}/nirdf_{RC}.dat
/bin/mv clsm/visdf.dat   land/{GRIDNAME}/visdf_{RC}.dat
/bin/mv clsm/lai.dat     land/{GRIDNAME}/lai_clim_{RC}.data
/bin/mv clsm/green.dat   land/{GRIDNAME}/green_clim_{RC}.data
/bin/mv clsm/lnfm.dat    land/{GRIDNAME}/lnfm_clim_{RC}.data
/bin/mv clsm/ndvi.dat    land/{GRIDNAME}/ndvi_clim_{RC}.data

# vegdyn_{RC}.dat file is nc4; for clarification, create link with proper file name extension
cd land/{GRIDNAME}
ln -s vegdyn_{RC}.dat vegdyn_{RC}.nc4
cd ../../

/bin/mv clsm/ar.new \\
        clsm/bf.dat \\
        clsm/ts.dat \\
        clsm/catchment.def \\
        clsm/cti_stats.dat \\
        clsm/tau_param.dat \\
        clsm/soil_param.dat \\
        clsm/mosaic_veg_typs_fracs \\
        clsm/soil_param.first \\
        clsm/bad_sat_param.tiles \\
        clsm/README \\
        clsm/lai.* \\
        clsm/AlbMap* \\
        clsm/g5fmt \\
        clsm/vegetation.hst2 \\
        clsm/pfaf_fractions.dat \\
        clsm/plots \\
        clsm/CLM_veg_typs_fracs \\
        clsm/Grid2Catch_TransferData.nc \\
        clsm/CLM_NDep_SoilAlb_T2m \\
        clsm/CLM4.5_abm_peatf_gdp_hdm_fc \\
        clsm/catch_params.nc4 \\
        clsm/catchcn_params.nc4 \\
        clsm/country_and_state_code.data \\
        land/{GRIDNAME}/clsm/

""" 
   mv_template = mv_template + get_change_til_file(grid_type)
   mv_template = mv_template + """

/bin/mv clsm/mkCatchParam.log ../logs/{GRIDNAME}/mkCatchParam.log

# move output into final directory tree (layout as of July 2023)

mkdir -p ../../geometry ../../land/shared ../../logs

echo "-----------------------------" 
echo "make_bcs ends date/time" 
echo `date` 
echo "-----------------------------" 

/bin/mv ../logs/{GRIDNAME}  ../../logs/.

/bin/mv geometry/{GRIDNAME} ../../geometry/.
/bin/mv land/{GRIDNAME}     ../../land/.

# cd out of and clean up the temporary directory

cd ../..

/bin/rm -r {TMP_DIR}

# if necessary, copy resolution-independent CO2 file from MAKE_BCS_INPUT_DIR to bcs dir 

if(-f land/shared/CO2_MonthlyMean_DiurnalCycle.nc4) then
    echo "CO2_MonthlyMean_DiurnalCycle.nc4 already present in bcs dir."
else
    /bin/cp -p {MAKE_BCS_INPUT_DIR}/land/CO2/v1/CO2_MonthlyMean_DiurnalCycle.nc4 land/shared/CO2_MonthlyMean_DiurnalCycle.nc4
    echo "Successfully copied CO2_MonthlyMean_DiurnalCycle.nc4 to bcs dir."
endif

# adjust permissions

chmod +rX -R geometry land logs

"""

   return mv_template

def check_script(expdir, script):

  if glob.glob(os.path.join(expdir, '**', script), recursive=True):
       print("-------------------------------------------------------------------")
       print("                 Abort:                                            ")
       print(" This BCS run is configured to create output that already exists:  ")
       print(expdir, ' has this grid ', script)
       print(" Please delete run dir and same resolution BCS files and resubmit. ")
       print("-------------------------------------------------------------------")
       exit()

