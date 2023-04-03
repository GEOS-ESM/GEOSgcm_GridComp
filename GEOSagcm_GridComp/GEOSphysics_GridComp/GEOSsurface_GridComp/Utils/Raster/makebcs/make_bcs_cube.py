#!/usr/bin/env python3
#
# source install/bin/g5_modules

import os
from make_bcs_questionnarie import *

cube_template = """#!/bin/csh -x

#SBATCH --output={EXPDIR}/{OUTDIR}/logs/{BCNAME}/{BCNAME2}.log
#SBATCH --error={EXPDIR}/{OUTDIR}/logs/{BCNAME}/{BCNAME2}.err
#SBATCH --account={account}
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --job-name={BCNAME2}.j
#SBATCH --constraint=sky|cas

cd {BCDIR}

if ( {STEP1} == True ) then
  /bin/ln -s {bin_dir}

  mkdir -p geometry land til rst data/MOM5 data/MOM6 clsm/plots

  ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/360x200 data/MOM5/360x200
  ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/720x410 data/MOM5/720x410
  ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/1440x1080 data/MOM5/1440x1080
  ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM6/72x36 data/MOM6/72x36
  ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM6/540x458 data/MOM6/540x458
  ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM6/1440x1080 data/MOM6/1440x1080

  if( -e CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}.stdout ) /bin/rm -f CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}.stdout

endif 

source bin/g5_modules
setenv MASKFILE {MASKFILE}
setenv MAKE_BCS_INPUT_DIR {MAKE_BCS_INPUT_DIR}
limit stacksize unlimited

if ( {STEP1} == True ) then
  bin/mkCubeFVRaster.x -x {NX} -y {NY} {NC} >/dev/null 
  bin/mkLandRaster.x -x {NX} -y {NY} -v -t {NT}
endif

if( {LATLON_OCEAN} == True ) then

   if ( {STEP1} == True ) then 
      bin/mkLatLonRaster.x -x {NX} -y {NY} -b DE -p PE -t 0 {IMO} {JMO} >/dev/null
      bin/CombineRasters.x -f 0 -t {NT} DE{IMO}xPE{JMO} Pfafstetter >/dev/null
      bin/CombineRasters.x -t {NT} CF{NC}x6C DE{IMO}xPE{JMO}-Pfafstetter
      setenv OMP_NUM_THREADS 1
      if ( {SKIPLAND} != True ) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
   endif

   if ( {STEP2} == True ) then 
      setenv OMP_NUM_THREADS {NCPUS}
      if ( {SKIPLAND} != True ) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
      chmod 755 bin/create_README.csh
      bin/create_README.csh
   endif
endif

if( {TRIPOL_OCEAN} == True ) then
   if ( {STEP1} == True ) then 
      bin/mkMOMAquaRaster.x -x {NX} -y {NY}  data/{MOM_VERSION}/{imo}x{jmo}/MAPL_Tripolar.nc > /dev/null
      /bin/cp til/Pfafstetter.til til/Pfafstetter-ORIG.til
      /bin/cp rst/Pfafstetter.rst rst/Pfafstetter-ORIG.rst
      bin/FillMomGrid.x -f 0 -g Pfafstetter-M {DATENAME}{IMO}x{POLENAME}{JMO} Pfafstetter data/{MOM_VERSION}/{imo}x{jmo}/MAPL_Tripolar.nc 
      /bin/mv til/Pfafstetter-M.til til/Pfafstetter.til
      /bin/mv rst/Pfafstetter-M.rst rst/Pfafstetter.rst
      bin/CombineRasters.x -f 0 -t {NT} {DATENAME}{IMO}x{POLENAME}{JMO} Pfafstetter >/dev/null
      bin/CombineRasters.x -t {NT} CF{NC}x6C {DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter
      bin/mk_runofftbl.x CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter
      setenv OMP_NUM_THREADS 1
      if ({SKIPLAND} != True) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter -v {lbcsv}
   endif

   if ( {STEP2} == True ) then 
      setenv OMP_NUM_THREADS {NCPUS}
      if ({SKIPLAND} != True) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter -v {lbcsv}
      chmod 755 bin/create_README.csh
      bin/create_README.csh
   endif
endif

if( {CUBED_SPHERE_OCEAN} == True ) then
   if ( {STEP1} == True ) then 
      bin/CombineRasters.x -f 0 -t {NT} CF{NC}x6C Pfafstetter >/dev/null
      bin/CombineRasters.x -t {NT} CF{NC}x6C CF{NC}x6C-Pfafstetter
      setenv OMP_NUM_THREADS 1
      if ({SKIPLAND} != True) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C_CF{NC}x6C-Pfafstetter -v {lbcsv}
   endif

   if ( {STEP2} == True ) then 
      setenv OMP_NUM_THREADS {NCPUS}
      if ({SKIPLAND} != True) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C_CF{NC}x6C-Pfafstetter -v {lbcsv} 
      chmod 755 bin/create_README.csh
      bin/create_README.csh
   endif
endif

if ( {STEP2} == True ) then 
/bin/mv clsm  clsm.C{NC}
/bin/cp til/CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.til clsm.C{NC}
if( {TRIPOL_OCEAN} == True ) /bin/cp til/CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.TRN clsm.C{NC}
/bin/rm clsm.C{NC}/CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.file

cd clsm.C{NC}
   /bin/mv irrig.dat irrigation_{RC}.dat
   /bin/mv vegdyn.data   vegdyn_{RC}.dat
   /bin/mv nirdf.dat      nirdf_{RC}.dat
   /bin/mv visdf.dat      visdf_{RC}.dat
   /bin/mv   lai.dat   lai_clim_{RC}.data
   /bin/mv green.dat green_clim_{RC}.data
   /bin/mv lnfm.dat   lnfm_clim_{RC}.data
   /bin/mv ndvi.dat   ndvi_clim_{RC}.data

/bin/rm -f sedfile
if( {CUBED_SPHERE_OCEAN} == True ) then
cat > sedfile << EOF
s/{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter/OC{nc}x{nc6}-CF/g
s/CF{NC}x6C/PE{nc}x{nc6}-CF/g
EOF
sed -f sedfile       CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.til > tile.file
/bin/mv -f tile.file CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.til
/bin/rm -f sedfile
else
cat > sedfile << EOF
s/CF{NC}x6C/PE{nc}x{nc6}-CF/g
s/{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter/PE{imo}x{jmo}-{DATENAME}/g
EOF
sed -f sedfile       CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.til > tile.file
/bin/mv -f tile.file CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.til
/bin/rm -f sedfile
endif

cd ../

/bin/rm -rf         CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}
/bin/mv clsm.C{NC} CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}
cd  CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}
mkdir clsm
/bin/mv ar.new \
        bf.dat \
        ts.dat \
        catchment.def \
        cti_stats.dat \
        tau_param.dat \
        soil_param.dat \
        mosaic_veg_typs_fracs \
        soil_param.first \
        README \
        bad_sat_param.tiles \
        lai.* \
        AlbMap* \
        plots \
        CLM_veg_typs_fracs \
        mkCatchParam.log \
        CLM_NDep_SoilAlb_T2m \
        CLM4.5_abm_peatf_gdp_hdm_fc \
        catch_params.nc4 \
        catchcn_params.nc4 \
        country_and_state_code.data \
        clsm
cd  ../ 

/bin/mv rst CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}
/bin/mv til CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}
/bin/mv {BCJOB} CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}/.

cd ../../

/bin/mv    {BCDIR}/{BCNAME} .
/bin/mv    {BCJOB}                {BCNAME}
/bin/mv    {EXPDIR}/{OUTDIR}/logs {BCNAME}/.
/bin/mv    {BCNAME}/clsm/mkCatchParam.log {BCNAME}/logs/{BCNAME}/mkCatchParam.log

mkdir -p {BCNAME}/geometry/CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}
mkdir -p {BCNAME}/land/CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}

/bin/mv  {BCNAME}/CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.til  {BCNAME}/geometry/CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}/.


/bin/mv  {BCNAME}/clsm               {BCNAME}/land/CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}/.
/bin/mv  {BCNAME}/irrigation_{RC}.dat  \
         {BCNAME}/vegdyn_{RC}.dat      \
         {BCNAME}/nirdf_{RC}.dat       \
         {BCNAME}/visdf_{RC}.dat       \
         {BCNAME}/lai_clim_{RC}.data   \
         {BCNAME}/green_clim_{RC}.data \
         {BCNAME}/lnfm_clim_{RC}.data  \
         {BCNAME}/ndvi_clim_{RC}.data  \
         {BCNAME}/land/CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}/.

/bin/mv {BCNAME}/rst {BCNAME}/til {BCNAME}/geometry/CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}/.
/bin/mv  {BCNAME}/CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}.j {BCNAME}/geometry/CF{NC}x6C_{DATENAME}{IMO}x{POLENAME}{JMO}/.
/bin/rm -r {OUTDIR}

endif  # STEP2

"""

def make_bcs_cube(config):

  bin_dir = os.getcwd()
  if 'install/bin' not in bin_dir:
    print(" please run this program in installed bin directory")
    return

  grid_type  = config['grid_type']
  if 'Cubed' not in grid_type :
     print('This is not a Cubed-Sphere grid')
     return

  resolution = config['resolution']

  account = get_account()
  IMO = '%04d'%config['imo']
  JMO = '%04d'%config['jmo']
  NC  = '%04d'%config['im']
  imo = config['imo']
  jmo = config['jmo']
  nc  = config['im']
  nc6 = nc*6
  RC = str(nc)+'x'+str(nc6)

  DATENAME = config['DATENAME']
  POLENAME = config['POLENAME']
  bcname = 'CF'+NC+'x6C_'+DATENAME+IMO+'x'+POLENAME+JMO
  SKIPLAND = config['skipland']

  if config['CUBED_SPHERE_OCEAN'] :
    bcname =  'CF'+ NC+'x6_CF'+NC+'x6C'
    DATENAME = 'CF'
    POLENAME = ''
    IMO = NC
    JMO = '6C'

  now   = datetime.now()
  tmp_dir =now.strftime("%Y%m%d%H%M%S") 
  expdir = config['expdir']
  scratch_dir = expdir+ tmp_dir+'/'+bcname+'.scratch/'
  log_dir     = expdir+'/'+tmp_dir+'/logs/'+ bcname
  bcjob       = scratch_dir+'/'+bcname+'.j'

  if os.path.exists(bcjob):
    print('please remove the run temprory directory: ' + expdir+'/'+ tmp_dir) 
    return

  os.makedirs(scratch_dir)
  if not os.path.exists(log_dir):
    os.makedirs(log_dir)

  STEP1 = True
  STEP2 = True
  BCNAME2 = bcname
  if resolution in ['c2880', 'c3072', 'c5760'] :
     STEP1 = True
     STEP2 = False
  script_string = cube_template.format(\
           account = account, \
           EXPDIR = config['expdir'], \
           OUTDIR = tmp_dir, \
           BCNAME  = bcname, \
           BCNAME2 = BCNAME2, \
           STEP1   = STEP1, \
           STEP2   = STEP2, \
           BCDIR = scratch_dir, \
           bin_dir = bin_dir, \
           MAKE_BCS_INPUT_DIR = config['inputdir'], \
           BCJOB =  bcjob, \
           DATENAME = DATENAME, \
           POLENAME = POLENAME, \
           SKIPLAND = SKIPLAND, \
           MOM_VERSION = config['MOM_VERSION'], \
           LATLON_OCEAN= config['LATLON_OCEAN'], \
           TRIPOL_OCEAN= config['TRIPOL_OCEAN'], \
           CUBED_SPHERE_OCEAN = config['CUBED_SPHERE_OCEAN'], \
           nc  = nc, \
           nc6 = nc6, \
           imo = config['imo'], \
           jmo = config['jmo'], \
           IRRIGTHRES = 2, \
           IMO = IMO, \
           JMO = JMO, \
           NC  = NC, \
           MASKFILE = config['MASKFILE'], \
           lbcsv    = config['lbcsv'], \
           NX = config['NX'], \
           NY = config['NY'], \
           NT = config['NT'], \
           RC = RC,\
           NCPUS = config['NCPUS'])

  cube_job = open(bcjob,'wt')
  cube_job.write(script_string)
  cube_job.close()

  if resolution in ['c2880', 'c3072', 'c5760'] :
     STEP1 = False
     STEP2 = True
     BCNAME2 = bcname+'-2'
     script_string = cube_template.format(\
           account = account, \
           EXPDIR = config['expdir'], \
           OUTDIR = tmp_dir, \
           BCNAME  = bcname, \
           BCNAME2 = BCNAME2, \
           STEP1   = STEP1, \
           STEP2   = STEP2, \
           BCDIR = scratch_dir, \
           bin_dir = bin_dir, \
           MAKE_BCS_INPUT_DIR = config['inputdir'], \
           BCJOB =  bcjob+'-2', \
           DATENAME = DATENAME, \
           POLENAME = POLENAME, \
           SKIPLAND = SKIPLAND, \
           MOM_VERSION = config['MOM_VERSION'], \
           LATLON_OCEAN= config['LATLON_OCEAN'], \
           TRIPOL_OCEAN= config['TRIPOL_OCEAN'], \
           CUBED_SPHERE_OCEAN = config['CUBED_SPHERE_OCEAN'], \
           nc  = nc, \
           nc6 = nc6, \
           imo = config['imo'], \
           jmo = config['jmo'], \
           IRRIGTHRES = 2, \
           IMO = IMO, \
           JMO = JMO, \
           NC  = NC, \
           MASKFILE = config['MASKFILE'], \
           lbcsv    = config['lbcsv'], \
           NX = config['NX'], \
           NY = config['NY'], \
           NT = config['NT'], \
           RC = RC,\
           NCPUS = config['NCPUS'])

     cube_job = open(bcjob+'-2','wt')
     cube_job.write(script_string)
     cube_job.close()

  interactive = os.getenv('SLURM_JOB_ID', default = None)
  if ( interactive ) :
     print('interactive mode\n')
     ntasks = os.getenv('SLURM_NTASKS', default = None)
     if ( not ntasks):
        nnodes = int(os.getenv('SLURM_NNODES', default = '1'))
        ncpus  = int(os.getenv('SLURM_CPUS_ON_NODE', default = '28'))
        subprocess.call(['chmod', '755', bcjob])
        log_name = bcjob+'.log'
        print(bcjob+  '  1>' + log_name  + '  2>&1')
        os.system(bcjob + ' 1>' + log_name+ ' 2>&1')
  else:
    print("sbatch " + bcjob +"\n")
    out = subprocess.check_output(['sbatch', bcjob])
    jobid = str(int(out.split()[3]))
    print( "Submitted batch job " + jobid)
    if resolution in ['c2880', 'c3072', 'c5760']:
      subprocess.call(['sbatch', '--dependency=afterok:'+jobid, bcjob+'-2'])
       
  print( "cd " + bin_dir)
  os.chdir(bin_dir)
 
  #print(script_string)

if __name__ == "__main__":

   answers = ask_questions()
   configs = get_configs_from_answers(answers)
   for config in configs:
      if 'Cubed-Sphere' in config['grid_type']:
         make_bcs_cube(config)

