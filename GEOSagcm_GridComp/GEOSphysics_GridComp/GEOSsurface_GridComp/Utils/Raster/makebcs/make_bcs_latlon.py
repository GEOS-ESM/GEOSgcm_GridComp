#!/usr/bin/env python3
#
# source install/bin/g5_modules

import os
from make_bcs_questionnarie import *


latlon_template = """#!/bin/csh -x

#SBATCH --output={EXPDIR}/{OUTDIR}/logs/{BCNAME}/{BCNAME}.log
#SBATCH --error={EXPDIR}/{OUTDIR}/logs/{BCNAME}/{BCNAME}.err
#SBATCH --account={account}
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --job-name={BCNAME}.j
#SBATCH --constraint=sky|cas

cd {BCDIR}
/bin/ln -s {bin_dir}
source bin/g5_modules
mkdir -p  geometry land til rst data/MOM5 data/MOM6 clsm/plots
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/360x200 data/MOM5/360x200
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/720x410 data/MOM5/720x410
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/1440x1080 data/MOM5/1440x1080
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM6/72x36 data/MOM6/72x36
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM6/540x458 data/MOM6/540x458
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM6/1440x1080 data/MOM6/1440x1080

cd data 
cd ../

if( -e DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}.stdout ) /bin/rm -f DC{IM}xPC{JM}_{DATENAME}{IMO}{POLENAME}{JMO}.stdout
setenv MASKFILE {MASKFILE}
setenv MAKE_BCS_INPUT_DIR {MAKE_BCS_INPUT_DIR}
limit stacksize unlimited
bin/mkLatLonRaster.x -x {NX} -y {NY}  -t -1 {IM} {JM} >/dev/null
bin/mkLandRaster.x -x {NX} -y {NY} -v -t {NT}

if( {LATLON_OCEAN} == True ) then
    bin/mkLatLonRaster.x -x {NX} -y {NY} -b DE -p PE -t 0 {IMO} {JMO} >/dev/null
    bin/CombineRasters.x -f 0 -t {NT} DE{IMO}xPE{JMO} Pfafstetter >/dev/null
    bin/CombineRasters.x -t {NT} DC{IM}xPC{JM} DE{IMO}xPE{JMO}-Pfafstetter
    setenv OMP_NUM_THREADS 1
    if ( {SKIPLAND} != True) bin/mkCatchParam.x -x {NX} -y {NY} -g DC{IM}xPC{JM}_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
    setenv OMP_NUM_THREADS {NCPUS}
    if ( {SKIPLAND} != True) bin/mkCatchParam.x -x {NX} -y {NY} -g DC{IM}xPC{JM}_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
    chmod 755 bin/create_README.csh
    bin/create_README.csh
endif
if( {TRIPOL_OCEAN} == True ) then
   bin/mkMOMAquaRaster.x -x {NX} -y {NY}  data/{MOM_VERSION}/{imo}x{jmo}/MAPL_Tripolar.nc > /dev/null
    /bin/cp til/Pfafstetter.til til/Pfafstetter-ORIG.til
    /bin/cp rst/Pfafstetter.rst rst/Pfafstetter-ORIG.rst
    bin/FillMomGrid.x -f 0 -g Pfafstetter-M {DATENAME}{IMO}x{POLENAME}{JMO} Pfafstetter data/{MOM_VERSION}/{imo}x{jmo}/MAPL_Tripolar.nc 
    /bin/mv til/Pfafstetter-M.til til/Pfafstetter.til
    /bin/mv rst/Pfafstetter-M.rst rst/Pfafstetter.rst
    bin/CombineRasters.x -f 0 -t {NT} {DATENAME}{IMO}x{POLENAME}{JMO} Pfafstetter >/dev/null
    bin/CombineRasters.x -t {NT} DC{IM}xPC{JM} {DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter
    bin/mk_runofftbl.x DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter
    setenv OMP_NUM_THREADS 1
    if ( {SKIPLAND} != True ) bin/mkCatchParam.x -x {NX} -y {NY} -g DE{IMO}xPE{JMO}_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
    setenv OMP_NUM_THREADS {NCPUS}
    if ( {SKIPLAND} != True ) bin/mkCatchParam.x -x {NX} -y {NY} -g DE{IMO}xPE{JMO}_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
    chmod 755 bin/create_README.csh
    bin/create_README.csh    
endif

/bin/mv clsm  clsm.{IM}x{JM}
/bin/cp til/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.til clsm.{IM}x{JM}
if( {TRIPOL_OCEAN} == True ) /bin/cp til/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.TRN clsm.{IM}x{JM}
/bin/rm clsm.{IM}x{JM}/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.file

cd clsm.{IM}x{JM}
   /bin/mv irrig.dat irrigation_{RS}_DC.dat
   /bin/mv vegdyn.data   vegdyn_{RS}_DC.dat
   /bin/mv nirdf.dat      nirdf_{RS}_DC.dat
   /bin/mv visdf.dat      visdf_{RS}_DC.dat
   /bin/mv   lai.dat   lai_clim_{RS}_DC.data
   /bin/mv green.dat green_clim_{RS}_DC.data
   /bin/mv lnfm.dat   lnfm_clim_{RS}_DC.data
   /bin/mv ndvi.dat   ndvi_clim_{RS}_DC.data
/bin/rm -f sedfile
cat > sedfile << EOF
s/DC{IM}xPC{JM}/PC{im}x{jm}-DC/g
s/{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter/PE{imo}x{jmo}-{DATENAME}/g
EOF
sed -f sedfile       DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.til > tile.file
/bin/mv -f tile.file DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.til
/bin/rm -f sedfile

cd ../

/bin/rm -rf              DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}
/bin/mv clsm.{IM}x{JM} DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}
cd  DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}
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

/bin/mv rst DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}
/bin/mv til DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}

cd ../../
/bin/mv    {BCDIR}/{BCNAME} .
/bin/mv    {BCJOB}                {BCNAME}
/bin/mv    {EXPDIR}/{OUTDIR}/logs  {BCNAME}/.
/bin/mv    {BCNAME}/clsm/mkCatchParam.log {BCNAME}/logs/{BCNAME}/mkCatchParam.log

mkdir -p {BCNAME}/geometry/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}
mkdir -p {BCNAME}/land/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}

/bin/mv  {BCNAME}/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.til  {BCNAME}/geometry/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}/.


/bin/mv  {BCNAME}/clsm               {BCNAME}/land/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}/.
/bin/mv  {BCNAME}/irrigation_{RS}_DC.dat \
         {BCNAME}/vegdyn_{RS}_DC.dat \
         {BCNAME}/nirdf_{RS}_DC.dat \
         {BCNAME}/visdf_{RS}_DC.dat \
         {BCNAME}/lai_clim_{RS}_DC.data \
         {BCNAME}/green_clim_{RS}_DC.data \
         {BCNAME}/lnfm_clim_{RS}_DC.data \
         {BCNAME}/ndvi_clim_{RS}_DC.data \
         {BCNAME}/land/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}/.

/bin/mv {BCNAME}/rst {BCNAME}/til {BCNAME}/geometry/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}/.
/bin/mv   {BCNAME}/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}.j {BCNAME}/geometry/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}/.
/bin/rm -r {OUTDIR}

"""

def make_bcs_latlon(config):
  bin_dir = os.getcwd()
  if 'install/bin' not in bin_dir:
    print(" please run this program in installed bin directory")
    return

  resolution = config['resolution']

  account = get_account()
  IMO = '%04d'%config['imo']
  JMO = '%04d'%config['jmo']
  IM  = '%04d'%config['im']
  JM  = '%04d'%config['jm']

  RS = str(config['im']) +'_' +str(config['jm'])

  DATENAME = config['DATENAME']
  POLENAME = config['POLENAME']

  bcname = 'DC'+IM+'xPC'+JM+'_'+ DATENAME+IMO+'x'+POLENAME+JMO 

  SKIPLAND = config['skipland'] 

  now   = datetime.now()
  tmp_dir =now.strftime("%Y%m%d%H%M%S") 
  expdir = config['expdir']
  scratch_dir = expdir+ tmp_dir+'/'+bcname+'.scratch/'
  log_dir     = expdir+'/'+tmp_dir+'/logs/' + bcname
  bcjob       = scratch_dir+'/'+bcname+'.j'

  if os.path.exists(bcjob):
    print('please remove the run temprory directory: ' +  scratch_dir) 
    return

  os.makedirs(scratch_dir)
  if not os.path.exists(log_dir):
    os.makedirs(log_dir)

  job_script = latlon_template.format(\
           account = account, \
           EXPDIR = config['expdir'], \
           OUTDIR = tmp_dir, \
           BCNAME = bcname, \
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
           im = config['im'], \
           jm = config['jm'], \
           imo = config['imo'], \
           jmo = config['jmo'], \
           IRRIGTHRES = 2, \
           IM = IM, \
           JM = JM, \
           IMO = IMO, \
           JMO = JMO, \
           MASKFILE = config['MASKFILE'], \
           lbcsv    = config['lbcsv'], \
           NX = config['NX'], \
           NY = config['NY'], \
           NT = config['NT'], \
           RS = RS,\
           NCPUS = config['NCPUS'])


  latlon_job = open(bcjob,'wt')
  latlon_job.write(job_script)
  latlon_job.close()

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
    subprocess.call(['sbatch', bcjob])

  print( "cd " + bin_dir)
  os.chdir(bin_dir)
 
  #print(job_script)

if __name__ == "__main__":

   answers = ask_questions(default_grid="Lat-Lon")
   configs = get_configs_from_answers(answers)
   for config in configs:
      if 'Lat-Lon' in config['grid_type']:
         make_bcs_latlon(config)

