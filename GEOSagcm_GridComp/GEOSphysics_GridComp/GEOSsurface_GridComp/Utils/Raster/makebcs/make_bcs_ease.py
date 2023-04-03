#!/usr/bin/env python3
#
# source install/bin/g5_modules

import os
from make_bcs_questionnarie import *


ease_template = """#!/bin/csh -x

#SBATCH --output={EXPDIR}/{OUTDIR}/logs/{GRIDNAME}/{GRIDNAME}.log
#SBATCH --error={EXPDIR}/{OUTDIR}/logs/{GRIDNAME}/{GRIDNAME}.err
#SBATCH --account={account}
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --job-name={GRIDNAME}.j
#SBATCH --constraint=sky|cas

cd {BCDIR}

/bin/ln -s {bin_dir}
source bin/g5_modules
mkdir -p geometry land til rst data/MOM5 data/MOM6 clsm/plots
limit stacksize unlimited

setenv MASKFILE {MASKFILE}
setenv MAKE_BCS_INPUT_DIR {MAKE_BCS_INPUT_DIR}
setenv OMP_NUM_THREADS 1
bin/mkEASETilesParam.x -ease_label {GRIDNAME} 
setenv OMP_NUM_THREADS 1
bin/mkCatchParam.x -g {GRIDNAME} -v {lbcsv} -x {NX} -y {NY}
setenv OMP_NUM_THREADS {NCPUS}
bin/mkCatchParam.x -g {GRIDNAME} -v {lbcsv} -x {NX} -y {NY}
chmod 755 bin/create_README.csh
bin/create_README.csh

/bin/mv clsm  clsm.{IM}x{JM}
/bin/cp til/{EASEVERSION}_{RES}_{RS}.til clsm.{IM}x{JM}

cd clsm.{IM}x{JM}
   /bin/mv irrig.dat irrigation_{RS}_DE.dat
   /bin/mv vegdyn.data   vegdyn_{RS}_DE.dat
   /bin/mv nirdf.dat      nirdf_{RS}_DE.dat
   /bin/mv visdf.dat      visdf_{RS}_DE.dat
   /bin/mv   lai.dat   lai_clim_{RS}_DE.data
   /bin/mv green.dat green_clim_{RS}_DE.data
   /bin/mv lnfm.dat   lnfm_clim_{RS}_DE.data
   /bin/mv ndvi.dat   ndvi_clim_{RS}_DE.data


cd ../
/bin/rm -rf              {EASEVERSION}_{RES}
/bin/mv clsm.{IM}x{JM} {EASEVERSION}_{RES}
cd  {EASEVERSION}_{RES} 
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
        bad_sat_param.tiles \
        README \
        lai.* \
        AlbMap* \
        g5fmt \
        vegetation.hst2 \
        pfaf_fractions.dat \
        plots \
        CLM_veg_typs_fracs \
        mkCatchParam.log \
        Grid2Catch_TransferData.nc \
        CLM_NDep_SoilAlb_T2m \
        CLM4.5_abm_peatf_gdp_hdm_fc \
        catch_params.nc4 \
        catchcn_params.nc4 \
        country_and_state_code.data \
        clsm
cd  ../ 
/bin/mv rst  {EASEVERSION}_{RES}
/bin/mv til  {EASEVERSION}_{RES}

cd ../../
/bin/mv    {BCDIR}/{GRIDNAME} .
/bin/mv    {BCJOB}                {GRIDNAME}
/bin/mv    {EXPDIR}/{OUTDIR}/logs   {GRIDNAME}/.
/bin/mv    {GRIDNAME}/clsm/mkCatchParam.log {GRIDNAME}/logs/{GRIDNAME}/mkCatchParam.log

mkdir -p {GRIDNAME}/geometry/{EASEVERSION}_{RES}
mkdir -p {GRIDNAME}/land/{EASEVERSION}_{RES}

/bin/mv  {GRIDNAME}/{EASEVERSION}_{RES}_{RS}.til  {GRIDNAME}/geometry/{EASEVERSION}_{RES}/.


/bin/mv  {GRIDNAME}/clsm               {GRIDNAME}/land/{EASEVERSION}_{RES}/.
/bin/mv  {GRIDNAME}/irrigation_{RS}_DE.dat  \
         {GRIDNAME}/vegdyn_{RS}_DE.dat      \
         {GRIDNAME}/nirdf_{RS}_DE.dat       \
         {GRIDNAME}/visdf_{RS}_DE.dat       \
         {GRIDNAME}/lai_clim_{RS}_DE.data   \
         {GRIDNAME}/green_clim_{RS}_DE.data \
         {GRIDNAME}/lnfm_clim_{RS}_DE.data  \
         {GRIDNAME}/ndvi_clim_{RS}_DE.data  \
         {GRIDNAME}/land/{EASEVERSION}_{RES}/.

/bin/mv {GRIDNAME}/rst {GRIDNAME}/til {GRIDNAME}/geometry/{EASEVERSION}_{RES}/.
/bin/mv {GRIDNAME}/{GRIDNAME}.j {GRIDNAME}/geometry/{EASEVERSION}_{RES}/.
/bin/rm -r {OUTDIR}

"""

def make_bcs_ease(config):
  bin_dir = os.getcwd()
  if 'install/bin' not in bin_dir:
    print(" please run this program in installed bin directory")
    return

  grid_type  = config['grid_type']
  if 'EASEv' not in grid_type :
     print('This is not a EASE grid')
     return

  resolution = config['resolution']

  gridname  = grid_type+'_'+ resolution
  now   = datetime.now()
  tmp_dir =now.strftime("%Y%m%d%H%M%S") 
  tmp_dir=f"{resolution}_{tmp_dir}"
  expdir = config['expdir']
  scratch_dir = expdir+ tmp_dir+'/'+gridname+'.scratch/'
  log_dir     = expdir+'/'+tmp_dir+'/logs/' + gridname
  job_script       = scratch_dir+'/'+gridname+'.j'
  if os.path.exists(job_script):
    print('please remove the run temprory directory: ' + expdir+'/'+ tmp_dir) 
    return

  os.makedirs(scratch_dir)
  if not os.path.exists(log_dir):
    os.makedirs(log_dir)

  account = get_account()
  ims = '%04d'%config['im']
  jms = '%04d'%config['jm']
  RS = str(config['im'])+'x'+ str(config['jm'])

  script_string = ease_template.format(\
           account = account, \
           EXPDIR = config['expdir'], \
           OUTDIR = tmp_dir, \
           GRIDNAME = gridname, \
           bin_dir = bin_dir, \
           MAKE_BCS_INPUT_DIR = config['inputdir'], \
           BCJOB =  job_script, \
           EASEVERSION = grid_type, \
           RES = resolution, \
           IM = ims, \
           JM = jms, \
           MASKFILE = config['MASKFILE'], \
           lbcsv    = config['lbcsv'], \
           NX = config['NX'], \
           NY = config['NY'], \
           RS = RS,\
           BCDIR = scratch_dir, \
           NCPUS = config['NCPUS'])


  ease_job = open(job_script,'wt')
  ease_job.write(script_string)
  ease_job.close()

  interactive = os.getenv('SLURM_JOB_ID', default = None)
  if ( interactive ) :
     print('interactive mode\n')
     ntasks = os.getenv('SLURM_NTASKS', default = None)
     if ( not ntasks):
        nnodes = int(os.getenv('SLURM_NNODES', default = '1'))
        ncpus  = int(os.getenv('SLURM_CPUS_ON_NODE', default = '28'))
        subprocess.call(['chmod', '755', job_script])
        log_name = job_script +'.log'
        print(job_script+  '  1>' + log_name  + '  2>&1')
        os.system(job_script + ' 1>' + log_name+ ' 2>&1')
  else:
    print("sbatch " + job_script +"\n")
    subprocess.call(['sbatch', job_script])

  print( "cd " + bin_dir)
  os.chdir(bin_dir)
 
  #print(script_string)

if __name__ == "__main__":

   answers = ask_questions(default_grid="EASEv2")
   configs = get_configs_from_answers(answers)
   for config in configs:
      if 'EASEv' in config['grid_type']:
         make_bcs_ease(config)

