#!/usr/bin/env python3
#
# source install/bin/g5_modules

import os
from questionnarie_bcs import *


ease_template = """#!/bin/csh -x

#SBATCH --output={EXPDIR}/{OUTDIR}/logs/{BCNAME}.log
#SBATCH --error={EXPDIR}/{OUTDIR}/logs/{BCNAME}.err
#SBATCH --account={account}
#SBATCH --time=12:00:00
#SBATCH --ntasks=28
#SBATCH --job-name={BCNAME}.j
#SBATCH --constraint=sky

cd {BCDIR}

/bin/ln -s {bin_dir}
source bin/g5_modules
mkdir -p til rst data/MOM5 data/MOM6 clsm/plots
cd data 
cd ../
limit stacksize unlimited

setenv MASKFILE {MASKFILE}
setenv MAKE_BCS_INPUT_DIR {MAKE_BCS_INPUT_DIR}
setenv OMP_NUM_THREADS 1
bin/mkEASETilesParam.x -ease_label {BCNAME} 
setenv OMP_NUM_THREADS 1
bin/mkCatchParam.x -g {BCNAME} -v {lbcsv} -x {NX} -y {NY}
setenv OMP_NUM_THREADS {NCPUS}
bin/mkCatchParam.x -g {BCNAME} -v {lbcsv} -x {NX} -y {NY}
chmod 755 bin/create_README.csh
bin/create_README.csh

/bin/mv clsm  clsm.{IM}x{JM}
/bin/cp til/SMAP_{EASEVERSION}_{HRCODE}_{RS}.til clsm.{IM}x{JM}

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
/bin/rm -rf              SMAP_{EASEVERSION}_{HRCODE}
/bin/mv clsm.{IM}x{JM} SMAP_{EASEVERSION}_{HRCODE}
                     cd  SMAP_{EASEVERSION}_{HRCODE} 
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
/bin/mv rst  SMAP_{EASEVERSION}_{HRCODE}
/bin/mv til  SMAP_{EASEVERSION}_{HRCODE}

cd ../../
/bin/mv    {BCDIR}/{BCNAME} .
/bin/mv    {BCJOB}                {BCNAME}
/bin/mv    {EXPDIR}/{OUTDIR}/logs   {BCNAME}/.
/bin/mv    {BCNAME}/clsm/mkCatchParam.log {BCNAME}/logs/mkCatchParam.log
/bin/rm -r {OUTDIR}

"""

def make_ease_bcs(config):
  bin_dir = os.getcwd()
  if 'install/bin' not in bin_dir:
    print(" please run this program in installed bin directory")
    return

  grid_type  = config['grid_type']
  if 'EASEv' not in grid_type :
     print('This is not a EASE grid')
     return

  resolution = config['resolution']

  EASElabel  = 'SMAP_'+grid_type+'_'+ resolution
  now   = datetime.now()
  tmp_dir =now.strftime("%Y%m%d%H%M%S") 
  expdir = config['expdir']
  scratch_dir = expdir+ tmp_dir+'/'+EASElabel+'.scratch/'
  log_dir     = expdir+'/'+tmp_dir+'/logs/'
  bcjob       = scratch_dir+'/'+EASElabel+'.j'
  if os.path.exists(bcjob):
    print('please remove the run temprory directory: ' + expdir+'/'+ tmp_dir) 
    return

  os.makedirs(scratch_dir)
  os.makedirs(log_dir)

  account = get_account()
  ims = '%04d'%config['im']
  jms = '%04d'%config['jm']
  RS = ims+'x'+jms
  job_script = ease_template.format(\
           account = account, \
           EXPDIR = config['expdir'], \
           OUTDIR = tmp_dir, \
           BCNAME = EASElabel, \
           bin_dir = bin_dir, \
           MAKE_BCS_INPUT_DIR = config['inputdir'], \
           BCJOB =  bcjob, \
           EASEVERSION = grid_type, \
           HRCODE = resolution, \
           IM = config['im'], \
           JM = config['jm'], \
           MASKFILE = config['MASKFILE'], \
           lbcsv    = config['lbcsv'], \
           NX = config['NX'], \
           NY = config['NY'], \
           RS = RS,\
           BCDIR = scratch_dir, \
           NCPUS = config['NCPUS'])


  ease_job = open(bcjob,'wt')
  ease_job.write(job_script)
  ease_job.close()

  interactive = os.getenv('SLURM_JOB_ID', default = None)
  if ( interactive ) :
     print('interactive mode\n')
     ntasks = os.getenv('SLURM_NTASKS', default = None)
     if ( not ntasks):
        nnodes = int(os.getenv('SLURM_NNODES', default = '1'))
        ncpus  = int(os.getenv('SLURM_CPUS_ON_NODE', default = '28'))
        subprocess.call(['chmod', '755', bcjob])
        print(bcjob+  '  1>' + log_name  + '  2>&1')
        os.system(bcjob + ' 1>' + log_name+ ' 2>&1')
  else:
    print("sbatch " + bcjob +"\n")
    subprocess.call(['sbatch', bcjob])

  print( "cd " + bin_dir)
  os.chdir(bin_dir)
 
  print(job_script)

if __name__ == "__main__":

   answers = ask_questions()
   config = get_config_from_answers(answers)
   print("make_ease_bcs")
   make_ease_bcs(config)

