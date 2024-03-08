#!/usr/bin/env python3
#
# source install/bin/g5_modules

import os
from make_bcs_questionary import *
from make_bcs_shared  import *

ease_template = """

setenv OMP_NUM_THREADS 1
bin/mkEASETilesParam.x -ease_label {GRIDNAME} 
setenv OMP_NUM_THREADS 1
bin/mkCatchParam.x -g {GRIDNAME} -v {lbcsv} -x {NX} -y {NY}
setenv OMP_NUM_THREADS {NCPUS}
bin/mkCatchParam.x -g {GRIDNAME} -v {lbcsv} -x {NX} -y {NY}
chmod 755 bin/create_README.csh
bin/create_README.csh

"""

def make_bcs_ease(config):
  bin_dir = os.getcwd()
  if '/bin' not in bin_dir:
     print("                                                                  ")
     print("******************************************************************")
     print("ERROR: Must run this program in installed bin directory.          ")
     print("******************************************************************")
     return

  grid_type  = config['grid_type']
  if 'EASEv' not in grid_type :
     print("                                                                  ")
     print("******************************************************************")
     print("ERROR: " + GRIDNAME + " is not an EASE grid.                      ")
     print("******************************************************************")
     return

  resolution = config['resolution']

  # if resolution is M01 in EASEv1 or EASEv2, do not execute
  # Instead, exist and print an on-screen message to the user
  if resolution == "M01":
     print("                                                                  ")
     print("******************************************************************")
     print("WARNING: Job for " + GRIDNAME + " not submitted!                  ")
     print("         EASEv[x]_M01 (~1 km) resolution requires custom process. ")
     print("         Contact GMAO Land Group for assistance.                  ")     
     print("******************************************************************")
     print("                                                                  ")
     return 

  GRIDNAME  = grid_type+'_'+ resolution
  now     = datetime.now()
  tmp_dir = now.strftime("%Y%m%d%H%M%S") 
  tmp_dir = f"{resolution}_{tmp_dir}"
  expdir  = config['expdir']
  scratch_dir = expdir+ tmp_dir+'/'+GRIDNAME+'.scratch/'
  log_dir     = expdir+'/'+tmp_dir+'/logs/' + GRIDNAME
  bcjob       = scratch_dir+'/'+GRIDNAME+'.j'

  check_script(expdir, GRIDNAME+'.j')

  os.makedirs(scratch_dir)
  if not os.path.exists(log_dir):
    os.makedirs(log_dir)

  account = get_account()
  ims = '%04d'%config['im']
  jms = '%04d'%config['jm']
  RS = str(config['im'])+'x'+ str(config['jm'])

  script_template = get_script_head() + ease_template + get_script_mv(config['grid_type'])
 
  script_string = script_template.format(\
           account = account, \
           EXPDIR = config['expdir'], \
           TMP_DIR = tmp_dir, \
           GRIDNAME  = GRIDNAME, \
           GRIDNAME2 = GRIDNAME, \
           bin_dir = bin_dir, \
           MAKE_BCS_INPUT_DIR = config['inputdir'], \
           IM = ims, \
           JM = jms, \
           MASKFILE = config['MASKFILE'], \
           lbcsv    = config['lbcsv'], \
           TRIPOL_OCEAN = False, \
           NX = config['NX'], \
           NY = config['NY'], \
           RS = '_'+RS,\
           RC = RS+'_DE',\
           SCRATCH_DIR = scratch_dir, \
           NCPUS = config['NCPUS'])

  ease_job = open(bcjob,'wt')
  ease_job.write(script_string)
  ease_job.close()

  interactive = os.getenv('SLURM_JOB_ID', default = None)
  if ( interactive ) :
     print('interactive mode\n')
     ntasks = os.getenv('SLURM_NTASKS', default = None)
     if ( not ntasks):
        nnodes = int(os.getenv('SLURM_NNODES', default = '1'))
        ncpus  = int(os.getenv('SLURM_CPUS_ON_NODE', default = '28'))
        subprocess.call(['chmod', '755', bcjob])
        log_name = bcjob +'.log'
        print(bcjob+  '  1>' + log_name  + '  2>&1')
        os.system(bcjob + ' 1>' + log_name+ ' 2>&1')
  else:
    print("sbatch " + bcjob +"\n")
    subprocess.call(['sbatch', bcjob])

  print( "cd " + bin_dir)
  os.chdir(bin_dir)
 
  #print(script_string)

if __name__ == "__main__":

   answers = ask_questions(default_grid="EASEv2")
   configs = get_configs_from_answers(answers)
   for config in configs:
      if 'EASEv' in config['grid_type']:
         make_bcs_ease(config)

