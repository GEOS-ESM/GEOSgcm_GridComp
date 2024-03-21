#!/usr/bin/env python3
#
# source install/bin/g5_modules

import os
from make_bcs_questionary import *
from make_bcs_shared import *


latlon_template = """

ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/360x200 data/MOM5/360x200
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/720x410 data/MOM5/720x410
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/1440x1080 data/MOM5/1440x1080
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM6/72x36 data/MOM6/72x36
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM6/540x458 data/MOM6/540x458
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM6/1440x1080 data/MOM6/1440x1080

if( -e {GRIDNAME}.stdout ) /bin/rm -f {GRIDNAME}.stdout

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

"""

def make_bcs_latlon(config):
  bin_dir = os.getcwd()
  if '/bin' not in bin_dir:
    print(" please run this program in installed bin directory")
    return

  resolution = config['resolution']
  orslv      = config['orslvs']

  account = get_account()
  IMO = '%04d'%config['imo']
  JMO = '%04d'%config['jmo']
  IM  = '%04d'%config['im']
  JM  = '%04d'%config['jm']

  RC = str(config['im']) +'x' +str(config['jm'])+'_DC'

  DATENAME = config['DATENAME']
  POLENAME = config['POLENAME']

  GRIDNAME = 'DC'+IM+'xPC'+JM+'_'+ DATENAME+IMO+'x'+POLENAME+JMO 

  SKIPLAND = config['skipland'] 

  now   = datetime.now()
  tmp_dir =now.strftime("%Y%m%d%H%M%S") 
  tmp_dir=f"{resolution}_{orslv}_{tmp_dir}"
  expdir = config['expdir']
  scratch_dir = expdir+ tmp_dir+'/'+GRIDNAME+'.scratch/'
  log_dir     = expdir+'/'+tmp_dir+'/logs/' + GRIDNAME
  bcjob       = scratch_dir+'/'+GRIDNAME+'.j'

  check_script(expdir, GRIDNAME+'.j')

  os.makedirs(scratch_dir)
  if not os.path.exists(log_dir):
    os.makedirs(log_dir)

  script_template = get_script_head() + latlon_template + get_script_mv(config['grid_type']) 
  script_string = script_template.format(\
           account = account, \
           EXPDIR = config['expdir'], \
           TMP_DIR= tmp_dir, \
           GRIDNAME  = GRIDNAME, \
           GRIDNAME2 = GRIDNAME, \
           SCRATCH_DIR = scratch_dir, \
           bin_dir = bin_dir, \
           MAKE_BCS_INPUT_DIR = config['inputdir'], \
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
           RC = RC,\
           RS = "-Pfafstetter",\
           NCPUS = config['NCPUS'])


  latlon_job = open(bcjob,'wt')
  latlon_job.write(script_string)
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
 
  #print(script_string)

if __name__ == "__main__":

   answers = ask_questions(default_grid="Lat-Lon")
   configs = get_configs_from_answers(answers)
   for config in configs:
      if 'Lat-Lon' in config['grid_type']:
         make_bcs_latlon(config)

