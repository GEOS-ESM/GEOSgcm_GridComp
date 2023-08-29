#!/usr/bin/env python3
#
# source install/bin/g5_modules

import os
from make_bcs_questionary import *
from make_bcs_shared import * 

stretched_cube_template = """

if ( {STEP1} == True ) then
  ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/360x200 data/MOM5/360x200
  ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/720x410 data/MOM5/720x410
  ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/1440x1080 data/MOM5/1440x1080
  ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM6/72x36 data/MOM6/72x36
  ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM6/540x458 data/MOM6/540x458
  ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM6/1440x1080 data/MOM6/1440x1080

  if( -e CF{NC}x6C-{SG}_{DATENAME}{IMO}x{POLENAME}{JMO}.stdout ) /bin/rm -f CF{NC}x6C-{SG}_{DATENAME}{IMO}x{POLENAME}{JMO}.stdout

endif 

if ( {STEP1} == True ) then
  bin/mkCubeFVRaster.x -x {NX} -y {NY} -s {SG} {STRETCH} {NC} >/dev/null 
  bin/mkLandRaster.x -x {NX} -y {NY} -v -t {NT}
endif

if( {LATLON_OCEAN} == True ) then

   if ( {STEP1} == True ) then 
      bin/mkLatLonRaster.x -x {NX} -y {NY} -b DE -p PE -t 0 {IMO} {JMO} >/dev/null
      bin/CombineRasters.x -f 0 -t {NT} DE{IMO}xPE{JMO} Pfafstetter >/dev/null
      bin/CombineRasters.x -t {NT} CF{NC}x6C-{SG} DE{IMO}xPE{JMO}-Pfafstetter
      setenv OMP_NUM_THREADS 1
      if ( {SKIPLAND} != True ) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C-{SG}_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
   endif

   if ( {STEP2} == True ) then 
      setenv OMP_NUM_THREADS {NCPUS}
      if ( {SKIPLAND} != True ) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C-{SG}_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
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
      bin/CombineRasters.x -t {NT} CF{NC}x6C-{SG} {DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter
      bin/mk_runofftbl.x CF{NC}x6C-{SG}_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter
      setenv OMP_NUM_THREADS 1
      if ({SKIPLAND} != True) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C-{SG}_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter -v {lbcsv}
   endif

   if ( {STEP2} == True ) then 
      setenv OMP_NUM_THREADS {NCPUS}
      if ({SKIPLAND} != True) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C-{SG}_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter -v {lbcsv}
      chmod 755 bin/create_README.csh
      bin/create_README.csh
   endif
endif

if( {CUBED_SPHERE_OCEAN} == True ) then
   if ( {STEP1} == True ) then 
      bin/mkCubeFVRaster.x -x {NX} -y {NY} {STRETCH} {NC} >/dev/null 
      bin/CombineRasters.x -f 0 -t {NT} CF{NC}x6C Pfafstetter >/dev/null
      bin/CombineRasters.x -t {NT} -s {SG} CF{NC}x6C CF{NC}x6C-Pfafstetter
      setenv OMP_NUM_THREADS 1
      if ({SKIPLAND} != True) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C-{SG}_CF{NC}x6C-Pfafstetter -v {lbcsv}
   endif

   if ( {STEP2} == True ) then 
      setenv OMP_NUM_THREADS {NCPUS}
      if ({SKIPLAND} != True) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C-{SG}_CF{NC}x6C-Pfafstetter -v {lbcsv} 
      chmod 755 bin/create_README.csh
      bin/create_README.csh
   endif
endif

"""
def make_bcs_stretched_cube(config):

  bin_dir = os.getcwd()
  if 'install/bin' not in bin_dir:
    print(" please run this program in installed bin directory")
    return

  grid_type  = config['grid_type']
  if 'Stretched_Cubed' not in grid_type :
     print('This is not a Stretched_Cubed-Sphere grid')
     return

  resolution = config['resolution']
  orslv      = config['orslvs']

  account = get_account()
  IMO = '%04d'%config['imo']
  JMO = '%04d'%config['jmo']
  NC  = '%04d'%config['im']
  imo = config['imo']
  jmo = config['jmo']
  nc  = config['im']
  nc6 = nc*6
  RC = str(nc)+'x'+str(nc6)

  if resolution in ['c270','c540','c1080', 'c2160'] :
      STRETCH  = '-F 2.5 -X -98.35 -Y 39.5' 
  if resolution in ['c1536'] :
      STRETCH  = '-F 3.0 -X -98.35 -Y 39.5'

  DATENAME = config['DATENAME']
  POLENAME = config['POLENAME']
  SKIPLAND = config['skipland']
      
  if resolution in ['c270', 'c540', 'c1080', 'c2160'] :
      SG   = 'SG001'
  if resolution in ['c1536'] :
      SG   = 'SG002'

  GRIDNAME = 'CF'+NC+'x6C-'+SG+'_'+DATENAME+IMO+'x'+POLENAME+JMO

  if config['CUBED_SPHERE_OCEAN'] :
    IMO = NC
    JMO = '6C'
    GRIDNAME   =  'CF'+ NC+'x6C-'+SG+'_CF'+NC+'x6C' 

  now   = datetime.now()
  tmp_dir =now.strftime("%Y%m%d%H%M%S") 
  tmp_dir=f"{resolution}_{orslv}_{tmp_dir}"
  expdir = config['expdir']
  scratch_dir = expdir+ tmp_dir+'/'+GRIDNAME+'.scratch/'
  log_dir     = expdir+'/'+tmp_dir+'/logs/'+ GRIDNAME
  bcjob       = scratch_dir+'/'+GRIDNAME+'.j'

  check_script(expdir, GRIDNAME+'.j')

  os.makedirs(scratch_dir)
  if not os.path.exists(log_dir):
    os.makedirs(log_dir)

  STEP1 = True
  STEP2 = True
  GRIDNAME2 = GRIDNAME
  script_template = get_script_head() + stretched_cube_template + get_script_mv(config['grid_type'])
  if resolution in ['c1080' ,'c1536', 'c2160'] :
     STEP1 = True
     STEP2 = False
     script_template = get_script_head() + stretched_cube_template 

  stretched_script_string = script_template.format(\
           account = account, \
           EXPDIR = config['expdir'], \
           TMP_DIR = tmp_dir, \
           GRIDNAME  = GRIDNAME, \
           GRIDNAME2 = GRIDNAME2, \
           STEP1   = STEP1, \
           STEP2   = STEP2, \
           SCRATCH_DIR = scratch_dir, \
           bin_dir = bin_dir, \
           MAKE_BCS_INPUT_DIR = config['inputdir'], \
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
           RS = '-Pfafstetter',\
           RC = RC,\
           SG = SG,\
           STRETCH = STRETCH, \
           NCPUS = config['NCPUS'])

  stretched_cube_job = open(bcjob,'wt')
  stretched_cube_job.write(stretched_script_string)
  stretched_cube_job.close()

  if resolution in ['c1080' ,'c1536', 'c2160'] :
     STEP1 = False
     STEP2 = True
     GRIDNAME2 = GRIDNAME+'-2'
     script_template = get_script_head() + stretched_cube_template + get_script_mv(config['grid_type'])
     stretched_script_string = script_template.format(\
           account = account, \
           EXPDIR = config['expdir'], \
           TMP_DIR = tmp_dir, \
           GRIDNAME  = GRIDNAME, \
           GRIDNAME2 = GRIDNAME2, \
           STEP1   = STEP1, \
           STEP2   = STEP2, \
           SCRATCH_DIR = scratch_dir, \
           bin_dir = bin_dir, \
           MAKE_BCS_INPUT_DIR = config['inputdir'], \
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
           SG = SG,\
           STRETCH = STRETCH, \
           RS = '-Pfafstetter',\
           NCPUS = config['NCPUS'])

     stretched_cube_job = open(bcjob+'-2','wt')
     stretched_cube_job.write(stretched_script_string)
     stretched_cube_job.close()

  interactive = os.getenv('SLURM_JOB_ID', default = None)
  if ( interactive ) :
     if resolution in ['c1080' ,'c1536', 'c2160'] :
        exit("resolution is too high for interactive mode")
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
    if resolution in ['c1080' ,'c1536', 'c2160'] :
      print("sbatch " + bcjob+'-2' + " depending on " + bcjob + "\n")
      subprocess.call(['sbatch', '--dependency=afterok:'+jobid, bcjob+'-2'])
      print()
       
  print( "cd " + bin_dir)
  os.chdir(bin_dir)
 
  #print(stretched_script_string)

if __name__ == "__main__":

   answers = ask_questions()
   configs = get_configs_from_answers(answers)
   for config in configs:
      if 'Stretched_Cubed-Sphere' in config['grid_type']:
         make_bcs_stretched_cube(config)

