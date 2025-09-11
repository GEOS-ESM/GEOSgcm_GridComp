#!/usr/bin/env python3
#
# source install/bin/g5_modules

import os
from make_bcs_questionary import *
from make_bcs_shared import * 
import re
from datetime import datetime
import subprocess

def _resolve_versions(bcs_version: str):
    """
    Returns (mom6_bathy_version, topo_version) as 'v1' or 'v2'
    Rule: v14+ -> v2; else -> v1. Non-numeric (NL3, ICA, GM4...) -> v1.
    """
    v = (bcs_version or "").strip()
    m = re.match(r'[vV]?(\d+)', v)
    vnum = int(m.group(1)) if m else None
    use_v2 = vnum is not None and vnum >= 14
    return ("v2" if use_v2 else "v1", "v2" if use_v2 else "v1")

def _build_mom6_link_lines(inputdir: str, preferred_version: str) -> str:
    """
    Build csh lines to symlink MOM6 datasets.
    Prefer `preferred_version` (v1 or v2) per size; if missing, fall back to the other version.
    Emits an echo note on fallback; warns if a size is missing in both.
    """
    sizes = ["72x36", "540x458", "1440x1080"]
    other = "v1" if preferred_version == "v2" else "v2"

    lines = ['if ( ! -d data/MOM6 ) mkdir -p data/MOM6']
    for size in sizes:
        src_pref = os.path.join(inputdir, "ocean", "MOM6", preferred_version, size)
        src_other = os.path.join(inputdir, "ocean", "MOM6", other,            size)

        if os.path.isdir(src_pref):
            src = src_pref
            note = ""
        elif os.path.isdir(src_other):
            src = src_other
            note = f'echo "NOTE: MOM6 {preferred_version}/{size} not found; using {other}/{size}"'
        else:
            lines.append(f'echo "WARNING: MOM6 {size} missing in both v1 and v2; skipping"')
            continue

        lines.append(f'if ( -e data/MOM6/{size} ) /bin/rm -f data/MOM6/{size}')
        if note:
            lines.append(note)
        lines.append(f'ln -s {src} data/MOM6/{size}')
    return "\n".join(lines)

cube_template = """

ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/360x200 data/MOM5/360x200
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/720x410 data/MOM5/720x410
ln -s {MAKE_BCS_INPUT_DIR}/ocean/MOM5/1440x1080 data/MOM5/1440x1080
{MOM6_LINK_LINES}


if( -e CF{NC}x6C{SGNAME}_{DATENAME}{IMO}x{POLENAME}{JMO}.stdout ) /bin/rm -f CF{NC}x6C{SGNAME}_{DATENAME}{IMO}x{POLENAME}{JMO}.stdout

bin/mkCubeFVRaster.x -x {NX} -y {NY} {SGPARAM} {STRETCH} {NC} >/dev/null 
bin/mkLandRaster.x -x {NX} -y {NY} -v -t {NT}

if( {LATLON_OCEAN} == True ) then
    bin/mkLatLonRaster.x -x {NX} -y {NY} -b DE -p PE -t 0 {IMO} {JMO} >/dev/null
    bin/CombineRasters.x -f 0 -t {NT} DE{IMO}xPE{JMO} Pfafstetter >/dev/null
    bin/CombineRasters.x -t {NT} CF{NC}x6C{SGNAME} DE{IMO}xPE{JMO}-Pfafstetter
    setenv OMP_NUM_THREADS {NCPUS}
    if ( {SKIPLAND} != True ) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C{SGNAME}_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
    setenv OMP_NUM_THREADS 1
    chmod 755 bin/create_README.csh
    bin/create_README.csh
endif

if( {TRIPOL_OCEAN} == True ) then
    bin/mkMOMAquaRaster.x -x {NX} -y {NY} -w {OCEAN_VERSION} data/{MOM_VERSION}/{imo}x{jmo}/MAPL_Tripolar.nc > /dev/null
    /bin/cp til/Pfafstetter.til til/Pfafstetter-ORIG.til
    /bin/cp rst/Pfafstetter.rst rst/Pfafstetter-ORIG.rst
    bin/FillMomGrid.x -f 0 -g Pfafstetter-M {OCEAN_VERSION}{DATENAME}{IMO}x{POLENAME}{JMO} Pfafstetter data/{MOM_VERSION}/{imo}x{jmo}/MAPL_Tripolar.nc 
    /bin/mv til/Pfafstetter-M.til til/Pfafstetter.til
    /bin/mv rst/Pfafstetter-M.rst rst/Pfafstetter.rst
    bin/CombineRasters.x -f 0 -t {NT} {OCEAN_VERSION}{DATENAME}{IMO}x{POLENAME}{JMO} Pfafstetter >/dev/null
    bin/CombineRasters.x -t {NT} CF{NC}x6C{SGNAME} {OCEAN_VERSION}{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter
    bin/mk_runofftbl.x -g CF{NC}x6C{SGNAME}_{OCEAN_VERSION}{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter -v {lbcsv}

    if ({SKIPLAND} != True) then
      bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C{SGNAME}_{OCEAN_VERSION}{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter -v {lbcsv} -p no
      bin/ExtractBCsFromOrig.py {BCS_DIR}  {lbcsv} CF{NC}x6C{SGNAME} {OCEAN_VERSION}{DATENAME}{IMO}x{POLENAME}{JMO}
    endif
 
    chmod 755 bin/create_README.csh
    bin/create_README.csh
endif

if( {CUBED_SPHERE_OCEAN} == True ) then
    if ( {IS_STRETCHED} == True ) then
       bin/mkCubeFVRaster.x -x {NX} -y {NY} {STRETCH} {NC} >/dev/null 
    endif
    bin/CombineRasters.x -f 0 -t {NT} CF{NC}x6C Pfafstetter >/dev/null
    bin/CombineRasters.x -t {NT} {SGPARAM} CF{NC}x6C CF{NC}x6C-Pfafstetter
    setenv OMP_NUM_THREADS {NCPUS}
    if ({SKIPLAND} != True) bin/mkCatchParam.x -x {NX} -y {NY} -g CF{NC}x6C{SGNAME}_CF{NC}x6C-Pfafstetter -v {lbcsv}
    setenv OMP_NUM_THREADS 1
    chmod 755 bin/create_README.csh
    bin/create_README.csh
endif

"""
def make_bcs_cube(config):

  bin_dir = os.getcwd()
  if '/bin' not in bin_dir:
    print(" please run this program in installed bin directory")
    return
  grid_type  = config['grid_type']
  
  if grid_type not in ["Stretched_CS", "Cubed-Sphere"] :
     print('This should be a Cubed-Sphere or Stretched Cubed-Sphere (Stretched_CS) grid')
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

  IS_STRETCHED = False
  STRETCH = ""
  SG = ""
  SGNAME = ""
  SGPARAM = ""
  if 'Stretched_CS' in grid_type:
     IS_STRETCHED = True
     SG       = config['SG']
     if SG == 'SG001' :
        STRETCH  = '-F 2.5 -X -98.35 -Y 39.5' 
     if SG == 'SG002' :
        STRETCH  = '-F 3.0 -X -98.35 -Y 39.5'
     SGNAME = '-'+''.join(SG)
     SGPARAM = '-s '+''.join(SG)


  DATENAME = config['DATENAME']
  POLENAME = config['POLENAME']
  OCEAN_VERSION = config['OCEAN_VERSION']
  SKIPLAND = config['skipland']
      
  GRIDNAME = 'CF'+NC+'x6C'+SGNAME+'_'+DATENAME+IMO+'x'+POLENAME+JMO

  if config['CUBED_SPHERE_OCEAN'] :
    DATENAME = 'CF'
    POLENAME = ''
    IMO = NC
    JMO = '6C'
    GRIDNAME   =  'CF'+ NC+'x6C'+SGNAME+'_CF'+NC+'x6C' 

  if config['TRIPOL_OCEAN'] :
      GRIDNAME = 'CF'+ NC+'x6C'+SGNAME+'_'+OCEAN_VERSION+IMO+'x'+JMO

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

  MOM6_BATHY_VERSION, TOPO_VERSION = _resolve_versions(config['lbcsv'])
  MOM6_LINK_LINES = _build_mom6_link_lines(config['inputdir'], MOM6_BATHY_VERSION)

  script_template = get_script_head() + cube_template + get_script_mv(config['grid_type'])

  script_string = script_template.format(\
           account = account, \
           EXPDIR = config['expdir'], \
           TMP_DIR = tmp_dir, \
           GRIDNAME  = GRIDNAME, \
           SCRATCH_DIR = scratch_dir, \
           bin_dir = bin_dir, \
           MAKE_BCS_INPUT_DIR = config['inputdir'], \
           BCS_DIR  = config['bcs_dir'], \
           DATENAME = DATENAME, \
           POLENAME = POLENAME, \
           OCEAN_VERSION = OCEAN_VERSION, \
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
           SGNAME  = SGNAME, \
           SGPARAM = SGPARAM, \
           TOPO_VERSION = TOPO_VERSION,
           mom6_bathy_version = MOM6_BATHY_VERSION,
           MOM6_LINK_LINES = MOM6_LINK_LINES,                              
           IS_STRETCHED = IS_STRETCHED, \
           NCPUS = config['NCPUS'])

  cube_job = open(bcjob,'wt')
  cube_job.write(script_string)
  cube_job.close()

  interactive = os.getenv('SLURM_JOB_ID', default = None)
  if ( interactive ) :
     if resolution in ['c1080' ,'c1536', 'c2160', 'c2880', 'c3072','c5760'] :
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
       
  print( "cd " + bin_dir)
  os.chdir(bin_dir)
 
  #print(stretched_script_string)

if __name__ == "__main__":

   answers = ask_questions()
   configs = get_configs_from_answers(answers)
   for config in configs:
       if config['grid_type'] in ["Stretched_CS", "Cubed-Sphere"]:
           make_bcs_cube(config)   
