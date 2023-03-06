#!/usr/bin/env python3
#
# source install/bin/g5_modules

import os
from bcs_utils import *


ease_template = """#!/bin/csh -x

#SBATCH --output={EXPDIR}/{OUTDIR}/logs/{BCNAME}.log
#SBATCH --error={EXPDIR}/{OUTDIR}/logs/{BCNAME}.err
#SBATCH --account={account}
#SBATCH --time=12:00:00
#SBATCH --ntasks=28
#SBATCH --job-name={BCNAME}.j
#SBATCH --constraint=sky

cd {BCDIR}

source bin/g5_modules
/bin/ln -s {bin_dir}
mkdir -p til rst data/MOM5 data/MOM6 clsm/plots
ln -s /discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM5/360x200 data/MOM5/360x200
ln -s /discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM5/720x410 data/MOM5/720x410
ln -s /discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM5/1440x1080 data/MOM5/1440x1080
ln -s /discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM6/72x36 data/MOM6/72x36
ln -s /discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM6/1440x1080 data/MOM6/1440x1080

cd data 
ln -s {input_dir} CATCH
cd ../

if( -e DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}.stdout ) /bin/rm -f DC{IM}xPC{JM}_{DATENAME}{IMO}{POLENAME}{JMO}.stdout
setenv MASKFILE {MASKFILE}
limit stacksize unlimited
bin/mkLatLonRaster.x -x {NX} -y {NY}  -t -1 {IM} {JM} >/dev/null
bin/mkLandRaster.x -x {NX} -y {NY} -v -t {NT}

if( {LATLON_OCEAN} == TRUE ) then
    bin/mkLatLonRaster.x -x {NX} -y {NY} -b DE -p PE -t 0 {IMO} {JMO} >/dev/null
    bin/CombineRasters.x -f 0 -t {NT} DE{IMO}xPE{JMO} Pfafstetter >/dev/null
    bin/CombineRasters.x -t {NT} DC{IM}xPC{JM} DE{IMO}xPE{JMO}-Pfafstetter
    setenv OMP_NUM_THREADS 1
    if ({SKIPLAND} != YES) bin/mkCatchParam.x -x {NX} -y {NY} -g DC{IM}xPC{JM}_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
    setenv OMP_NUM_THREADS {NCPUS}
    if ({SKIPLAND} != YES) bin/mkCatchParam.x -x {NX} -y {NY} -g DC{IM}xPC{JM}_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
    chmod 755 bin/create_README.csh
    bin/create_README.csh
endif
if( {TRIPOL_OCEAN} == TRUE ) then
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
    if ({SKIPLAND} != YES) bin/mkCatchParam.x -x {NX} -y {NY} -g DE{IMO}xPE{JMO}_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
    setenv OMP_NUM_THREADS {NCPUS}
    if ({SKIPLAND} != YES) bin/mkCatchParam.x -x {NX} -y {NY} -g DE{IMO}xPE{JMO}_DE{IMO}xPE{JMO}-Pfafstetter -v {lbcsv}
    chmod 755 bin/create_README.csh
    bin/create_README.csh    
endif

/bin/mv clsm  clsm.{IM}x{JM}
/bin/cp til/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.til clsm.{IM}x{JM}
if( {TRIPOL_OCEAN} == TRUE ) /bin/cp til/DC{IM}xPC{JM}_{DATENAME}{IMO}x{POLENAME}{JMO}-Pfafstetter.TRN clsm.{IM}x{JM}
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
/bin/mv    {BCNAME}/clsm/mkCatchParam.log {BCNAME}/logs/mkCatchParam.log
/bin/rm -r {OUTDIR}

"""

def make_cube_bcs(config):
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
  imo = '%04d'%config['imo']
  jmo = '%04d'%config['jmo']
  NC  = '%04d'%config['im']

  bcname = 'CF'+NC+'x6C_'+DATENAME+imo+'x'+POLENAME+jmo
  if config['CUBED_SPHERE_OCEAN'] :
    bcname =  'CF'+NC'x6_CF'+NC+'6C'
    DATENAME = 'CF'
    POLENAME = ''
    imo = NC
    jmo = '6C'
  
  

  now   = datetime.now()
  tmp_dir =now.strftime("%Y%m%d%H%M%S") 
  expdir = config['expdir']
  scratch_dir = expdir+ tmp_dir+'/'+bcname+'.scratch/'
  log_dir     = expdir+'/'+tmp_dir+'/logs/'
  bcjob       = scratch_dir+'/'+bcname+'.j'

  if os.path.exists(bcjob):
    print('please remove the run temprory directory: ' + expdir+'/'+ tmp_dir) 
    return

  os.makedirs(scratch_dir)
  os.makedirs(log_dir)

  job_script = ease_template.format(\
           account = account, \
           EXPDIR = config['expdir'], \
           OUTDIR = tmp_dir, \
           BCNAME = bcname, \
           bin_dir = bin_dir, \
           input_dir = config['inputdir'], \
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

