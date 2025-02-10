#!/usr/bin/env python3
#
#
import os
import subprocess
import shlex
import ruamel.yaml
import shutil
import questionary
import pathlib

smoothmap ={ 'C12'  : 773.91,
             'C24'  : 386.52,
             'C48'  : 193.07,
             'C90'  : 102.91,
             'C180' : 51.44,
             'C270' : 40,
             'C360' : 25.71,
             'C540' : 19,
             'C720' : 12.86,
             'C1080': 10,
             'C1120': 8.26,
             'C1440': 6.43,
             'C1539': 6,
             'C2160': 5,
             'C2880': 3.21
            }

def get_script_topo(answers) :
  head =  """#!/bin/csh -x

#SBATCH --output=topo.log
#SBATCH --error=topo.err
#SBATCH --account={account}
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --job-name=topo.j
"""

  constraint = '#SBATCH --constraint="[mil|cas]"'
  #if 'TRUE' not in BUILT_ON_SLES15:
  #   constraint = "#SBATCH --constraint=sky"

  topo_template = head + constraint + """

echo "-----------------------------" 
echo "make_topo starts date/time" 
echo `date` 
echo "-----------------------------" 

if( ! -d bin ) then
  /bin/ln -s {bin_dir}
endif

source bin/g5_modules  

if ( ! -e landm_coslat.nc ) then
  /bin/ln -s bin/landm_coslat.nc landm_coslat.nc
endif 

set source_topo = gmted_intel
set cutoff = 50
set smooths = {SMOOTHMAP}
set resolutions = {RESOLUTIONS}
set lowres  = 0
set highres = 0
set SG001 = ( 270 540 1080 2160 )
set SG002 = ( 1539 )

foreach res ($resolutions)
  if ($res < $cutoff) then
    set lowres  = 1
  endif
  if ($res > $cutoff) then
    set highres = 1
  endif
end

if ( $lowres == 1 ) then
cat << _EOF_ > bin_to_cube.nl
&binparams
  raw_latlon_data_file='{raw_latlon_data}'
  output_file='c360.gmted_fixedanarticasuperior.nc'
  ncube=360
/
_EOF_
   bin/bin_to_cube.x
endif

if ( $highres == 1 ) then
cat << _EOF_ > bin_to_cube.nl
&binparams
  raw_latlon_data_file='{raw_latlon_data}'
  output_file='c3000.gmted_fixedanarticasuperior.nc'
  ncube=3000
/
_EOF_
   bin/bin_to_cube.x
endif

@ count = 1
foreach im ($resolutions)
  @ jm = ($im * 6)
  set output_dir = output_$im

  if ( ! -e $output_dir ) then
    mkdir $output_dir
  endif

  set DO_SCHMIDT = ''
  set TARGET_LON = ''
  set TARGET_LAT = ''
  set STRETCH_FACTOR = ''

  foreach sg1 ($SG001)
     if ($im == $sg1) then
       set DO_SCHMIDT = 'DO_SCHMIDT:  .true.'
       set TARGET_LON = 'TARGET_LON:  -98.35 '
       set TARGET_LAT = 'TARGET_LAT:  39.5 '
       set STRETCH_FACTOR = 'STRETCH_FACTOR: 2.5 '
     endif
  end

  foreach sg2 ($SG002)
     if ($im == $sg2) then
       set DO_SCHMIDT = 'DO_SCHMIDT:  .true.'
       set TARGET_LON = 'TARGET_LON:   -98.35'
       set TARGET_LAT = 'TARGET_LAT:   39.5'
       set STRETCH_FACTOR = 'STRETCH_FACTOR: 3.0'
     endif
  end

  set config_file = GenScrip.rc
  set output_grid = PE${{im}}x${{jm}}-CF
  set scriptfile  = ${{output_grid}}.nc4

cat << _EOF_ > ${{config_file}}
CUBE_DIM: $im
output_scrip: ${{scriptfile}}
output_geos: c${{im}}_coords.nc4
${{DO_SCHMIDT}}
${{TARGET_LON}}
${{TARGET_LAT}}
${{STRETCH_FACTOR}}
_EOF_

   cat ${{config_file}}
   mpirun -np 6 bin/generate_scrip_cube_topo.x

   rm ${{config_file}}

   if ($im < $cutoff) then
     set intermediate_cube = c360.gmted_fixedanarticasuperior.nc
   else
     set intermediate_cube = c3000.gmted_fixedanarticasuperior.nc
   endif
   set jmax_segments=''
   if ($im == 2880) then
      set jmax_segments = --jmax_segments=32
   endif
   set output_grid = PE${{im}}x${{jm}}-CF
   bin/cube_to_target.x --grid_descriptor_file=$scriptfile --intermediate_cs_name=$intermediate_cube --output_data_directory=$output_dir --smoothing_scale=${{smooths[$count]}} --name_email_of_creator='gmao' --fine_radius=0 --output_grid=$output_grid --source_data_identifier=$source_topo $jmax_segments

   rm $scriptfile

   #convert to gmao
   cd $output_dir
   set arr = `ls *.nc`
   echo $arr
   ../bin/scrip_to_restart_topo.py -i $arr -o gwd_internal_rst
   ../bin/convert_to_gmao_output_topo.x -i $arr --im $im
   cd ..
   @ count = $count + 1
end
"""
  account = get_account()
  SMOOTHMAP = '( ' 
  RESOLUTIONS = '( '
  for res in answers['resolutions']:
      SMOOTHMAP = SMOOTHMAP + str(smoothmap[res]) + ' '
      RESOLUTIONS = RESOLUTIONS + str(res)[1:] + ' '

  SMOOTHMAP   = SMOOTHMAP + ' )'
  RESOLUTIONS = RESOLUTIONS + ' )'

  script_string = topo_template.format(\
       account = account, \
       bin_dir = answers['bin_dir'], \
       raw_latlon_data = answers['path_latlon']+ "/gmted_fixed_anartica_superior_caspian.nc4", \
       SMOOTHMAP   = SMOOTHMAP, \
       RESOLUTIONS = RESOLUTIONS )
  out_dir = answers['out_dir']
  pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
  topojob = out_dir+'/topo.j'
  topo_job = open(topojob,'wt')
  topo_job.write(script_string)
  topo_job.close()
  subprocess.call(['chmod', '755', topojob])

  print("\nJob script topo.j has been generated in "  + out_dir + "\n")
  
def get_user():
   cmd = 'whoami'
   p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
   (user, err) = p.communicate()
   p_status = p.wait()
   user = user.decode().split()
   return user[0]

def get_account():
   cmd = 'id -gn'
   p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
   (accounts, err) = p.communicate()
   p_status = p.wait()
   accounts = accounts.decode().split()
   return accounts[0]

def ask_questions():

   # See remap_utils.py for definitions of "choices", "message" strings, and "validate" lists
   # that are used multiple times.
   user_name = get_user()

   questions = [
       {
            "type": "path",
            "name": "bin_dir",
            "message": "Enter the root path of the bin:\n",
            "default": "./"
        },

       {
            "type": "path",
            "name": "out_dir",
            "message": "Enter the root path of the bin:\n",
            "default": "/discover/nobackup/"+user_name+"/BCS_TOPO/"
        },

       {
            "type": "path",
            "name": "path_latlon",
            "message": "Enter the path contains gmted_fixed_anartica_superior_caspian.nc4 :\n",
            "default": "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/topo/v1/"
        },

                
        {   
            "type": "checkbox",
            "name": "resolutions",
            "message": "Select resolutions: \n",
            "choices": ["C12","C24", "C48", "C90", "C180", "C360", "C720", "C1120", "C1440", "C2880", "SG001","SG002"]
        },

        {
            "type": "checkbox",
            "name": "SG001",
            "message": "Select resolution of SG001 grid: \n",
            "choices": ['C270', 'C540', 'C1080', 'C2160'], 
            "when": lambda x : 'SG001' in x.get('resolutions'),
        },
        {
            "type": "checkbox",
            "name": "SG002",
            "message": "Select resolution of SG002 grid: \n",
            "choices": ['C1539'],
            "when": lambda x : 'SG002' in x.get('resolutions'),
         },

        ]
   answers = questionary.prompt(questions)
   answers['bin_dir'] = os.path.abspath(answers['bin_dir'])
   if 'SG001' in answers:
      answers['resolutions'].remove('SG001')
      answers['resolutions'] = answers['resolutions'] + answers['SG001']
   if 'SG002' in answers:
      answers['resolutions'].remove('SG002')
      answers['resolutions'] = answers['resolutions'] + answers['SG002']
   print(answers['resolutions'])
   return answers

    
if __name__ == '__main__' :

   answers = ask_questions()
   get_script_topo(answers)
