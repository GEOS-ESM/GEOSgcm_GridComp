#!/usr/bin/env python3
#
# run as ./make_topo.py
#
# NOTE run time setups:
# c180, c360, c720, c1080, c1440, c2880, c1536, c1120, c2160 all need 1 node : 1h
# c5760, c540, c270 and c48 require 2 nodes :  4h 
# c90 1 node  : 3h
# c24 2 nodes : 8h
# c12 2 nodes : 19h - has to be run with qos=long
#
import os
import subprocess
import shlex
import questionary
import pathlib

#Bigger smoothing_scale ⇒ stronger smoothing ⇒ smaller GWD (less roughness).
smoothmap ={ 
#  # uniform
             'C12'  : 512.0,
             'C24'  : 305.0,
             'C48'  : 166.0,
             'C90'  : 96.2,
             'C180' : 51.21,
             'C360' : 28.95,
             'C720' : 19.5,
             'C1120': 8.26,
             'C1440': 12.0,
             'C2880': 3.285,
             'C5760': 3.0, 
#  # stretched (SG001/SG002)
             'C270' : 100.0,
             'C540' : 53.3,
             'C1080': 17.0,
             'C1536': 26.8,
             'C2160': 2.98
            }
## Tuning parameter for GWD amplitude for stretched grids this affects step in laplacian iterations.
alpha ={ 'C270' : 2.9,
         'C540' : 2.73,
         'C1080': 7.0,
         'C1536': 12.1,
         'C2160': 14.25,
         # regular grids do not use alpha!!!
        }

def get_script_topo(answers) :
  head =  """#!/bin/csh -x

#SBATCH --output=topo.log
#SBATCH --error=topo.err
#SBATCH --account={account}
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --job-name=topo_{res_tag}.j
"""

  constraint = '#SBATCH --constraint="[mil]"'

  if len(answers['resolutions']) == 1:
      res_tag = answers['resolutions'][0].lower()
  else:
      res_tag = "_".join(r.lower() for r in answers['resolutions'])

  topo_template = head + constraint + """

echo "-----------------------------"
echo "make_topo starts date/time"
echo `date`
echo "-----------------------------"

if( ! -d bin ) then
  /bin/ln -s {bin_dir}
endif

source bin/g5_modules
module load nco
module load cdo

if ( ! -e landm_coslat.nc ) then
  /bin/ln -s bin/landm_coslat.nc landm_coslat.nc
endif

set source_topo = gmted_intel
set smooths = {SMOOTHMAP}
set resolutions = {RESOLUTIONS}
set SG001 = ( 270 540 1080 2160 )
set SG002 = ( 1536 )

# Generate a single high-resolution intermediate cube (3000) for ALL resolutions
if ( ! -e c3000.gmted_fixedanarticasuperior.nc ) then
cat << _EOF_ > bin_to_cube.nl
&binparams
  raw_latlon_data_file='{raw_latlon_data}'
  output_file='c3000.gmted_fixedanarticasuperior.nc'
  ncube=3000
/
_EOF_
  bin/bin_to_cube.x
else
  echo "Reusing existing c3000.gmted_fixedanarticasuperior.nc"
endif

# Use the same intermediate for all resolutions
set intermediate_cube = c3000.gmted_fixedanarticasuperior.nc
if ( ! -e $intermediate_cube ) then
  echo "ERROR: Missing $intermediate_cube after generation"; exit 2
endif

@ count = 1
foreach im ($resolutions)
  @ jm = ($im * 6)
  set output_dir = output_${{im}}

  if ( ! -e $output_dir ) then
    mkdir $output_dir
  endif
  set alpha     = ""
  set ALPHALINE = ""

  set DO_SCHMIDT = ''
  set TARGET_LON = ''
  set TARGET_LAT = ''
  set STRETCH_FACTOR = ''
  set grid_type = ''

  foreach sg1 ($SG001)
     if ($im == $sg1) then
       set DO_SCHMIDT = 'DO_SCHMIDT:  true'
       set TARGET_LON = 'TARGET_LON:  -98.35 '
       set TARGET_LAT = 'TARGET_LAT:  39.5 '
       set STRETCH_FACTOR = 'STRETCH_FACTOR: 2.5 '
       set grid_type = sg001
     endif
  end

  foreach sg2 ($SG002)
     if ($im == $sg2) then
       set DO_SCHMIDT = 'DO_SCHMIDT:  true'
       set TARGET_LON = 'TARGET_LON:   -98.35'
       set TARGET_LAT = 'TARGET_LAT:   39.5'
       set STRETCH_FACTOR = 'STRETCH_FACTOR: 3.0'
       set grid_type = sg002
     endif
  end

  set config_file = GenScrip.yaml
  set output_grid = PE${{im}}x${{jm}}-CF
  set scriptfile  = ${{output_grid}}.nc4
  set smoothing_scale = ${{smooths[$count]}}

   if ( $im == 5760 ) then
        set extra_cli = "-l 13"   # run the Laplacian for 13 cycles
   else
        set extra_cli = ""
   endif

   # fill alpha only for stretched cases (injected switch from Python)
   {alpha_switch}

   # after your alpha_switch and smoothing_scale logic
   if ( "$alpha" == "" ) then
     set ALPHALINE = ''
   else
     set ALPHALINE = "ALPHA: $alpha"
   endif

cat << _EOF_ > ${{config_file}}
CUBE_DIM: $im
output_scrip: ${{scriptfile}}
output_geos: c${{im}}_coords.nc4
${{DO_SCHMIDT}}
${{TARGET_LON}}
${{TARGET_LAT}}
${{STRETCH_FACTOR}}
${{ALPHALINE}}
_EOF_

   cat ${{config_file}}
   mpirun -np 6 bin/generate_scrip_cube_topo.x

   # --- Add error-checking here ---
   if ( $status != 0 ) then
       echo "ERROR: generate_scrip_cube_topo.x failed (exit $status)"
       exit 1
   endif
   if ( ! -e ${{scriptfile}} ) then
       echo "ERROR: descriptor ${{scriptfile}} not created"
       exit 1
   endif

   echo "IM=$im  -> using intermediate: $intermediate_cube"
  
   #--------------------------------------------------------
   # Build jmax/rrfac flags
   #--------------------------------------------------------

       # --- rrfac_max  = ceil( max(rrfac) ) ------------------
       set rr = `cdo -s infon $scriptfile | \
                 awk '/rrfac/ {{v=$(NF-2); printf("%d",(v>int(v)?int(v)+1:int(v)));}}'`
       if ( "$rr" != "" ) then
        set rrfac = "--rrfac_max=$rr"
      else
          set rrfac = "--rrfac_max=1"
      endif

        bin/cube_to_target.x \
            --grid_descriptor_file=$scriptfile \
            --intermediate_cs_name=$intermediate_cube \
            --output_data_directory=$output_dir \
            --smoothing_scale=$smoothing_scale \
            --name_email_of_creator=gmao \
            --fine_radius=0 \
            --output_grid=$output_grid \
            --source_data_identifier=$source_topo \
            $rrfac $extra_cli

      # Safety check after cube_to_target.x
      if ( $status != 0 ) then
          echo "ERROR: cube_to_target.x failed (exit $status). Check stdout above."
          exit 1
      endif

      ls $output_dir/*.nc >& /dev/null
      if ( $status != 0 ) then
          echo "ERROR: cube_to_target.x returned 0, but wrote no *.nc files."
          exit 1
      endif


      #rm $scriptfile
      rm ${{config_file}}

     # convert to gmao
     cd $output_dir
     
     # choose exactly the PE* file (ignore any topo_smooth*.nc written by -z)
     set pe = `ls -1t PE${{im}}x${{jm}}*.nc | head -1`
     if ( "$pe" == "" ) then
       echo "ERROR: no PE file found for IM=$im JM=$jm"
       exit 1
     endif
     
     if ( "$grid_type" != "" ) then
       ../bin/scrip_to_restart_topo.py -i $pe -o gwd_internal_rst -g $grid_type
     else
       ../bin/scrip_to_restart_topo.py -i $pe -o gwd_internal_rst
     endif
     
     ../bin/convert_to_gmao_output_topo.x -i $pe --im $im
     cd ..
   @ count = $count + 1
end
"""
  account = get_account()
  SMOOTHMAP = '( '
  RESOLUTIONS = '( '
  for res in answers['resolutions']:
      SMOOTHMAP += str(smoothmap[res]) + ' '
      RESOLUTIONS += str(res)[1:] + ' '
  SMOOTHMAP   = SMOOTHMAP + ' )'
  RESOLUTIONS = RESOLUTIONS + ' )'

  # Explicit csh switch-case for alpha per resolution
  alpha_switch = "switch ($im)\n"
  for res in answers['resolutions']:
      numeric_res = res[1:]  # Remove leading "C"
      if res in alpha:
          alpha_switch += f"  case {numeric_res}:\n"
          alpha_switch += f"    set alpha = {alpha[res]}\n"
          alpha_switch += "    breaksw\n"
  alpha_switch += "  default:\n"
  alpha_switch += "    set alpha = ''\n"
  alpha_switch += "endsw\n"  
  

  script_string = topo_template.format(\
       account = account, \
       bin_dir = answers['bin_dir'], \
       raw_latlon_data = answers['path_latlon']+ "/gmted_fixed_anartica_superior_caspian.nc4", \
       SMOOTHMAP   = SMOOTHMAP, \
       alpha_switch = alpha_switch, \
       RESOLUTIONS = RESOLUTIONS, \
       res_tag = res_tag )
  out_dir = answers['out_dir']
  pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
  if len(answers['resolutions']) == 1:
      res_tag = answers['resolutions'][0].lower()
  else:
      res_tag = "_".join(r.lower() for r in answers['resolutions'])

  topojob = f"{out_dir}/topo_{res_tag}.j"

  topo_job = open(topojob,'wt')
  topo_job.write(script_string)
  topo_job.close()
  subprocess.call(['chmod', '755', topojob])

  print(f"\nJob script {os.path.basename(topojob)} has been generated in {out_dir}\n")


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
            "message": "Enter the path of the output directory:\n",
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
            "choices": ["C12","C24", "C48", "C90", "C180", "C360", "C720", "C1120", "C1440", "C2880", "C5760", "SG001","SG002"]
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
            "choices": ['C1536'],
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
