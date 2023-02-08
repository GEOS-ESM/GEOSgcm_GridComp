#!/usr/bin/env python3
#
# source install/bin/g5_modules
#
# Newer GEOS code should load a module with GEOSpyD Python3 if not run:
#   module load python/GEOSpyD/Min4.10.3_py3.9
#

import os
import socket
import subprocess
import shlex
import ruamel.yaml
import shutil
import questionary
import glob
from datetime import datetime

def get_account():
   cmd = 'id -gn'
   p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
   (accounts, err) = p.communicate()
   p_status = p.wait()
   accounts = accounts.decode().split()
   return accounts[0]

def get_config_from_answers(answers):

   grid_type  = answers['grid_type']
   lbcsv      = answers['bcs_version']
   resolution = answers['resolution']
   orslvs     = answers['orslvs']

   im = {}
   imo = {}
   jm = {}
   jmo = {}
   NX = 8640 
   NY = 4320
   NT = 232000000

   hostname = socket.gethostname()
   input_dir = ''
   if 'discover' in hostname:
      input_dir = '/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/'
   else:
      input_dir = '/nobackup/gmao_SIteam/ModelData/l_data/LandBCs_files_for_mkCatchParam/V001/' 

   maskfile = ''

   if orslvs in['O1','T2','T3','T4','T1MOM6','T2MOM6','T4MOM6']:
      maskfile = 'shared/mask/GEOS5_10arcsec_mask_freshwater-lakes.nc'
      if lbcsv in ['F25', 'GM4', 'ICA']:
         maskfile = 'global.cat_id.catch.DL'

   if orslvs in['O2','O3','CS']:
      maskfile = 'shared/mask/EOS5_10arcsec_mask.nc'
      if lbcsv in ['F25', 'GM4', 'ICA']:
         maskfile = 'global.cat_id.catch.GreatLakesCaspian_Updated.DL'

   if grid_type in ['EASEv1', 'EASEv2']:
      maskfile = 'shared/mask/EOS5_10arcsec_mask.nc'

   imo['O1'] = 360
   jmo['O1'] = 180

   imo['O2'] = 1440
   jmo['O2'] = 720

   imo['O3'] = 2880
   jmo['O3'] = 1440
   
   imo['T2'] = 360
   jmo['T2'] = 200
   
   imo['T3'] = 720
   jmo['T3'] = 410
   
   imo['T4'] = 1440
   jmo['T4'] = 1080
   
   imo['T1MOM6'] = 72
   jmo['T1MOM6'] = 36
   
   imo['T2MOM6'] = 360
   jmo['T2MOM6'] = 210
   
   imo['T4MOM6'] = 1440
   jmo['T4MOM6'] = 1080
   
   im['b'] = 144
   jm['b'] = 91
   
   im['c'] = 288
   jm['c'] = 181
   
   im['d'] = 576
   jm['d'] = 361
   
   im['e'] = 1152
   jm['e'] = 721
   
   cubed = [12,24,48, 90, 180,360,720, 768,1000,1152, 1440,1536,2880,3072,5760]
   for i in cubed:
      key = 'c'+str(i)
      im[key] = i
      jm[key] = i*6
   
   if grid_type == 'EASEv2' :
      im['M01'] = 34704 
      jm['M01'] = 14616
   
      im['M03'] = 11568
      jm['M03'] = 4872
      
      im['M09'] = 3856
      jm['M09'] = 1624

      im['M25'] = 1388  
      jm['M25'] = 584
      
      im['M36'] = 964
      jm['M36'] = 406
      
   if grid_type == 'EASEv1' :
      im['M01'] = 34668 
      jm['M01'] = 14688
   
      im['M03'] = 11556
      jm['M03'] = 4896
      
      im['M09'] = 3852
      jm['M09'] = 1632

      im['M25'] = 1383  
      jm['M25'] = 586

      im['M36'] = 963  
      jm['M36'] = 408

   if resolution in ['c768','c1440'] : 
     NX = 17280
     NY = 8640
   if resolution == 'c2800': 
     NX = 21600
     NY = 10800
   if resolution in ['c1536', 'c3072','c5760'] : 
     NX = 43200
     NY = 21600

   if 'GEOS5_10arcsec_mask' in maskfile:
      NX = 43200
      NY = 21600

   config = {}
   config['grid_type'] = grid_type
   config['lbcsv']     = lbcsv
   config['resolution']= resolution
   config['orslvs']    = orslvs
   config ['im']  = im[resolution]
   config ['jm']  = jm[resolution]
   config ['imo'] = imo[orslvs]
   config ['jmo'] = jmo[orslvs]
   config ['NX']  = NX
   config ['NY']  = NY
   config ['NT']  = NT
   config ['MASKFILE']  = maskfile
   user = os.getlogin()
   config ['expdir'] = '/discover/nobackup/'+user+'/BCS_PACKAGE/'+lbcsv+'/'
   now   = datetime.now()
   config ['outdir'] =now.strftime("%Y%m%d%H%M%S")
   config ['inputdir'] = input_dir
   config ['NCPUS'] = 20

   return config

def ask_questions():

   questions = [

        {
            "type": "select",
            "name": "bcs_version",
            "message": "Select land BCS version: \n \
    *BCs produced by this code will differ from BCs in archived directories!!! \n \
    These differences are caused by compiler changes, code improvements and bug \n \
    fixes that were implemented since the archived BCs in the above-mentioned \n \
    directories were originally created.  The impact of these differences on \n \
    science is insignificant, and the parameter files produced by current \n \
    code are scientifically equivalent to the corresponding archived BCs. \n",
            "choices": [ \
                  "F25 : Fortuna-2_5   (archived*: n/a)", \
   "GM4 : Ganymed-4_0   (archived*: /discover/nobackup/ltakacs/bcs/Ganymed-4_0/)", \
   "ICA : Icarus        (archived*: /discover/nobackup/ltakacs/bcs/Icarus/)", \
   "NL3 : Icarus-NLv3   (archived*: /discover/nobackup/ltakacs/bcs/Icarus-NLv3/)", \
   "NL4 : NLv4 [SMAPL4] (archived*: /discover/nobackup/projects/gmao/smap/bcs_NLv4/NLv4/)", \
   "NL5 : NLv5 [SMAPL4] (archived*: /discover/nobackup/projects/gmao/smap/SMAP_L4/L4_SM/bcs/CLSM_params/Icarus-NLv5_EASE/)", \
   "DEV : Development version"],
            "default": "NL3 : Icarus-NLv3   (archived*: /discover/nobackup/ltakacs/bcs/Icarus-NLv3/)",
        },

       {
            "type": "select",
            "name": "grid_type",
            "message": "Select grid type: \n ",
            "choices": ["Lat-Lon", "Cubed-Sphere", "EASEv2", "EASEv1"],
            "default": "Cubed-Sphere",
        },

       {
            "type": "select",
            "name": "resolution",
            "message": "Select lat-lon resolution: \n ",
            "choices": ["b -- 2   deg", "c -- 1  deg", "d -- 1/2 deg","e --  1/4 deg"],
            "default": "d -- 1/2 deg",
            "when": lambda x: x['grid_type'] == "Lat-Lon",
        },

       {
            "type": "select",
            "name": "resolution",
            "message": "Select cubed-sphere resolution: \n ",
            "choices": [ \
                 "c12   -- 8    deg", \
                 "c24   -- 4    deg", \
                 "c48   -- 2    deg", \
                 "c90   -- 1    deg", \
                 "c180  -- 1/2  deg ( 56   km)", \
                 "c360  -- 1/4  deg ( 28   km)", \
                 "c720  -- 1/8  deg ( 14   km)", \
                 "c768  -- 1/10 deg ( 12   km)", \
                 "c1000 -- 1/10 deg ( 10   km)", \
                 "c1152 -- 1/10 deg (  8   km)", \
                 "c1440 -- 1/16 deg (  7   km)", \
                 "c1536 -- 1/16 deg (  7   km)", \
                 "c2880 -- 1/32 deg (  3   km)", \
                 "c3072 -- 1/32 deg (  3   km)", \
                 "c5760 -- 1/64 deg (  1.5 km)"],
            "default": "c360  -- 1/4  deg ( 28   km)",
            "when": lambda x: x['grid_type'] == "Cubed-Sphere",
        },

       {
            "type": "select",
            "name": "resolution",
            "message": "Select EASE grid resolution: \n ",
            "choices": [ \
                 "M01  --  1km", \
                 "M03  --  3km", \
                 "M09  --  9km", \
                 "M25  -- 25km", \
                 "M36  -- 36km"],
            "default": "M09  --  9km",
            "when": lambda x: x['grid_type'] == "EASEv2" or x['grid_type'] == "EASEv1",
        },

       {
            "type": "select",
            "name": "orslvs",
            "message": "Select ocean resolution: \n ",
            "choices": [ \
                 "O1     --  Low-Resolution  Reynolds 1   deg (Lon/Lat Data-Ocean:    360x180 )", \
                 "O2     --  Med-Resolution  Reynolds 1/4 deg (Lon/Lat Data-Ocean:   1440x720 )", \
                 "O3     --  High-Resolution    OSTIA 1/8 deg (Lon/Lat Data-Ocean:   2880x1440)", \
                 "T2     --  Med-Resolution  Tripolar 1   deg (MOM-Tripolar-Ocean:    360x200 )", \
                 "T3     --  High-Resolution Tripolar 1/2 deg (MOM-Tripolar-Ocean:    720x410 )", \
                 "T4     --  High-Resolution Tripolar 1/4 deg (MOM-Tripolar-Ocean:   1440x1080)", \
                 "T1MOM6 --  Low-Resolution  Tripolar 5   deg (MOM6-Tripolar-Ocean:    72x36  )", \
                 "T2MOM6 --  Med-Resolution  Tripolar 1   deg (MOM6-Tripolar-Ocean:   360x210 )", \
                 "T4MOM6 --  High-Resolution Tripolar 1/4 deg (MOM6-Tripolar-Ocean:  1440x1080)", \
                 "CS     --  Cubed-Sphere Ocean               (Cubed-Sphere Data-Ocean        )"],
            "when": lambda x: x['grid_type'] == "Lat-Lon" or x['grid_type'] == "Cubed-Sphere",
        },
       

   ]
   answers_ = questionary.prompt(questions)
   answers  = {}
   for key, value in answers_.items():
      answers[key] = value.split()[0]
   if 'EASE' in answers['grid_type'] : answers['orslvs'] = 'O1'

   return answers

def print_config( config, indent = 0 ):
   for k, v in config.items():
     if isinstance(v, dict):
        print("   " * indent, f"{k}:")
        print_config(v, indent+1)
     else:
        print("   " * indent, f"{k}: {v}")

if __name__ == "__main__":

   answers = ask_questions()
   for key, value in answers.items():
      print(key, value)

   config = get_config_from_answers(answers)
   print('\n print config file:\n')
   print_config(config)   
