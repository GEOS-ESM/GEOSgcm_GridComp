#!/usr/bin/env python3
#
# source install/bin/g5_modules
#

import sys
import argparse
import textwrap
import ruamel.yaml
from make_bcs_questionary import *
from make_bcs_ease import *
from make_bcs_latlon import *
from make_bcs_cube import *

# Define the argument parser
def parse_args():

    program_description = textwrap.dedent(f'''
      Usage: 
                                                                              
      Boundary Conditions (BCs) Package:                                      
         Creates surface tile and other model parameter input files           
         (.til [tiles], .rst [raster], and land parameters) for               
         combinations of atmospheric resolution, ocean resolution,            
         and land BCs version.                                                
                                                                              
         STEP 1: Build the model.      (GCM or GEOSldas)                      
         STEP 2: cd [install-path]/bin                                        
         STEP 3: source g5_modules     (for bash or zsh use g5_modules.[z]sh) 
         STEP 4: ./make_bcs.py                                                   
                 Answer the following interactive questions:                  
                 a) Select Land BCs version.                                  
                 b) Select atmospheric resolution(s).                         
                 c) Select ocean resolution(s).                               
                      (Not relevant for land-only EASE-grid BCs.)             
                 d) Enter BCs output directory.                               
                 e) Enter sponsor code for computing account.                 
                                                                              
         To skip the generation of land parameter files (ie, mkCatchParam.x), 
         use:  ./make_bcsi.pc -noland                                             
         This option saves time when additional bcs are created that have     
         the exact same land parameters as an existing set of bcs because     
         the only difference between the two sets of bcs is the [non-tripolar]
         ocean resolution.     ''')
 
    parser = argparse.ArgumentParser(description='make bcs',epilog=program_description,formatter_class=argparse.RawDescriptionHelpFormatter)
    # define a mutually exclusive group of arguments
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-c', '--config_file',    help='YAML config file')

    # Parse using parse_known_args so we can pass the rest to the remap scripts
    # If config_file is used, then extra_args will be empty
    # If flattened_yaml is used, then extra_args will be populated
    args, extra_args = parser.parse_known_args()
    return args, extra_args

def main():

  question_flag = False

  # Parse the command line arguments from parse_args() capturing the arguments and the rest
  command_line_args, extra_args = parse_args()
  print(f'command_line_args: {command_line_args}')
  config_yaml = command_line_args.config_file

  configs = []
  if config_yaml:
      config = yaml_to_config(config_yaml)
      configs = [config]
  else:
      answers = ask_questions()
      configs = get_configs_from_answers(answers)
  for config in configs :
      if 'EASE' in config['grid_type']:
         make_bcs_ease(config)    
      if 'Lat-Lon' in config['grid_type']:
         make_bcs_latlon(config)    
      if 'Cubed-Sphere' in config['grid_type'] or 'Stretched_CS' in config['grid_type']:
         make_bcs_cube(config)    

if __name__ == '__main__' :
  #exit("The python version of make_bcs is not yet ready for general use.  Until further notice, please use csh script make_bcs")
  main()

