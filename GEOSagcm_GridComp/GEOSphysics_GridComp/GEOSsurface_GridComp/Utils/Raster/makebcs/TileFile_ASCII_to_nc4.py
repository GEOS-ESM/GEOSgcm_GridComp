#!/usr/bin/env python3
#
# source install/bin/g5_modules
#
# script to create nc4-formatted tile files from ASCII files for *existing* bcs;
# *not* used by make_bcs

import os
import glob
import subprocess

if __name__ == "__main__":

   BC_base = "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/"
   BC_versions = ["NL3", "GM4"]
   PWD = os.getcwd()
   print(PWD)
   
   for version in BC_versions:
      new_dir = BC_base+version
      if (not os.path.isdir(new_dir)):
         continue
      os.chdir(new_dir ) 
      tile_files = glob.glob("geometry/*/*-Pfafstetter.til")
      for til in tile_files:
         if not os.path.exists(til) :
            print(os.path.abspath(til) + " does not exists")
            continue
         nc4 = til.replace('-Pfafstetter.til','-Pfafstetter.nc4')
         if os.path.exists(nc4) :
            continue
         catch_file = til.replace('geometry', 'land')
         fname = os.path.basename(til)
         catch_file = catch_file.replace(fname, 'clsm/catchment.def')
         result = sp.run["./TileFile_ASCII_to_nc4.x", til, catch_file, capture_output=True, text=True)
         print(result.stdout)

   os.chdir(PWD)
