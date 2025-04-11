#!/usr/bin/env python3

from netCDF4 import Dataset
from scipy.io import FortranFile
import numpy as np
import sys
import os 
import glob

def get_i_j_pf(til_file):
""" This is just for cubed-sphere text tile."""
   ii   = []
   jj   = []
   pfaf = []
   with open(til_file, 'r') as f:
       for line in f:
          data = line.split()
          if data[0] == '100' :
             ii.append(int(data[4]))
             jj.append(int(data[5]))
             pfaf.append(int(data[8]))
   return ii, jj, pfaf

def create_reduced_nc4(input_file, output_file, indices):
    """
    Create a new nc4 file with with indices of for the input file.
    
    Parameters:
    input_file (str): Path to the input nc4 file
    output_file (str): Path to the output nc4 file
    indices(int array): 
    """
    
    # Open the input file in read mode
    try:
        src = Dataset(input_file, 'r')
        print(f"Successfully opened input file: {input_file}")
    except Exception as e:
        print(f"Error opening input file: {e}")
        return

    # Create the output file in write mode
    try:
        dst = Dataset(output_file, 'w', format='NETCDF4')
        print(f"Creating output file: {output_file}")
    except Exception as e:
        print(f"Error creating output file: {e}")
        src.close()
        return

    # Copy global attributes
    dst.setncatts(src.__dict__)

    # Copy dimensions with reduced size for unlimited/time dimensions
    for name, dimension in src.dimensions.items():
       if name == 'tile':
          dst.createDimension(name, len(indices))
       else:
          dst.createDimension(name, len(dimension))
    # Copy variables with reduced data
    for name, variable in src.variables.items():
        # Create variable with same attributes but possibly reduced dimensions
        dims = variable.dimensions
        new_dims = []
        for dim in dims:
            new_dims.append(dim)
        # Create the variable in the output file
        dst_var = dst.createVariable(
            name, 
            variable.datatype, 
            new_dims, 
        )
        
        # Copy variable attributes
        dst[name].setncatts(variable.__dict__)
        
        # Get the data and reduce it
        if len(dims) == 2 :
          dst_var[:,:] = variable[:,indices]
        if len(dims) == 1 :
          dst_var[:] = variable[indices]
        

    # Close both files
    src.close()
    dst.close()
    print(f"Successfully created reduced file: {output_file}")

def create_reduced_bin(input_file, output_file, indices):
   """
    Create a new ibinary file with with indices of for the input file.
    
    Parameters:
    input_file (str): Path to the input binary file
    output_file (str): Path to the output binaryfile
    indices(int array): 
   """
   nland = len(indices)
   fout = FortranFile(output_file, 'w')
   with FortranFile(input_file, 'r') as f:
     while True :
       try:
         a = f.read_reals(dtype=np.float32)
         b = f.read_reals(dtype=np.float32)
         a[12] = np.float32(nland)
         fout.write_record(a)
         b = b[indices]
         fout.write_record(b)
       except :
         break
   fout.close()
   print(f"Successfully created reduced file: {output_file}")

# Example usage
if __name__ == "__main__":
    # Replace these with your actual file paths
    # input_nc4 = "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v13/land/CF1440x6C/clsm/catchcn_params.nc4"
    # output_nc4 = "output_file_reduced.nc"
    bcs_dir = sys.argv[1]
    bcs_ver = sys.argv[2]
    air_res = sys.argv[3]
    ocn_res = sys.argv[4]
    agname = air_res + '_' + air_res
    orig_tile = bcs_dir + '/' + bcs_ver + '/geometry/' + agname + '/' + agname +'-Pfafstetter.til'
    if not os.path.exists(orig_tile) :
      print( "The original tile file must exist " + orig_tile)
      exit()
    
    mom_tile = 'til/' + air_res + '_' + ocn_res + '-Pfafstetter.til'
    if not os.path.exists(mom_tile) :
      print( "The MOM tile file must exist " + mom_tile)
      exit()
   
    ii1, jj1, pf1 = get_i_j_pf(orig_tile)
    ii2, jj2, pf2 = get_i_j_pf(mom_tile)
    i1 = 0
    indices =[]
    for i2 in range(len(ii2)):
       match =  ii1[i1] == ii2[i2] and jj1[i1] == jj2[i2] and pf1[i1] == pf2[i2]
       while not match :
          i1 = i1 +1
          match =  ii1[i1] == ii2[i2] and jj1[i1] == jj2[i2] and pf1[i1] == pf2[i2]
       indices.append(i1)
       i1 = i1 + 1 
    
    catch_params_file   = bcs_dir + '/' + bcs_ver + '/land/' + air_res + '/clsm/catch_params.nc4'
    catchcn_params_file = bcs_dir + '/' + bcs_ver + '/land/' + air_res + '/clsm/catchcn_params.nc4'
    mom_catch_params_file   = 'clsm/catch_params.nc4'
    mom_catchcn_params_file = 'clsm/catchcn_params.nc4'

    create_reduced_nc4(catch_params_file,   mom_catch_params_file, indices)
    create_reduced_nc4(catchcn_params_file, mom_catchcn_params_file, indices)

    files  = glob.glob(bcs_dir + '/' + bcs_ver + '/land/' + air_res + '/*.da*')
    for file in files:
      fname = os.path.basename(file)
      if 'vegdyn' in fname:  
         create_reduced_nc4( file,'clsm/vegdyn.data', indices)
      else:
         shortname = 'clsm/' + fname.split('_')[0] + '.dat'
         create_reduced_bin(file, shortname, indices)
