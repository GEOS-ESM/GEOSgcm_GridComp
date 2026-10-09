#!/usr/bin/env python3
import os
import glob
import subprocess
import numpy as np
from netCDF4 import Dataset
import multiprocessing as mp

# Define directory paths
IN_DIR = "upscaled"
OUT_DIR = "MODIS_8-DayClim"

YEAR_START = 2003
YEAR_END = 2025

def setup_directories():
    """
    Create the output directory and copy 2025 files to serve as templates.
    """
    os.makedirs(OUT_DIR, exist_ok=True)
    
    template_src = f"{IN_DIR}/2025"
    print(f"Syncing templates from {template_src} to {OUT_DIR}...")
    
    # Use rsync for robust copying, preserving file attributes and avoiding I/O truncation
    cmd = ["rsync", "-a", f"{template_src}/", f"{OUT_DIR}/"]
    subprocess.run(cmd, check=True)
    print("Template files copied successfully.")

def process_tile(filename):
    """
    Calculate the 2003-2025 climatological mean for a single LAI tile.
    """
    out_file = os.path.join(OUT_DIR, filename)
    
    # Initialize accumulators. 
    # int32 safely holds max short (32767) * 23 years without overflow.
    sum_array = np.zeros((46, 1200, 1200), dtype=np.int32)
    count_array = np.zeros((46, 1200, 1200), dtype=np.int16)
    
    for year in range(YEAR_START, YEAR_END + 1):
        in_file = os.path.join(IN_DIR, str(year), filename)
        
        if not os.path.exists(in_file):
            continue
            
        with Dataset(in_file, 'r') as nc:
            # Disable auto-scaling and auto-masking for max speed and raw integer math
            nc.set_auto_mask(False)
            nc.set_auto_scale(False)
            data = nc.variables['LAI'][:]
            
            # Mask valid pixels (exclude -9999)
            valid = data != -9999
            sum_array[valid] += data[valid]
            count_array[valid] += 1
            
    # Calculate climatological mean
    mean_array = np.full((46, 1200, 1200), -9999, dtype=np.int16)
    valid_counts = count_array > 0
    
    # Round the float mean to nearest integer to preserve original ScaleFactor logic
    mean_array[valid_counts] = np.round(
        sum_array[valid_counts] / count_array[valid_counts]
    ).astype(np.int16)
    
    # Write directly back to the template in update mode ('r+')
    with Dataset(out_file, 'r+') as nc_out:
        nc_out.set_auto_mask(False)
        nc_out.set_auto_scale(False)
        nc_out.variables['LAI'][:] = mean_array
        
    return filename

if __name__ == '__main__':
    # 1. Setup folders and clone templates
    setup_directories()
    
    # 2. Get list of all NC files in the output directory
    nc_files = [os.path.basename(f) for f in glob.glob(f"{OUT_DIR}/*.nc")]
    
    print(f"Starting climatology calculation for {len(nc_files)} tiles...")
    
    # 3. Process tiles in parallel
    # 20 processes is a safe and highly efficient number for a Discover compute node
    num_cores = min(20, mp.cpu_count()) 
    with mp.Pool(num_cores) as pool:
        for i, completed_file in enumerate(pool.imap_unordered(process_tile, nc_files), 1):
            print(f"[{i}/{len(nc_files)}] Completed {completed_file}")
            
    print("All climatology calculations completed successfully!")
