import os
import sys
import numpy as np
import netCDF4 as nc
import multiprocessing as mp

# ==============================================================================
# Configuration
# ==============================================================================
# Original 500m 4-day LAI data directory
IN_DIR = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/lai/v5/"
# Target directory
TPL_DIR_BASE = "upscaled/"
# Number of CPU cores for parallel processing
NUM_CORES = 18

# ==============================================================================
# Processing Logic
# ==============================================================================
def process_single_tile(in_file):
    base_name = os.path.basename(in_file)
    # Extract tile ID, e.g., "H01V05"
    tile_id = base_name.split('.')[1] 

    # Check if source file exists before processing
    if not os.path.exists(in_file):
        return

    print(f"Worker started processing tile: {tile_id}", flush=True)
    
    try:
        # Open source file in read mode
        with nc.Dataset(in_file, 'r') as src:
            src_var = src.variables['Lai_500m']
            
            # Disable auto-scaling to manage memory and math manually
            src_var.set_auto_maskandscale(False)
            
            # Loop through each year from 2003 to 2025
            for year_idx, year in enumerate(range(2003, 2026)):
                tpl_file = os.path.join(TPL_DIR_BASE, str(year), f"MODIS_lai_clim.{tile_id}.nc")
                
                if not os.path.exists(tpl_file):
                    print(f"Template missing for {year} {tile_id}, skipping...", flush=True)
                    continue
                
                # Each year has 92 time steps (4-day interval in the source data)
                t_start = year_idx * 92
                t_end = t_start + 92
                
                # 1. Read raw source data for one specific year (shape: 92, 2400, 2400)
                raw_src = src_var[t_start:t_end, :, :]
                
                # 2. Apply original scaling factor (0.1) and create masked array (FillValue = 255)
                mask = (raw_src == 255)
                phys_src = np.ma.array(raw_src.astype(np.float32) * 0.1, mask=mask)
                
                # 3. Spatiotemporal averaging (downscale to 46 steps, 1200x1200)
                # Reshape to create 2x2x2 windows: (t_out, 2, lat_out, 2, lon_out, 2)
                reshaped = phys_src.reshape(46, 2, 1200, 2, 1200, 2)
                chunk_mean = reshaped.mean(axis=(1, 3, 5))
                
                # 4. Convert physical LAI to template's raw short format (ScaleFactor = 0.01)
                # Formula: Raw_target = Physical * 100
                new_raw_data = np.ma.round(chunk_mean * 100.0).astype(np.int16)
                
                # Replace our calculated missing values with template's UNDEF (-9999)
                new_raw_data_filled = new_raw_data.filled(-9999)
                
                # 5. Overwrite the template (Mode 'r+' allows in-place modification)
                with nc.Dataset(tpl_file, 'r+') as tpl:
                    tpl_var = tpl.variables['LAI']
                    
                    # Access raw data to correctly read and write the -9999 value
                    tpl_var.set_auto_maskandscale(False) 
                    
                    # Read the original raw data from the template
                    orig_tpl_data = tpl_var[:]
                    
                    # Identify where the template originally had -9999 (water bodies/invalid areas)
                    keep_mask = (orig_tpl_data == -9999)
                    
                    # Restore -9999 in those exact grid cells, overriding the newly processed data
                    new_raw_data_filled[keep_mask] = -9999
                    
                    # Write the final assembled array back into the template
                    tpl_var[:] = new_raw_data_filled
                    
        print(f"Successfully updated {tile_id} across all years (2003-2025)", flush=True)
        
    except Exception as e:
        print(f"Error processing {base_name}: {e}", flush=True)

# ==============================================================================
# Main Execution
# ==============================================================================
if __name__ == '__main__':
    # Ensure the H index is passed as a command-line argument (via Slurm Array)
    if len(sys.argv) < 2:
        print("Usage: python upscale_global_array.py <H_INDEX>", flush=True)
        sys.exit(1)
        
    h_idx = int(sys.argv[1])
    h_str = f"H{h_idx:02d}"
    
    # Generate file list for all V tiles (V01 to V18) for the specific H tile
    file_list = []
    for v in range(1, 19):
        tile_name = f"{h_str}V{v:02d}"
        file_path = os.path.join(IN_DIR, f"MCD15A3H_lai_2003-2025.{tile_name}.nc")
        file_list.append(file_path)
        
    print(f"Node assigned to {h_str}. Submitting {len(file_list)} tasks to multiprocessing pool...", flush=True)
    
    # Process the 18 'V' tiles in parallel
    with mp.Pool(processes=NUM_CORES) as pool:
        pool.map(process_single_tile, file_list)
        
    print(f"All updates for {h_str} completed successfully!", flush=True)
