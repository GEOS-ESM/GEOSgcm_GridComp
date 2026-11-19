import xarray as xr
import os
import sys
import argparse
from dask.diagnostics import ProgressBar

if __name__ == "__main__":
    # ==========================================================================
    # 1. Parse Arguments (base_lat, base_lon)
    # ==========================================================================
    parser = argparse.ArgumentParser(description="Stitch 16 sub-tiles into a 10x10 master tile.")
    parser.add_argument('--base_lat', type=float, required=True, help="Base latitude (bottom-left corner)")
    parser.add_argument('--base_lon', type=float, required=True, help="Base longitude (bottom-left corner)")
    args = parser.parse_args()

    base_lat = args.base_lat
    base_lon = args.base_lon
    step = 2.5
    
    print(f"Initiating stitching for Master Tile starting at Lat: {base_lat}, Lon: {base_lon}", flush=True)
    
    # ==========================================================================
    # 2. Strict File Validation & Size Checking
    # ==========================================================================
    expected_files = []
    
    # Generate the expected 16 filenames based on the 4x4 grid math
    for row in range(4):
        for col in range(4):
            lat_min = base_lat + row * step
            lon_min = base_lon + col * step
            expected_files.append(f"temporary/MCD15A3H_2003_2025_lat_{lat_min}_lon_{lon_min}.nc")
            
    # Set the minimum size threshold to 2.8 GB (2.8 * 1024^3 bytes)
    MIN_SIZE_BYTES = 2.8 * 1024 * 1024 * 1024 
    
    print(f"Verifying existence and size (>= 2.8 GB) of {len(expected_files)} required sub-tiles...", flush=True)
    
    for f in expected_files:
        # Check 1: Verify file existence
        if not os.path.exists(f):
            print(f"\n FATAL ERROR: Required file is missing: {f}")
            print("Aborting stitch process to prevent incomplete master tile.")
            sys.exit(1)
            
        # Check 2: Verify file size meets the threshold
        f_size_bytes = os.path.getsize(f)
        f_size_gb = f_size_bytes / (1024 ** 3)
        if f_size_gb < 2.8:
            print(f"\n FATAL ERROR: File too small! {f}")
            print(f"Expected >= 2.8 GB, but found {f_size_gb:.2f} GB.")
            print("This usually indicates a killed or incomplete Slurm task. Aborting.")
            sys.exit(1)
            
    print("All 16 sub-tiles exist and pass the 2.8GB size check!", flush=True)
    
    # ==========================================================================
    # 3. Memory Loading & Stitching
    # ==========================================================================
    print("Loading data entirely into RAM (Bypassing Dask Write Locks)...", flush=True)
    ds_merged = xr.open_mfdataset(
        expected_files, 
        combine='by_coords'
    )
    
    ds_merged.load()
    
    # ==========================================================================
    # Calculate Grid Indices (H01-H36, V01-V18) and Format Filename
    # ==========================================================================
    h_idx = int((base_lon + 180.0) / 10.0) + 1
    v_idx = int((base_lat + 90.0) / 10.0) + 1
    
    # Format with leading zeros (e.g., 1 -> 01)
    out_file = f"output/MCD15A3H_lai_2003-2025.H{h_idx:02d}V{v_idx:02d}.nc"
    
    ds_merged.attrs['Description'] = f"Merged 10x10 degree MODIS LAI smoothed tile (Tile: H{h_idx:02d}V{v_idx:02d}, Base: Lat {base_lat}, Lon {base_lon})."
    
    # ==========================================================================
    # 4. Extreme Compression Encoding (ubyte + zlib)
    # ==========================================================================
    encoding_dict = {
        'Lai_500m': {
            'dtype': 'uint8',            # Maps exactly to netCDF 'ubyte'
            'scale_factor': 0.1,         
            'add_offset': 0.0,           
            '_FillValue': 255,           
            'zlib': True,                
            'complevel': 1               # Level 1 for ultra-fast run-length compression
        },
        'time': {
            'dtype': 'int32'
        },
        'lat': {'dtype': 'float32'},
        'lon': {'dtype': 'float32'}
    }
    
    print(f"Saving highly compressed (ubyte + zlib) dataset to {out_file}...", flush=True)
    
    #with ProgressBar():
    ds_merged.to_netcdf(out_file, encoding=encoding_dict)
        
    print(f"Stitching Complete! Master tile saved as: {out_file}", flush=True)
