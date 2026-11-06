import xarray as xr
import numpy as np
import pandas as pd
import glob
import os
import argparse
from scipy.signal import savgol_filter
import sys

# ==============================================================================
# 1. Core Algorithm & Preprocess Functions
# ==============================================================================

def timesat_upper_envelope_sg_nd(da, window_size=15, poly_order=2, iterations=3):
    """
    Applies an upper-envelope Savitzky-Golay filter to time-series data.
    This mimics the TIMESAT algorithm (Jönsson and Eklundh 2004), 
    which assumes negative biases in NDVI/LAI (due to clouds/aerosols) 
    and iteratively fits the upper envelope to recover the true vegetation peak.

    Jönsson, P., & Eklundh, L. (2004). TIMESAT—a program for analyzing 
    time-series of satellite sensor data. Computers & geosciences, 30(8), 833-845.
    """
    y = da.values 
    y_fit = np.copy(y)
    
    # Iteratively apply SG filter and force the fitted curve towards the upper envelope
    for _ in range(iterations):
        y_smooth = savgol_filter(y_fit, window_size, poly_order, axis=0)
        y_fit = np.maximum(y, y_smooth)
        
    # Final smoothing step
    y_final = savgol_filter(y_fit, window_size, poly_order, axis=0)
    
    # Clip physically impossible values (LAI bounded between 0.0 and 8.0 for this context)
    y_clipped = np.clip(y_final, 0.0, 8.0).astype(np.float32)
    
    return xr.DataArray(y_clipped, dims=da.dims, coords=da.coords)

def preprocess_modis(ds):
    """
    Extracts the date from the MODIS HDF/NetCDF filename and adds a 'time' dimension.
    Expected filename format: MCD15A3H.061_YYYYDDD_*.nc4 (where DDD is day of year).
    """
    fname = ds.encoding["source"]
    date_str = os.path.basename(fname).split('_')[1].split('.')[0]
    time_index = pd.to_datetime([date_str], format='%Y%j')
    ds = ds.expand_dims(time=time_index)
    return ds

def process_modis_chunk(ds):
    """
    Parses MODIS Quality Control (QC) flags via bitwise operations and applies masking.
    Retains only 'perfect' pixels or conditionally fills snow-covered pixels.
    """
    # Scale LAI (MODIS scale factor is 0.1) and exclude fill values (255)
    lai = ds['Lai_500m'].where(ds['Lai_500m'] != 255) * 0.1
    
    qc = ds['FparLai_QC'].where(ds['FparLai_QC'] != 255)
    exqc = ds['FparExtra_QC'].where(ds['FparExtra_QC'] != 255)
    
    qc_int = qc.fillna(255).astype(int)
    exqc_int = exqc.fillna(255).astype(int)
    
    # Bitwise decoding of MODIS QC flags
    is_good_main = (qc_int & 1) == 0                   # Bit 0: Good quality
    is_snow = ((exqc_int >> 2) & 1) == 1               # Bit 2: Snow/Ice
    is_internal_cloud = ((exqc_int >> 5) & 1) == 1     # Bit 5: Internal Cloud
    is_shadow = ((exqc_int >> 6) & 1) == 1             # Bit 6: Cloud Shadow
    
    # A pixel is 'perfect' if it has good quality and no clouds, snow, or shadows
    is_perfect = is_good_main & (~is_internal_cloud) & (~is_snow) & (~is_shadow)
    
    # Calculate a baseline for filling snow pixels (minimum valid LAI over time)
    perfect_lai = lai.where(is_perfect & (lai > 0.05))
    p_fillbase = perfect_lai.min(dim='time')
    p_fillbase = p_fillbase.fillna(0.0) # Ensure water bodies get 0.0 base, not 0.1
    
    # Apply baseline to snow pixels, keep original LAI otherwise
    lai_filtered = xr.where(is_snow, p_fillbase, lai)
    
    # Mask out bad pixels (not perfect and not snow)
    bad_mask = (~is_perfect) & (~is_snow)
    lai_filtered = lai_filtered.where(~bad_mask)
    
    # Cap anomalously high LAI values
    lai_filtered = lai_filtered.where(lai_filtered <= 25.0)
    
    return lai_filtered

# ==============================================================================
# 2. Main Execution Block with Dynamic Arguments
# ==============================================================================
if __name__ == "__main__":
    from dask.diagnostics import ProgressBar

    # Parse the task ID and base coordinates passed by the Bash script
    parser = argparse.ArgumentParser(description="Process a 2.5x2.5 degree MODIS tile.")
    parser.add_argument('--task_id', type=int, required=True, help="Slurm Array Task ID (0-15)")
    parser.add_argument('--base_lat', type=float, required=True, help="Base latitude (bottom-left corner) of the 10x10 tile")
    parser.add_argument('--base_lon', type=float, required=True, help="Base longitude (bottom-left corner) of the 10x10 tile")
    parser.add_argument('--data_dir', type=str, required=True, help="Path to the MODIS MCD15A3H.061 raw data directory")
    args = parser.parse_args()
    
    # --------------------------------------------------------------------------
    # Master 10x10 Grid Definition & Sub-tile Mapping
    # --------------------------------------------------------------------------
    base_lat = args.base_lat
    base_lon = args.base_lon
    step = 2.5
    
    # Calculate row (lat) and col (lon) index based on task_id (0-15)
    # Maps a 1D Slurm Array ID to a 4x4 2D spatial grid: row = id // 4, col = id % 4
    row = args.task_id // 4
    col = args.task_id % 4
    
    lat_min = base_lat + row * step
    lat_max = lat_min + step
    lon_min = base_lon + col * step
    lon_max = lon_min + step
    
    out_file = f"temporary/MCD15A3H_2003_2025_lat_{lat_min}_lon_{lon_min}.nc"
    
    print(f"--- Task ID {args.task_id} Started ---", flush=True)
    print(f"Master Tile Base: Lat={base_lat}, Lon={base_lon}", flush=True)
    print(f"Targeting Sub-tile: Lat [{lat_min} to {lat_max}], Lon [{lon_min} to {lon_max}]", flush=True)
    print(f"Output will be saved to: {out_file}", flush=True)
    
    # --------------------------------------------------------------------------
    # Data Loading
    # --------------------------------------------------------------------------
    data_dir = args.data_dir
    
    files = []
    for year in range(2002, 2027):
        file_pattern = os.path.join(data_dir, str(year), f"MCD15A3H.061_{year}*.nc4")
        files.extend(sorted(glob.glob(file_pattern)))
        
    ds = xr.open_mfdataset(
        files, 
        preprocess=preprocess_modis,
        combine='by_coords', 
        chunks="auto",
        parallel=True
    )
    
    # Slice target region (handle both ascending and descending latitude coords safely)
    lat_arr = ds['lat'].values
    if lat_arr[0] > lat_arr[-1]:
        ds_region = ds.sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max))
    else:
        ds_region = ds.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
        
    # ==========================================================================
    # Dask Workflow: QC -> Interpolation -> Smoothing
    # ==========================================================================
    print("Applying QC Filtering...", flush=True)
    lai_qc = process_modis_chunk(ds_region)
    
    print("Interpolating missing values...", flush=True)
    # Rechunk across space to ensure full time-series is available per pixel for interpolation
    lai_qc = lai_qc.chunk({'time': -1, 'lat': 'auto', 'lon': 'auto'})
    
    # Temporal interpolation and edge filling
    lai_filled = lai_qc.interpolate_na(dim='time', method='linear')
    lai_filled = lai_filled.bfill(dim='time').ffill(dim='time').fillna(0.0) 
    
    print("Applying S-G Smoothing...", flush=True)
    # Map the TIMESAT algorithm block-by-block via Dask
    lai_smoothed = lai_filled.map_blocks(
        timesat_upper_envelope_sg_nd,
        kwargs={'window_size': 15, 'poly_order': 2, 'iterations': 3},
        template=lai_filled
    )

    core_time_slice = slice('2003-01-01', '2025-12-31')    
    # Truncate to the core study period
    lai_smoothed_core = lai_smoothed.sel(time=core_time_slice)
    
    lai_smoothed_core.name = 'Lai_500m'
    lai_smoothed_core.attrs['units'] = "m^2/m^2"
    lai_smoothed_core.attrs['long_name'] = f"QC-Filtered Upper-Envelope Smoothed LAI ({lat_min} to {lat_max}, {lon_min} to {lon_max})"
    
    # ==========================================================================
    # I/O Execution & Teardown Hang Prevention
    # ==========================================================================
    print("Triggering lazy computation into RAM...", flush=True)
    # .load() forces Dask to execute the graph into memory BEFORE I/O.
    # This prevents Lustre filesystem lock conflicts during concurrent writes.
    lai_smoothed_core = lai_smoothed_core.load() 
    
    print(f"Writing physical data to {out_file}...", flush=True)
    lai_smoothed_core.to_dataset().to_netcdf(out_file)

    print(f"Task {args.task_id} successfully saved! Forcing immediate exit to prevent Dask hang...", flush=True)
    
    # Force OS-level exit to prevent Dask thread-pool from hanging during garbage collection
    import os
    os._exit(0)
