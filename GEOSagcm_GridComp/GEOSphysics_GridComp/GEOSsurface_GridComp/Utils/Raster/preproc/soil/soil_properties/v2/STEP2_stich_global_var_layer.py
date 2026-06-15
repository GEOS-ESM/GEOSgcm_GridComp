#!/usr/bin/env python3
import os
import re
import numpy as np
import xarray as xr

DX = 1.0 / 120.0   # 0.008333333333... degrees (30 arc-sec)
LON_SIZE = 43200   # 360 / DX
LAT_SIZE = 21600   # 180 / DX
FILL = -9999.0

# Parse lon/lat bounds from filename
PAT = re.compile(r"Lon_(-?\d+)_(-?\d+)_Lat_(-?\d+)_(-?\d+)_([A-Z0-9]+)\.nc$")

def lon_to_col(lon):
    # lon in [-180, 180), map to [0, LON_SIZE)
    return int(round((lon + 180.0) / DX))

def lat_to_row(lat):
    # lat in [-90, 90), map to [0, LAT_SIZE)
    return int(round((lat + 90.0) / DX))

def mosaic_global(variable, layer, tiles_dir, out_nc4):
    # Global grid centers (optional, but nice to write)
    lon = np.linspace(-180.0 + DX/2, 180.0 - DX/2, LON_SIZE, dtype=np.float32)
    lat = np.linspace(-90.0 + DX/2,  90.0 - DX/2,  LAT_SIZE, dtype=np.float32)

    global_data = np.full((LAT_SIZE, LON_SIZE), np.nan, dtype=np.float32)

    tiles = sorted(f for f in os.listdir(tiles_dir) if f.endswith(f"_{layer}.nc"))
    if not tiles:
        raise FileNotFoundError(f"No tiles for layer {layer} in {tiles_dir}")

    for fn in tiles:
        m = PAT.search(fn)
        if not m:
            # skip unexpected filenames
            continue

        lon0, lon1, lat0, lat1, lyr = m.groups()
        lon0 = float(lon0); lon1 = float(lon1)
        lat0 = float(lat0); lat1 = float(lat1)

        # Compute target slice in global array
        # Example tile: Lon_-180_-170 Lat_70_80 => col0..col1, row0..row1
        col0 = lon_to_col(lon0)
        col1 = lon_to_col(lon1)
        row0 = lat_to_row(lat0)
        row1 = lat_to_row(lat1)

        path = os.path.join(tiles_dir, fn)
        with xr.open_dataset(path) as ds:
            if variable not in ds:
                continue

            # Ensure tile data is (lat,lon)
            tile = ds[variable].transpose("lat", "lon").astype("float32").values

            # Safety: expected tile size
            if tile.shape != (1200, 1200):
                raise ValueError(f"Unexpected tile shape {tile.shape} in {fn}")

            global_data[row0:row1, col0:col1] = tile

    # Replace NaN with fill
    global_filled = np.where(np.isnan(global_data), FILL, global_data).astype(np.float32)

    out = xr.Dataset(
        {variable: (("lat", "lon"), global_filled)},
        coords={"lat": lat, "lon": lon},
        attrs={"description": f"Global {variable} mosaicked from HWSD2 tiles", "layer": layer}
    )

    # Set fill value attributes (helps CDO/NCO)
    out[variable].attrs["_FillValue"] = np.float32(FILL)
    out[variable].attrs["missing_value"] = np.float32(FILL)

    out.to_netcdf(out_nc4, format="NETCDF4")
    print(f"Saved: {out_nc4}")

if __name__ == "__main__":
    tiles_dir = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/output_STEP1/"
#Sand
#Silt
#Clay
#Coarse
#SOC
#Organic_Matter
#Share
#Porosity
#Bulk_Density
#Ref_Bulk_Density
#PH
#Sequence_Used    
    layer = "D2"
    variable = "Organic_Matter"  #"Sand"
    out_nc4 = f"Global_{variable}_{layer}.nc4"

    mosaic_global(variable, layer, tiles_dir, out_nc4)

