import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# -----------------------------
# User settings
# -----------------------------
V3_DIR  = "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/soil/SOIL-DATA/soil_properties/v3/"
#NEW_DIR = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/out_HWSDv2_NGDC_STATSGO_noMASK/"
NEW_DIR = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v4/"
OUTDIR  = "./soil_global_plots"
os.makedirs(OUTDIR, exist_ok=True)

# Variables we care about (Silt is derived)
PLOT_VARS = ["Mask", "Clay0_30", "Sand0_30", "Silt0_30", "OC0_30"] # Silt is derived from Clay+Sand

# Coarse global plotting resolution (degrees)
# 0.25 is nicer but heavier. 0.5 is fast.
GLOBAL_RES_DEG = 0.5

# Tile list (H01..H36, V01..V18)
HS = range(1, 37)
VS = range(1, 19)

# Optional zoom region (set to None to skip)
# Example: Caribbean-ish
ZOOM = None
#ZOOM = dict(name="Caribbean", lon_min=-85, lon_max=-55, lat_min=-25, lat_max=5)
#ZOOM = dict(name="Indonesia", lon_min=90, lon_max=150, lat_min=-15, lat_max=15)
#ZOOM  = dict(name="CONUS_tile_H09V13", lon_min=-100, lon_max=-90, lat_min=30, lat_max=40)
ZOOM   = dict(name="Africa_tile_H20V09", lon_min=10,   lon_max=20,  lat_min=-10, lat_max=0)

# Native zoom sampling step (1 = full, 2 = every other pixel, etc.)
ZOOM_STEP = 1

# -----------------------------
# Helpers
# -----------------------------
def tile_name(h, v):
    return f"SoilProperties_H{h:02d}V{v:02d}.nc"

def get_undef_sf(da):
    undef = da.attrs.get("UNDEF", -9999)
    sf = float(da.attrs.get("ScaleFactor", 1.0))
    return undef, sf

def read_percent(ds, varname):
    """
    Return var in physical percent units (float64), with undef->nan.
    """
    da = ds[varname]
    a = da.values.astype(np.float64)
    undef, sf = get_undef_sf(da)
    out = np.where(a != undef, a * sf, np.nan)
    return out

def read_mask(ds):
    da = ds["Mask"]
    a = da.values.astype(np.float64)
    undef = da.attrs.get("UNDEF", -9999)
    # keep 0/1, undef->nan
    return np.where(a != undef, a, np.nan)

def compute_silt(clay_pct, sand_pct):
    silt = 100.0 - clay_pct - sand_pct
    # clamp physically a bit; keep NaNs
    silt = np.where(np.isfinite(silt), np.clip(silt, 0.0, 100.0), np.nan)
    return silt

def land_mask_from(ds):
    clay = read_percent(ds, "Clay0_30")
    sand = read_percent(ds, "Sand0_30")
    return np.isfinite(clay) | np.isfinite(sand)

def accumulate_to_grid(lat2d, lon2d, val2d, lat_bins, lon_bins, sum_grid, cnt_grid):
    """
    Vectorized binning accumulation: adds val2d into (lat,lon) coarse bins.
    """
    m = np.isfinite(val2d)
    if not np.any(m):
        return

    latv = lat2d[m]
    lonv = lon2d[m]
    vv   = val2d[m]

    # bin indices
    yi = np.digitize(latv, lat_bins) - 1
    xi = np.digitize(lonv, lon_bins) - 1

    ok = (yi >= 0) & (yi < sum_grid.shape[0]) & (xi >= 0) & (xi < sum_grid.shape[1])
    yi = yi[ok]; xi = xi[ok]; vv = vv[ok]
    if yi.size == 0:
        return

    np.add.at(sum_grid, (yi, xi), vv)
    np.add.at(cnt_grid, (yi, xi), 1)

def finalize_mean(sum_grid, cnt_grid):
    mean = np.full_like(sum_grid, np.nan, dtype=np.float64)
    ok = cnt_grid > 0
    mean[ok] = sum_grid[ok] / cnt_grid[ok]
    return mean

def plot_global(field2d, lat_bins, lon_bins, title, out_png, vmin=None, vmax=None, cmap="viridis"):
    plt.figure(figsize=(11, 5))
    extent = [lon_bins[0], lon_bins[-1], lat_bins[0], lat_bins[-1]]
    im = plt.imshow(field2d, origin="lower", extent=extent, vmin=vmin, vmax=vmax, cmap=cmap, aspect="auto")
    plt.colorbar(im, label="percent")
    plt.title(title)
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def plot_global_diff(diff2d, lat_bins, lon_bins, title, out_png, vmax=5.0):
    plt.figure(figsize=(11, 5))
    extent = [lon_bins[0], lon_bins[-1], lat_bins[0], lat_bins[-1]]
    im = plt.imshow(diff2d, origin="lower", extent=extent, vmin=-vmax, vmax=vmax,
                    cmap="RdBu_r", aspect="auto")
    plt.colorbar(im, label="NEW - v3 (percentage points)")
    plt.title(title)
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def build_global_coarse_map(root_dir, varname, res_deg=0.5):
    """
    Build coarse global mean map for a variable from 36x18 tiles.
    varname can be Clay0_30, Sand0_30, OC0_30, or derived Silt0_30.
    """
    lat_bins = np.arange(-90.0, 90.0 + res_deg, res_deg)
    lon_bins = np.arange(-180.0, 180.0 + res_deg, res_deg)

    sum_grid = np.zeros((len(lat_bins)-1, len(lon_bins)-1), dtype=np.float64)
    cnt_grid = np.zeros((len(lat_bins)-1, len(lon_bins)-1), dtype=np.int64)
    max_grid = None
    if varname == "Mask":
        max_grid = np.zeros((len(lat_bins)-1, len(lon_bins)-1), dtype=np.float64)    

    for h in HS:
        for v in VS:
            f = os.path.join(root_dir, tile_name(h, v))
            if not os.path.exists(f):
                continue

            with xr.open_dataset(f, decode_cf=False, mask_and_scale=False) as ds:
                lat1 = ds["latitude"].values.astype(np.float64)
                lon1 = ds["longitude"].values.astype(np.float64)

                # build 2D grids (lat varies by row, lon by col)
                lat2d = np.repeat(lat1[:, None], lon1.size, axis=1)
                lon2d = np.repeat(lon1[None, :], lat1.size, axis=0)

                if varname == "Silt0_30":
                    clay = read_percent(ds, "Clay0_30")
                    sand = read_percent(ds, "Sand0_30")
                    val = compute_silt(clay, sand)
                
                    accumulate_to_grid(lat2d, lon2d, val, lat_bins, lon_bins, sum_grid, cnt_grid)
                
                elif varname == "Mask":
                    val = read_mask(ds)
                    m = np.isfinite(val)
                    if np.any(m):
                        latv = lat2d[m]
                        lonv = lon2d[m]
                        vv   = val[m]
                
                        yi = np.digitize(latv, lat_bins) - 1
                        xi = np.digitize(lonv, lon_bins) - 1
                
                        ok = (yi >= 0) & (yi < max_grid.shape[0]) & (xi >= 0) & (xi < max_grid.shape[1])
                        yi = yi[ok]; xi = xi[ok]; vv = vv[ok]
                
                        np.maximum.at(max_grid, (yi, xi), vv)
                
                else:
                    val = read_percent(ds, varname)
                    accumulate_to_grid(lat2d, lon2d, val, lat_bins, lon_bins, sum_grid, cnt_grid)

    if varname == "Mask":
        return max_grid, lat_bins, lon_bins
    
    mean = finalize_mean(sum_grid, cnt_grid)
    return mean, lat_bins, lon_bins

def build_zoom_native(root_dir, varname, lon_min, lon_max, lat_min, lat_max, step=1):
    """
    Build a native-resolution (tile-res) zoom mosaic for a lon/lat box.
    Stitches only tiles that intersect the region and samples with step.
    Returns (field, lon2d, lat2d) for plotting with pcolormesh/imshow.
    """
    chunks = []

    for h in HS:
        for v in VS:
            f = os.path.join(root_dir, tile_name(h, v))
            if not os.path.exists(f):
                continue

            with xr.open_dataset(f, decode_cf=False, mask_and_scale=False) as ds:
                lat = ds["latitude"].values.astype(np.float64)
                lon = ds["longitude"].values.astype(np.float64)

                # quick tile bbox reject
                if (lat.max() < lat_min) or (lat.min() > lat_max) or (lon.max() < lon_min) or (lon.min() > lon_max):
                    continue

                # index window in this tile
                ilat = np.where((lat >= lat_min) & (lat <= lat_max))[0]
                ilon = np.where((lon >= lon_min) & (lon <= lon_max))[0]
                if ilat.size == 0 or ilon.size == 0:
                    continue

                ilat = ilat[::step]
                ilon = ilon[::step]

                if varname == "Silt0_30":
                    clay = read_percent(ds, "Clay0_30")[np.ix_(ilat, ilon)]
                    sand = read_percent(ds, "Sand0_30")[np.ix_(ilat, ilon)]
                    val = compute_silt(clay, sand)
                elif varname == "Mask":
                    val = read_mask(ds)[np.ix_(ilat, ilon)]
                else:
                    val = read_percent(ds, varname)[np.ix_(ilat, ilon)]

                lat_sub = lat[ilat]
                lon_sub = lon[ilon]

                chunks.append((lat_sub, lon_sub, val))

    if not chunks:
        return None, None, None

    # We’ll mosaic by binning onto a regular lat/lon grid at native sample spacing:
    # infer approximate dx/dy from the first chunk
    lat0, lon0, _ = chunks[0]
    dlat = np.median(np.abs(np.diff(lat0))) if lat0.size > 1 else 0.008333
    dlon = np.median(np.abs(np.diff(lon0))) if lon0.size > 1 else 0.008333

    lat_bins = np.arange(lat_min, lat_max + dlat, dlat)
    lon_bins = np.arange(lon_min, lon_max + dlon, dlon)

    sum_grid = np.zeros((lat_bins.size, lon_bins.size), dtype=np.float64)
    cnt_grid = np.zeros((lat_bins.size, lon_bins.size), dtype=np.int64)

    for lat_sub, lon_sub, val in chunks:
        lat2d = np.repeat(lat_sub[:, None], lon_sub.size, axis=1)
        lon2d = np.repeat(lon_sub[None, :], lat_sub.size, axis=0)
        accumulate_to_grid(lat2d, lon2d, val, lat_bins, lon_bins, sum_grid, cnt_grid)

    mean = finalize_mean(sum_grid, cnt_grid)

    lon2d = np.repeat(lon_bins[None, :], lat_bins.size, axis=0)
    lat2d = np.repeat(lat_bins[:, None], lon_bins.size, axis=1)
    return mean, lon2d, lat2d

def plot_zoom(field, lon2d, lat2d, title, out_png, vmin=None, vmax=None, cmap="viridis"):
    plt.figure(figsize=(9, 6))
    extent = [lon2d.min(), lon2d.max(), lat2d.min(), lat2d.max()]
    im = plt.imshow(field, origin="lower", extent=extent, vmin=vmin, vmax=vmax, cmap=cmap, aspect="auto")
    plt.colorbar(im, label="percent")
    plt.title(title)
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def plot_zoom_diff(diff, lon2d, lat2d, title, out_png, vmax=5.0):
    plt.figure(figsize=(9, 6))
    extent = [lon2d.min(), lon2d.max(), lat2d.min(), lat2d.max()]
    im = plt.imshow(diff, origin="lower", extent=extent, vmin=-vmax, vmax=vmax,
                    cmap="RdBu_r", aspect="auto")
    plt.colorbar(im, label="NEW - v3 (percentage points)")
    plt.title(title)
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

# -----------------------------
# Main
# -----------------------------
def main():
    # Global coarse maps
    for v in PLOT_VARS:
        print("Building global coarse:", v)

        new_map, lat_bins, lon_bins = build_global_coarse_map(
            NEW_DIR, v, res_deg=GLOBAL_RES_DEG
        )
        v3_map, _, _ = build_global_coarse_map(
            V3_DIR, v, res_deg=GLOBAL_RES_DEG
        )

        # ---- variable-specific plotting settings ----
        if v == "Mask":
            vmin, vmax = 0, 1
            cmap = "gray"
            title_new = (
                f"NEW Mask (coarse {GLOBAL_RES_DEG}°)\n"
                "1 = soil properties exist, 0 = no soil"
            )
            title_v3 = (
                f"v3 Mask (coarse {GLOBAL_RES_DEG}°)\n"
                "1 = soil properties exist, 0 = no soil"
            )
        else:
            vmin, vmax = 0, 100
            cmap = "jet"
            title_new = f"NEW {v} (coarse {GLOBAL_RES_DEG}°)"
            title_v3  = f"v3 {v} (coarse {GLOBAL_RES_DEG}°)"

        # ---- Plot absolute maps ----
        plot_global(
            new_map, lat_bins, lon_bins,
            title_new,
            os.path.join(OUTDIR, f"global_NEW_{v}_{GLOBAL_RES_DEG}deg.png"),
            vmin=vmin, vmax=vmax, cmap=cmap
        )

        plot_global(
            v3_map, lat_bins, lon_bins,
            title_v3,
            os.path.join(OUTDIR, f"global_v3_{v}_{GLOBAL_RES_DEG}deg.png"),
            vmin=vmin, vmax=vmax, cmap=cmap
        )

        # ---- Plot difference (skip for Mask) ----
        if v != "Mask":
            diff = new_map - v3_map
            vmax_diff = 2.0 if "OC" in v else 10.0

            plot_global_diff(
                diff, lat_bins, lon_bins,
                f"NEW − v3 {v} (coarse {GLOBAL_RES_DEG}°)",
                os.path.join(OUTDIR, f"global_DIFF_{v}_{GLOBAL_RES_DEG}deg.png"),
                vmax=vmax_diff
            )


    # Optional zoom
    if ZOOM is not None:
        name = ZOOM["name"]
        lon_min, lon_max = ZOOM["lon_min"], ZOOM["lon_max"]
        lat_min, lat_max = ZOOM["lat_min"], ZOOM["lat_max"]

        for v in PLOT_VARS:
            print("Building zoom:", name, v)

            new_zoom, lon2d, lat2d = build_zoom_native(
                NEW_DIR, v, lon_min, lon_max, lat_min, lat_max, step=ZOOM_STEP
            )
            v3_zoom, lon2d2, lat2d2 = build_zoom_native(
                V3_DIR, v, lon_min, lon_max, lat_min, lat_max, step=ZOOM_STEP
            )

            if new_zoom is None or v3_zoom is None:
                print("  zoom region had no data for", v)
                continue

            if v == "Mask":
                vmin, vmax = 0, 1
                cmap = "gray"
                title_new = f"NEW Mask zoom {name}\n1 = soil properties exist, 0 = no soil"
                title_v3  = f"v3 Mask zoom {name}\n1 = soil properties exist, 0 = no soil"
            else:
                vmin, vmax = 0, 100
                cmap = "jet"
                title_new = f"NEW {v} zoom {name} (step={ZOOM_STEP})"
                title_v3  = f"v3 {v} zoom {name} (step={ZOOM_STEP})"

            plot_zoom(
                new_zoom, lon2d, lat2d,
                title_new,
                os.path.join(OUTDIR, f"zoom_{name}_NEW_{v}_step{ZOOM_STEP}.png"),
                vmin=vmin, vmax=vmax, cmap=cmap
            )

            plot_zoom(
                v3_zoom, lon2d2, lat2d2,
                title_v3,
                os.path.join(OUTDIR, f"zoom_{name}_v3_{v}_step{ZOOM_STEP}.png"),
                vmin=vmin, vmax=vmax, cmap=cmap
            )

            # Difference (skip for Mask)
            if v != "Mask":
                diffz = new_zoom - v3_zoom
                vmax_diff = 2.0 if "OC" in v else 10.0
                plot_zoom_diff(
                    diffz, lon2d, lat2d,
                    f"NEW − v3 {v} zoom {name} (step={ZOOM_STEP})",
                    os.path.join(OUTDIR, f"zoom_{name}_DIFF_{v}_step{ZOOM_STEP}.png"),
                    vmax=vmax_diff
                )

    print("Wrote plots under:", OUTDIR)

if __name__ == "__main__":
    main()

