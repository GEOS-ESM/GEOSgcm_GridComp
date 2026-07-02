#!/usr/bin/env python3
import os
import sys
import time
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
###
#python3 plot_stich_global_var_layer.py Sand D2
#python3 plot_stich_global_var_layer.py Silt D2
#python3 plot_stich_global_var_layer.py Clay D2

nc_dir = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/output_STEP1/"

VAR_RANGES = {
    "Sand": (0, 100),
    "Silt": (0, 100),
    "Clay": (0, 100),
    "Coarse": (0, 100),
    "SOC": (0, 100),
    "Organic_Matter": (0, 100),
    "Share": (0, 100),
    "Porosity": (0, 1),
    "Bulk_Density": (0, 2.5),
    "Ref_Bulk_Density": (0.5, 2.5),
    "PH": (0, 10),
    "Sequence_Used": (0, 8),
}

# Zoom control:
# - set to None for global
# - or set to [lonmin, lonmax, latmin, latmax] for zoom
zoom_extent = None
# zoom_extent = [-110, -100, 20, 30]  # example zoom


def stitch_plot(variable, layer, extent=None, heartbeat=25):
    t0 = time.time()

    files = sorted(f for f in os.listdir(nc_dir) if f.endswith(f"_{layer}.nc"))
    if not files:
        print(f"[WARN] No files found for layer {layer} in {nc_dir}")
        return 1

    print(f"\n[INFO] Start plot: variable={variable} layer={layer} tiles={len(files)} extent={extent}")

    vmin, vmax = VAR_RANGES.get(variable, (None, None))
    norm = Normalize(vmin=vmin, vmax=vmax) if (vmin is not None and vmax is not None) else None
    cmap = plt.colormaps["YlGnBu"]

    fig, ax = plt.subplots(figsize=(14, 7), subplot_kw={"projection": ccrs.PlateCarree()})
    plotted = False

    for i, fn in enumerate(files, start=1):
        path = os.path.join(nc_dir, fn)

        if i == 1 or (i % heartbeat == 0) or (i == len(files)):
            print(f"[INFO]   tile {i:3d}/{len(files)}: {fn}")

        # Use context manager to ensure file handles close even on errors
        with xr.open_dataset(path) as ds:
            if variable not in ds:
                continue

            da = ds[variable].transpose("lat", "lon")
            lat = ds["lat"].values
            lon = ds["lon"].values

            data = np.ma.masked_invalid(da.values)
            lon_grid, lat_grid = np.meshgrid(lon, lat)

            ax.pcolormesh(
                lon_grid, lat_grid, data,
                cmap=cmap, norm=norm,
                shading="auto",
                transform=ccrs.PlateCarree()
            )
            plotted = True

    if not plotted:
        print(f"[WARN] Variable '{variable}' not found in any tiles for layer {layer}")
        plt.close()
        return 2

    ax.coastlines(resolution="10m" if extent else "110m")
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color="gray", alpha=0.5, linestyle="--")
    gl.top_labels = False
    gl.right_labels = False

    if extent:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        tag = f"zoom_lon{extent[0]}_{extent[1]}_lat{extent[2]}_{extent[3]}"
    else:
        ax.set_global()
        tag = "global"

    ax.set_title(f"{variable} ({layer}) — {tag}", pad=20)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, orientation="vertical", label=variable)

    out = f"{variable}_{layer}_{tag}.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()

    dt = time.time() - t0
    print(f"[DONE] Saved: {out}   ({dt:.1f} s)")
    return 0


if __name__ == "__main__":
    # Usage:
    #   python3 plot_stich_global_var_layer.py Sand D2
    # If not provided, defaults:
    variable = sys.argv[1] if len(sys.argv) > 1 else "Sand"
    layer = sys.argv[2] if len(sys.argv) > 2 else "D2"

    print("\n========== HWSD2 Stitch Plotter ==========")
    print(f"[INFO] nc_dir:   {nc_dir}")
    print(f"[INFO] variable: {variable}")
    print(f"[INFO] layer:    {layer}")
    print(f"[INFO] zoom:     {zoom_extent}")
    print("=========================================\n")

    raise SystemExit(stitch_plot(variable, layer, extent=zoom_extent))

