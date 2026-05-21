import glob, os
import numpy as np
import xarray as xr

#d = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/out_HWSDv2_NGDC_STATSGO_noMASK"
d = "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/soil/SOIL-DATA/soil_properties/v3/"

vars_ = ["Clay0_30","Sand0_30","OC0_30","Clay30_100","Sand30_100","OC30_100","Gravel"]
files = sorted(glob.glob(os.path.join(d, "SoilProperties_H??V??.nc")))
print("Tiles:", len(files))

def get_undef_and_sf(da):
    # STEP5 tiles use these non-CF attribute names
    undef = da.attrs.get("UNDEF", -9999)
    sf = float(da.attrs.get("ScaleFactor", 1.0))
    return undef, sf

gmin = {v: np.inf for v in vars_}
gmax = {v: -np.inf for v in vars_}
gsum = {v: 0.0 for v in vars_}
gcnt = {v: 0 for v in vars_}

for f in files:
    with xr.open_dataset(f, decode_cf=False, mask_and_scale=False) as ds:
        for v in vars_:
            da = ds[v]
            a = da.values
            undef, sf = get_undef_and_sf(da)

            m = (a != undef)
            n = int(m.sum())
            if n == 0:
                continue

            x = a[m].astype(np.float64) * sf  # percent
            gmin[v] = min(gmin[v], float(x.min()))
            gmax[v] = max(gmax[v], float(x.max()))
            gsum[v] += float(x.sum())
            gcnt[v] += n

print("\nGlobal stats (percent):")
for v in vars_:
    if gcnt[v] == 0:
        print(v, "NO VALID DATA")
    else:
        print(f"{v:12s} min={gmin[v]:8.3f} mean={gsum[v]/gcnt[v]:8.3f} max={gmax[v]:8.3f} n={gcnt[v]}")

# Range checks in scaled units
for v in ["Sand0_30", "Gravel"]:
    neg = 0
    gt100 = 0
    bad_tiles = []

    for f in files:
        with xr.open_dataset(f, decode_cf=False, mask_and_scale=False) as ds:
            da = ds[v]
            a = da.values
            undef, sf = get_undef_and_sf(da)

            m = (a != undef)
            if not np.any(m):
                continue

            x = a[m].astype(np.float64) * sf
            neg += int((x < 0).sum())
            gt100 += int((x > 100).sum())

            if np.any((x < 0) | (x > 100)):
                bad_tiles.append(os.path.basename(f))

    print(f"\n{v}: count(x<0)={neg} count(x>100)={gt100} out of n={gcnt[v]}")
    print(f"{v}: tiles with out-of-range values:", bad_tiles[:20], "count=", len(bad_tiles))

