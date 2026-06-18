import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os


# ---- USER SETTINGS ----
# discover base current
DIR_OLD = "/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/soil/SOIL-DATA/soil_properties/v3/"
# discover what I produce
#DIR_NEW = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/out_HWSDv2_NGDC_STATSGO_AFSIS_noMASK/"   # w afsis
#DIR_NEW = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/out_HWSDv2_NGDC_STATSGO_noMASK/"   # no afsis
DIR_NEW = "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v4/"   # no afsis
TILE    = "SoilProperties_H20V09.nc" # africa tile
#TILE    = "SoilProperties_H10V17.nc"  #"SoilProperties_H13V06.nc"   # change per region
#TILE    =  "SoilProperties_H13V06.nc"   # change per region
VAR     = "Clay0_30"
# ----------------------


f_old = os.path.join(DIR_OLD, TILE)
f_new = os.path.join(DIR_NEW, TILE)

def get_undef_sf(da):
    undef = da.attrs.get("UNDEF", -9999)
    sf = float(da.attrs.get("ScaleFactor", 1.0))
    return undef, sf

with xr.open_dataset(f_old, decode_cf=False, mask_and_scale=False) as a, \
     xr.open_dataset(f_new, decode_cf=False, mask_and_scale=False) as b:

    xa = a[VAR].values.astype(np.float64)
    yb = b[VAR].values.astype(np.float64)

    undef_a, sf_a = get_undef_sf(a[VAR])
    undef_b, sf_b = get_undef_sf(b[VAR])

    # Convert each to physical percent first
    x = np.where(xa != undef_a, xa * sf_a, np.nan)
    y = np.where(yb != undef_b, yb * sf_b, np.nan)

diff = y - x  # percent (percentage points)

plt.figure(figsize=(6, 5))
im = plt.imshow(diff, origin="lower", cmap="RdBu_r", vmin=-20, vmax=20)
plt.colorbar(im, label=f"{VAR} diff (percentage points)")
plt.title(f"{VAR} NEW - OLD\n{TILE}")
plt.tight_layout()
plt.show()

old_undef = ~np.isfinite(x)
new_undef = ~np.isfinite(y)

print("old undef %:", 100*old_undef.mean())
print("new undef %:", 100*new_undef.mean())
print("old undef -> new filled %:", 100*np.mean(old_undef & ~new_undef))
print("old filled -> new undef %:", 100*np.mean(~old_undef & new_undef))

both = (~old_undef) & (~new_undef)
d = diff[both]
print("both-filled count:", both.sum())
print("diff stats (pp): min/mean/max =", np.min(d), np.mean(d), np.max(d))
print("fraction |diff|>1pp:", np.mean(np.abs(d) > 1.0))
print("fraction |diff|>5pp:", np.mean(np.abs(d) > 5.0))

