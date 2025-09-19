import sys 
import numpy as np
from netCDF4 import Dataset
#Main purpose: Processes lake data to be used in the river routing model.

file_lat1m, file_lon1m, file_catmap, file_lake_mantag, file_lakecat_manfix = sys.argv[1:6]

# Define constants
nlat = 10800
nlon = 21600

# Read data from files
lats = np.loadtxt("temp/outlet_lat.txt", dtype=float)  # Latitude of outlets
lons = np.loadtxt("temp/outlet_lon.txt", dtype=float)  # Longitude of outlets
lat1m = np.loadtxt(file_lat1m, dtype=float)  # Latitude grid
lon1m = np.loadtxt(file_lon1m, dtype=float)  # Longitude grid

# Function to find the nearest index in a coordinate array
def ind_nearest_coord(coord_array1, coord_array2):
    """
    Find the index of the nearest value in coord_array2 for each value in coord_array1.
    """
    indices = []
    for coord in coord_array1:
        index = np.argmin(np.abs(coord_array2 - coord))
        indices.append(index)
    return np.array(indices)

# Find nearest indices for latitudes and longitudes
lati = ind_nearest_coord(lats, lat1m)+1
loni = ind_nearest_coord(lons, lon1m)+1

#------------------------------------------------------------------------------------------------------
ns = 3917

# Allocate array
catchind = np.zeros((nlat, nlon), dtype=int)

# Read NetCDF file
def read_ncfile_int2d(filepath, varname, shape):
    # Open the NetCDF file and read the specified variable
    with Dataset(filepath, 'r') as nc:
        data = nc.variables[varname][:].reshape(shape)
        # Check for missing values and replace them with a default value
        fill_value = nc.variables[varname]._FillValue if hasattr(nc.variables[varname], '_FillValue') else None
        if fill_value is not None:
            data = np.where(data == fill_value, -9999, data)  # Replace missing values with 0
    return data

catchind = read_ncfile_int2d(file_catmap, "CatchIndex", (nlat, nlon))

# Calculate catid
catid = np.zeros(ns, dtype=int)
for i in range(ns):
    # Ensure indices are within bounds
    if 0 < loni[i] <= nlon and 0 < lati[i] <= nlat:
        catid[i] = catchind[lati[i] - 1, loni[i] - 1]  # Adjust for 0-based indexing in Python
    else:
        catid[i] = -1  # Assign a default value for out-of-bounds indices

#------------------------------------------------------------------------------------------------------
# Constants
nall = 291284
nv = 1782
nv3 = 2097

# Read input data
aca_all = np.loadtxt("temp/Pfaf_acar.txt")

# Initialize aca_model array
aca_model = np.full(ns, -9999.0)

# Map aca_model using catid
for i in range(ns):
    if catid[i] != -9999:
        cid = catid[i]
        aca_model[i] = aca_all[cid - 1]

# Read observation data
aca_obs = np.loadtxt("temp/outlet_lakeacaOBS.txt")
outid_INCON = np.zeros(nv, dtype=int)

# Filter inconsistent data
k = 0
for i in range(ns):
    if not (0.7 * aca_model[i] <= aca_obs[i] <= 1.3 * aca_model[i]):
        outid_INCON[k] = i + 1
        k += 1

#print(k)

#------------------------------------------------------------------------------------------------------
# Read tags
tag_INCON = np.loadtxt(file_lake_mantag, dtype=int)

# Update catid and aca_model based on tags
for i in range(nv):
    oid = outid_INCON[i]
    tag = tag_INCON[i]
    if tag >= 1:
        cid = tag
        catid[oid - 1] = cid
        aca_model[oid - 1] = aca_all[cid - 1]
    else:
        catid[oid - 1] = -9999
        aca_model[oid - 1] = -9999.0

# Compute flag_out
flag_out = np.where(aca_model != -9999, 1, 0)
#print(np.sum(flag_out))

#------------------------------------------------------------------------------------------------------
# Read lakeid_out and compute absolute differences
lakeid_out = np.loadtxt("temp/outlet_lakeid.txt", dtype=int)
acaABSDIF_out = np.abs(aca_model - aca_obs)

# Initialize collections
lakeid_collect = np.zeros(nv3, dtype=int)
outletid_collect = np.zeros(nv3, dtype=int)
acaABSDIF_collect = np.full(nv3, 1e10)
flag_2097_out = np.zeros(ns, dtype=int)
k = 0

# Collect valid outlets
for i in range(ns):
    if flag_out[i] == 1:
        lid = lakeid_out[i]
        flag = 1
        if k >= 1:
            for j in range(k):
                if lid == lakeid_collect[j]:
                    flag = 0
                    if acaABSDIF_out[i] < acaABSDIF_collect[j]:
                        flag_2097_out[outletid_collect[j]] = 0
                        outletid_collect[j] = i
                        acaABSDIF_collect[j] = acaABSDIF_out[i]
                        flag_2097_out[i] = 1
        if flag == 1:
            lakeid_collect[k] = lid
            outletid_collect[k] = i
            acaABSDIF_collect[k] = acaABSDIF_out[i]
            flag_2097_out[i] = 1
            k += 1

#print(np.sum(flag_2097_out))
np.savetxt("output/lake_outlet_flag_valid_2097.txt", flag_2097_out, fmt="%d")

#------------------------------------------------------------------------------------------------------
# Update catid with valid flags
catid = np.where(flag_2097_out == 0, -9999, catid)

# Collect valid outlet IDs
outidV = np.zeros(nv3, dtype=int)
k = 0
for i in range(ns):
    if flag_2097_out[i] == 1:
        outidV[k] = i
        k += 1

outidV += 1

# Fix multiple outlets in same catchment
catid_outfix_2097 = np.loadtxt(file_lakecat_manfix, dtype=int)
catid_outfix_out = np.full(ns, -9999, dtype=int)

for i in range(nv3):
    oid = outidV[i]
    catid_outfix_out[oid - 1] = catid_outfix_2097[i]

catid = np.where((catid_outfix_out != 0) & (catid_outfix_out != -9999), catid_outfix_out, catid)
np.savetxt("output/lake_outlet_catid.txt", catid, fmt="%d")


