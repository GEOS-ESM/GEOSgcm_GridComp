import numpy as np
from netCDF4 import Dataset
import routing_model_constants

nc_len = routing_model_constants.nc
nres = routing_model_constants.nres
nlake = routing_model_constants.nlake

def read_ascii(filename, count, dtype=float):
    return np.loadtxt(filename, dtype=dtype, max_rows=count)

# ---------- Catchment datasets ----------
Kstr_catchment = read_ascii("temp/Pfaf_Kstr_PR_fac1_0p35_0p45_0p2_n0p2.txt", nc_len)
Kv_catchment = read_ascii("temp/Pfaf_Kv_PR_0p35_0p45_0p2_n0p2.txt", nc_len)
area_catchment = read_ascii("temp/Pfaf_area.txt", nc_len)
lriv_catchment = read_ascii("temp/Pfaf_lriv_PR.txt", nc_len)
lstr_catchment = read_ascii("temp/Pfaf_lstr_PR.txt", nc_len)
Qin_catchment = read_ascii("temp/Pfaf_qin.txt", nc_len)
Qri_catchment = read_ascii("temp/Pfaf_qri.txt", nc_len)
Qstr_catchment = read_ascii("temp/Pfaf_qstr.txt", nc_len)
tosink_catchment = read_ascii("temp/Pfaf_tosink.txt", nc_len, dtype=int)
dnid_catchment = read_ascii("temp/downstream_1D_new_noadj.txt", nc_len, dtype=int)
area_catchment = area_catchment * 1.e6
lriv_catchment = lriv_catchment * 1.e3
lstr_catchment = lstr_catchment * 1.e3

# ---------- Reservoir datasets ----------
area_grand = read_ascii("temp/area_skm_grand.txt", nres)
cap_grand = read_ascii("temp/cap_max_grand.txt", nres)
catid_grand = read_ascii("temp/catid_dam_corr_aca_grand5000.txt", nres, dtype=int)
flag_grand = read_ascii("temp/flag_all_res.txt", nres, dtype=int)
fld_grand = read_ascii("temp/fldmainsec_grand.txt", nres, dtype=int)
elec_grand = read_ascii("temp/hydroelec_grand.txt", nres, dtype=int)
irr_grand = read_ascii("temp/irr_grand.txt", nres, dtype=int)
nav_grand = read_ascii("temp/nav_grand.txt", nres, dtype=int)
other_grand = read_ascii("temp/other_grand.txt", nres, dtype=int)
rec_grand = read_ascii("temp/rec_grand.txt", nres, dtype=int)
supply_grand = read_ascii("temp/watersupply_grand.txt", nres, dtype=int)

# ---------- Lake datasets ----------
catid_lake = read_ascii("temp/lake_outlet_catid.txt", nlake, dtype=int)
flag_lake = read_ascii("temp/lake_outlet_flag_valid_2097.txt", nlake, dtype=int)
area_lake = read_ascii("temp/lake_outlet_lakearea.txt", nlake)


# -------------------------------------------------------------------------
# Unit Conversions
# -------------------------------------------------------------------------
# Convert capacity from million cubic meters (MCM) to m3
cap_grand = cap_grand * 1.e6
area_grand= area_grand * 1.e6
area_lake = area_lake * 1.e6

# -------------------------------------------------------------------------
# Initialize Reservoir Properties (Arrays of size nc_len)
# -------------------------------------------------------------------------
cap_res = np.zeros(nc_len, dtype=np.float32)
area_res = np.zeros(nc_len, dtype=np.float32)
area_max_res = np.zeros(nc_len, dtype=np.float32)
wid_res = np.zeros(nc_len, dtype=np.float32)
    
type_res = np.zeros(nc_len, dtype=int)
fld_res = np.zeros(nc_len, dtype=int)
active_res = np.zeros(nc_len, dtype=int)

# -------------------------------------------------------------------------
# Aggregation Logic
# -------------------------------------------------------------------------
#print("Processing reservoir aggregation...")

# Filter for active reservoirs (flag_grand == 1)
active_mask = flag_grand == 1
    
# Get indices (adjusting Fortran 1-based index to Python 0-based index)
# Using catid - 1 to map to Python array indices
curr_cids = catid_grand[active_mask] - 1
    
# Verify indices are within bounds
valid_idx_mask = (curr_cids >= 0) & (curr_cids < nc_len)
curr_cids = curr_cids[valid_idx_mask]
    
# Slice the data arrays with the same masks
# We need to filter the original arrays by active_mask, then by valid_idx_mask
filtered_cap = cap_grand[active_mask][valid_idx_mask]
filtered_area = area_grand[active_mask][valid_idx_mask]
    
# Sum up capacities and areas for reservoirs in the same catchment
# np.add.at is equivalent to doing a loop and adding values to specific indices
np.add.at(cap_res, curr_cids, filtered_cap)
np.add.at(area_res, curr_cids, filtered_area)
    
# Handle flood control flag
# If any active reservoir in the catchment has flood control, set fld_res to 1
# We iterate to simulate the "if(fld_grand(i) == 1)" logic
fld_indices = np.where((flag_grand == 1) & (fld_grand == 1))[0]
fld_cids = catid_grand[fld_indices] - 1
# Filter valid bounds
fld_cids = fld_cids[(fld_cids >= 0) & (fld_cids < nc_len)]
fld_res[fld_cids] = 1

# Handle missing capacity data
cap_res[cap_res == 0.0] = -1.0

# -------------------------------------------------------------------------
# Assign Reservoir Types
# -------------------------------------------------------------------------
# Other(6) -> Rec(5) -> Nav(4) -> Supply(3) -> Elec(2) -> Irr(1).
# Later loops overwrite earlier ones if (area >= area_max).
    
#print("Assigning reservoir types...")

# Define the order and arrays corresponding to types 6 down to 1
type_configs = [
    (6, other_grand),
    (5, rec_grand),
    (4, nav_grand),
    (3, supply_grand),
    (2, elec_grand),
    (1, irr_grand)
]

for type_id, type_array in type_configs:
    # Find all active reservoirs of this type
    # Condition: flag_grand == 1 AND type_array == 1
    mask = (flag_grand == 1) & (type_array == 1)
    indices = np.where(mask)[0]
        
    for idx in indices:
        cid = catid_grand[idx] - 1 # 0-based index
            
        if 0 <= cid < nc_len:
            # Check if this reservoir is the largest processed so far for this catchment
            if area_grand[idx] >= area_max_res[cid]:
                type_res[cid] = type_id
                area_max_res[cid] = area_grand[idx]

# -------------------------------------------------------------------------
# Set Up Natural Lakes
# -------------------------------------------------------------------------
#print("Processing natural lakes...")
    
# Filter active lakes with valid catid
lake_mask = (flag_lake == 1) & (catid_lake > 0)
lake_indices = np.where(lake_mask)[0]

for idx in lake_indices:
    cid = catid_lake[idx] - 1 # 0-based index
        
    if 0 <= cid < nc_len:
        # If no reservoir type assigned yet AND no flood control
        if type_res[cid] == 0 and fld_res[cid] == 0:
            type_res[cid] = -1 # Type for lake
            area_res[cid] = area_lake[idx]

# -------------------------------------------------------------------------
# Final Calculations and Cleanup
# -------------------------------------------------------------------------
# Compute width (sqrt of area)
wid_res = np.sqrt(area_res)
    
# Mark active reservoirs based on type or flood control status
active_res[:] = 0 # Ensure clean state
# Condition: type_res != 0 OR fld_res == 1
active_condition = (type_res != 0) | (fld_res == 1)
active_res[active_condition] = 1

#print("Processing complete.")

# ---------- Create NetCDF ----------
fout = Dataset("output/route_parameters.nc", "w", format="NETCDF4")

# define dimensions
fout.createDimension("tile", nc_len)

# ---------- Create variables ----------
def create_var(name, data, dims, units, long_name):
    var = fout.createVariable(name, data.dtype, dims)
    var[:] = data
    var.units = units
    var.long_name = long_name
    return var

# Catchment variables
create_var("KSTR", np.float32(Kstr_catchment), ("tile",), "1", "K_parameter_for_local_streams")
create_var("K", np.float32(Kv_catchment), ("tile",), "1", "K_parameter_for_main_rivers")
create_var("LENGSC", np.float32(lriv_catchment), ("tile",), "m", "main_river_length_scale")
create_var("LSTR", np.float32(lstr_catchment), ("tile",), "m", "local_streams_length_scale")
create_var("QIN_CLMT", np.float32(Qin_catchment), ("tile",), "m+3 s-1", "climatology_of_catchment_inflow")
create_var("QRI_CLMT", np.float32(Qri_catchment), ("tile",), "m+3 s-1", "climatology_of_catchment_outflow")
create_var("QSTR_CLMT", np.float32(Qstr_catchment), ("tile",), "m+3 s-1", "climatology_of_catchment_stream_flow")
create_var("DOWNID", np.float32(dnid_catchment), ("tile",), "1", "downstream_catchment_id")
create_var("AREA_CATCH", np.float32(area_catchment), ("tile",), "m+2", "catchment_area")

# Reservoir variables
create_var("ACTIVE_RES", np.float32(active_res), ("tile",), "1", "active_reservoirs")
create_var("CAP_RES", np.float32(cap_res), ("tile",), "m+3", "max_capacity_of_reservoirs")
create_var("FLD_RES", np.float32(fld_res), ("tile",), "1", "flood_control_flag_of_reservoirs")
create_var("TYPE_RES", np.float32(type_res), ("tile",), "1", "type_of_reservoirs")
create_var("WID_RES", np.float32(wid_res), ("tile",), "m", "length_scale_of_reservoirs")


# Close file
fout.close()
print("route_parameters.nc created successfully!")
