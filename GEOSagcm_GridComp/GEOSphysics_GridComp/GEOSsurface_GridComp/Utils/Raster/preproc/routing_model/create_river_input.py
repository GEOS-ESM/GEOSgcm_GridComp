import numpy as np
from netCDF4 import Dataset

# Dimensions
nc = 291284
nr = 7250
nl = 3917

def read_ascii(filename, count, dtype=float):
    return np.loadtxt(filename, dtype=dtype, max_rows=count)

# ---------- Catchment datasets ----------
Kstr_catchment = read_ascii("temp/Pfaf_Kstr_PR_fac1_0p35_0p45_0p2_n0p2.txt", nc)
Kv_catchment = read_ascii("temp/Pfaf_Kv_PR_0p35_0p45_0p2_n0p2.txt", nc)
area_catchment = read_ascii("temp/Pfaf_area.txt", nc)
lriv_catchment = read_ascii("temp/Pfaf_lriv_PR.txt", nc)
lstr_catchment = read_ascii("temp/Pfaf_lstr_PR.txt", nc)
Qin_catchment = read_ascii("temp/Pfaf_qin.txt", nc)
Qri_catchment = read_ascii("temp/Pfaf_qri.txt", nc)
Qstr_catchment = read_ascii("temp/Pfaf_qstr.txt", nc)
tosink_catchment = read_ascii("temp/Pfaf_tosink.txt", nc, dtype=int)
dnid_catchment = read_ascii("temp/downstream_1D_new_noadj.txt", nc, dtype=int)

# ---------- Reservoir datasets ----------
area_reservoir = read_ascii("temp/area_skm_grand.txt", nr)
capmax_reservoir = read_ascii("temp/cap_max_grand.txt", nr)
catid_reservoir = read_ascii("temp/catid_dam_corr_aca_grand5000.txt", nr, dtype=int)
flag_reservoir = read_ascii("temp/flag_all_res.txt", nr, dtype=int)
fldmainsec_reservoir = read_ascii("temp/fldmainsec_grand.txt", nr, dtype=int)
elec_reservoir = read_ascii("temp/hydroelec_grand.txt", nr, dtype=int)
irr_reservoir = read_ascii("temp/irr_grand.txt", nr, dtype=int)
nav_reservoir = read_ascii("temp/nav_grand.txt", nr, dtype=int)
other_reservoir = read_ascii("temp/other_grand.txt", nr, dtype=int)
rec_reservoir = read_ascii("temp/rec_grand.txt", nr, dtype=int)
supply_reservoir = read_ascii("temp/watersupply_grand.txt", nr, dtype=int)

# ---------- Lake datasets ----------
catid_lake = read_ascii("temp/lake_outlet_catid.txt", nl, dtype=int)
flag_lake = read_ascii("temp/lake_outlet_flag_valid_2097.txt", nl, dtype=int)
area_lake = read_ascii("temp/lake_outlet_lakearea.txt", nl)

# ---------- Create NetCDF ----------
fout = Dataset("output/river_input.nc", "w", format="NETCDF4")

# define dimensions
fout.createDimension("pfaf", nc)
fout.createDimension("reservoir", nr)
fout.createDimension("lake", nl)

# ---------- Create variables ----------
def create_var(name, data, dims, units, long_name):
    var = fout.createVariable(name, data.dtype, dims)
    var[:] = data
    var.units = units
    var.long_name = long_name
    return var

# Catchment variables
create_var("Kstr", Kstr_catchment, ("pfaf",), "1", "K parameter for local streams")
create_var("K", Kv_catchment, ("pfaf",), "1", "K parameter for main rivers")
create_var("lengsc", lriv_catchment, ("pfaf",), "km", "main river length scale")
create_var("lstr", lstr_catchment, ("pfaf",), "km", "local streams length scale")
create_var("qin_clmt", Qin_catchment, ("pfaf",), "m3/s", "climatology of catchment inflow")
create_var("qri_clmt", Qri_catchment, ("pfaf",), "m3/s", "climatology of catchment outflow")
create_var("qstr_clmt", Qstr_catchment, ("pfaf",), "m3/s", "climatology of catchment stream flow")
create_var("downid", dnid_catchment, ("pfaf",), "1", "downstream catchment id")
create_var("area_catch", area_catchment, ("pfaf",), "km2", "catchment area")

# Reservoir variables
create_var("area_grand", area_reservoir, ("reservoir",), "km2", "reservoir area")
create_var("cap_grand", capmax_reservoir, ("reservoir",), "million m3", "reservoir maximum capacity")
create_var("catid_grand", catid_reservoir, ("reservoir",), "1", "catchment id for reservoir")
create_var("flag_grand", flag_reservoir, ("reservoir",), "1", "active flag for reservoir")
create_var("fld_grand", fldmainsec_reservoir, ("reservoir",), "1", "flood-control flag")
create_var("elec_grand", elec_reservoir, ("reservoir",), "1", "hydro-power flag")
create_var("irr_grand", irr_reservoir, ("reservoir",), "1", "irrigation flag")
create_var("nav_grand", nav_reservoir, ("reservoir",), "1", "navigation flag")
create_var("other_grand", other_reservoir, ("reservoir",), "1", "other use flag")
create_var("rec_grand", rec_reservoir, ("reservoir",), "1", "recreation flag")
create_var("supply_grand", supply_reservoir, ("reservoir",), "1", "water supply flag")

# Lake variables
create_var("catid_lake", catid_lake, ("lake",), "1", "catchment id for lake outlet")
create_var("flag_lake", flag_lake, ("lake",), "1", "active flag for lake")
create_var("area_lake", area_lake, ("lake",), "km2", "lake area")

# Close file
fout.close()
print("river_input.nc created successfully!")
