import sys  
import numpy as np
from netCDF4 import Dataset
import os
import glob

#Main purpose: Processes reservoir (dam) data: reads dam locations and usage information from GRanD database.

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

if __name__ == '__main__':

    file_latdam, file_londam, file_lat1m, file_lon1m, file_catmap, file_acadam, file_damcat_manfix, file_dam_manflag, file_dam_use, file_flood = sys.argv[1:11]
    # Parameter settings
    ns = 7250 #nr
    nlat = 10800
    nlon = 21600
    nc = 291284
    thres = 5000.0    

#----get dam lat lon ind----

    # Read data from ASCII files
    lats = np.loadtxt(file_latdam, dtype=np.float64, max_rows=ns)
    lons = np.loadtxt(file_londam, dtype=np.float64, max_rows=ns)
    lat1m = np.loadtxt(file_lat1m, dtype=np.float64, max_rows=nlat)
    lon1m = np.loadtxt(file_lon1m, dtype=np.float64, max_rows=nlon)

    # For each target coordinate, find the nearest index in the reference array.
    lati = ind_nearest_coord(lats, lat1m) 
    loni = ind_nearest_coord(lons, lon1m) 

#----get dam catchemnt indices----
    # Read the catchment indices file
    nc_file = file_catmap
    with Dataset(nc_file, 'r') as ncf:
        catchind = ncf.variables["CatchIndex"][:]
        if np.ma.is_masked(catchind):
            catchind = catchind.filled(-9999)

    # Initialize the array to store the output catchment ID values
    catid = np.empty(ns, dtype=int)

    # Loop over each index
    for i in range(ns):
        # For each index, retrieve the value from catchind.
        catid[i] = catchind[ lati[i], loni[i] ]

#----get dam drainage area--
    # Read full dataset for acar(drainage area) and catchment area from ASCII files
    acar_all = np.loadtxt("temp/Pfaf_acar.txt", dtype=float, max_rows=nc)
    area_all = np.loadtxt("output/Pfaf_area.txt", dtype=float, max_rows=nc)

    # Initialize arrays to store the selected values for each catchment
    acar = np.empty(ns, dtype=float)
    area = np.empty(ns, dtype=float)

    # Loop over each catchment index and assign values based on catid
    for i in range(ns):
        cid = catid[i]
        if cid != -9999:
            # Subtract 1 from cid to convert 1-based index to 0-based index for Python
            acar[i] = acar_all[cid - 1]
            area[i] = area_all[cid - 1]
        else:
            acar[i] = -9999.0
            area[i] = -9999.0

#----look for station: model drainage area is too small------
    # Read drainage area from GRAND database
    grand = np.loadtxt(file_acadam, dtype=float, max_rows=ns)

    # Initialize lists to store error information
    id_error = []

    # Loop over each catchment index
    for i in range(ns):
        # At here only care about large-scale dams
        if grand[i] > thres:
            if acar[i] < 0.8 * grand[i]:
                # Append error information; add 1 to i for 1-based indexing
                id_error.append(i + 1)

    # Get the number of errors found
    ne = len(id_error)

#----get corrected catid for above station--------------------
    # Read manually corrected catid array from ASCII files.
    catid_error = np.loadtxt(file_damcat_manfix, dtype=int, max_rows=ne)

    # Loop over each error index and update catid.
    # Note: We subtract 1 from resid_error values to convert from 1-based to 0-based indexing.
    for i in range(ne):
        rid = id_error[i]
        catid[rid - 1] = catid_error[i]

#----get dam drainage area with corrected catid--------------------
    # Initialize arrays to store the selected acar and area values for each catchment
    acar = np.empty(ns, dtype=float)
    area = np.empty(ns, dtype=float)

    # Loop over each catchment index
    for i in range(ns):
        cid = catid[i]
        if cid != -9999:
            # Adjust for 1-based indexing: subtract 1 when accessing the full dataset arrays
            acar[i] = acar_all[cid - 1]
            area[i] = area_all[cid - 1]
        else:
            acar[i] = -9999.0
            area[i] = -9999.0

#----look for station: model drainage area is too large------
    model = acar
    # we use a list to collect error indices.
    id_error = []
    
    # Loop over each catchment index
    for i in range(ns):
        # Check if the model value is greater than the threshold and grand is less than 80% of model
        if model[i] > thres:
            if grand[i] < 0.8 * model[i]:
                # Append 1-based index (i+1) to the error list
                id_error.append(i + 1)
                
    ne = len(id_error)

#----create flag for all dams------

    # Three more manual adjustment dams
    resid_man = np.array([5179, 289, 7070], dtype=int)
    catid_man = np.array([46616, 142851, 199281], dtype=int)
    nman = resid_man.size

    # Update specific indices in catid_all with manual adjustments.
    for i in range(nman):
        catid[resid_man[i] - 1] = catid_man[i]

    # Write the updated catid_all to an ASCII file
    np.savetxt("output/catid_dam_corr_aca_grand5000.txt", catid, fmt='%d')

    # Read dams flag (whether we still need it in the model) for the above uncorrect dams from a manually checked flag file
    flag_error = np.loadtxt(file_dam_manflag, dtype=int, max_rows=ne)

    # Initialize flag_all array with ones (default flag value)
    flag_all = np.ones(ns, dtype=int)

    # For each error entry, update flag_all at the specified index.
    # Adjust id from 1-based to 0-based indexing.
    for i in range(ne):
        id_val = id_error[i]
        flag_all[id_val - 1] = flag_error[i]

    # If drainage area in GRAND is small and also less than 0.5 times model drainage area, we do not need the dam.
    for i in range(ns):
        if grand[i] < 1.e3:
            if grand[i] < 0.5 * acar[i]:
                flag_all[i] = 0

    # If model drainage area is negative, set flag_all to 0 for that dam.
    for i in range(ns):
        if acar[i] < 0.:
            flag_all[i] = 0

    # Write the final flag_all array to an ASCII file.
    np.savetxt("output/flag_all_res.txt", flag_all, fmt='%d')

#----get dam main use---------------
    # Define category strings and corresponding output tags
    use_string = ["Irrigation", "Hydroelectricity", "Water supply", "Navigation", "Recreation"]
    use_out = ["irr", "hydroelec", "watersupply", "nav", "rec"]
    nu = len(use_string)

    # Read the main use data as strings from the GRAND file
    with open(file_dam_use, "r") as f:
        use = [line.strip() for line in f]
    if len(use) != ns:
        print(f"Warning: expected {ns} lines, but got {len(use)} lines.")

    # For each category in use_string, create a flag array and output the result
    for j in range(nu):
        flag = np.zeros(ns, dtype=int)
        # Set flag to 1 where the use value matches the current category
        for i in range(ns):
            if use[i] == use_string[j]:
                flag[i] = 1

        # Write the flag array
        out_filename = os.path.join("output", use_out[j] + "_grand.txt")
        np.savetxt(out_filename, flag, fmt='%d')    
#----flood use--------------------
    # Read the use_irr strings from the GRAND file
    with open(file_flood, "r") as f:
        use_irr = [line.strip() for line in f]

    # Initialize the flag array with zeros
    flag = np.zeros(ns, dtype=int)

    # Loop over each entry and set flag to 1 if the entry is not "NA"
    for i in range(ns):
        if use_irr[i] != "NA":
            flag[i] = 1

    # Write the flag array
    np.savetxt("output/fldmainsec_grand.txt", flag, fmt='%d')

#----other use--------------------
    use_out = "other"

    # Read the main use data from the GRAND file
    with open(file_dam_use, "r") as f:
        use = [line.strip() for line in f]
    if len(use) != ns:
        print(f"Warning: expected {ns} entries, but got {len(use)} entries.")

    # Initialize the flag array with zeros
    flag = np.zeros(ns, dtype=int)

    # Loop over each entry and set flag to 1 if the entry matches the specified categories
    for i in range(ns):
        if use[i] == "Fisheries" or use[i] == "NA" or use[i] == "Other":
            flag[i] = 1

    # Write the flag array
    np.savetxt("output/" + use_out + "_grand.txt", flag, fmt='%d')

