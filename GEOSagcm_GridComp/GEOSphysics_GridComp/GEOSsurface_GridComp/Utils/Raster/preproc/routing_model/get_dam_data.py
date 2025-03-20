import numpy as np
from netCDF4 import Dataset
import os
import glob

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

#----get dam lat lon ind----
    # Parameter settings
    ns = 7250 #nr
    nlat = 10800
    nlon = 21600

    # Read data from ASCII files; each file contains one number per line.
    lats = np.loadtxt("input/lat_dam_grand.txt", dtype=np.float64, max_rows=ns)
    lons = np.loadtxt("input/lon_dam_grand.txt", dtype=np.float64, max_rows=ns)
    lat1m = np.loadtxt("input/lat_1m.txt", dtype=np.float64, max_rows=nlat)
    lon1m = np.loadtxt("input/lon_1m.txt", dtype=np.float64, max_rows=nlon)

    # For each target coordinate, find the nearest index in the reference array.
    lati = ind_nearest_coord(lats, lat1m) 
    loni = ind_nearest_coord(lons, lon1m) 

    # Since NCL uses 1-based indexing, add 1 when writing the results.
    #np.savetxt("data/lati_dam_PR.txt", lati, fmt='%d')
    #np.savetxt("data/loni_dam_PR.txt", loni, fmt='%d')


#----get dam cat ind----
    # Read the NetCDF file
    nc_file = "input/CatchIndex.nc"
    with Dataset(nc_file, 'r') as nc:
        # Read the integer 2D variable "data"
        # We assume that the data is stored with shape (nlon, nlat)
        catchind = nc.variables["data"][:]
        # If catchind is a masked array, fill masked values with -9999
        if np.ma.is_masked(catchind):
            catchind = catchind.filled(-9999)
    # Read ASCII files containing indices for latitude and longitude.
    # It is assumed that the files contain one integer per line and use 1-based indexing.
    #lati = np.loadtxt("data/lati_dam.txt", dtype=int)  # shape: (ns,)
    #loni = np.loadtxt("data/loni_dam.txt", dtype=int)  # shape: (ns,)

    # Initialize the array to store the output catchment ID values
    catid = np.empty(ns, dtype=int)

    # Loop over each index (Fortran indices start at 1, so we subtract 1 for 0-based indexing in Python)
    for i in range(ns):
        # For each index, retrieve the value from catchind using (loni, lati) as indices.
        # Subtract 1 from the read indices to convert from 1-based to 0-based indexing.
        catid[i] = catchind[ lati[i], loni[i] ]

    # Write the output catid array to an ASCII file with one number per line.
    #np.savetxt("data/catid_dam_PR.txt", catid, fmt='%d')

#----get dam drainage area--
    # Define the number of catchments and the total number of entries in the full dataset
    nc = 291284

    # Alternative file (commented out in the original NCL script):
    # catid = np.loadtxt("data/catid_dam_corr_aca_grand5000.txt", dtype=int, max_rows=ns)

    # Read full dataset for acar and area from ASCII files
    acar_all = np.loadtxt("temp/Pfaf_acar.txt", dtype=float, max_rows=nc)
    area_all = np.loadtxt("output/Pfaf_area.txt", dtype=float, max_rows=nc)

    # Initialize arrays to store the selected values for each catchment
    acar = np.empty(ns, dtype=float)
    area = np.empty(ns, dtype=float)

    # Loop over each catchment index and assign values based on catid
    for i in range(ns):
        cid = catid[i]
        if cid != -9999:
            # Subtract 1 from cid to convert 1-based index (from ASCII file) to 0-based index for Python
            acar[i] = acar_all[cid - 1]
            area[i] = area_all[cid - 1]
        else:
            acar[i] = -9999.0
            area[i] = -9999.0

    # Write the output arrays to ASCII files, one number per line
    #np.savetxt("data/catch_aca_model_PR.txt", acar, fmt="%.6f")
    #np.savetxt("data/catch_area_model_PR.txt", area, fmt="%.6f")

#----look for uncorrect station------
    thres = 5000.0
    # Read data from ASCII files
    grand = np.loadtxt("input/catch_aca_grand.txt", dtype=float, max_rows=ns)

    # Initialize lists to store error information
    id_error = []

    # Loop over each catchment index
    for i in range(ns):
        if grand[i] > thres:
            if acar[i] < 0.8 * grand[i]:
                # Append error information; add 1 to i for 1-based indexing
                id_error.append(i + 1)

    # Get the number of errors found
    ne = len(id_error)

    # Write the error IDs to an ASCII file, one number per line
    #np.savetxt("data/id_error_aca_grand5000_PR.txt", np.array(id_error), fmt='%d')   

#----get corrected catid for above station--------------------

    # Read error arrays and the full catid array from ASCII files.
    # It is assumed that the files contain one number per line.
    catid_error = np.loadtxt("input/newcatid_error_aca_grand5000.txt", dtype=int, max_rows=ne)

    # Loop over each error index and update catid_all.
    # Note: We subtract 1 from resid_error values to convert from 1-based to 0-based indexing.
    for i in range(ne):
        rid = id_error[i]
        catid[rid - 1] = catid_error[i]

    # Write the updated catid_all array to an ASCII file.
    #np.savetxt("data/catid_dam_corr_aca_grand5000_noman_PR.txt", catid, fmt='%d') 

#----get dam drainage area after correction--------------------
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

    # Write the output arrays to ASCII files with one number per line
    #np.savetxt("data/catch_aca_model_corr_aca_grand5000_PR.txt", acar, fmt="%.6f")
    #np.savetxt("data/catch_area_model_corr_aca_grand5000_PR.txt", area, fmt="%.6f")

#----look for uncorrect station------
    # Define threshold and total number of catchments
    thres = 5000.0

    # The following files are read for completeness, though not used in the logic below.
    #area = np.loadtxt("data/catch_area_model_corr_aca_grand5000.txt", dtype=float, max_rows=ns)
    model = acar
    # Instead of preallocating an array with a fixed size (np in NCL),
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
    
    # Write the error indices to an ASCII file (one number per line)
    #np.savetxt("data/id_error_aca_model5000_PR.txt", np.array(id_error), fmt='%d')

#----create flag for stations------

    # Define manual adjustment arrays (1D arrays)
    resid_man = np.array([5179, 289, 7070], dtype=int)
    catid_man = np.array([46616, 142851, 199281], dtype=int)
    nman = resid_man.size

    # Update specific indices in catid_all with manual adjustments.
    # Convert from 1-based indexing (in NCL) to 0-based indexing (in Python)
    for i in range(nman):
        catid[resid_man[i] - 1] = catid_man[i]

    # Write the updated catid_all to an ASCII file
    np.savetxt("output/catid_dam_corr_aca_grand5000.txt", catid, fmt='%d')

    # Read flag_error and id_error arrays from ASCII files
    flag_error = np.loadtxt("input/flag_model5000.txt", dtype=int, max_rows=ne)
    #id_error = np.loadtxt("data/id_error_aca_model5000.txt", dtype=int, max_rows=ne)

    # Initialize flag_all array with ones (default flag value)
    flag_all = np.ones(ns, dtype=int)

    # For each error entry, update flag_all at the specified index.
    # Adjust id from 1-based to 0-based indexing.
    for i in range(ne):
        id_val = id_error[i]
        flag_all[id_val - 1] = flag_error[i]

    # Read the catchment area data from ASCII files
    #aca_grand = np.loadtxt("data/catch_aca_grand.txt", dtype=float, max_rows=ns)
    #aca_model = np.loadtxt("data/catch_aca_model_corr_aca_grand5000.txt", dtype=float, max_rows=ns)

    # Update flag_all based on conditions related to aca_grand and aca_model.
    # If aca_grand is less than 1.e3 and also less than 0.5 times aca_model, set flag to 0.
    for i in range(ns):
        if grand[i] < 1.e3:
            if grand[i] < 0.5 * acar[i]:
                flag_all[i] = 0

    # If aca_model is negative, set flag_all to 0 for that catchment.
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

    # Read the main use data as strings from the ASCII file (one entry per line)
    with open("input/main_use_grand.txt", "r") as f:
        use = [line.strip() for line in f]
    if len(use) != ns:
        print(f"Warning: expected {ns} lines, but got {len(use)} lines.")

    # For each category in use_string, create a flag array and output the result
    for j in range(nu):
        # Initialize the flag array with zeros
        flag = np.zeros(ns, dtype=int)
        # Set flag to 1 where the use value matches the current category
        for i in range(ns):
            if use[i] == use_string[j]:
                flag[i] = 1

        # Print the sum of the flag array (i.e., the count of matched entries)
       # print(np.sum(flag))
        
        # Write the flag array to an ASCII file, one number per line
        out_filename = os.path.join("output", use_out[j] + "_grand.txt")
        np.savetxt(out_filename, flag, fmt='%d')    
#----flood use--------------------
    # Read the use_irr strings from the ASCII file (one entry per line)
    with open("input/flood_use_grand.txt", "r") as f:
        use_irr = [line.strip() for line in f]

    # Initialize the flag array with zeros
    flag = np.zeros(ns, dtype=int)

    # Loop over each entry and set flag to 1 if the entry is not "NA"
    for i in range(ns):
        if use_irr[i] != "NA":
            flag[i] = 1

    # Print the sum of the flag array (i.e., count of non-"NA" entries)
    #print(np.sum(flag))

    # Write the flag array to an ASCII file, one number per line
    np.savetxt("output/fldmainsec_grand.txt", flag, fmt='%d')

#----other use--------------------
    use_out = "other"

    # Read the main use data from the ASCII file (assumed one entry per line)
    with open("input/main_use_grand.txt", "r") as f:
        use = [line.strip() for line in f]
    if len(use) != ns:
        print(f"Warning: expected {ns} entries, but got {len(use)} entries.")

    # Initialize the flag array with zeros
    flag = np.zeros(ns, dtype=int)

    # Loop over each entry and set flag to 1 if the entry matches the specified categories
    for i in range(ns):
        if use[i] == "Fisheries" or use[i] == "NA" or use[i] == "Other":
            flag[i] = 1

    # Print the sum of the flag array (i.e., count of matching entries)
    #print(np.sum(flag))

    # Write the flag array to an ASCII file, one number per line
    np.savetxt("output/" + use_out + "_grand.txt", flag, fmt='%d')

#----flood threshold-------------
if 1 == 0:

    thres_per = 1.0

    nday = 1827
    day_start = 276
    day_end = 2102

    # List files matching the pattern (assumed to be sorted in the same order as ls)
    files = sorted(glob.glob("/Volumes/PASSPORT5T/work/river/river_OL7000_Kv/*Qr*.txt"))
    nd = len(files)
    #print(nd)
    
    # Initialize data_ori array with shape (nc, nday)
    data_ori = np.empty((nc, nday), dtype=float)
    
    # Loop over each day from day_start to day_end (inclusive)
    for i in range(day_start, day_end + 1):
        if i % 10 == 0:
            print(i)
        # Read nc float values from file corresponding to day i and assign to the proper column
        data_ori[:, i - day_start] = np.loadtxt(files[i], dtype=float, max_rows=nc)
    
    # Sort each row (daily values for each grid cell) in descending order
    # Note: Sorting along axis=1 (the day dimension)
    data_sorted = np.sort(data_ori, axis=1)[:, ::-1]
    
    # Calculate threshold index based on thres_per percentage of days
    # For example, with thres_per=1, idx_thres becomes int(1/100 * 1827) = 18
    idx_thres = int(thres_per / 100.0 * nday)
    #print(idx_thres)
    
    # For each grid cell (each row), select the value at rank idx_thres-1 (i.e. the 18th highest value)
    output_data = data_sorted[:, idx_thres - 1]
    
    # Construct output filename and write the output_data to an ASCII file
    filename = "output/Pfaf_flood_qr_thres0" + str(int(thres_per)) + ".txt"
    np.savetxt(filename, output_data, fmt="%.6f")
