file1 = '/l_data/GEOLAND2-LAI/LAI-CLIM/netCDF/geoland2_lai_clim.H01V08.nc'
file2 = '/l_data/GEOLAND2-LAI/LAI-CLIM/netCDF2/geoland2_lai_clim.H01V08.nc'

id1 = NCDF_OPEN(file1)
id2 = NCDF_OPEN(file2)

NCDF_VARGET, id1, 'MMDD', MMDD1
NCDF_VARGET, id2, 'MMDD', MMDD2

NCDF_VARGET, id1, 'LAI', lai1
NCDF_VARGET, id2, 'LAI', lai2

NCDF_CLOSE, id1
NCDF_CLOSE, id2

end
