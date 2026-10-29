filename ='/l_data/MODIS/LAI/2011.04.23/MOD15A2.A2011113.h14v01.005.2011123123507.hdf'

fid =0l
fid=hdf_sd_start(filename, /read)
hdf_sd_fileinfo, fid, datasets, attributes
sds_id = hdf_sd_select(fid, 1)
hdf_sd_getinfo, sds_id, dims=dims, name=name, type=type, unit=unit, natts=natts
hdf_sd_getdata, sds_id,adum

HDF_SD_ENDACCESS, sds_id

save,adum,file =strtrim(filename,2)+'.idl'
end
