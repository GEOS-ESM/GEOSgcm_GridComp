data = intarr (1200,1200)
long_data = fltarr (1200)
lat_data  = fltarr (1200)
data (*,*) = 1
write_netcdf_file,filename = 'test.nc', nc_global = 43200l, nr_global = 21600l, nc_local = 1200,   $
       nr_local = 1200, n_timesteps = 2, scale_factor = 100., timesteps =['0115','0615'] ,       $
       varname = 'albedo', data = data, units = '-',undef=-9999, i_offset = 1, j_offset = 1,     $
       long_data= long_data, lat_data = lat_data
end
