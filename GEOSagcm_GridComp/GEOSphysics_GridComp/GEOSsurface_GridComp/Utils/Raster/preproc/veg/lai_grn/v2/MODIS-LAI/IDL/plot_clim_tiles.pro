n_levels = 21
loadct2,0
colors = indgen(n_levels)+ 80
upval = 7.
lwval = 0.
levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]
      !P.position=0
      !P.Multi = [0, 2, 1, 0, 1]
;restore,'/l_data/model_parameters/LAI/MODIS-LAI/YYYY.06.02/MOD15A2.h18v02_new.idl'
restore,'/l_data/model_parameters/LAI/MODIS-LAI/YYYY.06.02/MOD15A2.h11v05_new.idl'
 Erase, Color= 255
clim_lai = reverse(clim_lai,2,/overwrite)
contour, clim_lai,levels = levels,c_colors=colors,/isotropic

jx = 13
ix = 10

 filename = '/l_data/model_parameters/LAI/MODIS-LAI/MODIS_8-DayClim/netCDF/MODIS_lai_clim.H'+string(ix,'(i2.2)')+'V'+string(jx,'(i2.2)') + '.nc'
ncid = NCDF_OPEN(filename)
            ncdf_attget,ncid,'i_ind_offset_LL',iLL,/global
            ncdf_attget,ncid,'j_ind_offset_LL',jLL,/global
            ncdf_varget,ncid,'LAI',var1
            ncdf_attget,ncid,'LAI','ScaleFactor',sf
            ncdf_attget,ncid,'LAI','UNDEF',undef
            NCDF_CLOSE, ncid 
var1 (where (var1 eq undef)) = !VALUES.F_NAN
var1 = var1*sf

contour,var1(*,*,19),levels = levels,c_colors=colors,/isotropic

write_to_jpeg,file_name= 'MODIS_CLIM_TILE.jpg'

end
