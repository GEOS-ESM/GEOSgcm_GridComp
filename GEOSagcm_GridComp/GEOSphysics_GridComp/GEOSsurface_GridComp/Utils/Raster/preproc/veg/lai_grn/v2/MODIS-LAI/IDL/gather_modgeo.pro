;pro gather_netcdf,layer=layer,global_var1=global_var1,global_var2=global_var2,global_var3=global_var3,dir_name = dir_name
limit_global = [-60,-180,90,180]
limit_conus  = [29.,-125.,54.,-82.]
limit_americas = [-60,-180,90,-25]
limits = limit_global
file_dom ='_conus'
;goto,jump
filename = '/l_data/model_parameters/LAI/MODIS-LAI/MODIS_8-DayClim/netCDF/MODIS_lai_clim.H33V05.nc'
plot_int  = 1
ct_number = 25

ncid = NCDF_OPEN(filename)
ncdf_attget,ncid,'N_lon_global',nc_global,/global
ncdf_attget,ncid,'N_lat_global',nr_global,/global
ncdf_attget,ncid,'i_ind_offset_LL',iLL,/global
ncdf_attget,ncid,'j_ind_offset_LL',jLL,/global
ncdf_varget,ncid,'MMDD',mmdd

mmdd = string (mmdd)

N_lon = nc_global
N_lat = nr_global

dx    = 360./N_lon
dy    = 180./N_lat

x_high = fltarr(N_lon) 
y_high = fltarr(N_lat) 

for i = 0l,N_lon -1l do x_high(i) = i*dx -180. +  dx/2.
for i = 0l,N_lat -1l do y_high(i) = i*dx -90.  +  dy/2.

N_lon_plot = N_lon/(plot_int)
N_lat_plot = N_lat/(plot_int)

ii = indgen(N_lon_plot)*(long(plot_int)) + plot_int/2
jj = indgen(N_lat_plot)*(plot_int) + plot_int/2

x = x_high([ii])
y = y_high([jj])

if(plot_int eq 1) then begin
   
   x = x_high
   y = y_high

endif

vname1 = 'LAI'

ncdf_attget,ncid,vname1,'units',units
units=string(units)

ncdf_attget,ncid,vname1,'ScaleFactor',sf_lai
ncdf_attget,ncid,vname1,'UNDEF',undef
ncdf_varget,ncid,vname1,var

nc = n_elements(var(*,0))
nr = n_elements(var(0,*))

NCDF_CLOSE, ncid 


iLL = 0l
jLL = 0l


n_levels = 21
loadct2,0
colors = indgen(n_levels)+ 80
upval = 10.
lwval = 0.
levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]

for t = 0,n_elements(mmdd) - 1 do begin	

   global_var1 = fltarr (N_lon_plot,N_lat_plot)
   global_var1(*,*) = !VALUES.F_NAN

   for jx =1,18 do begin
      for ix =1,36 do begin
         filename = '/l_data/model_parameters/LAI/MODIS-LAI/MODIS_8-DayClim/netCDF/MODIS_lai_clim.H'+string(ix,'(i2.2)')+'V'+string(jx,'(i2.2)') + '.nc'
         if file_test (filename) then begin
            ncid = NCDF_OPEN(filename)
            ncdf_attget,ncid,'i_ind_offset_LL',iLL,/global
            ncdf_attget,ncid,'j_ind_offset_LL',jLL,/global
            ncdf_varget,ncid,vname1,var1
            ncdf_attget,ncid,vname1,'ScaleFactor',sf
            ncdf_attget,ncid,vname1,'UNDEF',undef
            NCDF_CLOSE, ncid 
            for j = long(jLL-1)+ plot_int/2,long(jLL) - 1l + nr -1, plot_int do begin
               for i = long(iLL-1)+ plot_int/2, long(iLL-1l) + nc -1, plot_int do begin                  
                  if(var1(i -(iLL-1l) ,j -(jLL-1l),t) ne undef) then global_var1(i/plot_int,j/plot_int) =$
                     sf_lai*var1(i -(iLL-1) ,j -(jLL-1),t)
               endfor	               
            endfor
         endif
      endfor
   endfor
      !P.position=0
      !P.Multi = [0, 1, 1, 0, 1]
      Erase, Color= 255
      MAP_SET,/ROBINSON,/hires,color= 0,/NoErase,limit=limits,title =mmdd(t)
      contour, global_var1,x,y,levels = levels,c_colors=colors,/overplot
      write_to_jpeg,file_name= 'MODIS_' + string (t+1,'(i2.2)') + '.jpg'
endfor

end
