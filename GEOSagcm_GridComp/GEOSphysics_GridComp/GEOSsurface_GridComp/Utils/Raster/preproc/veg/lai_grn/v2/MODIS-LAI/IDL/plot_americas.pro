limit_global = [-60,-180,90,180]
limit_conus  = [29.,-125.,54.,-82.]
limit_americas = [-60,-170,90,-25]
limit_east   = [29,-84.,36.,-77]
limit_east2   = [30,-90.,40.,-80] ; 
limit_east3   = [35.,-82.,42.,-73] ; MD area
limit_scan   = [55.,0.,75.,55.]
limits = limit_east3 ;limit_global
file_dom ='_conu2'
DEVICE, DECOMPOSED = 0
 loadct,8
TVLCT, newRed, newGreen, newBlue,/GET
newRed   = reverse (newRed)
newGreen = reverse (newGreen)
newBlue  = reverse (newBlue)

newRed(0) = 0
newGreen(0) = 0
newBlue(0) = 0
newRed(255) = 255
newGreen(255) = 255
newBlue(255) = 255
TVLCT, newRed, newGreen, newBlue

;loadct2,0
;colors = indgen(21) + 80
;n_levels = 21
;mg_loadct, 0, file='/home/sarith/IDL/lib/brewer/brewer.tbl',rgb_table=ct	
;newRed = ct(*,0)
;newGreen = ct(*,1)
;newBlue = ct(*,2)
;newRed(0) = 0
;newGreen(0) = 0
;newBlue(0) = 0
;newRed(255) = 255
;newGreen(255) = 255
;newBlue(255) = 255
;TVLCT, newRed, newGreen, newBlue

n_levels = 25
upval = 6.
lwval = 0.
levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]

colors = indgen(25)*10 + 5
;goto,jump
filename = '/l_data/model_parameters/LAI/MODIS-LAI/MODIS_8-DayClim/netCDF-with-gaps/MODIS_lai_clim.H33V05.nc'
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


;n_levels = 256
;mg_loadct, 0, file='/home/sarith/IDL/lib/brewer/brewer.tbl',rgb_table=ct
;      RED   = CT(0:N_levels-1,0)
;      GREEN = CT(0:N_levels-1,1)
;      BLUE  = CT(0:N_levels-1,2)
;
;colors = indgen(n_levels)


jump:

for t = 19,19 do begin ;0,n_elements(mmdd) - 1 do begin	
;	goto, jump2
   global_var1 = fltarr (N_lon_plot,N_lat_plot)
   global_var1(*,*) = !VALUES.F_NAN

   for jx =1,18 do begin
      for ix =1,36 do begin
         filename = '/l_data/model_parameters/LAI/MODIS-LAI/MODIS_8-DayClim/netCDF-with-gaps/MODIS_lai_clim.H'+string(ix,'(i2.2)')+'V'+string(jx,'(i2.2)') + '.nc'
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
	jump2:
	
      !P.position=0
      !P.Multi = [0, 1, 1, 0, 1]
      Erase, Color= 255
      MAP_SET,/ROBINSON,/hires,color= 0,/NoErase,limit=limits,title =mmdd(t),/isotropic
      contour, global_var1,x,y,levels = levels,c_colors=colors,/overplot;,/cell_fill
   alpha=fltarr(n_levels,2)
   alpha(*,0)=levels
   alpha(*,1)=levels
   h=[0,1]
   !P.position=[0.30, 0.0+0.005, 0.70, 0.015+0.005]
   clev = levels
   clev (*) = 1
  contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xticks=10
  contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(levels) -1 do xyouts,levels(k),1.1, string(levels(k),'(f4.2)'),orientation=90,color=0

thisDevice = !D.NAME
SET_PLOT, 'Z'
SET_PLOT, thisDevice

Device, Get_Decomposed=theDecomposedState, Get_Visual_Depth=theDepth
Device, Decomposed=1
color24 = 1

image = TVRD(true=color24)
TVLCT, red, green, blue,/get
image3d = image

write_jpeg,'MODIS_GLOBAL_' + string(t+1,'(i2.2)') + '.jpg',image3d,quality=100,true=1

;TVLCT, red, green, blue, /GET
;thisDevice = !D.NAME
;SET_PLOT, 'Z'
;SET_PLOT, thisDevice
;
;image = TVRD()
;thisImage = BYTSCL(image)
;s = SIZE(thisImage)
;image3d = BYTARR(s(1), s(2),3)
;red = reverse(red)
;green = reverse(green)
;blue = reverse(blue) 
;image3d( *, *,0) = red(thisImage)
;image3d( *, *,1) = green(thisImage)
;image3d( *, *,2) = blue(thisImage)
;write_jpeg,'MODIS_GLOBAL_' + string(t+1,'(i2.2)') + '.jpg',image3d,quality=100,true=3
;
;image = TVRD()
;thisImage = BYTSCL(image)
;s = SIZE(thisImage)
;image3d = BYTARR(s(1), s(2),3)
;image3d( *, *,0) = red(thisImage)
;image3d( *, *,1) = green(thisImage)
;image3d( *, *,2) = blue(thisImage)
;;WRITE_JPEG,strcompress(Region,/remove_all) + '.jpg', TVRD(), QUALITY=75
;write_jpeg,'MODIS_GLOBAL_' + string(t+1,'(i2.2)') + '.jpg',image3d,quality=100,true=3

;wdelete
endfor

end
