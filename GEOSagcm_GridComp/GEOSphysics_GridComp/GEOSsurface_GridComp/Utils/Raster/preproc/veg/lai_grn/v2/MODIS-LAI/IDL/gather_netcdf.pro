filename = '../MODIS_8-DayClim/netCDF/MODIS_lai_clim.H33V16.nc'
vname = 'LAI'
 
ncid = NCDF_OPEN(filename)
ncdf_attget,ncid,'N_lon_global',nc_global,/global
ncdf_attget,ncid,'N_lat_global',nr_global,/global
ncdf_attget,ncid,'i_ind_offset_LL',iLL,/global
ncdf_attget,ncid,'j_ind_offset_LL',jLL,/global

ncdf_attget,ncid,vname,'units',units
units=string(units)
ncdf_attget,ncid,vname,'ScaleFactor',sf
ncdf_attget,ncid,vname,'UNDEF',undef
ncdf_varget,ncid,vname,var
ncdf_varget,ncid,'MMDD',str
MMDD = string(str)

nc = 1200
nr = 1200 
NCDF_CLOSE, ncid 

global_var = fltarr(nc_global,nr_global)
global_var1 = fltarr(nc_global,nr_global)
global_var2 = fltarr(nc_global,nr_global)

global_var (*,*) = !VALUES.F_NAN
global_var1(*,*) = !VALUES.F_NAN
global_var2(*,*) = !VALUES.F_NAN
iLL = 0l
jLL = 0l
tstep = 32
for jx =1,18 do begin
    for ix =1,36 do begin
        filename = '/l_data/model_parameters/LAI/MODIS-LAI/MODIS_8-DayClim/netCDF/MODIS_lai_clim.H'+string(ix,'(i2.2)')+'V'+string(jx,'(i2.2)') + '.nc'
        if file_test (filename) then begin
            print, 'Processing :','MODIS_lai_clim.H'+string(ix,'(i2.2)')+'V'+string(jx,'(i2.2)') + '.nc'
            ncid = NCDF_OPEN(filename)
            ncdf_attget,ncid,'i_ind_offset_LL',iLL,/global
            ncdf_attget,ncid,'j_ind_offset_LL',jLL,/global
            ncdf_varget,ncid,vname ,var
             ncdf_attget,ncid,vname,'ScaleFactor',sf
            ncdf_attget,ncid,vname,'UNDEF',undef
            NCDF_CLOSE, ncid 
            for j = long(jLL-1),long(jLL) - 1l + nr -1 do begin
                for i = long(iLL-1), long(iLL-1l) + nc -1 do begin
                    if(var(i -(iLL-1l) ,j -(jLL-1l)) ne undef) then global_var (i,j) = sf*var (i -(iLL-1) ,j -(jLL-1),tstep)
  
                endfor
            endfor
        endif
    endfor
endfor
;global_var(where (global_var lt 1333.)) = 0

nc_low=720
nr_low=360
ncb = nc_low
nrb = nr_low
dx25=0.5
dy25=0.5
pgrid = fltarr (nc_low,nr_low)
pgrid1= fltarr (nc_low,nr_low)
pgrid2= fltarr (nc_low,nr_low)
k1 = nc_global/nc_low
k2 = k1 -1 

x1b=-180. + dx25/2.
xeb= x1b +(ncb-1)*dx25
y1b= -90. + dy25/2.
yeb= y1b +(nrb-1)*dy25

for j = 0,nr_low -1 do begin
    for i =0,nc_low-1 do begin
        pgrid  (i,j) = mean (global_var (i*long(k1):i*long(k1)+k2,j*long(k1):j*long(k1)+k2),/nan)
        pgrid1 (i,j) = mean (global_var1(i*long(k1):i*long(k1)+k2,j*long(k1):j*long(k1)+k2),/nan)
        pgrid2 (i,j) = mean (global_var2(i*long(k1):i*long(k1)+k2,j*long(k1):j*long(k1)+k2),/nan)
    endfor
endfor


sp,1,1,/f,/h,/color,landscape=1

loadct2,0
colors=90+indgen(11)
upval =10.0
lwval =0.
levp = [lwval,lwval+(upval-lwval)/10+indgen(10)*(upval-lwval)/10.]

contour,pgrid ,title=vname+' '+MMDD(tstep) ,c_colors=colors,levels=levp,xstyle=4,ystyle=4,/CELL_FILL
;contour,pgrid1,title=vname1,c_colors=colors,levels=levp,xstyle=4,ystyle=4,/CELL_FILL
;contour,pgrid2,title=vname2,c_colors=colors,levels=levp,xstyle=4,ystyle=4,/CELL_FILL

!P.position=[0.30, 0.0, 0.70, 0.015]
    alpha=fltarr(n_elements(levp),2)
    alpha(*,0)=levp
    alpha(*,1)=levp
    h=[0,1]
    contour,alpha,levp,h,levels=levp,c_colors=colors,/fill,/xstyle,/ystyle, $
      /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levp),max(levp)], $
      xtitle=' ', color=0, $
      C_charsize=1.0, charsize=1. ,xticks=10
    contour,alpha,levp,h,levels=levp,color=0,/overplot

!P.position=0
closeps,vname+'.'+MMDD(tstep)+'.ps'
end
