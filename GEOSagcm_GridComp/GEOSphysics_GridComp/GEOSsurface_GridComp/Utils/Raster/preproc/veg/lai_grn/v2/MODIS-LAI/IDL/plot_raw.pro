limit_global = [-60,-180,90,180]
limit_conus  = [29.,-125.,54.,-82.]
limit_americas = [-60,-180,90,-25]
limits = limit_global

project =1
ncols = 1200
nrows = 1200
nc=43200l
nr=21600l
lai = fltarr(nc,nr)

raster_pi = 3.1415926535898

nx = nc/36
ny = nr/18
dxy = 360./nc
x25 = fltarr(nc)
y25 = fltarr(nr)
for i = 0l,nc -1l do x25(i) = i*0.5/60. -180.
for i = 0l,nr -1l do y25(i) = i*0.5/60. -90.

x_high = fltarr(nc)
y_high = fltarr(nr)

for i = 0l,NC -1l do x_high(i) = i*dxy -180. +  dxy/2.
for i = 0l,NR -1l do y_high(i) = i*dxy -90.  +  dxy/2.

lai(*,*) =  !VALUES.F_NAN

n_levels = 256
mg_loadct, 0, file='/home/sarith/IDL/lib/brewer/brewer.tbl',rgb_table=ct
      RED   = CT(0:N_levels-1,0)
      GREEN = CT(0:N_levels-1,1)
      BLUE  = CT(0:N_levels-1,2)

colors = indgen(n_levels)
upval = 7.
lwval = 0.
levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]
spawn, 'ls /l_data/model_parameters/LAI/MODIS-LAI/2008.06.01/MOD15A2.*hdf',flist

for k = 0,n_elements(flist) - 1 do begin
   jx=1*strmid(flist(k),71,2)
   ix=1*strmid(flist(k),68,2)

   fid =0l
   fid=hdf_sd_start(flist(k), /read)
   hdf_sd_fileinfo, fid, datasets, attributes
   sds_id = hdf_sd_select(fid, 1)
   hdf_sd_getinfo, sds_id, dims=dims, name=name, type=type, unit=unit, natts=natts
   print,flist(k),name,ix,jx
   hdf_sd_getdata, sds_id,adum

   if (project eq 1 ) then begin
         UL_x = 0.
         UL_y = 0.
         LR_x = 0.
         LR_y = 0.
         xgrid = fltarr(ncols,nrows)
         ygrid = fltarr(ncols,nrows)

         HDF_SD_AttrInfo, FID, 1, name = attr_name, data = attr_data
   ;      hdf_sd_end, fid
         cmpstr = strcompress(attr_data)
;          print,cmpstr
;          wait,5
         pos1=strpos(cmpstr,"UpperLeftPointMtrs")
         pos2=strpos(cmpstr,"(",pos1)
         pos3=strpos(cmpstr,",",pos1)
         pos4=strpos(cmpstr,")",pos1)
         
         UL_x = float(strmid(cmpstr,pos2+1,pos3-pos2-1))
         UL_y = float(strmid(cmpstr,pos3+1,pos4-pos3-1))
         
         pos1=strpos(cmpstr,"LowerRightMtrs")
         pos2=strpos(cmpstr,"(",pos1)
         pos3=strpos(cmpstr,",",pos1)
         pos4=strpos(cmpstr,")",pos1)
         
         LR_x = float(strmid(cmpstr,pos2+1,pos3-pos2-1))
         LR_y = float(strmid(cmpstr,pos3+1,pos4-pos3-1))
         
         pos1=strpos(cmpstr,"ProjParams")
         pos2=strpos(cmpstr,"(",pos1)
         pos3=strpos(cmpstr,")",pos1)
         
;         print, strmid(cmpstr,pos2+1,pos3-pos2-1)
         project_param=fltarr(13)
         
         spos2=pos2
         for n =0,N_elements(project_param) -1 do begin
            spos3=strpos(cmpstr,",",spos2)
            project_param(n) = float(strmid(cmpstr,spos2+1,spos3-spos2-1))
            spos2=strpos(cmpstr,",",spos3)
         endfor
         
         dx = (LR_x - UL_x)/1200.
         dy = (UL_y - LR_y)/1200.
         x1 = UL_x + dx/2.
         y1 = UL_y - dy/2.
;         print,x1,y1,dx,dy,LR_x, UL_x,UL_y,LR_y
;         wait,3
      endif
   
   HDF_SD_ENDACCESS, sds_id
   HDF_SD_END, fid
   for j =0,1199 do begin
      for i =0,1199 do begin
         if ((adum(i,j) ge 0) and (adum(i,j) le 100)) then begin
            if(project eq 0) then lai(ix*1200l+i,jx*1200l+j) =  0.1*adum(i,j)
            if(project eq 1.) then begin
               xc = x1 + i*dx
               yc = y1 - j*dy
               ygrid(i,j) = yc/project_param(0)
               xgrid(i,j) = xc/project_param(0)/cos(ygrid(i,j))
               ygrid(i,j) = ygrid(i,j)*180./raster_pi
               xgrid(i,j) = xgrid(i,j)*180./raster_pi
               ii = where(xgrid(i,j) ge x25)
               jj = where(ygrid(i,j) ge y25)
               IX = floor(float(n_elements(ii))/float(ncols)) + 1
               JX = floor(float(n_elements(jj))/float(nrows)) + 1
               i1 = (ix -1)*1200l + 1l
               j1 = (jx -1)*1200l + 1l
               ik = n_elements(ii)-1 ;-i1 + 1
               jk = n_elements(jj)-1 ;-j1 + 1
              
               lai(ik,jk) =  0.1*adum(i,j)
 ;              print,ik,jk,x_high(ik),y_high(jk),lai(ik,jk),ix,jx,i1,j1 
            endif
         endif
      endfor
   endfor     
endfor
      !P.position=0
      !P.Multi = [0, 1, 1, 0, 1]
      Erase, Color= 0
      if(project eq 0) then $
      contour, lai,levels = levels,c_colors=colors
      if(project eq 1) then begin
         MAP_SET,/ROBINSON,/hires,color= 255,/NoErase,limit=limits
         contour, lai,x_high,y_high,levels = levels,c_colors=colors,/overplot
      endif
;      write_to_jpeg,file_name= 'MODIS_Americas_' + string
;      (t+1,'(i2.2)') + '.jpg'
thisDevice = !D.NAME
SET_PLOT, 'Z'

TVLCT, red, green, blue
SET_PLOT, thisDevice

image = TVRD()
thisImage = BYTSCL(image)
s = SIZE(thisImage)
image3d = BYTARR(s(1), s(2),3)
image3d( *, *,0) = red(thisImage)
image3d( *, *,1) = green(thisImage)
image3d( *, *,2) = blue(thisImage)
;WRITE_JPEG,strcompress(Region,/remove_all) + '.jpg', TVRD(), QUALITY=75
write_jpeg,'MODIS_RAW_PROJ.jpg',image3d,quality=75,true=3

wdelete

end
