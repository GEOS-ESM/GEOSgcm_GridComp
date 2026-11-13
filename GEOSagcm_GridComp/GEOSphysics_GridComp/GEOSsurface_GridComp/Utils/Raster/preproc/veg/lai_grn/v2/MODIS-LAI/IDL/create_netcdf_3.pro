ndays = 46
mmdd_dates = strarr(ndays)
t =0
for dofm =0,365,8 do begin  
;    lai = fltarr (i_sib,j_sib)
 ;   lai(*,*) =0.
 ;   lai_count = fltarr (i_sib,j_sib)
 ;   lai_count(*,*) = 0. 
 
    mytime=timegen(start=julday(1,1,2006),final=julday(1,1,2006)+dofm,unit='d')
    caldat,mytime,month,day,yyyy
    ii=n_elements(mytime)-1
    date_stamp=string(yyyy[ii],month[ii],day[ii],format='(i4.4,".",i2.2,".",i2.2)')
    MMDD= string(month[ii],day[ii],format='(i2.2,i2.2)')   ;format='(i2.2,".",i2.2)')   
;    print,MMDD
    mmdd_dates (t) = MMDD
    t = t + 1
endfor
ncols = 1200
nrows = 1200
undef = -9999
scale_factor = 0.01
xgrid = fltarr(ncols,nrows)
ygrid = fltarr(ncols,nrows)
raster_pi = 3.1415926535898
nc=43200l
nr=21600l
nx = nc/36
ny = nr/18
dxy = 360./nc
x25 = fltarr(nc)
y25 = fltarr(nr)
for i = 0l,nc -1l do x25(i) = i*0.5/60. -180.
for i = 0l,nr -1l do y25(i) = i*0.5/60. -90.

units='-'
;stop
; Loop through MODIS data files
for  j_modis = 8,9 do begin;0,17 do begin
    for i_modis = 0,35 do begin;0,35 do begin
        gridcell = string(i_modis,j_modis,format='("h",i2.2,"v",i2.2)')
        ; check availabiliy of MODIS tile
        fname = '/l_data/model_parameters/LAI/MODIS-LAI/IDL/tmp_new/'+gridcell+'.txt' 
        if file_test(fname) then begin
            n_files = 0
            ixjx = intarr (2,20)
            lai = intarr(ncols,nrows,ndays)
            lai (*,*,*) = undef  
            real_lai = intarr(ncols,nrows,ndays)
            real_lai (*,*,*) = undef  
            UL_x = 0.
            UL_y = 0.
            LR_x = 0.
            LR_y = 0.

            ; reading tile geography etc

            openr,1,fname
            readf,1,nfiles
            char_var ='                '
            readf,1,char_var
            close,1
            print,'Processing ..:',fname,'',strtrim(char_var,2)
            
            FileID=hdf_sd_start('/l_data/model_parameters/LAI/MODIS-LAI/'+char_var, /read)
            HDF_SD_AttrInfo, FileID, 1, name = attr_name, data = attr_data
            sds_id = hdf_sd_select(fileID, 1)
            hdf_sd_getdata, sds_id,adum
            hdf_sd_end, fileID
            cmpstr = strcompress(attr_data)
            pos1=strpos(cmpstr,"UpperLeftPointMtrs")
            pos2=strpos(cmpstr,"(",pos1)
            pos3=strpos(cmpstr,",",pos1)
            pos4=strpos(cmpstr,")",pos1)
            
            UL_x = float(strmid(cmpstr,pos2+1,pos3-pos2-1))
            UL_y = float(strmid(cmpstr,pos3+1,pos4-pos3-1))
            
;           UL(i_modis,j_modis,0) = UL_x
;           UL(i_modis,j_modis,1) = UL_y
            
            print, strmid(cmpstr,pos2+1,pos4-pos2-1)
            print, strmid(cmpstr,pos2+1,pos3-pos2-1),' ', strmid(cmpstr,pos3+1,pos4-pos3-1)
            
            pos1=strpos(cmpstr,"LowerRightMtrs")
            pos2=strpos(cmpstr,"(",pos1)
            pos3=strpos(cmpstr,",",pos1)
            pos4=strpos(cmpstr,")",pos1)
            
            LR_x = float(strmid(cmpstr,pos2+1,pos3-pos2-1))
            LR_y = float(strmid(cmpstr,pos3+1,pos4-pos3-1))
            
;                LR(i_modis,j_modis,0) = LR_x
;                LR(i_modis,j_modis,1) = LR_y
            
            print, strmid(cmpstr,pos2+1,pos4-pos2-1)
            print, strmid(cmpstr,pos2+1,pos3-pos2-1),' ', strmid(cmpstr,pos3+1,pos4-pos3-1)
            
            pos1=strpos(cmpstr,"ProjParams")
            pos2=strpos(cmpstr,"(",pos1)
            pos3=strpos(cmpstr,")",pos1)
            
            print, strmid(cmpstr,pos2+1,pos3-pos2-1)
            project_param=fltarr(13)
            
            spos2=pos2
            for n =0,N_elements(project_param) -1 do begin
                spos3=strpos(cmpstr,",",spos2)
                project_param(n) = float(strmid(cmpstr,spos2+1,spos3-spos2-1))
                spos2=strpos(cmpstr,",",spos3)
            endfor
            
            dx = (LR_x - UL_x)/ncols
            dy = (UL_y - LR_y)/nrows
            
            x1 = UL_x + dx/2.
            y1 = UL_y - dy/2.

            tstep = 0

            ; loop through 8-day climatological data

            for dofm =0,365,8 do begin  

                mytime=timegen(start=julday(1,1,2006),final=julday(1,1,2006)+dofm,unit='d')
                caldat,mytime,month,day,yyyy
                ii=n_elements(mytime)-1
                date_stamp=string(yyyy[ii],month[ii],day[ii],format='(i4.4,".",i2.2,".",i2.2)')
                MMDD= string(month[ii],day[ii],format='(i2.2,".",i2.2)')   
                filename = '/l_data/model_parameters/LAI/MODIS-LAI/YYYY.'+ MMDD +'/MOD15A2.' + gridcell+'_new.idl'   
                print,filename
                restore,filename

                for j = 0,nrows - 1 do begin
                    for i = 0,ncols - 1 do begin
                        if((clim_lai (i,j) ge 0.) and (clim_lai (i,j) le 10.)) then begin
                            lai (i,j,tstep) = fix(clim_lai(i,j)/scale_factor)
                        endif
                    endfor
                endfor


                ; horizontal gap filling
                lai_save = lai (*,*,tstep)
                for j = 0,nrows - 1 do begin
                    for i = 0,ncols - 1 do begin
                        if(lai (i,j,tstep) eq undef) then begin

                            if(adum(i,j) eq 253) then lai (i,j,tstep) = 1 ;253 = land cover assigned as barren, sparse vegetation (rock, tundra, desert.)
                            if(adum(i,j) eq 250) then lai (i,j,tstep) = 50 ;250 = land cover assigned as urban/built-up.
                            if(adum(i,j) eq 251)  then begin ; marshland 
                     ;       if((adum(i,j) eq 249) or (adum(i,j) eq 251) or  (adum(i,j) eq 252) or (adum(i,j) eq 255)) then begin ; unclassified, marshland, snow/ice, nodata)
                                if(max (lai_save ge 0)) then begin ; there must be some data within the MODIS tile to use aas neighbor
                                    l = 0
                                    subset2=undef
                                    
                                    while (max(subset2) lt 0)  do begin
                                        
                                        imx=i+l
                                        imn=i-l
                                        jmn=j-l
                                        jmx=j+l
                                        imn=MAX([imn,0])
                                        jmn=MAX([jmn,0])
                                        imx=MIN([imx,ncols-1])
                                        jmx=MIN([jmx,nrows-1])
                                       ; if((j_modis ge 2) and (j_modis le 13)) then begin 
                                           subset2= lai_save (imn:imx,jmn:jmx)              
                                       ; endif else begin
                                       ;    subset2= lai (imn:imx,jmn:jmx,tstep)
                                       ; endelse
                                        l = l + 1
                                        
                                    endwhile
                                    subset2 (where (subset2 eq undef)) = !VALUES.F_NAN
                                    lai (i,j,tstep) = min(subset2,/nan)
                                endif
                            endif

                        endif
                    endfor
                endfor

                tstep = tstep + 1
            endfor    

            ; loop through MODIS data pixels
            for j = 0,nrows - 1 do begin
                for i = 0,ncols - 1 do begin
                    xc = x1 + i*dx
                    yc = y1 - j*dy
                    ygrid(i,j) = yc/project_param(0)
                    xgrid(i,j) = xc/project_param(0)/cos(ygrid(i,j))
                    ygrid(i,j) = ygrid(i,j)*180./raster_pi
                    xgrid(i,j) = xgrid(i,j)*180./raster_pi
                    ii = where(xgrid(i,j) ge x25)
                    jj = where(ygrid(i,j) ge y25)
                    
                    if(max(lai (i,j,*)) gt 0.) then begin
                        ;HH and VV  output file
                        IX = floor(float(n_elements(ii))/float(ncols)) + 1
                        JX = floor(float(n_elements(jj))/float(nrows)) + 1
                        i1 = (ix -1)*1200l + 1l
                        j1 = (jx -1)*1200l + 1l

                        ; i,j indices of 0.5-minute netCDF data pixel
                        ik = n_elements(ii) -i1 + 1
                        jk = n_elements(jj) -j1 + 1

                        true = 0
                                ; n_files = 0 at the beginning of the
                                ; source data tile
                        if (n_files eq 20) then begin
                            print, 'increase n_files stopping'
                            stop
                        endif

                        if (n_files gt 0) then begin
                            for f = 0, n_files -1 do begin   
                                ; already have output HH VV file created
                                if((ixjx(0,f) eq ix) and (ixjx(1,f) eq jx)) then true = 1
                            endfor
                        endif

                        if (true eq 0) then begin
                            nc_filename = '../MODIS_8-DayClim/netCDF/MODIS_lai_clim.'+ 'H'+string(ix,'(i2.2)')+'V'+string(jx,'(i2.2)') + '.nc'
                            if (file_test(nc_filename)) then begin
                                print, 'Shouldnt we have visted here?',true ,n_files,'H'+string(ix,'(i2.2)')+'V'+string(jx,'(i2.2)') + '.nc'
 ;                               stop
                                ixjx(0,n_files) = ix
                                ixjx(1,n_files) = jx

                                id   = NCDF_OPEN(nc_filename)
                                if (n_files eq 0) then begin 
                                    NCDF_VARGET, id,'LAI',save_lai1
                                    save_lai1 (ik,jk,*) = lai (i,j,*)
                                endif
                                if (n_files eq 1) then begin 
                                    NCDF_VARGET, id,'LAI',save_lai2
                                    save_lai2 (ik,jk,*) = lai (i,j,*)
                                endif
                                
                                if (n_files eq 2) then begin 
                                    NCDF_VARGET, id,'LAI',save_lai3
                                    save_lai3 (ik,jk,*) = lai (i,j,*)
                                endif
                                
                                if (n_files eq 3) then begin 
                                    NCDF_VARGET, id,'LAI',save_lai4
                                    save_lai4 (ik,jk,*) = lai (i,j,*)
                                endif
                                
                                if (n_files eq 4) then begin 
                                    NCDF_VARGET, id,'LAI',save_lai5
                                    save_lai5 (ik,jk,*) = lai (i,j,*)
                                endif
                                
                                if (n_files eq 5) then begin 
                                    NCDF_VARGET, id,'LAI',save_lai6
                                    save_lai6 (ik,jk,*) = lai (i,j,*)
                                endif
                                 if (n_files eq 6) then begin 
                                    NCDF_VARGET, id,'LAI',save_lai7
                                    save_lai7 (ik,jk,*) = lai (i,j,*)
                                endif
                                ;           
                                if (n_files eq 7) then begin 
                                    NCDF_VARGET, id,'LAI',save_lai8
                                    save_lai8 (ik,jk,*) = lai (i,j,*)
                                endif
                                ;           
                                if (n_files eq 8) then begin 
                                    NCDF_VARGET, id,'LAI',save_lai9
                                    save_lai9 (ik,jk,*) = lai (i,j,*)
                                endif
                                ;                                          ;           
                                NCDF_CLOSE, id

                                n_files = n_files + 1

                            endif else begin
                                ixjx(0,n_files) = ix
                                ixjx(1,n_files) = jx
                                long_data = indgen(nx)*dxy + i1*dxy - dxy/2. -180.
                                lat_data  = indgen(ny)*dxy + j1*dxy - dxy/2. -90.
                                
                                write_netcdf_file,filename = nc_filename, $
                                  nc_global = nc, nr_global = nr, nc_local = nx, nr_local = ny, i_offset = i1, j_offset = j1, $
                                  long_data= long_data, lat_data = lat_data, n_timesteps = ndays, timesteps = mmdd_dates,cellsize=dxy*3600.,  $
                                  scale_factor = [scale_factor], units = [units], undef=undef, varname = ['LAI'], data1 = real_lai
                                
                                if (n_files eq 0) then begin 
                                    save_lai1 = intarr(ncols,nrows,ndays)
                                    save_lai1 (*,*,*)   = undef  
                                    save_lai1 (ik,jk,*) = lai (i,j,*)
                                endif
                                
                                if (n_files eq 1) then begin 
                                    save_lai2 = intarr(ncols,nrows,ndays)
                                    save_lai2 (*,*,*)   = undef  
                                    save_lai2 (ik,jk,*) = lai (i,j,*)
                                endif
                                
                                if (n_files eq 2) then begin 
                                    save_lai3 = intarr(ncols,nrows,ndays)
                                    save_lai3 (*,*,*)   = undef  
                                    save_lai3 (ik,jk,*) = lai (i,j,*)
                                endif
                                
                                if (n_files eq 3) then begin 
                                    save_lai4 = intarr(ncols,nrows,ndays)
                                    save_lai4 (*,*,*)   = undef  
                                    save_lai4 (ik,jk,*) = lai (i,j,*)
                                endif
                                
                                if (n_files eq 4) then begin 
                                    save_lai5 = intarr(ncols,nrows,ndays)
                                    save_lai5 (*,*,*)   = undef  
                                    save_lai5 (ik,jk,*) = lai (i,j,*)
                                endif
                                
                                if (n_files eq 5) then begin 
                                    save_lai6 = intarr(ncols,nrows,ndays)
                                    save_lai6 (*,*,*)   = undef  
                                    save_lai6 (ik,jk,*) = lai (i,j,*)
                                endif                              
 
                                if (n_files eq 6) then begin 
                                    save_lai7 = intarr(ncols,nrows,ndays)
                                    save_lai7 (*,*,*)   = undef  
                                    save_lai7 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 7) then begin 
                                    save_lai8 = intarr(ncols,nrows,ndays)
                                    save_lai8 (*,*,*)   = undef  
                                    save_lai8 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 8) then begin 
                                    save_lai9 = intarr(ncols,nrows,ndays)
                                    save_lai9 (*,*,*)   = undef  
                                    save_lai9 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 9) then begin 
                                    save_lai10 = intarr(ncols,nrows,ndays)
                                    save_lai10 (*,*,*)   = undef  
                                    save_lai10 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 10) then begin 
                                    save_lai11 = intarr(ncols,nrows,ndays)
                                    save_lai11 (*,*,*)   = undef  
                                    save_lai11 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 11) then begin 
                                    save_lai12 = intarr(ncols,nrows,ndays)
                                    save_lai12 (*,*,*)   = undef  
                                    save_lai12 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 12) then begin 
                                    save_lai13 = intarr(ncols,nrows,ndays)
                                    save_lai13 (*,*,*)   = undef  
                                    save_lai13 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 13) then begin 
                                    save_lai14 = intarr(ncols,nrows,ndays)
                                    save_lai14 (*,*,*)   = undef  
                                    save_lai14 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 14) then begin 
                                    save_lai15 = intarr(ncols,nrows,ndays)
                                    save_lai15 (*,*,*)   = undef  
                                    save_lai15 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 15) then begin 
                                    save_lai16 = intarr(ncols,nrows,ndays)
                                    save_lai16 (*,*,*)   = undef  
                                    save_lai16 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 16) then begin 
                                    save_lai17 = intarr(ncols,nrows,ndays)
                                    save_lai17 (*,*,*)   = undef  
                                    save_lai17 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 17) then begin 
                                    save_lai18 = intarr(ncols,nrows,ndays)
                                    save_lai18 (*,*,*)   = undef  
                                    save_lai18 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 18) then begin 
                                    save_lai19 = intarr(ncols,nrows,ndays)
                                    save_lai19 (*,*,*)   = undef  
                                    save_lai19 (ik,jk,*) = lai (i,j,*)
                                endif 

                                if (n_files eq 19) then begin 
                                    save_lai20 = intarr(ncols,nrows,ndays)
                                    save_lai20 (*,*,*)   = undef  
                                    save_lai20 (ik,jk,*) = lai (i,j,*)
                                endif 

;                                id   = NCDF_OPEN(nc_filename, /write) 
;                              NCDF_VARGET, id,'LAI',read_lai
;                              read_lai (ik,jk,*) = lai (i,j,*)
;                              NCDF_VARPUT, id,'LAI',read_lai
;                              NCDF_CLOSE, id

                              n_files = n_files + 1

                          endelse

                      endif else begin

                          for f = 0, n_files -1 do begin   
                                ; already have output HH VV file created
                                if((ixjx(0,f) eq ix) and (ixjx(1,f) eq jx)) then nf2 = f
                            endfor
                          ;nf2 = n_files -1

                          if (nf2 eq 0) then  save_lai1 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 1) then  save_lai2 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 2) then  save_lai3 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 3) then  save_lai4 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 4) then  save_lai5 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 5) then  save_lai6 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 6) then  save_lai7 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 7) then  save_lai8 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 8) then  save_lai9 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 9) then  save_lai10 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 10) then  save_lai11 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 11) then  save_lai12 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 12) then  save_lai13 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 13) then  save_lai14 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 14) then  save_lai15 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 15) then  save_lai16 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 16) then  save_lai17 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 17) then  save_lai18 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 18) then  save_lai19 (ik,jk,*) = lai (i,j,*)
                          if (nf2 eq 19) then  save_lai20 (ik,jk,*) = lai (i,j,*)
 
                      endelse



     ;                        if((clim_lai(i,j) ge 0.) or (clim_lai(i,j) le 10.)) then begin
     ;                       lai      (n_elements(ii)-1,n_elements(jj)-1) = lai      (n_elements(ii)-1,n_elements(jj)-1) + clim_lai(i,j)
     ;                       lai_count(n_elements(ii)-1,n_elements(jj)-1) = lai_count(n_elements(ii)-1,n_elements(jj)-1) + 1.
                  endif
              endfor
          endfor
          if (n_files gt 0) then begin
              for f = 0, n_files -1 do begin
                  ix = ixjx(0,f) 
                  jx = ixjx(1,f) 
                  nc_filename = '../MODIS_8-DayClim/netCDF/MODIS_lai_clim.'+ 'H'+string(ix,'(i2.2)')+'V'+string(jx,'(i2.2)') + '.nc'

                  id   = NCDF_OPEN(nc_filename,/write)

                  if (f eq 0) then  NCDF_VARPUT, id,'LAI',save_lai1 
                  if (f eq 1) then  NCDF_VARPUT, id,'LAI',save_lai2 
                  if (f eq 2) then  NCDF_VARPUT, id,'LAI',save_lai3 
                  if (f eq 3) then  NCDF_VARPUT, id,'LAI',save_lai4 
                  if (f eq 4) then  NCDF_VARPUT, id,'LAI',save_lai5 
                  if (f eq 5) then  NCDF_VARPUT, id,'LAI',save_lai6 
                  if (f eq 6) then  NCDF_VARPUT, id,'LAI',save_lai7 
                  if (f eq 7) then  NCDF_VARPUT, id,'LAI',save_lai8 
                  if (f eq 8) then  NCDF_VARPUT, id,'LAI',save_lai9 
                  if (f eq 9) then  NCDF_VARPUT, id,'LAI',save_lai10
                  if (f eq 10) then  NCDF_VARPUT, id,'LAI',save_lai11
                  if (f eq 11) then  NCDF_VARPUT, id,'LAI',save_lai12
                  if (f eq 12) then  NCDF_VARPUT, id,'LAI',save_lai13
                  if (f eq 13) then  NCDF_VARPUT, id,'LAI',save_lai14
                  if (f eq 14) then  NCDF_VARPUT, id,'LAI',save_lai15
                  if (f eq 15) then  NCDF_VARPUT, id,'LAI',save_lai16
                  if (f eq 16) then  NCDF_VARPUT, id,'LAI',save_lai17
                  if (f eq 17) then  NCDF_VARPUT, id,'LAI',save_lai18
                  if (f eq 18) then  NCDF_VARPUT, id,'LAI',save_lai19
                  if (f eq 19) then  NCDF_VARPUT, id,'LAI',save_lai20

                  NCDF_CLOSE, id
              endfor
          endif
      endif
  endfor	
endfor	

end
