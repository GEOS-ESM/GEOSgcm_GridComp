pro spatial_interpolation4
;sp,1,1,/f,/h,/color,/land
;UL=fltarr(36,18,2)
;LR=fltarr(36,18,2)
;loadct2,0

ncols = 1200
nrows = 1200

xgrid = fltarr(ncols,nrows)
ygrid = fltarr(ncols,nrows)
raster_pi = 3.1415926535898
;colors=80+indgen(21)
colors=90+indgen(11)
upval =10.0
lwval =0.
;levp = [lwval,lwval+(upval-lwval)/20+indgen(20)*(upval-lwval)/20.]
levp = [lwval,lwval+(upval-lwval)/10+indgen(10)*(upval-lwval)/10.]

i_sib=8640
j_sib=4320

x25=indgen(i_sib)*2.5/60. -180.
y25=indgen(j_sib)*2.5/60. -90.


;xxg = fltarr (43200,21600)
;yyg = fltarr (43200,21600)
;lai(*,*) = 10000.
;goto,jump

for dofm =280,365,8 do begin  
    lai = fltarr (i_sib,j_sib)
    lai(*,*) =0.
    lai_count = fltarr (i_sib,j_sib)
    lai_count(*,*) = 0. 
 
    mytime=timegen(start=julday(1,1,2006),final=julday(1,1,2006)+dofm,unit='d')
    caldat,mytime,month,day,yyyy
    ii=n_elements(mytime)-1
    date_stamp=string(yyyy[ii],month[ii],day[ii],format='(i4.4,".",i2.2,".",i2.2)')
    MMDD= string(month[ii],day[ii],format='(i2.2,".",i2.2)')   
    for  j_modis = 0,17 do begin
        for i_modis = 0,35 do begin
            gridcell = string(i_modis,j_modis,format='("h",i2.2,"v",i2.2)')
            fname = '/l_data/MODIS/LAI/IDL/tmp/'+gridcell+'.txt' 
            if file_test(fname) then begin
                
                UL_x = 0.
                UL_y = 0.
                LR_x = 0.
                LR_y = 0.
                
                openr,1,fname
                readf,1,nfiles
                char_var ='                '
                readf,1,char_var
                close,1
                print,'Processing ..:',fname,'',strtrim(char_var,2)
                
                FileID=hdf_sd_start('/l_data/MODIS/LAI/'+char_var, /read)
                HDF_SD_AttrInfo, FileID, 1, name = attr_name, data = attr_data
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
                
                filename = '/l_data/MODIS/LAI/YYYY.'+ MMDD +'/MOD15A2.' + gridcell+'.idl'   
                print,filename
                restore,filename

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
                        if((clim_lai(i,j) ge 0.) or (clim_lai(i,j) le 10.)) then begin
                            lai      (n_elements(ii)-1,n_elements(jj)-1) = lai      (n_elements(ii)-1,n_elements(jj)-1) + clim_lai(i,j)
                            lai_count(n_elements(ii)-1,n_elements(jj)-1) = lai_count(n_elements(ii)-1,n_elements(jj)-1) + 1.
                        endif
                    endfor
                endfor
;           file = '/l_data/MODIS/LAI/YYYY.'+ string(month[ii], '(i2.2)') $
;                  +'.'+string(day[ii], '(i2.2)')+'/MOD15A2.' + gridcell+'.idl'           
;       filename = '/l_data/MODIS/LAI/YYYY.03.22/MOD15A2.' + gridcell+'.idl'   
;       print,filename
;       restore,filename
;       plai =fltarr(1200,1200)
;       plai(*,*)=10000.
;       for j = 0,nrows - 1 do begin
;           for i = 0,ncols - 1 do begin
;               xxg= xgrid(i,j)
;               yyg= ygrid(i,j)
;               if((clim_lai(i,j) ge 0.) or (clim_lai(i,j) le 10.)) then begin
;                   plai(i,j)= clim_lai(i,j)
;               endif
;           endfor
;       endfor
;
            endif
        endfor
    endfor
    
    lai = lai/lai_count
    filename = '/l_data/MODIS/LAI/YYYY.' + MMDD +'/MOD15A2.global_2.5.idl'   
    print,filename
    save,lai,file=filename
;jump:
;    restore,'lai.idl'
;    for j = 0,4320 - 1 do begin
;        for i = 0,8640 - 1 do begin
;            if((lai(i,j) ge 0.) or (lai(i,j) le 10.)) then begin
;            endif else begin
;                lai(i,j)= 10000.
;            endelse
;        endfor
;    endfor
endfor


;contour,lai,x25,y25,c_colors=colors,levels=levp,c_label=levp,/fill,title=gridcell
;closeps,'lai.ps'
;save,UL,LR,project_param(n)


end
