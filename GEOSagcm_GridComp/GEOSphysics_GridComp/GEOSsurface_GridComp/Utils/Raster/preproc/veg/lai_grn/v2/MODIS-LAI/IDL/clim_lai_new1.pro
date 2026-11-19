yearB = 2000
yearE = 2013

; Looping through MODIS tiles

for  j_modis = 8,8 do begin
    for i_modis = 0,35 do begin
       gridcell = string(i_modis,j_modis,format='("h",i2.2,"v",i2.2)')
       fname = '/l_data/model_parameters/LAI/MODIS-LAI/IDL/tmp_new/'+gridcell+'.txt'

       ; if data files available for the tile

       if file_test(fname) then begin
           print,'Processing ..:',fname

           ; looping through 8-day composites

           for dofm =0,365,8 do begin                  
               openr,1,fname
               readf,1,nfiles
               nf_count =0
               char_var ='                '

               ; looping through years
               clim_lai = fltarr (1200,1200)
               px_count = fltarr (1200,1200)
               adum     = fltarr (1200,1200) 
               clim_lai(*,*) = !VALUES.F_NAN
               px_count(*,*) = 0.
               
               for year = yearB,yearE do begin
                   mytime=timegen(start=julday(1,1,year),final=julday(1,1,year)+dofm,unit='d')
                   caldat,mytime,month,day,yyyy
                   ii=n_elements(mytime)-1
                   date_stamp=string(yyyy[ii],month[ii],day[ii],format='(i4.4,".",i2.2,".",i2.2)')        

                   for n = nf_count,nfiles -1 do begin
                       readf,1,char_var
                       nf_count=n+1
                       if(strtrim(date_stamp,2) eq strtrim(strmid(char_var,0,10),2)) then begin
                         n =   nfiles+1
                         print,date_stamp,' ',char_var
                         if(string(char_var,'(a10)') eq '2011.04.23' or strtrim(char_var,2) $
                           eq '2005.02.18/MOD15A2.A2005049.h11v02.005.2008199065834.hdf' or $
                           strtrim(char_var,2) eq '2003.02.18/MOD15A2.A2003049.h12v04.005.2007279212643.hdf' $
                           or strtrim(char_var,2) eq '2010.10.24/MOD15A2.A2010297.h33v11.005.2010318080800.hdf' $
                           or strtrim(char_var,2) eq '2000.10.31/MOD15A2.A2000305.h33v11.005.2008211064833.hdf' $
                           or strtrim(char_var,2) eq '2011.02.18/MOD15A2.A2011049.h11v08.005.2011060002920.hdf' $
                           ) then begin ;eq $
;                            '2011.04.23/MOD15A2.A2011113.h12v01.005.2011123074311.hdf' or  $
;                            strtrim(char_var,2) eq $
;                            '2011.04.23/MOD15A2.A2011113.h13v01.005.2011123093522.hdf') then begin
;                             restore,'/l_data/MODIS/LAI/'+strtrim(char_var,2)+'.idl'
                         print,'skipping ..',strtrim(char_var,2)
                         endif else begin
                         ; reading HDF file
                         fid =0l
                         fid=hdf_sd_start('/l_data/model_parameters/LAI/MODIS-LAI/'+char_var, /read)
                         hdf_sd_fileinfo, fid, datasets, attributes
                         sds_id = hdf_sd_select(fid, 1)
                         hdf_sd_getinfo, sds_id, dims=dims, name=name, type=type, unit=unit, natts=natts
                         hdf_sd_getdata, sds_id,adum

                         HDF_SD_ENDACCESS, sds_id
                         HDF_SD_END, fid
                     endelse
                         for j =0,1199 do begin
                             for i =0,1199 do begin
                                 if ((adum(i,j) ge 0) and (adum(i,j) le 100)) then begin
                                     if (px_count (i,j) eq 0.) then clim_lai(i,j) = 0.
                                     clim_lai(i,j) =  clim_lai(i,j) + 0.1*adum(i,j)
                                     px_count (i,j) = px_count (i,j) + 1.
                                 endif
                             endfor
                         endfor                                                  
                      endif                     
                  endfor

;                  print,year,' ',date_stamp,n,nf_count,nfiles -1
                   if(nf_count eq nfiles) then begin
                       close,1
                       openr,1,fname
                       readf,1,nfiles
                       nf_count =0
                    endif
                endfor                              
                close,1

                clim_lai =  clim_lai/px_count
                print, max(clim_lai,/nan),min(clim_lai,/nan)
                save_file = '/l_data/model_parameters/LAI/MODIS-LAI/YYYY.'+ string(month[ii], '(i2.2)') $
                  +'.'+string(day[ii], '(i2.2)')+'/MOD15A2.' + gridcell+'_new.idl'
                save,clim_lai,file=save_file
            endfor
        endif
    endfor
endfor

end
