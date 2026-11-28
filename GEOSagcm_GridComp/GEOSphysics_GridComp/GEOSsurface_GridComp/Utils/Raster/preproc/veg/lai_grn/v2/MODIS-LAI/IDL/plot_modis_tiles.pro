n_levels = 21
loadct2,0
colors = indgen(n_levels)+ 80
upval = 10.
lwval = 0.
levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]

NLEN = 1200

for jx =1,18 do begin
    for ix =1,36 do begin
        filename = '/l_data/model_parameters/LAI/MODIS-LAI/MODIS_8-DayClim/MODIS_lai_clim.H'+string(ix,'(i2.2)')+'V'+string(jx,'(i2.2)') + '.nc'
        if file_test (filename) then begin
           IXJX = 'H'+string(ix,'(i2.2)')+'V'+string(jx,'(i2.2)')
            ncid = NCDF_OPEN(filename)
            ncdf_attget,ncid,'i_ind_offset_LL',iLL,/global
            ncdf_attget,ncid,'j_ind_offset_LL',jLL,/global
            ncdf_varget,ncid,'latitude',y
            ncdf_varget,ncid,'longitude',x
            ncdf_varget,ncid,'LAI',var1
            ncdf_attget,ncid,'LAI','ScaleFactor',sf
            ncdf_attget,ncid,'LAI','UNDEF',undef
            ncdf_varget,ncid,'MMDD', mmdd
            mmdd = string(mmdd)
            NCDF_CLOSE, ncid 
         !P.position=0
         !P.Multi = [0, 2, 2, 0, 1]   
          Erase, Color= 255
            for t = 21,24 do begin

               plot_var = fltarr (nlen, nlen)
               
               for j = 0,NLEN -1 do begin
                  for i = 0, NLEN -1 do begin 
                      if(var1(i,j,t) ne undef) then plot_var (i,j) = sf*var1(i,j,t)
                   endfor
               endfor
                contour, plot_var,x,y,levels = levels,c_colors=colors,title=IXJX + ' - ' + mmdd(t)
             endfor
            write_to_jpeg,file_name= 'MODIS_' + IXJX + '.jpg'
        endif
    endfor
endfor

end
