
loadct2,0
colors=90+indgen(11)
upval =10.0
lwval =0.
levp = [lwval,lwval+(upval-lwval)/10+indgen(10)*(upval-lwval)/10.]

i_sib=8640
j_sib=4320

x25=indgen(i_sib)*2.5/60. -180.
y25=indgen(j_sib)*2.5/60. -90.

;for dofm =0,365,8 do begin  
;for dofm =0,365,8 do begin  
for dofm =272,296,8 do begin
    mytime=timegen(start=julday(1,1,2006),final=julday(1,1,2006)+dofm,unit='d')
    caldat,mytime,month,day,yyyy
    ii=n_elements(mytime)-1
    date_stamp=string(yyyy[ii],month[ii],day[ii],format='(i4.4,".",i2.2,".",i2.2)')
    MMDD= string(month[ii],day[ii],format='(i2.2,".",i2.2)')   

    fname = '/l_data/model_parameters/LAI/MODIS-LAI/YYYY.' + MMDD +'/MOD15A2.global_2.5.idl'   
    fout = '/l_data/model_parameters/LAI/MODIS-LAI/MODIS_8-DayClim/MOD15A2.YYYY.' + MMDD + '.global_2.5min.data'
    openw, lun,fout,/get_lun,/F77_UNFORMATTED

    if file_test(fname) then begin
;        sp,1,1,/f,/h,/color,/land
        restore,fname
        for j = 0,4320 - 1 do begin
            for i = 0,8640 - 1 do begin
                if((lai(i,j) ge 0.) or (lai(i,j) le 10.)) then begin
                endif else begin
                    lai(i,j)= 9999.
                endelse
            endfor
        endfor
        contour,lai,x25,y25,c_colors=colors,levels=levp,c_label=levp,/fill,title=MMDD
        for j = 0,4320 - 1 do begin
            strip = lai (*,j)
            writeu,lun,strip
        endfor
        closeps,strtrim(fname)+'.ps'
        close,lun
        free_lun, lun

    endif
endfor
end
