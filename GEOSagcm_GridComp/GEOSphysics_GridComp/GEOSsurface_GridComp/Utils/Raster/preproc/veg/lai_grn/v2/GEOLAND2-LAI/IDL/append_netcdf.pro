 MMDD1 = ['0103', '0113', '0124', '0203', '0213', '0221', '0303', '0313', '0324', '0403', '0413', '0423', $
          '0503', '0513', '0524', '0603', '0613', '0623', '0703', '0713', '0724', '0803', '0813', '0824', $
          '0903', '0913', '0923', '1003', '1013', '1024', '1103', '1113', '1123', '1203', '1213', '1224' ]

MMDD = strarr(36)

for n = 1,n_elements(MMDD1) -1 do begin

    m = fix(strmid(mmdd1(n),0,2))
    d = fix(strmid(mmdd1(n),2,2))
    mytime=timegen(final=julday(m,d,2001)-5,unit='d');final=julday(m,d,2001)+5,
    caldat,mytime,month,day,year
    i=n_elements(mytime)-1
    mmdd(n-1) = string(month[i],day[i],format='(i2.2,i2.2)')
    print,mmdd1(n),' ',mmdd(n)

endfor

n = 0
    m = fix(strmid(mmdd1(n),0,2))
    d = fix(strmid(mmdd1(n),2,2))
    mytime=timegen(final=julday(m,d,2001)-5,unit='d');final=julday(m,d,2001)+5,
    caldat,mytime,month,day,year
    i=n_elements(mytime)-1
    mmdd(35) = string(month[i],day[i],format='(i2.2,i2.2)')
    print,mmdd1(n),' ',mmdd(n)

spawn, 'ls /l_data/GEOLAND2-LAI/LAI-CLIM/netCDF2/geoland2_lai_clim.*.nc',flist



for n = 0,n_elements(flist) -1 do begin

    print,flist(n)
    id = NCDF_OPEN(flist(n), /WRITE)
    NCDF_VARPUT, id, 'MMDD', MMDD
    NCDF_VARGET, id, 'LAI', lai_old
    LAI = lai_old

    for t = 0,34 do begin
        lai(*,*,t) = lai_old(*,*,t+1)
    endfor
    lai(*,*,35) = lai_old(*,*,0)
    NCDF_VARPUT, id, 'LAI', lai,offset=[0,0,0]
    NCDF_CLOSE, id

endfor

end
