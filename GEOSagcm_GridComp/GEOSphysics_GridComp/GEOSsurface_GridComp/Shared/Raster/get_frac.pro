pro get_frac, til=til, ntiles=ntiles, ii=ii, jj=jj, fr=fr
print, 'Reading '+til
start=8  ;; start index + header
tempstuff=read_ascii(til, data_start=start)
land=where(abs(tempstuff.field01[0,*]-100) lt 1.0E-5, ntiles)
itemp=tempstuff.field01[4,*]
jtemp=tempstuff.field01[5,*]
ftemp=tempstuff.field01[6,*]
ii=itemp[land]
print, 'Range of tile I Coordinates', max(ii), min(ii)
jj=jtemp[land]
print, 'Range of tile J Coordinates', max(jj), min(jj)
fr=ftemp[land]
print, 'Range of tile FR values', max(fr), min(fr)
return
end
