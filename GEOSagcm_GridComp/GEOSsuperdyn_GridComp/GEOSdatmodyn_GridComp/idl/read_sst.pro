pro read_sst,f=f,runthru=runthru

get_lun,lun
openr,lun,f,/f77

sst=fltarr(360,180)

date1=fltarr(6)
date2=fltarr(6)
res  =fltarr(2)


for i=0,10000 do begin

;;readu,lun,year1,month1,day1,hour1,min1,sec1

readu,lun,DATE1,DATE2,RES
readu,lun,sst

print,DATE1
print,DATE2
print,RES

if not keyword_set(runthru) then STOP

endfor

close,lun
STOP



return
end
