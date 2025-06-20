pro read_lnd_albedos,f=f,runthru=runthru,alb=alb

get_lun,lun
openr,lun,f,/f77

date1=fltarr(6)
date2=fltarr(6)
res  =fltarr(2)
readu,lun,DATE1,DATE2,RES
point_lun,lun,0
albd=fltarr( long( res(0)+.0001 ) )

albo={origfile:f, date1:date1*0.-9999. , date2:date2*0.-9999. , res:res, albedo:albd }

alb = replicate( albo , 14 )

for i=0,13 do begin

;;readu,lun,year1,month1,day1,hour1,min1,sec1

readu,lun,DATE1,DATE2,RES
readu,lun,albd

alb(i).date1  = date1
alb(i).date2  = date2
alb(i).res    = res
alb(i).albedo = albd

;;comment

if not keyword_set(runthru) then STOP

endfor

close,lun

return
end
