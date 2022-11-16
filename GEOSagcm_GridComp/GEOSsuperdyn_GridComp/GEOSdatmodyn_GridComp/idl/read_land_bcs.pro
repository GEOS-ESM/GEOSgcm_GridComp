pro read_land_bcs,f=f,runthru=runthru,dataset=dataset,serieslength=serieslength

get_lun,lun
openr,lun,f,/f77

date1=fltarr(6)
date2=fltarr(6)
res  =fltarr(2)
readu,lun,DATE1,DATE2,RES
point_lun,lun,0
valuetmp=fltarr( long( res(0)+.0001 ) )

record={origfile:f, date1:date1*0.-9999. , date2:date2*0.-9999. , res:res, value:valuetmp }

;dataset = replicate( record , 14 )
dataset = replicate( record , serieslength+2 )

;for i=0,13 do begin
for i=0,serieslength+1 do begin

;;readu,lun,year1,month1,day1,hour1,min1,sec1

readu,lun,DATE1,DATE2,RES
readu,lun,valuetmp

print,f
print,date1
print,date2
print,res


dataset(i).date1  = date1
dataset(i).date2  = date2
dataset(i).res    = res
dataset(i).value = valuetmp

;;comment

if not keyword_set(runthru) then STOP

endfor

close,lun

return
end
