function int2d,arrx,fromg=fromg,tog=tog,noshifts=noshifts

tolat=tog.lat
tolon=tog.lon
frlon=fromg.lon
frlat=fromg.lat

arr=arrx
shift_result=0

nxf=n_elements(frlon) & nyf=n_elements(frlat)
nx2=n_elements(tolon) & ny2=n_elements(tolat)

if not keyword_set(noshifts) then begin 
if min(tolon) lt -179./50. then begin
   tolon= shift( tolon, nx2/2 )
   tolon(where(tolon lt 0))=tolon(where(tolon lt 0))+360.
   shift_result=1
endif
if min(frlon) lt -179./50. then begin
   frlon= shift( frlon, nxf/2 )
   frlon(where(frlon lt 0))=frlon(where(frlon lt 0))+360.
   arr=shift(arr,nxf/2,0)
    
endif
endif

iarr=fltarr(ny2,nx2)
tmpa=fltarr(nx2,nyf)

for j=0,nyf-1 do begin
    d=interpol( arr(*,j) , frlon, tolon )
    tmpa(*,j)=d
endfor
tmpa=transpose(tmpa)
for j=0,nx2-1 do begin
    d=interpol( tmpa(*,j) , frlat, tolat )
    iarr(*,j)=d
endfor
iarr=transpose(iarr)
   
if shift_result eq 1 then begin
  iarr= shift( iarr, -nx2/2 , 0)
endif

for j=0,nx2-1 do begin
    extra=where( tolat gt max(frlat) )
    if min(extra) gt -1 then iarr( j , extra )=-9.99e9
endfor
for j=0,nx2-1 do begin
    extra=where( tolat lt min(frlat) )
    if min(extra) gt -1 then iarr( j , extra )=-9.99e9
endfor


return,iarr
end
