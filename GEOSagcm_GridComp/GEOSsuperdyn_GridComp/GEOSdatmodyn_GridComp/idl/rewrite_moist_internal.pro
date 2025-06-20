pro rewrite_moist_internal,f=f,writefile=writefile,odir=odir,xdir=xdir $
                       ,gridname=grid$,ogr=ogr,xgr=xgr  $
                       ,im=im,jm=jm,lm=lm,xim=xim,xjm=xjm

fx=f+'_'+grid$

close,1,2
openr,1,/f77,odir+f

 if keyword_set( writefile) then openw,2,/f77,xdir+fx


; Order of VARS in restart
; QV   (:,:,1:LM)
; QLLS (:,:,1:LM)
; QLCN (:,:,1:LM)
; CLLS (:,:,1:LM)
; CLCN (:,:,1:LM)
; QILS (:,:,1:LM)
; QICN (:,:,1:LM)

if not keyword_set(ogr) then ogr={lon:findgen(im)*360. /(1.*im)  , lat:findgen(jm)*180. /(1.*jm) }
if not keyword_set(xgr) then xgr={lon:findgen(xim)*360./(1.*xim) , lat:findgen(xjm)*180./(1.*xjm) }

a=fltarr(     n_elements(ogr.lon) ,  n_elements(ogr.lat)   )

;QV
for l=0,lm-1 do begin
    readu,1,a ; & a(*,0)=0.0
    ax = int2d( a , from=ogr, to=xgr, /noshifts ) ; & ax(*,0)=bad
    cleanunder,ax
    if keyword_set( writefile) then writeu,2,ax
endfor

;QLLS
for l=0,lm-1 do begin
    readu,1,a 
    ax = int2d( a , from=ogr, to=xgr, /noshifts ) ; & ax(*,0)=bad
    cleanunder,ax
    if keyword_set( writefile) then writeu,2,ax
endfor

;QLCN
for l=0,lm-1 do begin
    readu,1,a ; & a(*,0)=0.0
    ax = int2d( a , from=ogr, to=xgr, /noshifts ) ; & ax(*,0)=bad
    cleanunder,ax
    if keyword_set( writefile) then writeu,2,ax
endfor

;CLLS
for l=0,lm-1 do begin
    readu,1,a ; & a(*,0)=0.0
    ax = int2d( a , from=ogr, to=xgr, /noshifts ) ; & ax(*,0)=bad
    cleanunder,ax
    if keyword_set( writefile) then writeu,2,ax
endfor

;CLLS
for l=0,lm-1 do begin
    readu,1,a ; & a(*,0)=0.0
    ax = int2d( a , from=ogr, to=xgr, /noshifts ) ; & ax(*,0)=bad
    cleanunder,ax
    if keyword_set( writefile) then writeu,2,ax
endfor

;QILS
for l=0,lm-1 do begin
    readu,1,a ; & a(*,0)=0.0
    ax = int2d( a , from=ogr, to=xgr, /noshifts ) ; & ax(*,0)=bad
    cleanunder,ax
    if keyword_set( writefile) then writeu,2,ax
endfor

;QILS
for l=0,lm-1 do begin
    readu,1,a ; & a(*,0)=0.0
    ax = int2d( a , from=ogr, to=xgr, /noshifts ) ; & ax(*,0)=bad
    cleanunder,ax
    if keyword_set( writefile) then writeu,2,ax
endfor


close,1,2
return
end
