pro write_datmodyn_internal,f=f,writefile=writefile,odir=odir,xdir=xdir $
                       ,gridname=grid$,ogr=ogr,xgr=xgr  $
                       ,im=im,jm=jm,xim=xim,xjm=xjm,lm=lm


fx='datmodyn_internal_rst_'+grid$

close,1,2
openr,1,/f77,odir+f

 if keyword_set( writefile) then openw,2,/f77,xdir+fx


i6=lonarr(6)
i5=lonarr(5)


readu,1,i6
readu,1,i5



LM=i5(2)
aks=dblarr(i5(2)+1)
bks=dblarr(i5(2)+1)

readu,1,aks
readu,1,bks

ird=0

pref =  float( bks*100000.0d +aks )
;STOP

; Order of VARS in restart post-header
; U  (:,:,1:LM)
; V  (:,:,1:LM)
; PT (:,:,1:LM)
; PE (:,:,0:LM)
; PKZ(:,:,1:LM)


if not keyword_set(ogr) then ogr={lon:findgen(im)*360. /(1.*im)  , lat:findgen(jm)*180. /(1.*jm) }
if not keyword_set(xgr) then xgr={lon:findgen(xim)*360./(1.*xim) , lat:findgen(xjm)*180./(1.*xjm) }


a   =dblarr(     n_elements(ogr.lon) ,  n_elements(ogr.lat)   )
ux  =fltarr(     n_elements(xgr.lon) ,  n_elements(xgr.lat) , lm   )
vx  =fltarr(     n_elements(xgr.lon) ,  n_elements(xgr.lat) , lm   )
tx  =fltarr(     n_elements(xgr.lon) ,  n_elements(xgr.lat) , lm   )
plex=fltarr(     n_elements(xgr.lon) ,  n_elements(xgr.lat) , lm+1 )
pkx =fltarr(     n_elements(xgr.lon) ,  n_elements(xgr.lat) , lm   )
omx =fltarr(     n_elements(xgr.lon) ,  n_elements(xgr.lat) , lm+1 )

bad = double( 9.99999999999e+14 )

;U
for l=0,lm-1 do begin
    readu,1,a & a(*,0)=0.0
    ax = double( int2d( a , from=ogr, to=xgr, /noshifts ) ) & ax(*,0)=bad
    ux(*,*,l) = float( ax )
endfor

;V
for l=0,lm-1 do begin
    readu,1,a 
    ax = double( int2d( a , from=ogr, to=xgr, /noshifts ) ) ; & ax(*,0)=bad
    vx(*,*,l) = float( ax )
endfor

;PT
for l=0,lm-1 do begin
    readu,1,a ; & a(*,0)=0.0
    ax = double( int2d( a , from=ogr, to=xgr, /noshifts ) ) ; & ax(*,0)=bad
    tx(*,*,l) = float( ax )
endfor

;PE
for l=0,lm do begin
    readu,1,a ; & a(*,0)=0.0
    ax = double( int2d( a , from=ogr, to=xgr, /noshifts ) ) ; & ax(*,0)=bad
    plex(*,*,l) = float( ax )
endfor

;PKZ
for l=0,lm-1 do begin
    readu,1,a ; & a(*,0)=0.0
    ax = double( int2d( a , from=ogr, to=xgr, /noshifts ) ) ; & ax(*,0)=bad
    pkx(*,*,l) = float( ax )
endfor

tex=tx * pkx ; FV therm var to Temp

if keyword_set(writefile) then begin

;for l=0,lm do begin
;    writeu,2,pref(l)
;endfor
writeu,2,pref
for l=0,lm do begin
    writeu,2,plex(*,*,l)
endfor
for l=0,lm-1 do begin
    writeu,2,tex(*,*,l)
endfor
for l=0,lm-1 do begin
    writeu,2,ux(*,*,l)
endfor
for l=0,lm-1 do begin
    writeu,2,vx(*,*,l)
endfor
for l=0,lm do begin
    writeu,2,omx(*,*,l)
endfor

endif
close,1,2

return
end
