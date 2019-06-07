pro rewrite_fv_internal,f=f,writefile=writefile,odir=odir,xdir=xdir $
                       ,gridname=grid$,ogr=ogr,xgr=xgr  $
                       ,im=im,jm=jm,xim=xim,xjm=xjm,lm=lm


fx=f+'_'+grid$

close,1,2
openr,1,/f77,odir+f

 if keyword_set( writefile) then openw,2,/f77,xdir+fx


i6=lonarr(6)
i5=lonarr(5)


readu,1,i6
readu,1,i5

 if keyword_set( writefile) then writeu,2,i6
 if keyword_set( writefile) then writeu,2,i5


LM=i5(2)
aks=dblarr(i5(2)+1)
bks=dblarr(i5(2)+1)

readu,1,aks
readu,1,bks

 if keyword_set( writefile) then writeu,2,aks
 if keyword_set( writefile) then writeu,2,bks

ird=0

; Order of VARS in restart post-header
; U  (:,:,1:LM)
; V  (:,:,1:LM)
; PT (:,:,1:LM)
; PE (:,:,0:LM)
; PKZ(:,:,1:LM)


if not keyword_set(ogr) then ogr={lon:findgen(im)*360. /(1.*im)  , lat:findgen(jm)*180. /(1.*jm) }
if not keyword_set(xgr) then xgr={lon:findgen(xim)*360./(1.*xim) , lat:findgen(xjm)*180./(1.*xjm) }


a=dblarr(     n_elements(ogr.lon) ,  n_elements(ogr.lat)   )

bad = double( 9.99999999999e+14 )

;U
for l=0,lm-1 do begin
    readu,1,a & a(*,0)=0.0
    ax = double( int2d( a , from=ogr, to=xgr, /noshifts ) ) & ax(*,0)=bad
    if keyword_set( writefile) then writeu,2,ax
endfor

;V
for l=0,lm-1 do begin
    readu,1,a 
    ax = double( int2d( a , from=ogr, to=xgr, /noshifts ) ) ; & ax(*,0)=bad
    if keyword_set( writefile) then writeu,2,ax
endfor

;PT
for l=0,lm-1 do begin
    readu,1,a ; & a(*,0)=0.0
    ax = double( int2d( a , from=ogr, to=xgr, /noshifts ) ) ; & ax(*,0)=bad
    if keyword_set( writefile) then writeu,2,ax
endfor

;PE
for l=0,lm do begin
    readu,1,a ; & a(*,0)=0.0
    ax = double( int2d( a , from=ogr, to=xgr, /noshifts ) ) ; & ax(*,0)=bad
    if keyword_set( writefile) then writeu,2,ax
endfor

;PKZ
for l=0,lm-1 do begin
    readu,1,a ; & a(*,0)=0.0
    ax = double( int2d( a , from=ogr, to=xgr, /noshifts ) ) ; & ax(*,0)=bad
    if keyword_set( writefile) then writeu,2,ax
endfor


close,1,2

return
end
