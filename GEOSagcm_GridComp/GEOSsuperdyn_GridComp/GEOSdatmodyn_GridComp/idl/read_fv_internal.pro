pro read_fv_internal,f=f,im=im,jm=jm


close,1
openr,1,/f77,f


i6=lonarr(6)
i5=lonarr(5)


readu,1,i6
readu,1,i5


LM=i5(2)
aks=dblarr(i5(2)+1)
bks=dblarr(i5(2)+1)
a=dblarr(im,jm)


stop

readu,1,aks
readu,1,bks

ird=0

; Order of VARS in restart post-header
; U  (:,:,1:LM)
; V  (:,:,1:LM)
; PT (:,:,1:LM)
; PE (:,:,0:LM)
; PKZ(:,:,1:LM)


u  = fltarr(im,jm,lm)
v  = fltarr(im,jm,lm)
pt = fltarr(im,jm,lm)
pe = fltarr(im,jm,lm+1)
pkz= fltarr(im,jm,lm)


for l=0,lm-1 do begin
    readu,1,a
    u(*,*,l) = a
endfor
for l=0,lm-1 do begin
    readu,1,a
    v(*,*,l) = a
endfor
for l=0,lm-1 do begin
    readu,1,a
    pt(*,*,l) = a
endfor
for l=0,lm do begin
    readu,1,a
    pe(*,*,l) = a
endfor
for l=0,lm-1 do begin
    readu,1,a
    pkz(*,*,l) = a
endfor



stop
return
end
