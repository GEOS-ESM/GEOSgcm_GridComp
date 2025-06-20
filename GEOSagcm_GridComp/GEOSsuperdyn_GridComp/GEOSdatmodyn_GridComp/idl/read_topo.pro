pro read_topo,f=f,im=im,jm=jm,a=a

a=fltarr(im,jm)

close,1
openr,1,/f77,f

readu,1,a


close,1
return
end
