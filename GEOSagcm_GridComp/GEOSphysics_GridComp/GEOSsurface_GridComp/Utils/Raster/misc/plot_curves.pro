pro plot_curves

arw1=0.
arw2=0.
arw3=0.
arw4=0.
ars1=0.
ars2=0.
ars3=0.
cti_mean=0. 
cti_std =0.
cti_min =0.
cti_max =0.
cti_skew=0.
BEE  =0. 
PSIS =0.
POROS=0.
COND =0.
WPWET=0.
soildepth=0.
nbdep=0
nbdepl=0
wmin0=0.
cdcr1=0.
cdcr2=0.

file3='file.0000001'
openr,12,file3
readf,12,cti_mean, cti_std,cti_min, cti_max, cti_skew
readf,12,BEE, PSIS,POROS,COND,WPWET,soildepth
readf,12,nbdep,nbdepl,wmin0,cdcr1,cdcr2

catdef = fltarr(nbdep)
ar1    = fltarr(nbdep)
wmin   = fltarr(nbdep)

readf,12,catdef
readf,12,ar1
readf,12,wmin
readf,12,ars1,ars2,ars3
readf,12,arw1,arw2,arw3,arw4
close,12

ntot=nbdep
x=indgen(ntot)
y  = fltarr(ntot)

plot,x,y,xrange=[0.,max(catdef)],yrange=[0.,1.],linestyle=0,title='WMIN and AR1'
oplot,catdef,wmin,thick=2
oplot,[cdcr1,cdcr1],[0.,1]
oplot,[cdcr2,cdcr2],[0.,1]

;ntot=fix(catdef(nbdep-1))+ 1.
ntot=fix(cdcr1)+ 1.
x=indgen(ntot)
y=fltarr(ntot)
y2=fltarr(ntot)
for n =0,ntot-1 do begin

        y(n) =  arw4 + (1.- arw4)*(1.+ arw1*x(n))/(1.+ arw2*x(n)+  arw3*x(n)*x(n))
        y2(n) =  (1.+ ars1*x(n))/(1.+ ars2*x(n)+  ars3*x(n)*x(n))
endfor

oplot,x,y,color = 98
oplot,catdef,ar1,thick=2
oplot,x,y2,color = 98

thisDevice = !D.NAME
SET_PLOT, 'Z'

TVLCT, red, green, blue, /GET
SET_PLOT, thisDevice

image = TVRD()
thisImage = BYTSCL(image)
s = SIZE(thisImage)
image3d = BYTARR(s(1), s(2),3)
image3d( *, *,0) = red(thisImage)
image3d( *, *,1) = green(thisImage)
image3d( *, *,2) = blue(thisImage)

write_jpeg,'img.0000001.jpg',image3d,quality=75,true=3

wdelete
end
