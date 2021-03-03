;
; This IDL code is good to plot global maps of vegetation classes
; Contact - Sarith (sarith.p.mahanama@nasa.gov)
;
; Usage -Set path to rst and til files & output path
gfile = '../DC144_91/DC0144xPC0091_DE0360xPE0180-Pfaf.notiny'
pathout = '../DC144_91/'
; End user defined variables

sp,1,1,/f,/h,/color,/land
loadct2,0

ncb=2880
nrb=1440
dx25=7.5/60.
dy25=7.5/60.

x1b=-180. + dx25/2.
xeb= x1b +(ncb-1)*dx25
y1b= -90. + dy25/2.
yeb= y1b +(nrb-1)*dy25

;
xlim=intarr(2)
ylim=intarr(2)
;
nbm=10
levp=findgen(nbm+1)+1
pfstat=lonarr(1)
;
xylim=[-90.,-179.999,90.,180.]
xlen=xylim(3)-xylim(1)
ylen=xylim(2)-xylim(0)
print,'xlen=',xlen
print,'ylen=',ylen
print,xylim
init=replicate(0.,xlen,ylen)
x = indgen(xlen)+xlim(0)
y = indgen(ylen)+ylim(0)
;end
;
i=0l
pfc=long(0.)
tind=0l

cc=0
dum=0.
ncatg=0l
ncatt=0l

file15=gfile + '.til'
openr,15,file15
readf,15,ncatg
close,15

tilefile=strtrim(gfile,2)+'.til'
get_frac, til=tilefile, ntiles=ntiles ,ii=ii, jj=jj,fr=fr

ncatt = ntiles

print,' # of global tiles and land in the set :',ncatg,ncatt

veg1=intarr(ncatt)
veg2=intarr(ncatt)
frc1=fltarr(ncatt)
frc2=fltarr(ncatt)
frc3=fltarr(ncatt)
v1=0
v2=0
f1=0.
f2=0.
f3=0.

mask=replicate(-99.,ncatg)
file15=pathout + 'mosaic_veg_typs_fracs'
  openr,15,file15
  for i=0l,ncatt -1 do begin
	readf,15,'(i8,i8,2(2x,i3),2(2x,f6.4))',tind,pfc,v1,v2,f1,f2
;        print,tind,pfc,v1,v2,f1,f2,f3
	mask(tind)=i
	veg1(i)=v1
        veg2(i)=v2
        frc1(i)=f1   
	frc2(i)=f2
;        frc3(i)=f3
    endfor
close,15
print,'min max veg type 1',min(veg1),max(veg1)
print,'min max veg type 2',min(veg2),max(veg2)

;; BEGINNING THE PLOT
;; ------------------
tit=['PRIMARY VEGETATION TYPE', 'SECONDARY VEGETATION TYPE','PRIMARY FRACTION', 'SECONDARY FRACTION', 'BARESOIL FRACTION']

for iplot=0,3 do begin

if(iplot le 1) then begin
levp=[1,2,3,4,5,6,7,8,9]
colors=[98,96,72,70,9,88,84,81,79,255]
endif

if(iplot gt 1) then begin
upval=1.
lwval=0.
levp=[lwval,lwval+(upval-lwval)/10+indgen(9)*(upval-lwval)/10.,upval]
colors=[16,18,20,21,22,23,24,25,26,28,30]
endif

gplot, init,x,y, /shade, lev=levp, c_color=colors, title=tit(iplot), oceanmask=1,limit=xylim,oc_color=255

;; READING AND PLOTTING STEP BY STEP
;; ---------------------------------

cc=0
pfc=long(0.)
pfc1=long(0.)


  cumar=0.	
  cat25=lonarr(ncb)

  file15=gfile + '.7.5.rst'
  openr,15,file15
  idum=0l
  for i=0,nrb-1 do begin

	readu,15,idum,cat25,idum
	yu=y1b+i*dy25+dy25/2.
	yl=y1b+i*dy25-dy25/2.

	for j=0,ncb-1 do begin

	pfc=cat25(j)
	pfc1=cat25(j)

	if(pfc gt 0) then begin
	if (mask(pfc1) ne -99.) then begin
            
        if(iplot eq 0) then n=veg1(mask(pfc1))-1
        if(iplot eq 1) then n=veg2(mask(pfc1))-1

        if(iplot eq 2) then n=floor((frc1(mask(pfc1))-lwval)*10./(upval-lwval))
        if(iplot eq 3) then n=floor((frc2(mask(pfc1))-lwval)*10./(upval-lwval))
        if(iplot eq 4) then n=floor((frc3(mask(pfc1))-lwval)*10./(upval-lwval))

	xl=x1b+j*dx25-dx25/2.
	xr=x1b+j*dx25+dx25/2.
	xx=fltarr(5)
	yy=fltarr(5)
	xx=[xl,xl,xr,xr,xl]
	yy=[yu,yl,yl,yu,yu]

		if (n lt 0) then n=9
		polyfill,xx,yy,color=colors(n)
		cc=cc+1

	endif
	endif
	endfor

  endfor
  close,15

xleg=-160.
yleg1=-80.
close,12
subt1='1 = Broadleaf evergreen  2 = Broadleaf decidous  3 = Needleleaf'
subt2='4 = Grassland            5 = Broadleaf shrubs    6 = Dwarf '
subt3='7 = Baresoil             8 = Desert soil'
if(iplot le 1) then begin
xyouts,xleg,-60.,subt1,charsize=0.8,charthick=1.2
xyouts,xleg,-70.,subt2,charsize=0.8,charthick=1.2
xyouts,xleg,-80.,subt3,charsize=0.8,charthick=1.2
endif
endfor

fileps=pathout + 'mosaic_veg_typs_fracs.ps'

wait,30
closeps,fileps
end

