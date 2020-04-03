;
; This IDL code is good to plot global maps of soil types and soil
; soil hydraulic properties
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

file15=gfile +'.til'
openr,15,file15
readf,15,ncatg
close,15

tilefile=strtrim(gfile,2)+'.til'
get_frac, til=tilefile, ntiles=ntiles ,ii=ii, jj=jj,fr=fr

ncatt = ntiles

print,' # of global tiles and land in the set :',ncatg,ncatt

sol1=intarr(ncatt)
sol2=intarr(ncatt)
bee=fltarr(ncatt)
psi=fltarr(ncatt)
por=fltarr(ncatt)
con=fltarr(ncatt)
wil=fltarr(ncatt)
dep=fltarr(ncatt)

v1=0
v2=0
f1=0.
f2=0.
f3=0.
f4=0.
f5=0.
f6=0.
f7=0.

mask=replicate(-99.,ncatg)
file15=pathout + 'soil_param.dat'
  openr,15,file15
  for i=0l,ncatt -1 do begin
	readf,15,'(i8,i8,i3,i3,2f7.2,f6.3,f12.8,f7.4,f10.3)',tind,pfc,v1,v2,f1,f2,f3,f4,f5,f6
;        print,tind,pfc,v1,v2,f1,f2,f3
	mask(tind)=i
	sol1(i)=v1
        sol2(i)=v2
        bee(i)=f1   
	psi(i)=f2
        por(i)=f3
        con(i)=f4
        wil(i)=f5*f3
        dep(i)=f6
    endfor
close,15

;; BEGINNING THE PLOT
;; ------------------
tit=['SOIL CLASS (0-30cm)', 'SOIL CLASSS (0-100cm)' ,'B ','PSIS (m)','POROSITY (v/v)', 'Saturated hydraulic conductivity at the surface (m/s)', 'WILTING POINT (v/v)','SOIL DEPTH (mm)']

for iplot=0,7 do begin

if(iplot eq 2) then begin
upval=13.
lwval=3.
levp=[lwval,lwval+(upval-lwval)/10+indgen(9)*(upval-lwval)/10.,upval]
colors=[16,18,20,21,22,23,24,25,26,28,30]
endif

if(iplot eq 3) then begin
upval=-0.05
lwval=-0.85
levp=[lwval,lwval+(upval-lwval)/10+indgen(9)*(upval-lwval)/10.,upval]
colors=[16,18,20,21,22,23,24,25,26,28,30]
endif

if(iplot eq 4) then begin
upval=0.5
lwval=0.3
levp=[lwval,lwval+(upval-lwval)/10+indgen(9)*(upval-lwval)/10.,upval]
colors=[16,18,20,21,22,23,24,25,26,28,30]
endif

if(iplot eq 5) then begin
upval=0.03
lwval=0.001
levp=[lwval,lwval+(upval-lwval)/10+indgen(9)*(upval-lwval)/10.,upval]
colors=[16,18,20,21,22,23,24,25,26,28,30]
endif

if(iplot eq 6) then begin
upval=0.3
lwval=0.03
levp=[lwval,lwval+(upval-lwval)/10+indgen(9)*(upval-lwval)/10.,upval]
colors=[16,18,20,21,22,23,24,25,26,28,30]
endif

if(iplot eq 7) then begin
upval=8000.
lwval=0.
levp=[lwval,lwval+(upval-lwval)/10+indgen(9)*(upval-lwval)/10.,upval]
colors=[16,18,20,21,22,23,24,25,26,28,30]
endif

if(iplot le 1) then begin
colors=[16,17,18,20,21,23,24,25,27,28,29,30,31,32,33,34,36,37,39,40]
levp=indgen(13)+1.
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
            
        if(iplot eq 0) then n=sol1(mask(pfc1))-1
        if(iplot eq 1) then n=sol2(mask(pfc1))-1

        if(iplot eq 2) then n=floor((bee(mask(pfc1))-lwval)*10./(upval-lwval))
        if(iplot eq 3) then n=floor((psi(mask(pfc1))-lwval)*10./(upval-lwval))
        if(iplot eq 4) then n=floor((por(mask(pfc1))-lwval)*10./(upval-lwval))
        if(iplot eq 5) then n=floor((con(mask(pfc1))-lwval)*10./(upval-lwval))
        if(iplot eq 6) then n=floor((wil(mask(pfc1))-lwval)*10./(upval-lwval))
        if(iplot eq 7) then n=floor((dep(mask(pfc1))-lwval)*10./(upval-lwval))

	xl=x1b+j*dx25-dx25/2.
	xr=x1b+j*dx25+dx25/2.
	xx=fltarr(5)
	yy=fltarr(5)
	xx=[xl,xl,xr,xr,xl]
	yy=[yu,yl,yl,yu,yu]

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
subt1='1 = Sand            2 = Loamy Sand       3 = Sandy Loam'
subt2='4 = Silt Loam       5 = Silt             6 = Loam      '
subt3='7 = Sandy Clay Loam 8 = Silty Clay Loam  9 = Clay Loam'
subt4='10= Sandy Clay      11= Silty Clay       12= Clay' 
if(iplot le 1) then begin
xyouts,xleg,-55.,subt1,charsize=0.8,charthick=1.2
xyouts,xleg,-65.,subt2,charsize=0.8,charthick=1.2
xyouts,xleg,-75.,subt3,charsize=0.8,charthick=1.2
xyouts,xleg,-85.,subt4,charsize=0.8,charthick=1.2
endif
if(iplot eq 1) then begin
sp,1,3,/f,/h,/color
endif
endfor

fileps=pathout +'soil_param.ps'
closeps,fileps

end

