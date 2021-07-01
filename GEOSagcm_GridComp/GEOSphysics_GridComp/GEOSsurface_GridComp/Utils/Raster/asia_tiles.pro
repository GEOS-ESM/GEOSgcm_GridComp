;
; This IDL code is good to plot land-tile boundaries for the lat/lon
; extermeties bounded by xylim=[20.,110.,50.,150.]
; Contact - Sarith (sarith.p.mahanama@nasa.gov)
;
; Usage -Set path to rst and til files
gfile = '../clsm/DC0144xPC0091_DE0360xPE0180-Pfaf.notiny'
; End user defined variables

sp,1,1,/f,/h,/color,/land
;!p.multi=[0,1,2]
loadct2,0
ncb=8640
nrb=4320
dx25=2.5/60.
dy25=2.5/60.
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
; ZOOM IN HERE
xylim=[20.,110.,50.,150.]    ;Asia
xlen=xylim(3)-xylim(1)
ylen=xylim(2)-xylim(0)
print,'xlen=',xlen
print,'ylen=',ylen
print,xylim
init=replicate(0.,xlen,ylen)
x = indgen(xlen)+xylim(1)
y = indgen(ylen)+xylim(0)

;; BEGINNING THE PLOT
;; ------------------

tit='CATCHMENT DEFINITIONS'
fmt='(i3)'
levp=[0,99999,199999,299999,399999,499999,599999,699999,799999,899999,999999]
colors=[16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,128,130,132]
;
sp,1,1,/f,/h,/color,/land
contour,init,x,y,title=tit,xrange=[x(0),x(xlen-1)+1.],yrange=[y(0),y(ylen-1)+1.],xstyle=1,ystyle=1

;; READING AND PLOTTING STEP BY STEP
;; ---------------------------------

cc=0

pfc=long(0.)
pfc1=long(0.)
pfcl=long(0.)
pfcr=long(0.)

  cumar=0.	
  cat25=lonarr(ncb)
  cat25p=lonarr(ncb)

  file15=gfile+'.rst'
 idum=0l
 openr,15,file15;,/xdr
  for i=0,nrb-1 do begin

	readu,15,idum,cat25,idum
	yu=y1b+i*dy25+dy25/2.
	yl=y1b+i*dy25-dy25/2.

	for j=0,8160 do begin ;ASIA
	pfc=cat25(j)
	pfc1=cat25(j)
	if(j ne 0)then pfcl=cat25(j-1)
	if(j ne ncb-1)then pfcr=cat25(j+1)

	if(pfc ne 0) then begin

	if(i eq 0)then cat25p(j)=pfc
	xl=x1b+j*dx25-dx25/2.
	xr=x1b+j*dx25+dx25/2.
	xx=fltarr(5)
	yy=fltarr(5)
	xx=[xl,xl,xr,xr,xl]
	yy=[yu,yl,yl,yu,yu]
	if((xl ge min(x)) and (xr le max(x)+1.) and (yl ge min(y)) and (yu le max(y)+1.)) then begin

		ccid=strtrim((string(pfc)),2)
	if(pfc le 9999) then begin
		n=fix(strmid(ccid,3,1))
	endif
	if(pfc gt 9999) then begin
		n=fix(strmid(ccid,4,1))
	endif
	        n = pfc mod 18

		polyfill,xx,yy,color=colors(n)
		cc=cc+1

        if(pfc ne cat25p(j)) then oplot,[xl,xr],[yl,yl]
        if(pfc ne pfcl) then oplot,[xl,xl],[yl,yu]
        if(pfc ne pfcr) then oplot,[xr,xr],[yl,yu]
	endif
	endif

	cat25p(j)=pfc
	endfor

  endfor
  close,15

xleg=-80.
yleg1=12.
close,12

fileps1=gfile+'-asia.ps'

closeps,fileps1
end

