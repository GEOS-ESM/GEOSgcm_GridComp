loadct2,0
colors=80+indgen(21)
upval =8.0
lwval =0.
levp = [lwval,lwval+(upval-lwval)/20+indgen(20)*(upval-lwval)/20.]

i_sib=8640
j_sib=4320

x25=indgen(i_sib)*2.5/60. -180.
y25=indgen(j_sib)*2.5/60. -90.

i_gswp2=360
j_gswp2=180

x25=indgen(i_sib)*2.5/60. -180. + 1.25/60.
y25=indgen(j_sib)*2.5/60. -90. + 1.25/60.

xgswp2=indgen(i_gswp2) -180. + 0.5
ygswp2=indgen(j_gswp2) -90. + 0.5

soil_high = fltarr(i_sib,j_sib)
soil_low  = fltarr(i_gswp2,j_gswp2)

openr,1,'/l_data/clsm_para_generation/data/soil_depth_2.5.rst',/f77_unformatted

rdum = fltarr(i_sib)
for j=0,j_sib -1 do begin
    readu,1,rdum
    soil_high(*,j) = rdum(*)
endfor

rdum = fltarr(i_gswp2)
for j=0,j_gswp2 -1 do begin
    readu,1,rdum
    soil_low(*,j) = rdum(*)
endfor
close,1
sp,1,1,/f,/h,/color,/land

laxis=x25(60*24:100*24+24)
yhigh=soil_high(60*24:100*24+24,125*24+12)
plot,laxis,yhigh,yrange=[0.,5.]
for i =-12,12,2 do begin
oplot,laxis,soil_high(60*24:100*24+24,125*24+12+i),color=90+i
endfor
xlow=xgswp2(60:100)
ylow=soil_low(60:100,125)
oplot,xlow,ylow,psym=4
ylow=soil_low(60:100,126)
oplot,xlow,ylow,psym=5
;for i =60,100 do begin
;    ylow(
;endfor

sp,1,1,/f,/h,/color,/land

contour,soil_high,x25,y25,c_colors=colors,levels=levp,c_label=levp,/fill,title='High'

contour,soil_low,xgswp2,ygswp2,c_colors=colors,levels=levp,c_label=levp,/fill,title='Low'

closeps,'soil.ps'

end
