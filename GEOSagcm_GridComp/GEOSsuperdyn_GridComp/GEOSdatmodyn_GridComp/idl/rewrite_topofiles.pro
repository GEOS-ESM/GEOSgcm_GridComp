pro rewrite_topofiles,f=f,im=im,jm=jm,xim=xim,xjm=xjm,writefile=writefile,scalefactor=scale,aquaplanet=flatten,odir=odir,xdir=xdir $
                     ,gridname=grid$,ogr=ogr,xgr=xgr


fx=f+'_'+grid$
if keyword_set(flatten) then fx=fx+'_aqua'
close,1,2
openr,1,/f77,odir+f

 if keyword_set( writefile) then openw,2,/f77,xdir+fx

if keyword_set(scale) then begin
   if scale le 1.e-5 then begin 
       ScaleTopo=0.0
   endif else begin
       ScaleTopo=scale
   endelse
endif else begin
   ScaleTopo = 1.0
endelse

if keyword_set(flatten) then ScaleTopo=0.0

if not keyword_set(ogr) then ogr={lon:findgen(im)*360. /(1.*im)  , lat:findgen(jm)*180. /(1.*jm) }
if not keyword_set(xgr) then xgr={lon:findgen(xim)*360./(1.*xim) , lat:findgen(xjm)*180./(1.*xjm) }

a = fltarr( n_elements(ogr.lon) ,  n_elements(ogr.lat) )

    readu,1,a ; & a(*,0)=0.0
    ax = ScaleTopo * int2d( a , from=ogr, to=xgr, /noshifts )
    if keyword_set( writefile) then writeu,2,ax


close,1,2

return
end
