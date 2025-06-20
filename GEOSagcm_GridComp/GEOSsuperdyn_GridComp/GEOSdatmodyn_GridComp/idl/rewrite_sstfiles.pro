pro makeup_sstfiles,im=im,jm=jm,year0=year0,year1=year1,lon=lon,lat=lat,gridname=grid$,griddir=gdir$

sst=fltarr(im,jm)
        
year1$=strtrim( string(year1), 2) 
year0$=strtrim( string(year0), 2) 

        im$=strtrim( string(im), 2) ;& im$=strmid( im$ , strlen(im$)-3,3)
        jm$=strtrim( string(jm), 2) ;& jm$=strmid( jm$ , strlen(jm$)-3,3)

if not keyword_set(lon) or not keyword_set(lat) then begin
  lat=fltarr(jm)
  lon=fltarr(im)

  DX = float( fix( 100 * 360./im ) / 100. )
  DY = float( fix( 100 * 180./(jm-1) ) /100. )

  for j=0,jm-1 do begin
    lat(j)=j*Dy-90.
  endfor
 
  for i=0,im-1 do begin
    lon(i)=i*Dx-180.
  endfor
endif else begin
  print,' ===== Using input lon and lat r'
endelse

seq='reynolds_'
if not keyword_set(grid$) then begin
   f1='sst_'+seq+year0$+'_'+year1$+'_PC'+im$+'x'+jm$+'-DC.dat'
   f2='sstsi_'+seq+year0$+'_'+year1$+'_PC'+im$+'x'+jm$+'-DC.dat'
   f3='fraci_'+seq+year0$+'_'+year1$+'_PC'+im$+'x'+jm$+'-DC.dat'
   f4='SEAWIFS_KPAR_mon_clim_PC'+im$+'x'+jm$+'-DC.dat'
endif else begin
   f1='sst_'+seq+year0$+'_'+year1$+'_'+grid$+'.dat'
   f2='sstsi_'+seq+year0$+'_'+year1$+'_'+grid$+'.dat'
   f3='fraci_'+seq+year0$+'_'+year1$+'_'+grid$+'.dat'
   f4='SEAWIFS_KPAR_mon_clim_'+grid$+'.dat'
endelse

for j=0,jm-1 do begin

    sst(*,j) = 270.+30.*cos( 1.25*lat(j)*!pi/180. )

endfor


sstsi=sst



fraci=sstsi*0.

is_ice = where(sstsi lt 272.)
if min(is_ice ge 0) then fraci( is_ice )=1.00

kpar = sst*0.


;;STOP

close,1,2,3,4

openw,1,gdir$+f1,/f77
openw,2,gdir$+f2,/f77
openw,3,gdir$+f3,/f77
openw,4,gdir$+f4,/f77


for yy=year0,year1 do begin
for mm=1,12 do begin

    date0=[ yy*1.0 ,       mm*1.0 ,       1.0,  0.0 , 0.0, 0.0 ]
    date1=[ yy*1.0 ,       (mm+1)*1.0 ,   1.0,  0.0 , 0.0, 0.0 ]
    res  =[ im*1.0 , jm*1.0 ]

    if ( mm eq 12 ) then begin
       date1 = [ (yy+1)*1.0 , 1.0 , 1.0, 0.0, 0.0, 0.0 ]
    endif
 


    writeu,1,date0,date1,res
    writeu,1,sst

    writeu,2,date0,date1,res
    writeu,2,sstsi

    writeu,3,date0,date1,res
    writeu,3,fraci



endfor
endfor


yy=0
for mm=12,12 do begin
    date0=[ yy*1.0 ,       mm*1.0 ,       1.0,  0.0 , 0.0, 0.0 ]
    date1=[ yy*1.0 ,       (mm+1)*1.0 ,   1.0,  0.0 , 0.0, 0.0 ]
    res  =[ im*1.0 , jm*1.0 ]

    if ( mm eq 12 ) then begin
       date1 = [ (yy+1)*1.0 , 1.0 , 1.0, 0.0, 0.0, 0.0 ]
    endif
 
    writeu,4,date0,date1,res
    writeu,4,kpar
endfor

yy=1
for mm=1,12 do begin
    date0=[ yy*1.0 ,       mm*1.0 ,       1.0,  0.0 , 0.0, 0.0 ]
    date1=[ yy*1.0 ,       (mm+1)*1.0 ,   1.0,  0.0 , 0.0, 0.0 ]
    res  =[ im*1.0 , jm*1.0 ]

    if ( mm eq 12 ) then begin
       date1 = [ (yy+1)*1.0 , 1.0 , 1.0, 0.0, 0.0, 0.0 ]
    endif
 
    writeu,4,date0,date1,res
    writeu,4,kpar
endfor

yy=2
for mm=1,1 do begin
    date0=[ yy*1.0 ,       mm*1.0 ,       1.0,  0.0 , 0.0, 0.0 ]
    date1=[ yy*1.0 ,       (mm+1)*1.0 ,   1.0,  0.0 , 0.0, 0.0 ]
    res  =[ im*1.0 , jm*1.0 ]

    if ( mm eq 12 ) then begin
       date1 = [ (yy+1)*1.0 , 1.0 , 1.0, 0.0, 0.0, 0.0 ]
    endif
 
    writeu,4,date0,date1,res
    writeu,4,kpar
endfor




close,1,2,3,4




return
end






