pro makeup_sstfiles,im=im,jm=jm,year0=year0,year1=year1,lon=lon,lat=lat,corners=cr,gridname=grid$,griddir=griddir
;odir='/discover/nobackup/amolod/bcs/144x91/' & oim=144 & ojm=91
odir='/discover/nobackup/projects/gmao/share/dao_ops/fvInput/g5gcm/bcs/realtime/SST/360x180/' & oim=144 & ojm=91
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
  print,' ===== Using input lon and lat m'
endelse


 seq='simple1_'
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






close,1,2,3,4

openw,1,griddir+f1,/f77
openw,2,griddir+f2,/f77
openw,3,griddir+f3,/f77
openw,4,griddir+f4,/f77

;******************* nunber of points in grid*********************
latbegin=cr(0)
latend=cr(2)
lonbegin=cr(1)
lonend=cr(3)
dlat =1.
      dlon =1.
      ibegin = ROUND(  abs(-180.-lonbegin)/dlon )
      iend = ROUND(  abs(-180.-lonend)/dlon )
      jbegin = ROUND(  abs(-90.-latbegin)/dlat )
      jend = ROUND(  abs(-90.-latend)/dlat )
      itot = (iend-ibegin+1) * (jend-jbegin+1)

;*******************************************
fsst=odir+'dataoceanfile_MERRA_sst_1971-current.360x180.LE' 

get_lun,lun
openr,lun,fsst,/f77

sst=fltarr(360,180)

date1=fltarr(6)
date2=fltarr(6)
res  =fltarr(2)
res1  =fltarr(2)
kk=1
sstav=0
for i=0,10000 do begin
readu,lun,DATE1,DATE2,RES
readu,lun,sst
if (DATE1(0) lt year0)then goto,MC
if (DATE2(0) gt year1) then goto, ME
sst1=0. 
 for ii=ibegin,iend do begin
      for jj=jbegin,jend do begin
       sst1 = sst1 + sst(ii,jj)/itot
      endfor
      endfor

sstav=sstav+sst1
 

if (DATE1(1) eq DATE2(1))then begin
 kk=kk+1


goto, MC 
endif
 sstav=sstav/kk
date1(2)=1.00000
date2(2)=1.00000
res1(0)=1.00000
res1(1)=1.00000
;print,date1,date2,res1,sstav
   writeu,1,date1,date2,res1
    writeu,1,sstav

    sstsi=sstav
    writeu,2,date1,date2,res1
    writeu,2,sstsi

kk=1
sstav=0
MC:
   endfor
ME:
    close,lun
;************************************************************
;goto, MMM
get_lun,lun
ffraci=odir+'dataoceanfile_MERRA_fraci_1971-current.360x180.LE' 
openr,lun,ffraci,/f77 
sst=fltarr(360,180)

date1=fltarr(6)
date2=fltarr(6)
res  =fltarr(2)
res1  =fltarr(2)
fraci=fltarr(360,180)
kkf=1
fracav=0
for i=0,10000 do begin
readu,lun,DATE1,DATE2,RES
readu,lun,fraci
if (DATE1(0) lt year0)then goto,MCF
if (DATE2(0) gt year1) then goto, MEF
fraci1=0. 
 for ii=ibegin,iend do begin
      for jj=jbegin,jend do begin
       fraci1 = fraci1 + fraci(ii,jj)/itot
      endfor
      endfor

fracav=fracav+fraci1
 

if (DATE1(1) eq DATE2(1))then begin
 kkf=kkf+1


goto, MCF 
endif
 fracav=fracav/kkf
date1(2)=1.00000
date2(2)=1.00000
res1(0)=1.00000
res1(1)=1.00000
   writeu,3,date1,date2,res1
    writeu,3,fracav
   


 kpar=0.
        writeu,4,date1,date2,res1
    writeu,4,kpar  
   
;    stop 
kkf=1
fracav=0
MCF:
   endfor
MEF:
    close,lun
;*******************************************

MMM:
close,1,2,3,4

return
end
