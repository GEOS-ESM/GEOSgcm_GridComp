days1 = ['20030103/','20080103/','19990103/','20000103/','20010103/','20020103/','20090103/','20100103/','20110103/']         
days2 = ['20030113/','20080113/','19990113/','20000113/','20010113/','20020113/','20090113/','20100113/','20110113/'] 
days3 = ['20030124/','20080124/','19990124/','20000124/','20010124/','20020124/','20090124/','20100124/','20110124/'] 
days4 = ['20030203/','20080203/','19990203/','20000203/','20010203/','20020203/','20090203/','20100203/','20110203/'] 
days5 = ['20030213/','20080213/','19990213/','20000213/','20010213/','20020213/','20090213/','20100213/','20110213/'] 
days6 = ['20030221/','20080222/','19990221/','20000222/','20010221/','20020221/','20090221/','20100221/','20110221/'] 
days7 = ['20030303/','20080303/','19990303/','20000303/','20010303/','20020303/','20090303/','20100303/','20110303/'] 
days8 = ['20030313/','20080313/','19990313/','20000313/','20010313/','20020313/','20090313/','20100313/','20110313/'] 
days9 = ['20030324/','20080324/','19990324/','20000324/','20010324/','20020324/','20090324/','20100324/','20110324/'] 
days10= ['20030403/','20080403/','19990403/','20000403/','20010403/','20020403/','20090403/','20100403/','20110403/'] 
days11= ['20030413/','20080413/','19990413/','20000413/','20010413/','20020413/','20090413/','20100413/','20110413/'] 
days12= ['20080423/','19990423/','20000423/','20010423/','20020423/','20090423/','20100423/','20110423/'] 
days13= ['20080503/','19990503/','20000503/','20010503/','20020503/','20090503/','20100503/','20110503/'] 
days14= ['20080513/','19990513/','20000513/','20010513/','20020513/','20090513/','20100513/','20110513/'] 
days15= ['20030524/','20080524/','19990524/','20000524/','20010524/','20020524/','20090524/','20100524/','20110524/'] 
days16= ['20080603/','19990603/','20000603/','20010603/','20020603/','20090603/','20100603/','20110603/'] 
days17= ['20030613/','20080613/','19990613/','20000613/','20010613/','20020613/','20090613/','20110613/'] 
days18= ['20030623/','20080623/','19990623/','20000623/','20010623/','20020623/','20090623/','20100623/','20110623/'] 
days19= ['20080703/','19990703/','20000703/','20010703/','20020703/','20090703/','20100703/','20110703/'] 
days20= ['20030713/','20080713/','19990713/','20000713/','20010713/','20020713/','20090713/','20100713/','20110713/']
days21= ['20030724/','20080724/','19990724/','20000724/','20010724/','20020724/','20090724/','20100724/','20110724/'] 
days22= ['20030803/','20080803/','19990803/','20000803/','20010803/','20020803/','20090803/','20100803/','20110803/'] 
days23= ['20030813/','20080813/','19990813/','20000813/','20010813/','20020813/','20090813/','20100813/','20110813/'] 
days24= ['20030824/','20080824/','19990824/','20000824/','20010824/','20020824/','20090824/','20100824/','20110824/']
days25= ['20030903/','20080903/','19990903/','20000903/','20010903/','20020903/','20090903/','20100903/','20110903/'] 
days26= ['20030913/','20080913/','19990913/','20000913/','20010913/','20020913/','20090913/','20100913/','20110913/'] 
days27= ['20030923/','20080923/','19990923/','20000923/','20010923/','20020923/','20090923/','20100923/','20110923/'] 
days28= ['20031003/','20081003/','19991003/','20001003/','20011003/','20021003/','20091003/','20101003/','20111003/'] 
days29= ['20031013/','20081013/','19991013/','20001013/','20021013/','20091013/','20101013/','20111013/'] 
days30= ['20031024/','20081024/','19991024/','20001024/','20011024/','20021024/','20091024/','20101024/','20111024/']
days31= ['20031103/','20081103/','19991103/','20001103/','20021103/','20091103/','20101103/','20111103/'] 
days32= ['20031113/','20081113/','19991113/','20001113/','20021113/','20091113/','20101113/','20111113/'] 
days33= ['20031123/','20081123/','19991123/','20001123/','20021123/','20091123/','20101123/','20111123/'] 
days34= ['20031203/','20081203/','19991203/','20001203/','20021203/','20091203/','20101203/','20111203/'] 
days35= ['20031213/','20081213/','19991213/','20001213/','20021213/','20091213/','20101213/','20111213/'] 
days36= ['20071224/','20021224/','19981224/','19991224/','20001224/','20011224/','20091224/','20101224/']

mdays = [31,28,31,30,31,30,31,31,30,31,30,31]

modis_dates = ['01.01','01.09','01.25','02.02','02.10','02.18','03.06','03.14','03.22','04.07','04.15','04.23', $
               '05.01','05.09','05.25','06.02','06.10','06.26','07.04','07.12','07.20','08.05','08.13','08.21', $
               '09.06','09.14','09.22','09.30','10.16','10.24','11.01','11.09','11.25','12.03','12.11','12.27']

;    mytime=timegen(start=julday(1,1,1950),final=julday(1,1,1950)+reftime(t),unit='d')
;    caldat,mytime,month,day,year
;    i=n_elements(mytime)-1
;    start_date(t)=string(month[i],day[i],year[i],format='(i2.2,".",i2.2,".",i4.4)')

stop
                      
nc = 40320l
nr = 20160l
nx = nc/36
ny = nr/18
scale_factor = 0.01
undef = -9999

ncb=360
nrb=180
dx25=1.
dy25=1.
nc_low=ncb
nr_low=nrb
resol=strtrim((string(ncb)),2) +'x'+strtrim((string(nrb)),2)
x1b=-180. + dx25/2.
xeb= x1b +(ncb-1)*dx25
y1b= -90. + dy25/2.
yeb= y1b +(nrb-1)*dy25

dateline='DE'

pathout='/land/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1/'
tilefile=pathout + 'FV_'+resol +'_'+dateline+'_360x180_DE_NO_TINY.til'
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


for d = 1,36 do begin

    if (d eq 1 ) then days = days1 
    if (d eq 2 ) then days = days2 
    if (d eq 3 ) then days = days3 
    if (d eq 4 ) then days = days4 
    if (d eq 5 ) then days = days5 
    if (d eq 6 ) then days = days6 
    if (d eq 7 ) then days = days7 
    if (d eq 8 ) then days = days8 
    if (d eq 9 ) then days = days9 
    if (d eq 10) then days = days10
    if (d eq 11) then days = days11
    if (d eq 12) then days = days12
    if (d eq 13) then days = days13
    if (d eq 14) then days = days14
    if (d eq 15) then days = days15
    if (d eq 16) then days = days16
    if (d eq 17) then days = days17
    if (d eq 18) then days = days18
    if (d eq 19) then days = days19
    if (d eq 20) then days = days20
    if (d eq 21) then days = days21
    if (d eq 22) then days = days22
    if (d eq 23) then days = days23
    if (d eq 24) then days = days24
    if (d eq 25) then days = days25
    if (d eq 26) then days = days26
    if (d eq 27) then days = days27
    if (d eq 28) then days = days28
    if (d eq 29) then days = days29
    if (d eq 30) then days = days30
    if (d eq 31) then days = days31
    if (d eq 32) then days = days32
    if (d eq 33) then days = days33
    if (d eq 34) then days = days34
    if (d eq 35) then days = days35
    if (d eq 36) then days = days36

    mmdd = strmid(days(0),4,4)
    modis_mmdd = modis_dates (d -1)
    restore,  '/l_data/MODIS/LAI/YYYY.' + MODIS_MMDD +'/MOD15A2.global_2.5.idl'
    lai_modis = lai
    lai = fltarr (nc,nr)
    
    for jx =1,36 do begin

        j1 = (jx -1)*1120l
        j2 = jx*1120l - 1l

        for ix =1,36 do begin
            fname = '/l_data/GEOLAND2-LAI/LAI-CLIM/' + 'geoland2_lai_clim.'+'H'+string(ix,'(i2.2)')+'V'+string(jx,'(i2.2)') + '.' + mmdd
            if file_test(fname) then begin
;                print,fname
                openr, lun,fname,/get_lun,/F77_UNFORMATTED   
                ii1 =0l
                ii2 =0l
                jj1 =0l
                jj2 =0l
                sf  =0.
                uf  =0
                NX =0L
                NY=0L
                NC=0L
                NR=0L
                
                readu,lun,ii1,ii2,jj1,jj2,sf,uf,nx,ny,nc,nr
;                print, ii1,ii2,jj1,jj2,sf,uf,nx,ny,nc,nr
                ii1 = ii1 - 1
                ii2 = ii2 - 1
                jj1 = jj1 - 1
                jj2 = jj2 - 1

                for j = jj1,jj2 do begin
                    strip = intarr (nx)
                    readu,lun,strip
                    if (min(strip) eq uf) then strip (where (strip eq undef)) = !VALUES.F_NAN
                    strip = strip*sf
                    lai(ii1:ii2,j) = strip(*) 
                endfor 
                close,lun
                free_lun, lun                                             
            endif
        endfor
    endfor

sp,1,3,/f,/h,/color,landscape=0
loadct2,0

k1 = nc/nc_low
k2 = k1 -1 

M1 = 8640/nc_low
M2 = M1 -1 

lai_low_geo = fltarr (nc_low,nr_low)
lai_low_mod = fltarr (nc_low,nr_low)
lai_low_geo (*,*)  = !VALUES.F_NAN
lai_low_mod (*,*)  = !VALUES.F_NAN

for j = 0,nr_low -1 do begin
    for i =0,nc_low-1 do begin
        if(mean (lai(i*long(k1):i*long(k1)+k2,j*long(k1):j*long(k1)+k2),/nan) gt 0.) then $
        lai_low_geo (i,j) = mean (lai(i*long(k1):i*long(k1)+k2,j*long(k1):j*long(k1)+k2),/nan)

        lai_low_mod (i,j) = mean (lai_modis(i*long(M1):i*long(M1)+M2,j*long(M1):j*long(M1)+M2),/nan)
        
    endfor
endfor
upval=10.
lwval=0.
colors=90+indgen(11)
upval =10.0
lwval =0.
levp = [lwval,lwval+(upval-lwval)/10+indgen(10)*(upval-lwval)/10.]

contour,lai_low_mod,title='MODIS on ' + modis_mmdd,c_colors=colors,levels=levp,xstyle=4,ystyle=4,/CELL_FILL
contour,lai_low_geo,title='GEOSLAND2 on ' + mmdd,c_colors=colors,levels=levp,xstyle=4,ystyle=4,/CELL_FILL

!P.position=[0.15, 0.35+0.01, 0.95, 0.365+0.01]

alpha=fltarr(n_elements(levp),2)
  alpha(*,0)=levp
  alpha(*,1)=levp
  h=[0,1]
contour,alpha,levp,h,levels=levp,c_colors=colors,/fill,/xstyle,/ystyle, $
    /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levp),max(levp)], $
    xtitle=' ', color=0, $
    C_charsize=0.7, charsize=2. ,xticks=10
contour,alpha,levp,h,levels=levp,color=0,/overplot



!P.position=0
colors2=indgen(21) +80
upval=5.
lwval=-5.

levp=[lwval,lwval+(upval-lwval)/20+indgen(19)*(upval-lwval)/20.,upval]

lai_low_mod = lai_low_mod - lai_low_geo
contour,lai_low_mod,title='MODIS - GEOLAND2 ' ,c_colors=colors2,levels=levp,xstyle=4,ystyle=4,/CELL_FILL

!P.position=[0.15, 0.+0.01, 0.95, 0.015+0.01]

alpha=fltarr(n_elements(levp),2)
  alpha(*,0)=levp
  alpha(*,1)=levp
  h=[0,1]
contour,alpha,levp,h,levels=levp,c_colors=colors2,/fill,/xstyle,/ystyle, $
    /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levp),max(levp)], $
    xtitle=' ', color=0, $
    C_charsize=0.7, charsize=2. ,xticks=10
contour,alpha,levp,h,levels=levp,color=0,/overplot



!P.position=0

closeps,'/l_data/GEOLAND2-LAI/LAI-CLIM/geoland2_lai_clim.'+ mmdd + '.ps'
endfor

end
