
;odir='/discover/nobackup/aeichman/bcs/f2_1-2deg-amolod-20120224/' & oim=144 & ojm=91
;odir='/discover/nobackup/aeichman/barrow/bcs/144x91/' & oim=144 & ojm=91
odir='/discover/nobackup/aeichman/barrow/bcs_20141028/' & oim=144 & ojm=91
;gdirbase='/discover/nobackup/aeichman/laigrn/' & oim=144 & ojm=91
gdirbase='/discover/nobackup/aeichman/barrow/' & oim=144 & ojm=91

im=1
jm=1
;sst_impose=300.

gtyp$='XY'
;gtyp$='PC'


; SCSMEX -MERRA
;cr=[18,115,23,120] 
;casename='merra_scmx'
;im=1
;jm=1
;ntype=0
; SCSMEX
;cr=[70,-158,72,-156]
;casename='barrow2013'
im=1
jm=1
ntype=0

;kan-l-merra

cr=[66,-52,68,-50]
casename='kan-l-merra'


; MERRA_NAMMA
;cr=[15,-25,17,-23]
;casename='merra_namma'
;im=1
;jm=1
;ntype=0

; ARM SGP
;cr=[34,-100,38,-95] 
;casename='arm_97jul'
;im=1
;jm=1
;ntype=100

; COARE
;cr=[-5,152,2,158]
;casename='COARE'
;im=1
;jm=1
;sst_impose=300.
;ntype=0

; TWP-SVETA
;cr=[-14,129,-10,133] 
;casename='merra_twp'
;im=1
;jm=1
;sst_impose=300.
;ntype=100

; TRMM-SVETA
;cr=[-12,-63,-10,-61] 
;casename='merra_trmm'
;im=1
;jm=1
;ntype=100
; ARMSGP
;cr=[34,-100,39,-95] 
;casename='merra_armsgp'
;im=1
;jm=1
;ntype=100
; MERRA-  KWAJEX
;cr=[7,165,12,170] 
;casename='merra_kwjx'
;im=1
;jm=1
;ntype=0


;year0=1990 & year1=2006
year0=1990 & year1=2013

griddir=gdirbase+casename+'/'

mktile,im=im,jm=jm,gridtype=gtyp$,corners=cr,tilefile=tile$,gdirbase=gdirbase,griddir=griddir,gridname=grid$,lat=lat,lon=lon,ntype=ntype

print,'makeup_sstfiles' 
makeup_sstfiles,im=im,jm=jm,year0=year0,year1=year1,lon=lon,lat=lat,corners=cr,gridname=grid$,griddir=griddir

ogr={lon:-180.0+findgen(oim)*360. /(1.*oim)  , lat: -90.0+findgen(ojm)*180./(1.*ojm-1.) }
xgr={lon:lon                                 , lat:lat }

rewrite_topofiles,f='topo_DYN_ave_144x91_DC.data',/writefile,scalefactor=1.0,odir=odir,xdir=griddir  $
                     ,gridname=grid$,ogr=ogr,xgr=xgr
rewrite_topofiles,f='topo_DYN_ave_144x91_DC.data',/writefile,odir=odir,xdir=griddir  $
                     ,gridname=grid$,ogr=ogr,xgr=xgr,/aqua
rewrite_topofiles,f='topo_GWD_var_144x91_DC.data',/writefile,scalefactor=1.0,odir=odir,xdir=griddir  $
                     ,gridname=grid$,ogr=ogr,xgr=xgr
rewrite_topofiles,f='topo_TRB_var_144x91_DC.data',/writefile,scalefactor=1.0,odir=odir,xdir=griddir  $
                     ,gridname=grid$,ogr=ogr,xgr=xgr


rewrite_fv_internal,f='fvcore_internal_rst',/writefile,odir=odir,xdir=griddir $
                       ,gridname=grid$,ogr=ogr,xgr=xgr,lm=lm

rewrite_moist_internal,f='moist_internal_rst',/writefile,odir=odir,xdir=griddir $
                       ,gridname=grid$,ogr=ogr,xgr=xgr,lm=lm

write_datmodyn_internal,f='fvcore_internal_rst',/writefile,odir=odir,xdir=griddir $
                       ,gridname=grid$,ogr=ogr,xgr=xgr,lm=lm

end
