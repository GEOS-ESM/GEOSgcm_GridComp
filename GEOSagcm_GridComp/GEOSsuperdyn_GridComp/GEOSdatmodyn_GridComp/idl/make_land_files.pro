
;bcsodir='/discover/nobackup/aeichman/bcs/g1_0-2deg-ltakacs-20120224/'
bcsodir='/discover/nobackup/aeichman/barrow/bcs_20141028/'
;xdir='/discover/nobackup/aeichman/laigrn/'
xdir='/discover/nobackup/aeichman/barrow/'

;tilefile=bcsodir+'/FV_144x91_DC_360x180_DE.til'
tilefile=bcsodir+'/DC0144xPC0091_DE0360xPE0180-Pfafstetter.til'

readtilefile,f=tilefile,til=til

land=where(til(0,*) eq 100.) 
ntl=n_elements(land)

; COARE
;casename='COARE'
;itile=[ 30846] ; Ocean point at about 155 E 2 S
;tile$=strtrim( string(itile),2)
;catchname='catch_internal_rst.e19911023_21z'

;casename='newlaigreen'
;casename='barrow2013'
casename='kan-l-merra'
;itile=[11505]
;itile=[918]
itile=[2793]
 
tile$=strtrim( string(itile),2)
catchname='catch_internal_rst'

bcsxdir=xdir+casename+'/Landfiles/'
spawn,'mkdir -p '+bcsxdir

read_catch_internal,f=bcsodir+catchname,nt=ntl,cas=cas

read_vegdyn_internal,f=bcsodir+'/vegdyn_144x91_DC.dat',nt=ntl,veg=veg

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




read_land_bcs,f=bcsodir+'/green_clim_144x91_DC.data',/runthru,dataset=green,serieslength=12
;read_land_bcs,f=bcsodir+'/lai_clim_144x91_DC.data',/runthru,dataset=lai,serieslength=46
read_land_bcs,f=bcsodir+'/lai_clim_144x91_DC.data',/runthru,dataset=lai,serieslength=12


get_lun,lun
openw,lun,bcsxdir+'green.dat.TILE_'+ tile$,/f77

for n=0,13 do begin
    writeu,lun,green(n).date1,green(n).date2, [ 1.000*n_elements(itile) , 1.000 ] ; Ditto
    writeu,lun, green(n).value( itile )
endfor
close,lun

get_lun,lun
openw,lun,bcsxdir+'lai.dat.TILE_'+ tile$,/f77

;for n=0,47 do begin
for n=0,13 do begin
    writeu,lun,lai(n).date1,lai(n).date2, [ 1.000*n_elements(itile) , 1.000 ] ; Ditto
    writeu,lun, lai(n).value( itile )
endfor
close,lun



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; SST file-style
;  ... 
; date1(6),date2(6),ntl
; data(ntl)
;  ...
;       x 14


;read_land_bcs,f=bcsodir+'/visdf_144x91_DC.dat',/runthru,dataset=visdf,serieslength=12
read_land_bcs,f=bcsodir+'/visdf_144x91_DC.dat',/runthru,dataset=visdf,serieslength=46
;read_land_bcs,f=bcsodir+'/nirdf_144x91_DC.dat',/runthru,dataset=nirdf,serieslength=12
read_land_bcs,f=bcsodir+'/nirdf_144x91_DC.dat',/runthru,dataset=nirdf,serieslength=46



print , '  ---- visdf.dat ---- '
get_lun,lun
openw,lun,bcsxdir+'visdf.dat.TILE_'+ tile$  ,/f77

;for n=0,13 do begin
for n=0,47 do begin
    writeu,lun,visdf(n).date1,visdf(n).date2, [ 1.000*n_elements(itile) , 1.000 ] ; Need Dummy second element for
    writeu,lun,visdf(n).value( itile )                                           ; resolution
endfor
close,lun

get_lun,lun
openw,lun,bcsxdir+'nirdf.dat.TILE_'+ tile$,/f77

;for n=0,13 do begin
for n=0,47 do begin
    writeu,lun,nirdf(n).date1,nirdf(n).date2, [ 1.000*n_elements(itile) , 1.000 ] ; Ditto
    writeu,lun, nirdf(n).value( itile )
endfor
close,lun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print , '  ---- catch_internal_rst ---- '
name=tag_names( cas )
ntag=n_tags( cas )

get_lun,lun
openw,lun,bcsxdir+'catch_internal_rst.'+ tile$,/f77_u ; +'.ascii'

for n=1,ntag-1 do begin

    a=cas.(n)

    ss=size( a )

    if ss(0) eq 1 then begin
      writeu,lun,a(itile)        
    endif
    if ss(0) eq 2 then begin
      writeu,lun,reform( a(itile,*) )          
    endif

endfor

close,lun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print , '  ---- vegdyn_internal_rst ---- '
name=tag_names( veg )
ntag=n_tags( veg )

get_lun,lun
openw,lun,bcsxdir+'vegdyn_internal_rst.'+ tile$,/f77_u ; +'.ascii'

for n=1,ntag-1 do begin


    a=veg.(n)
    ss=size( a )

    if ss(0) eq 1 then begin
       writeu,lun,a(itile)  
    endif
    if ss(0) eq 2 then begin
      for j=0,ss(2)-1 do begin
       writeu,lun,a(itile,j)   
      endfor
    endif

endfor

close,lun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end
