pro get_lat_lon4tils

; code to derive lat/lons for each MOD10A1 grid point. The code 
; loops over MODIS files (hereinafter refered to as tiles), its 
; grid boxes and calculates lat/lon values using the same approach
; as given in C function provided by MODIS team (publically available). 
; The lat/lon values for each gridbox are written out in 10-degx10-deg files
; (to match MODIS files; each 2400x2400 elements), as well as in two 
; global files holding lat and lon each. The code will print 
; out even the tiles (and gridboxes) with no valid data (miss val:-999.)

; to run in terminal window: 
; type 'idl', then 'r. get_lat_lon4tils', then 'get_lat_lon4tils'
; requirements: a subdirectory 'MODIS_lat_lon' in the local space 
; created: March 2022 Biljana Orescanin SSAI@NASA

mis_val=-999.

lat_stiched=make_array(2400l*36, 2400l*18, value=mis_val,/float) 
lon_stiched=make_array(2400l*36, 2400l*18, value=mis_val,/float) ; Array[864, 432]

; Set Map projection 
sinusMap = MAP_PROJ_INIT('Sinusoidal', DATUM=8, CENTER_LAT=0., CENTER_LON=0,SPHERE_RADIUS=6371007.181)

; loop over MOIDS tiles
for ivtil=0,18-1 do begin
  for ihtil=0,36-1 do begin

    ; get 2-digit strings for vert and horiz tile/file numbering
    v_str=strmid('0'+strtrim(ivtil,2),1,2,/reverse)
    h_str=strmid('0'+strtrim(ihtil,2),1,2,/reverse)

    ; declare 2D arays to store lat/lons for the current file
    lat2d=make_array(2400l, 2400l, value=mis_val,/float) ;  fltarr(2400,2400)
    lon2d=make_array(2400l, 2400l, value=mis_val,/float) ; 

    ; loop over grid boxes within the current file
    for ivgrid=0,2400-1 do begin
      for ihgrid=0,2400-1 do begin

         ; get lats and lons from UV coordinates (requires the map projection loaded above)
         xOrigin=-20015109.354d +((ihtil)*2400l+ihgrid)*463.31271653d
         yOrigin= 10007554.677d -((ivtil)*2400l+ivgrid)*463.31271653d
         lonlat = MAP_PROJ_INVERSE([xOrigin], [yOrigin], MAP_STRUCTURE=sinusMap)      
         lat    =lonlat[1]
         lon    =lonlat[0]

         if lat eq lat then begin
           lat2d(ihgrid,ivgrid)=lat
           lon2d(ihgrid,ivgrid)=lon
           i_ind=(ihtil)*2400l+ihgrid
           j_ind=(ivtil)*2400l+ivgrid
           lat_stiched(i_ind,j_ind)=lat
           lon_stiched(i_ind,j_ind)=lon
         endif

      endfor ; ihgrid
    endfor ; ivgrid

    ;write 2D lat/lons into binary file
    openw, 10, 'MODIS_lat_lon/MODIS_hdres_lat_lon_v'+v_str+'_h'+h_str+'.bin.gz', /compress
    writeu,10, lat2D,lon2D
    close, 10

    ; place the current tile in global gird
    lat_stiched((ihtil)*2400l:(ihtil)*2400l+2399,(ivtil)*2400l:(ivtil)*2400l+2399l)=lat2d
    lon_stiched((ihtil)*2400l:(ihtil)*2400l+2399,(ivtil)*2400l:(ivtil)*2400l+2399l)=lon2d

  endfor ; ihtil
endfor ; ivril

; write stiched lat/lons into a text file
openw , 20, 'MODIS_lat_lon/MODIS_lat_stiched_hdres_global.txt.gz', /compress
printf, 20, lat_stiched
close , 20

openw , 30, 'MODIS_lat_lon/MODIS_lon_stiched_hdres_global.txt.gz', /compress
printf, 30, lon_stiched
close , 30

end
