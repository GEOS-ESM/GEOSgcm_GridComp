pro mktile,im=im,jm=jm,gridtype=gtyp$,corners=cr,tilefile=tile$,gdirbase=gdirbase,griddir=griddir,gridname=grid$,lat=lat,lon=lon,ntype=ntype

;   113307
;   2
;PC144x91-DC
;   144
;   91
;PE360x180-DE
;   360
;   180

case gtyp$ of

'PC': begin
    lat=fltarr(jm)
    lon=fltarr(im)
    DX = float( fix( 100 * 360./im ) / 100. )
    DY = float( fix( 100 * 180./(jm-1) ) /100. )

    for j=0,jm-1 do begin
       lat(j)=j*DY-90.
    endfor

    for i=0,im-1 do begin
        lon(i)=i*DX-180.
    endfor
    im$=strtrim( string(im), 2) 
    jm$=strtrim( string(jm), 2)
    grid$ = 'PC'+im$+'x'+jm$+'-DC'
    spawn,'mkdir -p '+griddir

end
'XY': begin
                        ; A-grid locations

    lat=fltarr(jm)
    lon=fltarr(im)

     cr = fix( cr )

     ns$=['S','N']
     ew$=['W','E']
     crs$=strtrim( string(abs(cr(0))) , 2) + ns$( fix( cr(0) gt 0 ) )
     crn$=strtrim( string(abs(cr(2))) , 2) + ns$( fix( cr(2) gt 0 ) )
     crw$=strtrim( string(abs(cr(1))) , 2)  + ew$( fix( cr(1) gt 0 ) )
     cre$=strtrim( string(abs(cr(3))) , 2)  + ew$( fix( cr(3) gt 0 ) )
     
    ; Corners cr = [ LAT_SOUTH, LON_WEST,LAT_NORTH, LON_EAST ] 
    DOMY= cr(2)-cr(0)
    DOMX= cr(3)-cr(1)

    DX = float( fix( 100 * DOMX/(im) ) / 100. )
    DY = float( fix( 100 * DOMY/(jm) ) /100. )

    for j=0,jm-1 do begin
        lat(j) = cr(0) + 0.5*DY + j*DY
    endfor

    for i=0,im-1 do begin
        lon(i) = cr(1) + 0.5*DX + i*DX
    endfor
    im$=strtrim( string(im), 2) 
    jm$=strtrim( string(jm), 2)
    grid$ = 'XY'+im$+'x'+jm$+'-C_'+crs$+'_'+crw$+'_'+crn$+'_'+cre$
    spawn,'mkdir -p '+griddir

end

endcase

;2345678912345678912345678901234567890123451234512345678901234561234567812345678123456781234567890123456123456781234567890123
;       0 -1800360  179.5000   89.5000    1   91  0.400000000000       1     360     180  1.000000000000       1     107.8957
;   i9       i9      f10.4      f10.4   i4  i5        f16.10        i8     i8       i8       f16.10          i8        f13.4
;  type    ???        lon       lat      i1   j1     frac1            nt1     i2     j2      frac2            nt2       ????   


;;      printf,1,form='("PC144x91-DC")'

frm='( i5 , 2x,a7,    f10.4 ,    f10.4 ,  i4,  i5,     f16.10 ,      i8,    i8,      i8,      f16.10,         i8,      f13.4 )' 

ntile=1
close,1

tile$='tile.data_'+'simple1_'+ grid$
openw,1,griddir+tile$

printf,1,'  '+strtrim(string(im*jm*1L),2)
printf,1,form='(i3)',2
printf,1,grid$  ; Atmos grid name
printf,1,form='(i6)',im
printf,1,form='(i6)',jm
printf,1,grid$  ; Ocean grid name
printf,1,form='(i6)',im
printf,1,form='(i6)',jm

for j=0,jm-1 do begin
    for i=0,im-1 do begin
      
        i$='0000'+strtrim( string(i+1), 2) & i$=strmid( i$ , strlen(i$)-3,3)
        j$='0000'+strtrim( string(j+1), 2) & j$=strmid( j$ , strlen(j$)-3,3)

        ij$='-'+i$+j$        
       
        printf,1,form=frm, ntype,ij$,lon(i),lat(j),i+1,j+1,1.000,ntile,i+1,j+1,1.000,ntile,123.456 

        ntile=ntile+1

    endfor 
endfor
close,1

return
end
