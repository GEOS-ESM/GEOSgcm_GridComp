;==============================================================================
; STEP5: Build 30-arcsec (nc=43200 x nr=21600) global soil-property arrays,
;        fill local holes horizontally, and write 36x18 BCS tiles to NetCDF.
;
; INPUTS:
;   - HWSDv2/NGDC/STATSGO soil fields (F77_UNFORMATTED row strips) from `path`:
;       clay_top_30sec.dat, sand_top_30sec.dat, oc_top_30sec.dat,
;       clay_sub_30sec.dat, sand_sub_30sec.dat, oc_sub_30sec.dat
;   - Gravel fraction (top/sub) from `path3`:
;       gravel_top_30sec.dat, gravel_sub_30sec.dat
;   - Source mask from `path`:
;       data_sources.msk
;   - Legacy tile-domain index from path2: tile_id.rst (ONLY file used from path2)
;
; TILE_ID USAGE (legacy):
;   tileid(i,j) is used ONLY to:
;     (1) identify pixels outside any valid 36×18 BCS tile (tileid==0),
;     (2) restrict hole-filling to valid tile domain (1..maxcat),
;     (3) suppress writing completely empty tiles (skip tiles with max(tileid)<1).
;   Soil validity and merging are driven by data presence (HWSDv2 values),
;   not by tileid.
;
; OUTPUT:
;   - NetCDF tile files:
;       SoilProperties_HxxVyy.nc
;     containing Mask + Clay/Sand/OC (top/sub) + Gravel
; NOTE:
; This workflow does NOT use a legacy HWSD land/sea mask. Validity and gap-filling
; are driven by data presence (clay/sand) plus tile-domain gating via tile_id.rst.
;============================================================================== 

nc=43200l
nr=21600l
nx = nc/36
ny = nr/18
dxy = 360./nc
sf  = .01
un='%'
uf = -9999

  ;HWSDv2+NGDC+STATSGO2+AFSIS
   ;path='/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/HWSDv2-NGDC-STATSGO-AFSIS/'

  ;HWSDv2+NGDC+STATSGO2
  path='/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/HWSDv2-NGDC-STATSGO/'

; for files I don't have on discover I need path2 so I can run IDL code
   path2 ='/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/'
   path3 = '/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/MERGED-DATA/'

 file0 = path2+'tile_id.rst'   

;HWSDv2+NGDC+STATSGO2
;path='/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/HWSDv2-NGDC-STATSGO/'
 file1 = 'clay_top_30sec.dat' 
 file2 = 'sand_top_30sec.dat' 
 file3 = 'oc_top_30sec.dat' 
 file4 = 'clay_sub_30sec.dat' 
 file5 = 'sand_sub_30sec.dat' 
 file6 = 'oc_sub_30sec.dat' 
;------------------------------------------------------------------------------
; Coarse fragments / Gravel handling (HWSDv2)
;------------------------------------------------------------------------------
; HWSDv2 provides a native coarse fragments (Gravel) field (coarse fragments, D2).
; This replaces the legacy HWSDv1.2 gravel proxy that was derived from textures.
;
; Inputs used here:
;   - We ingest coarse fragments via prebuilt 30-arcsec intermediates (derived
;     from HWSDv2 D2):
;       gravel_top_30sec.dat and gravel_sub_30sec.dat
;
; Implementation notes (legacy STEP5 compatibility):
;   - We read two gravel inputs (top/sub) as integer percent (0..100) with uf=-9999.
;   - To match legacy STEP5 behavior, we combine top/sub into one Gravel value using
;     the historical weighting: 0.3*top + 0.7*sub.
;   - The output Gravel is stored as short with ScaleFactor=0.01, so we convert
;     integer percent -> integer (percent*100) before writing.
;
; Historical (HWSDv1.2) inputs (no longer used):
;   ;file7 = path2+'HWSDv1.21/top_grav.bin'
;   ;file8 = path2+'HWSDv1.21/sub_grav.bin'
;------------------------------------------------------------------------------

  file7 = path3+'gravel_top_30sec.dat'
  file8 = path3+'gravel_sub_30sec.dat'

strip0 = lonarr(nc)
strip1 = fltarr(nc)
strip2 = fltarr(nc)
strip3 = fltarr(nc)
strip4 = fltarr(nc)
strip5 = fltarr(nc)
strip6 = fltarr(nc)
;  HWSDv1.2 ;strip7 = fltarr(nc)
;  HWSDv1.2 ;strip8 = fltarr(nc)
strip7 = intarr(nc)
strip8 = intarr(nc)

;------------------------------------------------------------------------------
; Mask concepts used in this workflow
;------------------------------------------------------------------------------
; 1) Legacy HWSD land mask (e.g., hwsd_mask.bin):
;    Used in older pipelines. NOT used in HWSDv2 processing here.
;
; 2) Data-driven validity mask (this script: `mask`):
;    mask(i,j)=1 when sufficient soil data are present (top or sub texture valid).
;    mask(i,j)=0 for holes inside valid tiles (eligible for gap filling).
;    mask(i,j)=-9999 for out-of-domain pixels (tileid==0 and no data).
;
; 3) NetCDF "Mask" variable written to tiles:
;    This is the final `mask` subset written to each NetCDF tile file.
;    It encodes data validity / fill state, NOT a scientific land/sea mask.
;------------------------------------------------------------------------------
mask  = intarr(nc,nr)
tileid= lonarr(nc,nr)

data1 = intarr(nc,nr)
data2 = intarr(nc,nr)
data3 = intarr(nc,nr)
data4 = intarr(nc,nr)
data5 = intarr(nc,nr)
data6 = intarr(nc,nr)
data7 = intarr(nc,nr)

openr, lun0 ,file0     , /get_lun,/F77_UNFORMATTED
openr, lun7 ,file7     , /get_lun,/F77_UNFORMATTED
openr, lun8 ,file8     , /get_lun,/F77_UNFORMATTED

 ;OPTION 1
openr, lun1 ,path+file1, /get_lun,/F77_UNFORMATTED
openr, lun2 ,path+file2, /get_lun,/F77_UNFORMATTED
openr, lun3 ,path+file3, /get_lun,/F77_UNFORMATTED
openr, lun4 ,path+file4, /get_lun,/F77_UNFORMATTED
openr, lun5 ,path+file5, /get_lun,/F77_UNFORMATTED
openr, lun6 ,path+file6, /get_lun,/F77_UNFORMATTED

for j = 0L, nr-1L do begin 

    readu,lun1  ,strip1 
    readu,lun2  ,strip2 
    readu,lun3  ,strip3
    readu,lun4  ,strip4 
    readu,lun5  ,strip5 
    readu,lun6  ,strip6
    readu,lun7  ,strip7
    readu,lun8  ,strip8
    readu,lun0  ,strip0
    strip3 = strip3 ; HWSDv2 alredy OC;/1.72 
    strip6 = strip6 ; HWSDv2 alredy OC;;/1.72  ; OM=> OC
    mask(*,j) = 0
    tileid(*,J) = strip0
    data1 (*,j) = fix(100.*strip1)
    data2 (*,j) = fix(100.*strip2) 
    data3 (*,j) = fix(100.*strip3)
    data4 (*,j) = fix(100.*strip4) 
    data5 (*,j) = fix(100.*strip5) 
    data6 (*,j) = fix(100.*strip6) 
    ;  HWSDv1.2 ;data7 (*,j) = fix(100.*(0.3*strip7 + 0.7*strip8)) 
    good7 = (strip7 ne uf)
    good8 = (strip8 ne uf)
    
    out = intarr(nc) + uf
    
    idx = where(good7 AND good8, nidx)
    if (nidx gt 0) then out[idx] = long(0.3*float(strip7[idx]) + 0.7*float(strip8[idx]) + 0.5)
    
    idx = where(good7 AND (NOT good8), nidx)
    if (nidx gt 0) then out[idx] = strip7[idx]
    
    idx = where(good8 AND (NOT good7), nidx)
    if (nidx gt 0) then out[idx] = strip8[idx]
    ; Convert out from integer percent -> integer (percent*100) for ScaleFactor=0.01
    idx = where(out ne uf, nidx2)
    if (nidx2 gt 0) then out[idx] = out[idx] * 100L

    
    data7(*,j) = out

endfor
; Read full src_mask (same shape as data1..data7)
src_mask = lonarr(nc,nr)
tmp = lonarr(nc)

openr, lunsrc2, path+'data_sources.msk', /get_lun, /F77_UNFORMATTED
for j=0L, nr-1L do begin
  readu, lunsrc2, tmp
  src_mask(*,j) = tmp
endfor
close, lunsrc2


m2 = src_mask mod 100L
is_statsgo_mask = (m2 eq 22L) OR (m2 eq 82L) OR (m2 eq 28L)

;tileid (where (tileid eq 0)) = uf
;print,max(tileid),min(tileid)
maxcat = 64770L
undef = uf

for j = 0L, nr-1L do begin
  for i = 0L, nc-1L do begin

    top_tex = ( (data1(i,j) gt 0) AND (data2(i,j) gt 0) )   ; clay_top & sand_top
    sub_tex = ( (data4(i,j) gt 0) AND (data5(i,j) gt 0) )   ; clay_sub & sand_sub

    if (top_tex OR sub_tex) then mask(i,j) = 1 $
    else if ((mask(i,j) eq 0) AND (tileid(i,j) eq 0)) then mask(i,j) = -9999

  endfor
endfor


data1 (where (mask le 0)) = undef
data2 (where (mask le 0)) = undef
data3 (where (mask le 0)) = undef
data4 (where (mask le 0)) = undef
data5 (where (mask le 0)) = undef
data6 (where (mask le 0)) = undef
data7 (where (mask le 0)) = undef

idx = where(is_statsgo_mask AND (data1 eq uf), n1)  ; clay_top missing in statsgo footprint
idx = where(is_statsgo_mask AND (data2 eq uf), n2)  ; sand_top missing
idx = where(is_statsgo_mask AND (data3 eq uf), n3)  ; oc_top missing
idx = where(is_statsgo_mask AND (data4 eq uf), n4)  ; clay_sub missing
idx = where(is_statsgo_mask AND (data5 eq uf), n5)  ; sand_sub missing
idx = where(is_statsgo_mask AND (data6 eq uf), n6)  ; oc_sub missing

print, 'PRE-FILL missing inside STAT footprint: clayT=',n1,' sandT=',n2,' ocT=',n3,' clayS=',n4,' sandS=',n5,' ocS=',n6
n1_pre = n1 & n2_pre = n2 & n3_pre = n3
n4_pre = n4 & n5_pre = n5 & n6_pre = n6

close,lun0
close,lun1
close,lun2
close,lun3
close,lun4
close,lun5
close,lun6
close,lun7
close,lun8

print, 'Done reading. Starting horizontal filling'

;stop
;goto,jump2
mask2 =mask
t0 = systime(/seconds)
DEBUG = 0L ;set to 1L when you want verbose prints
n_try     = 0L   ; times we enter fill block
n_fill    = 0L   ; times we actually fill (copy donor)
n_try_statsgo  = 0L
n_fill_statsgo = 0L
ms = 0L & mt = 0L & mval = 0L

; Also count how many STATSGO-footprint pixels remain missing AFTER fill
n_undef_statsgo_pre  = 0L
n_undef_statsgo_post = 0L
n_hitcap  = 0L   ; times we fail to find donor within lmax
n_bigL    = 0L   ; times l got "large" (diagnostic)
n_hole_statsgo = 0L
n_hole_statsgo_tile = 0L
n_hole_statsgo_mask1 = 0L

for j = 0L, nr-1L do begin
  try0  = n_try
  fill0 = n_fill
  hit0  = n_hitcap

  for i = 0L, nc-1L do begin

    ; STATSGO footprint tag (for diagnostics only)
    m2 = src_mask(i,j) mod 100L
    is_statsgo = (m2 eq 22L) OR (m2 eq 82L) OR (m2 eq 28L)

    ; define "hole" once (used for both debug + fill)
    hole_top = ( (data1(i,j) eq uf) OR (data2(i,j) eq uf) OR $
                 ((data1(i,j) eq 0L) AND (data2(i,j) eq 0L)) )

    hole_sub = ( (data4(i,j) eq uf) OR (data5(i,j) eq uf) OR $
                 ((data4(i,j) eq 0L) AND (data5(i,j) eq 0L)) )

    ; Debug counters: STAT footprint holes (global)
    if (is_statsgo AND (hole_top OR hole_sub)) then begin
      n_hole_statsgo = n_hole_statsgo + 1L
      if ((tileid(i,j) ge 1L) AND (tileid(i,j) le maxcat)) then $
        n_hole_statsgo_tile = n_hole_statsgo_tile + 1L
      if (mask(i,j) eq 1) then $
        n_hole_statsgo_mask1 = n_hole_statsgo_mask1 + 1L
    endif

    ; Only fill inside valid tile domain
    if ((tileid(i,j) ge 1L) AND (tileid(i,j) le maxcat)) then begin

       ; old (too strict it blocks all holes)           
       ; if ((mask(i,j) eq 1) AND (hole_top OR hole_sub)) then begin
       ; new (allows filling holes that are mask=0 but inside valid tiles)
       ;If this pixel is inside a valid tile and is not ocean (mask ≥ 0), and it has missing soil data, then try to fill it.
       ; Here mask means : Do I currently have enough soil data at this pixel to consider it valid?
         if ((mask(i,j) ge 0) AND (hole_top OR hole_sub)) then begin

                n_try = n_try + 1L
                if (is_statsgo) then n_try_statsgo = n_try_statsgo + 1L

                if (DEBUG AND (j ge 3995L) AND (j le 4005L) AND (n_try le 20L)) then begin
                  print, 'HOLE j=', j, ' i=', i, ' mask=', mask(i,j), ' data5=', data5(i,j), $
                         ' top0=', ((data1(i,j) eq 0L) AND (data2(i,j) eq 0L)), $
                         ' sub0=', ((data4(i,j) eq 0L) AND (data5(i,j) eq 0L))
                endif
                
                ; 50-pixel donor existence check (reuse for debug + skip)
                imn0 = MAX([i-50L,0L]) & imx0 = MIN([i+50L,nc-1L])
                jmn0 = MAX([j-50L,0L]) & jmx0 = MIN([j+50L,nr-1L])
                
                msk50 = max(mask2(imn0:imx0, jmn0:jmx0))
                dat50 = max(data5(imn0:imx0, jmn0:jmx0))
                
                if (DEBUG AND (j ge 3400L) AND (j le 4200L) AND ((n_try-try0) le 5L)) then begin
                  print, 'HOLE quickcheck: j=', j, ' i=', i, ' max(mask2@50)=', msk50, ' max(data5@50)=', dat50
                endif                
                ; if no land donor within 50 pixels, skip immediately (avoid expensive l search)
                if (msk50 eq 0) then begin
                  n_hitcap = n_hitcap + 1L
                  continue
                endif
                
                l = 3
                lmax = 500L   ; SCIENCE: do not borrow donor beyond ~2.5°
               ; lmax=300 ≈ 300×0.00833° ≈ 2.5° (~275 km lat)
               ; lmax=500 ≈ 4.17° (~460 km lat)
                subset =0.
                subset2=undef
               ; DEBUG only for first few holes at j=4000
               ;dbg = (j eq 4000L) AND (n_try le 5L)
               dbg = ((j ge 3995L) AND (j le 4005L) AND (n_try le 5L))
               if (dbg) then print, 'DEBUG start i=', i, ' mask=', mask(i,j), ' tileid=', tileid(i,j), ' data5=', data5(i,j), ' lmax=', lmax
                
              ;  while ((max(subset) eq 0) or (max(subset2) le 0))  do begin
                while ( ((max(subset) eq 0) OR (max(subset2) le 0)) AND (l le lmax) ) do begin
                    imx=i+l
                    imn=i-l
                    jmn=j-l
                    jmx=j+l
                    imn=MAX([imn,0])
                    jmn=MAX([jmn,0])
                    imx=MIN([imx,nc-1])
                    jmx=MIN([jmx,nr-1])
                    subset = mask2 (imn:imx,jmn:jmx)
                    subset2= data5 (imn:imx,jmn:jmx)
                    if (dbg AND ((l eq 3L) OR (l eq 10L) OR (l eq 50L) OR (l eq 100L) OR (l eq lmax))) then begin
                      ms = max(subset)
                      mt = max(subset2)
                      print, 'DEBUG i=', i, ' l=', l, ' max(subset)=', ms, ' max(subset2)=', mt
                    endif
                
                    mval = max(subset2,max_subscript)
;	            print,l, max(subset) ,max(subset2) 
                    ;l = l + 1
                    if (l lt 50L) then l = l + 1L else l = MIN([l*2L, lmax+1L])

	   
                endwhile
;                if (l gt lmax) then continue   ; SCIENCE: leave missing rather than far donor
                 if (l gt lmax) then begin
                   n_hitcap = n_hitcap + 1L
                   if ((n_hitcap mod 10000L) eq 0) then print, 'hitcap count=', n_hitcap, ' at j=', j
                   continue
                 endif
                
                l = l -1
                if (l gt 50L) then n_bigL = n_bigL + 1L
                d1=imx-imn+1
                d2=jmx-jmn+1

                ii0 = where((subset2 gt 0.) AND (subset2 ne undef), nsub)
                if (nsub le 0) then continue   ; no usable donors in window
                
                sube = float(subset2[ii0])
                med  = median(sube)
                
                ; prefer an exact (rounded) median if present, else closest-to-median
                imed = long(med + 0.5)
                midx = where(sube eq imed, nm)
                
                if (nm gt 0) then begin
                  max_subscript = ii0(midx(0))
                endif else begin
                  del = abs(float(sube) - med)
                  junk = min(del, kmin0)
                  max_subscript = ii0(kmin0)
                endelse                

                ii = max_subscript - d1*(max_subscript/d1) + imn ;i -l
                jj = max_subscript/d1 + jmn;j - l

                mask(i,j)  = 1
                mask2(i,j) = 1
                data1(i,j) = data1(ii,jj)
                data2(i,j) = data2(ii,jj)
                data3(i,j) = data3(ii,jj)
                data4(i,j) = data4(ii,jj)
                data5(i,j) = data5(ii,jj)
                data6(i,j) = data6(ii,jj)
                data7(i,j) = data7(ii,jj)

                n_fill = n_fill + 1L
                if (is_statsgo) then n_fill_statsgo = n_fill_statsgo + 1L
;                stop
            endif
        endif
    endfor
       dtry = (n_try-try0) & dfill = (n_fill-fill0) & dhit=(n_hitcap-hit0)  
       if (DEBUG AND (j ge 3400L) AND (j le 4200L) AND ((j mod 10L) eq 0) AND (dtry gt 0)) then begin
         print, 'j=', j, ' Δtry=', dtry, ' Δfill=', dfill, ' Δhit=', dhit
       endif       
endfor

idx = where(is_statsgo_mask AND (data1 eq uf), n1)
idx = where(is_statsgo_mask AND (data2 eq uf), n2)
idx = where(is_statsgo_mask AND (data3 eq uf), n3)
idx = where(is_statsgo_mask AND (data4 eq uf), n4)
idx = where(is_statsgo_mask AND (data5 eq uf), n5)
idx = where(is_statsgo_mask AND (data6 eq uf), n6)

print, 'FILL DONE elapsed_s=', systime(/seconds)-t0, $
       ' try=', n_try, ' fill=', n_fill, ' hitcap=', n_hitcap, ' bigL=', n_bigL

print, 'FILL efficiency: fill/try=', float(n_fill)/float(n_try), $
       ' hitcap/try=', float(n_hitcap)/float(n_try)

print, 'STAT footprint holes: total=', n_hole_statsgo, ' tile_ok=', n_hole_statsgo_tile
print, 'STAT footprint filled: try=', n_try_statsgo, ' fill=', n_fill_statsgo

print, 'PRE-FILL missing in STAT footprint: clayT=',n1_pre,' sandT=',n2_pre,' ocT=',n3_pre,' clayS=',n4_pre,' sandS=',n5_pre,' ocS=',n6_pre
print, 'POST-FILL missing in STAT footprint: clayT=',n1,' sandT=',n2,' ocT=',n3,' clayS=',n4,' sandS=',n5,' ocS=',n6


mask2=0
;jump2:
;stop  
for jx =1,18 do begin
    j1 = long((jx -1)*ny)
    j2 = long(jx*ny) - 1L
    for ix =1,36 do begin
        i1 = long((ix -1)*nx)
        i2 = long(ix*nx) - 1L       
        subset_tileid = tileid(i1:i2, j1:j2)
        if (max(subset_tileid) ge 1) then begin  
            ; OPTION 3 
            filename = '/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/out_HWSDv2_NGDC_STATSGO_noMASK/' + $
            ;filename = '/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/out_HWSDv2_NGDC_STATSGO_AFSIS_noMASK/' + $
                       'SoilProperties_H'+string(ix,'(i2.2)')+'V'+string(jx,'(i2.2)') + '.nc'

            subset_data1 = data1(i1:i2,j1:j2)
            subset_data2 = data2(i1:i2,j1:j2)
            subset_data3 = data3(i1:i2,j1:j2)
            subset_data4 = data4(i1:i2,j1:j2)
            subset_data5 = data5(i1:i2,j1:j2)
            subset_data6 = data6(i1:i2,j1:j2)
            subset_data7 = data7(i1:i2,j1:j2)

            ; but write the ORIGINAL mask as the NetCDF Mask variable:
            subset_mask_out = mask(i1:i2, j1:j2)

            long_data = indgen(nx)*dxy + (i1+1)*dxy - dxy/2. -180.
            lat_data  = indgen(ny)*dxy + (j1+1)*dxy - dxy/2. -90.

            STEP5_write_netcdf_file,filename = filename, $
                  nc_global = nc, nr_global = nr, nc_local = nx, nr_local = ny, i_offset = i1+1, j_offset = j1+1, $
                  long_data= long_data, lat_data = lat_data, cellsize=dxy*3600.,scale_factor = [1.,sf,sf,sf,sf,sf,sf,sf], $
                  units = ['_',un,un,un,un,un,un,un], undef=uf, varname = ['Mask','Clay0_30','Sand0_30','OC0_30','Clay30_100','Sand30_100','OC30_100','Gravel'], $
                  data1 = subset_mask_out, data2=subset_data1,data3=subset_data2,data4=subset_data3,data5=subset_data4,data6=subset_data5, $
                  data7 = subset_data6,data8=subset_data7

        endif        
    endfor
endfor

end
