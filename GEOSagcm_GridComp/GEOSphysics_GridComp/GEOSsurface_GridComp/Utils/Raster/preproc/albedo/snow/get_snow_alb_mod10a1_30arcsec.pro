pro get_snow_alb_mod10a1_30arcsec, year=year, h_s=h_s, h_e=h_e

; Code to stitch MODIS files, such as: MOD10A1.A2022002.h35v08.006.2022004050029.hdf
; into a static global 30arcsec grid and produce inputs for GEOS make_bcs package, such as:
; snow_alb_MOD10A1.061_30arcsec_H36V05.nc
; The original MODIS input files are on a sinusoidal (~500-m) grid. The final snow albedo
; product is on a 30arcsec regular lat/lon grid (date-line-edge and pole-edge)

; It first reads the MODIS-tile files to generate snow cover PDFs over a period of time.
; These PDFs are made based on grid cells that exceed a snow cover value (e.g. 60%). This is done
; so a decent amount of snow over a grid box ensures a better estimate of albedo. 
; The PDFs are then used to locate a user-defined cutoff for a %ile of snow coverage
; (e.g. top 10 %tile) and filter snow albedo values. 

; Input argument 'year' in the call specifies for what year(s) to run the code (to allow multiple
; serial runs). Only the PDF building can be processed using multiple serial runs, the remainder 
; of the code has to execute on a single CPU. If no year(s) is provided on the input, the default 
; year is used. Input argument can be an integer or an array of integers indicating year(s) to 
; process (e.g. year=[2019,2020,2021])

; Input arguments 'h_s' and 'h_e' are starting and ending horizontal tiles. These are optional
; and if provided define the range of horizontal MODIS tiles to be processed (it saves time if
; certain tiles need re-processing; otherwise should be ignored)

; dependencies:
; grid.pro (performes gridding)
; read_hdf_sd.pro (reads hdf files)
; read_mod10a1_hdf.pro (reads MODIS files)
; get_lat_lon4tils.pro  (perform conversion of original MODIS sinusoidal projection into lat/lon space)
; 
; [sub]directories to be created prior executing this code:
; data/           - output directory (user needs to create)
; MODIS_lat_lon/  - holds previousely created lat/lon info (user needs to create this with get_lat_lon4tils.pro listed in dependencies)
; MOD10A1_data/   - holds MODIS input data (user has to bring it in; follow the dir structure)

; There are five steps (runs) to process. Steps have to be processed in order. The current step must
; complete before the next is initiated. Use the 'goto' commands to control which step is exectued.
; When the current step is completed, uncomment the next-step goto command on lines 86-92

; To execture the code in terminal window type (without quotations):
; 'idl', then '.r get_snow_alb_mod10a1_30arcsec', then 'get_snow_alb_mod10a1_30arcsec'

; Created April 2022 Biljana Orescanin SSIA@NASA

; path to MODIS data (file name example: MOD10A1.A2022001.h24v06.061.2022003035827.hdf)

path_in='./MOD10A1_data/' ; here is only one file to provide an example

print, "You must make sure the 'path_in' is of the expected structure"

; set parameters to loop over
if keyword_set(year) then begin
  year=year[sort(year)] ; set an increasing order
  years=strtrim(year,2)
  print, ' Processing year(s): ', years
endif else begin
  print,'No year to run specified on the input. Processing for default year - 2020!'
  years        =['2020']
endelse

snow_cvr_min =   0. ; minimum snow cover to consider [%]
cvr_cutoff   =   0  ; snow cover cutoff in %s (to consider only tiles gt this value)
top_alb_limit=  99  ; the top %ile for snow albedo (I choose this value) to be considered (e.g. top 10% of the PDF)
top_limit    =  10  ; the top %ile for snow cover (I choose this value)  to be considered (e.g. top 10% of the PDF)
mis_val      =-999. ; has to be in range [0,255] b/c snow_alb_til is byte array
snow_min_str=strtrim(snow_cvr_min,2)
year_start   =years(0)
year_end     =years(-1)
print, 'Start/End Years:',year_start,'/',year_end

if keyword_set(h_s) and keyword_set(h_e) then begin
  h_start=h_s 
  h_end  =h_e 
endif else begin
  h_start=0  
  h_end  =35 
endelse
print, 'h_start:',h_start,' h_end:',h_end

v_start=0  
v_end  =17 

dim_1d=2400l*2400l
dim_2d=[2400l,2400l]

; ********* STEPS ***********
; [un]comment out the following goto's depending on what part of the code is to be executed
;goto, skip_extracting_snow_cvr   ; *** STEP 2
;goto, skip_calculating_snow_cvr  ; *** STEP 3
;goto, skip_reading_snow_alb      ; *** STEP 4
;goto, skip2ncoutput              ; *** STEP 5
; ********* STEPS ***********

; loop over years
for iyear=0,n_elements(years)-1 do begin

  ; loop over horiontal tiles
  for ih=h_start,h_end do begin
    ih_str=strmid('0'+strtrim(ih,2),1,2,/reverse)
  
    ; loop over vertical tiles
    for iv=v_start,v_end do begin
      iv_str=strmid('0'+strtrim(iv,2),1,2,/reverse)
  
      ; get all files for the current  h/v tile
      filename=file_search(path_in+years(iyear)+ $
                           '/n5eil01u.ecs.nsidc.org/DP4/MOST/MOD10A1.061/'+ $
                           years(iyear)+'.'+'*.*/MOD10A1.A'+years(iyear)+   $
                           '*.h'+ih_str+'v'+iv_str+'.*.hdf',count=nfiles)
  
      ; if no files available for this tile, skip to the next one (note: ~ 1/3 of tiles will be empty)
      if nfiles eq 0 then begin
        print, 'No data for tile: h'+ih_str+'v'+iv_str, ' skipping'
        continue
      endif
  
      ; declare arrays to hold all days info
      snow_cvr_pdf=make_array(dim_1d,101,/long,value=0l) ; to hold pdf of snow cover values
  
      ; loop over files (i.e. days) for this tile to stitch
      for ifile=0,nfiles-1 do begin
      
        ; read file
        data=read_mod10a1_hdf(filename(ifile))
      
        ; keep only what you need
                                                   ; See bottom of the code for details
        snow_cvr_til=data.NDSI_Snow_Cover          ; values [0,100] 
        snow_cvr_qc =data.NDSI_Snow_Cover_Basic_QA ; values [0,255] 0=best, 1=good, 2=ok, 3=poor-not used
      
        ; reform all arrays into 1D arrays so I can apply indexes from where command
        snow_cvr_til=reform(snow_cvr_til,dim_2d[0]*dim_2d[1],/overwrite)
        snow_cvr_qc =reform(snow_cvr_qc ,dim_2d[0]*dim_2d[1],/overwrite)
  
        ; If data present then accumulate non-missing values and their counts
        ind_val=where(snow_cvr_til gt snow_cvr_min and snow_cvr_til le 100. and $
                      snow_cvr_qc  eq 0,c_val)
  
        ; If no valida data of albedo and snow cover in this tile, skip. 
        if c_val eq 0 then begin
          continue
        endif
  
        ; Get PDFs of snow fraction and albedo
        for ival=0,c_val-1 do begin
          snow_cvr_pdf(ind_val(ival),fix(snow_cvr_til[ind_val(ival)]))++
        endfor
  
      endfor ; ifile
  
      openw, lun, 'data/data_out/snow_cov_pdfs_titch08_OD10A1.A.'+years(iyear)+    $
                  '.h'+ih_str+'v'+iv_str+'.bin.gz', /get_lun,/compress
      writeu,lun, snow_cvr_pdf 
      free_lun,lun
  
    endfor ; iv
  endfor ; ih
endfor ; iyear

print, 'Done creating PDFs of snow cover and albedo'
stop

skip_extracting_snow_cvr:
print, 'Starting snow fraction top percentile calculations'

  ; -- get snow fractions for the top n-percentile of PDFs at each pixel
  ; This will take a while since it has to be done pixel by pixel

  ; Few notes:
  ; - Needed is a mean albedo corresponding to the top 10% (by count) snow fraction values
  ;   but only among those vaues that correspond to 60% or more of snow over a grid box.
  ; - So, used is a snow cover cutoff value of 60 to get the top 90the percentile of snow cover values

; loop over horiontal tiles
for ih=h_start,h_end do begin
  ih_str=strmid('0'+strtrim(ih,2),1,2,/reverse)

  ; loop over vertical tiles
  for iv=v_start,v_end do begin
    iv_str=strmid('0'+strtrim(iv,2),1,2,/reverse)

    snow_cvr_pdf_tmp=make_array(dim_1d,101,/long,value=0l) ; to hold current pdf of snow cover values
    snow_cvr_pdf_all=make_array(dim_1d,101,/long,value=0l) ; to hold cummulatinve pdfs of snow cover values
    
    ; First read PDFs of all the years for each tile. Then, find the top %ile and write out the cutoffs

    ; get all files for the current h/v tile (all available years)
    filename_cvr=file_search('data/data_out/snow_cov_pdfs_titch08_OD10A1.A.*.h'+ih_str+'v'+iv_str+'.bin.gz',count=nfiles2)

    ; if no files available for this tile, skip to the next one (note: ~ 1/3 of tiles will be empty)
    if nfiles2 eq 0 then begin
      print, 'No data for tile: h'+ih_str+'v'+iv_str, ' skipping'
      continue
    endif
    
    for ifile=0,nfiles2-1 do begin
      openr,   lun, filename_cvr(ifile), /get_lun,/compress
      readu,   lun, snow_cvr_pdf_tmp 
      free_lun,lun

      snow_cvr_pdf_all=snow_cvr_pdf_all+snow_cvr_pdf_tmp
    endfor

    ; find the snow cover cutoff
    sf_lim_tail_mean=make_array(dim_1d,/float,value=mis_val) ; array to store the top 90th snow fraction %ile cutoff for each pixel

    for ipix=0,dim_1d-1 do begin
    
      tot_pix=total(snow_cvr_pdf_all(ipix,*)) ; tot # of valid snow fraction obs at this pixel

      if tot_pix lt 10 then continue ; if lt 10 valid elements leave as missing (no enough data)

      tot_sf_pix=0l
      for ibin=100,30,-1 do begin ; accumulate from the top bin down till you get enough data.
                                  ; Yet, don't go below 30th bin (i.e. snow fraction lt 30%)
                                  ; b/c there is not enough snow to make it a "reliable albedo estimate"
        tot_sf_pix=tot_sf_pix+snow_cvr_pdf_all(ipix,ibin)*ibin
        if total(snow_cvr_pdf_all(ipix,ibin:100))/tot_pix*100. gt top_limit then begin
          sf_lim_tail_mean(ipix)=tot_sf_pix/total(snow_cvr_pdf_all(ipix,ibin:100))
          break ; ibin loop
        endif
      endfor ; ibin

    endfor ; ipix

    ; write out the snow cover cutoff %ile for this MODIS tile
    openw, lun, 'data/data_out/snow_cov_cutoff_titch08_MOD10A1.A.h'+ih_str+'v'+iv_str+'_'+year_start+'_'+year_end+'.bin.gz', /get_lun,/compress
    writeu,lun, sf_lim_tail_mean ;,sf_max ; both are [dim_1d] float arrays
    free_lun,lun

  endfor ; iv
endfor ;ih

print, 'Done snow fraction top percentile calculations'
stop

skip_calculating_snow_cvr:
print, 'Starting reading in the albedo data to get mean value for those grids that have above the limit snow cover'

; loop over years
for iyear=0,n_elements(years)-1 do begin

  ; loop over horiontal tiles
  for ih=h_start,h_end do begin
    ih_str=strmid('0'+strtrim(ih,2),1,2,/reverse)
  
    ; loop over vertical tiles
    for iv=v_start,v_end do begin
      iv_str=strmid('0'+strtrim(iv,2),1,2,/reverse)
  
      ; get all files for the current h/v tile
      filename=file_search(path_in+years(iyear)+'/n5eil01u.ecs.nsidc.org/DP4/MOST/MOD10A1.061/'+years(iyear)+'.'+ $
                           '*.*/MOD10A1.A'+years(iyear)+'*.h'+ih_str+'v'+iv_str+'.*.hdf',count=nfiles3)
  
      ; if no files available for this tile, skip to the next one (note: ~ 1/3 of tiles will be empty)
      if nfiles3 eq 0 then begin
        print, 'No data for tile: h'+ih_str+'v'+iv_str, ' skipping'
        continue
      endif
  
      ; read in the snow cover cutoff %ile for this MODIS tile
      sf_lim_tail_mean=make_array(dim_1d,/float,value=mis_val) ; to store the top snow fraction %ile cutoff for each pixel
      openr, lun, 'data/data_out/snow_cov_cutoff_titch08_MOD10A1.A.h'+ih_str+'v'+iv_str+'_'+year_start+'_'+year_end+'.bin.gz', /get_lun,/compress
      readu,lun, sf_lim_tail_mean 
      free_lun,lun
  
      ; declare arrays to hold all days info
      accu_alb=make_array(dim_1d,/float,value=0l) ; to hold accumulated snow albedo values
      accu_cnt=make_array(dim_1d,/float,value=0l) ; to hold couhnts of accumulated values
  
      ; loop over files (i.e. days) for this tile to stitch
      for ifile=0,nfiles3-1 do begin
      
        ; read file
        data=read_mod10a1_hdf(filename(ifile))
      
        ; keep only what you need                  ; see bottom of the code for details
        snow_alb_til=data.SNOW_ALBEDO_DAILY_TILE   ; values [0,100] 
        snow_cvr_til=data.NDSI_Snow_Cover          ; values [0,100] 
        snow_cvr_qc =data.NDSI_Snow_Cover_Basic_QA ; values [0,255] 0=best, 1=good, 2=ok, 3=poor-not used
      
        ; reform all arrays into 1D arrays so I can apply indexes from where command
        snow_alb_til=reform(snow_alb_til,dim_2d[0]*dim_2d[1],/overwrite)
        snow_cvr_til=reform(snow_cvr_til,dim_2d[0]*dim_2d[1],/overwrite)
        snow_cvr_qc =reform(snow_cvr_qc ,dim_2d[0]*dim_2d[1],/overwrite)
  
        ; if no valid data of abledo and snow cover in this tile, skip. 
        ; If data present then accumulate non-missing values and their counts
        ind_val=where(snow_alb_til     gt 0 and snow_alb_til le 100 and $
                      snow_cvr_til     gt 0 and snow_cvr_til le 100 and $
                      sf_lim_tail_mean gt 0 and snow_cvr_qc  eq 0,c_val)

        if c_val eq 0 then begin
          continue
        endif
    
        ; get accumulated albedo
        for ipix=0,c_val-1 do begin
          if snow_cvr_til(ind_val[ipix]) ge sf_lim_tail_mean(ind_val[ipix]) then begin
             accu_alb(ind_val[ipix])=accu_alb(ind_val[ipix])+snow_alb_til(ind_val[ipix])
             accu_cnt(ind_val[ipix])=accu_cnt(ind_val[ipix])+1.
          endif
        endfor ; ipix
          
      endfor ; ifile

      ; write out the counts, cumulative and mean values of albedo and PDFs of albedo
      ; for this tile (these will be stitched once all completed)
      openw, lun, 'data/data_out/snow_alb_pdfs_08_MOD10A1.A.'+years(iyear)+           $
                  '.h'+ih_str+'v'+iv_str+'_'+strtrim(100-top_limit,2)+'%ile_cover2_gt_'+ $
                  strtrim(cvr_cutoff,2)+'.bin.gz', /get_lun,/compress
      writeu,lun, accu_alb,accu_cnt
      free_lun,lun
  
    endfor ; iv
  endfor ;ih
endfor ; iyear
print, 'Done with calculating mean snow albedo'
stop
    
skip_reading_snow_alb:

; -- All snow albedo PDFs are in, only in yearly files. Read them all 
;    and make the cumulative stats for mean albedo at a chosen %ile
;    This is to be done on a single CPU as a single sbatch job
print, 'Reading in all Snow Albedo accum and counts to form the overall mean albedo values'

; need to read in all the years for each tile, find top %ile and write out the mean albedo
; loop over horizontal tiles
for ih=h_start,h_end do begin
  ih_str=strmid('0'+strtrim(ih,2),1,2,/reverse)

  ; loop over vertical tiles
  for iv=v_start,v_end do begin
    iv_str=strmid('0'+strtrim(iv,2),1,2,/reverse)

    ; get all files for the current h/v tile (all available years)
    filename_alb=file_search('data/data_out/snow_alb_pdfs_08_MOD10A1.A.*.h'+ $
                             ih_str+'v'+iv_str+'_'+strtrim(100-top_limit,2)   + $
                             '%ile_cover2_gt_'+strtrim(cvr_cutoff,2)+'.bin.gz',count=nfiles4)

    ; if no files available for this tile, skip to the next one (note: ~ 1/3 of tiles will be empty)
    if nfiles4 eq 0 then begin
      print, 'No data for tile: h'+ih_str+'v'+iv_str, ' skipping'
      continue
    endif

    print, 'tile: h'+ih_str+'v'+iv_str
    
    ; declare arrays to hold all days info
    accu_alb_tmp=make_array(dim_1d,/float,value=0l) ; to hold accumulated snow albedo values
    accu_cnt_tmp=make_array(dim_1d,/float,value=0l) ; to hold counts of accumulated values
    snow_alb_accu_all=0.
    snow_alb_cnt_all =0.

    for ifile=0,nfiles4-1 do begin
      openr,   lun, filename_alb(ifile), /get_lun,/compress
      readu,   lun, accu_alb_tmp,accu_cnt_tmp
      free_lun,lun
      snow_alb_accu_all=snow_alb_accu_all+accu_alb_tmp 
      snow_alb_cnt_all =snow_alb_cnt_all +accu_cnt_tmp 
    endfor

    ; get mean albedo
    mean_alb=snow_alb_accu_all/snow_alb_cnt_all
    mean_alb[where(snow_alb_cnt_all lt 10,/null)]=mis_val

    ; write out the snow albedo max and cutoff %ile for this MODIS tile
    openw, lun, 'data/data_out/snow_alb_08_'+strtrim(top_alb_limit,2)+'_cutoff_MOD10A1.A.h'+ $
                ih_str+'v'+iv_str+'_'+year_start+'_'+year_end+'.bin.gz', /get_lun,/compress
    writeu,lun, mean_alb
    free_lun,lun

    ; reset arrays for accumulating snow albedo
    snow_alb_accu_all(*)=0l
    snow_alb_cnt_all (*)=0l

  endfor ; iv
endfor ;ih

print, 'Done with snow albedo cutoff %ile calculations'
stop

skip2ncoutput:
; -- Now all the albedo vaules are in. Output them on 
;    GEOS-friendly 10x10 deg tiles using nc format

; set new missing value
mis_val= 1.e15

print, 'Starting stitching'

; array to store the top albedo %ile cutoff for each pixel at a given MODIS tile
alb_lim_tail_mean=make_array(dim_1d,/float,value=mis_val) 

; create 36x18 tiles at 30sec arc resolution (each tile is 10x10deg with 1200x1200 grid boxes)
; loop over horiontal tiles
for ih=h_start,h_end do begin
  ih_str=strmid('0'+strtrim(ih+1,2),1,2,/reverse)

  ; loop over vertical tiles
  for iv=v_start,v_end do begin
    iv_str=strmid('0'+strtrim(iv+1,2),1,2,/reverse)

    ; declare array to store the output for this tile; fill with missing 
    alb_30sec_grid=make_array(1200l, 1200l, value=mis_val,/float)

    print, 'Creating tile: h'+ih_str+'v'+iv_str

    ; get min/max and all lat/lon values for this tile
    minlat       = iv   *10.-90.   
    maxlat       =(iv+1)*10.-90.   
    minlon       = ih   *10.-180.0
    maxlon       =(ih+1)*10.-180.0
    lat_positions=minlat+indgen(1200)*10./1200. 
    lon_positions=minlon+indgen(1200)*10./1200. 

    ; have to read +/- 1 MODIS tile in vertical direction
    ; and as many as needed to cover all valid lat/lons in horizontal direction
    
    alb_lim_tail_mean_all=[] ; To store all 500m resolution albedo values
    lat_all              =[] ; corresponding to this tile to be used in 
    lon_all              =[] ; gridding

    ; loop over vertical tiles
    for iiv=17-iv-1,17-iv+1 do begin
      iiv_str=strmid('0'+strtrim(iiv,2),1,2,/reverse)

      ; loop over horizontal tiles
      for iih=0,35 do begin
        iih_str=strmid('0'+strtrim(iih,2),1,2,/reverse)

        ; read in the cumulative and mean values for this tile (these will be stitched once all completed)
        filename_alb=file_search('data/data_out/snow_alb_08_'+strtrim(top_alb_limit,2)+  $
                                 '_cutoff_MOD10A1.A.h'+iih_str+'v'+iiv_str+'_'+year_start+  $
                                 '_'+year_end+'.bin.gz',count=n_files5)

        if n_files5 ne 1 then continue

        openr, lun, filename_alb, /get_lun,/compress
        readu,lun, alb_lim_tail_mean 
        free_lun,lun

        ; Reform the arrays back to 2 dimensions
        alb_lim_tail_mean=reform(alb_lim_tail_mean,dim_2d,/overwrite)

        ; read lat and lon for this tile
        lat2d        = fltarr(2400l,2400l)
        lon2d        = fltarr(2400l,2400l)
        file_lat_lon = 'MODIS_lat_lon/MODIS_hdres_lat_lon_v'+iiv_str+'_h'+iih_str+'.bin.gz'
        openr   , lun, file_lat_lon, /get_lun,/compress 
        readu   , lun, lat2d,lon2d
        free_lun, lun

        ; if no valid values for snow albedo on this tile, write all missing vals
        ind_fit=where(lon2d ge minlon and lon2d le maxlon and $
                      lat2d ge minlat and lat2d le maxlat and $
                      alb_lim_tail_mean ge 0 and alb_lim_tail_mean le 100 , c_fit)

        ; if no grids with corresponding lat/lons to this tile, skip it
        if c_fit eq 0 then continue 

        ; if valid grids corresponding to this tile exist then remember them
        alb_lim_tail_mean_all=[alb_lim_tail_mean_all,alb_lim_tail_mean[ind_fit]]
        lon_all              =[lon_all              ,lon2d            [ind_fit]]
        lat_all              =[lat_all              ,lat2d            [ind_fit]]

      endfor ; iih
    endfor ; iiv

    ; if there are any valid values, grid them to populate the tile
    if n_elements(alb_lim_tail_mean_all) gt 0 then begin

      print, 'Start gridding'

      ; regrid from 2400x2400 sinusoidal grid to 1200x1200 (30 arc sec) regular lat/lon grid (dateline-on-edge, pole-on-edge)
      alb_30sec_grid=grid(alb_lim_tail_mean_all,lat_all,lon_all,nlat=1200,nlon=1200,  $
                          region=[minlat,maxlat,minlon,maxlon],mis_val=mis_val)

      ; scale to range [0.0,1.0]
      ind_val=where(alb_30sec_grid ne mis_val,c_val)
      if c_val gt 0 then alb_30sec_grid[ind_val]=alb_30sec_grid[ind_val]/100.

    endif

    ; write snow albedo in NetCDF format for the current MODIS tile
    ;   *** Create a NCDF file 

    ; Set up the file & handler
    nc_file='data/data_out/snow_alb_all_08_Top'+strtrim(top_alb_limit,2)+ $
            'th_percentile_MOD10A1.A_30arcsec_'+year_start+'_'+year_end+     $
            '_H'+ih_str+'V'+iv_str+'.nc' 

    ; create the file
    fid=ncdf_create(nc_file,/CLOBBER,/NETCDF4_FORMAT) ; erese existing file and make a new one

    ; write global attributes 
    NCDF_ATTPUT, fid, 'N_lon_global'   , '43200'              , /GLOBAL,/char ; Total number of grids in i-direction
    NCDF_ATTPUT, fid, 'N_lat_global'   , '21600'              , /GLOBAL,/char ; Total number of grids in j-direction
    NCDF_ATTPUT, fid, 'i_ind_offset_LL', strtrim(ih*1200l+1,2), /GLOBAL,/char ; GEOS-friendly H grid box 
    NCDF_ATTPUT, fid, 'j_ind_offset_LL', strtrim(iv*1200l+1,2), /GLOBAL,/char ; GEOS-friendly V grid box
    NCDF_ATTPUT, fid, 'CellSize='      , '30arcsec'           , /GLOBAL,/char ; grid size
    NCDF_ATTPUT, fid, 'CreatedBy'      , 'borescan'           , /GLOBAL,/char ; user info
    spawn,'date',date1 ; get the time form the system
    NCDF_ATTPUT, fid, 'Date'           , strtrim(date1,2)     , /GLOBAL,/char ; time and date
    NCDF_ATTPUT, fid, 'Data_Resolution',' 1200 x 1200 '       , /GLOBAL,/char ; resolution of the run
    NCDF_ATTPUT, fid, 'Region:'        ,'MODIS tile H'+ih_str+'V'+iv_str, /GLOBAL,/char ; covearge

    ;Set Dimensions (there will be two dimentsions: lat and lon, i.e. x and y)
    d=indgen(2,/LONG) ; number of dimensions that will be created
    d[0]=ncdf_dimdef(fid,'N_lon',1200l)   ; number of longitudes   
    d[1]=ncdf_dimdef(fid,'N_lat',1200l)   ; number of latitudes  

    ; Define variables to be stored
    var_id1=ncdf_vardef(fid,'lon'        ,d[0],/FLOAT)
    var_id2=ncdf_vardef(fid,'lat'        ,d[1],/FLOAT)
    var_id3=ncdf_vardef(fid,'Snow_Albedo',d   ,/FLOAT)

    ; Change modes (to enter data section/group)
    ncdf_control,fid,/ENDEF

    ; Write the data
    ncdf_varput,fid,'lat', lat_positions
    ncdf_varput,fid,'lon', lon_positions

    ; missing value for MAPL is 1.e+15f
    ncdf_varput,fid,'Snow_Albedo', alb_30sec_grid 
    NCDF_ATTPUT,fid,'Snow_Albedo', "scale_factor", 1.0
    NCDF_ATTPUT,fid,'Snow_Albedo', "offset", 0.0
    NCDF_ATTPUT,fid,'Snow_Albedo', "missing_value",  1.e15 
    NCDF_ATTPUT,fid,'Snow_Albedo', "long_name", string("Snow Albedo"),/char
    NCDF_ATTPUT,fid,'Snow_Albedo', "units", "-",/char

    ; Close the file & release the handler   
    ncdf_close,fid

    print, 'Done creating NCDF file'
    
  endfor ; iv
endfor ;ih

; values in the SNOW_ALBEDO_DAILY_TILE are as following: 
 ; 1-100: snow albedo 
 ; 101: no decision
 ; 111: night
 ; 125: land
 ; 137: inland water
 ; 139: ocean
 ; 150: cloud
 ; 151: cloud detected as snow 250: missing
 ; 251: self-shadowing
 ; 252: land mask mismatch 253: BRDF failure
 ; 254: non-production mask

; NDSI snow cover general quality value   ; values in the SNOW_ALBEDO_DAILY_TILE are as following: 
 ; 0=best,                                ; 1–100: snow albedo 
 ; 1=good,                                ; 101: no decision
 ; 2=ok,                                  ; 111: night
 ; 3=poor-not used,                       ; 125: land
 ; 4=other-not used,                      ; 137: inland water
 ; 211=night,                             ; 139: ocean
 ; 239=ocean,                             ; 150: cloud
 ; 255=unusable L1B data or no data       ; 151: cloud detected as snow 250: missing
                                          ; 251: self-shadowing
                                          ; 252: land mask mismatch 253: BRDF failure
                                          ; 254: non-production mask

stop 
end
