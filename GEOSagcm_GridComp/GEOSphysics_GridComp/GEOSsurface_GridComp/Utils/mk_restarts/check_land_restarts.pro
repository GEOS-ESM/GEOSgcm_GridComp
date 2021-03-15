; =========================================================================
; USAGE :
; Edit lines to 44-47 specify paths to BCs dir and yjecatch{cn}_internal_rst file
; =========================================================================
;_____________________________________________________________________
;_____________________________________________________________________

FUNCTION NCDF_ISNCDF, FILENAME

;- Set return values

false = 0B
true = 1B

;- Establish error handler

catch, error_status
if error_status ne 0 then begin
  catch, /cancel
  return, false
endif

;- Try opening the file

cdfid = ncdf_open( filename )

;- If we get this far, open must have worked

ncdf_close, cdfid
catch, /cancel
return, true

END

;_____________________________________________________________________
;_____________________________________________________________________

pro plot_rst

; **********************************************************************************************************
; STEP (1) Specify below:
; -----------------------

BCSDIR   = '/discover/nobackup/smahanam/bcs/Heracles-4_3/Heracles-4_3_MERRA-3/CF0180x6C_DE1440xPE0720/'
GFILE    = 'CF0180x6C_DE1440xPE0720-Pfafstetter'
OutDir   = 'OutData2/'
int_rst  = 'catchcn_internal_rst'

; STEP (2) save :
; ---------------
; On dali : (a) module load tool/idl-8.5, (b) idl (c) .compile chk_restarts
;               and (d) plot_rst

; **********************************************************************************************************

; Setting up and select variables for plotting
; --------------------------------------------

TILFILE  = BCSDIR + 'til/' + GFILE + '.til'
RSTFILE  = BCSDIR + 'rst/' + GFILE + '.rst'

NTILES   = 0l
NG       = 0l
NC       = 0l
NR       = 0l

openr,1,BCSDIR + 'clsm/catchment.def'
readf,1,NTILES
close,1

openr,1,TILFILE
readf,1,NG,NC,NR
close,1

Var_Names = [ $
            'CDCR2'   , $       ; 0
            'BEE'     , $       ; 1    
            'POROS'   , $       ; 2
            'ITY1'    , $       ; 3
            'ITY2'    , $       ; 4
            'ITY3'    , $       ; 5
            'ITY4'    , $       ; 6
            'TC1'     , $       ; 7
            'TC2'     , $       ; 8
            'TC3'     , $       ; 9
            'TC4'     , $       ;10
            'CATDEF'  , $       ;11
            'RZEXC'   , $       ;12
            'SFEXC'   ]

N_VARS    = N_ELEMENTS (Var_Names)
PLOT_VARS = fltarr (NTILES,N_VARS) 
TMP_VAR1  = fltarr (NTILES)
TMP_VAR2  = fltarr (NTILES,4)


; Get file information : (1) model, (2) file format
; -------------------------------------------------

catch_model = boolean (strcmp(int_rst,'catchcn',7,/fold_case) eq 0)
ncdf_file   = boolean (ncdf_isncdf(OutDir + int_rst))

; Set up vector to grid for plotting
; ----------------------------------

NC_plot = 4320
NR_plot = 2160

tileid_plot = lonarr (NC_plot,NR_plot)

dx = NC/NC_plot
dy = NR/NR_plot

catrow = lonarr(nc)
cat    = lonarr(nc,dy)

openr,1,RSTFILE,/F77_UNFORMATTED

for j = 0l, NR_plot -1 do begin
   
   for i=0,dy -1 do begin
      readu,1,catrow
      cat (*,i) = catrow
   endfor
   
   for i = 0, NC_plot -1 do begin
      subset = cat (i*dx: (i+1)*dx -1,*)
      if (min (subset) le NTILES) then begin
         min1 = min(subset)
         subset(where (subset gt NTILES)) = 0
         hh = histogram(subset,bin=1,min = min1, locations=loc_val)
         dom_tile = max(hh,loc)	
         tileid_plot[i,j] = loc_val(loc)
      endif
   endfor

endfor

close,1

; Reading catch*_internal_rst
; ---------------------------

if (ncdf_file) then begin

   ncid = NCDF_OPEN(OutDir  + int_rst,/NOWRITE)
   result = ncdf_inquire( ncid)
   if(result.nvars gt 60) then catch_model = boolean (result.nvars lt 60)
   NCDF_VARGET, ncid,'CDCR2'  ,TMP_VAR1
   PLOT_VARS (*,0) = TMP_VAR1
   NCDF_VARGET, ncid,'BEE'    ,TMP_VAR1
   PLOT_VARS (*,1) = TMP_VAR1
   NCDF_VARGET, ncid,'POROS'    ,TMP_VAR1
   PLOT_VARS (*,2) = TMP_VAR1
   NCDF_VARGET, ncid,'TC'    ,TMP_VAR2 
   PLOT_VARS (*,7) = TMP_VAR2(*,0)
   PLOT_VARS (*,8) = TMP_VAR2(*,1)
   PLOT_VARS (*,9) = TMP_VAR2(*,2)
   PLOT_VARS (*,10)= TMP_VAR2(*,3)
   NCDF_VARGET, ncid,'CATDEF'    ,TMP_VAR1
   PLOT_VARS (*,11) = TMP_VAR1
   NCDF_VARGET, ncid,'RZEXC'    ,TMP_VAR1
   PLOT_VARS (*,12) = TMP_VAR1
   NCDF_VARGET, ncid,'SRFEXC'    ,TMP_VAR1
   PLOT_VARS (*,13) = TMP_VAR1

   if(catch_model) then begin

      NCDF_VARGET, ncid,'OLD_ITY'    ,TMP_VAR1
      PLOT_VARS (*,3) = TMP_VAR1
 
   endif else begin

      NCDF_VARGET, ncid,'ITY'    ,TMP_VAR2      
      PLOT_VARS (*,3) = TMP_VAR2(*,0)
      PLOT_VARS (*,4) = TMP_VAR2(*,1)
      PLOT_VARS (*,5) = TMP_VAR2(*,2)
      PLOT_VARS (*,6) = TMP_VAR2(*,3)

   endelse
   
   NCDF_CLOSE, ncid

endif else begin

   openr,1,OutDir  + int_rst, /F77_UNFORMATTED

   if(catch_model) then begin

      for i = 1,30 do begin
         readu,1,TMP_VAR1
         if (i eq  6) then PLOT_VARS (*,0) = TMP_VAR1
         if (i eq  8) then PLOT_VARS (*,1) = TMP_VAR1
         if (i eq  9) then PLOT_VARS (*,2) = TMP_VAR1
         if (i eq 30) then PLOT_VARS (*,3) = TMP_VAR1
      endfor

      readu,1,TMP_VAR2
      PLOT_VARS (*,7) = TMP_VAR2(*,0)
      PLOT_VARS (*,8) = TMP_VAR2(*,1)
      PLOT_VARS (*,9) = TMP_VAR2(*,2)
      PLOT_VARS (*,10)= TMP_VAR2(*,3)

      readu,1,TMP_VAR2
      readu,1,TMP_VAR1
      readu,1,TMP_VAR1
      PLOT_VARS (*,11) = TMP_VAR1
      readu,1,TMP_VAR1
      PLOT_VARS (*,12) = TMP_VAR1
      readu,1,TMP_VAR1
      PLOT_VARS (*,13) = TMP_VAR1

   endif else begin

      for i = 1,37 do begin
         readu,1,TMP_VAR1
         if (i eq  6) then PLOT_VARS (*,0) = TMP_VAR1
         if (i eq  8) then PLOT_VARS (*,1) = TMP_VAR1
         if (i eq  9) then PLOT_VARS (*,2) = TMP_VAR1
         if (i eq 30) then PLOT_VARS (*,3) = TMP_VAR1
         if (i eq 31) then PLOT_VARS (*,4) = TMP_VAR1
         if (i eq 32) then PLOT_VARS (*,5) = TMP_VAR1
         if (i eq 33) then PLOT_VARS (*,6) = TMP_VAR1
      endfor 
      readu,1,TMP_VAR2
      PLOT_VARS (*,7) = TMP_VAR2(*,0)
      PLOT_VARS (*,8) = TMP_VAR2(*,1)
      PLOT_VARS (*,9) = TMP_VAR2(*,2)
      PLOT_VARS (*,10)= TMP_VAR2(*,3)

      readu,1,TMP_VAR2
      readu,1,TMP_VAR2
      readu,1,TMP_VAR1
      readu,1,TMP_VAR1
      PLOT_VARS (*,11) = TMP_VAR1
      readu,1,TMP_VAR1
      PLOT_VARS (*,12) = TMP_VAR1
      readu,1,TMP_VAR1
      PLOT_VARS (*,13) = TMP_VAR1

   endelse

   close,1

endelse

; Plotting
; --------

spawn, 'mkdir -p ' + OutDir + 'plots'
load_colors

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[720,800], Z_Buffer=0
Erase,255
!p.background = 255

!P.position=0
!P.Multi = [0, 2, 3, 0, 0]

plot_6maps, ntiles, tileid_plot, PLOT_VARS(*,0), [min(PLOT_VARS(*,0)), max(PLOT_VARS(*,0))]  , Var_Names  (0)
plot_6maps, ntiles, tileid_plot, PLOT_VARS(*,1), [min(PLOT_VARS(*,1)), max(PLOT_VARS(*,1))]  , Var_Names  (1),advance =1 
plot_6maps, ntiles, tileid_plot, PLOT_VARS(*,2), [0.37,0.8]                                  , Var_Names  (2),advance =1 
plot_6maps, ntiles, tileid_plot, PLOT_VARS(*,11),[min(PLOT_VARS(*,11)), max(PLOT_VARS(*,11))], Var_Names (11),advance =1 
plot_6maps, ntiles, tileid_plot, PLOT_VARS(*,12),[min(PLOT_VARS(*,12)), max(PLOT_VARS(*,12))], Var_Names (12),advance =1 
plot_6maps, ntiles, tileid_plot, PLOT_VARS(*,13),[min(PLOT_VARS(*,13)), max(PLOT_VARS(*,13))], Var_Names (13),advance =1 

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 720, 800)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, OutDir + 'plots/soil_var.jpg', image24, True=1, Quality=100

plot_tc, NTILES, tileid_plot,OutDir + 'plots/', plot_vars (*,7), plot_vars (*,8), plot_vars (*,9), plot_vars (*,10)

if(catch_model) then begin
   plot_mosaic, ntiles, OutDir + 'plots/', tileid_plot, fix(plot_vars (*,3))
endif else begin
   plot_carbon, ntiles, OutDir + 'plots/', tileid_plot, fix(plot_vars (*,3:6))
endelse

end

;_____________________________________________________________________
;_____________________________________________________________________

pro check_regrid_carbon

; **********************************************************************************************************
; STEP (1) Specify below:
; -----------------------

BCSDIR1   = '/discover/nobackup/smahanam/bcs/Heracles-4_3/Heracles-4_3_MERRA-3/SMAP_EASEv2_M09/'
GFILE1    = 'SMAP_EASEv2_M09_3856x1624'
OutDir1   = '/discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/CatchCN/M09/20151231/'
int_rst1  = 'catchcn_internal_rst'

BCSDIR2   = '/discover/nobackup/smahanam/bcs/Heracles-4_3/Heracles-4_3_MERRA-3/CF0180x6C_DE1440xPE0720/'
GFILE2    = 'CF0180x6C_DE1440xPE0720-Pfafstetter'
OutDir2   = ''
int_rst2  = 'catchcn_internal_rst'

; STEP (2) save :
; ---------------
; On dali : (a) module load tool/idl-8.5, (b) idl (c) .compile chk_restarts
;               and (d) plot_rst

; **********************************************************************************************************

; Setting up and select variables for plotting
; --------------------------------------------
Var_Names = [ $
            'CDCR2'   , $       ; 0
            'BEE'     , $       ; 1    
            'POROS'   , $       ; 2
            'ITY1'    , $       ; 3
            'ITY2'    , $       ; 4
            'ITY3'    , $       ; 5
            'ITY4'    , $       ; 6
            'TC1'     , $       ; 7
            'TC2'     , $       ; 8
            'TC3'     , $       ; 9
            'TC4'     , $       ;10
            'CATDEF'  , $       ;11
            'RZEXC'   , $       ;12
            'SFEXC'   ]
NC_plot = 4320
NR_plot = 2160

;goto, jump

for resol = 1,2 do begin

if(resol eq 1) then begin
   BCSDIR   = BCSDIR1	
   TILFILE  = BCSDIR1 + 'til/' + GFILE1 + '.til'
   RSTFILE  = BCSDIR1 + 'rst/' + GFILE1 + '.rst'
endif else begin
   BCSDIR   = BCSDIR2	
   TILFILE  = BCSDIR2 + 'til/' + GFILE2 + '.til'
   RSTFILE  = BCSDIR2 + 'rst/' + GFILE2 + '.rst'
endelse	


NTILES   = 0l
NG       = 0l
NC       = 0l
NR       = 0l

openr,1,BCSDIR + 'clsm/catchment.def'
readf,1,NTILES
close,1

openr,1,TILFILE
readf,1,NG,NC,NR
close,1
; Set up vector to grid for plotting
; ----------------------------------



tileid_plot = lonarr (NC_plot,NR_plot)

dx = NC/NC_plot
dy = NR/NR_plot

catrow = lonarr(nc)
cat    = lonarr(nc,dy)

openr,1,RSTFILE,/F77_UNFORMATTED

for j = 0l, NR_plot -1 do begin
   
   for i=0,dy -1 do begin
      readu,1,catrow
      cat (*,i) = catrow
   endfor
   
   for i = 0, NC_plot -1 do begin
      subset = cat (i*dx: (i+1)*dx -1,*)
      if (min (subset) le NTILES) then begin
         min1 = min(subset)
         subset(where (subset gt NTILES)) = 0
         hh = histogram(subset,bin=1,min = min1, locations=loc_val)
         dom_tile = max(hh,loc)	
         tileid_plot[i,j] = loc_val(loc)
      endif
   endfor

endfor

close,1
if (resol eq 1) then begin
   tileid_plot1 = tileid_plot
   NTILES1 = NTILES
endif else begin
   tileid_plot2 = tileid_plot
   NTILES2 = NTILES
endelse
endfor

cnpft1 = fltarr (ntiles1, 888)
cnpft2 = fltarr (ntiles2, 888)
fvg1   = fltarr (ntiles1, 4)
fvg2   = fltarr (ntiles2, 4)
ncid = NCDF_OPEN(OutDir1  + int_rst1,/NOWRITE)
NCDF_VARGET, ncid,'TILE_ID'  ,TILE_ID
NCDF_VARGET, ncid,'CNPFT'  ,CNPFT1
NCDF_VARGET, ncid,'FVG'  ,fvg1
TILE_ID = long (TILE_ID) - 1l

CNPFT=CNPFT1
FVG  =FVG1

for k =0l,n_elements (CNPFT1(*,0)) -1l do CNPFT1(TILE_ID(k),*) = CNPFT(k,*)
for k =0l,n_elements (FVG1  (*,0)) -1l do FVG1  (TILE_ID(k),*) = FVG  (k,*)

CNPFT=0.
FVG  =0.
NCDF_CLOSE, ncid

ncid = NCDF_OPEN(OutDir2  + int_rst2,/NOWRITE)
NCDF_VARGET, ncid,'CNPFT'  ,CNPFT2
NCDF_VARGET, ncid,'FVG'  ,fvg2
NCDF_CLOSE, ncid
save,NTILES1,NTILES2,tileid_plot1,tileid_plot2,CNPFT1,CNPFT2, fvg1, fvg2,file = 'temp_file.idl'
;stop

jump:

restore,'temp_file.idl'

; Plotting
; --------

spawn, 'mkdir -p plots'
load_colors
limits = [-60,-180,90,180]

plot_varid = 14
cnpft1 = reform ( cnpft1,[ntiles1,3,4,74],/overwrite)
cnpft2 = reform ( cnpft2,[ntiles2,3,4,74],/overwrite)

for iv = 1,4 do begin

plot_vars1 = cnpft1(*,0,iv - 1,plot_varid-1)
plot_vars2 = cnpft2(*,0,iv - 1,plot_varid-1)


thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[720,1000], Z_Buffer=0
Erase,255
!p.background = 255

!P.position=0
!P.Multi = [0, 1, 2, 0, 1]

plot_2maps, ntiles1, tileid_plot1, plot_vars1(*), [min([PLOT_VARS1,plot_vars2],/nan),max([PLOT_VARS1,plot_vars2],/nan)], string(plot_varid,'(i2.2)')+'_v' + string(iv,'(i1.1)')
plot_2maps, ntiles2, tileid_plot2, plot_vars2(*), [min([PLOT_VARS1,plot_vars2],/nan),max([PLOT_VARS1,plot_vars2],/nan)], string(plot_varid,'(i2.2)')+'_v' + string(iv,'(i1.1)'),advance =1 

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 720, 1000)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'plots/pft_'+ string(plot_varid,'(i2.2)')+'_v' + string(iv,'(i1.1)') +'.jpg', image24, True=1, Quality=100
endfor
fvg1(where (fvg1 le 1.e-4)) = !VALUES.F_NAN
fvg2(where (fvg2 le 1.e-4)) = !VALUES.F_NAN

plot_fr, NTILES1, tileid_plot1,'plots/offl_', fvg1 (*,0), fvg1 (*,1), fvg1 (*,2), fvg1 (*,3)

plot_fr, NTILES2, tileid_plot2,'plots/agcm_', fvg2 (*,0), fvg2 (*,1), fvg2 (*,2), fvg2 (*,3)

end
;_____________________________________________________________________
;_____________________________________________________________________

PRO plot_2maps, ncat, tile_id, data, vlim, vname,advance = advance

lwval = vlim(0)
upval = vlim(1)

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

data_grid = fltarr (im,jm)
data_grid (*,*) =  !VALUES.F_NAN

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(tile_id[i,j] gt 0) then data_grid(i,j) =  data(tile_id[i,j] -1)
   endfor
endfor

limits = [-60,-180,90,180]

colors = [27,26,25,24,23,22,21,20,40,41,42,43,44,45,46,47,48]
n_levels = n_elements (colors)

levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]

if keyword_set (advance) then begin

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ADVANCE,/ISOTROPIC,/NOBORDER, title =vname
endif else begin
MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ISOTROPIC,/NOBORDER, title =vname
endelse

contour, data_grid,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot

levels_x = levels

alpha=fltarr(n_levels,2)
alpha(*,0)=levels
alpha(*,1)=levels
h=[0,1]

dx = (240.)/(n_levels-1)

clev = levels
clev (*) = 1
n=0
k = 0
fmt_string = '(f7.2)'
!P.position=[0.30, 0.0+0.005, 0.70, 0.015+0.005]

contour,alpha,levels_x,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels_x,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(colors) -1 do xyouts,levels_x[k],1.1,string(levels[k],format=fmt_string) ,orientation=90,color=0,charsize =0.8

!P.position=0

END

;_____________________________________________________________________
;_____________________________________________________________________

pro plot_fr, ncat, tile_id,out_path, VISDR, VISDF, NIRDR, NIRDF

load_colors
thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[720,500], Z_Buffer=0
Erase,255
!p.background = 255

!P.position=0
!P.Multi = [0, 2, 2, 0, 0]
limits = [-60,-180,90,180]

lwval = 0.
upval = 1.

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

colors = [27,26,25,24,23,22,21,20,40,41,42,43,44,45,46,47,48]
n_levels = n_elements (colors)

for map = 1,4 do begin

   if (map eq 1) then data = VISDR
   if (map eq 2) then data = VISDF	
   if (map eq 3) then data = NIRDR
   if (map eq 4) then data = NIRDF

   if (map eq 1) then ctitle = 'PF1'
   if (map eq 2) then ctitle = 'PF2'	
   if (map eq 3) then ctitle = 'SF1'
   if (map eq 4) then ctitle = 'SF2'

   levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]

   data_grid = fltarr (im,jm)
   data_grid (*,*) =  !VALUES.F_NAN

   for j = 0l, jm -1l do begin
       for i = 0l, im -1 do begin
           if(tile_id[i,j] gt 0) then data_grid(i,j) =  data(tile_id[i,j] -1)
       endfor
   endfor

   if(map eq 1) then begin
      MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ISOTROPIC,/NOBORDER, title = ctitle
   endif else begin
      MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ADVANCE,/ISOTROPIC,/NOBORDER, title = ctitle
   endelse

   contour, data_grid,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
   if(map eq 3) then begin
      !P.position=[0.25, 0.05, 0.75, 0.075]

      alpha=fltarr(n_levels,2)
      alpha(*,0)=levels
      alpha(*,1)=levels
      h=[0,1]		
      clev = levels
      clev (*) = 1
      n=0
      k = 0
      fmt_string = '(f6.2)'
      contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
      contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(colors) -1 do xyouts,levels[k],1.1,string(levels[k],format=fmt_string) ,orientation=90,color=0,charsize =0.8
   !P.position=0
   endif
endfor

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 720, 500)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, out_path +'FR.jpg', image24, True=1, Quality=100

end

;_____________________________________________________________________
;_____________________________________________________________________

pro load_colors

R = intarr (256)
G = intarr (256)
B = intarr (256)

R (*) = 255
G (*) = 255
B (*) = 255

r_drought = [0,   0,   0,   0,  47, 200, 255, 255, 255, 255, 249, 197]
g_drought = [0, 115, 159, 210, 255, 255, 255, 255, 219, 157,   0,   0]
b_drought = [0,   0,   0,   0,  67, 130, 255,   0,   0,   0,   0,   0]

colors = indgen (11) + 1
R (0:11) = r_drought
G (0:11) = g_drought
B (0:11) = b_drought

r_green = [200, 150,  47,  60,   0,   0,   0,   0]
g_green = [255, 255, 255, 230, 219, 187, 159, 131]
b_green = [200, 150,  67,  15,   0,   0,   0,   0]

r_blue  = [ 55,   0,   0,   0,   0,   0,   0,   0,   0,   0]
g_blue  = [255, 255, 227, 195, 167, 115,  83,   0,   0,   0]
b_blue  = [199, 255, 255, 255, 255, 255, 255, 255, 200, 130]

r_red   = [255, 240, 255, 255, 255, 255, 255, 233, 197]
g_red   = [255, 255, 219, 187, 159, 131,  51,  23,   0]
b_red   = [153,  15,   0,   0,   0,   0,   0,   0,   0]

r_grey  = [245, 225, 205, 185, 165, 145, 125, 105,  85]
g_grey  = [245, 225, 205, 185, 165, 145, 125, 105,  85]
b_grey  = [245, 225, 205, 185, 165, 145, 125, 105,  85]

r_type  = [255,106,202,251,  0, 29, 77,109,142,233,255,255,255,127,164,164,217,217,204,104,  0]
g_type  = [245, 91,178,154, 85,115,145,165,185, 23,131,131,191, 39, 53, 53, 72, 72,204,104, 70]
b_type  = [215,154,214,153,  0,  0,  0,  0, 13,  0,  0,200,  0,  4,  3,200,  1,200,204,200,200]

R (20:27) = r_green
G (20:27) = g_green
B (20:27) = b_green

R (30:39) = r_blue
G (30:39) = g_blue
B (30:39) = b_blue

R (40:48) = r_red
G (40:48) = g_red
B (40:48) = b_red

R (50:58) = r_grey
G (50:58) = g_grey
B (50:58) = b_grey

R (60:80) = r_type
G (60:80) = g_type
B (60:80) = b_type

TVLCT,R ,G ,B

end

;_____________________________________________________________________
;_____________________________________________________________________

PRO plot_6maps, ncat, tile_id, data, vlim, vname,advance = advance

lwval = vlim(0)
upval = vlim(1)
if (vname eq 'SOILDEPTH') then upval = 5000.

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

data_grid = fltarr (im,jm)
data_grid (*,*) =  !VALUES.F_NAN

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(tile_id[i,j] gt 0) then data_grid(i,j) =  data(tile_id[i,j] -1)
   endfor
endfor

limits = [-60,-180,90,180]
if file_test ('limits.idl') then restore,'limits.idl'

colors = [27,26,25,24,23,22,21,20,40,41,42,43,44,45,46,47,48]
n_levels = n_elements (colors)

levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]

if(vname eq 'POROS') then $
levels = [lwval,lwval+(0.57-lwval)/(n_levels -2) +indgen(n_levels -2)*(0.57-lwval)/(n_levels -2),upval]

if keyword_set (advance) then begin

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ADVANCE,/ISOTROPIC,/NOBORDER, title =vname
endif else begin
MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ISOTROPIC,/NOBORDER, title =vname
endelse

contour, data_grid,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot

levels_x = levels

if(vname eq 'POROS') then  begin
dxp = (0.8-0.37)/16.
levels_x = indgen(17)*dxp+ 0.37
endif

alpha=fltarr(n_levels,2)
alpha(*,0)=levels
alpha(*,1)=levels
h=[0,1]

dx = (240.)/(n_levels-1)

clev = levels
clev (*) = 1
n=0
k = 0
fmt_string = '(f7.2)'

if(vname eq 'CDCR2' ) then !P.position=[0.064, 0.675, 0.41, 0.69]
if(vname eq 'BEE'   ) then !P.position=[0.58, 0.675, 0.92, 0.69]
if(vname eq 'POROS' ) then !P.position=[0.064, 0.345, 0.41, 0.36]
if(vname eq 'CATDEF') then !P.position=[0.58, 0.345, 0.92, 0.36]
if(vname eq 'RZEXC' ) then !P.position=[0.064, 0.015, 0.41, 0.03]
if(vname eq 'SFEXC' ) then !P.position=[0.58, 0.015, 0.92, 0.03]

;!P.position=[0.064, 0.675, 0.41, 0.69]
;!P.position=[0.58, 0.0+0.005, 0.92, 0.015+0.005]

contour,alpha,levels_x,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels_x,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(colors) -1 do xyouts,levels_x[k],1.1,string(levels[k],format=fmt_string) ,orientation=90,color=0,charsize =0.8

;for l = 0,n_levels -2 do begin
;    k = l
;    xbox = [-120. + k*dx,-120. + k*dx, -120. + (k+1)*dx, -120. + (k+1)*dx,-120. + k*dx]
;    ybox = [-65., -55.,-55.,-65.,-65.]
;    polyfill, xbox,ybox,color=colors [k]
;
;    xyouts,xbox[1],ybox[2]+0.05,string(levels[l],format=fmt_string),color =0, orientation =90,charsize =0.8
;    k = k + 1     
;endfor
;
;l = n_levels -1
;xyouts,-120. + l*dx,ybox[2]+0.05,string(levels[l],format=fmt_string),color =0, orientation =90,charsize =0.8
!P.position=0

END
;_____________________________________________________________________
;_____________________________________________________________________

pro plot_tc, ncat, tile_id,out_path, VISDR, VISDF, NIRDR, NIRDF

load_colors
thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[720,500], Z_Buffer=0
Erase,255
!p.background = 255

!P.position=0
!P.Multi = [0, 2, 2, 0, 0]
limits = [-60,-180,90,180]

lwval = min ([min(VISDR), min(VISDF), min(NIRDR), min(NIRDF)])
upval = max ([max(VISDR), max(VISDF), max(NIRDR), max(NIRDF)])

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

colors = [27,26,25,24,23,22,21,20,40,41,42,43,44,45,46,47,48]
n_levels = n_elements (colors)

for map = 1,4 do begin

   if (map eq 1) then data = VISDR
   if (map eq 2) then data = VISDF	
   if (map eq 3) then data = NIRDR
   if (map eq 4) then data = NIRDF

   if (map eq 1) then ctitle = 'TC1'
   if (map eq 2) then ctitle = 'TC2'	
   if (map eq 3) then ctitle = 'TC3'
   if (map eq 4) then ctitle = 'TC4'

   levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]

   data_grid = fltarr (im,jm)
   data_grid (*,*) =  !VALUES.F_NAN

   for j = 0l, jm -1l do begin
       for i = 0l, im -1 do begin
           if(tile_id[i,j] gt 0) then data_grid(i,j) =  data(tile_id[i,j] -1)
       endfor
   endfor

   if(map eq 1) then begin
      MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ISOTROPIC,/NOBORDER, title = ctitle
   endif else begin
      MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ADVANCE,/ISOTROPIC,/NOBORDER, title = ctitle
   endelse

   contour, data_grid,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
   if(map eq 3) then begin
      !P.position=[0.25, 0.05, 0.75, 0.075]

      alpha=fltarr(n_levels,2)
      alpha(*,0)=levels
      alpha(*,1)=levels
      h=[0,1]		
      clev = levels
      clev (*) = 1
      n=0
      k = 0
      fmt_string = '(f6.2)'
      contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
      contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(colors) -1 do xyouts,levels[k],1.1,string(levels[k],format=fmt_string) ,orientation=90,color=0,charsize =0.8
   !P.position=0
   endif
endfor

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 720, 500)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, out_path +'TC.jpg', image24, True=1, Quality=100

end
; ==============================================================================
;                                  Mosaic classes    
; ==============================================================================

PRO plot_mosaic, ncat, outdir, tile_id, mos_type

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

mos_grid = intarr (im,jm)
mos_grid (*,*) =  !VALUES.F_NAN

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(tile_id[i,j] gt 0) then mos_grid(i,j) =  mos_type(tile_id[i,j] -1)
   endfor
endfor

limits = [-60,-180,90,180]
if file_test ('limits.idl') then restore,'limits.idl'

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[700,500], Z_Buffer=0

r_in  = [233,255,255,255,210,  0,  0,  0,204,170,255,220,205,  0,  0,170,  0, 40,120,140,190,150,255,255,  0,  0,  0,195,255,  0,255,  0]
g_in  = [ 23,131,191,255,255,255,155,  0,204,240,255,240,205,100,160,200, 60,100,130,160,150,100,180,235,120,150,220, 20,245, 70,255,  0]
b_in  = [  0,  0,  0,178,255,255,255,200,204,240,100,100,102,  0,  0,  0,  0,  0,  0,  0,  0,  0, 50,175, 90,120,130,  0,215,200,255,  0]
vtypes =[  1,  2,  3,  4,  5,  6,  7,  8, 10, 11, 14, 20, 30, 40, 50, 60, 70, 90,100,110,120,130,140,150,160,170,180,190,200,210,220,230]

red  = intarr (256)
green= intarr (256)
blue = intarr (256)
red  (255) = 255
green(255) = 255
blue (255) = 255

for k = 0, n_elements(vtypes) -1 do begin 
	red  (vtypes(k)) = r_in (k)
	green(vtypes(k)) = g_in (k)
	blue (vtypes(k)) = b_in (k)
endfor

TVLCT,red,green,blue

colors = vtypes
levels = vtypes

Erase,255
!p.background = 255
!P.position=0
!P.Multi = [0, 1, 1, 0, 1]
MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits
contour, mos_grid,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot

mos_name = strarr(6)
mos_name( 0)  = 'BL Evergreen' 
mos_name( 1)  = 'BL Deciduous' 
mos_name( 2)  = 'Needleleaf' 
mos_name( 3)  = 'Grassland' 
mos_name( 4)  = 'BL Shrubs' 
mos_name( 5)  = 'Dwarf' 

n_levels = 6;n_elements(vtypes)
alpha=fltarr(n_levels,2)
alpha(*,0)=levels [0:n_levels-1]
alpha(*,1)=levels [0:n_levels-1]
h=[0,1]
!P.position=[0.30, 0.0+0.005, 0.70, 0.015+0.005]
clev = levels
clev (*) = 1
contour,alpha,levels[0:5],h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
        /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[1,7], $
        xtitle=' ', color=0,xtickv=levels, $
        C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels[0:5],h,levels=levels,color=0,/overplot,c_label=clev
for k = 0,5 do xyouts,levels[k]+0.5,1.2,mos_name[k] ,orientation=90,color=0

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 700, 500)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, outdir +  '/mosaic_prim.jpg', image24, True=1, Quality=100


END
; ==============================================================================
;                                  CLM-Carbon classes    
; ==============================================================================

PRO plot_carbon,ncat, OutDir, tile_id, clm_type

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

clm_grid = intarr (im,jm,4)

clm_grid (*,*,*) =  !VALUES.F_NAN

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(tile_id[i,j] gt 0) then begin
	clm_grid(i,j,0) =  clm_type(tile_id[i,j] -1,0)
	clm_grid(i,j,1) =  clm_type(tile_id[i,j] -1,1)
	clm_grid(i,j,2) =  clm_type(tile_id[i,j] -1,2)
	clm_grid(i,j,3) =  clm_type(tile_id[i,j] -1,3)
      endif	
   endfor
endfor

clm_type = 0

limits = [-60,-180,90,180]
if file_test ('limits.idl') then restore,'limits.idl'

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[700,1000], Z_Buffer=0
;types= [  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,11a, 12, 13, 14,14a, 15,15a, 16,16a, 17]
r_in  = [106,202,251,  0, 29, 77,109,142,233,255,255,255,127,164,164,217,217,204,104,  0]
g_in  = [ 91,178,154, 85,115,145,165,185, 23,131,131,191, 39, 53, 53, 72, 72,204,104, 70]
b_in  = [154,214,153,  0,  0,  0,  0, 13,  0,  0,200,  0,  4,  3,200,  1,200,204,200,200]
vtypes= [  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

red  = intarr (256)
green= intarr (256)
blue = intarr (256)

red  (255) = 255
green(255) = 255
blue (255) = 255

for k = 0, n_elements(vtypes) -1 do begin 
	red  (vtypes(k)) = r_in (k)
	green(vtypes(k)) = g_in (k)
	blue (vtypes(k)) = b_in (k)
endfor

TVLCT,red,green,blue

colors = vtypes
levels = vtypes

clm_name = strarr(19)
clm_name( 0)  = 'NLEt'  ;  1 	needleleaf evergreen temperate tree
clm_name( 1)  = 'NLEB'  ;  2 	needleleaf evergreen boreal tree
clm_name( 2)  = 'NLDB'  ;  3  	needleleaf deciduous boreal tree
clm_name( 3)  = 'BLET'  ;  4 	broadleaf evergreen tropical tree
clm_name( 4)  = 'BLEt'  ;  5 	broadleaf evergreen temperate tree
clm_name( 5)  = 'BLDT'  ;  6 	broadleaf deciduous tropical tree
clm_name( 6)  = 'BLDt'  ;  7 	broadleaf deciduous temperate tree
clm_name( 7)  = 'BLDB'  ;  8 	broadleaf deciduous boreal tree
clm_name( 8)  = 'BLEtS' ;  9 	broadleaf evergreen temperate shrub
clm_name( 9)  = 'BLDtS' ; 10 	broadleaf deciduous temperate shrub [moisture +  deciduous]
clm_name(10)  = 'BLDtSm'; 11 	broadleaf deciduous temperate shrub [moisture stress only]
clm_name(11)  = 'BLDBS' ; 12 	broadleaf deciduous boreal shrub
clm_name(12)  = 'AC3G'  ; 13 	arctic c3 grass
clm_name(13)  = 'CC3G'  ; 14 	cool c3 grass [moisture +  deciduous]
clm_name(14)  = 'CC3Gm' ; 15 	cool c3 grass [moisture stress only]
clm_name(15)  = 'WC4G'  ; 16 	warm c4 grass [moisture +  deciduous]
clm_name(16)  = 'WC4Gm' ; 17 	warm c4 grass [moisture stress only]
clm_name(17)  = 'CROP'  ; 18 	crop [moisture +  deciduous]
clm_name(18)  = 'CROPm' ; 19 	crop [moisture stress only]

Erase,255
!p.background = 255
!P.position=0
!P.Multi = [0, 1, 2, 0, 1]

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits
contour, clm_grid[*,*,0],x,y,levels = levels,c_colors=colors,/cell_fill,/overplot

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/advance
contour, clm_grid[*,*,1],x,y,levels = levels,c_colors=colors,/cell_fill,/overplot

n_levels = n_elements(vtypes)
alpha=fltarr(n_levels,2)
alpha(*,0)=levels (0:n_levels-1)
alpha(*,1)=levels (0:n_levels-1)
h=[0,1]

!P.position=[0.30, 0.0+0.005, 0.70, 0.015+0.005]
clev = levels
clev (*) = 1
contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(vtypes) -2 do xyouts,levels[k]+0.5,1.1,clm_name(k) ,orientation=90,color=0
snapshot = TVRD()

TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 700, 1000)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, OutDir + '/CLM-Carbon_PRIM_veg_typs.jpg', image24, True=1, Quality=100

; now plotting secondary
thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[700,1000], Z_Buffer=0

Erase,255
!p.background = 255
!P.position=0
!P.Multi = [0, 1, 2, 0, 1]

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits
contour, clm_grid[*,*,2],x,y,levels = levels,c_colors=colors,/cell_fill,/overplot

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/advance
contour, clm_grid[*,*,3],x,y,levels = levels,c_colors=colors,/cell_fill,/overplot

n_levels = n_elements(vtypes)
alpha=fltarr(n_levels,2)
alpha(*,0)=levels (0:n_levels-1)
alpha(*,1)=levels (0:n_levels-1)
h=[0,1]

!P.position=[0.30, 0.0+0.005, 0.70, 0.015+0.005]
clev = levels
clev (*) = 1
contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(vtypes) -2 do xyouts,levels[k]+0.5,1.1,clm_name(k) ,orientation=90,color=0
snapshot = TVRD()

TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 700, 1000)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, OutDir + '/CLM-Carbon_SEC_veg_typs.jpg', image24, True=1, Quality=100

END
