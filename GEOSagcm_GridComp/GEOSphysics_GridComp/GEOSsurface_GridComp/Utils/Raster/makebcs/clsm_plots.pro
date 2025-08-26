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

; ----------------

function nint, x, LONG = long             ;Nearest Integer Function
;+
; NAME:
;	NINT
; PURPOSE:
;	Nearest integer function.
; EXPLANATION:   
;	NINT() is similar to the intrinsic ROUND function, with the following
;	two differences:
;	(1) if no absolute value exceeds 32767, then the array is returned as
;		as a type INTEGER instead of LONG
;	(2) NINT will work on strings, e.g. print,nint(['3.4','-0.9']) will
;		give [3,-1], whereas ROUND() gives an error message
;
; CALLING SEQUENCE:
;	result = nint( x, [ /LONG] )
;
; INPUT:
;	X - An IDL variable, scalar or vector, usually floating or double
;		Unless the LONG keyword is set, X must be between -32767.5 and 
;		32767.5 to avoid integer overflow
;
; OUTPUT
;	RESULT - Nearest integer to X
;
; OPTIONAL KEYWORD INPUT:
;	LONG - If this keyword is set and non-zero, then the result of NINT
;		is of type LONG.   Otherwise, the result is of type LONG if
;		any absolute values exceed 32767, and type INTEGER if all
;		all absolute values are less than 32767.
; EXAMPLE:
;	If X = [-0.9,-0.1,0.1,0.9] then NINT(X) = [-1,0,0,1]
;
; PROCEDURE CALL:
;	None:
; REVISION HISTORY:
;	Written W. Landsman        January 1989
;	Added LONG keyword         November 1991
;	Use ROUND if since V3.1.0  June 1993
;	Always start with ROUND function    April 1995
;	Return LONG values, if some input value exceed 32767
;		and accept string values   February 1998 
;       Use size(/TNAME) instead of DATATYPE()      October 2001
;-
xmax = max(x,min=xmin)
 xmax = abs(xmax) > abs(xmin)
 if (xmax gt 32767) or keyword_set(long) then begin
    if size(x,/TNAME) eq 'STRING' then b = round(float(x)) else b = round(x)
 end else begin
    if size(x,/TNAME) eq 'STRING' then b = fix(round(float(x))) else	$
	    b = fix(round(x))
 endelse

  return, b 
  end

; -------------------------------------------------------------------------------------------

FUNCTION IS_IN_DOMAIN, xylim, x,y

if (((x ge xylim(1)) and (x le xylim(3))) and  $
    ((y ge xylim(0)) and (y le xylim(2)))) then begin
 
   return_value = boolean(1) 

endif else begin

   return_value = boolean(0) 

endelse

return,return_value

END 

; #######################################################

FUNCTION Z0_VALUE, Z2CH, lai, SCALE4Z0

MIN_VEG_HEIGHT = 0.01
Z0_BY_ZVEG     = 0.13

if (SCALE4Z0  eq 2.) then begin
   return_value = SCALE4Z0 * Z0_BY_ZVEG * (Z2CH - (Z2CH - MIN_VEG_HEIGHT) * exp(-1.*LAI))
endif else begin
   return_value = Z0_BY_ZVEG * (Z2CH - SCALE4Z0 * (Z2CH - MIN_VEG_HEIGHT) * exp(-1.*LAI))
endelse

return,return_value

END

;===========================================================
;+
; NAME:
;        SHUFFLE
;
; PURPOSE:
;        This function returns the uniformly-shuffled elements of an array.
;
; CATEGORY:
;        Math.
;
; CALLING SEQUENCE:
;
;        Result = SHUFFLE( A [, Num])
;
; INPUTS:
;        A:        Array containing the elements to shuffle (e.g. INDGEN(100))
;
; OPTIONAL INPUTS:
;        Num:      Number of shuffled elements to return. Must be < N_ELEMENTS(A)+1
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;        SEED:     Number used to seed the random number generator, RANDOMU.
;
; OUTPUTS:
;        Returns the Num shuffled elements of the A array.
;
; OPTIONAL OUTPUT KEYWORD PARAMETERS:
;
;        INDICES:  Array of indices pointing to the shuffled elements of A.
;
; EXAMPLE:
;        Pick 10 unique random integers between the numbers 1..100:
;
;        i = INDGEN(100)
;        j = SHUFFLE(i,10)
;
; MODIFICATION HISTORY:
;        Written by:    Han Wen, January 1997.
;-
function SHUFFLE, A, Num, INDICES=Indices, SEED=Seed

         NP   = N_PARAMS()
         N    = N_ELEMENTS(A)
         if (N eq 0) then message, $
              'Must be called with 1-2 parameters: A [,Num]'
         if (NP eq 1) then Num = N

         r         = RANDOMU(Seed, N)
         Indices   = SORT(r)
         return, A(Indices(0:Num-1))
end

; #########################################################

PRO clsm_plots

; ##########################################################
; Calling Sequence:
; (1) get environment variables
; (2) reading catchment.def and setting map limits
; (3) generating NC_plot x NR_plot mask for plotting maps
; (4) plotting catchment-tiles in the Eastern United States
; (5) processing JPL Height
; (6) plotting CTI statistics
; (7) plotting vegetation types
; (8) plotting soil hydraulic properties
; (9) plotting elevation
; (10)plot LAI monthly climatology
; (11)generating NC_plot x NR_plot mask for plotting maps 
; (12)making  movies of Seasonal data
;
; Miscellaneous Routines
; (a1) check_satparam - Check ars and arw parameters 
; (a2) create_vec_file - For LIS/GSWP-2 type applications
; ##########################################################

; (1) Reading in Enviornment variables
; --------------------------------

gfile=GETENV('gfile')
path =GETENV('workdir')
NC   =1l*GETENV('NC')
NR   =1l*GETENV('NR')


; (2) Reading number of catchments
;---------------------------------

openr,1,'../catchment.def'
ncat = 0l
readf,1,ncat

if((stregex (gfile,'Pfafstetter') ge 0) or (stregex (gfile,'SMAP') ge 0)) then begin
; global plots
endif else begin

min_lon =  180.
max_lon = -180.
min_lat =   90.
max_lat =  -90.
a1 = 0.
a2 = 0.
a3 = 0.
a4 = 0.
k  = 0

for i = 0l,ncat -1l  do begin
   readf,1,k,k,a1, a2, a3, a4
   if (a1 lt min_lon) then min_lon = a1
   if (a2 gt max_lon) then max_lon = a2	
   if (a3 lt min_lat) then min_lat = a3
   if (a4 gt max_lat) then max_lat = a4
endfor 

limits = [floor(min_lat), floor(min_lon),ceil(max_lat),ceil(max_lon)]
if((ceil(max_lon) - floor(min_lon)) lt 180.) then save,limits,file ='limits.idl'

endelse

close,1

; (3) generating NC_plot x NR_plot mask for plotting maps
;--------------------------------------------------------

NC_plot = 4320
NR_plot = 2160

tile_id = lonarr (NC_plot,NR_plot)

dx = NC/NC_plot
dy = NR/NR_plot

catrow = lonarr(nc)
cat    = lonarr(nc,dy)

rst_file=path + '/rst/' + gfile+'*.rst'
openr,1,rst_file,/F77_UNFORMATTED

for j = 0l, NR_plot -1 do begin

   for i=0,dy -1 do begin
      readu,1,catrow
      cat (*,i) = catrow
   endfor

   for i = 0, NC_plot -1 do begin
      subset = cat (i*dx: (i+1)*dx -1,*)
      if (min (subset) le ncat) then begin
         min1 = min(subset)
         subset(where (subset gt ncat)) = 0
         hh = histogram(subset,bin=1,min = min1, locations=loc_val)
         dom_tile = max(hh,loc)	
         tile_id[i,j] = loc_val(loc)
      endif
   endfor

endfor

close,1

 
; (4) plotting catchment-tiles in the Eastern United States
;----------------------------------------------------------

plot_tiles,nc,nr,ncat,gfile,path

; plot countr_codes
country_codes, tile_id

; (5) Plot canopy height
; ----------------------

;canop_Height, nc,nr, tile_id, gfile, path

; (6) plotting CTI statistics
;----------------------------

filename = '../cti_stats.dat'
cti_mean = fltarr (ncat)
cti_std  = fltarr (ncat)
cti_skew = fltarr (ncat)

a1 = 0.
a2 = 0.
a3 = 0.
a4 = 0.
a5 = 0.
k =  0
openr,1,filename

readf,1,k

for i = 0l,ncat -1l  do begin
   readf,1,k,k,a1, a2, a3, a4, a5
   cti_mean (i) = a1
   cti_std  (i) = a2
   cti_skew (i) = a5
endfor

close,1
clm_file = '../CLM_veg_typs_fracs'
if (file_test (clm_file))  then begin
endif else begin
cti_mean = 0.961*cti_mean - 1.957  
endelse


;plot_vars2, ncat, tile_id, cti_mean, 'cti_mean'
;plot_vars2, ncat, tile_id, cti_std , 'cti_std'
;plot_vars2, ncat, tile_id, cti_skew, 'cti_skew'

plot_three_vars1, ncat, tile_id, cti_mean, cti_std, cti_skew

cti_mean = 0.
cti_std  = 0.
cti_skew = 0.

; (7) plotting vegetation types
;------------------------------

plot_mosaic, ncat, tile_id
clm_file = '../CLM_veg_typs_fracs'

if (file_test (clm_file))  then begin

  plot_clm   , ncat, tile_id
  plot_carbon, ncat, tile_id

; Now plot Ndep, T2m and SoilAlb
; ------------------------------

  filename = '../CLM_NDep_SoilAlb_T2m'    
  ndep  = fltarr (ncat)
  visdr = fltarr (ncat)
  visdf = fltarr (ncat)
  nirdr = fltarr (ncat)
  nirdf = fltarr (ncat)
  t2mm  = fltarr (ncat)
  t2mp  = fltarr (ncat)

  a1 = 0.
  a2 = 0.
  a3 = 0.
  a4 = 0.
  a5 = 0.
  a6 = 0.
  a7 = 0.
  
  openr,1,filename

  for i = 0l,ncat -1l  do begin
    readf,1,a1, a2, a3, a4, a5, a6, a7
    ndep (i) = a1
    visdr(i) = a2
    visdf(i) = a3
    nirdr(i) = a4
    nirdf(i) = a5
    t2mm (i) = a6
    t2mp (i) = a7

  endfor

  close,1
  plot_three_vars2, ncat, tile_id, ndep, t2mm, t2mp
  plot_soilalb, ncat, tile_id,VISDR, VISDF, NIRDR, NIRDF

  ndep  = 0.
  visdr = 0.
  visdf = 0.
  nirdr = 0.
  nirdf = 0.
  t2mm  = 0.
  t2mp  = 0.

endif	

; (8) plotting soil hydraulic properties
;---------------------------------------

filename = '../soil_param.dat'
bee  = fltarr (ncat)
psis = fltarr (ncat)
poro = fltarr (ncat)
cond = fltarr (ncat)
wwet = fltarr (ncat)
sdep = fltarr (ncat)
a1 = 0.
a2 = 0.
a3 = 0.
a4 = 0.
a5 = 0.
a6 = 0.
k =  0

openr,1,filename

for i = 0l,ncat -1l  do begin
   readf,1,k,k,k,k,a1, a2, a3, a4, a5, a6
   bee  (i) = a1
   psis (i) = a2
   poro (i) = a3
   cond (i) = a4
   wwet (i) = a5
   sdep (i) = a6
endfor

close,1

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[720,800], Z_Buffer=0
load_colors
Erase,255

!p.background = 255

!P.position=0
!P.Multi = [0, 2, 3, 0, 0]

plot_vars, ncat, tile_id, bee, [1.   ,8.], 'BEE'
plot_vars, ncat, tile_id, psis,[-1.85,-0.1],'PSIS',advance =1 
plot_vars, ncat, tile_id, poro,[0.37,0.8],'POROS',advance =1 
plot_vars, ncat, tile_id, cond,[2.37e-6,2.845e-4],'COND',advance =1 
plot_vars, ncat, tile_id, wwet,[0.01,0.45],'WPWET',advance =1 
plot_vars, ncat, tile_id, sdep,[1334.,5000.],'SOILDEPTH',advance =1 

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 720, 800)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'soil_param.jpg', image24, True=1, Quality=100

;plot_vars2, ncat, tile_id, bee, 'BEE'
;plot_vars2, ncat, tile_id, psis,'PSIS' 
;plot_vars2, ncat, tile_id, poro,'POROS'
;plot_vars2, ncat, tile_id, cond,'COND'
;plot_vars2, ncat, tile_id, wwet,'WPWET'
;plot_vars2, ncat, tile_id, sdep,'SOILDEPTH'

bee  = 0.
psis = 0.
poro = 0.
cond = 0.
wwet = 0.
sdep = 0.

; (9) plotting elevation
;-----------------------

filename = '../catchment.def'
elevation  = fltarr (ncat)

a1 = 0.
a2 = 0.
k =  0

openr,1,filename

readf,1,k

for i = 0l,ncat -1l  do begin
   readf,1,k,k,a1, a1, a1, a1, a2
   elevation (i) = a2
endfor

close,1

plot_vars2, ncat, tile_id, elevation, 'ELEVATION'

elevation = 0.

; (10) plot LAI monthly climatology
; --------------------------------

plot_lai, ncat, tile_id

; (11) vegetation height and roughness length
filename = '../vegdyn.data'
ncdf_file   = boolean (ncdf_isncdf(filename))

if (ncdf_file) then begin

   ncid = NCDF_OPEN(filename,/NOWRITE)
   NCDF_VARGET, ncid,'ITY', ITYP
   NCDF_VARGET, ncid,'Z2CH', Z2
   NCDF_VARGET, ncid,'ASCATZ0', ASZ0
   NCDF_CLOSE, ncid

endif else begin

  openr,1,filename,/F77_UNFORMATTED
  ityp  = fltarr (ncat)
  z2    = fltarr (ncat)
  asz0  = fltarr (ncat)

  readu,1,ITYP
  readu,1,Z2
  readu,1,ASZ0
  close,1

endelse

ASZ0 = ASZ0 * 1000.


plot_canoph, z2, tile_id

if (file_test ( '../CLM_veg_typs_fracs'))  then begin
SCALE4Z0  = 0.5
endif else begin
SCALE4Z0  = 2.
endelse

compute_zo,'ascat' , SCALE4Z0, ASZ0, Z2, tile_id
compute_zo,'icarus', SCALE4Z0, ASZ0, Z2, tile_id
compute_zo,'merged', SCALE4Z0, ASZ0, Z2, tile_id

; plotting irrigation parameters
; ------------------------------

;plot_crop_times, ncat, tile_id
;irrig_method, ncat, tile_id 
;plot_lai_minmax, ncat, tile_id
;irrig_fractions, ncat, tile_id


; (12) generating NC_plot x NR_plot mask for plotting maps
;--------------------------------------------------------

NC_movie = 720
NR_movie = 360
Ntiles_per_cell = 30
if(NC gt 8640) then Ntiles_per_cell = 800

vec_map = {NT:0, TID: lonarr (Ntiles_per_cell), TFrac : fltarr (Ntiles_per_cell)}
vec2grid = REPLICATE (vec_map,NC_movie,NR_movie)

dx  = NC/NC_movie
dy  = NR/NR_movie
cat = lonarr(nc,dy)
catrow = lonarr(nc)

rst_file=path + '/rst/' + gfile+'*.rst'
openr,1,rst_file,/F77_UNFORMATTED

for j = 0l, NR_movie -1 do begin
   for i=0,dy -1 do begin
      readu,1,catrow
      cat (*,i) = catrow
   endfor

   for i = 0, NC_movie -1 do begin
      subset = cat (i*dx: (i+1)*dx -1,*)
      cat_unq = Subset[uniq(Subset,sort(Subset))]
      k_land = where ((cat_unq ge 1l) and (cat_unq le ncat))
      if (max(k_land) ne -1) then begin
         hh      = intarr(n_elements (cat_unq)) 
         for k = 0,n_elements (cat_unq) -1 do hh[k] =  $
            n_elements(where (Subset eq  cat_unq[k]))         
         NCOUNT =  0

         for k = 0, n_elements (hh) -1 do begin
            if((cat_unq (k) ge 1) and (cat_unq (k) le ncat)) then begin
               if (NCOUNT eq Ntiles_per_cell) then begin
                  print, 'Increase Ntiles_per_cell'
                  stop
               endif
               vec2grid[i,j].NT = vec2grid[i,j].NT + 1
               vec2grid[i,j].TID  (NCOUNT) = cat_unq (k)
               vec2grid[i,j].TFrac(NCOUNT) = 1.*hh(k)/total(hh)
               NCOUNT =  NCOUNT + 1	
            endif
         endfor	
      endif
   endfor
endfor	

close,1

; (12) Making  movies of Seasonal data
;------------------------------------

make_movies, ncat, vec2grid, 'LAI'
make_movies, ncat, vec2grid, 'GREEN'
make_movies, ncat, vec2grid, 'VISDF'
make_movies, ncat, vec2grid, 'NIRDF'

END

; ==============================================================================
;                                  Catchment-CN classes    
; ==============================================================================

PRO plot_carbon,ncat, tile_id

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

clm_type = intarr (ncat,4)
clm_grid = intarr (im,jm,4)

filename = '../CLM_veg_typs_fracs'
openr,1,filename
k = 0
v = 0
fr= 0.
v1= 0
v2= 0
v3= 0
v4 =0

for i = 0l,ncat -1l  do begin
   readf,1,k,k,v1,v2,v3,v4,fr,fr,fr,fr,v,v
   clm_type(i,0) = v1
   clm_type(i,1) = v2	
   clm_type(i,2) = v3
   clm_type(i,3) = v4
endfor

close,1

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
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/advance
contour, clm_grid[*,*,1],x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

n_levels = n_elements(vtypes)
alpha=fltarr(n_levels,2)
alpha(*,0)=levels (0:n_levels-1)
alpha(*,1)=levels (0:n_levels-1)
h=[0,1]

!P.position=[0.15, 0.0+0.005, 0.85, 0.015+0.005]
clev = levels
clev (*) = 1
contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(clm_name) -1 do xyouts,levels[k]+0.5,1.1,clm_name(k) ,orientation=90,color=0
snapshot = TVRD()

TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 700, 1000)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'CatchmentCN_PRIM_veg_typs.jpg', image24, True=1, Quality=100

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
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/advance
contour, clm_grid[*,*,3],x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

n_levels = n_elements(vtypes)
alpha=fltarr(n_levels,2)
alpha(*,0)=levels (0:n_levels-1)
alpha(*,1)=levels (0:n_levels-1)
h=[0,1]

!P.position=[0.15, 0.0+0.005, 0.85, 0.015+0.005]
clev = levels
clev (*) = 1
contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(clm_name) -1 do xyouts,levels[k]+0.5,1.1,clm_name(k) ,orientation=90,color=0
snapshot = TVRD()

TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 700, 1000)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'CatchmentCN_SEC_veg_typs.jpg', image24, True=1, Quality=100

END

; ==============================================================================
;                                  CLM classes    
; ==============================================================================

PRO plot_clm,ncat, tile_id

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

clm_type = intarr (ncat,2)
clm_grid = intarr (im,jm,2)

filename = '../CLM_veg_typs_fracs'
openr,1,filename
k = 0
v = 0
fr= 0.
v1= 0
v2= 0

for i = 0l,ncat -1l  do begin
   readf,1,k,k,v,v,v,v,fr,fr,fr,fr,v1,v2
   clm_type(i,0) = v1
   clm_type(i,1) = v2	
endfor

close,1

clm_grid (*,*,*) =  !VALUES.F_NAN

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(tile_id[i,j] gt 0) then begin
	clm_grid(i,j,0) =  clm_type(tile_id[i,j] -1,0)
	clm_grid(i,j,1) =  clm_type(tile_id[i,j] -1,1)
      endif	
   endfor
endfor

clm_type = 0

limits = [-60,-180,90,180]
if file_test ('limits.idl') then restore,'limits.idl'

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[700,500], Z_Buffer=0

r_in  = [255,106,202,251,  0, 29, 77,109,142,233,255,255,127,164,217,204,  0]
g_in  = [245, 91,178,154, 85,115,145,165,185, 23,131,191, 39, 53, 72,204, 70]
b_in  = [215,154,214,153,  0,  0,  0,  0, 13,  0,  0,  0,  4,  3,  1,204,200]
vtypes= [  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17]

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

clm_name = strarr(16)
clm_name( 0)  = 'BARE' ;  1  	bare
clm_name( 1)  = 'NLEt' ;  2 	needleleaf evergreen temperate tree
clm_name( 2)  = 'NLEB' ;  3 	needleleaf evergreen boreal tree
clm_name( 3)  = 'NLDB' ;  4  	needleleaf deciduous boreal tree
clm_name( 4)  = 'BLET' ;  5 	broadleaf evergreen tropical tree
clm_name( 5)  = 'BLEt' ;  6 	broadleaf evergreen temperate tree
clm_name( 6)  = 'BLDT' ;  7 	broadleaf deciduous tropical tree
clm_name( 7)  = 'BLDt' ;  8 	broadleaf deciduous temperate tree
clm_name( 8)  = 'BLDB' ;  9 	broadleaf deciduous boreal tree
clm_name( 9)  = 'BLEtS'; 10 	broadleaf evergreen temperate shrub
clm_name(10)  = 'BLDtS'; 11 	broadleaf deciduous temperate shrub
clm_name(11)  = 'BLDBS'; 12 	broadleaf deciduous boreal shrub
clm_name(12)  = 'AC3G' ; 13 	arctic c3 grass
clm_name(13)  = 'CC3G' ; 14 	cool c3 grass
clm_name(14)  = 'WC4G' ; 15 	warm c4 grass
clm_name(15)  = 'CROP' ; 16 	crop

Erase,255
!p.background = 255
!P.position=0
!P.Multi = [0, 1, 1, 0, 1]

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits
contour, clm_grid[*,*,0],x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

n_levels = n_elements(vtypes)
alpha=fltarr(n_levels,2)
alpha(*,0)=levels (0:n_levels-1)
alpha(*,1)=levels (0:n_levels-1)
h=[0,1]

!P.position=[0.15, 0.0+0.005, 0.85, 0.015+0.005]
clev = levels
clev (*) = 1
contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(clm_name) -1 do xyouts,levels[k]+0.5,1.1,clm_name(k) ,orientation=90,color=0
snapshot = TVRD()

TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 700, 500)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'CLM_PRIM_veg_typs.jpg', image24, True=1, Quality=100

; now plotting secondary
thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[700,500], Z_Buffer=0

Erase,255
!p.background = 255
!P.position=0
!P.Multi = [0, 1, 1, 0, 1]

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits
contour, clm_grid[*,*,1],x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

n_levels = n_elements(vtypes)
alpha=fltarr(n_levels,2)
alpha(*,0)=levels (0:n_levels-1)
alpha(*,1)=levels (0:n_levels-1)
h=[0,1]

!P.position=[0.15, 0.0+0.005, 0.85, 0.015+0.005]
clev = levels
clev (*) = 1
contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(clm_name) -1 do xyouts,levels[k]+0.5,1.1,clm_name(k) ,orientation=90,color=0
snapshot = TVRD()

TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 700, 500)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'CLM_SEC_veg_typs.jpg', image24, True=1, Quality=100

END

; ==============================================================================
;                                  Make movies    
; ==============================================================================

PRO make_movies, ncat, vec2grid, vname

upval = 1.
lwval = 0.

if (vname eq 'LAI') then upval = 6.

if (vname eq 'LAI')   then filename = '../lai.dat'
if (vname eq 'GREEN') then filename = '../green.dat' 
if (vname eq 'VISDF') then filename = '../AlbMap.WS.8-day.tile.0.3_0.7.dat'
if (vname eq 'NIRDF') then filename = '../AlbMap.WS.8-day.tile.0.7_5.0.dat'

im = n_elements(vec2grid[*,0].NT)
jm = n_elements(vec2grid[0,*].NT)

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

limits = [-60,-180,90,180]
if file_test ('limits.idl') then restore,'limits.idl'

DEVICE, DECOMPOSED = 0

r_in  = [253,224,255,238,205,193,152,  0,124,  0,  0,  0,  0,  0,  0, 48,110, 85]
g_in  = [253,238,255,238,205,255,251,255,252,255,238,205,139,128,100,128,139,107]
b_in  = [253,224,  0,  0,  0,193,152,127,  0,  0,  0,  0,  0,  0,  0, 20, 61, 47]

n_levels = n_elements (r_in)

if (vname eq 'LAI')   then begin
   levels=[0.,0.25,0.5,0.75,1.,1.25,1.5,10. * indgen(11)*0.05+2.]
endif else begin
   levels=[0.,0.025,0.05,0.075,0.1,0.125,0.15,0.2,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0]
endelse

red  = intarr (256)
green= intarr (256)
blue = intarr (256)

red  (255) = 255
green(255) = 255
blue (255) = 255

for k = 0, N_levels -1 do begin 
	red  (k+1) = r_in (k)
	green(k+1) = g_in (k)
	blue (k+1) = b_in (k)
endfor

TVLCT,red,green,blue

colors = indgen (N_levels) + 1

compile_opt idl2

openr,1,filename,/F77_UNFORMATTED

yr=0.
mn =0.
dy =0.
dum =0.
yr1 =0.
mn1 =0.
dy1 =0.
yrg=0.
mng =0.
lai     = fltarr (im,jm)
lai1    = fltarr (im,jm)
lai2    = fltarr (im,jm)
lai_vec = fltarr (ncat)
lai1 [*,*] = !VALUES.F_NAN
lai2 [*,*] = !VALUES.F_NAN

alpha = fltarr(n_levels,2)
alpha [*,0] = levels
alpha [*,1] = levels
h = [0,1]
m_days = [31,28,31,30,31,30,31,31,30,31,30,31]

readu,1,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
readu,1,lai_vec 

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(vec2grid[i,j].NT gt 0) then begin
         for k = 0, vec2grid[i,j].NT -1 do begin 
            if (k eq 0) then lai1[i,j] = 0.
		      lai1[i,j] =  lai1[i,j]  + lai_vec[vec2grid[i,j].TID[k] -1]*vec2grid[i,j].TFrac[k]
         endfor	
      endif
   endfor
endfor

dofyr_b4 = ((float(julday(mn1,dy1,2001+yr1)-julday(12,31,2000))) - (float(julday(mn,dy,2001+yr)-julday(12,31,2000))))/2 + $
  float(julday(mn,dy,2001+yr)-julday(12,31,2000))

readu,1,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
dofyr_nxt = ((float(julday(mn1,dy1,2001+yr1)-julday(12,31,2000))) - (float(julday(mn,dy,2001+yr)-julday(12,31,2000))))/2 + $
  float(julday(mn,dy,2001+yr)-julday(12,31,2000))

readu,1,lai_vec 

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(vec2grid[i,j].NT gt 0) then begin
         for k = 0, vec2grid[i,j].NT -1 do begin 
            if (k eq 0) then lai2[i,j] = 0.
		      lai2[i,j] =  lai2[i,j]  + lai_vec[vec2grid[i,j].TID[k] -1]*vec2grid[i,j].TFrac[k]
         endfor	
      endif
   endfor
endfor

compile_opt idl2

video_file = vname+'.mp4'
video = idlffvideowrite(video_file)
framerate = 10
framedims = [750,512]
stream = video.addvideostream(framedims[0], framedims[1], framerate)
set_plot, 'z', /copy
device, set_resolution=framedims, set_pixel_depth=24, decomposed=0

for month = 1,12 do begin
    for day =1,m_days[month -1] do begin
       !P.position=0 
       dofyr_now = float(julday(month,day,2001+yr)-julday(12,31,2000))
       fac1 = (dofyr_now - dofyr_b4 )/(dofyr_nxt - dofyr_b4)
       fac2 = (dofyr_nxt - dofyr_now)/(dofyr_nxt - dofyr_b4)
       lai = fac1*lai2 + fac2*lai1

       dstamp =string(2001+fix(yr),'(i4.4)')+string(fix(month),'(i2.2)')+string(fix(day),'(i2.2)')
       Erase, Color= 255
       MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,title=vname + ':' + dstamp
       contour, lai,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
       MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2
 
       !P.position=[0.15, 0.0+0.005, 0.85, 0.015+0.005]
       contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
               /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
               xtitle=' ', color=255,xtickv=levels, $
               C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
       
       contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=levels
       if (vname eq 'LAI')   then begin
          for k = 0,n_elements(colors) -1 do xyouts,levels[k],1.1,string(levels[k],format='(f4.2)') ,orientation=90,color=0,charsize =0.8
       endif else begin
          for k = 0,n_elements(colors) -1 do xyouts,levels[k],1.1,string(levels[k],format='(f5.3)') ,orientation=90,color=0,charsize =0.8
       endelse

       timestamp = video.put(stream, tvrd(true=1))
       
       if(dofyr_now + 0.5 ge dofyr_nxt) then begin
          lai1 = lai2
          dofyr_b4 = dofyr_nxt 
          readu,1,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
          
          dofyr_nxt = ((float(julday(mn1,dy1,2001+yr1)-julday(12,31,2000))) - (float(julday(mn,dy,2001+yr)-julday(12,31,2000))))/2 + $
                      float(julday(mn,dy,2001+yr)-julday(12,31,2000))
          if((month eq 12) and (yr eq 2)) then yr = yr -1
          readu,1,lai_vec 

          for j = 0l, jm -1l do begin
             for i = 0l, im -1 do begin
                if(vec2grid[i,j].NT gt 0) then begin
                   for k = 0, vec2grid[i,j].NT -1 do begin 
                      if (k eq 0) then lai2[i,j] = 0.
                      lai2[i,j] =  lai2[i,j]  + lai_vec[vec2grid[i,j].TID[k] -1]*vec2grid[i,j].TFrac[k]
                   endfor	
                endif
             endfor
          endfor
       endif       
    endfor
 endfor    

close,1

device, /close
set_plot, strlowcase(!version.os_family) eq 'windows' ? 'win' : 'x'
video.cleanup

END

; ==============================================================================
;                                  Mosaic classes    
; ==============================================================================

PRO plot_mosaic, ncat, tile_id

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

mos_type = intarr (ncat)
mos_grid = intarr (im,jm)

filename = '../mosaic_veg_typs_fracs'
openr,1,filename
k = 0
v = 0
for i = 0l,ncat -1l  do begin
   readf,1,k,k,v
   mos_type(i) = v
endfor

close,1

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
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

mos_name = strarr(6)
mos_name( 0)  = 'BL Evergreen' 
mos_name( 1)  = 'BL Deciduous' 
mos_name( 2)  = 'Needleleaf' 
mos_name( 3)  = 'Grassland' 
mos_name( 4)  = 'BL Shrubs' 
mos_name( 5)  = 'Dwarf' 

n_levels = 6;n_elements(vtypes)
alpha=fltarr(n_levels+1,2)
alpha[*,0]=levels [0:n_levels]
alpha[*,1]=levels [0:n_levels]
h=[0,1]
!P.position=[0.30, 0.0+0.005, 0.70, 0.015+0.005]
clev = levels
clev (*) = 1
contour,alpha,levels[0:6],h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
        /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[1,7], $
        xtitle=' ', color=0,xtickv=levels, $
        C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels[0:6],h,levels=levels,color=0,/overplot,c_label=clev
for k = 0,5 do xyouts,levels[k]+0.5,1.2,mos_name[k] ,orientation=90,color=0

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 700, 500)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'mosaic_prim.jpg', image24, True=1, Quality=100


END

; ==============================================================================
;                    Catchment-tiles in the Eastern United States  
; ==============================================================================

PRO plot_tiles,nc,nr,ncat,gfile,path

dx=360./nc
dy=180./nr

glon = fltarr (nc)
glat = fltarr (nr)

for i = 0l,nc -1l do glon(i) = -180. + dx/2. + i*dx
for i = 0l,nr -1l do glat(i) =  -90. + dy/2. + i*dy

xylim = [35.,-82.,42.,-73]

if file_test ('limits.idl') then begin
 restore,'limits.idl'
 xylim = limits
endif

i1 = where((glon ge xylim(1)) and (glon lt xylim(1) + dx))
i2 = where((glon ge xylim(3)) and (glon lt xylim(3) + dx))
j1 = where((glat ge xylim(0)) and (glat lt xylim(0) + dy))
j2 = where((glat ge xylim(2)) and (glat lt xylim(2) + dy))

xlen=xylim(3)-xylim(1)
ylen=xylim(2)-xylim(0)
init=replicate(0.,xlen,ylen)
x = indgen(xlen)+xylim(1)
y = indgen(ylen)+xylim(0)

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[700,500], Z_Buffer=0
n_levels = 30
colors = indgen(30) + 90
load_colors
Erase,255
!p.background = 255

!P.position=0
!P.Multi = [0, 1, 1, 0, 1]

contour,init,x,y,title=tit,xrange=[x(0),x(xlen-1)+1.],yrange=[y(0),y(ylen-1)+1.],xstyle=1,ystyle=1,color =0

pfc =0l 
pfc1=0l 
pfcl=0l
pfcr=0l

cat =lonarr(nc)
catp=lonarr(nc)
rst_file=path + '/rst/' + gfile+'*.rst'
idum=0l
openr,1,rst_file,/F77_UNFORMATTED

for j= 0l,j2(0) do begin
 
   readu,1,cat
   if(j ge j1) then begin
      yu = -90. + j*dy + dy
      yl = -90. + j*dy 
      for i = i1(0),i2(0) do begin 
         pfc =cat(i)
         pfc1=cat(i)
         if((pfc ge 1) and (pfc le ncat)) then begin
            if(i ne 0)    then pfcl = cat(i-1)
            if(i ne nc-1) then pfcr = cat(i+1)
            if(j eq 0)then catp(i)=pfc
            xl= -180. + i*dx 
            xr= -180. + i*dx +dx
            xx=fltarr(5)
            yy=fltarr(5)
            xx=[xl,xl,xr,xr,xl]
            yy=[yu,yl,yl,yu,yu]
            n = pfc mod n_levels
            polyfill,xx,yy,color=colors(n)
            if(pfc ne catp(i)) then oplot,[xl,xr],[yl,yl],color =0
            if(pfc ne pfcl) then oplot,[xl,xl],[yl,yu],color =0
            if(pfc ne pfcr) then oplot,[xr,xr],[yl,yu],color =0
         endif
         catp(i)=pfc
      endfor
   endif
endfor
close,1
snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 700, 500)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'US-east.jpg', image24, True=1, Quality=100

END

;========================================================================
;                            Global maps
;========================================================================

PRO plot_vars, ncat, tile_id, data, vlim, vname,advance = advance

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

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ADVANCE,/ISOTROPIC,/NOBORDER
endif else begin
MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ISOTROPIC,/NOBORDER
endelse

contour, data_grid,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

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
fmt_string = '(f4.2)'
if(vname eq 'COND') then fmt_string = '(e8.2)'
if(vname eq 'SOILDEPTH') then fmt_string = '(i4)'
if(vname eq 'PSIS') then fmt_string = '(f5.2)'

if(vname eq 'BEE') then !P.position=[0.064, 0.675, 0.41, 0.69]
if(vname eq 'PSIS') then !P.position=[0.58, 0.675, 0.92, 0.69]
if(vname eq 'POROS') then !P.position=[0.064, 0.345, 0.41, 0.36]
if(vname eq 'COND') then !P.position=[0.58, 0.345, 0.92, 0.36]
if(vname eq 'WPWET') then !P.position=[0.064, 0.015, 0.41, 0.03]
if(vname eq 'SOILDEPTH') then  !P.position=[0.58, 0.015, 0.92, 0.03]

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

;________________________________________________________
;________________________________________________________
;________________________________________________________


PRO plot_vars2, ncat, tile_id, data, vname

lwval = min(data)
upval = max(data)
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

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[700,500], Z_Buffer=0

load_colors
colors = [27,26,25,24,23,22,21,20,40,41,42,43,44,45,46,47,48]
n_levels = n_elements (colors)
levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]

Erase,255
!p.background = 255

!P.position=0
!P.Multi = [0, 1, 1, 0, 1]

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,title = vname
contour, data_grid,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

alpha=fltarr(n_levels,2)
alpha(*,0)=levels
alpha(*,1)=levels
h=[0,1]
!P.position=[0.15, 0.0+0.005, 0.85, 0.015+0.005]
clev = levels
clev (*) = 1
contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
        /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
        xtitle=' ', color=0,xtickv=levels, $
        C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
if(vname eq 'COND') then begin
   for k = 0,n_elements(levels) -1 do xyouts,levels[k],1.1, string(levels[k],'(e10.3)'),orientation=90,color=0
endif else begin
   if ((vname eq 'SOILDEPTH') or (vname eq 'ELEVATION')) then begin
      for k = 0,n_elements(levels) -1 do xyouts,levels[k],1.1, string(levels[k],'(f5.0)'),orientation=90,color=0
   endif else begin
      for k = 0,n_elements(levels) -1 do xyouts,levels[k],1.1, string(levels[k],'(f5.2)'),orientation=90,color=0
   endelse
endelse

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 700, 500)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, vname + '.jpg', image24, True=1, Quality=100

END

;========================================================================
;                            Check ars and arw parameters
;========================================================================

PRO check_satparam

arw1=0.
arw2=0.
arw3=0.
arw4=0.
ars1=0.
ars2=0.
ars3=0.
cti_mean=0. 
cti_std =0.
cti_min =0.
cti_max =0.
cti_skew=0.
BEE  =0. 
PSIS =0.
POROS=0.
COND =0.
WPWET=0.
soildepth=0.
nbdep=0
nbdepl=0
wmin0=0.
cdcr1=0.
cdcr2=0.

file3='file.0000001'
openr,12,file3
readf,12,cti_mean, cti_std,cti_min, cti_max, cti_skew
readf,12,BEE, PSIS,POROS,COND,WPWET,soildepth
readf,12,nbdep,nbdepl,wmin0,cdcr1,cdcr2

catdef = fltarr(nbdep)
ar1    = fltarr(nbdep)
wmin   = fltarr(nbdep)

readf,12,catdef
readf,12,ar1
readf,12,wmin
readf,12,ars1,ars2,ars3
readf,12,arw1,arw2,arw3,arw4
close,12

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[700,500], Z_Buffer=0

load_colors
colors = [27,26,25,24,23,22,21,20,40,41,42,43,44,45,46,47,48]
n_levels = n_elements (colors)

Erase,255
!p.background = 255

!P.position=0
!P.Multi = [0, 1, 1, 0, 1]

ntot=nbdep
x=indgen(ntot)
y  = fltarr(ntot)

plot,x,y,xrange=[0.,max(catdef)],yrange=[0.,1.],linestyle=0,title='WMIN and AR1', color =0
oplot,catdef,wmin,color = 30
oplot,[cdcr1,cdcr1],[0.,1], color = 100
oplot,[cdcr2,cdcr2],[0.,1], color = 100

;ntot=fix(catdef(nbdep-1))+ 1.
ntot =fix(cdcr1)+ 1.
ntot2=fix(cdcr2)+ 1.

x=indgen(ntot2)
y=fltarr(ntot)
y2=fltarr(ntot2)
for n =0,ntot2-1 do begin

        if (n lt ntot) then y(n) =  arw4 + (1.- arw4)*(1.+ arw1*x(n))/(1.+ arw2*x(n)+  arw3*x(n)*x(n))
        y2(n) =  (1.+ ars1*x(n))/(1.+ ars2*x(n)+  ars3*x(n)*x(n))
endfor

oplot,x(0:ntot-1),y,color = 0, linestyle = 1
oplot,catdef,ar1,color = 220
oplot,x,y2,color = 0, linestyle = 1

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 700, 500)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG,'img.0000001.jpg', image24, True=1, Quality=100

end

;========================================================================
;                          Process JPL Canopy Height
;========================================================================

PRO canop_Height, nc, nr, tileid_plot, gfile, path

CanopH=read_tiff('/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/veg/veg_height/v1/Simard_Pinto_3DGlobalVeg_JGR.tif')
im=n_elements(CanopH(*,0))
jm=n_elements(CanopH(0,*))
CanopH = reverse(CanopH,2,/overwrite)

yh= dblarr(jm)
for i = 0l,jm -1l do yh(i) = i*1./120 -90.  + 1./240.
xh = dblarr(im)
for i = 0l,im -1l do xh(i) = i*1./120 -180. + 1./240.

N_tiles = 0l

openr,1,'../catchment.def'
readf,1,N_tiles
close,1 

canop_tiles = fltarr (N_tiles)
count_pix   = fltarr (N_tiles)

canop_tiles (*) = 0.01
count_pix   (*) = 0.

dx = IM/NC
dy = JM/NR

catrow  = lonarr (nc)
tile_id = lonarr (NC, nr)
rst_file= path + '/rst/' + gfile+'*.rst'

openr,1,rst_file,/F77_UNFORMATTED

for j = 0l, NR -1l do begin
   readu,1,catrow
   tile_id(*,j) = catrow(*)
   for i=0l, nc-1l do begin
      subset = CanopH (i*dx: (i+1)*dx -1,j*dy: (j+1)*dy -1)
      if ((catrow(i) ge 1) and (catrow(i) le N_tiles)) then begin 
        canop_tiles (catrow(i) -1) =  canop_tiles (catrow(i) -1) + mean (subset)
        count_pix   (catrow(i) -1) =  count_pix   (catrow(i) -1) + 1.
      endif
   endfor
endfor

close,1

canop_tiles (where (count_pix gt 0.)) = canop_tiles/count_pix
canop_tiles (where (canop_tiles lt 0.01)) = 0.01

openw,1,'Simard_Pinto_3DGlobalVeg_JGR.dat'

for k = 0l,n_tiles -1l do begin
   printf,1,format='(f7.3)',canop_tiles(k)
endfor

close,1

;openw,1,'../Simard_Pinto_3DGlobalVeg_JGR.bin',/F77_UNFORMATTED
;writeu,1,canop_tiles
;close,1
 
lwval = min(0.)

ip = n_elements(tileid_plot[*,0])
jp = n_elements(tileid_plot[0,*])

dx = 360. / ip
dy = 180. / jp

x = indgen(ip)*dx -180. +  dx/2.
y = indgen(jp)*dy -90.  +  dy/2.

data_grid = fltarr (ip,jp)
data_grid (*,*) =  !VALUES.F_NAN

for j = 0l, jp -1l do begin
   for i = 0l, ip -1 do begin
      if((tileid_plot[i,j] gt 0) and (tileid_plot[i,j] le N_tiles)) then data_grid(i,j) =  canop_tiles(tileid_plot[i,j] -1)
   endfor
endfor

upval = max(canop_tiles)

limits = [-60,-180,90,180]

colors = [27,26,25,24,23,22,21,20,40,41,42,43,44,45,46,47,48]
colors = reverse (colors)
n_levels = n_elements (colors)

levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[800,400], Z_Buffer=0
load_colors
Erase,255
!p.background = 255
!P.position=0
!P.Multi = [0, 1, 1, 0, 0]

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ISOTROPIC,/NOBORDER
contour, data_grid, x, y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

alpha=fltarr(n_levels,2)
alpha(*,0)=levels
alpha(*,1)=levels
h=[0,1]

dx = (240.)/(n_levels-1)

clev = levels
clev (*) = 1

  k = 0
   for l = 0,n_levels -2 do begin
      k = l
      xbox = [-120. + k*dx,-120. + k*dx, -120. + (k+1)*dx, -120. + (k+1)*dx,-120. + k*dx]
      ybox = [-65., -55.,-55.,-65.,-65.]
      polyfill, xbox,ybox,color=colors [k]
      xyouts,xbox[1],ybox[2]+0.05,string(levels[l],format='(f5.2)'),color =0, orientation =90,charsize =0.8
      k = k + 1     
   endfor

   for l = n_levels -1,n_levels -1 do begin
      xyouts,-120. + l*dx,ybox[2]+0.05,string(levels[l],format='(f5.2)'),color =0, orientation =90,charsize =0.8

   endfor	

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 800,400)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'Canopy_Height_onTiles.jpg', image24, True=1, Quality=100

spawn, "paste ../mosaic_veg_typs_fracs Simard_Pinto_3DGlobalVeg_JGR.dat > new_mos"
spawn, "/bin/mv new_mos ../mosaic_veg_typs_fracs"
spawn, "/bin/rm Simard_Pinto_3DGlobalVeg_JGR.dat"

tmp_data = read_ascii ("../mosaic_veg_typs_fracs")
openw,1,'../vegdyn.data',/F77_UNFORMATTED
writeu,1,tmp_data.field1(2,*)
writeu,1,tmp_data.field1(6,*)
close,1

END

;_____________________________________________________________________
;_____________________________________________________________________

PRO plot_three_vars2, ncat, tile_id, data1, data2, data3

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.
;stop
data_grid1 = fltarr (im,jm)
data_grid1 (*,*) =  !VALUES.F_NAN
data_grid2 = fltarr (im,jm)
data_grid2 (*,*) =  !VALUES.F_NAN
data_grid3 = fltarr (im,jm)
data_grid3 (*,*) =  !VALUES.F_NAN

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(tile_id[i,j] gt 0) then data_grid1(i,j) =  data1(tile_id[i,j] -1)
      if(tile_id[i,j] gt 0) then data_grid2(i,j) =  data2(tile_id[i,j] -1)
      if(tile_id[i,j] gt 0) then data_grid3(i,j) =  data3(tile_id[i,j] -1)
   endfor
endfor

limits = [-60,-180,90,180]
if file_test ('limits.idl') then restore,'limits.idl'

colors = [27,26,25,24,23,22,21,20,40,41,42,43,44,45,46,47,48]
n_levels = n_elements (colors)

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[720,900], Z_Buffer=0
load_colors
Erase,255
!p.background = 255

!P.position=0
!P.Multi = [0, 1, 3, 0, 1]

for n = 0,2 do begin

if (n eq 0) then begin
   upval = 350.
   lwval = 0.
   data = data_grid1
endif 

if (n eq 1) then begin
   upval = 300.
   lwval = 250.
   data = data_grid2
endif 

if (n eq 2) then begin
   upval =  300.
   lwval =  250.
   data = data_grid3
endif 

levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]
if(n eq 0) then levels = [indgen(15)*4.,65.,350.]
if(n eq 0) then MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/NOBORDER
if(n gt 0) then MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ADVANCE,/NOBORDER
contour, data,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

alpha=fltarr(n_levels,2)
alpha(*,0)=levels
alpha(*,1)=levels
h=[0,1]

dx = (240.)/(n_levels-1)

clev = levels
clev (*) = 1

  k = 0
   for l = 0,n_levels -2 do begin
      k = l
      xbox = [-120. + k*dx,-120. + k*dx, -120. + (k+1)*dx, -120. + (k+1)*dx,-120. + k*dx]
      ybox = [-65., -55.,-55.,-65.,-65.]
      polyfill, xbox,ybox,color=colors [k]
      xyouts,xbox[1],ybox[2]+0.05,string(levels[l],format='(f5.1)'),color =0, orientation =90,charsize =0.8
      k = k + 1     
   endfor

   for l = n_levels -1,n_levels -1 do begin
      xyouts,-120. + l*dx,ybox[2]+0.05,string(levels[l],format='(f5.1)'),color =0, orientation =90,charsize =0.8
   endfor	
	
endfor

;stop
snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 720, 900)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'CLM_Ndep_T2m.jpg', image24, True=1, Quality=100

END

;_____________________________________________________________________
;_____________________________________________________________________

pro plot_soilalb, ncat, tile_id,VISDR, VISDF, NIRDR, NIRDF


thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[720,500], Z_Buffer=0
load_colors
Erase,255
!p.background = 255

!P.position=0
!P.Multi = [0, 2, 2, 0, 0]
limits = [-60,-180,90,180]

lwval = 0.
upval = 0.65

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

   if (map eq 1) then ctitle = 'VISDR'
   if (map eq 2) then ctitle = 'VISDF'	
   if (map eq 3) then ctitle = 'NIRDR'
   if (map eq 4) then ctitle = 'NIRDF'

   if (map ge 3) then upval = 1.

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
   MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

   if((map eq 1) or (map eq 3)) then begin
      if(map eq 1) then !P.position=[0.25, 0.55, 0.75, 0.575]	
      if(map eq 3) then !P.position=[0.25, 0.05, 0.75, 0.075]

      alpha=fltarr(n_levels,2)
      alpha(*,0)=levels
      alpha(*,1)=levels
      h=[0,1]		
      clev = levels
      clev (*) = 1
      n=0
      k = 0
      fmt_string = '(f4.2)'
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
Write_JPEG, 'SoilAlb.jpg', image24, True=1, Quality=100

end

;_____________________________________________________________________
;_____________________________________________________________________

pro create_vec_file

dx         = 1.d0/12.
dy         = 1.d0/12.
DATELINE   = 1
global_bcs = 0
WORKDIR    = ''

nc = long(360./dx)
nr = long(180./dy)

openw,1,workdir + 'clsm/NLDAS-5arcmin_vec.data'

if(NOT (boolean (global_bcs))) then begin
  xylim = [35., -180., 80., -55.]
  x = indgen (nc)*dx -180. + dx/2.
  y = indgen (nr)*dy  -90. + dy/2.
  i1 = value_locate (x,  xylim(1)) + 1
  i2 = value_locate (x,  xylim(3)) 
  j1 = value_locate (y,  xylim(0)) + 1
  j2 = value_locate (y,  xylim(2)) 
  i_offset = i1
  j_offset = j1
  nc_domain = i2 - i1 + 1
  nr_domain = j2 - j1 + 1
  printf,1,format ='(2f8.4, i3, 4i5)', dx,dy, dateline, nc_domain,nr_domain,i_offset,j_offset	
endif else begin
  printf,1,dx,dy, dateline
endelse

SRTM_maxcat = 291284

nc_esa = 129600l
nr_esa = 64800l

nx = nc_esa/nc
ny = nr_esa/nr

ncid = NCDF_OPEN('/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/shared/mask/GEOS5_10arcsec_mask.nc')
;NCDF_VARGET, ncid,0, y
;NCDF_VARGET, ncid,1, x

n = 1l

subset = lonarr (nc_esa,ny)
for j = 0l,nr -1l do begin
    NCDF_VARGET, ncid,'CatchIndex', offset = [0,j*ny], count = [nc_esa,ny], SubSet
    for i = 0l,nc -1l do begin
;	NCDF_VARGET, ncid,'CatchIndex', offset = [i*nx,j*ny], count = [nx,ny], CatchIndex
	CatchIndex = SubSet(i*nx:(i+1)*nx -1,*)
	if(max(CatchIndex) gt SRTM_maxcat) then CatchIndex (where (CatchIndex gt SRTM_maxcat)) = 0
	if (max (CatchIndex) ge 1) then begin
          if(boolean (global_bcs)) then begin
	     printf,1,format ='(i7,2(1x,f10.5),2(1x,I5))',n,j*dy -90. + dy/2.,i*dx -180. + dx/2.,I+1,J+1
             n = n + 1
	  endif else begin
             if(IS_IN_DOMAIN(xylim, i*dx -180. + dx/2.,j*dy -90. + dy/2.)) then begin
                printf,1,format ='(i7,2(1x,f10.5),2(1x,I5))',n,j*dy -90. + dy/2.,i*dx -180. + dx/2.,I+1 - i_offset,J+1 - j_offset
	        n = n + 1
             endif
          endelse
         endif
    endfor
endfor
close,1
ncdf_close,ncid


end
;_________________________________________________________________
;_________________________________________________________________


PRO plot_three_vars1, ncat, tile_id, data1, data2, data3

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.
;stop
data_grid1 = fltarr (im,jm)
data_grid1 (*,*) =  !VALUES.F_NAN
data_grid2 = fltarr (im,jm)
data_grid2 (*,*) =  !VALUES.F_NAN
data_grid3 = fltarr (im,jm)
data_grid3 (*,*) =  !VALUES.F_NAN

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(tile_id[i,j] gt 0) then data_grid1(i,j) =  data1(tile_id[i,j] -1)
      if(tile_id[i,j] gt 0) then data_grid2(i,j) =  data2(tile_id[i,j] -1)
      if(tile_id[i,j] gt 0) then data_grid3(i,j) =  data3(tile_id[i,j] -1)
   endfor
endfor

limits = [-60,-180,90,180]
if file_test ('limits.idl') then restore,'limits.idl'

colors = [27,26,25,24,23,22,21,20,40,41,42,43,44,45,46,47,48]
n_levels = n_elements (colors)

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[720,900], Z_Buffer=0
Erase,255
!p.background = 255

!P.position=0
!P.Multi = [0, 1, 3, 0, 1]

for n = 0,2 do begin

if (n eq 0) then begin
   upval = 14.
   lwval = 6.
   data = data_grid1
endif 

if (n eq 1) then begin
   upval = 4.
   lwval = 0.
   data = data_grid2
endif 

if (n eq 2) then begin
   upval =  2.5
   lwval = -2.5
   data = data_grid3
endif 

levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]
if(n eq 0) then MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/NOBORDER
if(n gt 0) then MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ADVANCE,/NOBORDER
contour, data,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

alpha=fltarr(n_levels,2)
alpha(*,0)=levels
alpha(*,1)=levels
h=[0,1]

dx = (240.)/(n_levels-1)

clev = levels
clev (*) = 1

  k = 0
   for l = 0,n_levels -2 do begin
      k = l
      xbox = [-120. + k*dx,-120. + k*dx, -120. + (k+1)*dx, -120. + (k+1)*dx,-120. + k*dx]
      ybox = [-65., -55.,-55.,-65.,-65.]
      polyfill, xbox,ybox,color=colors [k]
      if (n eq 0) then xyouts,xbox[1],ybox[2]+0.05,string(levels[l],format='(f5.2)'),color =0, orientation =90,charsize =0.8
      if (n eq 1) then xyouts,xbox[1],ybox[2]+0.05,string(levels[l],format='(f5.2)'),color =0, orientation =90,charsize =0.8
      if (n eq 2) then xyouts,xbox[1],ybox[2]+0.05,string(levels[l],format='(f6.2)'),color =0, orientation =90,charsize =0.8
      k = k + 1     
   endfor

   for l = n_levels -1,n_levels -1 do begin
      if (n eq 0) then xyouts,-120. + l*dx,ybox[2]+0.05,string(levels[l],format='(f5.2)'),color =0, orientation =90,charsize =0.8
      if (n eq 1) then xyouts,-120. + l*dx,ybox[2]+0.05,string(levels[l],format='(f5.2)'),color =0, orientation =90,charsize =0.8
      if (n eq 2) then xyouts,-120. + l*dx,ybox[2]+0.05,string(levels[l],format='(f6.2)'),color =0, orientation =90,charsize =0.8
   endfor	
	
endfor

;stop
snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 720, 900)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'cti.jpg', image24, True=1, Quality=100

END
;_____________________________________________________________________
;_____________________________________________________________________

pro plot_lai, ncat, tile_id

lwval = 0.
upval = 7.

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

limits = [-60,-180,90,180]
if file_test ('limits.idl') then restore,'limits.idl'

r_in  = [253,224,255,238,205,193,152,  0,124,  0,  0,  0,  0,  0,  0, 48,110, 85]
g_in  = [253,238,255,238,205,255,251,255,252,255,238,205,139,128,100,128,139,107]
b_in  = [253,224,  0,  0,  0,193,152,127,  0,  0,  0,  0,  0,  0,  0, 20, 61, 47]

n_levels = n_elements (r_in)
levels=[0.,0.25,0.5,0.75,1.,1.25,1.5,10. * indgen(11)*0.05+2.]
red  = intarr (256)
green= intarr (256)
blue = intarr (256)

red  (255) = 255
green(255) = 255
blue (255) = 255

for k = 0, N_levels -1 do begin 
	red  (k+1) = r_in (k)
	green(k+1) = g_in (k)
	blue (k+1) = b_in (k)
endfor

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[1080,600], Z_Buffer=0

TVLCT,red,green,blue

colors = indgen (N_levels) + 1


Erase,255
!p.background = 255

!P.position=0
!P.Multi = [0, 3, 4, 0, 0]

file = '../lai.dat'

yr  = 0.
mn  = 0.
dy  = 0.
dum = 0.
yr1 = 0.
mn1 = 0.
dy1 = 0.
yrg = 0.
mng = 0.
lai  = fltarr (ncat)
lai1 = fltarr (ncat)
lai2 = fltarr (ncat)
mdays = [31,28,31,30,31,30,31,31,30,31,30,31]
mname = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

openr,1,file,/F77_UNFORMATTED
readu,1,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1

readu,1,lai1
dofyr_b4 = ((float(julday(mn1,dy1,2001+yr1)-julday(12,31,2000))) - (float(julday(mn,dy,2001+yr)-julday(12,31,2000))))/2 + $
  float(julday(mn,dy,2001+yr)-julday(12,31,2000))

readu,1,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1

dofyr_nxt = ((float(julday(mn1,dy1,2001+yr1)-julday(12,31,2000))) - (float(julday(mn,dy,2001+yr)-julday(12,31,2000))))/2 + $
  float(julday(mn,dy,2001+yr)-julday(12,31,2000))
readu,1,lai2

for month = 1,12 do begin
    lai_month  = fltarr (ncat)	
    for day =1,mdays[month -1] do begin
       dofyr_now = float(julday(month,day,2001+yr)-julday(12,31,2000))
       fac1 = (dofyr_now - dofyr_b4 )/(dofyr_nxt - dofyr_b4)
       fac2 = (dofyr_nxt - dofyr_now)/(dofyr_nxt - dofyr_b4)

       lai = fac1*lai2 + fac2*lai1
       lai_month(*) = lai_month(*) + lai (*)/mdays[month -1]

       if(dofyr_now + 0.5 ge dofyr_nxt) then begin
          lai1 = lai2
          dofyr_b4 = dofyr_nxt 
          readu,1,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1

          dofyr_nxt = ((float(julday(mn1,dy1,2001+yr1)-julday(12,31,2000))) - (float(julday(mn,dy,2001+yr)-julday(12,31,2000))))/2 + $
                      float(julday(mn,dy,2001+yr)-julday(12,31,2000))
          if((month eq 12) and (yr eq 2)) then yr = yr -1
          readu,1,lai2
       endif  	
    endfor

    data_grid = fltarr (im,jm)
    data_grid (*,*) =  !VALUES.F_NAN

    for j = 0l, jm -1l do begin
        for i = 0l, im -1 do begin
            if(tile_id[i,j] gt 0) then data_grid(i,j) =  lai_month(tile_id[i,j] -1)
        endfor
    endfor
    if(month gt 1) then begin
    MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ADVANCE,/NOBORDER,title=mname(month-1)
	endif else begin
    MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/NOBORDER,title=mname(month-1)
    endelse

    contour, data_grid,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
    MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2
endfor    

close,1	

!P.position=[0.15, 0.005, 0.85, 0.025]

alpha=fltarr(n_levels,2)
alpha(*,0)=levels
alpha(*,1)=levels
h=[0,1]

dx = (240.)/(n_levels-1)

clev = levels
clev (*) = 1
contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(colors) -1 do xyouts,levels[k],1.1,string(levels[k],format='(f4.2)') ,orientation=90,color=0,charsize =0.8
!P.position=0
snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 1080, 600)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'lai.jpg', image24, True=1, Quality=100

!P.Multi = 0
!P.position=0

end
 
;_____________________________________________________________________

pro load_random_colors

R = shuffle(indgen(256))
G = shuffle(indgen(256))
B = shuffle(indgen(256))

R (255) = 255
G (255) = 255
B (255) = 255
R (0) = 0
G (0) = 0
B (0) = 0

TVLCT,R ,G ,B

end



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

r_lct2  = [  0,   0,   0,   0,   0,   0,   0,   0,   0,  55, 120, 190, 240, 255, 255, 255, 255, 255, 233, 197, 158]
g_lct2  = [  0,   0,   0,  83, 115, 167, 195, 227, 255, 255, 255, 255, 255, 219, 187, 159, 131,  51,  23,   0,   0]
b_lct2  = [130, 200, 255, 255, 255, 255, 255, 255, 255, 199, 135,  67,  15,   0,   0,   0,   0 ,  0,   0,   0,   0]

r_veg  = [233,255,255,255,210,  0,  0,  0,204,170,255,220,205,  0,  0,170,  0, 40,120,140,190,150,255,255,  0,  0,  0,195,255,  0]
g_veg  = [ 23,131,191,255,255,255,155,  0,204,240,255,240,205,100,160,200, 60,100,130,160,150,100,180,235,120,150,220, 20,245, 70]
b_veg  = [  0,  0,  0,178,255,255,255,200,204,240,100,100,102,  0,  0,  0,  0,  0,  0,  0,  0,  0, 50,175, 90,120,130,  0,215,200]

r_grads_rb = [160, 110, 30,   0,   0,   0,   0, 160, 230, 230, 240, 250, 240]
g_grads_rb = [  0,  0,  60, 150, 200, 210, 220, 230, 220, 175, 130,  60,   0] 
b_grads_rb = [200, 220, 255, 255, 200, 140,  0,  50,  50,  45,  40,  60, 130] 

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

R (90:119) = r_veg
G (90:119) = g_veg
B (90:119) = b_veg

R (120:132) = r_grads_rb
G (120:132) = g_grads_rb
B (120:132) = b_grads_rb

R (140:160) = r_lct2
G (140:160) = g_lct2
B (140:160) = b_lct2
TVLCT,R ,G ,B

end

; -----------------------------------------------------------------------

pro jpl_tif2nc4

CanopH=read_tiff('/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/veg/veg_height/v1/Simard_Pinto_3DGlobalVeg_JGR.tif')
im=n_elements(CanopH(*,0))
jm=n_elements(CanopH(0,*))
CanopH = reverse(CanopH,2,/overwrite)

yh= dblarr(jm)
for i = 0l,jm -1l do yh(i) = i*1./120 -90.  + 1./240.
xh = dblarr(im)
for i = 0l,im -1l do xh(i) = i*1./120 -180. + 1./240.

id   = NCDF_CREATE('/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/veg/veg_height/v1/Simard_Pinto_3DGlobalVeg_JGR.nc4', /clobber, /NETCDF4_FORMAT) 
xid  = NCDF_DIMDEF(id, 'N_lon' , im)     ;Define x-dimension
yid  = NCDF_DIMDEF(id, 'N_lat' , jm)     ;Define y-dimension
NCDF_ATTPUT,id, 'CreatedBy', 'NASA GSFC GMAO Land Group',/global
NCDF_ATTPUT,id, 'Contact',   'NASA GSFC GMAO Land Group',/global
   
str_date=systime()
NCDF_ATTPUT,id, 'Date', str_date,/global 
vid  = NCDF_VARDEF(id,'longitude' , [xid], /DOUBLE)
vid  = NCDF_VARDEF(id,'latitude' ,  [yid], /DOUBLE)
vid  = NCDF_VARDEF(id,'CanopyHeight',[xid,yid], /SHORT)

NCDF_CONTROL, id, /ENDEF

NCDF_VARPUT, id,'longitude',xh
NCDF_VARPUT, id,'latitude', yh

for j = 0, jm -1 do begin
 NCDF_VARPUT, id,'CanopyHeight',offset=[0,j],count=[im,1],CanopH(*,j)
endfor

NCDF_CLOSE, id

end

; -----------------------------------------------------------------------

pro plot_canoph, z2, tileid_plot
N_tiles = 0l

openr,1,'../catchment.def'
readf,1,N_tiles
close,1 
ip = n_elements(tileid_plot[*,0])
jp = n_elements(tileid_plot[0,*])

dx = 360. / ip
dy = 180. / jp

x = indgen(ip)*dx -180. +  dx/2.
y = indgen(jp)*dy -90.  +  dy/2.

data_grid = fltarr (ip,jp)
data_grid (*,*) =  !VALUES.F_NAN

for j = 0l, jp -1l do begin
   for i = 0l, ip -1 do begin
      if((tileid_plot[i,j] gt 0) and (tileid_plot[i,j] le N_tiles)) then data_grid(i,j) =  z2(tileid_plot[i,j] -1)
   endfor
endfor

upval = max(z2)
lwval = min(z2)

limits = [-60,-180,90,180]

colors = [27,26,25,24,23,22,21,20,40,41,42,43,44,45,46,47,48]
colors = reverse (colors)
n_levels = n_elements (colors)

levels = [lwval,lwval+(upval-lwval)/(n_levels -1) +indgen(n_levels -1)*(upval-lwval)/(n_levels -1)]

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[800,400], Z_Buffer=0
load_colors
Erase,255
!p.background = 255
!P.position=0
!P.Multi = [0, 1, 1, 0, 0]

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ISOTROPIC,/NOBORDER
contour, data_grid, x, y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

alpha=fltarr(n_levels,2)
alpha(*,0)=levels
alpha(*,1)=levels
h=[0,1]

dx = (240.)/(n_levels-1)

clev = levels
clev (*) = 1

  k = 0
   for l = 0,n_levels -2 do begin
      k = l
      xbox = [-120. + k*dx,-120. + k*dx, -120. + (k+1)*dx, -120. + (k+1)*dx,-120. + k*dx]
      ybox = [-65., -55.,-55.,-65.,-65.]
      polyfill, xbox,ybox,color=colors [k]
      xyouts,xbox[1],ybox[2]+0.05,string(levels[l],format='(f5.2)'),color =0, orientation =90,charsize =0.8
      k = k + 1     
   endfor

   for l = n_levels -1,n_levels -1 do begin
      xyouts,-120. + l*dx,ybox[2]+0.05,string(levels[l],format='(f5.2)'),color =0, orientation =90,charsize =0.8

   endfor	

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 800,400)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'Canopy_Height_onTiles.jpg', image24, True=1, Quality=100

end

; ====================================================================================
pro country_codes, tile_id

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.
cnt_grid = intarr (im,jm)
cnt_grid (*,*) =  !VALUES.F_NAN

N_tiles = 0l
openr,1,'../catchment.def'
readf,1,N_tiles
close,1 
cnt_code = intarr (N_tiles)
st_code  = intarr (N_tiles)

openr,1,"../country_and_state_code.data"
k = 0l
i1 = 0
i2 = 0

for n = 0l, N_Tiles -1l do begin
readf,1,k,i1,i2
cnt_code (n) = i1
st_code  (n) = i2
endfor
close,1
us_ind   = where (cnt_code eq 243)
cnt_code (us_ind) = st_code (us_ind)
cnt_code (where (cnt_code eq 257)) =  !VALUES.F_NAN
tmp_data = 0

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(tile_id[i,j] gt 0) then begin
	cnt_grid(i,j) =  cnt_code(tile_id[i,j] -1) + 1
      endif	
   endfor
endfor

colors = indgen(256)
levels = colors
limits = [-60,-180,90,180]
thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[800,400], Z_Buffer=0
load_random_colors
Erase,255
!p.background = 255
!P.position=0
!P.Multi = [0, 1, 1, 0, 0]

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ISOTROPIC,/NOBORDER
for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if((cnt_grid[i,j] gt 0) and (cnt_grid[i,j] le 255)) then begin
	   yu = -90. + j*dy + dy
           yl = -90. + j*dy 
	   xl= -180. + i*dx 
           xr= -180. + i*dx +dx
           xx=fltarr(5)
           yy=fltarr(5)
           xx=[xl,xl,xr,xr,xl]
           yy=[yu,yl,yl,yu,yu]
           oplot,[xl,xr],[yl,yl],color =cnt_grid[i,j]           
      endif
   endfor
endfor

MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2
snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 800,400)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, 'Country_codes.jpg', image24, True=1, Quality=100
end
; ====================================================================================

pro compute_zo, pname, SCALE4Z0, ASZ0, Z2CH, tile_id

ncat = n_elements (Z2CH)

; Reading LAI and computing Z0
; ----------------------------

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

mdays     = [31,28,31,30,31,30,31,31,30,31,30,31]

zo_vec    = fltarr (ncat,   4)
ndvi_vec  = fltarr (ncat,   4)
ZOT       = fltarr (NCAT)

if (pname eq 'ascat') then goto, skip_lai 

lai_file  = '../lai.dat'
ndvi_file = '../ndvi.dat'

yr  = 0.
mn  = 0.
dy  = 0.
dum = 0.
yr1 = 0.
mn1 = 0.
dy1 = 0.
yrg = 0.
mng = 0.
lai  = fltarr (ncat)
lai1 = fltarr (ncat)
lai2 = fltarr (ncat)
ndvi  = fltarr (ncat)
ndvi1 = fltarr (ncat)
ndvi2 = fltarr (ncat)

openr,1,lai_file,/F77_UNFORMATTED
openr,2,ndvi_file,/F77_UNFORMATTED

readu,1,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
readu,1,lai1
dofyr_b4 = ((float(julday(mn1,dy1,2001+yr1)-julday(12,31,2000))) - (float(julday(mn,dy,2001+yr)-julday(12,31,2000))))/2 + $
  float(julday(mn,dy,2001+yr)-julday(12,31,2000))

readu,1,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
dofyr_nxt = ((float(julday(mn1,dy1,2001+yr1)-julday(12,31,2000))) - (float(julday(mn,dy,2001+yr)-julday(12,31,2000))))/2 + $
  float(julday(mn,dy,2001+yr)-julday(12,31,2000))
readu,1,lai2

readu,2,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
readu,2,ndvi1
ndvi_b4 = ((float(julday(mn1,dy1,2001+yr1)-julday(12,31,2000))) - (float(julday(mn,dy,2001+yr)-julday(12,31,2000))))/2 + $
  float(julday(mn,dy,2001+yr)-julday(12,31,2000))

readu,2,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
ndvi_nxt = ((float(julday(mn1,dy1,2001+yr1)-julday(12,31,2000))) - (float(julday(mn,dy,2001+yr)-julday(12,31,2000))))/2 + $
  float(julday(mn,dy,2001+yr)-julday(12,31,2000))
readu,2,ndvi2


for month = 1,12 do begin
    lai_month  = fltarr (ncat)	
    for day =1,mdays[month -1] do begin
       dofyr_now = float(julday(month,day,2001+yr)-julday(12,31,2000))
       fac1 = (dofyr_now - dofyr_b4 )/(dofyr_nxt - dofyr_b4)
       fac2 = (dofyr_nxt - dofyr_now)/(dofyr_nxt - dofyr_b4)
       lai = fac1*lai2 + fac2*lai1
       
       fac1 = (dofyr_now - ndvi_b4 )/(ndvi_nxt - ndvi_b4)
       fac2 = (ndvi_nxt - dofyr_now)/(ndvi_nxt - ndvi_b4)
       ndvi = fac1*ndvi2 + fac2*ndvi1	

; ==========================================================================================
; Here is the roughness length parameterization
; ==========================================================================================

	for n = 0l,ncat -1l do ZOT(n)  = Z0_VALUE(Z2CH(n), lai(n), SCALE4Z0)

       if((month eq  1) or (month eq  2) or (month eq 12)) then zo_vec (*,0) = zo_vec (*,0) + zot (*)/total (mdays ([11, 0, 1]))
       if((month eq  3) or (month eq  4) or (month eq  5)) then zo_vec (*,1) = zo_vec (*,1) + zot (*)/total (mdays ([ 2, 3, 4]))
       if((month eq  6) or (month eq  7) or (month eq  8)) then zo_vec (*,2) = zo_vec (*,2) + zot (*)/total (mdays ([ 5, 6, 7]))
       if((month eq  9) or (month eq 10) or (month eq 11)) then zo_vec (*,3) = zo_vec (*,3) + zot (*)/total (mdays ([ 8, 9,10]))

       if((month eq  1) or (month eq  2) or (month eq 12)) then ndvi_vec (*,0) = ndvi_vec (*,0) + ndvi (*)/total (mdays ([11, 0, 1]))
       if((month eq  3) or (month eq  4) or (month eq  5)) then ndvi_vec (*,1) = ndvi_vec (*,1) + ndvi (*)/total (mdays ([ 2, 3, 4]))
       if((month eq  6) or (month eq  7) or (month eq  8)) then ndvi_vec (*,2) = ndvi_vec (*,2) + ndvi (*)/total (mdays ([ 5, 6, 7]))
       if((month eq  9) or (month eq 10) or (month eq 11)) then ndvi_vec (*,3) = ndvi_vec (*,3) + ndvi (*)/total (mdays ([ 8, 9,10]))

       if(dofyr_now + 0.5 ge dofyr_nxt) then begin
          lai1 = lai2
          dofyr_b4 = dofyr_nxt 
          readu,1,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
          dofyr_nxt = ((float(julday(mn1,dy1,2001+yr1)-julday(12,31,2000))) - (float(julday(mn,dy,2001+yr)-julday(12,31,2000))))/2 + $
                      float(julday(mn,dy,2001+yr)-julday(12,31,2000))
          if((month eq 12) and (yr eq 2)) then yr = yr -1
          readu,1,lai2
       endif 

       if(dofyr_now + 0.5 ge ndvi_nxt) then begin
          ndvi1 = ndvi2
          ndvi_b4 = ndvi_nxt 
          readu,2,yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
          ndvi_nxt = ((float(julday(mn1,dy1,2001+yr1)-julday(12,31,2000))) - (float(julday(mn,dy,2001+yr)-julday(12,31,2000))))/2 + $
                      float(julday(mn,dy,2001+yr)-julday(12,31,2000))
          if((month eq 12) and (yr eq 2)) then yr = yr -1
          readu,2,ndvi2
       endif 
 	
    endfor	
endfor    

close,1	
close,2

skip_lai:

; now plotting
;-------------

sea_label = ['DJF','MAM','JJA','SON']

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[720,600], Z_Buffer=0
;Device, Set_Resolution=[720,500], Z_Buffer=0
load_colors
Erase,255
!p.background = 255

!P.position=0
!P.Multi = [0, 1, 2, 0, 0]
;!P.Multi = [0, 1, 1, 0, 0]
;!P.Multi = [0, 2, 2, 0, 0]
limits = [-60,-180,90,180]

colors = [74,77,35,34,33,32,25,24,23,22,21,20,41,42,43,44,45,46,47,48]
levels = [0.02,0.05,0.07,0.1,0.3,0.5,1,2,4,6,8,10,50,100,500,1000,2000,3000,4000,5000]
n_levels = n_elements (levels)

for season = 0,3,2 do begin

   data_grid = fltarr (im,jm)
   data_grid (*,*) =  !VALUES.F_NAN
   
   for j = 0l, jm -1l do begin
      for i = 0l, im -1 do begin
         if(tile_id[i,j] gt 0) then begin
            if (pname eq 'ascat')  then data_grid(i,j) =  asz0(tile_id[i,j] -1)
            if (pname eq 'icarus') then data_grid(i,j) =  1000.*zo_vec(tile_id[i,j] -1, season)
            if (pname eq 'merged')  then begin 
               data_grid(i,j) =  1000.*zo_vec(tile_id[i,j] -1, season)
               if(ndvi_vec(tile_id[i,j]-1, season) le 0.2) then data_grid(i,j) =  asz0(tile_id[i,j] -1)
            endif
         endif
      endfor
   endfor

    if(season ge 1) then begin
        MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/ADVANCE,/NOBORDER,title=pname +' : '+ sea_label (season)
    endif else begin
        MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits,/NOBORDER,title=pname +' : '+ sea_label (season) 
    endelse

    contour, data_grid,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot

    alpha=fltarr(n_levels,2)
    alpha(*,0)=levels
    alpha(*,1)=levels
    h=[0,1]		
    clev = levels
    clev (*) = 1
    dx = (240.)/(n_levels-1)

    clev = levels
    clev (*) = 1

    k = 0
    for l = 0,n_levels -2 do begin
        k = l
        xbox = [-120. + k*dx,-120. + k*dx, -120. + (k+1)*dx, -120. + (k+1)*dx,-120. + k*dx]
        ybox = [-65., -55.,-55.,-65.,-65.]
        polyfill, xbox,ybox,color=colors [k]
        if (l le 5) then xyouts,xbox[1],ybox[2]+0.05,string(levels[l],format='(f4.2)'),color =0, orientation =90,charsize =0.8
        if (l gt 5) then xyouts,xbox[1],ybox[2]+0.05,string(levels[l],format='(i4)'),color =0, orientation =90,charsize =0.8
        k = k + 1     
    endfor

    for l = n_levels -1,n_levels -1 do begin
        xyouts,-120. + l*dx,ybox[2]+0.05,string(levels[l],format='(i4)'),color =0, orientation =90,charsize =0.8
    endfor

endfor

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 720, 600)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_JPEG, pname + '_Z0.jpg', image24, True=1, Quality=100

end

; ------------------------------------

;; This program is being deprecated, we don't have IDATA valid path 
;pro proc_glass
;
;
;;IDATA = '/gpfsm/dnb43/projects/p03/RS_DATA/GLASS/LAI/AVHRR/V4/HDF/'
;;ODATA = '/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/veg/lai_grn/v4/GLASS-LAI/AVHRR.v4/'
;;LABEL = 'GLASS01B02.V04.A'
;;yearb = 1981
;;YEARe = 2017
;
;IDATA = '/gpfsm/dnb43/projects/p03/RS_DATA/GLASS/LAI/MODIS/V4/HDF/'
;ODATA = '/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/veg/lai_grn/v4/GLASS-LAI/MODIS.v4/'
;LABEL = 'GLASS01B01.V04.A'
;yearb = 2000
;YEARe = 2017
;
;nc = 7200
;nr = 3600
;nyrs = YEARe - yearb + 1
;
;lwval = 0.
;upval = 7.
;
;im = nc
;jm = nr
;
;dx = 360. / im
;dy = 180. / jm
;
;x = indgen(im)*dx -180. +  dx/2.
;y = indgen(jm)*dy -90.  +  dy/2.
;r_in  = [253,224,255,238,205,193,152,  0,124,  0,  0,  0,  0,  0,  0, 48,110, 85]
;g_in  = [253,238,255,238,205,255,251,255,252,255,238,205,139,128,100,128,139,107]
;b_in  = [253,224,  0,  0,  0,193,152,127,  0,  0,  0,  0,  0,  0,  0, 20, 61, 47]
;
;n_levels = n_elements (r_in)
;levels=[0.,0.25,0.5,0.75,1.,1.25,1.5,10. * indgen(11)*0.05+2.]
;red  = intarr (256)
;green= intarr (256)
;blue = intarr (256)
;
;red  (255) = 255
;green(255) = 255
;blue (255) = 255
;
;for k = 0, N_levels -1 do begin 
;	red  (k+1) = r_in (k)
;	green(k+1) = g_in (k)
;	blue (k+1) = b_in (k)
;endfor
;thisDevice = !D.Name
;set_plot,'Z'
;Device, Set_Resolution=[800,500], Z_Buffer=0
;TVLCT,red,green,blue
;colors = indgen (N_levels) + 1
;
;limits = [-60,-180,90,180]
;
;Erase,255
;!p.background = 255
;
;;file1 = '/gpfsm/dnb43/projects/p03/RS_DATA/GLASS/LAI/MODIS/V4/HDF/2008/GLASS01B01.V04.A2008185.hdf'
;file1 = '/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/veg/lai_grn/v4/GLASS-LAI/MODIS.v4/GLASS01B01.V04.AYYYY105.nc4'
;ncid = ncdf_open (file1)
;NCDF_VARGET, ncid,'LAI', adum
;ncdf_close,ncid
;
;;FileID=HDF_SD_Start(file1, /read)
;;sds_id = hdf_sd_select(FileID, 0)
;;hdf_sd_getdata, sds_id,adum
;;HDF_SD_END, FileID
;adum (where (adum eq  2550)) =  !VALUES.F_NAN
;adum = adum /100.
;
;
;MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits
;MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2,/USA
;;contour, adum,x,y,levels = levels,c_colors=colors,/cell_fill
;
;snapshot = TVRD()
;TVLCT, r, g, b, /Get
;Device, Z_Buffer=1
;Set_Plot, thisDevice
;image24 = BytArr(3, 800, 500)
;image24[0,*,*] = r[snapshot]
;image24[1,*,*] = g[snapshot]
;image24[2,*,*] = b[snapshot]
;Write_JPEG, 'global_map.jpg', image24, True=1, Quality=100
;stop
;
;for DOY = 1,361,8 do begin
;    LAI = intarr (nc,nr,nyrs)
;    LAI (*,*,*) = 2550
;    DDD = string (DOY,'(i3.3)')
;    for year = yearB, yearE do begin
;	YYYY = string (year, '(i4.4)')
;	filename = IDATA + YYYY + '/' + LABEL + yyyy + DDD + '.hdf'
;	FileID=HDF_SD_Start(filename, /read)
;	sds_id = hdf_sd_select(FileID, 0)
;	hdf_sd_getdata, sds_id,adum
;        adum = adum *10
;	HDF_SD_END, FileID
;	print, year,min(adum), max(adum)
;	LAI (*,*,year - yearB) = adum
;	
;    endfor
;
;indata = intarr (nc,nr)
;indata (*,*) = 2550
;
;for j = 0, nr -1 do begin
;   for i = 0, nc -1 do begin
;       if(min (LAI (i,j,*)) lt 2550) then begin
;	  syears = where (LAI (i,j,*) lt 2550)
;	  indata (i,j) = mean (LAI (i,j,syears))
;         ; if(mean (LAI (i,j,syears)) gt 500.) then  stop
;	 ; print, n_elements (syears), mean (LAI (i,j,syears))
;       endif
;   endfor
;endfor
;
;print, min (indata), max(indata)
;
;ofile = ODATA + LABEL + 'YYYY' + DDD + '.nc4'
;write_glass_output, indata, ofile
;
;endfor
;
;end
;
;; ----------------------------------------------------------------
;
;
;; This program is being deprecated, since it needs "pro proc_glass" 

;pro write_glass_output,indata,ofile
;
;nc = 7200
;nr = 3600
;
;id   = NCDF_CREATE(ofile, /clobber) ;Create netCDF output file
;xid  = NCDF_DIMDEF(id, 'N_lon', nc)                                              ;Define x-dimension
;yid  = NCDF_DIMDEF(id, 'N_lat', nr)                                              ;Define y-dimension
;NCDF_ATTPUT,id, 'CellSize_arcmin' , 3,/global
;NCDF_ATTPUT,id, 'CreatedBy', 'Sarith Mahanama GSFC/NASA',/global
;NCDF_ATTPUT,id, 'Contact', 'Anyone from GMAO Land Group',/global
;str_date=systime()
;NCDF_ATTPUT,id, 'Date', str_date,/global
;vid  = NCDF_VARDEF(id, 'lat',  yid, /DOUBLE)         ;Define latitude variable
;vid  = NCDF_VARDEF(id, 'lon', xid, /DOUBLE)         ;Define longitude variable
;vid  = NCDF_VARDEF(id, 'LAI', [xid, yid], /SHORT)
;NCDF_ATTPUT, id, vid, 'LongName','Leaf Area Index 8-Day 0.05-degrees GEO Grid climatology'
;NCDF_ATTPUT, id, vid, 'units', 'm^2/m^2'
;NCDF_ATTPUT, id, vid, 'scale_factor',0.01
;NCDF_ATTPUT, id, vid, 'valid_range','0 1000'
;NCDF_ATTPUT, id, vid, '_FillValue', 2550
;
;NCDF_CONTROL, id, /ENDEF
;
;dxy = 360.d/7200.d
;
;x = indgen (nc)*dxy -180. + dxy/2.d
;y = indgen (nr)*dxy  -90. + dxy/2.d
;
;NCDF_VARPUT, id,'lat', y
;NCDF_VARPUT, id,'lon', x
;for j =0, nr -1 do begin
;NCDF_VARPUT, id, 'LAI',offset=[0,nr-1 -j],count=[nc,1]   , nint(indata [*,j])
;endfor
;NCDF_CLOSE, id 
;
;end

; -------------------------------------------------------------------
 pro irrig_method, ncat, tile_id

; Rerad 0.25-degree GRIPC and MIRCA data

filename = '../irrig.dat'
limits = [-60.,-180.,90.,-180.]
if file_test ('limits.idl') then restore,'limits.idl'

id   = NCDF_OPEN (filename, /NOWRITE)

NCDF_VARGET, id,'IRRG_IRRIGFRAC_SPR',SPRINKLERV
NCDF_VARGET, id,'IRRG_IRRIGFRAC_DRP',DRIPV
NCDF_VARGET, id,'IRRG_IRRIGFRAC_FRW',FLOODV

NCDF_CLOSE, id

SPRINKLERV(where (SPRINKLERV gt 1.)) = !VALUES.F_NAN
DRIPV     (where (DRIPV      gt 1.)) = !VALUES.F_NAN
FLOODV    (where (FLOODV     gt 1.)) = !VALUES.F_NAN

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

SPRINKLER = REPLICATE (!VALUES.F_NAN,IM, JM)
DRIP      = REPLICATE (!VALUES.F_NAN,IM, JM)
FLOOD     = REPLICATE (!VALUES.F_NAN,IM, JM)

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(tile_id[i,j] gt 0) then begin
	SPRINKLER(i,j) = SPRINKLERV (tile_id[i,j] -1)
	DRIP     (i,j) = DRIPV      (tile_id[i,j] -1)
	FLOOD    (i,j) = FLOODV     (tile_id[i,j] -1)
      endif	
   endfor
endfor

SPRINKLER(where (SPRINKLER eq 0.)) = !VALUES.F_NAN
DRIP     (where (DRIP      eq 0.)) = !VALUES.F_NAN
FLOOD    (where (FLOOD     eq 0.)) = !VALUES.F_NAN

colors = indgen (21) + 140
levels = indgen (21)*0.1/2.

position_row1 = [0.02, 0.70, 0.98, 0.95]
position_row2 = [0.02, 0.40, 0.98, 0.65]
position_row3 = [0.02, 0.10, 0.98, 0.35]

position_col   = [0.20, 0.02, 0.80, 0.05] 

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[850,1100], Z_Buffer=0,SET_FONT='Helvetica Bold', /TT_FONT
load_colors
;set_plot,'PS'
!P.font=0
;Device, FILENAME= plotdir + 'IrrigMethod.ps',/color,/PORTRAIT,xsi=0.9*8.2, ysi=0.9*11.7, xoff=.7, yoff=.5, _extra=_extra,/INCHES
Erase,255
!p.background = 255 
!P.position=0

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits , charsize = .9, title = 'SPRINKLER FRACTION',  /noborder,/isotropic, position = position_row1
contour,sprinkler,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2   

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits , charsize = .9, title = 'DRIP FRACTION',  /noborder,/isotropic, position = position_row2
contour,drip,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2   

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits , charsize = .9, title = 'FLOOD FRACTION',  /noborder,/isotropic, position = position_row3
contour,flood,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2   

colorbar,n_levels = 21, levels = levels, colors = colors, labels = levels, position = position_col
;DEVICE, /CLOSE
;Set_Plot, thisDevice
snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 850, 1100)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_PNG, 'IrrigMethod.png' , image24

end
; ---------------------------------------------------------------------------------------------

 pro plot_lai_minmax, ncat, tile_id

; Rerad 0.25-degree GRIPC and MIRCA data

filename = '../irrig.dat'
limits = [-60.,-180.,90.,-180.]
if file_test ('limits.idl') then restore,'limits.idl'

id   = NCDF_OPEN (filename, /NOWRITE)

NCDF_VARGET, id,'IRRG_LAIMIN',LAI_MNV
NCDF_VARGET, id,'IRRG_LAIMAX',LAI_MXV

NCDF_CLOSE, id

LAI_MNV (where (LAI_MNV gt 100.)) = !VALUES.F_NAN
LAI_MXV (where (LAI_MXV gt 100.)) = !VALUES.F_NAN

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

LAI_MN = REPLICATE (!VALUES.F_NAN,IM, JM)
LAI_MX = REPLICATE (!VALUES.F_NAN,IM, JM)
for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(tile_id[i,j] gt 0) then begin
	LAI_MN(i,j) = LAI_MNV (tile_id[i,j] -1)
	LAI_MX(i,j) = LAI_MXV (tile_id[i,j] -1)
      endif	
   endfor
endfor
LAI_MN (where (LAI_MN eq 0.)) = !VALUES.F_NAN
LAI_MX (where (LAI_MX eq 0.)) = !VALUES.F_NAN

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[850,1100], Z_Buffer=0,SET_FONT='Helvetica Bold', /TT_FONT

load_colors
r_in  = [253,224,255,238,205,193,152,  0,124,  0,  0,  0,  0,  0,  0, 48,110, 85]
g_in  = [253,238,255,238,205,255,251,255,252,255,238,205,139,128,100,128,139,107]
b_in  = [253,224,  0,  0,  0,193,152,127,  0,  0,  0,  0,  0,  0,  0, 20, 61, 47]

n_levels = n_elements (r_in)
levels=[0.,0.25,0.5,0.75,1.,1.25,1.5,10. * indgen(11)*0.05+2.]
red  = intarr (256)
green= intarr (256)
blue = intarr (256)

red  (255) = 255
green(255) = 255
blue (255) = 255

for k = 0, N_levels -1 do begin 
	red  (k+1) = r_in (k)
	green(k+1) = g_in (k)
	blue (k+1) = b_in (k)
endfor

TVLCT,red,green,blue

colors = indgen (N_levels) + 1

position_row1 = [0.02, 0.50, 0.98, 0.95]
position_row2 = [0.02, 0.10, 0.98, 0.45]

position_col   = [0.20, 0.02, 0.80, 0.05] 


Erase,255
!p.background = 255 
!P.position=0

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits , charsize = 1.5, title = 'LAI Minimum',  /noborder,/isotropic, position = position_row1
contour,LAI_MN,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2   

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits , charsize = 1.5, title = 'LAI Maximum',  /noborder,/isotropic, position = position_row2
contour,LAI_MX,x,y,levels = levels, c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2
!P.position= position_col 
alpha=fltarr(n_levels,2)
alpha(*,0)=levels
alpha(*,1)=levels
h=[0,1]

dx = (240.)/(n_levels-1)

clev = levels
clev (*) = 1
contour,alpha,levels,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
           /noerase,yticks=1,ytickname=[' ',' '] ,xrange=[min(levels),max(levels)], $
           xtitle=' ', color=0,xtickv=levels, $
           C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)"
contour,alpha,levels,h,levels=levels,color=0,/overplot,c_label=clev
      for k = 0,n_elements(colors) -1 do xyouts,levels[k],1.1,string(levels[k],format='(f4.2)') ,orientation=90,color=0,charsize =0.8


snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 850, 1100)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_PNG, 'LAI_minmax.png' , image24

end
; ---------------------------------------------------------------------------------------------

 pro irrig_fractions, ncat, tile_id

; Rerad 0.25-degree GRIPC and MIRCA data
filename = '../irrig.dat'
limits = [-60.,-180.,90.,-180.]
if file_test ('limits.idl') then restore,'limits.idl'

id   = NCDF_OPEN (filename, /NOWRITE)

NCDF_VARGET, id,'IRRG_IRRIGFRAC',IRRIGFRAC_V
NCDF_VARGET, id,'IRRG_PADDYFRAC',IRRG_PADDYFRACV
NCDF_VARGET, id,'RAINFEDFRAC',RAINFEDFRACV

NCDF_CLOSE, id

IRRIGFRAC_V  (where (IRRIGFRAC_V   gt 1.)) = !VALUES.F_NAN
IRRG_PADDYFRACV  (where (IRRG_PADDYFRACV   gt 1.)) = !VALUES.F_NAN
RAINFEDFRACV(where (RAINFEDFRACV gt 1.)) = !VALUES.F_NAN

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

x = indgen(im)*dx -180. +  dx/2.
y = indgen(jm)*dy -90.  +  dy/2.

IRRG_IRRIGFRAC   = REPLICATE (!VALUES.F_NAN,IM, JM)
IRRG_PADDYFRAC   = REPLICATE (!VALUES.F_NAN,IM, JM)
RAINFEDFRAC = REPLICATE (!VALUES.F_NAN,IM, JM)

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      if(tile_id[i,j] gt 0) then begin
        IRRG_IRRIGFRAC  (i,j) = IRRIGFRAC_V   (tile_id[i,j] -1)
	IRRG_PADDYFRAC  (i,j) = IRRG_PADDYFRACV   (tile_id[i,j] -1)
	RAINFEDFRAC(i,j) = RAINFEDFRACV (tile_id[i,j] -1)
      endif	
   endfor
endfor

IRRG_IRRIGFRAC  (where (IRRG_IRRIGFRAC   eq 0.)) = !VALUES.F_NAN
IRRG_PADDYFRAC  (where (IRRG_PADDYFRAC   eq 0.)) = !VALUES.F_NAN
RAINFEDFRAC(where (RAINFEDFRAC eq 0.)) = !VALUES.F_NAN

colors = indgen (21) + 140
levels = indgen (21)*0.05/2.

position_row1 = [0.02, 0.70, 0.98, 0.95]
position_row2 = [0.02, 0.40, 0.98, 0.65]
position_row3 = [0.02, 0.10, 0.98, 0.35]

position_col   = [0.20, 0.02, 0.80, 0.05] 

thisDevice = !D.Name
set_plot,'Z'
Device, Set_Resolution=[850,1100], Z_Buffer=0,SET_FONT='Helvetica Bold', /TT_FONT
Erase,255
!p.background = 255 
!P.position=0

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits , charsize = 1.5, title = 'IRRIGATED CROP FRACTION',  /noborder,/isotropic, position = position_row1
contour,irrigfrac,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2   

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits , charsize = 1.5, title = 'PADDY FRACTION',  /noborder,/isotropic, position = position_row2
contour,paddyfrac,x,y,levels = levels, c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2

MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits , charsize = 1.5, title = 'RAINFED FRACTION',  /noborder,/isotropic, position = position_row3
contour,rainfedfrac,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2   


colorbar,n_levels = 21, levels = levels, colors = colors, labels = levels, position = position_col

snapshot = TVRD()
TVLCT, r, g, b, /Get
Device, Z_Buffer=1
Set_Plot, thisDevice
image24 = BytArr(3, 850, 1100)
image24[0,*,*] = r[snapshot]
image24[1,*,*] = g[snapshot]
image24[2,*,*] = b[snapshot]
Write_PNG, 'GIA-Hybrid_IrrigFracs.png' , image24

end

; -------------------------------------------------------------------------------------------

pro plot_crop_times, ncat, tile_id
  

filename = '../irrig.dat'
limits = [-60.,-180.,90.,-180.]
if file_test ('limits.idl') then restore,'limits.idl'

id   = NCDF_OPEN (filename, /NOWRITE)

NCDF_VARGET, id,'IRRG_DOY_PLANT',plantv
NCDF_VARGET, id,'IRRG_DOY_HARVEST',harvestv
NCDF_VARGET, id,'IRRG_CROPIRRIGFRAC',Fracv
NCDF_VARGET, id,'IRRG_TYPE',IRRG_TYPE_V
NCDF_VARGET, id,'CROPCLASSNAME',cropname
NCDF_CLOSE, id

cropname = string (cropname)

im = n_elements(tile_id[*,0])
jm = n_elements(tile_id[0,*])

dx = 360. / im
dy = 180. / jm

xx = indgen(im)*dx -180. +  dx/2.
yy = indgen(jm)*dy -90.  +  dy/2.
x = fltarr (IM,JM)
y = fltarr (IM,JM)
PLANT     = REPLICATE (!VALUES.F_NAN,IM, JM, 2, 26)
HARVEST   = REPLICATE (!VALUES.F_NAN,IM, JM, 2, 26)
FRAC      = REPLICATE (!VALUES.F_NAN,IM, JM, 26)
IRRG_TYPE = REPLICATE (!VALUES.F_NAN,IM, JM, 26)

for j = 0l, jm -1l do begin
   for i = 0l, im -1 do begin
      x (i,j) = xx(i)
      y (i,j) = yy(j)
      if(tile_id[i,j] gt 0) then begin
	PLANT    (i,j,*,*) = PLANTV     (tile_id[i,j] -1,*,*)
        HARVEST  (i,j,*,*) = HARVESTV   (tile_id[i,j] -1,*,*)
        FRAC     (i,j,*) = FRACV      (tile_id[i,j] -1,*)
        IRRG_TYPE(i,j,*) = IRRG_TYPE_V (tile_id[i,j] -1,*)
      endif	
   endfor
endfor
;for n = 0, 25 do begin
;  data1 = frac    (*,*,n)
;  data2 = plant   (*,*,0,n)
;  data3 = harvest (*,*,0,n)
;   for j = 0, nr - 1 do begin
;      for i = 0, nc -1 do begin
;         if (mask (i,j) gt 0.5) then begin
;            
;            if((Plant  (I,J,0,n) gt 400) and (frac (i,j,n) gt 0.)) then begin
;               print , i,j, n,Plant  (I,J,0,0:3), frac (i,j,0:3)
;            if((data2  (I,J) gt 400) and (data1 (i,j) gt 0.)) then begin
;               print , i,j, n,data2  (I,J), data1 (i,j)
;               
;            endif 
;         endif
;      endfor
;   endfor
;endfor

;fmask = where (frac gt  0.)
;frac    (where (frac    eq  0.)) = !VALUES.F_NAN
;plant   (where (plant   gt 400)) = !VALUES.F_NAN
;harvest (where (harvest gt 400)) = !VALUES.F_NAN
;stop

colors = indgen (21) + 140
levels = indgen (21)*0.05/2.
DOY   = [ 1, 32, 60,  91, 121, 152, 182, 213, 244, 274, 305, 335, 366, 370]
DOYM    = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349, 366, 370]
DOYL   = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349, 366, 370]
ITYP   = [1,2,3,4] 
colors2= [69, 145, 64, 66, 70, 71,73,75,76, 78, 80, 113, 114, 116, 117]   
colors3 = [69, 64,80]
row_dims1 = [0.10, 0.32, 0.54, 0.76]
row_dims2 = [0.27, 0.49, 0.71, 0.96]
;col_dmis1 = [0.02, 0.34, 0.66]
;col_dmis2 = [0.34, 0.66, 0.98]
col_dmis1 = [0.02, 0.25, 0.50, 0.75]
col_dmis2 = [0.25, 0.50, 0.75, 0.98]

;position_col1   = [0.04, 0.02, 0.32, 0.05] 
;position_col2   = [0.38, 0.02, 0.92, 0.05] 
position_col1   = [0.04, 0.02, 0.32, 0.05] 
position_col2   = [0.38, 0.02, 0.68, 0.05] 
position_col3   = [0.70, 0.02, 0.92, 0.05] 

page = 1
Row = 1 
A=findgen(16)*(!PI*2/16.)
usersym,0.1*cos(a),0.1*sin(a),/fill
thisDevice = !D.Name

set_plot,'Z'
Device, Set_Resolution=[850,1100], Z_Buffer=0,SET_FONT='Helvetica Bold', /TT_FONT
load_colors
;set_plot,'PS'
;!P.font=0
;Device, FILENAME= plotdir + 'gia_irrig_params.ps',/color,/PORTRAIT,xsi=0.9*8.2, ysi=0.9*11.7, xoff=.7, yoff=.5, _extra=_extra,/INCHES
Erase,255
!p.background = 255 
!P.position=0

for n = 0, 25 do begin

   if (row eq 1) then begin
      thisDevice = !D.Name
      set_plot,'Z'
      Device, Set_Resolution=[850,1100], Z_Buffer=0,SET_FONT='Helvetica Bold', /TT_FONT
      Erase,255
      !p.background = 255 
      !P.position=0
   endif
   
   data1 = frac    (*,*,n)
   fmask = where (data1 gt 0.)
   this = plant   (*,*,0,n)
   data2 = this (fmask)
   this = harvest (*,*,0,n)
   data3 = this (fmask) 
   this = irrigtype (*,*,n)
   data4 = this (fmask)     
   lons = x (fmask)
   lats = y (fmask)

   data1 (where (data1 le 0)) = !VALUES.F_NAN
   
   for col = 0,3 do begin
      print, col
      if (col eq 0) then begin
         ptitle = ' : frac'
         data_grid = data1
      endif
      
      if (col eq 1) then begin
         ptitle = ' : DOY plant' 
         data_grid = data2
      endif
      
      if (col eq 2) then begin
         ptitle = ' : DOY harvest'
         data_grid = data3
      endif
 
      if (col eq 3) then begin
         ptitle = ' : IRRG_TYPE'
         data_grid = data4
      endif

      plot_position = [col_dmis1(col),row_dims1(row-1), col_dmis2(col),row_dims2(row-1)]
      print, col, row, plot_position
      MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits , charsize = 0.8, title = cropname(n) + ptitle,  /noborder, position = plot_position
      
      if(col eq 0) then begin
         contour,data_grid,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
      endif else if (col eq 3) then begin
	for i = 0l, n_elements (data_grid) -1l do begin
	    if(data_grid(i) gt 0.)  then oplot,[lons(i), lons(i)], [lats(i), lats(i)],psym=8,color= colors3(value_locate (ITYP,data_grid(i)))	
	endfor
      endif else begin
	for i = 0l, n_elements (data_grid) -1l do begin
	    if((data_grid(i) gt 0.) and (data_grid(i) le 366.))  then oplot,[lons(i), lons(i)], [lats(i), lats(i)],psym=8,color= colors2(value_locate (doy,data_grid(i)))	
	endfor

      endelse
      
      MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2   
      
   endfor
   
   row = row + 1
   
   if (row eq 5) then begin
      
      colorbar,n_levels = 21, levels = levels, colors = colors, labels = levels, position = position_col1
      colorbar,n_levels = n_elements (DOY), levels = DOY, colors = colors2, labels = DOY, position = position_col2
      colorbar,n_levels = 4, levels = [1,2,3,4], colors = colors3, labels = ['Sprinkler', 'Drip', 'Flood',''], position = position_col3
;      ERASE
      snapshot = TVRD()
      TVLCT, r, g, b, /Get
      Device, Z_Buffer=1
      Set_Plot, thisDevice
      image24 = BytArr(3, 850, 1100)
      image24[0,*,*] = r[snapshot]
      image24[1,*,*] = g[snapshot]
      image24[2,*,*] = b[snapshot]
      Write_PNG, 'gia_irrig_params_' + string (page, '(i2.2)') + '.png' , image24
      
      row = 1 
      page = page + 1
      
   endif

   this = plant   (*,*,1,n)
   data2 = this (fmask)
	
   if(max (data2) gt 0) then begin
      if (row eq 1) then begin
         thisDevice = !D.Name
         set_plot,'Z'
         Device, Set_Resolution=[850,1100], Z_Buffer=0,SET_FONT='Helvetica Bold', /TT_FONT
         Erase,255
         !p.background = 255 
         !P.position=0
      endif
      
      this = plant   (*,*,1,n)
      data2 = this (fmask)
      this = harvest (*,*,1,n)
      data3 = this (fmask) 
  
      for col = 0,2 do begin
         if (col eq 0) then begin
            ptitle = ' : frac'
            data_grid = data1
         endif
         
         if (col eq 1) then begin
            ptitle = ' : DOY plant' 
            data_grid = data2
         endif
         
         if (col eq 2) then begin
            ptitle = ' : DOY harvest'
            data_grid = data3
         endif
         
         plot_position = [col_dmis1(col),row_dims1(row-1), col_dmis2(col),row_dims2(row-1)]
         print, col, row, plot_position
         MAP_SET,/CYLINDRICAL,/hires,color= 0,/NoErase,limit=limits , charsize = 0.8, title = cropname(n) + ptitle,  /noborder, position = plot_position
         if(col eq 0) then begin
            contour,data_grid,x,y,levels = levels,c_colors=colors,/cell_fill,/overplot
        endif else if (col eq 3) then begin
	     for i = 0l, n_elements (data_grid) -1l do begin
	         if(data_grid(i) gt 0.) then oplot,[lons(i), lons(i)], [lats(i), lats(i)],psym=8,color= colors3(value_locate (ITYP,data_grid(i)))	
	     endfor
         endif else begin
	     for i = 0l, n_elements (data_grid) -1l do begin
	         if((data_grid(i) gt 0.) and (data_grid(i) le 366.)) then oplot,[lons(i), lons(i)], [lats(i), lats(i)],psym=8,color= colors2(value_locate (doy,data_grid(i)))	
	     endfor
         endelse
         
         MAP_CONTINENTS,/COASTS,color=0,MLINETHICK=2   
         
      endfor
      row = row + 1
      
      if (row eq 5) then begin
         
         colorbar,n_levels = 21, levels = levels, colors = colors, labels = levels, position = position_col1
         colorbar,n_levels = n_elements (DOY), levels = DOY, colors = colors2, labels = DOY, position = position_col2
         colorbar,n_levels = 4, levels = [1,2,3,4], colors = colors3, labels = ['Sprinkler', 'Drip', 'Flood',''], position = position_col3
;	 ERASE
         snapshot = TVRD()
         TVLCT, r, g, b, /Get
         Device, Z_Buffer=1
         Set_Plot, thisDevice
         image24 = BytArr(3, 850, 1100)
         image24[0,*,*] = r[snapshot]
         image24[1,*,*] = g[snapshot]
         image24[2,*,*] = b[snapshot]
         Write_PNG, 'gia_irrig_params_' + string (page, '(i2.2)') + '.png' , image24
         
         row = 1 
         page = page + 1

      endif      
   endif   
endfor

;DEVICE, /CLOSE
;Set_Plot, thisDevice
end
; =========================================================================

pro colorbar,n_levels = n_levels, levels = levels, colors = colors, labels = labels,$
             position = position, vertical= vertical, horizontal = horizontal

  IF KEYWORD_SET(vertical) THEN BEGIN
  
     !P.position=[position(2) + 0.02, position(1) ,position(2) + 0.05, position(3)] 
     
     alpha=fltarr(2,n_levels)
     alpha(0,*)=levels
     alpha(1,*)=levels
     h=[-1,1]
     clev = levels
     clev (*) = 1
     k = 0
     
     levelsx = indgen (N_levels)

     contour,alpha,h,levelsx,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
             /noerase,xticks=1, xtickname=[' ',' '] ,yrange=[min(levelsx),max(levelsx)], $
             ytitle=' ', color=0, $
             C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)",ytickformat = "(A1)"
     contour,alpha,h,levelsx,levels=levels,color=0,/overplot,c_label=clev
     for k = 0,n_elements(levels) -1 do xyouts,1.1, levelsx[k],labels(k) ,color=0,charsize =1.2

  endif else begin

;     !P.position=[position(0) +0.1, position(1)-0.06 ,position(2) - 0.1, position(1)-0.04] 
     !P.position=[0.15, 0.01, 0.85, 0.04]  
     IF KEYWORD_SET(position) then !P.position=position ;

     alpha=fltarr(n_levels,2)
     alpha(*,0)=levels
     alpha(*,1)=levels
     h=[0,1]
     clev = levels
     clev (*) = 1
     k = 0
     
     levelsx = indgen (N_levels)

     contour,alpha,levelsx,h,levels=levels,c_colors=colors,/fill,/xstyle,/ystyle, $
             /noerase,yticks=1, ytickname=[' ',' '] ,xrange=[min(levelsx),max(levelsx)], $
             xtitle=' ', color=0, $
             C_charsize=1.0, charsize=0.5 ,xtickformat = "(A1)",ytickformat = "(A1)"
     contour,alpha,levelsx,h,levels=levels,color=0,/overplot,c_label=clev

     if(n_levels eq 4) then begin
 
	 for k = 0,n_elements(levels) -1 do xyouts, levelsx[k],1.1, strtrim(labels(k)), color=0,charsize =0.8,orientation=90    
    endif else begin
     if(max (labels) le 1.) then begin
     for k = 0,n_elements(levels) -1 do xyouts, levelsx[k],1.1,string(labels(k),'(f5.2)') ,color=0,charsize =0.8,orientation=90    
 
       endif else begin
       for k = 0,n_elements(levels) -1 do xyouts, levelsx[k],1.1,string(labels(k),'(i3.3)') ,color=0,charsize =0.8,orientation=90    
      endelse
     endelse
  endelse

     !P.position=0

end
