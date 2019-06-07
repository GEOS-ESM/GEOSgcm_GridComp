#include "Raster.h"
PROGRAM mkLISTilesPara
 
  use process_hres_data
  use MAPL_ConstantsMod
  use process_hres_data
  use rmTinyCatchParaMod
  use MAPL_SortMod

  implicit none
  integer  , parameter :: nc_esa = 129600, nr_esa = 64800, SRTM_maxcat = 291284
  real     , parameter :: pi= RASTER_PI
  integer  , parameter :: HYDRO1k_maxcat = 6000000
  integer, parameter   :: nc_gswp2 = 360, nr_gswp2 = 180, n_gswp2 =15238 
  integer, parameter :: max_pfaf_smap = 100
  character(40) :: arg
  integer       ::  i, N_args, iargc, status
  character*300 :: latlon_vector_file
  integer       :: nc, nr
  character*200 :: gfile
  character*200 :: tmpstring, tmpstring1, tmpstring2
  character*128 :: MaskFile
  character*100 ::  char_string
  real :: dx, dy
  integer :: ncells, dateline, nc_domain,nr_domain,i_offset,j_offset
  
  N_args = iargc()

  if(N_args /= 2) then
     print *,'USAGE : bin/mkLISTilesPara -vfile filename'
     stop
  end if

  i=0      
  
  do while ( i < N_args )
     
     i = i+1
     
     call getarg(i,arg)
     
     if     ( trim(arg) == '-vfile' ) then
        i = i+1
        call getarg(i,latlon_vector_file)
	
     else ! stop for any other arguments

        print *,'USAGE : bin/mkLISTilesPara -vfile filename'
        stop        

     endif     
end do

call system('mkdir -p data/ ; mkdir -p til/ ; mkdir -p rst/ ; mkdir -p clsm/plots')
call system('cd data/ ; ln -s /discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/ CATCH')  
call system('cd ..')
     
!   Check for the 10 arc-sec MaskFile
! -----------------------------------

call getenv ("MASKFILE"        ,MaskFile        )

if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then       
   ! Use new ESA based MaskFile
   print *, 'Using MaskFile ', trim(MaskFile)
   call create_files_esa (nc, nr, gfile,latlon_vector_file)
 
else
   ! Use old mask file (Ganymed-4_0 and before
   call create_files (nc,nr,gfile,latlon_vector_file)
endif

 open (10,file= trim(latlon_vector_file), form = 'formatted',action= 'read', status = 'old')
  
 read (10,'(a)') char_string
  
 read (char_string,*,IOSTAT=status) ncells, dx,dy, dateline,nc_domain,nr_domain,i_offset,j_offset
 
 close (10, status='keep')

 ! create Grid2Catch transfer file
 ! -------------------------------
 
 CALL CREATE_ROUT_PARA_FILE (NC, NR, trim(gfile), deltaXY=dx)  

tmpstring1 = '-e EASE -g '//trim(gfile) 
write(tmpstring2,'(2(a2,x,i5,x))')'-x',nc,'-y',nr
tmpstring = 'bin/mkCatchParam_openmp '//trim(tmpstring2)//' '//trim(tmpstring1)
print *,trim(tmpstring)

call system(tmpstring)

contains



!
! --------------------------------------------------------------------------------------------
!

SUBROUTINE create_files (nc,nr,gfile,filename)

  implicit none

  real, allocatable, dimension (:,:) :: &
     cti_mean, cti_std, cti_min, cti_max, cti_skew
  character (*), intent (in) :: filename
  character*200, intent (out) :: gfile
  integer, intent(out) :: nc,nr
  integer :: ncells,status,i,j,k,ix1,ix2,iy1,iy2,im,jm,ii,jj, dateline, CatCount
  real :: dx, dy, dxh, dyh,d2r,lats
  integer :: n,nc_domain,nr_domain,i_offset,j_offset, nc_global,nr_global
  integer :: catid_index, NBINS, NPLUS, catNo
  character*11 :: glabel1
  character*9  :: glabel2 
  character*100 ::  gtopo30,char_string
  real :: x0, y0
  real :: lat,lon,da, x1,x2,x3,x4,x5, ox1,ox2,ox3,ox4,ox5
  integer, dimension(:,:), allocatable :: tileid
  real*4, dimension (:,:), allocatable ::    q0
  real, dimension (:), allocatable :: tile_ele,tile_area_land
  integer, allocatable :: catid(:,:)
  integer, allocatable :: subset_catid(:,:)
  REAL,    allocatable, DIMENSION (:) :: loc_val
  INTEGER, ALLOCATABLE, DIMENSION (:) :: density, loc_int
  logical, dimension (:), allocatable :: unq_mask      
  logical :: counted

  nc = 8640
  nr = 4320

  allocate(catid       (1:nc,1:nr))
  nc_domain = 0
  nr_domain = 0
  i_offset  = 0
  j_offset  = 0
 
  open (10,file='data/CATCH/global.cat_id.catch.DL', form='formatted', &
       action='read', status='old')!
  
  do j=1,nr
     read(10,*)(catid(i,j),i=1,nc)
  end do

  print *,'CATID MIN/MAX : ',minval(catid),maxval(catid)

  close (10,status='keep') 
 !
 ! opening LIS latlon vector file
 !
  dateline = 1
  open (10,file= trim(filename), form = 'formatted',action= 'read', status = 'old')
  
  read (10,'(a)') char_string

  ! Reading #of grid cells in the domain or the global file, dx, dy, dateline, nc_domain,nr_domain,i_offset,j_offset 
  ! added dateline - (1) SET DATELINE = 1 for GSWP2, NLDAS if  the dateline is along the western edge, (2)
  ! SET DATELINE = 0 for MERRA GEOS5 grids when the dateline lies along the center of the 1st column
  ! nc_domain,nr_domain,i_offset,j_offset are optional, they are only used when working on a rectangular subset of the globe 
  !(for e.g. NLDAS), where nc_domain = # of columns in the domain, nr_domain  = # of rows in the domain, 
  ! i_offset = i index of LL grid cell of the domain minus 1, when the domain grid is overlayed on the similar global grid, 
  !   counting columns from the dateline to the dateline, 
  ! j_offset is the same counting rows from the south pole to the north pole. Note that the global grid has the 
  ! dateline and the south pole as the western and sothern edges 

  read (char_string,*,IOSTAT=status) ncells, dx,dy, dateline, nc_domain,nr_domain,i_offset,j_offset
 
  da = dx*dy*MAPL_RADIUS*MAPL_RADIUS*pi*pi/180./180./1000000. 

  if (dateline == 0) then

     !For MERRA, GEOS5 type
     ! --------------------
     x0 = -180.           ! DC
     y0 = -90. + dy       ! PC
     
  else

     ! For GSWP-2, NLDAS
     ! -----------------
     x0 = -180. + dx/2.
     y0 = -90.  + dy/2.
     
  endif

  dxh= 360._8/nc
  dyh= 180._8/nr
  d2r= PI/180._8

  catCount = 0

  nc_global = nint(360./dx)
  nr_global = nint(180./dy)

  IM = nc/nc_global      ! i-Pixels in GCM box
  JM = nr/nr_global

  if (nc_domain == 0) nc_domain = nc_global
  if (nr_domain == 0) nr_domain = nr_global
  
  write(glabel2,'(i4.4,1a,i4.4)')nc_domain,'x',nr_domain 
  write(glabel1,'(i5.5,1a,i5.5)')nc_global,'x',nr_global

  if (dateline == 0) then
      gfile = 'DC_'//glabel1//'_PC_'//glabel2
   else
      gfile = 'DE_'//glabel1//'_PE_'//glabel2
   endif

! reading in other input data files
! ---------------------------------
  
  allocate (cti_mean(1:nc_gswp2,1:nr_gswp2))
  allocate (cti_std (1:nc_gswp2,1:nr_gswp2))
  allocate (cti_min (1:nc_gswp2,1:nr_gswp2))
  allocate (cti_max (1:nc_gswp2,1:nr_gswp2))
  allocate (cti_skew(1:nc_gswp2,1:nr_gswp2))
  allocate (q0      (1:nc      ,1:nr      ))


  cti_mean = 0.
  cti_std  = 0.
  cti_min  = 0.
  cti_max  = 0.
  cti_skew = 0.
  q0       = 0.


! opening and writing  .til and .rst cti_stat and catchment.def files
! -------------------------------------------------------------------

  open (11,file='/discover/nobackup/rreichle/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1/FV_360x180_DE_360x180_DE.til', &
       form='formatted',action='READ',status='OLD')
  open (12,file='/discover/nobackup/rreichle/l_data/geos5/bcs/SiB2_V2/DE/GSWP2_1by1/cti_stats.dat',                &
       form='formatted',action='READ',status='OLD')

  read (11,*) i
  read (11,*) i
  read (11,'(a)') 
  read (11,*) i
  read (11,*) i
  read (11,'(a)') 
  read (11,*) i
  read (11,*) i
  read (12,*) i
  
  do n =1,n_gswp2

     read(11,'(i10,i9,2f10.4,2i5,f16.12,3i8,f16.12,i8,f13.4)')k,k,lon,lat,i,j 
     read(12,'(i8,i8,5(1x,f8.4))')k,k, x1,x2,x3,x4,x5
     cti_mean(i,j) = x1
     cti_std (i,j) = x2
     cti_min (i,j) = x3
     cti_max (i,j) = x4
     cti_skew(i,j) = x5

  end do

  close(11,status='keep')
  close(12,status='keep')

  gtopo30 = 'data/CATCH/srtm30_withKMS_2.5x2.5min.data'
  open (13,file=trim(gtopo30),form='unformatted',status='old',convert='little_endian')
  read (13) q0
  close(13,status='keep')
  
  open  (20,file ='til/'//trim(gfile)//'.til-first',form = 'formatted'  ,action= 'write', status = 'unknown') 
  open  (21,file ='rst/'//trim(gfile)//'.rst',form = 'unformatted',action= 'write', status = 'unknown')
  open  (22,file ='clsm/cti_stats.dat-first',       form = 'formatted'  ,action= 'write', status = 'unknown')
  open  (23,file ='clsm/catchment.def',       form = 'formatted'  ,action= 'write', status = 'unknown')

  allocate (tileid  (1:nc      ,1:nr      ))

  tileid   = 0

  write (20,*) ncells
  write (20,*) 2
  write (20,'(a11)') glabel1
  write (20,*) nc_global
 if (dateline == 0) then
     write (20,*) nr_global + 1
  else
     write (20,*) nr_global
  endif
  write (20,'(a9)') glabel2
  write (20,*) nc_domain
  if (dateline == 0) then
     write (20,*) nr_domain + 1
  else
     write (20,*) nr_domain
  endif

  write(22,*) ncells

  allocate(subset_catid (1:IM, 1:JM))

  do CatNo = 1, ncells 

     read (10,*) k,lat,lon

     i = nint((lon - x0)/dx) + 1
     j = nint((lat - y0)/dy) + 1

     ix1 = NINT ((lon -dx/2. + 180.) / dxh) + 1 
     ix2 = NINT ((lon +dx/2. + 180.) / dxh) 

     if(ix1 < 0) ix1 = 1 ! a vertical strip of 0.5 dx that lies on the Eastern hemesphere in DC grid is discarded

     iy1 = NINT ((lat - dy/2. +90.) / dyh) + 1
     iy2 = NINT ((lat + dy/2. +90.) / dyh) 

     counted = .false.
     subset_catid = -9999
     do jj = iy1,iy2
        do ii = ix1,ix2 
           if ((catid (ii,jj)  >= 10000) .and. (catid (ii,jj) <= HYDRO1k_maxcat)) then 
              if (.not. counted) catCount = catCount + 1
              subset_catid (ii - ix1 +1, jj -iy1 +1) = catid (ii,jj)
              tileid (ii,jj) = catCount
              counted = .true.
           endif
        end do
     end do
    
     IF(COUNTED) THEN

        NPLUS = count(subset_catid >= 1 .and. subset_catid <= HYDRO1k_maxcat)
        allocate (loc_int (1:NPLUS))
        allocate (unq_mask(1:NPLUS))
        loc_int = pack(subset_catid,mask = (subset_catid >= 1 .and. subset_catid <= HYDRO1k_maxcat)) ! loc_int contains catch_indices of non-ocean ESA pixels 
        call MAPL_Sort (loc_int)
        unq_mask = .true.
        do n = 2,NPLUS 
           unq_mask(n) = .not.(loc_int(n) == loc_int(n-1)) ! count number of unique numbers in loc_int for binning
        end do
        NBINS = count(unq_mask)
           
        if (NBINS >= 1) then
           allocate(loc_val (1:NBINS))
           allocate(density (1:NBINS))
           loc_val = 1.*pack(loc_int,mask =unq_mask)                               ! loc_val contains available non-ocean catch_indices within the i,j grid cell,
           ! Those numbers will be used as bin values
           call histogram (im*jm, NBINS, density, loc_val, real(subset_catid))   ! density is the pixel count for each bin value
           catid_index = loc_val (maxloc(density,1))  
           deallocate (loc_val, density)
        else
           print *,'Check suspicous'
           catid_index = loc_int (1)          
        endif

        deallocate (loc_int, unq_mask)
        x1 =cti_mean(ceiling(real(i)/(real(nc_global)/360.)),ceiling(real(j)/(real(nr_global)/180.))) 
        x2 =cti_std (ceiling(real(i)/(real(nc_global)/360.)),ceiling(real(j)/(real(nr_global)/180.)))
        x3 =cti_min (ceiling(real(i)/(real(nc_global)/360.)),ceiling(real(j)/(real(nr_global)/180.)))
        x4 =cti_max (ceiling(real(i)/(real(nc_global)/360.)),ceiling(real(j)/(real(nr_global)/180.)))
        x5 =cti_skew(ceiling(real(i)/(real(nc_global)/360.)),ceiling(real(j)/(real(nr_global)/180.)))
        
        if((CatCount == 1).and.(x1==0.)) then
           x1 = 12.5830
           x2 = 2.6039 
           x3 = 7.8878 
           x4 = 23.1378 
           x5 = 0.8388    
        endif
        
        if(x1.gt.0.) then
           ox1 = x1 
           ox2 = x2 
           ox3 = x3 
           ox4 = x4 
           ox5 = x5 
        endif

        if (nr_global < 1000) then
           write(20,'(i10,i9,2f10.4,2i5,f19.12,i10,e13.4,i15,i8)') 100,catid_index,lon,lat,i-i_offset, &
                j-j_offset,real(count(tileid (ix1:ix2,iy1:iy2) > 0))/real((ix2-ix1+1)*(iy2-iy1+1)), &
                CatNo,da*cos(lat*pi/180.)
           if(x1.gt.0.) then 
              write(22,'(i8,i8,5(1x,f8.4))')catcount,catid_index, x1,x2,x3,x4,x5
           else
              write(22,'(i8,i8,5(1x,f8.4))')catcount,catid_index, ox1,ox2,ox3,ox4,ox5 
           endif
        else
           write(20,'(i10,i9,2f10.4,2i5,f19.12,i10,e13.4,i15,i8)') 100,catid_index,lon,lat,i-i_offset,j-j_offset, &
                real(count(tileid (ix1:ix2,iy1:iy2) > 0))/real((ix2-ix1+1)*(iy2-iy1+1)),CatNo,da*cos(lat*pi/180.) 
           if(x1.gt.0.) then 
              write(22,'(i8,i8,5(1x,f8.4))')catcount,catid_index, x1,x2,x3,x4,x5
           else
              write(22,'(i8,i8,5(1x,f8.4))')catcount,catid_index, ox1,ox2,ox3,ox4,ox5 
           endif
        endif
     else
        ncells = ncells - 1        
     ENDIF
     end do

     allocate (tile_ele(1:ncells))
     allocate (tile_area_land (1:ncells))
     tile_ele = 0.
     tile_area_land = 0.
    
     do j = 1,nr
        write(21) tileid(:,j)
        lats = -90._8 + (j - 0.5_8)*dyh
        
        do i = 1,nc
           if((tileid(i,j) > 0).and.(tileid(i,j) <= ncells))then
              tile_ele(tileid(i,j)) = tile_ele(tileid(i,j)) + q0(i,j)*   &	
                   (sin(d2r*(lats+0.5*dyh)) -sin(d2r*(lats-0.5*dyh)))*(dxh*d2r)
              tile_area_land(tileid(i,j)) = tile_area_land(tileid(i,j)) +  &	
                   (sin(d2r*(lats+0.5*dyh)) -sin(d2r*(lats-0.5*dyh)))*(dxh*d2r)
           endif
        enddo
     end do
     
     tile_ele = tile_ele/tile_area_land
     
     ! adjustment Global Mean Topography to 614.649 (615.662 GTOPO 30) m
    ! ---------------------------
    tile_ele =tile_ele*(614.649D0 /656.7712) 

    close (10,status = 'keep')
    close (20,status = 'keep')
    open  (19,file ='til/'//trim(gfile)//'.til-first',form = 'formatted'  ,action= 'read', status = 'old')   
    open  (20,file ='til/'//trim(gfile)//'.til',form = 'formatted'  ,action= 'write', status = 'unknown') 
    open  (21,file ='clsm/cti_stats.dat-first',       form = 'formatted'  ,action= 'read', status = 'old')
    open  (22,file ='clsm/cti_stats.dat',       form = 'formatted'  ,action= 'write', status = 'unknown')

    read (19,*) i
    read (21,*) i 
    write(20,*) ncells
    write(22,*) ncells
    write(23,*) ncells

    read (19,*) i
    write(20,*) i
    read (19,'(a11)') glabel1
    write (20,'(a11)') glabel1
    read (19,*) nc_global
    write (20,*) nc_global
    read  (19,*) nr_global
    write (20,*) nr_global
    read  (19,'(a9)') glabel2
    write (20,'(a9)') glabel2
    read (19,*) nc_domain
    write (20,*) nc_domain
    read (19,*) nr_domain
    write (20,*) nr_domain

   do n =1,ncells

     read (21,'(i8,i8,5(1x,f8.4))') k,ix1,x1,x2,x3,x4,x5
     write (22,'(i8,i8,5(1x,f8.4))') k,ix1,x1,x2,x3,x4,x5
     read (19,'(i10,i9,2f10.4,2i5,f19.12,i10,e13.4,i15,i8)')k,ix1,lon,lat,i,j, x1, catid_index,x2 

     if (dateline == 0) then
        write(20,'(i10,i9,2f10.4,2i5,f19.12,i10,e13.4,i15,i8)')k,ix1,lon,lat,i,j+1, x1, catid_index,x2
     else
        write(20,'(i10,i9,2f10.4,2i5,f19.12,i10,e13.4,i15,i8)')k,ix1,lon,lat,i,j, x1, catid_index,x2  
     endif
     write(23,'(i8,i8,5(2x,f9.4))')n,ix1,lon-dx/2.,lon+dx/2.,lat-dy/2.,lat+dy/2.,tile_ele(n)
   end do
   close (19,status = 'keep')
   close (20,status = 'keep')
   close (21,status = 'keep')
   close (22,status = 'keep')
   close (23,status = 'keep')
   close (10,status = 'keep')

END SUBROUTINE create_files

!
! --------------------------------------------------------------------------------------------
!

SUBROUTINE create_files_esa (nc, nr, gfile,filename)

  implicit none
  integer, intent (out) :: nc, nr
  ! For regridding
  integer :: LakeType =10, IceType =11, OceanType =0
  integer, allocatable, target, dimension (:,:) :: geos_msk, high_msk 
  REAL,    allocatable, DIMENSION (:) :: loc_val
  INTEGER, ALLOCATABLE, DIMENSION (:) :: density, loc_int
  logical, dimension (:), allocatable :: unq_mask      
  integer :: dx_esa, dy_esa, NBINS, NPLUS
  integer :: catid_index, catNo,catCount
  integer*8, allocatable, dimension (:) ::  SRTM_catid 
  integer*8 :: pfaf_dbl
  real,    dimension (max_pfaf_smap) :: pfaf_area
  integer, dimension (max_pfaf_smap) :: pfaf_index
  character (*), intent (in) :: filename
  character*200, intent (out) :: gfile
  integer :: ncells,status,i,j,k,ix1,ix2,iy1,iy2,im,jm, dateline, ix,jx
  integer :: ncid,ii,jj
  real    :: dx, dy
  real (kind =8) :: dxh, dyh,d2r,lats, dyvh, dxvh
  integer :: n,nc_domain,nr_domain,i_offset,j_offset, nc_global,nr_global,msk2rst
  character*11 :: glabel1
  character*9  :: glabel2 
  character*100 ::  gtopo30,char_string
  real :: x0, y0
  real :: lat,lon,da, x1,x2,x3,x4,x5, ox1,ox2,ox3,ox4,ox5
  integer, dimension(:,:), allocatable :: tileid
  real*4, dimension (:,:), allocatable ::    q0, raster
  real, dimension (:), allocatable     :: tile_ele,tile_area_land
  logical :: regrid, counted

  include 'netcdf.inc'

  nc_domain = 0
  nr_domain = 0
  i_offset  = 0
  j_offset  = 0
  dateline = 1 ! Default DE
 !
 ! opening LIS latlon vector file
 !-------------------------------

  open (10,file= trim(filename), form = 'formatted',action= 'read', status = 'old')
  
  read (10,'(a)') char_string

  ! Reading #of grid cells in the domain or the global file, dx, dy, ,nc_domain,nr_domain,i_offset,j_offset 
  ! added dateline - (1) SET DATELINE = 1 for GSWP2, NLDAS if  the dateline is along the western edge, (2)
  ! SET DATELINE = 0 for MERRA GEOS5 grids when the dateline lies along the center of the 1st column
  ! nc_domain,nr_domain,i_offset,j_offset are optional, they are only used when working on a rectangular subset of the globe 
  !(for e.g. NLDAS), where nc_domain = # of columns in the domain, nr_domain  = # of rows in the domain, 
  ! i_offset = i index of LL grid cell of the domain minus 1, when the domain grid is overlayed on the similar global grid, 
  !   counting columns from the dateline to the dateline, 
  ! j_offset is the same counting rows from the south pole to the north pole. Note that the global grid has the 
  ! dateline and the south pole as the western and sothern edges 

  read (char_string,*,IOSTAT=status) ncells, dx,dy, dateline,nc_domain,nr_domain,i_offset,j_offset
 
  da = MAPL_RADIUS*MAPL_RADIUS/1000000. 

  if (dateline == 0) then
     ! For MERRA, GEOS5 type
     x0 = -180.           ! DC
     y0 = -90. + dy       ! PC
  else    
     ! For GSWP-2, NLDAS
     x0 = -180. + dx/2.
     y0 = -90.  + dy/2.
  endif
  
  ! Output raster file resolution
  ! -----------------------------
  
  nc = 43200
  nr = 21600
  
  nc_global = nint(360./dx)
  nr_global = nint(180./dy)
  
  if (nc_domain == 0) nc_domain = nc_global
  if (nr_domain == 0) nr_domain = nr_global
  
  write(glabel2,'(i4.4,1a,i4.4)')nc_domain,'x',nr_domain 
  write(glabel1,'(i5.5,1a,i5.5)')nc_global,'x',nr_global

  if (dateline == 0) then
      gfile = 'DC_'//glabel1//'_PC_'//glabel2
   else
      gfile = 'DE_'//glabel1//'_PE_'//glabel2
   endif
  !
  ! reading in other input data files
  ! ---------------------------------
  
  allocate (raster  (1:8640, 1:4320 ))
  allocate (q0      (1:nc      ,1:nr))
  allocate (tile_ele(1:ncells))
  allocate (tile_area_land (1:ncells))
  raster   = 0.
  q0       = 0.
  tile_ele = 0.
  tile_area_land = 0.
  
  ! opening and writing  .til and .rst and catchment.def files
  ! ----------------------------------------------------------
  
  gtopo30 = 'data/CATCH/srtm30_withKMS_2.5x2.5min.data'
  open (13,file=trim(gtopo30),form='unformatted',status='old',convert='little_endian')
  read (13) raster
  close(13,status='keep')
  
  call RegridRasterReal(raster,q0)
  deallocate (raster)
  
  open  (20,file ='til/'//trim(gfile)//'.til-first',form = 'formatted'  ,action= 'write', status = 'unknown') 
  open  (21,file ='rst/'//trim(gfile)//'.rst',form = 'unformatted',action= 'write', status = 'unknown')
  open  (23,file ='clsm/catchment.def',       form = 'formatted'  ,action= 'write', status = 'unknown')
  
  allocate (tileid  (1:nc      ,1:nr      ))
  
  tileid   = 0
  
  write (20,*) ncells
  write (20,*) 2
  write (20,'(a11)') glabel1
  write (20,*) nc_global
  if (dateline == 0) then
     write (20,*) nr_global + 1
  else
     write (20,*) nr_global
  endif
  write (20,'(a9)') glabel2
  write (20,*) nc_domain
  if (dateline == 0) then
     write (20,*) nr_domain + 1
  else
     write (20,*) nr_domain
  endif
 
  regrid = .true.
  dx_esa = ceiling(real(nc_esa) / real(nc_global)) ! x-dimension (or # of ESA columns within the raster grid cell)
  dy_esa = ceiling(real(nr_esa) / real(nr_global)) ! y-dimension (or # of ESA rows within the raster grid cell)
  msk2rst = nc_esa/nc
  
  allocate(SRTM_catid  (1:SRTM_maxcat+2))
  allocate(geos_msk    (1:dx_esa,1:dy_esa))
  allocate(high_msk    (1:msk2rst,1:msk2rst))
   
  status    = NF_OPEN ('data/CATCH/GEOS5_10arcsec_mask.nc', NF_NOWRITE, ncid)
  status    = NF_GET_VARA_INT64 (ncid,3,(/1/),(/SRTM_maxcat/),SRTM_catid(1:SRTM_maxcat))  ! Read pfafstetter IDs
  if(status /=0) then
     PRINT *, NF_STRERROR(STATUS)
     print *, 'Problem with NF_OPEN : GEOS5_10arcsec_mask.nc'
  endif
  
  SRTM_catid (SRTM_maxcat + 1) = 190000000
  SRTM_catid (SRTM_maxcat + 2) = 200000000 
  dxvh= 360._8/nc_esa
  dyvh= 180._8/nr_esa
  d2r= PI/180._8

  dxh = 360._8/nc
  dyh = 180._8/nr

  catCount = 0

  do catNo = 1, ncells 
     
     pfaf_index = 0
     pfaf_area  = 0.   

     read (10,*) k,lat,lon
     
     i = nint((lon - x0)/dx) + 1
     j = nint((lat - y0)/dy) + 1
     
     ix1 = NINT ((lon -dx/2. + 180.)/dxh) + 1 
     ix2 = NINT ((lon +dx/2. + 180.)/dxh) 

     if(ix1 < 0) ix1 = 1 ! a vertical strip of 0.5 dx that lies on the Eastern hemesphere in DC grid is discarded

     iy1 = NINT ((lat - dy/2. +90.)/dyh) + 1
     iy2 = NINT((lat + dy/2. +90.)/dyh) 
 
     counted = .false.

     do jj = iy1,iy2
        do ii = ix1,ix2
          status  = NF_GET_VARA_INT (ncid,4,(/(ii-1)*msk2rst+1,(jj-1)*msk2rst +1/),  &
               (/msk2rst,msk2rst/),high_msk) 
          if(count(high_msk >= 1 .and. high_msk <= SRTM_maxcat) > 0) then
	    if (.not. counted) catCount       = catCount + 1
	    tileid (ii,jj) = catCount
            counted = .true.
          endif
        end do
     end do

     if (counted) then
        status  = NF_GET_VARA_INT (ncid,4,(/(ix1-1)*msk2rst+1,(iy1-1)*msk2rst +1/), &
             (/dx_esa,dy_esa/),geos_msk) ! Read 10-arcsec rows that lie within the raster row 'j' 
        if(status /=0) then
           PRINT *, NF_STRERROR(STATUS)
           print *, 'Problem with NF_GET_VARA_INT GEOS5_10arcsec_mask.nc',status
        endif
    
        catid_index  = 0
     
        if(maxval (geos_msk) > SRTM_maxcat) then
           if(maxval (geos_msk) == 200000000) where (geos_msk == 200000000) geos_msk = SRTM_maxcat + 2 
           if(maxval (geos_msk) == 190000000) where (geos_msk == 190000000) geos_msk = SRTM_maxcat + 1
        endif
        
        NPLUS = count(geos_msk >= 1 .and. geos_msk <= SRTM_maxcat) ! Count non-ocean ESA pixels within 
        
        if (NPLUS > 0) then ! check whether there are Non-ocean ESA pixels 
           ! catID of the largest Pfafstetter catchment within the grid cell                
           allocate (loc_int (1:NPLUS))
           allocate (unq_mask(1:NPLUS))
           loc_int = pack(geos_msk,mask = (geos_msk >= 1 .and. geos_msk <= SRTM_maxcat)) ! loc_int contains catch_indices of non-ocean ESA pixels 
           call MAPL_Sort (loc_int)
           unq_mask = .true.
           do n = 2,NPLUS 
              unq_mask(n) = .not.(loc_int(n) == loc_int(n-1)) ! count number of unique numbers in loc_int for binning
           end do
           NBINS = count(unq_mask)
           
           if(NBINS > max_pfaf_smap) then
              print *,' # of Pfaf within exceeded : ', NBINS
              stop
           endif
           
           if (NBINS >= 1) then
              allocate(loc_val (1:NBINS))
              allocate(density (1:NBINS))
              loc_val = 1.*pack(loc_int,mask =unq_mask)                               ! loc_val contains available non-ocean catch_indices within the i,j grid cell,
              ! Those numbers will be used as bin values
              call histogram (dx_esa*dy_esa, NBINS, density, loc_val, real(geos_msk))   ! density is the pixel count for each bin value
              catid_index = loc_val (maxloc(density,1))                         ! picks maximum density as the dominant catchment_index at i,j
              
              do n = 1, nbins
                 
                 pfaf_index(n) = nint (loc_val (n))
                 
                 do jx = iy1,iy2 
                    lats = -90._8 + (jx - 0.5_8)*dyvh
                    do ix = ix1, ix2 
                       if(geos_msk (ix -ix1 + 1,jx - iy1 + 1) == nint (loc_val (n))) then 
                          pfaf_area(n) = pfaf_area(n) + (sin(d2r*(lats+0.5*dyvh)) -sin(d2r*(lats-0.5*dyvh)))*(dxvh*d2r)
                       endif
                    end do
                 end do
              end do
              
              deallocate (loc_val, density)
           else
              print *,'Check suspicous'
              catid_index = loc_int (1)
           endif
           deallocate (loc_int, unq_mask)
        
        endif
       
        write(20,'(i10,i9,2f10.4,2i5,f19.12,i10,e13.4,i15,i8)') 100,catid_index,lon,lat,i-i_offset,j-j_offset,    &
             real(count(tileid (ix1:ix2,iy1:iy2) > 0))/real((ix2-ix1+1)*(iy2-iy1+1)),catid_index,da*sum(pfaf_area (1:nbins)),  &
             SRTM_catid(catid_index), catNo 
     else
        ncells = ncells - 1    
     endif
  end do
  close (20, status = 'keep')
  status    = NF_CLOSE (ncid)  
  deallocate (geos_msk)

! rst file creating raster file
! -----------------------------

 
  do j = 1,nr
     write(21) tileid(:,j)
     lats = -90._8 + (j - 0.5_8)*dyh

     do i = 1,nc
        if((tileid(i,j) > 0).and.(tileid(i,j) <= ncells))then
           tile_ele(tileid(i,j)) = tile_ele(tileid(i,j)) + q0(i,j)*   &	
                  (sin(d2r*(lats+0.5*dyh)) -sin(d2r*(lats-0.5*dyh)))*(dxh*d2r)
           tile_area_land(tileid(i,j)) = tile_area_land(tileid(i,j)) +  &	
                  (sin(d2r*(lats+0.5*dyh)) -sin(d2r*(lats-0.5*dyh)))*(dxh*d2r)
        endif
     enddo
  end do

  tile_ele = tile_ele/(tile_area_land + 1.e-20)
  
  ! adjustment Global Mean Topography to 614.649 (615.662 GTOPO 30) m
  ! ---------------------------
  tile_ele =tile_ele*(614.649D0 /656.7712) 
  open  (19,file ='til/'//trim(gfile)//'.til-first',form = 'formatted'  ,action= 'read', status = 'old')   
  open  (20,file ='til/'//trim(gfile)//'.til',form = 'formatted'  ,action= 'write', status = 'unknown') 
  read (19,*) i
  write(20,*) ncells
  read (19,*) i
  write(20,*) i
  read (19,'(a11)') glabel1
  write (20,'(a11)') glabel1
  read (19,*) nc_global
  write (20,*) nc_global
  read  (19,*) nr_global
  write (20,*) nr_global
  read  (19,'(a9)') glabel2
  write (20,'(a9)') glabel2
  read (19,*) nc_domain
  write (20,*) nc_domain
  read (19,*) nr_domain
  write (20,*) nr_domain
  
  write(23,*) ncells
  
  do n =1,ncells
     read (19,'(i10,i9,2f10.4,2i5,f19.12,i10,e13.4,i15,i8)')k,ix1,lon,lat,i,j, x1, catid_index,x2,pfaf_dbl,catNo 
     if (dateline == 0) then
        write(20,'(i10,i9,2f10.4,2i5,f19.12,i10,e13.4,i15,i8)')k,ix1,lon,lat,i,j+1, x1, catid_index,x2,pfaf_dbl,catNo 
     else
        write(20,'(i10,i9,2f10.4,2i5,f19.12,i10,e13.4,i15,i8)')k,ix1,lon,lat,i,j  , x1, catid_index,x2,pfaf_dbl,catNo 
     endif
     write(23,'(i8,i8,5(2x,f9.4))')n,ix1,lon-dx/2.,lon+dx/2.,lat-dy/2.,lat+dy/2.,tile_ele(n)
  end do
  close (19,status = 'keep')  
  close (20,status = 'keep')
  close (21,status = 'keep')
  close (23,status = 'keep')
  close (10,status = 'keep')
  
END SUBROUTINE create_files_esa

END PROGRAM mkLISTilesPara
