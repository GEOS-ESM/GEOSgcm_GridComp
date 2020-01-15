#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

PROGRAM create_vegdyn_ndvi

! USAGE : ./create_vegdyn_ndvi BCSDIR GFILE, IMxJM JPLH OUTDIR
! EXAMPLES
! bin/create_vegdyn_ndvi /discover/nobackup/ltakacs/bcs/Ganymed-4_0/SMAP_EASEv2/SMAP_EASEv2_M36/ SMAP_EASEv2_M36_964x406 964x406_DE 0 M2/M36/
! bin/create_vegdyn_ndvi /discover/nobackup/ltakacs/bcs/Ganymed-4_0/SMAP_EASEv2/SMAP_EASEv2_M09/ SMAP_EASEv2_M09_3856x1624 3856x1624_DE 0 M2/M09/
! bin/create_vegdyn_ndvi /discover/nobackup/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_MERRA-2/DC0144xPC0091_DE1440xPE0720/ DC0144xPC0091_DE1440xPE0720-Pfafstetter  144x91_DC 0 M2/DC144/
! bin/create_vegdyn_ndvi /discover/nobackup/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_MERRA-2/DC0288xPC0181_DE1440xPE0720/ DC0288xPC0181_DE1440xPE0720-Pfafstetter 288x181_DC 0 M2/DC288/
! bin/create_vegdyn_ndvi /discover/nobackup/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_MERRA-2/DC0576xPC0361_DE1440xPE0720/ DC0576xPC0361_DE1440xPE0720-Pfafstetter 576x361_DC 0 M2/DC576/
! bin/create_vegdyn_ndvi /discover/nobackup/ltakacs/bcs/Ganymed-4_0/Ganymed-4_0_MERRA-2/DC1152xPC0721_DE1440xPE0720/ DC1152xPC0721_DE1440xPE0720-Pfafstetter 1152x721_DC 0 M2/DC1152/
! bin/create_vegdyn_ndvi /discover/nobackup/projects/gmao/ssd/land/l_data/geos5/bcs/CLSM_params/mkCatchParam_SMAP_L4SM_v001/SMAP_EASEv2_M36/ SMAP_EASEv2_M36_964x406 964x406_DE 0 smapv1/M36/ 
! bin/create_vegdyn_ndvi /discover/nobackup/projects/gmao/ssd/land/l_data/geos5/bcs/CLSM_params/mkCatchParam_SMAP_L4SM_v001/SMAP_EASEv2_M09/ SMAP_EASEv2_M09_3856x1624 3856x1624_DE 0 smapv1/M09/
! bin/create_vegdyn_ndvi /discover/nobackup/projects/gmao/ssd/land/l_data/geos5/bcs/CLSM_params/mkCatchParam_SMAP_L4SM_v001/DC0576xPC0361_DE0360xPE0180/  DC0576xPC0361_DE0360xPE0180-Pfafstetter 576x361_DC 0 smapv1/DC576/
! bin/create_vegdyn_ndvi /discover/nobackup/projects/gmao/ssd/land/l_data/geos5/bcs/CLSM_params/mkCatchParam_SMAP_L4SM_v002/SMAP_EASEv2_M36/ SMAP_EASEv2_M36_964x406 964x406_DE 1 smapv2/M36/  
! bin/create_vegdyn_ndvi /discover/nobackup/projects/gmao/ssd/land/l_data/geos5/bcs/CLSM_params/mkCatchParam_SMAP_L4SM_v002/SMAP_EASEv2_M09/ SMAP_EASEv2_M09_3856x1624 3856x1624_DE 1 smapv2/M09/  

  use rmTinyCatchParaMod

  implicit none
  
  integer            :: NTILES, N, I, k, iargc, JPLH
  integer, parameter :: nx = 8640, ny = 4320

  character*400 :: BCSDIR, OUTDIR, GFILE, IMxJM, arg(5)
  real, dimension (6) :: VGZ2 = (/35.0, 20.0, 17.0, 0.6, 0.5, 0.6/) ! Dorman and Sellers (1989)
  logical :: file_exists
  real, pointer, dimension (:)  :: z2, z0, ityp
  include 'netcdf.inc'	

  I = iargc()
  IMxJM =''
  GFILE =''

  if(iargc() /= 5) then
     print *, "Wrong Number of arguments: ", iargc()
     print *, "Usage : ./create_vegdyn_ndvi BCSDIR GFILE, IMxJM JPLH OUTDIR"
     stop
  endif
  
  do n=1,5
     call getarg(n,arg(n))
  enddo
  
  read(arg(1),'(a)') BCSDIR
  read(arg(2),'(a)') GFILE
  read(arg(3),'(a)') IMxJM
  read(arg(4),*    ) JPLH
  read(arg(5),'(a)') OUTDIR

  ! create dirs/links
  ! -----------------

  call system('mkdir -p data ; cd data/ ; ln -s /discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/ CATCH')
  call system('mkdir -p '//trim(OUTDIR)) 

  open (10,file = trim(BCSDIR)//'/clsm/catchment.def', &
       form= 'formatted', action = 'read', status = 'old')
  read (10,*) NTILES
  close (10, status = 'keep')
    
  print '(a)', trim(BCSDIR)//'/'//trim(GFILE)
   
  allocate (ityp (1:NTILES))
  allocate (z2   (1:NTILES))
 
  inquire(file=trim(BCSDIR)//'/vegdyn_'//trim(IMxJM)//'.dat',exist=file_exists)
  if(file_exists) then
     open (10,file = trim(BCSDIR)//'/vegdyn_'//trim(IMxJM)//'.dat',&
          form= 'unformatted', action = 'read', status = 'old')
     read (10) ityp     
     close (10)
  else
     print '(a)', 'MISSING VEGDYN FILE reading from mosaic_veg_typs_fracs'
     open (10,file = trim(BCSDIR)//'/clsm/mosaic_veg_typs_fracs',&
       form= 'formatted', action = 'read', status = 'old')
     do n =1, ntiles
        read (10, *) k,k,ityp(n)
     end do
  end if

  if(JPLH == 1) then
     call jpl_canoph_this (ntiles, nx,ny,trim(BCSDIR)//'/rst/'//trim(GFILE), z2)
  else
     Z2  = VGZ2(NINT(ITYP)) 
  endif

  open (20,file=trim(OUTDIR)//'vegdyn_'//trim(IMxJM)//'.dat',status='unknown',action='write',form='unformatted', &
           convert='little_endian')

  write (20) ityp
  write (20) Z2
  
  print *,'ITYP : ', minval(ityp), maxval (ityp)
  print *,'Z2   : ', minval(z2  ), maxval (z2  )
  
  call ascat_r0_this (ntiles, nx,ny,trim(BCSDIR)//'/rst/'//trim(GFILE), z0)
  write (20) Z0

  close (20, status = 'keep')

  print *,'Z0   : ', minval(z0  ), maxval (z0  )
  
  call  gimms_clim_ndvi (ntiles, nx,ny,trim(BCSDIR)//'/rst/'//trim(GFILE), trim(OUTDIR), trim(IMxJM))
contains


! -----------------------------------------------------------------------------------
    
    SUBROUTINE ascat_r0_this (ntiles, nc,nr,gfiler, z0)
      
      implicit none

      ! 1) ASCAT roughness 
      ! /discover/nobackup/adarmeno/projects/k14/arlems-roughness.x3600_y1800_t1.nc4
     
      integer, intent (in)               :: ntiles, nc, nr
      real, pointer, dimension (:), intent (inout) :: z0
      character(*), intent (in)          :: gfiler
      integer  , parameter               :: N_lon_ascat = 3600, N_lat_ascat = 1800
      integer                            :: i,j, status, varid, ncid
      integer                            :: tid, cid
      REAL, ALLOCATABLE, dimension (:)   :: count_pix
      REAL, ALLOCATABLE, dimension (:,:) :: z0_grid, data_grid
      INTEGER, ALLOCATABLE, dimension (:,:) :: tile_id
      character*100                      :: fout

      ! READ CLM4.5 source data files and regrid
      ! ----------------------------------------

      status  = NF_OPEN ('data/CATCH/arlems-roughness.x3600_y1800_t1.nc4', NF_NOWRITE, ncid)
            
      allocate (z0_grid   (1 : NC         , 1 : NR))
      allocate (data_grid (1 : N_lon_ascat, 1 : N_lat_ascat)) 

      status  = NF_INQ_VARID (ncid,'roughness',VarID) ; VERIFY_(STATUS)
      status  = NF_GET_VARA_REAL (ncid,VarID, (/1,1,1/),(/N_lon_ascat, N_lat_ascat,1/), data_grid) ; VERIFY_(STATUS)

      call RegridRasterReal(data_grid, z0_grid)

      status = NF_CLOSE(ncid)

      ! Grid to tile
      ! ------------
                
      ! Reading tile-id raster file

      allocate(tile_id(1:nc,1:nr))
      
      open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
           form='unformatted',convert='little_endian')
      
      do j=1,nr
         read(10)tile_id(:,j)
      end do       

      close (10,status='keep')     
      
      allocate (z0        (1:NTILES))
      allocate (count_pix (1:NTILES))
      
      z0        = 0.
      count_pix = 0.

      do j = 1,nr
         do i = 1, nc
            if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then

               ! z0 0. < 0.1
               if((z0_grid(i,j) >= 2.0e-6).and.(z0_grid(i,j) <= 0.1)) then 
                 z0 (tile_id(i,j)) = z0 (tile_id(i,j)) + z0_grid(i,j)
                 count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1. 
               endif

            endif
         end do
      end do
      
      where (count_pix > 0.) z0 = z0/count_pix
      where (z0 == 0.)       z0 = 2.0e-6
      
      deallocate (count_pix)
      deallocate (z0_grid)
      deallocate (tile_id)
      
  END SUBROUTINE ascat_r0_this
     ! ----------------------------------------------------------------------------------------------------------------------------

    SUBROUTINE gimms_clim_ndvi (ntiles, nc,nr,gfiler, ThisDir, IMxJM)
      
      implicit none
      ! Producing :  GIMMS NDVI 15-day climatology from 5 arcmin data
      !  24 values per tile
      
      integer, intent (in)                   :: NTILES, nc, nr 
      character(*), intent (in)              :: gfiler, ThisDir, IMxJM
      integer  , parameter                   :: N_lon_gimms = 4320, N_lat_gimms = 2160
      integer                                :: status, varid, ncid1, ncid2,ncid
      real, dimension (:,:), allocatable     :: ndvi_grid, data_grid 
      integer, dimension (:,:), allocatable     ::data_grid2
      REAL, ALLOCATABLE, dimension (:)       :: ndvi, count_pix
      INTEGER, ALLOCATABLE, dimension (:,:)  :: tile_id
      integer                                :: yr,mn,yr1,mn1, k,t,i,j,l
      integer, parameter :: scale_fac = 10000
      real,    parameter :: val_min = -0.3, val_max = 1.

      ! Grid to tile
      ! ------------
                
      ! Reading tile-id raster file
      
      allocate(tile_id(1:nc,1:nr))
      
      open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
           form='unformatted',convert='little_endian')
      
      do j=1,nr
         read(10)tile_id(:,j)
      end do
      
      close (10,status='keep')     
      
      ! READ GIMMS NDVI source data files and regrid
      ! ----------------------------------------
      
      status  = NF_OPEN ('data/CATCH/ndvi3g_geo_v1_YYYY_0106.nc4', NF_NOWRITE, ncid1) ; VERIFY_(STATUS)
      status  = NF_OPEN ('data/CATCH/ndvi3g_geo_v1_YYYY_0712.nc4', NF_NOWRITE, ncid2) ; VERIFY_(STATUS)
      status  = NF_INQ_VARID (ncid2,'ndvi',VarID) ; VERIFY_(STATUS)
      
      allocate (ndvi_grid   (1:NC,1:NR))
      allocate (data_grid (1 : N_lon_gimms, 1 : N_lat_gimms)) 
      allocate (data_grid2(1 : N_lon_gimms, 1 : N_lat_gimms)) 
      allocate (ndvi (1:NTILES))
      allocate (count_pix (1:NTILES))
      
      ! writing tile-spaced output
      ! --------------------------
      
      open (31,file=trim(ThisDir)//'ndvi_clim_'//trim(IMxJM)//'.data',status='unknown',action='write',form='unformatted', &
           convert='little_endian')
      
      do K=0,13
         yr = (k+11)/12
         mn = mod(k+11,12)+1
         yr1= (k+12)/12
         mn1= mod(k+12,12)+1

         ndvi = 0.
         count_pix = 0.
         t = k
         if (k == 0 ) then 
            t = 12
            ncid = ncid2
            write(31) float((/yr,mn,16,0,0,0,yr1,mn1,1,0,0,0,NTILES,1/))

            status  = NF_GET_VARA_INT (ncid,VarID, (/1,1,t/),(/N_lon_gimms, N_lat_gimms,1/), data_grid2) ; VERIFY_(STATUS)

            do j = 1,  N_lat_gimms
               data_grid (:,j) =   data_grid2 (:,N_lat_gimms - (j-1)) / real(scale_fac)
            end do
            
            call RegridRasterReal(data_grid, ndvi_grid)
            
            do j = 1,nr
               do i = 1, nc
                  if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then        
                     if((ndvi_grid(i,j) >= val_min).and.(ndvi_grid(i,j) <= val_max)) then 
                        ndvi (tile_id(i,j)) = ndvi (tile_id(i,j)) + ndvi_grid(i,j)
                        count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1.                     
                     endif
                  endif
               end do
            end do
            
            where (count_pix > 0.) ndvi = ndvi /count_pix
            write(31) ndvi
            
         elseif (k == 13) then

            t = 1
            ncid = ncid1
            write(31) float((/yr,mn,1,0,0,0,yr,mn,16,0,0,0,NTILES,1/))

            status  = NF_GET_VARA_INT (ncid,VarID, (/1,1,t/),(/N_lon_gimms, N_lat_gimms,1/), data_grid2) ; VERIFY_(STATUS)

            do j = 1, N_lat_gimms
               data_grid (:,j) = data_grid2 (:,N_lat_gimms - (j-1)) / real(scale_fac)
            end do
            
            call RegridRasterReal(data_grid, ndvi_grid)
            
            do j = 1,nr
               do i = 1, nc
                  if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then        
                     if((ndvi_grid(i,j) >= val_min).and.(ndvi_grid(i,j) <= val_max)) then 
                        ndvi (tile_id(i,j)) = ndvi (tile_id(i,j)) + ndvi_grid(i,j)
                        count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1.                     
                     endif
                  endif
               end do
            end do
            
            where (count_pix > 0.) ndvi = ndvi /count_pix
            write(31) ndvi

         else

            do l = 1, 0 , -1
               t = k*2 - l
               if (k <= 6) ncid = ncid1
               if (k >= 7) ncid = ncid2
               if (k >= 7) t = t - 12
               if(l == 1) write(31) float((/yr,mn,1,0,0,0,yr,mn,16,0,0,0,NTILES,1/))
               if(l == 0) write(31) float((/yr,mn,16,0,0,0,yr1,mn1,1,0,0,0,NTILES,1/))

               ndvi = 0.
               count_pix = 0.

               status  = NF_GET_VARA_INT (ncid,VarID, (/1,1,t/),(/N_lon_gimms, N_lat_gimms,1/), data_grid2) ; VERIFY_(STATUS)

               do j = 1,  N_lat_gimms
                  data_grid (:,j) =   data_grid2 (:,N_lat_gimms - (j-1)) / real(scale_fac)
               end do

               call RegridRasterReal(data_grid, ndvi_grid)
               
               do j = 1,nr
                  do i = 1, nc
                     if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then        
                        if((ndvi_grid(i,j) >= val_min).and.(ndvi_grid(i,j) <= val_max)) then 
                           ndvi (tile_id(i,j)) = ndvi (tile_id(i,j)) + ndvi_grid(i,j)
                           count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1.                     
                        endif
                     endif
                  end do
               end do
               
               where (count_pix > 0.) ndvi = ndvi /count_pix
               write(31) ndvi
               
            end do            
         endif
      end do

      close(31,status='keep')    

    END SUBROUTINE gimms_clim_ndvi

  ! ----------------------------------------------------------------------------------------------------------------------------

    SUBROUTINE jpl_canoph_this (ntiles, nc,nr,gfiler, z2)

      implicit none

      ! 1) JPL Canopy Height 
      ! /discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001//Simard_Pinto_3DGlobalVeg_JGR.nc4
     
      integer, intent (in)               :: nc, nr, ntiles
      real, pointer, dimension (:), intent (inout) :: z2
      character(*), intent (in)          :: gfiler
      integer  , parameter               :: N_lon_jpl = 43200, N_lat_jpl = 21600
      integer                            :: i,j, status, varid, ncid
      REAL, ALLOCATABLE, dimension (:)   :: count_pix
      INTEGER, ALLOCATABLE, dimension (:,:) :: data_grid, z2_grid
      INTEGER, ALLOCATABLE, dimension (:,:) :: tile_id
      character*100                      :: fout

      ! READ JPL source data files and regrid
      ! -------------------------------------

      status  = NF_OPEN ('data/CATCH/Simard_Pinto_3DGlobalVeg_JGR.nc4', NF_NOWRITE, ncid)
            
      allocate (z2_grid   (1 : NC         , 1 : NR))
      allocate (data_grid (1 : N_lon_jpl, 1 : N_lat_jpl)) 

      status  = NF_INQ_VARID (ncid,'CanopyHeight',VarID) ; VERIFY_(STATUS)
      status  = NF_GET_VARA_INT (ncid,VarID, (/1,1/),(/N_lon_jpl, N_lat_jpl/), data_grid) ; VERIFY_(STATUS)
      
      call RegridRaster(data_grid, z2_grid)

      status = NF_CLOSE(ncid)

      ! Grid to tile
      ! ------------
                
      ! Reading tile-id raster file

      allocate(tile_id(1:nc,1:nr))
      
      open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
           form='unformatted',convert='little_endian')
      
      do j=1,nr
         read(10)tile_id(:,j)
      end do       

      close (10,status='keep')     
      
      if(.not.associated(z2)) allocate (z2        (1:NTILES))
      allocate (count_pix (1:NTILES))
      
      z2        = 0.
      count_pix = 0.

      do j = 1,nr
         do i = 1, nc
            if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then

               if(z2_grid(i,j) >= 0.) then 
                 z2 (tile_id(i,j)) = z2 (tile_id(i,j)) + real (z2_grid(i,j))
                 count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1. 
               endif

            endif
         end do
      end do
      
      where (count_pix >   0.) z2 = z2/count_pix
      where (z2        < 0.01) z2 = 0.01            ! to ensure Z2 >= MIN_VEG_HEIGHT

      deallocate (count_pix)
      deallocate (z2_grid)
      deallocate (tile_id)
      
    END SUBROUTINE jpl_canoph_this

!--------------------------------------

  END PROGRAM create_vegdyn_ndvi
