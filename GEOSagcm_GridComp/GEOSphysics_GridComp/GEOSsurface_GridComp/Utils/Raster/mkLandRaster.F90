
#include "Raster.h"

  Program MakeLandRaster

    use LogRectRasterizeMod
    use MAPL_HashMod
    use process_hres_data
    use MAPL_SortMod
    use rmTinyCatchParaMod, ONLY: SRTM_maxcat

! Program to create a surface raster file at 2.5' that has
! the ocean divided with a regular lat-lon DE grid. Its inputs
! are Sarith's formatted 2.5' raster of the Pfafstetter catchments with 
! ocean points set to 0. Output is an unformatted integer*4 2.5" raster
! with the ocean points set to negative numbers.

    implicit none


    integer, parameter     :: InUNIT   = 20
    integer, parameter     :: nx0 = 8640
    integer, parameter     :: ny0 = 4320

    integer                :: command_argument_count
    integer                :: i,j,k,n, status,ncid, ncid2
    integer                :: ip, nxt
    integer                :: type, maxtiles, nx, ny
    integer                :: count0,count1,count_rate

    REAL_                  :: dx, dy, d2r   ! Grid spacing of raster grid
    REAL_                  :: xmin, ymin, xmax, ymax, xs, ys, da

    REAL_,     allocatable :: cc(:), ss(:)
    REAL_ ,    allocatable :: rTable(:,:)

    integer,       pointer :: Raster(:,:)
    integer,   allocatable, target :: Raster0(:,:)
    integer,   allocatable :: iTable(:,:)

    integer :: Pix, IceID, LakeID
    integer :: CreateHash, IncrementHash, Hash

    logical                :: DoZip
    logical                :: Verb
    logical                :: regrid, reynolds_sst = .false.

    REAL_                  :: VV(4)
    REAL_                  :: PI=RASTER_PI
    
    ! ESA/SRTM ocean/land/ice/lake mask parameters
    ! --------------------------------------------

    integer, parameter :: nc_esa = 129600, nr_esa = 64800
    integer, allocatable, target, dimension (:,:) :: geos_msk, geos_msk2 
    REAL,    allocatable, DIMENSION (:) :: loc_val
    INTEGER, ALLOCATABLE, DIMENSION (:) :: density, loc_int
    REAL,    allocatable, DIMENSION (:) :: loc_val2
    INTEGER, ALLOCATABLE, DIMENSION (:) :: density2, loc_int2
    logical, dimension (:), allocatable :: unq_mask, unq_mask2      
    integer, pointer    , dimension (:,:) :: subset, subset2
    integer :: dx_esa, dy_esa, NBINS, NPLUS,NPLUS2, NBINS2, NonReynold
    
    integer*8, allocatable, dimension (:) ::  SRTM_catid

    character*1            :: Opt
    character*128          :: arg
    character*128          :: TilFile, RstFile
    character*128          :: Tildir,  Rstdir
    character*128          :: GridName
    character*128          :: InputFile
    character*128          :: MaskFile
    character*128          :: &
    Usage = "mkLandRaster -x nx -y ny -v -h -z -t maxtiles -l LandFile -g GridName"
    include 'netcdf.inc'
    call execute_command_line('cd data/ ; ln -s /discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/ CATCH')
    call execute_command_line('cd ..')

! Process Arguments
!------------------

    nx        = nx0
    ny        = ny0
    GridName  = "Pfafstetter"
    DoZip     = .false.  ! No zipping
    Verb      = .false.  ! Run quiet
    tildir    = 'til/'   ! Write in current dir
    rstdir    = 'rst/'   ! Write in current dir
    maxtiles  = 50000
    InputFile = &
          "data/CATCH/global.cat_id.catch.DL"
     
    I = command_argument_count()

    if(I > 13) then
       print *, "Wrong Number of arguments: ", i
       print *, trim(Usage)
       call exit(1)
    end if

    nxt = 1

    do while(nxt<=I)
       call get_command_argument(nxt,arg)
       if(arg(1:1)/='-') then
          print *, trim(Usage)
          call exit(1)
       end if
       opt=arg(2:2)
       if(len(trim(arg))==2) then
          if(scan(opt,'zvh')==0) then
             nxt = nxt + 1
             call get_command_argument(nxt,arg)
          end if
       else
          arg = arg(3:)
       end if
       select case (opt)
       case ('x')
          read(arg,'(i6)') nx
       case ('y')
          read(arg,'(i6)') ny
       case ('l')
          InputFile = arg
       case ('g')
          Gridname = arg
       case ('v')
          Verb = .true.
       case ('z')
          DoZip = .true.
       case ('h')
          tildir = ''
          rstdir = ''
       case ('t')
          read(arg,*) maxtiles
       case default
          print *, trim(Usage)
          call exit(1)
       end select

       nxt = nxt + 1
    end do
    
!   Check for the 10 arc-sec MaskFile (SM)

    call get_environment_variable ("MASKFILE"        ,MaskFile        )

    print *,'Using Mask file : ',trim(MaskFile)

    call system_clock(count0,count_rate)

    if(DoZip) then
       TilFile = trim(tildir)//trim(Gridname)//'.til.gz'
       RstFile = trim(rstdir)//trim(Gridname)//'.rst.gz'
    else
       TilFile = trim(tildir)//trim(Gridname)//'.til'
       RstFile = trim(rstdir)//trim(Gridname)//'.rst'
    end if

    allocate(cc(nx),       stat=STATUS); VERIFY_(STATUS)
    allocate(ss(nx),       stat=STATUS); VERIFY_(STATUS)

    allocate(iTable(0:3,maxtiles),stat=status)
    VERIFY_(STATUS)
    allocate(rTable(1:4,maxtiles),stat=status)
    VERIFY_(STATUS)

! Compute ys and xs at centers of raster cell

    xmin = -180.0_8
    xmax =  180.0_8
    ymin =  -90.0_8
    ymax =   90.0_8

    dx    = (xmax-xmin)/nx
    dy    = (ymax-ymin)/ny
    d2r   = (2._8*PI)/360.0_8
 
    do i = 1,nx
       xs = (xmin + dx*(i-0.5_8))*d2r
       cc(i) = cos(xs)
       ss(i) = sin(xs)
    enddo

    InputFile = 'data/CATCH/'//trim(MaskFile)
  
    if (index(trim(MaskFile),'GEOS5_10arcsec_mask')/=0) then
       ! 10 arcsec new mask
       LakeID = 190000000
       IceID  = 200000000
       dx_esa = nc_esa / nx
       dy_esa = nr_esa / ny
       regrid = .false.
       regrid = nx/=nc_esa .or. ny/=nr_esa
       allocate(SRTM_catid  (1:SRTM_maxcat))
       allocate(geos_msk    (1:nc_esa,1:dy_esa))
       allocate (raster (1:nx, 1:ny)) 

       InputFile = 'data/CATCH/'//trim(MaskFile)

       status    = NF_OPEN (InputFile, NF_NOWRITE, ncid)

       if(status /=0) then
          PRINT *, NF_STRERROR(STATUS)
          print *, 'Problem with NF_OPEN',InputFile
       endif

       status    = NF_GET_VARA_INT64 (ncid,3,(/1/),(/SRTM_maxcat/),SRTM_catid   )  

       if(index(trim(MaskFile),'GEOS5_10arcsec_mask_freshwater-lakes')/=0) then

          ! Special case for Reynolds SST - i.e. Great lakes and the Caspian 
          ! sea are treated as freshwater bodies. We make sure both land-tiles don't 
          ! when the 30 arc-sec mask (raster file) is created from the 10 arc-sec 
          ! GEOS5_10arcsec_mask_freshwater-lakes.nc 
 
          print *, 'Using Reynolds SSTs MASKFILE',trim(MaskFile)
          reynolds_sst = .true.
  
          InputFile = 'data/CATCH/GEOS5_10arcsec_mask.nc'
          status    = NF_OPEN (InputFile, NF_NOWRITE, ncid2)
          allocate(geos_msk2    (1:nc_esa,1:dy_esa))
       endif

       if(Verb) then
          call system_clock(count1)
          print *, 'Opened land file. Time = ', (count1-count0)/float(count_rate)
          count0 = count1
       end if

       do j=1,ny
          status  = NF_GET_VARA_INT  (ncid,4,(/1,(j-1)*dy_esa +1/),(/nc_esa,dy_esa/),geos_msk)  ! Read 10-arcsec rows that lie within the raster row 'j'  

          if(reynolds_sst) &
          status  = NF_GET_VARA_INT  (ncid2,4,(/1,(j-1)*dy_esa +1/),(/nc_esa,dy_esa/),geos_msk2)  ! Read 10-arcsec MERRA2/OSTIA mask file

          if (regrid) then
             if(status /=0) then
                PRINT *, NF_STRERROR(STATUS)
                print *, 'Problem with NF_GET_VARA_INT',InputFile,status
             endif
             
             do i = 1,nx
                if (associated (subset)) NULLIFY (subset)
                subset => geos_msk ((i-1)*dx_esa + 1 : i*dx_esa, 1:dy_esa)          ! rectangular array contains 10-arcsec pixels that lie within the raster grid cell at i,j
                if(maxval (subset) > SRTM_maxcat) then
                   where (subset == LakeID)subset = SRTM_maxcat + 1
                   where (subset == IceID) subset = SRTM_maxcat + 2
                endif
                
                if (maxval (subset) > 0) then  ! check whether there are Non-ocean 10-arcsec pixels 
                   
                   NPLUS = count(subset >= 1 .and. subset <= SRTM_maxcat + 2)                  ! Count non-ocean 10-arcsec pixels within  
                   allocate (loc_int (1:NPLUS))
                   allocate (unq_mask(1:NPLUS))
                   loc_int = pack(subset,mask = (subset >= 1 .and. subset <= SRTM_maxcat + 2)) ! loc_int contains catch_indices of non-ocean 10-arcsec pixels 
                   call MAPL_Sort (loc_int)
                   unq_mask = .true.
                   do n = 2,NPLUS 
                      unq_mask(n) = .not.(loc_int(n) == loc_int(n-1))                          ! count number of unique numbers in loc_int for binning 
                   end do
                   NBINS = count(unq_mask)

                   if(reynolds_sst) then
                      if (associated (subset2)) NULLIFY (subset2)
                      subset2 => geos_msk2 ((i-1)*dx_esa + 1 : i*dx_esa, 1:dy_esa)  
                      if(maxval (subset2) > SRTM_maxcat) then
                         where (subset2 == LakeID)subset2 = SRTM_maxcat + 1
                         where (subset2 == IceID) subset2 = SRTM_maxcat + 2
                      endif
                      NPLUS2 = count(subset2 >= 1 .and. subset2 <= SRTM_maxcat + 2)                  ! Count non-ocean 10-arcsec pixels within  
                      allocate (loc_int2 (1:NPLUS2))
                      allocate (unq_mask2(1:NPLUS2))
                      loc_int2 = pack(subset2,mask = (subset2 >= 1 .and. subset2 <= SRTM_maxcat + 2)) ! loc_int contains catch_indices of non-ocean 10-arcsec pixels 
                      call MAPL_Sort (loc_int2)
                      unq_mask2 = .true.
                      do n = 2,NPLUS2 
                         unq_mask2(n) = .not.(loc_int2(n) == loc_int2(n-1))                          ! count number of unique numbers in loc_int for binning 
                      end do
                      NBINS2 = count(unq_mask2)
                   endif

                   if (NBINS > 1) then
                      allocate(loc_val (1:NBINS))
                      allocate(density (1:NBINS))
                      loc_val = 1.*pack(loc_int,mask =unq_mask)                              ! loc_val contains available non-ocean catch_indices within the i,j grid cell,
                      ! Those numbers will be used as bin values
                      call histogram (dx_esa*dy_esa, NBINS, density, loc_val, real(subset))  ! density is the pixel count for each bin value 
                      raster (i,j) = loc_val (maxloc(density,1))                             ! picks maximum density as the dominant catchment_index at i,j

                      if(reynolds_sst) then
                         allocate(loc_val2 (1:NBINS2))
                         allocate(density2 (1:NBINS2))
                         loc_val2 = 1.*pack(loc_int2,mask =unq_mask2) 
                         call histogram (dx_esa*dy_esa, NBINS2, density2, loc_val2, real(subset2))
                         NonReynold = loc_val2 (maxloc(density2,1))
                         if(NonReynold /= raster (i,j)) then
                            if(NonReynold >= 1 .and. NonReynold <= SRTM_maxcat) raster (i,j) = NonReynold ! reset to NonReynold land
                         endif
                      endif

                      deallocate (loc_val, density)
                      if(reynolds_sst) deallocate (loc_val2, density2)
                   else
                      raster (i,j) = loc_int (1)
                   endif
                   deallocate (loc_int, unq_mask)
                   if(reynolds_sst) deallocate (loc_int2, unq_mask2)
                   
                   if(raster (i,j) == SRTM_maxcat + 1) raster (i,j) = LakeID
                   if(raster (i,j) == SRTM_maxcat + 2) raster (i,j) = IceID
                endif
             end do
          else
             raster (:,j) = geos_msk (:,1)
          endif
       end do

       status    = NF_CLOSE (ncid)
       if(reynolds_sst) status    = NF_CLOSE (ncid2)
    else
       LakeID = 6190000
       IceID  = 6200000
       ! Input file: (formatted for now)

       allocate(raster0(nx0,ny0),stat=STATUS); VERIFY_(STATUS)

       regrid = nx/=nx0 .or. ny/=ny0
       
       if(regrid) then
          allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
       else
          raster => raster0
       end if

       open (InUnit,file=InputFile, form='formatted',status='old')

       if(Verb) then
          call system_clock(count1)
          print *, 'Opened land file. Time = ', (count1-count0)/float(count_rate)
          count0 = count1
       end if

       do j=1,ny0
          read (InUnit,*) raster0(:,j)
       enddo

!    call ReadRaster(InputFile,raster) ! reead unformatted

       if(regrid) then
          call RegridRaster(Raster0,raster)
       endif
       
    endif

    if(Verb) then
       call system_clock(count1)
       print *, 'Read land file. Time = ', (count1-count0)/float(count_rate)
       count0 = count1
    end if

! Loop over ocean latitudes

     if(Verb) then
        write (6, '(A)', advance='NO') 'Processing catchment data'
     end if

    ip = 0
    Hash  = CreateHash(4*1024)

    LATITUDES: do j=1,ny

       where(raster(:,j)>=LakeID)
          raster(:,j) = (raster(:,j)/10000)*10000
       end where

       ys  = ymin + dy*(j-0.5_8)
       da  = (sin(d2r*(ys+0.5*dy)) -   &
              sin(d2r*(ys-0.5*dy)) )*(dx*d2r)

       vv(3) = da*ys
       vv(4) = da

       LONGITUDES: do i=1,nx

          vv(1) = ss(i)*da
          vv(2) = cc(i)*da

          Pix = raster(i,j)

          k = MAPL_HashIncrement(Hash,Pix);

          raster(i,j) = k

          if(k<=ip) then    ! Bump the counter and the lon, lat, and area sums
             iTable(1,k) = iTable(1,k) + 1
             rTable(1,k) = rTable(1,k) + vv(1)
             rTable(2,k) = rTable(2,k) + vv(2)
             rTable(3,k) = rTable(3,k) + vv(3)
             rTable(4,k) = rTable(4,k) + vv(4)
          else              ! We have a new tile in the exchange grid
             ip = ip + 1
             if(ip>maxtiles) then
                print *, "Exceeded maxtiles = ", maxtiles, &
                         ". Use -t option to increase."
                call exit(1)
             end if
             
             select case (Pix)
             case (0)
                type = 0
             case(190000000,6190000)
                type = 19
             case(200000000,6200000)
                type = 20
             case default
                type = 100
             end select

             iTable(0 ,ip) = type
             iTable(1 ,ip) = 1
             rTable(:4,ip) = vv

             iTable( 2,ip) = Pix
             iTable( 3,ip) = 1
          end if

     end do LONGITUDES

     if(Verb) then
        if(mod(j,200)==0) then
           write (6, '(A)', advance='NO') '.'
        end if
     end if

  end do LATITUDES

  close(InUnit)
  
  if(Verb) then
     call system_clock(count1)
     print *
     print *, "Found ",ip," unique tiles."
     print *, "Tiling Time = ",(count1-count0)/float(count_rate)
     count0 = count1
     type = sum(iTable(1,:ip))
     print *, type, nx*ny, nx*ny-type
  end if

  call DestroyHash(Hash)

! Compute proper longitude and latitude in degrees and compress
! the real table for WriteTiling.

  if(Verb) print *, "Computing weighted lons and lats..."
  
  do k=1,ip
     rTable(1,k) = atan2(rTable(1,k),rTable(2,k))/d2r
     rTable(2,k) =       rTable(3,k)/rTable(4,k)     
     rTable(3,k) = rTable(4,k)
  end do

! Sort table by the first grid type in descending order.

  if(Verb) print *, "Sorting..."
  call SortTiling(Raster,rTable(:,:ip),iTable(:,:ip))

  if(Verb) then
     call system_clock(count1)
     print *,  "Done Sorting. Time = ", (count1-count0)/float(count_rate)
     count0 = count1
  end if

  do k=1,ip
     iTable( 3,k) = 1
  end do

  if(Verb) print *, 'Writing til file...'
  call WriteTiling(TilFile, (/Gridname/), (/ip/), (/1/), (/ip/),      &
                   nx, ny, iTable(:,:ip), rTable(:4,:ip), Dozip, Verb )

  if(Verb) print *, 'Writing raster file...'
  call WriteRaster( RstFile, Raster, DoZip)
  
  if(Verb) then
     call system_clock(count1)
     print *,  "Done Writing. Time = ", (count1-count0)/float(count_rate)
     count0 = count1
  end if

! All done

  if(Verb) print * , 'Terminated Normally'
  call exit(0)

contains

subroutine RegridRaster(Rin,Rout)

  integer, intent(IN)  :: Rin(:,:)
  integer, intent(OUT) :: Rout(:,:)

  REAL_  :: xx, yy
  integer :: i,j,ii,jj

  xx = size(Rin ,1)/float(size(Rout,1))
  yy = size(Rin ,2)/float(size(Rout,2))

  do j=1,size(Rout,2)
     jj = (j-1)*yy + 1
     do i=1,size(Rout,1)
        ii = (i-1)*xx + 1
        Rout(i,j) = Rin(ii,jj)
     end do
  end do

end subroutine RegridRaster


end Program MakeLandRaster

