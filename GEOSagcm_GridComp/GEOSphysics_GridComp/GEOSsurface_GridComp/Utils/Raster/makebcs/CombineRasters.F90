#define I_AM_MAIN
#include "MAPL_ErrLog.h"

program mkOverlaySimple

  use LogRectRasterizeMod,     ONLY: WriteRaster, WriteTiling, SortTiling
  use MAPL_SortMod
  use MAPL_HashMod
  use MAPL_ExceptionHandling
  use MAPL_Constants

! Overlay atmosphere, land, and ocean rasters, creating a .idx file.
! The ocean raster should be defined everywhere, or at least, everywhere
! that the land Raster is non zero. The land raster uses the actual
! continental outlines to define land and has zero else where.
! The land and ocean files produce a surface raster by masking the
! ocean with the land. The resulting Raster is the overlaid with the
! atmosphere producing a unique index for each pair of atmosphere and
! surface rstare values. This pair together with the assigned index
! is written in the .idx file. This file is ascii and is sorted by
! surface types. This puts the ocean points first, since they are given 
! negative values.

  implicit none

  integer, parameter     :: RSTUNIT1  = 20
  integer, parameter     :: RSTUNIT2  = 21
  integer, parameter     :: TILUNIT1  = 22
  integer, parameter     :: TILUNIT2  = 23

  real(REAL64),   parameter     :: PI        = MAPL_PI_R8

  integer                :: command_argument_count
  integer                :: nxt, argl, fill
  integer                :: i, j, k, l, ip
  integer                :: STATUS, i1, i2, nvars, rvars
  integer                :: ip1, ip2, nf1, nf2
  integer                :: nx1, nx2, ny1, ny2, nx, ny
  integer                :: maxtiles, hash
  integer                :: count0,count1,count_rate
  integer                :: j_start, j_end   ! Latitude range for subset processing
  integer                :: subset_mode      ! 0=full, 1=fast, 2=full_hash

  integer,   allocatable :: RST1(:,:)
  integer,   allocatable :: RST2(:  )
  integer,   allocatable :: iTable(:,:)

  real(REAL64) ,    allocatable :: Table1(:,:) 
  real(REAL64) ,    allocatable :: Table2(:,:) 
  real(REAL64) ,    allocatable :: rTable(:,:)
  real(REAL64) ,    allocatable :: cc(:), ss(:)
  real(REAL64)                  :: dx, dy, area, xc, yc, d2r, vv(4)
  real(REAL64)                  :: lats, lons, da

  logical                :: DoZip
  logical                :: Verb
  logical                :: found
  logical                :: Merge
  logical                :: s_flag=.false.
  logical                :: process_this_lat
  logical                :: accumulate_this_lat
                         
  character*4            :: tildir, rstdir
  character*14           :: sg
  character*1            :: Opt 
  character*128          :: arg
  character*128          :: Overlay=''
  character*128          :: GridName1, GridName2
  character*128          :: Grid1, Grid2
  character*128          :: TilFile, RstFile
  character*128          :: &
      Usage = "CombineRasters -v -h -z -t MT -g GF -f TYPE -s SG -j J1:J2 -m MODE BOTTOMRASTER TOPRASTER"

  character*128          :: Iam = "CombineRasters"
  integer :: Pix1, Pix2

! Argument defaults

    DoZip = .false. ! No zipping
    Verb  = .false. ! Run quiet
    fill  = -1      ! Fill all
    tildir='til/'   ! Write in current dir
    rstdir='rst/'   ! Write in current dir
    maxtiles=4000000
    j_start = -1    ! Process all latitudes
    j_end = -1
    subset_mode = 0 ! 0=full, 1=fast testing, 2=full hash

    I = command_argument_count()

    if(I < 2 .or. I > 17) then
       print *, trim(Usage)
       print *, ""
       print *, "Subset options:"
       print *, "  -j J1:J2   Process only latitudes J1 to J2"
       print *, "  -m MODE    Subset mode: 0=full(default), 1=fast, 2=full_hash"
       print *, "             Mode 1: Fast testing - skip unselected lats (relaxed bitwise)"
       print *, "             Mode 2: Full hash - read all, accumulate subset (strict bitwise)"
       call exit(1)
    end if

    nxt = 1
    call get_command_argument(nxt,arg)

    do while(arg(1:1)=='-')
       opt=arg(2:2)
       argl = len(trim(arg))

       if(argl==2) then
          if(scan(opt,'zvh')==0) then
             nxt = nxt + 1
             call get_command_argument(nxt,arg)
          end if
       else
          arg = arg(3:)
       end if

       select case (opt)
       case ('v')
          Verb  = .true.
       case ('z')
          DoZip = .true.
       case ('h')
          tildir = ''
          rstdir = ''
       case ('f')
          read(arg,*) fill
       case ('t')
          read(arg,*) maxtiles
       case ('g')
          Overlay = trim(arg)
       case ('s')
          sg = trim(arg)
          s_flag=.true.
       case ('j')
          ! Parse latitude range J1:J2
          i = index(arg, ':')
          if (i > 0) then
             read(arg(1:i-1),*) j_start
             read(arg(i+1:),*) j_end
          else
             print *, "Error: -j option requires format J1:J2"
             call exit(1)
          end if
       case ('m')
          read(arg,*) subset_mode
          if (subset_mode < 0 .or. subset_mode > 2) then
             print *, "Error: -m MODE must be 0, 1, or 2"
             call exit(1)
          end if
       case default
          print *, trim(Usage)
          call exit(1)
       end select

       nxt = nxt + 1
       call get_command_argument(nxt,arg)
    end do

    Grid1 = ARG

    nxt = nxt + 1
    call get_command_argument(nxt,arg)

    Grid2 = ARG

    Merge = fill/=-1


    if(trim(Overlay)=='') then
       if(Merge) then
         if (s_flag) then
           Overlay = trim(Grid1)//"-"//trim(sg)//"-"//adjustl(Grid2)
         else
           Overlay = trim(Grid1)//"-"//adjustl(Grid2)
         endif
       else
         if (s_flag) then
           Overlay = trim(Grid1)//"-"//trim(sg)//"_"//adjustl(Grid2)
         else
           Overlay = trim(Grid1)//"_"//adjustl(Grid2)
         endif
       end if
    end if

    if(DoZip) then
       TilFile = trim(tildir)//trim(Overlay)//'.til.gz'
       RstFile = trim(rstdir)//trim(Overlay)//'.rst.gz'
    else
       TilFile = trim(tildir)//trim(Overlay)//'.til'
       RstFile = trim(rstdir)//trim(Overlay)//'.rst'
    end if

! Input files:
     
    if (s_flag) then
      open (TILUNIT1,file=trim(tildir)//trim(adjustl(Grid1))//"-"//trim(sg)//'.til', & 
           form=  'formatted',                          status='old')
      open (TILUNIT2,file=trim(tildir)//trim(adjustl(Grid2))//'.til', & 
           form=  'formatted',                          status='old')

      open (RSTUNIT1,file=trim(rstdir)//trim(adjustl(Grid1))//"-"//trim(sg)//'.rst', & 
           form=  'unformatted', convert='little_endian',  status='old')
      open (RSTUNIT2,file=trim(rstdir)//trim(adjustl(Grid2))//'.rst', & 
           form=  'unformatted', convert='little_endian',  status='old')
    else
      open (TILUNIT1,file=trim(tildir)//trim(adjustl(Grid1))//'.til', & 
           form=  'formatted',                          status='old')
      open (TILUNIT2,file=trim(tildir)//trim(adjustl(Grid2))//'.til', & 
           form=  'formatted',                          status='old')

      open (RSTUNIT1,file=trim(rstdir)//trim(adjustl(Grid1))//'.rst', & 
           form=  'unformatted', convert='little_endian',  status='old')
      open (RSTUNIT2,file=trim(rstdir)//trim(adjustl(Grid2))//'.rst', & 
           form=  'unformatted', convert='little_endian',  status='old')
    endif

    call system_clock(count0,count_rate)

! Read raster sizes info from .til headers

    read(TILUNIT1,*) ip1, nf1, nx1, ny1
    read(TILUNIT2,*) ip2, nf2, nx,  ny

! Both grids must be based on same shape rasters

    _ASSERT(NX1==NX,'needs informative message')
    _ASSERT(NY1==NY,'needs informative message')

! Validate subset range
    if (j_start > 0) then
       _ASSERT(j_start >= 1 .and. j_start <= ny, 'j_start out of range')
       _ASSERT(j_end >= j_start .and. j_end <= ny, 'j_end out of range')
       if (Verb) then
          print *, 'Subset mode: ', subset_mode
          print *, 'Processing latitudes: ', j_start, ' to ', j_end, ' of ', ny
       end if
    else
       j_start = 1
       j_end = ny
    end if

! allocate space

    if(Merge) then
       rvars = 5
       nvars = 3
    else
       rvars = 6
       nvars = 7
    endif

    allocate(iTable(0:nvars,maxtiles),stat=status)
    VERIFY_(STATUS)
    allocate(rTable(1:rvars,maxtiles),stat=status)
    VERIFY_(STATUS)

    iTable = 0
    rTable = 0.0_REAL64

    allocate(rst1(nx,ny),             stat=status)
    VERIFY_(STATUS)

    allocate(cc(nx),ss(nx), rst2(nx), stat=status)
    VERIFY_(STATUS)

    allocate(Table1(6,ip1),         stat=status)
    VERIFY_(STATUS)
    allocate(Table2(6,ip2),           stat=status)
    VERIFY_(STATUS)

! The rasters cover the sphere in lat-lon 

    dx  = 360.0_REAL64/nx
    dy  = 180.0_REAL64/ny
    d2r = PI/180.0_REAL64

! sine and cosine of raster longitudes

    do i=1,nx
       lons  = d2r*(-180.0_REAL64 + (i-0.5_REAL64)*dx)
       cc(i) = cos(lons)
       ss(i) = sin(lons)
    enddo

! Read input tables

    ! Try to open binary versions first, fall back to text
    logical :: use_binary1, use_binary2
    integer :: ios
    character*256 :: BinFile1, BinFile2

    if (s_flag) then
       BinFile1 = trim(tildir)//trim(adjustl(Grid1))//"-"//trim(sg)//'.til.bin'
       BinFile2 = trim(tildir)//trim(adjustl(Grid2))//'.til.bin'
    else
       BinFile1 = trim(tildir)//trim(adjustl(Grid1))//'.til.bin'
       BinFile2 = trim(tildir)//trim(adjustl(Grid2))//'.til.bin'
    end if

    ! Check if binary files exist
    inquire(file=trim(BinFile1), exist=use_binary1)
    inquire(file=trim(BinFile2), exist=use_binary2)
    use_binary1 = use_binary1 .and. use_binary2  ! Both must exist

    if (use_binary1) then
       if(Verb) print *, 'Reading binary .til files...'
       open (TILUNIT1, file=trim(BinFile1), form='unformatted', &
             convert='little_endian', status='old')
       open (TILUNIT2, file=trim(BinFile2), form='unformatted', &
             convert='little_endian', status='old')
       call ReadTilTable(TILUNIT1, .false., ip1, Table1, Grid1, nx1, ny1, Verb)
       call ReadTilTable(TILUNIT2, .false., ip2, Table2, Grid2, nx2, ny2, Verb)
    else
       if(Verb) print *, 'Reading text .til files...'
       ! Files already open from earlier
       read(TILUNIT1,*) k
       read(TILUNIT1,*) Grid1
       read(TILUNIT1,*) nx1
       read(TILUNIT1,*) ny1
       do j=2,k
          read(TILUNIT1,*)
          read(TILUNIT1,*)
          read(TILUNIT1,*)
       end do

       if(Verb) print *, 'First  input grid: ',trim(Grid1), ip1, nx1, ny1

       read(TILUNIT2,*) k
       read(TILUNIT2,*) Grid2
       read(TILUNIT2,*) nx2
       read(TILUNIT2,*) ny2
       do j=2,k
          read(TILUNIT2,*)
          read(TILUNIT2,*)
          read(TILUNIT2,*)
       end do

       if(Verb) print *, 'Second input grid: ',trim(Grid2), ip2, nx2, ny2

       do k=1,ip1
          read(TILUNIT1,*) Table1(1:2,k),lons,lons,Table1(3:6,k)
       enddo

       do k=1,ip2
          read(TILUNIT2,*) Table2(1:2,k),lons,lons,Table2(3:6,k)
       enddo
    end if

    if(Verb) then
       call system_clock(count1)
       print *, 'Read Tables. Time = ', (count1-count0)/float(count_rate)
       count0 = count1
    end if

    ip    = 0
    Hash  = MAPL_HashCreate(8*1024)

    if(Verb) write (6, '(A)', advance='NO') ' Started Overlay'

    LATITUDES: do j=1,ny
       
       ! Determine if we process this latitude
       ! Mode 1: skip file reads outside range
       ! Mode 2: read all but accumulate only in range
       ! Mode 0: process all
       if (subset_mode == 1) then
          accumulate_this_lat = (j >= j_start .and. j <= j_end)
          if (.not. accumulate_this_lat) cycle
       else
          accumulate_this_lat = (j >= j_start .and. j <= j_end)
       end if

       lats = -90.0_REAL64 + (j - 0.5_REAL64)*dy
       da   = (sin(d2r*(lats+0.5_REAL64*dy)) - &
               sin(d2r*(lats-0.5_REAL64*dy))   )*(dx*d2r)

       vv(3)   = da*lats
       vv(4)   = da

       read(RSTUNIT1) Rst1(:,j)
       read(RSTUNIT2) Rst2

       LONGITUDES: do i=1,nx
          vv(1) = da*ss(i)
          vv(2) = da*cc(i)

          Pix1  = Rst1(i,j)
          Pix2  = Rst2(i)

! If it is an overlay or we are at the tile type of the first grid
! that is being overlayed during a merge, hash the tile indeces from
! the two grids. If it is a merge and we are at the tile type of the
! first grid that is being retained, hash it into a single pseudo bin
! of the second grid.

          if(.not.Merge) then
             k = MAPL_HashIncrement(Hash,Pix1,Pix2)
          else if (nint(Table2(1,Pix2))==fill ) then
             k = MAPL_HashIncrement(Hash,Pix1,   0)
          else
             k = MAPL_HashIncrement(Hash,   0,Pix2)
          end if

          if (k <= 0 .or. k > maxtiles) then
             write (6, '(A,I0,A,2(I0,A),2(I0,A))') 'Bad hash index k=', k, &
                  ' Pix1=', Pix1, ' Pix2=', Pix2, ' i=', i, ' j=', j
             call exit(2)
          end if

          found     = (k > 0 .and. k <= ip) ! The grid indices were in hash
          Rst1(i,j) = k

! Table variables
!
!  iTable(0)    :: Surface type
!  iTable(1)    :: tile count 
!  iTable(2)    :: I_1 I of first grid
!  iTable(3)    :: J_1
!  iTable(4)    :: I_2 I of 2nd   grid
!  iTable(5)    :: J_2
!
!  rTable(1)    :: sum of lons
!  rTable(2)    :: sum of lats
!  rTable(3)    :: area
!  rTable(4)    :: of first grid box area
!  rTable(5)    :: of 2nd   grid box area


          if(found) then    ! Bump the counter and the lon, lat, and area sums
             if (accumulate_this_lat) then
                iTable( 1,k) = iTable( 1,k) + 1
                if(.not.Merge) then
                   rTable(:4,k) = rTable(:4,k) + vv(:4)
                else
                   rTable(:3,k) = rTable(:3,k) + vv(:3)
                end if
             end if
          else if (accumulate_this_lat) then ! New tile: only create in range
             ip = ip + 1
             if(ip>maxtiles) then
                print *, "Exceeded maxtiles = ", maxtiles," at j= ",  j, &
                          " ny=", ny,  " i=", i, " nx=", nx
                print *, "Use -t option to increase."
                call exit(1)
             end if

             iTable(0 ,ip) = nint(Table2(1,Pix2))
             iTable(1 ,ip) = 1

             if(.not.Merge) then
                rTable(:4,ip) = vv(:4)
             else
                rTable(:3,ip) = vv(:3)
             end if

             if(.not.Merge) then ! This is a normal overlay of two grids
                iTable(2,ip) = nint(Table1(3,Pix1))
                iTable(3,ip) = nint(Table1(4,Pix1))
                rTable(5,ip) =      Table1(2,Pix1)
                iTable(4,ip) = nint(Table2(3,Pix2))
                iTable(5,ip) = nint(Table2(4,Pix2))
                rTable(6,ip) =      Table2(2,Pix2)
                iTable(6,ip) = nint(Table1(6,Pix1))
                iTable(7,ip) = nint(Table2(6,Pix2))
             elseif(nint(Table2(1,Pix2))==fill) then ! A use 1st grid area of a merge
                iTable(2,ip) = nint(Table1(3,Pix1))
                iTable(3,ip) = nint(Table1(4,Pix1))
                rTable(4,ip) =      Table1(2,Pix1)
                rTable(5,ip) =      Table1(2,Pix1)
             else                                    ! A retain 2nd grid area of a merge
                iTable(2,ip) = nint(Table2(3,Pix2))
                iTable(3,ip) = nint(Table2(4,Pix2))
                rTable(4,ip) =      Table2(2,Pix2)
                rTable(5,ip) =      Table2(2,Pix2)
             end if

          end if

       end do LONGITUDES    !  End raster I-loop
       if(Verb) then
          if(mod(j,200)==0) then
             write (6, '(A)', advance='NO') '.'
          end if
       end if

    enddo LATITUDES  !  End raster J-loop
    
    if(Verb) then
       call system_clock(count1)
       print *,  "Found ",ip," unique tiles. Time = ", (count1-count0)/float(count_rate)
       count0 = count1
    end if

! Done with some space and with hash

!    deallocate(Table1,Table2,Rst2,ss,cc,stat=STATUS)
    VERIFY_(STATUS)

!    call MAPL_HashDestroy(Hash)

! Compute proper longitude and latitude in degrees and compress
! the real table for WriteTiling.

    if(Verb) print *, "Computing weighted lons and lats..."

    do k=1,ip
       rTable(1,k) = atan2(rTable(1,k),rTable(2,k))/d2r
       if(rTable(4,k) < 1.0e-15_REAL64) then 
          rTable(2,k) =       rTable(3,k)/(rTable(4,k) + 1.0e-15_REAL64)
       else
          rTable(2,k) =       rTable(3,k)/rTable(4,k)
       endif
       rTable(3,k) = rTable(4,k)
       rTable(4,k) = rTable(5,k)
    end do

    if(rvars==6) then
       rTable(5,:ip) = rTable(6,:ip)
    end if

! Sort table by the first grid type in descending order.

    if(Verb) print *, "Sorting..."
    call SortTiling(Rst1,rTable(:,:ip),iTable(:,:ip))


    if(Verb) then
       call system_clock(count1)
       print *,  "Done Sorting. Time = ", (count1-count0)/float(count_rate)
       count0 = count1
    end if

! Write .til and .rst files

    if(Verb) print *, 'Writing til file...'
    if(.not.Merge) then
       call WriteTiling(TilFile,(/Grid1,Grid2/), (/nx1,nx2/),(/ny1,ny2/),(/ip1,ip2/), &
                     nx, ny,    iTable(:,:ip),rTable(:5,:ip),Dozip, Verb)
    else
       call WriteTiling(TilFile,(/Overlay/),(/nx1/), (/ny1/), (/ip1/),    &
                     nx, ny,    iTable(:,:ip),rTable(:4,:ip),Dozip, Verb)
    end if

    if(Verb) then
       call system_clock(count1)
       print *,  "Done Writing. Time = ", (count1-count0)/float(count_rate)
       count0 = count1
    end if

    if(Verb) print *, 'Writing raster file...'
    call WriteRaster(RstFile,Rst1,DoZip)

    if(Verb) then
       call system_clock(count1)
       print *,  "Done Wrting. Time = ", (count1-count0)/float(count_rate)
       count0 = count1
    end if

! Clean up

    deallocate(Rst1, iTable, rTable)

! All done

    if(Verb) print * , 'Terminated Normally'
    call exit(0)


  end program mkOverlaySimple

!===================================================================

contains

  subroutine ReadTilTable(Unit, IsText, ip, Table, GridName, nx_out, ny_out, Verb)
    ! Read a .til file (text or binary format)
    ! If IsText=.true., reads formatted ASCII
    ! If IsText=.false., reads unformatted binary
    
    integer, intent(IN)        :: Unit
    logical, intent(IN)        :: IsText
    integer, intent(IN)        :: ip
    real(REAL64), intent(OUT)  :: Table(:,:)
    character*(*), intent(OUT) :: GridName
    integer, intent(OUT)       :: nx_out, ny_out
    logical, intent(IN)        :: Verb
    
    integer :: k, num_grids
    real(REAL64) :: dummy_lon, dummy_lat
    real(REAL64) :: row(1:8)
    
    if (IsText) then
       ! Text format (original)
       read(Unit,*) k
       read(Unit,*) GridName
       read(Unit,*) nx_out
       read(Unit,*) ny_out
       do k=2,k
          read(Unit,*)
          read(Unit,*)
          read(Unit,*)
       end do
       
       do k=1,ip
          read(Unit,*) Table(1:2,k), dummy_lon, dummy_lon, Table(3:6,k)
       end do
    else
       ! Binary format
       read(Unit) num_grids
       read(Unit) GridName
       read(Unit) nx_out, ny_out
       
       do k=1,ip
          read(Unit) row(1:8)
          Table(1:2,k) = row(1:2)
          Table(3:6,k) = row(3:8)
       end do
    end if
    
  end subroutine ReadTilTable

!===================================================================
