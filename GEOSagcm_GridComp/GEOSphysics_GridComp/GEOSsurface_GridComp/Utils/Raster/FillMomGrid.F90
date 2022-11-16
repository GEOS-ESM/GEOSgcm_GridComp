!#include "Raster.h"
#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif


program FillMomGrid

  use LogRectRasterizeMod
  use MAPL_SortMod
  use MAPL_HashMod
  use MAPL_ConstantsMod
  use, intrinsic :: iso_fortran_env, only: REAL64

! Modifies Pfafstetter.rst such that for every pixel within a MOM ocean
! grid cell, its value is redirected pointing to ocean (surf type 0) in 
! the corresponding .til file
! in the resulting tile file, MOM ocean grid wil be covered 100% by ocean
! tiles. Curretly this is done only for 30 deg poleward where sea ice
! could be present 

  implicit none

  integer, parameter     :: RSTUNIT1  = 20
  integer, parameter     :: RSTUNIT2  = 21
  integer, parameter     :: TILUNIT1  = 22
  integer, parameter     :: TILUNIT2  = 23

  REAL(KIND=REAL64),   parameter     :: PI        = MAPL_PI

  integer                :: command_argument_count
  integer                :: nxt, argl, fill
  integer                :: i, j, k, l, ip
  integer                :: STATUS, i1, i2, nvars, rvars
  integer                :: ip1, ip2, nf1, nf2, ip3
  integer                :: io, jo
  integer                :: im, jm
  integer                :: nx1, nx2, ny1, ny2, nx, ny
  integer                :: maxtiles, hash
  integer                :: LineOcn 
  integer                :: count0,count1,count_rate

  REAL(KIND=REAL64)      :: xmin, ymin, xmax, ymax, xs, ys

  REAL(KIND=REAL64),     pointer     :: MOMLAT(:,:)          ! Lats of MOM's T-cell centers
  REAL(KIND=REAL64),     pointer     :: MOMWET(:,:)          ! TMASK of MOM's grid cells

  integer,   allocatable :: RST1(:,:), RST3(:,:), RST0(:,:)
  integer,   allocatable :: RST2(:  ) 
  integer,   allocatable :: iTable(:,:)

  REAL(KIND=REAL64),    allocatable :: Table1(:,:) 
  REAL(KIND=REAL64),    allocatable :: Table2(:,:) 
  REAL(KIND=REAL64),    allocatable :: rTable(:,:)
  REAL(KIND=REAL64),    allocatable :: cc(:), ss(:)
  REAL(KIND=REAL64)                 :: dx, dy, area, xc, yc, d2r, vv(4)
  REAL(KIND=REAL64)                 :: lats, lons, da

  logical                :: DoZip
  logical                :: Verb
  logical                :: found
  logical                :: Merge
                         
  character*4            :: tildir, rstdir
  character*1            :: Opt 
  character*128          :: arg
  character*128          :: Overlay=''
  character*128          :: OverlayO=''
  character*128          :: GridName1, GridName2
  character*128          :: Grid1, Grid2
  character*128          :: TilFile, RstFile
  character*128          :: TilFile1
  character*128          :: GridFile
  character*128          :: &
      Usage = "FillMomGrid -v -z -t MT -g GF -f TYPE BOTTOMRASTER TOPRASTER MOM_GRIDSPEC"

  integer :: Pix1, Pix2
  integer,   allocatable, dimension (:) :: pfaf_yes, pfaf_rem, pfaf_new
  
INCLUDE "netcdf.inc"
! Argument defaults

    DoZip = .false. ! No zipping
    Verb  = .false. ! Run quiet
    fill  = -1      ! Fill all
    tildir='til/'   ! Write in current dir
    rstdir='rst/'   ! Write in current dir
    maxtiles=4000000

    I = command_argument_count()

    if(I < 2 .or. I > 11) then
       print *, trim(Usage)
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
       !case ('o')
       !   OverlayO = trim(arg)
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

    nxt = nxt + 1
    call get_command_argument(nxt,arg)

    GridFile = arg

    if(trim(Overlay)=='') then
      print*, 'Must Provide Overlay'
      call exit(0)
    end if
    !if(trim(OverlayO)=='') then
    !  print*, 'Must Provide OverlayO'
    !  call exit(0)
    !end if

    call ReadGridFile(GridFile, MOMLAT, MOMWET)

    print*, 'MOM grid dims'
    print*, size(MOMWET,dim=1), size(MOMWET,dim=2)
    print*, 'grid1: ', trim(adjustl(Grid1))
    print*, 'grid2: ', trim(adjustl(Grid2))


    if(DoZip) then
       TilFile = trim(tildir)//trim(Overlay)//'.til.gz'
       RstFile = trim(rstdir)//trim(Overlay)//'.rst.gz'
    else
       TilFile = trim(tildir)//trim(Overlay)//'.til'
       RstFile = trim(rstdir)//trim(Overlay)//'.rst'
       TilFile1 = trim(tildir)//trim(OverlayO)//'.til'
    end if

! Input files:
     
    open (TILUNIT1,file=trim(tildir)//trim(adjustl(Grid1))//'.til', & 
         form=  'formatted',                          status='old')
    open (TILUNIT2,file=trim(tildir)//trim(adjustl(Grid2))//'.til', & 
         form=  'formatted',                          status='old')

    open (RSTUNIT1,file=trim(rstdir)//trim(adjustl(Grid1))//'.rst', & 
         form=  'unformatted', convert='little_endian',  status='old')
    open (RSTUNIT2,file=trim(rstdir)//trim(adjustl(Grid2))//'.rst', & 
         form=  'unformatted', convert='little_endian',  status='old')

    call system_clock(count0,count_rate)

! Read raster sizes info from .til headers

    read(TILUNIT1,*) ip1, nf1, nx1, ny1
    read(TILUNIT2,*) ip2, nf2, nx,  ny

    print*, 'nx1 = ', nx1, ' nx = ', nx
    print*, 'ny1 = ', ny1, ' ny = ', ny

! Both grids must be based on same shape rasters

    ASSERT_(NX1==NX)
    ASSERT_(NY1==NY)

! allocate space

    allocate(rst1(nx,ny), stat=status)
    VERIFY_(STATUS)

    allocate(rst2(nx), stat=status)
    VERIFY_(STATUS)

    allocate(Table1(8,ip1),         stat=status)
    VERIFY_(STATUS)
    allocate(Table2(8,ip2),           stat=status)
    VERIFY_(STATUS)

    allocate(iTable(0:3,maxtiles),stat=status)
    VERIFY_(STATUS)
    allocate(rTable(1:4,maxtiles),stat=status)
    VERIFY_(STATUS)
    
    allocate(rst0(nx,ny), stat=status)
    VERIFY_(STATUS)
    allocate(rst3(nx,ny), stat=status)
    VERIFY_(STATUS)

! Read input tables

    read(TILUNIT1,*) k
    read(TILUNIT1,*) Grid1
    read(TILUNIT1,*) im
    read(TILUNIT1,*) jm
    do j=2,k
       read(TILUNIT1,*)
       read(TILUNIT1,*)
       read(TILUNIT1,*)
    end do

    print *, 'First  input grid: ',trim(Grid1), ip1, im, jm

    read(TILUNIT2,*) k
    read(TILUNIT2,*) Grid2
    read(TILUNIT2,*) nx2
    read(TILUNIT2,*) ny2
    do j=2,k
       read(TILUNIT2,*)
       read(TILUNIT2,*)
       read(TILUNIT2,*)
    end do

    print *, 'Second input grid: ',trim(Grid2), ip2, nx2, ny2

    do k=1,ip1
       read(TILUNIT1,*) Table1(:,k)
    enddo

    do k=1,ip2
       read(TILUNIT2,*) Table2(:,k)
       if(nint(Table2(1,k)) == 0) then
          LineOcn = nint(Table2(8,k))
       endif
    enddo

    xmin = -180.0_8
    xmax =  180.0_8
    ymin =  -90.0_8
    ymax =   90.0_8

    dx    = (xmax-xmin)/nx
    dy    = (ymax-ymin)/ny
    d2r   = (2._8*PI)/360.0_8

    if(Verb) then
       call system_clock(count1)
       print *, 'Read Tables. Time = ', (count1-count0)/float(count_rate)
       count0 = count1
    end if


    if(Verb) write (6, '(A)', advance='NO') ' Started Overlay'

    LATITUDES: do j=1,ny

       lats = -90._8 + (j - 0.5_8)*dy
       da   = (sin(d2r*(lats+0.5*dy)) - &
               sin(d2r*(lats-0.5*dy))   )*(dx*d2r)

       area  = da

       read(RSTUNIT2) Rst1(:,j)
       read(RSTUNIT1) Rst2
       Rst0 (:,j) = Rst1(:,j)
       LONGITUDES: do i=1,nx

          Pix1  = Rst1(i,j)
          Pix2  = Rst2(i)
          if(Pix2 <= 0) cycle 
          io = nint(Table1(5,Pix2))
          jo = nint(Table1(6,Pix2))
          if(MOMLAT(io,jo) > 37.5 .or. MOMLAT(io,jo) < -60.0) then  ! if at higher latitudes
             if(MOMWET(io,jo) > 0.5) then ! if this is a MOM ocean point
               if(Pix1 /= LineOcn) then
                    Rst1(i,j) = LineOcn
                    Table2(2,Pix1) = Table2(2,Pix1) - area
                    Table2(2,LineOcn) = Table2(2,LineOcn) + area
               endif
             endif
          endif

       end do LONGITUDES    !  End raster I-loop
       if(Verb) then
          if(mod(j,200)==0) then
             write (6, '(A)', advance='NO') '.'
          end if
       end if

    enddo LATITUDES  !  End raster J-loop

    
    ! reindex Pfafstetter
    ! -------------------
    allocate (pfaf_yes (1:ip2))
    allocate (pfaf_new (1:ip2))

    pfaf_yes = 0
    pfaf_new = 0
    
    do j = 1, ny
       do i = 1, nx
          pfaf_yes(Rst1(i,j))=Rst1(i,j) 
       end do
    end do
    k = count (pfaf_yes == 0)
    print *, ' # of submerged cats', k
    ip3 = ip2 - k
    
    allocate (pfaf_rem (1: ip3))
    pfaf_rem = pack (pfaf_yes, mask = (pfaf_yes > 0))

    ! new pfaf table
    
    do k=1,ip3
       iTable(0,k) = nint(Table2(1,pfaf_rem(k)))
       iTable(2:3,k) = nint(Table2(5:6,pfaf_rem(k)))
       rTable(1,k) = Table2(3,pfaf_rem(k))
       rTable(2,k) = Table2(4,pfaf_rem(k))
       rTable(3,k) = Table2(2,pfaf_rem(k))
       rTable(4,k) = rTable(3,pfaf_rem(k))
       pfaf_new (pfaf_rem(k)) =k 
    enddo

    ! new pfaf raster

    do j = 1, ny
       do i = 1, nx
          Rst3 (i,j) = pfaf_new (Rst1(i,j)) 
       end do
    end do

    ! --------------------------------
    !            SUMMARY
    ! --------------------------------
    print *, '                           Original       FillMOM         Re-indexed'
    print *, ' '
    print *, ' # of ocean pixels  ', count (Rst0 == ip2-2), count (Rst1 == ip2-2), count (Rst3 == ip3-2)
    print *, ' # of lake pixels   ', count (Rst0 == ip2-1), count (Rst1 == ip2-1), count (Rst3 == ip3-1)
    print *, ' # of landice pixels', count (Rst0 == ip2  ), count (Rst1 == ip2  ), count (Rst3 == ip3  )
    print *, ' # of land pixels   ', count (Rst0 <= ip2-3), count (Rst1 <= ip2-3), count (Rst3 <= ip3-3)
    print *, ' MINVAL             ', minval(Rst0), minval (Rst1), minval (Rst3)
    print *, ' MAXVAL             ', maxval(Rst0), maxval (Rst1), maxval (Rst3)   

    
! Write .til and .rst files
    print *, 'Writing land til file...'
    call WriteTiling(TilFile, (/Grid2/), (/ip3/), (/1/), (/ip3/),      &
                   nx, ny, iTable(:,:ip3), rTable(:4,:ip3), Dozip, Verb )

    !do k=1,ip1
    !   iTable(0,k) = nint(Table1(1,k))
    !   iTable(2:3,k) = nint(Table1(5:6,k))
    !   rTable(1,k) = Table1(3,k)
    !   rTable(2,k) = Table1(4,k)
    !   rTable(3,k) = Table1(2,k)
    !   rTable(4,k) = rTable(3,k)
    !enddo

    !print *, 'Writing ocn til file...'
    !call WriteTiling(TilFile1, (/Grid1/), (/im/), (/jm/), (/ip1/),      &
    !               nx, ny, iTable(:,:ip1), rTable(:4,:ip1), Dozip, Verb )

    print *, 'Writing raster file...'
    call WriteRaster(RstFile,Rst3,DoZip)

    if(Verb) then
       call system_clock(count1)
       print *,  "Done Wrting. Time = ", (count1-count0)/float(count_rate)
       count0 = count1
    end if

! Clean up

    deallocate(Rst1)

! All done

    print * , 'Terminated Normally'
    call exit(0)

contains

  subroutine FieldSize(NCID,name,XY,nn)
    integer, intent(IN) :: NCID
    character(*), intent(IN) :: name
    integer, intent(out) :: XY
    integer, intent(in ) :: nn

    integer :: ID, ITMP, ndims, dimid(3)

    ITMP = NF_INQ_VARID    (NCID,  NAME, ID )
!    print *, name
    ASSERT_(ITMP==NF_NOERR)

    ITMP = NF_INQ_VARNDIMS (NCID, ID, ndims)
!    print *, ndims
    ASSERT_(ITMP==NF_NOERR)
    !ASSERT_(ndims==2)

    itmp = NF_INQ_VARDIMID (NCID, ID, diMId)
!    print *, dimid
    ASSERT_(ITMP==NF_NOERR)

    itmp = NF_INQ_DIMLEN   (NCID, DIMID(nn),XY)
!    print *, Xy
    ASSERT_(ITMP==NF_NOERR)

    return
  end subroutine FieldSize

  subroutine ReadGridFile(FILE, LAT, WET)
    
    character*(*),      intent(IN ) :: FILE
    REAL(KIND=REAL64), pointer                  :: LAT(:,:)
    REAL(KIND=REAL64), pointer                  :: WET(:,:)

    integer :: STATUS, NCID, VARID
    integer :: SIZ_XVERT_X
    integer :: SIZ_YVERT_Y 
    logical :: newstyle
    integer :: ID, ITMP

    Status=NF_OPEN(FILE,NF_NOWRITE,NCID)
    ASSERT_(STATUS==NF_NOERR)

    call fieldSize(NCID,'lon_centers',SIZ_XVERT_X,1)
    call fieldSize(NCID,'lat_centers',SIZ_YVERT_Y,2)

    allocate(LAT(SIZ_XVERT_X,SIZ_YVERT_Y),stat=STATUS)
    ASSERT_(STATUS==0)
    allocate(WET(SIZ_XVERT_X,SIZ_YVERT_Y),stat=STATUS)
    ASSERT_(STATUS==0)

    STATUS = NF_INQ_VARID     (NCID,  'lat_centers', VARID )
    ASSERT_(STATUS==NF_NOERR)
    status = NF_GET_VAR_DOUBLE(NCID, VARID, LAT)
    ASSERT_(STATUS==NF_NOERR)

    STATUS = NF_INQ_VARID     (NCID,  'mask', VARID )
    ASSERT_(STATUS==NF_NOERR)
    STATUS = NF_GET_VAR_DOUBLE(NCID, VARID, WET)
    ASSERT_(STATUS==NF_NOERR)

  end subroutine READGRIDFILE


  end program FillMomGrid

!===================================================================

