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

  REAL(KIND=8),   parameter     :: PI        = MAPL_PI

  integer                :: IARGC
  integer                :: nxt, argl, fill
  integer                :: i, j, k, l, ip
  integer                :: STATUS, i1, i2, nvars, rvars
  integer                :: ip1, ip2
  integer                :: io, jo
  integer                :: nx1, nx2, ny1, ny2, nx, ny
  integer                :: maxtiles, hash
  integer                :: LineOcn 
  integer                :: count0,count1,count_rate

  REAL(KIND=8),     pointer     :: MOMLAT(:,:)          ! Lats of MOM's T-cell centers
  REAL(KIND=8),     pointer     :: MOMWET(:,:)          ! TMASK of MOM's grid cells

  integer,   allocatable :: RST1(:,:)
  integer,   allocatable :: RST2(:  )
  integer,   allocatable :: iTable(:,:)

  REAL(KIND=8) ,    allocatable :: Table1(:,:) 
  REAL(KIND=8) ,    allocatable :: Table2(:,:) 
  REAL(KIND=8) ,    allocatable :: rTable(:,:)
  REAL(KIND=8) ,    allocatable :: cc(:), ss(:)
  REAL(KIND=8)                  :: dx, dy, area, xc, yc, d2r, vv(4)
  REAL(KIND=8)                  :: lats, lons, da

  logical                :: DoZip
  logical                :: Verb
  logical                :: found
  logical                :: Merge
                         
  character*4            :: tildir, rstdir
  character*1            :: Opt 
  character*128          :: arg
  character*128          :: Overlay=''
  character*128          :: GridName1, GridName2
  character*128          :: Grid1, Grid2
  character*128          :: TilFile, RstFile
  character*128          :: GridFile
  character*128          :: &
      Usage = "FillMomGrid -v -h -z -t MT -g GF -f TYPE BOTTOMRASTER TOPRASTER MOM_GRIDSPEC"

  integer :: Pix1, Pix2

INCLUDE "netcdf.inc"
! Argument defaults

    DoZip = .false. ! No zipping
    Verb  = .false. ! Run quiet
    fill  = -1      ! Fill all
    tildir='til/'   ! Write in current dir
    rstdir='rst/'   ! Write in current dir
    maxtiles=4000000

    I = iargc()

    if(I < 2 .or. I > 11) then
       print *, trim(Usage)
       call exit(1)
    end if

    nxt = 1
    call getarg(nxt,arg)

    do while(arg(1:1)=='-')
       opt=arg(2:2)
       argl = len(trim(arg))

       if(argl==2) then
          if(scan(opt,'zvh')==0) then
             nxt = nxt + 1
             call getarg(nxt,arg)
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
       case default
          print *, trim(Usage)
          call exit(1)
       end select

       nxt = nxt + 1
       call getarg(nxt,arg)
    end do

    Grid1 = ARG

    nxt = nxt + 1
    call getarg(nxt,arg)

    Grid2 = ARG

    nxt = nxt + 1
    call getarg(nxt,arg)

    GridFile = arg

    if(trim(Overlay)=='') then
      print*, 'Must Provide Overlay'
      call exit(0)
    end if

    call ReadGridFile(GridFile, MOMLAT, MOMWET)

    print*, 'MOM grid dims'
    print*, size(MOMWET,dim=1), size(MOMWET,dim=2)


    if(DoZip) then
       TilFile = trim(tildir)//trim(Overlay)//'.til.gz'
       RstFile = trim(rstdir)//trim(Overlay)//'.rst.gz'
    else
       TilFile = trim(tildir)//trim(Overlay)//'.til'
       RstFile = trim(rstdir)//trim(Overlay)//'.rst'
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

    read(TILUNIT1,*) ip1, nx1, ny1
    read(TILUNIT2,*) ip2, nx,  ny

! Both grids must be based on same shape rasters

    ASSERT_(NX1==NX)
    ASSERT_(NY1==NY)

! allocate space

    allocate(rst1(nx,ny), stat=status)
    VERIFY_(STATUS)

    allocate(rst2(nx), stat=status)
    VERIFY_(STATUS)

    allocate(Table1(6,ip1),         stat=status)
    VERIFY_(STATUS)
    allocate(Table2(6,ip2),           stat=status)
    VERIFY_(STATUS)

! Read input tables

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
       if(nint(Table2(1,k)) == 0) then
          LineOcn = nint(Table2(6,k))
       endif
    enddo


    if(Verb) then
       call system_clock(count1)
       print *, 'Read Tables. Time = ', (count1-count0)/float(count_rate)
       count0 = count1
    end if


    if(Verb) write (6, '(A)', advance='NO') ' Started Overlay'

    LATITUDES: do j=1,ny

       read(RSTUNIT2) Rst1(:,j)
       read(RSTUNIT1) Rst2

       LONGITUDES: do i=1,nx

          Pix1  = Rst1(i,j)
          Pix2  = Rst2(i)
          if(Pix2 <= 0) cycle 
          io = nint(Table1(3,Pix2))
          jo = nint(Table1(4,Pix2))
          if(MOMLAT(io,jo) > 30.0 .or. MOMLAT(io,jo) < -30.0) then  ! if at higher latitudes
             if(MOMWET(io,jo) > 0.5) then ! if this is a MOM ocean point
               if(Pix1 /= LineOcn) Rst1(i,j) = LineOcn
             endif
          endif

       end do LONGITUDES    !  End raster I-loop
       if(Verb) then
          if(mod(j,200)==0) then
             write (6, '(A)', advance='NO') '.'
          end if
       end if

    enddo LATITUDES  !  End raster J-loop
    

! Write .til and .rst files

    if(Verb) print *, 'Writing raster file...'
    call WriteRaster(RstFile,Rst1,DoZip)

    if(Verb) then
       call system_clock(count1)
       print *,  "Done Wrting. Time = ", (count1-count0)/float(count_rate)
       count0 = count1
    end if

! Clean up

    deallocate(Rst1)

! All done

    if(Verb) print * , 'Terminated Normally'
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
    REAL(KIND=8), pointer                  :: LAT(:,:)
    REAL(KIND=8), pointer                  :: WET(:,:)

    integer :: STATUS, NCID, VARID, j
    integer :: SIZ_XVERT_X, SIZ_XVERT_Y
    integer :: SIZ_YVERT_X, SIZ_YVERT_Y 
    logical :: newstyle
    integer :: ID, ITMP

    Status=NF_OPEN(FILE,NF_NOWRITE,NCID)
    ASSERT_(STATUS==NF_NOERR)

    ITMP = NF_INQ_VARID    (NCID, 'x_vert_T', ID )
    newstyle = ITMP==NF_NOERR


    if( NEWSTYLE) then

       call fieldSize(NCID,'x_vert_T',SIZ_XVERT_X,1)
       call fieldSize(NCID,'y_vert_T',SIZ_YVERT_Y,2)

       allocate(LAT(SIZ_XVERT_X,SIZ_YVERT_Y),stat=STATUS)
       ASSERT_(STATUS==0)
       allocate(WET(SIZ_XVERT_X,SIZ_YVERT_Y),stat=STATUS)
       ASSERT_(STATUS==0)

       STATUS = NF_INQ_VARID     (NCID,  'y_T', VARID )
       ASSERT_(STATUS==NF_NOERR)
       status = NF_GET_VAR_DOUBLE(NCID, VARID, LAT)
       ASSERT_(STATUS==NF_NOERR)

       STATUS = NF_INQ_VARID     (NCID,  'wet', VARID )
       ASSERT_(STATUS==NF_NOERR)
       STATUS = NF_GET_VAR_DOUBLE(NCID, VARID, WET)
       ASSERT_(STATUS==NF_NOERR)

!!$       print *, 'Newstyle'
!!$       print *, 'xs: ',xvert(1,1,:)
!!$       print *, 'ys: ',yvert(1,1,:)

    else

       call fieldSize(NCID,'geolon_vert_t',SIZ_XVERT_X,1)
       call fieldSize(NCID,'geolat_vert_t',SIZ_YVERT_Y,2)


       SIZ_XVERT_X = SIZ_XVERT_X-1
       SIZ_YVERT_Y = SIZ_YVERT_Y-1
!       print *, SIZ_XVERT_X,SIZ_YVERT_Y

       allocate(LAT(SIZ_XVERT_X,SIZ_YVERT_Y),stat=STATUS)
       ASSERT_(STATUS==0)
       allocate(WET(SIZ_XVERT_X,SIZ_YVERT_Y),stat=STATUS)
       ASSERT_(STATUS==0)

       STATUS = NF_INQ_VARID     (NCID,  'geolat_t', VARID )
       ASSERT_(STATUS==NF_NOERR)
       status = NF_GET_VAR_DOUBLE(NCID, VARID, LAT)
       ASSERT_(STATUS==NF_NOERR)

       STATUS = NF_INQ_VARID     (NCID,  'wet', VARID )
       ASSERT_(STATUS==NF_NOERR)
       STATUS = NF_GET_VAR_DOUBLE(NCID, VARID, WET)
       ASSERT_(STATUS==NF_NOERR)

!!$       print *, 'Oldstyle'
!!$       print *, 'xs: ',vertx(1,1),vertx(2,1),vertx(2,2),vertx(1,2)
!!$       print *, 'ys: ',verty(1,1),verty(2,1),verty(2,2),verty(1,2)


    endif

  end subroutine READGRIDFILE

  end program FillMomGrid

!===================================================================

