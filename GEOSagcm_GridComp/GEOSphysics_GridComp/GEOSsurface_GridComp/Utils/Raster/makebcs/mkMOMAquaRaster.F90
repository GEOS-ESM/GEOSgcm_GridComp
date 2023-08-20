#define I_AM_MAIN
#include "MAPL_ErrLog.h"
program MOMraster

  use LogRectRasterizeMod
  use MAPL_ExceptionHandling
  implicit none

! this program builds a rasterized grid whose cells are 2.5 by 2.5 minutes
! one degree will contain 24 cell (60/2.5=24)
! User needs to provide inputs: x_size_deg,y_size_deg,init_lat, init_lon,cell_size
! the meaning of input is explained below, for convenience inputs can be given
! via namelist hence can be changed at runtime

  integer                :: im, jm                ! dimensions of MOM grid
  real(kind=8),     pointer     :: xvert(:,:,:)          ! Lons of MOM's vertices
  real(kind=8),     pointer     :: yvert(:,:,:)          ! Lats of MOM's vertices
  real(kind=8)                  :: xmin, xmax
  integer                :: i, j, nxt,k
  integer                :: status, command_argument_count
  character*(128)        :: GridFile
  character*(128)        :: GridName=''
  character*(128)        :: arg
  character*(128)        :: OCEAN_VERSION
  character*(2)          :: opt
  character*(128)        :: &
      Usage = "mkMOMAquaRaster -x rx -y ry -z -v -g GridName -h GridSpecFile -w OCEAN_VERSION"
  character*(128)        :: Iam = "mkMOMAquaRaster"

! argument defaults

  logical                :: DoZip =.false.
  logical                :: Verb  =.false.
  logical                :: Here  =.false.
  integer                :: Nc    = 8640
  integer                :: NR    = 4320

  real(kind=8)                  :: tol
INCLUDE "netcdf.inc"

! Process Arguments
!------------------

    I = command_argument_count()

    if(I < 1 .or. I > 8) then
       print *, "Wrong Number of arguments: ", i
       print *, trim(Usage)
       call exit(1)
    end if

    nxt = 1
    call get_command_argument(nxt,arg)
    do while(arg(1:1)=='-')
       opt=arg(2:2)
       if(len(trim(arg))==2) then
          if(scan(opt,'zvh')==0) then
             nxt = nxt + 1
             call get_command_argument(nxt,arg)
          endif
       else
          arg = arg(3:)
       end if
       select case (opt)
       case ('x')
          read(arg,'(i6)') nc
       case ('y')
          read(arg,'(i6)') nr
       case ('z')
          DoZip = .true.
       case ('v')
          Verb  = .true.
       case ('h')
          Here  = .true.
       case ('g')
          GridName = trim(arg) 
       case ('w')
          OCEAN_VERSION = trim(arg) 
       case default
          print *, trim(Usage)
          call exit(1)
       end select
       nxt = nxt + 1
       call get_command_argument(nxt,arg)
    end do

    GridFile = arg

    call ReadGridFile(GridFile,xvert,yvert)

    IM = size(xvert,1)
    JM = size(xvert,2)

    if(trim(GridName)=='')write(Gridname,'(A4,A1,I4.4,A1,I4.4)')trim(OCEAN_VERSION),"-",im,"x",jm

    if(DoZip) GridName = trim(Gridname)//'.gz'

! Report

    if(Verb) then
       print *, 'Read MOM grid from file ', GridFile
       print *, '   Raster grid sizes     : ', nc, nr
       print *, '   MOM grid sizes        : ', im, jm
    endif

! Get surface types from integer raster file

    tol = 1.0e-12 ! the default for the rasterization routine
    !Atanas:  some heuristics for very coarse MOM grids
    if (nr  > 300*jm) tol = 1.0e-5

    call LRRasterize(GridName,xvert,yvert,nc=nc,nr=nr,&
                     SurfaceType=0,Verb=Verb,Here=Here,tol=tol)

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
    _ASSERT(ITMP==NF_NOERR,'needs informative message')

    ITMP = NF_INQ_VARNDIMS (NCID, ID, ndims)
!    print *, ndims
    _ASSERT(ITMP==NF_NOERR,'needs informative message')
    !_ASSERT(ndims==2,'needs informative message')

    itmp = NF_INQ_VARDIMID (NCID, ID, diMId)
!    print *, dimid
    _ASSERT(ITMP==NF_NOERR,'needs informative message')

    itmp = NF_INQ_DIMLEN   (NCID, DIMID(nn),XY)
!    print *, Xy
    _ASSERT(ITMP==NF_NOERR,'needs informative message')

    return
  end subroutine FieldSize

  subroutine ReadGridFile(FILE,XVERT,YVERT)
    
    character*(*),      intent(IN ) :: FILE
    real(kind=8), pointer                  :: XVERT(:,:,:)
    real(kind=8), pointer                  :: YVERT(:,:,:)

    integer :: STATUS, NCID, VARID
    integer :: SIZ_XVERT_X, SIZ_XVERT_Y
    integer :: SIZ_YVERT_X, SIZ_YVERT_Y 
    real(kind=8), pointer :: VERTX(:,:),VERTY(:,:)

    Status=NF_OPEN(FILE,NF_NOWRITE,NCID)
    _ASSERT(STATUS==NF_NOERR,'needs informative message')


    call fieldSize(NCID,'lon_corners',SIZ_XVERT_X,1)
    call fieldSize(NCID,'lat_corners',SIZ_YVERT_Y,2)


    allocate(VERTX(SIZ_XVERT_X,SIZ_YVERT_Y),stat=STATUS)
    _ASSERT(STATUS==0,'needs informative message')
    allocate(VERTY(SIZ_XVERT_X,SIZ_YVERT_Y),stat=STATUS)
    _ASSERT(STATUS==0,'needs informative message')

!       print *, SIZ_XVERT_X,SIZ_YVERT_Y

    SIZ_XVERT_X = SIZ_XVERT_X-1
    SIZ_YVERT_Y = SIZ_YVERT_Y-1

    allocate(XVERT(SIZ_XVERT_X,SIZ_YVERT_Y,4),stat=STATUS)
    _ASSERT(STATUS==0,'needs informative message')
    allocate(YVERT(SIZ_XVERT_X,SIZ_YVERT_Y,4),stat=STATUS)
    _ASSERT(STATUS==0,'needs informative message')

    STATUS = NF_INQ_VARID     (NCID,  'lon_corners', VARID )
    _ASSERT(STATUS==NF_NOERR,'needs informative message')
    status = NF_GET_VAR_DOUBLE(NCID, VARID, VERTX)
    _ASSERT(STATUS==NF_NOERR,'needs informative message')

    STATUS = NF_INQ_VARID     (NCID,  'lat_corners', VARID )
    _ASSERT(STATUS==NF_NOERR,'needs informative message')
    STATUS = NF_GET_VAR_DOUBLE(NCID, VARID, VERTY)
    _ASSERT(STATUS==NF_NOERR,'needs informative message')

!!$       print *, 'Oldstyle'
!!$       print *, 'xs: ',vertx(1,1),vertx(2,1),vertx(2,2),vertx(1,2)
!!$       print *, 'ys: ',verty(1,1),verty(2,1),verty(2,2),verty(1,2)

    XVERT(:,:,1) = VERTX(1:siz_xvert_x  ,1:siz_yvert_y  )
    XVERT(:,:,2) = VERTX(2:siz_xvert_x+1,1:siz_yvert_y  )
    XVERT(:,:,3) = VERTX(2:siz_xvert_x+1,2:siz_yvert_y+1)
    XVERT(:,:,4) = VERTX(1:siz_xvert_x  ,2:siz_yvert_y+1)

    yVERT(:,:,1) = VERTy(1:siz_xvert_x  ,1:siz_yvert_y  )
    yVERT(:,:,2) = VERTy(2:siz_xvert_x+1,1:siz_yvert_y  )
    yVERT(:,:,3) = VERTy(2:siz_xvert_x+1,2:siz_yvert_y+1)
    yVERT(:,:,4) = VERTy(1:siz_xvert_x  ,2:siz_yvert_y+1)


  end subroutine READGRIDFILE

end program MOMraster
