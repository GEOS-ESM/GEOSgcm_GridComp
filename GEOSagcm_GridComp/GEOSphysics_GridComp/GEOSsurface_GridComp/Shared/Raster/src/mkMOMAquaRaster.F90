
!   $Id$

#include "Raster.h"

program MOMraster

  use LogRectRasterizeMod

  implicit none

! this program builds a rasterized grid whose cells are 2.5 by 2.5 minutes
! one degree will contain 24 cell (60/2.5=24)
! User needs to provide inputs: x_size_deg,y_size_deg,init_lat, init_lon,cell_size
! the meaning of input is explained below, for convenience inputs can be given
! via namelist hence can be changed at runtime

  integer                :: im, jm                ! dimensions of MOM grid
  REAL_,     pointer     :: xvert(:,:,:)          ! Lons of MOM's vertices
  REAL_,     pointer     :: yvert(:,:,:)          ! Lats of MOM's vertices
  REAL_                  :: xmin, xmax
  integer                :: i, j, nxt,k
  integer                :: status, iargc
  character*(128)        :: GridFile
  character*(128)        :: GridName=''
  character*(128)        :: arg
  character*(2)          :: opt
  character*(128)        :: &
      Usage = "mkMOMAquaRaster -x rx -y ry -z -v -g GridName -h GridSpecFile"

! argument defaults

  logical                :: DoZip =.false.
  logical                :: Verb  =.false.
  logical                :: Here  =.false.
  integer                :: Nc    = 8640
  integer                :: NR    = 4320

INCLUDE "netcdf.inc"

! Process Arguments
!------------------

    I = iargc()

    if(I < 1 .or. I > 8) then
       print *, "Wrong Number of arguments: ", i
       print *, trim(Usage)
       call exit(1)
    end if

    nxt = 1
    call getarg(nxt,arg)
    do while(arg(1:1)=='-')
       opt=arg(2:2)
       if(len(trim(arg))==2) then
          if(scan(opt,'zvh')==0) then
             nxt = nxt + 1
             call getarg(nxt,arg)
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
       case default
          print *, trim(Usage)
          call exit(1)
       end select
       nxt = nxt + 1
       call getarg(nxt,arg)
    end do

    GridFile = arg

    call ReadGridFile(GridFile,xvert,yvert)

    IM = size(xvert,1)
    JM = size(xvert,2)

    if(trim(GridName)=='')write(Gridname,'(A2,I4.4,A3,I4.4)') "TM",im,"xTM",jm

    if(DoZip) GridName = trim(Gridname)//'.gz'

! Report

    if(Verb) then
       print *, 'Read MOM grid from file ', GridFile
       print *, '   Raster grid sizes     : ', nc, nr
       print *, '   MOM grid sizes        : ', im, jm
    endif

! Get surface types from integer raster file

    call LRRasterize(GridName,xvert,yvert,nc=nc,nr=nr,&
                     SurfaceType=0,Verb=Verb,Here=Here)

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

  subroutine ReadGridFile(FILE,XVERT,YVERT)
    
    character*(*),      intent(IN ) :: FILE
    REAL_, pointer                  :: XVERT(:,:,:)
    REAL_, pointer                  :: YVERT(:,:,:)

    integer :: STATUS, NCID, VARID, j
    integer :: SIZ_XVERT_X, SIZ_XVERT_Y
    integer :: SIZ_YVERT_X, SIZ_YVERT_Y 
    REAL_, pointer :: VERTX(:,:),VERTY(:,:)
    logical :: newstyle
    integer :: ID, ITMP

    Status=NF_OPEN(FILE,NF_NOWRITE,NCID)
    ASSERT_(STATUS==NF_NOERR)

    ITMP = NF_INQ_VARID    (NCID, 'x_vert_T', ID )
    newstyle = ITMP==NF_NOERR


    if( NEWSTYLE) then

       call fieldSize(NCID,'x_vert_T',SIZ_XVERT_X,1)
       call fieldSize(NCID,'y_vert_T',SIZ_YVERT_Y,2)

       allocate(XVERT(SIZ_XVERT_X,SIZ_YVERT_Y,4),stat=STATUS)
       ASSERT_(STATUS==0)
       allocate(YVERT(SIZ_XVERT_X,SIZ_YVERT_Y,4),stat=STATUS)
       ASSERT_(STATUS==0)

       STATUS = NF_INQ_VARID     (NCID,  'x_vert_T', VARID )
       ASSERT_(STATUS==NF_NOERR)
       status = NF_GET_VAR_DOUBLE(NCID, VARID, XVERT)
       ASSERT_(STATUS==NF_NOERR)

       STATUS = NF_INQ_VARID     (NCID,  'y_vert_T', VARID )
       ASSERT_(STATUS==NF_NOERR)
       STATUS = NF_GET_VAR_DOUBLE(NCID, VARID, YVERT)
       ASSERT_(STATUS==NF_NOERR)

!!$       print *, 'Newstyle'
!!$       print *, 'xs: ',xvert(1,1,:)
!!$       print *, 'ys: ',yvert(1,1,:)

    else

       call fieldSize(NCID,'geolon_vert_t',SIZ_XVERT_X,1)
       call fieldSize(NCID,'geolat_vert_t',SIZ_YVERT_Y,2)

       allocate(VERTX(SIZ_XVERT_X,SIZ_YVERT_Y),stat=STATUS)
       ASSERT_(STATUS==0)
       allocate(VERTY(SIZ_XVERT_X,SIZ_YVERT_Y),stat=STATUS)
       ASSERT_(STATUS==0)

!       print *, SIZ_XVERT_X,SIZ_YVERT_Y

       SIZ_XVERT_X = SIZ_XVERT_X-1
       SIZ_YVERT_Y = SIZ_YVERT_Y-1

       allocate(XVERT(SIZ_XVERT_X,SIZ_YVERT_Y,4),stat=STATUS)
       ASSERT_(STATUS==0)
       allocate(YVERT(SIZ_XVERT_X,SIZ_YVERT_Y,4),stat=STATUS)
       ASSERT_(STATUS==0)

       STATUS = NF_INQ_VARID     (NCID,  'geolon_vert_t', VARID )
       ASSERT_(STATUS==NF_NOERR)
       status = NF_GET_VAR_DOUBLE(NCID, VARID, VERTX)
       ASSERT_(STATUS==NF_NOERR)

       STATUS = NF_INQ_VARID     (NCID,  'geolat_vert_t', VARID )
       ASSERT_(STATUS==NF_NOERR)
       STATUS = NF_GET_VAR_DOUBLE(NCID, VARID, VERTY)
       ASSERT_(STATUS==NF_NOERR)

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

    endif

  end subroutine READGRIDFILE

end program MOMraster
