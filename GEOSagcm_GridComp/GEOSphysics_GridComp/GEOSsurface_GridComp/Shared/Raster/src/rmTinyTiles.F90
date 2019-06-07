#include "Raster.h"

 PROGRAM rmTinyTiles

! !INTERFACE:
!
!     rmTinyTiles -x nx -y ny -b lon0 -p pos -T type  im jm
!
! !ARGUMENTS:
!
!     -x: Size of longitude (1st) dimension of raster. default: 8640
!     -y: Size of latitude  (2nd) dimension of raster. default: 4320
!     -b: longitude of left edge of first box. This can be a 
!         floating point number in degrees or one of the following 
!         two character strings denoting the standard positions:
!         DE, DC, GE, GC, for Dateline or Greenwich, Center or Edge.
!         default: 'DC'
!     -p: Position of poes relative to grid.This can be a 
!         one of the following two character strings: PE or PC, for
!         Pole on the Edge, or Pole on the Center of the first box.
!         PC is the FV convention, with the firat and last boxes being
!         half as big as the others.
!         default: 'PC'
!     -t: Grid type used in .til file. Usually only Atmophere(=-1)
!         or Ocean(=0); default is the default of LLrasterize, currently
!         Atmosphere.
!     im: Size of longitude (1st) dimension of grid.
!     jm: Size of latitude  (2nd) dimension of grid.
!
! This program removes tiny tiles in
! DC${IM}xPC${JM}_DE0360xPE0180-Pfafstetter.til/
! DC${IM}xPC${JM}_DE0360xPE0180-Pfafstetter.rst, and 
! creates a new pair, 
! DC${IM}xPC${JM}_DE0360xPE0180-Pfaf.update.til/ 
! DC${IM}xPC${JM}_DE0360xPE0180-Pfaf.update.rst  
!  
! Sarith Mahanama - Feb 6, 2010
! Email: sarith.p.mahanama@nasa.gov

   use rmTinyCatchParaMod

  implicit none

    integer              :: NC = 8640, NR = 4320
    character*20         :: DQ = 'DC'
    character*2          :: PT = 'PC'
    integer              :: PUSHLEFT = 10000
    character*128        :: GridName = ''
    character*128        :: ARG
    character*1          :: opt
    integer              :: II, JJ, Type
    logical              :: UseType = .false., DoZip=.false., Verb=.false.
    logical              :: Here = .false.
    integer              :: I, J, status, iargc, nxt
    real*8               :: dx, dy, lon0
    real*8,  allocatable :: xs(:), ys(:)
    character*128        :: &
      Usage = "rmTinyTiles -x nx -y ny -v -h -z -g Gridname -b lon0 -p pos -t type  im jm"
    character*2 :: dateline,poles	 

     CHARACTER*128 :: Grid2

    print *,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    print *,'Begin removing tiny tiles.................................:'
    print *,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

 
! Process Arguments
!------------------

    I = iargc()

    if(I < 2 .or. I > 17) then
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
       case ('b')
          dq = trim(arg)
       case ('p')
          pt = arg(1:2)
       case ('z')
          DoZip = .true.
       case ('v')
          Verb = .true.
       case ('h')
          Here = .true.
       case ('g')
          GridName = trim(arg)
       case ('t')
          read(arg,*) type
          UseType = .true.
       case default
          print *, trim(Usage)
          call exit(1)
       end select
       nxt = nxt + 1
       call getarg(nxt,arg)
    end do

    read(arg,'(i5)') ii

    nxt = nxt + 1
    call getarg(nxt,arg)

    read(arg,'(i5)') jj

    nxt = nxt + 1
    call getarg(nxt,arg)
    Grid2 = ARG
 
    if(trim(Gridname) == '') then
       write(Gridname,'(A2,I4.4,A1,A2,I4.4)') dq(1:2),ii,"x",pt(1:2),jj
    endif
 
    write (dateline,'(A2)')dq(1:2)
    write (poles,   '(A2)')pt(1:2)

! Allocate and define the Cell vertices
!--------------------------------------

    allocate(xs(ii+1),ys(jj+1),stat=STATUS)
    VERIFY_(STATUS)

    dx = 360.0_8/float(ii)
   
    select case (dq(1:2))
    case ('DC')
       lon0 = -180. - dx*0.5
    case ('DE')
       lon0 = -180.
    case ('GC')
       lon0 = - dx*0.5
    case ('GE')
       lon0 = 0.0
    case default
       read(dq,*) lon0
       dq = 'XX'
    end select

    do i=1,ii+1
       xs(i) = lon0 + (i-1)*dx
    enddo

    ys(   1) = -90._8
    ys(jj+1) =  90._8

    select case (pt)
    case ('PC')
       dy = 180._8 / (jj-1)
       ys(2) = ys(1) + 0.5*dy
    case ('PE')
       dy = 180._8 / (jj  )
       ys(2) = ys(1) +     dy
    case default
       print *, " Bad pole grid type. Must be PE or PC:", PT
       print *, trim(Usage)
       call exit(1)
    end select

    Gridname =trim(Gridname)//'_'//adjustl(Grid2)
    print *,'Input tile space....:','til/'//trim(Gridname)//'-Pfafstetter'

    call remove_tiny_tiles (                         &
       dateline,poles,Gridname)

    print *,'============================================================'
    print *,'DONE removing tiny tiles....................................'
    print *,'============================================================'


END PROGRAM rmTinyTiles

