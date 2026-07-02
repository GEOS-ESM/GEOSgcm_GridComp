#define I_AM_MAIN
#include "MAPL_ErrLog.h"

  program mkEASERaster

! BOP

! !PROGRAM:  mkEASERaster -- Rasterizes a EASE grid

! !INTERFACE:
!
!     mkEASERaster -x nx -y ny -g EASELabel -T type m
!
! !ARGUMENTS:
!
!     -x: Size of longitude (1st) dimension of raster. default: 8640
!     -y: Size of latitude  (2nd) dimension of raster. default: 4320
!     -g: EASElabel, which is the format EASEv2_Mxx
!     -t: Grid type used in .til file. Usually only Atmophere(=-1)
!         or Ocean(=0); default is the default of LLrasterize, currently
!         Atmosphere.
!
    use LogRectRasterizeMod,         ONLY: LRRasterize
    use MAPL

    implicit none

    character*128        :: EASELabel = ''
    character*128        :: ARG
    character*1          :: opt
    integer              :: nc_ease, nr_ease, surfaceType, nc, nr
    logical              :: UseType = .false., DoZip=.false., Verb=.false.
    logical              :: Here = .false.
    integer              :: I, J, status, nxt
    real*8,  allocatable :: xs(:,:), ys(:,:)
    real          :: x,y, xout, yout
    character*128        :: &
      Usage = "mkEASERaster -x nx -y ny -v -h -z -g EASELabel -t type"
    character*128        :: Iam = "mkEASERaster"

! Process Arguments
!------------------

    I = command_argument_count()

    if(I < 2 .or. I > 13) then
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
          Verb = .true.
       case ('h')
          Here = .true.
       case ('g')
          EASELabel = trim(arg)
       case ('t')
          read(arg,*) surfaceType
          UseType = .true.
       case default
          print *, trim(Usage)
          call exit(1)
       end select
       nxt = nxt + 1
       call get_command_argument(nxt,arg)
    end do

    call MAPL_ease_extent( EASELabel, nc_ease, nr_ease, _RC)

! Allocate and define the Cell vertices
!--------------------------------------

    allocate (xs ( nc_ease+1, nr_ease+1))
    allocate (ys ( nc_ease+1, nr_ease+1))

    do  j = 1, nr_ease+1
       do i = 1, nc_ease+1
          x = real(i-1)        -0.5
          y = real(nr_ease - j)+0.5
          call MAPL_ease_inverse(EASELabel, x, y, yout, xout)
          ys (i,j) = dble(yout)
          xs (i,j) = dble(xout)
        end do
     end do

! Produce the .rst and .til files
!--------------------------------

    if(UseType) then
       call  LRRasterize(EASELabel,xs,ys,nc=nc,nr=nr,Here=Here,Verb=Verb,   &
                         SurfaceType=surfaceType                            )
    else
       call  LRRasterize(EASELabel,xs,ys,nc=nc,nr=nr,Here=Here,Verb=Verb    )
    end if
 
   deallocate(xs,ys)

  end program mkEASERaster

