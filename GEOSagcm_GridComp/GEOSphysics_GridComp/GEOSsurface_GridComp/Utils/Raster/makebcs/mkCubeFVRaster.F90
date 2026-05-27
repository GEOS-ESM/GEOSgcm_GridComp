#define I_AM_MAIN
#include "MAPL_ErrLog.h"

  program mkCubeFVRaster

!BOP

! !ROUTINE:  mkCubeFVraster -- Rasterizes FV cubed sphere grid
!
! !INTERFACE:
!
!     mkCubeFVraster -x RX -y RY -z -h -v -g GN ncells
!
! !DESCRIPTION:
!
! !USES:
!
    use CubedSphere_GridMod
    use LogRectRasterizeMod
    use MAPL_ExceptionHandling
!EOP

    implicit none

    integer              :: NC = 8640, NR = 4320
    character*128        :: ARG, SG, GridName=''
    character*1          :: opt
    integer              :: ncells !  Cells on edges of cubed faces
    integer              :: I, J, N, status, command_argument_count, nxt
    integer              :: js,jv
    real*8               :: dx, dy
    real*8,  allocatable :: xs(:,:), ys(:,:), xv(:,:,:), yv(:,:,:)
    integer              :: RasterUnit=20, GridUnit=30
    logical              :: g_case=.false.
    logical              :: s_case=.false.
    logical              :: dozip=.false.
    logical              :: Here=.false.
    logical              :: Verb=.false.
    character*128        :: Usage="mkCubeFVraster -x RX -y RY -z -h -v -g GN -s SG ncells" 
    character*128        :: Iam ="mkCubeFVraster"
    type(stretch) :: stg
    logical :: tlon=.false.
    logical :: tlat=.false.
    logical :: tfac=.false.

! Process Arguments
!------------------

    I = command_argument_count()

    if(I < 1 .or. I > 16) then
       print *, Usage
       call exit(66)
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
      case ('X')
          read(arg,'(F12.4)') stg%target_lon
          tlon = .true.
       case ('Y')    
          read(arg,'(F12.4)') stg%target_lat
          tlat = .true.
       case ('F') 
          read(arg,'(F12.4)') stg%stretch_factor
          tfac = .true.
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
       case ('s')
          SG = trim(arg)
          s_case=.true.
       case ('g')
          GridName = trim(arg)
          g_case=.true. 
       case default
          print *, Usage
          call exit(1)
       end select

       nxt = nxt + 1
       call get_command_argument(nxt,arg)
    end do

    read(arg,'(i6)') ncells

! Allocate and define the Cell vertices
!--------------------------------------
    allocate(xs(ncells+1,(ncells+1)*6), ys(ncells+1,(ncells+1)*6),stat=STATUS)
    VERIFY_(STATUS)

    if (tlon .and. tlat .and. tfac) then
       call Get_CubedSphere_Grid(ncells+1, (ncells+1)*6, xs, ys, 0, stg=stg)
    else
       call Get_CubedSphere_Grid(ncells+1, (ncells+1)*6, xs, ys, 0, .true.)
    end if

    print *, 'Finished Cube Grid.'

! Allocate and define 2-D array of quadrilaterals
!------------------------------------------------

!!$    allocate(xv(ncells,ncells*6,4), yv(ncells,ncells*6,4),stat=STATUS)
!!$    VERIFY_(STATUS)
!!$
!!$    do n=0,5
!!$       do j=1,ncells
!!$          do i=1,ncells
!!$             jv         =  ncells   *n + j
!!$             js         = (ncells+1)*n + j
!!$
!!$             xv(i,jv,1) = xs(i  ,js  )
!!$             xv(i,jv,2) = xs(i+1,js  )
!!$             xv(i,jv,3) = xs(i+1,js+1)
!!$             xv(i,jv,4) = xs(i  ,js+1)
!!$             yv(i,jv,1) = ys(i  ,js  )
!!$             yv(i,jv,2) = ys(i+1,js  )
!!$             yv(i,jv,3) = ys(i+1,js+1)
!!$             yv(i,jv,4) = ys(i  ,js+1)
!!$          end do
!!$       end do
!!$    end do
!!$
!!$! Free the space for the vertices
!!$!--------------------------------
!!$
!!$    deallocate(xs,ys)

! Set the GridName
!-----------------
    if(trim(GridName)=='') write(GridName,'(A,I4.4,A)') 'CF',ncells,'x6C'

    if (.not.(g_case) .and. s_case) then
        write(GridName,'(A2,I4.4,A4,A5)') 'CF',ncells,'x6C-',trim(SG)
    endif

    if(DoZip) GridName = trim(Gridname)//'.gz'

! Produce the .rst and .til files
!--------------------------------


    call  LRRasterize(GridName,xs,ys,nc=nc,nr=nr,Here=Here,jseg=6, Verb=Verb)

! All done
!---------

    call exit(0)

  end program mkCubeFVRaster
