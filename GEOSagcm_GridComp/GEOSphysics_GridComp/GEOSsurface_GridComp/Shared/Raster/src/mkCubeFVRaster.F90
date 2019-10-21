!   $Id$

#include "Raster.h"

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

!EOP

    implicit none

    integer              :: NC = 8640, NR = 4320
    character*128        :: ARG, GridName=''
    character*1          :: opt
    integer              :: ncells !  Cells on edges of cubed faces
    integer              :: I, J, N, status, iargc, nxt
    integer              :: js,jv
    real*8               :: dx, dy
    real*8,  allocatable :: xs(:,:), ys(:,:), xv(:,:,:), yv(:,:,:)
    integer              :: RasterUnit=20, GridUnit=30
    logical              :: dozip=.false.
    logical              :: Here=.false.
    logical              :: Verb=.false.
    character*128        :: Usage="mkCubeFVraster -x RX -y RY -z -h -v -g GN ncells" 

! Process Arguments
!------------------

    I = iargc()

    if(I < 1 .or. I > 10) then
       print *, Usage
       call exit(66)
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
          Verb = .true.
       case ('h')
          Here = .true.
       case ('g')
          GridName = trim(arg)
       case default
          print *, Usage
          call exit(1)
       end select

       nxt = nxt + 1
       call getarg(nxt,arg)
    end do

    read(arg,'(i6)') ncells

! Allocate and define the Cell vertices
!--------------------------------------

    allocate(xs(ncells+1,(ncells+1)*6), ys(ncells+1,(ncells+1)*6),stat=STATUS)
    VERIFY_(STATUS)


    call Get_CubedSphere_Grid(ncells+1, (ncells+1)*6, xs, ys, 0, .true.)

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

    if(DoZip) GridName = trim(Gridname)//'.gz'

! Produce the .rst and .til files
!--------------------------------


    call  LRRasterize(GridName,xs,ys,nc=nc,nr=nr,Here=Here,jseg=6, Verb=Verb)

! All done
!---------

    call exit(0)

  end program mkCubeFVRaster
