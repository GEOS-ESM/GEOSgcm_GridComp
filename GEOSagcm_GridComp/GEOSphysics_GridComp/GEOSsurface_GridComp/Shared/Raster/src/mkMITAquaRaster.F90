

!  $Id$

#include "Raster.h"

      program MAIN

      use LogRectRasterizeMod

      implicit none

      integer, parameter      :: NFACES = 6
      integer, parameter      :: IUNIT  = 11,   OUNIT = 12
      integer, parameter      :: RKIND  = 8
      integer, parameter      :: LL = 1, LR = 2, UR = 3, UL = 4

      integer                 :: iargc
      integer                 :: Nc    = 8640, NR   = 4320

      type Ptr2
         real(kind=RKIND), pointer  :: V(:,:)
      end type Ptr2

      type(Ptr2)                    :: X(4), Y(4)
      real(kind=RKIND),     pointer :: XG(:,:  ), YG(:,:  )
      real(kind=RKIND), allocatable :: XV(:,:,:), YV(:,:,:)

      integer                 :: STATUS
      integer                 :: STATARRAY(12)
      integer                 :: IBEG, IEND
      integer                 :: K, N, i, j, NX, NY, Length, NXT
      integer, allocatable    :: RASTER(:,:)

      character*(128)         :: RasterFile, arg
      character*80            :: FACEFILE(NFACES) = (/               &
           '/tile001.mitgrid', '/tile002.mitgrid', '/tile003.mitgrid',  &
           '/tile004.mitgrid', '/tile005.mitgrid', '/tile006.mitgrid'  /)

      logical         :: DoZip =.false.
      logical         :: Verb  =.false.
      character*(128) :: GridDir
      character*(2)   :: OPT
      character*(128) :: &
          Usage = "mkMITAquaRaster -x rx -y ry -z -v GridDir"

! Get source grid directory and destination raster file names
!------------------------------------------------------------

   i = iargc()

   if(I < 2 .or. i > 7) then
      print *, "Wrong Number of arguments: ", i
      print *, trim(Usage)
      call exit(1)
   end if
   
   nxt = 1
   call getarg(nxt,arg)
   do while(arg(1:1)=='-')
      opt=arg(2:2)
      if(len(trim(arg))==2) then
         if(scan(opt,'zv')==0) then
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
      case default
         print *, trim(Usage)
         call exit(1)
      end select
      nxt = nxt + 1
      call getarg(nxt,arg)
   end do

   GridDir = arg

! Find size of a grid face
!-------------------------

   open (20,file=trim(GridDir)//'/tile001.mitgrid', status='old')
   call fstat(20,statarray)
   close (20)

   LENGTH = statarray(8)/rkind

   do k=16,20
      if(mod(length,k)==0) exit
   enddo

   if(k==21) then
      print *, 'Bad GMIT grid file', trim(GridDir)//'/tile001.mitgrid'
      call exit(1)
   end if

   nx  = nint(sqrt(length/real(k,kind=rkind)))
   ny  = nx

   LENGTH = nx*ny*rkind

   if(Verb) then
      print *
      print *, 'Processing cube ', nx-1, ' points on a side.'
   end if

    write(RasterFile,'(A2,I4.4,A3,I4.4)') "CM",nx-1,"xCM",6*(nx-1)
    if(DoZip) RasterFile = trim(RasterFile)//'.gz'

!  Initialize the raster array
!-----------------------------

   allocate(RASTER(nc,nr)  ,stat=STATUS); VERIFY_(STATUS)
   RASTER = RASTERUNDEF

! Allocate space for symmetric version of coordinates
!----------------------------------------------------

   allocate(XG(NX    ,NY    ),stat=STATUS); VERIFY_(STATUS)
   allocate(YG(NX    ,NY    ),stat=STATUS); VERIFY_(STATUS)

! Set pointers to coordinates of four corners (asymmetric)
!---------------------------------------------------------

   X(LL)%V => XG(1:NX-1,1:NY-1)
   X(LR)%V => XG(2:NX  ,1:NY-1)
   X(UR)%V => XG(2:NX  ,2:NY  )
   X(UL)%V => XG(1:NX-1,2:NY  )

   Y(LL)%V => YG(1:NX-1,1:NY-1)
   Y(LR)%V => YG(2:NX  ,1:NY-1)
   Y(UR)%V => YG(2:NX  ,2:NY  )
   Y(UL)%V => YG(1:NX-1,2:NY  )

! Allocate space for contigous vertices of one face
!--------------------------------------------------

   allocate(Xv(6*(NX-1),NY-1,4),stat=STATUS); VERIFY_(STATUS)
   allocate(Yv(6*(NX-1),NY-1,4),stat=STATUS); VERIFY_(STATUS)

   FACES: do K=1,NFACES

! Read vertcies for each face
!----------------------------

      open (IUNIT, FILE=trim(GridDir)//trim(FACEFILE(k)), &
            ACCESS='DIRECT', RECL=LENGTH, STATUS='OLD')

      read (IUNIT,REC=6) XG
      read (IUNIT,REC=7) YG

      close(IUNIT)

      if(Verb) then
         print *
         print *, 'Processing face ', K
         print *, 'MAXVAL(XG), MINVAL(XG) = ', maxval(XG), minval(XG)
         print *, 'MAXVAL(YG), MINVAL(YG) = ', maxval(YG), minval(YG)
      endif

      IBEG = 1 + (NX-1)*(K-1)
      IEND =     (NX-1)*(K  )

! Copy face corners into contiguous array expected by LRrasterize
!----------------------------------------------------------------

      CORNERS: do N=1,4
         XV(IBEG:IEND,:,N) = X(N)%V
         YV(IBEG:IEND,:,N) = Y(N)%V
      end do CORNERS

   end do FACES

   call LRRasterize(RASTERFILE,XV,YV,nc=nc,nr=nr,&
                    SurfaceType=0,Verb=Verb)

   call exit(0)

 end program MAIN
