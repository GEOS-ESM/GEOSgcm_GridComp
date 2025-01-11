#define I_AM_MAIN
#include "MAPL_ErrLog.h"

      program MAIN

      use LogRectRasterizeMod,     ONLY: LRRasterize
      use MAPL_ExceptionHandling

      implicit none

      integer, parameter      :: IUNIT  = 11,   OUNIT = 12
      integer                 :: iargc

      integer, parameter      :: RKIND  = 8
      INTEGER                 :: NC
      INTEGER                 :: NX, NY

      integer                 :: STATARRAY(12)
      integer(kind=8)         :: filesize
      integer(kind=8)         :: Length
      integer                 :: K
      integer                 :: i, j
      integer                 :: KF, L, NF
      integer                 :: KP
      integer                 :: IBEG, IEND
      integer                 :: I1,IN
      integer                 :: J1,JN
      integer                 :: N
      integer                 :: nxt
      REAL(kind=RKIND), POINTER      :: XG(:,:), YG(:,:)
      real(kind=RKIND), pointer      :: XT(:,:  ), YT(:,:  )
      real(kind=RKIND), allocatable  :: XV(:,:,:), YV(:,:,:)

      character(len=128)      :: GridDir
      character(len=128)      :: &
          Usage = "mkMITAquaRaster -x rx -y ry -z -v GridDir"
      character*(128)         :: RasterFile, arg
      CHARACTER (LEN=128)     :: FACEFILE
      CHARACTER (LEN=2)       :: opt
      logical                 :: haveit
      logical                 :: isLLC
      logical                 :: DoZip =.false.
      logical                 :: Verb  =.false.
      integer                 :: pNX, pNY
      integer                 :: nFaces
      integer                 :: nTotal
      integer                 :: nCPU
      integer                 :: status
      integer, ALLOCATABLE    :: scaleX(:), scaleY(:)

      integer, parameter      :: LL = 1, LR = 2, UR = 3, UL = 4
      integer                 :: Ncol    = 8640, NRow   = 4320

      type Ptr2
         real(kind=RKIND), pointer  :: V(:,:)
      end type Ptr2

      type(Ptr2)              :: X(4), Y(4)


   ! specific values for LLC90
!      integer, parameter :: BLNKSZ = 4
!      integer, parameter :: sNX = 45
!      integer, parameter :: sNY = 45
!      integer, parameter :: BLNKSZ = 21
!      integer, parameter :: sNX = 30
!      integer, parameter :: sNY = 30
!      integer, parameter :: BLNKSZ = 42
!      integer, parameter :: sNX = 15
!      integer, parameter :: sNY = 30
!      integer, parameter :: BLNKSZ = 108
!      integer, parameter :: sNX = 15
!      integer, parameter :: sNY = 15

!      integer, dimension(*), parameter :: blankList = & ! intel 2008 feature!
!      integer, dimension(BLNKSZ), parameter :: blankList = &
! sizes for 45x45
!           [2,13,14,23]
! sizes for 30x30
!           [1,2,3,5,6,28,29,30,31,32,33,49,50,52,53,72,81,90,99,108,117]
!#15x30   nprocs = 192
!#  blankList(1:42)=
!           [1,2,3,4,5,6,9,10,11,12,55,56,57,58,59,60,61,62,63,64,65,66, &
!            97,98,99,100,103,104,105,106,143,144,&
!            161,162,179,180,197,198,215,216,233,234]
!#15x15   nprocs = 360
!#  blankList(1:108)
! [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,21,22,23,24,& 
!  65,71,75,76,90,95,96,101,102,109,110,111,112,113,114,115,116,117,118,119,&
!  120,121,122,123,124,125,126,127,128,129,130,131,132,&
!  188,189,190,193,194,195,196,199,&
!  200,201,202,203,205,206,207,208,209,211,212,213,214,215,216,242,247,253,&
!  267,268,269,270,287,288,305,306,323,324,341,342,359,360,362,376,377,378,&
!  380,381,382,395,396,400,412,413,414,430]

      integer :: sNX
      integer :: sNY
      integer, parameter :: MAXBLNKSZ = 25000
      integer            :: BLNKSZ
      integer, dimension(MAXBLNKSZ) :: blankList

      real(kind=RKIND) :: areamin, xc, yc
      character(len=128)      :: Iam = "mkMITAquaRaster"
      NAMELIST /W2_EXCH2_PARM01/   sNx, SNy,  blankList

! Get source grid directory and destination raster file names
!------------------------------------------------------------

   i = command_argument_count()

   if(I < 2 .or. i > 7) then
      print *, "Wrong Number of arguments: ", i
      print *, trim(Usage)
      call exit(1)
   end if
   
   nxt = 1
   call get_command_argument(nxt,arg)
   do while(arg(1:1)=='-')
      opt=arg(2:2)
      if(len(trim(arg))==2) then
         if(scan(opt,'zv')==0) then
            nxt = nxt + 1
            call get_command_argument(nxt,arg)
         endif
      else
         arg = arg(3:)
      end if

      select case (opt)
      case ('x')
         read(arg,'(i6)') ncol
      case ('y')
         read(arg,'(i6)') nrow
      case ('z')
         DoZip = .true.
      case ('v')
         Verb  = .true.
      case default
         print *, trim(Usage)
         call exit(1)
      end select
      nxt = nxt + 1
      call get_command_argument(nxt,arg)
   end do

   GridDir = arg

   blanklist = 0
   OPEN(UNIT=iUnit, file='data.exch2', form='formatted', status='old')
   READ(UNIT=iUnit,NML=W2_EXCH2_PARM01)
   CLOSE(iUnit)
   BLNKSZ =  count(blanklist /= 0)

      ! Open Facet 3 first. It is always a square (CS or LLC)
      open (IUNIT,file=trim(GridDir)//'/tile003.mitgrid', status='old')
      call fstat(IUNIT,statarray)
      close (IUNIT)
      filesize = statarray(8)

      !ALT: Kludge for LLC4320
      if (filesize <= 0) filesize = 2389893248
!      print *,'file size=',filesize

      LENGTH = filesize/rkind

      do k=16,20
         if(mod(length,k)==0) exit
      enddo

      if(k==21) then
         print *, 'Bad GMIT grid file', trim(GridDir)//'/tile003.mitgrid'
         call exit(1)
      end if

      nc  = nint(sqrt(length/real(k,kind=rkind)))
 !     nc = 4321

      nx  = nc-1
      ny  = nc-1

      LENGTH = nx*ny*rkind

      ! Open Facet 1 to check sizes CS or LLC)
      open (IUNIT,file=trim(GridDir)//'/tile001.mitgrid', status='old')
      call fstat(IUNIT,statarray)
      close (IUNIT)

      filesize = statarray(8)

      !ALT: Kludge for LLC4320
      if (filesize <= 0) filesize = 7168573568
!      print *,'file size=',filesize

      LENGTH = filesize/(rkind * k)

      if (LENGTH == NC*NC) then ! cubed-sphere
         isLLC = .false.
      elseif (LENGTH == (3*(NC-1)+1)*NC) then ! LLC
         isLLC = .true.
      else
         print *, 'ERROR: Unknown grid type'
         call exit(1)
      end if

      if(Verb) then
         print *
         if (isLLC) then
            print *, 'Processing LLC ', nx
         else
            print *, 'Processing CS ', nx
         end if
      end if


! Find size of a grid face
!-------------------------

      if (isLLC) then
         nFaces = 5
      else
         ! CS
         nFaces = 6
      end if
      allocate(scaleX(nFaces), scaleY(nFaces), stat=status)
      scaleX = 1
      scaleY = 1
      if (isLLC) then
         scaleX = [1,1,1,3,3]
         scaleY = [3,3,1,1,1]
      end if


      pNX = NX/sNX
      pNY = NY/sNY

      nTotal = dot_product(scaleX*pNX,scaleY*pNY)
      nCPU = ntotal-BLNKSZ

      write(RasterFile,'(A2,I4.4,A3,I4.4)') "LLC",nCPU*sNX,"xLLC",sNY
      if(DoZip) RasterFile = trim(RasterFile)//'.gz'

!  Initialize the raster array
!-----------------------------

!   allocate(RASTER(nc,nr)  ,stat=STATUS); VERIFY_(STATUS)
!   RASTER = RASTERUNDEF

! Allocate space for symmetric version of coordinates
!----------------------------------------------------

   allocate(XT(sNX+1, sNY+1 ),stat=STATUS); VERIFY_(STATUS)
   allocate(YT(sNX+1, sNY+1 ),stat=STATUS); VERIFY_(STATUS)

! Set pointers to coordinates of four corners (asymmetric)
!---------------------------------------------------------

   X(LL)%V => XT(1:sNX  ,1:sNY  )
   X(LR)%V => XT(2:sNX+1,1:sNY  )
   X(UR)%V => XT(2:sNX+1,2:sNY+1)
   X(UL)%V => XT(1:sNX  ,2:sNY+1)

   Y(LL)%V => YT(1:sNX  ,1:sNY  )
   Y(LR)%V => YT(2:sNX+1,1:sNY  )
   Y(UR)%V => YT(2:sNX+1,2:sNY+1)
   Y(UL)%V => YT(1:sNX  ,2:sNY+1)

! Allocate space for contigous vertices of one face
!--------------------------------------------------

   allocate(Xv(nCPU*sNX,sNY,4),stat=STATUS); VERIFY_(STATUS)
   allocate(Yv(nCPU*sNX,sNY,4),stat=STATUS); VERIFY_(STATUS)

   XV = 0.0
   YV = 0.0

   KF = 0
   KP = 0
   FACES: do K=1,NFACES
      write(facefile,"(A,I3.3,A)") 'tile',K,'.mitgrid'
      if(Verb) then
         print *, trim(facefile)
      end if



! Allocate space for symmetric version of coordinates
!----------------------------------------------------

      allocate(XG(scaleX(k)*NX+1, scaleY(k)*NY+1),stat=STATUS)
      VERIFY_(STATUS)
      allocate(YG(scaleX(k)*NX+1, scaleY(k)*NY+1),stat=STATUS)
      VERIFY_(STATUS)

      xg=0.0
      yg=0.0

      LENGTH = size(XG)*rkind
!      print *,'DEBUG:length=',length, rkind

! Read vertcies for each face
!----------------------------

      open (IUNIT, FILE=trim(GridDir)//trim(FACEFILE), &
            ACCESS='DIRECT', RECL=LENGTH, STATUS='OLD',convert='big_endian')

!      read (IUNIT,REC=5) rA 
      read (IUNIT,REC=6) XG
      read (IUNIT,REC=7) YG


      close(IUNIT)

      if(Verb) then
         PRINT *, 'MIN/MAXVAL(XG) = ', MINVAL(XG), MAXVAL(XG)
         PRINT *, 'MIN/MAXVAL(YG) = ', MINVAL(YG), MAXVAL(YG)
      end if

      DO J = 1, scaleY(k)*pNY
         DO I = 1, scaleX(k)*pNX
            KF = KF + 1
            ! check if KF is in the blackList. If YES, skip
            haveIt = .false.
            DO L = 1, BLNKSZ
               if (blanklist(L) > KF) exit
               if (blanklist(L) == KF) then
                  haveIt = .true.
                  exit
               end if
            END DO
            if (haveIt) then
!               if(Verb) then
!                  print *, 'skipping ',kf
!               end if
               cycle
            end if
            KP = KP + 1

            I1 = (I-1)*sNX + 1
            IN =  I   *sNX + 1
            J1 = (J-1)*sNY + 1
            JN =  J   *sNY + 1
            XT = XG(I1:IN,J1:JN)
            YT = YG(I1:IN,J1:JN)

            IBEG = 1 + sNX*(KP-1)
            IEND =     sNX*(KP  )

            ! Copy face corners into contiguous array expected by LRrasterize
            !----------------------------------------------------------------

            CORNERS: do N=1,4
               XV(IBEG:IEND,:,N) = X(N)%V
               YV(IBEG:IEND,:,N) = Y(N)%V
            end do CORNERS

         END DO
      END DO
      deallocate(YG)
      deallocate(XG)
   end do FACES
   if(Verb) then
      print *,'KP=',kp
   end if

   open(unit=9, file='mit-llc90.ascii', form='formatted', status='new')
   write(9,*) XV(:,:,1)
   write(9,*) YV(:,:,1)
   close(9)
   !stop 'writing coordinates only, exiting...'

   call LRRasterize(RASTERFILE,XV,YV,nc=ncol,nr=nrow,&
                    SurfaceType=0,Verb=Verb)

   call exit(0)

 end program MAIN
