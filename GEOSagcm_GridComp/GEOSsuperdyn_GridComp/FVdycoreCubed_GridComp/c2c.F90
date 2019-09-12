
#define VERIFY_(A) if(MAPL_VRFY(A,Iam,__LINE__))call exit(-1)
#define ASSERT_(A) if(MAPL_ASRT(A,Iam,__LINE__))call exit(-1)

 
#if defined(__INTEL_COMPILER)
# define _FTELL ftelli8
#elif defined(__PGI)
# define _FTELL ftell64
#else
# define _FTELL ftell
#endif

program gmao_regrid

  use ESMF
  use MAPL_Mod
  use CubedSphereGridFactoryMod
  implicit none

  integer, parameter :: GridType_Unknown = 0
  integer, parameter :: GridType_LatLon = 1
  integer, parameter :: GridType_CubedSphere = 2

  character(ESMF_MAXSTR) :: str
  character(ESMF_MAXSTR) :: f_in
  type(ESMF_VM)          :: VM
  integer                :: filetype
  logical                :: changeResolution

#ifndef __GFORTRAN__
  integer*4              :: iargc
  external               :: iargc

  integer(kind=8)         :: _FTELL
  external      :: _FTELL
#endif

! ErrLog variables
!-----------------

  integer                      :: STATUS
  character(len=ESMF_MAXSTR)   :: Iam="GMAO_Regrid"

  integer :: comm
  
  type Regrid_GridInfo
     character(len=ESMF_MAXSTR) :: filename
     integer                    :: IM
     integer                    :: JM
     integer                    :: gridtype
  end type Regrid_GridInfo

  type(Regrid_GridInfo)         :: gi
  type(Regrid_GridInfo)         :: gout
  type(ESMF_Config)             :: config
  class (AbstractRegridder), pointer :: regridder
  integer                       :: nargs

  integer                       :: unit_r, unit_w
  integer                       :: i,j,k,n
  integer                       :: ic
  integer                       :: im
  integer                       :: ndes, myid

  integer(kind=8)               :: RecStart, RecEnd

  integer, parameter            :: LM = 72
  real, allocatable             :: var_in(:,:)
  real, allocatable             :: var_out(:,:)
  real                          :: pref(LM+1)

  integer, parameter            :: iA=ichar('a')
  integer, parameter            :: mA=ichar('A')
  integer, parameter            :: iZ=ichar('z')
  integer, parameter            :: i0=ichar('0')
  integer, parameter            :: i9=ichar('9')

  integer, parameter :: LatLonRes(2,8) = RESHAPE( SOURCE =&
       [  72,  46,   & ! A - 4 degree
         144,  91,   & ! B - 2 degree
         288, 181,   & ! C - 1 degree
         540, 361,   & ! D - 1/2 degree MERRA
         576, 361,   & ! D - 1/2 degree
        1152, 721,   & ! E - 1/4 degree
        1440, 721,   & ! F - 1/4 degree rectangular
        2880, 1441], & ! G - 1/8 degree rectangular  
       SHAPE = [2,8] )

  type(MAPL_NCIO) :: inNCIO,outNCIO
  integer           :: nDims, dimSizes(4),nSpatialDims
  type (ESMF_Grid) :: gridIn, gridOut

! Begin
   
  nargs = iargc() ! get command line argument info

  if (nargs /= 3) then
     call getarg(0,str)
     write(*,*) "Usage:",trim(str)," <file_in> <file_out> <Resolution(i.e. C180)>"
     call exit(2)
  end if

  call getarg(1,f_in)

  call ESMF_Initialize (vm=vm, logKindFlag=ESMF_LOGKIND_NONE, rc=status)
  VERIFY_(STATUS)

  call ESMF_VmGet(VM, localPet=myid, petCount=ndes, rc=status)
  VERIFY_(STATUS)

  if (ndes /= 1) then
     print *,''
     print *,'ERROR: currently PARALLEL jobs not supported'
     print *,''
     ASSERT_(.false.)
  end if

  if (MAPL_AM_I_Root(vm)) then
     call GuessFileType(f_in, filetype, rc=status)
     VERIFY_(STATUS)
  end if

  call MAPL_CommsBcast(vm, DATA=filetype, N=1, ROOT=0, RC=status)
  VERIFY_(STATUS)

!  print *, filetype
  gi%filename = f_in

  ! determine grid type, and compute/guess im,im
  if (filetype ==0) then
     InNCIO = MAPL_NCIOOpen(f_in,rc=status)
     VERIFY_(STATUS)
     call GetGridInfo(gi, filetype, ncinfo=InNCIO, rc=status)
     VERIFY_(STATUS)   
  else
     call GetGridInfo(gi, filetype, rc=status)
     VERIFY_(STATUS)
  end if

  call getarg(2,str)
  gout%filename = str

  call getarg(3,str)

  ! convert STR to upprecase
  do i = 1, len_trim(str)
     ic=ichar(str(i:i))
     if(ic >= iA .and. ic <= iZ) ic=ic+(mA-iA)
     str(i:i)=char(ic)
  end do
  
  if (str(1:1) == 'B') then
     gout%gridtype = GridType_LatLon
     gout%IM = 144
     gout%JM = 91
  else if (str(1:1) == 'C') then
     if (str(2:2) /= ' ') then 
        ic = ichar(str(2:2))
        ASSERT_(ic >= i0 .and. ic <= i9) ! Must be a number
        read(str(2:),*) im
        gout%gridtype = GridType_CubedSphere
        gout%IM = im
        gout%JM = 6*im
     else
        gout%gridtype = GridType_LatLon
        gout%IM = 288
        gout%JM = 181
     end if
  else if (str(1:1) == 'D') then
     gout%gridtype = GridType_LatLon
     gout%IM = 576
     gout%JM = 361
  else if (str(1:1) == 'E') then
     gout%gridtype = GridType_LatLon
     gout%IM = 1152
     gout%JM = 721
  else if (str(1:1) == 'F') then
     if (str(2:2) /= ' ') then 
        ic = ichar(str(2:2))
        ASSERT_(ic >= i0 .and. ic <= i9) ! Must be a number
        read(str(2:),*) im
        gout%gridtype = GridType_CubedSphere
        gout%IM = im
        gout%JM = 6*im
     else
        gout%gridtype = GridType_LatLon
        gout%IM = 1440
        gout%JM = 721
     end if
  else if (str(1:1) == 'G') then
     if (str(2:2) /= ' ') then 
        ic = ichar(str(2:2))
        ASSERT_(ic >= i0 .and. ic <= i9) ! Must be a number
        read(str(2:),*) im
        gout%gridtype = GridType_CubedSphere
        gout%IM = im
        gout%JM = 6*im
     else
        gout%gridtype = GridType_LatLon
        gout%IM = 2880
        gout%JM = 1441
     end if
  else
     VERIFY_(999)
  end if

  changeResolution =  gi%im /= gout%im .or. gi%jm /= gout%jm

  if (changeResolution) then

     ! create horz regridder

     gridIn = grid_manager%make_grid(CubedSphereGridFactory(im_world=gi%im,lm=1,nx=1,ny=1))
     gridOut = grid_manager%make_grid(LatLonGridFactory(im_world=gout%im, jm_world=gout%jm, lm=1, nx=1, ny=1))

     regridder => regridder_manager%make_regridder(gridIn, gridOut, REGRID_METHOD_BILINEAR, rc=status)
     VERIFY_(STATUS)

     ! allocate buffers
     allocate(var_in(gi%im, gi%jm), stat=status)
     VERIFY_(STATUS)
     allocate(var_out(gout%im, gout%jm), stat=status)
     VERIFY_(STATUS)

     if (filetype ==0) then

        call MAPL_NCIOChangeRes(InNCIO,OutNCIO,latSize=gout%jm,lonSize=gout%im,rc=status)
        VERIFY_(STATUS)
        call MAPL_NCIOSet(OutNCIO,filename=gout%filename)
        call MAPL_NCIOCreateFile(OutNCIO)
        do n=1,InNCIO%nVars
           call MAPL_NCIOVarGetDims(InNCIO,InNCIO%vars(n)%name,nDims,dimSizes,nSpatialDims=nSpatialDims)
 
           if (nSpatialDims == 2) then
              call MAPL_VarRead(InNCIO,InNCIO%vars(n)%name,var_in)
              call regridder%regrid(var_in, var_out, rc=status)
              VERIFY_(STATUS)
              call MAPL_VarWrite(OutNCIO,InNCIO%vars(n)%name,var_out)
           else if (nSpatialDims ==3) then
              do i=1,dimSizes(3) 
                 call MAPL_VarRead(InNCIO,InNCIO%vars(n)%name,var_in,lev=i)
                 call regridder%regrid(var_in, var_out, rc=status)
                 VERIFY_(STATUS)
                 call MAPL_VarWrite(OutNCIO,InNCIO%vars(n)%name,var_out,lev=i)
              end do
           end if
        enddo
        call MAPL_NCIOClose(OutNCIO)
     else
        ! open files
        UNIT_R = GetFile(gi%filename, rc=status)
        VERIFY_(STATUS)

        UNIT_W = GetFile(gout%filename, rc=status)
        VERIFY_(STATUS)

        i=0
        ! do until EOF
        do while (.true.)

        !  read record (level, slice, etc)
           read(unit_r, err=100, end=200) var_in
           i = i+1
        !  if not VertOnly
        !    transform
        !  write
           call regridder%regrid(var_in, var_out, rc=status)
           VERIFY_(STATUS)
           write(unit_w) var_out
   !        print *,'record ',i
           cycle

   100     continue
           RecEnd = _FTELL(unit_r)
           backspace(unit_r)
           RecStart = _FTELL(unit_r)
           ASSERT_(4*((LM+1)+2) == RecEnd-RecStart)
           print *,'WARNING: encoutered shorter record, assuming PREF'
           read (unit_r) pref
           write(unit_w) pref
        end do
   200  continue
        ! end do
        ! close files
        call FREE_FILE(UNIT_R)
        call FREE_FILE(UNIT_W)

     end if

     deallocate(var_out, var_in)

  else
     print *, 'No change in resolution! Nothing to be done. Copy input to output yourself!'
  end if
  call ESMF_Finalize (RC=status)
  VERIFY_(STATUS)

#undef VERIFY_
#undef ASSERT_

#include "MAPL_Generic.h"

  contains

! ================================================================

    subroutine GuessFileType(filename, filetype, rc)
      use ESMF
      use MAPL_Mod

      implicit none

      ! Arguments
      !----------
      character(len=*),  intent(IN   ) :: filename
      integer,           intent(INOUT) :: filetype
      integer, optional, intent(  OUT) :: RC

      ! ErrLog variables
      !-----------------
      
      integer                      :: STATUS
      character(len=ESMF_MAXSTR)   :: Iam="GMAO_Regrid"

      character(len=1)             :: word(4)
      character(len=1)             :: TwoWords(8)
      integer, parameter           :: hdf5(8) = (/137, 72, 68, 70, &
                                                 13, 10, 26, 10 /)
      integer                      :: irec
      integer                      :: unit
      integer                      :: i, nx, cwrd
      logical                      :: typehdf5


      UNIT = GETFILE(FILENAME, DO_OPEN=0, ALL_PES=.false., RC=STATUS)
      VERIFY_(STATUS)
      
      INQUIRE(IOLENGTH=IREC) WORD
      open (UNIT=UNIT, FILE=FILENAME, FORM='unformatted', ACCESS='DIRECT', RECL=IREC, IOSTAT=status)
      VERIFY_(STATUS)
      
! Read first 8 characters and compare with HDF5 signature
      read (UNIT, REC=1, ERR=100) TwoWords(1:4)
      read (UNIT, REC=2, ERR=100) TwoWords(5:8)
      call FREE_FILE(UNIT)

      typehdf5 = .true.
      filetype = -1 ! Unknown

      do i = 1, 8
         if (iachar(TwoWords(i)) /= hdf5(i)) then
            typehdf5 = .false.
            exit
         end if
      end do
      if (typehdf5) then
         print *, 'HDF5 file'
         filetype = 0 ! HDF5
         RETURN_(ESMF_SUCCESS)

      end if

      ! Attempt to identify as fortran binary
      cwrd = transfer(TwoWords(1:4), irec)
      ! check if divisible by 4 
      irec = cwrd/4
      filetype = irec
      if (cwrd /= 4*irec) then
         print *, "ERROR: not a Fortran binary"
         RETURN_(ESMF_FAILURE)
      end if

      RETURN_(ESMF_SUCCESS)

100   continue
      RETURN_(ESMF_FAILURE)

    end subroutine GuessFileType

    subroutine GetGridInfo(gi, filetype, ncinfo, rc)
      use ESMF
      use MAPL_Mod
      use MAPL_IOMod

      implicit none

      type(Regrid_GridInfo)         :: gi
      integer                       :: filetype
      type(MAPL_NCIO), optional, intent(in) :: ncinfo
      integer, optional, intent(OUT):: RC

      integer :: i6, im, jm
      integer :: i
      logical :: found

      gi%gridtype = GridType_Unknown
      if (filetype == 0) then
         gi%im=-1
         gi%jm=-1
         do i=1,ncinfo%ndims
            if ( trim(ncinfo%dims(i)%name) == "lon" ) then
               gi%im = ncinfo%dims(i)%len
            else if (trim(ncinfo%dims(i)%name) == "lat" ) then
               gi%jm = ncinfo%dims(i)%len
            end if
         enddo
         ASSERT_(gi%im /= -1)
         ASSERT_(gi%jm /= -1)
         if (gi%jm == gi%im*6) then
            gi%gridtype = GridType_CubedSphere
         else
            gi%gridtype = GridType_LatLon
         end if
         RETURN_(ESMF_SUCCESS)
      end if

      if (filetype == 6) then
         print *, ""
         print *,"FV binary not supported yet"
         print *, ""
         RETURN_(ESMF_FAILURE)
      end if

      ! check for cubed-sphere
      i6 = filetype/6
      if (filetype == 6*i6) then
         im = nint(sqrt(real(i6)))
         if (i6 == im*im) then
            ! cubed-sphere
            print *, 'cubed sphere C',im
            gi%gridtype = GridType_CubedSphere
            gi%im = im
            gi%jm = 6*im
            RETURN_(ESMF_SUCCESS)
         end if
      end if

      ! check for "known" lat-lon gridsizes

      found = .false.
      do i = 1, size(LatLonRes,2)
         im = LatLonRes(1,i)
         jm =  LatLonRes(2,i)
         if (filetype == im*jm) then
            gi%gridtype = GridType_LatLon
            gi%im = im
            gi%jm = jm
            found = .true.
            exit
         end if
      end do
      if (.not. found) then
         gi%gridtype = GridType_Unknown
         print *, ""
         print *, 'Unknown input grid type, assumming TILE, currently not supported'
         print *, ""
         RETURN_(ESMF_FAILURE)
      end if

      RETURN_(ESMF_SUCCESS)

    end subroutine GetGridInfo

end program gmao_regrid


