
#define _VERIFY(A) if(MAPL_VRFY(A,Iam,__LINE__))call exit(-1)
#define _ASSERT(A) if(MAPL_ASRT(A,Iam,__LINE__))call exit(-1,'needs informative message')

 
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
  use, intrinsic :: iso_fortran_env

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
  type(MAPL_HorzTransform)      :: Trans
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
  real(real64), allocatable     :: var64_in(:,:)
  real(real64), allocatable     :: var64_out(:,:)
  real                          :: pref(LM+1)
  real(real64)                  :: pref64(LM+1)

  integer, parameter            :: iA=ichar('a')
  integer, parameter            :: mA=ichar('A')
  integer, parameter            :: iZ=ichar('z')
  integer, parameter            :: i0=ichar('0')
  integer, parameter            :: i9=ichar('9')

  integer, parameter :: LatLonRes(2,6) = RESHAPE( SOURCE =&
       [  72,  46,   & ! A - 4 degree
         144,  91,   & ! B - 2 degree
         288, 181,   & ! C - 1 degree
         540, 361,   & ! D - 1/2 degree MERRA
         576, 361,   & ! D - 1/2 degree
        1152, 721] , &  ! E - 1/4 degree
       SHAPE = [2,6] )

  type(MAPL_NCIO) :: inNCIO,outNCIO
  integer           :: nDims, dimSizes(3)
  character(len=ESMF_MAXSTR) :: gridnamef, gridnameo, tileFile
  integer :: isFV3
  integer :: headr1(6)
  integer :: headr2(5)

! Begin
   
  nargs = iargc() ! get command line argument info

  if (nargs /= 4) then
     call getarg(0,str)
     write(*,*) "Usage:",trim(str)," <file_in> <file_out> <Resolution(i.e. C180)> <isFV3>"
     call exit(2)
  end if

  call getarg(1,f_in)

  call ESMF_Initialize (vm=vm, logKindFlag=ESMF_LOGKIND_NONE, rc=status)
  _VERIFY(STATUS)

  call ESMF_VmGet(VM, localPet=myid, petCount=ndes, rc=status)
  _VERIFY(STATUS)

  if (ndes /= 1) then
     print *,''
     print *,'ERROR: currently PARALLEL jobs not supported'
     print *,''
     _ASSERT(.false.,'needs informative message')
  end if

  if (MAPL_AM_I_Root(vm)) then
     call GuessFileType(f_in, filetype, rc=status)
     _VERIFY(STATUS)
  end if

  call MAPL_CommsBcast(vm, DATA=filetype, N=1, ROOT=0, RC=status)
  _VERIFY(STATUS)

!  print *, filetype
  gi%filename = f_in

  ! determine grid type, and compute/guess im,im
  if (filetype ==0) then
     InNCIO = MAPL_NCIOOpen(f_in,rc=status)
     _VERIFY(STATUS)
     call GetGridInfo(gi, filetype, f_in, ncinfo=InNCIO, rc=status)
     _VERIFY(STATUS)   
  else
     call GetGridInfo(gi, filetype, f_in, rc=status)
     _VERIFY(STATUS)
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
        _ASSERT(ic >= i0 .and. ic <= i9,'needs informative message') ! Must be a number
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
  else
     _VERIFY(999)
  end if

  changeResolution =  gi%im /= gout%im .or. gi%jm /= gout%jm
  if (gout%gridtype /= GridType_CubedSphere .and. gi%gridtype /= GridType_CubedSphere) then
     _ASSERT(.false.,'needs informative message')
  end if
  call MAPL_GenGridName(gi%im, gi%jm, gridname=gridnamef, geos_style=.false.)
  call MAPL_GenGridName(gout%im, gout%jm, gridname=gridnameo, geos_style=.false.)
  tileFile=trim(adjustl(gridnamef)) // '_' //trim(adjustl(gridnameo))  // '.bin'

  call getarg(4,str)
  read(str,'(I10)'),isFV3

  if (changeResolution) then

     ! create horz transform
     call MAPL_HorzTransformCreate (Trans, tileFile, gridnamef, gridnameo, rc=STATUS)
     _VERIFY(STATUS)

     ! allocate buffers
     if (isFV3 == 0) then
        allocate(var_in(gi%im, gi%jm), stat=status)
        _VERIFY(STATUS)
        allocate(var_out(gout%im, gout%jm), stat=status)
        _VERIFY(STATUS)
     else
        allocate(var_in(gi%im, gi%jm), stat=status)
        _VERIFY(STATUS)
        allocate(var_out(gout%im, gout%jm), stat=status)
        _VERIFY(STATUS)
        allocate(var64_in(gi%im, gi%jm), stat=status)
        _VERIFY(STATUS)
        allocate(var64_out(gout%im, gout%jm), stat=status)
        _VERIFY(STATUS)
     end if

     if (isFV3 /=0) then

        if (filetype ==0) then

           call MAPL_NCIOChangeRes(InNCIO,OutNCIO,latSize=gout%jm,lonSize=gout%im,rc=status)
           _VERIFY(STATUS)
           call MAPL_NCIOSet(OutNCIO,filename=gout%filename)
           call MAPL_NCIOCreateFile(OutNCIO)
           do n=1,InNCIO%nVars
              call MAPL_NCIOVarGetDims(InNCIO,InNCIO%vars(n)%name,nDims,dimSizes)
              if (ndims ==1) then
                 call MAPL_VarRead(inNCIO,InNCIO%vars(n)%name,pref64)
                 call MAPL_VarWrite(OutNCIO,InNCIO%vars(n)%name,pref64)
              else if (ndims ==2) then
                 call MAPL_VarRead(InNCIO,InNCIO%vars(n)%name,var64_in)
                 var_in=var64_in
                 call MAPL_HorzTransformRun(Trans, var_in, var_out, RC=status)
                 _VERIFY(STATUS)
                 var64_out=var_out
                 call MAPL_VarWrite(OutNCIO,InNCIO%vars(n)%name,var64_out)
              else if (ndims ==3) then
                 do i=1,dimSizes(3) 
                    call MAPL_VarRead(InNCIO,InNCIO%vars(n)%name,var64_in,lev=i)
                    var_in=var64_in
                    call MAPL_HorzTransformRun(Trans, var_in, var_out, RC=status)
                    _VERIFY(STATUS)
                    var64_out=var_out
                    call MAPL_VarWrite(OutNCIO,InNCIO%vars(n)%name,var64_out,lev=i)
                 end do
              end if
           enddo
           call MAPL_NCIOClose(OutNCIO)
        else
           ! open files
           UNIT_R = GetFile(gi%filename, rc=status)
           _VERIFY(STATUS)

           UNIT_W = GetFile(gout%filename, rc=status)
           _VERIFY(STATUS)

           i=0
           ! do until EOF
           read(unit_r)headr1
           write(unit_w)headr1
           read(unit_r)headr2
           headr2(1) = gout%im
           headr2(2) = gout%jm
           write(unit_w)headr2
           read(unit_r)pref64 
           write(unit_w)pref64
           read(unit_r)pref64 
           write(unit_w)pref64
           do while (.true.)

           !  read record (level, slice, etc)
              read(unit_r, end=400) var64_in
              var_in=var64_in
              i = i+1
              call MAPL_HorzTransformRun(Trans, var_in, var_out, RC=status)
              _VERIFY(STATUS)
              var64_out=var_out
              write(unit_w) var64_out
              cycle

           end do
      400  continue
           ! end do
           ! close files
           call FREE_FILE(UNIT_R)
           call FREE_FILE(UNIT_W)

        end if

        deallocate(var_out, var_in, var64_out, var64_in)

     else

        if (filetype ==0) then

           call MAPL_NCIOChangeRes(InNCIO,OutNCIO,latSize=gout%jm,lonSize=gout%im,rc=status)
           _VERIFY(STATUS)
           call MAPL_NCIOSet(OutNCIO,filename=gout%filename)
           call MAPL_NCIOCreateFile(OutNCIO)
           do n=1,InNCIO%nVars
              call MAPL_NCIOVarGetDims(InNCIO,InNCIO%vars(n)%name,nDims,dimSizes)
              if (ndims ==1) then

              else if (ndims ==2) then
                 call MAPL_VarRead(InNCIO,InNCIO%vars(n)%name,var_in)
                 call MAPL_HorzTransformRun(Trans, var_in, var_out, RC=status)
                 _VERIFY(STATUS)
                 call MAPL_VarWrite(OutNCIO,InNCIO%vars(n)%name,var_out)
              else if (ndims ==3) then
                 do i=1,dimSizes(3) 
                    call MAPL_VarRead(InNCIO,InNCIO%vars(n)%name,var_in,lev=i)
                    call MAPL_HorzTransformRun(Trans, var_in, var_out, RC=status)
                    _VERIFY(STATUS)
                    call MAPL_VarWrite(OutNCIO,InNCIO%vars(n)%name,var_out,lev=i)
                 end do
              end if
           enddo
           call MAPL_NCIOClose(OutNCIO)
        else
           ! open files
           UNIT_R = GetFile(gi%filename, rc=status)
           _VERIFY(STATUS)

           UNIT_W = GetFile(gout%filename, rc=status)
           _VERIFY(STATUS)

           i=0
           ! do until EOF
           do while (.true.)

           !  read record (level, slice, etc)
              read(unit_r, err=100, end=200) var_in
              i = i+1
           !  if not VertOnly
           !    transform
           !  write
              call MAPL_HorzTransformRun(Trans, var_in, var_out, RC=status)
              _VERIFY(STATUS)
              write(unit_w) var_out
      !        print *,'record ',i
              cycle

      100     continue
              RecEnd = _FTELL(unit_r)
              backspace(unit_r)
              RecStart = _FTELL(unit_r)
              _ASSERT(4*((LM+1)+2) == RecEnd-RecStart,'needs informative message')
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
 
     end if

     call MAPL_HorzTransformDestroy(Trans,rc=STATUS)
     _VERIFY(STATUS)
  else
     print *, 'No change in resolution! Nothing to be done. Copy input to output yourself!'
  end if
  call ESMF_Finalize (RC=status)
  _VERIFY(STATUS)

#undef _VERIFY
#undef _ASSERT

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
      _VERIFY(STATUS)
     
      INQUIRE(IOLENGTH=IREC) WORD
      open (UNIT=UNIT, FILE=FILENAME, FORM='unformatted', ACCESS='DIRECT', RECL=IREC, IOSTAT=status)
      _VERIFY(STATUS)
      
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
         _RETURN(ESMF_SUCCESS)

      end if

      ! Attempt to identify as fortran binary
      cwrd = transfer(TwoWords(1:4), irec)
      ! check if divisible by 4 
      irec = cwrd/4
      filetype = irec
      if (cwrd /= 4*irec) then
         print *, "ERROR: not a Fortran binary"
         _RETURN(ESMF_FAILURE)
      end if

      _RETURN(ESMF_SUCCESS)

100   continue
      _RETURN(ESMF_FAILURE)

    end subroutine GuessFileType

    subroutine GetGridInfo(gi, filetype, filename, ncinfo, rc)
      use ESMF
      use MAPL_Mod
      use MAPL_IOMod

      implicit none

      type(Regrid_GridInfo)         :: gi
      integer                       :: filetype
      character(len=*)              :: filename
      type(MAPL_NCIO), optional, intent(in) :: ncinfo
      integer, optional, intent(OUT):: RC

      integer :: i6, im, jm
      integer :: i
      logical :: found
      integer :: headr1(6)
      integer :: headr2(5)
      integer :: unit

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
         _ASSERT(gi%im /= -1,'needs informative message')
         _ASSERT(gi%jm /= -1,'needs informative message')
         if (gi%jm == gi%im*6) then
            gi%gridtype = GridType_CubedSphere
         else
            gi%gridtype = GridType_LatLon
         end if
         _RETURN(ESMF_SUCCESS)
      end if

      if (filetype == 6) then
         UNIT = GETFILE(FILENAME, DO_OPEN=0, ALL_PES=.false., RC=STATUS)
         _VERIFY(STATUS)
         open (UNIT=UNIT, FILE=FILENAME, FORM='unformatted',  IOSTAT=status)
         _VERIFY(STATUS)
         read(unit)headr1
         read(unit)headr2
         close(unit)
         gi%im=headr2(1)
         gi%jm=headr2(2)
         _RETURN(ESMF_SUCCESS)
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
            _RETURN(ESMF_SUCCESS)
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
         _RETURN(ESMF_FAILURE)
      end if

      _RETURN(ESMF_SUCCESS)

    end subroutine GetGridInfo

end program gmao_regrid


