! c2l_CFIO_offline - A variant of Example 2 using exception handling.
!
! 1. Initialize the ESMF
! 2. Create a grid
! 3. Use the time manager to create clock and time objects
! 3. Read bundles form a file
! 4. Writes a bundle to a file, changing resolution
!
! Needs 6 processors to run as cube 2 latlon interpolator
!
! % mprirun -n 6 c2l_CFIO_offline.x
!
! Arlindo da Silva <arlindo.dasilva@nasa.gov>, October 2008
! Update for c2l: William.M.Putman@nasa.gov, July 2009
!----------------------------------------------------------------------------

!                                 The exception handling macros ae here
#  include "MAPL_Generic.h"

   Program c2l_CFIO_offline

   use ESMF
   use MAPL_Mod

   implicit NONE

   type(ESMF_Grid), EXTERNAL             :: AppGridCreateF

!  Basic ESMF objects being used in this example
!  ---------------------------------------------
   type(ESMF_Grid)         :: grid3d         ! Cubed-Sphere Grid
   type(ESMF_Clock)        :: clock          ! used for output file
   type(ESMF_Time)         :: Time(4)        ! Time objects
   type(ESMF_TimeInterval) :: TimeStep       ! used to define a clock
   type(ESMF_VM)           :: vm             ! ESMF Virtual Machine

!  Bundles to hold data to be read in
!  ----------------------------------
   type(ESMF_FieldBundle) :: Bundle3D

!  The CFIO object associated with a disk file
!  -------------------------------------------
   type(MAPL_CFIO) :: cfio

!  Hardwired file names
!  --------------------
   character(len=120) ::  in_Filename
   character(len=120) :: out_Filename 

!  Basic information about the parallel environment
!         PET = Persistent Execution Threads
!  In the current implementation, a PET is equivalent 
!  to an MPI process
!  ------------------------------------------------
   integer :: myPET   ! The local PET number
   integer :: nPET    ! The total number of PETs you are running on

   integer :: rc
   integer :: i, j, k

   integer :: Nx = 1, Ny=6                             ! Layout
  !integer :: IM_World=360, JM_World=2160, LM_World=0   ! Grid dimensions
   integer :: IM_World, JM_World
   integer :: nlon, nlat
   integer :: LM_World=1

   integer :: yr, mon, day, hr, mn

   integer, pointer :: resolution(:)  ! for output fie

   character(len=120) :: str_arg
#ifndef __GFORTRAN__
   external :: getarg, iargc
   integer iargc
#endif
   integer nargs

!  Defined required local variables for MAPL_Exceptions
!  ----------------------------------------------------
   __Iam__('c2l_CFIO_offline')

!                             -----

    call Main()

CONTAINS

    subroutine Main()

    nargs = IARGC()
    if ((nargs /= 11) .and. (nargs /= 12)) then
       __raise__(MAPL_RC_ERROR,"ABORT: need 11 or 12 arguments LM_World(if 3D data) IM_World,JM_World, nlon,nlat, in_filename, out_filename, year, month, day, hour, minute")
    endif

    if (nargs==11) then
    CALL GETARG(1, str_arg)
    read (str_arg,'(I10)') IM_World
    CALL GETARG(2, str_arg)
    read (str_arg,'(I10)') JM_World
    CALL GETARG(3, str_arg)
    read (str_arg,'(I10)') nlon
    CALL GETARG(4, str_arg)
    read (str_arg,'(I10)') nlat
    CALL GETARG(5, str_arg)
    in_Filename = str_arg
    CALL GETARG(6, str_arg)
    out_Filename = str_arg
    CALL GETARG(7, str_arg)
    read (str_arg,'(I10)') yr
    CALL GETARG(8, str_arg)
    read (str_arg,'(I10)') mon
    CALL GETARG(9, str_arg)
    read (str_arg,'(I10)') day
    CALL GETARG(10, str_arg)
    read (str_arg,'(I10)') hr
    CALL GETARG(11, str_arg)
    read (str_arg,'(I10)') mn
    else
    CALL GETARG(1, str_arg)
    read (str_arg,'(I10)') LM_World
    CALL GETARG(2, str_arg)
    read (str_arg,'(I10)') IM_World
    CALL GETARG(3, str_arg)
    read (str_arg,'(I10)') JM_World
    CALL GETARG(4, str_arg)
    read (str_arg,'(I10)') nlon
    CALL GETARG(5, str_arg)
    read (str_arg,'(I10)') nlat
    CALL GETARG(6, str_arg)
    in_Filename = str_arg
    CALL GETARG(7, str_arg)
    out_Filename = str_arg
    CALL GETARG(8, str_arg)
    read (str_arg,'(I10)') yr
    CALL GETARG(9, str_arg)
    read (str_arg,'(I10)') mon
    CALL GETARG(10, str_arg)
    read (str_arg,'(I10)') day
    CALL GETARG(11, str_arg)
    read (str_arg,'(I10)') hr
    CALL GETARG(12, str_arg)
    read (str_arg,'(I10)') mn
    endif


!   Initialize the ESMF. For performance reasons, it is important
!    to turn OFF ESMF's automatic logging feature
!   -------------------------------------------------------------
    call ESMF_Initialize (LogKindFlag=ESMF_LOGKIND_NONE, vm=vm, __RC__ )

!   Check the number of processors
!   ------------------------------
    call ESMF_VMGet(vm, localPET=myPET, PETcount=nPET)  
    if ( nPET /= 6 ) then
       __raise__(MAPL_RC_ERROR,"Invalid number of PETS; needs to be 6 PEs")
    end if

    if ( MAPL_am_I_root() ) then
         print *
         print *, 'Starting ' // Iam // ' with ', nPET, ' PETs ...'
         print *
    end if

!   Create a global Lat-Lon grids on a 2x1 layout
!   Uses origin defaults: date-line and South Pole
!   -----------------------------------------------
    __try__

!      3D grid: 
!      ---------------------------------------------------
       grid3d = AppGridCreateF(IM_World, JM_World, LM_World, Nx, Ny, STATUS )

    __except__

        __raise__(ESMF_RC_OBJ_NOT_CREATED,"cannot create grids")

    __endtry__

!   Set the time as the one on the hardwired file name
!   --------------------------------------------------
    call ESMF_CalendarSetDefault ( ESMF_CALKIND_GREGORIAN, __RC__ )
    call ESMF_TimeSet(Time(1), yy=yr, mm=mon, dd=day,  h=hr,  m=mn, s=0, __RC__ )

!   Create empty Bundles
!   --------------------
    Bundle3D = ESMF_FieldBundleCreate ( name='c2l', __RC__ )
    call ESMF_FieldBundleSet ( bundle3D, grid=grid3d, __RC__ )

!                               ---------------
!                               Reading Bundles
!                               ---------------

__try__

!   Read given time on file into a 3D Bundle: all variables
!   -------------------------------------------------------
    call MAPL_CFIORead  ( in_Filename, Time(1), Bundle3d, &
                          verbose=.true., __rc__ )

__except__

    __raise__(ESMF_RC_FILE_READ,"Cannot read bundles")

__endtry__

!                               ---------------
!                               Writing Bundles
!                               ---------------

!   You can write a Bundle to a file at the same or at another resolution.
!   The data will be interpolated horizontally, and possibly thinned
!   vertically.
!   ---------------------------------------------------------------------
    allocate ( resolution(2), __STAT__ )
    resolution = (/ nlon, nlat /)

!   We first create a clock with a initial time on file
!   ---------------------------------------------------
    call ESMF_TimeIntervalSet( TimeStep, h=6, m=0, s=0, __RC__ )
    Clock = ESMF_ClockCreate ( name="c2l", timeStep=TimeStep, &
                               startTime=Time(1), __RC__ )

__try__

!   Before writing to the file, we create a CFIO object 
!   Here we define a file with reduced resolution
!   ---------------------------------------------------
    call MAPL_CFIOCreate ( cfio, out_Filename, clock, Bundle3d,  &
                           resolution=resolution, &
                           descr='Bundle Write Test',            &
                           __rc__)

!   Then we write to it
!   -------------------
    call MAPL_cfioWrite ( cfio, Clock, Bundle3d, verbose = .true., &
                          __rc__ ) 

!   and finally close the file
!   --------------------------
    call MAPL_cfioDestroy ( cfio, __rc__ )

__except__

    __raise__(ESMF_RC_FILE_WRITE,"cannot write bundle")

__endtry__

!   All done
!   --------
    call ESMF_Finalize ( __RC__ )

  end subroutine Main

end Program c2l_CFIO_offline

