
#include "MAPL_Generic.h"

!=============================================================================
!BOP

!  !MODULE: GEOS_WgcmGridCompMod -- A Module to compute wave properties via the
!          wavewatch3 wave model

!  !INTERFACE:

#define RUN_COUPLED

module GEOS_WaveWatchGridCompMod

!  !USES:

    use ESMF
    use MAPL_Mod
    use NUOPC
    use NUOPC_Model, only: label_Advance, label_DataInitialize

!   WW3 modules
    use WMMAPLMD, only : WW3_SetServices => SetServices

    use, intrinsic :: ISO_FORTRAN_ENV

    implicit none
    private


!   Private state
!   -------------

    type WaveModel_State
        private

        type(ESMF_Config) :: CF          ! Private Config
 
        logical:: verbose = .false.      ! verbose messages
 
        real    :: dt     = 0.0          ! wave model time step, s
 
        logical :: stokes = .true.       ! output Stokes drift velocity fields
        real, pointer, dimension(:) :: depths => null()  ! depths for Stokes diagnostics
 
    end type WaveModel_State


!   Hook for the ESMF
!   -----------------
    type WaveModel_Wrap
        type (WaveModel_State), pointer :: ptr => null()
    end type WaveModel_Wrap


!  !PUBLIC MEMBER FUNCTIONS:

    public SetServices

    integer :: WW3GC
!=============================================================================

!  !DESCRIPTION:
! 
!

!EOP

contains

!BOP

! ! IROUTINE: SetServices -- Sets ESMF services for this component

! ! INTERFACE:

    subroutine SetServices(GC, RC)

! ! ARGUMENTS:

        type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
        integer, optional                  :: RC  ! return code

! ! DESCRIPTION: This version uses the MAPL_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!                our instance of a generic state and putting it in the 
!                gridded component (GC). Here we only need to set the run method and
!                add the state variable specifications (also generic) to our instance
!                of the generic state. This is the way our true state variables get into
!                the ESMF_State INTERNAL, which is in the MAPL_MetaComp.

!EOP

!=============================================================================
!
! ErrLog Variables

        character(len=ESMF_MAXSTR)      :: Iam
        integer                         :: STATUS
        character(len=ESMF_MAXSTR)      :: COMP_NAME

! Local derived type aliases
        type(MAPL_MetaComp),  pointer   :: MAPL
        type(ESMF_Config)               :: CF     ! global config

        type (WaveModel_State), pointer :: self   ! private internal state
        type (WaveModel_Wrap)           :: wrap


!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

        Iam = 'SetServices'
        call ESMF_GridCompGet(GC, NAME=COMP_NAME, __RC__)
        Iam = trim(COMP_NAME) // Iam

! Wrap the private internal state for storing in GC
! -------------------------------------------------

        allocate(self, __STAT__)
        wrap%ptr => self
 
! Load private Config Attributes
! ------------------------------

        self%CF = ESMF_ConfigCreate(__RC__)

!       call ESMF_ConfigLoadFile(self%CF, UMWM_CONFIG_FILE, __RC__)

! process the config ....

! add a child component (the NUOPC wrapped ww3)
        ! this internally executes the SetServices method
        ! of the child
        WW3GC = MAPL_AddChild(GC, NAME='WW3', SS=WW3_SetServices, RC=STATUS)
        VERIFY_(STATUS)

! add import, export (and possibly internal) state variables


! Set the Initialize, Run entry point
! -----------------------------------

!ALT        call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_INITIALIZE, Initialize, __RC__)
!ALT        call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN,        Run,        __RC__)
!ALT        call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_FINALIZE,   Finalize,   __RC__)


! Store private internal state in GC
! ----------------------------------
        call ESMF_UserCompSetInternalState(GC, 'WaveModel_State', wrap, STATUS)
        VERIFY_(STATUS)


! Set the state variable specs.
! -----------------------------

! !INTERNAL STATE:


! !IMPORT STATE:

!ALT: we need to zero these 3
!        impFieldName(i)     = 'seahgt'
!        impFieldStdName(i)  = 'sea_surface_height_above_sea_level'
!        impFieldName(i)     = 'uucurr'
!        impFieldStdName(i)  = 'surface_eastward_sea_water_velocity'
!        impFieldName(i)     = 'vvcurr'
!        impFieldStdName(i)  = 'surface_northward_sea_water_velocity'

! the next 3 are needed but have different names: 'uutrue','vvtrue','seaice'
        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'U10M',                                 &
            LONG_NAME      = '10-meter_eastward_wind',               &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            RESTART        = MAPL_RestartSkip,        __RC__) 
  
        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'V10M',                                 &
            LONG_NAME      = '10-meter_northward_wind',              &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            RESTART        = MAPL_RestartSkip,        __RC__)

        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'FRACI',                                &
            LONG_NAME      = 'ice_covered_fraction_of_tile',         &
            UNITS          = '1',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            RESTART        = MAPL_RestartSkip,        __RC__)


!  AGCM -> WGCM


!  OGCM -> WGCM
  


! !EXPORT STATE:
        ! Name change: 'charno', 'z0rlen' 
        call MAPL_AddExportSpec(GC,                                  &
             SHORT_NAME    = 'CHARNOCK',                             &
             LONG_NAME     = 'wave_model_charnock_coefficient',      &
             UNITS         = '1',                                    &
             DIMS          = MAPL_DimsHorzOnly,                      &
             VLOCATION     = MAPL_VLocationNone,     __RC__) 
     
        call MAPL_AddExportSpec(GC,                                  &
             SHORT_NAME    = 'Z0',                                   &
             LONG_NAME     = 'surface_roughness',                    &
             UNITS         = 'm',                                    &
             DIMS          = MAPL_DimsHorzOnly,                      &
             VLOCATION     = MAPL_VLocationNone,     __RC__) 


        call MAPL_AddExportSpec(GC,                                  &
             SHORT_NAME    = 'USTAR',                                &
             CHILD_ID      = WW3GC,                  __RC__)

        call MAPL_AddExportSpec(GC,                                  &
             SHORT_NAME    = 'DCP',                                  &
             CHILD_ID      = WW3GC,                  __RC__)

        call MAPL_AddExportSpec(GC,                                  &
             SHORT_NAME    = 'SWH',                                  &
             CHILD_ID      = WW3GC,                  __RC__)

        call MAPL_AddExportSpec(GC,                                  &
             SHORT_NAME    = 'EDF',                                  &
             CHILD_ID      = WW3GC,                  __RC__)

 
 

! Set the Profiling timers
! ------------------------

        call MAPL_TimerAdd(GC, name='TOTAL'        , __RC__)
        call MAPL_TimerAdd(GC, name='INITIALIZE'   , __RC__)
        call MAPL_TimerAdd(GC, name='RUN'          , __RC__)
        call MAPL_TimerAdd(GC, name='FINALIZE'     , __RC__)



!ALT: we need to terminate child's import so they do not "bubble up". We will fill them explicitly

! this should be irrelevant here because the children are not MAPL components
!        call MAPL_TerminateImport()

! Set generic init and final methods
! ----------------------------------

        call MAPL_GenericSetServices(GC, __RC__)

        RETURN_(ESMF_SUCCESS)

    end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0
- execute all of required NUOPC phases (0,1,3)
- phase 0 : general init
- phase 1: ww3 init, advertise the import and export fields
- phase 3: create grids, realize the fields
- create regridding route handles (import and export)
- regrid imports
#endif

#if 0
!BOP

! ! IROUTINE: INITIALIZE -- Initialize method for the WM component

! !INTERFACE:

   subroutine Initialize(GC, IMPORT, EXPORT, CLOCK, RC)

! !ARGUMENTS:

      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
      integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: 

!EOP


! ErrLog Variables

      character(len=ESMF_MAXSTR)   :: Iam
      integer                      :: STATUS
      character(len=ESMF_MAXSTR)   :: COMP_NAME

! Local derived type aliases

      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_Grid)              :: GRID
      type(ESMF_VM)                :: VM

      type(ESMF_Alarm)             :: run_alarm
      type(ESMF_TimeInterval)      :: ring_interval
      real(ESMF_KIND_R8)           :: time_step

      type (WaveModel_State), pointer :: self   ! private internal state
      type (WaveModel_Wrap)           :: wrap

! Local Variables

      integer :: COMM ! MPI communicator from VM
      integer :: myPE
      integer :: nPEs
      integer :: IM, JM, LM
      integer :: IM_world, JM_world
      integer :: I, J, K, L

      real, pointer, dimension(:,:) :: LATS => NULL()
      real, pointer, dimension(:,:) :: LONS => NULL()

      integer :: COUNTS(ESMF_MAXDIM)
      type (ESMF_GridComp),      pointer  :: GCS(:) => null()
      type (ESMF_State),         pointer  :: GIM(:) => null()
      type (ESMF_State),         pointer  :: GEX(:) => null()

      integer :: phases(3), p
      character(len=ESMF_MAXSTR) :: lbl
      character(len=ESMF_MAXSTR), parameter :: phase_lbl(2) = (/'IPDv03p1', 'IPDv03p3'/)
    
      type(ESMF_Time) :: currTime
      type(ESMF_Field) :: field

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

      Iam = 'Initialize'
      call ESMF_GridCompGet(GC, name=COMP_NAME, __RC__)
      Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
! ---------------------------------

      call MAPL_GetObjectFromGC(GC, MAPL, __RC__)


! Get my internal private state
! -----------------------------
      call ESMF_UserCompGetInternalState(GC, 'WaveModel_State', wrap, STATUS)
      VERIFY_(STATUS)

      self => wrap%ptr


! Get layout from the grid
! ------------------------

      call ESMF_VMGetCurrent(VM, __RC__)

      call ESMF_VMGet(VM, mpiCommunicator=COMM, localPet=myPE, petCount=nPEs, __RC__)

! Start the timers
! ----------------

      call MAPL_TimerOn(MAPL, 'TOTAL',      __RC__)
      call MAPL_TimerOn(MAPL, 'INITIALIZE', __RC__)

! Get parameters from generic state.
! ----------------------------------

      call MAPL_Get(MAPL, GCS=GCS, GIM=GIM, GEX=GEX, RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_Set(MAPL, ChildInit=.false., RC=STATUS)
      VERIFY_(STATUS)

      phases = [0,6,7]
!      DO P=1,7
!         lbl=''
!         k=p
!         call NUOPC_CompSearchRevPhaseMap(GCS(WW3GC), &
!              methodFlag=ESMF_METHOD_INITIALIZE, &
!              phaseIndex=k, &
!              phaseLabel=lbl, &
!              rc=status)
!         VERIFY_(STATUS)
!         ASSERT_(lbl == phase_lbl(p))
!         print *,'WW3GC:phase ',p,k,trim(lbl)
!      END DO

      !ALT at this point PHASES should be [0,6,7]
      DO P = 1,3
         call ESMF_GridCompInitialize(GCS(WW3GC), importState=GIM(WW3GC), &
              exportState=GEX(WW3GC), clock=clock, phase=PHASES(P), &
              userRC=status )
         if (status /= 0) print *,'WW3GC_INIT: status ', status,' phase',P
         VERIFY_(STATUS)
      END DO
      call ESMF_GridCompSet(GCS(WW3GC), clock=clock, __RC__)

! time stamping
      call ESMF_ClockGet(clock, currTime=currTime, rc=status)
      VERIFY_(STATUS)
#if 0
      call ESMF_StateGet(GIM(WW3GC), "uutrue", field, rc=status)
      VERIFY_(STATUS)
      call NUOPC_SetTimestamp(field, currTime, rc=status)
      VERIFY_(STATUS)
      call ESMF_StateGet(GIM(WW3GC), "vvtrue", field, rc=status)
      VERIFY_(STATUS)
      call NUOPC_SetTimestamp(field, currTime, rc=status)
      VERIFY_(STATUS)
      call ESMF_StateGet(GIM(WW3GC), "seaice", field, rc=status)
      VERIFY_(STATUS)
      call NUOPC_SetTimestamp(field, currTime, rc=status)
      VERIFY_(STATUS)
!      print *,'calling DataInitialize'
      call ESMF_MethodExecute (GCS(WW3GC), label=label_DataInitialize, RC=STATUS)
      VERIFY_(STATUS)
#endif
! Get the grid
! ------------
! this might be involved. If the call below does not work, we need to get a
! field for let say import state and then extract the grid

#ifdef RUN_COUPLED
      call ESMF_GridCompGet( GCS(WW3GC), grid=GRID, __RC__)
      ! make this grid my own
      call ESMF_GridCompSet( GC, grid=GRID, __RC__)
#endif


!ALT this is bad!!!
      call MAPL_GridCompSetEntryPoint ( GCS(WW3GC), ESMF_METHOD_WRITERESTART, MAPL_GenericRecord, RC=STATUS)
      VERIFY_(STATUS)

! Call GenericInitialize
! ----------------------
      call MAPL_GenericInitialize(GC, IMPORT, EXPORT, CLOCK, __RC__)


! Get the time step
! -----------------
      call MAPL_Get(MAPL, RunAlarm=run_alarm, __RC__)
      call ESMF_AlarmGet(run_alarm, ringInterval=ring_interval, __RC__)

      call ESMF_TimeIntervalGet(ring_interval, s_r8=time_step, __RC__)
      self%dt = real(time_step)


! Stop the timers
! ---------------

      call MAPL_TimerOff(MAPL, 'INITIALIZE', __RC__)
      call MAPL_TimerOff(MAPL, 'TOTAL',      __RC__)


! All Done
! --------

      RETURN_(ESMF_SUCCESS)

   end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! ! IROUTINE: RUN -- Run method for the WM component

! !INTERFACE:

   subroutine Run(GC, IMPORT, EXPORT, CLOCK, RC)

! !ARGUMENTS:

      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
      integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: 

!EOP

! ErrLog Variables

      character(len=ESMF_MAXSTR)    :: Iam
      integer                       :: STATUS
      integer                       :: iSTAT
      character(len=ESMF_MAXSTR)    :: COMP_NAME


! Local derived type aliases

      type(WaveModel_State), pointer  :: self => null()
      type(WaveModel_Wrap)            :: wrap

      type(MAPL_MetaComp), pointer :: MAPL
      type (ESMF_State)            :: INTERNAL
      type(ESMF_Grid)              :: GRID
      type(ESMF_VM)                :: VM

! Local Variables

      integer :: isd, ied, jsd, jed  ! halo indexes
      integer :: isc, iec, jsc, jec  ! comp indexes
      integer :: hwn
      integer :: hwe
      integer :: hws
      integer :: hww
      integer :: i1, in, j1, jn

! Pointers from Import state
      real, pointer :: U10M(:,:) => null()
      real, pointer :: V10M(:,:) => null()
      real, pointer :: fraci(:,:) => null()

! Pointers to child's Import state

      real, pointer, dimension(:,:) :: uutrue => null()
      real, pointer, dimension(:,:) :: vvtrue => null()
      real, pointer, dimension(:,:) :: seaice => null()

! Pointers to my Export state

      real, pointer, dimension(:,:) :: z0 => null()
      real, pointer, dimension(:,:) :: charnock => null()


! Pointers to child's Export state

      real, pointer, dimension(:,:) :: z0rlen => null()
      real, pointer, dimension(:,:) :: charno => null()

      type (ESMF_GridComp), pointer :: GCS(:)
      type (ESMF_State),    pointer :: GIM(:)
      type (ESMF_State),    pointer :: GEX(:)

      type (ESMF_Clock) :: myclock
      type (ESMF_Time)  :: currTime
      type (ESMF_Field) :: field

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

      Iam = 'Run'

      call ESMF_GridCompGet(GC, name=COMP_NAME, GRID=GRID, __RC__)
      Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
! ---------------------------------

      call MAPL_GetObjectFromGC(GC, MAPL, __RC__)

! Start the timers
! ----------------

      call MAPL_TimerOn(MAPL, 'TOTAL', __RC__)
      call MAPL_TimerOn(MAPL, 'RUN',   __RC__)

! Get parameters from generic state.
! ----------------------------------
      call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_GRID_INTERIOR(GRID,I1,IN,J1,JN)
      i1 = i1-1
      j1 = j1-1

! Get my internal private state
! -----------------------------
      call ESMF_UserCompGetInternalState(GC, 'WaveModel_State', wrap, STATUS)
      VERIFY_(STATUS)

      self => wrap%ptr

! Get pointers to inputs
! ----------------------

#ifdef RUN_COUPLED
      call MAPL_GetPointer(IMPORT, U10M,  'U10M',  __RC__)
      call MAPL_GetPointer(IMPORT, V10M,  'V10M',  __RC__)
      call MAPL_GetPointer(IMPORT, FRACI, 'FRACI', __RC__)

! Get pointers to child's imports
      call MAPL_GetPointer(GIM(WW3GC), UUTRUE, 'uutrue', __RC__)
      call MAPL_GetPointer(GIM(WW3GC), VVTRUE, 'vvtrue', __RC__)
      call MAPL_GetPointer(GIM(WW3GC), SEAICE, 'seaice', __RC__)

      isd = lbound(uutrue,1)
      ied = ubound(uutrue,1)
      jsd = lbound(uutrue,2)
      jed = ubound(uutrue,2)
      ! (isc:iec,jsc:jec) is the "computational" domain.
      ! these should be used _only_ for WW3GC vars

! Copy in. We need to be careful about the haloed vars
      hww = lbound(u10m,1) + i1 - isd
      hwe = ied - ubound(u10m,1) - i1 
      hws = lbound(u10m,2) + j1 - jsd
      hwn = jed - ubound(u10m,2) - j1 

      isc = isd + hwe
      iec = ied - hww
      jsc = jsd + hws
      jec = jed - hwn
 
      uutrue = 0.0
      vvtrue = 0.0
      seaice = 0.0
      uutrue(isc:iec, jsc:jec) = u10m
      vvtrue(isc:iec, jsc:jec) = v10m
      seaice(isc:iec, jsc:jec) = fraci
#endif

! Call WW3
! -----------------
      if (MAPL_AM_I_Root()) write (OUTPUT_UNIT,*) 'DEBUG::WW3GC  WW3 run...'

!@      call ESMF_GridCompRun (GCS(WW3GC), importState=GIM(WW3GC), &
!@           exportState=GEX(WW3GC), clock=CLOCK, userRC=STATUS, RC = istat )
!@      VERIFY_(STATUS)
!      print *,'calling modelAdvance'

      myclock = ESMF_ClockCreate(clock, rc=status) 
      VERIFY_(STATUS)
      call ESMF_ClockAdvance(myclock, rc=status)
      VERIFY_(STATUS)

      call NUOPC_SetTimestamp(GEX(WW3GC), myclock, rc=status)
      VERIFY_(STATUS)

!ALT: i am not sure if we need to time stamp import state
      call ESMF_ClockGet(myclock, currTime=currTime, rc=status)
      VERIFY_(STATUS)
      call ESMF_StateGet(GIM(WW3GC), "uutrue", field, rc=status)
      VERIFY_(STATUS)
      call NUOPC_SetTimestamp(field, currTime, rc=status)
      VERIFY_(STATUS)
      call ESMF_StateGet(GIM(WW3GC), "vvtrue", field, rc=status)
      VERIFY_(STATUS)
      call NUOPC_SetTimestamp(field, currTime, rc=status)
      VERIFY_(STATUS)
      call ESMF_StateGet(GIM(WW3GC), "seaice", field, rc=status)
      VERIFY_(STATUS)
      call NUOPC_SetTimestamp(field, currTime, rc=status)
      VERIFY_(STATUS)

      call ESMF_MethodExecute (GCS(WW3GC), label=label_Advance, RC=STATUS)
      VERIFY_(STATUS)

! Get pointers from export state
! ------------------------------
#ifdef ENABLE_EXPORTS 
      call MAPL_GetPointer(EXPORT, Z0,        'Z0',       alloc=.true., __RC__)
      call MAPL_GetPointer(EXPORT, CHARNOCK,  'CHARNOCK', alloc=.true., __RC__)
      call MAPL_GetPointer(GEX(WW3GC), z0rlen, 'z0rlen', __RC__)
      call MAPL_GetPointer(GEX(WW3GC), charno, 'charno', __RC__)

      ! Copy out. Again take only "computational" domain
      if(associated(z0)) z0 = z0rlen(isc:iec, jsc:jec)
      if(associated(charnock)) charnock = charno(isc:iec, jsc:jec)
#endif


! Stop the timers
! ---------------

      call MAPL_TimerOff(MAPL, 'RUN',   __RC__)
      call MAPL_TimerOff(MAPL, 'TOTAL', __RC__)

! All Done
! --------

      RETURN_(ESMF_SUCCESS)

   end subroutine Run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! ! IROUTINE: FINALIZE -- Finalize method for the WM component

! !INTERFACE:

   subroutine Finalize(GC, IMPORT, EXPORT, CLOCK, RC)

! !ARGUMENTS:

      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
      integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: 

!EOP

! ErrLog Variables

      character(len=ESMF_MAXSTR)    :: Iam
      integer                       :: STATUS
      character(len=ESMF_MAXSTR)    :: COMP_NAME

! Local derived type aliases

      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_Grid)              :: GRID
      type(ESMF_VM)                :: VM

! Local global variables

      integer :: COUNTS(ESMF_MAXDIM)

! Local Variables

      integer               :: IM, JM, LM
      integer               :: IM_world, JM_world
      integer               :: COMM ! MPI communicator from VM
      integer               :: myPE
      integer               :: nPEs

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

      Iam = 'Finalize'
      call ESMF_GridCompGet(GC, name=COMP_NAME, GRID=GRID, __RC__)
      Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
! ---------------------------------

      call MAPL_GetObjectFromGC(GC, MAPL, __RC__)


! Start the timers
! ----------------

      call MAPL_TimerOn(MAPL, 'TOTAL',    __RC__)
      call MAPL_TimerOn(MAPL, 'FINALIZE', __RC__)


! Get parameters from generic state.
! ----------------------------------
!     None


! Stop the timers
! ---------------

      call MAPL_TimerOff(MAPL, 'FINALIZE', __RC__)
      call MAPL_TimerOff(MAPL, 'TOTAL',    __RC__)


! Call GenericFinalize
! ----------------------
      call MAPL_GenericFinalize(GC, IMPORT, EXPORT, CLOCK, __RC__)
      VERIFY_(STATUS)


! All Done
! --------

      RETURN_(ESMF_SUCCESS)

   end subroutine Finalize

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_WaveWatchGridCompMod

