 !------------------------------------------------------------------------------
#include "MAPL_Generic.h"
!
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: AdvCore_GridCompMod
!
! !DESCRIPTION: 
!    This a MAPL component that can be used in
!    either with offline or online applications to advect an arbitrary set
!    of constituents.
!
! \paragraph{Scientific Description:}
!
!   The advection scheme used is that from the FVdycore grid-point
!   dynamical core.  It runs on a sphere and uses finite-volume
!   discretization techniques. The advection is time split into a
!   horizontal phase that is assumed to be vertically Lagrangian and a
!   vertical remap phase. A complete description of the core from
!   which this component is taken may be found in:
!
!   \begin{quote}
!   Lin, S.-J. 2004, A vertically Lagrangian Finite-Volume Dynamical 
!   Core for Global Models. {\em Mon. Wea. Rev.}, {\bf 132}, 2293-2307.
!   \end{quote}
!
!  \paragraph{Code Implementation:}
!
!    It code uses the MAPL (http://MAPLCode.org/maplwiki/) to
!    encapsulate the FV advection scheme as an ESMF gridded component
!    using the ESMF paradigm of initialize, run and finalize methods,
!    and their SetServices routine. As in all ESMF codes, only
!    SetServices is public and the interface consists of of a Clock
!    and Import and Export states.  The import state includes a
!    specialized description of the motion field in terms of C-grid
!    winds and mass fluxes. These are assumed to have been accumulated
!    over the time interval specified in the resource file. The
!    default of this interval is 1800 seconds. The layer pressure
!    thicknesses in the import state are assumed to be the
!    instantaneous values valid at the beginning of this interval.  If
!    these thicknesses are friendly they will be updated to values
!    valid at the end of the interval, consistent with the given
!    motion field.  Mixing ratios of the constituents to be advected
!    are placed ESMF Fields within an ESMF Bundle in the Import
!    state. Each Field in the Bundle is tested for ``Friendliness'' to
!    advection; if friendly it is advected and its values updated.
!
!    Currently no Export capability is implemented. 
!
! !INTERFACE:

module AdvCore_GridCompMod

!
! !USES:

      use ESMF
      use MAPL_Mod
      use m_set_eta,       only: set_eta
      use fv_arrays_mod,   only: fv_atmos_type, FVPRC, REAL4, REAL8
      use fms_mod,         only: fms_init, set_domain, nullify_domain
      use fv_control_mod,  only: fv_init1, fv_init2, fv_end
      use fv_tracer2d_mod, only: offline_tracer_advection
      use fv_mp_mod,       only: is,ie, js,je, is_master, tile
      use fv_grid_utils_mod, only: g_sum

      USE FV_StateMod,     only: AdvCoreTracers => T_TRACERS
      USE FV_StateMod,     only: FV_Atm
      use CubeGridPrototype, only: register_grid_and_regridders

      implicit none
      private

      integer     :: nx, ny
      integer     :: npes_x, npes_y
      real(FVPRC) :: dt
      logical     :: FV3_DynCoreIsRunning=.false.
      integer     :: AdvCore_Advection=1
      logical     :: chk_mass=.false.

      integer,  parameter :: ntiles_per_pe = 1

! Tracer I/O History stuff
! -------------------------------------
      integer, parameter         :: ntracers=11
      integer                    :: ntracer
      character(len=ESMF_MAXSTR) :: myTracer
      character(len=ESMF_MAXSTR) :: tMassStr
      real(FVPRC), SAVE          :: TMASS0(ntracers)
      real(REAL8), SAVE          ::  MASS0
      logical    , SAVE          :: firstRun=.true.

! !PUBLIC MEMBER FUNCTIONS:

      public SetServices
      logical, allocatable, save :: grids_on_my_pe(:)

!EOP

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: SetServices - Externally visible registration routine
!
! !INTERFACE:
!
      subroutine SetServices(GC, rc)
!
! !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: GC
      integer, optional,   intent(  out) :: RC
!
! !DESCRIPTION:
!
!     User-supplied setservices routine.
!     The register routine sets the subroutines to be called
!     as the init, run, and finalize routines.  Note that those are
!     private to the module.
!
!EOP

      character(len=ESMF_MAXSTR)              :: IAm
      integer                                 :: STATUS
      character(len=ESMF_MAXSTR)              :: COMP_NAME
      type (MAPL_MetaComp),      pointer      :: MAPL
      character(len=ESMF_MAXSTR)              :: DYCORE
      type(ESMF_VM)                           :: VM
      integer                                 :: comm, ndt
      integer                                 :: p_split=1

!=============================================================================

! Begin...

      ! Get my name and set-up traceback handle
      ! ---------------------------------------
    
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, vm=vm, RC=STATUS )
      _VERIFY(STATUS)
      Iam = trim(COMP_NAME) // 'SetServices'

!BOS

! !IMPORT STATE:
!
    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'MFX',                                       &
         LONG_NAME  = 'pressure_weighted_eastward_mass_flux',      &
         UNITS      = 'Pa m+2 s-1',                                &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     _VERIFY(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'MFY',                                       &
         LONG_NAME  = 'pressure_weighted_northward_mass_flux',     &
         UNITS      = 'Pa m+2 s-1',                                &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     _VERIFY(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'CX',                                        &
         LONG_NAME  = 'eastward_accumulated_courant_number',       &
         UNITS      = '',                                          &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     _VERIFY(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'CY',                                        &
         LONG_NAME  = 'northward_accumulated_courant_number',      &
         UNITS      = '',                                          &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     _VERIFY(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'PLE0',                                      &
         LONG_NAME  = 'pressure_at_layer_edges_before_advection',  &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
     _VERIFY(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'PLE1',                                      &
         LONG_NAME  = 'pressure_at_layer_edges_after_advection',   &                
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
     _VERIFY(STATUS)

    call MAPL_AddImportSpec( gc,                              &
        SHORT_NAME = 'TRADV',                                        &
        LONG_NAME  = 'advected_quantities',                        &
        UNITS      = 'unknown',                                    &
        DATATYPE   = MAPL_BundleItem,               &
        RC=STATUS  )
    _VERIFY(STATUS)

  !EXPORT STATE:
     call MAPL_AddExportSpec ( gc,                                  &
          SHORT_NAME = 'AREA',                                      &
          LONG_NAME  = 'agrid_cell_area',                           &
          UNITS      = 'm+2'  ,                                     &
          DIMS       = MAPL_DimsHorzOnly,                           &
          VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     _VERIFY(STATUS)
 

! 3D Tracers
     do ntracer=1,ntracers
        write(myTracer, "('TEST_TRACER',i1.1)") ntracer-1
        call MAPL_AddExportSpec ( gc,                             &
             SHORT_NAME = TRIM(myTracer),                         &
             LONG_NAME  = TRIM(myTracer),                         &
             UNITS      = '1',                                    &
             DIMS       = MAPL_DimsHorzVert,                      &
             VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
        _VERIFY(STATUS)
     enddo

!EOS

      ! Set the Profiling timers
      !-------------------------
      call MAPL_TimerAdd(GC,    name="INITIALIZE"  ,RC=STATUS)
      _VERIFY(STATUS)
      call MAPL_TimerAdd(GC,    name="RUN"         ,RC=STATUS)
      _VERIFY(STATUS)
      call MAPL_TimerAdd(GC,    name="FINALIZE"    ,RC=STATUS)
      _VERIFY(STATUS)
      call MAPL_TimerAdd(GC,    name="TOTAL"       ,RC=STATUS)
      _VERIFY(STATUS)


! Register methods with MAPL
! --------------------------

      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize, RC=status )
      _VERIFY(STATUS)
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,          Run,       RC=status )
      _VERIFY(STATUS)
      call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_FINALIZE,     Finalize,  RC=status)
      _VERIFY(STATUS)

      ! Check if AdvCore is running without FV3_DynCoreIsRunning, if yes then setup the MAPL Grid 
      ! ----------------------------------------------------------------------------
      call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
      _VERIFY(STATUS)
      call MAPL_GetResource(MAPL, DYCORE, 'DYCORE:', default="", RC=STATUS )
      _VERIFY(STATUS)
      call MAPL_GetResource(MAPL, AdvCore_Advection , label='AdvCore_Advection:', &
                                  default=AdvCore_Advection, RC=STATUS )
      _VERIFY(STATUS)
      if(adjustl(DYCORE)=="FV3") FV3_DynCoreIsRunning = .true.
      if(adjustl(DYCORE)=="FV3+ADV") FV3_DynCoreIsRunning = .true.

      ! Start up FMS/MPP
      !-------------------------------------------
      call ESMF_VMGet(VM,mpiCommunicator=comm,rc=STATUS)
      _VERIFY(STATUS)
      call fms_init(comm)
      _VERIFY(STATUS)

      if (.NOT. FV3_DynCoreIsRunning) then
      ! Make sure FV3 is setup
      ! -----------------------
         call register_grid_and_regridders()
         call fv_init1(FV_Atm, dt, grids_on_my_pe, p_split)
      ! Get Domain decomposition
      !-------------------------
         call MAPL_GetResource( MAPL, nx, 'NX:', default=0, RC=STATUS )
         _VERIFY(STATUS)
         FV_Atm(1)%layout(1) = nx
         call MAPL_GetResource( MAPL, ny, 'NY:', default=0, RC=STATUS )
         _VERIFY(STATUS)
         if (FV_Atm(1)%flagstruct%grid_type == 4) then
            FV_Atm(1)%layout(2) = ny
         else
            FV_Atm(1)%layout(2) = ny / 6
         end if
      ! Get Resolution Information
      !---------------------------
      ! FV grid dimensions setup from MAPL
         call MAPL_GetResource( MAPL, FV_Atm(1)%flagstruct%npx, 'IM:', default= 32, RC=STATUS )
         _VERIFY(STATUS)
         call MAPL_GetResource( MAPL, FV_Atm(1)%flagstruct%npy, 'JM:', default=192, RC=STATUS )
         _VERIFY(STATUS)
         call MAPL_GetResource( MAPL, FV_Atm(1)%flagstruct%npz, 'LM:', default= 72, RC=STATUS )
         _VERIFY(STATUS)
      ! FV likes npx;npy in terms of cell vertices
         if (FV_Atm(1)%flagstruct%npy == 6*FV_Atm(1)%flagstruct%npx) then
            FV_Atm(1)%flagstruct%ntiles = 6
            FV_Atm(1)%flagstruct%npy    = FV_Atm(1)%flagstruct%npx+1
            FV_Atm(1)%flagstruct%npx    = FV_Atm(1)%flagstruct%npx+1
         else
            FV_Atm(1)%flagstruct%ntiles = 1
            FV_Atm(1)%flagstruct%npy    = FV_Atm(1)%flagstruct%npy+1
            FV_Atm(1)%flagstruct%npx    = FV_Atm(1)%flagstruct%npx+1
         endif
      endif

      call MAPL_GetResource( MAPL, ndt, 'RUN_DT:', default=0, RC=STATUS )
      _VERIFY(STATUS)
      DT = ndt

      ! Start up FV if AdvCore is running without FV3_DynCoreIsRunning
      !--------------------------------------------------
      if (.NOT. FV3_DynCoreIsRunning) then
         call fv_init2(FV_Atm, dt, grids_on_my_pe, p_split)
      end if

      ! Ending with a Generic SetServices call is a MAPL requirement 
      !-------------------------------------------------------------
      call MAPL_GenericSetServices    ( GC, rc=STATUS)
      _VERIFY(STATUS)

      _RETURN(ESMF_SUCCESS)

      end subroutine SetServices

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize - initialization routine
!
! !INTERFACE:
!
  subroutine Initialize(GC, IMPORT, EXPORT, CLOCK, RC)
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION:
!     This initialization routine creates the import and export states,
!     as well as the internal state, which is attached to the component.
!     It also determines the distribution (and therefore the grid) 
!     and performs allocations of persistent data, 
!
!EOP
!=============================================================================
!BOC

      character(len=ESMF_MAXSTR)         :: IAm
      integer                            :: STATUS
      character(len=ESMF_MAXSTR)         :: COMP_NAME
      type(ESMF_Config)                  :: CF
      type (MAPL_MetaComp),      pointer :: MAPL
      type (ESMF_VM)                     :: VM
      real, pointer                      :: temp2d(:,:)
      integer                            :: IS, IE, JS, JE
      logical                            :: gridCreated
      type(ESMF_Grid)                    :: grid

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

      Iam = "Initialize"
      call ESMF_GridCompGet ( GC, name=COMP_NAME, config=CF, vm=VM, RC=STATUS )
      _VERIFY(STATUS)
      Iam = trim(COMP_NAME) // trim(Iam)

      ! Retrieve the pointer to the state
      ! ---------------------------------
      call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
      _VERIFY(STATUS)

      call MAPL_TimerOn(MAPL,"TOTAL")
      call MAPL_TimerOn(MAPL,"INITIALIZE")

      gridCreated=.false.
      call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
      _VERIFY(STATUS)
      call ESMF_GridCompGet(GC,grid=grid,rc=status)
      if (status == ESMF_SUCCESS) then
         call ESMF_GridValidate(grid,rc=status)
         if (status==ESMF_SUCCESS) GridCreated = .true.
      end if

      if (.not.GridCreated) then
         call MAPL_GridCreate(GC, rc=status)
         _VERIFY(STATUS)
      endif

      call MAPL_GenericInitialize(GC, IMPORT, EXPORT, CLOCK, RC=STATUS)
      _VERIFY(STATUS)
      ! Compute Grid-Cell Area
      ! ----------------------
      if (.NOT. FV3_DynCoreIsRunning) then
         IS = FV_Atm(1)%bd%isc
         IE = FV_Atm(1)%bd%iec
         JS = FV_Atm(1)%bd%jsc
         JE = FV_Atm(1)%bd%jec
         call MAPL_GetPointer(EXPORT, temp2d, 'AREA', ALLOC=.TRUE., rc=status)
         _VERIFY(STATUS)
         temp2d = FV_Atm(1)%gridstruct%area(IS:IE,JS:JE)
      endif

      call MAPL_TimerOff(MAPL,"INITIALIZE")
      call MAPL_TimerOff(MAPL,"TOTAL")

      _RETURN(ESMF_SUCCESS)

      end subroutine Initialize
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run - run routine
!
! !INTERFACE:
!
      subroutine Run(GC, IMPORT, EXPORT, CLOCK, RC)
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION:
! 
! The Run method advanced the advection one long time step, as
! specified in the configuration.  This may be broken down int a
! number of internal, small steps, also configurable.
!
!EOP
!=============================================================================
!BOC
! !LOCAL VARIABLES:
      character(len=ESMF_MAXSTR)    :: IAm
      integer                       :: STATUS
      character(len=ESMF_MAXSTR)    :: COMP_NAME
      type (ESMF_Grid)              :: ESMFGRID
      type (MAPL_MetaComp), pointer :: MAPL
      type (ESMF_Alarm)             :: ALARM

! Imports
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iCX
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iCY
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iMFX
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iMFY
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iPLE0
      REAL(REAL8), POINTER, DIMENSION(:,:,:)   :: iPLE1

! Locals
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: CX
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: CY
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: MFX
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: MFY
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: PLE0
      REAL(FVPRC), POINTER, DIMENSION(:,:,:)   :: PLE1
      REAL(FVPRC), POINTER, DIMENSION(:)       :: AK
      REAL(FVPRC), POINTER, DIMENSION(:)       :: BK
      REAL(REAL8), allocatable :: ak_r8(:),bk_r8(:)
      REAL(FVPRC), POINTER, DIMENSION(:,:,:,:) :: TRACERS
      REAL(FVPRC) :: MASS1, TMASS1(ntracers)
      TYPE(AdvCoreTracers), POINTER :: advTracers(:)
      type(ESMF_FieldBundle) :: TRADV
      type(ESMF_Field)       :: field
      type(ESMF_Array)       :: array
      INTEGER :: IM, JM, LM, N, NQ, LS
      REAL(FVPRC) :: PTOP, PINT
      REAL(REAL8) :: ptop_r8,pint_r8
! Temporaries for exports/tracers
      REAL, POINTER :: temp3D(:,:,:)
      real(REAL4),        pointer     :: tracer_r4 (:,:,:)
      real(REAL8),        pointer     :: tracer_r8 (:,:,:)
      character(len=ESMF_MAXSTR)    :: fieldName
      type(ESMF_TypeKind_Flag)      :: kind
      character(len=ESMF_MAXSTR)    :: STRING
! Excluding tracers
      type(ESMF_FieldBundle), save        :: bundleAdv
      type (ESMF_Config)                  :: CF
      logical                             :: adjustTracers
      type(ESMF_Alarm)                    :: predictorAlarm
      type(ESMF_Grid)                     :: bgrid
      integer, save                       :: nq_saved = 0
      integer                             :: i,j
      integer                             :: nqt
      logical                             :: tend
      logical                             :: exclude
      character(len=ESMF_MAXSTR)          :: tmpstring
      character(len=ESMF_MAXSTR)          :: adjustTracerMode
      character(len=ESMF_MAXSTR), allocatable :: xlist(:)
      character(len=ESMF_MAXSTR), allocatable :: biggerlist(:)
      integer, parameter                  :: XLIST_MAX = 60

! Get my name and set-up traceback handle
! ---------------------------------------

      Iam = 'Run'
      call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, grid=ESMFGRID, RC=STATUS )
      _VERIFY(STATUS)
      Iam = trim(COMP_NAME) // Iam

!WMP  if (AdvCore_Advection>0) then

! Get parameters from generic state.
!-----------------------------------
      call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
      _VERIFY(STATUS)
      call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,   &
                                RUNALARM = ALARM, &
                                      RC = STATUS )
      _VERIFY(STATUS)

      call MAPL_TimerOn(MAPL,"TOTAL")
      call MAPL_TimerOn(MAPL,"RUN")

! Get AKs and BKs for vertical grid
!----------------------------------
      AllOCATE( AK(LM+1) ,stat=STATUS )
      _VERIFY(STATUS)
      AllOCATE( BK(LM+1) ,stat=STATUS )
      _VERIFY(STATUS)
      AllOCATE( AK_r8(LM+1) ,stat=STATUS )
      _VERIFY(STATUS)
      AllOCATE( BK_r8(LM+1) ,stat=STATUS )
      _VERIFY(STATUS)
      call set_eta(LM,LS,ptop_r8,pint_r8,ak_r8,bk_r8)
      ptop=ptop_r8
      pint=pint_r8
      ak=ak_r8
      bk=bk_r8

      CALL MAPL_GetPointer(IMPORT, iPLE0, 'PLE0', ALLOC = .TRUE., RC=STATUS)
      _VERIFY(STATUS)
      CALL MAPL_GetPointer(IMPORT, iPLE1, 'PLE1', ALLOC = .TRUE., RC=STATUS)
      _VERIFY(STATUS)
      CALL MAPL_GetPointer(IMPORT, iMFX,   'MFX', ALLOC = .TRUE., RC=STATUS)
      _VERIFY(STATUS)
      CALL MAPL_GetPointer(IMPORT, iMFY,   'MFY', ALLOC = .TRUE., RC=STATUS)
      _VERIFY(STATUS)
      CALL MAPL_GetPointer(IMPORT, iCX,     'CX', ALLOC = .TRUE., RC=STATUS)
      _VERIFY(STATUS)
      CALL MAPL_GetPointer(IMPORT, iCY,     'CY', ALLOC = .TRUE., RC=STATUS)
      _VERIFY(STATUS)

      ALLOCATE( PLE0(IM,JM,LM+1) )
      ALLOCATE( PLE1(IM,JM,LM+1) )
      ALLOCATE(  MFX(IM,JM,LM  ) )
      ALLOCATE(  MFY(IM,JM,LM  ) )
      ALLOCATE(   CX(IM,JM,LM  ) )
      ALLOCATE(   CY(IM,JM,LM  ) )

      PLE0 = iPLE0
      PLE1 = iPLE1
       MFX = iMFX
       MFY = iMFY
        CX = iCX
        CY = iCY

      ! The quantities to be advected come as friendlies in a bundle
      !  in the import state.
      !--------------------------------------------------------------

      call ESMF_StateGet(IMPORT, "TRADV", TRADV, rc=STATUS)
      _VERIFY(STATUS)

      !-------------------------------------------------------------------
      ! ALT: this section attempts to limit the amount of advected tracers
      !-------------------------------------------------------------------
      adjustTracers = .false.
      call MAPL_GetResource ( MAPL, adjustTracerMode, &
           'EXCLUDE_ADVECTION_TRACERS:', &
           default='ALWAYS', rc=status )
      _VERIFY(STATUS)
      if (adjustTracerMode == 'ALWAYS') then
         adjustTracers = .true.
      else if (adjustTracerMode == 'PREDICTOR') then
         !get PredictorAlarm from clock
         call ESMF_ClockGetAlarm(clock, alarmName='PredictorAlarm', &
              alarm=PredictorAlarm, rc=status)
         if (status == ESMF_SUCCESS) then
            !check if ringing
            if (ESMF_AlarmIsRinging(predictorAlarm)) then
               adjustTracers = .true.
            end if
         end if
      else
         call WRITE_PARALLEL('Invalid option, ignored')
         adjustTracers = .false.
      end if
      if (adjustTracers) then
         if (firstRun) then
            ! get the list of excluded tracers from resource
            n = 0
            call ESMF_ConfigFindLabel ( CF,'EXCLUDE_ADVECTION_TRACERS_LIST:',rc=STATUS )
            if(STATUS==ESMF_SUCCESS) then

               tend  = .false.
               allocate(xlist(XLIST_MAX), stat=status)
               _VERIFY(STATUS)
               do while (.not.tend)
                  call ESMF_ConfigGetAttribute (CF,value=tmpstring,default='',rc=STATUS) !ALT: we don't check return status!!!
                  if (tmpstring /= '')  then
                     n = n + 1
                     if (n > size(xlist)) then
                        allocate( biggerlist(2*n), stat=status )
                        _VERIFY(STATUS)
                        biggerlist(1:n-1)=xlist
                        call move_alloc(from=biggerlist, to=xlist)
                     end if
                     xlist(n) = tmpstring
                  end if
                  call ESMF_ConfigNextLine(CF,tableEnd=tend,rc=STATUS )
                  _VERIFY(STATUS)
               enddo
            end if

            ! Count the number of tracers
            !---------------------
            call ESMF_FieldBundleGet(TRADV, grid=bgrid,fieldCount=nqt,  RC=STATUS)
            _VERIFY(STATUS)
            BundleAdv = ESMF_FieldBundleCreate ( name='xTRADV', rc=STATUS )
            _VERIFY(STATUS)
            call ESMF_FieldBundleSet ( BundleAdv, grid=bgrid, rc=STATUS )
            _VERIFY(STATUS)
            !loop over NQ in TRADV
            do i = 1, nqt
               !get field from TRADV and its name
               call ESMF_FieldBundleGet(TRADV, fieldIndex=i, field=field, rc=status)
               _VERIFY(STATUS)
               call ESMF_FieldGet(FIELD, name=fieldname, RC=STATUS)
               _VERIFY(STATUS)
               !exclude everything that is not cloud/water species
               if ( (FV3_DynCoreIsRunning      ) .and. &
                   ( (TRIM(fieldname) == 'Q'       ) .or. &
                     (TRIM(fieldname) == 'QLCN'    ) .or. &
                     (TRIM(fieldname) == 'QLLS'    ) .or. &
                     (TRIM(fieldname) == 'QICN'    ) .or. &
                     (TRIM(fieldname) == 'QILS'    ) .or. &
                     (TRIM(fieldname) == 'CLCN'    ) .or. &
                     (TRIM(fieldname) == 'CLLS'    ) .or. &
                     (TRIM(fieldname) == 'NCPL'    ) .or. &
                     (TRIM(fieldname) == 'NCPI'    ) .or. &
                     (TRIM(fieldname) == 'QRAIN'   ) .or. &
                     (TRIM(fieldname) == 'QSNOW'   ) .or. &
                     (TRIM(fieldname) == 'QGRAUPEL') ) ) then
                   ! write(STRING,'(A,A)') "ADV is excluding ", TRIM(fieldname)
                   ! call WRITE_PARALLEL( trim(STRING)   )
                     n = n + 1
                     if (n > size(xlist)) then
                        allocate( biggerlist(2*n), stat=status )
                        _VERIFY(STATUS)
                        biggerlist(1:n-1)=xlist
                        call move_alloc(from=biggerlist, to=xlist)
                     end if
                     xlist(n) = TRIM(fieldname)
               end if
               !loop over exclude_list
               exclude = .false.
               do j = 1, n
                  if (fieldname == xlist(j)) then
                     exclude = .true.
                     exit
                  end if
               end do
               if (.not. exclude) then
                  call MAPL_FieldBundleAdd(BundleAdv, FIELD, RC=STATUS)
                  _VERIFY(STATUS)
               end if
            end do

            if (allocated(xlist)) then
           !   ! Just in case xlist was allocated, but nothing was in it, could have garbage
           !   if (n > 0) then
           !      call ESMF_FieldBundleRemove(TRADV, fieldNameList=xlist, &
           !         relaxedFlag=.true., rc=status)
           !      _VERIFY(STATUS)
           !   end if
               deallocate(xlist)
            end if

         end if ! firstRun
         TRADV = bundleAdv
      end if ! adjustTracers

      call ESMF_FieldBundleGet(TRADV, fieldCount=NQ,    rc=STATUS)
      _VERIFY(STATUS)

      if (NQ > 0) then
         ! We allocate a list of tracers big enough to hold all items in the bundle
         !-------------------------------------------------------------------------
         ALLOCATE( TRACERS(IM,JM,LM,NQ),stat=STATUS )
         _VERIFY(STATUS)
         ALLOCATE( advTracers(NQ),stat=STATUS )
         _VERIFY(STATUS)

         if (NQ /= NQ_SAVED) then
            write(STRING,'(A,I5,A)') "AdvCore is Advecting the following ", nq, " tracers in FV3:"
            call WRITE_PARALLEL( trim(STRING)   )
         end if

         ! Go through the bundle copying the friendlies into the tracer list.
         !-------------------------------------------------------------------------
         do N=1,NQ
            call ESMF_FieldBundleGet (TRADV, fieldIndex=N, field=FIELD, RC=STATUS)
            _VERIFY(STATUS)
            call ESMF_FieldGet  (field, array=array, name=fieldName, RC=STATUS)
            _VERIFY(STATUS)
            call ESMF_ArrayGet(array,typekind=kind, rc=status )
            _VERIFY(STATUS)
            advTracers(N)%is_r4 = (kind == ESMF_TYPEKIND_R4)   ! Is real*4?
            advTracers(N)%tName = fieldName

            if (NQ /= NQ_SAVED) then
               call WRITE_PARALLEL( trim(fieldName) )
            endif

            if (advTracers(N)%is_r4) then
               call ESMF_ArrayGet(array,farrayptr=tracer_r4, rc=status )
               _VERIFY(STATUS)
               advTracers(N)%content_r4 => tracer_r4
               TRACERS(:,:,:,N) = advTracers(N)%content_r4
            else
               call ESMF_ArrayGet(array,farrayptr=tracer_r8, rc=status )
               _VERIFY(STATUS)
               advTracers(N)%content => tracer_r8
               TRACERS(:,:,:,N) = advTracers(N)%content
            end if
         end do

         if (NQ /= NQ_SAVED) then
            NQ_SAVED = NQ
         end if

         if (chk_mass) then
        ! Check Mass conservation
            if (firstRun .and. AdvCore_Advection>0) then
               MASS0 = g_sum(FV_Atm(1)%domain, PLE0(:,:,LM), is,ie, js,je, FV_Atm(1)%ng, FV_Atm(1)%gridstruct%area_64, 1, .true.)
               call global_integral(TMASS0, TRACERS, PLE0, IM,JM,LM,NQ)
               if (MASS0 /= 0.0) TMASS0=TMASS0/MASS0
            elseif (firstRun) then
               MASS0 = g_sum(FV_Atm(1)%domain, PLE1(:,:,LM), is,ie, js,je, FV_Atm(1)%ng, FV_Atm(1)%gridstruct%area_64, 1, .true.)
               call global_integral(TMASS0, TRACERS, PLE1, IM,JM,LM,NQ)
               if (MASS0 /= 0.0) TMASS0=TMASS0/MASS0
            endif
         endif
         firstRun=.false.

         ! Run FV3 advection
         !------------------
         if (AdvCore_Advection>0) then
         call offline_tracer_advection(TRACERS, PLE0, PLE1, MFX, MFY, CX, CY, &
                                       FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, FV_Atm(1)%bd, &
                                       FV_Atm(1)%domain, AK, BK, PTOP, FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz,   &
                                       NQ, dt)
         endif

         ! Update tracer mass conservation
         !-------------------------------------------------------------------------
         if (chk_mass) then 
            MASS1 = g_sum(FV_Atm(1)%domain, PLE1(:,:,LM), is,ie, js,je, FV_Atm(1)%ng, FV_Atm(1)%gridstruct%area_64, 1, .true.)
            call global_integral(TMASS1, TRACERS, PLE1, IM,JM,LM,NQ)
            if (MASS1 /= 0.0) TMASS1=TMASS1/MASS1
         endif

         if (chk_mass .and. is_master()) then
#ifdef PRINT_MASS
            write(6,100)  MASS0   , &
                         TMASS0(2), &
                         TMASS0(3), &
                         TMASS0(4), &
                         TMASS0(5)
            write(6,102)  MASS1   , &
                         TMASS1(2), &
                         TMASS1(3), &
                         TMASS1(4), &
                         TMASS1(5)
#endif
            write(6,103) ( MASS1   - MASS0   )/ MASS0   , &
                         (TMASS1(2)-TMASS0(2))/TMASS0(2), &
                         (TMASS1(3)-TMASS0(3))/TMASS0(3), &
                         (TMASS1(4)-TMASS0(4))/TMASS0(4), &
                         (TMASS1(5)-TMASS0(5))/TMASS0(5)
 100        format('Tracer M0  : ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14)
 101        format('Tracer Ma  : ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14)
 102        format('Tracer M1  : ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14)
 103        format('Tracer Mdif: ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14,' ',e21.14)
         endif

         ! Go through the bundle copying tracers back to the bundle.
         !-------------------------------------------------------------------------
         do N=1,NQ
            if (advTracers(N)%is_r4) then
               advTracers(N)%content_r4 = TRACERS(:,:,:,N)
            else
               advTracers(N)%content    = TRACERS(:,:,:,N)
            end if
! Fill Export States
            write(myTracer, "('TEST_TRACER',i1.1)") N-1
            call MAPL_GetPointer(EXPORT, temp3D, TRIM(myTracer), rc=status)
            _VERIFY(STATUS)
            if ((associated(temp3D)) .and. (N<=ntracers)) then
               temp3D = TRACERS(:,:,:,N)
            endif
         enddo

         ! Deallocate the list of tracers
         !-------------------------------------------------------------------------
         DEALLOCATE( TRACERS,stat=STATUS )
         _VERIFY(STATUS)

      end if ! NQ > 0

      deallocate( advTracers, stat=STATUS )
      _VERIFY(STATUS)
      DEALLOCATE( AK ,stat=STATUS )
      _VERIFY(STATUS)
      DEALLOCATE( BK ,stat=STATUS )
      _VERIFY(STATUS)

      DEALLOCATE( PLE0 )
      DEALLOCATE( PLE1 )
      DEALLOCATE(  MFX )
      DEALLOCATE(  MFY )
      DEALLOCATE(   CX )
      DEALLOCATE(   CY )

      call MAPL_TimerOff(MAPL,"RUN")
      call MAPL_TimerOff(MAPL,"TOTAL")

!WMP  end if ! AdvCore_Advection

      _RETURN(ESMF_SUCCESS)

      end subroutine Run
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize - user supplied finalize routine
!
! !INTERFACE:
!
  subroutine Finalize(GC, IMPORT, EXPORT, CLOCK, RC)
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION:
!    Finalize merely destroys the FVadv object that was created in Initialize
!    and releases the space for the persistent data .
!
!EOP
!=============================================================================
!BOC
! !LOCAL VARIABLES:

      character(len=ESMF_MAXSTR)    :: IAm
      integer                       :: STATUS
      character(len=ESMF_MAXSTR)    :: COMP_NAME

! Get my name and set-up traceback handle
! ---------------------------------------

      Iam = 'Finalize'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
      _VERIFY(STATUS)
      Iam = trim(COMP_NAME) // TRIM(Iam)

      ! Clean up FV if AdvCore is running without FV3_DynCoreIsRunning
      !--------------------------------------------------
      if (.NOT. FV3_DynCoreIsRunning) then
         call fv_end(FV_Atm, grids_on_my_pe, .false.)
      endif

      call MAPL_GenericFinalize(GC, IMPORT, EXPORT, CLOCK, RC)
      _VERIFY(STATUS)

      _RETURN(ESMF_SUCCESS)
      end subroutine Finalize


subroutine global_integral (QG,Q,PLE,IM,JM,KM,NQ)

      real(FVPRC), intent(OUT)   :: QG(NQ)
      real(FVPRC), intent(IN)    :: Q(IM,JM,KM,NQ)
      real(FVPRC), intent(IN)    :: PLE(IM,JM,KM+1)
      integer,     intent(IN)    :: IM,JM,KM,NQ
! Locals
      integer   :: k,n
      real(REAL8), allocatable ::    dp(:,:,:)
      real(FVPRC), allocatable :: qsum1(:,:)

      allocate(    dp(im,jm,km) )
      allocate( qsum1(im,jm)    )

! Compute Pressure Thickness
! --------------------------
      do k=1,KM
         dp(:,:,k) = PLE(:,:,k+1)-PLE(:,:,k)
      enddo

! Loop over Tracers
! -----------------
     do n=1,NQ
        qsum1(:,:) = 0.d0
        do k=1,KM
           qsum1(:,:) = qsum1(:,:) + Q(:,:,k,n)*dp(:,:,k)
        enddo
        qg(n) = g_sum(FV_Atm(1)%domain, qsum1, is,ie, js,je, FV_Atm(1)%ng, FV_Atm(1)%gridstruct%area_64, 1, .true.)
     enddo

     deallocate( dp )
     deallocate( qsum1 )

end subroutine global_integral

!EOC
!------------------------------------------------------------------------------

end module AdvCore_GridCompMod
