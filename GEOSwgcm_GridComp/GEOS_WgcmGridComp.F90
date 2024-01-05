
#include "MAPL_Generic.h"

!=============================================================================
!BOP

!  !MODULE: GEOS_WgcmGridCompMod -- A Module to compute wave properties

!  !INTERFACE:

module GEOS_WgcmGridCompMod

!  !USES:

    use ESMF
    use MAPL_Mod

    use GEOS_UMWMGridCompMod,      only : UMWM_SetServices => SetServices 
    use GEOS_WaveWatchGridCompMod, only :  WW3_SetServices => SetServices

    use bl_seaspray_mod, only : mabl_sea_spray => online_spray

    use, intrinsic :: ISO_FORTRAN_ENV

    implicit none
    private

    character(len=*), parameter :: WGCM_CONFIG_FILE = 'WGCM.rc'

    character(len=*), parameter :: wave_model_ww3   = 'WW3'
    character(len=*), parameter :: wave_model_umwm  = 'UMWM'
    character(len=*), parameter :: wave_model_data  = 'WM.data'
    character(len=*), parameter :: wave_model_idealized = 'JONSWAP'

    integer :: WM

!   Private state
!   -------------

    type WaveModel_State
        private

        type(ESMF_Config) :: CF                     ! Private Config
 
        logical :: verbose = .false.                ! verbose messages

        real    :: dt = 0.0                         ! time step, s 

        character(len=ESMF_MAXSTR) :: wave_model    ! name of the wave model 
    end type WaveModel_State


!   Hook for the ESMF
!   -----------------
    type WaveModel_Wrap
        type (WaveModel_State), pointer :: ptr => null()
    end type WaveModel_Wrap


!  !PUBLIC MEMBER FUNCTIONS:

    public SetServices

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

        call ESMF_ConfigLoadFile(self%CF, WGCM_CONFIG_FILE, __RC__)

        call ESMF_ConfigGetAttribute(self%CF, self%verbose, label='verbose:', default=.false., __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%wave_model, label='wave_model:', default='__unknown__', __RC__)

        self%wave_model = ESMF_UtilStringUpperCase(self%wave_model, __RC__)

        
        

! Add a child component (the NUOPC wrapped ww3)
! ---------------------------------------------
        WM = -1
        select case (self%wave_model)
            case (wave_model_umwm)
                WM = MAPL_AddChild(GC, NAME='UMWM', SS=UMWM_SetServices, __RC__)
            case (wave_model_ww3)
                WM = MAPL_AddChild(GC, NAME='WW3plug', SS=WW3_SetServices,  __RC__)
            case default
                __raise__(MAPL_RC_ERROR, 'Unrecognized wave model name in ' // trim(WGCM_CONFIG_FILE))
        end select


! Set the Initialize, Run entry point
! -----------------------------------

        call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_INITIALIZE, Initialize, __RC__)
        call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN,        Run,        __RC__)
        call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_FINALIZE,   Finalize,   __RC__)


! Store private internal state in GC
! ----------------------------------
        call ESMF_UserCompSetInternalState(GC, 'WaveModel_State', wrap, STATUS)
        VERIFY_(STATUS)


! Set the state variable specs
! -----------------------------

! !INTERNAL STATE:
       



! !IMPORT STATE:

!  AGCM -> WGCM

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
            SHORT_NAME     = 'U10N',                                 &
            LONG_NAME      = 'equivalent_neutral_10-meter_eastward_wind', &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            RESTART        = MAPL_RestartSkip,        __RC__) 
  
        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'V10N',                                 &
            LONG_NAME      = 'equivalent_neutral_10-meter_northward_wind', &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            RESTART        = MAPL_RestartSkip,        __RC__)


        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'RHOS',                                 &
            LONG_NAME      = 'air_density_at_surface',               &
            UNITS          = 'kg m-3',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 1.28,                                   &
            RESTART        = MAPL_RestartOptional,    __RC__)

        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'TSKINW',                               &
            LONG_NAME      = 'open_water_skin_temperature',          &
            UNITS          = 'K',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 280.0,                                  &
            RESTART        = MAPL_RestartOptional,    __RC__)

        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'TS',                                   &
            LONG_NAME      = 'skin_temperature',                     &
            UNITS          = 'K',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 280.0,                                  &
            RESTART        = MAPL_RestartOptional,    __RC__)

        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'FRLAND',                               &
            LONG_NAME      = 'fraction_of_land',                     &
            UNITS          = '1',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            RESTART        = MAPL_RestartSkip,        __RC__)

        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'FROCEAN',                              &
            LONG_NAME      = 'fraction_of_ocean',                    &
            UNITS          = '1',                                    &
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

        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'PS',                                   &
            LONG_NAME      = 'surface_pressure',                     &
            UNITS          = 'Pa',                                   &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 101325.0,                               &
            RESTART        = MAPL_RestartSkip,       __RC__)

       call MAPL_AddImportSpec(GC,                                   &
            SHORT_NAME     = 'Q10M',                                 &
            LONG_NAME      = '10-meter_specific_humidity',           &
            UNITS          = 'kg kg',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 0.8,                                    &
            RESTART        = MAPL_RestartSkip,       __RC__)

       call MAPL_AddImportSpec(GC,                                   &
            SHORT_NAME     = 'RH2M',                                 &
            LONG_NAME      = 'near-surface_relative_humidity',       &
            UNITS          = '%',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 0.8,                                    &
            RESTART        = MAPL_RestartSkip,       __RC__)

       call MAPL_AddImportSpec(GC,                                   &
            SHORT_NAME     = 'T10M',                                 &
            LONG_NAME      = '10-meter_air_temperature',             &
            UNITS          = 'K',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 288.0,                                  &
            RESTART        = MAPL_RestartSkip,       __RC__)

        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'LHFX',                                 &
            LONG_NAME      = 'total_latent_energy_flux',             &
            UNITS          = 'W m-2',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 0.0,                                    &
            RESTART        = MAPL_RestartSkip,       __RC__)

        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'SH',                                   &
            LONG_NAME      = 'sensible_heat_flux_from_turbulence',   &
            UNITS          = 'W m-2',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 0.0,                                    &
            RESTART        = MAPL_RestartSkip,       __RC__)


!  OGCM -> WGCM
  
        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'TW',                                   &
            LONG_NAME      = 'temperature',                          &
            UNITS          = 'K',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            RESTART        = MAPL_RestartSkip,                       &
            VLOCATION      = MAPL_VLocationNone,      __RC__)

        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'UW',                                   &
            LONG_NAME      = 'zonal_velocity_of_surface_water',      &
            UNITS          = 'm s-1 ',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 0.0,                                    &
            RESTART        = MAPL_RestartOptional,    __RC__)

        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'VW',                                   &
            LONG_NAME      = 'meridional_velocity_of_surface_water', &
            UNITS          = 'm s-1 ',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 0.0,                                    &
            RESTART        = MAPL_RestartOptional,    __RC__)

        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'DW_WGCM',                              &
            LONG_NAME      = 'sea_floor_depth',                      &
            UNITS          = 'm',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            RESTART        = MAPL_RestartSkip,        __RC__)


! !EXPORT STATE:

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'CHARNOCK',                             &
            CHILD_ID       = WM,                     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'USTAR',                                &
            CHILD_ID       = WM,                     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'SWH',                                  &
            CHILD_ID       = WM,                     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'DCP',                                  &
            CHILD_ID       = WM,                     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
             SHORT_NAME    = 'EDF',                                  &
             CHILD_ID      = WM,                     __RC__)



        !
        ! Sea spray diagnostics
        !
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'EDFP',                                 &
            LONG_NAME      = 'wave_energy_dissipation_flux_parameterized', &
            UNITS          = 'kg s-3',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            LONG_NAME      = 'sensible_heat_flux_from_turbulence',   &
            UNITS          = 'W m-2',                                &
            SHORT_NAME     = 'SHFX',                                 &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__) 

        call MAPL_AddExportSpec(GC,                                  &
            LONG_NAME      = 'total_latent_energy_flux',             &
            UNITS          = 'W m-2',                                &
            SHORT_NAME     = 'LHFX',                                 &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME      = 'SHFX_TURB',                            &
           LONG_NAME       = 'sensible_heat_carried_by_turbulence',  &
           UNITS           = 'W m-2',                                &
           DIMS            = MAPL_DimsHorzOnly,                      &
           VLOCATION       = MAPL_VLocationNone,     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME      = 'LHFX_TURB',                            &
           LONG_NAME       = 'latent_heat_carried_by_turbulence',    &
           UNITS           = 'W m-2',                                &
           DIMS            = MAPL_DimsHorzOnly,                      &
           VLOCATION       = MAPL_VLocationNone,     __RC__)
 
        call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME      = 'SHFX_TOT',                             &
           LONG_NAME       = 'sensible_heat_medited_by_seaspray',    &
           UNITS           = 'W m-2',                                &
           DIMS            = MAPL_DimsHorzOnly,                      &
           VLOCATION       = MAPL_VLocationNone,     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME      = 'LHFX_TOT',                             &
           LONG_NAME       = 'latent_heat_mediated_by_seaspray',     &
           UNITS           = 'W m-2',                                &
           DIMS            = MAPL_DimsHorzOnly,                      &
           VLOCATION       = MAPL_VLocationNone,     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME      = 'SHFX_SPRAY',                           &
           LONG_NAME       = 'sensible_heat_contribution_from_sea_spray', &
           UNITS           = 'W m-2',                                &
           DIMS            = MAPL_DimsHorzOnly,                      &
           VLOCATION       = MAPL_VLocationNone,     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME      = 'LHFX_SPRAY',                           &
           LONG_NAME       = 'latent_heat_contribution_from_sea_spray',   &
           UNITS           = 'W m-2',                                &
           DIMS            = MAPL_DimsHorzOnly,                      &
           VLOCATION       = MAPL_VLocationNone,     __RC__)


! Set the Profiling timers
! ------------------------

        call MAPL_TimerAdd(GC, name='TOTAL'       , __RC__)
        call MAPL_TimerAdd(GC, name='INITIALIZE'  , __RC__)
        call MAPL_TimerAdd(GC, name='RUN'         , __RC__)
        call MAPL_TimerAdd(GC, name='-SEA_SPRAY'  , __RC__)
        call MAPL_TimerAdd(GC, name='FINALIZE'    , __RC__)


!ALT: we need to terminate child's import so they do not "bubble up". We will fill them explicitly

! this should be irrelevant here because the children are not MAPL components
!        call MAPL_TerminateImport()

! Set generic init and final methods
! ----------------------------------

        call MAPL_GenericSetServices(GC, __RC__)

        RETURN_(ESMF_SUCCESS)

    end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! ! IROUTINE: INITIALIZE -- Initialize method for the UMWM component

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

      character(len=ESMF_MAXSTR) :: Iam
      integer                    :: STATUS
      character(len=ESMF_MAXSTR) :: COMP_NAME

! Local derived type aliases

      type(MAPL_MetaComp), pointer    :: MAPL
      type(ESMF_Grid)                 :: GRID

      type(ESMF_Alarm)                :: run_alarm
      type(ESMF_TimeInterval)         :: ring_interval
      real(ESMF_KIND_R8)              :: time_step

      type (WaveModel_State), pointer :: self   ! private internal state
      type (WaveModel_Wrap)           :: wrap

! Local Variables

      type (ESMF_GridComp),      pointer  :: GCS(:) => null()
      type (ESMF_State),         pointer  :: GIM(:) => null()
      type (ESMF_State),         pointer  :: GEX(:) => null()
 

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


! Start the timers
! ----------------

      call MAPL_TimerOn(MAPL, 'TOTAL',        __RC__)
      call MAPL_TimerOn(MAPL, 'INITIALIZE',  __RC__)


! Set the grid explicitly if the WM instance is WW3
! -------------------------------------------------

      ! this section needs to be executed only for WW3:
      ! propagate the WW3 grid up to the WGCM.
#if 0
      if (self%wave_model == wave_model_ww3) then 
          call MAPL_Get(MAPL, GCS=GCS, GIM=GIM, GEX=GEX, __RC__)
          call MAPL_Set(MAPL, ChildInit=.false., __RC__)

          call ESMF_GridCompInitialize(GCS(WM), importState=GIM(WM), &
                   exportState=GEX(WM), clock=CLOCK, userRC=STATUS)
          VERIFY_(STATUS)

          call ESMF_GridCompGet(GCS(WM), grid=GRID, __RC__)
          call ESMF_GridCompSet(GC, grid=GRID, __RC__)
      end if
#endif

! Get the grid
! ------------

      call ESMF_GridCompGet(GC, grid=GRID, __RC__)


! Generic initialize
! ------------------
      call MAPL_GenericInitialize(GC, IMPORT, EXPORT, CLOCK, __RC__)


! Get parameters from generic state
! ---------------------------------

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

! ! IROUTINE: RUN -- Run method for the UMWM component

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

      type(WaveModel_State), pointer :: self => null()
      type(WaveModel_Wrap)           :: wrap

      type(MAPL_MetaComp), pointer   :: MAPL
      type(ESMF_State)               :: INTERNAL
      type(ESMF_Grid)                :: GRID
      type(ESMF_VM)                  :: VM

      type (ESMF_GridComp), pointer  :: GCS(:)
      type (ESMF_State),    pointer  :: GIM(:)
      type (ESMF_State),    pointer  :: GEX(:)


! Pointers from Import state       
      real, pointer, dimension(:,:) :: U10M => null()
      real, pointer, dimension(:,:) :: V10M => null()
      real, pointer, dimension(:,:) :: RHOS => null()
      real, pointer, dimension(:,:) :: PS   => null()
      real, pointer, dimension(:,:) :: TS   => null()
      real, pointer, dimension(:,:) :: T10M => null()
      real, pointer, dimension(:,:) :: RH2M => null()
      real, pointer, dimension(:,:) :: FRACI=> null()

      real, pointer, dimension(:,:) :: LHFX => null()
      real, pointer, dimension(:,:) :: SHFX => null()

! Pointers to my Export state

      real, pointer, dimension(:,:) :: WM_USTAR
      real, pointer, dimension(:,:) :: WM_SWH
      real, pointer, dimension(:,:) :: WM_DCP
      real, pointer, dimension(:,:) :: WM_EDF

      !
      ! Sea spray diagnostics
      !
      real, pointer, dimension(:,:) :: WM_EDFP      => null()

      real, pointer, dimension(:,:) :: WM_LHFX      => null()
      real, pointer, dimension(:,:) :: WM_SHFX      => null()

      real, pointer, dimension(:,:) :: SHFX_TURB    => null()
      real, pointer, dimension(:,:) :: LHFX_TURB    => null()
      real, pointer, dimension(:,:) :: SHFX_TOT     => null()
      real, pointer, dimension(:,:) :: LHFX_TOT     => null()
      real, pointer, dimension(:,:) :: SHFX_SPRAY   => null()
      real, pointer, dimension(:,:) :: LHFX_SPRAY   => null()

! Pointers to child's Export state
!     N/A

! Local variables
      integer :: IM, JM
      integer :: i, j

! MABL sea spray
      real :: spray_edf_factor
      real :: spray_hss
      real :: spray_hll
      real :: spray_hwave
      real :: spray_cwave
      real :: spray_p
      real :: spray_usr
      real :: spray_w10m
      real :: spray_massf
      real :: spray_hs_tot
      real :: spray_hl_tot
      real :: spray_usr_new
      real :: spray_S_bar1
      real :: spray_z_r
      real :: spray_omega
      real :: spray_alpha
      real :: spray_vfm

      integer :: DO_SEA_SPRAY

      real, parameter :: SPRAY_SOURCE_STRENGTH = 0.4
      real, parameter :: SPRAY_FEEDBACK        = 0.2
 
! TODO: need to be consistent across the wave models; move to a config file
      real, parameter :: FRACTION_ICE_SUPPRESS_WAVES = 0.8



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
      call MAPL_TimerOn(MAPL, 'RUN',  __RC__)


! Get parameters from generic state
! ---------------------------------
      call MAPL_Get (MAPL, IM=IM, JM=JM, GCS=GCS, GIM=GIM, GEX=GEX, __RC__)


! Get my internal private state
! -----------------------------
      call ESMF_UserCompGetInternalState(GC, 'WaveModel_State', wrap, STATUS)
      VERIFY_(STATUS)

      self => wrap%ptr


! Run children (specific WM)
! --------------------------
      call MAPL_GenericRunChildren (GC, IMPORT, EXPORT, CLOCK, RC=STATUS)
      VERIFY_(STATUS)

! MABL sea spray parameterization, Bao et al, 2011
! ------------------------------------------------
    call MAPL_TimerOn(MAPL, '-SEA_SPRAY')

    call MAPL_GetResource( MAPL, DO_SEA_SPRAY, Label="USE_SEA_SPRAY:", DEFAULT=1, __RC__)

! Get pointers to inputs
! ----------------------
      call MAPL_GetPointer(IMPORT, U10M,   'U10M',    __RC__)
      call MAPL_GetPointer(IMPORT, V10M,   'V10M',    __RC__)
      call MAPL_GetPointer(IMPORT, T10M,   'T10M',    __RC__)
      call MAPL_GetPointer(IMPORT, RH2M,   'RH2M',    __RC__)
      call MAPL_GetPointer(IMPORT, RHOS,   'RHOS',    __RC__)
      call MAPL_GetPointer(IMPORT, TS,     'TS',      __RC__)
      call MAPL_GetPointer(IMPORT, PS,     'PS',      __RC__)
      call MAPL_GetPointer(IMPORT, FRACI,  'FRACI',   __RC__)
      call MAPL_GetPointer(IMPORT, LHFX,   'LHFX',    __RC__)
      call MAPL_GetPointer(IMPORT, SHFX,   'SH',      __RC__)


      call MAPL_GetPointer(EXPORT, WM_USTAR,   'USTAR',     __RC__)
      call MAPL_GetPointer(EXPORT, WM_SWH,     'SWH',       __RC__)
      call MAPL_GetPointer(EXPORT, WM_DCP,     'DCP',       __RC__)
      call MAPL_GetPointer(EXPORT, WM_EDF,     'EDF',       __RC__)
      
      call MAPL_GetPointer(EXPORT, WM_LHFX,    'LHFX',      __RC__)
      call MAPL_GetPointer(EXPORT, WM_SHFX,    'SHFX',      __RC__)

      ! sea spray diagnostics
      call MAPL_GetPointer(EXPORT, WM_EDFP,    'EDFP',      __RC__)

      call MAPL_GetPointer(EXPORT, LHFX_TURB,  'LHFX_TURB' , alloc=(DO_SEA_SPRAY/=0), __RC__)
      call MAPL_GetPointer(EXPORT, SHFX_TURB,  'SHFX_TURB' , alloc=(DO_SEA_SPRAY/=0), __RC__)
      call MAPL_GetPointer(EXPORT, LHFX_TOT,   'LHFX_TOT'  , alloc=(DO_SEA_SPRAY/=0), __RC__)
      call MAPL_GetPointer(EXPORT, SHFX_TOT,   'SHFX_TOT'  , alloc=(DO_SEA_SPRAY/=0), __RC__)
      call MAPL_GetPointer(EXPORT, LHFX_SPRAY, 'LHFX_SPRAY', alloc=(DO_SEA_SPRAY/=0), __RC__)
      call MAPL_GetPointer(EXPORT, SHFX_SPRAY, 'SHFX_SPRAY', alloc=(DO_SEA_SPRAY/=0), __RC__)

! Sanity diagnostics
! ------------------

      ! heat fluxes
      if (associated(WM_LHFX))    WM_LHFX = LHFX
      if (associated(WM_SHFX))    WM_SHFX = SHFX

    
    PARAMETERIZED_ENERGY_DISSIPATION_FLUX: if (associated(WM_EDFP)) then
         
         WM_EDFP = 1.0*(-0.4+0.25*sqrt(U10M**2 + V10M**2))  ! cEwave
         where (WM_USTAR /= MAPL_UNDEF)
             WM_EDFP = (RHOS/MAPL_RHOWTR) * WM_EDFP*WM_USTAR**2
         elsewhere
             WM_EDFP = 0.0
         end where

         where (WM_EDFP < 0.0 .and. WM_USTAR /= MAPL_UNDEF)
             WM_EDFP = (RHOS/MAPL_RHOWTR) * 3.5 * WM_USTAR**3.5            
         end where
         
         WM_EDFP = WM_EDFP  * MAPL_RHOWTR  ! convert units from 'm3 s-3' to match the units of EDF 'kg s-3'

    end if PARAMETERIZED_ENERGY_DISSIPATION_FLUX



    DIAGNOSTICS_SPRAY_FLUXES: if ( associated(SHFX_SPRAY) .or. &
                                   associated(LHFX_SPRAY) ) then
      
      ASSERT_(associated(SHFX_SPRAY))
      ASSERT_(associated(LHFX_SPRAY))

      ASSERT_(associated(SHFX_TURB))
      ASSERT_(associated(LHFX_TURB))

      ASSERT_(associated(SHFX_TOT))
      ASSERT_(associated(LHFX_TOT))

      SHFX_TOT   = SHFX
      LHFX_TOT   = LHFX
      SHFX_TURB  = SHFX
      LHFX_TURB  = LHFX
      SHFX_SPRAY = 0.0
      LHFX_SPRAY = 0.0

      if (self%wave_model == wave_model_umwm) then
          ! scale up UMWM:EDF to agree with EDFP and WW3
          spray_edf_factor = 3.0
      else
          spray_edf_factor = 1.0
      end if 

      do j = 1, JM
          do i = 1, IM
          
          spray_w10m = sqrt(U10M(i,j)**2 + V10M(i,j)**2)

          if ( (spray_w10m > 7.0) .and. &
               (WM_USTAR(i,j) /= MAPL_UNDEF) .and. &
               (WM_DCP(i,j) /= MAPL_UNDEF) .and. &
               (WM_SWH(i,j) > 0.5) .and. &
               (FRACI(i,j) < FRACTION_ICE_SUPPRESS_WAVES) ) then

              spray_hss   = SHFX(i,j)     ! SHWTR
              spray_hll   = LHFX(i,j)     ! HLATWTR
              spray_hwave = WM_SWH(i,j)
              spray_cwave = WM_DCP(i,j)
#if(1)
              spray_p = spray_edf_factor * WM_EDF(i,j)  / MAPL_RHOWTR   ! the factor 3 is to agree with EDFP

              ! protect against negative energy dissipation flux 
              if (spray_p < 0) then
                  spray_p = RHOS(i,j) * 3.5 * WM_USTAR(i,j)**3.5
              end if

              if (spray_p < tiny(spray_p)) then
                  cycle
              end if
#else
              spray_p     = WM_EDFP(i,j) / MAPL_RHOWTR
#endif
              spray_usr   = WM_USTAR(i,j)
              
              
              call mabl_sea_spray(SPRAY_SOURCE_STRENGTH, &
                                  SPRAY_FEEDBACK,        &
                                  spray_w10m,            &
                                  10.0,                  &
                                  TS(i,j)   - 273.15,    &  ! T:=(TSKINW!=MAPL_UNDEF)?TSKINW:TS
                                  T10M(i,j) - 273.15,    &
                                  min(RH2M(i,j) * 0.01, 0.99),   &  ! s = ...?
                                  PS(i,j)   * 0.01,      &
                                  spray_hss,             &
                                  spray_hll,             &
                                  spray_hwave,           &
                                  spray_cwave,           &
                                  spray_p,               &
                                  spray_usr,             &
                                  spray_massf,           &
                                  spray_hs_tot,          &
                                  spray_hl_tot,          &
                                  spray_usr_new,         &
                                  spray_S_bar1,          &
                                  spray_z_r,             &
                                  spray_omega,           &
                                  spray_alpha,           &
                                  spray_vfm)

              if (abs(spray_hss - SHFX(i,j)) > tiny(spray_hss) .or. &
                  abs(spray_hll - LHFX(i,j)) > tiny(spray_hll)) then
                  print *, 'DEBUG::WAVES_PHYS  ***Heat Fluxes do not match after SPRAY()'
              end if    
   
              SHFX_TURB(i,j)  = spray_hss
              LHFX_TURB(i,j)  = spray_hll

              SHFX_SPRAY(i,j) = (spray_hs_tot - spray_hss) * (1 - FRACI(i,j))
              LHFX_SPRAY(i,j) = (spray_hl_tot - spray_hll) * (1 - FRACI(i,j))

              SHFX_TOT(i,j)   = spray_hss + (spray_hs_tot - spray_hss) * (1 - FRACI(i,j))
              LHFX_TOT(i,j)   = spray_hll + (spray_hl_tot - spray_hll) * (1 - FRACI(i,j))
          end if


          if (abs(SHFX_SPRAY(i,j)) > 2e2 .or. abs(LHFX_SPRAY(i,j)) > 2e2) then
              print *, SHFX_SPRAY(i,j),   &
                       LHFX_SPRAY(i,j),   &
                       WM_USTAR(i,j),     &
                       sqrt(U10M(i,j)**2 + V10M(i,j)**2), &
                       TS(i,j)  - 273.15, & 
                       T10M(i,j)- 273.15, &
                       RH2M(i,j) * 0.01,  &  ! s = ...?
                       PS(i,j)   * 0.01,  &
                       WM_SWH(i,j),       &
                       WM_DCP(i,j),       &
#if(1)
                       WM_EDF(i,j)
#else
                       WM_EDFP(i,j)
#endif
          end if

          end do
      end do
       
#ifdef DEBUG
      print *, ' *** DEBUG   WM:SH_SPRAY = ', minval(SHFX_SPRAY), maxval(SHFX_SPRAY)
      print *, ' *** DEBUG   WM:LH_SPRAY = ', minval(LHFX_SPRAY), maxval(LHFX_SPRAY)
#endif


   end if DIAGNOSTICS_SPRAY_FLUXES

   call MAPL_TimerOff(MAPL, '-SEA_SPRAY')



! Stop the timers
! ---------------

      call MAPL_TimerOff(MAPL, 'RUN',  __RC__)
      call MAPL_TimerOff(MAPL, 'TOTAL', __RC__)

! All Done
! --------

      RETURN_(ESMF_SUCCESS)

   end subroutine Run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! ! IROUTINE: FINALIZE -- Finalize method for the UMWM component

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

      call MAPL_TimerOn(MAPL, 'TOTAL',     __RC__)
      call MAPL_TimerOn(MAPL, 'FINALIZE', __RC__)

! Get parameters from generic state.
! ----------------------------------

      call MAPL_Get(MAPL, IM=IM, JM=JM, LM=LM, __RC__)

      call MAPL_GridGet(GRID, globalCellCountPerDim=COUNTS, __RC__)

      IM_world = COUNTS(1)
      JM_world = COUNTS(2)

! Get layout from the grid
! ------------------------

      call ESMF_VMGetCurrent(VM, __RC__)

      call ESMF_VMGet(VM, mpiCommunicator=COMM, localPet=myPE, petCount=nPEs, __RC__)

! Get parameters
!---------------






! Stop the timers
! ---------------

      call MAPL_TimerOff(MAPL, 'FINALIZE', __RC__)
      call MAPL_TimerOff(MAPL, 'TOTAL',     __RC__)

! Call GenericFinalize
! ----------------------
      call MAPL_MemUtilsWrite(VM, 'GEOSUMWM_GridComp:BeforeGenericFinalize', __RC__)

      call MAPL_GenericFinalize(GC, IMPORT, EXPORT, CLOCK, __RC__)
      VERIFY_(STATUS)

      call MAPL_MemUtilsWrite(VM, 'GEOSUMWM_GridComp:AfterGenericFinalize', __RC__)

! All Done
! --------

      RETURN_(ESMF_SUCCESS)

   end subroutine Finalize





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_WgcmGridCompMod

