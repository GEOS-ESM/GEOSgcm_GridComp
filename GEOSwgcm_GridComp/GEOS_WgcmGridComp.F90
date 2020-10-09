
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
            SHORT_NAME     = 'DW',                                   &
            LONG_NAME      = 'bathymetry',                           &
            UNITS          = 'm',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 0.0,                                    &
            RESTART        = MAPL_RestartOptional,    __RC__)

!       call MAPL_AddImportSpec(GC,                                  &
!           SHORT_NAME     = 'SSKINW',                               &
!           LONG_NAME      = 'water_skin_salinity',                  &
!           UNITS          = 'psu',                                  &
!           DIMS           = MAPL_DimsHorzOnly,                      &
!           VLOCATION      = MAPL_VLocationNone,                     &
!           DEFAULT        = 30.0,                                   &
!           RESTART        = MAPL_RestartSkip,       __RC__)
!
!       call MAPL_AddImportSpec(GC,                                  &
!           SHORT_NAME     = 'HLATN',                                &
!           LONG_NAME      = 'total_latent_energy_flux',             &
!           UNITS          = 'W m-2',                                &
!           DIMS           = MAPL_DimsHorzOnly,                      &
!           VLOCATION      = MAPL_VLocationNone,                     &
!           DEFAULT        = 0.0,                                    &
!           RESTART        = MAPL_RestartSkip,       __RC__)
!
!
!       call MAPL_AddImportSpec(GC,                                  &
!           SHORT_NAME     = 'SHOUT',                                &
!           LONG_NAME      = 'upward_sensible_heat_flux',            &
!           UNITS          = 'W m-2',                                &
!           DIMS           = MAPL_DimsHorzOnly,                      &
!           VLOCATION      = MAPL_VLocationNone,                     &
!           DEFAULT        = 0.0,                                    &
!           RESTART        = MAPL_RestartSkip,       __RC__)
!


! !EXPORT STATE:
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'U10M',                                 &
            LONG_NAME      = '10-meter_eastward_wind',               &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,      __RC__) 
  
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'V10M',                                 &
            LONG_NAME      = '10-meter_northward_wind',              &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,      __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'U10N',                                 &
            LONG_NAME      = 'equivalent_neutral_10-meter_eastward_wind', &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,      __RC__) 
  
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'V10N',                                 &
            LONG_NAME      = 'equivalent_neutral_10-meter_northward_wind', &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,      __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'W10N',                                 &
            LONG_NAME      = 'equivalent_neutral_10-meter_wind',     &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,      __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'W10M',                                 &
            LONG_NAME      = '10-meter_wind',                        &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,      __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'RHOS',                                 &
            LONG_NAME      = 'air_density_at_surface',               &
            UNITS          = 'kg m-3',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,      __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'UW',                                   &
            LONG_NAME      = 'zonal_velocity_of_surface_water',      &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,      __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'VW',                                   &
            LONG_NAME      = 'meridional_velocity_of_surface_water', &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,      __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'DW',                                   &
            LONG_NAME      = 'bathymetry',                           &
            UNITS          = 'm',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,      __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'NUW',                                  &
            LONG_NAME      = 'sea_water_kinematic_viscosity',        &
            UNITS          = 'm2 s-1',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'FRACI',                                &
            LONG_NAME      = 'ice_covered_fraction_of_tile',         &
            UNITS          = '1',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,      __RC__)


        !
        ! WM diagnostics
        !

#if (1)
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'CHARNOCK',                             &
            CHILD_ID       = WM,                     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'SHFX_SPRAY',                             &
            CHILD_ID       = WM,                     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'LHFX_SPRAY',                             &
            CHILD_ID       = WM,                     __RC__)
#else
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'SWH',                                  &
            LONG_NAME      = 'sea_surface_wave_significant_height',  &
            UNITS          = 'm',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'SWHW',                                 &
            LONG_NAME      = 'sea_surface_wind_wave_significant_height', &
            UNITS          = 'm',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'SWHS',                                 &
            LONG_NAME      = 'sea_surface_swell_significant_height', &
            UNITS          = 'm',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'MWP',                                  &
            LONG_NAME      = 'mean_wave_period',                     &
            UNITS          = 's',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'MWD',                                  &
            LONG_NAME      = 'mean_wave_direction',                  &
            UNITS          = 'rad',                                  &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'MSS',                                  &
            LONG_NAME      = 'mean_squared_slope',                   &
            UNITS          = 'rad',                                  &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'MWL',                                  &
            LONG_NAME      = 'mean_wavelength',                      &
            UNITS          = 'm',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'DWD',                                  &
            LONG_NAME      = 'dominant_wave_direction',              &
            UNITS          = 'rad',                                  &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'DWL',                                  &
            LONG_NAME      = 'dominant_wave_length',                 &
            UNITS          = 'm',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'DWP',                                  &
            LONG_NAME      = 'dominant_wave_period',                 &
            UNITS          = 's',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'DCP0',                                 &
            LONG_NAME      = 'dominant_phase_speed_intrinsic',       &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'DCG0',                                 &
            LONG_NAME      = 'dominant_group_speed_intrinsic',       &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'DCP',                                  &
            LONG_NAME      = 'dominant_phase_speed',                 &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'DCG',                                  &
            LONG_NAME      = 'dominant_group_speed',                 &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'EDFX',                                 &
            LONG_NAME      = 'wave_energy_dissipation_flux_x_component', &
            UNITS          = 'kg s-3',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'EDFY',                                 &
            LONG_NAME      = 'wave_energy_dissipation_flux_y_component', &
            UNITS          = 'kg s-3',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'EDF',                                  &
            LONG_NAME      = 'wave_energy_dissipation_flux',         &
            UNITS          = 'kg s-3',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'EGFX',                                 &
            LONG_NAME      = 'wave_energy_growth_flux_x_component',  &
            UNITS          = 'kg s-3',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'EGFY',                                 &
            LONG_NAME      = 'wave_energy_growth_flux_y_component',  &
            UNITS          = 'kg s-3',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'EGF',                                  &
            LONG_NAME      = 'wave_energy_growth_flux',              &
            UNITS          = 'kg s-3',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'Z0',                                   &
            LONG_NAME      = 'surface_roughness',                    &
            UNITS          = 'm',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__) 

        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'CD',                                   &
            LONG_NAME      = 'drag_coefficient_of_air',              &
            UNITS          = '1',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__)
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'CHARNOCK',                             &
            LONG_NAME      = 'wave_model_charnock_coefficient',      &
            UNITS          = '1',                                    &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__) 
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'TAU',                                  &
            LONG_NAME      = 'total drag',                           &
            UNITS          = 'N m-2',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__) 
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'TAUX',                                 &
            LONG_NAME      = 'total drag, x-component',              &
            UNITS          = 'N m-2',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__) 
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'TAUY',                                 &
            LONG_NAME      = 'total drag, y-component',              &
            UNITS          = 'N m-2',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__) 
        
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'TAU_FORM',                             &
            LONG_NAME      = 'form drag',                            &
            UNITS          = 'N m-2',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__) 
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'TAUX_FORM',                            &
            LONG_NAME      = 'form drag, x-component',               &
            UNITS          = 'N m-2',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__) 
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'TAUY_FORM',                            &
            LONG_NAME      = 'form drag, y-component',               &
            UNITS          = 'N m-2',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__) 
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'TAU_SKIN',                             &
            LONG_NAME      = 'skin drag',                            &
            UNITS          = 'N m-2',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__) 
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'TAUX_SKIN',                            &
            LONG_NAME      = 'skin drag, x-component',               &
            UNITS          = 'N m-2',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__) 
    
        call MAPL_AddExportSpec(GC,                                  &
            SHORT_NAME     = 'TAUY_SKIN',                            &
            LONG_NAME      = 'skin drag, y-component',               &
            UNITS          = 'N m-2',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,     __RC__) 


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
#endif

! Set the Profiling timers
! ------------------------

        call MAPL_TimerAdd(GC, name='TOTAL'        , __RC__)
        call MAPL_TimerAdd(GC, name='INITIALIZE'  , __RC__)
        call MAPL_TimerAdd(GC, name='RUN'         , __RC__)
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

! this section needs to be executed only for WW3
      call MAPL_Get(MAPL, GCS=GCS, GIM=GIM, GEX=GEX, __RC__)
      call MAPL_Set(MAPL, ChildInit=.false., __RC__)

      call ESMF_GridCompInitialize(GCS(WM), importState=GIM(WM), &
           exportState=GEX(WM), clock=CLOCK, userRC=status )
      VERIFY_(STATUS)

      call ESMF_GridCompGet(GCS(WM), grid=grid, __RC__)
      call ESMF_GridCompSet(GC, grid=grid, __RC__)
! end of WW3 section

! Get the grid
! ------------

      call ESMF_GridCompGet( GC, grid=GRID, __RC__)


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

! Local Variables


! Pointers to my Export state

      real, pointer, dimension(:,:) :: z0 => null()
      real, pointer, dimension(:,:) :: charnock => null()

! Pointers to child's Export state

      real, pointer, dimension(:,:) :: z0rlen => null()
      real, pointer, dimension(:,:) :: charno => null()

      type (ESMF_GridComp),      pointer  :: GCS(:)
      type (ESMF_State),         pointer  :: GIM(:)
      type (ESMF_State),         pointer  :: GEX(:)
      

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

! Get parameters from generic state.
! ----------------------------------
      call MAPL_Get (MAPL, GCS=GCS, GIM=GIM, GEX=GEX, __RC__)


! Get my internal private state
! -----------------------------
      call ESMF_UserCompGetInternalState(GC, 'WaveModel_State', wrap, STATUS)
      VERIFY_(STATUS)

      self => wrap%ptr

! Get pointers to inputs
! ----------------------
#if (0)
      call MAPL_GetPointer(IMPORT, U10M,   'U10M',    __RC__)
      call MAPL_GetPointer(IMPORT, V10M,   'V10M',    __RC__)
      call MAPL_GetPointer(IMPORT, FRACI, 'FRACI',    __RC__)
#endif

      if (MAPL_AM_I_Root()) write (OUTPUT_UNIT,*) 'DEBUG::WGCM  Run...'

      call MAPL_GenericRunChildren (GC, IMPORT, EXPORT, CLOCK, RC=STATUS )

      if (MAPL_AM_I_Root()) write (OUTPUT_UNIT,*) 'DEBUG::WGCM  ...done.'

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

