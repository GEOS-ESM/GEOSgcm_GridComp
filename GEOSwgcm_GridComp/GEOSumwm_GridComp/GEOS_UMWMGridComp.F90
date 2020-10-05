#include "MAPL_Generic.h"

!=============================================================================
!BOP

!  !MODULE: GEOS_WgcmGridCompMod -- A Module to compute wave properties via the
!          UMWM wave model

!  !INTERFACE:

module GEOS_UMWMGridCompMod

!  !USES:

    use ESMF
    use MAPL_Mod

!   UMWM modules
    use UMWM_module,        only: umwm_version       => version

    use UMWM_module,        only: umwm_isGlobal      => isGlobal,      &
                                  umwm_restart       => restart,       &
                                  umwm_fillEstuaries => fillEstuaries, &
                                  umwm_fillLakes     => fillLakes

    use UMWM_module,       only:  umwm_mm            => mm,            &
                                  umwm_nm            => nm,            &
                                  umwm_lat           => lat,           &
                                  umwm_lon           => lon,           &
                                  umwm_d_2d          => d_2d,          &
                                  umwm_mask          => mask

    use UMWM_module,        only: umwm_om            => om,            &
                                  umwm_pm            => pm,            &
                                  umwm_fmin          => fmin,          &
                                  umwm_fmax          => fmax,          &
                                  umwm_fprog         => fprog

    use UMWM_module,        only: umwm_wspd          => wspd,          &
                                  umwm_wdir          => wdir,          &
                                  umwm_rhow0         => rhow0,         &
                                  umwm_rhow          => rhow,          &
                                  umwm_rhoa0         => rhoa0,         &
                                  umwm_rhoa          => rhoa,          &
                                  umwm_rhorat        => rhorat,        &
                                  umwm_nu_water_     => nu_water_,     &
                                  umwm_fice          => fice,          &
                                  umwm_uc            => uc,            &
                                  umwm_vc            => vc,            &
                                  umwm_d             => d

    use UMWM_module,        only: umwm_nu_air        => nu_air,        &
                                  umwm_nu_water      => nu_water,      & 
                                  umwm_g             => g,             &
                                  umwm_kappa         => kappa,         &
                                  umwm_gustiness     => gustiness,     &
                                  umwm_sfct          => sfct,          &   
                                  umwm_z             => z,             &  
                                  umwm_dmin          => dmin,          &
                                  umwm_explim        => explim,        &  
                                  umwm_sin_fac       => sin_fac,       &  
                                  umwm_sin_diss1     => sin_diss1,     &
                                  umwm_sin_diss2     => sin_diss2,     &
                                  umwm_sds_fac       => sds_fac,       &
                                  umwm_sds_power     => sds_power,     &
                                  umwm_mss_fac       => mss_fac,       &
                                  umwm_snl_fac       => snl_fac,       &
                                  umwm_sdt_fac       => sdt_fac,       &
                                  umwm_sbf_fac       => sbf_fac,       &
                                  umwm_sbp_fac       => sbp_fac,       &
                                  umwm_fice_lth      => fice_lth,      &
                                  umwm_fice_uth      => fice_uth

    use UMWM_module,        only: umwm_e             => e,             &
                                  umwm_ef            => ef,            &
                                  umwm_k             => k,             &
                                  umwm_cp0           => cp0,           &
                                  umwm_cg0           => cg0

    use UMWM_module,        only: umwm_im            => im,            &
                                  umwm_imm           => imm,           &
                                  umwm_istart        => istart,        &
                                  umwm_iend          => iend,          &
                                  umwm_iistart       => iistart,       &
                                  umwm_iiend         => iiend,         &
                                  umwm_ni            => ni,            &
                                  umwm_mi            => mi

    use UMWM_module,        only: umwm_sumt          => sumt,          &
                                  umwm_dtg           => dtg,           &
                                  umwm_dta           => dta,           &
                                  umwm_dtamin        => dtamin

    use UMWM_module,        only: umwm_ustar         => ustar,         &
                                  umwm_cd            => cd,            &
                                  umwm_ht            => ht,            &
                                  umwm_hts           => hts,           &
                                  umwm_htw           => htw,           &
                                  umwm_mwp           => mwp,           &
                                  umwm_mwd           => mwd,           &
                                  umwm_mss           => mss,           &
                                  umwm_mwl           => mwl,           &
                                  umwm_dwd           => dwd,           &
                                  umwm_dwl           => dwl,           &
                                  umwm_dwp           => dwp,           &
                                  umwm_dcp0          => dcp0,          &
                                  umwm_dcg0          => dcg0,          &
                                  umwm_dcp           => dcp,           &
                                  umwm_dcg           => dcg,           &
                                  umwm_momx          => momx,          &
                                  umwm_momy          => momy,          &
                                  umwm_cgmxx         => cgmxx,         &
                                  umwm_cgmxy         => cgmxy,         &
                                  umwm_cgmyy         => cgmyy,         &
                                  umwm_epsx_ocn      => epsx_ocn,      &
                                  umwm_epsy_ocn      => epsy_ocn,      &
                                  umwm_epsx_atm      => epsx_atm,      &
                                  umwm_epsy_atm      => epsy_atm,      &
                                  umwm_taux          => taux,          &
                                  umwm_tauy          => tauy,          &
                                  umwm_taux_form     => taux_form,     &
                                  umwm_tauy_form     => tauy_form,     &
                                  umwm_taux_skin     => taux_skin,     &
                                  umwm_tauy_skin     => tauy_skin      

    use UMWM_init,          only: umwm_environment   => environment,   &
                                  umwm_alloc         => alloc,         &
                                  umwm_grid          => grid,          &
                                  umwm_masks         => masks,         &
                                  umwm_partition     => partition,     &
                                  umwm_remap         => remap,         &
                                  umwm_initialize    => init


    use UMWM_source_functions, only: umwm_s_in       => sin_d12,       &
                                     umwm_s_ds       => sds_d12,       &
                                     umwm_s_nl       => snl_d12,       &
                                     umwm_s_ice      => s_ice
                        
    use UMWM_stress,        only: umwm_stress_       => stress

    use UMWM_physics,       only: umwm_source        => source,        &
                                  umwm_diag          => diag

    use UMWM_advection,     only: umwm_propagation   => propagation,   &
                                  umwm_refraction    => refraction


    use UMWM_stokes,        only: umwm_stokes_drift  => stokes_drift

    use UMWM_util,          only: umwm_dealloc       => dealloc,       &
                                  umwm_remap_mn2i    => remap_mn2i,    &  
                                  umwm_remap_i2mn    => remap_i2mn

    use UMWM_util,          only: sigWaveHeight, meanWavePeriod


    use bl_seaspray_mod,    only: mabl_sea_spray     => online_spray


    use, intrinsic :: ISO_FORTRAN_ENV

    implicit none
    private

    character(len=*), parameter :: UMWM_CONFIG_FILE = 'UMWM.rc'

    real, parameter :: FRACTION_ICE_SUPPRESS_WAVES = 0.8
    real, parameter :: NORTH_POLE_CAP_LATITUDE = 89.0


!   Private state
!   -------------

    type WaveModel_State
        private

        type(ESMF_Config) :: CF             ! Private Config
 
        logical:: verbose   = .false.       ! verbose messages
 
        real    :: dt        = 0.0          ! wave model time step, s
 
        integer :: n_split      = 0         ! number of substeps for time-splitting
        integer :: max_substeps = 0         ! max number of sub-steps
        
 
        integer :: om        = 0            ! number of frequency/wavenumber bins
        integer :: pm        = 0            ! number of direction bins
 
        real    :: fmin      = 0.0          ! lowest  frequency bin, Hz
        real    :: fmax      = 0.0          ! highest frequency bin, Hz
        real    :: fprog     = 0.0          ! highest prognostic frequency bin, Hz
 
        real    :: nu_air    = 0.0          ! kinematic viscosity of air, m2 s-1
        real    :: nu_water  = 0.0          ! kinematic viscosity of water, m2 s-1
        real    :: sfct      = 0.0          ! surface tension, N m-1
        real    :: gustiness = 0.0          ! random wind gustiness factor (should be between 0 and 0.2)
        real    :: dmin      = 0.0          ! depth limiter, m
        real    :: explim    = 0.0          ! exponent limiter (0.69 ~ 100% growth)
        real    :: sin_fac   = 0.0          ! input factor from following winds
        real    :: sin_diss1 = 0.0          ! damping factor from opposing winds
        real    :: sin_diss2 = 0.0          ! damping factor from swell overrunning wind
        real    :: sds_fac   = 0.0          ! breaking dissipation factor
        real    :: sds_power = 0.0          ! saturation spectrum power
        real    :: mss_fac   = 0.0          ! mean-square-slope adjustment to Sds
        real    :: snl_fac   = 0.0          ! wave energy downshifting factor
        real    :: sdt_fac   = 0.0          ! dissipation due to turbulence factor
        real    :: sbf_fac   = 0.0          ! bottom friction coefficient, m s-1
        real    :: sbp_fac   = 0.0          ! bottom percolation coefficient, m s-1

        real    :: fice_lth  = 0.0          ! sea ice fraction - lower threshold for attenuation
        real    :: fice_uth  = 0.0          ! sea ice fraction - upper threshold for attenuation

        real    :: charnock_sf=0.0          ! value of Charnock when wave supported stress is 0

        logical :: stokes    = .true.       ! output Stokes drift velocity fields
        real, pointer, dimension(:) :: depths => null()  ! depths for Stokes diagnostics
 
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

        integer                         :: NUM_FREQUENCY_BINS  ! number of frequency/wavenumber bins
        integer                         :: NUM_DIRECTIONS      ! number of descrete directions 

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

        call ESMF_ConfigLoadFile(self%CF, UMWM_CONFIG_FILE, __RC__)

        call ESMF_ConfigGetAttribute(self%CF, self%verbose,   label='verbose:', default=.false., __RC__)

        call ESMF_ConfigGetAttribute(self%CF, self%om,        label='FREQUENCIES:',    __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%pm,        label='DIRECTIONS:' ,    __RC__)
    
        call ESMF_ConfigGetAttribute(self%CF, self%fmin,      label='MIN_FREQUENCY:',  __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%fmax,      label='MAX_FREQUENCY:',  __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%fprog,     label='MAX_PROGNFREQ:',  __RC__)
    
        call ESMF_ConfigGetAttribute(self%CF, self%n_split,   label='N_SPLIT:',        __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%max_substeps, label='MAX_SUBSTEPS:',__RC__)
    
        call ESMF_ConfigGetAttribute(self%CF, self%fice_lth,  label='SEAICE_LTH:',     __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%fice_uth,  label='SEAICE_UTH:',     __RC__)

        call ESMF_ConfigGetAttribute(self%CF, self%charnock_sf, label='CHARNOCK_SF:',  __RC__)

        call ESMF_ConfigGetAttribute(self%CF, self%nu_air,    label='nu_air:',         __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%nu_water,  label='nu_water:',       __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%sfct,      label='sfct:',           __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%gustiness, label='gustiness:',      __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%dmin,      label='dmin:',           __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%explim,    label='explim:',         __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%sin_fac,   label='sin_fac:',        __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%sin_diss1, label='sin_diss1:',      __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%sin_diss2, label='sin_diss2:',      __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%sds_fac,   label='sds_fac:',        __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%sds_power, label='sds_power:',      __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%mss_fac,   label='mss_fac:',        __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%snl_fac,   label='snl_fac:',        __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%sdt_fac,   label='sdt_fac:',        __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%sbf_fac,   label='sbf_fac:',        __RC__)
        call ESMF_ConfigGetAttribute(self%CF, self%sbp_fac,   label='sbp_fac:',        __RC__)


        ASSERT_(self%om > 0)
        ASSERT_(self%pm > 0) 

!!!     ASSERT_(mod(self%pm,8) /= 0)

        ASSERT_(self%fmin  >  0.0)
        ASSERT_(self%fmax  >  0.0)
        ASSERT_(self%fmax  >  self%fmin)
        ASSERT_(self%fprog >= self%fmin)
        ASSERT_(self%fprog <= self%fmax)



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
        NUM_FREQUENCY_BINS = self%om
        NUM_DIRECTIONS     = self%pm

! Set the state variable specs.
! -----------------------------

!  !INTERNAL STATE:

        call MAPL_AddInternalSpec(GC,                                &
            SHORT_NAME     = 'WM_E',                                 &
            LONG_NAME      = 'sea_surface_wave_energy_spectrum',     &
            UNITS          = 'm4 rad-1',                             &
            DIMS           = MAPL_DimsHorzOnly,                      &
            UNGRIDDED_DIMS = (/NUM_FREQUENCY_BINS, NUM_DIRECTIONS/), &
            VLOCATION      = MAPL_VLocationNone,                     &
            ADD2EXPORT     = .true.,                                 &
            RESTART        = MAPL_RestartOptional,                   &
            DEFAULT        = tiny(0.0),               __RC__)

!       call MAPL_AddInternalSpec(GC,                                &
!           SHORT_NAME     = 'WM_K',                                 &
!           LONG_NAME      = 'sea_surface_wave_wavenumber',          &
!           UNITS          = 'rad m-1',                              &
!           DIMS           = MAPL_DimsHorzOnly,                      &
!           UNGRIDDED_DIMS = (/NUM_FREQUENCY_BINS/),                 &
!           VLOCATION      = MAPL_VLocationNone,                     &
!           ADD2EXPORT     = .true.,                                 &
!           RESTART        = MAPL_RestartOptional,                   &
!           DEFAULT        = 1.0,                    __RC__)

        call MAPL_AddInternalSpec(GC,                                &
            SHORT_NAME     = 'WM_USTAR',                             &
            LONG_NAME      = 'friction_velocity_of_air',             &
            UNITS          = 'm s-1',                                &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            ADD2EXPORT     = .true.,                                 &
            RESTART        = MAPL_RestartOptional,                   &
            DEFAULT        = 0.2,                     __RC__)

! !IMPORT STATE:

!  AGCM -> WM

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

!       call MAPL_AddImportSpec(GC,                                  &
!           SHORT_NAME     = 'UA',                                   &
!           LONG_NAME      = 'surface_eastward_wind',                &
!           UNITS          = 'm s-1',                                &
!           DIMS           = MAPL_DimsHorzOnly,                      &
!           VLOCATION      = MAPL_VLocationNone,                     &
!           RESTART        = MAPL_RestartSkip,       __RC__) 
!  
!       call MAPL_AddImportSpec(GC,                                  &
!           SHORT_NAME     = 'VA',                                   &
!           LONG_NAME      = 'surface_northward_wind',               &
!           UNITS          = 'm s-1',                                &
!           DIMS           = MAPL_DimsHorzOnly,                      &
!           VLOCATION      = MAPL_VLocationNone,                     &
!           RESTART        = MAPL_RestartSkip,       __RC__)

        call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME     = 'RHOS',                                 &
            LONG_NAME      = 'air_density_at_surface',               &
            UNITS          = 'kg m-3',                               &
            DIMS           = MAPL_DimsHorzOnly,                      &
            VLOCATION      = MAPL_VLocationNone,                     &
            DEFAULT        = 1.28,                                   &
            RESTART        = MAPL_RestartOptional,    __RC__)

!        call MAPL_AddImportSpec(GC,                                  &
!            SHORT_NAME     = 'DZ',                                   &
!            LONG_NAME      = 'surface_layer_height',                 &
!            UNITS          = 'm',                                    &
!            DIMS           = MAPL_DimsHorzOnly,                      &
!            VLOCATION      = MAPL_VLocationNone,                     &
!            RESTART        = MAPL_RestartOptional,  __RC__)        

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


!  OGCM -> WM
  
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
        ! UMWM diagnostics
        !

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


! Set the Profiling timers
! ------------------------

        call MAPL_TimerAdd(GC, name='TOTAL'        , __RC__)
        call MAPL_TimerAdd(GC, name='INITIALIZE'   , __RC__)
        call MAPL_TimerAdd(GC, name='RUN'          , __RC__)
        call MAPL_TimerAdd(GC, name='-WM_SET'      , __RC__)
        call MAPL_TimerAdd(GC, name='-WM_INIT'     , __RC__)
        call MAPL_TimerAdd(GC, name='-WM_RUN'      , __RC__)
        call MAPL_TimerAdd(GC, name='--WM_PHYSICS' , __RC__)
        call MAPL_TimerAdd(GC, name='--WM_DYNAMICS', __RC__)
        call MAPL_TimerAdd(GC, name='--WM_ADVECT'  , __RC__)
        call MAPL_TimerAdd(GC, name='--WM_REFRACT' , __RC__)
        call MAPL_TimerAdd(GC, name='--WM_S_IN'    , __RC__)
        call MAPL_TimerAdd(GC, name='--WM_S_DS'    , __RC__)
        call MAPL_TimerAdd(GC, name='--WM_S_NL'    , __RC__)
        call MAPL_TimerAdd(GC, name='--WM_SOURCE'  , __RC__)
        call MAPL_TimerAdd(GC, name='--WM_STRESS'  , __RC__)
        call MAPL_TimerAdd(GC, name='-WM_DIAG'     , __RC__)
        call MAPL_TimerAdd(GC, name='-WM_GET'      , __RC__)
        call MAPL_TimerAdd(GC, name='-SEA_SPRAY' , __RC__)
        call MAPL_TimerAdd(GC, name='FINALIZE'     , __RC__)
        call MAPL_TimerAdd(GC, name='-WM_FINALIZE' , __RC__)

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

      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_Grid)              :: GRID
      type(ESMF_VM)                :: VM

      type(ESMF_Alarm)               :: run_alarm
      type(ESMF_TimeInterval)        :: ring_interval
      real(ESMF_KIND_R8)             :: time_step

      type (WaveModel_State), pointer :: self   ! private internal state
      type (WaveModel_Wrap)           :: wrap

! Local Variables



      integer               :: COMM ! MPI communicator from VM
      integer               :: myPE
      integer               :: nPEs
      integer               :: IM, JM, LM
      integer               :: IM_world, JM_world
      integer               :: I, J, K, L

      real, pointer, dimension(:,:) :: LATS => NULL()
      real, pointer, dimension(:,:) :: LONS => NULL()

      integer :: COUNTS(ESMF_MAXDIM)


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

! Get the grid
! ------------

      call ESMF_GridCompGet( GC, grid=GRID, __RC__)

! Call GenericInitialize
! ----------------------
      call MAPL_MemUtilsWrite(VM, 'GEOSUMWM_GridComp:BeforeGenericInitialize', __RC__)

      call MAPL_GenericInitialize(GC, IMPORT, EXPORT, CLOCK, __RC__)

      call MAPL_MemUtilsWrite(VM, 'GEOSUMWM_GridComp:AfterGenericInitialize', __RC__)

! Get parameters from generic state.
! ----------------------------------

      call MAPL_Get(MAPL, IM=IM, JM=JM, LM=LM, LATS=LATS, LONS=LONS, __RC__)

      call MAPL_GridGet(GRID, globalCellCountPerDim=COUNTS, __RC__)

      IM_world = COUNTS(1)
      JM_world = COUNTS(2)


! Initialize UMWM
! ---------------

      !call MAPL_MemUtilsWrite(VM, 'GEOSUMWM_GridComp:BeforeGEOS_UMWM_INIT', __RC__)

      !call GEOS_UMWM_Initialize( COMM, ... , __RC__)

! Get the time step
! -----------------
      call MAPL_Get(MAPL, RunAlarm=run_alarm, __RC__)
      call ESMF_AlarmGet(run_alarm, ringInterval=ring_interval, __RC__)

      call ESMF_TimeIntervalGet(ring_interval, s_r8=time_step, __RC__)
      self%dt = real(time_step)

      if (MAPL_AM_I_ROOT()) then
          write (*, '(A, I)') 'Wave model (dynamics) time step      = ', nint(self%dt)
          write (*, '(A, I)') 'Wave model (dynamics) frequency bins = ', self%om
          write (*, '(A, I)') 'Wave model (dynamics) directions     = ', self%pm
      end if

      !call MAPL_MemUtilsWrite(VM, 'GEOSUMWM_GridComp:AfterGEOS_UMWM_INIT', __RC__)

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
      character(len=ESMF_MAXSTR)    :: COMP_NAME

! Pointers from Internal state
      real, pointer, dimension(:,:,:,:) :: WM_E     => null()
      real, pointer, dimension(:,:)     :: WM_USTAR => null()

! Pointers from Import state
      real, pointer, dimension(:,:) :: U10N => null()
      real, pointer, dimension(:,:) :: V10N => null()
      real, pointer, dimension(:,:) :: U10M => null()
      real, pointer, dimension(:,:) :: V10M => null()
      real, pointer, dimension(:,:) :: RHOS => null()

      real, pointer, dimension(:,:) :: FRACICE => null()

      real, pointer, dimension(:,:) :: UW => null()
      real, pointer, dimension(:,:) :: VW => null()
      real, pointer, dimension(:,:) :: DW => null()

      real, pointer, dimension(:,:) :: TS     => null()
      real, pointer, dimension(:,:) :: TSKINW => null()

      real, pointer, dimension(:,:) :: PS   => null()
      real, pointer, dimension(:,:) :: Q10M => null()
      real, pointer, dimension(:,:) :: T10M => null()
      real, pointer, dimension(:,:) :: RH2M => null()

      real, pointer, dimension(:,:) :: LHFX => null()
      real, pointer, dimension(:,:) :: SHFX => null()

! Pointers to Export state
      real, pointer, dimension(:,:) :: WM_WIND_10N => null()
      real, pointer, dimension(:,:) :: WM_WIND_10M => null()

      real, pointer, dimension(:,:) :: WM_U10N => null()
      real, pointer, dimension(:,:) :: WM_V10N => null()
      real, pointer, dimension(:,:) :: WM_U10M => null()
      real, pointer, dimension(:,:) :: WM_V10M => null()

      real, pointer, dimension(:,:) :: WM_RHOS => null()

      real, pointer, dimension(:,:) :: WM_UW => null()
      real, pointer, dimension(:,:) :: WM_VW => null()
      real, pointer, dimension(:,:) :: WM_DW => null()

      real, pointer, dimension(:,:) :: WM_FRACICE => null()

      !
      ! UMWM diagnostics
      !

      real, pointer, dimension(:,:) :: WM_NUW  => null()

      real, pointer, dimension(:,:) :: WM_SWH  => null()
      real, pointer, dimension(:,:) :: WM_SWHS => null()
      real, pointer, dimension(:,:) :: WM_SWHW => null()

      real, pointer, dimension(:,:) :: WM_MWP  => null()
      real, pointer, dimension(:,:) :: WM_MWD  => null()
      real, pointer, dimension(:,:) :: WM_MSS  => null()
      real, pointer, dimension(:,:) :: WM_MWL  => null()
      real, pointer, dimension(:,:) :: WM_DWD  => null()
      real, pointer, dimension(:,:) :: WM_DWL  => null()
      real, pointer, dimension(:,:) :: WM_DWP  => null()
      real, pointer, dimension(:,:) :: WM_DCP0 => null()
      real, pointer, dimension(:,:) :: WM_DCG0 => null()
      real, pointer, dimension(:,:) :: WM_DCP  => null()
      real, pointer, dimension(:,:) :: WM_DCG  => null()
  
  
      real, pointer, dimension(:,:) :: WM_Z0    => null()
      real, pointer, dimension(:,:) :: WM_CD    => null()
      real, pointer, dimension(:,:) :: WM_CHARNOCK => null()
  
  
      real, pointer, dimension(:,:) :: WM_TAU   => null()
      real, pointer, dimension(:,:) :: WM_TAUX  => null()
      real, pointer, dimension(:,:) :: WM_TAUY  => null()
  
      real, pointer, dimension(:,:) :: WM_TAU_FORM  => null()
      real, pointer, dimension(:,:) :: WM_TAUX_FORM => null()
      real, pointer, dimension(:,:) :: WM_TAUY_FORM => null()
  
      real, pointer, dimension(:,:) :: WM_TAU_SKIN  => null()
      real, pointer, dimension(:,:) :: WM_TAUX_SKIN => null()
      real, pointer, dimension(:,:) :: WM_TAUY_SKIN => null()

      real, pointer, dimension(:,:) :: WM_EDF       => null()
      real, pointer, dimension(:,:) :: WM_EDFX      => null()
      real, pointer, dimension(:,:) :: WM_EDFY      => null()

      real, pointer, dimension(:,:) :: WM_EGF       => null()
      real, pointer, dimension(:,:) :: WM_EGFX      => null()
      real, pointer, dimension(:,:) :: WM_EGFY      => null()

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


! Local derived type aliases

      type(WaveModel_State), pointer  :: self => null()
      type(WaveModel_Wrap)            :: wrap

      type(MAPL_MetaComp), pointer :: MAPL
      type (ESMF_State)            :: INTERNAL
      type(ESMF_Grid)              :: GRID
      type(ESMF_VM)                :: VM

! Local global variables

      integer :: COUNTS(ESMF_MAXDIM)

! Local Variables

      integer :: IM, JM, LM
      integer :: IM_world, JM_world
      integer :: COMM ! MPI communicator from VM
      integer :: myPE
      integer :: nPEs

      integer :: time_substeps

      integer :: ii, NT
      integer :: i, j
      integer :: o, p
      integer :: isd, ied, jsd, jed  ! halo indexes
      integer :: isc, iec, jsc, jec  ! comp indexes
      real, allocatable, dimension(:,:) :: hvar    ! 2d halo buffer

      real    :: tau_, tau_form_, tau_skin_

      real, pointer, dimension(:,:) :: LATS => NULL()
      real, pointer, dimension(:,:) :: LONS => NULL()

      ! MABL sea spray
      real    :: spray_hss
      real    :: spray_hll
      real    :: spray_hwave
      real    :: spray_cwave
      real    :: spray_p
      real    :: spray_usr
      real    :: spray_massf
      real    :: spray_hs_tot
      real    :: spray_hl_tot
      real    :: spray_usr_new
      real    :: spray_S_bar1
      real    :: spray_z_r
      real    :: spray_omega
      real    :: spray_alpha
      real    :: spray_vfm

      real, parameter :: SPRAY_SOURCE_STRENGTH = 0.4
      real, parameter :: SPRAY_FEEDBACK        = 0.2
 
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

      call MAPL_Get(MAPL, IM=IM, JM=JM, LM=LM, &
                          LATS=LATS, LONS=LONS, &
                          INTERNAL_ESMF_STATE=INTERNAL, __RC__)

      call MAPL_GridGet(GRID, globalCellCountPerDim=COUNTS, __RC__)

      IM_world = COUNTS(1)
      JM_world = COUNTS(2)

! Get layout from the grid
! ------------------------

      call ESMF_VMGetCurrent(VM, __RC__)

      call ESMF_VMGet(VM, mpiCommunicator=COMM, localPet=myPE, petCount=nPEs,  __RC__)

! Get my internal private state
! -----------------------------
      call ESMF_UserCompGetInternalState(GC, 'WaveModel_State', wrap, STATUS)
      VERIFY_(STATUS)

      self => wrap%ptr

! Get pointers to inputs
! ----------------------
      call MAPL_GetPointer(INTERNAL, WM_E,     'WM_E',        __RC__)
      call MAPL_GetPointer(INTERNAL, WM_USTAR, 'WM_USTAR',    __RC__)

      ! AGCM     
      call MAPL_GetPointer(IMPORT, U10N,    'U10N',    __RC__)
      call MAPL_GetPointer(IMPORT, V10N,    'V10N',    __RC__)

      call MAPL_GetPointer(IMPORT, U10M,    'U10M',    __RC__)
      call MAPL_GetPointer(IMPORT, V10M,    'V10M',    __RC__)

      call MAPL_GetPointer(IMPORT, RHOS,    'RHOS',    __RC__)

      call MAPL_GetPointer(IMPORT, FRACICE, 'FRACI',   __RC__)

      call MAPL_GetPointer(IMPORT, TS,      'TS',      __RC__)
      call MAPL_GetPointer(IMPORT, TSKINW,  'TSKINW',  __RC__)

      call MAPL_GetPointer(IMPORT, PS,      'PS',      __RC__)
      call MAPL_GetPointer(IMPORT, Q10M,    'Q10M',    __RC__)
      call MAPL_GetPointer(IMPORT, T10M,    'T10M',    __RC__)
      call MAPL_GetPointer(IMPORT, RH2M,    'RH2M',    __RC__)

      call MAPL_GetPointer(IMPORT, LHFX,    'LHFX',    __RC__)
      call MAPL_GetPointer(IMPORT, SHFX,    'SH',      __RC__)


      ! OGCM
      call MAPL_GetPointer(IMPORT, UW,      'UW',      __RC__)
      call MAPL_GetPointer(IMPORT, VW,      'VW',      __RC__)
      call MAPL_GetPointer(IMPORT, DW,      'DW',      __RC__)


      if (MAPL_AM_I_ROOT() .and. self%verbose) then
          print *, 'DEBUG::UMWM  DW      = ', minval(DW), maxval(DW)
          print *, 'DEBUG::UMWM  VW      = ', minval(VW), maxval(VW)
          print *, 'DEBUG::UMWM  W10N    = ', minval(sqrt(U10N*U10N + V10N*V10N)), maxval(sqrt(U10N*U10N + V10N*V10N))
          print *, 'DEBUG::UMWM  FRACICE = ', minval(FRACICE), maxval(FRACICE)

          print *, 'DEBUG::UMWM   E      = ', minval(WM_E), maxval(WM_E)
          print *, 'DEBUG::UMWM   UST    = ', minval(WM_USTAR), maxval(WM_USTAR)
      end if



! Get pointers from export state
! ------------------------------
      call MAPL_GetPointer(EXPORT, WM_WIND_10N,  'W10N', alloc=.true., __RC__)
      call MAPL_GetPointer(EXPORT, WM_WIND_10M,  'W10M', alloc=.true., __RC__)

      call MAPL_GetPointer(EXPORT, WM_U10N,      'U10N', __RC__)
      call MAPL_GetPointer(EXPORT, WM_V10N,      'V10N', __RC__)
      call MAPL_GetPointer(EXPORT, WM_U10M,      'U10M', __RC__)
      call MAPL_GetPointer(EXPORT, WM_V10M,      'V10M', __RC__)

      call MAPL_GetPointer(EXPORT, WM_RHOS,      'RHOS',     __RC__)

      call MAPL_GetPointer(EXPORT, WM_UW,        'UW',       __RC__)
      call MAPL_GetPointer(EXPORT, WM_VW,        'VW',       __RC__)

      call MAPL_GetPointer(EXPORT, WM_DW,        'DW',        __RC__)

      call MAPL_GetPointer(EXPORT, WM_NUW,       'NUW',       __RC__)
     
      call MAPL_GetPointer(EXPORT, WM_FRACICE,   'FRACI',     __RC__)


      ! wave model diagnostics
      call MAPL_GetPointer(EXPORT, WM_SWH,       'SWH', alloc=.true.,  __RC__)
      call MAPL_GetPointer(EXPORT, WM_SWHS,      'SWHS',      __RC__)
      call MAPL_GetPointer(EXPORT, WM_SWHW,      'SWHW',      __RC__)
   
      call MAPL_GetPointer(EXPORT, WM_MWP,       'MWP',       __RC__)
      call MAPL_GetPointer(EXPORT, WM_MWD,       'MWD',       __RC__)
      call MAPL_GetPointer(EXPORT, WM_MSS,       'MSS',       __RC__)
      call MAPL_GetPointer(EXPORT, WM_MWL,       'MWL',       __RC__)
      call MAPL_GetPointer(EXPORT, WM_DWD,       'DWD',       __RC__)
      call MAPL_GetPointer(EXPORT, WM_DWL,       'DWL',       __RC__)
      call MAPL_GetPointer(EXPORT, WM_DWP,       'DWP',       __RC__)
      call MAPL_GetPointer(EXPORT, WM_DCP0,      'DCP0',      __RC__)
      call MAPL_GetPointer(EXPORT, WM_DCG0,      'DCG0',      __RC__)
      call MAPL_GetPointer(EXPORT, WM_DCP,       'DCP',    alloc=.true.,  __RC__)
      call MAPL_GetPointer(EXPORT, WM_DCG,       'DCG',       __RC__)
  
  
      call MAPL_GetPointer(EXPORT, WM_Z0,        'Z0',        __RC__)
      call MAPL_GetPointer(EXPORT, WM_CD,        'CD',        __RC__)
  
      call MAPL_GetPointer(EXPORT, WM_TAU,       'TAU',       __RC__)
      call MAPL_GetPointer(EXPORT, WM_TAUX,      'TAUX',      __RC__)
      call MAPL_GetPointer(EXPORT, WM_TAUY,      'TAUY',      __RC__)
      call MAPL_GetPointer(EXPORT, WM_TAU_FORM,  'TAU_FORM',  __RC__)
      call MAPL_GetPointer(EXPORT, WM_TAUX_FORM, 'TAUX_FORM', __RC__)
      call MAPL_GetPointer(EXPORT, WM_TAUY_FORM, 'TAUY_FORM', __RC__)
      call MAPL_GetPointer(EXPORT, WM_TAU_SKIN,  'TAU_SKIN',  __RC__)
      call MAPL_GetPointer(EXPORT, WM_TAUX_SKIN, 'TAUX_SKIN', __RC__)
      call MAPL_GetPointer(EXPORT, WM_TAUY_SKIN, 'TAUY_SKIN', __RC__)

      call MAPL_GetPointer(EXPORT, WM_EDF,       'EDF',    alloc=.true., __RC__)
      call MAPL_GetPointer(EXPORT, WM_EDFX,      'EDFX',      __RC__)
      call MAPL_GetPointer(EXPORT, WM_EDFY,      'EDFY',      __RC__)

      call MAPL_GetPointer(EXPORT, WM_EGF,       'EGF',       __RC__)
      call MAPL_GetPointer(EXPORT, WM_EGFX,      'EGFX',      __RC__)
      call MAPL_GetPointer(EXPORT, WM_EGFY,      'EGFY',      __RC__)

      call MAPL_GetPointer(EXPORT, WM_CHARNOCK,  'CHARNOCK', alloc=.true., __RC__)

      ! sea spray diagnostics
      call MAPL_GetPointer(EXPORT, WM_EDFP,      'EDFP',        alloc=.true., __RC__)

      call MAPL_GetPointer(EXPORT, WM_LHFX,      'LHFX',      __RC__)
      call MAPL_GetPointer(EXPORT, WM_SHFX,      'SHFX',      __RC__)

      call MAPL_GetPointer(EXPORT, LHFX_TURB,    'LHFX_TURB' , alloc=.true.,  __RC__)
      call MAPL_GetPointer(EXPORT, SHFX_TURB,    'SHFX_TURB' , alloc=.true.,  __RC__)
      call MAPL_GetPointer(EXPORT, LHFX_TOT,     'LHFX_TOT'  , alloc=.true.,  __RC__)
      call MAPL_GetPointer(EXPORT, SHFX_TOT,     'SHFX_TOT'  , alloc=.true.,  __RC__)
      call MAPL_GetPointer(EXPORT, LHFX_SPRAY,   'LHFX_SPRAY', alloc=.true.,  __RC__)
      call MAPL_GetPointer(EXPORT, SHFX_SPRAY,   'SHFX_SPRAY', alloc=.true.,  __RC__)


! Sanity diagnostics
! ------------------

      ! wind
      if (associated(WM_WIND_10N))   WM_WIND_10N = sqrt(U10N*U10N + V10N*V10N)
      if (associated(WM_WIND_10M))   WM_WIND_10M = sqrt(U10M*U10M + V10M*V10M)

      if (associated(WM_U10N))       WM_U10N = U10N
      if (associated(WM_V10N))       WM_V10N = V10N
      if (associated(WM_U10M))       WM_U10M = U10M
      if (associated(WM_V10M))       WM_V10M = V10M

      ! air density
      if (associated(WM_RHOS))       WM_RHOS = RHOS

      ! surface currents
      if (associated(WM_UW  ))       WM_UW = UW
      if (associated(WM_VW  ))       WM_VW = VW

      ! ocean depth
      if (associated(WM_DW  ))       WM_DW = DW

      ! heat fluxes
      if (associated(WM_LHFX))       WM_LHFX = LHFX
      if (associated(WM_SHFX))       WM_SHFX = SHFX

      ! sea-ice
      if (associated(WM_FRACICE))    WM_FRACICE = FRACICE


! Call UMWM PHYSICS
! -----------------
      if (MAPL_AM_I_Root()) write (OUTPUT_UNIT,*) 'DEBUG::UMWM  WM physics...'


! Wave model physics 
! -----------------------------------------------------
!   0. get ocean currents, etc...

!   ** For tiles not covered by sea-ice:
!   1. Initialize UMWM
!   2. Set UMWM state variables from the INTERNAL state
!   3. Update UMWM forcings from IMPORT state 
!   4. Integrate sources in time
!      4.1 physics
!   5. Wave propagation
!      5.1 advection 
!      5.2 refraction
!   6. Update UMWM diagnostics
!   7. Update INTERNAL state       
!   8. Update EXPORTs
!   9. Free UMWM memory
!
!  99. interface with AGCM


      call MAPL_TimerOn(MAPL, '--WM_SET' )


! 1. Initialize UMWM
! ------------------------------------------------------
      umwm_isGlobal  = .true.
      umwm_restart   = .true.

      umwm_fillEstuaries = .false.
      umwm_fillLakes     = .false.

      umwm_g         = MAPL_GRAV
      umwm_kappa     = MAPL_KARMAN

      umwm_z         = 10.0

      umwm_dtg       = self%dt

      umwm_om        = self%om 
      umwm_pm        = self%pm
      umwm_fmin      = self%fmin
      umwm_fmax      = self%fmax
      umwm_fprog     = self%fprog

      umwm_fice_lth  = self%fice_lth
      umwm_fice_uth  = self%fice_uth

      umwm_nu_air    = self%nu_air
      umwm_nu_water  = self%nu_water
      umwm_sfct      = self%sfct
      umwm_gustiness = self%gustiness
      umwm_dmin      = self%dmin
      umwm_explim    = self%explim
      umwm_sin_fac   = self%sin_fac
      umwm_sin_diss1 = self%sin_diss1
      umwm_sin_diss2 = self%sin_diss2
      umwm_sds_fac   = self%sds_fac
      umwm_sds_power = self%sds_power
      umwm_mss_fac   = self%mss_fac
      umwm_snl_fac   = self%snl_fac
      umwm_sdt_fac   = self%sdt_fac
      umwm_sbf_fac   = self%sbf_fac
      umwm_sbp_fac   = self%sbp_fac


      ! prep halo
      umwm_mm = IM + 2    ! domain size in x
      umwm_nm = JM + 2    ! domain size in y

      ! comp indexes
      isc = lbound(LATS, 1)
      iec = ubound(LATS, 1)
      jsc = lbound(LATS, 2) 
      jec = ubound(LATS, 2)

      ! halo indexes
      isd = isc - 1
      ied = iec + 1
      jsd = jsc - 1
      jed = jec + 1

      NT = (ied-isd+1)*(jed-jsd+1)

      ! 2d buffer for halo transform
      allocate(hvar(isd:ied,jsd:jed), __STAT__)


      call umwm_environment('INIT')

      call umwm_alloc(1)

    
      ! coordinates
      hvar = MAPL_UNDEF
      hvar(isc:iec,jsc:jec) = LONS(:,:)
      call ESMFL_Halo(GRID, hvar, __RC__)
      if (hvar(isd,jsd) == MAPL_UNDEF) hvar(isd,jsd) = 0.5*(hvar(isd+1,jsd) + hvar(isd,jsd+1))
      if (hvar(ied,jed) == MAPL_UNDEF) hvar(ied,jed) = 0.5*(hvar(ied-1,jed) + hvar(ied,jed-1))
      if (hvar(isd,jed) == MAPL_UNDEF) hvar(isd,jed) = 0.5*(hvar(isd,jed-1) + hvar(isd+1,jed))
      if (hvar(ied,jsd) == MAPL_UNDEF) hvar(ied,jsd) = 0.5*(hvar(ied-1,jsd) + hvar(ied,jsd+1))
      umwm_lon = (180.0/MAPL_PI) * hvar

      hvar = MAPL_UNDEF
      hvar(isc:iec,jsc:jec) = LATS(:,:)
      call ESMFL_Halo(GRID, hvar, __RC__)
      if (hvar(isd,jsd) == MAPL_UNDEF) hvar(isd,jsd) = 0.5*(hvar(isd+1,jsd) + hvar(isd,jsd+1))
      if (hvar(ied,jed) == MAPL_UNDEF) hvar(ied,jed) = 0.5*(hvar(ied-1,jed) + hvar(ied,jed-1))
      if (hvar(isd,jed) == MAPL_UNDEF) hvar(isd,jed) = 0.5*(hvar(isd,jed-1) + hvar(isd+1,jed))
      if (hvar(ied,jsd) == MAPL_UNDEF) hvar(ied,jsd) = 0.5*(hvar(ied-1,jsd) + hvar(ied,jsd+1))
      umwm_lat = (180.0/MAPL_PI) * hvar
 

      ! ocean depth
      hvar = MAPL_UNDEF
      hvar(isc:iec,jsc:jec) = DW(:,:)
      call ESMFL_Halo(GRID, hvar, __RC__)
      if (hvar(isd,jsd) == MAPL_UNDEF) hvar(isd,jsd) = 0.5*(hvar(isd+1,jsd) + hvar(isd,jsd+1))
      if (hvar(ied,jed) == MAPL_UNDEF) hvar(ied,jed) = 0.5*(hvar(ied-1,jed) + hvar(ied,jed-1))
      if (hvar(isd,jed) == MAPL_UNDEF) hvar(isd,jed) = 0.5*(hvar(isd,jed-1) + hvar(isd+1,jed))
      if (hvar(ied,jsd) == MAPL_UNDEF) hvar(ied,jsd) = 0.5*(hvar(ied-1,jsd) + hvar(ied,jsd+1))
      umwm_d_2d = hvar


      ! maskout the North pole by treating it as land pints 
      where (umwm_lat > NORTH_POLE_CAP_LATITUDE) umwm_d_2d = tiny(umwm_d_2d)

      ! ...also mask out sea ice
      hvar = MAPL_UNDEF
      hvar(isc:iec,jsc:jec) = FRACICE(:,:)
      call ESMFL_Halo(GRID, hvar, __RC__)

      if (hvar(isd,jsd) == MAPL_UNDEF) hvar(isd,jsd) = 0.5*(hvar(isd+1,jsd) + hvar(isd,jsd+1))
      if (hvar(ied,jed) == MAPL_UNDEF) hvar(ied,jed) = 0.5*(hvar(ied-1,jed) + hvar(ied,jed-1))
      if (hvar(isd,jed) == MAPL_UNDEF) hvar(isd,jed) = 0.5*(hvar(isd,jed-1) + hvar(isd+1,jed))
      if (hvar(ied,jsd) == MAPL_UNDEF) hvar(ied,jsd) = 0.5*(hvar(ied-1,jsd) + hvar(ied,jsd+1))

      where (hvar >=  FRACTION_ICE_SUPPRESS_WAVES) umwm_d_2d = tiny(umwm_d_2d)


      call umwm_grid()             !!! TODO: all arrays are time invariant, no need to do math at every time step
      call umwm_masks()

      call umwm_partition()

      call umwm_alloc(2)
      call umwm_remap()

      ! NOTE - treatment of estuaries and isolated sea points is not
      !        an issue for the WM physics, also it is not active
      !        in the offline runs, so we do not have to do anything
      !        here

      call MAPL_TimerOff(MAPL, '-WM_SET' )



      UMWM_WAVE_MODEL: if (.true.) then !!! umwm_im > 0 .or. .true.) then

          call MAPL_TimerOn(MAPL, '-WM_SET' )

          ! wind speed
          hvar = MAPL_UNDEF
          hvar(isc:iec,jsc:jec) = WM_WIND_10N(:,:)
          call ESMFL_Halo(GRID, hvar, __RC__)
          if (hvar(isd,jsd) == MAPL_UNDEF) hvar(isd,jsd) = 0.5*(hvar(isd+1,jsd) + hvar(isd,jsd+1))
          if (hvar(ied,jed) == MAPL_UNDEF) hvar(ied,jed) = 0.5*(hvar(ied-1,jed) + hvar(ied,jed-1))
          if (hvar(isd,jed) == MAPL_UNDEF) hvar(isd,jed) = 0.5*(hvar(isd,jed-1) + hvar(isd+1,jed))
          if (hvar(ied,jsd) == MAPL_UNDEF) hvar(ied,jsd) = 0.5*(hvar(ied-1,jsd) + hvar(ied,jsd+1))
          umwm_wspd = umwm_remap_mn2i(hvar)


          ! wind direction
          hvar = MAPL_UNDEF
          hvar(isc:iec,jsc:jec) = atan2(V10N, U10N)
          call ESMFL_Halo(GRID, hvar, __RC__)
          if (hvar(isd,jsd) == MAPL_UNDEF) hvar(isd,jsd) = 0.5*(hvar(isd+1,jsd) + hvar(isd,jsd+1))
          if (hvar(ied,jed) == MAPL_UNDEF) hvar(ied,jed) = 0.5*(hvar(ied-1,jed) + hvar(ied,jed-1))
          if (hvar(isd,jed) == MAPL_UNDEF) hvar(isd,jed) = 0.5*(hvar(isd,jed-1) + hvar(isd+1,jed))
          if (hvar(ied,jsd) == MAPL_UNDEF) hvar(ied,jsd) = 0.5*(hvar(ied-1,jsd) + hvar(ied,jsd+1))
          umwm_wdir = umwm_remap_mn2i(hvar)

          ! air density
          hvar = MAPL_UNDEF
          hvar(isc:iec,jsc:jec) = RHOS
          call ESMFL_Halo(GRID, hvar, __RC__)
          if (hvar(isd,jsd) == MAPL_UNDEF) hvar(isd,jsd) = 0.5*(hvar(isd+1,jsd) + hvar(isd,jsd+1))
          if (hvar(ied,jed) == MAPL_UNDEF) hvar(ied,jed) = 0.5*(hvar(ied-1,jed) + hvar(ied,jed-1))
          if (hvar(isd,jed) == MAPL_UNDEF) hvar(isd,jed) = 0.5*(hvar(isd,jed-1) + hvar(isd+1,jed))
          if (hvar(ied,jsd) == MAPL_UNDEF) hvar(ied,jsd) = 0.5*(hvar(ied-1,jsd) + hvar(ied,jsd+1))
          umwm_rhoa = umwm_remap_mn2i(hvar)

          ! water viscosity
          hvar = MAPL_UNDEF
          hvar(isc:iec,jsc:jec) = TS
          call ESMFL_Halo(GRID, hvar, __RC__)
          if (hvar(isd,jsd) == MAPL_UNDEF) hvar(isd,jsd) = 0.5*(hvar(isd+1,jsd) + hvar(isd,jsd+1))
          if (hvar(ied,jed) == MAPL_UNDEF) hvar(ied,jed) = 0.5*(hvar(ied-1,jed) + hvar(ied,jed-1))
          if (hvar(isd,jed) == MAPL_UNDEF) hvar(isd,jed) = 0.5*(hvar(isd,jed-1) + hvar(isd+1,jed))
          if (hvar(ied,jsd) == MAPL_UNDEF) hvar(ied,jsd) = 0.5*(hvar(ied-1,jsd) + hvar(ied,jsd+1))
          ! temperature dependent water viscosity assuming salinity of ...
          hvar = ( -1.04166667e-05*(hvar-273.15)**3 + &
                    1.19142857e-03*(hvar-273.15)**2 + &
                   -5.99654762e-02*(hvar-273.15)    + &
                    1.85338571e+00 ) * 1e-6
          umwm_nu_water_ = umwm_remap_mn2i(hvar)


          ! ocean currents
          hvar = MAPL_UNDEF
          hvar(isc:iec,jsc:jec) = UW(:,:)
          call ESMFL_Halo(GRID, hvar, __RC__)
          where (hvar == MAPL_UNDEF) hvar = 0.0   ! TODO: this is fine for land points, but what about the halo corners?
          umwm_uc = umwm_remap_mn2i(hvar)         !       use the averaged of the two physical-values, similarly to 2D ocean depth?

          hvar = MAPL_UNDEF
          hvar(isc:iec,jsc:jec) = VW(:,:)
          call ESMFL_Halo(GRID, hvar, __RC__)
          where (hvar == MAPL_UNDEF) hvar = 0.0   ! TODO: this is fine for land points, but what about the halo corners?
          umwm_vc = umwm_remap_mn2i(hvar)         !       use the averaged of the two physical-values, similarly to 2D ocean depth?

          ! sea ice
          hvar = MAPL_UNDEF
          hvar(isc:iec,jsc:jec) = FRACICE(:,:)
          call ESMFL_Halo(GRID, hvar, __RC__)
          if (hvar(isd,jsd) == MAPL_UNDEF) hvar(isd,jsd) = 0.5*(hvar(isd+1,jsd) + hvar(isd,jsd+1))
          if (hvar(ied,jed) == MAPL_UNDEF) hvar(ied,jed) = 0.5*(hvar(ied-1,jed) + hvar(ied,jed-1))
          if (hvar(isd,jed) == MAPL_UNDEF) hvar(isd,jed) = 0.5*(hvar(isd,jed-1) + hvar(isd+1,jed))
          if (hvar(ied,jsd) == MAPL_UNDEF) hvar(ied,jsd) = 0.5*(hvar(ied-1,jsd) + hvar(ied,jsd+1))
          umwm_fice = umwm_remap_mn2i(hvar)

          ! water density
          umwm_rhow   = MAPL_RHO_SEAWATER
          umwm_rhow0  = MAPL_RHO_SEAWATER
          umwm_rhorat = umwm_rhoa / umwm_rhow
    
          call MAPL_TimerOff(MAPL, '-WM_SET' )
    
    
          call MAPL_TimerOn(MAPL, '-WM_INIT' )
    
          call umwm_initialize()
    
          call MAPL_TimerOff(MAPL, '-WM_INIT' )


          call MAPL_TimerOn(MAPL, '-WM_SET' )
          ! wave emergy
          do p = 1, umwm_pm
              do o = 1, umwm_om
    
                  hvar = MAPL_UNDEF
                  hvar(isc:iec,jsc:jec) = WM_E(:,:,o,p)
                  call ESMFL_Halo(GRID, hvar, __RC__)
    
                  if (hvar(isd,jsd) == MAPL_UNDEF) hvar(isd,jsd) = 0.5*(hvar(isd+1,jsd) + hvar(isd,jsd+1))
                  if (hvar(ied,jed) == MAPL_UNDEF) hvar(ied,jed) = 0.5*(hvar(ied-1,jed) + hvar(ied,jed-1))
                  if (hvar(isd,jed) == MAPL_UNDEF) hvar(isd,jed) = 0.5*(hvar(isd,jed-1) + hvar(isd+1,jed))
                  if (hvar(ied,jsd) == MAPL_UNDEF) hvar(ied,jsd) = 0.5*(hvar(ied-1,jsd) + hvar(ied,jsd+1))
    
                  ii = 1
                  do j = jsd, jed
                      do i = isd, ied
                          if (umwm_mask(i-isd+1,j-jsd+1) == 1) then
                              umwm_e(o,p,ii) = hvar(i,j)
                              ii = ii + 1
                          end if
                      end do
                  end do
                  
                  ASSERT_(ii-1 == umwm_iend-umwm_istart+1)

              end do
          end do

          ! friction velocity
          hvar = MAPL_UNDEF
          hvar(isc:iec,jsc:jec) = WM_USTAR
          call ESMFL_Halo(GRID, hvar, __RC__)

          if (hvar(isd,jsd) == MAPL_UNDEF) hvar(isd,jsd) = 0.5*(hvar(isd+1,jsd) + hvar(isd,jsd+1))
          if (hvar(ied,jed) == MAPL_UNDEF) hvar(ied,jed) = 0.5*(hvar(ied-1,jed) + hvar(ied,jed-1))
          if (hvar(isd,jed) == MAPL_UNDEF) hvar(isd,jed) = 0.5*(hvar(isd,jed-1) + hvar(isd+1,jed))
          if (hvar(ied,jsd) == MAPL_UNDEF) hvar(ied,jsd) = 0.5*(hvar(ied-1,jsd) + hvar(ied,jsd+1))

          ii = 1
          do j = jsd, jed
              do i = isd, ied
                  if (umwm_mask(i-isd+1,j-jsd+1) == 1) then
                      umwm_ustar(ii) = hvar(i,j)
                      ii = ii + 1
                  end if
              end do
          end do
          
          ASSERT_(ii-1 == umwm_iend-umwm_istart+1)

          call MAPL_TimerOff(MAPL, '-WM_SET' )

    
    
          if (MAPL_AM_I_ROOT() .and. self%verbose) then
              write (*, '(A)'   ) 'UMWM is initialized'
          end if 
    
      
          call MAPL_TimerOn(MAPL, '-WM_RUN' )
          umwm_sumt = 0.0
          umwm_dtamin = self%dt/self%n_split   ! set explicitly the time step
    
          time_substeps = 0

          call MAPL_TimerOn(MAPL,  '--WM_PHYSICS' )

          PHYSICS_TIME_INTEGRATION: do while ((umwm_sumt < umwm_dtg) .and. (time_substeps < self%max_substeps))

#ifdef DEBUG
              if (MAPL_AM_I_ROOT()) write (*, '(A)'   ) 'UMWM integrate source functions...'
#endif

              call MAPL_TimerOn(MAPL,  '--WM_S_IN' )
              call umwm_s_in()           ! compute source input term Sin
              call MAPL_TimerOff(MAPL, '--WM_S_IN' )
    
              call MAPL_TimerOn(MAPL,  '--WM_S_DS' )
              call umwm_s_ds()           ! compute source dissipation term Sds
              call MAPL_TimerOff(MAPL, '--WM_S_DS' )
    
              call MAPL_TimerOn(MAPL,  '--WM_S_NL' )
              call umwm_s_nl()           ! compute non-linear source term Snl
              call MAPL_TimerOff(MAPL, '--WM_S_NL' ) 
   
              call MAPL_TimerOn(MAPL,  '--WM_S_ICE')
              call umwm_s_ice()          ! compute sea ice attenuation term Sice
              call MAPL_TimerOff(MAPL, '--WM_S_ICE') 
   

              call MAPL_TimerOn(MAPL,  '--WM_SOURCE' )
              call umwm_source()         ! integrate source functions
              call MAPL_TimerOff(MAPL, '--WM_SOURCE' ) 

              ! non-negative energy
!!!           where (umwm_e < 0.0) umwm_e = tiny(umwm_e)

              call MAPL_TimerOn(MAPL,  '--WM_STRESS' )
              call umwm_stress_('atm')   ! compute wind stress and drag coefficient
              call umwm_stress_('ocn')   ! compute stress into ocean top and bottom
              call MAPL_TimerOff(MAPL, '--WM_STRESS' )

              time_substeps = time_substeps + 1

#ifdef DEBUG
              if (MAPL_AM_I_ROOT()) write (*, '(F5.3)'   ) umwm_sumt/umwm_dtg
#endif

          end do PHYSICS_TIME_INTEGRATION

          call MAPL_TimerOff(MAPL,  '--WM_PHYSICS' )

#ifdef DEBUG
          if (self%verbose) print *, 'DEBUG::UMWM  time substeps =', time_substeps
#endif

          if (MAPL_AM_I_ROOT() .and. self%verbose) then
              write (*, '(A)'   ) 'UMWM time integration is done for this time step.'
          end if
   

          call MAPL_TimerOn(MAPL, '--WM_DYNAMICS')

          umwm_dtg = self%dt
          umwm_dta = self%dt                      ! TODO: double check
           

!!!       TODO: time integration based on CFL
!!!             use umwm_advection_init::dtamin to split the model time step
!!!             keep in mind that the loop needs to include setting the 1D energy array
!!!       N_timesteps = max(1, self%dt / umwm_dtamin)
!!!       TIME_INTEGRATION: do steps = 1, N_timesteps

#if (0)
          ! update UMWM ef, note that the spatial indexes slice is needed
          umwm_ef(:,:,umwm_istart:umwm_iend) = umwm_e(:,:,umwm_istart:umwm_iend)
#endif   
          ! non-negative energy
!!!       where (umwm_e  < 0.0) umwm_e  = tiny(umwm_e)
!!!       where (umwm_ef < 0.0) umwm_ef = tiny(umwm_ef)

          call umwm_propagation() 
          umwm_e(:,:,umwm_istart:umwm_iend) = umwm_ef(:,:,umwm_istart:umwm_iend)

!!!       where (umwm_e  < 0.0) umwm_e  = tiny(umwm_e)
!!!       where (umwm_ef < 0.0) umwm_ef = tiny(umwm_ef)

          call umwm_refraction() 
          umwm_e(:,:,umwm_istart:umwm_iend) = umwm_ef(:,:,umwm_istart:umwm_iend)

!!!       where (umwm_e  < 0.0) umwm_e  = tiny(umwm_e)
!!!       where (umwm_ef < 0.0) umwm_ef = tiny(umwm_ef)

!!!       TODO
!!!       end do TIME_INTEGRATION

          call MAPL_TimerOff(MAPL, '--WM_DYNAMICS') 



          call MAPL_TimerOn(MAPL,  '--WM_STRESS' )
          call umwm_stress_('atm')   ! compute wind stress and drag coefficient
!!!       call umwm_stress_('ocn')   ! compute stress into ocean top and bottom
          call MAPL_TimerOff(MAPL, '--WM_STRESS' )


          call MAPL_TimerOn(MAPL, '-WM_DIAG')

    !!!   call umwm_stokes_drift()
          call umwm_diag()
    
          call MAPL_TimerOff(MAPL, '-WM_DIAG')


          if (MAPL_AM_I_ROOT() .and. self%verbose) then
              write (*, '(A)'   ) 'UMWM diagnostics is done for this time step.'
          end if


          call MAPL_TimerOn(MAPL, '-WM_GET' )
          ! copy out state variables and diagnostics
        
          do p = 1, umwm_pm
              do o = 1, umwm_om

                  hvar = MAPL_UNDEF
                  ii   = 1
  
                  do j = jsd, jed
                      do i = isd, ied
                          if (umwm_mask(i-isd+1,j-jsd+1) == 1) then
                              hvar(i,j) = umwm_e(o,p,ii)
                              ii = ii + 1
                          else
                              hvar(i,j) = 0.0 
                          end if
                      end do
                  end do
  
                  WM_E(:,:,o,p) = hvar(isc:iec,jsc:jec)
              end do
          end do
  
  
          hvar = MAPL_UNDEF
          ii   = 1
  
          do j = jsd, jed
              do i = isd, ied
                  if (umwm_mask(i-isd+1,j-jsd+1) == 1) then
                      hvar(i,j) = umwm_ustar(ii)
                      ii = ii + 1
                  else
                      hvar(i,j) = 0.0
                  end if
              end do
          end do
  
          WM_USTAR(:,:) = hvar(isc:iec,jsc:jec)
  
      
          call MAPL_TimerOff(MAPL, '-WM_GET' )
 

          if (MAPL_AM_I_ROOT() .and. self%verbose) then
              print *, 'DEBUG::UMWM  _E      = ', minval(WM_E), maxval(WM_E)
              print *, 'DEBUG::UMWM  _UST    = ', minval(WM_USTAR), maxval(WM_USTAR)
          end if



      end if UMWM_WAVE_MODEL


      if (associated(WM_NUW)) then
          call diagnostics(WM_NUW,  umwm_nu_water_, RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_SWH)) then
          call diagnostics(WM_SWH,  umwm_ht,        RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_SWHS)) then
          call diagnostics(WM_SWHS, umwm_hts,       RC=STATUS)
          VERIFY_(STATUS)
      end if
 
      if (associated(WM_SWHW)) then
          call diagnostics(WM_SWHW, umwm_htw,       RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_MWP)) then
          call diagnostics(WM_MWP,  umwm_mwp,       RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_MWD)) then
          call diagnostics(WM_MWD,  umwm_mwd,       RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_MSS)) then
          call diagnostics(WM_MSS,  umwm_mss,       RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_MWL)) then
          call diagnostics(WM_MWL,  umwm_mwl,       RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_DWD)) then
          call diagnostics(WM_DWD,  umwm_dwd,       RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_DWL)) then
          call diagnostics(WM_DWL,  umwm_dwl,       RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_DWP)) then
          call diagnostics(WM_DWP,  umwm_dwp,       RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_DCP0)) then
          call diagnostics(WM_DCP0, umwm_dcp0,      RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_DCG0)) then
          call diagnostics(WM_DCG0, umwm_dcg0,      RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_DCP)) then
          call diagnostics(WM_DCP,  umwm_dcp,       RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_DCG)) then
          call diagnostics(WM_DCG,  umwm_dcg,       RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_CD)) then
          call diagnostics(WM_CD,   umwm_cd,        RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_TAU)) then
          call diagnostics(WM_TAU, sqrt(umwm_taux**2 + umwm_tauy**2), RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_TAU_FORM)) then
          call diagnostics(WM_TAU_FORM, sqrt(umwm_taux_form**2 + &
                                             umwm_tauy_form**2), RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_TAU_SKIN)) then
          call diagnostics(WM_TAU_SKIN, sqrt(umwm_taux_skin**2 + &
                                             umwm_tauy_skin**2), RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_EDF)) then
          call diagnostics(WM_EDF, sqrt(umwm_epsx_ocn**2 + umwm_epsy_ocn**2), RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_EDFX)) then
          call diagnostics(WM_EDFX, umwm_epsx_ocn, RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_EDFY)) then
          call diagnostics(WM_EDFY, umwm_epsy_ocn, RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_EGF)) then
          call diagnostics(WM_EGF, sqrt(umwm_epsx_atm**2 + umwm_epsy_atm**2), RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_EGFX)) then
          call diagnostics(WM_EGFX, umwm_epsx_atm, RC=STATUS)
          VERIFY_(STATUS)
      end if

      if (associated(WM_EGFY)) then
          call diagnostics(WM_EGFY, umwm_epsy_atm, RC=STATUS)
          VERIFY_(STATUS)
      end if


      DIAGNOSTICS_CHARNOCK: if (associated(WM_CHARNOCK)) then
          hvar = MAPL_UNDEF
          ii   = 1

          do j = jsd, jed
              do i = isd, ied
                  if (umwm_mask(i-isd+1,j-jsd+1) == 1) then

                      tau_      = sqrt(umwm_taux(ii)**2 + umwm_tauy(ii)**2)
                      tau_form_ = sqrt(umwm_taux_form(ii)**2 + umwm_tauy_form(ii)**2)
                      tau_skin_ = sqrt(umwm_taux_skin(ii)**2 + umwm_tauy_skin(ii)**2)

                      ! TODO: need to make sure that this behaves well numerically
                      if (tau_ > tiny(tau_)) then
                          if (tau_form_ > 0.9999*tau_) tau_form_ = 0.9999*tau_
                          hvar(i,j) = self%charnock_sf / sqrt(1 - tau_form_/tau_)
                      else
                          hvar(i,j) = 0.0185
                      end if

                      ii = ii + 1
                  else
                      hvar(i,j) = MAPL_UNDEF
                  end if
              end do
          end do

          WM_CHARNOCK(:,:) = hvar(isc:iec,jsc:jec)
          where (FRACICE > FRACTION_ICE_SUPPRESS_WAVES) WM_CHARNOCK = 0.0

          !TODO : temporary workaround until W2A recognizes MAPL_UNDEF
          where (WM_CHARNOCK == MAPL_UNDEF) WM_CHARNOCK = 0.0

      end if DIAGNOSTICS_CHARNOCK


      DIAGNOSTICS_Z0: if (associated(WM_Z0)) then

          WM_Z0 = MAPL_UNDEF

          where (WM_CHARNOCK /= MAPL_UNDEF .and.  WM_USTAR /= MAPL_UNDEF)
!             WM_Z0 = WM_CHARNOCK * WM_USTAR**2 / MAPL_GRAV
              WM_Z0 = (0.11*MAPL_NUAIR)/max(1e-6, WM_USTAR) + WM_CHARNOCK * WM_USTAR**2 / MAPL_GRAV
          end where

          ! TODO: for now default to GEOS:OCEANICEZ0=1e-3 m over sea ice
          where (FRACICE > FRACTION_ICE_SUPPRESS_WAVES) WM_Z0 = 1.0e-3

      end if DIAGNOSTICS_Z0


! Free the memory used by UMWM
! -----------------------------
      call umwm_dealloc() 



! MABL sea spray parameterization, Bao et al, 2011
! ------------------------------------------------
    call MAPL_TimerOn(MAPL, '-SEA_SPRAY')


    PARAMETERIZED_ENERGY_DISSIPATION_FLUX: if (associated(WM_EDFP)) then
         
         WM_EDFP = 1.0*(-0.4+0.25*WM_WIND_10M)  ! cEwave
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

      SHFX_TOT   = MAPL_UNDEF
      LHFX_TOT   = MAPL_UNDEF
      SHFX_TURB  = MAPL_UNDEF
      LHFX_TURB  = MAPL_UNDEF
      SHFX_SPRAY = MAPL_UNDEF
      LHFX_SPRAY = MAPL_UNDEF

      do j = jsc, jec
          do i = isc, iec

          if ((WM_USTAR(i,j) > 1e-2)   .and. &
              (WM_WIND_10M(i,j) > 1.0) .and. &
              (WM_EDF(i,j) > 1e-3)     .and. &
              (TS(i,j) > MAPL_TICE)    .and. &
              (FRACICE(i,j) < FRACTION_ICE_SUPPRESS_WAVES)) then

              spray_hss   = SHFX(i,j)     ! SHWTR
              spray_hll   = LHFX(i,j)     ! HLATWTR
              spray_hwave = max(WM_SWH(i,j), 0.01)
              spray_cwave = max(WM_DCP(i,j), 0.01)
#if(1)
              spray_p     = 3 * WM_EDF(i,j)  / MAPL_RHOWTR   ! the factor 3 is to agree with EDFP
#else
              spray_p     =     WM_EDFP(i,j) / MAPL_RHOWTR
#endif
              spray_usr   = WM_USTAR(i,j)

              
              if (spray_hwave == MAPL_UNDEF) spray_hwave = tiny(spray_hwave)
              if (spray_cwave == MAPL_UNDEF) spray_cwave = tiny(spray_cwave)
              if (spray_p     == MAPL_UNDEF) spray_p     = tiny(spray_p)
              if (spray_usr   == MAPL_UNDEF) spray_usr   = tiny(spray_usr)
              
              call mabl_sea_spray(SPRAY_SOURCE_STRENGTH, &
                                  SPRAY_FEEDBACK,        &
                                  max(WM_WIND_10M(i,j),0.01), &
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

              SHFX_SPRAY(i,j) = (spray_hs_tot - spray_hss) * (1 - FRACICE(i,j))
              LHFX_SPRAY(i,j) = (spray_hl_tot - spray_hll) * (1 - FRACICE(i,j))

              SHFX_TOT(i,j)   = spray_hss + (spray_hs_tot - spray_hss) * (1 - FRACICE(i,j))
              LHFX_TOT(i,j)   = spray_hll + (spray_hl_tot - spray_hll) * (1 - FRACICE(i,j))
          else
              SHFX_TOT(i,j)   = SHFX(i,j)
              LHFX_TOT(i,j)   = LHFX(i,j)

              SHFX_TURB(i,j)  = SHFX(i,j)
              LHFX_TURB(i,j)  = LHFX(i,j)

              SHFX_SPRAY(i,j) = 0.0
              LHFX_SPRAY(i,j) = 0.0
          end if


          if (abs(SHFX_SPRAY(i,j)) > 1e2 .or. abs(LHFX_SPRAY(i,j)) > 1e2) then
              print *, SHFX_SPRAY(i,j),   &
                       LHFX_SPRAY(i,j),   &
                       WM_WIND_10M(i,j),  &
                       TS(i,j)  - 273.15, & 
                       T10M(i,j)- 273.15, &
                       RH2M(i,j) * 0.01,  &  ! s = ...?
                       PS(i,j)   * 0.01,  &
                       WM_SWH(i,j),       &
                       WM_DCP(i,j),       &
                       WM_EDF(i,j)
          end if

          end do
      end do
       
#ifdef DEBUG
      print *, ' *** DEBUG   WM:SH_SPRAY = ', minval(SHFX_SPRAY), maxval(SHFX_SPRAY)
      print *, ' *** DEBUG   WM:LH_SPRAY = ', minval(LHFX_SPRAY), maxval(LHFX_SPRAY)
#endif


   end if DIAGNOSTICS_SPRAY_FLUXES

   call MAPL_TimerOn(MAPL, '-SEA_SPRAY')





! Stop the timers
! ---------------

      call MAPL_TimerOff(MAPL, '-RUN',  __RC__)
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
      call MAPL_TimerOn(MAPL, '-FINALIZE', __RC__)

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

! Call UMWM_FINALIZE
! ------------

      call MAPL_MemUtilsWrite(VM, 'GEOSUMWM_GridComp:BeforeGEOS_UMWM_FINALIZE', __RC__)

      call MAPL_TimerOn(MAPL, '--UMWM_FINALIZE', __RC__)

      !call GEOS_UMWM_FINALIZE( ... , __RC__)

      call MAPL_TimerOff(MAPL, '--UMWM_FINALIZE', __RC__)

      call MAPL_MemUtilsWrite(VM, 'GEOSUMWM_GridComp:AfterGEOS_UMWM_FINALIZE', __RC__)

! Stop the timers
! ---------------

      call MAPL_TimerOff(MAPL, 'FINALIZE', __RC__)
      call MAPL_TimerOff(MAPL, 'TOTAL',    __RC__)

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

!BOP

! ! IROUTINE: DIAGNOSTICS -- maps internal UMWM 1D variables into 2D arrays,
! !                          can be used to update the GEOS diagnostics

! !INTERFACE:

   subroutine diagnostics(var2d, var1d, RC)

      use UMWM_module, only: umwm_mask => mask

      implicit none

! !ARGUMENTS:

      real, dimension(:), intent(in ) :: var1d
      real, dimension(:,:), intent(out) :: var2d

      integer, intent(out) :: RC

! ! DESCRIPTION: 

!EOP

! ErrLog Variables

      character(len=ESMF_MAXSTR) :: Iam
      integer :: STATUS


! Local Variables
      integer :: isc, iec, jsc, jec
      integer :: isd, ied, jsd, jed

      integer :: i, j, ii

      real, allocatable, dimension(:,:) :: hvar


      Iam = 'UMWM::diagnostics()'
      
      ! comp indexes
      isc = lbound(var2d, 1)
      iec = ubound(var2d, 1)
      jsc = lbound(var2d, 2) 
      jec = ubound(var2d, 2)

      ! halo indexes
      isd = isc - 1
      ied = iec + 1
      jsd = jsc - 1
      jed = jec + 1

      allocate(hvar(isd:ied, jsd:jed), __STAT__)

      hvar = MAPL_UNDEF
      ii   = 1

      do j = jsd, jed
          do i = isd, ied
              if (umwm_mask(i-isd+1,j-jsd+1) == 1) then
                  hvar(i,j) = var1d(ii)
                  ii = ii + 1
              else
                  hvar(i,j) = MAPL_UNDEF
              end if
          end do
      end do

      var2d(:,:) = hvar(isc:iec,jsc:jec)

      RETURN_(ESMF_SUCCESS)

  end subroutine diagnostics

end module GEOS_UMWMGridCompMod









#if (0)
   subroutine Put_UMWM_Input (INPUT_NAME, INPUT, GRID, IM_world, JM_world, RC)

      implicit none

      character(len=*), intent(in)        :: INPUT_NAME
      real, pointer, dimension(:,:)       :: INPUT
      type(ESMF_Grid), intent(in)         :: GRID
      integer, intent(in)                 :: IM_world
      integer, intent(in)                 :: JM_world 

      integer, optional, intent(  out)    :: RC     ! Error code:

      ! This contains the INPUT from every PE on a global array
      ! that every PE sees
      real, dimension(IM_world, JM_world) :: tmp_global

      ! These are tmp_global masked with zeroes for UMWM
      real, dimension(IM_world, JM_world) :: tmp_masked

      tmp_global = 0.0
      tmp_masked = 0.0

      ! Gather distributed import into global variables...
      ! --------------------------------------------------

      call ArrayGather(INPUT, tmp_global, GRID, __RC__)

      ! ... and gathered variables only on root, so broadcast
      ! -----------------------------------------------------

      call MAPL_CommsBcast(VM, DATA=tmp_global, &
                           N=size(tmp_global),  &
                           ROOT=0, __RC__)

      ! Wavewatch can have instabilities if we send MAPL_UNDEF
      ! and it is used. So, set MAPL_UNDEF to 0.0 in a masked
      ! variable (we keep the global version set to MAPL_UNDEF
      ! so we can "remask" the output of UMWM)
      ! ------------------------------------------------------

      where (tmp_global == MAPL_UNDEF) 
         tmp_masked = 0.0
      elsewhere 
         tmp_masked = tmp_global
      endwhere

      ! Now put the masked array into the Wavewatch data
      ! ------------------------------------------------

      call GEOS_UMWM_Put(INPUT_NAME, tmp_masked, RC=STATUS)
      VERIFY_(STATUS)

      ! All Done
      ! --------

      RETURN_(ESMF_SUCCESS)

   end subroutine Put_UMWM_Input

   subroutine Get_UMWM_Output (INPUT, OUTPUT, INPUTMASK, GRID, IM_world, JM_world, RC)

      implicit none

      character(len=*), intent(in)        :: INPUT
      real, pointer, dimension(:,:)       :: OUTPUT
      logical, dimension(:,:), intent(in) :: INPUTMASK
      type(ESMF_Grid), intent(in)         :: GRID
      integer, intent(in)                 :: IM_world
      integer, intent(in)                 :: JM_world 


      integer, optional, intent(  out) :: RC     ! Error code:

      ! This contains output from UMWM on the global grid, but only
      ! with that processor's section of the world filled
      real, dimension(IM_world, JM_world) :: tmp_globallocal

      ! This is tmp_globallocal gathered together filling the world
      real, dimension(IM_world, JM_world) :: tmp_global

      tmp_globallocal = 0.0
      tmp_global = 0.0

      ! Fill the exports
      ! ----------------

      call GEOS_UMWM_Get(INPUT, tmp_globallocal, RC=STATUS)
      VERIFY_(STATUS)

      ! W3S2XY only returns each processor's data on the global grid, 
      ! it does not accumulate the data globally. Here we reduce the
      ! data using a MAPL call.
      ! NOTE: For some outputs, this might need to be a Sum rather 
      !       than a Max. Max works here as Z0 is never less than 0.
      !       But if a field could be negative, you'd want to sum
      !       against a 0 general fill.
      ! -------------------------------------------------------------

      call MAPL_CommsAllReduceMax(VM, tmp_globallocal, tmp_global, size(tmp_globallocal), __RC__)

      ! Now scatter global output to each processor correctly
      ! -----------------------------------------------------

      call ArrayScatter(OUTPUT, tmp_global, GRID, __RC__)

      ! Now that the data is scattered, mask out the parts with UNDEF that were
      ! undefined on the input arrays. INPUTMASK is set against the winds at
      ! present.
      ! -----------------------------------------------------------------------

      where (.not. INPUTMASK) OUTPUT = MAPL_UNDEF

      ! All Done
      ! --------

      RETURN_(ESMF_SUCCESS)

   end subroutine Get_UMWM_Output
#endif

