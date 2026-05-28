
!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
module GEOS_SurfaceGridCompMod

!BOP
! !MODULE: GEOS_Surface --- A composite component for the surface components.

! !DESCRIPTION:
!
!   {\tt GEOS\_Surface} is a light-weight gridded component that implements the
!      interface to the tiled surface components. The surface computational components
!      (LAND, LAKE, OCEAN, LANDICE) are its children. All of {\tt GEOS\_Surface}'s imports and exports
!      are in the atmospheric model's grid. In {\tt GEOS\_Surface} these are transformed to the
!      exchange grid, and the relevant portions of the exchange grid are passed to
!      each of the children. The children's results are them replaced in
!      the full exchange grid and transformed back to the atmospheric grid.
!
!      {\tt GEOS\_Surface} has two run stages, as do its children. These are meant
!      to interface with the two stages of {\tt GEOS\_Turbulence}. During the first run
!      stage, the children all produce surface exchange coefficients, and during the
!      second, they update the surface state and produce final values of the fluxes.
!
!      {\tt GEOS\_Surface} keeps a Private Internal State called 'SURF\_state' in the
!      component object. In this state it saves the tranforms between the atmospheric
!      grid and each of the children's exchange grids. This should be done more
!      elegantly once ESMF has exchange grid support. It also has a Internal State
!      that is used to communicate between the two run methods. These internal states
!      do not need to be saved in restarts.
!
!      The four children of {\tt GEOS\_Surface} are given the names:
!      'LAKE', which treats inland freshwater bodies; 'LANDICE', which treats permanent
!      glaciers; 'LAND', which treats all other land surface types, both bare and vegetated,
!      as well as vegetated wetlands not considered freshwater bodies; and  'SALTWATER', which
!      performs the surface calculations for all ocean areas. All four operate in lists
!      of tiles that are nonoverlapping subsets of the exhange grid, and their union---the
!      full exchange grid---tiles the entire sphere.
!
!      By default MAPL\_Generic tries to resolve Imports and Exports among
!      the children; but the children of {\tt GEOS\_Surface} do not talk directly to each other,
!      and all communication between them would need to be performed by {\tt GEOS\_Surface} manipulating
!      their Import and Export states.
!
! !USES:

  use ESMF
  use MAPL
  use MAPL_ESMFFieldBundleRead
  use GEOS_UtilsMod

  use GEOS_LakeGridCompMod,      only : LakeSetServices     => SetServices
  use GEOS_LandiceGridCompMod,   only : LandiceSetServices  => SetServices
  use GEOS_LandGridCompMod,      only : LandSetServices     => SetServices
  use GEOS_SaltwaterGridCompMod, only : OceanSetServices    => SetServices

  use m_mpif90, only: MP_INTEGER, MP_REAL, MP_STATUS_SIZE
  use StieglitzSnow, only : NUM_DUDP, NUM_DUSV, NUM_DUWT, NUM_DUSD, &
                            NUM_BCDP, NUM_BCSV, NUM_BCWT, NUM_BCSD, &
                            NUM_OCDP, NUM_OCSV, NUM_OCWT, NUM_OCSD, &
                            NUM_SUDP, NUM_SUSV, NUM_SUWT, NUM_SUSD, &
                            NUM_SSDP, NUM_SSSV, NUM_SSWT, NUM_SSSD
  use SurfParams,    only : SurfParams_init
  USE CATCH_CONSTANTS, ONLY :                 &
       N_SNOW_LAND      => CATCH_N_SNOW

  implicit none
  private

  real(kind=ESMF_KIND_R8), parameter :: pi = 3.14159265358979323846
  real(kind=ESMF_KIND_R8), parameter :: rad_to_deg = 180.0/pi  ! degree-radian conversion

  type( ESMF_VM ) :: VMG


  ! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

  integer ::        LAKE
  integer ::     LANDICE
  integer ::       OCEAN
  integer ::        LAND

#ifdef AQUA_PLANET
  integer, parameter :: NUM_CHILDREN = 1
#else
  integer, parameter :: NUM_CHILDREN = 4
#endif

  INTEGER            :: catchswim,landicegoswim

  character(len=ESMF_MAXSTR), pointer :: GCNames(:)
  integer                    :: CHILD_MASK(NUM_CHILDREN)
  integer :: DO_OBIO, ATM_CO2
  integer :: DO_WAVES
  integer :: CHOOSEMOSFC
  logical :: DO_GOSWIM
  logical :: DO_FIRE_DANGER
  logical :: DO_DATA_ATM4OCN

! used only when DO_OBIO==1 or ATM_CO2 == ATM_CO2_FOUR
  integer, parameter :: NB_CHOU_UV   = 5 ! Number of UV bands
  integer, parameter :: NB_CHOU_NIR  = 3 ! Number of near-IR bands
  integer, parameter :: NB_CHOU      = NB_CHOU_UV + NB_CHOU_NIR ! Total number of bands
  integer, parameter :: NB_OBIO = 33    !total number of bands for OradBio
  integer, parameter :: ATM_CO2_FOUR = 4
!

  character(len=ESMF_MAXSTR) :: LAND_PARAMS ! land parameter option

  type T_Routing
     integer :: srcTileID, dstTileID,     &
                srcIndex=-1, dstIndex=-1, &
                srcPE=-1, dstPE=-1, SeqIdx=-1
     real    :: weight
  end type T_Routing

  type T_RiverRouting
     type(T_Routing), pointer :: LocalRoutings(:) => NULL()
     integer, pointer :: karray(:)
     integer, pointer :: kdx(:)
     integer, pointer :: BlockSizes(:), displ(:)
  end type T_RiverRouting

! Internal state and its wrapper
! ------------------------------

  type T_SURFACE_STATE
     private
     type (MAPL_LocStreamXFORM)  :: XFORM_IN (NUM_CHILDREN)
     type (MAPL_LocStreamXFORM)  :: XFORM_OUT(NUM_CHILDREN)
     type (T_RiverRouting), pointer   :: RoutingType => NULL()
  end type T_SURFACE_STATE

  type SURF_WRAP
     type (T_SURFACE_STATE), pointer :: PTR
  end type SURF_WRAP

  interface FILLIN_TILE
     module procedure FILLIN_TILE1D
     module procedure FILLIN_TILE2D
  end interface

  interface FILLOUT_TILE
     module procedure FILLOUT_TILE1D
     module procedure FILLOUT_UNGRIDDED
  end interface

  interface MKTILE
     module procedure MKTILE_1D
     module procedure MKTILE_UNGRIDDED
  end interface

   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer,             intent(  OUT) :: RC  ! return code

! !DESCRIPTION: This version uses the GEOS\_GenericSetServices, which in addition
!                to setting default IRF methods, also allocates
!   our instance of a generic state and puts it in the
!   gridded component (GC). Here we override the Initialize and Run methods.
!   The Run method is a two-stage method that implemets the interaction
!   between the 2-stage children representing the various surface types and the 2-stage
!   turbulence run methods.\\
!
!
!   Note that, in addition to its explicit exports,
!   the entire internal state, which is used to communicate between the two run stages,
!   is exported using the ``friendly-to-self'' mechanism.\\
!
! Imports are read-only quantities computed by other gridded components.\\
!
! Note that the turbulence fluxes appearing in the import state are
! the values computed by the first run stage of turbulence using fixed
! surface conditions. The Export versions of these fluxes are the final
! values actually used in the surface budgets. The same applies to some
! of the radiative fluxes, for which the values exported here are those
! actually used in the budget.

!EOP

!=============================================================================
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Locals

    integer                                 :: I
    type (T_SURFACE_STATE), pointer         :: SURF_INTERNAL_STATE
    type (SURF_wrap)                        :: WRAP
    type (MAPL_MetaComp    ), pointer       :: MAPL
    INTEGER                                 :: LSM_CHOICE, DO_CICE_THERMO
    character(len=ESMF_MAXSTR)              :: SURFRC
    type(ESMF_Config)                       :: SCF        ! info from Surface Config File

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    ! Are we running DataAtm?
    !------------------------
    call MAPL_GetResource ( MAPL, DO_DATA_ATM4OCN, Label="USE_DATA_ATM4OCN:" , DEFAULT=.FALSE., RC=STATUS)
    VERIFY_(STATUS)

! Create Surface Config
    call MAPL_GetResource (MAPL, SURFRC, label = 'SURFRC:', default = 'GEOS_SurfaceGridComp.rc', RC=STATUS) ; VERIFY_(STATUS)
    SCF = ESMF_ConfigCreate(rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigLoadFile(SCF,SURFRC,rc=status) ; VERIFY_(STATUS)

    ! Get CHOICE OF Land Surface Model (1:Catch, 2:Catch-CN4.0, 3:Catch-CN4.5)
    ! -------------------------------------------------------
    call MAPL_GetResource    (MAPL, LSM_CHOICE,    label="LSM_CHOICE:",             DEFAULT=1, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource    (MAPL, DO_OBIO,       label="USE_OCEANOBIOGEOCHEM:",   DEFAULT=0, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource    (MAPL, DO_WAVES,      label="USE_WAVES:",              DEFAULT=0, RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource    (SCF,  ATM_CO2,       label='ATM_CO2:',                DEFAULT=2,          __RC__ )
    call MAPL_GetResource    (SCF,  catchswim,     label='N_CONST_LAND4SNWALB:',    DEFAULT=0,          __RC__ )
    call MAPL_GetResource    (SCF,  landicegoswim, label='N_CONST_LANDICE4SNWALB:', DEFAULT=0,          __RC__ )
    if     (LSM_CHOICE.eq.1) then
       call MAPL_GetResource (SCF,  LAND_PARAMS,   label='LAND_PARAMS:',            DEFAULT="NRv7.2",   __RC__ )
    elseif (LSM_CHOICE.eq.2) then
       call MAPL_GetResource (SCF,  LAND_PARAMS,   label='LAND_PARAMS:',            DEFAULT="CN_CLM40",  __RC__ )
!    elseif (LSM_CHOICE.eq.3) then
!       call MAPL_GetResource (SCF,  LAND_PARAMS,   label='LAND_PARAMS:',            DEFAULT="CN_CLM45", __RC__ )
    else
       _ASSERT(.FALSE.,'unknown LSM_CHOICE')
    end if
    call MAPL_GetResource    (SCF,  CHOOSEMOSFC,   label='CHOOSEMOSFC:',            DEFAULT=1,          __RC__ )
    call MAPL_GetResource    (SCF,  DO_FIRE_DANGER,label='FIRE_DANGER:',            DEFAULT=.false.,    __RC__ )

    call ESMF_ConfigDestroy(SCF, __RC__ )

    if ((catchswim/=0) .or. (landicegoswim/=0) .or. (DO_OBIO/=0)) then
       DO_GOSWIM=.true.
    else
       DO_GOSWIM=.false.
    endif

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run1, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2, RC=STATUS )
    VERIFY_(STATUS)

! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( SURF_INTERNAL_STATE, stat=STATUS )
    VERIFY_(STATUS)
    WRAP%PTR => SURF_INTERNAL_STATE

! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'SURF_state',wrap,status )
    VERIFY_(STATUS)


! Set the state variable specs.
! -----------------------------
!BOS
!  !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_pressure',                  &
        UNITS              = 'Pa',                                &
        SHORT_NAME         = 'PS',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_air_temperature',           &
        UNITS              = 'K',                                 &
        SHORT_NAME         = 'TA',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_air_specific_humidity',     &
        UNITS              = 'kg kg-1',                           &
        SHORT_NAME         = 'QA',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_wind_speed',                &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'SPEED',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'UA',                                        &
        LONG_NAME  = 'eastward_wind_bottom_level',                &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'VA',                                        &
        LONG_NAME  = 'northward_wind_bottom_level',               &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_layer_height',              &
        UNITS              = 'm',                                 &
        SHORT_NAME         = 'DZ',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface geopotential height',       &
        UNITS              = 'm+2 s-2',                           &
        SHORT_NAME         = 'PHIS',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'sensible_heat_flux',                &
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'SH',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'eastward_surface_stress_on_air',    &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUX',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'northward_surface_stress_on_air',   &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUY',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'evaporation',                       &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'EVAP',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'dewfall',                           &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'DEWL',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'frostfall',                         &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'FRSL',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        LONG_NAME          = 'derivative_of_sensible_heat_wrt_dry_static_energy',&
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'DSH',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_eastward_surface_stress_wrt_Us', &
        UNITS              = 'N s m-3',                           &
        SHORT_NAME         = 'DFU',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_northward_surface_stress_wrt_Us', &
        UNITS              = 'N s m-3',                           &
        SHORT_NAME         = 'DFV',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_evaporation_wrt_QS',  &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'DEVAP',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_dewfall_wrt_QS',  &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'DDEWL',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_frostfall_wrt_QS',  &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'DFRSL',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

!   Total precip from MOIST (for backwards compatibility when not using PRECIP_FILE)
!   --------------------------------------------------------------------------------
    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME='TPREC',                                       &
         LONG_NAME ='total_precipitation',                         &
         UNITS     ='kg m-2 s-1',                                  &
         DEFAULT   = MAPL_UNDEF,                                   &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
         RESTART   = MAPL_RestartSkip,                  RC=STATUS  )
     VERIFY_(STATUS)

!   Convective precip from MOIST (for backwards compatibility when not using PRECIP_FILE)
!   -------------------------------------------------------------------------------------
    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME='CN_PRCP',                                     &
         LONG_NAME ='convective_precipitation',                    &
         UNITS     ='kg m-2 s-1',                                  &
         DEFAULT   = MAPL_UNDEF,                                   &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
         RESTART   = MAPL_RestartSkip,                  RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'liquid_water_convective_precipitation', &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'PCU',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'liquid_water_large_scale_precipitation', &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'PLS',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'snowfall',                          &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'SNO',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'icefall',                           &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'ICE',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'freezing_rain_fall',                &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'FRZR',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME         = 'DRPARN',                            &
        LONG_NAME          = 'normalized_surface_downwelling_PAR_beam_flux', &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME         = 'DFPARN',                            &
        LONG_NAME          = 'normalized_surface_downwelling_PAR_diffuse_flux', &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         SHORT_NAME         = 'DRNIRN'                      ,&
         LONG_NAME          = 'normalized_surface_downwelling_nir_beam_flux',&
         UNITS              = '1'                           ,&
         DIMS               = MAPL_DimsHorzOnly,             &
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         SHORT_NAME         = 'DFNIRN'                      ,&
         LONG_NAME          = 'normalized_surface_downwelling_nir_diffuse_flux',&
         UNITS              = '1'                           ,&
         DIMS               = MAPL_DimsHorzOnly,             &
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         SHORT_NAME         = 'DRUVRN'                      ,&
         LONG_NAME          = 'normalized_surface_downwelling_uvr_beam_flux',&
         UNITS              = '1'                           ,&
         DIMS               = MAPL_DimsHorzOnly,             &
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         SHORT_NAME         = 'DFUVRN'                      ,&
         LONG_NAME          = 'normalized_surface_downwelling_uvr_diffuse_flux',&
         UNITS              = '1'                           ,&
         DIMS               = MAPL_DimsHorzOnly,             &
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME         = 'LWDNSRF',                           &
        LONG_NAME          = 'surface_absorbed_longwave_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME         = 'ALW',                               &
        LONG_NAME          = 'linearization_of_surface_emitted_longwave_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME         = 'BLW',                               &
        LONG_NAME          = 'linearization_of_surface_emitted_longwave_flux', &
        UNITS              = 'W m-2 K-1',                         &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    if((DO_OBIO/=0).OR. (ATM_CO2 == ATM_CO2_FOUR)) call OBIO_setServices(NB_CHOU, RC)

    if (DO_GOSWIM) then

       call MAPL_AddImportSpec(GC,                              &
            SHORT_NAME         = 'AERO_DP',                           &
            LONG_NAME          = 'aerosol_depositions',                &
            UNITS              = 'kg m-2 s-1',                        &
            DIMS               = MAPL_DimsHorzOnly,                   &
            VLOCATION          = MAPL_VLocationNone,                  &
            DATATYPE           = MAPL_BundleItem,                     &
            RC=STATUS  )
       VERIFY_(STATUS)

    end if

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DTSDT',                                     &
         LONG_NAME  = 'skin_temperature_analysis_tendency',        &
         UNITS      = 'K s-1',                                     &
         RESTART    = MAPL_RestartSkip,                            &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    if (DO_WAVES /= 0) then
       call MAPL_AddImportSpec(GC,                                    &
            SHORT_NAME         = 'CHARNOCK',                          &
            LONG_NAME          = 'charnock_coefficient',              &
            UNITS              = '1',                                 &
            RESTART            = MAPL_RestartOptional,                &
            DEFAULT            = 0.0185,                              &
            DIMS               = MAPL_DimsHorzOnly,                   &
            VLOCATION          = MAPL_VLocationNone,                  &
            RC=STATUS  )
       VERIFY_(STATUS)
    end if

    if (DO_DATA_ATM4OCN) then
       call MAPL_AddImportSpec ( gc,                               &
         SHORT_NAME = 'DISCHARGE',                                 &
         LONG_NAME  = 'river_discharge_at_ocean_points',           &
         UNITS      = 'kg m-2 s-1',                                &
         RESTART    = MAPL_RestartSkip,                            &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
       VERIFY_(STATUS)
    end if

!  !EXPORT STATE:

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_reflectivity_for_visible_beam',   &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVR',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_reflectivity_for_visible_diffuse',&
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVF',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_reflectivity_for_near_infrared_beam', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNR',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_reflectivity_for_near_infrared_diffuse', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNF',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'EMIS',                              &
        LONG_NAME          = 'surface_emissivity',                &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'Z0',                                &
        LONG_NAME          = 'surface_roughness',                 &
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'Z0H',                               &
        LONG_NAME          = 'surface_roughness_for_heat',        &
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'RI',                                &
        LONG_NAME          = 'surface_bulk_richardson_number',    &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'RE',                                &
        LONG_NAME          = 'surface_reynolds_number',           &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRACI',                             &
        LONG_NAME          = 'ice_covered_fraction_of_grid_cell',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'OFRACI',                            &
        LONG_NAME          = 'ice_covered_fraction_of_ocean_area',&
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'QDWL',                              &
        LONG_NAME          = 'surface_liquid_condensate',         &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'QFRL',                              &
        LONG_NAME          = 'surface_ice_condensate',            &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'SHAT',                              &
        LONG_NAME          = 'effective_surface_dry_static_energy',&
        UNITS              = 'm+2 s-2',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELUS',                             &
        LONG_NAME          = 'change_of_surface_eastward_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        FIELD_TYPE         = MAPL_VectorField,                    &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELVS',                             &
        LONG_NAME          = 'change_of_surface_northward_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        FIELD_TYPE         = MAPL_VectorField,                    &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELSS',                             &
        LONG_NAME          = 'change_of_surface_dry_static_energy',&
        UNITS              = 'm+2 s-2',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELTS',                             &
        LONG_NAME          = 'change_of_surface_skin_temperature',&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELQS',                             &
        LONG_NAME          = 'change_of_surface_specific_humidity',&
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DLQLL',                             &
        LONG_NAME          = 'change_of_surface_liquid_condensate',&
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DLQIL',                             &
        LONG_NAME          = 'change_of_surface_frozen_condensate',&
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRLAND',                            &
        LONG_NAME          = 'fraction_of_land',                  &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRLANDICE',                         &
        LONG_NAME          = 'fraction_of_land_ice',              &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRLAKE',                            &
        LONG_NAME          = 'fraction_of_lake',                  &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FROCEAN',                           &
        LONG_NAME          = 'fraction_of_ocean',                 &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'USTAR',                             &
        LONG_NAME          = 'surface_velocity_scale',            &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TSTAR',                             &
        LONG_NAME          = 'surface_temperature_scale',         &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'QSTAR',                             &
        LONG_NAME          = 'surface_moisture_scale',            &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'BSTAR',                             &
        LONG_NAME          = 'surface_buoyancy_scale',            &
        UNITS              = 'm s-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TSOIL1',                            &
        LONG_NAME          = 'soil_temperature_layer_1'          ,&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TSOIL2',                            &
        LONG_NAME          = 'soil_temperature_layer_2'          ,&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TSOIL3',                            &
        LONG_NAME          = 'soil_temperature_layer_3'          ,&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TSOIL4',                            &
        LONG_NAME          = 'soil_temperature_layer_4'          ,&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TSOIL5',                            &
        LONG_NAME          = 'soil_temperature_layer_5'          ,&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TSOIL6',                            &
        LONG_NAME          = 'soil_temperature_layer_6'          ,&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_snow_on_land',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ASNOW'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'downward_heat_flux_into_snow',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'SHSNOW'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'averaged_snow_temperature' ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'AVETSNOW'                  ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_snow_on_land',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSNOW'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_saturated_zone',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSAT'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_unsaturated_zone',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPUNST'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_wilting_zone',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPWLT'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_land_incl_snow',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSURF'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_saturated_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FRSAT'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_unsaturated_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FRUST'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_wilting_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FRWLT'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'snow_mass'                         ,&
        UNITS              = 'kg m-2'                            ,&
        SHORT_NAME         = 'SNOMAS'                            ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'soil_wetness_surface'              ,&
        UNITS              = '1'                                 ,&
        SHORT_NAME         = 'WET1'                              ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'soil_wetness_surface_for_chem'     ,&
        UNITS              = '1'                                 ,&
        SHORT_NAME         = 'WET1_FOR_CHEM'                     ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'soil_wetness_rootzone'             ,&
        UNITS              = '1'                                 ,&
        SHORT_NAME         = 'WET2'                              ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'soil_wetness_profile'              ,&
        UNITS              = '1'                                 ,&
        SHORT_NAME         = 'WET3'                              ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'soil_moisture_surface'             ,&
        UNITS              = 'm3 m-3'                            ,&
        SHORT_NAME         = 'WCSF'                              ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'soil_moisture_rootzone'            ,&
        UNITS              = 'm3 m-3'                            ,&
        SHORT_NAME         = 'WCRZ'                              ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'soil_moisture_profile'             ,&
        UNITS              = 'm3 m-3'                            ,&
        SHORT_NAME         = 'WCPR'                              ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'leaf_area_index'                   ,&
        UNITS              = '1'                                 ,&
        SHORT_NAME         = 'LAI'                               ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'vegetation_greenness_fraction'     ,&
        UNITS              = '1'                                 ,&
        SHORT_NAME         = 'GRN'                               ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'canopy_height'                     ,&
        UNITS              = 'm'                                 ,&
        SHORT_NAME         = 'Z2CH'                              ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'root_length'                       ,&
        UNITS              = 'm m-3'                             ,&
        SHORT_NAME         = 'ROOTL'                             ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'sensible_heat_flux_from_turbulence',&
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'SH',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'eastward_surface_stress',           &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUX',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        FIELD_TYPE         = MAPL_VectorField,                    &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'northward_surface_stress',          &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUY',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        FIELD_TYPE         = MAPL_VectorField,                    &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'evaporation_from_turbulence',       &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'EVAP',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '10-meter_eastward_wind',                                &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'U10M',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       FIELD_TYPE = MAPL_VectorField,                                        &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '10-meter_northward_wind',                               &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'V10M',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       FIELD_TYPE = MAPL_VectorField,                                        &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'equivalent_neutral_10-meter_eastward_wind',             &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'U10N',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       FIELD_TYPE = MAPL_VectorField,                                        &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'equivalent_neutral_10-meter_northward_wind',            &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'V10N',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       FIELD_TYPE = MAPL_VectorField,                                        &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '50-meter_eastward_wind',                                &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'U50M',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       FIELD_TYPE = MAPL_VectorField,                                        &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '50-meter_northward_wind',                               &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'V50M',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       FIELD_TYPE = MAPL_VectorField,                                        &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '10-meter_air_temperature',                              &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'T10M',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '10-meter_specific_humidity',                            &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'Q10M',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '2-meter_eastward_wind',                                 &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'U2M',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       FIELD_TYPE = MAPL_VectorField,                                        &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '2-meter_northward_wind',                                &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'V2M',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       FIELD_TYPE = MAPL_VectorField,                                        &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '2-meter_air_temperature',                               &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'T2M',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '2-meter_specific_humidity',                             &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'Q2M',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_air_temperature',                               &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'TA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_air_specific_humidity',                         &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'QA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_eastward_wind',                                 &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'UA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       FIELD_TYPE = MAPL_VectorField,                                        &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_northward_wind',                                &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'VA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       FIELD_TYPE = MAPL_VectorField,                                        &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       SHORT_NAME = 'GUST',                                                  &
       LONG_NAME  = 'gustiness',                                             &
       UNITS      = 'm s-1',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       SHORT_NAME = 'VENT',                                                  &
       LONG_NAME  = 'surface_ventilation_velocity',                          &
       UNITS      = 'm s-1',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
       LONG_NAME          = 'land(1)_water(0)_ice(2)_flag',      &
       UNITS              = '1',                                 &
       SHORT_NAME         = 'LWI',                               &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'snow_depth_within_snow_covered_area_fraction_on_land'   ,&
        UNITS              = 'm'                                 ,&
        SHORT_NAME         = 'SNOWDP'                            ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'eastward_stress_over_water',&
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUXW'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        FIELD_TYPE         = MAPL_VectorField            ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'northward_stress_over_water',&
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUYW'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        FIELD_TYPE         = MAPL_VectorField            ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'eastward_stress_over_ice',  &
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUXI'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        FIELD_TYPE         = MAPL_VectorField            ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'northward_stress_over_ice',  &
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUYI'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        FIELD_TYPE         = MAPL_VectorField            ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'open_water_upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHWTR'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'sea_ice_upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHICE'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'open_water_latent_energy_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLATWTR'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'sea_ice_latent_energy_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLATICE'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)
     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'open_water_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWNDWTR'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'sea_ice_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWNDICE'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'open_water_net_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWNDWTR'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'sea_ice_net_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWNDICE'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'ocean_snowfall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'SNOWOCN'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'ocean_icefall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'ICEFOCN'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'ocean_snow_and_icefall'    ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'SPTOTOCN'                  ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'ocean_rainfall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'RAINOCN'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  )
     VERIFY_(STATUS)



  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'evaporation'               ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'EVAPOUT'                   ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'sublimation'          ,&
    UNITS              = 'kg m-2 s-1'                     ,&
    SHORT_NAME         = 'SUBLIM'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'upward_sensible_heat_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'SHOUT'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'land_surface_skin_temperature' ,&
    UNITS              = 'K'                             ,&
    SHORT_NAME         = 'LST'                           ,&
    DIMS               = MAPL_DimsHorzOnly               ,&
    VLOCATION          = MAPL_VLocationNone              ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'river_discharge_at_ocean_points',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'DISCHARGE'                 ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'runoff_total_flux'         ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'RUNOFF'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'river_drainage_at_ocean_points',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'DRAINAGE'                  ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'interception_loss_latent_heat_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPINT'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'baresoil_evaporation_latent_heat_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPSOI'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'transpiration_latent_heat_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPVEG'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snowpack_evaporation_latent_heat_flux_on_land',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPICE'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snowpack_evaporation_latent_heat_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPSNO'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil moisture in Upper 10cm'    ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WAT10CM'                   ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'total soil moisture'       ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WATSOI'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil frozen water content' ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'ICESOI'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'baseflow_flux_land'        ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'BASEFLOW'                  ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'overland_runoff_including_throughflow'  ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'RUNSURF'                   ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'EVLAND',                    &
    LONG_NAME          = 'total_evapotranspiration_land',          &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'PRLAND',                    &
    LONG_NAME          = 'Total_precipitation_land',  &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SNOLAND',                   &
    LONG_NAME          = 'snowfall_land',             &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DRPARLAND',                 &
    LONG_NAME          = 'surface_downwelling_PAR_beam_flux', &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DFPARLAND',                 &
    LONG_NAME          = 'surface_downwelling_PAR_diffuse_flux', &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                    &
     SHORT_NAME         = 'LHSNOW',                    &
     LONG_NAME          = 'Latent_heat_flux_snow',     &
     UNITS              = 'W m-2',                     &
     DIMS               = MAPL_DimsHorzOnly,           &
     VLOCATION          = MAPL_VLocationNone,          &
                                            RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                    &
     SHORT_NAME         = 'TCSORIG',                   &
     LONG_NAME          = 'Input_tc_for_snow',         &
     UNITS              = 'K',                         &
     DIMS               = MAPL_DimsHorzOnly,           &
     VLOCATION          = MAPL_VLocationNone,          &
                                            RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                    &
     SHORT_NAME         = 'TPSN1IN',                   &
     LONG_NAME          = 'Input_temp_of_top_snow_lev',&
     UNITS              = 'K',                         &
     DIMS               = MAPL_DimsHorzOnly,           &
     VLOCATION          = MAPL_VLocationNone,          &
                                            RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                    &
     SHORT_NAME         = 'TPSN1OUT',                  &
     LONG_NAME          = 'Output_temp_of_top_snow_lev',&
     UNITS              = 'K',                         &
     DIMS               = MAPL_DimsHorzOnly,           &
     VLOCATION          = MAPL_VLocationNone,          &
                                            RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                    &
     SHORT_NAME         = 'SWNETSNOW',                    &
     LONG_NAME          = 'Net_shortwave_flux_snow',        &
     UNITS              = 'W m-2',                     &
     DIMS               = MAPL_DimsHorzOnly,           &
     VLOCATION          = MAPL_VLocationNone,          &
                                            RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                    &
     SHORT_NAME         = 'LWUPSNOW',                    &
     LONG_NAME          = 'surface_emitted_longwave_flux_snow',         &
     UNITS              = 'W m-2',                     &
     DIMS               = MAPL_DimsHorzOnly,           &
     VLOCATION          = MAPL_VLocationNone,          &
                                            RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                    &
     SHORT_NAME         = 'LWDNSNOW',                    &
     LONG_NAME          = 'surface_absorbed_longwave_flux_snow',         &
     UNITS              = 'W m-2',                     &
     DIMS               = MAPL_DimsHorzOnly,           &
     VLOCATION          = MAPL_VLocationNone,          &
                                            RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                    &
     SHORT_NAME         = 'GHSNOW',                    &
     LONG_NAME          = 'Ground_heating_flux_snow',       &
     UNITS              = 'W m-2',                     &
     DIMS               = MAPL_DimsHorzOnly,           &
     VLOCATION          = MAPL_VLocationNone,          &
                                            RC=STATUS  )
   VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LHLAND',                    &
    LONG_NAME          = 'Latent_heat_flux_land',     &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SHLAND',                    &
    LONG_NAME          = 'Sensible_heat_flux_land',   &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SWLAND',                    &
    LONG_NAME          = 'Net_shortwave_flux_land',   &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SWDOWNLAND',                &
    LONG_NAME          = 'Incident_shortwave_flux_land',   &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LWLAND',                    &
    LONG_NAME          = 'Net_longwave_flux_land',         &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHLAND',                    &
    LONG_NAME          = 'Ground_heating_flux_land',       &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHTSKIN',                   &
    LONG_NAME          = 'Ground_heating_flux_for_skin_temp',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SMLAND',                    &
    LONG_NAME          = 'Snowmelt_flux_land',        &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'QINFIL',                    &
    LONG_NAME          = 'Soil_water_infiltration_rate',        &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TWLAND',                    &
    LONG_NAME          = 'total_water_storage_land',  &
    UNITS              = 'kg m-2',                    &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TELAND',                    &
    LONG_NAME          = 'Total_energy_storage_land', &
    UNITS              = 'J m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TSLAND',                    &
    LONG_NAME          = 'Total_snow_storage_land',   &
    UNITS              = 'kg m-2',                    &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DWLAND',                    &
    LONG_NAME          = 'rate_of_change_of_total_land_water',&
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DHLAND',                    &
    LONG_NAME          = 'rate_of_change_of_total_land_energy',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPLAND',                    &              ! a.k.a. SPSHLAND
    LONG_NAME          = 'Spurious_sensible_heat_flux_land',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPLH',                      &
    LONG_NAME          = 'Spurious_latent_heat_flux_land',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPWATR',                    &              ! a.k.a. SPEVLAND
    LONG_NAME          = 'Spurious_evapotranspiration_flux_land',&
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPSNOW',                    &
    LONG_NAME          = 'Spurious_snow_energy_flux_land',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snowmelt_flux'             ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'SMELT'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_emitted_longwave_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'HLWUP'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    LONG_NAME          = 'surface_net_downward_longwave_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'LWNDSRF'                   ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
    VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    LONG_NAME          = 'surface_net_downward_shortwave_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'SWNDSRF'                   ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
    VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'total_latent_heat_flux_consistent_with_evaporation_from_turbulence'  ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'LHFX'                      ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'ACCUM',                     &
    LONG_NAME          = 'net_ice_accumulation_rate', &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                         &
    SHORT_NAME         = 'ITY',                       &
    LONG_NAME          = 'vegetation_type',           &
    UNITS              = '1',                         &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                         &
    SHORT_NAME         = 'NITY',                      &
    LONG_NAME          = 'NCEP_vegetation_type',      &
    UNITS              = '1',                         &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
       LONG_NAME          = 'liquid_water_convective_precipitation', &
       UNITS              = 'kg m-2 s-1',                        &
       SHORT_NAME         = 'PCU',                               &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
       LONG_NAME          = 'liquid_water_large_scale_precipitation', &
       UNITS              = 'kg m-2 s-1',                        &
       SHORT_NAME         = 'PLS',                               &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'total_precipitation',               &
       UNITS              = 'kg m-2 s-1',                        &
       SHORT_NAME         = 'PRECTOT',                           &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'convective_precipitation',          &
       UNITS              = 'kg m-2 s-1',                        &
       SHORT_NAME         = 'CN_PRCP',                           &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'snowfall',                          &
       UNITS              = 'kg m-2 s-1',                        &
       SHORT_NAME         = 'SNO',                               &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'icefall',                           &
       UNITS              = 'kg m-2 s-1',                        &
       SHORT_NAME         = 'ICE',                               &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'freezing_rain_fall',                &
       UNITS              = 'kg m-2 s-1',                        &
       SHORT_NAME         = 'FRZR',                              &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'TSKINW',                    &
    LONG_NAME          = 'open_water_skin_temperature',&
    UNITS              = 'K'                         ,&
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
    RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'TSKINICE',                    &
    LONG_NAME          = 'sea_ice_skin_temperature',&
    UNITS              = 'K'                         ,&
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
    RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'HICE',                     &
    LONG_NAME          = 'grid_cell_mean_ice_thickness',&
    UNITS              = 'm'                         ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'HSNO',                     &
    LONG_NAME          = 'grid_cell_mean_snow_thickness',&
    UNITS              = 'm'                         ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'FRZMLT',                     &
    LONG_NAME          = 'freezing_melting_potential',&
    UNITS              = 'W m-2'                         ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'TSKINWCICE',                    &
    LONG_NAME          = 'CICE_water_skin_temperature',&
    UNITS              = 'K'                         ,&
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
    RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'ISTSFC',                    &
    LONG_NAME          = 'snow_or_ice_surface_temperature',&
    UNITS              = 'C'                         ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'SSKINW',                   &
    LONG_NAME          = 'sea_skin_layer_salinity',   &
    UNITS              = 'psu'                       ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'MELTT',                     &
    LONG_NAME          = 'top_ice_melt'              ,&
    UNITS              = 'm s-1'                     ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'MELTB',                     &
    LONG_NAME          = 'basal_ice_melt'            ,&
    UNITS              = 'm s-1'                     ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'MELTL',                     &
    LONG_NAME          = 'lateral_ice_melt'          ,&
    UNITS              = 'm s-1'                     ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'MELTS',                     &
    LONG_NAME          = 'snow_melt'            ,&
    UNITS              = 'm s-1'                     ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'FRAZIL',                    &
    LONG_NAME          = 'frazil_ice_growth'         ,&
    UNITS              = 'm s-1'           ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'CONGEL',                    &
    LONG_NAME          = 'congelation_ice_growth'    ,&
    UNITS              = 'm s-1'                     ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'SNOICE',                    &
    LONG_NAME          = 'snow-ice_formation'        ,&
    UNITS              = 'm s-1'                     ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'DAIDTT',                                &
    LONG_NAME          = 'ice_area_tendency_dueto_thermodynamics', &
    UNITS              = '\% day-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'DVIDTT',                                &
    LONG_NAME          = 'ice_volume_tendency_dueto_thermodynamics', &
    UNITS              = 'cm day-1',                                 &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'DAIDTD',                                &
    LONG_NAME          = 'ice_area_tendency_dueto_dynamics', &
    UNITS              = '\% day-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'DVIDTD',                                &
    LONG_NAME          = 'ice_volume_tendency_dueto_dynamics', &
    UNITS              = 'cm day-1',                                 &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'FBOT',                                &
    LONG_NAME          = 'net_downward_heat_flux_from_ice_to_ocean', &
    UNITS              = 'W m-2',                                 &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'HFLUX',                                &
    LONG_NAME          = 'heat_flux_bw_saltwater_ocean',         &
    UNITS              = 'W m-2',                                 &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
       SHORT_NAME         = 'WATERFLUX',                      &
       LONG_NAME          = 'FRESHWATER_flux_bw_saltwater_ocean',    &
       UNITS              = 'kg m-2 s-1',                        &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
       SHORT_NAME         = 'SALTFLUX',                      &
       LONG_NAME          = 'salt_flux_bw_saltwater_ocean',    &
       UNITS              = 'kg m-2 s-1',                        &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'FSWTHRU',                                &
    LONG_NAME          = 'SW_flux_thru_ice_to_ocean' ,&
    UNITS              = 'W m-2',                                 &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'FSWABS',                                &
    LONG_NAME          = 'SW_flux_absorbed_by_skin_layer',        &
    UNITS              = 'W m-2',                                 &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'USTARI',                                &
    LONG_NAME          = 'ice_ocean_friction_velocity',         &
    UNITS              = 'm s-1',                                 &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'FHOCN',                     &
    LONG_NAME          = 'actual_ocean_ice_flux'     ,&
    UNITS              = 'W m-2'                     ,&
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
       LONG_NAME          = 'snow_mass_layer_1',                 &
       UNITS              = 'kg m-2',                            &
       SHORT_NAME         = 'WESNN1',                            &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
       LONG_NAME          = 'snow_mass_layer_2',                 &
       UNITS              = 'kg m-2',                            &
       SHORT_NAME         = 'WESNN2',                            &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
       LONG_NAME          = 'snow_mass_layer_3',                 &
       UNITS              = 'kg m-2',                            &
       SHORT_NAME         = 'WESNN3',                            &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
       LONG_NAME          = 'vegetation_interception_water_storage',      &
       UNITS              = 'kg m-2',                            &
       SHORT_NAME         = 'CAPAC',                             &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
       LONG_NAME          = 'dew_point_temperature_at_2_m',      &
       UNITS              = 'K',                                 &
       SHORT_NAME         = 'T2MDEW',                            &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
       LONG_NAME          = 'wet_bulb_temperature_at_2_m',       &
       UNITS              = 'K',                                 &
       SHORT_NAME         = 'T2MWET',                            &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'near-surface_relative_humidity',    &
       UNITS              = '%',                                 &
       SHORT_NAME         = 'RH2M',                              &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'near-surface_wind_speed',           &
       UNITS              = 'm s-1',                             &
       SHORT_NAME         = 'UU10M',                             &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'depth_of_cool_layer',               &
       UNITS              = 'm',                                 &
       SHORT_NAME         = 'DCOOL',                             &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'depth_at_base_of_warm_layer',       &
       UNITS              = 'm',                                 &
       SHORT_NAME         = 'DWARM',                             &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'temperature_drop_across_cool_layer',&
       UNITS              = 'K',                                 &
       SHORT_NAME         = 'TDROP',                             &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'net_cooling_in_cool_layer',         &
       UNITS              = 'W m-2',                             &
       SHORT_NAME         = 'QCOOL',                             &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'solar_heating_in_cool_layer',       &
       UNITS              = 'W m-2',                             &
       SHORT_NAME         = 'SWCOOL',                            &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'ustar_over_water_layer',            &
       UNITS              = 'm s-1',                             &
       SHORT_NAME         = 'USTARW',                            &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'mean_temperature_of_interface_layer', &
       UNITS              = 'K',                                 &
       SHORT_NAME         = 'TBAR',                              &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'Saunders_parameter',                &
       UNITS              = '1',                                 &
       SHORT_NAME         = 'LCOOL',                             &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'buoyancy_generation_in_cool_layer', &
       UNITS              = 'm+2 s-3',                           &
       SHORT_NAME         = 'BCOOL',                             &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'temperature_at_base_of_cool_layer', &
       UNITS              = 'K',                                 &
       SHORT_NAME         = 'TDEL',                              &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'foundation_temperature_for_interface_layer', &
       UNITS              = 'K',                                 &
       SHORT_NAME         = 'TS_FOUND',                          &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'foundation_salinity_for_interface_layer', &
       UNITS              = 'PSU',                                 &
       SHORT_NAME         = 'SS_FOUND',                          &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'net_heating_in_warm_layer',         &
       UNITS              = 'W m-2',                             &
       SHORT_NAME         = 'QWARM',                             &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'solar_heating_in_warm_layer',       &
       UNITS              = 'W m-2',                             &
       SHORT_NAME         = 'SWWARM',                            &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'Langmuir_number',                   &
       UNITS              = '1',                                 &
       SHORT_NAME         = 'LANGM',                             &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'Similarity_function_in_warm_layer', &
       UNITS              = '1',                                 &
       SHORT_NAME         = 'PHIW',                              &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'relaxation_time_of_TW_to_TS_FOUND', &
       UNITS              = 's',                                 &
       SHORT_NAME         = 'TAUTW',                             &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'Stability_parameter_in_Warm_Layer', &
       UNITS              = '1',                                 &
       SHORT_NAME         = 'ZETA_W',                            &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                    &
       LONG_NAME          = 'departure_of_mean_interface_temperature_from_foundation_temperature', &
       UNITS              = 'K',                                 &
       SHORT_NAME         = 'TWMTF',                             &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

! GOSWIM EXPORTS (from land snow - catchment/catchmentCN)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'dust_mass_in_snow_bin_1'   ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RDU001'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW_LAND/)             ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'dust_mass_in_snow_bin_2'   ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RDU002'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW_LAND/)             ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'dust_mass_in_snow_bin_3'   ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RDU003'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW_LAND/)             ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'dust_mass_in_snow_bin_4'   ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RDU004'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW_LAND/)             ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'dust_mass_in_snow_bin_5'   ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RDU005'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW_LAND/)             ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'hydrophobic_black_carbon_mass_in_snow_bin_1',&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RBC001'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW_LAND/)             ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'hydrophilic_black_carbon_mass_in_snow_bin_2',&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RBC002'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW_LAND/)             ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'hydrophobic_organic_carbon_mass_in_snow_bin_1',&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'ROC001'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW_LAND/)             ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'hydrophilic_organic_carbon_mass_in_snow_bin_2',&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'ROC002'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW_LAND/)             ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_1',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU001'                ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_2',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU002'                ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_3',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU003'                ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_4',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU004'                ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_5',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU005'                ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_black_carbon_mass_flux_from_the_bottom_layer_bin_1',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTBC001'                ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_black_carbon_mass_flux_from_the_bottom_layer_bin_2',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTBC002'                ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_organic_carbon_mass_flux_from_the_bottom_layer_bin_1',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTOC001'                ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_organic_carbon_mass_flux_from_the_bottom_layer_bin_2',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTOC002'                ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'depth_to_water_table_from_surface_in_peat',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'PEATCLSM_WATERLEVEL'               ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'change_in_free_surface_water_reservoir_on_peat',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'PEATCLSM_FSWCHANGE'        ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSOIL1',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZGT1'                     ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSOIL2',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZGT2'                     ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSOIL3',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZGT3'                     ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSOIL4',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZGT4'                     ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSOIL5',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZGT5'                     ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSOIL6',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZGT6'                     ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_PRMC_and_GWETPROF',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZPR'                      ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_RZMC_and_GWETROOT',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZRZ'                      ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_SFMC_and_GWETTOP',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZSF'                      ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSATLAND_TUNSTLAND_and_TWLTLAND',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZTS'                      ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'soil_wilting_point_in_degree_of_saturation_units'  ,&
       UNITS              = '1'                         ,&
       SHORT_NAME         = 'WPWET'                     ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'soil_wilting_point_in_equivalent_mass_of_total_profile_water'  ,&
       UNITS              = 'kg m-2'                         ,&
       SHORT_NAME         = 'WPEMW'                     ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'soil_wilting_point_in_volumetric_units'  ,&
       UNITS              = 'm3 m-3'                    ,&
       SHORT_NAME         = 'WPMC'                      ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'maximum_soil_water_content_above_wilting_point'  ,&
       UNITS              = 'kg m-2'                    ,&
       SHORT_NAME         = 'CDCR2'                     ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'soil_porosity'             ,&
       UNITS              = 'm3 m-3'                    ,&
       SHORT_NAME         = 'POROS'                     ,&
       DIMS               = MAPL_DimsHorzOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
     VERIFY_(STATUS)

  IF(LSM_CHOICE > 1) THEN

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_exposed_leaf-area_index',&
          UNITS              = '1'                         ,&
          SHORT_NAME         = 'CNLAI'                     ,&
          DIMS               = MAPL_DimsHorzOnly,           &
          VLOCATION          = MAPL_VLocationNone,          &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_total_leaf-area_index'  ,&
          UNITS              = '1'                         ,&
          SHORT_NAME         = 'CNTLAI'                    ,&
          DIMS               = MAPL_DimsHorzOnly,           &
          VLOCATION          = MAPL_VLocationNone,          &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_exposed_stem-area_index',&
          UNITS              = '1'                         ,&
          SHORT_NAME         = 'CNSAI'                     ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_total_carbon'           ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'CNTOTC'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_total_vegetation_carbon',&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'CNVEGC'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_total_root_carbon'      ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'CNROOT'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     if (LSM_CHOICE == 3) then
        call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_fine_root_carbon'       ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'CNFROOTC'                  ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
        VERIFY_(STATUS)
     endif

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_net_primary_production' ,&
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'CNNPP'                     ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_gross_primary_production',&
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'CNGPP'                     ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_total_soil_respiration' ,&
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'CNSR'                      ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_net_ecosystem_exchange' ,&
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'CNNEE'                     ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'abstract_C_pool_to_meet_excess_MR_demand' ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'CNXSMR'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_added_to_maintain_positive_C' ,&
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'CNADD'                     ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_carbon_loss_to_fire'    ,&
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'CNLOSS'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'CN_fractional_area_burn_rate' ,&
          UNITS              = 's-1'                       ,&
          SHORT_NAME         = 'CNBURN'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'absorbed_PAR'              ,&
          UNITS              = 'W m-2'                     ,&
          SHORT_NAME         = 'PARABS'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'incident_PAR'              ,&
          UNITS              = 'W m-2'                     ,&
          SHORT_NAME         = 'PARINC'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'saturated_stomatal_conductance' ,&
          UNITS              = 'm s-1'                     ,&
          SHORT_NAME         = 'SCSAT'                     ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'unstressed_stomatal_conductance' ,&
          UNITS              = 'm s-1'                     ,&
          SHORT_NAME         = 'SCUNS'                     ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'transpiration coefficient' ,&
          UNITS              = '1'                         ,&
          SHORT_NAME         = 'BTRANT'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'solar induced fluorescence',&
          UNITS              = 'umol m-2 sm s-1'           ,&
          SHORT_NAME         = 'SIF'                       ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'fire season length'        ,&
          UNITS              = 'days'                      ,&
          SHORT_NAME         = 'CNFSEL'                    ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

  END IF


  if (DO_WAVES /= 0) then
    call MAPL_AddExportSpec(GC,                             &
       LONG_NAME  = 'surface_pressure',                     &
       UNITS      = 'Pa',                                   &
       SHORT_NAME = 'PS',                                   &
       DIMS       = MAPL_DimsHorzOnly,                      &
       VLOCATION  = MAPL_VLocationNone,                     &
       RC=STATUS  )
    VERIFY_(STATUS)
  end if

  if (DO_FIRE_DANGER) then

    ! hourly

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'FFMC',                      &
         LONG_NAME  = 'fine fuel moisture code',   &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'GFMC',                      &
         LONG_NAME  = 'grass fuel moisture code',  &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DMC',                       &
         LONG_NAME  = 'duff moisture code',        &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DC',                        &
         LONG_NAME  = 'drought code',              &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'FWI',                       &
         LONG_NAME  = 'fire weather index',        &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'BUI',                       &
         LONG_NAME  = 'buildup index',             &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ISI',                       &
         LONG_NAME  = 'initial spread index',      &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DSR',                       &
         LONG_NAME  = 'daily severity rating',     &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    ! daily

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'FFMC_DAILY',                &
         LONG_NAME  = 'fine fuel moisture code (daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DMC_DAILY',                 &
         LONG_NAME  = 'duff moisture code (daily)',&
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DC_DAILY',                  &
         LONG_NAME  = 'drought code (daily)',      &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'FWI_DAILY',                 &
         LONG_NAME  = 'fire weather index (daily)',&
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'BUI_DAILY',                 &
         LONG_NAME  = 'buildup index (daily)',     &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ISI_DAILY',                 &
         LONG_NAME  = 'initial spread index (daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DSR_DAILY',                 &
         LONG_NAME  = 'daily severity rating (daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'FFMC_DAILY_',               &
         LONG_NAME  = 'fine fuel moisture code (daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DMC_DAILY_',                &
         LONG_NAME  = 'duff moisture code (daily)',&
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DC_DAILY_',                 &
         LONG_NAME  = 'drought code (daily)',      &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'FWI_DAILY_',                &
         LONG_NAME  = 'fire weather index (daily)',&
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'BUI_DAILY_',                &
         LONG_NAME  = 'buildup index (daily)',     &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ISI_DAILY_',                &
         LONG_NAME  = 'initial spread index (daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DSR_DAILY_',                &
         LONG_NAME  = 'daily severity rating (daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    ! flammability and ignition sources

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'VPD',                       &
         LONG_NAME  = 'vapor pressure deficit',    &
         UNITS      = 'Pa',                        &
         DIMS       = MAPL_DimsHorzOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)
  end if

! !INTERNAL STATE:

!  These are here only because they are passed between run1 and run2.
!  They don't need to be saved in restarts. Note they are all exported
!  by being made friendly to self.
!  Some may be needed by turbulence, but not in a Friendly way; others
!  are only diagnostics.


     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'QS',                                &
        LONG_NAME          = 'surface_specific_humidity',         &
        UNITS              = 'kg kg-1',                           &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'TS',                                &
        LONG_NAME          = 'surface_temperature',          &
        UNITS              = 'K',                                 &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CT',                                &
        LONG_NAME          = 'surface_exchange_coefficient_for_heat', &
        UNITS              = 'kg m-2 s-1',                        &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CQ',                                &
        LONG_NAME          = 'surface_exchange_coefficient_for_moisture', &
        UNITS              = 'kg m-2 s-1',                        &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CM',                                &
        LONG_NAME          = 'surface_exchange_coefficient_for_momentum', &
        UNITS              = 'kg m-2 s-1',                        &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CN',                                &
        LONG_NAME          = 'surface_neutral_drag_coefficient',  &
        UNITS              = '1',                                 &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'THAT',                              &
        LONG_NAME          = 'effective_surface_skin_temperature',&
        UNITS              = 'K',                                 &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'QHAT',                              &
        LONG_NAME          = 'effective_surface_specific_humidity',&
        UNITS              = 'kg kg-1',                           &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'UHAT',                              &
        LONG_NAME          = 'effective_surface_eastward_velocity',&
        UNITS              = 'm s-1',                             &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'VHAT',                              &
        LONG_NAME          = 'effective_surface_northward_velocity',&
        UNITS              = 'm s-1',                             &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        LONG_NAME          = 'air_density_at_surface',            &
        UNITS              = 'kg m-3',                            &
        SHORT_NAME         = 'RHOS',                              &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC                           ,&
        LONG_NAME          = 'zero_plane_displacement_height'    ,&
        UNITS              = 'm'                                 ,&
        SHORT_NAME         = 'D0'                                ,&
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )

     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC                           ,&
        LONG_NAME          = 'discharge_adjustment_factor'       ,&
        UNITS              = '1'                                 ,&
        SHORT_NAME         = 'DISCHARGE_ADJUST'                  ,&
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsTileOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )

     VERIFY_(STATUS)
!EOS

    OCEAN    = MAPL_AddChild(GC, NAME='SALTWATER', SS=OceanSetServices, RC=STATUS)
    VERIFY_(STATUS)
#ifndef AQUA_PLANET
    LAKE     = MAPL_AddChild(GC, NAME='LAKE', SS=LakeSetServices, RC=STATUS)
    VERIFY_(STATUS)
    LANDICE  = MAPL_AddChild(GC, NAME='LANDICE', SS=LandiceSetServices, RC=STATUS)
    VERIFY_(STATUS)
    LAND     = MAPL_AddChild(GC, NAME='LAND', SS=LandSetServices, RC=STATUS)
    VERIFY_(STATUS)
#endif

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_Get(MAPL, GCNAMES = GCNames, RC=STATUS)
    VERIFY_(STATUS)
    _ASSERT(NUM_CHILDREN == size(GCNames),'needs informative message')

    CHILD_MASK(OCEAN  ) = MAPL_OCEAN
#ifndef AQUA_PLANET
    CHILD_MASK(LAKE   ) = MAPL_LAKE
    CHILD_MASK(LANDICE) = MAPL_LANDICE
    CHILD_MASK(LAND   ) = MAPL_LAND
#endif

! By default MAPL_Generic tries to resolve Imports and Exports among
! the children; but our children do not talk to each other, only to us
! --------------------------------------------------------------------

    ! Note; SURFSTATE is only connected between AGCM and OGCM if USE_CICE_Thermo is set
    ! to 2. Otherwise it is not and indeed, SURFSTATE does not seem to exist. As such,
    ! we use the old TerminateImport of all of ocean for less than 2

    call MAPL_GetResource ( MAPL, DO_CICE_THERMO, Label="USE_CICE_Thermo:" , DEFAULT=0, _RC)
    if (DO_CICE_THERMO == 2) then
       call MAPL_TerminateImport    ( GC, SHORT_NAMES=['SURFSTATE'],    &
                                      CHILD_IDS=[OCEAN],  RC=STATUS  )
       VERIFY_(STATUS)
    else
       call MAPL_TerminateImport    ( GC, CHILD = OCEAN,   RC=STATUS  )
       VERIFY_(STATUS)
    endif
#ifndef AQUA_PLANET
    call MAPL_TerminateImport    ( GC, CHILD = LAKE,    RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_TerminateImport    ( GC, CHILD = LANDICE, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_TerminateImport    ( GC, CHILD = LAND,    RC=STATUS  )
    VERIFY_(STATUS)
#endif

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="InitChild"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="LocStreamCreate"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="LocStreamXForm"    ,RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerAdd(GC,    name="-RUN1"   ,RC=STATUS)
    VERIFY_(STATUS)

    do I=1,NUM_CHILDREN
       call MAPL_TimerAdd(GC,    name="--RUN1_"//trim(GCNames(I))  ,RC=STATUS)
       VERIFY_(STATUS)
    end do

    call MAPL_TimerAdd(GC,    name="-RUN2"   ,RC=STATUS)
    VERIFY_(STATUS)

    do I=1,NUM_CHILDREN
       call MAPL_TimerAdd(GC,    name="--RUN2_"//trim(GCNames(I))  ,RC=STATUS)
       VERIFY_(STATUS)
    end do

! Call SetServices for children
!------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

    contains

    subroutine OBIO_setServices(NB_CHOU, RC)

      integer,           intent(   IN) ::  NB_CHOU
      integer, optional, intent(  OUT) ::  RC

      character(len=ESMF_MAXSTR), parameter     :: IAm="OBIO_setServices"
      integer                                   :: STATUS

      call MAPL_AddImportSpec(GC,                              &
           SHORT_NAME         = 'CO2SC',                             &
           LONG_NAME          = 'CO2 Surface Concentration Bin 001', &
           UNITS              = '1e-6',                              &
           DIMS               = MAPL_DimsHorzOnly,                   &
           VLOCATION          = MAPL_VLocationNone,                  &
           RC=STATUS  )
      VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                              &
            SHORT_NAME         = 'DROBIO',                            &
            LONG_NAME          = 'surface_downwelling_shortwave_beam_flux_per_OBIO_band', &
            UNITS              = 'W m-2',                             &
            DIMS               = MAPL_DimsHorzOnly,                   &
            UNGRIDDED_DIMS     = (/NB_OBIO/),                         &
            VLOCATION          = MAPL_VLocationNone,                  &
            RC=STATUS  )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                              &
            SHORT_NAME         = 'DFOBIO',                            &
            LONG_NAME          = 'surface_downwelling_shortwave_diffuse_flux_per_OBIO_band', &
            UNITS              = 'W m-2',                             &
            DIMS               = MAPL_DimsHorzOnly,                   &
            UNGRIDDED_DIMS     = (/NB_OBIO/),                         &
            VLOCATION          = MAPL_VLocationNone,                  &
            RC=STATUS  )
       VERIFY_(STATUS)

      RETURN_(ESMF_SUCCESS)
    end subroutine OBIO_setServices

  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Initialize -- Initialize method for the GEOS Surface component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Initialize method of the Surface Composite Gridded Component.
!   It reads the tiling file that defines the exchange grid and sets-up the
!   location streams for its children. It then does a Generic_Initialize

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp    ), pointer   :: MAPL
    type (MAPL_MetaComp    ), pointer   :: CHILD_MAPL
    type (MAPL_LocStream       )            :: LOCSTREAM
    type (MAPL_LocStream       )            :: EXCH
    type (MAPL_LocStream       )            :: CHILD_LS
    type (ESMF_Grid            )            :: GRID
    type (ESMF_GridComp        ), pointer   :: GCS(:)
    type (ESMF_State           ), pointer   :: GIM(:), GEX(:)
    character(len=ESMF_MAXSTR)              :: TILEFILE
    character(len=ESMF_MAXSTR)              :: ROUTINGFILE
    character(len=ESMF_MAXSTR)              :: DischargeAdjustFile

    type (T_SURFACE_STATE), pointer         :: SURF_INTERNAL_STATE
    type (SURF_wrap)                        :: WRAP
    integer                                 :: I
    real, pointer                           :: FRLAND   (:,:) => NULL()
    real, pointer                           :: FRLANDICE(:,:) => NULL()
    real, pointer                           :: FRLAKE   (:,:) => NULL()
    real, pointer                           :: FROCEAN  (:,:) => NULL()

    real, pointer                           :: PUME(:,:)
    real, pointer                           :: PCME(:,:)
    real, pointer                           :: EVAP(:,:)
    real, pointer                           :: DISCHARGE_ADJUST(:)
    real, allocatable                       :: PUMETILE(:)
    real, allocatable                       :: PCMETILE(:)
    real, allocatable                       :: EVAPTILE(:)
    real, allocatable                       :: PUMEDISTILE(:)
    real, allocatable                       :: PCMEDISTILE(:)

    type (ESMF_Field)                       :: FIELD
    type (ESMF_FieldBundle)                 :: Bundle
    type (ESMF_Time)                        :: ClimateTime
    type (ESMF_State)                       :: INTERNAL
    integer, pointer, dimension(:)          :: TILETYPE       => NULL()
    integer                                 :: NT
    integer                                 :: K
    integer                                 :: userRC, NumInitPhases
    INTEGER                                 :: LSM_CHOICE

    character(len=ESMF_MAXSTR), parameter   :: INITIALIZED_EXPORTS(4) = (/'FROCEAN  ', &
                             'FRLAKE   ', &
                             'FRLAND   ', &
                             'FRLANDICE' /)


!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Initialize"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"INITIALIZE")

! Get my internal private state. This contains the transforms
!  between the exchange grid and the atmos grid.
!-------------------------------------------------------------

    call ESMF_UserCompGetInternalState(gc, 'SURF_state', wrap, status)
    VERIFY_(STATUS)

    SURF_INTERNAL_STATE => WRAP%PTR

! Get the grid
! ------------

    call ESMF_GridCompGet( GC, grid=GRID,  RC=STATUS )
    VERIFY_(STATUS)

! Create the LocStream for the full exchange grid and put it in the state
! -----------------------------------------------------------------------

    call MAPL_TimerOn(MAPL,"LocStreamCreate")

    call MAPL_Get(MAPL, ExchangeGrid=exch, rc=status)
    VERIFY_(STATUS)

    LOCSTREAM = EXCH

    call MAPL_Set (MAPL, LocStream=LOCSTREAM, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOff(MAPL,"LocStreamCreate")

! Get the children's GCS
!-----------------------

    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, RC=STATUS )
    VERIFY_(STATUS)

! Create the children's location streams as subsets of the exhange grid
!----------------------------------------------------------------------

    call MAPL_TimerOn(MAPL,"LocStreamCreate")

    do I = 1, NUM_CHILDREN

       call MAPL_LocStreamCreate(CHILD_LS, LOCSTREAM,                  &
                                 NAME = GCNAMES(I) ,                   &
                                 MASK = (/CHILD_MASK(I)/),             &
                                 RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetObjectFromGC ( GCS(I) ,   CHILD_MAPL,   RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_Set (CHILD_MAPL, LOCSTREAM=CHILD_LS, RC=STATUS )
       VERIFY_(STATUS)

    end do
    call MAPL_TimerOff(MAPL,"LocStreamCreate")

! Call Initialize for every Child
!--------------------------------

    call MAPL_TimerOn(MAPL,"InitChild")
    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerOff(MAPL,"InitChild")

! Create LocStreams Surface_to_child and child_to_surface transforms
!-------------------------------------------------------------------

    call MAPL_Get(MAPL,             &
         TILETYPES = TILETYPE,                   &
         LOCSTREAM = LOCSTREAM,                  &
         INTERNAL_ESMF_STATE = INTERNAL,         &
                                       RC=STATUS )
    VERIFY_(STATUS)


! Static grid exports
!--------------------

    call MAPL_GetPointer(EXPORT,    FRLAND,     'FRLAND', ALLOC=.true.,  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,    FRLAKE,     'FRLAKE', ALLOC=.true.,  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FRLANDICE,  'FRLANDICE', ALLOC=.true.,  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   FROCEAN,    'FROCEAN', ALLOC=.true.,  RC=STATUS)
    VERIFY_(STATUS)

! Fractional areas of each type onthe atmospheric grid, which is the grid
!   attached to the surface locstream
!------------------------------------------------------------------------

    call MAPL_LocStreamFracArea( LOCSTREAM, MAPL_OCEAN  ,  FROCEAN  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamFracArea( LOCSTREAM, MAPL_LAND   ,  FRLAND   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamFracArea( LOCSTREAM, MAPL_LAKE   ,  FRLAKE   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamFracArea( LOCSTREAM, MAPL_LANDICE,  FRLANDICE, RC=STATUS)
    VERIFY_(STATUS)

    FRLANDICE = max(min(FRLANDICE,1.0),0.0)
    FRLAND    = max(min(FRLAND   ,1.0),0.0)
    FRLAKE    = max(min(FRLAKE   ,1.0),0.0)
    FROCEAN   = max(min(FROCEAN  ,1.0),0.0)

! Create transforms to and from the child streams and the surface stream
!   and save them in the surface internal state.
!-----------------------------------------------------------------------

    call MAPL_TimerOn(MAPL,"LocStreamXForm")
    do I = 1, NUM_CHILDREN
       call MAPL_GetObjectFromGC ( GCS(I) ,   CHILD_MAPL,   RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_Get (CHILD_MAPL, LOCSTREAM=CHILD_LS, RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_LocStreamCreateXform ( XFORM=SURF_INTERNAL_STATE%XFORM_IN(I), &
                                        LocStreamOut=CHILD_LS, &
                                        LocStreamIn=LOCSTREAM, &
                                        NAME=GCNAMES(I), &
                                        RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_LocStreamCreateXform ( XFORM=SURF_INTERNAL_STATE%XFORM_OUT(I), &
                                        LocStreamOut=LOCSTREAM, &
                                        LocStreamIn=CHILD_LS, &
                                        NAME=GCNAMES(I), &
                                        MASK_OUT=TILETYPE == CHILD_MASK(I), &
                                        RC=STATUS )
       VERIFY_(STATUS)
    end do
    call MAPL_TimerOff(MAPL,"LocStreamXForm")


! ======================================================================
!ALT: the next section addresses the problem when export variables have been
!     assigned values during Initialize. To prevent "connected" exports
!     being overwritten by DEFAULT in the Import spec in the other component
!     we label them as being "initailized by restart". A better solution
!     would be to move the computation to phase 2 of Initialize and
!     eliminate this section alltogether
! ======================================================================
    DO I = 1, size(INITIALIZED_EXPORTS)
       call ESMF_StateGet(EXPORT,INITIALIZED_EXPORTS(I), FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", &
                              VALUE=MAPL_InitialRestart, RC=STATUS)
       VERIFY_(STATUS)
    END DO

! Init land and snow constants, currently different in Icarus and GEOSldas
!-----------------------------------------------------------------------
    call MAPL_GetResource ( MAPL, LSM_CHOICE, Label="LSM_CHOICE:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    call SurfParams_init(LAND_PARAMS,LSM_CHOICE,RC=STATUS)
    VERIFY_(STATUS)

! Handle river routing (if required)
!-----------------------------------
    call MAPL_GetResource ( MAPL, RoutingFile, Label="ROUTING_FILE:", &
         DEFAULT="", RC=STATUS)
    VERIFY_(STATUS)

    if (RoutingFile /= "") then
       call InitializeRiverRouting(SURF_INTERNAL_STATE%RoutingType, &
            RoutingFile, LocStream, rc=STATUS)
       VERIFY_(STATUS)

       call MAPL_GetResource ( MAPL, DischargeAdjustFile, Label="DISCHARGE_ADJUST_FILE:", &
            DEFAULT="null", RC=STATUS)
       VERIFY_(STATUS)

       if (DischargeAdjustFile /= "null") then

          bundle = ESMF_FieldBundleCreate (NAME='DISCHARGE_ADJUST', RC=STATUS)
          VERIFY_(STATUS)
          call ESMF_FieldBundleSet(bundle, GRID=GRID, RC=STATUS)
          VERIFY_(STATUS)


          call ESMF_TimeSet (ClimateTime, YY=1997, MM=12, DD=15,          &
                                           H=23,    M=30,   S=0, RC=STATUS )

          call MAPL_CFIORead(DischargeAdjustFile, ClimateTime, Bundle, RC=STATUS)
          VERIFY_(STATUS)

          call ESMFL_BundleGetPointerToData(Bundle,'PRECTOT'    , PUME, RC=STATUS)
          VERIFY_(STATUS)
          call ESMFL_BundleGetPointerToData(Bundle,'PRECTOTCORR', PCME, RC=STATUS)
          VERIFY_(STATUS)
          call ESMFL_BundleGetPointerToData(Bundle,'EVLAND',      EVAP, RC=STATUS)
          VERIFY_(STATUS)

          PCME = PCME - EVAP
          PUME = PUME - EVAP

          NT = size(TileType)

          allocate(PUMETILE   (NT), PCMETILE   (NT))
          allocate(PUMEDISTILE(NT), PCMEDISTILE(NT))


          call MAPL_LocStreamTransform( LOCSTREAM, PUMETILE, PUME, RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_LocStreamTransform( LOCSTREAM, PCMETILE, PCME, RC=STATUS)
          VERIFY_(STATUS)

          call RouteRunoff(SURF_INTERNAL_STATE%RoutingType, PUMETILE, PUMEDISTILE, RC=STATUS)
          VERIFY_(STATUS)
          call RouteRunoff(SURF_INTERNAL_STATE%RoutingType, PCMETILE, PCMEDISTILE, RC=STATUS)
          VERIFY_(STATUS)

          call MAPL_GetPointer(INTERNAL, DISCHARGE_ADJUST, 'DISCHARGE_ADJUST',  RC=STATUS)
          VERIFY_(STATUS)

          DISCHARGE_ADJUST = 1.0


          where(PCMEDISTILE /= 0.0 .and. PUMEDISTILE /= 0.0 .and. TileType == MAPL_OCEAN) &
             DISCHARGE_ADJUST = PUMEDISTILE/PCMEDISTILE

          deallocate(PUMETILE   , PCMETILE   )
          deallocate(PUMEDISTILE, PCMEDISTILE)

          ! Destroy bundle
          call MAPL_FieldBundleDestroy(bundle, rc=status)
          VERIFY_(STATUS)

       end if
    end if



! All Done
!---------

    call MAPL_TimerOff(MAPL,"INITIALIZE")

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

  subroutine InitializeRiverRouting(RoutingType, RoutingFile, Stream, rc)
    type(T_RiverRouting), pointer    :: RoutingType
    character(len=*),        intent(IN) :: RoutingFile
    type(MAPL_LocStream), intent(IN) :: Stream
    integer, optional,    intent(OUT):: rc

    type(T_Routing), pointer         :: LocalRoutings(:) => NULL()
    integer, pointer :: karray(:)
    integer, pointer :: kdx(:)
    integer, pointer :: BlockSizes(:), displ(:)

    type(ESMF_VM)            :: VM
    integer                  :: comm, nDEs, myPE, i, numRoutings
    type(T_Routing), pointer :: Routing
    type(T_Routing), pointer :: tmpLocalRoutings(:)
    integer, pointer         :: TmpActive(:,:)
    integer, pointer         :: Active(:,:)
    integer, pointer         :: ActiveGlobal(:,:)
    integer                  :: numActive, numLocalRoutings
    integer, pointer         :: Local_Id(:)
    integer                  :: unit
    integer :: ksum, n, nsdx
    integer :: ntotal, k
    integer, allocatable :: kn(:), kseq(:), tmparray(:)
#ifdef DEBUG
    character(len=ESMF_MAXSTR)  :: routefile
#endif

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm="InitializeRiverRouting"
    integer                                 :: STATUS

    call ESMF_VMGetCurrent(VM,                                RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_VMGet       (VM,       mpiCommunicator =comm,   RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_VMGet       (VM, localpet=MYPE, petcount=nDEs,  RC=STATUS)
    VERIFY_(STATUS)

! Open the trn file and read the number of "Routings"
!  or land-ocean tile pairs that exchange runoff

    UNIT = GETFILE(RoutingFile, RC=status)
    VERIFY_(STATUS)

    if ( MAPL_am_I_root(vm) ) then
       read(unit, iostat=status) numRoutings
       VERIFY_(STATUS)
    end if

    call MAPL_CommsBcast(vm, DATA=numRoutings, N=1, ROOT=0, RC=status)
    VERIFY_(STATUS)

! Allocate a list that will hold the maximum possible number
!  of local routings. The list will compacted later
!  to the exact number.

    allocate(tmpLocalRoutings(numRoutings))
    allocate(tmpActive(2,numRoutings))

! Read routings from the file

!ALT: this section should be done by READ_PARALLEL when overloaded
    if ( MAPL_am_I_root(vm) ) read(unit) tmpLocalRoutings%srcTileID
    call MAPL_CommsBcast(vm, DATA=tmpLocalRoutings%srcTileID, N=size(tmpLocalRoutings%srcTileID), ROOT=0, RC=status)
    VERIFY_(STATUS)
    if ( MAPL_am_I_root(vm) ) read(unit) tmpLocalRoutings%dstTileID
    call MAPL_CommsBcast(vm, DATA=tmpLocalRoutings%dstTileID, N=size(tmpLocalRoutings%dstTileID), ROOT=0, RC=status)
    VERIFY_(STATUS)
    if ( MAPL_am_I_root(vm) ) read(unit) tmpLocalRoutings%weight
    call MAPL_CommsBcast(vm, DATA=tmpLocalRoutings%weight, N=size(tmpLocalRoutings%weight), ROOT=0, RC=status)
    VERIFY_(STATUS)

    call FREE_FILE(unit, RC=STATUS)
    VERIFY_(STATUS)

    where(tmpLocalRoutings%dstTileID==0) tmpLocalRoutings%dstTileID = tmpLocalRoutings%srcTileID

!  and do some initial processing.

    numLocalRoutings= 0
    numActive       = 0

    call MAPL_LocStreamGet(Stream, LOCAL_ID=Local_Id, rc=status)
    VERIFY_(STATUS)
    do i=1,numRoutings

! For routings with a src or dst tile in the local PE,
!  convert the tile's IDs to local indeces.
!  Tiles that are not in the local processor are
!  assigned an index of -1.

       Routing => tmpLocalRoutings(i)
       Routing%seqIdx = i

       call Tile2Index(Routing, Local_Id)

       if(Routing%srcIndex>0 .and. Routing%dstIndex>0) then
          Routing%srcPE = myPE
          Routing%dstPE = myPE
       endif

! If either the routing's source or destination tile
!  was found to be in the local processor, add the routing
!  to the list of local routings.

       if(Routing%srcIndex>0 .or. Routing%dstIndex>0) then
          numLocalRoutings = numLocalRoutings + 1
          TmpLocalRoutings(numLocalRoutings) = Routing

! Tiles in the local PE that will be involved in communication are
!  places in an "active" list.

          if(Routing%srcIndex<0 .or. Routing%dstIndex<0) then
             numActive = numActive + 1
             if(Routing%srcIndex>0) then
                TmpActive(1,numActive) = Routing%srcTileId
                TmpActive(2,numActive) = Routing%srcIndex
             else
                TmpActive(1,numActive) = Routing%dstTileId
                TmpActive(2,numActive) = Routing%dstIndex
             end if
          end if
       endif
    enddo

! Compact the lists of local and active routings

    if(associated(LocalRoutings)) deallocate(LocalRoutings)
    allocate(LocalRoutings(numLocalRoutings))
    LocalRoutings = tmpLocalRoutings(:numLocalRoutings)
    deallocate(tmpLocalRoutings)

    allocate(Active(3,numActive))
    Active(1:2,:) = tmpActive(:,:numActive)
    Active(3,:)   = myPE
    deallocate(tmpActive)

! We now create a global version of the Active list to
!  figure out the communication pattern

    allocate(BlockSizes(0:nDEs-1),displ(0:nDEs))

    call MPI_AllGather(numActive, 1, MP_INTEGER, &
                       BlockSizes,1, MP_INTEGER, &
                       comm,status)
    VERIFY_(STATUS)

    BlockSizes = BlockSizes*3
    displ(0)=0
    do i=1,nDEs
       displ(i) = displ(i-1) + BlockSizes(i-1)
    enddo

    allocate(ActiveGlobal(3,displ(nDEs)/3))

    call MPI_AllGatherV(Active      ,size(Active), MP_INTEGER, &
                        ActiveGlobal,blocksizes, displ, MP_INTEGER, &
                        comm,status)
    VERIFY_(STATUS)

! Using the global list, we now visit all the local tiles that will be
!   active in communication and find the PE and Index they
!   are sending to or receiving from. This is sufficient to do the
!   routing if we do an mpi send/recv pair for every Active routing.

    do i=1,numLocalRoutings
       if(LocalRoutings(i)%srcIndex>0 .neqv. LocalRoutings(i)%dstIndex>0) then
          if    (LocalRoutings(i)%srcIndex>0) then
             call FindTile(ActiveGlobal, LocalRoutings(i)%dstTileID,  &
                           LocalRoutings(i)%dstPE, LocalRoutings(i)%dstIndex)
             LocalRoutings(i)%srcPE = myPE
          else
             call FindTile(ActiveGlobal, LocalRoutings(i)%srcTileID,  &
                           LocalRoutings(i)%srcPE, LocalRoutings(i)%srcIndex)
             LocalRoutings(i)%dstPE = myPE
          end if
       elseif(LocalRoutings(i)%srcIndex>0 .and. LocalRoutings(i)%dstIndex>0) then
          LocalRoutings(i)%srcPE = myPE
          LocalRoutings(i)%dstPE = myPE
       else
          _ASSERT(.FALSE.,'needs informative message')
       end if
    end do

    deallocate(BlockSizes,displ,Active,ActiveGlobal)

#ifdef DEBUG
    write(routefile,"(A6,I2.2)") 'route.', mype
    unit = 7
    open(unit=unit, file=routeFile, form='formatted', iostat=status)
    VERIFY_(STATUS)
    do i=1,min(numLocalRoutings,100)
       write(unit, '(6I10,F10.3)')  &
            LocalRoutings(i)%srcTileID, &
            LocalRoutings(i)%dstTileID, &
            LocalRoutings(i)%srcIndex,  &
            LocalRoutings(i)%dstIndex,  &
            LocalRoutings(i)%srcPE,     &
            LocalRoutings(i)%dstPE,     &
            LocalRoutings(i)%weight
    end do
    close(unit)

#endif

    !ALT NEW ROUTING to make communication more effective

    nsdx=0
    do i=1,numLocalRoutings
       !notneeded if (mype == Routing(i)%dstPE) nddx = nddx+1
       if (mype == LocalRoutings(i)%srcPE) nsdx = nsdx+1
    end do
    allocate(kdx(nsdx), blocksizes(nDEs), _STAT)
    blocksizes=0

    ! exchange with everybody else
    call MPI_AllGather(nsdx, 1, MP_Integer, &
         blocksizes, 1, MP_Integer, comm, status)
    _VERIFY(status)

    ! now everybody has blocksizes(nDEs)

    ntotal = sum(blocksizes) ! should be same as # of paired sources and sinks (npairs)
    _ASSERT(ntotal==numRoutings, 'Number source/sinks does not match')
    allocate (karray(numRoutings), _STAT) !declare as target!!!
    karray = 0
    allocate (displ(0:nDEs), _STAT) !declare as target!!!

    ksum = 0
    displ(0)=ksum
    do n=1,nDEs
       ksum = ksum + blocksizes(n)
       displ(n)=ksum
    end do
    ! as another sanity check: ksum should be the same as npairs
    _ASSERT(displ(nDEs)==ntotal, 'Displ source/sinks does not match')

    allocate(kseq(nsdx), _STAT)
    allocate(tmparray(numRoutings), _STAT)
    ! local k index
    k=0
    do i=1,size(LocalRoutings)
       if (mype==LocalRoutings(i)%srcPE) then
          k=k+1
          kseq(k) = LocalRoutings(i)%seqIdx
          kdx(k) = i
       end if
    end do

    call MPI_AllGatherV(kseq, nsdx, MP_Integer, &
         tmparray, blocksizes, displ, MP_Integer, comm, status)
    _VERIFY(STATUS)

    deallocate(kseq)
    do n=1,nDEs
       do k=1,blocksizes(n)
          i=tmparray(displ(n-1)+k)
          karray(i)=k
       end do
    end do
    deallocate(tmparray)

    allocate(RoutingType, _STAT)
    RoutingType%LocalRoutings => LocalRoutings
    RoutingType%karray => karray
    RoutingType%kdx => kdx
    RoutingType%BlockSizes => BlockSizes
    RoutingType%displ => displ

    return

  contains

    subroutine FindTile(Table, TileID, PE, Index)
      integer, intent(IN) :: Table(:,:), TileID
      integer, intent(OUT) :: PE, Index

      integer :: i

      do i=1,size(Table,2)
         if(TileId/=Table(1,i)) cycle
         PE = Table(3,i)
         Index = Table(2,i)
         exit
      end do

      return
    end subroutine FindTile

    subroutine Tile2Index(Routing, IdList)
      type(T_Routing), intent(INOUT) :: Routing
      integer,         intent(IN   ) :: IdList(:)

      integer :: i

      Routing%srcIndex = -1
      Routing%dstIndex = -1

      do i=1,size(IdList)
         if   (Routing%srcTileID==IdList(i)) then
            Routing%srcIndex = i
            if(Routing%dstIndex>0) exit
         elseif(Routing%dstTileID==IdList(i)) then
            Routing%dstIndex = i
            if(Routing%srcIndex>0) exit
         endif
      end do

      return
    end subroutine Tile2Index

  end subroutine InitializeRiverRouting


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!BOP

! !IROUTINE: RUN1 -- First stage Run method for the Surface component

! !INTERFACE:

  subroutine RUN1 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Interfaces to the children RUN1 methods, which compute
!   the surface exchange coefficients. In addition to exchange coefficients
!   for heat, moisture, and momentum, it also computes effective
!   surface values of the diffused quantities on the atmospheric grid.
!   These are exchange-coefficient-weighted averages of the tile values
!   within an atmospheric grid box.

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

    type (MAPL_MetaComp),  pointer  :: MAPL
    type (ESMF_GridComp),      pointer  :: GCS(:)
    type (ESMF_State),         pointer  :: GIM(:)
    type (ESMF_State),         pointer  :: GEX(:)
    type (ESMF_State)                   :: INTERNAL
    type (MAPL_LocStream)               :: LOCSTREAM
    character(len=ESMF_MAXSTR), pointer :: GCNames(:)
    integer                             :: NT
    integer, pointer, dimension(:)      :: TYPE

    type (T_SURFACE_STATE), pointer     :: SURF_INTERNAL_STATE
    type (SURF_wrap)                    :: WRAP

    type (ESMF_Time)                    :: CurrentTime

! Pointers to imports

    real, pointer, dimension(:,:) :: PS      => NULL()
    real, pointer, dimension(:,:) :: TA      => NULL()
    real, pointer, dimension(:,:) :: QA      => NULL()
    real, pointer, dimension(:,:) :: DZ      => NULL()
    real, pointer, dimension(:,:) :: UU      => NULL()
    real, pointer, dimension(:,:) :: UWINDLM => NULL()
    real, pointer, dimension(:,:) :: VWINDLM => NULL()
    real, pointer, dimension(:,:) :: PCU     => NULL()
    real, pointer, dimension(:,:) :: PHIS    => NULL()

    real, pointer, dimension(:,:) :: CHARNOCK => NULL()


! Pointers to gridded internals

    real, pointer, dimension(:,:) :: CT   => NULL()
    real, pointer, dimension(:,:) :: CQ   => NULL()
    real, pointer, dimension(:,:) :: CM   => NULL()
    real, pointer, dimension(:,:) :: CN   => NULL()
    real, pointer, dimension(:,:) :: TH   => NULL()
    real, pointer, dimension(:,:) :: QH   => NULL()
    real, pointer, dimension(:,:) :: SH   => NULL()
    real, pointer, dimension(:,:) :: UH   => NULL()
    real, pointer, dimension(:,:) :: VH   => NULL()
    real, pointer, dimension(:,:) :: TS   => NULL()
    real, pointer, dimension(:,:) :: QS   => NULL()
    real, pointer, dimension(:,:) :: RHOS => NULL()
    real, pointer, dimension(:,:) :: D0   => NULL()

! Pointers to gridded exports

    real, pointer, dimension(:,:) :: RI     => NULL()
    real, pointer, dimension(:,:) :: RE     => NULL()
    real, pointer, dimension(:,:) :: QDWL   => NULL()
    real, pointer, dimension(:,:) :: QFRL   => NULL()
    real, pointer, dimension(:,:) :: USTAR  => NULL()
    real, pointer, dimension(:,:) :: BSTAR  => NULL()
    real, pointer, dimension(:,:) :: LAI    => NULL()
    real, pointer, dimension(:,:) :: GRN    => NULL()
    real, pointer, dimension(:,:) :: ROOTL  => NULL()
    real, pointer, dimension(:,:) :: Z2CH   => NULL()
    real, pointer, dimension(:,:) :: VNT    => NULL()
    real, pointer, dimension(:,:) :: GST    => NULL()
    real, pointer, dimension(:,:) :: Z0     => NULL()
    real, pointer, dimension(:,:) :: MOT2M  => NULL()
    real, pointer, dimension(:,:) :: MOQ2M  => NULL()
    real, pointer, dimension(:,:) :: MOU2M  => NULL()
    real, pointer, dimension(:,:) :: MOV2M  => NULL()
    real, pointer, dimension(:,:) :: MOT10M => NULL()
    real, pointer, dimension(:,:) :: MOQ10M => NULL()
    real, pointer, dimension(:,:) :: MOU10M => NULL()
    real, pointer, dimension(:,:) :: MOV10M => NULL()
    real, pointer, dimension(:,:) :: MOU50M => NULL()
    real, pointer, dimension(:,:) :: MOV50M => NULL()
    real, pointer, dimension(:,:) :: ITY    => NULL()
    real, pointer, dimension(:,:) :: NITY   => NULL()
    real, pointer, dimension(:,:) :: Z0H    => NULL()

    real, pointer  :: T2MDEW   (:,:)      => NULL()
    real, pointer  :: T2MWET   (:,:)      => NULL()
    real, pointer  :: RH2M     (:,:)      => NULL()
    real, pointer  :: UU10M    (:,:)      => NULL()

! Pointers to tile versions of imports

    real, pointer, dimension(:) :: PSTILE      => NULL()
    real, pointer, dimension(:) :: TATILE      => NULL()
    real, pointer, dimension(:) :: QATILE      => NULL()
    real, pointer, dimension(:) :: DZTILE      => NULL()
    real, pointer, dimension(:) :: UUTILE      => NULL()
    real, pointer, dimension(:) :: UWINDLMTILE => NULL()
    real, pointer, dimension(:) :: VWINDLMTILE => NULL()
    real, pointer, dimension(:) :: PCUTILE     => NULL()

    real, pointer, dimension(:) :: CHARNOCKTILE=> NULL()

! Pointers to tiled versions of internals

    real, pointer, dimension(:) :: CTTILE => NULL()
    real, pointer, dimension(:) :: CMTILE => NULL()
    real, pointer, dimension(:) :: CQTILE => NULL()
    real, pointer, dimension(:) :: CNTILE => NULL()
    real, pointer, dimension(:) :: RETILE => NULL()
    real, pointer, dimension(:) :: RITILE => NULL()
    real, pointer, dimension(:) :: THTILE => NULL()
    real, pointer, dimension(:) :: QHTILE => NULL()
    real, pointer, dimension(:) :: UHTILE => NULL()
    real, pointer, dimension(:) :: VHTILE => NULL()
    real, pointer, dimension(:) :: TSTILE => NULL()
    real, pointer, dimension(:) :: QSTILE => NULL()
    real, pointer, dimension(:) :: D0TILE => NULL()

! Pointers to tiled versions of exports

    real, pointer, dimension(:) :: LAITILE     => NULL()
    real, pointer, dimension(:) :: GRNTILE     => NULL()
    real, pointer, dimension(:) :: ROOTLTILE   => NULL()
    real, pointer, dimension(:) :: Z2CHTILE    => NULL()
    real, pointer, dimension(:) :: VNTTILE     => NULL()
    real, pointer, dimension(:) :: GSTTILE     => NULL()
    real, pointer, dimension(:) :: Z0HTILE     => NULL()
    real, pointer, dimension(:) :: Z0TILE      => NULL()
    real, pointer, dimension(:) :: MOT2MTILE   => NULL()
    real, pointer, dimension(:) :: MOQ2MTILE   => NULL()
    real, pointer, dimension(:) :: MOU2MTILE   => NULL()
    real, pointer, dimension(:) :: MOV2MTILE   => NULL()
    real, pointer, dimension(:) :: MOT10MTILE  => NULL()
    real, pointer, dimension(:) :: MOQ10MTILE  => NULL()
    real, pointer, dimension(:) :: MOU10MTILE  => NULL()
    real, pointer, dimension(:) :: MOV10MTILE  => NULL()
    real, pointer, dimension(:) :: MOU50MTILE  => NULL()
    real, pointer, dimension(:) :: MOV50MTILE  => NULL()
    real, pointer, dimension(:) :: ITYTILE     => NULL()

! Pointers to other stuff

    real, pointer, dimension(:,:) :: RH2MO => NULL()
    real, pointer, dimension(:,:) :: ALHX  => NULL()

    type (MAPL_MetaComp), pointer :: CHILD_MAPL
    integer                           :: I

    integer    :: iUseInterp
    logical    :: UseInterp

    integer    :: IM, JM, YEAR, MONTH, DAY, HR, SE, MN

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = 'Run1'
    call ESMF_GridCompGet( GC, name=COMP_NAME, VM=VMG, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_UserCompGetInternalState(GC, 'SURF_state', wrap, status)
    VERIFY_(STATUS)

    SURF_INTERNAL_STATE => WRAP%PTR

! Start Total timer
!------------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"-RUN1" )

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,             &
         LOCSTREAM = LOCSTREAM,                  &
         GIM       = GIM,                        &
         GEX       = GEX,                        &
         TILETYPES = TYPE,                       &
         GCS       = GCS,                        &
         GCNAMES   = GCNAMES,                    &
         INTERNAL_ESMF_STATE = INTERNAL,         &
                                       RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_ClockGet(CLOCK, currTime=CurrentTime, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_TimeGet (currentTime,               &
                       YY=YEAR, MM=MONTH, DD=DAY, &
                       H=HR,    M=MN,     S=SE,   &
                                        RC=STATUS )

! Pointers to imports
!--------------------

    call MAPL_GetPointer(IMPORT  , PS    , 'PS'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , DZ    , 'DZ'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , UU    , 'SPEED' ,  RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(IMPORT  , UWINDLM , 'UA'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , VWINDLM , 'VA'  ,  RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(IMPORT  , TA    , 'TA'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , QA    , 'QA'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , PCU   , 'PCU'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , PHIS  , 'PHIS'  ,  RC=STATUS); VERIFY_(STATUS)

    if (DO_WAVES /= 0) then
      call MAPL_GetPointer(IMPORT  , CHARNOCK , 'CHARNOCK',  RC=STATUS); VERIFY_(STATUS)
    end if

! Pointers to grid outputs
!-------------------------

! These are computed by the children in tile space and transformed

    call MAPL_GetPointer(EXPORT  , RI    , 'RI'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RE    , 'RE'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , GRN   , 'GRN'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ROOTL , 'ROOTL' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , Z2CH  , 'Z2CH'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , VNT   , 'VENT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , GST   , 'GUST'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , Z0    , 'Z0'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , Z0H   , 'Z0H'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , NITY  , 'NITY'  ,  RC=STATUS); VERIFY_(STATUS)

! if we are running MO sfc layer, get these exports from the MO values
  if (CHOOSEMOSFC.eq.1) then
!!AMM alloc=true for mot2m and moq2m because a bunch of other diagnostics rely on them
    call MAPL_GetPointer(EXPORT  , MOT2M , 'T2M' ,  ALLOC = .true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MOQ2M , 'Q2M' ,  ALLOC = .true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MOU2M , 'U2M' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MOV2M , 'V2M' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MOT10M , 'T10M' ,ALLOC = .true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MOQ10M , 'Q10M' ,ALLOC = .true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MOU10M , 'U10M' ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MOV10M , 'V10M' ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MOU50M , 'U50M' ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MOV50M , 'V50M' ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , UU10M   , 'UU10M'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RH2M    , 'RH2M'   ,  ALLOC = .true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , T2MDEW , 'T2MDEW' ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , T2MWET , 'T2MWET' ,   RC=STATUS); VERIFY_(STATUS)
  else
    NULLIFY(MOU2M)
    NULLIFY(MOV2M)
    NULLIFY(MOT2M)
    NULLIFY(MOQ2M)
    NULLIFY(MOU10M)
    NULLIFY(MOV10M)
    NULLIFY(MOT10M)
    NULLIFY(MOQ10M)
    NULLIFY(MOU50M)
    NULLIFY(MOV50M)
  endif

! Need to force LAI if GRN is required.

    call MAPL_GetPointer(EXPORT  , LAI   , 'LAI'   ,  alloc=associated(GRN), RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ITY   , 'ITY'   ,  alloc=associated(NITY), RC=STATUS)
    VERIFY_(STATUS)

! These are computed by SURFACE in grid space and have no tile versions

    call MAPL_GetPointer(EXPORT  , SH    , 'SHAT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , QDWL  , 'QDWL'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , QFRL  , 'QFRL'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , BSTAR , 'BSTAR' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , USTAR , 'USTAR' ,  RC=STATUS); VERIFY_(STATUS)

! These are force-allocated because run2 needs them or their space

    call MAPL_GetPointer(INTERNAL, TS    , 'TS'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QS    , 'QS'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CT    , 'CT'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CQ    , 'CQ'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CM    , 'CM'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CN    , 'CN'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QH    , 'QHAT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, TH    , 'THAT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, UH    , 'UHAT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, VH    , 'VHAT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, RHOS  , 'RHOS'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, D0    , 'D0'    ,  RC=STATUS); VERIFY_(STATUS)

! Size of exchange grid
!----------------------

    NT = size(TYPE)

!  Allocate tile versions of imports
!-----------------------------------

    allocate(   PSTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   DZTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   UUTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   UWINDLMTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   VWINDLMTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   TATILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   QATILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  PCUTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    if (DO_WAVES /= 0) then
      allocate(CHARNOCKTILE(NT), STAT=STATUS)
      VERIFY_(STATUS)
    end if

! Imports at the tiles
!---------------------
     call MAPL_GetResource(MAPL, iUseInterp, 'INTERPOLATE_ATMLM:', &
         default=0, RC=STATUS )
    VERIFY_(STATUS)
    useInterp = (iUseInterp /= 0)

    call MAPL_LocStreamTransform( LOCSTREAM,  PSTILE,  PS, INTERP=useInterp, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM,  TATILE,  TA, INTERP=useInterp, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM,  QATILE,  QA, INTERP=useInterp, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM,  DZTILE,  DZ, INTERP=useInterp, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM,  UUTILE,  UU, INTERP=useInterp, RC=STATUS); VERIFY_(STATUS)

    call MAPL_LocStreamTransform( LOCSTREAM,  UWINDLMTILE,  UWINDLM, INTERP=useInterp, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM,  VWINDLMTILE,  VWINDLM, INTERP=useInterp, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM,  PCUTILE,      PCU,     RC=STATUS); VERIFY_(STATUS)

    if (DO_WAVES /= 0) then
      call MAPL_LocStreamTransform( LOCSTREAM,  CHARNOCKTILE,  CHARNOCK, RC=STATUS); VERIFY_(STATUS)
    end if

! Allocate tile versions of internal
!------------------------------------

!  We do not need a tile version of RHOS

    allocate(   CTTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   CQTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   CMTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   CNTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   TSTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   QSTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   THTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   QHTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   UHTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   VHTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   D0TILE(NT), STAT=STATUS)
    VERIFY_(STATUS)

! Allocate tile versions of needed exports that are filled by children
!---------------------------------------------------------------------

    call MKTILE(LAI   , LAITILE     , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(GRN   , GRNTILE     , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ROOTL , ROOTLTILE   , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(Z2CH  , Z2CHTILE    , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(VNT   , VNTTILE     , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(GST   , GSTTILE     , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(Z0H   , Z0HTILE     , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(Z0    , Z0TILE      , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MOT2M , MOT2MTILE  , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MOQ2M , MOQ2MTILE  , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MOU2M , MOU2MTILE , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MOV2M , MOV2MTILE , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MOT10M , MOT10MTILE  , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MOQ10M , MOQ10MTILE  , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MOU10M , MOU10MTILE , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MOV10M , MOV10MTILE , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MOU50M , MOU50MTILE , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MOV50M , MOV50MTILE , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ITY    , ITYTILE    , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RI    , RITILE      , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RE    , RETILE      , NT, RC=STATUS); VERIFY_(STATUS)

! If the child does not produce them, we want these zeroed.
!---------------------------------------------------------

    UHTILE = 0.0
    VHTILE = 0.0
    D0TILE = 0.0

! Do the run1 (surface layer calculations) for each child.
!--------------------------------------------------------

    do I = 1, NUM_CHILDREN
       call DOCDS(I, NT, RC=STATUS)
       VERIFY_(STATUS)
    end do

! Grid exports
!-------------

    if(associated(GRNTILE)) then
       where(GRNTILE /= MAPL_UNDEF) GRNTILE =  GRNTILE*LAITILE
    endif

    if(associated(    RI)) then
       call MAPL_LocStreamTransform( LOCSTREAM,     RI,     RITILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(    RE)) then
       call MAPL_LocStreamTransform( LOCSTREAM,     RE,     RETILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(   LAI)) then
       call MAPL_LocStreamTransform( LOCSTREAM,    LAI,    LAITILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(   GRN)) then
       call MAPL_LocStreamTransform( LOCSTREAM,    GRN,    GRNTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( ROOTL)) then
       call MAPL_LocStreamTransform( LOCSTREAM,  ROOTL,  ROOTLTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  Z2CH)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   Z2CH,   Z2CHTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(   VNT)) then
       call MAPL_LocStreamTransform( LOCSTREAM,    VNT,    VNTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(   GST)) then
       call MAPL_LocStreamTransform( LOCSTREAM,    GST,    GSTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(   Z0H)) then
       call MAPL_LocStreamTransform( LOCSTREAM,    Z0H,    Z0HTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(    Z0)) then
       call MAPL_LocStreamTransform( LOCSTREAM,     Z0,     Z0TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(ITY)) then
       call MAPL_LocStreamTransform( LOCSTREAM,    ITY,    ITYTILE, SAMPLE=.true., RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(MOU50M)) then
       call MAPL_LocStreamTransform( LOCSTREAM, MOU50M, MOU50MTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(MOV50M)) then
       call MAPL_LocStreamTransform( LOCSTREAM, MOV50M, MOV50MTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(MOT10M)) then
       call MAPL_LocStreamTransform( LOCSTREAM,  MOT10M,  MOT10MTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(MOQ10M)) then
       call MAPL_LocStreamTransform( LOCSTREAM,  MOQ10M,  MOQ10MTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(MOU10M)) then
       call MAPL_LocStreamTransform( LOCSTREAM, MOU10M, MOU10MTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(MOV10M)) then
       call MAPL_LocStreamTransform( LOCSTREAM, MOV10M, MOV10MTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(MOT2M)) then
       call MAPL_LocStreamTransform( LOCSTREAM,  MOT2M,  MOT2MTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(MOQ2M)) then
       call MAPL_LocStreamTransform( LOCSTREAM,  MOQ2M,  MOQ2MTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(MOU2M)) then
       call MAPL_LocStreamTransform( LOCSTREAM, MOU2M, MOU2MTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(MOV2M)) then
       call MAPL_LocStreamTransform( LOCSTREAM, MOV2M, MOV2MTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(MOU10M) .and. associated(MOV10M)) then
       if(associated(UU10M)) UU10M = sqrt(MOU10M**2 + MOV10M**2)
    end if

    if(associated(GRN)) then
       where(GRN /= MAPL_UNDEF .and. LAI> 0.0)
          GRN =  GRN/LAI
       elsewhere
          GRN =  MAPL_UNDEF
       end where
    endif

      if( associated(NITY) ) then
         NITY = ITY
         where ( ITY==1 ) NITY=1
         where ( ITY==2 ) NITY=2
         where ( ITY==3 ) NITY=4
         where ( ITY==4 ) NITY=7
         where ( ITY==5 ) NITY=9
         where ( ITY==6 ) NITY=10
         where ( ITY==7 ) NITY=11
         where ( ITY==13) NITY=13
      endif


!    1  ...  broadleave-evergreen trees (tropical forest)  MAP TO ITYP=1
!    2  ...  broadleave-deciduous trees   MAP to ITYP=2
!    3  ...  broadleave and needle leave trees (mixed forest)
!                   (For this, we map 1/2 to ITYP=2 and 1/2 to ITYP=3)
!    4  ...  needleleave-evergreen trees   MAP to ITYP=3
!    5  ...  needleleave-deciduous trees (larch)   MAP to ITYP=3
!    6  ...  broadleave trees with groundcover (savanna)
!                   (For this, we map 1/10 to ITYP=2 and 9/10 to ITYP=4)
!    7  ...  groundcover only (perenial)    MAP to ITYP=4
!    8  ...  broadleave shrubs with perenial groundcover
!                   (For this, we map 0.25(?) to ITYP=5 and 0.75 to ITYP=4)
!    9  ...  broadleave shrubs with bare soil  MAP to ITYP=5
!   10  ...  dwarf trees and shrubs with ground cover (trunda) MAP to ITYP=6
!   11  ...  bare soil  MAP to ITYP=7
!   12  ...  cultivations (use parameters from type 7)  MAP to ITYP=4
!   13  ...  glacial



! Effective surface values on atmos grid. These and the ceoffs
!   are forced exports because run2 needs them.
!-------------------------------------------------------------

    THTILE = THTILE*CTTILE
    QHTILE = QHTILE*CQTILE
    UHTILE = UHTILE*CMTILE
    VHTILE = VHTILE*CMTILE

    call MAPL_LocStreamTransform( LOCSTREAM, CT, CTTILE, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, CM, CMTILE, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, CQ, CQTILE, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, CN, CNTILE, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, TH, THTILE, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, QH, QHTILE, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, UH, UHTILE, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, VH, VHTILE, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, D0, D0TILE, RC=STATUS)
    VERIFY_(STATUS)

! These are in the internal state

    QH    = QH/CQ
    TH    = TH/CT
    UH    = UH/CM
    VH    = VH/CM
    RHOS  = PS / ( MAPL_RGAS*TA*(1.+MAPL_VIREPS*QA) )

    if(associated(QDWL )) QDWL  = 0.
    if(associated(QFRL )) QFRL  = 0.
    if(associated(SH   )) SH    = MAPL_CP*TH + PHIS
    if(associated(USTAR)) USTAR = sqrt(CM*UU/RHOS)
    if(associated(BSTAR)) BSTAR = (MAPL_GRAV/(RHOS*sqrt(CM*max(UU,1.e-30)/RHOS))) *  &
        (CT*(TH-TA-(MAPL_GRAV/MAPL_CP)*DZ)/TA + MAPL_VIREPS*CQ*(QH-QA))

! if we are running MO sfc layer, get these exports from the MO values
  if (CHOOSEMOSFC.eq.1) then

     IM = size(PS,1)
     JM = size(PS,2)
     allocate(RH2MO(IM,JM),STAT=STATUS)
     VERIFY_(STATUS)

     if(associated(MOT2M).and.associated(MOQ2M)) then
        RH2MO = min(MOQ2M/GEOS_QSAT(MOT2M,PS,PASCALS=.true.),1.0)*100.
        if(associated(RH2M)) RH2M = RH2MO
     endif

     deallocate(RH2MO)

     if(associated(T2MDEW)) then
      T2MDEW = MOT2M
      do i = 1,4
       T2MDEW = T2MDEW + (MOQ2M-GEOS_QSAT(T2MDEW,PS,PASCALS=.true.))/GEOS_DQSAT(T2MDEW,PS,PASCALS=.true.)
      enddo
     endif

     if(associated(T2MWET)) then
      allocate(ALHX(IM,JM),STAT=STATUS)
      VERIFY_(STATUS)
      T2MWET = MOT2M
      do i = 1,10
       ALHX = (MAX(MIN(1.-(T2MWET-233.16)/40.,1.),0.))**4
       ALHX = (1.0-ALHX)*MAPL_ALHL + ALHX*MAPL_ALHS
       T2MWET = T2MWET + (((ALHX/MAPL_CP)/(1.+(ALHX/MAPL_CP)*GEOS_DQSAT(T2MWET,PS,PASCALS=.true.)))* &
                                  (MOQ2M-GEOS_QSAT(T2MWET,PS,PASCALS=.true.)))
      enddo
      deallocate(ALHX)
     endif

  endif          ! end of MO sfc layer if sequence

! Clean-up
!---------

    if(associated(  LAITILE)) deallocate(   LAITILE)
    if(associated(  GRNTILE)) deallocate(   GRNTILE)
    if(associated(ROOTLTILE)) deallocate( ROOTLTILE)
    if(associated( Z2CHTILE)) deallocate(  Z2CHTILE)
    if(associated(  VNTTILE)) deallocate(   VNTTILE)
    if(associated(  GSTTILE)) deallocate(   GSTTILE)
    if(associated(  Z0HTILE)) deallocate(   Z0HTILE)
    if(associated(   Z0TILE)) deallocate(    Z0TILE)
    if(associated(   ITYTILE))deallocate(   ITYTILE)
    if(associated(MOU50MTILE))deallocate(MOU50MTILE)
    if(associated(MOV50MTILE))deallocate(MOV50MTILE)
    if(associated(MOT2MTILE)) deallocate( MOT2MTILE)
    if(associated(MOQ2MTILE)) deallocate( MOQ2MTILE)
    if(associated(MOU2MTILE)) deallocate( MOU2MTILE)
    if(associated(MOV2MTILE)) deallocate( MOV2MTILE)
    if(associated(MOT10MTILE))deallocate(MOT10MTILE)
    if(associated(MOQ10MTILE))deallocate(MOQ10MTILE)
    if(associated(MOU10MTILE))deallocate(MOU10MTILE)
    if(associated(MOV10MTILE))deallocate(MOV10MTILE)
    if(associated(   RITILE)) deallocate(    RITILE)
    if(associated(   RETILE)) deallocate(    RETILE)
    if(associated(   CTTILE)) deallocate(    CTTILE)
    if(associated(   CQTILE)) deallocate(    CQTILE)
    if(associated(   CMTILE)) deallocate(    CMTILE)
    if(associated(   CNTILE)) deallocate(    CNTILE)
    if(associated(   TSTILE)) deallocate(    TSTILE)
    if(associated(   QSTILE)) deallocate(    QSTILE)
    if(associated(   THTILE)) deallocate(    THTILE)
    if(associated(   QHTILE)) deallocate(    QHTILE)
    if(associated(   UHTILE)) deallocate(    UHTILE)
    if(associated(   VHTILE)) deallocate(    VHTILE)
    if(associated(   D0TILE)) deallocate(    D0TILE)
    if(associated(   DZTILE)) deallocate(    DZTILE)
    if(associated(   PSTILE)) deallocate(    PSTILE)
    if(associated(  PCUTILE)) deallocate(   PCUTILE)
    if(associated(   QATILE)) deallocate(    QATILE)
    if(associated(   TATILE)) deallocate(    TATILE)
    if(associated(   UUTILE)) deallocate(    UUTILE)
    if(associated(   UWINDLMTILE)) deallocate(    UWINDLMTILE)
    if(associated(   VWINDLMTILE)) deallocate(    VWINDLMTILE)

    if(associated(CHARNOCKTILE))   deallocate(CHARNOCKTILE)

!  All done
!-----------

    call MAPL_TimerOff(MAPL,"-RUN1" )
    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)

    contains

      subroutine DOCDS(type, NT, RC)
        integer,           intent( IN) :: type
        integer,           intent( IN) :: NT
        integer, optional, intent(OUT) :: RC

!  Locals

        character(len=ESMF_MAXSTR)   :: IAm
        integer                      :: STATUS
        integer                      :: N

        type (MAPL_LocStreamXFORM)   :: XFORM
        real, pointer                :: DUM(:) => NULL()

!  Begin...
!----------

        IAM = trim(COMP_NAME) //  "DOCDS"

        call MAPL_TimerOn(MAPL,           trim(GCNames(type)))
        call MAPL_TimerOn(MAPL,"--RUN1_"//trim(GCNames(type)))

        call MAPL_Get(MAPL, GCNAMES = GCNAMES, RC=STATUS )
        VERIFY_(STATUS)

! Fill the child's locstream imports from the Surface exchange grid imports
!--------------------------------------------------------------------------

        XFORM = surf_internal_state%xform_in(type)

        call FILLIN_TILE(GIM(type),  'PS',  PSTILE, XFORM, RC=STATUS)
        VERIFY_(STATUS)
        call FILLIN_TILE(GIM(type),  'DZ',  DZTILE, XFORM, RC=STATUS)
        VERIFY_(STATUS)
        call FILLIN_TILE(GIM(type),  'UU',  UUTILE, XFORM, RC=STATUS)
        VERIFY_(STATUS)
        call FILLIN_TILE(GIM(type),  'UWINDLMTILE',  UWINDLMTILE, XFORM, RC=STATUS)
        VERIFY_(STATUS)
        call FILLIN_TILE(GIM(type),  'VWINDLMTILE',  VWINDLMTILE, XFORM, RC=STATUS)
        VERIFY_(STATUS)
        call FILLIN_TILE(GIM(type),  'TA',  TATILE, XFORM, RC=STATUS)
        VERIFY_(STATUS)
        call FILLIN_TILE(GIM(type),  'QA',  QATILE, XFORM, RC=STATUS)
        VERIFY_(STATUS)
        call FILLIN_TILE(GIM(type), 'PCU', PCUTILE, XFORM, RC=STATUS)
        VERIFY_(STATUS)
        if (DO_WAVES /= 0) then
          call FILLIN_TILE(GIM(type), 'CHARNOCK', CHARNOCKTILE, XFORM, RC=STATUS)
          VERIFY_(STATUS)
        end if


! Allocate the child's needed exports
!------------------------------------

!  Note that the first batch is really forced by the allocation in RUN1 proper.

        call MAPL_GetPointer(GEX(type), dum, 'TST', ALLOC=associated( TSTILE), RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum, 'QST', ALLOC=associated( QSTILE), RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,  'TH', ALLOC=associated( THTILE), RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,  'QH', ALLOC=associated( QHTILE), RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum, 'CHT', ALLOC=associated( CTTILE), RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum, 'CQT', ALLOC=associated( CQTILE), RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum, 'CMT', ALLOC=associated( CMTILE), RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum, 'CNT', ALLOC=associated( CNTILE), RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum, 'RIT', ALLOC=associated( RITILE), RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,  'Z0', ALLOC=associated( Z0TILE), RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,'MOU50M',ALLOC=associated(MOU50MTILE),RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,'MOV50M',ALLOC=associated(MOV50MTILE),RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,'MOT10M',ALLOC=associated(MOT10MTILE),RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,'MOQ10M',ALLOC=associated(MOQ10MTILE),RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,'MOU10M',ALLOC=associated(MOU10MTILE),RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,'MOV10M',ALLOC=associated(MOV10MTILE),RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,'MOT2M',ALLOC=associated(MOT2MTILE),RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,'MOQ2M',ALLOC=associated(MOQ2MTILE),RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,'MOU2M',ALLOC=associated(MOU2MTILE),RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,'MOV2M',ALLOC=associated(MOV2MTILE),RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum, 'Z0H', ALLOC=associated(Z0HTILE), RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,'VENT', ALLOC=associated(VNTTILE), RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,'GUST', ALLOC=associated(GSTTILE), RC=STATUS)
        VERIFY_(STATUS)

! These cannot be verified, because they dont exists in all children.
!-------------------------------------------------------------------

        call MAPL_GetPointer(GEX(type), dum,   'LAI', ALLOC=associated(  LAITILE), notFoundOK=.true., RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,   'GRN', ALLOC=associated(  GRNTILE), notFoundOK=.true., RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum, 'ROOTL', ALLOC=associated(ROOTLTILE), notFoundOK=.true., RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,  'Z2CH', ALLOC=associated( Z2CHTILE), notFoundOK=.true., RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,    'D0', ALLOC=associated(   D0TILE), notFoundOK=.true., RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,    'UH', ALLOC=associated(   UHTILE), notFoundOK=.true., RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,    'VH', ALLOC=associated(   VHTILE), notFoundOK=.true., RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,   'RET', ALLOC=associated(   RETILE), notFoundOK=.true., RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(GEX(type), dum,   'ITY', ALLOC=associated(  ITYTILE), notFoundOK=.true., RC=STATUS)
        VERIFY_(STATUS)

! Call Child
!-----------

        call ESMF_GridCompRun(GCS(type), &
             importState=GIM(type), exportState=GEX(type), &
             clock=CLOCK, PHASE=1, userRC=STATUS )
        VERIFY_(STATUS)

! Use childs exports to fill exchange grid exports.
!--------------------------------------------------

        XFORM = surf_internal_state%xform_out(type)

! Again the first batch is forced, but we test anyway since it does not hurt.

        if(associated(TSTILE)) then
           call FILLOUT_TILE(GEX(type),   'TST',   TSTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        endif
        if(associated(QSTILE)) then
           call FILLOUT_TILE(GEX(type),   'QST',   QSTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        endif
        if(associated(THTILE)) then
           call FILLOUT_TILE(GEX(type),    'TH',   THTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        endif
        if(associated(QHTILE)) then
           call FILLOUT_TILE(GEX(type),    'QH',   QHTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        endif
        if(associated(UHTILE)) then
           call FILLOUT_TILE(GEX(type),    'UH',   UHTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        endif
        if(associated(VHTILE)) then
           call FILLOUT_TILE(GEX(type),    'VH',   VHTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        endif
        if(associated(CTTILE)) then
           call FILLOUT_TILE(GEX(type),   'CHT',   CTTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        endif
        if(associated(CQTILE)) then
           call FILLOUT_TILE(GEX(type),   'CQT',   CQTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        endif
        if(associated(CMTILE)) then
           call FILLOUT_TILE(GEX(type),   'CMT',   CMTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        endif
        if(associated(CNTILE)) then
           call FILLOUT_TILE(GEX(type),   'CNT',   CNTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        endif
        if(associated(D0TILE)) then
           call FILLOUT_TILE(GEX(type),    'D0',    D0TILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if

       if(associated(RITILE)) then
           call FILLOUT_TILE(GEX(type),   'RIT',   RITILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        endif
       if(associated(RETILE)) then
           call FILLOUT_TILE(GEX(type),   'RET',   RETILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        endif
        if(associated(LAITILE)) then
           call FILLOUT_TILE(GEX(type),   'LAI',  LAITILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(GRNTILE)) then
           call FILLOUT_TILE(GEX(type),   'GRN',  GRNTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(ROOTLTILE)) then
           call FILLOUT_TILE(GEX(type), 'ROOTL', ROOTLTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(Z2CHTILE)) then
           call FILLOUT_TILE(GEX(type),  'Z2CH',  Z2CHTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(VNTTILE)) then
           call FILLOUT_TILE(GEX(type),  'VENT',   VNTTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(GSTTILE)) then
           call FILLOUT_TILE(GEX(type),  'GUST',   GSTTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(Z0HTILE)) then
           call FILLOUT_TILE(GEX(type),   'Z0H',   Z0HTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(Z0TILE)) then
           call FILLOUT_TILE(GEX(type),    'Z0',    Z0TILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(ITYTILE)) then
           call FILLOUT_TILE(GEX(type),   'ITY',   ITYTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(MOU50MTILE)) then
           call FILLOUT_TILE(GEX(type),'MOU50M',MOU50MTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(MOV50MTILE)) then
           call FILLOUT_TILE(GEX(type),'MOV50M',MOV50MTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(MOT10MTILE)) then
           call FILLOUT_TILE(GEX(type), 'MOT10M', MOT10MTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(MOQ10MTILE)) then
           call FILLOUT_TILE(GEX(type), 'MOQ10M', MOQ10MTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(MOU10MTILE)) then
           call FILLOUT_TILE(GEX(type),'MOU10M',MOU10MTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(MOV10MTILE)) then
           call FILLOUT_TILE(GEX(type),'MOV10M',MOV10MTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(MOT2MTILE)) then
           call FILLOUT_TILE(GEX(type), 'MOT2M', MOT2MTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(MOQ2MTILE)) then
           call FILLOUT_TILE(GEX(type), 'MOQ2M', MOQ2MTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(MOU2MTILE)) then
           call FILLOUT_TILE(GEX(type),'MOU2M',MOU2MTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if
        if(associated(MOV2MTILE)) then
           call FILLOUT_TILE(GEX(type),'MOV2M',MOV2MTILE, XFORM, RC=STATUS)
           VERIFY_(STATUS)
        end if

        call MAPL_TimerOff(MAPL,"--RUN1_"//trim(GCNames(type)))
        call MAPL_TimerOff(MAPL,           trim(GCNames(type)))
        RETURN_(ESMF_SUCCESS)

      end subroutine DOCDS

    end subroutine RUN1

!BOP

! !IROUTINE: RUN2 -- Second Run method for the Surface component

! !INTERFACE:

  subroutine RUN2 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION:

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),     pointer   :: MAPL
    type (MAPL_SunOrbit)                :: ORBIT
    type (ESMF_State)                   :: INTERNAL
    type (ESMF_Alarm)                   :: ALARM
    type (ESMF_GridComp),      pointer  :: GCS(:)
    type (ESMF_State),         pointer  :: GIM(:)
    type (ESMF_State),         pointer  :: GEX(:)
    type (ESMF_TimeInterval)            :: DELT
    type (MAPL_LocStream)               :: LOCSTREAM

    type (ESMF_FieldBundle)             :: Bundle
    type (ESMF_Field)                   :: Field
    type (ESMF_Grid)                    :: GRID
    type (ESMF_Time)                    :: CurrentTime
    character(len=ESMF_MAXSTR)          :: PRECIP_FILE

    type (T_SURFACE_STATE), pointer     :: surf_internal_state
    type (SURF_wrap)                    :: wrap
    type(ESMF_VM)                       :: VM
    character(len=ESMF_MAXSTR), pointer :: GCNames(:)
    character(len=ESMF_MAXSTR)          :: DischargeAdjustFile
    integer                             :: IM, JM,  NT, comm, myPE, nDEs
    real                                :: SC, MG, SB
    logical                             :: USE_NRLSSI2
    real,    pointer, dimension(:,:)    :: LATS     => NULL()
    real,    pointer, dimension(:,:)    :: LONS     => NULL()
    real,    pointer, dimension(:)      :: tileLATS => NULL()
    real,    pointer, dimension(:)      :: tileLONS => NULL()
    integer, pointer, dimension(:)      :: tiletype => NULL()
    real,    pointer, dimension(:)      :: dataptr  => NULL()
    real,    pointer, dimension(:)      :: DISCHARGE_ADJUST  => NULL()
    integer                             :: K
    character(len=ESMF_MAXSTR), pointer :: AERO_DP_FIELD_NAME(:)
    integer                             :: NUM_AERO_DP
    character(len=ESMF_MAXSTR)          :: FIELD_NAME
    integer                             :: N_DUDP, N_DUSV, N_DUWT, N_DUSD
    integer                             :: N_BCDP, N_BCSV, N_BCWT, N_BCSD
    integer                             :: N_OCDP, N_OCSV, N_OCWT, N_OCSD
    integer                             :: N_SUDP, N_SUSV, N_SUWT, N_SUSD
    integer                             :: N_SSDP, N_SSSV, N_SSWT, N_SSSD
    integer                             :: N_DUDP1, N_DUSV1, N_DUWT1, N_DUSD1
    integer                             :: N_BCDP1, N_BCSV1, N_BCWT1, N_BCSD1
    integer                             :: N_OCDP1, N_OCSV1, N_OCWT1, N_OCSD1
    integer                             :: N_SUDP1, N_SUSV1, N_SUWT1, N_SUSD1
    integer                             :: N_SSDP1, N_SSSV1, N_SSWT1, N_SSSD1

! Pointers to imports

    real, pointer, dimension(:,:)   :: PS        => NULL()
    real, pointer, dimension(:,:)   :: DZ        => NULL()
    real, pointer, dimension(:,:)   :: UU        => NULL()
    real, pointer, dimension(:,:)   :: EVAP      => NULL()
    real, pointer, dimension(:,:)   :: SH        => NULL()
    real, pointer, dimension(:,:)   :: DEVAP     => NULL()
    real, pointer, dimension(:,:)   :: DSH       => NULL()
    real, pointer, dimension(:,:)   :: PCU       => NULL()
    real, pointer, dimension(:,:)   :: PLS       => NULL()
    real, pointer, dimension(:,:)   :: SNOFL     => NULL()
    real, pointer, dimension(:,:)   :: ICEFL     => NULL()
    real, pointer, dimension(:,:)   :: FRZRFL    => NULL()
    real, pointer, dimension(:,:)   :: TAUX      => NULL()
    real, pointer, dimension(:,:)   :: TAUY      => NULL()
    real, pointer, dimension(:,:)   :: DRPARN    => NULL()
    real, pointer, dimension(:,:)   :: DFPARN    => NULL()
    real, pointer, dimension(:,:)   :: DRNIRN    => NULL()
    real, pointer, dimension(:,:)   :: DFNIRN    => NULL()
    real, pointer, dimension(:,:)   :: DRUVRN    => NULL()
    real, pointer, dimension(:,:)   :: DFUVRN    => NULL()
    real, pointer, dimension(:,:)   :: LWDNSRF   => NULL()
    real, pointer, dimension(:,:)   :: ALW       => NULL()
    real, pointer, dimension(:,:)   :: BLW       => NULL()
    real, pointer, dimension(:,:)   :: AERO_DP_I => NULL()
    real, pointer, dimension(:,:,:) :: AERO_DP   => NULL()
    real, pointer, dimension(:,:,:) :: DUDP      => NULL()
    real, pointer, dimension(:,:,:) :: DUSV      => NULL()
    real, pointer, dimension(:,:,:) :: DUWT      => NULL()
    real, pointer, dimension(:,:,:) :: DUSD      => NULL()
    real, pointer, dimension(:,:,:) :: BCDP      => NULL()
    real, pointer, dimension(:,:,:) :: BCSV      => NULL()
    real, pointer, dimension(:,:,:) :: BCWT      => NULL()
    real, pointer, dimension(:,:,:) :: BCSD      => NULL()
    real, pointer, dimension(:,:,:) :: OCDP      => NULL()
    real, pointer, dimension(:,:,:) :: OCSV      => NULL()
    real, pointer, dimension(:,:,:) :: OCWT      => NULL()
    real, pointer, dimension(:,:,:) :: OCSD      => NULL()
    real, pointer, dimension(:,:,:) :: SUDP      => NULL()
    real, pointer, dimension(:,:,:) :: SUSV      => NULL()
    real, pointer, dimension(:,:,:) :: SUWT      => NULL()
    real, pointer, dimension(:,:,:) :: SUSD      => NULL()
    real, pointer, dimension(:,:,:) :: SSDP      => NULL()
    real, pointer, dimension(:,:,:) :: SSSV      => NULL()
    real, pointer, dimension(:,:,:) :: SSWT      => NULL()
    real, pointer, dimension(:,:,:) :: SSSD      => NULL()
    real, pointer, dimension(:,:)   :: DTSDT     => NULL()

! Pointers to internals

    real, pointer, dimension(:,:) :: TS   => NULL()
    real, pointer, dimension(:,:) :: QS   => NULL()
    real, pointer, dimension(:,:) :: CM   => NULL()
    real, pointer, dimension(:,:) :: CT   => NULL()
    real, pointer, dimension(:,:) :: CQ   => NULL()
    real, pointer, dimension(:,:) :: CN   => NULL()
    real, pointer, dimension(:,:) :: TH   => NULL()
    real, pointer, dimension(:,:) :: QH   => NULL()
    real, pointer, dimension(:,:) :: UH   => NULL()
    real, pointer, dimension(:,:) :: VH   => NULL()
    real, pointer, dimension(:,:) :: RHOS => NULL()
    real, pointer, dimension(:,:) :: D0   => NULL()

! Pointers to exports

    real, pointer, dimension(:,:) :: LST       => NULL()
    real, pointer, dimension(:,:) :: FRI       => NULL()
    real, pointer, dimension(:,:) :: OFRI      => NULL()
    real, pointer, dimension(:,:) :: EMISS     => NULL()
    real, pointer, dimension(:,:) :: ALBVR     => NULL()
    real, pointer, dimension(:,:) :: ALBVF     => NULL()
    real, pointer, dimension(:,:) :: ALBNF     => NULL()
    real, pointer, dimension(:,:) :: ALBNR     => NULL()
    real, pointer, dimension(:,:) :: DELSS     => NULL()
    real, pointer, dimension(:,:) :: DELUS     => NULL()
    real, pointer, dimension(:,:) :: DELVS     => NULL()
    real, pointer, dimension(:,:) :: DELTS     => NULL()
    real, pointer, dimension(:,:) :: DELQS     => NULL()
    real, pointer, dimension(:,:) :: DLQLL     => NULL()
    real, pointer, dimension(:,:) :: DLQIL     => NULL()
    real, pointer, dimension(:,:) :: TSOIL1    => NULL()
    real, pointer, dimension(:,:) :: TSOIL2    => NULL()
    real, pointer, dimension(:,:) :: TSOIL3    => NULL()
    real, pointer, dimension(:,:) :: TSOIL4    => NULL()
    real, pointer, dimension(:,:) :: TSOIL5    => NULL()
    real, pointer, dimension(:,:) :: TSOIL6    => NULL()
    real, pointer, dimension(:,:) :: ASNOW     => NULL()
    real, pointer, dimension(:,:) :: SHSNOW    => NULL()
    real, pointer, dimension(:,:) :: AVETSNOW  => NULL()
    real, pointer, dimension(:,:) :: TPSNO     => NULL()
    real, pointer, dimension(:,:) :: TPUST     => NULL()
    real, pointer, dimension(:,:) :: TPSAT     => NULL()
    real, pointer, dimension(:,:) :: TPWLT     => NULL()
    real, pointer, dimension(:,:) :: TPSURF    => NULL()
    real, pointer, dimension(:,:) :: FRSAT     => NULL()
    real, pointer, dimension(:,:) :: FRUST     => NULL()
    real, pointer, dimension(:,:) :: FRWLT     => NULL()
    real, pointer, dimension(:,:) :: SNOMAS    => NULL()
    real, pointer, dimension(:,:) :: SNOWDP    => NULL()
    real, pointer, dimension(:,:) :: WET1      => NULL()
    real, pointer, dimension(:,:) :: WET1_FOR_CHEM => NULL()
    real, pointer, dimension(:,:) :: WET2      => NULL()
    real, pointer, dimension(:,:) :: WET3      => NULL()
    real, pointer, dimension(:,:) :: WCSF      => NULL()
    real, pointer, dimension(:,:) :: WCRZ      => NULL()
    real, pointer, dimension(:,:) :: WCPR      => NULL()
    real, pointer, dimension(:,:) :: WESNN1    => NULL()
    real, pointer, dimension(:,:) :: WESNN2    => NULL()
    real, pointer, dimension(:,:) :: WESNN3    => NULL()
    real, pointer, dimension(:,:) :: CAPAC     => NULL()
    real, pointer, dimension(:,:) :: TAUXO     => NULL()
    real, pointer, dimension(:,:) :: TAUYO     => NULL()
    real, pointer, dimension(:,:) :: EVAPO     => NULL()
    real, pointer, dimension(:,:) :: SHO       => NULL()
    real, pointer, dimension(:,:) :: USTAR     => NULL()
    real, pointer, dimension(:,:) :: TSTAR     => NULL()
    real, pointer, dimension(:,:) :: QSTAR     => NULL()
    real, pointer, dimension(:,:) :: U10M      => NULL()
    real, pointer, dimension(:,:) :: V10M      => NULL()
    real, pointer, dimension(:,:) :: U10N      => NULL()
    real, pointer, dimension(:,:) :: V10N      => NULL()
    real, pointer, dimension(:,:) :: U50M      => NULL()
    real, pointer, dimension(:,:) :: V50M      => NULL()
    real, pointer, dimension(:,:) :: T10M      => NULL()
    real, pointer, dimension(:,:) :: Q10M      => NULL()
    real, pointer, dimension(:,:) :: U2M       => NULL()
    real, pointer, dimension(:,:) :: V2M       => NULL()
    real, pointer, dimension(:,:) :: T2M       => NULL()
    real, pointer, dimension(:,:) :: Q2M       => NULL()
    real, pointer, dimension(:,:) :: UAX       => NULL()
    real, pointer, dimension(:,:) :: VAX       => NULL()
    real, pointer, dimension(:,:) :: TA        => NULL()
    real, pointer, dimension(:,:) :: QA        => NULL()
    real, pointer, dimension(:,:) :: LWI       => NULL()
    real, pointer, dimension(:,:) :: FROCEAN   => NULL()
    real, pointer, dimension(:,:) :: FRLAKE    => NULL()
    real, pointer, dimension(:,:) :: FRLAND    => NULL()
    real, pointer, dimension(:,:) :: FRLANDICE => NULL()
    real, pointer, dimension(:,:) :: HLATN     => NULL()
    real, pointer, dimension(:,:) :: PS_       => NULL()



    real, pointer, dimension(:,:) :: HLATWTR  => NULL()
    real, pointer, dimension(:,:) :: HLATICE  => NULL()
    real, pointer, dimension(:,:) :: SHWTR    => NULL()
    real, pointer, dimension(:,:) :: SHICE    => NULL()
    real, pointer, dimension(:,:) :: TAUXW    => NULL()
    real, pointer, dimension(:,:) :: TAUXI    => NULL()
    real, pointer, dimension(:,:) :: TAUYW    => NULL()
    real, pointer, dimension(:,:) :: TAUYI    => NULL()
    real, pointer, dimension(:,:) :: LWNDWTR  => NULL()
    real, pointer, dimension(:,:) :: SWNDWTR  => NULL()
    real, pointer, dimension(:,:) :: LWNDICE  => NULL()
    real, pointer, dimension(:,:) :: SWNDICE  => NULL()
    real, pointer, dimension(:,:) :: SNOWOCN  => NULL()
    real, pointer, dimension(:,:) :: ICEFOCN  => NULL()
    real, pointer, dimension(:,:) :: SPTOTOCN => NULL()
    real, pointer, dimension(:,:) :: RAINOCN  => NULL()
    real, pointer, dimension(:,:) :: TSKINW   => NULL()
    real, pointer, dimension(:,:) :: TSKINICE => NULL()

    real, pointer, dimension(:,:) :: DCOOL    => NULL()
    real, pointer, dimension(:,:) :: DWARM    => NULL()
    real, pointer, dimension(:,:) :: TDROP    => NULL()
    real, pointer, dimension(:,:) :: QCOOL    => NULL()
    real, pointer, dimension(:,:) :: SWCOOL   => NULL()
    real, pointer, dimension(:,:) :: USTARW   => NULL()
    real, pointer, dimension(:,:) :: TBAR     => NULL()
    real, pointer, dimension(:,:) :: LCOOL    => NULL()
    real, pointer, dimension(:,:) :: BCOOL    => NULL()
    real, pointer, dimension(:,:) :: TDEL     => NULL()
    real, pointer, dimension(:,:) :: TS_FOUND => NULL()
    real, pointer, dimension(:,:) :: SS_FOUND => NULL()
    real, pointer, dimension(:,:) :: QWARM    => NULL()
    real, pointer, dimension(:,:) :: SWWARM   => NULL()
    real, pointer, dimension(:,:) :: LANGM    => NULL()
    real, pointer, dimension(:,:) :: PHIW     => NULL()
    real, pointer, dimension(:,:) :: TAUTW    => NULL()
    real, pointer, dimension(:,:) :: ZETA_W   => NULL()
    real, pointer, dimension(:,:) :: TWMTF    => NULL()

    real, pointer, dimension(:,:) :: HICE       => NULL()
    real, pointer, dimension(:,:) :: HSNO       => NULL()
    real, pointer, dimension(:,:) :: FRZMLT     => NULL()
    real, pointer, dimension(:,:) :: TSKINWCICE => NULL()
    real, pointer, dimension(:,:) :: ISTSFC     => NULL()
    real, pointer, dimension(:,:) :: SSKINW     => NULL()
    real, pointer, dimension(:,:) :: MELTT      => NULL()
    real, pointer, dimension(:,:) :: MELTB      => NULL()
    real, pointer, dimension(:,:) :: MELTS      => NULL()
    real, pointer, dimension(:,:) :: MELTL      => NULL()
    real, pointer, dimension(:,:) :: FRAZIL     => NULL()
    real, pointer, dimension(:,:) :: CONGEL     => NULL()
    real, pointer, dimension(:,:) :: SNOICE     => NULL()
    real, pointer, dimension(:,:) :: DAIDTT     => NULL()
    real, pointer, dimension(:,:) :: DVIDTT     => NULL()
    real, pointer, dimension(:,:) :: DAIDTD     => NULL()
    real, pointer, dimension(:,:) :: DVIDTD     => NULL()
    real, pointer, dimension(:,:) :: FBOT       => NULL()
    real, pointer, dimension(:,:) :: HFLUX      => NULL()
    real, pointer, dimension(:,:) :: WFLUX      => NULL()
    real, pointer, dimension(:,:) :: SFLUX      => NULL()
    real, pointer, dimension(:,:) :: FSWTHRU    => NULL()
    real, pointer, dimension(:,:) :: FSWABS     => NULL()
    real, pointer, dimension(:,:) :: USTARI     => NULL()
    real, pointer, dimension(:,:) :: FHOCN      => NULL()

    real, pointer, dimension(:,:) :: EVAPOU       => NULL()
    real, pointer, dimension(:,:) :: SUBLIM       => NULL()
    real, pointer, dimension(:,:) :: SHOU         => NULL()
    real, pointer, dimension(:,:) :: HLWUP        => NULL()
    real, pointer, dimension(:,:) :: LWNDSRF      => NULL()
    real, pointer, dimension(:,:) :: SWNDSRF      => NULL()
    real, pointer, dimension(:,:) :: RUNOFF       => NULL()
    real, pointer, dimension(:,:) :: DISCHARGE    => NULL()
    real, pointer, dimension(:,:) :: RUNSURF      => NULL()
    real, pointer, dimension(:,:) :: BASEFLOW     => NULL()
    real, pointer, dimension(:,:) :: ACCUM        => NULL()
    real, pointer, dimension(:,:) :: SMELT        => NULL()
    real, pointer, dimension(:,:) :: EVEG         => NULL()
    real, pointer, dimension(:,:) :: EINT         => NULL()
    real, pointer, dimension(:,:) :: EICE         => NULL()
    real, pointer, dimension(:,:) :: ESNO         => NULL()
    real, pointer, dimension(:,:) :: ESOI         => NULL()
    real, pointer, dimension(:,:) :: WAT10CM      => NULL()
    real, pointer, dimension(:,:) :: WATSOI       => NULL()
    real, pointer, dimension(:,:) :: ICESOI       => NULL()
    real, pointer, dimension(:,:) :: EVLAND       => NULL()
    real, pointer, dimension(:,:) :: PRLAND       => NULL()
    real, pointer, dimension(:,:) :: SNOLAND      => NULL()
    real, pointer, dimension(:,:) :: DRPARLAND    => NULL()
    real, pointer, dimension(:,:) :: DFPARLAND    => NULL()
    real, pointer, dimension(:,:) :: LHSNOW       => NULL()
    real, pointer, dimension(:,:) :: TCSORIG      => NULL()
    real, pointer, dimension(:,:) :: TPSN1IN      => NULL()
    real, pointer, dimension(:,:) :: TPSN1OUT     => NULL()
    real, pointer, dimension(:,:) :: SWNETSNOW    => NULL()
    real, pointer, dimension(:,:) :: LWUPSNOW     => NULL()
    real, pointer, dimension(:,:) :: LWDNSNOW     => NULL()
    real, pointer, dimension(:,:) :: GHSNOW       => NULL()
    real, pointer, dimension(:,:) :: LHLAND       => NULL()
    real, pointer, dimension(:,:) :: SHLAND       => NULL()
    real, pointer, dimension(:,:) :: SWLAND       => NULL()
    real, pointer, dimension(:,:) :: SWDOWNLAND   => NULL()
    real, pointer, dimension(:,:) :: LWLAND       => NULL()
    real, pointer, dimension(:,:) :: GHLAND       => NULL()
    real, pointer, dimension(:,:) :: GHTSKIN      => NULL()
    real, pointer, dimension(:,:) :: SMLAND       => NULL()
    real, pointer, dimension(:,:) :: QINFIL       => NULL()
    real, pointer, dimension(:,:) :: TWLAND       => NULL()
    real, pointer, dimension(:,:) :: TELAND       => NULL()
    real, pointer, dimension(:,:) :: TSLAND       => NULL()
    real, pointer, dimension(:,:) :: DWLAND       => NULL()
    real, pointer, dimension(:,:) :: DHLAND       => NULL()
    real, pointer, dimension(:,:) :: SPLAND       => NULL()
    real, pointer, dimension(:,:) :: SPLH         => NULL()
    real, pointer, dimension(:,:) :: SPWATR       => NULL()
    real, pointer, dimension(:,:) :: SPSNOW       => NULL()
    real, pointer, dimension(:,:) :: RCU          => NULL()
    real, pointer, dimension(:,:) :: RLS          => NULL()
    real, pointer, dimension(:,:) :: SNO          => NULL()
    real, pointer, dimension(:,:) :: ICE          => NULL()
    real, pointer, dimension(:,:) :: FRZR         => NULL()
    real, pointer, dimension(:,:) :: TPREC        => NULL()
    real, pointer, dimension(:,:) :: CN_PRCP      => NULL()
    real, pointer, dimension(:,:) :: PRECTOT      => NULL()
    real, pointer, dimension(:,:) :: PRECCU       => NULL()
    real, pointer, dimension(:,:) :: T2MDEW       => NULL()
    real, pointer, dimension(:,:) :: T2MWET       => NULL()

! GOSWIM (internal/export variables from catch/catchcn)
    real, pointer, dimension(:,:,:) :: RDU001      => NULL()
    real, pointer, dimension(:,:,:) :: RDU002      => NULL()
    real, pointer, dimension(:,:,:) :: RDU003      => NULL()
    real, pointer, dimension(:,:,:) :: RDU004      => NULL()
    real, pointer, dimension(:,:,:) :: RDU005      => NULL()
    real, pointer, dimension(:,:,:) :: RBC001      => NULL()
    real, pointer, dimension(:,:,:) :: RBC002      => NULL()
    real, pointer, dimension(:,:,:) :: ROC001      => NULL()
    real, pointer, dimension(:,:,:) :: ROC002      => NULL()
    real, pointer, dimension(:,:)   :: RMELTDU001  => NULL()
    real, pointer, dimension(:,:)   :: RMELTDU002  => NULL()
    real, pointer, dimension(:,:)   :: RMELTDU003  => NULL()
    real, pointer, dimension(:,:)   :: RMELTDU004  => NULL()
    real, pointer, dimension(:,:)   :: RMELTDU005  => NULL()
    real, pointer, dimension(:,:)   :: RMELTBC001  => NULL()
    real, pointer, dimension(:,:)   :: RMELTBC002  => NULL()
    real, pointer, dimension(:,:)   :: RMELTOC001  => NULL()
    real, pointer, dimension(:,:)   :: RMELTOC002  => NULL()
    real, pointer, dimension(:,:)   :: PEATCLSM_WATERLEVEL => NULL()
    real, pointer, dimension(:,:)   :: PEATCLSM_FSWCHANGE  => NULL()
    real, pointer, dimension(:,:)   :: DZGT1       => NULL()
    real, pointer, dimension(:,:)   :: DZGT2       => NULL()
    real, pointer, dimension(:,:)   :: DZGT3       => NULL()
    real, pointer, dimension(:,:)   :: DZGT4       => NULL()
    real, pointer, dimension(:,:)   :: DZGT5       => NULL()
    real, pointer, dimension(:,:)   :: DZGT6       => NULL()
    real, pointer, dimension(:,:)   :: DZPR        => NULL()
    real, pointer, dimension(:,:)   :: DZRZ        => NULL()
    real, pointer, dimension(:,:)   :: DZSF        => NULL()
    real, pointer, dimension(:,:)   :: DZTS        => NULL()
    real, pointer, dimension(:,:)   :: WPWET       => NULL()
    real, pointer, dimension(:,:)   :: WPEMW       => NULL()
    real, pointer, dimension(:,:)   :: WPMC        => NULL()
    real, pointer, dimension(:,:)   :: CDCR2       => NULL()
    real, pointer, dimension(:,:)   :: POROS       => NULL()

! CN model
    real, pointer, dimension(:,:) :: CNLAI       => NULL()
    real, pointer, dimension(:,:) :: CNTLAI      => NULL()
    real, pointer, dimension(:,:) :: CNSAI       => NULL()
    real, pointer, dimension(:,:) :: CNTOTC      => NULL()
    real, pointer, dimension(:,:) :: CNVEGC      => NULL()
    real, pointer, dimension(:,:) :: CNROOT      => NULL()
    real, pointer, dimension(:,:) :: CNFROOTC    => NULL()
    real, pointer, dimension(:,:) :: CNNPP       => NULL()
    real, pointer, dimension(:,:) :: CNGPP       => NULL()
    real, pointer, dimension(:,:) :: CNSR        => NULL()
    real, pointer, dimension(:,:) :: CNNEE       => NULL()
    real, pointer, dimension(:,:) :: CNXSMR      => NULL()
    real, pointer, dimension(:,:) :: CNADD       => NULL()
    real, pointer, dimension(:,:) :: CNLOSS      => NULL()
    real, pointer, dimension(:,:) :: CNBURN      => NULL()
    real, pointer, dimension(:,:) :: PARABS      => NULL()
    real, pointer, dimension(:,:) :: PARINC      => NULL()
    real, pointer, dimension(:,:) :: SCSAT       => NULL()
    real, pointer, dimension(:,:) :: SCUNS       => NULL()
    real, pointer, dimension(:,:) :: BTRANT      => NULL()
    real, pointer, dimension(:,:) :: SIF         => NULL()
    real, pointer, dimension(:,:) :: CNFSEL      => NULL()

! Fire danger
    real, pointer, dimension(:,:) :: FFMC        => NULL()
    real, pointer, dimension(:,:) :: GFMC        => NULL()
    real, pointer, dimension(:,:) :: DMC         => NULL()
    real, pointer, dimension(:,:) :: DC          => NULL()
    real, pointer, dimension(:,:) :: ISI         => NULL()
    real, pointer, dimension(:,:) :: BUI         => NULL()
    real, pointer, dimension(:,:) :: FWI         => NULL()
    real, pointer, dimension(:,:) :: DSR         => NULL()

    real, pointer, dimension(:,:) :: FFMC_DAILY  => NULL()
    real, pointer, dimension(:,:) :: DMC_DAILY   => NULL()
    real, pointer, dimension(:,:) :: DC_DAILY    => NULL()
    real, pointer, dimension(:,:) :: ISI_DAILY   => NULL()
    real, pointer, dimension(:,:) :: BUI_DAILY   => NULL()
    real, pointer, dimension(:,:) :: FWI_DAILY   => NULL()
    real, pointer, dimension(:,:) :: DSR_DAILY   => NULL()

    real, pointer, dimension(:,:) :: FFMC_DAILY_ => NULL()
    real, pointer, dimension(:,:) :: DMC_DAILY_  => NULL()
    real, pointer, dimension(:,:) :: DC_DAILY_   => NULL()
    real, pointer, dimension(:,:) :: ISI_DAILY_  => NULL()
    real, pointer, dimension(:,:) :: BUI_DAILY_  => NULL()
    real, pointer, dimension(:,:) :: FWI_DAILY_  => NULL()
    real, pointer, dimension(:,:) :: DSR_DAILY_  => NULL()

    real, pointer, dimension(:,:) :: VPD         => NULL()

!   These are the tile versions of the imports

    real, pointer, dimension(:) :: PSTILE          => NULL()
    real, pointer, dimension(:) :: DZTILE          => NULL()
    real, pointer, dimension(:) :: UUTILE          => NULL()
    real, pointer, dimension(:) :: EVPTILE         => NULL()
    real, pointer, dimension(:) :: SHFTILE         => NULL()
    real, pointer, dimension(:) :: DEVTILE         => NULL()
    real, pointer, dimension(:) :: DSHTILE         => NULL()
    real, pointer, dimension(:) :: TMPTILE         => NULL()
    real, pointer, dimension(:) :: PCUTILE         => NULL()
    real, pointer, dimension(:) :: PLSTILE         => NULL()
    real, pointer, dimension(:) :: SNOFLTILE       => NULL()
    real, pointer, dimension(:) :: ICEFLTILE       => NULL()
    real, pointer, dimension(:) :: FRZRFLTILE      => NULL()
    real, pointer, dimension(:) :: TAUXTILE        => NULL()
    real, pointer, dimension(:) :: TAUYTILE        => NULL()
    real, pointer, dimension(:) :: DFPTILE         => NULL()
    real, pointer, dimension(:) :: DRPTILE         => NULL()
    real, pointer, dimension(:) :: DFNTILE         => NULL()
    real, pointer, dimension(:) :: DRNTILE         => NULL()
    real, pointer, dimension(:) :: DFUTILE         => NULL()
    real, pointer, dimension(:) :: DRUTILE         => NULL()
    real, pointer, dimension(:) :: LWBTILE         => NULL()
    real, pointer, dimension(:) :: ALWTILE         => NULL()
    real, pointer, dimension(:) :: BLWTILE         => NULL()
    real, pointer, dimension(:,:) :: DUDPTILE      => NULL()
    real, pointer, dimension(:,:) :: DUSVTILE      => NULL()
    real, pointer, dimension(:,:) :: DUWTTILE      => NULL()
    real, pointer, dimension(:,:) :: DUSDTILE      => NULL()
    real, pointer, dimension(:,:) :: BCDPTILE      => NULL()
    real, pointer, dimension(:,:) :: BCSVTILE      => NULL()
    real, pointer, dimension(:,:) :: BCWTTILE      => NULL()
    real, pointer, dimension(:,:) :: BCSDTILE      => NULL()
    real, pointer, dimension(:,:) :: OCDPTILE      => NULL()
    real, pointer, dimension(:,:) :: OCSVTILE      => NULL()
    real, pointer, dimension(:,:) :: OCWTTILE      => NULL()
    real, pointer, dimension(:,:) :: OCSDTILE      => NULL()
    real, pointer, dimension(:,:) :: SUDPTILE      => NULL()
    real, pointer, dimension(:,:) :: SUSVTILE      => NULL()
    real, pointer, dimension(:,:) :: SUWTTILE      => NULL()
    real, pointer, dimension(:,:) :: SUSDTILE      => NULL()
    real, pointer, dimension(:,:) :: SSDPTILE      => NULL()
    real, pointer, dimension(:,:) :: SSSVTILE      => NULL()
    real, pointer, dimension(:,:) :: SSWTTILE      => NULL()
    real, pointer, dimension(:,:) :: SSSDTILE      => NULL()
    real, pointer, dimension(:)   :: DTSDTTILE     => NULL()

! These are tile versions of internals

    real, pointer, dimension(:) :: TSTILE => NULL()
    real, pointer, dimension(:) :: QSTILE => NULL()
    real, pointer, dimension(:) :: THTILE => NULL()
    real, pointer, dimension(:) :: QHTILE => NULL()
    real, pointer, dimension(:) :: UHTILE => NULL()
    real, pointer, dimension(:) :: VHTILE => NULL()
    real, pointer, dimension(:) :: CTTILE => NULL()
    real, pointer, dimension(:) :: CQTILE => NULL()
    real, pointer, dimension(:) :: CMTILE => NULL()

! These are tile versions of exports

    real, pointer, dimension(:) :: LSTTILE      => NULL()
    real, pointer, dimension(:) :: FRTILE       => NULL()
    real, pointer, dimension(:) :: OFRTILE      => NULL()
    real, pointer, dimension(:) :: EMISSTILE    => NULL()
    real, pointer, dimension(:) :: ALBVRTILE    => NULL()
    real, pointer, dimension(:) :: ALBVFTILE    => NULL()
    real, pointer, dimension(:) :: ALBNFTILE    => NULL()
    real, pointer, dimension(:) :: ALBNRTILE    => NULL()
    real, pointer, dimension(:) :: DTSTILE      => NULL()
    real, pointer, dimension(:) :: DQSTILE      => NULL()
    real, pointer, dimension(:) :: TSOIL1TILE   => NULL()
    real, pointer, dimension(:) :: TSOIL2TILE   => NULL()
    real, pointer, dimension(:) :: TSOIL3TILE   => NULL()
    real, pointer, dimension(:) :: TSOIL4TILE   => NULL()
    real, pointer, dimension(:) :: TSOIL5TILE   => NULL()
    real, pointer, dimension(:) :: TSOIL6TILE   => NULL()
    real, pointer, dimension(:) :: ASNOWTILE    => NULL()
    real, pointer, dimension(:) :: SHSNOWTILE   => NULL()
    real, pointer, dimension(:) :: AVETSNOWTILE => NULL()
    real, pointer, dimension(:) :: TPSNOTILE    => NULL()
    real, pointer, dimension(:) :: TPUSTTILE    => NULL()
    real, pointer, dimension(:) :: TPSATTILE    => NULL()
    real, pointer, dimension(:) :: TPWLTTILE    => NULL()
    real, pointer, dimension(:) :: TPSURFTILE   => NULL()
    real, pointer, dimension(:) :: FRSATTILE    => NULL()
    real, pointer, dimension(:) :: FRUSTTILE    => NULL()
    real, pointer, dimension(:) :: FRWLTTILE    => NULL()
    real, pointer, dimension(:) :: SNOWTILE     => NULL()
    real, pointer, dimension(:) :: SNODTILE     => NULL()
    real, pointer, dimension(:) :: WET1TILE     => NULL()
    real, pointer, dimension(:) :: WET2TILE     => NULL()
    real, pointer, dimension(:) :: WET3TILE     => NULL()
    real, pointer, dimension(:) :: WCSFTILE     => NULL()
    real, pointer, dimension(:) :: WCRZTILE     => NULL()
    real, pointer, dimension(:) :: WCPRTILE     => NULL()
    real, pointer, dimension(:) :: WESNN1TILE   => NULL()
    real, pointer, dimension(:) :: WESNN2TILE   => NULL()
    real, pointer, dimension(:) :: WESNN3TILE   => NULL()
    real, pointer, dimension(:) :: CAPACTILE    => NULL()
    real, pointer, dimension(:) :: HLATNTILE    => NULL()

    real, pointer, dimension(:) :: HLATWTRTILE    => NULL()
    real, pointer, dimension(:) :: HLATICETILE    => NULL()
    real, pointer, dimension(:) ::   SHWTRTILE    => NULL()
    real, pointer, dimension(:) ::   SHICETILE    => NULL()
    real, pointer, dimension(:) ::   TAUXWTILE    => NULL()
    real, pointer, dimension(:) ::   TAUXITILE    => NULL()
    real, pointer, dimension(:) ::   TAUYWTILE    => NULL()
    real, pointer, dimension(:) ::   TAUYITILE    => NULL()
    real, pointer, dimension(:) :: LWNDWTRTILE    => NULL()
    real, pointer, dimension(:) :: SWNDWTRTILE    => NULL()
    real, pointer, dimension(:) :: LWNDICETILE    => NULL()
    real, pointer, dimension(:) :: SWNDICETILE    => NULL()
    real, pointer, dimension(:) :: SNOWOCNTILE    => NULL()
    real, pointer, dimension(:) :: ICEFOCNTILE    => NULL()
    real, pointer, dimension(:) ::SPTOTOCNTILE    => NULL()
    real, pointer, dimension(:) :: RAINOCNTILE    => NULL()
    real, pointer, dimension(:) ::  TSKINWTILE    => NULL()
    real, pointer, dimension(:) ::  TSKINICETILE  => NULL()

    real, pointer, dimension(:) :: DCOOL_TILE    => NULL()
    real, pointer, dimension(:) :: DWARM_TILE    => NULL()
    real, pointer, dimension(:) :: TDROP_TILE    => NULL()
    real, pointer, dimension(:) :: QCOOL_TILE    => NULL()
    real, pointer, dimension(:) :: SWCOOL_TILE   => NULL()
    real, pointer, dimension(:) :: USTARW_TILE   => NULL()
    real, pointer, dimension(:) :: TBAR_TILE     => NULL()
    real, pointer, dimension(:) :: LCOOL_TILE    => NULL()
    real, pointer, dimension(:) :: BCOOL_TILE    => NULL()
    real, pointer, dimension(:) :: TDEL_TILE     => NULL()
    real, pointer, dimension(:) :: TS_FOUND_TILE => NULL()
    real, pointer, dimension(:) :: SS_FOUND_TILE => NULL()
    real, pointer, dimension(:) :: QWARM_TILE    => NULL()
    real, pointer, dimension(:) :: SWWARM_TILE   => NULL()
    real, pointer, dimension(:) :: LANGM_TILE    => NULL()
    real, pointer, dimension(:) :: PHIW_TILE     => NULL()
    real, pointer, dimension(:) :: TAUTW_TILE    => NULL()
    real, pointer, dimension(:) :: ZETA_W_TILE   => NULL()
    real, pointer, dimension(:) :: TWMTF_TILE    => NULL()

    real, pointer, dimension(:) ::  HICETILE        => NULL()
    real, pointer, dimension(:) ::  HSNOTILE        => NULL()
    real, pointer, dimension(:) ::  FRZMLTTILE      => NULL()
    real, pointer, dimension(:) ::  TSKINWCICETILE  => NULL()
    real, pointer, dimension(:) ::  ISTSFCTILE      => NULL()
    real, pointer, dimension(:) ::  SSKINWTILE      => NULL()
    real, pointer, dimension(:) ::  MELTTTILE       => NULL()
    real, pointer, dimension(:) ::  MELTBTILE       => NULL()
    real, pointer, dimension(:) ::  MELTSTILE       => NULL()
    real, pointer, dimension(:) ::  MELTLTILE       => NULL()
    real, pointer, dimension(:) ::  FRAZILTILE      => NULL()
    real, pointer, dimension(:) ::  CONGELTILE      => NULL()
    real, pointer, dimension(:) ::  SNOICETILE      => NULL()
    real, pointer, dimension(:) ::  DAIDTTTILE      => NULL()
    real, pointer, dimension(:) ::  DVIDTTTILE      => NULL()
    real, pointer, dimension(:) ::  DAIDTDTILE      => NULL()
    real, pointer, dimension(:) ::  DVIDTDTILE      => NULL()
    real, pointer, dimension(:) ::  FBOTTILE        => NULL()
    real, pointer, dimension(:) ::  HFLUXTILE       => NULL()
    real, pointer, dimension(:) ::  WFLUXTILE       => NULL()
    real, pointer, dimension(:) ::  SFLUXTILE       => NULL()
    real, pointer, dimension(:) ::  FSWTHRUTILE     => NULL()
    real, pointer, dimension(:) ::  FSWABSTILE      => NULL()
    real, pointer, dimension(:) ::  USTARITILE      => NULL()
    real, pointer, dimension(:) ::  FHOCNTILE       => NULL()

    real, pointer, dimension(:) :: EVAPOUTILE       => NULL()
    real, pointer, dimension(:) :: SUBLIMTILE       => NULL()
    real, pointer, dimension(:) :: SHOUTILE         => NULL()
    real, pointer, dimension(:) :: HLWUPTILE        => NULL()
    real, pointer, dimension(:) :: LWNDSRFTILE      => NULL()
    real, pointer, dimension(:) :: SWNDSRFTILE      => NULL()
    real, pointer, dimension(:) :: RUNOFFTILE       => NULL()
    real, pointer, dimension(:) :: RUNSURFTILE      => NULL()
    real, pointer, dimension(:) :: DISCHARGETILE    => NULL()
    real, pointer, dimension(:) :: BASEFLOWTILE     => NULL()
    real, pointer, dimension(:) :: ACCUMTILE        => NULL()
    real, pointer, dimension(:) :: SMELTTILE        => NULL()
    real, pointer, dimension(:) :: EVEGTILE         => NULL()
    real, pointer, dimension(:) :: EINTTILE         => NULL()
    real, pointer, dimension(:) :: EICETILE         => NULL()
    real, pointer, dimension(:) :: ESNOTILE         => NULL()
    real, pointer, dimension(:) :: ESOITILE         => NULL()
    real, pointer, dimension(:) :: WAT10CMTILE      => NULL()
    real, pointer, dimension(:) :: WATSOITILE       => NULL()
    real, pointer, dimension(:) :: ICESOITILE       => NULL()
    real, pointer, dimension(:) :: EVLANDTILE       => NULL()
    real, pointer, dimension(:) :: PRLANDTILE       => NULL()
    real, pointer, dimension(:) :: SNOLANDTILE      => NULL()
    real, pointer, dimension(:) :: DRPARLANDTILE    => NULL()
    real, pointer, dimension(:) :: DFPARLANDTILE    => NULL()
    real, pointer, dimension(:) :: LHSNOWTILE       => NULL()
    real, pointer, dimension(:) :: SWNETSNOWTILE    => NULL()
    real, pointer, dimension(:) :: LWUPSNOWTILE     => NULL()
    real, pointer, dimension(:) :: LWDNSNOWTILE     => NULL()
    real, pointer, dimension(:) :: TCSORIGTILE      => NULL()
    real, pointer, dimension(:) :: TPSN1INTILE      => NULL()
    real, pointer, dimension(:) :: TPSN1OUTTILE     => NULL()
    real, pointer, dimension(:) :: GHSNOWTILE       => NULL()
    real, pointer, dimension(:) :: LHLANDTILE       => NULL()
    real, pointer, dimension(:) :: SHLANDTILE       => NULL()
    real, pointer, dimension(:) :: SWLANDTILE       => NULL()
    real, pointer, dimension(:) :: SWDOWNLANDTILE   => NULL()
    real, pointer, dimension(:) :: LWLANDTILE       => NULL()
    real, pointer, dimension(:) :: GHLANDTILE       => NULL()
    real, pointer, dimension(:) :: GHTSKINTILE      => NULL()
    real, pointer, dimension(:) :: SMLANDTILE       => NULL()
    real, pointer, dimension(:) :: QINFILTILE       => NULL()
    real, pointer, dimension(:) :: TWLANDTILE       => NULL()
    real, pointer, dimension(:) :: TELANDTILE       => NULL()
    real, pointer, dimension(:) :: TSLANDTILE       => NULL()
    real, pointer, dimension(:) :: DWLANDTILE       => NULL()
    real, pointer, dimension(:) :: DHLANDTILE       => NULL()
    real, pointer, dimension(:) :: SPLANDTILE       => NULL()
    real, pointer, dimension(:) :: SPLHTILE         => NULL()
    real, pointer, dimension(:) :: SPWATRTILE       => NULL()
    real, pointer, dimension(:) :: SPSNOWTILE       => NULL()
    real, pointer, dimension(:,:) :: RDU001TILE     => NULL()
    real, pointer, dimension(:,:) :: RDU002TILE     => NULL()
    real, pointer, dimension(:,:) :: RDU003TILE     => NULL()
    real, pointer, dimension(:,:) :: RDU004TILE     => NULL()
    real, pointer, dimension(:,:) :: RDU005TILE     => NULL()
    real, pointer, dimension(:,:) :: RBC001TILE     => NULL()
    real, pointer, dimension(:,:) :: RBC002TILE     => NULL()
    real, pointer, dimension(:,:) :: ROC001TILE     => NULL()
    real, pointer, dimension(:,:) :: ROC002TILE     => NULL()
    real, pointer, dimension(:) :: RMELTDU001TILE   => NULL()
    real, pointer, dimension(:) :: RMELTDU002TILE   => NULL()
    real, pointer, dimension(:) :: RMELTDU003TILE   => NULL()
    real, pointer, dimension(:) :: RMELTDU004TILE   => NULL()
    real, pointer, dimension(:) :: RMELTDU005TILE   => NULL()
    real, pointer, dimension(:) :: RMELTBC001TILE   => NULL()
    real, pointer, dimension(:) :: RMELTBC002TILE   => NULL()
    real, pointer, dimension(:) :: RMELTOC001TILE   => NULL()
    real, pointer, dimension(:) :: RMELTOC002TILE   => NULL()
    real, pointer, dimension(:) :: PEATCLSM_WATERLEVELTILE  => NULL()
    real, pointer, dimension(:) :: PEATCLSM_FSWCHANGETILE   => NULL()
    real, pointer, dimension(:) :: DZGT1TILE        => NULL()
    real, pointer, dimension(:) :: DZGT2TILE        => NULL()
    real, pointer, dimension(:) :: DZGT3TILE        => NULL()
    real, pointer, dimension(:) :: DZGT4TILE        => NULL()
    real, pointer, dimension(:) :: DZGT5TILE        => NULL()
    real, pointer, dimension(:) :: DZGT6TILE        => NULL()
    real, pointer, dimension(:) :: DZPRTILE         => NULL()
    real, pointer, dimension(:) :: DZRZTILE         => NULL()
    real, pointer, dimension(:) :: DZSFTILE         => NULL()
    real, pointer, dimension(:) :: DZTSTILE         => NULL()
    real, pointer, dimension(:) :: WPWETTILE        => NULL()
    real, pointer, dimension(:) :: WPEMWTILE        => NULL()
    real, pointer, dimension(:) :: WPMCTILE         => NULL()
    real, pointer, dimension(:) :: CDCR2TILE        => NULL()
    real, pointer, dimension(:) :: POROSTILE        => NULL()


    real, pointer, dimension(:) :: CNLAITILE        => NULL()
    real, pointer, dimension(:) :: CNTLAITILE       => NULL()
    real, pointer, dimension(:) :: CNSAITILE        => NULL()
    real, pointer, dimension(:) :: CNTOTCTILE       => NULL()
    real, pointer, dimension(:) :: CNVEGCTILE       => NULL()
    real, pointer, dimension(:) :: CNROOTTILE       => NULL()
    real, pointer, dimension(:) :: CNFROOTCTILE     => NULL()
    real, pointer, dimension(:) :: CNNPPTILE        => NULL()
    real, pointer, dimension(:) :: CNGPPTILE        => NULL()
    real, pointer, dimension(:) :: CNSRTILE         => NULL()
    real, pointer, dimension(:) :: CNNEETILE        => NULL()
    real, pointer, dimension(:) :: CNXSMRTILE       => NULL()
    real, pointer, dimension(:) :: CNADDTILE        => NULL()
    real, pointer, dimension(:) :: CNLOSSTILE       => NULL()
    real, pointer, dimension(:) :: CNBURNTILE       => NULL()
    real, pointer, dimension(:) :: PARABSTILE       => NULL()
    real, pointer, dimension(:) :: PARINCTILE       => NULL()
    real, pointer, dimension(:) :: SCSATTILE        => NULL()
    real, pointer, dimension(:) :: SCUNSTILE        => NULL()
    real, pointer, dimension(:) :: BTRANTTILE       => NULL()
    real, pointer, dimension(:) :: SIFTILE          => NULL()
    real, pointer, dimension(:) :: CNFSELTILE       => NULL()

    real, pointer, dimension(:) :: SLITILE      => NULL()
    real, pointer, dimension(:) :: ZTHTILE      => NULL()

! Fire danger
    real, pointer, dimension(:) :: FFMCTILE       => NULL()
    real, pointer, dimension(:) :: GFMCTILE       => NULL()
    real, pointer, dimension(:) :: DMCTILE        => NULL()
    real, pointer, dimension(:) :: DCTILE         => NULL()
    real, pointer, dimension(:) :: ISITILE        => NULL()
    real, pointer, dimension(:) :: BUITILE        => NULL()
    real, pointer, dimension(:) :: FWITILE        => NULL()
    real, pointer, dimension(:) :: DSRTILE        => NULL()

    real, pointer, dimension(:) :: FFMCDAILYTILE  => NULL()
    real, pointer, dimension(:) :: DMCDAILYTILE   => NULL()
    real, pointer, dimension(:) :: DCDAILYTILE    => NULL()
    real, pointer, dimension(:) :: ISIDAILYTILE   => NULL()
    real, pointer, dimension(:) :: BUIDAILYTILE   => NULL()
    real, pointer, dimension(:) :: FWIDAILYTILE   => NULL()
    real, pointer, dimension(:) :: DSRDAILYTILE   => NULL()

    real, pointer, dimension(:) :: FFMCDAILYTILE_ => NULL()
    real, pointer, dimension(:) :: DMCDAILYTILE_  => NULL()
    real, pointer, dimension(:) :: DCDAILYTILE_   => NULL()
    real, pointer, dimension(:) :: ISIDAILYTILE_  => NULL()
    real, pointer, dimension(:) :: BUIDAILYTILE_  => NULL()
    real, pointer, dimension(:) :: FWIDAILYTILE_  => NULL()
    real, pointer, dimension(:) :: DSRDAILYTILE_  => NULL()

    real, pointer, dimension(:) :: VPDTILE        => NULL()


    real, pointer, dimension(:,:) :: TMP => NULL()
    real, pointer, dimension(:,:) :: TTM => NULL()
    real, pointer, dimension(:,:) :: QTM => NULL()
    real, pointer, dimension(:,:) :: Z0  => NULL()
    real, pointer, dimension(:,:) :: FAC => NULL()
    real, pointer, dimension(:,:) :: TAU => NULL()
    real, pointer, dimension(:,:) :: DTS => NULL()
    real, pointer, dimension(:,:) :: DQS => NULL()

    real, pointer, dimension(:,:) :: DRPAR => NULL()
    real, pointer, dimension(:,:) :: DFPAR => NULL()
    real, pointer, dimension(:,:) :: DRNIR => NULL()
    real, pointer, dimension(:,:) :: DFNIR => NULL()
    real, pointer, dimension(:,:) :: DRUVR => NULL()
    real, pointer, dimension(:,:) :: DFUVR => NULL()
    real, pointer, dimension(:,:) :: ZTH   => NULL()
    real, pointer, dimension(:,:) :: SLR   => NULL()

! following three active only when DO_OBIO==1 or ATM_CO2 == ATM_CO2_FOUR (=4)
!   IMPORTS
    real, pointer, dimension(:,:)   :: CO2SC     => NULL()
    real, pointer, dimension(:,:,:) :: DRBAND    => NULL()
    real, pointer, dimension(:,:,:) :: DFBAND    => NULL()

!   tiled versio of IMPORTS
    real, pointer, dimension(:)   :: CO2SCTILE     => NULL()
    real, pointer, dimension(:,:) :: DRBANDTILE    => NULL()
    real, pointer, dimension(:,:) :: DFBANDTILE    => NULL()
!
    real, pointer, dimension(:,:) :: DISCHARGE_IM    => NULL()

! for reading "forced" precip
    real, pointer, dimension(:,:)           :: PTTe => NULL()
    Integer                                 :: fieldcount
    Type(esmf_field)                        :: bundle_field
    Character(len=ESMF_MAXSTR), allocatable :: fieldnames(:)

! interpolate wind for wind stress
    real, pointer, dimension(:,:) :: UUA     => NULL()
    real, pointer, dimension(:,:) :: VVA     => NULL()
    real, pointer, dimension(:  ) :: UUATILE => NULL()
    real, pointer, dimension(:  ) :: VVATILE => NULL()

    real, pointer                           :: UU10M    (:,:) => NULL()
    real, pointer                           :: RH2M     (:,:) => NULL()
    real, pointer                           :: ALHX     (:,:) => NULL()
    integer                       :: I, N
    integer :: YEAR, MONTH, DAY, HR, MN, SE

    integer  :: iUseInterp
    logical  :: UseInterp

    integer  :: USE_PP_TAPER
    real     :: PP_TAPER_LAT_LOW,PP_TAPER_LAT_HIGH, FACT
    INTEGER                                 :: LSM_CHOICE
    real, allocatable :: PCSCALE(:,:)
    real, allocatable :: PRECSUM(:,:)
    character(len=ESMF_MAXPATHLEN) :: SolCycFileName
    logical :: PersistSolar
    logical :: allocateRunoff

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, name=COMP_NAME, VM=VMG, Grid=GRID, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Run2"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Start Total timer
!------------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"-RUN2" )

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,             &
         LOCSTREAM = LOCSTREAM,                  &
         GCS       = GCS,                        &
         GCNAMES   = GCNAMES,                    &
         GIM       = GIM,                        &
         GEX       = GEX,                        &
         LATS      = LATS,                       &
         LONS      = LONS,                       &
         TILELATS  = tileLATS,                   &
         TILELONS  = tileLONS,                   &
         TILETYPES = TILETYPE,                   &
         ORBIT     = ORBIT,                      &
         INTERNAL_ESMF_STATE = INTERNAL,         &
                                       RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_UserCompGetInternalState(gc, 'SURF_state', wrap, status)
    VERIFY_(STATUS)

    SURF_INTERNAL_STATE => WRAP%PTR

! Get CHOICE OF  Land Surface Model (1:Catch, 2:Catch-CN)
! and Runoff Routing Model (0: OFF, 1: ON)
! -------------------------------------------------------

    call MAPL_GetResource ( MAPL, LSM_CHOICE, Label="LSM_CHOICE:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)

! Pointers to gridded imports
!----------------------------

    call MAPL_GetPointer(IMPORT  , PS      , 'PS'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , DZ      , 'DZ'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , UU      , 'SPEED'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , EVAP    , 'EVAP'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , SH      , 'SH'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , DEVAP   , 'DEVAP'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , DSH     , 'DSH'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , TAUX    , 'TAUX'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , TAUY    , 'TAUY'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , DRPARN  , 'DRPARN' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , DFPARN  , 'DFPARN' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , DRNIRN  , 'DRNIRN' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , DFNIRN  , 'DFNIRN' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , DRUVRN  , 'DRUVRN' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , DFUVRN  , 'DFUVRN' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , LWDNSRF , 'LWDNSRF',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , ALW     , 'ALW'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , BLW     , 'BLW'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , DTSDT   , 'DTSDT' ,   RC=STATUS); VERIFY_(STATUS)

! Horizontal dimensions needed to allocate local arrays
!------------------------------------------------------

    IM = size(PS,1)
    JM = size(PS,2)

    NUM_AERO_DP = 0

    ! Get contents of AERO_DP bundle, put them in an ungridded dimension
    if (DO_GOSWIM) then

       call ESMF_StateGet(      IMPORT, 'AERO_DP',   Bundle   , RC=STATUS); VERIFY_(STATUS)
       call ESMF_FieldBundleGet(Bundle, fieldCOUNT=NUM_AERO_DP, RC=STATUS); VERIFY_(STATUS)

       allocate(AERO_DP_FIELD_NAME(NUM_AERO_DP), STAT=STATUS); VERIFY_(STATUS)
       allocate(AERO_DP(     IM,JM,NUM_AERO_DP), STAT=STATUS); VERIFY_(STATUS)

       do K = 1, NUM_AERO_DP
          call ESMF_FieldBundleGet(Bundle, K, Field,              RC=STATUS); VERIFY_(STATUS)
          call ESMF_FieldGet(Field, NAME=AERO_DP_FIELD_NAME(K),   RC=STATUS); VERIFY_(STATUS)
          call ESMFL_BundleGetPointerToData(Bundle, K, AERO_DP_I, RC=STATUS); VERIFY_(STATUS)
          AERO_DP(:,:,K) = AERO_DP_I
       end do

       ! Sort AERO_DP into individual constituents and cross check
       N_DUDP = 0 ; N_DUSV = 0 ; N_DUWT = 0 ; N_DUSD = 0
       N_BCDP = 0 ; N_BCSV = 0 ; N_BCWT = 0 ; N_BCSD = 0
       N_OCDP = 0 ; N_OCSV = 0 ; N_OCWT = 0 ; N_OCSD = 0
       N_SUDP = 0 ; N_SUSV = 0 ; N_SUWT = 0 ; N_SUSD = 0
       N_SSDP = 0 ; N_SSSV = 0 ; N_SSWT = 0 ; N_SSSD = 0
       N_DUDP1 = 0 ; N_DUSV1 = 0 ; N_DUWT1 = 0 ; N_DUSD1 = 0
       N_BCDP1 = 0 ; N_BCSV1 = 0 ; N_BCWT1 = 0 ; N_BCSD1 = 0
       N_OCDP1 = 0 ; N_OCSV1 = 0 ; N_OCWT1 = 0 ; N_OCSD1 = 0
       N_SUDP1 = 0 ; N_SUSV1 = 0 ; N_SUWT1 = 0 ; N_SUSD1 = 0
       N_SSDP1 = 0 ; N_SSSV1 = 0 ; N_SSWT1 = 0 ; N_SSSD1 = 0

       call ESMF_VMGetCurrent(VM,                                RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_VMGet       (VM,       mpiCommunicator =comm,   RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_VMGet       (VM, localpet=MYPE, petcount=nDEs,  RC=STATUS)
       VERIFY_(STATUS)

       do K = 1, NUM_AERO_DP
          FIELD_NAME = trim(adjustl(AERO_DP_FIELD_NAME(K)))
          if ( index(FIELD_NAME, 'DUDP') /= 0 ) N_DUDP = N_DUDP + 1
          if ( index(FIELD_NAME, 'DUSV') /= 0 ) N_DUSV = N_DUSV + 1
          if ( index(FIELD_NAME, 'DUWT') /= 0 ) N_DUWT = N_DUWT + 1
          if ( index(FIELD_NAME, 'DUSD') /= 0 ) N_DUSD = N_DUSD + 1
          if ( index(FIELD_NAME, 'BCDP') /= 0 ) N_BCDP = N_BCDP + 1
          if ( index(FIELD_NAME, 'BCSV') /= 0 ) N_BCSV = N_BCSV + 1
          if ( index(FIELD_NAME, 'BCWT') /= 0 ) N_BCWT = N_BCWT + 1
          if ( index(FIELD_NAME, 'BCSD') /= 0 ) N_BCSD = N_BCSD + 1
          if ( index(FIELD_NAME, 'OCDP') /= 0 ) N_OCDP = N_OCDP + 1
          if ( index(FIELD_NAME, 'OCSV') /= 0 ) N_OCSV = N_OCSV + 1
          if ( index(FIELD_NAME, 'OCWT') /= 0 ) N_OCWT = N_OCWT + 1
          if ( index(FIELD_NAME, 'OCSD') /= 0 ) N_OCSD = N_OCSD + 1
          if ( index(FIELD_NAME, 'SUDP') /= 0 ) N_SUDP = N_SUDP + 1
          if ( index(FIELD_NAME, 'SUSV') /= 0 ) N_SUSV = N_SUSV + 1
          if ( index(FIELD_NAME, 'SUWT') /= 0 ) N_SUWT = N_SUWT + 1
          if ( index(FIELD_NAME, 'SUSD') /= 0 ) N_SUSD = N_SUSD + 1
          if ( index(FIELD_NAME, 'SSDP') /= 0 ) N_SSDP = N_SSDP + 1
          if ( index(FIELD_NAME, 'SSSV') /= 0 ) N_SSSV = N_SSSV + 1
          if ( index(FIELD_NAME, 'SSWT') /= 0 ) N_SSWT = N_SSWT + 1
          if ( index(FIELD_NAME, 'SSSD') /= 0 ) N_SSSD = N_SSSD + 1

          if ( index(FIELD_NAME, 'DUDP001') /= 0 ) N_DUDP1 = k
          if ( index(FIELD_NAME, 'DUSV001') /= 0 ) N_DUSV1 = k
          if ( index(FIELD_NAME, 'DUWT001') /= 0 ) N_DUWT1 = k
          if ( index(FIELD_NAME, 'DUSD001') /= 0 ) N_DUSD1 = k
          if ( index(FIELD_NAME, 'BCDP001') /= 0 ) N_BCDP1 = k
          if ( index(FIELD_NAME, 'BCSV001') /= 0 ) N_BCSV1 = k
          if ( index(FIELD_NAME, 'BCWT001') /= 0 ) N_BCWT1 = k
          if ( index(FIELD_NAME, 'BCSD001') /= 0 ) N_BCSD1 = k
          if ( index(FIELD_NAME, 'OCDP001') /= 0 ) N_OCDP1 = k
          if ( index(FIELD_NAME, 'OCSV001') /= 0 ) N_OCSV1 = k
          if ( index(FIELD_NAME, 'OCWT001') /= 0 ) N_OCWT1 = k
          if ( index(FIELD_NAME, 'OCSD001') /= 0 ) N_OCSD1 = k
          if ( index(FIELD_NAME, 'SUDP003') /= 0 ) N_SUDP1 = k
          if ( index(FIELD_NAME, 'SUSV003') /= 0 ) N_SUSV1 = k
          if ( index(FIELD_NAME, 'SUWT003') /= 0 ) N_SUWT1 = k
          if ( index(FIELD_NAME, 'SUSD003') /= 0 ) N_SUSD1 = k
          if ( index(FIELD_NAME, 'SSDP001') /= 0 ) N_SSDP1 = k
          if ( index(FIELD_NAME, 'SSSV001') /= 0 ) N_SSSV1 = k
          if ( index(FIELD_NAME, 'SSWT001') /= 0 ) N_SSWT1 = k
          if ( index(FIELD_NAME, 'SSSD001') /= 0 ) N_SSSD1 = k

       end do

       if (N_DUDP /= NUM_DUDP) STOP 'NUM_DUDP mismatch with AERO_DP'
       if (N_DUSV /= NUM_DUSV) STOP 'NUM_DUSV mismatch with AERO_DP'
       if (N_DUWT /= NUM_DUWT) STOP 'NUM_DUWT mismatch with AERO_DP'
       if (N_DUSD /= NUM_DUSD) STOP 'NUM_DUSD mismatch with AERO_DP'
       if (N_BCDP /= NUM_BCDP) STOP 'NUM_BCDP mismatch with AERO_DP'
       if (N_BCSV /= NUM_BCSV) STOP 'NUM_BCSV mismatch with AERO_DP'
       if (N_BCWT /= NUM_BCWT) STOP 'NUM_BCWT mismatch with AERO_DP'
       if (N_BCSD /= NUM_BCSD) STOP 'NUM_BCSD mismatch with AERO_DP'
       if (N_OCDP /= NUM_OCDP) STOP 'NUM_OCDP mismatch with AERO_DP'
       if (N_OCSV /= NUM_OCSV) STOP 'NUM_OCSV mismatch with AERO_DP'
       if (N_OCWT /= NUM_OCWT) STOP 'NUM_OCWT mismatch with AERO_DP'
       if (N_OCSD /= NUM_OCSD) STOP 'NUM_OCSD mismatch with AERO_DP'
       if (N_SUDP /= NUM_SUDP) STOP 'NUM_SUDP mismatch with AERO_DP'
       if (N_SUSV /= NUM_SUSV) STOP 'NUM_SUSV mismatch with AERO_DP'
       if (N_SUWT /= NUM_SUWT) STOP 'NUM_SUWT mismatch with AERO_DP'
       if (N_SUSD /= NUM_SUSD) STOP 'NUM_SUSD mismatch with AERO_DP'
       if (N_SSDP /= NUM_SSDP) STOP 'NUM_SSDP mismatch with AERO_DP'
       if (N_SSSV /= NUM_SSSV) STOP 'NUM_SSSV mismatch with AERO_DP'
       if (N_SSWT /= NUM_SSWT) STOP 'NUM_SSWT mismatch with AERO_DP'
       if (N_SSSD /= NUM_SSSD) STOP 'NUM_SSSD mismatch with AERO_DP'

       if ( NUM_DUDP /= 0 ) then
           K = N_DUDP1
           allocate( DUDP(IM,JM,NUM_DUDP), STAT=STATUS ); VERIFY_(STATUS)
           DUDP(:,:,1:NUM_DUDP) = AERO_DP(:,:,K:K+NUM_DUDP-1)
       endif

       if ( NUM_DUSV /= 0 ) then
           K = N_DUSV1
           allocate( DUSV(IM,JM,NUM_DUSV), STAT=STATUS ); VERIFY_(STATUS)
           DUSV(:,:,1:NUM_DUSV) = AERO_DP(:,:,K:K+NUM_DUSV-1)
       endif

       if ( NUM_DUWT /= 0 ) then
           K = N_DUWT1
           allocate( DUWT(IM,JM,NUM_DUWT), STAT=STATUS ); VERIFY_(STATUS)
           DUWT(:,:,1:NUM_DUWT) = AERO_DP(:,:,K:K+NUM_DUWT-1)
       endif

       if ( NUM_DUSD /= 0 ) then
           K = N_DUSD1
           allocate( DUSD(IM,JM,NUM_DUSD), STAT=STATUS ); VERIFY_(STATUS)
           DUSD(:,:,1:NUM_DUSD) = AERO_DP(:,:,K:K+NUM_DUSD-1)
       endif

   !---------- Black Carbon ----------
       if ( NUM_BCDP /= 0 ) then
           K = N_BCDP1
           allocate( BCDP(IM,JM,NUM_BCDP), STAT=STATUS ); VERIFY_(STATUS)
           BCDP(:,:,1:NUM_BCDP) = AERO_DP(:,:,K:K+NUM_BCDP-1)
       endif

       if ( NUM_BCSV /= 0 ) then
           K = N_BCSV1
           allocate( BCSV(IM,JM,NUM_BCSV), STAT=STATUS ); VERIFY_(STATUS)
           BCSV(:,:,1:NUM_BCSV) = AERO_DP(:,:,K:K+NUM_BCSV-1)
       endif

       if ( NUM_BCWT /= 0 ) then
           K = N_BCWT1
           allocate( BCWT(IM,JM,NUM_BCWT), STAT=STATUS ); VERIFY_(STATUS)
           BCWT(:,:,1:NUM_BCWT) = AERO_DP(:,:,K:K+NUM_BCWT-1)
       endif

       if ( NUM_BCSD /= 0 ) then
           K = N_BCSD1
           allocate( BCSD(IM,JM,NUM_BCSD), STAT=STATUS ); VERIFY_(STATUS)
           BCSD(:,:,1:NUM_BCSD) = AERO_DP(:,:,K:K+NUM_BCSD-1)
       endif

   !---------- Organic Carbon ----------
       if ( NUM_OCDP /= 0 ) then
           K = N_OCDP1
           allocate( OCDP(IM,JM,NUM_OCDP), STAT=STATUS ); VERIFY_(STATUS)
           OCDP(:,:,1:NUM_OCDP) = AERO_DP(:,:,K:K+NUM_OCDP-1)
       endif

       if ( NUM_OCSV /= 0 ) then
           K = N_OCSV1
           allocate( OCSV(IM,JM,NUM_OCSV), STAT=STATUS ); VERIFY_(STATUS)
           OCSV(:,:,1:NUM_OCSV) = AERO_DP(:,:,K:K+NUM_OCSV-1)
       endif

       if ( NUM_OCWT /= 0 ) then
           K = N_OCWT1
           allocate( OCWT(IM,JM,NUM_OCWT), STAT=STATUS ); VERIFY_(STATUS)
           OCWT(:,:,1:NUM_OCWT) = AERO_DP(:,:,K:K+NUM_OCWT-1)
       endif

       if ( NUM_OCSD /= 0 ) then
           K = N_OCSD1
           allocate( OCSD(IM,JM,NUM_OCSD), STAT=STATUS ); VERIFY_(STATUS)
           OCSD(:,:,1:NUM_OCSD) = AERO_DP(:,:,K:K+NUM_OCSD-1)
       endif

   !---------- Sulfate (aerosol component only; SO4) ----------
       if ( NUM_SUDP /= 0 ) then
           K = N_SUDP1
           allocate( SUDP(IM,JM,NUM_SUDP), STAT=STATUS ); VERIFY_(STATUS)
           SUDP(:,:,1:NUM_SUDP) = AERO_DP(:,:,K:K+NUM_SUDP-1)
       endif

       if ( NUM_SUSV /= 0 ) then
           K = N_SUSV1
           allocate( SUSV(IM,JM,NUM_SUSV), STAT=STATUS ); VERIFY_(STATUS)
           SUSV(:,:,1:NUM_SUSV) = AERO_DP(:,:,K:K+NUM_SUSV-1)
       endif

       if ( NUM_SUWT /= 0 ) then
           K = N_SUWT1
           allocate( SUWT(IM,JM,NUM_SUWT), STAT=STATUS ); VERIFY_(STATUS)
           SUWT(:,:,1:NUM_SUWT) = AERO_DP(:,:,K:K+NUM_SUWT-1)
       endif

       if ( NUM_SUSD /= 0 ) then
           K = N_SUSD1
           allocate( SUSD(IM,JM,NUM_SUSD), STAT=STATUS ); VERIFY_(STATUS)
           SUSD(:,:,1:NUM_SUSD) = AERO_DP(:,:,K:K+NUM_SUSD-1)
       endif

   !---------- Sea Salt ----------
       if ( NUM_SSDP /= 0 ) then
           K = N_SSDP1
           allocate( SSDP(IM,JM,NUM_SSDP), STAT=STATUS ); VERIFY_(STATUS)
           SSDP(:,:,1:NUM_SSDP) = AERO_DP(:,:,K:K+NUM_SSDP-1)
       endif

       if ( NUM_SSSV /= 0 ) then
           K = N_SSSV1
           allocate( SSSV(IM,JM,NUM_SSSV), STAT=STATUS ); VERIFY_(STATUS)
           SSSV(:,:,1:NUM_SSSV) = AERO_DP(:,:,K:K+NUM_SSSV-1)
       endif

       if ( NUM_SSWT /= 0 ) then
           K = N_SSWT1
           allocate( SSWT(IM,JM,NUM_SSWT), STAT=STATUS ); VERIFY_(STATUS)
           SSWT(:,:,1:NUM_SSWT) = AERO_DP(:,:,K:K+NUM_SSWT-1)
       endif

       if ( NUM_SSSD /= 0 ) then
           K = N_SSSD1
           allocate( SSSD(IM,JM,NUM_SSSD), STAT=STATUS ); VERIFY_(STATUS)
           SSSD(:,:,1:NUM_SSSD) = AERO_DP(:,:,K:K+NUM_SSSD-1)
       endif

       deallocate(AERO_DP_FIELD_NAME)
       deallocate(AERO_DP)

    end if

! Read in precip data. This is used in 'coupled' replay
!------------------------------------------------------

    call MAPL_GetResource(MAPL,PRECIP_FILE,LABEL="PRECIP_FILE:",default="null", RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_ClockGet(CLOCK, currTime=CurrentTime, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_TimeGet (currentTime,               &
                       YY=YEAR, MM=MONTH, DD=DAY, &
                       H=HR,    M=MN,     S=SE,   &
                                        RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_TimeSet (currentTime,               &
                       YY=YEAR, MM=MONTH, DD=DAY, &
                       H=HR,    M =30,    S = 0,  &
                                        RC=STATUS )
    VERIFY_(STATUS)


! These exports are the rainfalls and total snowfall that
!  the children of surface see. They can be the exports of
!  moist or can be read from a file.
!---------------------------------------------------------

    call MAPL_GetPointer(EXPORT, RCU , 'PCU' , alloc=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RLS , 'PLS' , alloc=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SNO , 'SNO' , alloc=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ICE , 'ICE' , alloc=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FRZR, 'FRZR', alloc=.true., RC=STATUS); VERIFY_(STATUS)

! These are the precips imported from moist
!------------------------------------------

    call MAPL_GetPointer(IMPORT, PCU     , 'PCU'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PLS     , 'PLS'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SNOFL   , 'SNO'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, ICEFL   , 'ICE'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, FRZRFL  , 'FRZR'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TA      , 'TA'      ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TPREC   , 'TPREC'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PRECCU  , 'CN_PRCP' ,  RC=STATUS); VERIFY_(STATUS)

! This is the default behavior, with all surface components seeing uncorrected precip
!------------------------------------------------------------------------------------

    RCU = PCU
    RLS = PLS
    SNO = SNOFL
    ICE = ICEFL
    FRZR= FRZRFL

! This code is used to have the surface components (land, salwater, etc.) see
!  a "corrected" precip according to Rolf et al. This was used in MERRA-2 and
!  would typically be done one during reanalysis or replay. Also, OGCM-coupled
!  replays require special treatment.
!-----------------------------------------------------------------------------

    REPLACE_PRECIP: if(PRECIP_FILE /= "null") then

       bundle = ESMF_FieldBundleCreate (NAME='PRECIP', RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_FieldBundleSet(bundle, GRID=GRID, RC=STATUS)
       VERIFY_(STATUS)

     ! call MAPL_CFIORead( PRECIP_FILE, CurrentTime, Bundle, RC=STATUS)
     ! VERIFY_(STATUS)
       call MAPL_read_bundle( Bundle, PRECIP_FILE, CurrentTime, regrid_method=REGRID_METHOD_CONSERVE, RC=status)
       VERIFY_(STATUS)
       call ESMFL_BundleGetPointerToData(Bundle,'PRECTOT',PTTe, RC=STATUS)
       VERIFY_(STATUS)


! Catchment required convective and large-scale rain and total snowfall,
!  but files were created for convective precip, total precip, and total
!  snowfall. The following scales model precips (including additions of
!  ICE and FRZR) by the total file precip.
!
! Per 05/2019 discussions with Rolf Reichle and Andrea Molod, the
! following treatment is applied:
!   In the case of tiny (< 0.1 mm/day) model precip, corrected precip
!   is parsed by the freezing point into large-scale rain or snow.
!   Future development using a ramp or more sophisticated approach is desired.
!   Note the original correction precip threshold of 1-4 mm/d was
!   deemed too small, and a 273.15 K temperature threshold instead of
!   MAPL_ICE ( = 273.16 K )
!-----------------------------------------------------------------------

       allocate( PCSCALE(IM,JM), stat=STATUS )
       VERIFY_(STATUS)
       allocate( PRECSUM(IM,JM), stat=STATUS )
       VERIFY_(STATUS)

       ! PRECSUM = uncorrected total precip
       ! PTTe    = total precip from file
       
       PRECSUM = RCU+RLS+SNO+ICE  ! do *not* add FRZR, which is liquid not solid and (probably) incl. in RCU+RLS
                                  ! see comment re. FRZR in GEOS_CatchGridComp.F90 by reichle, 6/6/2025
       
       where (PTTe == MAPL_UNDEF)
          RCU = PCU
          RLS = PLS
          SNO = SNOFL
          ICE = ICEFL
          FRZR= FRZRFL
       elsewhere (PTTe <= 0.0)         ! PTTe /= MAPL_UNDEF .AND. PTTe <= 0.0
          RCU = 0.
          RLS = 0.
          SNO = 0.
          ICE = 0.
          FRZR= 0.
       elsewhere (PRECSUM > 1.15741e-6)  ! Above is not true .AND. model precip > 0.1  mm/day. Note: 0.1mm/day = 0.1/86400 = 1.15741e-6
          PCSCALE = PTTe / PRECSUM
          RCU  = PCSCALE*RCU
          RLS  = PCSCALE*RLS
          SNO  = PCSCALE*SNO
          ICE  = PCSCALE*ICE
          FRZR = PCSCALE*FRZR
       elsewhere (TA > MAPL_TICE)    ! Above is not true .AND. model is warmer than freezing
          RCU = 0.
          RLS = PTTe
          SNO = 0.
          ICE = 0.
          FRZR= 0.
       elsewhere (TA <= MAPL_TICE)
          RCU = 0.
          RLS = 0.
          SNO = PTTe
          ICE = 0.
          FRZR= 0.
       endwhere

       where(RLS<0.0)
          RCU = RCU + RLS
          RLS = 0.0
       endwhere

       deallocate(PCSCALE,PRECSUM)
! Destroy the bundle and its fields
!----------------------------------

       Call ESMF_FieldBundleGet(bundle,fieldcount=fieldcount,rc=status)
       VERIFY_(STATUS)

       Allocate(fieldnames(fieldcount))
       Call ESMF_FieldBundleGet(bundle,fieldNameList=fieldnames,rc=status)
       VERIFY_(STATUS)

       Do I = 1,fieldCount
          Call ESMF_FieldBundleGet(bundle,trim(fieldnames(i)),field=bundle_field,rc=status)
          VERIFY_(STATUS)
          Call ESMF_FieldDestroy(bundle_field,noGarbage=.true.,rc=status)
          VERIFY_(STATUS)
       Enddo
       deAllocate(fieldnames)

       call ESMF_FieldBundleDestroy(bundle,noGarbage=.true.,rc=STATUS)
       VERIFY_(STATUS)

!      deallocate(PTTe)

! Apply latitude taper to replace the model-generated precip
!  only at low latitudes,
!-----------------------------------------------------------

       call MAPL_GetResource ( MAPL, USE_PP_TAPER, Label="USE_PP_TAPER:", DEFAULT=0, RC=STATUS)
       VERIFY_(STATUS)

       TAPER_PRECIP: if(USE_PP_TAPER/=0) then
          call MAPL_GetResource ( MAPL, PP_TAPER_LAT_LOW , Label="PP_TAPER_LAT_LOW:" , DEFAULT=50.0, RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_GetResource ( MAPL, PP_TAPER_LAT_HIGH, Label="PP_TAPER_LAT_HIGH:", DEFAULT=60.0, RC=STATUS)
          VERIFY_(STATUS)

          PP_TAPER_LAT_LOW  = PP_TAPER_LAT_LOW *(MAPL_PI/180.)
          PP_TAPER_LAT_HIGH = PP_TAPER_LAT_HIGH*(MAPL_PI/180.)

          _ASSERT(PP_TAPER_LAT_HIGH > PP_TAPER_LAT_LOW,'needs informative message')

          FACT = 1.0/(PP_TAPER_LAT_HIGH-PP_TAPER_LAT_LOW)

          where(abs(lats)>=PP_TAPER_LAT_LOW .and. abs(lats)<=PP_TAPER_LAT_HIGH)
             RCU = (PCU   *(abs(lats)-PP_TAPER_LAT_LOW) + RCU *(PP_TAPER_LAT_HIGH-abs(lats)))*FACT
             RLS = (PLS   *(abs(lats)-PP_TAPER_LAT_LOW) + RLS *(PP_TAPER_LAT_HIGH-abs(lats)))*FACT
             SNO = (SNOFL *(abs(lats)-PP_TAPER_LAT_LOW) + SNO *(PP_TAPER_LAT_HIGH-abs(lats)))*FACT
             ICE = (ICEFL *(abs(lats)-PP_TAPER_LAT_LOW) + ICE *(PP_TAPER_LAT_HIGH-abs(lats)))*FACT
             FRZR= (FRZRFL*(abs(lats)-PP_TAPER_LAT_LOW) + FRZR*(PP_TAPER_LAT_HIGH-abs(lats)))*FACT
          end where

          where(abs(lats)>PP_TAPER_LAT_HIGH)
             RCU = PCU
             RLS = PLS
             SNO = SNOFL
             ICE = ICEFL
             FRZR= FRZRFL
          end where

       end if TAPER_PRECIP

       where(RCU <0.0) RCU  = 0.0
       where(RLS <0.0) RLS  = 0.0
       where(SNO <0.0) SNO  = 0.0
       where(ICE <0.0) ICE  = 0.0
       where(FRZR<0.0) FRZR = 0.0

    end if REPLACE_PRECIP

! Pointers to gridded internals
!------------------------------

    call MAPL_GetPointer(INTERNAL, CM      , 'CM'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CT      , 'CT'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CQ      , 'CQ'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CN      , 'CN'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, TH      , 'THAT'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QH      , 'QHAT'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, UH      , 'UHAT'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, VH      , 'VHAT'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, TS      , 'TS'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QS      , 'QS'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, RHOS    , 'RHOS'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, D0      , 'D0'     ,  RC=STATUS); VERIFY_(STATUS)

! Pointers to gridded exports
!----------------------------

    call MAPL_GetPointer(EXPORT  , LST     , 'LST'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ALBVR   , 'ALBVR'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ALBVF   , 'ALBVF'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ALBNR   , 'ALBNR'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ALBNF   , 'ALBNF'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , EMISS   , 'EMIS'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DELSS   , 'DELSS'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DELUS   , 'DELUS'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DELVS   , 'DELVS'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DELTS   , 'DELTS'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DELQS   , 'DELQS'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DLQLL   , 'DLQLL'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DLQIL   , 'DLQIL'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , LWI     , 'LWI'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TSOIL1  , 'TSOIL1' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TSOIL2  , 'TSOIL2' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TSOIL3  , 'TSOIL3' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TSOIL4  , 'TSOIL4' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TSOIL5  , 'TSOIL5' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TSOIL6  , 'TSOIL6' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ASNOW   , 'ASNOW'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SHSNOW  , 'SHSNOW' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , AVETSNOW, 'AVETSNOW', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TPSNO   , 'TPSNOW' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TPUST   , 'TPUNST' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TPSAT   , 'TPSAT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TPWLT   , 'TPWLT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TPSURF  , 'TPSURF' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , FRSAT   , 'FRSAT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , FRUST   , 'FRUST'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , FRWLT   , 'FRWLT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SNOMAS  , 'SNOMAS' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SNOWDP  , 'SNOWDP' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WET1    , 'WET1'   ,  RC=STATUS); VERIFY_(STATUS)
    ! NOTE: GOCART's dust code expects WET1 to have all the cells with MAPL_UNDEF
    !       (aka not land) to be replaced with 1.0. We want WET1 to have
    !       MAPL_UNDEF over non-land points, so we need a separate export to pass
    !       to GOCART.
    call MAPL_GetPointer(EXPORT  , WET1_FOR_CHEM    , 'WET1_FOR_CHEM'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WET2    , 'WET2'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WET3    , 'WET3'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WCSF    , 'WCSF'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WCRZ    , 'WCRZ'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WCPR    , 'WCPR'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WESNN1  , 'WESNN1' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WESNN2  , 'WESNN2' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WESNN3  , 'WESNN3' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , CAPAC   , 'CAPAC'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TAUXO   , 'TAUX'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TAUYO   , 'TAUY'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , EVAPO   , 'EVAP'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SHO     , 'SH'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , USTAR   , 'USTAR'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TSTAR   , 'TSTAR'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , QSTAR   , 'QSTAR'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , U10N    , 'U10N'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , V10N    , 'V10N'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , UAX     , 'UA'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , VAX     , 'VA'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TA      , 'TA'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , QA      , 'QA'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , HLATN   , 'LHFX'   ,  RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT  , DCOOL   , 'DCOOL'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DWARM   , 'DWARM'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TDROP   , 'TDROP'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , QCOOL   , 'QCOOL'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SWCOOL  , 'SWCOOL' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , USTARW  , 'USTARW' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TBAR    , 'TBAR'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , LCOOL   , 'LCOOL'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , BCOOL   , 'BCOOL'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TDEL    , 'TDEL'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TS_FOUND, 'TS_FOUND', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SS_FOUND, 'SS_FOUND', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , QWARM   , 'QWARM'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SWWARM  , 'SWWARM' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , LANGM   , 'LANGM'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , PHIW    , 'PHIW'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TAUTW   , 'TAUTW'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ZETA_W  , 'ZETA_W' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TWMTF,    'TWMTF'  , RC=STATUS); VERIFY_(STATUS)


! if we are running Louis sfc layer, get these exports from the gridded Louis fluxes
  if (CHOOSEMOSFC.eq.0) then
    call MAPL_GetPointer(EXPORT  , U50M    , 'U50M'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , V50M    , 'V50M'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , U10M    , 'U10M'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , V10M    , 'V10M'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , T10M    , 'T10M'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , Q10M    , 'Q10M'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , U2M     , 'U2M'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , V2M     , 'V2M'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , UU10M   , 'UU10M'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RH2M    , 'RH2M'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , T2M     , 'T2M'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , Q2M     , 'Q2M'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , T2MDEW  , 'T2MDEW' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , T2MWET  , 'T2MWET' ,  RC=STATUS); VERIFY_(STATUS)
  endif


    call MAPL_GetPointer(EXPORT  , HLATWTR   , 'HLATWTR'   ,  alloc=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , HLATICE   , 'HLATICE'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  ,   SHWTR   , 'SHWTR'     ,  alloc=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  ,   SHICE   , 'SHICE'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  ,   TAUXW   , 'TAUXW'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  ,   TAUXI   , 'TAUXI'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  ,   TAUYW   , 'TAUYW'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  ,   TAUYI   , 'TAUYI'     ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , LWNDWTR   , 'LWNDWTR'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SWNDWTR   , 'SWNDWTR'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , LWNDICE   , 'LWNDICE'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SWNDICE   , 'SWNDICE'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SNOWOCN   , 'SNOWOCN'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ICEFOCN   , 'ICEFOCN'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SPTOTOCN  , 'SPTOTOCN'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RAINOCN   , 'RAINOCN'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TSKINW    , 'TSKINW'    ,  alloc=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TSKINICE  , 'TSKINICE'  ,  RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT  , HICE    , 'HICE'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , HSNO    , 'HSNO'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , FRZMLT  , 'FRZMLT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TSKINWCICE, 'TSKINWCICE',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ISTSFC  , 'ISTSFC'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SSKINW  , 'SSKINW'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MELTT   , 'MELTT'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MELTB   , 'MELTB'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MELTS   , 'MELTS'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , MELTL   , 'MELTL'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , FRAZIL  , 'FRAZIL'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , CONGEL  , 'CONGEL'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SNOICE  , 'SNOICE'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DAIDTT  , 'DAIDTT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DVIDTT  , 'DVIDTT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DAIDTD  , 'DAIDTD'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DVIDTD  , 'DVIDTD'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , FBOT    , 'FBOT'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , HFLUX   , 'HFLUX'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WFLUX   ,'WATERFLUX',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SFLUX   , 'SALTFLUX',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , FSWTHRU , 'FSWTHRU' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , FSWABS  , 'FSWABS'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , USTARI  , 'USTARI'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , FHOCN   , 'FHOCN'   ,  RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT  , EVAPOU  , 'EVAPOUT',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SUBLIM  , 'SUBLIM' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SHOU    , 'SHOUT'  ,  alloc=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , HLWUP   , 'HLWUP'  ,  alloc=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , LWNDSRF , 'LWNDSRF',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SWNDSRF , 'SWNDSRF',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RUNOFF  , 'RUNOFF' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DISCHARGE,'DISCHARGE',RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RUNSURF , 'RUNSURF',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , BASEFLOW, 'BASEFLOW', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ACCUM   , 'ACCUM'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SMELT   , 'SMELT'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , EVEG    , 'EVPVEG' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , EINT    , 'EVPINT' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , EICE    , 'EVPICE' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ESNO    , 'EVPSNO' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ESOI    , 'EVPSOI' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WAT10CM , 'WAT10CM',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WATSOI  , 'WATSOI' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ICESOI  , 'ICESOI' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , EVLAND  , 'EVLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , PRLAND  , 'PRLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SNOLAND  , 'SNOLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DRPARLAND  , 'DRPARLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DFPARLAND  , 'DFPARLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , LHSNOW  , 'LHSNOW' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SWNETSNOW  , 'SWNETSNOW' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , LWUPSNOW  , 'LWUPSNOW' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , LWDNSNOW  , 'LWDNSNOW' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TCSORIG  , 'TCSORIG' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TPSN1IN  , 'TPSN1IN' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TPSN1OUT  , 'TPSN1OUT' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , GHSNOW  , 'GHSNOW' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , LHLAND  , 'LHLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SHLAND  , 'SHLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SWLAND  , 'SWLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SWDOWNLAND  , 'SWDOWNLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , LWLAND  , 'LWLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , GHLAND  , 'GHLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , GHTSKIN  , 'GHTSKIN' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SMLAND  , 'SMLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , QINFIL  , 'QINFIL' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TWLAND  , 'TWLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TELAND  , 'TELAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TSLAND  , 'TSLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DWLAND  , 'DWLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DHLAND  , 'DHLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SPLAND  , 'SPLAND' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SPLH    , 'SPLH'   ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SPWATR  , 'SPWATR' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , SPSNOW  , 'SPSNOW' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RDU001  , 'RDU001' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RDU002  , 'RDU002' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RDU003  , 'RDU003' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RDU004  , 'RDU004' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RDU005  , 'RDU005' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RBC001  , 'RBC001' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RBC002  , 'RBC002' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ROC001  , 'ROC001' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , ROC002  , 'ROC002' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RMELTDU001 , 'RMELTDU001',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RMELTDU002 , 'RMELTDU002',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RMELTDU003 , 'RMELTDU003',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RMELTDU004 , 'RMELTDU004',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RMELTDU005 , 'RMELTDU005',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RMELTBC001 , 'RMELTBC001',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RMELTBC002 , 'RMELTBC002',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RMELTOC001 , 'RMELTOC001',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , RMELTOC002 , 'RMELTOC002',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , PEATCLSM_WATERLEVEL, 'PEATCLSM_WATERLEVEL', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , PEATCLSM_FSWCHANGE,  'PEATCLSM_FSWCHANGE',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DZGT1   ,  'DZGT1' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DZGT2   ,  'DZGT2' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DZGT3   ,  'DZGT3' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DZGT4   ,  'DZGT4' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DZGT5   ,  'DZGT5' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DZGT6   ,  'DZGT6' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DZPR    ,  'DZPR'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DZRZ    ,  'DZRZ'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DZSF    ,  'DZSF'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DZTS    ,  'DZTS'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WPWET   ,  'WPWET' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WPEMW   ,  'WPEMW' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , WPMC    ,  'WPMC'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , CDCR2   ,  'CDCR2' ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , POROS   ,  'POROS' ,  RC=STATUS); VERIFY_(STATUS)

    IF(LSM_CHOICE > 1) THEN
       call MAPL_GetPointer(EXPORT  , CNLAI   , 'CNLAI'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNTLAI  , 'CNTLAI' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNSAI   , 'CNSAI'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNTOTC  , 'CNTOTC' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNVEGC  , 'CNVEGC' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNROOT  , 'CNROOT' ,  RC=STATUS); VERIFY_(STATUS)
       if (LSM_CHOICE == 3) then
           call MAPL_GetPointer(EXPORT  , CNFROOTC, 'CNFROOTC' ,RC=STATUS); VERIFY_(STATUS)
       endif
       call MAPL_GetPointer(EXPORT  , CNNPP   , 'CNNPP'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNGPP   , 'CNGPP'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNSR    , 'CNSR'   ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNNEE   , 'CNNEE'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNXSMR  , 'CNXSMR' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNADD   , 'CNADD'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNLOSS  , 'CNLOSS' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNBURN  , 'CNBURN' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , PARABS  , 'PARABS' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , PARINC  , 'PARINC' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , SCSAT   , 'SCSAT'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , SCUNS   , 'SCUNS'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , BTRANT  , 'BTRANT' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , SIF     , 'SIF'    ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , CNFSEL  , 'CNFSEL' ,  RC=STATUS); VERIFY_(STATUS)
    END IF

    if (DO_WAVES /= 0) then
       call MAPL_GetPointer(EXPORT  , PS_     , 'PS',  alloc=.true., RC=STATUS); VERIFY_(STATUS)
       PS_ = PS
    end if

    if (DO_FIRE_DANGER) then
       call MAPL_GetPointer(EXPORT  , FFMC        , 'FFMC'       ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , GFMC        , 'GFMC'       ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , DMC         , 'DMC'        ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , DC          , 'DC'         ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , FWI         , 'FWI'        ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , BUI         , 'BUI'        ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , ISI         , 'ISI'        ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , DSR         , 'DSR'        ,  RC=STATUS); VERIFY_(STATUS)

       call MAPL_GetPointer(EXPORT  , FFMC_DAILY  , 'FFMC_DAILY' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , DMC_DAILY   , 'DMC_DAILY'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , DC_DAILY    , 'DC_DAILY'   ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , FWI_DAILY   , 'FWI_DAILY'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , BUI_DAILY   , 'BUI_DAILY'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , ISI_DAILY   , 'ISI_DAILY'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , DSR_DAILY   , 'DSR_DAILY'  ,  RC=STATUS); VERIFY_(STATUS)

       call MAPL_GetPointer(EXPORT  , FFMC_DAILY_ , 'FFMC_DAILY_',  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , DMC_DAILY_  , 'DMC_DAILY_' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , DC_DAILY_   , 'DC_DAILY_'  ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , FWI_DAILY_  , 'FWI_DAILY_' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , BUI_DAILY_  , 'BUI_DAILY_' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , ISI_DAILY_  , 'ISI_DAILY_' ,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT  , DSR_DAILY_  , 'DSR_DAILY_' ,  RC=STATUS); VERIFY_(STATUS)

       call MAPL_GetPointer(EXPORT  , VPD         , 'VPD'        ,  RC=STATUS); VERIFY_(STATUS)
    end if

! Force allocation for ice fraction for lwi mask

    call MAPL_GetPointer(EXPORT  , FRI     , 'FRACI'  ,  alloc=associated(LWI) , rC=STATUS)
    VERIFY_(STATUS)

! Not sure about why alloc of FRI depends on LWI, but copy the logic anyway 
    call MAPL_GetPointer(EXPORT  , OFRI    , 'OFRACI' ,  alloc=associated(LWI) , rC=STATUS)
    VERIFY_(STATUS)

!    FRI =  max(min(FRI,1.0),0.0)

! RiverRouting: force allocations of RUNOFF from continental components,
!   and make sure RoutingFile was specified.
!-----------------------------------------------------------------------

    if (associated(DISCHARGE)) then
       _ASSERT(associated(SURF_INTERNAL_STATE%RoutingType),'needs informative message')
!ALT       call MAPL_GetPointer(EXPORT, RUNOFF, 'RUNOFF', alloc=.true.,  RC=STATUS)
!ALT       VERIFY_(STATUS)
    end if

! Allocate some work arrays in grid and tile space
!-------------------------------------------------

    NT = size(TILETYPE)

    allocate(   Z0 (IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   TMP(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   TTM(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   QTM(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   FAC(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   TAU(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   DTS(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   DQS(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)

    allocate( DRPAR(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate( DFPAR(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate( DRNIR(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate( DFNIR(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate( DRUVR(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate( DFUVR(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate( ZTH  (IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    allocate( SLR  (IM,JM), STAT=STATUS)
    VERIFY_(STATUS)

    allocate(  DTSTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  DQSTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)

! Get the insolation and zenith angle on grid and tiles
!------------------------------------------------------

    call ESMF_ClockGet(CLOCK,     TIMESTEP=DELT, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_SunGetInsolation(LONS, LATS,      &
              ORBIT, ZTH, SLR, &
              INTV  = DELT,    &
              CLOCK = CLOCK,   &
              RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, SC, 'SOLAR_CONSTANT:', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, SolCycFileName, "SOLAR_CYCLE_FILE_NAME:", DEFAULT='/dev/null', RC=STATUS)
    VERIFY_(STATUS)

    if(SolCycFileName /= '/dev/null') THEN

       call MAPL_GetResource( MAPL, USE_NRLSSI2, "USE_NRLSSI2:", DEFAULT=.TRUE., RC=STATUS)
       VERIFY_(STATUS)

       if (USE_NRLSSI2) then
          call MAPL_GetResource( MAPL, PersistSolar, "PERSIST_SOLAR:", DEFAULT=.TRUE., RC=STATUS)
          VERIFY_(STATUS)

          call MAPL_SunGetSolarConstant(CLOCK,trim(SolCycFileName),SC,MG,SB,PersistSolar=PersistSolar,rc=STATUS)
          VERIFY_(STATUS)
       else
          call MAPL_SunGetSolarConstant(CLOCK,trim(SolCycFileName),SC,RC=STATUS)
          VERIFY_(STATUS)
       endif
    else if(SC<0.0) then
       call MAPL_SunGetSolarConstant(CURRENTTIME,SC,RC=STATUS)
       VERIFY_(STATUS)
    end if

    DRPAR = DRPARN * SLR * SC
    DFPAR = DFPARN * SLR * SC
    DRNIR = DRNIRN * SLR * SC
    DFNIR = DFNIRN * SLR * SC
    DRUVR = DRUVRN * SLR * SC
    DFUVR = DFUVRN * SLR * SC

    allocate(  SLITILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  ZTHTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)

    call MAPL_SunGetInsolation(     &
         tileLONS, tileLATS, ORBIT, &
         ZTHTILE, SLITILE,          &
         INTV  = DELT,              &
         CLOCK = CLOCK,             &
                          RC=STATUS )
    VERIFY_(STATUS)

    ZTHTILE = max(0.0,ZTHTILE)

! We need atmospheric version of the run1 outputs put back on tiles
!------------------------------------------------------------------

    allocate(   TSTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   QSTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   THTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   QHTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   UHTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   VHTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   CTTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   CQTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   CMTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)

    call MAPL_LocStreamTransform( LOCSTREAM, THTILE, TH, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, QHTILE, QH, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, UHTILE, UH, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, VHTILE, VH, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, CTTILE, CT, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, CQTILE, CQ, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, CMTILE, CM, RC=STATUS); VERIFY_(STATUS)

! Allocate tile versions of imports
!----------------------------------

    allocate(   PSTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   DZTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(   UUTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  EVPTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  SHFTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  DEVTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  DSHTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  TMPTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  PCUTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  PLSTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(SNOFLTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(ICEFLTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(FRZRFLTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate( TAUXTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate( TAUYTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  DRUTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  DFUTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  DRPTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  DFPTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  DRNTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  DFNTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  LWBTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  ALWTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  BLWTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(DUDPTILE(NT,NUM_DUDP), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(DUSVTILE(NT,NUM_DUSV), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(DUWTTILE(NT,NUM_DUWT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(DUSDTILE(NT,NUM_DUSD), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(BCDPTILE(NT,NUM_BCDP), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(BCSVTILE(NT,NUM_BCSV), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(BCWTTILE(NT,NUM_BCWT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(BCSDTILE(NT,NUM_BCSD), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(OCDPTILE(NT,NUM_OCDP), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(OCSVTILE(NT,NUM_OCSV), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(OCWTTILE(NT,NUM_OCWT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(OCSDTILE(NT,NUM_OCSD), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(SUDPTILE(NT,NUM_SUDP), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(SUSVTILE(NT,NUM_SUSV), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(SUWTTILE(NT,NUM_SUWT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(SUSDTILE(NT,NUM_SUSD), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(SSDPTILE(NT,NUM_SSDP), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(SSSVTILE(NT,NUM_SSSV), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(SSWTTILE(NT,NUM_SSWT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(SSSDTILE(NT,NUM_SSSD), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(  DTSDTTILE(NT), STAT=STATUS)
    VERIFY_(STATUS)

! Transform imports to the tiles
!-------------------------------

    call MAPL_LocStreamTransform( LOCSTREAM, PSTILE   , PS,      RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, DZTILE   , DZ,      RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, UUTILE   , UU,      RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, EVPTILE  , EVAP,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, SHFTILE  , SH,      RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, DEVTILE  , DEVAP,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, DSHTILE  , DSH,     RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, PCUTILE  , RCU,     RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, PLSTILE  , RLS,     RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, SNOFLTILE, SNO,     RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, ICEFLTILE, ICE,     RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, FRZRFLTILE, FRZR,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, TAUXTILE , TAUX,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, TAUYTILE , TAUY,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, DRUTILE  , DRUVR,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, DFUTILE  , DFUVR,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, DRPTILE  , DRPAR,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, DFPTILE  , DFPAR,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, DRNTILE  , DRNIR,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, DFNTILE  , DFNIR,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, LWBTILE  , LWDNSRF, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, ALWTILE  , ALW,     RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, BLWTILE  , BLW,     RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, DTSDTTILE, DTSDT,   RC=STATUS); VERIFY_(STATUS)

    if (DO_GOSWIM) then
       do K = 1, NUM_DUDP
          call MAPL_LocStreamTransform(LOCSTREAM, DUDPTILE(:,K), DUDP(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_DUSV
          call MAPL_LocStreamTransform(LOCSTREAM, DUSVTILE(:,K), DUSV(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_DUWT
          call MAPL_LocStreamTransform(LOCSTREAM, DUWTTILE(:,K), DUWT(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_DUSD
          call MAPL_LocStreamTransform(LOCSTREAM, DUSDTILE(:,K), DUSD(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_BCDP
          call MAPL_LocStreamTransform(LOCSTREAM, BCDPTILE(:,K), BCDP(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_BCSV
          call MAPL_LocStreamTransform(LOCSTREAM, BCSVTILE(:,K), BCSV(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_BCWT
          call MAPL_LocStreamTransform(LOCSTREAM, BCWTTILE(:,K), BCWT(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_BCSD
          call MAPL_LocStreamTransform(LOCSTREAM, BCSDTILE(:,K), BCSD(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_OCDP
          call MAPL_LocStreamTransform(LOCSTREAM, OCDPTILE(:,K), OCDP(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_OCSV
          call MAPL_LocStreamTransform(LOCSTREAM, OCSVTILE(:,K), OCSV(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_OCWT
          call MAPL_LocStreamTransform(LOCSTREAM, OCWTTILE(:,K), OCWT(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_OCSD
          call MAPL_LocStreamTransform(LOCSTREAM, OCSDTILE(:,K), OCSD(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_SUDP
          call MAPL_LocStreamTransform(LOCSTREAM, SUDPTILE(:,K), SUDP(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_SUSV
          call MAPL_LocStreamTransform(LOCSTREAM, SUSVTILE(:,K), SUSV(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_SUWT
          call MAPL_LocStreamTransform(LOCSTREAM, SUWTTILE(:,K), SUWT(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_SUSD
          call MAPL_LocStreamTransform(LOCSTREAM, SUSDTILE(:,K), SUSD(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_SSDP
          call MAPL_LocStreamTransform(LOCSTREAM, SSDPTILE(:,K), SSDP(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_SSSV
          call MAPL_LocStreamTransform(LOCSTREAM, SSSVTILE(:,K), SSSV(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_SSWT
          call MAPL_LocStreamTransform(LOCSTREAM, SSWTTILE(:,K), SSWT(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
       do K = 1, NUM_SSSD
          call MAPL_LocStreamTransform(LOCSTREAM, SSSDTILE(:,K), SSSD(:,:,K), RC=STATUS); VERIFY_(STATUS)
       end do
    end if

    ! option to interpolate effective wind vectors such that
    ! stresses computed in all children will have smooth curl/divergence
    call MAPL_GetResource(MAPL, iUseInterp, 'INTERPOLATE_ATMTAU:', &
         default=0, RC=STATUS )
    VERIFY_(STATUS)
    useInterp = (iUseInterp /= 0)

    allocate( UUA(IM,JM), STAT=STATUS ); VERIFY_(STATUS)
    allocate( VVA(IM,JM), STAT=STATUS ); VERIFY_(STATUS)

! The import taus come from turbulence and are stresses on the atmosphere
    UUA = (-TAUX/CM + UH)  ! is division safe here??
    VVA = (-TAUY/CM + VH)  ! is division safe here??

    call MKTILE(UUA , UUATILE,  NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(VVA , VVATILE,  NT,RC=STATUS); VERIFY_(STATUS)

    call MAPL_LocStreamTransform( LOCSTREAM, UUATILE, UUA, INTERP=useInterp, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM, VVATILE, VVA, INTERP=useInterp, RC=STATUS); VERIFY_(STATUS)

    deallocate(UUA)
    deallocate(VVA)

! The import taus come from turbulence and are stresses on the atmosphere

    TAUXTILE = -TAUXTILE
    TAUYTILE = -TAUYTILE
    DSHTILE  = DSHTILE * MAPL_CP ! ???

! If a grid export is required, allocate its tile version
!--------------------------------------------------------

    call MKTILE(LST     ,LSTTILE     ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ALBVR   ,ALBVRTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ALBVF   ,ALBVFTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ALBNR   ,ALBNRTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ALBNF   ,ALBNFTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(EMISS   ,EMISSTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(FRI     ,FRTILE      ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(OFRI    ,OFRTILE     ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TSOIL1  ,TSOIL1TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TSOIL2  ,TSOIL2TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TSOIL3  ,TSOIL3TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TSOIL4  ,TSOIL4TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TSOIL5  ,TSOIL5TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TSOIL6  ,TSOIL6TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ASNOW   ,ASNOWTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SHSNOW  ,SHSNOWTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(AVETSNOW,AVETSNOWTILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TPSNO   ,TPSNOTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TPUST   ,TPUSTTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TPSAT   ,TPSATTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TPWLT   ,TPWLTTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TPSURF  ,TPSURFTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(FRSAT   ,FRSATTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(FRUST   ,FRUSTTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(FRWLT   ,FRWLTTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SNOMAS  ,SNOWTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SNOWDP  ,SNODTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WET1    ,WET1TILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WET2    ,WET2TILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WET3    ,WET3TILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WCSF    ,WCSFTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WCRZ    ,WCRZTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WCPR    ,WCPRTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WESNN1  ,WESNN1TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WESNN2  ,WESNN2TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WESNN3  ,WESNN3TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(CAPAC   ,CAPACTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(HLATN   ,HLATNTILE   ,NT,RC=STATUS); VERIFY_(STATUS)

    call MKTILE(HLATWTR   ,HLATWTRTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(HLATICE   ,HLATICETILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(  SHWTR   ,  SHWTRTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(  SHICE   ,  SHICETILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(  TAUXW   ,  TAUXWTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(  TAUXI   ,  TAUXITILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(  TAUYW   ,  TAUYWTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(  TAUYI   ,  TAUYITILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(LWNDWTR   ,LWNDWTRTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SWNDWTR   ,SWNDWTRTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(LWNDICE   ,LWNDICETILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SWNDICE   ,SWNDICETILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SNOWOCN   ,SNOWOCNTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ICEFOCN   ,ICEFOCNTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SPTOTOCN  ,SPTOTOCNTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RAINOCN   ,RAINOCNTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TSKINW,   TSKINWTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TSKINICE, TSKINICETILE ,NT,RC=STATUS); VERIFY_(STATUS)

    call MKTILE(DCOOL,   DCOOL_TILE  ,   NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DWARM,   DWARM_TILE  ,   NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TDROP,   TDROP_TILE  ,   NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(QCOOL,   QCOOL_TILE  ,   NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SWCOOL,  SWCOOL_TILE ,   NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(USTARW,  USTARW_TILE  ,  NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TBAR,    TBAR_TILE   ,   NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(LCOOL,   LCOOL_TILE  ,   NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(BCOOL,   BCOOL_TILE  ,   NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TDEL,    TDEL_TILE     , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TS_FOUND,TS_FOUND_TILE , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SS_FOUND,SS_FOUND_TILE , NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(QWARM,   QWARM_TILE  ,   NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SWWARM,  SWWARM_TILE ,   NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(LANGM,   LANGM_TILE ,    NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(PHIW,    PHIW_TILE ,     NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TAUTW,   TAUTW_TILE ,    NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ZETA_W,  ZETA_W_TILE ,   NT, RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TWMTF,   TWMTF_TILE ,    NT, RC=STATUS); VERIFY_(STATUS)

    call MKTILE(HICE   ,      HICETILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(HSNO   ,      HSNOTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(FRZMLT ,      FRZMLTTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TSKINWCICE,   TSKINWCICETILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ISTSFC ,      ISTSFCTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SSKINW ,      SSKINWTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MELTT  ,       MELTTTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MELTB  ,       MELTBTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MELTS  ,       MELTSTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(MELTL  ,       MELTLTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(FRAZIL ,       FRAZILTILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(CONGEL ,       CONGELTILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SNOICE ,       SNOICETILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DAIDTT ,       DAIDTTTILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DVIDTT ,       DVIDTTTILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DAIDTD ,       DAIDTDTILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DVIDTD ,       DVIDTDTILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(FBOT   ,       FBOTTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(HFLUX  ,       HFLUXTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WFLUX  ,       WFLUXTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SFLUX  ,       SFLUXTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(FSWTHRU,       FSWTHRUTILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(FSWABS ,       FSWABSTILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(USTARI ,       USTARITILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(FHOCN  ,       FHOCNTILE ,NT,RC=STATUS); VERIFY_(STATUS)

    call MKTILE(EVAPOU  ,EVAPOUTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SUBLIM  ,SUBLIMTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SHOU    ,SHOUTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(HLWUP   ,HLWUPTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(LWNDSRF ,LWNDSRFTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SWNDSRF ,SWNDSRFTILE ,NT,RC=STATUS); VERIFY_(STATUS)

    allocateRunoff = .false.
    if (associated(RUNOFF)) allocateRunoff = .true.

    if (associated(SURF_INTERNAL_STATE%RoutingType) .or. DO_DATA_ATM4OCN) then ! routing file exists or we run DataAtm
       allocate(DISCHARGETILE(NT),stat=STATUS); VERIFY_(STATUS)
       DISCHARGETILE=MAPL_Undef
       allocateRunoff = .true.
    end if
    if (allocateRunoff) then
       allocate(RUNOFFTILE(NT),stat=STATUS); VERIFY_(STATUS)
       RUNOFFTILE = 0.0
    end if

    call MKTILE(RUNSURF ,RUNSURFTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(BASEFLOW,BASEFLOWTILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ACCUM   ,ACCUMTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SMELT   ,SMELTTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(EVEG    ,EVEGTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(EINT    ,EINTTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(EICE    ,EICETILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ESNO    ,ESNOTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ESOI    ,ESOITILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WAT10CM ,WAT10CMTILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WATSOI  ,WATSOITILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ICESOI  ,ICESOITILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(EVLAND  ,EVLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(PRLAND  ,PRLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SNOLAND  ,SNOLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DRPARLAND  ,DRPARLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DFPARLAND  ,DFPARLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(LHSNOW  ,LHSNOWTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SWNETSNOW  ,SWNETSNOWTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(LWUPSNOW  ,LWUPSNOWTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(LWDNSNOW  ,LWDNSNOWTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TCSORIG  ,TCSORIGTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TPSN1IN  ,TPSN1INTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TPSN1OUT  ,TPSN1OUTTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(GHSNOW  ,GHSNOWTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(LHLAND  ,LHLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SHLAND  ,SHLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SWLAND  ,SWLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SWDOWNLAND  ,SWDOWNLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(LWLAND  ,LWLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(GHLAND  ,GHLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(GHTSKIN  ,GHTSKINTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SMLAND  ,SMLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(QINFIL  ,QINFILTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TWLAND  ,TWLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TELAND  ,TELANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(TSLAND  ,TSLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DWLAND  ,DWLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DHLAND  ,DHLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SPLAND  ,SPLANDTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SPLH    ,SPLHTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SPWATR  ,SPWATRTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(SPSNOW  ,SPSNOWTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RDU001  ,RDU001TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RDU002  ,RDU002TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RDU003  ,RDU003TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RDU004  ,RDU004TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RDU005  ,RDU005TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RBC001  ,RBC001TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RBC002  ,RBC002TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ROC001  ,ROC001TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(ROC002  ,ROC002TILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RMELTDU001 ,RMELTDU001TILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RMELTDU002 ,RMELTDU002TILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RMELTDU003 ,RMELTDU003TILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RMELTDU004 ,RMELTDU004TILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RMELTDU005 ,RMELTDU005TILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RMELTBC001 ,RMELTBC001TILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RMELTBC002 ,RMELTBC002TILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RMELTOC001 ,RMELTOC001TILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(RMELTOC002 ,RMELTOC002TILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(PEATCLSM_WATERLEVEL,PEATCLSM_WATERLEVELTILE,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(PEATCLSM_FSWCHANGE ,PEATCLSM_FSWCHANGETILE ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DZGT1   ,DZGT1TILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DZGT2   ,DZGT2TILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DZGT3   ,DZGT3TILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DZGT4   ,DZGT4TILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DZGT5   ,DZGT5TILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DZGT6   ,DZGT6TILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DZPR    ,DZPRTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DZRZ    ,DZRZTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DZSF    ,DZSFTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(DZTS    ,DZTSTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WPWET   ,WPWETTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WPEMW   ,WPEMWTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(WPMC    ,WPMCTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(CDCR2   ,CDCR2TILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    call MKTILE(POROS   ,POROSTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
    
    IF (LSM_CHOICE > 1) THEN
       call MKTILE(CNLAI   ,CNLAITILE   ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNTLAI  ,CNTLAITILE  ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNSAI   ,CNSAITILE   ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNTOTC  ,CNTOTCTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNVEGC  ,CNVEGCTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNROOT  ,CNROOTTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
       if (LSM_CHOICE == 3) then
          call MKTILE(CNFROOTC,CNFROOTCTILE  ,NT,RC=STATUS);VERIFY_(STATUS)
       endif
       call MKTILE(CNNPP   ,CNNPPTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNGPP   ,CNGPPTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNSR    ,CNSRTILE    ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNNEE   ,CNNEETILE   ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNXSMR  ,CNXSMRTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNADD   ,CNADDTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNLOSS  ,CNLOSSTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNBURN  ,CNBURNTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(PARABS  ,PARABSTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(PARINC  ,PARINCTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(SCSAT   ,SCSATTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(SCUNS   ,SCUNSTILE   ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(BTRANT  ,BTRANTTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(SIF     ,SIFTILE     ,NT,RC=STATUS); VERIFY_(STATUS)
       call MKTILE(CNFSEL  ,CNFSELTILE  ,NT,RC=STATUS); VERIFY_(STATUS)
    END IF

    if (DO_FIRE_DANGER) then
       call MKTILE(FFMC,        FFMCTILE,          NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(GFMC,        GFMCTILE,          NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(DMC,         DMCTILE,           NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(DC,          DCTILE,            NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(FWI,         FWITILE,           NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(BUI,         BUITILE,           NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(ISI,         ISITILE,           NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(DSR,         DSRTILE,           NT,  RC=STATUS); VERIFY_(STATUS)

       call MKTILE(FFMC_DAILY,  FFMCDAILYTILE,     NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(DMC_DAILY,   DMCDAILYTILE,      NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(DC_DAILY,    DCDAILYTILE,       NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(FWI_DAILY,   FWIDAILYTILE,      NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(BUI_DAILY,   BUIDAILYTILE,      NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(ISI_DAILY,   ISIDAILYTILE,      NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(DSR_DAILY,   DSRDAILYTILE,      NT,  RC=STATUS); VERIFY_(STATUS)

       call MKTILE(FFMC_DAILY_, FFMCDAILYTILE_,    NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(DMC_DAILY_,  DMCDAILYTILE_,     NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(DC_DAILY_,   DCDAILYTILE_,      NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(FWI_DAILY_,  FWIDAILYTILE_,     NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(BUI_DAILY_,  BUIDAILYTILE_,     NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(ISI_DAILY_,  ISIDAILYTILE_,     NT,  RC=STATUS); VERIFY_(STATUS)
       call MKTILE(DSR_DAILY_,  DSRDAILYTILE_,     NT,  RC=STATUS); VERIFY_(STATUS)

       call MKTILE(VPD,         VPDTILE,           NT,  RC=STATUS); VERIFY_(STATUS)
    end if


    FRTILE  = 0.0
    OFRTILE = MAPL_UNDEF

!  Cycle through all continental children (skip ocean),
!   collecting RUNOFFTILE exports.

    do I = 1, NUM_CHILDREN
       if (I == OCEAN) cycle
       call DOTYPE(I,RC=STATUS)
       VERIFY_(STATUS)
    end do

! Create the Discharge for the ocean. This is an import of
!  Saltwater, which simply makes a copy to an export.
!  That export is what goes to GcmGridComp for coupling to the OGCM.
!  Because OGCM gets tile variables under the table from saltwater.
!  we have to go through an unecessary packing and unpacking of the
!  between the globally tiled discharge and the ocean only tiled discharge.
!--------------------------------------------------------------

    if(associated(DISCHARGETILE)) then

       ! Create discharge at exit tiles by routing runoff

       if (DO_DATA_ATM4OCN) then
          call MAPL_GetPointer(IMPORT  , DISCHARGE_IM, 'DISCHARGE',  RC=STATUS); VERIFY_(STATUS)
          call MAPL_LocStreamTransform( LOCSTREAM,  RUNOFFTILE, DISCHARGE_IM, RC=STATUS)
          VERIFY_(STATUS)
          ! it seems redundant to fill both DISCHARGETILE and RUNOFFTILE
          ! but this is done in case we need to output RUNOFF
          ! and not to change the existing code too much
          DISCHARGETILE = RUNOFFTILE

       else
          call RouteRunoff(SURF_INTERNAL_STATE%RoutingType, RUNOFFTILE, DISCHARGETILE, RC=STATUS)
          VERIFY_(STATUS)
       end if

       !-------------------------------------------------------------------------------------
       !  Special treatment for doing ocean-coupled atmospheric replays to an analysis
       !  that used "corrected" precips.

       !  This is done to prevent an ocean mass budget imbalance, since the atmosphere
       !  water budget is in balance only with the "uncorrected" precip, and the land needs
       !  corrected precips.

       !  It involves using the uncorrected precip over ocean tiles, which is done here.

       !  we also need to modify the river discharge that the the ocean sees
       !  to account, on average, for the difference between corrected and uncorrected precip
       !  over the continents.
       !-------------------------------------------------------------------------------------

       call MAPL_GetResource ( MAPL, DischargeAdjustFile, Label="DISCHARGE_ADJUST_FILE:", &
            DEFAULT="null", RC=STATUS)
       VERIFY_(STATUS)

       ! Do not correct the precip over ocean tiles
       if(Precip_File /= "null") then

          call MAPL_LocStreamTransform( LOCSTREAM, TMPTILE  , PCU,     RC=STATUS); VERIFY_(STATUS)
          where(tiletype == MAPL_OCEAN)  PCUTILE = TMPTILE

          call MAPL_LocStreamTransform( LOCSTREAM, TMPTILE  , PLS,     RC=STATUS); VERIFY_(STATUS)
          where(tiletype == MAPL_OCEAN)  PLSTILE = TMPTILE

          call MAPL_LocStreamTransform( LOCSTREAM, TMPTILE  , SNOFL,   RC=STATUS); VERIFY_(STATUS)
          where(tiletype == MAPL_OCEAN)  SNOFLTILE = TMPTILE

          call MAPL_LocStreamTransform( LOCSTREAM, TMPTILE  , ICEFL,   RC=STATUS); VERIFY_(STATUS)
          where(tiletype == MAPL_OCEAN)  ICEFLTILE = TMPTILE

          call MAPL_LocStreamTransform( LOCSTREAM, TMPTILE  , FRZRFL,   RC=STATUS); VERIFY_(STATUS)
          where(tiletype == MAPL_OCEAN)  FRZRFLTILE = TMPTILE

       end if

       ! Adjust the discharge going to the ocean
       if(Precip_File /= "null" .and. DischargeAdjustFile /= "null") then

          call MAPL_GetPointer(INTERNAL, DISCHARGE_ADJUST, 'DISCHARGE_ADJUST',  RC=STATUS)
          VERIFY_(STATUS)

          DISCHARGETILE = DISCHARGETILE * DISCHARGE_ADJUST

       end if

       ! Create a version of the discharge (may be adjusted) on atmospheric grid for export.

       !  The total dicharge into the world ocean in kg s$^{-1}$ is $\Sigma_{i,j} D^a_{i,j} A_{i,j} = \Sigma_k D^t_k A_k $.
       !  On the left, $D^a_{i,j}$ is the DISCHARGE (kg m$^{-2}$ s$^{-1}$) at atmospheric grid cell $i,j$, and
       !  $A_{i,j}$ the area of the grid cell; on the right, $D^t_k$ is DISCHARGETILE at tile $k$, and $A_k$ is
       !  the area of the tile. $D^t_k$ is non-zero only at ocean tiles that are river outlets.

       if (associated(DISCHARGE)) then
          call MAPL_LocStreamTransform(LOCSTREAM, DISCHARGE, DISCHARGETILE, RC=STATUS)
          VERIFY_(STATUS)
       end if

    endif

! Run the Ocean
!--------------

    call DOTYPE(TYPE=OCEAN, RC=STATUS)
    VERIFY_(STATUS)


! Total precipitation diagnostic from SurfaceGridComp,
!  including any correction. The uncorrected comes from moist.
!-------------------------------------------------------------

  ! Convective Precipitation
  ! ------------------------
    call MAPL_GetPointer(EXPORT, CN_PRCP, 'CN_PRCP', ALLOC=.true., RC=STATUS)
    VERIFY_(STATUS)

    if(PRECIP_FILE /= "null") then
       TMPTILE = PCUTILE
       call MAPL_LocStreamTransform( LOCSTREAM, CN_PRCP, TMPTILE, RC=STATUS)
       VERIFY_(STATUS)
    else
       CN_PRCP = PRECCU
    endif
    CN_PRCP = MAX(CN_PRCP, 0.0)

  ! Total Precipitation
  ! -------------------
    call MAPL_GetPointer(EXPORT, PRECTOT, 'PRECTOT', ALLOC=.true., RC=STATUS)
    VERIFY_(STATUS)

    if(PRECIP_FILE /= "null") then
       TMPTILE = PCUTILE + PLSTILE + SNOFLTILE + ICEFLTILE  ! do *not* add FRZR, which is liquid not solid and (probably) incl. in PCUTILE+PCSTILE
       call MAPL_LocStreamTransform( LOCSTREAM, PRECTOT, TMPTILE, RC=STATUS)
       VERIFY_(STATUS)
    else
       PRECTOT = TPREC
    endif
    PRECTOT = MAX(PRECTOT, 0.0)

! New effective temperature and humidity
!---------------------------------------

    call MAPL_LocStreamTransform( LOCSTREAM,  DTS,    DTSTILE, RC=STATUS); VERIFY_(STATUS)
    call MAPL_LocStreamTransform( LOCSTREAM,  DQS,    DQSTILE, RC=STATUS); VERIFY_(STATUS)

    TH = TH + DTS
    QH = QH + DQS

! Transform other  exports from exchange grid to agcm grid
!---------------------------------------------------------

    if(associated(  RUNOFF)) then
       call MAPL_LocStreamTransform( LOCSTREAM, RUNOFF, RUNOFFTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(    TS)) then
       call MAPL_LocStreamTransform( LOCSTREAM,     TS,     TSTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(    QS)) then
       call MAPL_LocStreamTransform( LOCSTREAM,     QS,     QSTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( ALBVR)) then
       call MAPL_LocStreamTransform( LOCSTREAM,  ALBVR,  ALBVRTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( ALBVF)) then
       call MAPL_LocStreamTransform( LOCSTREAM,  ALBVF,  ALBVFTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( ALBNR)) then
       call MAPL_LocStreamTransform( LOCSTREAM,  ALBNR,  ALBNRTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( ALBNF)) then
       call MAPL_LocStreamTransform( LOCSTREAM,  ALBNF,  ALBNFTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( EMISS)) then
       call MAPL_LocStreamTransform( LOCSTREAM,  EMISS,  EMISSTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( FRI  )) then
       call MAPL_LocStreamTransform( LOCSTREAM,  FRI  ,     FRTILE, RC=STATUS) 
       VERIFY_(STATUS)
    endif
    if(associated( OFRI )) then
       call MAPL_LocStreamTransform( LOCSTREAM,  OFRI ,    OFRTILE, RC=STATUS) 
       VERIFY_(STATUS)
    endif
    if(associated(TSOIL1)) then
       call MAPL_LocStreamTransform( LOCSTREAM, TSOIL1, TSOIL1TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TSOIL2)) then
       call MAPL_LocStreamTransform( LOCSTREAM, TSOIL2, TSOIL2TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TSOIL3)) then
       call MAPL_LocStreamTransform( LOCSTREAM, TSOIL3, TSOIL3TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TSOIL4)) then
       call MAPL_LocStreamTransform( LOCSTREAM, TSOIL4, TSOIL4TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TSOIL5)) then
       call MAPL_LocStreamTransform( LOCSTREAM, TSOIL5, TSOIL5TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TSOIL6)) then
       call MAPL_LocStreamTransform( LOCSTREAM, TSOIL6, TSOIL6TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  WET1)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   WET1,   WET1TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  WET2)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   WET2,   WET2TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  WET3)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   WET3,   WET3TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  WCSF)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   WCSF,   WCSFTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  WCRZ)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   WCRZ,   WCRZTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  WCPR)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   WCPR,   WCPRTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(WESNN1)) then
       call MAPL_LocStreamTransform( LOCSTREAM, WESNN1, WESNN1TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(WESNN2)) then
       call MAPL_LocStreamTransform( LOCSTREAM, WESNN2, WESNN2TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(WESNN3)) then
       call MAPL_LocStreamTransform( LOCSTREAM, WESNN3, WESNN3TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( CAPAC)) then
       call MAPL_LocStreamTransform( LOCSTREAM,  CAPAC,  CAPACTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(ASNOW )) then
       call MAPL_LocStreamTransform( LOCSTREAM,  ASNOW,  ASNOWTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SHSNOW )) then
       call MAPL_LocStreamTransform( LOCSTREAM, SHSNOW, SHSNOWTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(AVETSNOW )) then
       call MAPL_LocStreamTransform( LOCSTREAM, AVETSNOW, AVETSNOWTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TPSNO )) then
       call MAPL_LocStreamTransform( LOCSTREAM,  TPSNO,  TPSNOTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TPUST )) then
       call MAPL_LocStreamTransform( LOCSTREAM,  TPUST,  TPUSTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TPSAT )) then
       call MAPL_LocStreamTransform( LOCSTREAM,  TPSAT,  TPSATTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TPWLT )) then
       call MAPL_LocStreamTransform( LOCSTREAM,  TPWLT,  TPWLTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TPSURF )) then
       call MAPL_LocStreamTransform( LOCSTREAM,  TPSURF, TPSURFTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(FRSAT )) then
       call MAPL_LocStreamTransform( LOCSTREAM,  FRSAT,  FRSATTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(FRUST )) then
       call MAPL_LocStreamTransform( LOCSTREAM,  FRUST,  FRUSTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(FRWLT )) then
       call MAPL_LocStreamTransform( LOCSTREAM,  FRWLT,  FRWLTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SNOMAS)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SNOMAS,   SNOWTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SNOWDP)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SNOWDP,   SNODTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  HLATN)) then
       call MAPL_LocStreamTransform( LOCSTREAM, HLATN,   HLATNTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(  HLATWTR)) then
       call MAPL_LocStreamTransform( LOCSTREAM, HLATWTR,   HLATWTRTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(  HLATICE)) then
       call MAPL_LocStreamTransform( LOCSTREAM, HLATICE,   HLATICETILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(    SHWTR)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   SHWTR,     SHWTRTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(    SHICE)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   SHICE,     SHICETILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(    TAUXW)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   TAUXW,     TAUXWTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(    TAUXI)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   TAUXI,     TAUXITILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(    TAUYW)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   TAUYW,     TAUYWTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(    TAUYI)) then
       call MAPL_LocStreamTransform( LOCSTREAM,   TAUYI,     TAUYITILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(  LWNDWTR)) then
       call MAPL_LocStreamTransform( LOCSTREAM, LWNDWTR,   LWNDWTRTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(  SWNDWTR)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SWNDWTR,   SWNDWTRTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(  LWNDICE)) then
       call MAPL_LocStreamTransform( LOCSTREAM, LWNDICE,   LWNDICETILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(  SWNDICE)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SWNDICE,   SWNDICETILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(TSKINW)) then
       call MAPL_LocStreamTransform( LOCSTREAM,TSKINW,TSKINWTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(TSKINICE)) then
       call MAPL_LocStreamTransform( LOCSTREAM,TSKINICE,TSKINICETILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( DCOOL)) then
       call MAPL_LocStreamTransform( LOCSTREAM, DCOOL, DCOOL_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( DWARM)) then
       call MAPL_LocStreamTransform( LOCSTREAM, DWARM, DWARM_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( TDROP)) then
       call MAPL_LocStreamTransform( LOCSTREAM, TDROP, TDROP_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( QCOOL)) then
       call MAPL_LocStreamTransform( LOCSTREAM, QCOOL, QCOOL_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( SWCOOL)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SWCOOL, SWCOOL_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( USTARW)) then
       call MAPL_LocStreamTransform( LOCSTREAM, USTARW, USTARW_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( TBAR)) then
       call MAPL_LocStreamTransform( LOCSTREAM, TBAR,  TBAR_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated( LCOOL)) then
       call MAPL_LocStreamTransform( LOCSTREAM, LCOOL, LCOOL_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( BCOOL)) then
       call MAPL_LocStreamTransform( LOCSTREAM, BCOOL, BCOOL_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( TDEL)) then
       call MAPL_LocStreamTransform( LOCSTREAM, TDEL, TDEL_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( TS_FOUND)) then
       call MAPL_LocStreamTransform( LOCSTREAM, TS_FOUND, TS_FOUND_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( SS_FOUND)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SS_FOUND, SS_FOUND_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( QWARM)) then
       call MAPL_LocStreamTransform( LOCSTREAM, QWARM, QWARM_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( SWWARM)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SWWARM, SWWARM_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( LANGM)) then
       call MAPL_LocStreamTransform( LOCSTREAM, LANGM, LANGM_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( PHIW)) then
       call MAPL_LocStreamTransform( LOCSTREAM, PHIW, PHIW_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( TAUTW)) then
       call MAPL_LocStreamTransform( LOCSTREAM, TAUTW, TAUTW_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( ZETA_W)) then
       call MAPL_LocStreamTransform( LOCSTREAM, ZETA_W, ZETA_W_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( TWMTF)) then
       call MAPL_LocStreamTransform( LOCSTREAM, TWMTF, TWMTF_TILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(  HICE)) then
       call MAPL_LocStreamTransform( LOCSTREAM, HICE,   HICETILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(  HSNO)) then
       call MAPL_LocStreamTransform( LOCSTREAM, HSNO,   HSNOTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(FRZMLT)) then
       call MAPL_LocStreamTransform( LOCSTREAM,FRZMLT,FRZMLTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(TSKINWCICE)) then
       call MAPL_LocStreamTransform( LOCSTREAM,TSKINWCICE,TSKINWCICETILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(ISTSFC)) then
       call MAPL_LocStreamTransform( LOCSTREAM,ISTSFC,ISTSFCTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(SSKINW)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SSKINW,SSKINWTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(MELTT)) then
       call MAPL_LocStreamTransform( LOCSTREAM,MELTT,MELTTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(MELTB)) then
       call MAPL_LocStreamTransform( LOCSTREAM,MELTB,MELTBTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(MELTS)) then
       call MAPL_LocStreamTransform( LOCSTREAM,MELTS,MELTSTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(MELTL)) then
       call MAPL_LocStreamTransform( LOCSTREAM,MELTL,MELTLTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(FRAZIL)) then
       call MAPL_LocStreamTransform( LOCSTREAM,FRAZIL,FRAZILTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(CONGEL)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CONGEL,CONGELTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(SNOICE)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SNOICE,SNOICETILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(DAIDTT)) then
       call MAPL_LocStreamTransform( LOCSTREAM,DAIDTT,DAIDTTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(DVIDTT)) then
       call MAPL_LocStreamTransform( LOCSTREAM,DVIDTT,DVIDTTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(DAIDTD)) then
       call MAPL_LocStreamTransform( LOCSTREAM,DAIDTD,DAIDTDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(DVIDTD)) then
       call MAPL_LocStreamTransform( LOCSTREAM,DVIDTD,DVIDTDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(FBOT)) then
       call MAPL_LocStreamTransform( LOCSTREAM,FBOT  ,FBOTTILE,   RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(HFLUX)) then
       call MAPL_LocStreamTransform( LOCSTREAM,HFLUX ,HFLUXTILE,   RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(WFLUX)) then
       call MAPL_LocStreamTransform( LOCSTREAM,WFLUX ,WFLUXTILE,   RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(SFLUX)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SFLUX ,SFLUXTILE,   RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(FSWTHRU)) then
       call MAPL_LocStreamTransform( LOCSTREAM,FSWTHRU ,FSWTHRUTILE,   RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(FSWABS)) then
       call MAPL_LocStreamTransform( LOCSTREAM,FSWABS ,FSWABSTILE,   RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(USTARI)) then
       call MAPL_LocStreamTransform( LOCSTREAM,USTARI ,USTARITILE,   RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(FHOCN )) then
       call MAPL_LocStreamTransform( LOCSTREAM,FHOCN ,FHOCNTILE,   RC=STATUS)
       VERIFY_(STATUS)
    endif

!****************************************************************************

    if(associated(  RAINOCN)) then
       call MAPL_LocStreamTransform( LOCSTREAM, RAINOCN,   RAINOCNTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(  SNOWOCN)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SNOWOCN,   SNOWOCNTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(  ICEFOCN)) then
       call MAPL_LocStreamTransform( LOCSTREAM, ICEFOCN,   ICEFOCNTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated( SPTOTOCN)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SPTOTOCN, SPTOTOCNTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif


    if(associated(  EVAPOU)) then
       call MAPL_LocStreamTransform( LOCSTREAM, EVAPOU,   EVAPOUTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  SUBLIM)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SUBLIM,   SUBLIMTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  SHOU)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SHOU,   SHOUTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(   LST)) then
       call MAPL_LocStreamTransform( LOCSTREAM,  LST,   LSTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  HLWUP)) then
       call MAPL_LocStreamTransform( LOCSTREAM, HLWUP,   HLWUPTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  LWNDSRF)) then
       call MAPL_LocStreamTransform( LOCSTREAM, LWNDSRF,   LWNDSRFTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  SWNDSRF)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SWNDSRF,   SWNDSRFTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  RUNSURF)) then
       call MAPL_LocStreamTransform( LOCSTREAM, RUNSURF,   RUNSURFTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  BASEFLOW)) then
       call MAPL_LocStreamTransform( LOCSTREAM, BASEFLOW,   BASEFLOWTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  ACCUM)) then
       call MAPL_LocStreamTransform( LOCSTREAM, ACCUM,   ACCUMTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  SMELT)) then
       call MAPL_LocStreamTransform( LOCSTREAM, SMELT,   SMELTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  EVEG)) then
       call MAPL_LocStreamTransform( LOCSTREAM, EVEG,   EVEGTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  EINT)) then
       call MAPL_LocStreamTransform( LOCSTREAM, EINT,   EINTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  EICE)) then
       call MAPL_LocStreamTransform( LOCSTREAM, EICE,   EICETILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  ESOI)) then
       call MAPL_LocStreamTransform( LOCSTREAM, ESOI,   ESOITILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(  ESNO)) then
       call MAPL_LocStreamTransform( LOCSTREAM, ESNO,   ESNOTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(WAT10CM)) then
       call MAPL_LocStreamTransform( LOCSTREAM, WAT10CM,WAT10CMTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(WATSOI)) then
       call MAPL_LocStreamTransform( LOCSTREAM, WATSOI, WATSOITILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(ICESOI)) then
       call MAPL_LocStreamTransform( LOCSTREAM, ICESOI, ICESOITILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(EVLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,EVLAND,EVLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(PRLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,PRLAND,PRLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SNOLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SNOLAND,SNOLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(DRPARLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,DRPARLAND,DRPARLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(DFPARLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,DFPARLAND,DFPARLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
     if(associated(LHSNOW)) then
        call MAPL_LocStreamTransform( LOCSTREAM,LHSNOW,LHSNOWTILE, RC=STATUS)
        VERIFY_(STATUS)
     endif
     if(associated(SWNETSNOW)) then
        call MAPL_LocStreamTransform( LOCSTREAM,SWNETSNOW,SWNETSNOWTILE, RC=STATUS)
        VERIFY_(STATUS)
     endif
     if(associated(LWUPSNOW)) then
        call MAPL_LocStreamTransform( LOCSTREAM,LWUPSNOW,LWUPSNOWTILE, RC=STATUS)
        VERIFY_(STATUS)
     endif
     if(associated(LWDNSNOW)) then
        call MAPL_LocStreamTransform( LOCSTREAM,LWDNSNOW,LWDNSNOWTILE, RC=STATUS)
        VERIFY_(STATUS)
     endif
     if(associated(TCSORIG)) then
        call MAPL_LocStreamTransform( LOCSTREAM,TCSORIG,TCSORIGTILE, RC=STATUS)
        VERIFY_(STATUS)
     endif
     if(associated(TPSN1IN)) then
        call MAPL_LocStreamTransform( LOCSTREAM,TPSN1IN,TPSN1INTILE, RC=STATUS)
        VERIFY_(STATUS)
     endif
     if(associated(TPSN1OUT)) then
        call MAPL_LocStreamTransform( LOCSTREAM,TPSN1OUT,TPSN1OUTTILE, RC=STATUS)
        VERIFY_(STATUS)
     endif
     if(associated(GHSNOW)) then
        call MAPL_LocStreamTransform( LOCSTREAM,GHSNOW,GHSNOWTILE, RC=STATUS)
        VERIFY_(STATUS)
     endif
    if(associated(LHLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,LHLAND,LHLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SHLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SHLAND,SHLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SWLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SWLAND,SWLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SWDOWNLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SWDOWNLAND,SWDOWNLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(LWLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,LWLAND,LWLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(GHLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,GHLAND,GHLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(GHTSKIN)) then
       call MAPL_LocStreamTransform( LOCSTREAM,GHTSKIN,GHTSKINTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SMLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SMLAND,SMLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(QINFIL)) then
       call MAPL_LocStreamTransform( LOCSTREAM,QINFIL,QINFILTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TWLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,TWLAND,TWLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TELAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,TELAND,TELANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(TSLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,TSLAND,TSLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(DWLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,DWLAND,DWLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(DHLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,DHLAND,DHLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SPLAND)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SPLAND,SPLANDTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SPLH  )) then
       call MAPL_LocStreamTransform( LOCSTREAM,SPLH  ,SPLHTILE  , RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SPWATR)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SPWATR,SPWATRTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SPSNOW)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SPSNOW,SPSNOWTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

    DO N = 1, N_SNOW_LAND
       if(associated(RDU001)) call MAPL_LocStreamTransform( LOCSTREAM,RDU001(:,:,N) ,RDU001TILE(:,N), RC=STATUS); VERIFY_(STATUS)
       if(associated(RDU002)) call MAPL_LocStreamTransform( LOCSTREAM,RDU002(:,:,N) ,RDU002TILE(:,N), RC=STATUS); VERIFY_(STATUS)
       if(associated(RDU003)) call MAPL_LocStreamTransform( LOCSTREAM,RDU003(:,:,N) ,RDU003TILE(:,N), RC=STATUS); VERIFY_(STATUS)
       if(associated(RDU004)) call MAPL_LocStreamTransform( LOCSTREAM,RDU004(:,:,N) ,RDU004TILE(:,N), RC=STATUS); VERIFY_(STATUS)
       if(associated(RDU005)) call MAPL_LocStreamTransform( LOCSTREAM,RDU005(:,:,N) ,RDU005TILE(:,N), RC=STATUS); VERIFY_(STATUS)
       if(associated(RBC001)) call MAPL_LocStreamTransform( LOCSTREAM,RBC001(:,:,N) ,RBC001TILE(:,N), RC=STATUS); VERIFY_(STATUS)
       if(associated(RBC002)) call MAPL_LocStreamTransform( LOCSTREAM,RBC002(:,:,N) ,RBC002TILE(:,N), RC=STATUS); VERIFY_(STATUS)
       if(associated(ROC001)) call MAPL_LocStreamTransform( LOCSTREAM,ROC001(:,:,N) ,ROC001TILE(:,N), RC=STATUS); VERIFY_(STATUS)
       if(associated(ROC002)) call MAPL_LocStreamTransform( LOCSTREAM,ROC002(:,:,N) ,ROC002TILE(:,N), RC=STATUS); VERIFY_(STATUS)
    END DO

    if(associated(RMELTDU001 ))call MAPL_LocStreamTransform(LOCSTREAM,RMELTDU001 ,RMELTDU001TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(RMELTDU002 ))call MAPL_LocStreamTransform(LOCSTREAM,RMELTDU002 ,RMELTDU002TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(RMELTDU003 ))call MAPL_LocStreamTransform(LOCSTREAM,RMELTDU003 ,RMELTDU003TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(RMELTDU004 ))call MAPL_LocStreamTransform(LOCSTREAM,RMELTDU004 ,RMELTDU004TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(RMELTDU005 ))call MAPL_LocStreamTransform(LOCSTREAM,RMELTDU005 ,RMELTDU005TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(RMELTBC001 ))call MAPL_LocStreamTransform(LOCSTREAM,RMELTBC001 ,RMELTBC001TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(RMELTBC002 ))call MAPL_LocStreamTransform(LOCSTREAM,RMELTBC002 ,RMELTBC002TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(RMELTOC001 ))call MAPL_LocStreamTransform(LOCSTREAM,RMELTOC001 ,RMELTOC001TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(RMELTOC002 ))call MAPL_LocStreamTransform(LOCSTREAM,RMELTOC002 ,RMELTOC002TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(PEATCLSM_WATERLEVEL))call MAPL_LocStreamTransform(LOCSTREAM,PEATCLSM_WATERLEVEL,PEATCLSM_WATERLEVELTILE,RC=STATUS); VERIFY_(STATUS)
    if(associated(PEATCLSM_FSWCHANGE ))call MAPL_LocStreamTransform(LOCSTREAM,PEATCLSM_FSWCHANGE ,PEATCLSM_FSWCHANGETILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(DZGT1 ))call MAPL_LocStreamTransform(LOCSTREAM,DZGT1 ,DZGT1TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(DZGT2 ))call MAPL_LocStreamTransform(LOCSTREAM,DZGT2 ,DZGT2TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(DZGT3 ))call MAPL_LocStreamTransform(LOCSTREAM,DZGT3 ,DZGT3TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(DZGT4 ))call MAPL_LocStreamTransform(LOCSTREAM,DZGT4 ,DZGT4TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(DZGT5 ))call MAPL_LocStreamTransform(LOCSTREAM,DZGT5 ,DZGT5TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(DZGT6 ))call MAPL_LocStreamTransform(LOCSTREAM,DZGT6 ,DZGT6TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(DZPR  ))call MAPL_LocStreamTransform(LOCSTREAM,DZPR  ,DZPRTILE , RC=STATUS); VERIFY_(STATUS)
    if(associated(DZRZ  ))call MAPL_LocStreamTransform(LOCSTREAM,DZRZ  ,DZRZTILE , RC=STATUS); VERIFY_(STATUS)
    if(associated(DZSF  ))call MAPL_LocStreamTransform(LOCSTREAM,DZSF  ,DZSFTILE , RC=STATUS); VERIFY_(STATUS)
    if(associated(DZTS  ))call MAPL_LocStreamTransform(LOCSTREAM,DZTS  ,DZTSTILE , RC=STATUS); VERIFY_(STATUS)
    if(associated(WPWET ))call MAPL_LocStreamTransform(LOCSTREAM,WPWET ,WPWETTILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(WPEMW ))call MAPL_LocStreamTransform(LOCSTREAM,WPEMW ,WPEMWTILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(WPMC  ))call MAPL_LocStreamTransform(LOCSTREAM,WPMC  ,WPMCTILE , RC=STATUS); VERIFY_(STATUS)
    if(associated(CDCR2 ))call MAPL_LocStreamTransform(LOCSTREAM,CDCR2 ,CDCR2TILE, RC=STATUS); VERIFY_(STATUS)
    if(associated(POROS ))call MAPL_LocStreamTransform(LOCSTREAM,POROS ,POROSTILE, RC=STATUS); VERIFY_(STATUS)

    if(associated(CNLAI)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNLAI ,CNLAITILE , RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNTLAI)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNTLAI,CNTLAITILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNSAI)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNSAI ,CNSAITILE , RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNTOTC)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNTOTC,CNTOTCTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNVEGC)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNVEGC,CNVEGCTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNROOT)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNROOT,CNROOTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNFROOTC)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNFROOTC,CNFROOTCTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNNPP)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNNPP ,CNNPPTILE , RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNGPP)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNGPP ,CNGPPTILE , RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNSR)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNSR  ,CNSRTILE  , RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNNEE)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNNEE ,CNNEETILE , RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNXSMR)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNXSMR,CNXSMRTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNADD)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNADD ,CNADDTILE , RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNLOSS)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNLOSS,CNLOSSTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNBURN)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNBURN,CNBURNTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(PARABS)) then
       call MAPL_LocStreamTransform( LOCSTREAM,PARABS,PARABSTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(PARINC)) then
       call MAPL_LocStreamTransform( LOCSTREAM,PARINC,PARINCTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SCSAT)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SCSAT ,SCSATTILE , RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SCUNS)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SCUNS ,SCUNSTILE , RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(BTRANT)) then
       call MAPL_LocStreamTransform( LOCSTREAM,BTRANT,BTRANTTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(SIF)) then
       call MAPL_LocStreamTransform( LOCSTREAM,SIF   ,SIFTILE   , RC=STATUS)
       VERIFY_(STATUS)
    endif
    if(associated(CNFSEL)) then
       call MAPL_LocStreamTransform( LOCSTREAM,CNFSEL,CNFSELTILE, RC=STATUS)
       VERIFY_(STATUS)
    endif

! Fire danger
    if (associated(FFMC)) then
       call MAPL_LocStreamTransform(LOCSTREAM, FFMC, FFMCTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(GFMC)) then
       call MAPL_LocStreamTransform(LOCSTREAM, GFMC, GFMCTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(DMC)) then
       call MAPL_LocStreamTransform(LOCSTREAM, DMC, DMCTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(DC)) then
       call MAPL_LocStreamTransform(LOCSTREAM, DC, DCTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(ISI)) then
       call MAPL_LocStreamTransform(LOCSTREAM, ISI, ISITILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(BUI)) then
       call MAPL_LocStreamTransform(LOCSTREAM, BUI, BUITILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(FWI)) then
       call MAPL_LocStreamTransform(LOCSTREAM, FWI, FWITILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(DSR)) then
       call MAPL_LocStreamTransform(LOCSTREAM, DSR, DSRTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if

    if (associated(FFMC_DAILY)) then
       call MAPL_LocStreamTransform(LOCSTREAM, FFMC_DAILY, FFMCDAILYTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(DMC_DAILY)) then
       call MAPL_LocStreamTransform(LOCSTREAM, DMC_DAILY, DMCDAILYTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(DC_DAILY)) then
       call MAPL_LocStreamTransform(LOCSTREAM, DC_DAILY, DCDAILYTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(ISI_DAILY)) then
       call MAPL_LocStreamTransform(LOCSTREAM, ISI_DAILY, ISIDAILYTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(BUI_DAILY)) then
       call MAPL_LocStreamTransform(LOCSTREAM, BUI_DAILY, BUIDAILYTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(FWI_DAILY)) then
       call MAPL_LocStreamTransform(LOCSTREAM, FWI_DAILY, FWIDAILYTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(DSR_DAILY)) then
       call MAPL_LocStreamTransform(LOCSTREAM, DSR_DAILY, DSRDAILYTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if

    if (associated(FFMC_DAILY_)) then
       call MAPL_LocStreamTransform(LOCSTREAM, FFMC_DAILY_, FFMCDAILYTILE_, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(DMC_DAILY_)) then
       call MAPL_LocStreamTransform(LOCSTREAM, DMC_DAILY_, DMCDAILYTILE_, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(DC_DAILY_)) then
       call MAPL_LocStreamTransform(LOCSTREAM, DC_DAILY_, DCDAILYTILE_, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(ISI_DAILY_)) then
       call MAPL_LocStreamTransform(LOCSTREAM, ISI_DAILY_, ISIDAILYTILE_, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(BUI_DAILY_)) then
       call MAPL_LocStreamTransform(LOCSTREAM, BUI_DAILY_, BUIDAILYTILE_, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(FWI_DAILY_)) then
       call MAPL_LocStreamTransform(LOCSTREAM, FWI_DAILY_, FWIDAILYTILE_, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(DSR_DAILY_)) then
       call MAPL_LocStreamTransform(LOCSTREAM, DSR_DAILY_, DSRDAILYTILE_, RC=STATUS)
       VERIFY_(STATUS)
    end if
    if (associated(VPD)) then
       call MAPL_LocStreamTransform(LOCSTREAM, VPD, VPDTILE, RC=STATUS)
       VERIFY_(STATUS)
    end if


! Fill exports computed on agcm grid
!-----------------------------------

    if(associated( DELTS)) DELTS = DTS
    if(associated( DELQS)) DELQS = DQS
    if(associated( DELSS)) DELSS = MAPL_CP*DTS
    if(associated( DELUS)) DELUS = 0.0
    if(associated( DELVS)) DELVS = 0.0
    if(associated( TAUXO)) TAUXO = -TAUX
    if(associated( TAUYO)) TAUYO = -TAUY
    if(associated( EVAPO)) EVAPO = EVAP + DEVAP*DQS
    if(associated(   SHO))   SHO = SH   + DSH  *DTS *MAPL_CP
    if(associated( DLQLL)) DLQLL = 0.0
    if(associated( DLQIL)) DLQIL = 0.0

    TAU = sqrt(TAUX**2+TAUY**2)

    if(  associated(USTAR ) .or. associated(TSTAR ) .or.      &
         associated(QSTAR )                              ) then
       FAC   = sqrt(TAU/RHOS)
       if(associated( USTAR)) USTAR = FAC
       if(associated( TSTAR)) TSTAR = (SH/MAPL_CP + DSH  *DTS)/(RHOS*FAC)
       if(associated( QSTAR)) QSTAR = (EVAP       + DEVAP*DQS)/(RHOS*FAC)
    end if

    if (DO_DATA_ATM4OCN) then
       ! dataAtm operates only on "saltwater" tiles.
       ! we need to handle grid boxes withot any ocean
       ! and avoid division by 0
       where (CN == MAPL_Undef)
          CM = 0.01
          CT = 0.01
          CQ = 0.01
          CN = 0.01
          D0 = 0.0
       end where
    end if

    FAC = sqrt(CN)/MAPL_KARMAN
    Z0  = max((DZ-D0),10.)/(exp(1.0/FAC)-1.0)

! if we are running Louis sfc layer, get these exports from the gridded Louis fluxes
  if (CHOOSEMOSFC.eq.0) then

! 50m winds

    if(  associated(  U50M) .or. associated(  V50M)      ) then
       TMP = alog(1.0 + (50.)/Z0)*FAC
       if(associated(  U50M))  U50M = UH - TAUX*TMP/CM
       if(associated(  V50M))  V50M = VH - TAUY*TMP/CM
    end if

! 10m

    if(  associated(  U10M) .or. associated(  V10M) .or.      &
         associated(  T10M) .or. associated(  Q10M)      ) then
       TMP = alog(1.0 + (10.)/Z0)*FAC
       if(associated(  U10M) .and. associated(  V10M)) then
          U10M = UH - TAUX*TMP/CM
          V10M = VH - TAUY*TMP/CM
          if(associated(UU10M)) UU10M = sqrt(U10M**2 + V10M**2)
       else
          if(associated(  U10M))  U10M = UH - TAUX*TMP/CM
          if(associated(  V10M))  V10M = VH - TAUY*TMP/CM
       end if
       if(associated(T10M))  T10M =      TH - (SH/MAPL_CP + DSH  *DTS)*TMP/CT
       if(associated(Q10M))  Q10M = max( QH - (EVAP       + DEVAP*DQS)*TMP/CQ , 0.0 )
    end if
  endif

! For either sfc layer, compute neutral 10 m winds here
    if(  associated(  U10N) .or. associated(  V10N)      ) then
       TMP = alog(1.0 + (10.)/Z0)/(sqrt(max(TAU*RHOS,tiny(1.0))) * MAPL_KARMAN )
       if(associated(  U10N))  U10N = UH - TAUX*TMP
       if(associated(  V10N))  V10N = VH - TAUY*TMP
    end if

! if we are running Louis sfc layer, get these exports from the gridded Louis fluxes
  if (CHOOSEMOSFC.eq.0) then

! 2m; do this all the time for max/min

     if(  associated(U2M) .or. associated(V2M) .or.        &
          associated(T2M) .or. associated(Q2M)         ) then
       TMP = alog(1.0 + (2.0)/Z0)*FAC
       if(associated(  U2M) .or. associated(  V2M)) then
          if(associated(  U2M))  U2M = UH - TAUX*TMP/CM
          if(associated(  V2M))  V2M = VH - TAUY*TMP/CM
       end if
       TTM =       TH - (SH/MAPL_CP + DSH  *DTS)*TMP/CT
       QTM =  max( QH - (EVAP       + DEVAP*DQS)*TMP/CQ , 0.0 )
       if(associated(T2M))  T2M = TTM
       if(associated(Q2M))  Q2M = QTM
     endif

     if(associated(T2MDEW)) then
      T2MDEW = TTM
      do i = 1,4
       T2MDEW = T2MDEW + (QTM-GEOS_QSAT(T2MDEW,PS,PASCALS=.true.))/GEOS_DQSAT(T2MDEW,PS,PASCALS=.true.)
      enddo
     endif

     if(associated(T2MWET)) then
      T2MWET = TTM
      allocate(ALHX(IM,JM),STAT=STATUS)
      VERIFY_(STATUS)
      do i = 1,10
       ALHX = (MAX(MIN(1.-(T2MWET-233.16)/40.,1.),0.))**4
       ALHX = (1.0-ALHX)*MAPL_ALHL + ALHX*MAPL_ALHS
       T2MWET = T2MWET + ((ALHX/MAPL_CP)/(1.+(ALHX/MAPL_CP)*GEOS_DQSAT(T2MWET,PS,PASCALS=.true.)))* &
                                  (QTM-GEOS_QSAT(T2MWET,PS,PASCALS=.true.))
      enddo
      deallocate(ALHX)
     endif

     QTM = min(QTM/GEOS_QSAT(TTM,PS,PASCALS=.true.),1.0)*100.

     if(associated(RH2M))  RH2M = QTM

  endif   ! end choosemosfc sequence for louis scheme

! surface air values

    if(  associated(  UAX) .or. associated(  VAX) .or.      &
         associated(  TA) .or. associated(  QA)      ) then
       if(associated(  UAX) .or. associated(  VAX)) then
          if(associated(  UAX))  UAX = UH - TAUX/CM
          if(associated(  VAX))  VAX = VH - TAUY/CM
       end if
       if(associated(TA))  TA = TH - (SH/MAPL_CP + DSH  *DTS)/CT
       if(associated(QA))  QA = QH - (EVAP       + DEVAP*DQS)/CQ
    end if

! Set Integer LWI flag
!---------------------

      if( associated(LWI) ) then
          call MAPL_GetPointer( EXPORT, FROCEAN, 'FROCEAN', RC=STATUS )
          VERIFY_(STATUS)
          call MAPL_GetPointer( EXPORT, FRLAKE,  'FRLAKE' , RC=STATUS )
          VERIFY_(STATUS)
                                                      LWI = 1.0  ! Land
          where ( FROCEAN+FRLAKE >= 0.6             ) LWI = 0.0  ! Water
          where ( LWI==0 .and. FRI>0.5              ) LWI = 2.0  ! Ice
          where ( LWI==0 .and. TS<271.40            ) LWI = 2.0  ! Ice
      endif


! Fill WET1_FOR_CHEM over non-land points to 1.0
!-----------------------------------------------

! NOTE: GOCART's dust code expects WET1 to have all the cells with MAPL_UNDEF
!       (aka not land) to be replaced with 1.0. We want WET1 to have
!       MAPL_UNDEF over non-land points, so we need a separate export to pass
!       to GOCART.

      if( associated(WET1_FOR_CHEM) ) then
          WET1_FOR_CHEM = WET1
          where(WET1_FOR_CHEM == MAPL_UNDEF) WET1_FOR_CHEM = 1.0
      endif


! Fill imports/exports for OBIO
!-------------------------------
      if((DO_OBIO/=0) .OR. (ATM_CO2 == ATM_CO2_FOUR)) then
        call OBIO_fillExports(OCEAN, IMPORT,&
                              LOCSTREAM, GIM,&
                              surf_internal_state%xform_in(OCEAN), &
                              NT, NB_CHOU,&
                              CO2SC, DRBAND, DFBAND, &
                              CO2SCTILE, DRBANDTILE, DFBANDTILE, RC)
      else
        nullify(  CO2SCTILE   )
        nullify(  DRBANDTILE  )
        nullify(  DFBANDTILE  )
      endif

! Moved change of units for soil temperature export variables down to Catch[CN] Gridcomp.
! With this change, gridded TSOIL[n] exports from Surface and tile-space TP[n] exports
! from Catch are now consistently in units of Kelvin.
! - rreichle & borescan, 6 Nov 2020
!
!-----------------------
!      if( associated(TSOIL1) ) then
!             where ( TSOIL1 /= MAPL_Undef ) TSOIL1 = TSOIL1 + MAPL_TICE
!      endif
!      if( associated(TSOIL2) ) then
!             where ( TSOIL2 /= MAPL_Undef ) TSOIL2 = TSOIL2 + MAPL_TICE
!      endif
!      if( associated(TSOIL3) ) then
!             where ( TSOIL3 /= MAPL_Undef ) TSOIL3 = TSOIL3 + MAPL_TICE
!      endif
!      if( associated(TSOIL4) ) then
!             where ( TSOIL4 /= MAPL_Undef ) TSOIL4 = TSOIL4 + MAPL_TICE
!      endif
!      if( associated(TSOIL5) ) then
!             where ( TSOIL5 /= MAPL_Undef ) TSOIL5 = TSOIL5 + MAPL_TICE
!      endif
!      if( associated(TSOIL6) ) then
!             where ( TSOIL6 /= MAPL_Undef ) TSOIL6 = TSOIL6 + MAPL_TICE
!      endif

! Fill SNOMAS over Glaciers
!--------------------------

!      if( associated(SNOMAS) ) then
!          call MAPL_GetPointer( EXPORT, FRLANDICE, 'FRLANDICE', RC=STATUS )
!          VERIFY_(STATUS)
!         where ( SNOMAS /= MAPL_UNDEF ) SNOMAS = SNOMAS*1.E-3
!          where ( FRLANDICE > 0.9 ) SNOMAS = 4
!      endif

! Clean-up
!---------

    deallocate(TMP)
    deallocate(TTM)
    deallocate(QTM)
    deallocate(Z0 )
    deallocate(FAC)
    deallocate(TAU)
    deallocate(DTS)
    deallocate(DQS)

    deallocate(DRPAR)
    deallocate(DFPAR)
    deallocate(DRNIR)
    deallocate(DFNIR)
    deallocate(DRUVR)
    deallocate(DFUVR)
    deallocate(ZTH  )
    deallocate(SLR  )

    deallocate(   DTSTILE)
    deallocate(   DQSTILE)
    deallocate(   SLITILE)
    deallocate(   ZTHTILE)
    deallocate(    PSTILE)
    deallocate(    UUTILE)
    deallocate(    DZTILE)
    deallocate( SNOFLTILE)
    deallocate( ICEFLTILE)
    deallocate( FRZRFLTILE)
    deallocate(   TMPTILE)
    deallocate(   PCUTILE)
    deallocate(   PLSTILE)
    deallocate(   LWBTILE)
    deallocate(   DRPTILE)
    deallocate(   DFPTILE)
    deallocate(   SHFTILE)
    deallocate(   DSHTILE)
    deallocate(   EVPTILE)
    deallocate(   DEVTILE)
    deallocate(   ALWTILE)
    deallocate(   BLWTILE)
    deallocate(  TAUXTILE)
    deallocate(  TAUYTILE)
    deallocate(   DFNTILE)
    deallocate(   DRNTILE)
    deallocate(   DFUTILE)
    deallocate(   DRUTILE)

    if(associated( CO2SCTILE   )) deallocate( CO2SCTILE   )
    if(associated(  DUDPTILE   )) deallocate(  DUDPTILE   )
    if(associated(  DUSVTILE   )) deallocate(  DUSVTILE   )
    if(associated(  DUWTTILE   )) deallocate(  DUWTTILE   )
    if(associated(  DUSDTILE   )) deallocate(  DUSDTILE   )
    if(associated(  BCDPTILE   )) deallocate(  BCDPTILE   )
    if(associated(  BCSVTILE   )) deallocate(  BCSVTILE   )
    if(associated(  BCWTTILE   )) deallocate(  BCWTTILE   )
    if(associated(  BCSDTILE   )) deallocate(  BCSDTILE   )
    if(associated(  OCDPTILE   )) deallocate(  OCDPTILE   )
    if(associated(  OCSVTILE   )) deallocate(  OCSVTILE   )
    if(associated(  OCWTTILE   )) deallocate(  OCWTTILE   )
    if(associated(  OCSDTILE   )) deallocate(  OCSDTILE   )
    if(associated(  SUDPTILE   )) deallocate(  SUDPTILE   )
    if(associated(  SUSVTILE   )) deallocate(  SUSVTILE   )
    if(associated(  SUWTTILE   )) deallocate(  SUWTTILE   )
    if(associated(  SUSDTILE   )) deallocate(  SUSDTILE   )
    if(associated(  SSDPTILE   )) deallocate(  SSDPTILE   )
    if(associated(  SSSVTILE   )) deallocate(  SSSVTILE   )
    if(associated(  SSWTTILE   )) deallocate(  SSWTTILE   )
    if(associated(  SSSDTILE   )) deallocate(  SSSDTILE   )
    if(associated(  DRBANDTILE )) deallocate(  DRBANDTILE )
    if(associated(  DFBANDTILE )) deallocate(  DFBANDTILE )

    if(associated(     DTSDTTILE)) deallocate(      DTSDTTILE)

    if(associated(TSTILE      )) deallocate(TSTILE      )
    if(associated(QSTILE      )) deallocate(QSTILE      )
    if(associated(THTILE      )) deallocate(THTILE      )
    if(associated(QHTILE      )) deallocate(QHTILE      )
    if(associated(UHTILE      )) deallocate(UHTILE      )
    if(associated(VHTILE      )) deallocate(VHTILE      )
    if(associated(CTTILE      )) deallocate(CTTILE      )
    if(associated(CQTILE      )) deallocate(CQTILE      )
    if(associated(CMTILE      )) deallocate(CMTILE      )
    if(associated(TSOIL1TILE  )) deallocate(TSOIL1TILE  )
    if(associated(TSOIL2TILE  )) deallocate(TSOIL2TILE  )
    if(associated(TSOIL3TILE  )) deallocate(TSOIL3TILE  )
    if(associated(TSOIL4TILE  )) deallocate(TSOIL4TILE  )
    if(associated(TSOIL5TILE  )) deallocate(TSOIL5TILE  )
    if(associated(TSOIL6TILE  )) deallocate(TSOIL6TILE  )
    if(associated(WET1TILE    )) deallocate(WET1TILE    )
    if(associated(WET2TILE    )) deallocate(WET2TILE    )
    if(associated(WET3TILE    )) deallocate(WET3TILE    )
    if(associated(WCSFTILE    )) deallocate(WCSFTILE    )
    if(associated(WCRZTILE    )) deallocate(WCRZTILE    )
    if(associated(WCPRTILE    )) deallocate(WCPRTILE    )
    if(associated(WESNN1TILE  )) deallocate(WESNN1TILE  )
    if(associated(WESNN2TILE  )) deallocate(WESNN2TILE  )
    if(associated(WESNN3TILE  )) deallocate(WESNN3TILE  )
    if(associated(CAPACTILE   )) deallocate(CAPACTILE   )
    if(associated(ASNOWTILE   )) deallocate(ASNOWTILE   )
    if(associated(SHSNOWTILE  )) deallocate(SHSNOWTILE  )
    if(associated(AVETSNOWTILE)) deallocate(AVETSNOWTILE)
    if(associated(TPSNOTILE   )) deallocate(TPSNOTILE   )
    if(associated(TPUSTTILE   )) deallocate(TPUSTTILE   )
    if(associated(TPSATTILE   )) deallocate(TPSATTILE   )
    if(associated(TPWLTTILE   )) deallocate(TPWLTTILE   )
    if(associated(TPSURFTILE  )) deallocate(TPSURFTILE  )
    if(associated(FRSATTILE   )) deallocate(FRSATTILE   )
    if(associated(FRUSTTILE   )) deallocate(FRUSTTILE   )
    if(associated(FRWLTTILE   )) deallocate(FRWLTTILE   )
    if(associated(SNOWTILE    )) deallocate(SNOWTILE    )
    if(associated(SNODTILE    )) deallocate(SNODTILE    )
    if(associated(HLATNTILE   )) deallocate(HLATNTILE   )

    if(associated(HLATWTRTILE   )) deallocate(HLATWTRTILE   )
    if(associated(HLATICETILE   )) deallocate(HLATICETILE   )
    if(associated(  SHWTRTILE   )) deallocate(  SHWTRTILE   )
    if(associated(  SHICETILE   )) deallocate(  SHICETILE   )
    if(associated(  TAUXWTILE   )) deallocate(  TAUXWTILE   )
    if(associated(  TAUXITILE   )) deallocate(  TAUXITILE   )
    if(associated(  TAUYWTILE   )) deallocate(  TAUYWTILE   )
    if(associated(  TAUYITILE   )) deallocate(  TAUYITILE   )
    if(associated(LWNDWTRTILE   )) deallocate(LWNDWTRTILE   )
    if(associated(SWNDWTRTILE   )) deallocate(SWNDWTRTILE   )
    if(associated(LWNDICETILE   )) deallocate(LWNDICETILE   )
    if(associated(SWNDICETILE   )) deallocate(SWNDICETILE   )
    if(associated(SNOWOCNTILE   )) deallocate(SNOWOCNTILE   )
    if(associated(ICEFOCNTILE   )) deallocate(ICEFOCNTILE   )
    if(associated(SPTOTOCNTILE  )) deallocate(SPTOTOCNTILE  )
    if(associated(RAINOCNTILE   )) deallocate(RAINOCNTILE   )
    if(associated(TSKINWTILE  )) deallocate(TSKINWTILE  )
    if(associated(TSKINICETILE  )) deallocate(TSKINICETILE  )

    if(associated(DCOOL_TILE    )) deallocate(DCOOL_TILE     )
    if(associated(DWARM_TILE    )) deallocate(DWARM_TILE     )
    if(associated(TDROP_TILE    )) deallocate(TDROP_TILE     )
    if(associated(QCOOL_TILE    )) deallocate(QCOOL_TILE     )
    if(associated(SWCOOL_TILE   )) deallocate(SWCOOL_TILE    )
    if(associated(USTARW_TILE   )) deallocate(USTARW_TILE    )
    if(associated(TBAR_TILE     )) deallocate(TBAR_TILE      )
    if(associated(LCOOL_TILE    )) deallocate(LCOOL_TILE     )
    if(associated(BCOOL_TILE    )) deallocate(BCOOL_TILE     )
    if(associated(TDEL_TILE     )) deallocate(TDEL_TILE      )
    if(associated(TS_FOUND_TILE )) deallocate(TS_FOUND_TILE  )
    if(associated(SS_FOUND_TILE )) deallocate(SS_FOUND_TILE  )
    if(associated(QWARM_TILE    )) deallocate(QWARM_TILE     )
    if(associated(SWWARM_TILE   )) deallocate(SWWARM_TILE    )
    if(associated(LANGM_TILE    )) deallocate(LANGM_TILE     )
    if(associated(PHIW_TILE     )) deallocate(PHIW_TILE      )
    if(associated(TAUTW_TILE    )) deallocate(TAUTW_TILE     )
    if(associated(ZETA_W_TILE   )) deallocate(ZETA_W_TILE    )
    if(associated(TWMTF_TILE    )) deallocate(TWMTF_TILE     )

    if(associated(HICETILE    )) deallocate(HICETILE    )
    if(associated(HSNOTILE    )) deallocate(HSNOTILE    )
    if(associated(FRZMLTTILE  )) deallocate(FRZMLTTILE  )
    if(associated(TSKINWCICETILE  )) deallocate(TSKINWCICETILE  )
    if(associated(ISTSFCTILE  )) deallocate(ISTSFCTILE  )
    if(associated(SSKINWTILE  )) deallocate(SSKINWTILE  )
    if(associated(MELTTTILE   )) deallocate(MELTTTILE  )
    if(associated(MELTBTILE   )) deallocate(MELTBTILE  )
    if(associated(MELTSTILE   )) deallocate(MELTSTILE  )
    if(associated(MELTLTILE   )) deallocate(MELTLTILE  )
    if(associated(FRAZILTILE  )) deallocate(FRAZILTILE )
    if(associated(CONGELTILE  )) deallocate(CONGELTILE )
    if(associated(SNOICETILE  )) deallocate(SNOICETILE )
    if(associated(DAIDTTTILE  )) deallocate(DAIDTTTILE )
    if(associated(DVIDTTTILE  )) deallocate(DVIDTTTILE )
    if(associated(DAIDTDTILE  )) deallocate(DAIDTDTILE )
    if(associated(DVIDTDTILE  )) deallocate(DVIDTDTILE )
    if(associated(FBOTTILE    )) deallocate(FBOTTILE   )
    if(associated(HFLUXTILE   )) deallocate(HFLUXTILE  )
    if(associated(WFLUXTILE   )) deallocate(WFLUXTILE  )
    if(associated(SFLUXTILE   )) deallocate(SFLUXTILE  )
    if(associated(FSWTHRUTILE )) deallocate(FSWTHRUTILE)
    if(associated(FSWABSTILE  )) deallocate(FSWABSTILE )
    if(associated(USTARITILE  )) deallocate(USTARITILE )
    if(associated(FHOCNTILE   )) deallocate(FHOCNTILE  )

    if(associated(EVAPOUTILE  )) deallocate(EVAPOUTILE  )
    if(associated(SUBLIMTILE  )) deallocate(SUBLIMTILE  )
    if(associated(SHOUTILE    )) deallocate(SHOUTILE    )
    if(associated(LSTTILE     )) deallocate(LSTTILE     )
    if(associated(HLWUPTILE   )) deallocate(HLWUPTILE   )
    if(associated(LWNDSRFTILE )) deallocate(LWNDSRFTILE )
    if(associated(SWNDSRFTILE )) deallocate(SWNDSRFTILE )
    if(associated(RUNOFFTILE  )) deallocate(RUNOFFTILE  )
    if(associated(DISCHARGETILE))deallocate(DISCHARGETILE)
    if(associated(RUNSURFTILE )) deallocate(RUNSURFTILE )
    if(associated(BASEFLOWTILE)) deallocate(BASEFLOWTILE)
    if(associated(ACCUMTILE   )) deallocate(ACCUMTILE   )
    if(associated(SMELTTILE   )) deallocate(SMELTTILE   )
    if(associated(EVEGTILE    )) deallocate(EVEGTILE    )
    if(associated(EINTTILE    )) deallocate(EINTTILE    )
    if(associated(EICETILE    )) deallocate(EICETILE    )
    if(associated(ESNOTILE    )) deallocate(ESNOTILE    )
    if(associated(ESOITILE    )) deallocate(ESOITILE    )
    if(associated(WAT10CMTILE )) deallocate(WAT10CMTILE )
    if(associated(WATSOITILE  )) deallocate(WATSOITILE  )
    if(associated(ICESOITILE  )) deallocate(ICESOITILE  )
    if(associated(EVLANDTILE  )) deallocate(EVLANDTILE  )
    if(associated(PRLANDTILE  )) deallocate(PRLANDTILE  )
    if(associated(SNOLANDTILE  )) deallocate(SNOLANDTILE  )
    if(associated(DRPARLANDTILE  )) deallocate(DRPARLANDTILE  )
    if(associated(DFPARLANDTILE  )) deallocate(DFPARLANDTILE  )
    if(associated(LHSNOWTILE  )) deallocate(LHSNOWTILE  )
    if(associated(SWNETSNOWTILE  )) deallocate(SWNETSNOWTILE  )
    if(associated(LWUPSNOWTILE  )) deallocate(LWUPSNOWTILE  )
    if(associated(LWDNSNOWTILE  )) deallocate(LWDNSNOWTILE  )
    if(associated(TCSORIGTILE  )) deallocate(TCSORIGTILE  )
    if(associated(TPSN1INTILE  )) deallocate(TPSN1INTILE  )
    if(associated(TPSN1OUTTILE  )) deallocate(TPSN1OUTTILE  )
    if(associated(GHSNOWTILE  )) deallocate(GHSNOWTILE  )
    if(associated(LHLANDTILE  )) deallocate(LHLANDTILE  )
    if(associated(SHLANDTILE  )) deallocate(SHLANDTILE  )
    if(associated(SWLANDTILE  )) deallocate(SWLANDTILE  )
    if(associated(SWDOWNLANDTILE  )) deallocate(SWDOWNLANDTILE  )
    if(associated(LWLANDTILE  )) deallocate(LWLANDTILE  )
    if(associated(GHLANDTILE  )) deallocate(GHLANDTILE  )
    if(associated(GHTSKINTILE )) deallocate(GHTSKINTILE )
    if(associated(SMLANDTILE  )) deallocate(SMLANDTILE  )
    if(associated(QINFILTILE  )) deallocate(QINFILTILE  )
    if(associated(TWLANDTILE  )) deallocate(TWLANDTILE  )
    if(associated(TELANDTILE  )) deallocate(TELANDTILE  )
    if(associated(TSLANDTILE  )) deallocate(TSLANDTILE  )
    if(associated(DWLANDTILE  )) deallocate(DWLANDTILE  )
    if(associated(DHLANDTILE  )) deallocate(DHLANDTILE  )
    if(associated(SPLANDTILE  )) deallocate(SPLANDTILE  )
    if(associated(SPLHTILE    )) deallocate(SPLHTILE    )
    if(associated(SPWATRTILE  )) deallocate(SPWATRTILE  )
    if(associated(SPSNOWTILE  )) deallocate(SPSNOWTILE  )
    if(associated(RDU001TILE  )) deallocate(RDU001TILE  )
    if(associated(RDU002TILE  )) deallocate(RDU002TILE  )
    if(associated(RDU003TILE  )) deallocate(RDU003TILE  )
    if(associated(RDU004TILE  )) deallocate(RDU004TILE  )
    if(associated(RDU005TILE  )) deallocate(RDU005TILE  )
    if(associated(RBC001TILE  )) deallocate(RBC001TILE  )
    if(associated(RBC002TILE  )) deallocate(RBC002TILE  )
    if(associated(ROC001TILE  )) deallocate(ROC001TILE  )
    if(associated(ROC002TILE  )) deallocate(ROC002TILE  )
    if(associated(RMELTDU001TILE )) deallocate(RMELTDU001TILE )
    if(associated(RMELTDU002TILE )) deallocate(RMELTDU002TILE )
    if(associated(RMELTDU003TILE )) deallocate(RMELTDU003TILE )
    if(associated(RMELTDU004TILE )) deallocate(RMELTDU004TILE )
    if(associated(RMELTDU005TILE )) deallocate(RMELTDU005TILE )
    if(associated(RMELTBC001TILE )) deallocate(RMELTBC001TILE )
    if(associated(RMELTBC002TILE )) deallocate(RMELTBC002TILE )
    if(associated(RMELTOC001TILE )) deallocate(RMELTOC001TILE )
    if(associated(RMELTOC002TILE )) deallocate(RMELTOC002TILE )
    if(associated(PEATCLSM_WATERLEVELTILE)) deallocate(PEATCLSM_WATERLEVELTILE)
    if(associated(PEATCLSM_FSWCHANGETILE )) deallocate(PEATCLSM_FSWCHANGETILE )
    if(associated(DZGT1TILE   )) deallocate(DZGT1TILE   )
    if(associated(DZGT2TILE   )) deallocate(DZGT2TILE   )
    if(associated(DZGT3TILE   )) deallocate(DZGT3TILE   )
    if(associated(DZGT4TILE   )) deallocate(DZGT4TILE   )
    if(associated(DZGT5TILE   )) deallocate(DZGT5TILE   )
    if(associated(DZGT6TILE   )) deallocate(DZGT6TILE   )
    if(associated(DZPRTILE    )) deallocate(DZPRTILE    )
    if(associated(DZRZTILE    )) deallocate(DZRZTILE    )
    if(associated(DZSFTILE    )) deallocate(DZSFTILE    )
    if(associated(DZTSTILE    )) deallocate(DZTSTILE    )
    if(associated(WPWETTILE   )) deallocate(WPWETTILE   )
    if(associated(WPEMWTILE   )) deallocate(WPEMWTILE   )
    if(associated(WPMCTILE    )) deallocate(WPMCTILE    )
    if(associated(CDCR2TILE   )) deallocate(CDCR2TILE   )
    if(associated(POROSTILE   )) deallocate(POROSTILE   )

    if(associated(CNLAITILE   )) deallocate(CNLAITILE   )
    if(associated(CNTLAITILE  )) deallocate(CNTLAITILE  )
    if(associated(CNSAITILE   )) deallocate(CNSAITILE   )
    if(associated(CNTOTCTILE  )) deallocate(CNTOTCTILE  )
    if(associated(CNVEGCTILE  )) deallocate(CNVEGCTILE  )
    if(associated(CNROOTTILE  )) deallocate(CNROOTTILE  )
    if(associated(CNFROOTCTILE)) deallocate(CNFROOTCTILE)
    if(associated(CNNPPTILE   )) deallocate(CNNPPTILE   )
    if(associated(CNGPPTILE   )) deallocate(CNGPPTILE   )
    if(associated(CNSRTILE    )) deallocate(CNSRTILE    )
    if(associated(CNNEETILE   )) deallocate(CNNEETILE   )
    if(associated(CNXSMRTILE  )) deallocate(CNXSMRTILE  )
    if(associated(CNADDTILE   )) deallocate(CNADDTILE   )
    if(associated(CNLOSSTILE  )) deallocate(CNLOSSTILE  )
    if(associated(CNBURNTILE  )) deallocate(CNBURNTILE  )
    if(associated(PARABSTILE  )) deallocate(PARABSTILE  )
    if(associated(PARINCTILE  )) deallocate(PARINCTILE  )
    if(associated(SCSATTILE   )) deallocate(SCSATTILE   )
    if(associated(SCUNSTILE   )) deallocate(SCUNSTILE   )
    if(associated(BTRANTTILE  )) deallocate(BTRANTTILE  )
    if(associated(SIFTILE     )) deallocate(SIFTILE     )
    if(associated(CNFSELTILE  )) deallocate(CNFSELTILE  )
    if(associated(ALBNFTILE ))  deallocate(ALBNFTILE)
    if(associated(ALBNRTILE ))  deallocate(ALBNRTILE)
    if(associated(ALBVFTILE ))  deallocate(ALBVFTILE)
    if(associated(ALBVRTILE ))  deallocate(ALBVRTILE)
    if(associated(EMISSTILE ))  deallocate(EMISSTILE)
    if(associated(FRTILE    ))  deallocate(FRTILE   )
    if(associated(OFRTILE   ))  deallocate(OFRTILE  )

    if(associated(DUDP)) deallocate( DUDP )
    if(associated(DUWT)) deallocate( DUWT )
    if(associated(DUSD)) deallocate( DUSD )
    if(associated(BCDP)) deallocate( BCDP )
    if(associated(BCSV)) deallocate( BCSV )
    if(associated(BCWT)) deallocate( BCWT )
    if(associated(BCSD)) deallocate( BCSD )
    if(associated(OCDP)) deallocate( OCDP )
    if(associated(OCSV)) deallocate( OCSV )
    if(associated(OCWT)) deallocate( OCWT )
    if(associated(OCSD)) deallocate( OCSD )
    if(associated(SSDP)) deallocate( SSDP )
    if(associated(SSSV)) deallocate( SSSV )
    if(associated(SSWT)) deallocate( SSWT )
    if(associated(SSSD)) deallocate( SSSD )

    if(associated( UUATILE  )) deallocate( UUATILE   )
    if(associated( VVATILE  )) deallocate( VVATILE   )

    call MAPL_TimerOff(MAPL,"-RUN2" )
    call MAPL_TimerOff(MAPL,"TOTAL")

!  All done
!-----------

    RETURN_(ESMF_SUCCESS)

  contains


    subroutine DOTYPE(type,RC)
      integer,           intent(IN ) :: TYPE
      integer, optional, intent(OUT) :: RC

!  Locals

      character(len=ESMF_MAXSTR)   :: IAm
      integer                      :: STATUS
      real, pointer                :: PTR1(:)
      type (MAPL_LocStreamXFORM)   :: XFORM
      real, pointer                :: DUM(:)

      call MAPL_TimerOn(MAPL,           trim(GCNames(type)))
      call MAPL_TimerOn(MAPL,"--RUN2_"//trim(GCNames(type)))

      Iam = trim(COMP_NAME)//"RUN2_DOTYPE"

! Fill the child's import state on his location stream from
! variables on Surface's location stream.
!----------------------------------------------------------

      XFORM = surf_internal_state%xform_in(type)

      call FILLIN_TILE(GIM(type), 'PS',     PSTILE,  XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DZ',     DZTILE,  XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'UU',     UUTILE,  XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'EVAP',   EVPTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'SH',     SHFTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DEVAP',  DEVTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DSH',    DSHTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'SNO',    SNOFLTILE,XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'ICE',    ICEFLTILE,XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'FRZR',   FRZRFLTILE,XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'PCU',    PCUTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'PLS',    PLSTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'TAUX',   TAUXTILE,XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'TAUY',   TAUYTILE,XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'UUA',    UUATILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'VVA',    VVATILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DRPAR',  DRPTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DFPAR',  DFPTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DRNIR',  DRNTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DFNIR',  DFNTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DRUVR',  DRUTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DFUVR',  DFUTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'LWDNSRF',LWBTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'ALW',    ALWTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'BLW',    BLWTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DTSDT', DTSDTTILE,XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DUDP', DUDPTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DUSV', DUSVTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DUWT', DUWTTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DUSD', DUSDTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'BCDP', BCDPTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'BCSV', BCSVTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'BCWT', BCWTTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'BCSD', BCSDTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'OCDP', OCDPTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'OCSV', OCSVTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'OCWT', OCWTTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'OCSD', OCSDTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'SUDP', SUDPTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'SUSV', SUSVTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'SUWT', SUWTTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'SUSD', SUSDTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'SSDP', SSDPTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'SSSV', SSSVTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'SSWT', SSWTTILE, XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'SSSD', SSSDTILE, XFORM, RC=STATUS); VERIFY_(STATUS)

      if (associated(DISCHARGETILE)) then
         call FILLIN_TILE(GIM(type), 'DISCHARGE',  DISCHARGETILE,  XFORM, RC=STATUS); VERIFY_(STATUS)
      end if

      call FILLIN_TILE(GIM(type), 'ZTH',    ZTHTILE, XFORM, RC=STATUS); VERIFY_(STATUS)

      call FILLIN_TILE(GIM(type), 'THATM',  THTILE,  XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'QHATM',  QHTILE,  XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'UHATM',  UHTILE,  XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'VHATM',  VHTILE,  XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'CTATM',  CTTILE,  XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'CQATM',  CQTILE,  XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'CMATM',  CMTILE,  XFORM, RC=STATUS); VERIFY_(STATUS)

! Force allocation of the exports needed from the child by the FILLOUTs below
!----------------------------------------------------------------------------

! Some cannot be verified, because some children don't produce them

      call MAPL_GetPointer(GEX(type), dum, 'LST'     , ALLOC=associated(LSTTILE     ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TP1'     , ALLOC=associated(TSOIL1TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TP2'     , ALLOC=associated(TSOIL2TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TP3'     , ALLOC=associated(TSOIL3TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TP4'     , ALLOC=associated(TSOIL4TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TP5'     , ALLOC=associated(TSOIL5TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TP6'     , ALLOC=associated(TSOIL6TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'ASNOW'   , ALLOC=associated(ASNOWTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SHSNOW'  , ALLOC=associated(SHSNOWTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'AVETSNOW', ALLOC=associated(AVETSNOWTILE), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TPSNOW'  , ALLOC=associated(TPSNOTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TPUNST'  , ALLOC=associated(TPUSTTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TPSAT'   , ALLOC=associated(TPSATTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TPWLT'   , ALLOC=associated(TPWLTTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TPSURF'  , ALLOC=associated(TPSURFTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'FRSAT'   , ALLOC=associated(FRSATTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'FRUST'   , ALLOC=associated(FRUSTTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'FRWLT'   , ALLOC=associated(FRWLTTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SNOWMASS', ALLOC=associated(SNOWTILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SNOWDP'  , ALLOC=associated(SNODTILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WET1'    , ALLOC=associated(WET1TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WET2'    , ALLOC=associated(WET2TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WET3'    , ALLOC=associated(WET3TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WCSF'    , ALLOC=associated(WCSFTILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WCRZ'    , ALLOC=associated(WCRZTILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WCPR'    , ALLOC=associated(WCPRTILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WESNN1'  , ALLOC=associated(WESNN1TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WESNN2'  , ALLOC=associated(WESNN2TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WESNN3'  , ALLOC=associated(WESNN3TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'CAPAC'   , ALLOC=associated(CAPACTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RUNOFF'  , ALLOC=associated(RUNOFFTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RUNSURF' , ALLOC=associated(RUNSURFTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'BASEFLOW', ALLOC=associated(BASEFLOWTILE), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'ACCUM'   , ALLOC=associated(ACCUMTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SMELT'   , ALLOC=associated(SMELTTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'EVPVEG'  , ALLOC=associated(EVEGTILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'EVPINT'  , ALLOC=associated(EINTTILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'EVPICE'  , ALLOC=associated(EICETILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'EVPSNO'  , ALLOC=associated(ESNOTILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'EVPSOI'  , ALLOC=associated(ESOITILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WAT10CM' , ALLOC=associated(WAT10CMTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WATSOI'  , ALLOC=associated(WATSOITILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'ICESOI'  , ALLOC=associated(ICESOITILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'EVLAND'  , ALLOC=associated(EVLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'PRLAND'  , ALLOC=associated(PRLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SNOLAND'  , ALLOC=associated(SNOLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DRPARLAND'  , ALLOC=associated(DRPARLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DFPARLAND'  , ALLOC=associated(DFPARLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'LHSNOW'  , ALLOC=associated(LHSNOWTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SWNETSNOW'  , ALLOC=associated(SWNETSNOWTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'LWUPSNOW'  , ALLOC=associated(LWUPSNOWTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'LWDNSNOW'  , ALLOC=associated(LWDNSNOWTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TCSORIG'  , ALLOC=associated(TCSORIGTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TPSN1IN'  , ALLOC=associated(TPSN1INTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TPSN1OUT'  , ALLOC=associated(TPSN1OUTTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'GHSNOW'  , ALLOC=associated(GHSNOWTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'LHLAND'  , ALLOC=associated(LHLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SHLAND'  , ALLOC=associated(SHLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SWLAND'  , ALLOC=associated(SWLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SWDOWNLAND'  , ALLOC=associated(SWDOWNLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'LWLAND'  , ALLOC=associated(LWLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'GHLAND'  , ALLOC=associated(GHLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'GHTSKIN' , ALLOC=associated(GHTSKINTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SMLAND'  , ALLOC=associated(SMLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'QINFIL'  , ALLOC=associated(QINFILTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TWLAND'  , ALLOC=associated(TWLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TELAND'  , ALLOC=associated(TELANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TSLAND'  , ALLOC=associated(TSLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DWLAND'  , ALLOC=associated(DWLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DHLAND'  , ALLOC=associated(DHLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SPLAND'  , ALLOC=associated(SPLANDTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SPLH'    , ALLOC=associated(SPLHTILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SPWATR'  , ALLOC=associated(SPWATRTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SPSNOW'  , ALLOC=associated(SPSNOWTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'FRACI'   , ALLOC=associated(   FRTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
! in case FRACI removed in future
      call MAPL_GetPointer(GEX(type), dum, 'FRACI'   , ALLOC=associated(  OFRTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RDU001'  , ALLOC=associated(RDU001TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RDU002'  , ALLOC=associated(RDU002TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RDU003'  , ALLOC=associated(RDU003TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RDU004'  , ALLOC=associated(RDU004TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RDU005'  , ALLOC=associated(RDU005TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RBC001'  , ALLOC=associated(RBC001TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RBC002'  , ALLOC=associated(RBC002TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'ROC001'  , ALLOC=associated(ROC001TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'ROC002'  , ALLOC=associated(ROC002TILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RMELTDU001' , ALLOC=associated(RMELTDU001TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RMELTDU002' , ALLOC=associated(RMELTDU002TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RMELTDU003' , ALLOC=associated(RMELTDU003TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RMELTDU004' , ALLOC=associated(RMELTDU004TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RMELTDU005' , ALLOC=associated(RMELTDU005TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RMELTBC001' , ALLOC=associated(RMELTBC001TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RMELTBC002' , ALLOC=associated(RMELTBC002TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RMELTOC001' , ALLOC=associated(RMELTOC001TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RMELTOC002' , ALLOC=associated(RMELTOC002TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'PEATCLSM_WATERLEVEL', ALLOC=associated(PEATCLSM_WATERLEVELTILE), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'PEATCLSM_FSWCHANGE' , ALLOC=associated(PEATCLSM_FSWCHANGETILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DZGT1' , ALLOC=associated(DZGT1TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DZGT2' , ALLOC=associated(DZGT2TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DZGT3' , ALLOC=associated(DZGT3TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DZGT4' , ALLOC=associated(DZGT4TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DZGT5' , ALLOC=associated(DZGT5TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DZGT6' , ALLOC=associated(DZGT6TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DZPR'  , ALLOC=associated(DZPRTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DZRZ'  , ALLOC=associated(DZRZTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DZSF'  , ALLOC=associated(DZSFTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DZTS'  , ALLOC=associated(DZTSTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WPWET' , ALLOC=associated(WPWETTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WPEMW' , ALLOC=associated(WPEMWTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'WPMC'  , ALLOC=associated(WPMCTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'CDCR2' , ALLOC=associated(CDCR2TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'POROS' , ALLOC=associated(POROSTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)


      IF (LSM_CHOICE > 1) THEN
         call MAPL_GetPointer(GEX(type), dum, 'CNLAI'   , ALLOC=associated(CNLAITILE   ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNTLAI'  , ALLOC=associated(CNTLAITILE  ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNSAI'   , ALLOC=associated(CNSAITILE   ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNTOTC'  , ALLOC=associated(CNTOTCTILE  ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNVEGC'  , ALLOC=associated(CNVEGCTILE  ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNROOT'  , ALLOC=associated(CNROOTTILE  ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         if (LSM_CHOICE == 3) then
	    call MAPL_GetPointer(GEX(type), dum, 'CNFROOTC' , ALLOC=associated(CNFROOTCTILE), notFoundOK=.true., RC=STATUS)
         endif
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNNPP'   , ALLOC=associated(CNNPPTILE   ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNGPP'   , ALLOC=associated(CNGPPTILE   ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNSR'    , ALLOC=associated(CNSRTILE    ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNNEE'   , ALLOC=associated(CNNEETILE   ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNXSMR'  , ALLOC=associated(CNXSMRTILE  ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNADD'   , ALLOC=associated(CNADDTILE   ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNLOSS'  , ALLOC=associated(CNLOSSTILE  ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNBURN'  , ALLOC=associated(CNBURNTILE  ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'PARABS'  , ALLOC=associated(PARABSTILE  ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'PARINC'  , ALLOC=associated(PARINCTILE  ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'SCSAT'   , ALLOC=associated(SCSATTILE   ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'SCUNS'   , ALLOC=associated(SCUNSTILE   ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'BTRANT'  , ALLOC=associated(BTRANTTILE  ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'SIF'     , ALLOC=associated(SIFTILE     ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'CNFSEL'  , ALLOC=associated(CNFSELTILE  ), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
      END IF
      call MAPL_GetPointer(GEX(type), dum, 'HLATWTR' , ALLOC=associated(HLATWTRTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'HLATICE' , ALLOC=associated(HLATICETILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SHWTR'   , ALLOC=associated(  SHWTRTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SHICE'   , ALLOC=associated(  SHICETILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TAUXW'   , ALLOC=associated(  TAUXWTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TAUXI'   , ALLOC=associated(  TAUXITILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TAUYW'   , ALLOC=associated(  TAUYWTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TAUYI'   , ALLOC=associated(  TAUYITILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'LWNDWTR' , ALLOC=associated(LWNDWTRTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SWNDWTR' , ALLOC=associated(SWNDWTRTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'LWNDICE' , ALLOC=associated(LWNDICETILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SWNDICE' , ALLOC=associated(SWNDICETILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'RAINOCN' , ALLOC=associated(RAINOCNTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SNOWOCN' , ALLOC=associated(SNOWOCNTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'ICEFOCN' , ALLOC=associated(ICEFOCNTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SPTOTOCN', ALLOC=associated(SPTOTOCNTILE), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TSKINW', ALLOC=associated(TSKINWTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TSKINICE', ALLOC=associated(TSKINICETILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_GetPointer(GEX(type), dum, 'DCOOL' ,   ALLOC=associated(DCOOL_TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DWARM' ,   ALLOC=associated(DWARM_TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TDROP' ,   ALLOC=associated(TDROP_TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'QCOOL' ,   ALLOC=associated(QCOOL_TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SWCOOL',   ALLOC=associated(SWCOOL_TILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'USTARW' ,  ALLOC=associated(USTARW_TILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TBAR'  ,   ALLOC=associated(TBAR_TILE     ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'LCOOL' ,   ALLOC=associated(LCOOL_TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'BCOOL' ,   ALLOC=associated(BCOOL_TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TDEL'    , ALLOC=associated(TDEL_TILE     ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TS_FOUND', ALLOC=associated(TS_FOUND_TILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'QWARM',    ALLOC=associated(QWARM_TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SWWARM',   ALLOC=associated(SWWARM_TILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'LANGM',    ALLOC=associated(LANGM_TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'PHIW',     ALLOC=associated(PHIW_TILE     ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TAUTW',    ALLOC=associated(TAUTW_TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'ZETA_W',   ALLOC=associated(ZETA_W_TILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TWMTF',    ALLOC=associated(TWMTF_TILE    ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_GetPointer(GEX(type), dum, 'HICE'  , ALLOC=associated(HICETILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'HSNO'  , ALLOC=associated(HSNOTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'FRZMLT', ALLOC=associated(FRZMLTTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'TSKINWCICE', ALLOC=associated(TSKINWCICETILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'ISTSFC', ALLOC=associated(ISTSFCTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SSKINW', ALLOC=associated(SSKINWTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'MELTT', ALLOC=associated(MELTTTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'MELTB', ALLOC=associated(MELTBTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'MELTS', ALLOC=associated(MELTSTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'MELTL', ALLOC=associated(MELTLTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'FRAZIL',ALLOC=associated(FRAZILTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'CONGEL',ALLOC=associated(CONGELTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SNOICE',ALLOC=associated(SNOICETILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DAIDTT',ALLOC=associated(DAIDTTTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DVIDTT',ALLOC=associated(DVIDTTTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DAIDTD',ALLOC=associated(DAIDTDTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DVIDTD',ALLOC=associated(DVIDTDTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'FBOT'  ,ALLOC=associated(FBOTTILE   ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'HFLUX' ,ALLOC=associated(HFLUXTILE  ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum,'WATERFLUX',ALLOC=associated(WFLUXTILE), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SALTFLUX',ALLOC=associated(SFLUXTILE), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'FSWTHRU',ALLOC=associated(FSWTHRUTILE), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'FSWABS',ALLOC=associated(FSWABSTILE), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'USTARI',ALLOC=associated(USTARITILE), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'FHOCN' ,ALLOC=associated(FHOCNTILE ), notFoundOK=.true., RC=STATUS)
      VERIFY_(STATUS)

      if (DO_FIRE_DANGER) then
         call MAPL_GetPointer(GEX(type), dum, 'FFMC', ALLOC=associated(FFMCTILE), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'GFMC', ALLOC=associated(GFMCTILE), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'DMC',  ALLOC=associated(DMCTILE),  notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'DC',   ALLOC=associated(DCTILE),   notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'ISI',  ALLOC=associated(ISITILE),  notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'BUI',  ALLOC=associated(BUITILE),  notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'FWI',  ALLOC=associated(FWITILE),  notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'DSR',  ALLOC=associated(DSRTILE),  notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)

         call MAPL_GetPointer(GEX(type), dum, 'FFMC_DAILY',  ALLOC=associated(FFMCDAILYTILE),  notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'DMC_DAILY',   ALLOC=associated(DMCDAILYTILE),   notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'DC_DAILY',    ALLOC=associated(DCDAILYTILE),    notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'ISI_DAILY',   ALLOC=associated(ISIDAILYTILE),   notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'BUI_DAILY',   ALLOC=associated(BUIDAILYTILE),   notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'FWI_DAILY',   ALLOC=associated(FWIDAILYTILE),   notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'DSR_DAILY',   ALLOC=associated(DSRDAILYTILE),   notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)

         call MAPL_GetPointer(GEX(type), dum, 'FFMC_DAILY_', ALLOC=associated(FFMCDAILYTILE_), notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'DMC_DAILY_',  ALLOC=associated(DMCDAILYTILE_),  notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'DC_DAILY_',   ALLOC=associated(DCDAILYTILE_),   notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'ISI_DAILY_',  ALLOC=associated(ISIDAILYTILE_),  notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'BUI_DAILY_',  ALLOC=associated(BUIDAILYTILE_),  notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'FWI_DAILY_',  ALLOC=associated(FWIDAILYTILE_),  notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(GEX(type), dum, 'DSR_DAILY_',  ALLOC=associated(DSRDAILYTILE_),  notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)

         call MAPL_GetPointer(GEX(type), dum, 'VPD',         ALLOC=associated(VPDTILE),        notFoundOK=.true., RC=STATUS)
         VERIFY_(STATUS)
      end if


! All children can produce these

      call MAPL_GetPointer(GEX(type), dum, 'DELTS'  , ALLOC=associated(DTSTILE)    , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'DELQS'  , ALLOC=associated(DQSTILE)    , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'HLATN'  , ALLOC=associated(HLATNTILE)  , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'EVAPOUT', ALLOC=associated(EVAPOUTILE) , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SUBLIM' , ALLOC=associated(SUBLIMTILE) , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SHOUT'  , ALLOC=associated(SHOUTILE)   , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'HLWUP'  , ALLOC=associated(HLWUPTILE)  , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'LWNDSRF', ALLOC=associated(LWNDSRFTILE), RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'SWNDSRF', ALLOC=associated(SWNDSRFTILE), RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'ALBVR'  , ALLOC=associated(ALBVRTILE)  , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'ALBVF'  , ALLOC=associated(ALBVFTILE)  , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'ALBNR'  , ALLOC=associated(ALBNRTILE)  , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'ALBNF'  , ALLOC=associated(ALBNFTILE)  , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(type), dum, 'EMIS'   , ALLOC=associated(EMISSTILE)  , RC=STATUS)
      VERIFY_(STATUS)


! Run the child
!--------------

      call ESMF_GridCompRun(GCS(type), &
           importState=GIM(type), exportState=GEX(type), &
           clock=CLOCK, PHASE=2, userRC=STATUS )
      VERIFY_(STATUS)

! Fill variables on Surface's location stream from the child's
!  export state, which is on his location stream.
!-------------------------------------------------------------

      XFORM = SURF_INTERNAL_STATE%XFORM_OUT(type)

      if(associated(TSTILE)) then
         call FILLOUT_TILE(GEX(type), 'TST',        TSTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(QSTILE)) then
         call FILLOUT_TILE(GEX(type), 'QST',        QSTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(DTSTILE)) then
         call FILLOUT_TILE(GEX(type), 'DELTS',     DTSTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(DQSTILE)) then
         call FILLOUT_TILE(GEX(type), 'DELQS',     DQSTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(ALBVRTILE)) then
         call FILLOUT_TILE(GEX(type), 'ALBVR',   ALBVRTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(ALBVFTILE)) then
         call FILLOUT_TILE(GEX(type), 'ALBVF',   ALBVFTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(ALBNRTILE)) then
         call FILLOUT_TILE(GEX(type), 'ALBNR',   ALBNRTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(ALBNFTILE)) then
         call FILLOUT_TILE(GEX(type), 'ALBNF',   ALBNFTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(EMISSTILE)) then
         call FILLOUT_TILE(GEX(type), 'EMIS',    EMISSTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(   FRTILE)) then
         call FILLOUT_TILE(GEX(type), 'FRACI',      FRTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(  OFRTILE)) then
         call FILLOUT_TILE(GEX(type), 'FRACI',     OFRTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TSOIL1TILE)) then
         call FILLOUT_TILE(GEX(type), 'TP1',    TSOIL1TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TSOIL2TILE)) then
         call FILLOUT_TILE(GEX(type), 'TP2',    TSOIL2TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TSOIL3TILE)) then
         call FILLOUT_TILE(GEX(type), 'TP3',    TSOIL3TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TSOIL4TILE)) then
         call FILLOUT_TILE(GEX(type), 'TP4',    TSOIL4TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TSOIL5TILE)) then
         call FILLOUT_TILE(GEX(type), 'TP5',    TSOIL5TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TSOIL6TILE)) then
         call FILLOUT_TILE(GEX(type), 'TP6',    TSOIL6TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SNOWTILE)) then
         call FILLOUT_TILE(GEX(type), 'SNOWMASS', SNOWTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(ASNOWTILE)) then
         call FILLOUT_TILE(GEX(type), 'ASNOW',   ASNOWTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SHSNOWTILE)) then
         call FILLOUT_TILE(GEX(type), 'SHSNOW',  SHSNOWTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(AVETSNOWTILE)) then
         call FILLOUT_TILE(GEX(type), 'AVETSNOW',  AVETSNOWTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TPSNOTILE)) then
         call FILLOUT_TILE(GEX(type), 'TPSNOW',  TPSNOTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TPUSTTILE)) then
         call FILLOUT_TILE(GEX(type), 'TPUNST',  TPUSTTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TPSATTILE)) then
         call FILLOUT_TILE(GEX(type), 'TPSAT' ,  TPSATTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TPWLTTILE)) then
         call FILLOUT_TILE(GEX(type), 'TPWLT' ,  TPWLTTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TPSURFTILE)) then
         call FILLOUT_TILE(GEX(type), 'TPSURF',  TPSURFTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(FRSATTILE)) then
         call FILLOUT_TILE(GEX(type), 'FRSAT',   FRSATTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(FRUSTTILE)) then
         call FILLOUT_TILE(GEX(type), 'FRUST',   FRUSTTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(FRWLTTILE)) then
         call FILLOUT_TILE(GEX(type), 'FRWLT',   FRWLTTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SNODTILE)) then
         call FILLOUT_TILE(GEX(type), 'SNOWDP',   SNODTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(WET1TILE)) then
         call FILLOUT_TILE(GEX(type), 'WET1',     WET1TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(WET2TILE)) then
         call FILLOUT_TILE(GEX(type), 'WET2',     WET2TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(WET3TILE)) then
         call FILLOUT_TILE(GEX(type), 'WET3',     WET3TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(WCSFTILE)) then
         call FILLOUT_TILE(GEX(type), 'WCSF',     WCSFTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(WCRZTILE)) then
         call FILLOUT_TILE(GEX(type), 'WCRZ',     WCRZTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(WCPRTILE)) then
         call FILLOUT_TILE(GEX(type), 'WCPR',     WCPRTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(WESNN1TILE)) then
         call FILLOUT_TILE(GEX(type), 'WESNN1', WESNN1TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(WESNN2TILE)) then
         call FILLOUT_TILE(GEX(type), 'WESNN2', WESNN2TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(WESNN3TILE)) then
         call FILLOUT_TILE(GEX(type), 'WESNN3', WESNN3TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CAPACTILE)) then
         call FILLOUT_TILE(GEX(type), 'CAPAC',   CAPACTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(RUNOFFTILE)) then
         call FILLOUT_TILE(GEX(type), 'RUNOFF', RUNOFFTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(RUNSURFTILE)) then
         call FILLOUT_TILE(GEX(type), 'RUNSURF',RUNSURFTILE,XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(BASEFLOWTILE)) then
         call FILLOUT_TILE(GEX(type), 'BASEFLOW',BASEFLOWTILE,XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(ACCUMTILE)) then
         call FILLOUT_TILE(GEX(type), 'ACCUM',   ACCUMTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SMELTTILE)) then
         call FILLOUT_TILE(GEX(type), 'SMELT',   SMELTTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(EVEGTILE)) then
         call FILLOUT_TILE(GEX(type), 'EVPVEG',   EVEGTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(EINTTILE)) then
         call FILLOUT_TILE(GEX(type), 'EVPINT',   EINTTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(ESOITILE)) then
         call FILLOUT_TILE(GEX(type), 'EVPSOI',   ESOITILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(EICETILE)) then
         call FILLOUT_TILE(GEX(type), 'EVPICE',   EICETILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(ESNOTILE)) then
         call FILLOUT_TILE(GEX(type), 'EVPSNO',   ESNOTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(WAT10CMTILE)) then
         call FILLOUT_TILE(GEX(type), 'WAT10CM',WAT10CMTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(WATSOITILE))  then
         call FILLOUT_TILE(GEX(type), 'WATSOI', WATSOITILE,  XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(ICESOITILE))  then
         call FILLOUT_TILE(GEX(type), 'ICESOI', ICESOITILE,  XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(EVLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'EVLAND',   EVLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(PRLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'PRLAND',   PRLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SNOLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'SNOLAND',   SNOLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(DRPARLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'DRPARLAND',   DRPARLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(DFPARLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'DFPARLAND',   DFPARLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
       if(associated(LHSNOWTILE)) then
          call FILLOUT_TILE(GEX(type), 'LHSNOW',   LHSNOWTILE, XFORM, RC=STATUS)
          VERIFY_(STATUS)
       end if
       if(associated(TCSORIGTILE)) then
          call FILLOUT_TILE(GEX(type), 'TCSORIG',   TCSORIGTILE, XFORM, RC=STATUS)
          VERIFY_(STATUS)
       end if
       if(associated(TPSN1INTILE)) then
          call FILLOUT_TILE(GEX(type), 'TPSN1IN',   TPSN1INTILE, XFORM, RC=STATUS)
          VERIFY_(STATUS)
       end if
       if(associated(TPSN1OUTTILE)) then
          call FILLOUT_TILE(GEX(type), 'TPSN1OUT',   TPSN1OUTTILE, XFORM, RC=STATUS)
          VERIFY_(STATUS)
       end if
       if(associated(SWNETSNOWTILE)) then
          call FILLOUT_TILE(GEX(type), 'SWNETSNOW',   SWNETSNOWTILE, XFORM, RC=STATUS)
          VERIFY_(STATUS)
       end if
       if(associated(LWDNSNOWTILE)) then
          call FILLOUT_TILE(GEX(type), 'LWDNSNOW',   LWDNSNOWTILE, XFORM, RC=STATUS)
          VERIFY_(STATUS)
       end if
       if(associated(LWUPSNOWTILE)) then
          call FILLOUT_TILE(GEX(type), 'LWUPSNOW',   LWUPSNOWTILE, XFORM, RC=STATUS)
          VERIFY_(STATUS)
       end if
      if(associated(LHLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'LHLAND',   LHLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SHLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'SHLAND',   SHLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SWLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'SWLAND',   SWLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SWDOWNLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'SWDOWNLAND',   SWDOWNLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(LWLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'LWLAND',   LWLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(GHLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'GHLAND',   GHLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(GHSNOWTILE)) then
         call FILLOUT_TILE(GEX(type), 'GHSNOW',   GHSNOWTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(GHTSKINTILE)) then
         call FILLOUT_TILE(GEX(type), 'GHTSKIN',  GHTSKINTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SMLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'SMLAND',   SMLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(QINFILTILE)) then
         call FILLOUT_TILE(GEX(type), 'QINFIL',   QINFILTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TWLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'TWLAND',   TWLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TELANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'TELAND',   TELANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TSLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'TSLAND',   TSLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(DWLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'DWLAND',   DWLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(DHLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'DHLAND',   DHLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SPLANDTILE)) then
         call FILLOUT_TILE(GEX(type), 'SPLAND',   SPLANDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SPLHTILE  )) then
         call FILLOUT_TILE(GEX(type), 'SPLH'  ,   SPLHTILE,   XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SPWATRTILE)) then
         call FILLOUT_TILE(GEX(type), 'SPWATR',   SPWATRTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SPSNOWTILE)) then
         call FILLOUT_TILE(GEX(type), 'SPSNOW',   SPSNOWTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(RDU001TILE)) call FILLOUT_TILE(GEX(type), 'RDU001' , RDU001TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RDU002TILE)) call FILLOUT_TILE(GEX(type), 'RDU002' , RDU002TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RDU003TILE)) call FILLOUT_TILE(GEX(type), 'RDU003' , RDU003TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RDU004TILE)) call FILLOUT_TILE(GEX(type), 'RDU004' , RDU004TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RDU005TILE)) call FILLOUT_TILE(GEX(type), 'RDU005' , RDU005TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RBC001TILE)) call FILLOUT_TILE(GEX(type), 'RBC001' , RBC001TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RBC002TILE)) call FILLOUT_TILE(GEX(type), 'RBC002' , RBC002TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(ROC001TILE)) call FILLOUT_TILE(GEX(type), 'ROC001' , ROC001TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(ROC002TILE)) call FILLOUT_TILE(GEX(type), 'ROC002' , ROC002TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RMELTDU001TILE)) call FILLOUT_TILE(GEX(type), 'RMELTDU001' , RMELTDU001TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RMELTDU002TILE)) call FILLOUT_TILE(GEX(type), 'RMELTDU002' , RMELTDU002TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RMELTDU003TILE)) call FILLOUT_TILE(GEX(type), 'RMELTDU003' , RMELTDU003TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RMELTDU004TILE)) call FILLOUT_TILE(GEX(type), 'RMELTDU004' , RMELTDU004TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RMELTDU005TILE)) call FILLOUT_TILE(GEX(type), 'RMELTDU005' , RMELTDU005TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RMELTBC001TILE)) call FILLOUT_TILE(GEX(type), 'RMELTBC001' , RMELTBC001TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RMELTBC002TILE)) call FILLOUT_TILE(GEX(type), 'RMELTBC002' , RMELTBC002TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RMELTOC001TILE)) call FILLOUT_TILE(GEX(type), 'RMELTOC001' , RMELTOC001TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(RMELTOC002TILE)) call FILLOUT_TILE(GEX(type), 'RMELTOC002' , RMELTOC002TILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(PEATCLSM_WATERLEVELTILE)) call FILLOUT_TILE(GEX(type), 'PEATCLSM_WATERLEVEL', PEATCLSM_WATERLEVELTILE, XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(PEATCLSM_FSWCHANGETILE))  call FILLOUT_TILE(GEX(type), 'PEATCLSM_FSWCHANGE' , PEATCLSM_FSWCHANGETILE , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(DZGT1TILE))  call FILLOUT_TILE(GEX(type), 'DZGT1'  , DZGT1TILE  , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(DZGT2TILE))  call FILLOUT_TILE(GEX(type), 'DZGT2'  , DZGT2TILE  , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(DZGT3TILE))  call FILLOUT_TILE(GEX(type), 'DZGT3'  , DZGT3TILE  , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(DZGT4TILE))  call FILLOUT_TILE(GEX(type), 'DZGT4'  , DZGT4TILE  , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(DZGT5TILE))  call FILLOUT_TILE(GEX(type), 'DZGT5'  , DZGT5TILE  , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(DZGT6TILE))  call FILLOUT_TILE(GEX(type), 'DZGT6'  , DZGT6TILE  , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(DZPRTILE ))  call FILLOUT_TILE(GEX(type), 'DZPR'   , DZPRTILE   , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(DZRZTILE ))  call FILLOUT_TILE(GEX(type), 'DZRZ'   , DZRZTILE   , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(DZSFTILE ))  call FILLOUT_TILE(GEX(type), 'DZSF'   , DZSFTILE   , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(DZTSTILE ))  call FILLOUT_TILE(GEX(type), 'DZTS'   , DZTSTILE   , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(WPWETTILE))  call FILLOUT_TILE(GEX(type), 'WPWET'  , WPWETTILE  , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(WPEMWTILE))  call FILLOUT_TILE(GEX(type), 'WPEMW'  , WPEMWTILE  , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(WPMCTILE ))  call FILLOUT_TILE(GEX(type), 'WPMC'   , WPMCTILE   , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(CDCR2TILE))  call FILLOUT_TILE(GEX(type), 'CDCR2'  , CDCR2TILE  , XFORM, RC=STATUS);VERIFY_(STATUS)
      if(associated(POROSTILE))  call FILLOUT_TILE(GEX(type), 'POROS'  , POROSTILE  , XFORM, RC=STATUS);VERIFY_(STATUS)

      if(associated(CNLAITILE)) then
         call FILLOUT_TILE(GEX(type), 'CNLAI' ,   CNLAITILE , XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNTLAITILE)) then
         call FILLOUT_TILE(GEX(type), 'CNTLAI',   CNTLAITILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNSAITILE)) then
         call FILLOUT_TILE(GEX(type), 'CNSAI' ,   CNSAITILE , XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNTOTCTILE)) then
         call FILLOUT_TILE(GEX(type), 'CNTOTC',   CNTOTCTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNVEGCTILE)) then
         call FILLOUT_TILE(GEX(type), 'CNVEGC',   CNVEGCTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNROOTTILE)) then
         call FILLOUT_TILE(GEX(type), 'CNROOT',   CNROOTTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNFROOTCTILE)) then
         call FILLOUT_TILE(GEX(type), 'CNFROOTC', CNFROOTCTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNNPPTILE)) then
         call FILLOUT_TILE(GEX(type), 'CNNPP' ,   CNNPPTILE , XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNGPPTILE)) then
         call FILLOUT_TILE(GEX(type), 'CNGPP' ,   CNGPPTILE , XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNSRTILE)) then
         call FILLOUT_TILE(GEX(type), 'CNSR'  ,   CNSRTILE  , XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNNEETILE)) then
         call FILLOUT_TILE(GEX(type), 'CNNEE' ,   CNNEETILE , XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNXSMRTILE)) then
         call FILLOUT_TILE(GEX(type), 'CNXSMR',   CNXSMRTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNADDTILE)) then
         call FILLOUT_TILE(GEX(type), 'CNADD' ,   CNADDTILE , XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNLOSSTILE)) then
         call FILLOUT_TILE(GEX(type), 'CNLOSS',   CNLOSSTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNBURNTILE)) then
         call FILLOUT_TILE(GEX(type), 'CNBURN',   CNBURNTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(PARABSTILE)) then
         call FILLOUT_TILE(GEX(type), 'PARABS',   PARABSTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(PARINCTILE)) then
         call FILLOUT_TILE(GEX(type), 'PARINC',   PARINCTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SCSATTILE)) then
         call FILLOUT_TILE(GEX(type), 'SCSAT' ,   SCSATTILE , XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SCUNSTILE)) then
         call FILLOUT_TILE(GEX(type), 'SCUNS' ,   SCUNSTILE , XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(BTRANTTILE)) then
         call FILLOUT_TILE(GEX(type), 'BTRANT',   BTRANTTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SIFTILE)) then
         call FILLOUT_TILE(GEX(type), 'SIF'   ,   SIFTILE   , XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CNFSELTILE)) then
         call FILLOUT_TILE(GEX(type), 'CNFSEL',   CNFSELTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(EVAPOUTILE)) then
         call FILLOUT_TILE(GEX(type), 'EVAPOUT',  EVAPOUTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SUBLIMTILE)) then
         call FILLOUT_TILE(GEX(type), 'SUBLIM',   SUBLIMTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SHOUTILE)) then
         call FILLOUT_TILE(GEX(type), 'SHOUT',   SHOUTILE,  XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(LSTTILE)) then
         call FILLOUT_TILE(GEX(type), 'LST',      LSTTILE,  XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(HLATNTILE)) then
         call FILLOUT_TILE(GEX(type), 'HLATN',   HLATNTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(HLATWTRTILE)) then
         call FILLOUT_TILE(GEX(type), 'HLATWTR',  HLATWTRTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(HLATICETILE)) then
         call FILLOUT_TILE(GEX(type), 'HLATICE',  HLATICETILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SHWTRTILE)  ) then
         call FILLOUT_TILE(GEX(type), 'SHWTR',      SHWTRTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SHICETILE)  ) then
         call FILLOUT_TILE(GEX(type), 'SHICE',      SHICETILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TAUXWTILE)  ) then
         call FILLOUT_TILE(GEX(type), 'TAUXW',      TAUXWTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TAUXITILE)  ) then
         call FILLOUT_TILE(GEX(type), 'TAUXI',      TAUXITILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TAUYWTILE)  ) then
         call FILLOUT_TILE(GEX(type), 'TAUYW',      TAUYWTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TAUYITILE)  ) then
         call FILLOUT_TILE(GEX(type), 'TAUYI',      TAUYITILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(LWNDWTRTILE)) then
         call FILLOUT_TILE(GEX(type), 'LWNDWTR',  LWNDWTRTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SWNDWTRTILE)) then
         call FILLOUT_TILE(GEX(type), 'SWNDWTR',  SWNDWTRTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(LWNDICETILE)) then
         call FILLOUT_TILE(GEX(type), 'LWNDICE',  LWNDICETILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SWNDICETILE)) then
         call FILLOUT_TILE(GEX(type), 'SWNDICE',  SWNDICETILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(RAINOCNTILE)) then
         call FILLOUT_TILE(GEX(type), 'RAINOCN',  RAINOCNTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SNOWOCNTILE)) then
         call FILLOUT_TILE(GEX(type), 'SNOWOCN',  SNOWOCNTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(ICEFOCNTILE)) then
         call FILLOUT_TILE(GEX(type), 'ICEFOCN',  ICEFOCNTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SPTOTOCNTILE)) then
         call FILLOUT_TILE(GEX(type), 'SPTOTOCN', SPTOTOCNTILE,XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(HLWUPTILE)) then
         call FILLOUT_TILE(GEX(type), 'HLWUP',   HLWUPTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(LWNDSRFTILE)) then
         call FILLOUT_TILE(GEX(type), 'LWNDSRF',LWNDSRFTILE,XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SWNDSRFTILE)) then
         call FILLOUT_TILE(GEX(type), 'SWNDSRF',SWNDSRFTILE,XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TSKINWTILE)) then
         call FILLOUT_TILE(GEX(type), 'TSKINW',TSKINWTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TSKINICETILE)) then
         call FILLOUT_TILE(GEX(type), 'TSKINICE',TSKINICETILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(DCOOL_TILE)) then
         call FILLOUT_TILE(GEX(type), 'DCOOL', DCOOL_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(DWARM_TILE)) then
         call FILLOUT_TILE(GEX(type), 'DWARM', DWARM_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TDROP_TILE)) then
         call FILLOUT_TILE(GEX(type), 'TDROP', TDROP_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(QCOOL_TILE)) then
         call FILLOUT_TILE(GEX(type), 'QCOOL', QCOOL_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SWCOOL_TILE)) then
         call FILLOUT_TILE(GEX(type), 'SWCOOL',SWCOOL_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(USTARW_TILE)) then
         call FILLOUT_TILE(GEX(type), 'USTARW', USTARW_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TBAR_TILE)) then
         call FILLOUT_TILE(GEX(type), 'TBAR',  TBAR_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(LCOOL_TILE)) then
         call FILLOUT_TILE(GEX(type), 'LCOOL', LCOOL_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(BCOOL_TILE)) then
         call FILLOUT_TILE(GEX(type), 'BCOOL', BCOOL_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(TDEL_TILE)) then
         call FILLOUT_TILE(GEX(type), 'TDEL', TDEL_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(TS_FOUND_TILE)) then
         call FILLOUT_TILE(GEX(type), 'TS_FOUND', TS_FOUND_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(SS_FOUND_TILE)) then
         call FILLOUT_TILE(GEX(type), 'SS_FOUND', SS_FOUND_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(QWARM_TILE)) then
         call FILLOUT_TILE(GEX(type), 'QWARM', QWARM_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(SWWARM_TILE)) then
         call FILLOUT_TILE(GEX(type), 'SWWARM', SWWARM_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(LANGM_TILE)) then
         call FILLOUT_TILE(GEX(type), 'LANGM', LANGM_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(PHIW_TILE)) then
         call FILLOUT_TILE(GEX(type), 'PHIW', PHIW_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(TAUTW_TILE)) then
         call FILLOUT_TILE(GEX(type), 'TAUTW', TAUTW_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(ZETA_W_TILE)) then
         call FILLOUT_TILE(GEX(type), 'ZETA_W', ZETA_W_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(TWMTF_TILE)) then
         call FILLOUT_TILE(GEX(type), 'TWMTF', TWMTF_TILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if(associated(HICETILE)) then
         call FILLOUT_TILE(GEX(type), 'HICE',  HICETILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(HSNOTILE)) then
         call FILLOUT_TILE(GEX(type), 'HSNO',  HSNOTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(FRZMLTTILE)) then
         call FILLOUT_TILE(GEX(type), 'FRZMLT',FRZMLTTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(TSKINWCICETILE)) then
         call FILLOUT_TILE(GEX(type), 'TSKINWCICE',TSKINWCICETILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(ISTSFCTILE)) then
         call FILLOUT_TILE(GEX(type), 'ISTSFC',ISTSFCTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SSKINWTILE)) then
         call FILLOUT_TILE(GEX(type), 'SSKINW',SSKINWTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(MELTTTILE)) then
         call FILLOUT_TILE(GEX(type), 'MELTT',MELTTTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(MELTBTILE)) then
         call FILLOUT_TILE(GEX(type), 'MELTB',MELTBTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(MELTSTILE)) then
         call FILLOUT_TILE(GEX(type), 'MELTS',MELTSTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(MELTLTILE)) then
         call FILLOUT_TILE(GEX(type), 'MELTL',MELTLTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(FRAZILTILE)) then
         call FILLOUT_TILE(GEX(type), 'FRAZIL',FRAZILTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(CONGELTILE)) then
         call FILLOUT_TILE(GEX(type), 'CONGEL',CONGELTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SNOICETILE)) then
         call FILLOUT_TILE(GEX(type), 'SNOICE',SNOICETILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(DAIDTTTILE)) then
         call FILLOUT_TILE(GEX(type), 'DAIDTT',DAIDTTTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(DVIDTTTILE)) then
         call FILLOUT_TILE(GEX(type), 'DVIDTT',DVIDTTTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(DAIDTDTILE)) then
         call FILLOUT_TILE(GEX(type), 'DAIDTD',DAIDTDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(DVIDTDTILE)) then
         call FILLOUT_TILE(GEX(type), 'DVIDTD',DVIDTDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(FBOTTILE)) then
         call FILLOUT_TILE(GEX(type), 'FBOT'  ,FBOTTILE,   XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(HFLUXTILE)) then
         call FILLOUT_TILE(GEX(type), 'HFLUX' ,HFLUXTILE,   XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(WFLUXTILE)) then
         call FILLOUT_TILE(GEX(type), 'WATERFLUX',WFLUXTILE,   XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(SFLUXTILE)) then
         call FILLOUT_TILE(GEX(type), 'SALTFLUX',SFLUXTILE,   XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(FSWTHRUTILE)) then
         call FILLOUT_TILE(GEX(type), 'FSWTHRU' ,FSWTHRUTILE,   XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(FSWABSTILE)) then
         call FILLOUT_TILE(GEX(type), 'FSWABS' ,FSWABSTILE,   XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(USTARITILE)) then
         call FILLOUT_TILE(GEX(type), 'USTARI' ,USTARITILE,   XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(FHOCNTILE)) then
         call FILLOUT_TILE(GEX(type), 'FHOCN'  ,FHOCNTILE,   XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if(associated(DISCHARGETILE)) then
         call FILLOUT_TILE(GEX(type), 'DISCHARGE'  ,DISCHARGETILE,   XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

! Fire danger
      if (associated(FFMCTILE)) then
         call FILLOUT_TILE(GEX(type), 'FFMC', FFMCTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(GFMCTILE)) then
         call FILLOUT_TILE(GEX(type), 'GFMC', GFMCTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(DMCTILE)) then
         call FILLOUT_TILE(GEX(type), 'DMC', DMCTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(DCTILE)) then
         call FILLOUT_TILE(GEX(type), 'DC', DCTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(ISITILE)) then
         call FILLOUT_TILE(GEX(type), 'ISI', ISITILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(BUITILE)) then
         call FILLOUT_TILE(GEX(type), 'BUI', BUITILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(FWITILE)) then
         call FILLOUT_TILE(GEX(type), 'FWI', FWITILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(DSRTILE)) then
         call FILLOUT_TILE(GEX(type), 'DSR', DSRTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if (associated(FFMCDAILYTILE)) then
         call FILLOUT_TILE(GEX(type), 'FFMC_DAILY', FFMCDAILYTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(DMCDAILYTILE)) then
         call FILLOUT_TILE(GEX(type), 'DMC_DAILY', DMCDAILYTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(DCDAILYTILE)) then
         call FILLOUT_TILE(GEX(type), 'DC_DAILY', DCDAILYTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(ISIDAILYTILE)) then
         call FILLOUT_TILE(GEX(type), 'ISI_DAILY', ISIDAILYTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(BUIDAILYTILE)) then
         call FILLOUT_TILE(GEX(type), 'BUI_DAILY', BUIDAILYTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(FWIDAILYTILE)) then
         call FILLOUT_TILE(GEX(type), 'FWI_DAILY', FWIDAILYTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(DSRDAILYTILE)) then
         call FILLOUT_TILE(GEX(type), 'DSR_DAILY', DSRDAILYTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if

      if (associated(FFMCDAILYTILE_)) then
         call FILLOUT_TILE(GEX(type), 'FFMC_DAILY_', FFMCDAILYTILE_, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(DMCDAILYTILE_)) then
         call FILLOUT_TILE(GEX(type), 'DMC_DAILY_', DMCDAILYTILE_, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(DCDAILYTILE_)) then
         call FILLOUT_TILE(GEX(type), 'DC_DAILY_', DCDAILYTILE_, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(ISIDAILYTILE_)) then
         call FILLOUT_TILE(GEX(type), 'ISI_DAILY_', ISIDAILYTILE_, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(BUIDAILYTILE_)) then
         call FILLOUT_TILE(GEX(type), 'BUI_DAILY_', BUIDAILYTILE_, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(FWIDAILYTILE_)) then
         call FILLOUT_TILE(GEX(type), 'FWI_DAILY_', FWIDAILYTILE_, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(DSRDAILYTILE_)) then
         call FILLOUT_TILE(GEX(type), 'DSR_DAILY_', DSRDAILYTILE_, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if
      if (associated(VPDTILE)) then
         call FILLOUT_TILE(GEX(type), 'VPD', VPDTILE, XFORM, RC=STATUS)
         VERIFY_(STATUS)
      end if


      call MAPL_TimerOff(MAPL,"--RUN2_"//trim(GCNames(type)))
      call MAPL_TimerOff(MAPL,           trim(GCNames(type)))

      RETURN_(ESMF_SUCCESS)

    end subroutine DOTYPE

  end subroutine RUN2

  subroutine MKTILE_1D(VAR, TILEVAR, NT, RC)
    real, pointer                  :: VAR(:,:)
    real, pointer                  :: TILEVAR(:)
    integer,           intent(IN)  :: NT
    integer, optional, intent(OUT) :: RC

    character(len=ESMF_MAXSTR)   :: IAm='MKTILE_1D'
    integer                      :: STATUS

    if(associated(VAR) .and. .not.associated(TILEVAR)) then
       allocate(TILEVAR(NT), STAT=STATUS)
       VERIFY_(STATUS)
       TILEVAR = MAPL_Undef
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MKTILE_1D

  subroutine MKTILE_UNGRIDDED(VAR, TILEVAR, NT, RC)
    real, pointer                  :: VAR(:,:,:)
    real, pointer                  :: TILEVAR(:,:)
    integer,           intent(IN)  :: NT
    integer, optional, intent(OUT) :: RC

    character(len=ESMF_MAXSTR)   :: IAm='MKTILE_UNGRIDDED'
    integer                      :: STATUS

    if(associated(VAR) .and. .not.associated(TILEVAR)) then
       allocate(TILEVAR(NT,size(VAR,3)), STAT=STATUS)
       VERIFY_(STATUS)
       TILEVAR = MAPL_Undef
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MKTILE_UNGRIDDED

  subroutine FILLIN_TILE1D(STATE, NAME, TILE, XFORM, RC)
    type(ESMF_STATE),   intent(INOUT) :: STATE
    character(len=*)               :: NAME
    real                           :: TILE(:)
    type (MAPL_LocStreamXFORM)     :: XFORM
    integer, optional, intent(OUT) :: RC

!  Locals

    character(len=ESMF_MAXSTR)   :: IAm='FILLIN_TILE1D'
    integer                      :: STATUS
    real, pointer                :: PTR(:)
    type (ESMF_Field)            :: field
    type (ESMF_StateItem_Flag)   :: itemType

    call ESMF_StateGet(state, name, itemType=itemType, rc=status)
    VERIFY_(STATUS)

    if (itemType == ESMF_STATEITEM_NOTFOUND) then

! If the field is not in the state being filled, we do nothing.
!--------------------------------------------------------------
      RETURN_(ESMF_SUCCESS)

    else

! Get the pointer to the variable to be filled.
!----------------------------------------------

      call MAPL_GetPointer(STATE, PTR, NAME, RC=STATUS)
      VERIFY_(STATUS)

! Fill the variable from the provided stream variable.
!-----------------------------------------------------

      call MAPL_LocStreamTransform( PTR, XFORM, TILE, RC=STATUS )
      VERIFY_(STATUS)
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine FILLIN_TILE1D

  subroutine FILLIN_TILE2D(STATE, NAME, TILE, XFORM, RC)
    type(ESMF_STATE),   intent(INOUT) :: STATE
    character(len=*)               :: NAME
    real                           :: TILE(:,:)
    type (MAPL_LocStreamXFORM)     :: XFORM
    integer, optional, intent(OUT) :: RC

!  Locals

    character(len=ESMF_MAXSTR)   :: IAm='FILLIN_TILE2D'
    integer                      :: STATUS
    real, pointer                :: PTR(:,:)
    integer                      :: I
    type (ESMF_Field)            :: field
    type (ESMF_StateItem_Flag)   :: itemType

    call ESMF_StateGet(state, name, itemType=itemType, rc=status)
    VERIFY_(STATUS)

    if (itemType == ESMF_STATEITEM_NOTFOUND) then

! If the field is not in the state being filled, we do nothing.
!--------------------------------------------------------------

       RETURN_(ESMF_SUCCESS)

! Get the pointer to the variable to be filled.
!----------------------------------------------
    else

      call MAPL_GetPointer(STATE, PTR, NAME, RC=STATUS)
      VERIFY_(STATUS)

! Fill the variable from the provided stream variable.
!-----------------------------------------------------

      do I = 1, SIZE(PTR,2)
         call MAPL_LocStreamTransform( PTR(:,I), XFORM, TILE(:,I), RC=STATUS )
         VERIFY_(STATUS)
      end do
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine FILLIN_TILE2D

  subroutine FILLOUT_TILE1D(STATE, NAME, TILE, XFORM, RC)
    type(ESMF_STATE),   intent(INOUT) :: STATE
    character(len=*),   intent(IN   ) :: NAME
    real                              :: TILE(:)
    type (MAPL_LocStreamXFORM)        :: XFORM
    integer, optional,  intent(OUT  ) :: RC

!  Locals

    character(len=ESMF_MAXSTR)   :: IAm='FILLOUT_TILE1D'
    integer                      :: STATUS
    real, pointer                :: PTR(:)
    type (ESMF_Field)            :: field
    type (ESMF_StateItem_Flag)   :: itemType

    call ESMF_StateGet(state, name, itemType=itemType, rc=status)
    VERIFY_(STATUS)

    if (itemType == ESMF_STATEITEM_NOTFOUND) then

! If the field is not in the state being filled, we do nothing.
!--------------------------------------------------------------

       RETURN_(ESMF_SUCCESS)

    else

      call MAPL_GetPointer(STATE, PTR, NAME, RC=STATUS)
      VERIFY_(STATUS)

      _ASSERT(associated(PTR),'needs informative message')

      call MAPL_LocStreamTransform( TILE, XFORM, PTR, RC=STATUS )
      VERIFY_(STATUS)

    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine FILLOUT_TILE1D

  subroutine FILLOUT_UNGRIDDED(STATE, NAME, TILE, XFORM, RC)
    type(ESMF_STATE),   intent(INOUT) :: STATE
    character(len=*),   intent(IN   ) :: NAME
    real                              :: TILE(:,:)
    type (MAPL_LocStreamXFORM)        :: XFORM
    integer, optional,  intent(OUT  ) :: RC

!  Locals

    character(len=ESMF_MAXSTR)   :: IAm='FILLOUT_UNGRIDDED'
    integer                      :: STATUS
    real, pointer                :: PTR(:,:)
    integer                      :: I
    type (ESMF_Field)            :: field
    type (ESMF_StateItem_Flag)   :: itemType

    call ESMF_StateGet(state, name, itemType=itemType, rc=status)
    VERIFY_(STATUS)

    if (itemType == ESMF_STATEITEM_NOTFOUND) then

! If the field is not in the state being filled, we do nothing.
!--------------------------------------------------------------

       RETURN_(ESMF_SUCCESS)

    else

      call MAPL_GetPointer(STATE, PTR, NAME, RC=STATUS)
      VERIFY_(STATUS)

      _ASSERT(associated(PTR),'needs informative message')
      if (size(tile,2)==0) then
         RETURN_(ESMF_SUCCESS)
      else
         do I = 1, SIZE(PTR,2)
            call MAPL_LocStreamTransform( TILE(:,I), XFORM, PTR(:,I), RC=STATUS )
            VERIFY_(STATUS)
         enddo
      end if
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine FILLOUT_UNGRIDDED

    subroutine RouteRunoff(RoutingType, Runoff, Discharge, rc)
      type(T_RiverRouting),  intent(IN ) :: RoutingType
      real,             intent(IN ) :: Runoff(:)
      real,             intent(OUT) :: Discharge(:)
      integer, optional,intent(OUT) :: rc

      character(len=ESMF_MAXSTR)   :: IAm="RouteRunoff"
      integer                      :: STATUS

      type(T_Routing), pointer :: Routing(:)
      integer, pointer :: karray(:)
      integer, pointer :: kdx(:)
      integer, pointer :: BlockSizes(:), displ(:)

      type(ESMF_VM) :: VM
      integer       :: myPE, nDEs, comm
      integer       :: i
      real          :: TileDischarge
      integer       :: mpstatus(MP_STATUS_SIZE)
      integer :: n, k
      real, allocatable :: td(:), tarray(:)

      call ESMF_VMGetCurrent(VM,                                RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_VMGet       (VM,       mpiCommunicator =comm,   RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_VMGet       (VM, localpet=MYPE, petcount=nDEs,  RC=STATUS)
      VERIFY_(STATUS)

      Routing => RoutingType%LocalRoutings
      karray => RoutingType%karray
      kdx => RoutingType%kdx
      BlockSizes => RoutingType%BlockSizes
      displ => RoutingType%displ

      Discharge   = 0.0

      n=size(kdx)
      allocate(td(n), _STAT)
      allocate(tarray(displ(nDEs)), _STAT)
      do k=1,n
         i=kdx(k)

         TileDischarge = Runoff(Routing(i)%SrcIndex)*Routing(i)%weight
         TileDischarge = max(TileDischarge, 0.0)
         td(k) = TileDischarge
      end do
      call MPI_AllGatherV(td, n, MP_Real, &
           tarray, blocksizes, displ, MP_Real, comm, status)
      _VERIFY(STATUS)

      do i=1,size(Routing)
         if(Routing(i)%DstPE==myPE) then
            n=Routing(i)%srcPe
            k=karray(Routing(i)%seqIdx)
            TileDischarge=tarray(displ(n)+k)
            Discharge(Routing(i)%DstIndex) = Discharge(Routing(i)%DstIndex) + TileDischarge
         end if
      end do
      deallocate(td, _STAT)
      deallocate(tarray, _STAT)

      RETURN_(ESMF_SUCCESS)
    end subroutine RouteRunoff

    subroutine OBIO_fillExports(type, IMPORT, &
                                LOCSTREAM, GIM, &
                                XFORM, &
                                NT, NB_CHOU, &
                                CO2SC, DRBAND, DFBAND, &
                                CO2SCTILE, DRBANDTILE, DFBANDTILE, RC)

      integer,           intent(   IN )     :: TYPE
      type(ESMF_State),  intent(inout)      :: IMPORT ! Import state
      type(MAPL_LocStream), intent(in)      :: LOCSTREAM
      type(ESMF_State), pointer, intent(in) :: GIM(:)
      type(MAPL_LocStreamXFORM), intent(in) :: XFORM

      integer,           intent(   IN) ::  NT
      integer,           intent(   IN) ::  NB_CHOU
      integer, optional, intent(  OUT) ::  RC

      real, pointer ::  CO2SCTILE(:), CO2SC(:,:),&
                        DRBANDTILE(:,:), DRBAND(:,:,:),&
                        DFBANDTILE(:,:), DFBAND(:,:,:)

      character(len=ESMF_MAXSTR), parameter     :: IAm="OBIO_fillExports"
      integer                                   :: STATUS
      integer                                   :: K

!      type (ESMF_State),      pointer         :: GIM(:)
!      type (T_SURFACE_STATE), pointer         :: SURF_INTERNAL_STATE
!      type (MAPL_LocStream)                   :: LOCSTREAM
!      type (MAPL_LocStreamXFORM)              :: XFORM

      call MAPL_GetPointer(IMPORT  , CO2SC   , 'CO2SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, DRBAND    , 'DROBIO'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, DFBAND    , 'DFOBIO'   , RC=STATUS); VERIFY_(STATUS)

      allocate(CO2SCTILE(NT), STAT=STATUS)
      VERIFY_(STATUS)
      allocate(DRBANDTILE(  NT,NB_OBIO), STAT=STATUS)
      VERIFY_(STATUS)
      allocate(DFBANDTILE(  NT,NB_OBIO), STAT=STATUS)
      VERIFY_(STATUS)

      call MAPL_LocStreamTransform(LOCSTREAM, CO2SCTILE, CO2SC,    RC=STATUS); VERIFY_(STATUS)
      do K = 1, NB_OBIO
         call MAPL_LocStreamTransform(LOCSTREAM, DRBANDTILE(:,K), DRBAND(:,:,K), RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_LocStreamTransform(LOCSTREAM, DFBANDTILE(:,K), DFBAND(:,:,K), RC=STATUS)
         VERIFY_(STATUS)
      end do

!      XFORM = surf_internal_state%xform_in(type)
      call FILLIN_TILE(GIM(type), 'CO2SC',     CO2SCTILE,     XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DRBAND',    DRBANDTILE,    XFORM, RC=STATUS); VERIFY_(STATUS)
      call FILLIN_TILE(GIM(type), 'DFBAND',    DFBANDTILE,    XFORM, RC=STATUS); VERIFY_(STATUS)

      RETURN_(ESMF_SUCCESS)
    end subroutine OBIO_fillExports

end module GEOS_SurfaceGridCompMod


