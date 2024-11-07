!  $Id$

#include "MAPL_Generic.h"

module GEOS_LandiceGridCompMod

!=============================================================================
!BOP

! !MODULE: GEOS_LandiceGridCompMod -- Implements slab landice tiles.

! !==========================================================================
! An improved version over the slab landice 

! TODO :

! - Add multiple elevation classes support to account for ice sheet topo changes   

! - Add more layers for a more realistic treatment of ice energy budget

! - Compute surface mass balance (SMB) as an export for coupling to a dynamic
!   ice-sheet model in future
! !==========================================================================

! !INTERFACE:

! !USES:

  use sfclayer   !  use module that contains surface layer routines
  use StieglitzSnow, only:                       &
       snowrt      => StieglitzSnow_snowrt,      &
       SNOW_ALBEDO => StieglitzSnow_snow_albedo, &
       TRID        => StieglitzSnow_trid,        & 
       MINSWE      => StieglitzSnow_MINSWE,      &
       cpw         => StieglitzSnow_CPW,         &
       N_CONSTIT,                                &
       NUM_DUDP, NUM_DUSV, NUM_DUWT, NUM_DUSD,   &
       NUM_BCDP, NUM_BCSV, NUM_BCWT, NUM_BCSD,   &
       NUM_OCDP, NUM_OCSV, NUM_OCWT, NUM_OCSD,   &
       NUM_SUDP, NUM_SUSV, NUM_SUWT, NUM_SUSD,   &
       NUM_SSDP, NUM_SSSV, NUM_SSWT, NUM_SSSD
  use ESMF
  use MAPL
  use GEOS_UtilsMod
  use DragCoefficientsMod
  
  implicit none
  private

  integer, parameter :: ICE   = 1
  integer, parameter :: SNOW  = 2
  integer, parameter :: NUM_SUBTILES = 2
  integer, parameter :: NUM_SNOW_LAYERS = 15
  integer, parameter :: NUM_ICE_LAYERS  = 15
  integer, parameter :: NUM_SNOICE_LAYERS = NUM_SNOW_LAYERS+NUM_ICE_LAYERS
  real,    parameter :: rad_to_deg      = 180.0 / 3.1415926

 
  ! snowrt related constants
  ! will move these to a global module later 
  real,    parameter :: ALHE     = MAPL_ALHL   ! J/kg  @15C
  real,    parameter :: ALHM     = MAPL_ALHF   ! J/kg 
  real,    parameter :: TF       = MAPL_TICE   ! K
  real,    parameter :: RHOW     = MAPL_RHOWTR ! kg/m^3

  real,    parameter :: RHOFRESH = 300.        ! kg/m^3  density of fresh snow
  real,    parameter :: RHOICE   = 917.        ! kg/m^3  pure ice density
  real,    parameter :: MAXSNDZ  = 15.0        ! m
  real,    parameter :: BIG      = 1.e10
  real,    parameter :: condice  = 2.25        ! @ 0 C [W/m/K] 
  real,    parameter :: MINFRACSNO = 1.e-20    ! mininum sno/ice fraction for
                                               ! heat diffusion of ice layers to take effect
  real,    parameter :: LWCTOP     = 1.        ! top thickness to compute LWC. 1m taken from
                                               ! Fettweis et al 2011  
  real,    parameter :: VISMAX    = 0.96       ! parameter for snow_albedo
  real,    parameter :: NIRMAX    = 0.68       ! parameter for snow_albedo
  real,    parameter :: SLOPE     = 1.0        ! parameter for snow_albedo

  ! taken from CICE
   real,   parameter :: &                       ! currently used only
          AWTVDR = 0.00318, &! visible, direct  ! for history and
          AWTIDR = 0.00182, &! near IR, direct  ! diagnostics
          AWTVDF = 0.63282, &! visible, diffuse
          AWTIDF = 0.36218   ! near IR, diffuse  


  !real,    dimension(NUM_SNOW_LAYERS), parameter   :: DZMAX = (/0.08, 0.12, big/)
  real,    dimension(NUM_SNOW_LAYERS), parameter   :: DZMAX = (/0.08, 0.08, 0.08 &
             , 0.15, 0.25, big, big, big, big, big, big, big, big, big, big/)         
  real,    dimension(NUM_ICE_LAYERS), parameter   :: DZMAXI = (/0.08, 0.08, 0.08 &
             , 0.15, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0/)         
  


  integer,    parameter :: TAR_PE     = 43
  integer,    parameter :: TAR_TILE   = 1
  integer               :: N_CONST_LANDICE4SNWALB, AEROSOL_DEPOSITION, CHOOSEMOSFC

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !DESCRIPTION:
! 
!   {\tt GEOS\_Landice} is a light-weight gridded component that updates
!      the landice tiles
!

!EOP

   contains

!BOP
! !IROUTINE: SetServices -- Sets ESMF services for this component

!INTERFACE:

  subroutine SetServices ( GC, RC )

    !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: 
!                This version uses the MAPL\_GenericSetServices, which sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp. The import
!   and internal variables are allocated and initialized by generic.  Here
!   generic is used for tiles.

!EOP

!=============================================================================

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME
    character(len=ESMF_MAXSTR)              :: SURFRC
    type(ESMF_Config)                       :: SCF 

!=============================================================================

    type(MAPL_MetaComp), pointer            :: MAPL


! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run1, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2, RC=STATUS )
    VERIFY_(STATUS)

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource (MAPL, SURFRC, label = 'SURFRC:', default = 'GEOS_SurfaceGridComp.rc', RC=STATUS) ; VERIFY_(STATUS)
    SCF = ESMF_ConfigCreate(rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigLoadFile(SCF,SURFRC,rc=status) ; VERIFY_(STATUS)
    call MAPL_GetResource (SCF, N_CONST_LANDICE4SNWALB, label='N_CONST_LANDICE4SNWALB:', DEFAULT=0, __RC__ )
    call MAPL_GetResource (SCF, AEROSOL_DEPOSITION,     label='AEROSOL_DEPOSITION:',     DEFAULT=0, __RC__ )
    call MAPL_GetResource (SCF, CHOOSEMOSFC,            label='CHOOSEMOSFC:',            DEFAULT=1, __RC__ )
    call ESMF_ConfigDestroy      (SCF, __RC__)

! Set the state variable specs.
! -----------------------------

!BOS

!  !Export state:

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'EMIS',                              &
        LONG_NAME          = 'surface_emissivity',                &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_reflectivity_for_visible_beam',   &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVR',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_reflectivity_for_visible_diffuse',&
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVF',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_reflectivity_for_near_infrared_beam', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNR',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_reflectivity_for_near_infrared_diffuse', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNF',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TST',                               &
        LONG_NAME          = 'surface_temperature',          &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'LST',                               &
        LONG_NAME          = 'land_surface_skin_temperature',          &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'QST',                               &
        LONG_NAME          = 'surface_specific_humidity',         &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TH',                                &
        LONG_NAME          = 'turbulence_surface_skin_temperature', &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'QH',                                &
        LONG_NAME          = 'turbulence_surface_specific_humidity', &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELTS',                             &
        LONG_NAME          = 'change_of_surface_skin_temperature',&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELQS',                             &
        LONG_NAME          = 'change_of_surface_specific_humidity',&
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CHT',                               &
        LONG_NAME          = 'surface_heat_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CMT',                               &
        LONG_NAME          = 'surface_momentum_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CQT',                               &
        LONG_NAME          = 'surface_moisture_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CNT',                               &
        LONG_NAME          = 'neutral_drag_coefficient',          &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'RIT',                               &
        LONG_NAME          = 'surface_bulk_richardson_number',    &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'ACCUM',                             &
        LONG_NAME          = 'net_ice_accumulation_rate',         &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_ice_evaporation_energy_flux_over_glaciated_surface',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPICE_GL'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'sublimation'               ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'SUBLIM'                    ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_mass_over_glaciated_surface'                 ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'SNOMAS_GL'                 ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_mass_over_glaciated_surface'   ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'SNOWMASS'                  ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_depth_over_glaciated_surface'  ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNOWDP_GL'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_snow_covered_area_of_glaciated_surface',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ASNOW_GL'                  ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_layer_density',        &
    UNITS              = 'kg m-3'                    ,&
    SHORT_NAME         = 'RHOSNOW'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_layer_temperature',    &
    UNITS              = 'deg C'                     ,&
    SHORT_NAME         = 'TSNOW'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'aggregated_ice_layer_temperature',    &
    UNITS              = 'deg C'                     ,&
    SHORT_NAME         = 'TICE0'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    UNGRIDDED_DIMS     = (/NUM_ICE_LAYERS/)          ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_layer_water_content',  &
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WSNOW'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_layer_thickness',      &
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'ZSNOW'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_layer_density_change_due_to_densification',   &
    UNITS              = 'kg m-3'                    ,&
    SHORT_NAME         = 'DRHOS0'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_layer_mass_residual_due_to_densification',   &
    UNITS              = 'kg m-2 s-1'                    ,&
    SHORT_NAME         = 'WESNEX'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'total_snow_mass_residual_due_to_densification',   &
    UNITS              = 'kg m-2 s-1'                    ,&
    SHORT_NAME         = 'WESNEXT'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'top_snow_layer_mass_change_due_to_sublimation_and_condensation',   &
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'WESNSC'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'top_snow_layer_thickness_change_due_to_sub_con',   &
    UNITS              = 'm s-1'                     ,&
    SHORT_NAME         = 'SNDZSC'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'top_snow_layer_mass_change_due_to_precip',   &
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'WESNPREC'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'top_snow_layer_thickness_change_due_to_precip',   &
    UNITS              = 'm s-1'                     ,&
    SHORT_NAME         = 'SNDZPREC'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'top_snow_layer_thickness_change_due_to_percolation',   &
    UNITS              = 'm s-1'                     ,&
    SHORT_NAME         = 'SNDZ1PERC'                 ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_layer_mass_change_due_to_percolation',   &
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'WESNPERC'                  ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_layer_mass_change_due_to_densification',   &
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'WESNDENS'                  ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_layer_mass_change_due_to_repartition',   &
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'WESNREPAR'                 ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'frozen_runoff_due_to_fixed_max_depth',   &
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'WESNBOT'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'contribution_to_surface_mass_balance_from_rain_frozen_onto_bare_ice',   &
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'RAINRFZ'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_melt_flux'             ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'SMELT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'ice_melt_flux'             ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'IMELT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_broadband_albedo',   &
    UNITS              = '1'                ,&
    SHORT_NAME         = 'SNOWALB'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'aggregated_snow_ice_broadband_albedo',   &
    UNITS              = '1'                ,&
    SHORT_NAME         = 'SNICEALB'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'melt_water_production',   &
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'MELTWTR'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snowpack_meltwater_content',   &
    UNITS              = 'kg m-2'                ,&
    SHORT_NAME         = 'MELTWTRCONT'               ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'liquid_water_content_in_top_x_m',   &
    UNITS              = '1'                 ,&
    SHORT_NAME         = 'LWC'               ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'runoff_total_flux'         ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'RUNOFF'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GUST',                      &
    LONG_NAME          = 'gustiness',                 &
    UNITS              = 'm s-1',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'VENT',                      &
    LONG_NAME          = 'surface_ventilation_velocity',&
    UNITS              = 'm s-1',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_roughness'         ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'Z0'                        ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_roughness_for_heat',&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'Z0H'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOT2M',                     &
        LONG_NAME          = 'temperature 2m wind from MO sfc', &
        UNITS              = 'K',                         &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOQ2M',                     &
        LONG_NAME          = 'humidity 2m wind from MO sfc',    &
        UNITS              = 'kg kg-1',                   &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU2M',                    &
        LONG_NAME          = 'zonal 2m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV2M',                    &
        LONG_NAME          = 'meridional 2m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOT10M',                     &
        LONG_NAME          = 'temperature 10m wind from MO sfc', &
        UNITS              = 'K',                         &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOQ10M',                     &
        LONG_NAME          = 'humidity 10m wind from MO sfc',    &
        UNITS              = 'kg kg-1',                   &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU10M',                    &
        LONG_NAME          = 'zonal 10m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV10M',                    &
        LONG_NAME          = 'meridional 10m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU50M',                    &
        LONG_NAME          = 'zonal 50m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV50M',                    &
        LONG_NAME          = 'meridional 50m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'evaporation'               ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'EVAPOUT'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHOUT'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'surface_emitted_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLWUP'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'surface_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWNDSRF'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'surface_net_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWNDSRF'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'total_latent_energy_flux'  ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLATN'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'downward_heat_flux_in_ice' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'DNICFLX'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'Ground_heating_snow'       ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'GHSNOW'                    ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'glacier_ice_heating_flux'  ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'GHTSKIN'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
        LONG_NAME          = 'vegetation_type'           ,&
        UNITS              = '1'                         ,&
        SHORT_NAME         = 'ITY'                       ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_1',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU001'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_2',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU002'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_3',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU003'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_4',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU004'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_5',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU005'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_black_carbon_mass_flux_from_the_bottom_layer_bin_1',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTBC001'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_black_carbon_mass_flux_from_the_bottom_layer_bin_2',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTBC002'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_organic_carbon_mass_flux_from_the_bottom_layer_bin_1',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTOC001'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_organic_carbon_mass_flux_from_the_bottom_layer_bin_2',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTOC002'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)

!  !Internal state:

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'TS',                                &
        LONG_NAME          = 'surface_skin_temperature',          &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_SUBTILES/),                    &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 280.0,                               &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'QS',                                &
        LONG_NAME          = 'surface_specific_humidity',         &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_SUBTILES/),                    &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.01,                                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'FR',                                &
        LONG_NAME          = 'ice_fraction',                      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_SUBTILES/),                    &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 1.0/NUM_SUBTILES,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CH',                                &
        LONG_NAME          = 'surface_heat_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_SUBTILES/),                    &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 1.0e-4,                              &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CM',                                &
        LONG_NAME          = 'surface_momentum_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_SUBTILES/),                    &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 1.0e-4,                              &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CQ',                                &
        LONG_NAME          = 'surface_moisture_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_SUBTILES/),                    &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 1.0e-4,                              &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'WESN',                              &
        LONG_NAME          = 'snow_layer_mass',                   &
        UNITS              = 'kg m-2',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/),                 &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'HTSN',                              &
        LONG_NAME          = 'snow_layer_heat_content',           &
        UNITS              = 'J m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/),                 &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'SNDZ',                              &
        LONG_NAME          = 'snow_layer_depth',                  &
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/),                 &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'TICE',                              &
        LONG_NAME          = 'ice_layer_temperature',             &
        UNITS              = 'k',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_SUBTILES, NUM_ICE_LAYERS/),    &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     if (N_CONST_LANDICE4SNWALB /=0) then

        call MAPL_AddInternalSpec(GC,                                &
          SHORT_NAME         = 'IRDU001',                            &
          LONG_NAME          = 'dust_mass_in_snow_bin_1_in_GL',      &
          UNITS              = 'kg m-2',                             &
          DIMS               = MAPL_DimsTileOnly,                    &
          UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/),                  &
          VLOCATION          = MAPL_VLocationNone,                   &
          RESTART            = MAPL_RestartOptional,                 &
          FRIENDLYTO         = trim(COMP_NAME),                      &
                                                          RC=STATUS  ) 
        VERIFY_(STATUS)

        call MAPL_AddInternalSpec(GC,                                &
          SHORT_NAME         = 'IRDU002',                            &
          LONG_NAME          = 'dust_mass_in_snow_bin_2_in_GL',      &
          UNITS              = 'kg m-2',                             &
          DIMS               = MAPL_DimsTileOnly,                    &
          UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/),                  &
          RESTART            = MAPL_RestartOptional,                 &
          VLOCATION          = MAPL_VLocationNone,                   &
          FRIENDLYTO         = trim(COMP_NAME),                      &
                                                          RC=STATUS  ) 
        VERIFY_(STATUS)

        call MAPL_AddInternalSpec(GC,                                &
          SHORT_NAME         = 'IRDU003',                            &
          LONG_NAME          = 'dust_mass_in_snow_bin_3_in_GL',      &
          UNITS              = 'kg m-2',                             &
          DIMS               = MAPL_DimsTileOnly,                    &
          UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/),                  &
          VLOCATION          = MAPL_VLocationNone,                   &
          RESTART            = MAPL_RestartOptional,                 &
          FRIENDLYTO         = trim(COMP_NAME),                      &
                                                          RC=STATUS  ) 
        VERIFY_(STATUS)

        call MAPL_AddInternalSpec(GC,                                &
          SHORT_NAME         = 'IRDU004',                            &
          LONG_NAME          = 'dust_mass_in_snow_bin_4_in_GL',      &
          UNITS              = 'kg m-2',                             &
          DIMS               = MAPL_DimsTileOnly,                    &
          UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/),                  &
          VLOCATION          = MAPL_VLocationNone,                   &
          RESTART            = MAPL_RestartOptional,                 &
          FRIENDLYTO         = trim(COMP_NAME),                      &
                                                          RC=STATUS  ) 
        VERIFY_(STATUS)

        call MAPL_AddInternalSpec(GC,                                &
          SHORT_NAME         = 'IRDU005',                            &
          LONG_NAME          = 'dust_mass_in_snow_bin_5_in_GL',      &
          UNITS              = 'kg m-2',                             &
          DIMS               = MAPL_DimsTileOnly,                    &
          UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/),                  &
          VLOCATION          = MAPL_VLocationNone,                   &
          RESTART            = MAPL_RestartOptional,                 &
          FRIENDLYTO         = trim(COMP_NAME),                      &
                                                          RC=STATUS  ) 
        VERIFY_(STATUS)

        call MAPL_AddInternalSpec(GC,                                 &
          SHORT_NAME         = 'IRBC001',                             &
          LONG_NAME          = 'hydrophobic_black_carbon_mass_in_snow_bin_1_in_GL', &
          UNITS              = 'kg m-2',                              &
          DIMS               = MAPL_DimsTileOnly,                     &
          UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/),                   &
          VLOCATION          = MAPL_VLocationNone,                    &
          RESTART            = MAPL_RestartOptional,                  &
          FRIENDLYTO         = trim(COMP_NAME),                       &
                                                          RC=STATUS  ) 
        VERIFY_(STATUS)

        call MAPL_AddInternalSpec(GC,                                 &
          SHORT_NAME         = 'IRBC002',                             &
          LONG_NAME          = 'hydrophilic_black_carbon_mass_in_snow_bin_2_in_GL', &
          UNITS              = 'kg m-2',                              &
          DIMS               = MAPL_DimsTileOnly,                     &
          UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/),                   &
          VLOCATION          = MAPL_VLocationNone,                    &
          RESTART            = MAPL_RestartOptional,                  &
          FRIENDLYTO         = trim(COMP_NAME),                       &
                                                          RC=STATUS  ) 
        VERIFY_(STATUS)

        call MAPL_AddInternalSpec(GC,                                  &
          SHORT_NAME         = 'IROC001',                              &
          LONG_NAME          = 'hydrophobic_organic_carbon_mass_in_snow_in_GL', &
          UNITS              = 'kg m-2',                               &
          DIMS               = MAPL_DimsTileOnly,                      &
          UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/),                    &
          VLOCATION          = MAPL_VLocationNone,                     &
          RESTART            = MAPL_RestartOptional,                   &
          FRIENDLYTO         = trim(COMP_NAME),                        &
                                                          RC=STATUS  ) 
        VERIFY_(STATUS)

        call MAPL_AddInternalSpec(GC,                                  &
          SHORT_NAME         = 'IROC002',                              &
          LONG_NAME          = 'hydrophilic_organic_carbon_mass_in_snow_in_GL', &
          UNITS              = 'kg m-2',                               &
          DIMS               = MAPL_DimsTileOnly,                      &
          UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/),                    &
          VLOCATION          = MAPL_VLocationNone,                     &
          RESTART            = MAPL_RestartOptional,                   &
          FRIENDLYTO         = trim(COMP_NAME),                        &
                                                          RC=STATUS  ) 
        VERIFY_(STATUS)

     end if

!  !Import state:

! Flux forcings

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'ALW',                               &
        LONG_NAME          = 'linearization_of_surface_emitted_longwave_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'BLW',                               &
        LONG_NAME          = 'linearization_of_surface_emitted_longwave_flux', &
        UNITS              = 'W m-2 K-1',                         &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                             ,&
        LONG_NAME          = 'surface_downwelling_par_beam_flux' ,&
        UNITS              = 'W m-2'                             ,&
        SHORT_NAME         = 'DRPAR'                             ,&
        DIMS               = MAPL_DimsTileOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_par_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DFPAR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_nir_beam_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DRNIR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_nir_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DFNIR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_uvr_beam_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DRUVR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_uvr_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DFUVR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'LWDNSRF',                           &
        LONG_NAME          = 'surface_absorbed_longwave_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'evaporation',                       &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'EVAP ',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'upward_sensible_heat_flux',         &
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'SH',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_evaporation',         &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'DEVAP',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_upward_sensible_heat_flux', &
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'DSH',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'snowfall',                          &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'SNO',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

! Surface air quantities

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_air_temperature',           &
        UNITS              = 'K',                                 &
        SHORT_NAME         = 'TA',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_air_specific_humidity',     &
        UNITS              = 'kg kg-1',                           &
        SHORT_NAME         = 'QA',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_wind_speed',                &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'UU',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'levellm_uwind',                     &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'UWINDLMTILE',                       &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'levellm_vwind',                     &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'VWINDLMTILE',                       &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_layer_height',              &
        UNITS              = 'm',                                 &
        SHORT_NAME         = 'DZ',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_pressure',                  &
        UNITS              = 'Pa',                                &
        SHORT_NAME         = 'PS',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'liquid_water_convective_precipitation',&
        UNITS              = 'kg m-2 s-1'                        ,&
        SHORT_NAME         = 'PCU'                               ,&
        DIMS               = MAPL_DimsTileOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                            ,&
        LONG_NAME          = 'liquid_water_large_scale_precipitation',&
        UNITS              = 'kg m-2 s-1'                       ,&
        SHORT_NAME         = 'PLS'                              ,&
        DIMS               = MAPL_DimsTileOnly                  ,&
        VLOCATION          = MAPL_VLocationNone                 ,&
                                                      RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'THATM',                             &
        LONG_NAME          = 'effective_surface_skin_temperature',&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'QHATM',                             &
        LONG_NAME          = 'effective_surface_specific_humidity',&
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CTATM',                             &
        LONG_NAME          = 'surface_exchange_coefficient_for_heat', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CQATM',                             &
        LONG_NAME          = 'surface_exchange_coefficient_for_moisture', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

! GOSWIM IMPORTS FROM GOCART

    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'dust_dry_depos_all_bins',     &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'DUDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_DUDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'dust_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'DUSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_DUSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'dust_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'DUWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_DUWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'dust_gravity_sett_all_bins',  &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'DUSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_DUSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'black_carbon_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'BCDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_BCDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'black_carbon_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'BCSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_BCSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'black_carbon_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'BCWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_BCWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
 
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'black_carbon_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'BCSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_BCSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'organic_carbon_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'OCDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_OCDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'organic_carbon_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'OCSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_OCSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'organic_carbon_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'OCWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_OCWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'organic_carbon_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'OCSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_OCSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'sulfate_dry_depos_all_bins',  &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SUDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SUDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'sulfate_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SUSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SUSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'sulfate_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SUWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SUWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'sulfate_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SUSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SUSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'sea_salt_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SSDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SSDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'sea_salt_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SSSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SSSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'sea_salt_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SSWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SSWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                              &
         LONG_NAME          = 'sea_salt_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SSSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SSSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)

!EOS

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="RUN1"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN2"  ,RC=STATUS)
    VERIFY_(STATUS)
  
! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)


    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!BOP

! !IROUTINE: RUN1 -- First Run stage for the LandIce component

!INTERFACE:

subroutine RUN1 ( GC, IMPORT, EXPORT, CLOCK, RC )

  !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Periodically refreshes the ozone mixing ratios.

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp),   pointer   :: MAPL
  type (ESMF_State       )            :: INTERNAL
  type (ESMF_Alarm       )            :: ALARM
  type (ESMF_Config      )            :: CF

! pointers to export

   real, pointer, dimension(:  )  :: TH
   real, pointer, dimension(:  )  :: QH
   real, pointer, dimension(:  )  :: LST
   real, pointer, dimension(:  )  :: TST
   real, pointer, dimension(:  )  :: QST
   real, pointer, dimension(:  )  :: CHT
   real, pointer, dimension(:  )  :: CMT
   real, pointer, dimension(:  )  :: CQT
   real, pointer, dimension(:  )  :: CNT
   real, pointer, dimension(:  )  :: RIT
   real, pointer, dimension(:  )  :: VNT
   real, pointer, dimension(:  )  :: GST
   real, pointer, dimension(:  )  :: Z0EXP
   real, pointer, dimension(:  )  :: Z0HEXP
   real, pointer, dimension(:  )  :: MOT2M
   real, pointer, dimension(:  )  :: MOQ2M
   real, pointer, dimension(:  )  :: MOU2M
   real, pointer, dimension(:  )  :: MOV2M
   real, pointer, dimension(:  )  :: MOT10M
   real, pointer, dimension(:  )  :: MOQ10M
   real, pointer, dimension(:  )  :: MOU10M
   real, pointer, dimension(:  )  :: MOV10M
   real, pointer, dimension(:  )  :: MOU50M
   real, pointer, dimension(:  )  :: MOV50M
   real, pointer, dimension(:  )  :: ITY

! pointers to internal

   real, pointer, dimension(:,:)  :: TS
   real, pointer, dimension(:,:)  :: QS
   real, pointer, dimension(:,:)  :: FR
   real, pointer, dimension(:,:)  :: CH
   real, pointer, dimension(:,:)  :: CM
   real, pointer, dimension(:,:)  :: CQ

! pointers to import

   real, pointer, dimension(:)    :: UU
   real, pointer, dimension(:)    :: UWINDLMTILE
   real, pointer, dimension(:)    :: VWINDLMTILE
   real, pointer, dimension(:)    :: DZ     
   real, pointer, dimension(:)    :: TA
   real, pointer, dimension(:)    :: QA     
   real, pointer, dimension(:)    :: PS

   integer                        :: N
   integer                        :: NT
   integer                        :: niter


   real, allocatable              :: URA(:)
   real, allocatable              :: UUU(:)
   real, allocatable              :: LAI(:)
   real, allocatable              :: CN (:)
   real, allocatable              :: RE (:)
   real, allocatable              :: ZT (:)
   real, allocatable              :: ZQ (:)
   real, allocatable              :: PCU(:)
   real, allocatable              :: UCN(:)
   real, allocatable              :: Z0 (:,:)
   real, allocatable              :: WW (:,:)
   real, allocatable              :: U50M (:)
   real, allocatable              :: V50M (:)
   real, allocatable              :: T10M (:)
   real, allocatable              :: Q10M (:)
   real, allocatable              :: U10M (:)
   real, allocatable              :: V10M (:)
   real, allocatable              :: T2M (:)
   real, allocatable              :: Q2M (:)
   real, allocatable              :: U2M (:)
   real, allocatable              :: V2M (:)
   real, allocatable              :: RHO(:)
   real, allocatable              :: VKH(:)
   real, allocatable              :: fakelai(:)
   real, allocatable              :: VKM(:)
   real, allocatable              :: USTAR(:)
   real, allocatable              :: XX(:)
   real, allocatable              :: YY(:)
   real, allocatable              :: CU(:)
   real, allocatable              :: CT(:)
   real, allocatable              :: RIB(:)
   real, allocatable              :: ZETA(:)
   real, allocatable              :: WS(:)
   integer, allocatable           :: IWATER(:)
   real, allocatable              :: PSMB(:)
   real, allocatable              :: PSL(:)

   real, parameter :: LANDICEBAREZ0  = 0.005
   real, parameter :: LANDICESNOWZ0  = 0.001

   integer                        :: CHOOSEZ0
!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Run1"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN1" )

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,             &
         INTERNAL_ESMF_STATE = INTERNAL,         &
                                       RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, CHOOSEZ0, Label="CHOOSEZ0:", DEFAULT=3, RC=STATUS)
    VERIFY_(STATUS)

! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,UU     , 'UU'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,UWINDLMTILE     , 'UWINDLMTILE'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VWINDLMTILE     , 'VWINDLMTILE'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DZ     , 'DZ'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TA     , 'TA'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,QA     , 'QA'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PS     , 'PS'     ,    RC=STATUS)
   VERIFY_(STATUS)

! Pointers to internals
!----------------------

   call MAPL_GetPointer(INTERNAL,TS   , 'TS'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,QS   , 'QS'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,FR   , 'FR'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CH   , 'CH'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CM   , 'CM'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CQ   , 'CQ'     ,    RC=STATUS)
   VERIFY_(STATUS)

! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,QH    , 'QH'      , alloc=.true.,   RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TH    , 'TH'      , alloc=.true.,   RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,QST   , 'QST'     , alloc=.true.,   RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TST   , 'TST'     , alloc=.true.,   RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CHT   , 'CHT'     , alloc=.true.,   RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CMT   , 'CMT'     , alloc=.true.,   RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CQT   , 'CQT'     , alloc=.true.,   RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CNT   , 'CNT'     , alloc=.true.,   RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RIT   , 'RIT'     , alloc=.true.,   RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,LST   , 'LST'     ,                 RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,VNT   , 'VENT'    , RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,GST   , 'GUST'    , RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Z0EXP , 'Z0'      , RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Z0HEXP, 'Z0H'     , RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOT2M, 'MOT2M'   ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOQ2M, 'MOQ2M'   ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOU2M, 'MOU2M'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOV2M, 'MOV2M'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOT10M, 'MOT10M'   ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOQ10M, 'MOQ10M'   ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOU10M, 'MOU10M'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOV10M, 'MOV10M'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOU50M, 'MOU50M'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOV50M, 'MOV50M'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ITY   , 'ITY'     ,    RC=STATUS)
   VERIFY_(STATUS)

   NT = size(TA)

   allocate(URA(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(UUU(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(LAI(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(RE (NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CN (NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(ZT (NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(ZQ (NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(PCU(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(UCN(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(Z0 (size(CM,1),size(CM,2)),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(WW (size(CM,1),size(CM,2)),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(T2M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(Q2M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(U2M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(v2M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(T10M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(Q10M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(U10M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(v10M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(U50M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(v50M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(RHO(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(PSMB(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(PSL(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(VKH(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(fakelai(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(VKM(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(USTAR(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(XX(NT)   ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(YY(NT)   ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CU(NT)   ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CT(NT)   ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(RIB(NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(ZETA(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(WS(NT)   ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(IWATER(NT),STAT=STATUS)
   VERIFY_(STATUS)

!  Compute drag coefficient at tiles
!-----------------------------------
   CHT = 0.0
   CMT = 0.0
   CQT = 0.0
   CNT = 0.0
   RIT = 0.0
   if(associated(GST)) GST = 0.0
   if(associated(MOU50M)) MOU50M = 0.0
   if(associated(MOV50M)) MOV50M = 0.0
   if(associated(MOT10M)) MOT10M = 0.0
   if(associated(MOQ10M)) MOQ10M = 0.0
   if(associated(MOU10M)) MOU10M = 0.0
   if(associated(MOV10M)) MOV10M = 0.0
   if(associated( MOT2M))  MOT2M = 0.0
   if(associated( MOQ2M))  MOQ2M = 0.0
   if(associated( MOU2M))  MOU2M = 0.0
   if(associated( MOV2M))  MOV2M = 0.0
   if(associated(   ITY))    ITY = 13.

   TH  = 0.0
   QH  = 0.0
   TST = 0.0
   QST = 0.0

   do N=1,NUM_SUBTILES

   if(CHOOSEMOSFC.eq.0) then

    call louissurface(4,N,UU,WW,PS,TA,TS,QA,QS,PCU,LAI,Z0,DZ,CM,CN,RIB,ZT,ZQ,CH,CQ,UUU,UCN,RE)

   elseif (CHOOSEMOSFC.eq.1)then

      niter = 6   ! number of internal iterations in the helfand MO surface layer routine
      IWATER = 4
      ! roughness length scale set accroding to Ettema et al. (2010)
      if(N==ICE) then
         Z0(:,N)=LANDICEBAREZ0
      else
         Z0(:,N)=LANDICESNOWZ0
      endif

    PSMB = PS * 0.01            ! convert to MB
    fakelai = 1.e-4
! Approximate pressure at top of surface layer: hydrostatic, eqn of state using avg temp and press
    PSL = PSMB * (1. - (DZ*MAPL_GRAV)/(MAPL_RGAS*(TA+TS(:,N)) ) ) /   &
               (1. + (DZ*MAPL_GRAV)/(MAPL_RGAS*(TA+TS(:,N)) ) )

    call helfsurface( UWINDLMTILE,VWINDLMTILE,TA,TS(:,N),QA,QS(:,N),PSL,PSMB,Z0(:,N),fakelai,  &
                      IWATER,DZ,niter,nt,RHO,VKH,VKM,USTAR,XX,YY,CU,CT,RIB,ZETA,WS, &
                      t2m,q2m,u2m,v2m,t10m,q10m,u10m,v10m,u50m,v50m,CHOOSEZ0)

    CM(:,N)  = VKM
    CH(:,N)  = VKH
    CQ(:,N)  = VKH
      
    CN = (MAPL_KARMAN/ALOG(DZ/Z0(:,N) + 1.0)) * (MAPL_KARMAN/ALOG(DZ/Z0(:,N) + 1.0))
    ZT = Z0(:,N)
    ZQ = Z0(:,N)  
    RE = 0.
    UUU = UU
    UCN = 0.

!  Aggregate to tiles for MO only diagnostics
!--------------------------------------------
      if(associated(MOU50M))MOU50M = MOU50M + U50M(:)*FR(:,N)
      if(associated(MOV50M))MOV50M = MOV50M + V50M(:)*FR(:,N)
      if(associated(MOT10M))MOT10M = MOT10M + T10M(:)*FR(:,N)
      if(associated(MOQ10M))MOQ10M = MOQ10M + Q10M(:)*FR(:,N)
      if(associated(MOU10M))MOU10M = MOU10M + U10M(:)*FR(:,N)
      if(associated(MOV10M))MOV10M = MOV10M + V10M(:)*FR(:,N)
      if(associated(MOT2M))MOT2M = MOT2M + T2M(:)*FR(:,N)
      if(associated(MOQ2M))MOQ2M = MOQ2M + Q2M(:)*FR(:,N)
      if(associated(MOU2M))MOU2M = MOU2M + U2M(:)*FR(:,N)
      if(associated(MOV2M))MOV2M = MOV2M + V2M(:)*FR(:,N)

    endif

!  Aggregate to tiles
!--------------------

      CHT     = CHT + CH(:,N)*FR(:,N)
      CMT     = CMT + CM(:,N)*FR(:,N)
      CQT     = CQT + CQ(:,N)*FR(:,N)
      CNT     = CNT + CN(:  )*FR(:,N)
      RIT     = RIT + RIB(:  )*FR(:,N)

      TH      = TH  + CH(:,N)*TS(:,N)*FR(:,N)
      QH      = QH  + CQ(:,N)*QS(:,N)*FR(:,N)

      TST     = TST + TS(:,N)*FR(:,N)
      QST     = QST + QS(:,N)*FR(:,N)

   end do

   TH = TH /CHT
   QH = QH /CQT
   if(associated(Z0EXP)) Z0EXP = Z0(:,1)
   if(associated(Z0HEXP)) Z0HEXP = ZT
   if(associated(VNT)) VNT = UUU
   if(associated(LST)) LST = TST

   deallocate(URA)
   deallocate(UUU)
   deallocate(LAI)
   deallocate(RE )
   deallocate(CN )
   deallocate(Z0 )
   deallocate(ZT )
   deallocate(ZQ )
   deallocate(WW )
   deallocate(PCU)
   deallocate(UCN)
   deallocate(U50M )
   deallocate(V50M )
   deallocate(T10M )
   deallocate(Q10M )
   deallocate(U10M )
   deallocate(V10M )
   deallocate(T2M )
   deallocate(Q2M )
   deallocate(U2M )
   deallocate(V2M )
   deallocate(RHO)
   deallocate(VKH)
   deallocate(fakelai)
   deallocate(VKM)
   deallocate(USTAR)
   deallocate(XX)
   deallocate(YY)
   deallocate(CU)
   deallocate(CT)
   deallocate(RIB)
   deallocate(ZETA)
   deallocate(WS)
   deallocate(IWATER)
   deallocate(PSMB)
   deallocate(PSL)

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN1" )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)

 end subroutine RUN1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: RUN2 -- Second Run method for the LandIce component

!INTERFACE:

subroutine RUN2 ( GC, IMPORT, EXPORT, CLOCK, RC )

  !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

  !DESCRIPTION: 
!  Periodically refreshes the ozone mixing ratios.

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp), pointer   :: MAPL
  type (ESMF_State       )            :: INTERNAL
  type (ESMF_Alarm       )            :: ALARM
  type (ESMF_Config      )            :: CF

  real, pointer                       :: LATS(:)
  real, pointer                       :: LONS(:)
  !integer, pointer                    :: TILETYPES(:)
  type(MAPL_SunOrbit)                 :: ORBIT

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Run2"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN2")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,             &
         CF        = CF,                         &
         INTERNAL_ESMF_STATE = INTERNAL,         &
         ORBIT     = ORBIT,                      &
         TILELATS  = LATS,                       &
         TILELONS  = LONS,                       &
         !TILETYPES = TILETYPES,                  &     
         RUNALARM  = ALARM,                      &
                                       RC=STATUS )
    VERIFY_(STATUS)

! Do the calculations
!--------------------

    if ( ESMF_AlarmIsRinging(ALARM, RC=STATUS) ) then
       VERIFY_(STATUS)
       call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
       VERIFY_(STATUS)
       call LANDICECORE(RC=STATUS )
       VERIFY_(STATUS)
    end if
    VERIFY_(STATUS)

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN2")
   call MAPL_TimerOff(MAPL,"TOTAL")
   
   RETURN_(ESMF_SUCCESS)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine LANDICECORE(RC)
   integer, optional, intent(OUT) :: RC
     
!  Locals

   character(len=ESMF_MAXSTR)     :: IAm
   integer                        :: STATUS

! pointers to export

   real, pointer, dimension(:  )  :: EMISS
   real, pointer, dimension(:  )  :: ALBVF 
   real, pointer, dimension(:  )  :: ALBVR 
   real, pointer, dimension(:  )  :: ALBNF 
   real, pointer, dimension(:  )  :: ALBNR 
   real, pointer, dimension(:  )  :: DELTS
   real, pointer, dimension(:  )  :: DELQS
   real, pointer, dimension(:  )  :: TST
   real, pointer, dimension(:  )  :: QST
   real, pointer, dimension(:  )  :: EVPICE
   real, pointer, dimension(:  )  :: SUBLIM
   real, pointer, dimension(:  )  :: ACCUM
   real, pointer, dimension(:  )  :: SMELT
   real, pointer, dimension(:  )  :: IMELT
   real, pointer, dimension(:  )  :: RAINRFZ
   real, pointer, dimension(:  )  :: SNOWALB
   real, pointer, dimension(:  )  :: SNICEALB
   real, pointer, dimension(:  )  :: MELTWTR
   real, pointer, dimension(:  )  :: MELTWTRCONT
   real, pointer, dimension(:  )  :: LWC
   real, pointer, dimension(:  )  :: RUNOFF
   real, pointer, dimension(:  )  :: SNOMAS
   real, pointer, dimension(:  )  :: SNOWMASS
   real, pointer, dimension(:  )  :: SNOWDP
   real, pointer, dimension(:  )  :: ASNOW
   real, pointer, dimension(:,:)  :: RHOSNOW
   real, pointer, dimension(:,:)  :: TSNOW
   real, pointer, dimension(:,:)  :: TICE0
   real, pointer, dimension(:,:)  :: WSNOW
   real, pointer, dimension(:,:)  :: ZSNOW
   real, pointer, dimension(:,:)  :: DRHOS0
   real, pointer, dimension(:,:)  :: WESNEX
   real, pointer, dimension(:  )  :: WESNEXT
   real, pointer, dimension(:  )  :: WESC 
   real, pointer, dimension(:  )  :: SDSC 
   real, pointer, dimension(:  )  :: WEPRE 
   real, pointer, dimension(:  )  :: SDPRE
   real, pointer, dimension(:  )  :: SD1PC 
   real, pointer, dimension(:,:)  :: WEPERC
   real, pointer, dimension(:,:)  :: WEREP
   real, pointer, dimension(:  )  :: WEBOT
   real, pointer, dimension(:,:)  :: WEDENS
   real, pointer, dimension(:  )  :: EVAPOUT
   real, pointer, dimension(:  )  :: SHOUT
   real, pointer, dimension(:  )  :: HLATN
   real, pointer, dimension(:  )  :: HLWUP
   real, pointer, dimension(:  )  :: SWNDSRF
   real, pointer, dimension(:  )  :: LWNDSRF
   real, pointer, dimension(:  )  :: DNICFLX
   real, pointer, dimension(:  )  :: GHSNOW
   real, pointer, dimension(:  )  :: GHTSKIN
   real, pointer, dimension(:)   :: RMELTDU001
   real, pointer, dimension(:)   :: RMELTDU002
   real, pointer, dimension(:)   :: RMELTDU003
   real, pointer, dimension(:)   :: RMELTDU004
   real, pointer, dimension(:)   :: RMELTDU005
   real, pointer, dimension(:)   :: RMELTBC001
   real, pointer, dimension(:)   :: RMELTBC002
   real, pointer, dimension(:)   :: RMELTOC001
   real, pointer, dimension(:)   :: RMELTOC002

! pointers to internal

   real, pointer, dimension(:,:)  :: TS
   real, pointer, dimension(:,:)  :: QS
   real, pointer, dimension(:,:)  :: FR
   real, pointer, dimension(:,:)  :: CH
   real, pointer, dimension(:,:)  :: CM
   real, pointer, dimension(:,:)  :: CQ

   real, pointer, dimension(:,:)  :: WESN
   real, pointer, dimension(:,:)  :: HTSN
   real, pointer, dimension(:,:)  :: SNDZ
   real, pointer, dimension(:,:,:):: TICE
   real, pointer, dimension(:,:)  :: IRDU001
   real, pointer, dimension(:,:)  :: IRDU002
   real, pointer, dimension(:,:)  :: IRDU003
   real, pointer, dimension(:,:)  :: IRDU004
   real, pointer, dimension(:,:)  :: IRDU005
   real, pointer, dimension(:,:)  :: IRBC001
   real, pointer, dimension(:,:)  :: IRBC002
   real, pointer, dimension(:,:)  :: IROC001
   real, pointer, dimension(:,:)  :: IROC002

! pointers to import

   real, pointer, dimension(:)    :: ALW
   real, pointer, dimension(:)    :: BLW
   real, pointer, dimension(:)    :: LWDNSRF
   real, pointer, dimension(:)    :: EVAP
   real, pointer, dimension(:)    :: SH
   real, pointer, dimension(:)    :: DEV
   real, pointer, dimension(:)    :: DSH
   real, pointer, dimension(:)    :: SNO
   real, pointer, dimension(:)    :: PS
   real, pointer, dimension(:)    :: PCU
   real, pointer, dimension(:)    :: PLS
   real, pointer, dimension(:)    :: THATM
   real, pointer, dimension(:)    :: QHATM
   real, pointer, dimension(:)    :: CTATM
   real, pointer, dimension(:)    :: CQATM
   real, pointer, dimension(:)    :: DRPAR
   real, pointer, dimension(:)    :: DFPAR
   real, pointer, dimension(:)    :: DRNIR
   real, pointer, dimension(:)    :: DFNIR
   real, pointer, dimension(:)    :: DRUVR
   real, pointer, dimension(:)    :: DFUVR
   real, pointer, dimension(:)    :: TA
   real, pointer, dimension(:)    :: UU
   real, pointer, dimension(:,:)  :: DUDP
   real, pointer, dimension(:,:)  :: DUSV
   real, pointer, dimension(:,:)  :: DUWT
   real, pointer, dimension(:,:)  :: DUSD
   real, pointer, dimension(:,:)  :: BCDP
   real, pointer, dimension(:,:)  :: BCSV
   real, pointer, dimension(:,:)  :: BCWT
   real, pointer, dimension(:,:)  :: BCSD
   real, pointer, dimension(:,:)  :: OCDP
   real, pointer, dimension(:,:)  :: OCSV
   real, pointer, dimension(:,:)  :: OCWT
   real, pointer, dimension(:,:)  :: OCSD
   real, pointer, dimension(:,:)  :: SUDP
   real, pointer, dimension(:,:)  :: SUSV
   real, pointer, dimension(:,:)  :: SUWT
   real, pointer, dimension(:,:)  :: SUSD
   real, pointer, dimension(:,:)  :: SSDP
   real, pointer, dimension(:,:)  :: SSSV
   real, pointer, dimension(:,:)  :: SSWT
   real, pointer, dimension(:,:)  :: SSSD

   real,    allocatable           :: SHF(:)
   real,    allocatable           :: LHF(:)
   real,    allocatable           :: SHD(:)
   real,    allocatable           :: LHD(:)
   real,    allocatable           :: CFQ(:)
   real,    allocatable           :: CFT(:)
   real,    allocatable           :: MLT(:)
   real,    allocatable           :: DTS(:)
   real,    allocatable           :: DQS(:)
   real,    allocatable           :: SWN(:)
   real,    allocatable           :: DIF(:)
   real,    allocatable           :: ULW(:)

   real                           :: DT
   real                           :: LANDICECAP
   real                           :: LANDICEALB
   real                           :: LANDICEEMISS
   real                           :: LANDICEDEPTH
   real                           :: LANDICECOND
   real                           :: LANDICETDEEP

   ! snowrt related vars
   real,  dimension(1)            :: ZONEAREA
   real                           :: ZC1, ZDEP, ALPHA, ZKL
   real,  dimension(1)            :: TKGND
   real,  dimension(NUM_SNOW_LAYERS)    :: TKSNO

   real,    allocatable           :: PRECIP  (:)
   real,    allocatable           :: RAIN    (:)
   real,    allocatable           :: RAINRF  (:)
   real,    allocatable           :: PERC    (:)
   real,    allocatable           :: MELTI   (:)
   real,    allocatable           :: FROZFRAC(:,:)
   real,    allocatable           :: TPSN    (:,:)
   real,    allocatable           :: AREASC  (:)
   real,    allocatable           :: HCORR   (:)
   real,    allocatable           :: ghflxsno(:)
   real,    allocatable           :: ghflxice(:)
   real,    allocatable           :: HLWO    (:)
   real,    allocatable           :: EVAPO   (:)
   real,    allocatable           :: SHFO    (:)
   real,    allocatable           :: LHFO    (:)
   real,    allocatable           :: ZTH     (:)
   real,    allocatable           :: SLR     (:)
   real,    allocatable           :: EVAPI   (:)
   real,    allocatable           :: DEVAPDT (:)
   integer, allocatable           :: ITYPE   (:)
   real,    allocatable           :: LAI     (:)
   real,    allocatable           :: GRN     (:)
   real,    allocatable           :: MODISFAC(:)
   real,    allocatable           :: SNOVR(:), SNONR(:), SNOVF(:), SNONF(:) 
   real,    allocatable           :: LNDVR(:), LNDNR(:), LNDVF(:), LNDNF(:) 
   real,    allocatable           :: VSUVR   (:)
   real,    allocatable           :: VSUVF   (:)
   real,    allocatable           :: SWNETSNOW(:)
   real,    allocatable           :: RADDN   (:)
   real,    allocatable           :: FHGND   (:)
   real,    allocatable           :: DRHO0   (:,:)
   real,    allocatable           :: EXCS    (:,:)
   real,    allocatable           :: WESNSC(:), SNDZSC(:), WESNPREC(:),   & 
                                     SNDZPREC(:),  SNDZ1PERC(:) 
   real,    allocatable           :: WESNPERC(:,:), WESNDENS(:,:), WESNREPAR(:,:) 
   real,    allocatable           :: WESNBOT(:)
   real,    allocatable           :: LANDICELT(:)
   real,    allocatable           :: RCONSTIT(:,:,:)
   real,    allocatable           :: TOTDEPOS(:,:)
   real,    allocatable           :: RMELT(:,:)
   real,    allocatable           :: WESNN(:,:)
   real,    allocatable           :: HTSNN(:,:)
   real,    allocatable           :: SNDZN(:,:)


   ! snow albedo calculation stuff

   type(ESMF_Time)                :: CURRENT_TIME
   type(ESMF_TimeInterval)        :: DELT
   real                           :: DT_SOLAR
   type (ESMF_TimeInterval)       :: TINT
   type (ESMF_Time)               :: NOW
   type (ESMF_Time)               :: BEFORE
   type (ESMF_Time)               :: MODELSTART
   type (ESMF_Alarm)              :: SOLALARM
   logical                        :: solalarmison
   logical                        :: debugzth

   integer                        :: N, NT
   integer                        :: K, L, KL

   ! vars for debugging
   type(ESMF_VM)                  ::  VM
   real                           ::  LONSD, LATSD
   integer                        ::  comm, mype

   real, parameter :: LANDICEALB_   = 0.58
   real, parameter :: LANDICEEMISS_ = 0.99999
   real, parameter :: LANDICEDEPTH_ = 0.07 ! water equiv depth of top layer
   real, parameter :: LANDICECOND_  = 1.2  ! ice conductivity divided by depth of bottom layer
   real, parameter :: LANDICETDEEP_ = 230. ! deep ice temperature

!  Begin...
!----------

   IAm =  trim(COMP_NAME) // "LANDICECORE"

! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,ALW    , 'ALW'    , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,BLW    , 'BLW'    , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,LWDNSRF, 'LWDNSRF', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,EVAP   , 'EVAP'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SH     , 'SH'     , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DEV    , 'DEVAP'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DSH    , 'DSH'    , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PS     , 'PS'     , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PCU    , 'PCU'    , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PLS    , 'PLS'    , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SNO    , 'SNO'    , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,THATM  , 'THATM'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,QHATM  , 'QHATM'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,CTATM  , 'CTATM'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,CQATM  , 'CQATM'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DRPAR  , 'DRPAR'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DFPAR  , 'DFPAR'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DRNIR  , 'DRNIR'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DFNIR  , 'DFNIR'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DRUVR  , 'DRUVR'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DFUVR  , 'DFUVR'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TA     , 'TA'     , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,UU     , 'UU'     , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DUDP   , 'DUDP'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DUSV   , 'DUSV'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DUWT   , 'DUWT'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DUSD   , 'DUSD'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,BCDP   , 'BCDP'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,BCSV   , 'BCSV'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,BCWT   , 'BCWT'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,BCSD   , 'BCSD'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,OCDP   , 'OCDP'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,OCSV   , 'OCSV'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,OCWT   , 'OCWT'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,OCSD   , 'OCSD'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SUDP   , 'SUDP'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SUSV   , 'SUSV'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SUWT   , 'SUWT'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SUSD   , 'SUSD'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SSDP   , 'SSDP'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SSSV   , 'SSSV'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SSWT   , 'SSWT'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SSSD   , 'SSSD'   , RC=STATUS); VERIFY_(STATUS)

! Pointers to internals
!----------------------

   call MAPL_GetPointer(INTERNAL,TS   , 'TS'     , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,QS   , 'QS'     , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,FR   , 'FR'     , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CH   , 'CH'     , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CM   , 'CM'     , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CQ   , 'CQ'     , RC=STATUS); VERIFY_(STATUS)

   call MAPL_GetPointer(INTERNAL,WESN , 'WESN'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,HTSN , 'HTSN'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,SNDZ , 'SNDZ'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,TICE , 'TICE'   , RC=STATUS); VERIFY_(STATUS)
   if (N_CONST_LANDICE4SNWALB /= 0) then
      call MAPL_GetPointer(INTERNAL,IRDU001 ,'IRDU001', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL,IRDU002 ,'IRDU002', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL,IRDU003 ,'IRDU003', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL,IRDU004 ,'IRDU004', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL,IRDU005 ,'IRDU005', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL,IRBC001 ,'IRBC001', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL,IRBC002 ,'IRBC002', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL,IROC001 ,'IROC001', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL,IROC002 ,'IROC002', RC=STATUS); VERIFY_(STATUS)
   end if

! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,EMISS  , 'EMIS'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBVF  , 'ALBVF'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBVR  , 'ALBVR'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBNF  , 'ALBNF'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBNR  , 'ALBNR'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,DELTS  , 'DELTS'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,DELQS  , 'DELQS'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,EVPICE , 'EVPICE_GL' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SUBLIM , 'SUBLIM' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ACCUM  , 'ACCUM'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SMELT  , 'SMELT'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,IMELT  , 'IMELT'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RAINRFZ, 'RAINRFZ', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SNOWALB, 'SNOWALB', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SNICEALB,'SNICEALB',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MELTWTR, 'MELTWTR', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MELTWTRCONT, 'MELTWTRCONT', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,LWC    , 'LWC'    , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RUNOFF , 'RUNOFF' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SNOMAS , 'SNOMAS_GL' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SNOWMASS,'SNOWMASS',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SNOWDP , 'SNOWDP_GL' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ASNOW  , 'ASNOW_GL'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RHOSNOW, 'RHOSNOW', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TSNOW  , 'TSNOW'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TICE0  , 'TICE0'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,WSNOW  , 'WSNOW'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ZSNOW  , 'ZSNOW'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,DRHOS0 , 'DRHOS0' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,WESNEX , 'WESNEX' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,WESNEXT, 'WESNEXT', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,WESC   , 'WESNSC' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SDSC   , 'SNDZSC' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,WEPRE  , 'WESNPREC' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SDPRE  , 'SNDZPREC' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SD1PC  , 'SNDZ1PERC', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,WEPERC , 'WESNPERC' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,WEDENS , 'WESNDENS' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,WEREP  , 'WESNREPAR', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,WEBOT  , 'WESNBOT', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,QST    , 'QST'    , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TST    , 'TST'    , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,EVAPOUT, 'EVAPOUT' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SHOUT  , 'SHOUT'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLATN  , 'HLATN'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLWUP  , 'HLWUP'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,LWNDSRF, 'LWNDSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SWNDSRF, 'SWNDSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,DNICFLX, 'DNICFLX' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,GHSNOW,  'GHSNOW'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,GHTSKIN, 'GHTSKIN' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RMELTDU001,'RMELTDU001',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RMELTDU002,'RMELTDU002',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RMELTDU003,'RMELTDU003',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RMELTDU004,'RMELTDU004',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RMELTDU005,'RMELTDU005',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RMELTBC001,'RMELTBC001',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RMELTBC002,'RMELTBC002',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RMELTOC001,'RMELTOC001',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RMELTOC002,'RMELTOC002',  RC=STATUS); VERIFY_(STATUS)

! Get the time step
! -----------------

    call MAPL_Get(MAPL, HEARTBEAT = DT, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, DT, Label="DT:", DEFAULT=DT, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, LANDICEEMISS, Label="LANDICEEMISS:",  &
                            DEFAULT=LANDICEEMISS_, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, LANDICEALB,   Label="LANDICEALBEDO:", &
                            DEFAULT=LANDICEALB_,   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, LANDICEDEPTH, Label="LANDICEDEPTH:",  &
                            DEFAULT=LANDICEDEPTH_, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, LANDICECOND , Label="LANDICECOND:",   &
                            DEFAULT=LANDICECOND_, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, LANDICETDEEP, Label="LANDICETDEEP:",  &
                            DEFAULT=LANDICETDEEP_, RC=STATUS)
    VERIFY_(STATUS)

    NT = size(ALW)

    allocate(MLT (NT), STAT=STATUS)
    VERIFY_(STATUS)                
    allocate(DTS (NT), STAT=STATUS)
    VERIFY_(STATUS)                
    allocate(DQS (NT), STAT=STATUS)
    VERIFY_(STATUS)                
    allocate(SHF (NT), STAT=STATUS)
    VERIFY_(STATUS)                
    allocate(LHF (NT), STAT=STATUS)
    VERIFY_(STATUS)                
    allocate(SHD (NT), STAT=STATUS)
    VERIFY_(STATUS)                
    allocate(LHD (NT), STAT=STATUS)
    VERIFY_(STATUS)                
    allocate(CFT (NT), STAT=STATUS)
    VERIFY_(STATUS)                
    allocate(CFQ (NT), STAT=STATUS)
    VERIFY_(STATUS)                
    allocate(SWN (NT), STAT=STATUS)
    VERIFY_(STATUS)                
    allocate(DIF (NT), STAT=STATUS)
    VERIFY_(STATUS)                
    allocate(ULW (NT), STAT=STATUS)
    VERIFY_(STATUS)

    allocate(PRECIP(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(RAIN(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(RAINRF(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(PERC(NT)  , STAT=STATUS)
    VERIFY_(STATUS)
    allocate(MELTI(NT) , STAT=STATUS)
    VERIFY_(STATUS)
    allocate(FROZFRAC(NT,NUM_SNOW_LAYERS), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(TPSN(NT,NUM_SNOW_LAYERS)    , STAT=STATUS)
    VERIFY_(STATUS)
    allocate(AREASC(NT), STAT=STATUS)
    VERIFY_(STATUS)
    allocate(HCORR(NT) , STAT=STATUS)
    VERIFY_(STATUS)
    allocate(ghflxsno(NT) , STAT=STATUS)
    VERIFY_(STATUS)
    allocate(ghflxice(NT) , STAT=STATUS)
    VERIFY_(STATUS)
    allocate(HLWO(NT)  , STAT=STATUS)
    VERIFY_(STATUS)
    allocate(EVAPO(NT) , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(LHFO(NT)  , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(SHFO(NT)  , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(ZTH(NT)   , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(SLR(NT)   , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(EVAPI(NT) , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(DEVAPDT(NT), STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(ITYPE(NT) , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(LAI(NT)   , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(GRN(NT)   , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(MODISFAC(NT), STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(SNOVR(NT), SNONR(NT), SNOVF(NT), SNONF(NT) , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(LNDVR(NT), LNDNR(NT), LNDVF(NT), LNDNF(NT) , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(VSUVR(NT) , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(VSUVF(NT) , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(SWNETSNOW(NT) , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(RADDN(NT) , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(FHGND(NT) , STAT=STATUS)
    VERIFY_(STATUS) 
    allocate(DRHO0(NT,NUM_SNOW_LAYERS)    , STAT=STATUS)
    VERIFY_(STATUS)
    allocate(EXCS(NT,NUM_SNOW_LAYERS)     , STAT=STATUS)
    VERIFY_(STATUS)
    allocate(WESNSC(NT), SNDZSC(NT), WESNPREC(NT),   & 
             SNDZPREC(NT),  SNDZ1PERC(NT),           & 
             WESNBOT(NT),                            &
             STAT=STATUS) 
    VERIFY_(STATUS)
    allocate(WESNPERC(NT,NUM_SNOW_LAYERS),            &
             WESNDENS(NT,NUM_SNOW_LAYERS),            &
             WESNREPAR(NT,NUM_SNOW_LAYERS),           &
             STAT=STATUS) 
    VERIFY_(STATUS)
    allocate(LANDICELT(NT) , STAT=STATUS)
    VERIFY_(STATUS) 

    allocate(WESNN  (NUM_SNOW_LAYERS,NT))
    VERIFY_(STATUS)  
    allocate(HTSNN  (NUM_SNOW_LAYERS,NT))
    VERIFY_(STATUS)  
    allocate(SNDZN  (NUM_SNOW_LAYERS,NT))
    VERIFY_(STATUS)  
    allocate(RCONSTIT(NT, NUM_SNOW_LAYERS, N_CONSTIT) , STAT=STATUS)
    VERIFY_(STATUS)  
    allocate(TOTDEPOS(NT, N_CONSTIT) , STAT=STATUS)
    VERIFY_(STATUS)  
    allocate(RMELT(NT, N_CONSTIT) , STAT=STATUS)
    VERIFY_(STATUS)      

    call ESMF_VMGetCurrent(VM,                                RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_VMGet       (VM,       mpiCommunicator =comm,   RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_VMGet(VM, localPet=mype, rc=status)
    VERIFY_(STATUS)

    if(associated(EVAPOUT ))  EVAPOUT  = 0.0 
    if(associated(SUBLIM  ))  SUBLIM   = 0.0
    if(associated(SHOUT   ))  SHOUT    = 0.0 
    if(associated(HLATN   ))  HLATN    = 0.0 
    if(associated(DELTS   ))  DELTS    = 0.0 
    if(associated(DELQS   ))  DELQS    = 0.0 
    if(associated(SWNDSRF ))  SWNDSRF  = 0.0
    if(associated(LWNDSRF ))  LWNDSRF  = 0.0
    if(associated(DNICFLX ))  DNICFLX  = 0.0 
    if(associated(GHSNOW  ))  GHSNOW   = 0.0 
    if(associated(GHTSKIN ))  GHTSKIN  = 0.0 
    if(associated(IMELT   ))  IMELT    = 0.0 
    if(associated(RUNOFF  ))  RUNOFF   = 0.0 
    if(associated(EVPICE  ))  EVPICE   = 0.0 
    if(associated(HLWUP   ))  HLWUP    = 0.0
    if(associated(TICE0   ))  TICE0    = 0.0
    if(associated(ACCUM   ))  ACCUM    = 0.0 
    if(associated(MELTWTR ))  MELTWTR  = 0.0
    
    TOTDEPOS  = 0.0
    RCONSTIT  = 0.0
    RMELT     = 0.0

    ! Zero the light-absorbing aerosol (LAA) deposition rates from  GOCART:
    
    select case (AEROSOL_DEPOSITION)
    case (0)
       DUDP(:,:)=0.
       DUSV(:,:)=0.
       DUWT(:,:)=0.
       DUSD(:,:)=0.
       BCDP(:,:)=0.
       BCSV(:,:)=0.
       BCWT(:,:)=0.
       BCSD(:,:)=0.
       OCDP(:,:)=0.
       OCSV(:,:)=0.
       OCWT(:,:)=0.
       OCSD(:,:)=0.
       
    case (2)
       DUDP(:,:)=0.
       DUSV(:,:)=0.
       DUWT(:,:)=0.
       DUSD(:,:)=0.
       
    case (3)
       BCDP(:,:)=0.
       BCSV(:,:)=0.
       BCWT(:,:)=0.
       BCSD(:,:)=0.
       
    case (4)
       OCDP(:,:)=0.
       OCSV(:,:)=0.
       OCWT(:,:)=0.
       OCSD(:,:)=0.
       
    end select

! Convert the dimentions for LAAs from GEOS_SurfGridComp.F90 to GEOS_LandIceGridComp.F90
! Note: Explanations of each variable
! TOTDEPOS(:,1): Combined dust deposition from size bin 1 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,2): Combined dust deposition from size bin 2 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,3): Combined dust deposition from size bin 3 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,4): Combined dust deposition from size bin 4 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,5): Combined dust deposition from size bin 5 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,6): Combined hydrophobic BC deposition from size bin 1 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,7): Combined hydrophilic BC deposition from size bin 2 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,8): Combined hydrophobic OC deposition from size bin 1 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,9): Combined hydrophilic OC deposition from size bin 2 (dry, conv-scav, ls-scav, sed)
!============================= Possible future applications ====================================
! TOTDEPOS(:,10): Combined sulfate deposition from size bin 3 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,11): Combined sea salt deposition from size bin 1 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,12): Combined sea salt deposition from size bin 2 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,13): Combined sea salt deposition from size bin 3 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,14): Combined sea salt deposition from size bin 4 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,15): Combined sea salt deposition from size bin 5 (dry, conv-scav, ls-scav, sed)

        TOTDEPOS(:,1) = DUDP(:,1) + DUSV(:,1) + DUWT(:,1) + DUSD(:,1)
        TOTDEPOS(:,2) = DUDP(:,2) + DUSV(:,2) + DUWT(:,2) + DUSD(:,2)
        TOTDEPOS(:,3) = DUDP(:,3) + DUSV(:,3) + DUWT(:,3) + DUSD(:,3)
        TOTDEPOS(:,4) = DUDP(:,4) + DUSV(:,4) + DUWT(:,4) + DUSD(:,4)
        TOTDEPOS(:,5) = DUDP(:,5) + DUSV(:,5) + DUWT(:,5) + DUSD(:,5)
        TOTDEPOS(:,6) = BCDP(:,1) + BCSV(:,1) + BCWT(:,1) + BCSD(:,1)
        TOTDEPOS(:,7) = BCDP(:,2) + BCSV(:,2) + BCWT(:,2) + BCSD(:,2)
        TOTDEPOS(:,8) = OCDP(:,1) + OCSV(:,1) + OCWT(:,1) + OCSD(:,1)
        TOTDEPOS(:,9) = OCDP(:,2) + OCSV(:,2) + OCWT(:,2) + OCSD(:,2)
!============================= Possible future applications ====================================
!        TOTDEPOS(:,10) = SUDP(:,1) + SUSV(:,1) + SUWT(:,1) + SUSD(:,1)
!        TOTDEPOS(:,11) = SSDP(:,1) + SSSV(:,1) + SSWT(:,1) + SSSD(:,1)
!        TOTDEPOS(:,12) = SSDP(:,2) + SSSV(:,2) + SSWT(:,2) + SSSD(:,2)
!        TOTDEPOS(:,13) = SSDP(:,3) + SSSV(:,3) + SSWT(:,3) + SSSD(:,3)
!        TOTDEPOS(:,14) = SSDP(:,4) + SSSV(:,4) + SSWT(:,4) + SSSD(:,4)
!        TOTDEPOS(:,15) = SSDP(:,5) + SSSV(:,5) + SSWT(:,5) + SSSD(:,5)

! --------------- GOSWIM PROGRNOSTICS ---------------------------

! Conversion of the masses of the snow impurities
! Note: Explanations of each variable
! Number of snow layer is 15: N = 1-15
! RCONSTIT(NT,N,1): Dust mass from bin 1 in layer N
! RCONSTIT(NT,N,2): Dust mass from bin 2 in layer N
! RCONSTIT(NT,N,3): Dust mass from bin 3 in layer N
! RCONSTIT(NT,N,4): Dust mass from bin 4 in layer N
! RCONSTIT(NT,N,5): Dust mass from bin 5 in layer N
! RCONSTIT(NT,N,6): Hydrophobic BC mass from bin 1 in layer N
! RCONSTIT(NT,N,7): Hydrophilic BC mass from bin 2 in layer N
! RCONSTIT(NT,N,8): Hydrophobic OC mass from bin 1 in layer N
! RCONSTIT(NT,N,9): Hydrophilic OC mass from bin 2 in layer N
!============================= Possible future applications ====================================
! RCONSTIT(NT,N,10): Sulfate mass from size bin 3 in layer N
! RCONSTIT(NT,N,11): Sea salt mass from size bin 1 in layer N
! RCONSTIT(NT,N,12): Sea salt mass from size bin 2 in layer N
! RCONSTIT(NT,N,13): Sea salt mass from size bin 3 in layer N
! RCONSTIT(NT,N,14): Sea salt mass from size bin 4 in layer N
! RCONSTIT(NT,N,15): Sea salt mass from size bin 5 in layer N

    if (N_CONST_LANDICE4SNWALB /=0) then
        RCONSTIT(:,:,1) = IRDU001(:,:)
        RCONSTIT(:,:,2) = IRDU002(:,:)
        RCONSTIT(:,:,3) = IRDU003(:,:)
        RCONSTIT(:,:,4) = IRDU004(:,:)
        RCONSTIT(:,:,5) = IRDU005(:,:)
        RCONSTIT(:,:,6) = IRBC001(:,:)
        RCONSTIT(:,:,7) = IRBC002(:,:)
        RCONSTIT(:,:,8) = IROC001(:,:)
        RCONSTIT(:,:,9) = IROC002(:,:)

!============================= Possible future applications ====================================
!        RCONSTIT(:,:,10) = IRSU003(:,:)
!        RCONSTIT(:,:,11) = IRSS001(:,:)
!        RCONSTIT(:,:,12) = IRSS002(:,:)
!        RCONSTIT(:,:,13) = IRSS003(:,:)
!        RCONSTIT(:,:,14) = IRSS004(:,:)
!        RCONSTIT(:,:,15) = IRSS005(:,:)
    end if

    LANDICELT = 0.0 
    ZONEAREA  = 1.0 
    ! zc1 is not the actual thickness, but the vertical coordinate which is +ve upward
    ZC1       = -DZMAXI(1) * 0.5
    TKGND     = condice ! use value for ice at 0 degC
    PRECIP    = PCU + PLS + SNO
    RAIN      = PCU + PLS 
    PERC      = 0.0
    MELTI     = 0.0
    FROZFRAC  = 0.0
    TPSN      = 0.0  
    AREASC    = 0.0
    HCORR     = 0.0
    ghflxsno  = 0.0
    ghflxice  = 0.0

    HLWO      = 0.0
    EVAPO     = 0.0
    SHFO      = 0.0
    LHFO      = 0.0
    DTS       = 0.0

    ! just to be safe
    DRHO0     = 0.0
    EXCS      = 0.0
    WESNSC    = 0.0
    SNDZSC    = 0.0
    WESNPREC  = 0.0
    SNDZPREC  = 0.0  
    SNDZ1PERC = 0.0 
    WESNPERC  = 0.0 
    WESNDENS  = 0.0 
    WESNREPAR = 0.0 
    RAINRF    = 0.0 
    MLT       = 0.0 
    LNDVR     = 0.0
    LNDNR     = 0.0
    LNDVF     = 0.0
    LNDNF     = 0.0

    debugzth = .false.

    ! --------------------------------------------------------------------------
    ! Get the current time. 
    ! --------------------------------------------------------------------------

    call ESMF_ClockGet( CLOCK, currTime=CURRENT_TIME, startTime=MODELSTART, TIMESTEP=DELT,  RC=STATUS )
    VERIFY_(STATUS)
    if (MAPL_AM_I_Root(VM).and.debugzth) then
      print *,' start time of clock '
      CALL ESMF_TimePrint ( MODELSTART, OPTIONS="string", RC=STATUS )
    endif

        ! --------------------------------------------------------------------------
        ! retrieve the zenith angle
        ! --------------------------------------------------------------------------

!! The next sequence is to make sure that the albedo here and in solar are in sync
!!
! Need to know when Solar was called last, so first get the solar alarm
        call ESMF_ClockGetAlarm ( CLOCK, alarmname="SOLAR_Alarm", ALARM=SOLALARM, RC=STATUS )
        VERIFY_(STATUS)
! Get the interval of the solar alarm - first get it in seconds
        call ESMF_ConfigGetAttribute ( CF, DT_SOLAR, Label="SOLAR_DT:", DEFAULT=DT, RC=STATUS )
        VERIFY_(STATUS)
! Now make an ESMF interval from the increment in seconds
        CALL ESMF_TimeIntervalSet ( TINT, S=NINT(DT_SOLAR), RC=STATUS )
        VERIFY_(STATUS)
! Now print out the solar alarm interval
        if (MAPL_AM_I_Root(VM).and.debugzth) CALL ESMF_TimeIntervalPrint ( TINT, OPTIONS="string", RC=STATUS )
! Now find out if it is ringing now: if so, set "BEFORE" to last time it rang before now
         solalarmison = ESMF_AlarmIsRinging(SOLALARM,RC=STATUS)
         VERIFY_(STATUS)
         if (MAPL_AM_I_Root(VM).and.debugzth)print *,' logical for solar alarm ',solalarmison
!     if so, set "BEFORE" to last time it rang before now
        if(solalarmison) then
         if (MAPL_AM_I_Root(VM).and.debugzth)print *,' In catch, solar alarm is ringing '
         NOW = CURRENT_TIME
         BEFORE = NOW - TINT
! Now print out the last time solar alarm rang
         if (MAPL_AM_I_Root(VM).and.debugzth)CALL ESMF_TimePrint ( BEFORE, OPTIONS="string", RC=STATUS )
!     If alarm is not ringing now, find out when it rang last
        else
         if (MAPL_AM_I_Root(VM).and.debugzth)print *,' In catch, solar alarm is not ringing '
         call ESMF_AlarmGet ( SOLALARM, prevRingTime=BEFORE, RC=STATUS )
         VERIFY_(STATUS)
! PrevRingTime can lie: if alarm never went off yet it gives next alarm time, not prev.
         if(BEFORE > CURRENT_TIME) then
          BEFORE = BEFORE-TINT
          if (MAPL_AM_I_Root(VM).and.debugzth)print *,' In catch, solar alarm not ringing, prev time lied '
          if (MAPL_AM_I_Root(VM).and.debugzth)CALL ESMF_TimePrint ( BEFORE, OPTIONS="string", RC=STATUS )
         else
          if (MAPL_AM_I_Root(VM).and.debugzth)print *,' In catch, solar alarm not ringing, prev time okay '
          if (MAPL_AM_I_Root(VM).and.debugzth)CALL ESMF_TimePrint ( BEFORE, OPTIONS="string", RC=STATUS )
         endif
! Now print out the last time solar alarm rang
        endif

! Get the zenith angle at the center of the time between the last solar call and the next one
        call MAPL_SunGetInsolation(LONS, LATS,      &
            ORBIT, ZTH, SLR, &
            INTV = TINT,     &
            currTime=BEFORE+DELT,  &
            RC=STATUS )
        VERIFY_(STATUS)

        ZTH = max(0.0,ZTH)

    do N=1,NUM_SUBTILES  

    CFT = (CH(:,N)/CTATM)
    CFQ = (CQ(:,N)/CQATM)
    SHF = CFT*(SH   + DSH*(TS(:,N)-THATM))
    LHF = CFQ*(EVAP + DEV*(QS(:,N)-QHATM))*MAPL_ALHS
    SHD = CFT*DSH
    LHD = CFQ*DEV*MAPL_ALHS*GEOS_DQSAT(TS(:,N), PS, PASCALS=.TRUE., RAMP=0.0)
    SWN = ((DRUVR+DRPAR+DRNIR) + (DFUVR+DFPAR+DFNIR))*(1.0-LANDICEALB)
    DIF = 0.0
    ULW = ALW + BLW*TS(:,N)

    LANDICECAP= (MAPL_RHOWTR*MAPL_CAPICE*LANDICEDEPTH)

    EVAPI   = LHF / MAPL_ALHS 
    DEVAPDT = LHD / MAPL_ALHS
    RADDN   = LWDNSRF + SWN 

    PERC  = 0.0 
    MELTI = 0.0 


    if(N==SNOW) then

       ITYPE = 9
       LAI   = 0.0
       GRN   = 0.0 
       MODISFAC = 1.0

       !*** have to do a transpose of these internals since their dimensions in SNOW_ALBEDO
       !*** are reversed
       WESNN = transpose(WESN)
       HTSNN = transpose(HTSN)
       SNDZN = transpose(SNDZ)
       !*** call new/shared routine to compute albedo 

       call    SNOW_ALBEDO(NT, NUM_SNOW_LAYERS, N_CONST_LANDICE4SNWALB, ITYPE, LAI, ZTH, & 
                   RHOFRESH, VISMAX, NIRMAX, SLOPE, &     !0.96, 0.68, 1.0,  & !
                   WESNN, HTSNN, SNDZN,        & ! snow stuff
                   LNDVR, LNDNR, LNDVF, LNDNF, & ! instantaneous snow-free albedos on tiles
                   SNOVR, SNONR, SNOVF, SNONF, & ! instantaneous snow albedos on tiles
                   RCONSTIT, UU, TS(:,SNOW), DRPAR, DFPAR & ! When only N_constit > 0 (oprional)
                   )

       VSUVR     = DRPAR + DRUVR
       VSUVF     = DFPAR + DFUVR
       SWNETSNOW = (1.-SNOVR)*VSUVR + (1.-SNOVF)*VSUVF + (1.-SNONR)*DRNIR + (1.-SNONF)*DFNIR
       RADDN     = LWDNSRF + SWNETSNOW 
       SWN       = SWNETSNOW
       if(associated(SNOWALB)) then
           where(FR(:,N) > 0.0)
               SNOWALB = SNOVR*AWTVDR + SNOVF*AWTVDF + SNONR*AWTIDR + SNONF*AWTIDF
           elsewhere
               SNOWALB  = MAPL_UNDEF
           endwhere
           where(ZTH < 1.e-6)
               SNOWALB = MAPL_UNDEF
           endwhere
       endif
    endif

    if(N==ICE) then
       do k=1,NT
          if(FR(k,N) > MINFRACSNO) then     
             call SOLVEICELAYER(NUM_ICE_LAYERS, DT, TICE(k,N,:), DZMAXI, 0,   &
                              MELTI(k), DTSS=DTS(k),  RUNOFF=PERC(k),                 &
                              lhturb=LHF(k),hlwtc=ULW(k),hsturb=SHF(k),raddn=RADDN(k),        &
                              dlhdtc=LHD(k),dhsdtc=SHD(k),dhlwtc=BLW(k),rain=RAIN(k),    &
                              rainrf=RAINRF(k),                                          & 
                              lhflux=LHFO(k),shflux=SHFO(k),hlwout=HLWO(k),evapout=EVAPO(k), &
                              ghflxice=ghflxice(k))
          else
             TICE(k,N,:) =  TICE(k,SNOW,:)
          endif
       enddo 
       TS(:,N)   =  TICE(:,N,1)
       if(associated(RUNOFF))   RUNOFF   = RUNOFF + FR(:,N) * PERC
    endif

    if(N==SNOW) then 
       LANDICELT  =  TICE(:,N,1) - MAPL_TICE
       do k=1,NT
#if 0 
          LATSD=LATS(K)*rad_to_deg
          LONSD=LONS(K)*rad_to_deg
          !if(abs(LATSD-0.700003698112E+02) < 1.e-3 .and. &
          !   abs(LONSD-(-0.539905136947E+02)) < 1.e-3 ) then
          !if(abs(LATSD-0.605467530483E+02) < 1.e-3 .and. &
          !   abs(LONSD-(-0.433431029954E+02)) < 1.e-3 ) then
          if(abs(LATSD-0.807870232172E+02) < 1.e-3 .and. &
             abs(LONSD-(-0.154247429558E+02)) < 1.e-3 ) then
            print*, 'PE = ', mype, ' tile = ',k 
          endif  
#endif
          TKSNO = condice 

          call SNOWRT( LONS(k), LATS(k),                                       &  ! in     [radians]  !!!
                   1,NUM_SNOW_LAYERS,MAPL_LANDICE,                             &  ! in    
                   MAXSNDZ, RHOFRESH, DZMAX,                                   &  ! in    
                   LANDICELT(k),ZONEAREA,TKGND,PRECIP(k),SNO(k),TA(k),DT,      &  ! in    
                   EVAPI(k),DEVAPDT(k),SHF(k),SHD(k),ULW(k),BLW(k),            &  ! in    
                   RADDN(k),ZC1,TOTDEPOS(k,:),                                 &  ! in    
                   WESN(k,:),HTSN(k,:),SNDZ(k,:), RCONSTIT(k,:,:),             &  ! inout    
                   HLWO(k), FROZFRAC(k,:),TPSN(k,:), RMELT(k,:),               &  ! out    
                   AREASC(k),FR(K,N),PERC(k),FHGND(k),                         &  ! out   
                   EVAPO(k),SHFO(k),LHFO(k),HCORR(k),ghflxsno(k),              &  ! out   
                   SNDZSC(k), WESNPREC(k), SNDZPREC(k),SNDZ1PERC(k),           &  ! out    
                   WESNPERC(k,:), WESNDENS(k,:), WESNREPAR(k,:), MLT(k),       &  ! out      
                   EXCS(k,:), DRHO0(k,:), WESNBOT(k), TKSNO, DTS(k)       )       ! out   
                                                                                 

          ! Snow impurities update
           if (N_CONST_LANDICE4SNWALB /= 0) then
              if(associated(IRDU001)) IRDU001(k,:) = RCONSTIT(k,:,1) 
              if(associated(IRDU002)) IRDU002(k,:) = RCONSTIT(k,:,2) 
              if(associated(IRDU003)) IRDU003(k,:) = RCONSTIT(k,:,3) 
              if(associated(IRDU004)) IRDU004(k,:) = RCONSTIT(k,:,4) 
              if(associated(IRDU005)) IRDU005(k,:) = RCONSTIT(k,:,5) 
              if(associated(IRBC001)) IRBC001(k,:) = RCONSTIT(k,:,6) 
              if(associated(IRBC002)) IRBC002(k,:) = RCONSTIT(k,:,7) 
              if(associated(IROC001)) IROC001(k,:) = RCONSTIT(k,:,8) 
              if(associated(IROC002)) IROC002(k,:) = RCONSTIT(k,:,9) 
           end if
           if(associated(RMELTDU001)) RMELTDU001(k) = RMELT(k,1) 
           if(associated(RMELTDU002)) RMELTDU002(k) = RMELT(k,2) 
           if(associated(RMELTDU003)) RMELTDU003(k) = RMELT(k,3) 
           if(associated(RMELTDU004)) RMELTDU004(k) = RMELT(k,4) 
           if(associated(RMELTDU005)) RMELTDU005(k) = RMELT(k,5) 
           if(associated(RMELTBC001)) RMELTBC001(k) = RMELT(k,6) 
           if(associated(RMELTBC002)) RMELTBC002(k) = RMELT(k,7) 
           if(associated(RMELTOC001)) RMELTOC001(k) = RMELT(k,8) 
           if(associated(RMELTOC002)) RMELTOC002(k) = RMELT(k,9) 

          if(associated(LWC ))then
               ZDEP = sum(SNDZ(k,:))
               if(sum(WESN(k,:)) > MINSWE) then
                   if(ZDEP <= LWCTOP) then
                      LWC(k) = sum(WESN(k,:)*(1.-FROZFRAC(k,:)))/sum(WESN(k,:))
                   else
                      KL  = 0
                      ZKL = 0.0 
                      do l=1,NUM_SNOW_LAYERS
                         ZKL = ZKL + SNDZ(k,l) 
                         if(ZKL > LWCTOP) then
                           KL = l
                           exit
                         endif
                      enddo 
                      ALPHA = 1.0 - (ZKL-LWCTOP)/SNDZ(k,KL)
                      LWC(k) = (sum(WESN(k,1:KL-1)*(1.-FROZFRAC(k,1:KL-1)))+ &
                                ALPHA*WESN(k,KL)*(1.-FROZFRAC(k,KL))) / &
                               (sum(WESN(k,1:KL-1))+ALPHA*WESN(k,KL))  
                   endif
               else
                   LWC(k) = 0.0
               endif            
          endif
          if(FR(K,N) < MINFRACSNO) then
             TICE(k,N,:) =  TICE(k,ICE,:)
          else
              call SOLVEICELAYER(NUM_ICE_LAYERS, DT, TICE(k,N,:), DZMAXI, 1,   &
                              MELTI(k),                    &
                              condsno=TKSNO(NUM_SNOW_LAYERS),       & 
                              !tsn=TPSN(k,NUM_SNOW_LAYERS),          & 
                              fhgnd=FHGND(k),          & 
                              sndz=SNDZ(k,NUM_SNOW_LAYERS)          &
                              )
              if(associated(RUNOFF)) RUNOFF(K)   = RUNOFF(K) + FR(K,N) * MELTI(K)
          endif   
       enddo   
       WESNSC = EVAPO
       !PERC = PERC + MELTI
       if(associated(RUNOFF))   RUNOFF   = RUNOFF + PERC
       TS(:,N) = TPSN(:,1)+MAPL_TICE
       if(associated(MELTWTRCONT )) MELTWTRCONT = sum(WESN*(1.-FROZFRAC),dim=2)
    endif

    DQS       = GEOS_QSAT(TS(:,N), PS, PASCALS=.TRUE.,RAMP=0.0) - QS(:,N)
    QS(:,N)   = QS(:,N) + DQS  

    LHF = LHFO
    SHF = SHFO
    ULW = HLWO 

    if(associated(EVAPOUT)) EVAPOUT = EVAPOUT + FR(:,N)*EVAPO
    if(associated(SUBLIM )) SUBLIM  = SUBLIM  + FR(:,N)*EVAPO
    if(associated(SHOUT  )) SHOUT   = SHOUT   + FR(:,N)*SHF
    if(associated(HLATN  )) HLATN   = HLATN   + FR(:,N)*LHF

    if(associated(DELTS )) DELTS = DELTS + DTS*CFT*FR(:,N)
    if(associated(DELQS )) DELQS = DELQS + DQS*CFQ*FR(:,N)
    if(associated(EVPICE)) EVPICE = EVPICE + FR(:,N)*LHF

    !if(associated(RUNOFF))   RUNOFF   = RUNOFF + FR(:,N) * PERC
    if(associated(IMELT ))   IMELT    = IMELT  + FR(:,N) * MELTI

    if(associated(SWNDSRF )) SWNDSRF = SWNDSRF + SWN * FR(:,N)
    if(associated(LWNDSRF )) LWNDSRF = LWNDSRF + (LWDNSRF - ULW) * FR(:,N)
    if(associated(HLWUP   )) HLWUP   = HLWUP +   ULW * FR(:,N)
    if(associated(DNICFLX )) DNICFLX = DNICFLX + DIF * FR(:,N)
    if(associated(GHSNOW  )) GHSNOW  = ghflxsno
    if(associated(ACCUM   )) ACCUM   = ACCUM - FR(:,N) * EVAPO  
    if(associated(MELTWTR )) MELTWTR = MELTWTR + FR(:,N) * MELTI  

    if(associated(TICE0   )) then
       do k=1,NT
         TICE0(k,:) =  TICE0(k,:)  + TICE(k,N,:) * FR(k,N)
       enddo
    endif 

    enddo  

    FR(:,ICE) = max(1.0-FR(:,SNOW), 0.0)

    if(associated(GHTSKIN )) GHTSKIN = ghflxsno*FR(:,SNOW) + ghflxice*FR(:,ICE)
    if(associated(ACCUM )) ACCUM      = ACCUM + PRECIP   
    if(associated(EMISS )) EMISS      = LANDICEEMISS

    if(associated(SNOWMASS)) SNOWMASS = sum(WESN,dim=2)
    if(associated(SNOMAS))   SNOMAS   = sum(WESN,dim=2)
    if(associated(SNOWDP))   SNOWDP   = sum(SNDZ,dim=2)
    if(associated(ASNOW))    ASNOW    = FR(:,SNOW)
    if(associated(SMELT ))   SMELT    = PERC
    if(associated(RAINRFZ )) RAINRFZ  = FR(:,ICE)  * RAINRF
    if(associated(MELTWTR )) MELTWTR  = MELTWTR + MLT  

! Update snow and landice albedos to anticipate
!   next radiation calculation
!-----------------------------------------------
          call MAPL_SunGetInsolation(LONS, LATS,      &
            ORBIT, ZTH, SLR, &
            INTV = TINT,     &
            currTime=CURRENT_TIME+DELT,  &
            RC=STATUS )
          VERIFY_(STATUS)

          ZTH = max(0.0,ZTH)

       ITYPE = 9

       call    SNOW_ALBEDO(NT, NUM_SNOW_LAYERS, N_CONST_LANDICE4SNWALB, ITYPE, LAI, ZTH, & 
                   RHOFRESH, VISMAX, NIRMAX, SLOPE,  &   ! 0.96, 0.68, 1.0,  & !
                   WESNN, HTSNN, SNDZN,        & ! snow stuff
                   LNDVR, LNDNR, LNDVF, LNDNF, & ! instantaneous snow-free albedos on tiles
                   SNOVR, SNONR, SNOVF, SNONF, & ! instantaneous snow albedos on tiles
                   RCONSTIT, UU, TS(:,SNOW), DRPAR, DFPAR & ! When only N_constit > 0 (oprional)
                   )

    if(associated(ALBVF )) ALBVF = FR(:,ICE)*LANDICEALB+FR(:,SNOW)*SNOVF
    if(associated(ALBVR )) ALBVR = FR(:,ICE)*LANDICEALB+FR(:,SNOW)*SNOVR
    if(associated(ALBNF )) ALBNF = FR(:,ICE)*LANDICEALB+FR(:,SNOW)*SNONF
    if(associated(ALBNR )) ALBNR = FR(:,ICE)*LANDICEALB+FR(:,SNOW)*SNONR


    if(associated(SNICEALB ))  then
       SNICEALB = FR(:,ICE)*LANDICEALB +                   & 
                  FR(:,SNOW)*(SNOVR*AWTVDR + SNOVF*AWTVDF  &
                            + SNONR*AWTIDR + SNONF*AWTIDF)
       where(ZTH < 1.e-6)
          SNICEALB = MAPL_UNDEF
       endwhere
    endif

! Copies for export
!------------------

    if(associated(TST  )) then
       TST = 0.0
       do N=1,NUM_SUBTILES
          TST = TST + TS(:,N)*FR(:,N)
       enddo
    end if

    if(associated(QST  )) then
       QST = 0.0
       do N=1,NUM_SUBTILES
          QST = QST + QS(:,N)*FR(:,N)
       end do
    end if

    if(associated(RHOSNOW  )) then
       RHOSNOW = 0.0
       do N=1,NUM_SNOW_LAYERS 
         !where(FR(:,SNOW) > 0.0 .and. SNDZ(:,N) > 0.0)
         where(sum(WESN,dim=2) > MINSWE)
            RHOSNOW(:,N) = WESN(:,N) / FR(:,SNOW) / SNDZ(:,N)
         elsewhere
            RHOSNOW(:,N) = MAPL_UNDEF 
         endwhere
       enddo
    end if

    if(associated(TSNOW  )) then
       TSNOW = 0.0
       do N=1,NUM_SNOW_LAYERS 
         where(FR(:,SNOW) > 0.0 .and. SNDZ(:,N) > 0.0)
            TSNOW(:,N) = TPSN(:,N)
         elsewhere
            TSNOW(:,N) = MAPL_UNDEF 
         endwhere
       enddo
    end if

    if(associated(WSNOW  )) then
       WSNOW = WESN
    end if

    if(associated(ZSNOW  )) then
       ZSNOW = SNDZ
    end if

    if(associated(DRHOS0  )) then
       DRHOS0 = DRHO0
    end if

    if(associated(WESNEX  )) then
       WESNEX = EXCS
    end if

    if(associated(WESNEXT  )) then
       WESNEXT = sum(EXCS,dim=2)
    end if

    if(associated(WESC )) then
       WESC =  WESNSC 
    end if

    if(associated(SDSC )) then
       SDSC =  SNDZSC / DT
    end if

    if(associated(WEPRE )) then
       WEPRE =  WESNPREC / DT
    end if

    if(associated(SDPRE )) then
       SDPRE =  SNDZPREC / DT
    end if

    if(associated(SD1PC )) then
       SD1PC =  SNDZ1PERC / DT
    end if

    if(associated(WEPERC )) then
       WEPERC =  WESNPERC / DT
    end if

    if(associated(WEDENS )) then
       WEDENS =  WESNDENS / DT
    end if

    if(associated(WEREP )) then
       WEREP =  WESNREPAR / DT
    end if

    if(associated(WEBOT )) then
       WEBOT =  WESNBOT / DT
    end if

    if(allocated (MLT)) deallocate(MLT , STAT=STATUS); VERIFY_(STATUS)              
    if(allocated (DTS)) deallocate(DTS , STAT=STATUS); VERIFY_(STATUS)              
    if(allocated (DQS)) deallocate(DQS , STAT=STATUS); VERIFY_(STATUS)              
    if(allocated (SHF)) deallocate(SHF , STAT=STATUS); VERIFY_(STATUS)              
    if(allocated (LHF)) deallocate(LHF , STAT=STATUS); VERIFY_(STATUS)              
    if(allocated (SHD)) deallocate(SHD , STAT=STATUS); VERIFY_(STATUS)              
    if(allocated (LHD)) deallocate(LHD , STAT=STATUS); VERIFY_(STATUS)              
    if(allocated (CFT)) deallocate(CFT , STAT=STATUS); VERIFY_(STATUS)              
    if(allocated (CFQ)) deallocate(CFQ , STAT=STATUS); VERIFY_(STATUS)              
    if(allocated (SWN)) deallocate(SWN , STAT=STATUS); VERIFY_(STATUS)              
    if(allocated (DIF)) deallocate(DIF , STAT=STATUS); VERIFY_(STATUS)              
    if(allocated (ULW)) deallocate(ULW , STAT=STATUS); VERIFY_(STATUS)                  
    if(allocated (PRECIP  )) deallocate(PRECIP  , STAT=STATUS); VERIFY_(STATUS)                  
    if(allocated (RAIN    )) deallocate(RAIN    , STAT=STATUS); VERIFY_(STATUS)                  
    if(allocated (RAINRF  )) deallocate(RAINRF  , STAT=STATUS); VERIFY_(STATUS)                  
    if(allocated (PERC    )) deallocate(PERC    , STAT=STATUS); VERIFY_(STATUS)                  
    if(allocated (MELTI   )) deallocate(MELTI   , STAT=STATUS); VERIFY_(STATUS)                  
    if(allocated (FROZFRAC)) deallocate(FROZFRAC, STAT=STATUS); VERIFY_(STATUS)
    if(allocated (TPSN    )) deallocate(TPSN    , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (AREASC  )) deallocate(AREASC  , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (HCORR   )) deallocate(HCORR   , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (ghflxsno)) deallocate(ghflxsno, STAT=STATUS); VERIFY_(STATUS)
    if(allocated (ghflxice)) deallocate(ghflxice, STAT=STATUS); VERIFY_(STATUS)
    if(allocated (HLWO    )) deallocate(HLWO    , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (EVAPO   )) deallocate(EVAPO   , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (LHFO    )) deallocate(LHFO    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (SHFO    )) deallocate(SHFO    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (ZTH     )) deallocate(ZTH     , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (SLR     )) deallocate(SLR     , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (EVAPI   )) deallocate(EVAPI   , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (DEVAPDT )) deallocate(DEVAPDT , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (ITYPE   )) deallocate(ITYPE   , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (LAI     )) deallocate(LAI     , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (GRN     )) deallocate(GRN     , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (MODISFAC)) deallocate(MODISFAC, STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (SNOVR    )) deallocate(SNOVR    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (SNONR    )) deallocate(SNONR    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (SNOVF    )) deallocate(SNOVF    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (SNONF    )) deallocate(SNONF    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (LNDVR    )) deallocate(LNDVR    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (LNDNR    )) deallocate(LNDNR    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (LNDVF    )) deallocate(LNDVF    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (LNDNF    )) deallocate(LNDNF    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (VSUVR    )) deallocate(VSUVR    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (VSUVF    )) deallocate(VSUVF    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (SWNETSNOW)) deallocate(SWNETSNOW, STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (RADDN    )) deallocate(RADDN    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (FHGND    )) deallocate(FHGND    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (DRHO0    )) deallocate(DRHO0    , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (EXCS     )) deallocate(EXCS     , STAT=STATUS); VERIFY_(STATUS) 
    if(allocated (WESNSC   )) deallocate(WESNSC   , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (SNDZSC   )) deallocate(SNDZSC   , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (WESNPREC )) deallocate(WESNPREC , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (SNDZPREC )) deallocate(SNDZPREC , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (SNDZ1PERC)) deallocate(SNDZ1PERC, STAT=STATUS); VERIFY_(STATUS)
    if(allocated (WESNBOT  )) deallocate(WESNBOT  , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (LANDICELT)) deallocate(LANDICELT, STAT=STATUS); VERIFY_(STATUS)
    if(allocated (RCONSTIT )) deallocate(RCONSTIT , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (TOTDEPOS )) deallocate(TOTDEPOS , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (RMELT    )) deallocate(RMELT    , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (WESNN    )) deallocate(WESNN    , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (HTSNN    )) deallocate(HTSNN    , STAT=STATUS); VERIFY_(STATUS)
    if(allocated (SNDZN    )) deallocate(SNDZN    , STAT=STATUS); VERIFY_(STATUS)

!  All done                          
!-----------                         

    RETURN_(ESMF_SUCCESS)            

  end subroutine LANDICECORE

!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


     subroutine SOLVEICELAYER(NICE, dts, TICE, ICEDZ, UPPER_BND,     &
                              MELT, DTSS,  RUNOFF,                   &
                              lhturb,hlwtc,hsturb,raddn,             &
                              dlhdtc,dhsdtc,dhlwtc,rain,rainrf,      &
                              lhflux,shflux,hlwout,evapout,          &
                              condsno, fhgnd, sndz, ghflxice      )

     implicit none 

     integer, intent(in)  :: NICE
     real,    intent(in    ) :: DTS

     real,    intent(inout ) :: TICE(NICE)
     real,    intent(in    ) :: ICEDZ(NICE)
     integer, intent(in    ) :: UPPER_BND
     real,    intent(out   ) :: MELT
     !  UPPER_BND == 0
     real, optional, intent(out) ::  DTSS
     real, optional, intent(out) ::  RUNOFF
     real, optional, intent(in ) ::  lhturb,hlwtc,hsturb,raddn
     real, optional, intent(in ) ::  dlhdtc,dhsdtc,dhlwtc
     real, optional, intent(in ) ::  rain
     real, optional, intent(out) ::  rainrf
     real, optional, intent(out) ::  lhflux,shflux,hlwout,evapout 
     real, optional, intent(out) ::  ghflxice
     !  UPPER_BND == 1
     real, optional, intent(in ) ::  condsno, fhgnd, sndz 

     ! Locals
     real :: melti,frrain,dtr,tsx,mass,snowd,rainf,denom,alhv,hcorr,            &
           enew,eold,tdum,fnew,tnew,icedens,densfac,hnew                      
     integer  ::  i     
     real, dimension(size(TICE)  ) :: tpsn
     real, dimension(size(TICE)  ) :: dtc,q,cl,cd,cr
     real, dimension(size(TICE)+1) :: fhsn,df


     
      df     = 0.
      dtc    = 0.
      fhsn   = 0.

      MELT   = 0.

       
      if(UPPER_BND == 0) then
         rainrf = 0.0
         RUNOFF = 0.0
      endif

      tpsn = TICE - tf

      alhv   = alhe + alhm                            !randy

      fhsn(NICE+1) = 0.0  
      df(NICE+1)   = 0.0 

!**** Calculate heat fluxes between snow layers.

      do i=2,NICE
         df(i)   =  -condice/((ICEDZ(i-1)+ICEDZ(i))*0.5)
         fhsn(i) =  df(i)*(tpsn(i-1)-tpsn(i))
      enddo

      if(present(ghflxice)) ghflxice = fhsn(2)

      if(UPPER_BND == 0) then

!**** Initial estimate of net surface flux & its change with Tc
        fhsn(1) = lhturb + hsturb + hlwtc - raddn
        df(1)   = -(dlhdtc + dhsdtc + dhlwtc)

      else
        df(1)   = -sqrt(condice*condsno)/((ICEDZ(1)+sndz)*0.5)
        !fhsn(1) = df(1)*(TSN - tpsn(1)) 
        fhsn(1) = fhgnd
      endif 

!**** Prepare array elements for solution & coefficient matrices.
!**** Terms are as follows:  left (cl), central (cd) & right (cr)
!**** diagonal terms in coefficient matrix & solution (q) terms.

        do i=1,NICE

           cl(i) = df(i)
           cd(i) = cpw*RHOICE*ICEDZ(i)/dts - df(i) - df(i+1)
           cr(i) = df(i+1)
           q(i)  = fhsn(i+1)-fhsn(i)

        enddo

        cl(1)    = 0.
        cr(NICE) = 0.


!**** Solve the tri-diagonal matrix for implicit change in Tc.

        call TRID(dtc,cl,cd,cr,q,NICE)

!**** Check temperature changes for passages across critical points,i.e.
!**** If implicit change has taken layer past melting/freezing, correct.

       do i=1,NICE
          if(tpsn(i)+dtc(i) > 0.) then
             melti  = (tpsn(i)+dtc(i))*cpw*MAPL_RHOWTR*ICEDZ(i)/MAPL_ALHF  
             MELT   = MELT   + melti 
             if(UPPER_BND == 0) then
                RUNOFF = RUNOFF  + melti 
             endif
             dtc(i) = -tpsn(i)
             tpsn(i) = tpsn(i) + dtc(i)
             if(i == 1) then 
                if(UPPER_BND == 0) then
                   RUNOFF = RUNOFF + rain * dts
                endif
             endif
          elseif(tpsn(i)+dtc(i) == 0.0) then
             tpsn(i) = tpsn(i) + dtc(i)
             if(i == 1) then 
                if(UPPER_BND == 0) then
                   RUNOFF = RUNOFF + rain * dts
                endif
             endif
          else  ! temp < 0, refreeze rain if any
             tpsn(i) = tpsn(i) + dtc(i)
             if(i == 1) then 
                if(UPPER_BND == 0) then
                 !*** only latent heat of rain is used to raise ice temp.  
                 !*** since AGCM assumes rain has 0 heat content
                 dtr = rain*dts*alhm/(RHOICE*cpw*ICEDZ(i))   
                 if(tpsn(i)+dtr > 0.0) then
                   frrain = max(dtr-(-tpsn(i))/dtr, 1.)
                   dtr = -tpsn(i)
                 else  
                   frrain = 0.0
                 endif   
                 tpsn(i) = tpsn(i) + dtr
                 dtc(i) = dtc(i) + dtr
                 RUNOFF = RUNOFF + frrain * rain * dts
                 rainrf = rainrf + (1.-frrain) * rain * dts 
                endif
             endif
          endif
          if(UPPER_BND == 0) then
            if(i == 1) then
               lhflux = lhturb + dlhdtc*dtc(1)
               shflux = hsturb + dhsdtc*dtc(1)
               hlwout = hlwtc  + dhlwtc*dtc(1)
               evapout = lhflux/alhv
            endif
          endif
       enddo
 
       MELT   = MELT  / dts  

       if(present(RUNOFF)) RUNOFF = RUNOFF / dts  
       if(present(rainrf)) rainrf = rainrf / dts  

       if(present(dtss)) dtss = dtc(1) 

       TICE = tpsn + tf

     end subroutine SOLVEICELAYER  

end subroutine RUN2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_LandiceGridCompMod

