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
  use StieglitzSnow, only: snowrt, SNOW_ALBEDO, TRID, N_CONSTIT, &
       NUM_DUDP, NUM_DUSV, NUM_DUWT, NUM_DUSD, &
       NUM_BCDP, NUM_BCSV, NUM_BCWT, NUM_BCSD, &
       NUM_OCDP, NUM_OCSV, NUM_OCWT, NUM_OCSD, &
       NUM_SUDP, NUM_SUSV, NUM_SUWT, NUM_SUSD, &
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
  !real,    parameter :: RHOMA    = 500.        ! kg/m^3  maximum snow density
  real,    parameter :: RHOICE   = 917.        ! kg/m^3  pure ice density
  real,    parameter :: MINSWE   = 0.013       ! kg/m^2  min SWE to avoid immediate melt
  real,    parameter :: MAXSNDZ  = 15.0        ! m
  real,    parameter :: ZERO     = 0.
  real,    parameter :: ONE      = 1.
  real,    parameter :: BIG      = 1.e10
  real,    parameter :: cpw      = 2065.22    !  @ 0 C [J/kg/K]
  real,    parameter :: condice  = 2.25       !  @ 0 C [W/m/K] 
  real,    parameter :: MINFRACSNO = 1.e-20   !  mininum sno/ice fraction for
                                              !  heat diffusion of ice layers to take effect
  real,    parameter :: LWCTOP     = 1.       !  top thickness to compute LWC. 1m taken from
                                              !  Fettweis et al 2011  
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
        LONG_NAME          = 'surface_outgoing_longwave_flux',&
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
        LONG_NAME          = 'linearization_of_surface_upwelling_longwave_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'BLW',                               &
        LONG_NAME          = 'linearization_of_surface_upwelling_longwave_flux', &
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
        LONG_NAME          = 'surface_downwelling_longwave_flux', &
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
  integer, pointer                    :: TILETYPES(:)
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
         TILETYPES = TILETYPES,                  &     
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
   real                           :: DUM
   real,  dimension(NUM_SNOW_LAYERS)    :: TKSNO

   real                           :: WSS

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

    DUM       = 0.0

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
    WSS       = 0.0

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
                   RHOFRESH, 0.96, 0.68, 1.0,  & !
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

          call SNOWRT(1,NUM_SNOW_LAYERS,TILETYPES(k),                          &
                   LANDICELT(k),ZONEAREA,TKGND,PRECIP(k),SNO(k),TA(k),DT,      &
                   EVAPI(k),DEVAPDT(k),SHF(k),SHD(k),ULW(k),BLW(k),            &
                   DUM,HLWO(k),RADDN(k),ZC1,TOTDEPOS(k,:),WSS,                 &
                   WESN(k,:),HTSN(k,:),SNDZ(k,:),                              &
                   FROZFRAC(k,:),TPSN(k,:), RCONSTIT(k,:,:), RMELT(k,:),       & 
                   AREASC(k),FR(K,N),PERC(k),FHGND(k),                         &
                   EVAPO(k),SHFO(k),LHFO(k),HCORR(k),ghflxsno(k),              &
                   SNDZSC(k), WESNPREC(k), SNDZPREC(k),SNDZ1PERC(k),           & 
                   WESNPERC(k,:), WESNDENS(k,:), WESNREPAR(k,:), MLT(k),       &   
                   EXCS(k,:), DRHO0(k,:), WESNBOT(k), TKSNO, DTS(k),           &
                   MAXSNDZ, RHOFRESH, DZMAX)

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

    if(associated(GHTSKIN )) GHTSKIN = (-1.0)*(ghflxsno*FR(:,SNOW) + ghflxice*FR(:,ICE)) 
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
                   RHOFRESH, 0.96, 0.68, 1.0,  & !
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

#if 0
     subroutine snowrt(N_sm, N_snow,                           &
           t1,area,tkgnd,precip,snowf,ts,dts, eturb,dedtc,hsturb,dhsdtc,       &
           hlwtc,dhlwtc,desdtc,hlwout,raddn,zc1,                               &
           wesn,htsnn,sndz,fices,tpsn,                                         &
           areasc,areasc0,pre,fhgnd,evap,shflux,lhflux,hcorr                   &
            !, wesnsc, sndzsc, wesnprec, sndzprec,  sndz1perc                   &   
            , sndzsc, wesnprec, sndzprec,  sndz1perc                           &   
            , wesnperc, wesndens, wesnrepar, mltwtr                            &   
            , excs, drho0, wesnbot, tksno, dtss, maxsndepth, rhofs,   &
              pid, ktile)

!*********************************************************************
! AUTHORS:  M. Stieglitz, M. Suarez, R. Koster & S. Dery.
! VERSION:  2003b - This version last updated:  05/30/03.
!*********
! INPUTS:
!*********
!  t1     : Temperature of catchment zones  [C]
!  ts     : Air temperature [K]
!  area   : Fraction of snow-free area in each catchment zone [0-1]
!  precip : Precipitation (Rain+snowfall) [kg/m^2/s == mm/s]
!  snowf  : Snowfall per unit area of catchment [kg/m^2/s == mm/s]
!  dts    : Time step  [s]
!  eturb  : Evaporation per unit area of snow [kg/m^2/s == mm/s]
!  dedtc  : d(eturb)/d(ts) [kg/m^2/s/K]
!  hsturb : Sensible heat flux per unit area of snow  [W/m^2]
!  dhsdtc : d(hsturb)/d(ts)  [W/m^2/K]
!  hlwtc  : Emitted IR per unit area of snow  [W/m^2]
!  dhlwtc : d(hlwtc)/d(ts)  [W/m^2/K]
!  raddn  : Net solar + incident terrestrial per unit area of snow [W/m^2]
!  tkgnd  : Thermal diffusivity of soil in catchment zones [W/m/K]
!  zc1    : Half-thickness of top soil layer [m]
!***  Bin Zhao added *************************************************
!  maxsndepth :  Maximum snow depth beyond which snow gets thrown away
!*********
! UPDATES:
!*********
!  wesn   : Layer water contents per unit area of catchment [kg/m^2]
!  htsnn  : Layer heat contents relative to liquid water at 0 C [J/m^2]
!  sndz   : Layer depths [m]
!*********
! OUTPUTS: 
!*********
!  tpsn   : Layer temperatures [C]
!  fices  : Layer frozen fraction [0-1]
!  areasc : Areal snow coverage at beginning of step [0-1]
!  areasc0: Areal snow coverage at end of step [0-1]
!  pre    : Liquid water flow from snow base [kg/m^2/s]
!  fhgnd  : Heat flux at snow base at catchment zones  [W/m^2]
!  hlwout : Final emitted IR flux per unit area of snow [W/m^2]
!  lhflux : Final latent heat flux per unit area of snow [W/m^2]
!  shflux : Final sensible heat flux per unit area of snow   [W/m^2]
!  evap   : Final evaporation per unit area of snow   [kg/m^2/s]
!***  Bin Zhao added *************************************************
!  sndzsc    :  top layer thickness change due to sublimation/condensation
!  wesnprec  :  top layer water content change due to precip (different from precip itself)
!  sndzprec  :  top layer thickness change due to precip 
!  sndz1perc :  top layer thickness change due to percolation
!  wesnperc  :  layer water content change due to percolation
!  wesndens  :  layer water content change due to densification
!  wesnrepar :  layer water content change due to relayer
!  mltwtr    :  total melt water production rate
!  excs      :  frozen part of water content from densification excess
!  drho0     :  layer density change due to densification
!  wesnbot   :  excessive water content due to thickness exceeding maximum depth
!  tksno     :  layer conductivity
!  dtss      :  top layer temperature change
!*********************************************************************
! NOTA:  By convention, wesn is representative for a catchment area
! equal to 1 whereas sndz is relative to the area covered by snow only.
!*********************************************************************


      implicit none

!      real, parameter :: lhv    = 2.4548E6 !  2.5008e6   !  @ 0 C [J/kg]
!      real, parameter :: lhs    = 2.8368E6 !  2.8434e6 !  @ 0 C [J/kg]
!      real, parameter :: lhf    = (lhs-lhv)  !  @ 0 C [J/kg]
      
!rr      real, parameter :: cpw_liquid = 4185. ! [J/kg/K]

!      real, parameter :: tfrz   = 273.16     !  @ 0 C [K]
!      real, parameter :: rhofs  = 150.       !  [kg/m^3]
!      real, parameter :: rhoma  = 500.       !  [kg/m^3]
!      real, parameter :: rhow   = 1000.      !  [kg/m^3]
!      real, parameter :: wemin  = 13.        !  [kg/m^2]
      real, parameter :: snfr   = 0.01       !  holding capacity
      real, parameter :: small  = 1.e-6      !  small number
!      integer, parameter :: nlay = 3         !  number of layers
!      integer, parameter :: N_sm   = 3         !  number of zones
!      real   , parameter :: MIN_SNOW_MASS = .013 ! kg/M**2 equiv to 0.1% area
      integer, parameter ::  N_constit = 1

 
      integer, intent(in)  :: N_sm, N_snow

      real,    intent(in ) :: t1(N_sm),area(N_sm),tkgnd(N_sm)
      real,    intent(in ) :: ts,precip,snowf,dts,dedtc,raddn,hlwtc
      real,    intent(in ) :: dhsdtc,desdtc,dhlwtc,eturb,hsturb,zc1
      real,    intent(inout):: wesn(N_snow),htsnn(N_snow),sndz(N_snow)
      real,    intent(out) :: tpsn(N_snow),fices(N_snow),fhgnd(N_sm)
      real,    intent(out) :: hlwout,lhflux,shflux,areasc0,evap,areasc,pre

      real,    intent(out) ::  wesnprec
      !real,    intent(out) :: wesnsc, wesnprec
      real,    intent(out) :: sndzsc, sndzprec
      real,    intent(out) :: sndz1perc
      real,    intent(out) :: wesnperc(N_snow)
      real,    intent(out) :: wesndens(N_snow)
      real,    intent(out) :: wesnrepar(N_snow)
      real,    intent(out) :: mltwtr
      real,    intent(out) :: excs(N_snow)
      real,    intent(out) :: drho0(N_snow)
      real,    intent(out) :: wesnbot
      real,    intent(out) :: tksno(N_snow)
      real,    intent(out) :: dtss
      real,    intent(in)  :: maxsndepth
      real,    intent(in)  :: rhofs
      integer, intent(in)  :: pid, ktile

!Locals
      real :: tsx, mass,snowd,rainf,denom,alhv,lhturb,dlhdtc,hcorr,            &
           enew,eold,tdum,fnew,tnew,icedens,densfac,hnew,scale,t1ave,          &
           flxnet,fdum,dw,waterin,waterout,snowin,snowout, mtwt,              &
           waterbal,precision,flow,term,dz,w(0:N_snow),HTSPRIME
      real :: excsdz, excswe, sndzsum, melti
      real, dimension(size(wesn)  ) :: cmpc,dens
      real, dimension(size(wesn)  ) :: tksn
      real, dimension(size(wesn)  ) :: dtc,q,cl,cd,cr
      real, dimension(size(wesn)+1) :: fhsn,df
      real, dimension(size(wesn)  ) :: htest,ttest,ftest

      logical, dimension(size(wesn)  ) ::  ice1,tzero, ice10,tzero0
      real                          :: topthick
      real, dimension(size(wesn)-1) :: thickdist
      real,  dimension(size(wesn), N_constit):: rconstit
      integer :: i,izone
      logical :: logdum,kflag

       snowd = sum(wesn)
       snowin = snowd

!rr   correction for "cold" snow
       tsx   = min(ts-tf,0.)*cpw
       
!rr   correction for heat content of rain
!rr       tsx_rain = max(ts-tf,0.)*cpw_liquid
       
       df     = 0.
       dtc    = 0.
       tpsn   = 0.
       fices  = 0.
       areasc = 0.
       areasc0= 0.
       pre    = 0.
       fhgnd  = 0.
       hlwout = 0.
       shflux = 0.
       lhflux = 0.
       evap   = 0.
       excs   = 0.
       hcorr  = 0.
       dens   = rhofs
       rainf  = precip - snowf   ! [kg/m^2/s]

       !wesnsc = 0.
       sndzsc = 0.
       wesnprec = 0.
       sndzprec = 0.
       sndz1perc = 0.
       wesnperc = 0.
       wesndens = 0.
       wesnrepar = 0.
       wesnbot = 0.
       tksno   = condice
       dtss    = 0. 
       excswe  = 0.


       rconstit = 0.0 

       if(snowd <= MINSWE) then ! no snow
!         Assume initial (very small) snow water melts; new area is based
!         on new snowfall

          areasc = min(snowd/wemin,1.)
          areasc0 = 0.
          pre = snowd/dts + areasc*rainf
          wesn  = 0.
          hcorr = hcorr + raddn*areasc + sum(htsnn)/dts
          htsnn = 0.
          sndz  = 0.
          mltwtr = snowd/dts

          if(snowf > 0.) then   ! only initialize with non-liquid part of precip
                                ! liquid part runs off (above)

             wesn    = snowf*dts/float(N_snow)  
             htsnn   = (tsx-alhm)*wesn
             areasc0 = min((snowf*dts)/wemin,1.)
             sndz = wesn/rhofs
!             hcorr  = hcorr - (tsx-alhm)*snowf    ! randy
             hcorr  = hcorr - tsx*snowf           ! randy
             call FindTargetThickDist(N_snow, sndz, topthick, thickdist)
             !call relayer(N_snow, htsnn, wesn, sndz)
             call relayer2(N_snow, N_constit, topthick, thickdist, &
                            htsnn, wesn, sndz, rconstit, pid, ktile)
             call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)
             endif

        return ! if there was no snow at start of time step

       endif


      call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)
      mtwt = sum(wesn*(1.-fices)) 

!**** Determine the fractional snow coverage

       areasc = min(snowd/wemin,1.)

!**** Set the mean density & diffusivity of the layers

       do i=1,N_snow
         if(sndz(i) > 0) dens(i) = max(wesn(i)/(areasc*sndz(i)),rhofs)
       enddo
       tksn  = 3.2217e-06*dens**2
       tksno = tksn

!**** Determine temperature & frozen fraction of snow layers

         call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)

!**** Calculate the ground-snow energy flux at 3 zones

       denom = 1./(sndz(N_snow)*0.5-zc1)
       fhgnd = -sqrt(tkgnd*tksn(N_snow))*area*denom*(tpsn(N_snow)-t1)
       fhsn(N_snow+1) = sum(fhgnd)
       do i=1,N_sm
        df(N_snow+1)=df(N_snow+1)-sqrt(tkgnd(i)*tksn(N_snow))*area(i)*denom
       enddo


!**** Ensure against excessive heat flux between ground and snow:
!**** if heat flux found to cause the lowest snow layer temperature
!**** to "overshoot" (e.g. to become higher than the ground temperature
!**** when it had been lower), reduce the heat flux.  If the lowest 
!**** snow layer starts off at zero and the new temperature is greater
!**** than zero, reduce the heat flux to melt only half of the lowest
!**** layer snow.
!**** 
      t1ave=sum(t1*area)/sum(area)
      htest=htsnn
      htest(N_snow)=htest(N_snow)+fhsn(N_snow+1)*dts*areasc

      call get_tf_nd(N_snow, htest, wesn, ttest, ftest)

      scale=1.
      if((t1ave-tpsn(N_snow))*(t1ave-ttest(N_snow)) .lt. 0.) then
         scale=0.5*(tpsn(N_snow)-t1ave)/(tpsn(N_snow)-ttest(N_snow))
         endif
      if(tpsn(N_snow) .eq. 0. .and. ttest(N_snow) .gt. 0. .and.                &
                abs(fhsn(N_snow+1)) .gt. 1.e-10) then
         scale=(-0.5*htsnn(N_snow)/(dts*areasc))/fhsn(N_snow+1)
         endif

      fhsn(N_snow+1)=fhsn(N_snow+1)*scale
         df(N_snow+1)=df(N_snow+1)*scale
         fhgnd=fhgnd*scale


!**** Calculate heat fluxes between snow layers.

       do i=2,N_snow
         df(i) =  -sqrt(tksn(i-1)*tksn(i))/((sndz(i-1)+sndz(i))*0.5)
         fhsn(i)= df(i)*(tpsn(i-1)-tpsn(i))
       enddo
         
 
!**** Effective heat of vaporization includes bringing snow to 0 C

        alhv   = alhe + alhm                            !randy
!        alhv   = alhe + fices(1)*alhm + tpsn(1)*cpw    !randy

!**** Initial estimate of latent heat flux change with Tc

        lhturb = alhv*eturb
        dlhdtc = alhv*dedtc

!**** Initial estimate of net surface flux & its change with Tc

        fhsn(1) = lhturb + hsturb + hlwtc - raddn
        df(1)   = -(dlhdtc + dhsdtc + dhlwtc)

!**** Prepare array elements for solution & coefficient matrices.
!**** Terms are as follows:  left (cl), central (cd) & right (cr)
!**** diagonal terms in coefficient matrix & solution (q) terms.

        do i=1,N_snow

         call get_tf0d(htsnn(i),wesn(i),tdum,fdum, ice1(i),tzero(i))

         if(ice1(i)) then
           cl(i) = df(i)
           cd(i) = cpw*wesn(i)/dts - df(i) - df(i+1)
           cr(i) = df(i+1)
           q(i)  = fhsn(i+1)-fhsn(i)
         else
           cl(i) = 0.
           cd(i) = 1.
           cr(i) = 0.
           q(i)  = 0.
         endif

        enddo

        cl(1)    = 0.
        cr(N_snow) = 0.

        do i=1,N_snow-1
          if(.not.ice1(i)) cl(i+1) = 0.
        enddo

        do i=2,N_snow
          if(.not.ice1(i)) cr(i-1) = 0.
        enddo


!**** Solve the tri-diagonal matrix for implicit change in Tc.

        call TRID(dtc,cl,cd,cr,q,N_snow)

!**** Check temperature changes for passages across critical points,i.e.
!**** If implicit change has taken layer past melting/freezing, correct.

       do i=1,N_snow
          if(tpsn(i)+dtc(i) > 0. .or. htsnn(i)+wesn(i)*cpw*dtc(i) > 0.) then
             dtc(i)=-tpsn(i)
             endif
          if(.not.ice1(i)) dtc(i)=0.
          enddo
 
!**** Further adjustments; compute new values of h associated with
!**** all adjustments.

       eold=sum(htsnn)

       do i=1,N_snow

!**** Quick check for "impossible" condition:

          if(.not.tzero(i) .and. .not.ice1(i)) then
             write(*,*) 'bad snow condition: fice,tpsn =',fices(i),tpsn(i)
             stop
             endif

!****  Condition 1: layer starts fully frozen (temp < 0.)

          if(.not.tzero(i)) then
             tnew=tpsn(i)+dtc(i)
             fnew=1.

             endif

!****  Condition 2: layer starts with temp = 0, fices < 1.
!      Corrections for flxnet calculation: Koster, March 18, 2003.

          if(.not.ice1(i)) then
             tnew=0.
             if(i==1) flxnet= fhsn(i+1)+df(i+1)*(dtc(i)-dtc(i+1))              &
                   -fhsn(i)-df(i)*dtc(i)
             if(i > 1 .and. i < N_snow) flxnet=                                &
                    fhsn(i+1)+df(i+1)*(dtc(i)-dtc(i+1))                        &
                   -fhsn(i)-df(i)*(dtc(i-1)-dtc(i))
             if(i==N_snow) flxnet=fhsn(i+1)+df(i+1)*dtc(i)                     &
                   -fhsn(i)-df(i)*(dtc(i-1)-dtc(i))
             HTSPRIME=HTSNN(I)+AREASC*FLXNET*DTS
             call get_tf0d(HTSPRIME,wesn(i),  tdum,fnew,logdum,logdum)
             fnew=amax1(0.,  amin1(1.,  fnew))

             endif

!****  Condition 3: layer starts with temp = 0, fices = 1.
!      Corrections for flxnet calculation: Koster, March 18, 2003.

          if(ice1(i) .and. tzero(i)) then
             if(dtc(i) < 0.) then
                tnew=tpsn(i)+dtc(i)
                fnew=1.
                endif
             if(dtc(i) >= 0.) then
                tnew=0.
                if(i==1) flxnet=fhsn(i+1)+df(i+1)*(dtc(i)-dtc(i+1))            &
                       -fhsn(i)-df(i)*dtc(i)
                if(i > 1 .and. i < N_snow) flxnet=                             &
                    fhsn(i+1)+df(i+1)*(dtc(i)-dtc(i+1))                        &
                   -fhsn(i)-df(i)*(dtc(i-1)-dtc(i))
                if(i==N_snow) flxnet=fhsn(i+1)+df(i+1)*dtc(i)                  &
                   -fhsn(i)-df(i)*(dtc(i-1)-dtc(i))

                HTSPRIME=HTSNN(I)+AREASC*FLXNET*DTS
                call get_tf0d(HTSPRIME,wesn(i),   tdum,fnew,logdum,logdum)
                fnew=amax1(0.,  amin1(1.,  fnew))
                endif
             endif

!**** Now update heat fluxes & compute sublimation or deposition.

         if(i == 1) then
             dtss   = dtc(1)
             lhflux = lhturb + dlhdtc*dtc(1)
             shflux = hsturb + dhsdtc*dtc(1)
             hlwout = hlwtc  + dhlwtc*dtc(1)
             evap = lhflux/alhv
             dw = -evap*dts*areasc
             if(-dw > wesn(1) ) then
                dw = -wesn(1)
                evap = -dw/(dts*areasc)
!                shflux=shflux+(lhflux-evap*alhv)
                hcorr=hcorr+(lhflux-evap*alhv)*areasc
                lhflux=evap*alhv
                endif
             wesn(1)  = wesn(1) + dw
             denom = 1./dens(1)
             if(dw > 0.) denom = 1./rhoma
             sndz(1) = sndz(1) + dw*denom
             !wesnsc = dw
             sndzsc = dw*denom
             endif

         if(i == N_snow) then
             do izone=1,N_sm
                fhgnd(izone)=fhgnd(izone)+area(izone)*df(N_snow+1)*dtc(N_snow)
                enddo
             endif

!**** Now update thermodynamic quantities.

          htsnn(i)=(cpw*tnew-fnew*alhm)*wesn(i)
          tpsn(i) = tnew    
          fices(i)= fnew
        enddo

!**** Store excess heat in hcorr.

       enew=sum(htsnn)
       hcorr=hcorr-((enew-eold)/dts+areasc*(lhflux+shflux+hlwout-raddn)        &
                      -areasc*(fhsn(N_snow+1)+df(N_snow+1)*dtc(N_snow))        &
                      )

       call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)

       mltwtr = max(0., sum(wesn*(1.-fices)) - mtwt)
       mltwtr = mltwtr / dts

!rr!**** Add rainwater and snow at ts., bal. budget with shflux.
!rr   (tried and failed 19 Jun 2003, reichle)
!rr
!rr       wesn (1) = wesn (1) + (rainf*areasc+snowf)*dts
!rr       htsnn(1) = htsnn(1) + (tsx -alhm)*(snowf*dts) + tsx_rain*rainf*dts
!rr       sndz (1) = sndz (1) + (snowf/rhofs)*dts
!rr       !  shflux   = shflux   + tsx*snowf                   ! randy
!rr       hcorr   = hcorr   - (tsx-alhm)*snowf - tsx_rain*rainf ! randy


!**** Add rainwater at 0 C, snow at ts., bal. budget with shflux.

       wesn (1) = wesn (1) + (rainf*areasc+snowf)*dts
       htsnn(1) = htsnn(1) + (tsx -alhm)*(snowf*dts)
       sndz (1) = sndz (1) + (snowf/rhofs)*dts
       wesnprec = (rainf*areasc+snowf)*dts
       sndzprec = (snowf/rhofs)*dts
!       shflux  = shflux   + tsx*snowf          ! randy
!       hcorr   = hcorr   - (tsx-alhm)*snowf     ! randy
        hcorr   = hcorr   -  tsx*snowf          ! randy
       
       snowd=sum(wesn)

       call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)

!**** Move meltwater through the pack.
!**** Updated by Koster, August 27, 2002.

       pre = 0.
       flow = 0.

       wesnperc = wesn

       do i=1,N_snow

        if(flow > 0.) then
         wesn (i) =  wesn(i) + flow
         call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)  
      endif

        pre  = max((1.-fices(i))*wesn(i), 0.)
        flow = 0.

        if(snowd > wemin) then

          icedens=wesn(i)*fices(i)/(sndz(i)+1.e-20)
          densfac=amax1(0., amin1(1., icedens/rhofs))
          term=densfac*snfr*(sndz(i)*rhow-wesn(i)*fices(i))
      
          if(pre > term) then
            pre = min(pre - term, wesn(i))
            wesn(i) = wesn(i) - pre
            flow = pre
          endif
        else
          wesn(i) = wesn(i) - pre
          flow = pre
       endif

!**** Adjust top layer snow depth to get proper density values
!**** But limit this change for large throughflow (STEPH 06/19/03)

        if(i==1)then
          dz=min(flow/dens(i),0.5*sndz(i))
          sndz(i)=sndz(i)-dz
          sndz1perc = -dz
        endif
       enddo

       wesnperc = wesn - wesnperc

       pre = flow/dts
       snowd=sum(wesn)

!**** Update snow density by compaction (Pitman et al. 1991)

       excs = 0.
       mass = 0.
       w    = 0.
       drho0 = 0.

       wesndens = wesn

       if(snowd > wemin) then ! Compaction only after full coverage.

          do i=1,N_snow
             dens(i) = rhofs
             if(sndz(i)>0.) dens(i) = max(wesn(i)/(sndz(i)),rhofs)
          enddo

          drho0 = dens 
          
          cmpc    = exp(14.643 - (4000./min(tpsn+tf,tf))-.02*dens)

          do i=1,N_snow
             w(i) = wesn(i)
             mass = mass + 0.5*(w(i)+w(i-1))
             dens(i) = dens(i)*(1. + (dts*0.5e-7*9.81)*mass*cmpc(i))
             
!**** Clip densities below maximum value, adjust quantities accordingly
!**** while conserving heat & mass (STEPH 06/21/03).

             if(dens(i) > rhoma) then
                excs(i) = (dens(i)-rhoma)*sndz(i)
                wesn(i) = wesn(i) - excs(i)
                hnew = (cpw*tpsn(i)-fices(i)*alhm)*wesn(i)
                hcorr= hcorr+(htsnn(i)-hnew)/dts
                htsnn(i)= hnew
                dens(i) = rhoma
             endif
          enddo
          drho0 = dens - drho0
       endif
 
       wesndens = wesn - wesndens

       !pre  = pre + sum(excs)/dts
       pre  = pre + sum(excs*(1.-fices))/dts
       excs = excs * fices
       sndz = wesn/dens
                
       sndzsum = sum(sndz)
       if(sndzsum > maxsndepth) then
           excsdz  = sndzsum - maxsndepth
           excswe  = dens(N_snow) * excsdz 
           wesn(N_snow) = wesn(N_snow) - excswe
           hnew = (cpw*tpsn(N_snow)-fices(N_snow)*alhm)*wesn(N_snow)
           htsnn(N_snow)= hnew
           sndz(N_snow) = sndz(N_snow) - excsdz
           wesnbot = excswe
       endif

#if 0
      if(pid == TAR_PE .and. ktile == TAR_TILE)then
         write(*,*) 'BEFORE RELAYER *********************' 
         write(*,*) (wesn(i), i=1,N_snow)
         write(*,*) 'total swe = ', sum(wesn)
      endif 
#endif
       
!**** Restore layers to sigma values.
       
       wesnrepar = wesn

       do i=1,N_snow
         call get_tf0d(htsnn(i),wesn(i),tdum,fdum,ice10(i),tzero0(i))
         enddo

       call FindTargetThickDist(N_snow, sndz, topthick, thickdist)

       !call relayer(N_snow, htsnn, wesn, sndz)
       call relayer2(N_snow, N_constit, topthick, thickdist, &
                     htsnn, wesn, sndz, rconstit, pid, ktile)

       wesnrepar = wesn - wesnrepar

#if 0
      if(pid == TAR_PE .and. ktile == TAR_TILE)then
         write(*,*) 'AFTER RELAYER *********************' 
         write(*,*) (wesn(i), i=1,N_snow)
         write(*,*) 'total swe = ', sum(wesn)
      endif 
#endif

       call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)
       
!**** Check that (ice10,tzero) conditions are conserved through
!**** relayering process (or at least that (fices,tpsn) conditions don't 
!**** go through the (1,0) point); excess goes to hcorr.

       do i=1,N_snow
          kflag=.false.
          if(ice10(i).and.tzero0(i) .and.                                      &
             (fices(i) .ne. 1. .or. tpsn(i) .ne. 0.) ) kflag=.true.
          if(.not.ice10(i).and.tzero0(i) .and.                                 &
             (fices(i) .eq. 1. .and. tpsn(i) .lt. 0.) ) kflag=.true.
          if(ice10(i).and. .not.tzero0(i) .and.                                &
             (fices(i) .ne. 1. .and. tpsn(i) .eq. 0.) ) kflag=.true.

          if(kflag) then
             hnew=-alhm*wesn(i)
             hcorr=hcorr+(htsnn(i)-hnew)/dts
             htsnn(i)=hnew
             tpsn(i)=0.
             fices(i)=1.
             endif

          enddo


!**** Reset fractional area coverage.

       areasc0 = min(sum(wesn)/wemin,1.)
 
!**** Final check for water balance.

       waterin   = (rainf*areasc+snowf)*dts + max(dw,0.)
       waterout  = pre*dts - min(dw,0.)
       snowout   = sum(wesn) + sum(excs) + excswe
       waterbal  = snowin + waterin - waterout - snowout
       precision = snowout*small
       if((waterbal > precision).and.(waterbal > small)) then
#if 0
         write(*,*) 'Warning: Imbalance in snow water budget!', waterbal
         write(*,*) 'waterin   = ', waterin
         write(*,*) 'snowin    = ', snowin
         write(*,*) 'waterout  = ', waterout
         write(*,*) 'snowout   = ', snowout
         write(*,*) 'dw   = ', dw
         write(*,*) 'excswe  = ', excswe
         write(*,*) 'sum(excs)  = ', sum(excs)
         write(*,*) 'snowf*dts  = ', snowf*dts
         write(*,*) 'sum(wesn) = ', sum(wesn)
         write(*,*) (wesn(i), i=1,N_snow)
         write(*,*) 'pid = ', pid
         write(*,*) 'ktile = ', ktile
#endif
         stop
       endif
      
      return  !  end snow

      end subroutine snowrt

! **********************************************************************

      subroutine relayer(N_snow, htsnn, wesn, sndz)
      
      implicit none
      integer, intent(in) :: N_snow
      real, intent(inout) :: htsnn(N_snow),wesn(N_snow),sndz(N_snow)
      
      real, dimension(size(sndz),2)   :: ds
      real, dimension(size(sndz))     :: dovp
      real, dimension(size(sndz)+1)   :: sdold, sdnew
      real, dimension(size(sndz)+1)   :: sdoldr, sdnewr
      real, dimension(size(sndz)+1,2) :: h, s
      
      integer :: i, j, n
      integer :: kmax
      real    :: dzdiff, restthick, dzold
      integer, dimension(size(sndz)) :: mark

      logical :: lth_satisfy 
      
!      real, parameter :: dz1max = 0.08   ! [m]
!      real, parameter :: wemin  = 13.0   ! [kg/m2]
      real, parameter :: small  = 1.e-20 
      real :: areasc,dz
      
!**** Initialize some variables.
      
      h  = 0.
      s  = 0.
      ds = 0.
      dz = 0.
      
      areasc = min(sum(wesn)/wemin,1.)
      
!**** Compute specific heat & water contents of old layers.

      do i=1,N_snow
         if (sndz(i) > 0.) then
            h(i,1) = htsnn(i)/sndz(i)
            h(i,2) =  wesn(i)/sndz(i)
         endif
      enddo
      
!**** Obtain old & new layer thicknesses & boundaries.
      
      sdold = 0.
      sdnew = 0.
      sdoldr = 0.
      sdnewr = 0.
      
      do i=N_snow,1,-1
         sdold(i) = sdold(i+1) + sndz(i)
      enddo

      do i=1,N_snow
        sdoldr(i+1) = sdoldr(i) + sndz(i)
      enddo
      
      sndz = sdold(1)/float(N_snow)
      !sndz = sdold(N_snow+1)/float(N_snow)

      mark = 0
      do
        lth_satisfy = .true.
        do i=1,N_snow
          if(mark(i) == 0 .and. sndz(i) > dzmax(i)) then
            sndz(i)  = dzmax(i)
            mark(i)  = 1
            lth_satisfy = .false.
          endif
        enddo
        if(lth_satisfy) exit 
        dzdiff = 0.0
        do i=1,N_snow
            if(mark(i) == 1) then
                dzdiff = dzdiff + sndz(i)
            endif
        enddo
        restthick = (sdold(1)-dzdiff)/float(N_snow-sum(mark))
        !restthick = (sdold(N_snow+1)-dzdiff)/float(N_snow-sum(mark))
        do i=1,N_snow
           if(mark(i) == 0) then
               sndz(i) = restthick
           endif
        enddo
      enddo 

      !kmax = 0 
      !dzdiff = 0.0
      !do i=1,N_snow
      !   if(sndz(i) > dzmax(i)) then
      !      dzdiff   = dzdiff + sndz(i) - dzmax(i)
      !      sndz(i)  = dzmax(i)
      !      kmax = kmax + 1
      !   endif
      !enddo
      
      !sndz(kmax+1:) = (sdold(1)-sum(sndz(1:kmax)))/float(N_snow-kmax)
      
      do i=N_snow,1,-1
         sdnew(i) = sdnew(i+1) + sndz(i)
      enddo
      do i=1,N_snow
         sdnewr(i+1) = sdnewr(i) + sndz(i)
      enddo
      
!**** Since the snow boundary has moved, redistribute heat  
!     contents & water equivalents of old to new snow layers.

      dzold = 0.0 
                                                              
      do i=1,N_snow
         
         j = i
         dz=sdnew(i+1)-sdold(i+1)
         !dz=sdold(i+1)-sdnew(i+1)
         if(dz < 0.) j = i + 1
         if(dzold > sndz(i)) then
            call FindOverlap(N_snow, i, sdoldr, sdnewr, dz, dovp)
            do n=1,N_snow
              s(i+1,:) = s(i+1,:) + h(n,:) * dovp(n)
            enddo
         else
            s(i+1,:) = h(j,:)*dz
         endif
         dzold = dz
         ds(i,:)  = s(i,:) - s(i+1,:)
      enddo

      htsnn = htsnn + ds(:,1)
      wesn  = wesn  + ds(:,2) 
      
      if(sum(wesn) < wemin) sndz = sndz /(areasc + small)
      return
      
      end subroutine relayer

! **********************************************************************

      subroutine FindOverlap(N, m, sdoldr, sdnewr, dz, dovp)

      implicit none
      real, parameter     :: epsil = 1.e-6
      integer, intent(in) :: N, m
      real, intent(in), dimension(N+1) :: sdoldr, sdnewr
      real, intent(in)                 :: dz
      real, intent(out), dimension(N)  :: dovp

      real y1, y2
      integer i, j, k, ks, ke


      dovp = 0.0

      y1 = sdnewr(m+1)
      y2 = y1 + dz

      do k=1,N
         if((y1-sdoldr(k))>epsil) ks = k
         if((y2-sdoldr(k))>epsil) ke = k
      enddo

      dovp(ks) = sdoldr(ks+1) - y1
      y2 = y1
      do k=ks+1, ke
          y2 = y2+dovp(k-1)
          dovp(k) = sdoldr(k+1) - y2
      enddo
      !if(m==2) then
      !   write(*,*) (dovp(k), k=1,N)
      !endif
      end subroutine FindOverlap


! **********************************************************************

      subroutine FindTargetThickDist(N_snow, sndz, topthick, thickdist)

      integer, intent(in) :: N_snow
      real, intent(in)                       :: sndz(N_snow)
      real, intent(out)                      :: topthick
      real, intent(out), dimension(N_snow-1) :: thickdist

       
      real, dimension(N_snow)                :: sndzt
      real                                   :: totald, dzdiff, restthick
      integer                                :: i
      integer, dimension(N_snow)             :: mark
      logical                                :: lth_satisfy

      totald = sum(sndz)
      sndzt = totald/float(N_snow)

      mark = 0
      do
        lth_satisfy = .true.
        do i=1,N_snow
          if(mark(i) == 0 .and. sndzt(i) > dzmax(i)) then
            sndzt(i) = dzmax(i)
            mark(i)  = 1
            lth_satisfy = .false.
          endif
        enddo
        if(lth_satisfy) exit 
        dzdiff = 0.0
        do i=1,N_snow
            if(mark(i) == 1) then
                dzdiff = dzdiff + sndzt(i)
            endif
        enddo
        restthick = (totald-dzdiff)/float(N_snow-sum(mark))
        do i=1,N_snow
           if(mark(i) == 0) then
               sndzt(i) = restthick
           endif
        enddo
      enddo 

      topthick = sndzt(1)
      totald = totald - topthick
      do i=2,N_snow
        thickdist(i-1) = sndzt(i)/totald  
      enddo

      return

      end subroutine FindTargetThickDist  


! **********************************************************************

      subroutine relayer2(N_snow, N_constit, thick_toplayer, thickdist, & 
                           htsnn, wesn, sndz, rconstit, pid, ktile)
      
      implicit none
      integer, intent(in) :: N_snow, N_constit
      real, intent(in) :: thick_toplayer
      real, intent(in), dimension(N_snow-1) :: thickdist
      real, intent(inout) :: htsnn(N_snow),wesn(N_snow),sndz(N_snow)
      real, intent(inout) :: rconstit(N_snow,N_constit)
      integer, intent(in) :: pid, ktile 
      
      real, dimension(size(sndz),2+N_Constit) :: h, s
      
      integer :: i, j, k, ilow, ihigh
      
!      real, parameter :: dz1max = 0.08   ! [m]
!      real, parameter :: wemin  = 13.0   ! [kg/m2]
      real, parameter :: small  = 1.e-20 
      real :: areasc,dz

      real :: totalthick
      real, dimension(size(sndz)) :: thickness, tol_old, bol_old, tol_new,    &
                                         bol_new

!**** thick_toplayer: the assigned (final) thickness of the topmost layer (m)
!**** thickdist: the assigned (final) distribution of thickness in layers
!****        2 through N_snow, in terms of fraction
!**** h: array holding specific heat, water, and constituent contents
!**** s: array holding the total and final heat, water, and constit. contents
!**** ilow: first layer used in a particular relayering calculation
!**** ihigh: final layer used in a particular relayering calculation 
!**** totalthick: total thickness of layers 2 through N_snow
!**** thickness: array holding final thicknesses (m) of the snow layers
!**** tol_old(i): depth (from surface) of the top of layer i, before          &
!****               relayering
!**** bol_old(i): depth (from surface) of the bottom of layer i, before       &
!****               relayering
!**** tol_old(i): depth (from surface) of the top of layer i, after           &
!****               relayering
!**** bol_old(i): depth (from surface) of the bottom of layer i, after        &
!****               relayering


      thickness(1)=thick_toplayer

      totalthick=sum(sndz)-thick_toplayer
      do i=1,N_snow-1
        thickness(i+1)=thickdist(i)*totalthick
        enddo
      
#if 0
      if(pid == TAR_PE .and. ktile == TAR_TILE)then
         write(*,*) 'RELAYER2 A. *********************' 
         write(*,*) (sndz(i), i=1,N_snow)
          write(*,*) 'total sndz = ', sum(sndz)
         write(*,*) (thickness(i), i=1,N_snow)
          write(*,*) 'total thickness = ', sum(thickness)
      endif 
#endif
      
!**** Initialize some variables.
      
      h  = 0.
      s  = 0.
      dz = 0.
      
      areasc = min(sum(wesn)/wemin,1.)
      
!**** Compute specific heat & water contents of old layers.

      do i=1,N_snow
         if (sndz(i) > 0.) then
            h(i,1) = htsnn(i)/sndz(i)
            h(i,2) =  wesn(i)/sndz(i)
            do k=1,N_Constit
               h(i,2+k)=rconstit(i,k)/sndz(i)
               enddo
            endif
         enddo
      
!**** Determine old and new boundary depths (cumulative from top)
!**** (tol refers to "top of layer", bol refers to "bottom of layer"

      tol_old(1)=0.
      bol_old(1)=sndz(1)
      tol_new(1)=0.
      bol_new(1)=thickness(1)

      do i=2,N_snow
        tol_old(i)=bol_old(i-1)
        bol_old(i)=bol_old(i-1)+sndz(i)
        tol_new(i)=bol_new(i-1)
        bol_new(i)=bol_new(i-1)+thickness(i)
        enddo
        
!**** Redistribute quantities

!**** Step 1: Do top layer
      ihigh=1
      do k=1,N_snow
        if(bol_old(k) .lt. bol_new(1)) ihigh=k+1
        enddo

      do k=1,ihigh
         if(k .lt. ihigh) dz=sndz(k)
         if(k .eq. ihigh) dz=bol_new(1)-tol_old(k)
         !s(1,:)=s(1,:)+h(1,:)*dz
         s(1,:)=s(1,:)+h(k,:)*dz
#if 0
      if(pid == TAR_PE .and. ktile == TAR_TILE)then
         write(*,*) 'RELAYER2 B. *********************' 
         write(*,*) 'k = ', k, ' dz = ',dz, h(k,2)
      endif 
#endif
         enddo
#if 0
      if(pid == TAR_PE .and. ktile == TAR_TILE)then
         write(*,*) 'RELAYER2 C. *********************' 
         write(*,*) 'ihigh = ', ihigh
      endif 
#endif


!**** Step 2: Do remaining layers
      do i=2,N_snow

         ilow=ihigh
#if 0
      if(pid == TAR_PE .and. ktile == TAR_TILE)then
         write(*,*) 'i = ', i, 'ihigh = ', ihigh, ' ilow = ', ilow
      endif 
#endif
         do k=ilow,N_snow
            if(bol_old(k) .lt. bol_new(i)) ihigh=k+1
#if 0
      if(pid == TAR_PE .and. ktile == TAR_TILE)then
         write(*,*) 'i = ', i, ' k = ', k, 'ihigh = ', ihigh
         write(*,*) 'i = ', i, ' k = ', k, bol_old(k), bol_new(i)
      endif 
#endif
            enddo

         !if(ihigh .gt. N_snow) then
         !   print*, i, ihigh, ilow
         !   do k=1,N_snow
         !     print*, 'k = ', k, sndz(k), thickness(k), bol_old(k), bol_new(k) 
         !   enddo  
         !   do k=1,N_snow-1
         !     print*, 'k = ', k, thickdist(k)
         !   enddo  
         !   stop 
         !   print*, 'pid = ', pid, ' ktile = ',ktile
         !endif  
#if 0
      if(pid == TAR_PE .and. ktile == TAR_TILE)then
         write(*,*) 'i = ', i, 'ihigh = ', ihigh, ' ilow = ', ilow
      endif 
#endif
         if(ihigh .eq. N_snow+1)  ihigh=N_snow ! Account for potential truncation problem 
      
         do k=ilow,ihigh
           if(k .eq. ilow .and. k .lt. ihigh) dz=bol_old(k)-tol_new(i)
           if(k .eq. ilow .and. k .eq. ihigh) dz=bol_new(i)-tol_new(i)
           if(k .gt. ilow .and. k .lt. ihigh) dz=bol_old(k)-tol_old(k)
           if(k .gt. ilow .and. k .eq. ihigh) dz=bol_new(i)-tol_old(k)
           s(i,:)=s(i,:)+h(k,:)*dz
#if 0
      if(pid == TAR_PE .and. ktile == TAR_TILE)then
         write(*,*) 'i = ', i, ' k = ', k, ' dz = ',dz, h(k,2)
      endif 
#endif
           enddo

        enddo


      htsnn = s(:,1)
      wesn  = s(:,2)
      do k=1,N_Constit
        rconstit(:,k)=s(:,2+k)
        enddo
      sndz=thickness
      
      if(sum(wesn) < wemin) sndz = sndz /(areasc + small)
      return
      
      end subroutine relayer2


! **********************************************************************

      subroutine get_tf0d(h,w,t,f,ice1,tzero)

      implicit none
      
      real, parameter :: cpw    = 2065.22   !  @ 0 C [J/kg/K]
!      real, parameter :: lhv    = 2.4548E6 !  2.5008e6   !  @ 0 C [J/kg]
!      real, parameter :: lhs    = 2.8368E6 !  2.8434e6 !  @ 0 C [J/kg]
!rr   real, parameter :: lhv    = 2.5008e6  !  @ 0 C [J/kg]
!rr   real, parameter :: lhs    = 2.8434e6  !  @ 0 C [J/kg]
!      real, parameter :: lhf    = (lhs-lhv) !  @ 0 C [J/kg]
      real, parameter :: tfac=1./cpw
      real, parameter :: ffac=1./alhm
      
      real, intent(in )   :: w, h
      real, intent(out)   :: t, f

      logical, intent(out) :: ice1,tzero
      
      real :: hbw
      
      hbw=0.
      if(w > 0.) hbw = h/w
      
      if(hbw < -1.00001*alhm) then
            t = (hbw+alhm)*tfac
            f = 1.
            ice1=.true.
            tzero=.false.
         elseif(hbw > -0.99999*alhm) then
            t = 0.
            f =-hbw*ffac
            ice1=.false.
            tzero=.true.
         else
            t = 0.
            f = 1.
            ice1=.true.
            tzero=.true.
         endif
      
      if(f < 0.) then
         t = hbw*tfac
         f = 0.
         endif
      
      if(w == 0.) then
         t = 0.
         f = 0.
         endif

      return

      end subroutine get_tf0d
      
! **********************************************************************
      
      subroutine get_tf_nd(N,h,w,t,f)
      
!     n-dimensional version of get_tf
!     
!     avoid slow "where" statements
!     
!     can be called for any number of layers or catchments, for example

!     1.) call get_tf_nd( ncatm, htsnn1(1:ncatm), wesn1(1:ncatm),
!                         tpsn(1:ncatm),f(1:ncatm) )
!     
!     2.) call get_tf_nd(N_snow, h, w, t, f)

!     reichle, 22 Aug 2002
!     reichle, 29 Apr 2003 (updated parameter values)

      integer, intent(in) :: N
      
      real, dimension(n), intent(in)    :: h, w
      real, dimension(n), intent(out)   :: t, f
      
!     local variables
   
      real, parameter :: cpw    = 2065.22   !  @ 0 C [J/kg/K]
!      real, parameter :: lhv    = 2.4548E6 !  2.5008e6   !  @ 0 C [J/kg]
!      real, parameter :: lhs    = 2.8368E6 !  2.8434e6 !  @ 0 C [J/kg]
!rr   real, parameter :: lhv    = 2.5008e6  !  @ 0 C [J/kg]
!rr   real, parameter :: lhs    = 2.8434e6 !  @ 0 C [J/kg]
!      real, parameter :: lhf    = (lhs-lhv) !  @ 0 C [J/kg]
      real, parameter :: tfac=1./cpw
      real, parameter :: ffac=1./alhm
      
      integer :: i      
      
      real :: hbw
            
      do i=1,N
         
         if(w(i) .gt. 0.0) then
            hbw = h(i)/w(i)
         else
            hbw = 0.
         endif
         
         if(hbw .lt. -alhm) then
            t(i) = (hbw+alhm)*tfac
            f(i) = 1.
         elseif(hbw .gt. -alhm) then
            t(i) =  0.
            f(i) = -hbw*ffac
         else
            t(i) = 0.
            f(i) = 1.
         endif
         
         if(f(i) .lt. 0.) then
            t(i) = hbw*tfac
            f(i) = 0.
         endif
         
         if(w(i) .eq. 0.) then
            t(i) = 0.
            f(i) = 0.
         endif
         
      end do
      
      return
      
      end subroutine get_tf_nd
      
! **********************************************************************

      SUBROUTINE TRID(X,DD,D,RD,B,N)
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: N
      REAL*4, INTENT(IN), DIMENSION(N) :: DD, RD
      REAL*4, INTENT(INOUT), DIMENSION(N) :: D, B
      REAL*4, INTENT(OUT),DIMENSION(N) :: X

      integer I,J
      real*4  RSF
      RSF=0.
      DO 10 I=2,N
         J=N+1-I
         if(D(J+1).ne.0.) RSF=RD(J)/D(J+1)
         D(J)=D(J)-DD(J+1)*RSF
   10 B(J)=B(J)- B(J+1)*RSF
      if(D(1).ne.0.) X(1)=B(1)/D(1)
      DO 20 J=2,N
   20 if(D(J).ne.0.) X(J)=(B(J)-DD(J)*X(J-1))/D(J)
      RETURN
      END SUBROUTINE TRID

      SUBROUTINE SIBALB (                                                      &
                         NCH, ITYP, VLAI, VGRN, ZTH,                           &
                         SCALVDR,SCALVDF,SCALIDR,SCALIDF,                      &
                         WESN,SNDZ,                                            &
      			 AVISDR, ANIRDR, AVISDF, ANIRDF,                       &
                         ASNVDR, ASNNDR, ASNVDF, ASNNDF                        &
              		)

      IMPLICIT NONE

! OUTPUTS:
! AVISDR:   visible, direct albedo.
! ANIRDR:   near infra-red, direct albedo.
! AVISDF:   visible, diffuse albedo.
! ANIRDF:   near infra-red, diffuse albedo.

! INPUTS:
! SCALVDR:  MODIS scale factor for visible, direct.
! SCALVDF:  MODIS scale factor for visible, diffuse.
! SCALIDR:  MODIS scale factor for NIR, direct.
! SCALIDF:  MODIS scale factor for NIR, diffuse.
! VLAI:     the leaf area index.
!VGRN:     the greenness index.
! ZTH:      The cosine of the solar zenith angle.
! SNW:      Snow cover in meters water equivalent.


      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) ::  ITYP
      REAL, INTENT(IN), DIMENSION(NCH) :: SCALVDR,SCALVDF,SCALIDR,SCALIDF,     &
                                             VLAI, VGRN,  ZTH
      REAL, INTENT(IN), DIMENSION(:,:) :: WESN, SNDZ
      REAL, INTENT(OUT), DIMENSION(NCH) :: AVISDR, ANIRDR, AVISDF,             &
                        ANIRDF, ASNVDR, ASNNDR, ASNVDF, ASNNDF


      REAL, PARAMETER :: ALVDRS = 0.100
      REAL, PARAMETER :: ALIDRS = 0.200
      REAL, PARAMETER :: ALVDRD = 0.300
      REAL, PARAMETER :: ALIDRD = 0.350
      REAL, PARAMETER :: ALVDRI = 0.700
      REAL, PARAMETER :: ALIDRI = 0.700


!      REAL, PARAMETER :: WEMIN  = 13.0   ! [KG/M2]

! ALVDRS:  Albedo of soil for visible   direct  solar radiation.
! ALIDRS:  Albedo of soil for infra-red direct  solar radiation.
! ALVDFS:  Albedo of soil for visible   diffuse solar radiation.
! ALIDFS:  Albedo of soil for infra-red diffuse solar radiation.

      INTEGER, PARAMETER :: NLAI = 14

      REAL, PARAMETER :: EPSLN = 1.E-6
      REAL, PARAMETER :: BLAI = 0.5
      REAL, PARAMETER :: DLAI = 0.5

      REAL, PARAMETER :: ALATRM = BLAI + (NLAI - 1) * DLAI - EPSLN

      INTEGER, PARAMETER :: NTYPS_SIB=9

      REAL :: SWE, TOTDEP, AREASC, DENSITY, DENS_EXC, FRACV, SNWMASK,          &
                  AMASK, ASNVDR_VEG, ASNNDR_VEG, ASNVDF_VEG, ASNNDF_VEG
      REAL GK_B





! ITYP: Vegetation type as follows:
!                  1:  BROADLEAF EVERGREEN TREES
!                  2:  BROADLEAF DECIDUOUS TREES
!                  3:  NEEDLELEAF TREES
!                  4:  GROUND COVER
!                  5:  BROADLEAF SHRUBS
!                  6:  DWARF TREES (TUNDRA)
!                  7:  BARE SOIL
!                  8:  DESERT
!                  9:  ICE
!  NCH: Chip index
!

	INTEGER I, LAI
	REAL FAC, GAMMA, BETA, ALPHA, DX, DY, ALA, FVEG
        REAL, DIMENSION(2) :: GRN
        REAL, DIMENSION(4,NTYPS_SIB) :: SNWALB (4, NTYPS_SIB)
        REAL, DIMENSION(NTYPS_SIB) :: SNWMSK


      DATA GRN /0.33, 0.67/
      !REAL, PARAMETER :: SNWALB_VISMAX = 1.0 
      REAL, PARAMETER :: SNWALB_VISMAX = 0.921  ! matches GK94
      REAL, PARAMETER :: SNWALB_VISMIN = 0.5
      !REAL, PARAMETER :: SNWALB_NIRMAX = 0.8
      REAL, PARAMETER :: SNWALB_NIRMAX = 0.725  ! matches GK94
      REAL, PARAMETER :: SNWALB_NIRMIN = 0.3
!      REAL, PARAMETER :: RHOFS = 150.                ! DENSITY OF FRESH SNOW
      REAL, DIMENSION(NTYPS_SIB) :: SNWMID

!       DATA SNWALB/.85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50                                      &
!      		  /

!      DATA SNWMSK/25., 5., 10., 0.2, 0.5, 0.2, 0.1, 0.1, 0.1/

!***  grassland and tundra values arbitrarily increased.
      DATA SNWMID /50.,30.,45.,20.,30.,20.,2.,2.,2./




! [ Definition of Functions: ]
!
!	REAL COEFFSIB

! --------------------------------------------------



!   Constants used in albedo calculations:

      REAL ALVDR (NLAI, 2, NTYPS_SIB)
      REAL BTVDR (NLAI, 2, NTYPS_SIB)
      REAL GMVDR (NLAI, 2, NTYPS_SIB)
      REAL ALIDR (NLAI, 2, NTYPS_SIB)
      REAL BTIDR (NLAI, 2, NTYPS_SIB)
      REAL GMIDR (NLAI, 2, NTYPS_SIB)
      
!  (Data statements for ALVDR described in full; data statements for
!   other constants follow same framework.)


!    BROADLEAF EVERGREEN (ITYP=4); GREEN=0.33; LAI: .5-7
	DATA (ALVDR (I, 1, 1), I = 1, 14)                                      &
      	  /0.0808, 0.0796, 0.0792, 0.0790, 10*0.0789/

!    BROADLEAF EVERGREEN (ITYP=4); GREEN=0.67; LAI: .5-7
	DATA (ALVDR (I, 2, 1), I = 1, 14)                                      &
      	  /0.0788, 0.0775, 0.0771, 0.0769, 10*0.0768/

!    BROADLEAF DECIDUOUS (ITYP=1); GREEN=0.33; LAI: .5-7
	DATA (ALVDR (I, 1, 2), I = 1, 14)                                      &
      	  /0.0803, 0.0790, 0.0785, 0.0784, 3*0.0783, 7*0.0782/

!    BROADLEAF DECIDUOUS (ITYP=1); GREEN=0.67; LAI: .5-7
	DATA (ALVDR (I, 2, 2), I = 1, 14)                                      &
      	  /0.0782, 0.0770, 0.0765, 0.0763, 10*0.0762/

!    NEEDLELEAF (ITYP=3); GREEN=0.33; LAI=.5-7
	DATA (ALVDR (I, 1, 3), I = 1, 14)                                      &
      	  /0.0758, 0.0746, 0.0742, 0.0740, 10*0.0739/

!    NEEDLELEAF (ITYP=3); GREEN=0.67; LAI=.5-7
	DATA (ALVDR (I, 2, 3), I = 1, 14)                                      &
      	  /0.0683, 0.0672, 0.0667, 2*0.0665, 9*0.0664/

!    GROUNDCOVER (ITYP=2); GREEN=0.33; LAI=.5-7    
	DATA (ALVDR (I, 1, 4), I = 1, 14)                                      &
      	  /0.2436, 0.2470, 0.2486, 0.2494, 0.2498, 0.2500, 2*0.2501,           &
      		6*0.2502 /

!    GROUNDCOVER (ITYP=2); GREEN=0.67; LAI=.5-7
	DATA (ALVDR (I, 2, 4), I = 1, 14) /14*0.1637/

!    BROADLEAF SHRUBS (ITYP=5); GREEN=0.33,LAI=.5-7
        DATA (ALVDR (I, 1, 5), I = 1, 14)                                      &
          /0.0807, 0.0798, 0.0794, 0.0792, 0.0792, 9*0.0791/

!    BROADLEAF SHRUBS (ITYP=5); GREEN=0.67,LAI=.5-7
        DATA (ALVDR (I, 2, 5), I = 1, 14)                                      &
          /0.0787, 0.0777, 0.0772, 0.0771, 10*0.0770/

!    DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.33,LAI=.5-7
        DATA (ALVDR (I, 1, 6), I = 1, 14)                                      &
          /0.0802, 0.0791, 0.0787, 0.0786, 10*0.0785/

!    DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.67,LAI=.5-7
        DATA (ALVDR (I, 2, 6), I = 1, 14)                                      &
          /0.0781, 0.0771, 0.0767, 0.0765, 0.0765, 9*0.0764/


!    BARE SOIL
	DATA (ALVDR (I, 1, 7), I = 1, 14) /14*ALVDRS/
	DATA (ALVDR (I, 2, 7), I = 1, 14) /14*ALVDRS/

!    DESERT
	DATA (ALVDR (I, 1, 8), I = 1, 14) /14*ALVDRD/
	DATA (ALVDR (I, 2, 8), I = 1, 14) /14*ALVDRD/

!    ICE
	DATA (ALVDR (I, 1, 9), I = 1, 14) /14*ALVDRI/
	DATA (ALVDR (I, 2, 9), I = 1, 14) /14*ALVDRI/
!****
!**** -------------------------------------------------
	DATA (BTVDR (I, 1, 1), I = 1, 14)                                      &
        /0.0153, 0.0372, 0.0506, 0.0587, 0.0630, 0.0652, 0.0663,               &
      	0.0668, 0.0671, 0.0672, 4*0.0673 /
	DATA (BTVDR (I, 2, 1), I = 1, 14)                                      &
     	  /0.0135, 0.0354, 0.0487, 0.0568, 0.0611, 0.0633, 0.0644,             &
     	0.0650, 0.0652, 0.0654, 0.0654, 3*0.0655 /
	DATA (BTVDR (I, 1, 2), I = 1, 14)                                      &
      	  /0.0148, 0.0357, 0.0462, 0.0524, 0.0554, 0.0569, 0.0576,             &
      	0.0579, 0.0580, 0.0581, 0.0581, 3*0.0582 /
	DATA (BTVDR (I, 2, 2), I = 1, 14)                                      &
      	  /0.0131, 0.0342, 0.0446, 0.0508, 0.0539, 0.0554, 0.0560,             &
      	0.0564, 0.0565, 5*0.0566 /
	DATA (BTVDR (I, 1, 3), I = 1, 14)                                      &
      	  /0.0108, 0.0334, 0.0478, 0.0571, 0.0624, 0.0652, 0.0666,             &
      	0.0673, 0.0677, 0.0679, 4*0.0680 /
	DATA (BTVDR (I, 2, 3), I = 1, 14)                                      &
      	  /0.0034, 0.0272, 0.0408, 0.0501, 0.0554, 0.0582, 0.0597,             &
      		0.0604, 0.0608, 0.0610, 4*0.0611 /
	DATA (BTVDR (I, 1, 4), I = 1, 14)                                      &
      	  /0.2050, 0.2524, 0.2799, 0.2947, 0.3022, 0.3059, 0.3076,             &
      		0.3085, 0.3088, 0.3090, 4*0.3091 /
	DATA (BTVDR (I, 2, 4), I = 1, 14)                                      &
      	  /0.1084, 0.1404, 0.1617, 0.1754, 0.1837, 0.1887, 0.1915,             &
      		0.1931, 0.1940, 0.1946, 0.1948, 0.1950, 2*0.1951  /
        DATA (BTVDR (I, 1, 5), I = 1, 14)                                      &
          /0.0203, 0.0406, 0.0548, 0.0632, 0.0679, 0.0703, 0.0716,             &
           0.0722, 0.0726, 0.0727, 0.0728, 0.0728, 0.0728, 0.0729 /
        DATA (BTVDR (I, 2, 5), I = 1, 14)                                      &
          /0.0184, 0.0385, 0.0526, 0.0611,  0.0658, 0.0683, 0.0696,            &
           0.0702, 0.0705, 0.0707, 4*0.0708 /
        DATA (BTVDR (I, 1, 6), I = 1, 14)                                      &
          /0.0199, 0.0388, 0.0494,  0.0554, 0.0584, 0.0599, 0.0606,            &
           0.0609, 0.0611, 5*0.0612  /
        DATA (BTVDR (I, 2, 6), I = 1, 14)                                      &
          /0.0181, 0.0371, 0.0476, 0.0537,  0.0568, 0.0583, 0.0590,            &
           0.0593, 0.0595, 0.0595, 4*0.0596 /
	DATA (BTVDR (I, 1, 7), I = 1, 14) /14*0./
	DATA (BTVDR (I, 2, 7), I = 1, 14) /14*0./
	DATA (BTVDR (I, 1, 8), I = 1, 14) /14*0./
	DATA (BTVDR (I, 2, 8), I = 1, 14) /14*0./
	DATA (BTVDR (I, 1, 9), I = 1, 14) /14*0./
	DATA (BTVDR (I, 2, 9), I = 1, 14) /14*0./

!****
!**** -----------------------------------------------------------
	DATA (GMVDR (I, 1, 1), I = 1, 14)                                      &
       	  /0.0814, 0.1361, 0.2078, 0.2650, 0.2986, 0.3169,  0.3265,            &
        	   0.3313, 0.3337, 0.3348, 0.3354, 0.3357, 2*0.3358 /
	DATA (GMVDR (I, 2, 1), I = 1, 14)                                      &
       	  /0.0760, 0.1336, 0.2034, 0.2622, 0.2969, 0.3159,  0.3259,            &
       	   0.3309, 0.3333, 0.3346, 0.3352, 0.3354, 2*0.3356 /
	DATA (GMVDR (I, 1, 2), I = 1, 14)                                      &
       	  /0.0834, 0.1252, 0.1558, 0.1927, 0.2131,   0.2237, 0.2290,           &
       	   0.2315, 0.2327, 0.2332, 0.2335, 2*0.2336, 0.2337 /
	DATA (GMVDR (I, 2, 2), I = 1, 14)                                      &
       	  /0.0789, 0.1235, 0.1531, 0.1912, 0.2122, 0.2232,  0.2286,            &
      	   0.2312, 0.2324, 0.2330, 0.2333, 0.2334, 2*0.2335 /
	DATA (GMVDR (I, 1, 3), I = 1, 14)                                      &
       	  /0.0647, 0.1342, 0.2215, 0.2968, 0.3432, 0.3696, 0.3838,             &
       	   0.3912, 0.3950, 0.3968, 0.3978, 0.3982, 0.3984, 0.3985 /
	DATA (GMVDR (I, 2, 3), I = 1, 14)                                      &
       	  /0.0258, 0.1227, 0.1999, 0.2825, 0.3339, 0.3634, 0.3794,             &
       	   0.3877, 0.3919, 0.3940, 0.3950, 0.3956, 0.3958, 0.3959 /
	DATA (GMVDR (I, 1, 4), I = 1, 14)                                      &
       	  /0.3371, 0.5762, 0.7159, 0.7927, 0.8324, 0.8526,  0.8624,            &
       	   0.8671, 0.8693, 0.8704, 0.8709, 0.8710, 2*0.8712 /
	DATA (GMVDR (I, 2, 4), I = 1, 14)                                      &
       	  /0.2634, 0.4375, 0.5532, 0.6291, 0.6763, 0.7048, 0.7213,             &
       	   0.7310, 0.7363, 0.7395, 0.7411, 0.7420, 0.7426, 0.7428 /
        DATA (GMVDR (I, 1, 5), I = 1, 14)                                      &
           /0.0971, 0.1544, 0.2511, 0.3157, 0.3548, 0.3768, 0.3886,            &
            0.3948, 0.3978, 0.3994, 0.4001, 0.4006, 0.4007, 0.4008 /
        DATA (GMVDR (I, 2, 5), I = 1, 14)                                      &
           /0.0924, 0.1470, 0.2458, 0.3123, 0.3527, 0.3756, 0.3877,            &
            0.3942, 0.3974, 0.3990, 0.3998, 0.4002, 0.4004, 0.4005 /
        DATA (GMVDR (I, 1, 6), I = 1, 14)                                      &
           /0.0970, 0.1355, 0.1841, 0.2230, 0.2447,  0.2561, 0.2617,           &
            0.2645, 0.2658, 0.2664, 0.2667, 3*0.2669 /
        DATA (GMVDR (I, 2, 6), I = 1, 14)                                      &
           /0.0934, 0.1337, 0.1812, 0.2213, 0.2437, 0.2554, 0.2613,            &
            0.2642, 0.2656, 0.2662, 0.2665, 0.2667, 0.2667, 0.2668 /
	DATA (GMVDR (I, 1, 7), I = 1, 14) /14*1./
	DATA (GMVDR (I, 2, 7), I = 1, 14) /14*1./
	DATA (GMVDR (I, 1, 8), I = 1, 14) /14*1./
	DATA (GMVDR (I, 2, 8), I = 1, 14) /14*1./
	DATA (GMVDR (I, 1, 9), I = 1, 14) /14*1./
	DATA (GMVDR (I, 2, 9), I = 1, 14) /14*1./

!****
!****  -----------------------------------------------------------

	DATA (ALIDR (I, 1, 1), I = 1, 14)                                      &
       	  /0.2867,  0.2840, 0.2828, 0.2822, 0.2819, 0.2818, 2*0.2817,          &
       	   6*0.2816 /
	DATA (ALIDR (I, 2, 1), I = 1, 14)                                      &
        	  /0.3564, 0.3573, 0.3577, 0.3580, 2*0.3581, 8*0.3582 /
	DATA (ALIDR (I, 1, 2), I = 1, 14)                                      &
       	  /0.2848, 0.2819, 0.2804, 0.2798, 0.2795, 2*0.2793, 7*0.2792 /
	DATA (ALIDR (I, 2, 2), I = 1, 14)                                      &
       	  /0.3544, 0.3550, 0.3553, 2*0.3555, 9*0.3556 /
	DATA (ALIDR (I, 1, 3), I = 1, 14)                                      &
       	  /0.2350, 0.2311, 0.2293, 0.2285, 0.2281, 0.2280, 8*0.2279 /
	DATA (ALIDR (I, 2, 3), I = 1, 14)                                      &
       	  /0.2474, 0.2436, 0.2418, 0.2410, 0.2406, 0.2405, 3*0.2404,           &
       	   5*0.2403 /
	DATA (ALIDR (I, 1, 4), I = 1, 14)                                      &
       	  /0.5816, 0.6157, 0.6391, 0.6556, 0.6673, 0.6758, 0.6820,             &
       	   0.6866, 0.6899, 0.6924, 0.6943, 0.6956, 0.6966, 0.6974 /
	DATA (ALIDR (I, 2, 4), I = 1, 14)                                      &
       	  /0.5489, 0.5770, 0.5955, 0.6079, 0.6163, 0.6221, 0.6261,             &
       	   0.6288, 0.6308, 0.6321, 0.6330, 0.6337, 0.6341, 0.6344 /
        DATA (ALIDR (I, 1, 5), I = 1, 14)                                      &
           /0.2845, 0.2837, 0.2832, 0.2831, 0.2830, 9*0.2829 /
        DATA (ALIDR (I, 2, 5), I = 1, 14)                                      &
           /0.3532, 0.3562, 0.3578,  0.3586, 0.3590, 0.3592, 0.3594,           &
            0.3594, 0.3594, 5*0.3595 /
        DATA (ALIDR (I, 1, 6), I = 1, 14)                                      &
           /0.2825, 0.2812, 0.2806, 0.2803, 0.2802, 9*0.2801 /
        DATA (ALIDR (I, 2, 6), I = 1, 14)                                      &
           /0.3512, 0.3538,  0.3552, 0.3559, 0.3562, 0.3564, 0.3565,           &
            0.3565, 6*0.3566 /
	DATA (ALIDR (I, 1, 7), I = 1, 14) /14*ALIDRS/
	DATA (ALIDR (I, 2, 7), I = 1, 14) /14*ALIDRS/
	DATA (ALIDR (I, 1, 8), I = 1, 14) /14*ALIDRD/
	DATA (ALIDR (I, 2, 8), I = 1, 14) /14*ALIDRD/
	DATA (ALIDR (I, 1, 9), I = 1, 14) /14*ALIDRI/
	DATA (ALIDR (I, 2, 9), I = 1, 14) /14*ALIDRI/

!****
!**** -----------------------------------------------------------
	DATA (BTIDR (I, 1, 1), I = 1, 14)                                      &
       	  /0.1291, 0.1707, 0.1969, 0.2125, 0.2216,   0.2267, 0.2295,           &
       	   0.2311, 0.2319, 0.2323, 0.2326, 2*0.2327, 0.2328 /
	DATA (BTIDR (I, 2, 1), I = 1, 14)                                      &
       	  /0.1939, 0.2357, 0.2598, 0.2735, 0.2810,  0.2851, 0.2874,            &
       	   0.2885, 0.2892, 0.2895, 0.2897, 3*0.2898 /
	DATA (BTIDR (I, 1, 2), I = 1, 14)                                      &
       	  /0.1217, 0.1522, 0.1713, 0.1820,   0.1879,  0.1910, 0.1926,          &
      	   0.1935, 0.1939, 0.1942, 2*0.1943, 2*0.1944 /
	DATA (BTIDR (I, 2, 2), I = 1, 14)                                      &
       	  /0.1781, 0.2067, 0.2221, 0.2301,   0.2342,  0.2363, 0.2374,          &
       	   0.2379, 0.2382, 0.2383, 2*0.2384, 2*0.2385 /
	DATA (BTIDR (I, 1, 3), I = 1, 14)                                      &
       	  /0.0846, 0.1299, 0.1614, 0.1814, 0.1935,   0.2004, 0.2043,           &
           0.2064, 0.2076, 0.2082, 0.2085, 2*0.2087, 0.2088 /
	DATA (BTIDR (I, 2, 3), I = 1, 14)                                      &
       	  /0.0950, 0.1410, 0.1722, 0.1921, 0.2042, 0.2111,  0.2151,            &
       	   0.2172, 0.2184, 0.2191, 0.2194, 0.2196, 2*0.2197 /
	DATA (BTIDR (I, 1, 4), I = 1, 14)                                      &
       	  /0.5256, 0.7444, 0.9908, 1.2700, 1.5680, 1.8505, 2.0767,             &
       	   2.2211, 2.2808, 2.2774, 2.2362, 2.1779, 2.1160, 2.0564 /
	DATA (BTIDR (I, 2, 4), I = 1, 14)                                      &
       	  /0.4843, 0.6714, 0.8577, 1.0335, 1.1812, 1.2858, 1.3458,             &
       	   1.3688, 1.3685, 1.3546, 1.3360, 1.3168, 1.2989, 1.2838 /
	DATA (BTIDR (I, 1, 5), I = 1, 14)                                      &
           /0.1498, 0.1930, 0.2201, 0.2364, 0.2460, 0.2514, 0.2544,            &
            0.2560, 0.2569, 0.2574, 0.2577, 0.2578, 0.2579, 0.2579 /
        DATA (BTIDR (I, 2, 5), I = 1, 14)                                      &
           /0.2184, 0.2656, 0.2927, 0.3078, 0.3159,  0.3202, 0.3224,           &
            0.3235, 0.3241, 0.3244, 0.3245, 3*0.3246 /
        DATA (BTIDR (I, 1, 6), I = 1, 14)                                      &
           /0.1369, 0.1681, 0.1860, 0.1958, 0.2010,  0.2038, 0.2053,           &
            0.2060, 0.2064, 0.2066, 0.2067, 3*0.2068 /
        DATA (BTIDR (I, 2, 6), I = 1, 14)                                      &
           /0.1969, 0.2268, 0.2416,  0.2488, 0.2521, 0.2537, 0.2544,           &
            0.2547, 0.2548, 5*0.2549 / 
	DATA (BTIDR (I, 1, 7), I = 1, 14) /14*0./
	DATA (BTIDR (I, 2, 7), I = 1, 14) /14*0./
	DATA (BTIDR (I, 1, 8), I = 1, 14) /14*0./
	DATA (BTIDR (I, 2, 8), I = 1, 14) /14*0./
	DATA (BTIDR (I, 1, 9), I = 1, 14) /14*0./
	DATA (BTIDR (I, 2, 9), I = 1, 14) /14*0./

!****
!**** --------------------------------------------------------------
	DATA (GMIDR (I, 1, 1), I = 1, 14)                                      &
       	  /0.1582, 0.2581, 0.3227, 0.3635, 0.3882, 0.4026, 0.4108,             &
       	   0.4154, 0.4179, 0.4193, 0.4200, 0.4204, 0.4206, 0.4207 /
	DATA (GMIDR (I, 2, 1), I = 1, 14)                                      &
       	  /0.1934, 0.3141, 0.3818, 0.4200, 0.4415, 0.4533, 0.4598,             &
       	   0.4633, 0.4651, 0.4662, 0.4667, 0.4671, 2*0.4672 /
	DATA (GMIDR (I, 1, 2), I = 1, 14)                                      &
       	  /0.1347, 0.1871, 0.2277, 0.2515, 0.2651, 0.2727, 0.2768,             &
       	   0.2790, 0.2801, 0.2808, 0.2811, 0.2812, 0.2813, 0.2814 /
	DATA (GMIDR (I, 2, 2), I = 1, 14)                                      &
       	  /0.1440, 0.2217, 0.2629, 0.2839, 0.2947, 0.3003, 0.3031,             &
       	   0.3046, 0.3054, 0.3058, 0.3060, 2*0.3061, 0.3062 /
	DATA (GMIDR (I, 1, 3), I = 1, 14)                                      &
       	  /0.1372, 0.2368, 0.3235, 0.3839, 0.4229, 0.4465, 0.4602,             &
       	   0.4679, 0.4722, 0.4745, 0.4758, 0.4764, 0.4768, 0.4770 /
	DATA (GMIDR (I, 2, 3), I = 1, 14)                                      &
       	  /0.1435, 0.2524, 0.3370, 0.3955, 0.4332, 0.4563, 0.4697,             &
       	   0.4773, 0.4815, 0.4839, 0.4851, 0.4858, 0.4861, 0.4863 /
	DATA (GMIDR (I, 1, 4), I = 1, 14)                                      &
       	  /0.4298, 0.9651, 1.6189, 2.4084, 3.2992, 4.1928, 4.9611,             &
       	   5.5095, 5.8085, 5.9069, 5.8726, 5.7674, 5.6346, 5.4944 /
	DATA (GMIDR (I, 2, 4), I = 1, 14)                                      &
       	  /0.4167, 0.8974, 1.4160, 1.9414, 2.4147, 2.7803, 3.0202,             &
      	   3.1468, 3.1954, 3.1932, 3.1676, 3.1328, 3.0958, 3.0625 /
        DATA (GMIDR (I, 1, 5), I = 1, 14)                                      &
           /0.1959, 0.3203, 0.3985, 0.4472, 0.4766, 0.4937, 0.5034,            &
            0.5088, 0.5117, 0.5134, 0.5143, 0.5147, 0.5150, 0.5152 /
        DATA (GMIDR (I, 2, 5), I = 1, 14)                                      &
           /0.2328, 0.3859, 0.4734, 0.5227, 0.5498, 0.5644, 0.5720,            &
            0.5761, 0.5781, 0.5792, 0.5797, 0.5800, 0.5802, 0.5802 /
        DATA (GMIDR (I, 1, 6), I = 1, 14)                                      &
           /0.1447, 0.2244, 0.2698, 0.2953, 0.3094, 0.3170, 0.3211,            &
            0.3233, 0.3244, 0.3250, 0.3253, 0.3255, 0.3256, 0.3256 /
        DATA (GMIDR (I, 2, 6), I = 1, 14)                                      &
           /0.1643, 0.2624, 0.3110, 0.3347, 0.3461, 0.3517, 0.3543,            &
            0.3556, 0.3562, 0.3564, 0.3565, 0.3566, 0.3566, 0.3566 /
	DATA (GMIDR (I, 1, 7), I = 1, 14) /14*1./
	DATA (GMIDR (I, 2, 7), I = 1, 14) /14*1./
	DATA (GMIDR (I, 1, 8), I = 1, 14) /14*1./
	DATA (GMIDR (I, 2, 8), I = 1, 14) /14*1./
	DATA (GMIDR (I, 1, 9), I = 1, 14) /14*1./
	DATA (GMIDR (I, 2, 9), I = 1, 14) /14*1./

!**** -----------------------------------------------------------


!FPP$ EXPAND (COEFFSIB)

      DO I=1,NCH
        ALA = AMIN1 (AMAX1 (ZERO, VLAI(I)), ALATRM)
        LAI = 1 + MAX(0, INT((ALA-BLAI)/DLAI) )
        DX = (ALA - (BLAI+(LAI-1)*DLAI)) * (ONE/DLAI)
        DY = (VGRN(I)- GRN(1)) * (ONE/(GRN(2) - GRN(1)))

        ALPHA = COEFFSIB (ALVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
        BETA  = COEFFSIB (BTVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
        GAMMA = COEFFSIB (GMVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)

        GAMMA = MAX(GAMMA,0.01)

        AVISDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
        AVISDF(I) = ALPHA-BETA                                                 &
                 + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))

        ALPHA = COEFFSIB (ALIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
        BETA  = COEFFSIB (BTIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
        GAMMA = COEFFSIB (GMIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)

        GAMMA = MAX(GAMMA,0.01)

        ANIRDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
        ANIRDF(I) = ALPHA-BETA                                                 &
                 + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))

! SCALE TO MODIS VALUES (SNOW-FREE)

	  AVISDR(I) = AVISDR(I) * SCALVDR(I)
          ANIRDR(I) = ANIRDR(I) * SCALIDR(I)
	  AVISDF(I) = AVISDF(I) * SCALVDF(I)
          ANIRDF(I) = ANIRDF(I) * SCALIDF(I)

! PROTECT AGAINST BAD SCALING

	  AVISDR(I) = AMIN1( 1., AMAX1( 0., AVISDR(I) ) )
          ANIRDR(I) = AMIN1( 1., AMAX1( 0., ANIRDR(I) ) )
	  AVISDF(I) = AMIN1( 1., AMAX1( 0., AVISDF(I) ) )
          ANIRDF(I) = AMIN1( 1., AMAX1( 0., ANIRDF(I) ) )

! SNOW ALBEDOES

        SWE=sum(WESN(I,:))
        !TOTDEP=sum(SNDZ(I,:))
        TOTDEP=SNDZ(I,1)
        AREASC = MIN(SWE/WEMIN,1.)
        !DENSITY=(SWE/(AREASC+1.e-20)) / (TOTDEP+1.e-20)
        !*** only use top layer density to dentermine albedo 
        DENSITY=(WESN(I,1)/(AREASC+1.e-20)) / (TOTDEP+1.e-20)
        DENS_EXC=MAX(0., DENSITY-RHOFRESH)

        !ASNVDR(I) = MAX(SNWALB_VISMIN, SNWALB_VISMAX - DENS_EXC*.0006)
        !ASNNDR(I) = MAX(SNWALB_NIRMIN, SNWALB_NIRMAX - DENS_EXC*.0006)
        !ASNVDF(I) = MAX(SNWALB_VISMIN, SNWALB_VISMAX - DENS_EXC*.0006)
        !ASNNDF(I) = MAX(SNWALB_NIRMIN, SNWALB_NIRMAX - DENS_EXC*.0006)

        !*** Greuell & Konzelmann 94
        !*** features a higher albedo at high densities
        GK_B = (0.85-0.58)/(RHOFRESH-RHOMA)
        ASNVDR(I) = MAX(SNWALB_VISMIN, SNWALB_VISMAX + GK_B*DENS_EXC)
        ASNNDR(I) = MAX(SNWALB_NIRMIN, SNWALB_NIRMAX + GK_B*DENS_EXC)
        ASNVDF(I) = MAX(SNWALB_VISMIN, SNWALB_VISMAX + GK_B*DENS_EXC)
        ASNNDF(I) = MAX(SNWALB_NIRMIN, SNWALB_NIRMAX + GK_B*DENS_EXC)


! ACCOUNT FOR VEGETATION MASKING, FOR EACH COMPONENT

! A) FIRST DO MASKING IN VEGETATED FRACTION:
        FAC = SWE / (SWE + SNWMID(ITYP(I)))
        ASNVDR_VEG=AVISDR(I) + (ASNVDR(I)-AVISDR(I))*FAC
        ASNNDR_VEG=ANIRDR(I) + (ASNNDR(I)-ANIRDR(I))*FAC
        ASNVDF_VEG=AVISDF(I) + (ASNVDF(I)-AVISDF(I))*FAC
        ASNNDF_VEG=ANIRDF(I) + (ASNNDF(I)-ANIRDF(I))*FAC

! B) NOW ACCOUNT FOR SUBGRID VEGETATION FRACTION
        FVEG=AMIN1( 1., VLAI(I)/2. )
        ASNVDR(I)=ASNVDR(I)*(1.-FVEG)+ASNVDR_VEG*FVEG
        ASNNDR(I)=ASNNDR(I)*(1.-FVEG)+ASNNDR_VEG*FVEG
        ASNVDF(I)=ASNVDF(I)*(1.-FVEG)+ASNVDF_VEG*FVEG
        ASNNDF(I)=ASNNDF(I)*(1.-FVEG)+ASNNDF_VEG*FVEG


!        AMASK=FVEG*EXP(-TOTDEP/SNWMSK(ITYP(I)))

!
!	  AVISDR(I) = AVISDR(I) + (SNWALB(1,ITYP(I)) - AVISDR(I)) * FAC
!          ANIRDR(I) = ANIRDR(I) + (SNWALB(2,ITYP(I)) - ANIRDR(I)) * FAC
!	  AVISDF(I) = AVISDF(I) + (SNWALB(3,ITYP(I)) - AVISDF(I)) * FAC
!          ANIRDF(I) = ANIRDF(I) + (SNWALB(4,ITYP(I)) - ANIRDF(I)) * FAC
!	  ENDIF

!  (ORIGINAL NOTES FROM STIEGLITZ:)

!        ALBSNW = 0.913 - .0006*DENSITY
!        ALBSNW = AMIN1(1., AMAX1(ALBSNW,0.5))

!  a) formulation
!  While 3 basic land covers are allowed (bare soil, snow, and vegetation),
!  the grid cell is divided into 2 fractions: a vegetated fraction (A_v),
!  and a non-vegetated fraction (A_u).  In the absence of snow cover
!  A_v = m A_v0
!  A_u = 1 - A_v
!  where A_v0 is defined to be the grid cell specified value for the
!  vegetated fraction, A_u is the resultant bare soil fraction, and m, the
!  snow masking fraction, is unity.  In the absence of snow A_u is simply
!  the unvegetated portion of the grid cell.  In the presence of snow, the
!  green vegetation masked by snow is
!  m = exp( - d_s / d_m)
!  where d_s is the snow depth(not water equivalent) and d_m is the
!  vegetation specific masking depth.  A_u now becomes that portion of the
!  grid cell where vegetation is not visible
!  b) masking depths
!  DATA SMK/0.1,0.2,0.2,0.5,2.0,5.0,10.0,25.0/
!              SMK array types - desert, tundra, grass, shrub, woodland,
!  deciduous, evergreen, rain forest


!        FRACV=0.5
!        AMASK=FRACV*EXP(-TOTDEP/SNWMSK)
!        ALBSNW=ALBSNW*(1.-AMASK)+ALBAVE*AMASK
        
!        IF (SNW (I) .GT. ZERO) THEN
!	  FAC = SNW(I) / (SNW(I) + SNWMID(ITYP(I)))
!
!	  AVISDR(I) = AVISDR(I) + (SNWALB(1,ITYP(I)) - AVISDR(I)) * FAC
!          ANIRDR(I) = ANIRDR(I) + (SNWALB(2,ITYP(I)) - ANIRDR(I)) * FAC
!	  AVISDF(I) = AVISDF(I) + (SNWALB(3,ITYP(I)) - AVISDF(I)) * FAC
!          ANIRDF(I) = ANIRDF(I) + (SNWALB(4,ITYP(I)) - ANIRDF(I)) * FAC
!	  ENDIF

        ENDDO

      RETURN
      END SUBROUTINE SIBALB

      FUNCTION COEFFSIB(TABLE, NTABL, LAI ,DX, DY)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NTABL, LAI

      REAL, INTENT(IN) :: DX, DY
      REAL, INTENT(IN), DIMENSION(NTABL,2) :: TABLE
      REAL COEFFSIB

      COEFFSIB = (TABLE(LAI,  1)                                                  &
             + (TABLE(LAI  ,2) - TABLE(LAI  ,1)) * DY ) * (1.0-DX)             &
             + (TABLE(LAI+1,1)                                                 &
             + (TABLE(LAI+1,2) - TABLE(LAI+1,1)) * DY ) * DX

      RETURN
      END FUNCTION COEFFSIB     
#endif

end subroutine RUN2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_LandiceGridCompMod

