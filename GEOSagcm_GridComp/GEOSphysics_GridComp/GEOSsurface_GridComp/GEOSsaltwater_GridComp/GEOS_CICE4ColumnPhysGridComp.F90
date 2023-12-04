
!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
module GEOS_CICE4ColumnPhysGridComp

!BOP
! !MODULE: GEOS_CICE4ColumnPhysGridComp -- Implements CICE4 on ice tiles.

! !DESCRIPTION:
!
!   {\tt GEOS\_CICE4ColumnPhys} is a light-weight gridded component that updates
!      the skin sub-tiles at saltwater points, be they ocean, estuary, or salt
!      lake. Currently each tile can have multiple subtiles representing ice categories,
!      and includes implementation of LANL CICE thermodynamics.
!
!      The component is written with a two stage run method for use with
!      semi-implicit turbulence components. The first run stage computes
!      exchange coefficients for heat, moisture and momentum at each sub-tile
!      and combines these to tile space, accounting for sub tile variability
!      by passing back an effective surface value of the exchanged quantity.
!

! !USES:

  use sfclayer  ! using module that contains sfc layer code
  use ESMF
  use MAPL
  use GEOS_UtilsMod
  use DragCoefficientsMod

! LANL CICE Thermodynamics modules
  use ice_kinds_mod
  use ice_constants,      only: TFfresh, puny, c0
  use ice_constants,      only: rad_to_deg
  use ice_constants,      only: awtvdr, awtvdf, awtidr, awtidf
  use ice_constants,      only: m2_to_km2
  use ice_constants,      only: depressT
  use ice_constants,      only: Tocnfrz
  use ice_constants,      only: Lfresh, rhos, cp_ice
  use ice_domain_size,    only: init_column_physics
  use ice_itd,            only: init_itd, ilyr1, cleanup_itd
  use ice_therm_vertical, only: init_thermo_vertical,  &
                                init_vertical_profile, &
                                thermo_vertical,       &
                                diagnose_internal_ice_temp, &
                                frzmlt_bottom_lateral
  use ice_state,          only: nt_tsfc, nt_iage, nt_volpn, init_trcr_depend
  use ice_therm_itd,      only: linear_itd, add_new_ice, lateral_melt,    &
                                freeboard_ccsm
  use ice_init,           only: input_data, set_state_var, &
                                alloc_column_physics, dealloc_column_physics
  use ice_age,            only: iage_converter
  use ice_meltpond,       only: compute_ponds

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

  integer, parameter :: ICE = 1
  integer, parameter :: NUM_3D_ICE_TRACERS = 3
  integer, parameter :: NUM_SNOW_LAYERS    = 1

  logical ::      DUAL_OCEAN

  type cice_state
       integer:: CHOOSEMOSFC
  end type cice_state

  type cice_state_wrap
      type(cice_state), pointer :: ptr
  end type

#define PACKIT   1
#define UNPACKIT 2

  contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

    !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices, which sets
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

! Local derived type aliases

    type (MAPL_MetaComp),  pointer          :: MAPL
    type (ESMF_Config)                      :: CF

    integer                                 :: NUM_SUBTILES        ! = NUM_ICE_CATEGORIES
    integer                                 :: NUM_ICE_LAYERS      ! set via resource parameter
    integer                                 :: NUM_ICE_CATEGORIES  ! set via resource parameter
    integer ::      iDUAL_OCEAN

    type(cice_state_wrap) :: wrap
    type(cice_state), pointer :: mystate
    character(len=ESMF_MAXSTR)     :: SURFRC
    type(ESMF_Config)              :: SCF

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, _RC )
    Iam = trim(COMP_NAME) // 'SetServices'

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)


    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, _RC )

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run1, _RC )

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2, _RC )

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,  Finalize, _RC )

! Get constants from CF
! ---------------------

    call ESMF_ConfigGetAttribute(CF, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" , _RC)
    call ESMF_ConfigGetAttribute(CF, NUM_ICE_LAYERS    , Label="CICE_N_ICE_LAYERS:"     , _RC)
    NUM_SUBTILES  = NUM_ICE_CATEGORIES

    call MAPL_GetResource(MAPL, iDUAL_OCEAN, 'DUAL_OCEAN:', default=0, _RC )
    DUAL_OCEAN = iDUAL_OCEAN /= 0

! Set the state variable specs.
! -----------------------------

!BOS

!  !EXPORT STATE:

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'EMIS',                              &
        LONG_NAME          = 'surface_emissivity',                &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_visible_beam',   &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVR',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_visible_diffuse',&
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVF',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_near_infrared_beam', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNR',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_near_infrared_diffuse', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNF',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )


     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'evaporation'               ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'EVAPOUT'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'sublimation'               ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'SUBLIM'                    ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )


     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHOUT'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'sea_ice_upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHICE'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'surface_emitted_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLWUP'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'sea_ice_outgoing_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLWUPICE'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'sea_ice_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWNDICE'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'surface_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWNDSRF'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'sea_ice_net_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWNDICE'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'surface_net_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWNDSRF'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'total_latent_energy_flux'  ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLATN'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'sea_ice_latent_energy_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLATICE'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'total_surface_heat_flux_over_the_whole_tile' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'FSURF'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'total_surface_heat_flux_over_the_ice_tile' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'FSURFICE'                  ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TST',                               &
        LONG_NAME          = 'surface_skin_temperature',          &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'QST',                               &
        LONG_NAME          = 'surface_specific_humidity',         &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TH',                                &
        LONG_NAME          = 'turbulence_surface_temperature',    &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'QH',                                &
        LONG_NAME          = 'turbulence_surface_specific_humidity', &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'UH',                                &
        LONG_NAME          = 'turbulence_surface_zonal_velocity', &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'VH',                                &
        LONG_NAME          = 'turbulence_surface_meridional_velocity', &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELTS',                             &
        LONG_NAME          = 'change_of_surface_skin_temperature',&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELQS',                             &
        LONG_NAME          = 'change_of_surface_specific_humidity',&
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CHT',                               &
        LONG_NAME          = 'surface_heat_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CMT',                               &
        LONG_NAME          = 'surface_momentum_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CQT',                               &
        LONG_NAME          = 'surface_moisture_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CNT',                               &
        LONG_NAME          = 'neutral_drag_coefficient',          &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'RIT',                               &
        LONG_NAME          = 'surface_bulk_richardson_number',    &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'RET',                               &
        LONG_NAME          = 'surface_reynolds_number',           &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRACI',                             &
        LONG_NAME          = 'ice_covered_fraction_of_tile',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRACINEW',                             &
        LONG_NAME          = 'ice_covered_fraction_of_tile_after_update',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'GUST',                      &
        LONG_NAME          = 'gustiness',                 &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly,           &
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'VENT',                      &
        LONG_NAME          = 'surface_ventilation_velocity',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly,           &
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'surface_roughness'         ,&
        UNITS              = 'm'                         ,&
        SHORT_NAME         = 'Z0'                        ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'surface_roughness_for_heat',&
        UNITS              = 'm'                         ,&
        SHORT_NAME         = 'Z0H'                       ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOT2M',                     &
        LONG_NAME          = 'temperature 2m wind from MO sfc', &
        UNITS              = 'K',                         &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOQ2M',                     &
        LONG_NAME          = 'humidity 2m wind from MO sfc',    &
        UNITS              = 'kg kg-1',                   &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU2M',                    &
        LONG_NAME          = 'zonal 2m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV2M',                    &
        LONG_NAME          = 'meridional 2m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOT10M',                     &
        LONG_NAME          = 'temperature 10m wind from MO sfc', &
        UNITS              = 'K',                         &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOQ10M',                     &
        LONG_NAME          = 'humidity 10m wind from MO sfc',    &
        UNITS              = 'kg kg-1',                   &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU10M',                    &
        LONG_NAME          = 'zonal 10m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV10M',                    &
        LONG_NAME          = 'meridional 10m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU50M',                    &
        LONG_NAME          = 'zonal 50m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV50M',                    &
        LONG_NAME          = 'meridional 50m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'eastward_stress_over_ice',  &
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUXI'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'northward_stress_over_ice',  &
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUYI'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PENUVR',                             &
        LONG_NAME          = 'penetrative_uvr_direct_flux_through_sea_ice',           &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PENUVF',                             &
        LONG_NAME          = 'penetrative_uvr_diffuse_flux_through_sea_ice',           &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PENPAR',                             &
        LONG_NAME          = 'penetrative_par_direct_flux_through_sea_ice',           &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        _RC  )

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PENPAF',                             &
        LONG_NAME          = 'penetrative_par_diffuse_flux_through_sea_ice',           &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        _RC  )

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'TFREEZE',                        &
        LONG_NAME          = 'freezing_temperature_for_interface_layer',&
        UNITS              = 'K',                               &
        DIMS               = MAPL_DimsTileOnly,                 &
        VLOCATION          = MAPL_VLocationNone,                &
        _RC  )

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'surface_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWDNSRF'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'surface_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWDNSRF'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

!  !INTERNAL STATE:

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'HSKINI',                            &
        LONG_NAME          = 'ice_skin_layer_mass',               &
        UNITS              = 'kg m-2',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = 0.5*MAPL_RHOWTR,                     &
                                                       _RC  )

     call MAPL_AddInternalSpec(GC,                                &
         SHORT_NAME         = 'TSKINI',                            &
         LONG_NAME          = 'ice_skin_temperature',              &
         UNITS              = 'K',                                 &
         UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
         DIMS               = MAPL_DimsTileOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         FRIENDLYTO         = 'SEAICE',                            &
         DEFAULT            = MAPL_TICE-1.8,                       &
                                           _RC  )

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'SSKINI',                            &
        LONG_NAME          = 'ice_skin_salinity',                 &
        UNITS              = 'psu',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = 30.0,                                &
                                                       _RC  )

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'QS',                                &
        LONG_NAME          = 'surface_specific_humidity',         &
        UNITS              = 'kg kg-1',                           &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.01,                                &
                                                       _RC  )

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CH',                                &
        LONG_NAME          = 'surface_heat_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 1.0e-4,                              &
                                                       _RC  )

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CM',                                &
        LONG_NAME          = 'surface_momentum_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 1.0e-4,                              &
                                                       _RC  )

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CQ',                                &
        LONG_NAME          = 'surface_moisture_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 1.0e-4,                              &
                                                       _RC  )

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'Z0',                                &
        LONG_NAME          = 'aerodynamic_roughness',             &
        UNITS              = 'm',                                 &
        DEFAULT            = 0.00005,                             &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'WW',                                &
        LONG_NAME          = 'vertical_velocity_scale_squared',   &
        UNITS              = 'm+2 s-2',                           &
        DEFAULT            = 0.0,                                 &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

!  !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'ALW',                               &
        LONG_NAME          = 'linearization_of_surface_emitted_longwave_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'BLW',                               &
        LONG_NAME          = 'linearization_of_surface_emitted_longwave_flux', &
        UNITS              = 'W m-2 K-1',                         &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'LWDNSRF',                           &
        LONG_NAME          = 'surface_absorbed_longwave_flux',    &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC                             ,&
        LONG_NAME          = 'surface_downwelling_par_beam_flux' ,&
        UNITS              = 'W m-2'                             ,&
        SHORT_NAME         = 'DRPAR'                             ,&
        DIMS               = MAPL_DimsTileOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       _RC  )

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_par_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DFPAR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  _RC  )

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_nir_beam_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DRNIR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  _RC  )

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_nir_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DFNIR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  _RC  )

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_uvr_beam_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DRUVR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  _RC  )

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_uvr_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DFUVR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  _RC  )

    call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'evaporation',                       &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'EVAP ',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'upward_sensible_heat_flux',         &
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'SH',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'eastward_surface_stress',           &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUX',                              &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'northward_surface_stress',          &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUY',                              &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_evaporation',         &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'DEVAP',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_upward_sensible_heat_flux', &
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'DSH',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'snowfall',                          &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'SNO',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

! Surface air quantities

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_air_temperature',           &
        UNITS              = 'K',                                 &
        SHORT_NAME         = 'TA',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_air_specific_humidity',     &
        UNITS              = 'kg kg-1',                           &
        SHORT_NAME         = 'QA',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_wind_speed',                &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'UU',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'levellm_uwind',                     &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'UWINDLMTILE',                       &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'levellm_vwind',                     &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'VWINDLMTILE',                       &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_layer_height',              &
        UNITS              = 'm',                                 &
        SHORT_NAME         = 'DZ',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_pressure',                  &
        UNITS              = 'Pa',                                &
        SHORT_NAME         = 'PS',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'liquid_water_convective_precipitation',&
        UNITS              = 'kg m-2 s-1'                        ,&
        SHORT_NAME         = 'PCU'                               ,&
        DIMS               = MAPL_DimsTileOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       _RC  )

     call MAPL_AddImportSpec(GC                            ,&
        LONG_NAME          = 'liquid_water_large_scale_precipitation',&
        UNITS              = 'kg m-2 s-1'                       ,&
        SHORT_NAME         = 'PLS'                              ,&
        DIMS               = MAPL_DimsTileOnly                  ,&
        VLOCATION          = MAPL_VLocationNone                 ,&
                                                      _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'THATM',                             &
        LONG_NAME          = 'effective_surface_skin_temperature',&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'QHATM',                             &
        LONG_NAME          = 'effective_surface_specific_humidity',&
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'UHATM',                             &
        LONG_NAME          = 'effective_surface_zonal_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VHATM',                             &
        LONG_NAME          = 'effective_surface_meridional_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CTATM',                             &
        LONG_NAME          = 'surface_exchange_coefficient_for_heat', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CQATM',                             &
        LONG_NAME          = 'surface_exchange_coefficient_for_moisture', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CMATM',                             &
        LONG_NAME          = 'surface_exchange_coefficient_for_momentum', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'FRACICE',                           &
        LONG_NAME          = 'ice_covered_fraction_of_tile',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'UW',                                &
        LONG_NAME          = 'zonal_velocity_of_surface_water',   &
        UNITS              = 'm s-1 ',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'UI',                                &
        LONG_NAME          = 'zonal_velocity_of_surface_ice',     &
        UNITS              = 'm s-1 ',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VW',                                &
        LONG_NAME          = 'meridional_velocity_of_surface_water',   &
        UNITS              = 'm s-1 ',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VI',                                &
        LONG_NAME          = 'meridional_velocity_of_surface_ice',     &
        UNITS              = 'm s-1 ',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )


     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'SS_FOUND',                          &
        LONG_NAME          = 'foundation_salinity_for_interface_layer',               &
        UNITS              = 'psu',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 30.0,                                &

                                                       _RC  )

      call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TS_FOUND',                          &
        LONG_NAME          = 'foundation_temperature_for_interface_layer',            &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 280.0,                               &

                                                       _RC  )

      call MAPL_AddImportSpec(GC                         ,&
          SHORT_NAME         = 'FRZMLT'                    ,&
          LONG_NAME          = 'freeze_melt_potential',     &
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          DEFAULT            = 0.0,                         &
          _RC  )

      !call MAPL_AddImportSpec(GC,                                  &
      !    SHORT_NAME         = 'TFREEZE',                        &
      !    LONG_NAME          = 'freezing_temperature_for_interface_layer',&
      !    UNITS              = 'K',                               &
      !    DIMS               = MAPL_DimsTileOnly,                 &
      !    VLOCATION          = MAPL_VLocationNone,                &
      !    DEFAULT            = MAPL_TICE-1.8,                     &
      !    _RC  )

! Additions for LANL CICE Thermodynamics
!----------------------------------

!-------------------Exports---------------------------------------------------------------

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'FRAZIL',                    &
    LONG_NAME          = 'frazil_ice_growth'         ,&
    UNITS              = 'm s-1'           ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'FRESH',                     &
    LONG_NAME          = 'fresh_water_flux_to_ocean' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'FSALT',                     &
    LONG_NAME          = 'salt_flux_to_ocean'        ,&
    UNITS              = 'kg m-2 s-1'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'FHOCN',                     &
    LONG_NAME          = 'net_heat_flux_to_ocean'    ,&
    UNITS              = 'W m-2'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'PICE',                      &
    LONG_NAME          = 'sea_ice_pressure_loading'  ,&
    UNITS              = 'Pa'                        ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'FSWTHRU',                   &
    LONG_NAME          = 'SW_flux_thru_ice_to_ocean' ,&
    UNITS              = 'W m-2'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'FSWABS',                         &
    LONG_NAME          = 'SW_flux_absorbed_by_skin_layer' ,&
    UNITS              = 'W m-2'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'CONGEL',                    &
    LONG_NAME          = 'congelation_ice_growth'    ,&
    UNITS              = 'm s-1'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'SNOICE',                    &
    LONG_NAME          = 'snow_ice_formation'        ,&
    UNITS              = 'm s-1'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'MELTT',                     &
    LONG_NAME          = 'top_ice_melt'              ,&
    UNITS              = 'm s-1'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'MELTB',                     &
    LONG_NAME          = 'basal_ice_melt'            ,&
    UNITS              = 'm s-1'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'MELTL',                     &
    LONG_NAME          = 'lateral_ice_melt'          ,&
    UNITS              = 'm s-1'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'MELTS',                     &
    LONG_NAME          = 'snow_melt'                 ,&
    UNITS              = 'm s-1'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'HICE',                        &
    LONG_NAME          = 'grid_cell_mean_ice_thickness',&
    UNITS              = 'm'                         ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'HSNO',                         &
    LONG_NAME          = 'grid_cell_mean_snow_thickness',&
    UNITS              = 'm'                         ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'TSKINICE',                  &
    LONG_NAME          = 'snow_or_ice_surface_temperature',&
    UNITS              = 'K'                         ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'IAGE',                      &
    LONG_NAME          = 'sea_ice_age'               ,&
    UNITS              = 'years'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC,                           &
    SHORT_NAME         = 'DAIDTT',                                 &
    LONG_NAME          = 'ice_area_tendency_dueto_thermodynamics', &
    UNITS              = '% day-1',                                &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

  call MAPL_AddExportSpec(GC,                           &
    SHORT_NAME         = 'DVIDTT',                                   &
    LONG_NAME          = 'ice_volume_tendency_dueto_thermodynamics', &
    UNITS              = 'cm day-1',                                 &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

  call MAPL_AddExportSpec(GC,                                        &
    SHORT_NAME         = 'FBOT',                                     &
    LONG_NAME          = 'net_downward_heat_flux_from_ice_to_ocean', &
    UNITS              = 'W m-2',                                    &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

  call MAPL_AddExportSpec(GC,                                   &
    SHORT_NAME         = 'USTARI'                   ,           &
    LONG_NAME          = 'ice_ocean_friction_velocity',         &
    UNITS              = 'm s-1'                   ,            &
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC                    ,                         &
    SHORT_NAME         = 'HICEUNT',                                       &
    LONG_NAME          = 'grid_cell_mean_ice_thickness_untouched_by_run2',&
    UNITS              = 'm'                         ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC                    ,     &
    SHORT_NAME         = 'SNOONICE',                  &
    LONG_NAME          = 'snow_fall_on_top_of_ice',   &
    UNITS              = 'kg m-2 s-1'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'SIALB'                       ,&
    LONG_NAME          = 'broad_band_sea_ice_albedo'   ,&
    UNITS              = '1'                           ,&
    DIMS               = MAPL_DimsTileOnly             ,&
    VLOCATION          = MAPL_VLocationNone            ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHTSKIN',                   &
    LONG_NAME          = 'Ground_heating_for_skin_temp',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           _RC  )

  call MAPL_AddExportSpec(GC                         ,&
    SHORT_NAME         = 'FRZMLT'                    ,&
    LONG_NAME          = 'freeze_melt_potential',     &
    UNITS              = 'W m-2'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           _RC  )

  ! CMIP5 exports; this is only one part of the list, the rest are in CICEDyna

  call MAPL_AddExportSpec(GC,                          &
    SHORT_NAME         = 'evap_CMIP5'                ,&
    LONG_NAME          = 'water_evaporation_flux',    &
    UNITS              = 'kg m-2 s-1'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'pr_CMIP5'                  ,                                        &
    LONG_NAME          = 'surface_rainfall_rate_into_the_sea_ice_portion_of_the_grid_cell',   &
    UNITS              = 'kg m-2 s-1'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'prsn_CMIP5'                ,                                        &
    LONG_NAME          = 'surface_snowfall_rate_into_the_sea_ice_portion_of_the_grid_cell',   &
    UNITS              = 'kg m-2 s-1'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'grFrazil_CMIP5'            ,&
    LONG_NAME          = 'frazil_sea_ice_growth_rate',&
    UNITS              = 'kg m-2 s-1'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'grCongel_CMIP5'                    ,&
    LONG_NAME          = 'congelation_sea_ice_growth_rate',   &
    UNITS              = 'kg m-2 s-1'                        ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'grLateral_CMIP5'              ,&
        LONG_NAME          = 'lateral_sea_ice_growth_rate'  ,&
        UNITS              = 'kg m-2 s-1'                   ,&
        DIMS               = MAPL_DimsTileOnly              ,&
        VLOCATION          = MAPL_VLocationNone             ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'snoToIce_CMIP5'            ,&
        LONG_NAME          = 'snow_ice_formation_rate',   &
        UNITS              = 'kg m-2 s-1'                ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'snomelt_CMIP5'             ,&
        LONG_NAME          = 'snow_melt_rate',            &
        UNITS              = 'kg m-2 s-1'                ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'tmelt_CMIP5'               ,                 &
        LONG_NAME          = 'rate_of_melt_at_upper_surface_of_sea_ice',   &
        UNITS              = 'kg m-2 s-1'                                 ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'bmelt_CMIP5'               ,     &
        LONG_NAME          = 'rate_of_melt_at_sea_ice_base',   &
        UNITS              = 'kg m-2 s-1'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'sfdsi_CMIP5'               ,         &
        LONG_NAME          = 'downward_sea_ice_basal_salt_flux',   &
        UNITS              = 'kg m-2 s-1'                         ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'hfsifrazil_CMIP5'             ,                          &
        LONG_NAME          = 'heat_flux_into_sea_water_due_to_frazil_ice_formation',   &
        UNITS              = 'W m-2'                        ,&
        DIMS               = MAPL_DimsTileOnly              ,&
        VLOCATION          = MAPL_VLocationNone             ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME         = 'ialb_CMIP5'                   ,&
        LONG_NAME          = 'bare_sea_ice_albedo'          ,&
        UNITS              = '1'                            ,&
        DIMS               = MAPL_DimsTileOnly              ,&
        VLOCATION          = MAPL_VLocationNone             ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'rsdssi_CMIP5'                 ,             &
        LONG_NAME          = 'surface_downwelling_shortwave_flux_in_air' ,&
        UNITS              = 'W m-2'                        ,&
        DIMS               = MAPL_DimsTileOnly              ,&
        VLOCATION          = MAPL_VLocationNone             ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'rsussi_CMIP5'                 ,           &
        LONG_NAME          = 'surface_upwelling_shortwave_flux_in_air' ,&
        UNITS              = 'W m-2'                        ,&
        DIMS               = MAPL_DimsTileOnly              ,&
        VLOCATION          = MAPL_VLocationNone             ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'fsitherm_CMIP5'               ,                           &
        LONG_NAME          = 'water_flux_into_sea_water_due_to_sea_ice_thermodynamics' ,&
        UNITS              = 'kg m-2 s-1'                   ,&
        DIMS               = MAPL_DimsTileOnly              ,&
        VLOCATION          = MAPL_VLocationNone             ,&
                                               _RC  )

  call MAPL_AddExportSpec(GC,                    &
          SHORT_NAME         = 'FCONDTOP'                  ,             &
          LONG_NAME          = 'conductive_heat_flux_at_ice_top_surface',&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  )


  call MAPL_AddExportSpec(GC,                            &
         SHORT_NAME         = 'FCONDBOT'                  ,                &
          LONG_NAME          = 'conductive_heat_flux_at_ice_bottom_surface',&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  )

  call MAPL_AddExportSpec(GC,                            &
         SHORT_NAME         = 'NEWICEERG'                  ,                &
          LONG_NAME          = 'heat_flux_associated_with_new_ice_generation',&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  )

  call MAPL_AddExportSpec(GC,                            &
         SHORT_NAME          = 'SUBLIMFLX'                  ,                &
          LONG_NAME          = 'heat_flux_associated_with_sublimation_of_snow_ice',&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  )

!  Category dimensional exports

   call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME         = 'FCONDBOTN'                 ,                            &
          LONG_NAME          = 'conductive_heat_flux_at_ice_bottom_over_ice_categories',&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  )

   call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME         = 'FCONDTOPN'                 ,                            &
          LONG_NAME          = 'conductive_heat_flux_at_ice_snow_surface_over_ice_categories',&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  )

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'FSURFN'                    ,&
        LONG_NAME          = 'net_heat_flux_at_ice_snow_surface_over_ice_categories' ,&
        UNITS              = 'W m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        _RC  )

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'SHICEN'                    ,&
        LONG_NAME          = 'sea_ice_upward_sensible_heat_flux_over_ice_categories' ,&
        UNITS              = 'W m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        _RC  )

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'HLWUPN'                     ,&
        LONG_NAME          = 'outgoing_longwave_flux_at_ice_snow_surface_over_ice_categories',&
        UNITS              = 'W m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        _RC  )

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'LWNDSRFN'                     ,&
        LONG_NAME          = 'net_downward_longwave_flux_at_ice_snow_surface_over_ice_categories',&
        UNITS              = 'W m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        _RC  )

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'FSWSFCN'                    ,&
        LONG_NAME          = 'SW_absorbed_at_ice_snow_surface_over_ice_categories' ,&
        UNITS              = 'W m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        _RC  )

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'TSURFN'                    ,&
        LONG_NAME          = 'ice_snow_surface_temperature_over_ice_categories' ,&
        UNITS              = 'K'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        _RC  )

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'ALBIN'                    ,&
        LONG_NAME          = 'ice_surface_albedo_over_ice_categories' ,&
        UNITS              = '1'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        _RC  )

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'ALBSN'                    ,&
        LONG_NAME          = 'snow_surface_albedo_over_ice_categories' ,&
        UNITS              = '1'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        _RC  )

   call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'saturation_specific_humidity_using_geos_formula',&
        UNITS              = 'kg kg-1'                   ,&
        SHORT_NAME         = 'QSAT1'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        _RC  )

   call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'saturation_specific_humidity_using_bulk_formula',&
        UNITS              = 'kg kg-1'                   ,&
        SHORT_NAME         = 'QSAT2'                     ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        _RC  )

   ! this export actually has dimensions(nlayer,nice) but is collapsed
   ! into one for ease of history
   call MAPL_AddExportSpec(GC,                    &
          LONG_NAME          = 'internal_ice_temperature_over_ice_categories',&
          UNITS              = 'degC'                ,&
          SHORT_NAME         = 'TINZ'                 ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES*NUM_ICE_LAYERS/) ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  )

! budget terms

   call MAPL_AddExportSpec(GC, &
        SHORT_NAME         = 'DELTAVOL1',                         &
        LONG_NAME          = 'total_change_in_ice_volume_each_cat',&
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
        VLOCATION          = MAPL_VLocationNone,                  &
        RC=STATUS  )
   VERIFY_(STATUS)

!-------------------Internal--------------------------------------------------------------

    call MAPL_AddInternalSpec(GC,                                  &
        SHORT_NAME         = 'FR',                                &
        LONG_NAME          = 'subtile_fractions_of_grid_cell',    &
        UNITS              = '1',                                 &
         PRECISION          = MAPL_R8,                             &    ! Bin, Yury: Please listen to Matt and Atanas! Kindly work on interfacing
         DIMS               = MAPL_DimsTileOnly,                   &    ! all the R8 variables- internally, within CICE and doing GEOS computations in
         UNGRIDDED_DIMS     = (/NUM_SUBTILES/),                    &    ! R4. SA. Aug.2015
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = 'OCEAN:SEAICE',                      &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

   call MAPL_AddInternalSpec(GC,                                &
        SHORT_NAME         = 'VOLICE',                            &
        LONG_NAME          = 'ice_category_volume_per_unit_area_of_grid_cell',&
        UNITS              = 'm',                                 &
        PRECISION          = MAPL_R8,                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

   call MAPL_AddInternalSpec(GC,                                &
        SHORT_NAME         = 'VOLSNO',                            &
        LONG_NAME          = 'snow_category_volume_per_unit_area_of_grid_cell',&
        UNITS              = 'm',                                 &
        PRECISION          = MAPL_R8,                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

   call MAPL_AddInternalSpec(GC,                                &
        SHORT_NAME         = 'VOLPOND',                           &
        LONG_NAME          = 'pond_category_volume_per_unit_area_of_grid_cell',&
        UNITS              = 'm',                                 &
        PRECISION          = MAPL_R8,                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

   call MAPL_AddInternalSpec(GC,                                &
        SHORT_NAME         = 'APONDN',                            &
        LONG_NAME          = 'pond_concentration',                &
        UNITS              = '1',                                 &
        PRECISION          = MAPL_R8,                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

   call MAPL_AddInternalSpec(GC,                                &
        SHORT_NAME         = 'HPONDN',                            &
        LONG_NAME          = 'pond_depth',                        &
        UNITS              = 'm',                                 &
        PRECISION          = MAPL_R8,                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

   call MAPL_AddInternalSpec(GC,                                &
        SHORT_NAME         = 'ERGICE',                            &
        LONG_NAME          = 'ice_category_layer_internal_energy',&
        UNITS              = 'J m-2',                             &
        PRECISION          = MAPL_R8,                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        UNGRIDDED_DIMS     = (/NUM_ICE_LAYERS,NUM_ICE_CATEGORIES/),&
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

   call MAPL_AddInternalSpec(GC,                                &
        SHORT_NAME         = 'ERGSNO',                            &
        LONG_NAME          = 'snow_category_layer_internal_energy',&
        UNITS              = 'J m-2',                             &
        PRECISION          = MAPL_R8,                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS,NUM_ICE_CATEGORIES/),&
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

   call MAPL_AddInternalSpec(GC,                                &
        SHORT_NAME         = 'TAUAGE',                            &
        LONG_NAME          = 'volume_weighted_mean_ice_age',      &
        UNITS              = 's',                                 &
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

   call MAPL_AddInternalSpec(GC,                               &
        SHORT_NAME         = 'SLMASK',                           &
        LONG_NAME          = 'salt_water_lake_mask',             &
        UNITS              = '1',                                &
        DIMS               = MAPL_DimsTileOnly,                  &
        VLOCATION          = MAPL_VLocationNone,                 &
        DEFAULT            = 0.0,                                &
                                                       _RC  )

!-------------------Imports---------------------------------------------------------------

   call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TAUXBOT',                           &
        LONG_NAME          = 'eastward_stress_at_base_of_ice',    &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

   call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TAUYBOT',                           &
        LONG_NAME          = 'northward_stress_at_base_of_ice',   &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       _RC  )

   call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'UUA',                             &
        LONG_NAME          = 'interpolated_effective_surface_zonal_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
        _RC  )

   call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VVA',                             &
        LONG_NAME          = 'interpolated_effective_surface_meridional_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
        _RC  )


!EOS

    allocate(mystate,_STAT)
    call MAPL_GetResource (MAPL, SURFRC, label = 'SURFRC:', default = 'GEOS_SurfaceGridComp.rc', _RC)
    SCF = ESMF_ConfigCreate(_RC)
    call ESMF_ConfigLoadFile     (SCF,SURFRC,_RC)
    call MAPL_GetResource (SCF, mystate%CHOOSEMOSFC, label='CHOOSEMOSFC:', DEFAULT=1, _RC )
    call ESMF_ConfigDestroy      (SCF, _RC)
    wrap%ptr => mystate
    call ESMF_UserCompSetInternalState(gc, 'cice_private', wrap,_RC)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="INITIALIZE"   ,         _RC)

    call MAPL_TimerAdd(GC,    name="RUN1"   ,               _RC)
    call MAPL_TimerAdd(GC,    name="RUN2"  ,                _RC)

    call MAPL_TimerAdd(GC,    name="-Thermo1"    ,          _RC)
    call MAPL_TimerAdd(GC,    name="-Thermo2"    ,          _RC)
    call MAPL_TimerAdd(GC,    name="-Albedo"     ,          _RC)

    call MAPL_TimerAdd(GC,    name="FINALIZE"   ,         _RC)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC,  _RC )

! Set the Run entry point
! -----------------------

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Initialize -- Initialize method for the GEOS CICE Thermodynamic

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Initialize method of the CICE thermodynamics Gridded Component.
!   It then does a Generic\_Initialize and also CICE data structures

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

    integer                             :: NUM_SUBTILES        ! = NUM_ICE_CATEGORIES
    integer                             :: NUM_ICE_LAYERS      ! set via resource parameter
    integer                             :: NUM_ICE_CATEGORIES  ! set via resource parameter

! Local derived type aliases

    type (MAPL_MetaComp    ), pointer   :: MAPL => null()

    real                                :: DTI
    real                                :: ALBICEV, ALBSNOWV, ALBICEI, ALBSNOWI
    real                                :: USTAR_MIN, AHMAX
    real                                :: KSNO
    real                                :: ICE_REF_SALINITY
    real                                :: SNOWPATCH
    real                                :: DALB_MLT

    character(len=ESMF_MAXSTR)          :: CONDTYPE
    character(len=ESMF_MAXSTR)          :: SHORTWAVE

    integer                             :: DO_POND
    integer                             :: PRES_ICE

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, _RC )
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

    call MAPL_GetResource ( MAPL, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" ,     _RC)
    call MAPL_GetResource ( MAPL, NUM_ICE_LAYERS    , Label="CICE_N_ICE_LAYERS:"     ,     _RC)
    NUM_SUBTILES  = NUM_ICE_CATEGORIES

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"INITIALIZE")

    call MAPL_Get(MAPL, HEARTBEAT = DTI, _RC)

    call MAPL_GetResource ( MAPL, DTI,       Label="CICE_DT:",           DEFAULT=DTI,               _RC)
    call MAPL_GetResource ( MAPL, ALBICEV,   Label="ALBICEV:",           DEFAULT=0.73,              _RC)
    call MAPL_GetResource ( MAPL, ALBICEI,   Label="ALBICEI:",           DEFAULT=0.33,              _RC)
    call MAPL_GetResource ( MAPL, ALBSNOWV,  Label="ALBSNOWV:",          DEFAULT=0.96,              _RC)
    call MAPL_GetResource ( MAPL, ALBSNOWI,  Label="ALBSNOWI:",          DEFAULT=0.68,              _RC)
    call MAPL_GetResource ( MAPL, CONDTYPE,  Label="CICE_CONDUCTIVITY:", DEFAULT="bubbly",          _RC)
    call MAPL_GetResource ( MAPL, SHORTWAVE, Label="CICE_SHORTWAVE:" ,   DEFAULT="shortwave_ccsm" , _RC)
    call MAPL_GetResource ( MAPL, DO_POND,   Label="CICE_DO_POND:" ,     DEFAULT=0,                 _RC)
    call MAPL_GetResource ( MAPL, USTAR_MIN, Label="CICE_USTAR_MIN:",    DEFAULT=0.001,             _RC)
    call MAPL_GetResource ( MAPL, AHMAX,     Label="CICE_AH_MAX:",       DEFAULT=0.5,               _RC)
    call MAPL_GetResource ( MAPL, SNOWPATCH, Label="CICE_SNOW_PATCH:",   DEFAULT=0.02,              _RC)
    call MAPL_GetResource ( MAPL, DALB_MLT,  Label="CICE_DALB_MLT:",     DEFAULT=-0.075,            _RC)
    call MAPL_GetResource ( MAPL, ICE_REF_SALINITY,  Label="ICE_REF_SALINITY:" , DEFAULT=4.0,       _RC)

    ! It is desired to sometimes run the coupled model with prescribed ice.
    ! 1: prescribe ice, as in AMIP mode.
    call MAPL_GetResource ( MAPL, PRES_ICE,  Label="PRESCRIBED_ICE:" , DEFAULT=1, _RC)

    if (PRES_ICE == 1) then
       KSNO = 2.0  ! sea ice conductivity used in zero-layer ice param.
    else
       KSNO = 0.3  ! true snow conductivity
    endif

    if(MAPL_AM_I_ROOT()) then
          print*, 'Model time step = ', DTI
          print*, 'Sea ice albedo parameters:'
          print*, 'ALBICEV  = ', ALBICEV
          print*, 'ALBICEI  = ', ALBICEI
          print*, 'ALBSNOWV = ', ALBSNOWV
          print*, 'ALBSNOWI = ', ALBSNOWI
          print*, 'Sea ice conductivity parameterization:'
          print*, 'CONDTYPE = ', CONDTYPE

          if (DO_POND == 1) then
             print*, 'DO explicit melt ponding'
          else
             print*, 'DO NOT do any explicit melt ponding'
          endif

          print*, 'Sea ice shortwave parameterization:'
          print*, 'shortwave = ', SHORTWAVE
          print*, 'ustar_min = ', USTAR_MIN
          print*, 'ahmax     = ', AHMAX
    endif

    call init_column_physics(NUM_ICE_CATEGORIES,NUM_ICE_LAYERS)
    call alloc_column_physics( MAPL_AM_I_Root(), Iam )
    call input_data (DTI, ALBICEV, ALBSNOWV, ALBICEI, ALBSNOWI, CONDTYPE, USTAR_MIN, AHMAX, &
                     KSNO, ICE_REF_SALINITY, SNOWPATCH, DALB_MLT)
    call init_thermo_vertical
    call init_itd
    call init_trcr_depend(.true., (DO_POND==1))  ! 2nd argument must evaluate to a logical, e.g., TR_POND = DO_POND == 1

! Call Initialize for every Child
!--------------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  _RC)

! All Done
!---------

    call MAPL_TimerOff(MAPL,"INITIALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP
! !IROUTINE: RUN1 -- First Run stage for the Saltwater component
! !INTERFACE:

subroutine RUN1 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:
  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Periodically refreshes the sea-surface conditions

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)      :: IAm
  integer                         :: STATUS
  character(len=ESMF_MAXSTR)      :: COMP_NAME

! Locals

  integer                         :: NUM_SUBTILES        ! = NUM_ICE_CATEGORIES
  integer                         :: NUM_ICE_LAYERS      ! set via resource parameter
  integer                         :: NUM_ICE_CATEGORIES  ! set via resource parameter

  type (MAPL_MetaComp), pointer   :: MAPL => null()
  type (ESMF_Config)              :: CF
  type (ESMF_State   )            :: INTERNAL

! pointers to export

   real, pointer, dimension(:  )  :: TH     => null()
   real, pointer, dimension(:  )  :: QH     => null()
   real, pointer, dimension(:  )  :: UH     => null()
   real, pointer, dimension(:  )  :: VH     => null()
   real, pointer, dimension(:  )  :: TST    => null()
   real, pointer, dimension(:  )  :: QST    => null()
   real, pointer, dimension(:  )  :: CHT    => null()
   real, pointer, dimension(:  )  :: CMT    => null()
   real, pointer, dimension(:  )  :: CQT    => null()
   real, pointer, dimension(:  )  :: CNT    => null()
   real, pointer, dimension(:  )  :: RIT    => null()
   real, pointer, dimension(:  )  :: RET    => null()
   real, pointer, dimension(:  )  :: Z0O    => null()
   real, pointer, dimension(:  )  :: Z0H    => null()
   real, pointer, dimension(:  )  :: MOT2M  => null()
   real, pointer, dimension(:  )  :: MOQ2M  => null()
   real, pointer, dimension(:  )  :: MOU2M  => null()
   real, pointer, dimension(:  )  :: MOV2M  => null()
   real, pointer, dimension(:  )  :: MOT10M => null()
   real, pointer, dimension(:  )  :: MOQ10M => null()
   real, pointer, dimension(:  )  :: MOU10M => null()
   real, pointer, dimension(:  )  :: MOV10M => null()
   real, pointer, dimension(:  )  :: MOU50M => null()
   real, pointer, dimension(:  )  :: MOV50M => null()
   real, pointer, dimension(:  )  :: GST    => null()
   real, pointer, dimension(:  )  :: VNT    => null()
   real, pointer, dimension(:  )  :: TF     => null()
   real, pointer, dimension(:  )  :: FRACI  => null()

   real, pointer, dimension(:,:)  :: QSAT1  => null()
   real, pointer, dimension(:,:)  :: QSAT2  => null()

! pointers to internal

   real, pointer, dimension(:,:)  :: TI  => null()      ! ice skin temperature
   real, pointer, dimension(:,:)  :: QS  => null()
   real, pointer, dimension(:,:)  :: CH  => null()
   real, pointer, dimension(:,:)  :: CM  => null()
   real, pointer, dimension(:,:)  :: CQ  => null()
   real, pointer, dimension(:,:)  :: WW  => null()
   real, pointer, dimension(:,:)  :: Z0  => null()
   real, pointer, dimension(:,:)  :: TAUAGE  => null()
   real, pointer, dimension(:)    :: SLMASK  => null()
   real(kind=MAPL_R8), pointer, dimension(:,:)   :: FR      => null()
   real(kind=MAPL_R8), pointer, dimension(:,:)   :: VOLICE  => null()
   real(kind=MAPL_R8), pointer, dimension(:,:)   :: VOLSNO  => null()
   real(kind=MAPL_R8), pointer, dimension(:,:)   :: VOLPOND => null()
   real(kind=MAPL_R8), pointer, dimension(:,:,:) :: ERGICE  => null()
   real(kind=MAPL_R8), pointer, dimension(:,:,:) :: ERGSNO  => null()

! pointers to import

   real, pointer, dimension(:)    :: UU  => null()
   real, pointer, dimension(:)    :: UWINDLMTILE => null()
   real, pointer, dimension(:)    :: VWINDLMTILE => null()
   real, pointer, dimension(:)    :: UW  => null()
   real, pointer, dimension(:)    :: UI  => null()
   real, pointer, dimension(:)    :: VW  => null()
   real, pointer, dimension(:)    :: VI  => null()
   real, pointer, dimension(:)    :: DZ  => null()
   real, pointer, dimension(:)    :: TA  => null()
   real, pointer, dimension(:)    :: QA  => null()
   real, pointer, dimension(:)    :: PS  => null()
   real, pointer, dimension(:)    :: PCU => null()
   real, pointer, dimension(:)    :: FI  => null()
   real, pointer, dimension(:)    :: SW  => null()

   real, pointer, dimension(:)    :: AREA => null()
   real, pointer, dimension(:)    :: LATS_ORIGINAL => null()
   real, pointer, dimension(:)    :: LONS_ORIGINAL => null()

   integer                        :: N
   integer                        :: NT
   integer                        :: NC
   integer                        :: niter


   real, allocatable              :: TS (:,:)
   real, allocatable              :: US (:)
   real, allocatable              :: VS (:)
   real, allocatable              :: CN (:)
   real, allocatable              :: RE (:)
   real, allocatable              :: ZT (:)
   real, allocatable              :: ZQ (:)
   real, allocatable              :: UUU(:)
   real, allocatable              :: CHB(:)
   real, allocatable              :: CQB(:)
   real, allocatable              :: CMB(:)
   real, allocatable              :: UCN(:)
   real, allocatable              :: LAI(:)
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
   real, allocatable              :: PSMB(:)
   real, allocatable              :: PSL(:)
   integer, allocatable           :: IWATER(:)
   real(kind=MAPL_R8), allocatable:: FRI(:)

   real(kind=MAPL_R8), allocatable:: vice0(:)
   real(kind=MAPL_R8), allocatable:: vice1(:)
   real(kind=MAPL_R8)  ::  TOTALAREA, ALLTOTALAREA
   real(kind=MAPL_R8)  ::  TOTALAREA1, ALLTOTALAREA1
   real(kind=MAPL_R8)  ::  maxl
   real                ::  maxlat, maxlon
   type(ESMF_VM)                :: VMG
   integer                        :: i
   real                                :: LATSO, LONSO

   integer                        :: K
   logical                        :: L_STOP
   integer                        :: IDUM, JDUM
   real                           :: DT

   integer,            allocatable    :: TRCRTYPE(:)
   real(kind=MAPL_R8), allocatable    :: TRACERSDB2(:,:)
   real(kind=MAPL_R8)                 :: DTDB
   real(kind=MAPL_R8), dimension(1)   :: FRWATERDB,  FRCICEDB
   real(kind=MAPL_R8), dimension(1)   :: FHOCNLDB, FRESHLDB, FSALTLDB
   real                               :: MIN_FREEZE_SALINITY

   real(kind=MAPL_R8), allocatable    :: FR_TMP(:)
   real(kind=MAPL_R8), allocatable    :: VOLICE_TMP(:)
   real(kind=MAPL_R8), allocatable    :: VOLSNO_TMP(:)
   real(kind=MAPL_R8), allocatable    :: ERGICE_TMP(:,:)
   real(kind=MAPL_R8), allocatable    :: ERGSNO_TMP(:,:)


   real            :: ICEZ0
   real, parameter :: HPBL = 1000.

   integer         :: PRES_ICE
   integer         :: CHOOSEMOSFC
   integer         :: CHOOSEZ0
   type(cice_state_wrap) :: wrap
   type(cice_state), pointer :: mystate

! load balancing variables
   integer :: NUMMAX, pet, CICECOREBalanceHandle, L1, LN
   integer :: HorzDims, numIntSlices, numIntSlices8, numExpSlices
   real, target, allocatable :: BUFIMP(:), BUFINT(:), BUFEXP(:)
   real(kind=MAPL_R8), target, allocatable :: BUFINT8(:)
   real, pointer :: PTR1(:), PTR2(:,:), PTR3(:,:,:)
   real(kind=MAPL_R8), pointer :: PTR1R8(:), PTR2R8(:,:), PTR3R8(:,:,:)
   !integer   :: SLICESimp(100) ! increase size if more than 100 imports
   integer :: COMM
   type(ESMF_VM)                :: VM
   logical, allocatable, dimension(:) :: TILE_WITH_ICE
   logical :: loadBalance
   integer :: numUsedImp   ! number of imports actually used
   !character(len=ESMF_MAXSTR), dimension(29) :: NAMESimp
   real,               pointer    :: LATS(:)
   real,               pointer    :: LONS(:)
   integer :: NT_ORIGINAL

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run1"
    call ESMF_GridCompGet( GC, name=COMP_NAME, _RC )
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

    call MAPL_GetResource ( MAPL, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" ,     _RC)
    call MAPL_GetResource ( MAPL, NUM_ICE_LAYERS    , Label="CICE_N_ICE_LAYERS:"     ,     _RC)
    NUM_SUBTILES  = NUM_ICE_CATEGORIES

! Start Total timer
!------------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RUN1" )

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,                          &
         INTERNAL_ESMF_STATE = INTERNAL,         &
         TILEAREA = AREA,         &
         TILELATS = LATS_ORIGINAL,         &
         TILELONS = LONS_ORIGINAL,         &
         _RC )

! Get parameters (0:Louis, 1:Monin-Obukhov)
! -----------------------------------------
    call ESMF_UserCompGetInternalState(gc,'cice_private',wrap,_RC)
    mystate => wrap%ptr
    CHOOSEMOSFC = mystate%CHOOSEMOSFC

    call MAPL_GetResource ( MAPL, CHOOSEZ0,    Label="CHOOSEZ0:",    DEFAULT=3, _RC)

! Get roughness parameters with and without CICE Thermodynamics
! -------------------------------------------------------------
! icez0 value is based on literature (Ask Bin). It could be revisited, later.
    call MAPL_GetResource ( MAPL, ICEZ0,       Label="ICEZ0:" ,          DEFAULT=5.0e-4, _RC)
    call MAPL_GetResource ( MAPL, PRES_ICE,    Label="PRESCRIBED_ICE:" , DEFAULT=1,      _RC)
    call MAPL_GetResource ( MAPL, MIN_FREEZE_SALINITY, Label="MIN_FREEZE_SALINITY:" , DEFAULT=0.0,    _RC)

    call MAPL_Get(MAPL, HEARTBEAT = DT, _RC)
    call MAPL_GetResource ( MAPL, DT, Label="DT:", DEFAULT=DT, _RC)
    DTDB = REAL(DT, kind=MAPL_R8)

    call MAPL_GetResource ( MAPL, loadBalance    , Label="CICE_LOAD_BALANCE:", &
        DEFAULT=.TRUE., _RC)

   call ESMF_VMGetCurrent(VM, _RC)
   call ESMF_VMGet(VM, mpiCommunicator=COMM, localPet=pet, _RC)
   call ESMF_VMBarrier(VM, _RC)
   call MAPL_TimerOn(MAPL,    "-In_ReDist_RUN1")
   NT_ORIGINAL = size(LONS_ORIGINAL)
!load balance setup
   if(loadBalance) then

      allocate(TILE_WITH_ICE(NT_ORIGINAL), _STAT)
      TILE_WITH_ICE = .true.
      call MAPL_BalanceCreate(OrgLen=NT_ORIGINAL, Comm=COMM, Handle=CICECOREBalanceHandle, BalLen=NT, BufLen=NUMMAX, _RC)
     HorzDims = NT_ORIGINAL   ! Slice size for buffer packing

!****IMPORTANT****!!! Adjust the relevant buffer(s) and pointer assigments BufferPacking.h and BufferUnpacking.h if import/internal/export fields are added/deleted
#include "BufferPacking_RUN1.h"

   else  ! no load_balance

#include "GetPtr_RUN1.h"
      NT = NT_ORIGINAL
      LATS => LATS_ORIGINAL
      LONS => LONS_ORIGINAL

   end if
   call MAPL_TimerOff(MAPL,    "-In_ReDist_RUN1")

   NT = size(TA)
  ! if(NT == 0) then
  !    call MAPL_TimerOff(MAPL,"RUN1" )
  !    call MAPL_TimerOff(MAPL,"TOTAL")
  !    RETURN_(ESMF_SUCCESS)
  ! end if

   allocate(RE (NT)  ,   _STAT)
   allocate(CN (NT)  ,   _STAT)
   allocate(ZT (NT)  ,   _STAT)
   allocate(T2M (NT)  ,  _STAT)
   allocate(Q2M (NT)  ,  _STAT)
   allocate(U2M (NT)  ,  _STAT)
   allocate(V2M (NT)  ,  _STAT)
   allocate(T10M (NT)  , _STAT)
   allocate(Q10M (NT)  , _STAT)
   allocate(U10M (NT)  , _STAT)
   allocate(V10M (NT)  , _STAT)
   allocate(U50M (NT)  , _STAT)
   allocate(V50M (NT)  , _STAT)
   allocate(ZQ (NT)  ,   _STAT)
   allocate(UUU(NT)  ,   _STAT)
   allocate(RHO(NT) ,    _STAT)
   allocate(PSMB(NT) ,   _STAT)
   allocate(PSL(NT) ,    _STAT)
   allocate(VKH(NT) ,    _STAT)
   allocate(fakelai(NT) ,_STAT)
   allocate(VKM(NT) ,    _STAT)
   allocate(USTAR(NT) ,  _STAT)
   allocate(XX(NT)   ,   _STAT)
   allocate(YY(NT)   ,   _STAT)
   allocate(CU(NT)   ,   _STAT)
   allocate(CT(NT)   ,   _STAT)
   allocate(RIB(NT)  ,   _STAT)
   allocate(ZETA(NT) ,   _STAT)
   allocate(WS(NT)   ,   _STAT)
   allocate(IWATER(NT),  _STAT)
   allocate(LAI(NT)  ,   _STAT)
   allocate(CHB(NT)  ,   _STAT)
   allocate(CQB(NT)  ,   _STAT)
   allocate(CMB(NT)  ,   _STAT)
   allocate(UCN(NT)  ,   _STAT)
   allocate(US (NT)  ,   _STAT)
   allocate(VS (NT)  ,   _STAT)
   allocate(TS (NT,NUM_SUBTILES),   _STAT)
   allocate(FRI (NT) ,   _STAT)

#if 0
   call ESMF_GridCompGet( GC, VM=VMG, _RC )
   allocate(vice0 (NT) ,   _STAT)
   vice0 =  sum(VOLICE,dim=2)
   TOTALAREA = sum(vice0*AREA*(MAPL_RADIUS**2), mask=SLMASK<0.5)
   !deallocate(vice0)

   call ESMF_VMBarrier(VMG, _RC)
   call MAPL_CommsAllReduceSum(VMG, TOTALAREA, ALLTOTALAREA, 1, _RC)

    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run1 total ice0  = ', &
                                 ALLTOTALAREA

   call MAPL_GetResource ( MAPL, LATSO, Label="LATSO:", DEFAULT=70.0, _RC)
   call MAPL_GetResource ( MAPL, LONSO, Label="LONSO:", DEFAULT=70.0, _RC)
#endif

! do a cleanup here, in case transformation from tripolar to tile induces round-off errors
   allocate(TRCRTYPE   (NUM_3D_ICE_TRACERS),                      _STAT)
   allocate(TRACERSDB2 (NUM_3D_ICE_TRACERS , NUM_ICE_CATEGORIES), _STAT)
   allocate(FR_TMP (NUM_ICE_CATEGORIES),      _STAT)
   allocate(VOLICE_TMP (NUM_ICE_CATEGORIES),  _STAT)
   allocate(VOLSNO_TMP (NUM_ICE_CATEGORIES),  _STAT)
   allocate(ERGICE_TMP (NUM_ICE_LAYERS, NUM_ICE_CATEGORIES),   _STAT)
   allocate(ERGSNO_TMP (NUM_SNOW_LAYERS, NUM_ICE_CATEGORIES),  _STAT)

   TRCRTYPE(nt_tsfc)  = 0  ! ice/snow surface temperature
   TRCRTYPE(nt_iage)  = 1  ! volume-weighted ice age
   TRCRTYPE(nt_volpn) = 0  ! melt pond volume
   IDUM = 0
   JDUM = 0
   do k=1, NT
#if 0
       maxlat = LATS(k) * rad_to_deg
       maxlon = LONS(k) * rad_to_deg
       if((abs(maxlat-LATSO) < 1.e-3) .and. (abs(maxlon-LONSO) < 1.e-3)) then
          print*, 'before clean up'
          do i=1,NUM_ICE_CATEGORIES
             print*, FR(K,i), VOLICE(K,i)
             print*, VOLSNO(K,i), ERGSNO(K,1,i)
             if(VOLSNO(K,i) > 0.0_8) then
                print*, ERGSNO(K,1,i)/VOLSNO(K,i), (Lfresh+&
                  ERGSNO(K,1,i)/VOLSNO(K,i)/rhos)/cp_ice
             endif
          enddo
       endif
#endif
       TRACERSDB2(nt_tsfc,:) = REAL(TI(K,ICE:)-TFfresh, kind=MAPL_R8)
       TRACERSDB2(nt_iage,:) = REAL(TAUAGE(K,:),        kind=MAPL_R8)
       TRACERSDB2(nt_volpn,:)= VOLPOND(K,:)
       FRCICEDB              = sum(FR(K,ICE:))
       FRWATERDB             = REAL(1.0,kind=MAPL_R8) - FRCICEDB

       FHOCNLDB              = REAL(0.0,             kind=MAPL_R8)
       FRESHLDB              = REAL(0.0,             kind=MAPL_R8)
       FSALTLDB              = REAL(0.0,             kind=MAPL_R8)

       FR_TMP(:) = FR(K,ICE:)
       VOLICE_TMP(:) = VOLICE(K,:)
       VOLSNO_TMP(:) = VOLSNO(K,:)
       ERGICE_TMP(:,:) = ERGICE(K,:,:)
       ERGSNO_TMP(:,:) = ERGSNO(K,:,:)

       call cleanup_itd (1,1,1,1,1,1,DTDB, &
            FR_TMP,        TRACERSDB2,     &
            VOLICE_TMP,    VOLSNO_TMP,     &
            ERGICE_TMP,    ERGSNO_TMP,     &
            FRWATERDB,     FRCICEDB,       &
            TRCRTYPE,      FRESHLDB,       &
            FSALTLDB,      FHOCNLDB,       &
            .true.,        L_STOP,         &
            IDUM,            JDUM,         &
            limit_aice_in=.true.)

       _ASSERT(.not.L_STOP,'needs informative message')

       FR(K,ICE:)    = FR_TMP(:)
       VOLICE(K,:)   = VOLICE_TMP(:)
       VOLSNO(K,:)   = VOLSNO_TMP(:)
       ERGICE(K,:,:) = ERGICE_TMP(:,:)
       ERGSNO(K,:,:) = ERGSNO_TMP(:,:)

       TI(K,ICE:)   = REAL(TRACERSDB2(nt_tsfc, :)+TFfresh,   kind=MAPL_R4)
       TAUAGE(K,:)  = REAL(TRACERSDB2(nt_iage, :),           kind=MAPL_R4)
       VOLPOND(K,:) =      TRACERSDB2(nt_volpn,:)

       !if((abs(maxlat-LATSO) < 1.e-3) .and. (abs(maxlon-LONSO) < 1.e-3)) then
       !   print*, 'after  clean up'
       !   do i=1,NUM_ICE_CATEGORIES
       !      print*, FR(K,i),VOLICE(K,i)
       !   enddo
       !endif
   enddo
   deallocate(TRCRTYPE  , _STAT)
   deallocate(TRACERSDB2, _STAT)
   deallocate(FR_TMP,     _STAT)
   deallocate(VOLICE_TMP, _STAT)
   deallocate(VOLSNO_TMP, _STAT)
   deallocate(ERGICE_TMP, _STAT)
   deallocate(ERGSNO_TMP, _STAT)


#if 0
   allocate(vice1 (NT) ,   _STAT)
   vice1 =  sum(VOLICE,dim=2)
   TOTALAREA1 = sum(vice1*AREA*(MAPL_RADIUS**2), mask=SLMASK<0.5)
   !maxl = 0.0_8
   !maxlat = 0.0
   !maxlon = 0.0
   !do i=1,size(vice1)
   !  if(maxl < vice0(i)-vice1(i))then
   !        maxl = vice0(i)-vice1(i)
   !        maxlat = LATS(i) * rad_to_deg
   !        maxlon = LONS(i) * rad_to_deg
   !  endif
   !enddo
   !print*, maxl, maxlat, maxlon

   call ESMF_VMBarrier(VMG, _RC)
   call MAPL_CommsAllReduceSum(VMG, TOTALAREA1, ALLTOTALAREA1, 1, _RC)

   if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run1 total ice1  = ', &
                                 ALLTOTALAREA1

   TOTALAREA1 = sum((vice1-vice0)*AREA*(MAPL_RADIUS**2), mask=SLMASK<0.5)
   call MAPL_CommsAllReduceSum(VMG, TOTALAREA1, ALLTOTALAREA1, 1, _RC)

   if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run1 total ice1 dif = ', &
                                 ALLTOTALAREA1
   if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run1 total ice1 dif(%) = ', &
                               100*abs(ALLTOTALAREA1)/ALLTOTALAREA
   deallocate(vice0)
   deallocate(vice1)
#endif

   TS(:,ICE:) = TI           ! TI(in K): returned from CICEDyna
   US = UI
   VS = VI

   ! refresh QS based on the updated TS
   do N=1,NUM_SUBTILES
     QS(:,N) = GEOS_QSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.)
   enddo
   if(associated(QSAT1)) QSAT1 = QS

!  Clear the output tile accumulators
!------------------------------------

                       CHB = 0.0
                       CQB = 0.0
                       CMB = 0.0
   if(associated(CMT)) CMT = 0.0
   if(associated(TST)) TST = 0.0
   if(associated(QST)) QST = 0.0
   if(associated(CNT)) CNT = 0.0
   if(associated(RIT)) RIT = 0.0
   if(associated(RET)) RET = 0.0
   if(associated(TH )) TH  = 0.0
   if(associated(QH )) QH  = 0.0
   if(associated(UH )) UH  = 0.0
   if(associated(VH )) VH  = 0.0
   if(associated(Z0O)) Z0O = 0.0
   if(associated(Z0H)) Z0H = 0.0
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
   if(associated(VNT)) VNT = 0.0
   if(associated(GST)) GST = 0.0

   SUB_TILES: do N=1,NUM_SUBTILES

! Choose sfc layer: if CHOOSEMOSFC is 1 (default), choose helfand MO,
!                   if CHOOSEMOSFC is 0          , choose louis

      sfc_layer: if(CHOOSEMOSFC.eq.0) then
         call louissurface(1,N,UU,WW,PS,TA,TS,QA,QS,PCU,LAI,Z0,DZ,CM,CN,RIB,ZT,ZQ,CH,CQ,UUU,UCN,RE)

      elseif (CHOOSEMOSFC.eq.1) then

         niter = 6   ! number of internal iterations in the helfand MO surface layer routine
         IWATER= 5
         Z0(:,N)= ICEZ0

         PSMB = PS * 0.01            ! convert to MB
         fakelai  = 1.e-4
         ! Approximate pressure at top of surface layer: hydrostatic, eqn of state using avg temp and press
         PSL = PSMB * (1. - (DZ*MAPL_GRAV)/(MAPL_RGAS*(TA+TS(:,N)) ) ) /   &
                      (1. + (DZ*MAPL_GRAV)/(MAPL_RGAS*(TA+TS(:,N)) ) )

         call helfsurface( UWINDLMTILE,VWINDLMTILE,TA,TS(:,N),QA,QS(:,N),PSL,PSMB,Z0(:,N),        &
                           fakelai,IWATER,DZ,niter,nt,RHO,VKH,VKM,USTAR,XX,YY,CU,CT,RIB,ZETA,WS,  &
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
            if(associated(MOT2M ))MOT2M  = MOT2M  + T2M (:)*FR(:,N)
            if(associated(MOQ2M ))MOQ2M  = MOQ2M  + Q2M (:)*FR(:,N)
            if(associated(MOU2M ))MOU2M  = MOU2M  + U2M (:)*FR(:,N)
            if(associated(MOV2M ))MOV2M  = MOV2M  + V2M (:)*FR(:,N)

      endif sfc_layer

      !  Aggregate to tiles
      !--------------------
                             CHB     = CHB + CH(:,N)*FR(:,N)
                             CQB     = CQB + CQ(:,N)*FR(:,N)
                             CMB     = CMB + CM(:,N)*FR(:,N)
         if(associated(TST)) TST     = TST + TS(:,N)*FR(:,N)
         if(associated(QST)) QST     = QST + QS(:,N)*FR(:,N)
         if(associated(CNT)) CNT     = CNT + CN(:  )*FR(:,N)
         if(associated(RIT)) RIT     = RIT + RIB(: )*FR(:,N)
         if(associated(RET)) RET     = RET + RE(:  )*FR(:,N)
         if(associated(Z0O)) Z0O     = Z0O + Z0(:,N)*FR(:,N)
         if(associated(Z0H)) Z0H     = Z0H + ZT(:  )*FR(:,N)
         if(associated(VNT)) VNT     = VNT + UUU    *FR(:,N)

      !  Aggregate effective, CD-weighted, surface values of T and Q
      !-------------------------------------------------------------

         if(associated(TH)) TH      = TH  + CH(:,N)*TS(:,N)*FR(:,N)
         if(associated(QH)) QH      = QH  + CQ(:,N)*QS(:,N)*FR(:,N)
         if(associated(UH)) UH      = UH  + CM(:,N)*US(:  )*FR(:,N)
         if(associated(VH)) VH      = VH  + CM(:,N)*VS(:  )*FR(:,N)


      WW(:,N) = max(CH(:,N)*(TS(:,N)-TA-(MAPL_GRAV/MAPL_CP)*DZ)/TA + MAPL_VIREPS*CQ(:,N)*(QS(:,N)-QA),0.0)
      WW(:,N) = (HPBL*MAPL_GRAV*WW(:,N))**(2./3.)
      if(associated(GST)) GST     = GST + WW(:,N)*FR(:,N)
      if(associated(QSAT2)) QSAT2(:,N) = 1.0/RHO*11637800.0*exp(-5897.8/TS(:,N))

   end do SUB_TILES

   FRI = sum(FR, dim=2)

   if(associated(MOU50M)) call Normalize(MOU50M, FRI)
   if(associated(MOV50M)) call Normalize(MOV50M, FRI)
   if(associated(MOT10M)) call Normalize(MOT10M, FRI)
   if(associated(MOQ10M)) call Normalize(MOQ10M, FRI)
   if(associated(MOU10M)) call Normalize(MOU10M, FRI)
   if(associated(MOV10M)) call Normalize(MOV10M, FRI)
   if(associated(MOT2M )) call Normalize(MOT2M,  FRI)
   if(associated(MOQ2M )) call Normalize(MOQ2M,  FRI)
   if(associated(MOU2M )) call Normalize(MOU2M,  FRI)
   if(associated(MOV2M )) call Normalize(MOV2M,  FRI)

                          call Normalize(CHB,    FRI)
                          call Normalize(CQB,    FRI)
                          call Normalize(CMB,    FRI)
   if(associated(TST   )) call Normalize(TST,    FRI)
   if(associated(QST   )) call Normalize(QST,    FRI)
   if(associated(CNT   )) call Normalize(CNT,    FRI)
   if(associated(RIT   )) call Normalize(RIT,    FRI)
   if(associated(RET   )) call Normalize(RET,    FRI)
   if(associated(Z0O   )) call Normalize(Z0O,    FRI)
   if(associated(Z0H   )) call Normalize(Z0H,    FRI)
   if(associated(GST   )) call Normalize(GST,    FRI)
   if(associated(VNT   )) call Normalize(VNT,    FRI)

   if(associated(TH    )) call Normalize(TH,     FRI)
   if(associated(QH    )) call Normalize(QH,     FRI)
   if(associated(UH    )) call Normalize(UH,     FRI)
   if(associated(VH    )) call Normalize(VH,     FRI)

   if(associated(CHT   )) CHT = CHB
   if(associated(CQT   )) CQT = CQB
   if(associated(CMT   )) CMT = CMB

   if(associated(FRACI )) FRACI = REAL(FRI, KIND=MAPL_R4)

   if(associated(TF    )) call FreezingTemperature(TF, SW, MIN_FREEZE_SALINITY, PRES_ICE==1, kelvin=.true.)

    call ESMF_VMBarrier(VM, _RC)
    call MAPL_TimerOn(MAPL,    "-Out_ReDist_RUN1")
    if(loadBalance) then
#include "BufferUnpacking_RUN1.h"
       deallocate(BUFIMP,BUFINT,BUFINT8,BUFEXP,_STAT)
       deallocate(TILE_WITH_ICE, _STAT)

       call MAPL_BalanceDestroy(Handle=CICECOREBalanceHandle, _RC)
    endif
    call MAPL_TimerOff(MAPL,    "-Out_ReDist_RUN1")

   deallocate(UUU)
   deallocate(LAI)
   deallocate(CHB)
   deallocate(CQB)
   deallocate(CMB)
   deallocate(RE )
   deallocate(CN )
   deallocate(ZT )
   deallocate(ZQ )
   deallocate(UCN)
   deallocate(TS )
   deallocate(FRI)
   deallocate(VS )
   deallocate(US )
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
! !IROUTINE: RUN2 -- Second Run stage for the Saltwater component

! !INTERFACE:

subroutine RUN2 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Periodically refreshes the sea-surface conditions

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp), pointer       :: MAPL => null()
  type (ESMF_State       )            :: INTERNAL
  type (MAPL_SunOrbit)                :: ORBIT
  type (ESMF_Config      )            :: CF

  integer                             :: NUM_SUBTILES        ! = NUM_ICE_CATEGORIES
  integer                             :: NUM_ICE_LAYERS      ! set via resource parameter
  integer                             :: NUM_ICE_CATEGORIES  ! set via resource parameter

  real, pointer, dimension(:)         :: LATS_ORIGINAL => null()
  real, pointer, dimension(:)         :: LONS_ORIGINAL => null()

  real, pointer, dimension(:)         :: AREA => null()     ! needed to calculate TILEAREA in SaltWaterCore

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run2"
    call ESMF_GridCompGet( GC, name=COMP_NAME, _RC )
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

    call MAPL_GetResource ( MAPL, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" ,     _RC)
    call MAPL_GetResource ( MAPL, NUM_ICE_LAYERS    , Label="CICE_N_ICE_LAYERS:"     ,     _RC)
    NUM_SUBTILES  = NUM_ICE_CATEGORIES

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN2" )

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,             &
         TILELATS  = LATS_ORIGINAL ,                      &
         TILELONS  = LONS_ORIGINAL ,                      &
         TILEAREA  = AREA ,                      &
         ORBIT     = ORBIT,                      &
         INTERNAL_ESMF_STATE = INTERNAL,         &
         CF = CF,                                &
                                       _RC )

! Update the skin variables each step
!------------------------------------

    call CICECORE(NT_ORIGINAL=size(LONS_ORIGINAL), _RC )

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN2" )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine CICECORE(NT_ORIGINAL,RC)

   integer,           intent(IN ) :: NT_ORIGINAL
   integer, optional, intent(OUT) :: RC

!  Locals

   character(len=ESMF_MAXSTR)     :: IAm
   integer                        :: STATUS


! pointers to export

   real, pointer, dimension(:  )  :: EMISS   => null()
   real, pointer, dimension(:  )  :: ALBVF   => null()
   real, pointer, dimension(:  )  :: ALBVR   => null()
   real, pointer, dimension(:  )  :: ALBNF   => null()
   real, pointer, dimension(:  )  :: ALBNR   => null()
   real, pointer, dimension(:  )  :: EVAPOUT => null()
   real, pointer, dimension(:  )  :: SUBLIM  => null()
   real, pointer, dimension(:  )  :: SNOWOCN => null()
   real, pointer, dimension(:  )  :: RAINOCN => null()
   real, pointer, dimension(:  )  :: SHICE   => null()
   real, pointer, dimension(:  )  :: SHOUT   => null()
   real, pointer, dimension(:  )  :: HLATN   => null()
   real, pointer, dimension(:  )  :: HLATWTR => null()
   real, pointer, dimension(:  )  :: HLATICE => null()
   real, pointer, dimension(:  )  :: FSURFe  => null()
   real, pointer, dimension(:  )  :: FSURFICE=> null()
   real, pointer, dimension(:  )  :: HLWUP   => null()
   real, pointer, dimension(:  )  :: HLWUPe  => null()
   real, pointer, dimension(:  )  :: LWNDSRF => null()
   real, pointer, dimension(:  )  :: SWNDSRF => null()
   real, pointer, dimension(:  )  :: SWNDICE => null()
   real, pointer, dimension(:  )  :: LWNDICE => null()

   real, pointer, dimension(:  )  :: DELTS  => null()
   real, pointer, dimension(:  )  :: DELQS  => null()
   real, pointer, dimension(:  )  :: TST    => null()
   real, pointer, dimension(:  )  :: QST    => null()
   real, pointer, dimension(:  )  :: TAUXW  => null()
   real, pointer, dimension(:  )  :: TAUYW  => null()
   real, pointer, dimension(:  )  :: TAUXI  => null()
   real, pointer, dimension(:  )  :: TAUYI  => null()
   real, pointer, dimension(:  )  :: TAUXO  => null()
   real, pointer, dimension(:  )  :: TAUYO  => null()
   real, pointer, dimension(:  )  :: USTR3  => null()
   real, pointer, dimension(:  )  :: UUEX   => null()
   real, pointer, dimension(:  )  :: PSEX   => null()
   real, pointer, dimension(:  )  :: PENUVR => null()
   real, pointer, dimension(:  )  :: PENUVF => null()
   real, pointer, dimension(:  )  :: PENPAR => null()
   real, pointer, dimension(:  )  :: PENPAF => null()
   real, pointer, dimension(:  )  :: FRACI  => null()
   real, pointer, dimension(:  )  :: FRACINEW => null()

   real, pointer, dimension(:  )  :: FRAZIL     => null()          ! CICE related exports.
   real, pointer, dimension(:  )  :: CONGELO    => null()
   real, pointer, dimension(:  )  :: SNOICEO    => null()
   real, pointer, dimension(:  )  :: FRESH      => null()
   real, pointer, dimension(:  )  :: FSALT      => null()
   real, pointer, dimension(:  )  :: FHOCN      => null()
   real, pointer, dimension(:  )  :: PICE       => null()
   real, pointer, dimension(:  )  :: FSWTRUO    => null()
   real, pointer, dimension(:  )  :: FSWABSO    => null()
   real, pointer, dimension(:  )  :: MELTL      => null()
   real, pointer, dimension(:  )  :: MELTTL     => null()
   real, pointer, dimension(:  )  :: MELTBL     => null()
   real, pointer, dimension(:  )  :: MELTSL     => null()
   real, pointer, dimension(:  )  :: HICE       => null()
   real, pointer, dimension(:  )  :: HSNO       => null()
   real, pointer, dimension(:  )  :: HICEUNT    => null()
   real, pointer, dimension(:  )  :: SNOONICE   => null()
   real, pointer, dimension(:  )  :: TSKINICE   => null()
   real, pointer, dimension(:  )  :: IAGE       => null()
   real, pointer, dimension(:  )  :: DAIDTT     => null()
   real, pointer, dimension(:  )  :: DVIDTT     => null()
   real, pointer, dimension(:  )  :: FBOTL      => null()
   real, pointer, dimension(:  )  :: USTARI     => null()
   real, pointer, dimension(:  )  :: FCONDTOP   => null()
   real, pointer, dimension(:  )  :: FCONDB     => null()
   real, pointer, dimension(:  )  :: NIERG      => null()
   real, pointer, dimension(:  )  :: SBLXOUT    => null()
   real, pointer, dimension(:  )  :: SIALB      => null()
   real, pointer, dimension(:  )  :: GHTSKIN    => null()
   real, pointer, dimension(:)    :: FRZMLTe    => null()
   real, pointer, dimension(:)    :: LWDNSRFe   => null()
   real, pointer, dimension(:)    :: SWDNSRFe   => null()

   ! pointers to category dimensional exports (CICE)
   real, pointer, dimension(:,:)  :: SHICEN      => null()
   real, pointer, dimension(:,:)  :: HLWUPN      => null()
   real, pointer, dimension(:,:)  :: FSURFN      => null()
   real, pointer, dimension(:,:)  :: TSURFN      => null()
   real, pointer, dimension(:,:)  :: LWNDSRFN    => null()
   real, pointer, dimension(:,:)  :: FSWSFCN     => null()
   real, pointer, dimension(:,:)  :: ALBINe      => null()
   real, pointer, dimension(:,:)  :: ALBSNe      => null()
   real, pointer, dimension(:,:)  :: FCONDBOTN   => null()
   real, pointer, dimension(:,:)  :: FCONDTOPN   => null()
   real, pointer, dimension(:,:)  :: TINZ        => null()

   ! pointers to CMIP5 exports (CICE)
   real, pointer, dimension(:  )  :: EVAP_C5        => null()
   real, pointer, dimension(:  )  :: PR_C5          => null()
   real, pointer, dimension(:  )  :: PRSN_C5        => null()
   real, pointer, dimension(:  )  :: GRFRAZIL_C5    => null()
   real, pointer, dimension(:  )  :: GRCONGEL_C5    => null()
   real, pointer, dimension(:  )  :: GRLATERAL_C5   => null()
   real, pointer, dimension(:  )  :: SNOTOICE_C5    => null()
   real, pointer, dimension(:  )  :: SNOMELT_C5     => null()
   real, pointer, dimension(:  )  :: TMELT_C5       => null()
   real, pointer, dimension(:  )  :: BMELT_C5       => null()
   real, pointer, dimension(:  )  :: SFDSI_C5       => null()
   real, pointer, dimension(:  )  :: HFSIFRAZIL_C5  => null()
   real, pointer, dimension(:  )  :: IALB_C5        => null()
   real, pointer, dimension(:  )  :: RSDSSI_C5      => null()
   real, pointer, dimension(:  )  :: RSUSSI_C5      => null()
   real, pointer, dimension(:  )  :: FSITHERM_CMIP5 => null()

! pointers to internal

   real, pointer, dimension(:,:)  :: TI    => null()
   real, pointer, dimension(:  )  :: HI    => null()
   real, pointer, dimension(:  )  :: SI    => null()
   real, pointer, dimension(:,:)  :: QS    => null()
   real, pointer, dimension(:,:)  :: CH    => null()
   real, pointer, dimension(:,:)  :: CQ    => null()
   real, pointer, dimension(:,:)  :: CM    => null()

   real(kind=MAPL_R8), pointer, dimension(:,:)   :: FR8     => null()     ! CICE related
   real(kind=MAPL_R8), pointer, dimension(:,:)   :: VOLICE  => null()
   real(kind=MAPL_R8), pointer, dimension(:,:)   :: VOLSNO  => null()
   real(kind=MAPL_R8), pointer, dimension(:,:)   :: VOLPOND => null()
   real(kind=MAPL_R8), pointer, dimension(:,:)   :: APONDN  => null()
   real(kind=MAPL_R8), pointer, dimension(:,:)   :: HPONDN  => null()
   real(kind=MAPL_R8), pointer, dimension(:,:,:) :: ERGICE  => null()
   real(kind=MAPL_R8), pointer, dimension(:,:,:) :: ERGSNO  => null()

   real, pointer, dimension(:,:)   :: TAUAGE => null()
   real, pointer, dimension(:)     :: SLMASK => null()

! pointers to import

   real, pointer, dimension(:)    :: ALW => null()
   real, pointer, dimension(:)    :: BLW => null()
   real, pointer, dimension(:)    :: LWDNSRF => null()
   real, pointer, dimension(:)    :: DRPAR => null()
   real, pointer, dimension(:)    :: DFPAR => null()
   real, pointer, dimension(:)    :: DRNIR => null()
   real, pointer, dimension(:)    :: DFNIR => null()
   real, pointer, dimension(:)    :: DRUVR => null()
   real, pointer, dimension(:)    :: DFUVR => null()
   real, pointer, dimension(:)    :: EVAP  => null()
   real, pointer, dimension(:)    :: SH => null()
   real, pointer, dimension(:)    :: TAUX => null()
   real, pointer, dimension(:)    :: TAUY => null()
   real, pointer, dimension(:)    :: DEV => null()
   real, pointer, dimension(:)    :: DSH => null()
   real, pointer, dimension(:)    :: SNO => null()
   real, pointer, dimension(:)    :: PLS => null()
   real, pointer, dimension(:)    :: PCU => null()
   real, pointer, dimension(:)    :: PS => null()
   real, pointer, dimension(:)    :: UU => null()
!!$   real, pointer, dimension(:)    :: TF => null()
   real, pointer, dimension(:)    :: FI => null()
   real, pointer, dimension(:)    :: THATM => null()
   real, pointer, dimension(:)    :: QHATM => null()
   real, pointer, dimension(:)    :: UHATM => null()
   real, pointer, dimension(:)    :: VHATM => null()
   real, pointer, dimension(:)    :: UUA   => null()
   real, pointer, dimension(:)    :: VVA   => null()
   real, pointer, dimension(:)    :: CTATM => null()
   real, pointer, dimension(:)    :: CQATM => null()
   real, pointer, dimension(:)    :: CMATM => null()
   real, pointer, dimension(:)    :: UW => null()
   real, pointer, dimension(:)    :: UI => null()
   real, pointer, dimension(:)    :: VW => null()
   real, pointer, dimension(:)    :: VI => null()

   real, pointer, dimension(:)    :: TAUXBOT   => null()                  ! CICE related
   real, pointer, dimension(:)    :: TAUYBOT   => null()
   real, pointer, dimension(:)    :: TW        => null()
   real, pointer, dimension(:)    :: SW        => null()
   real, pointer, dimension(:)    :: FRZMLT    => null()

   real, pointer, dimension(:,:)       :: TS  => null()
   real, allocatable,    dimension(:)              :: SHF
   real, allocatable,    dimension(:)              :: EVP
   real, allocatable,    dimension(:)              :: SHD
   real, allocatable,    dimension(:)              :: EVD
   real, allocatable,    dimension(:)              :: CFQ
   real, allocatable,    dimension(:)              :: CFT
   !real, allocatable,    dimension(:)              :: UUA
   !real, allocatable,    dimension(:)              :: VVA
   real, allocatable,    dimension(:)              :: TXI
   real, allocatable,    dimension(:)              :: TYI
   real, allocatable,    dimension(:)              :: DQS
   real, allocatable,    dimension(:)              :: DTS
   real, allocatable,    dimension(:)              :: DTX
   real, allocatable,    dimension(:)              :: DTY
   real, allocatable,    dimension(:)              :: SWN
   real, allocatable,    dimension(:)              :: PEN
   real, allocatable,    dimension(:)              :: LHF
   real, allocatable,    dimension(:)              :: ZTH
   real, allocatable,    dimension(:)              :: SLR
   real, allocatable,    dimension(:)              :: VSUVR
   real, allocatable,    dimension(:)              :: VSUVF

   integer                             :: N
   real                                :: DT
   real                                :: MAXSALINITY
   real                                :: MINSALINITY

! following are related  to CICE

   integer                             :: DIAG_ICE_BUDGET           ! default (=0) is to not compute certain flux budgets over Sea Ice
   integer                             :: NSUB, I, K, L
   integer                             :: DO_POND

   real                                :: LATSO, LONSO
   real,    dimension(1)               :: LATSD, LONSD
   !real                                :: TOTALAREA, ALLTOTALAREA
   logical, dimension(1)               :: OBSERVE

   real,    dimension(1)               :: FRZ_ONSET, MLT_ONSET
   real,    dimension(1)               :: RDUM
   real                                :: FRZMLT_MAX
   real(kind=MAPL_R8)                  :: DTDB
   real(kind=MAPL_R8), dimension(1)    :: FRZMLTDB, TSCDB, TFDB, TAUXBOTDB, TAUYBOTDB, &
                                          TBOTDB, FBOTDB, RSIDEDB

   real, allocatable,   dimension(:)              :: FSWABS
   real                                           :: YDAY
   real, allocatable,  dimension(:)               :: ALBVRI
   real, allocatable,   dimension(:)              :: ALBVFI
   real, allocatable,   dimension(:)              :: ALBNRI
   real, allocatable,   dimension(:)              :: ALBNFI

   integer,            allocatable    :: TRCRTYPE      (:)
   real,               allocatable    :: TRACERS       (:,:)
   real,               allocatable    :: TF            (:)
   real,               allocatable    :: MELTLN        (:)
   real,               allocatable    :: FRAZLN        (:)
   real,               allocatable    :: FRESHN        (:)
   real,               allocatable    :: FRESHL        (:)
   real,               allocatable    :: FSALTN        (:)
   real,               allocatable    :: FSALTL        (:)
   real,               allocatable    :: FHOCNN        (:)
   real,               allocatable    :: FHOCNL        (:)
   real,               allocatable    :: RSIDE         (:)
   real,               allocatable    :: FSWTHRU       (:,:)        ! FSWTHRU is also an EXPORT
   real,               allocatable    :: FCOND         (:,:)
   real,               allocatable    :: FCONDBOT      (:,:)
   real,               allocatable    :: TBOT          (:)
   real,               allocatable    :: FBOT          (:)
   real,               allocatable    :: ALBVRN        (:,:)
   real,               allocatable    :: ALBNRN        (:,:)
   real,               allocatable    :: ALBVFN        (:,:)
   real,               allocatable    :: ALBNFN        (:,:)
   real,               allocatable    :: FSWSFC        (:,:)
   real,               allocatable    :: FSWINT        (:,:)
   real,               allocatable    :: ISWABS        (:,:,:)
   real,               allocatable    :: SSWABS        (:,:,:)
   real,               allocatable    :: MELTT         (:)
   real,               allocatable    :: MELTS         (:)
   real,               allocatable    :: MELTB         (:)
   real,               allocatable    :: CONGEL        (:)
   real,               allocatable    :: SNOICE        (:)

   real,               allocatable    :: TS_OLD        (:,:)

   real,               allocatable    :: ALBIN         (:,:)
   real,               allocatable    :: ALBSN         (:,:)
   real,               allocatable    :: ALBPND        (:,:)

   real,               allocatable    :: DRUVRTHRU     (:,:)
   real,               allocatable    :: DFUVRTHRU     (:,:)
   real,               allocatable    :: DRPARTHRU     (:,:)
   real,               allocatable    :: DFPARTHRU     (:,:)

   real,               allocatable    :: TOTALFLUX     (:)
   real,               allocatable    :: NEWICEERG     (:) ! newly generated ice energy <=0 (W m-2)
   real,               allocatable    :: SBLX          (:) !
   real,               allocatable    :: FSURF         (:) !

!  Following arrays have to be R8 for CICE
   real(kind=MAPL_R8), allocatable     :: AICENINIT    (:,:)
   real(kind=MAPL_R8), allocatable     :: VICENINIT    (:,:)
   real(kind=MAPL_R8), allocatable     :: FRCICE       (:)
   real(kind=MAPL_R8), allocatable     :: FR_OLD       (:)
   real(kind=MAPL_R8), allocatable     :: VOLSNO_OLD   (:)
   real(kind=MAPL_R8), allocatable     :: VOLICE_OLD   (:)
   real(kind=MAPL_R8), allocatable     :: VOLICE_DELTA (:,:)
   real(kind=MAPL_R8), allocatable     :: FR8TMP       (:,:)
   real(kind=MAPL_R8)                  :: ERGICE_TMP(NUM_ICE_LAYERS,  NUM_ICE_CATEGORIES)
   real(kind=MAPL_R8)                  :: ERGSNO_TMP(NUM_SNOW_LAYERS, NUM_ICE_CATEGORIES)

   real(kind=MAPL_R8), allocatable     :: TEMPVOLICE(:,:)
   real, pointer                       :: DELTAVOL1 (:,:) => null()

!  -------------------------------------------------------------------

   type (ESMF_TimeInterval)            :: DELT
   real                                :: DT_SOLAR
   type (ESMF_TimeInterval)            :: TINT
   type (ESMF_Time)                    :: CURRENT_TIME
   type (ESMF_Time)                    :: NOW
   type (ESMF_Time)                    :: BEFORE
   type (ESMF_Time)                    :: MODELSTART
   type (ESMF_Alarm)                   :: SOLALARM
   logical                             :: solalarmison
   type(ESMF_VM)                       :: VM
   type(ESMF_VM)                       :: VMG               ! for CICE
   integer                             :: MYPE              ! for CICE
   logical                             :: debugzth

   real(kind=MAPL_R8)                  :: TOTALAREAN, ALLTOTALAREAN
   real(kind=MAPL_R8)                  :: TOTALAREAS, ALLTOTALAREAS
   real(kind=MAPL_R8)                  :: TOTALAREAN1, ALLTOTALAREAN1
   real(kind=MAPL_R8)                  :: TOTALAREAS1, ALLTOTALAREAS1
   real(kind=MAPL_R8),       allocatable     :: vice0(:)
   real(kind=MAPL_R8),       allocatable     :: vice1(:)

   real                                :: EMSICE

   real, parameter                     :: SALTWATERCAP    = MAPL_CAPWTR
   real, parameter                     :: SALTWATERICECAP = MAPL_CAPICE

! load balancing variables
   integer :: NT, NUMMAX, pet, CICECOREBalanceHandle, L1, LN
   integer :: HorzDims, numIntSlices, numIntSlices8, numExpSlices
   real, target, allocatable :: BUFIMP(:), BUFINT(:), BUFEXP(:)
   real(kind=MAPL_R8), target, allocatable :: BUFINT8(:)
   real, pointer :: PTR1(:), PTR2(:,:), PTR3(:,:,:)
   real(kind=MAPL_R8), pointer :: PTR1R8(:), PTR2R8(:,:), PTR3R8(:,:,:)
   !integer   :: SLICESimp(100) ! increase size if more than 100 imports
   integer :: COMM
   logical, dimension(NT_ORIGINAL) :: TILE_WITH_ICE
   logical :: loadBalance
   integer :: numUsedImp   ! number of imports actually used
   !character(len=ESMF_MAXSTR), dimension(29) :: NAMESimp
   real,               pointer    :: LATS(:)
   real,               pointer    :: LONS(:)

!  Begin...
!----------

   IAm =  trim(COMP_NAME) // "CICECORE"

! Get the time step
! -----------------

    call MAPL_Get(MAPL, HEARTBEAT = DT, _RC)
    call MAPL_GetResource ( MAPL, DT, Label="DT:", DEFAULT=DT, _RC)

! Get parameters
! --------------

    call MAPL_GetResource ( MAPL, EMSICE,      Label="CICE_EMSICE:",   DEFAULT=0.99999, _RC)

    call MAPL_GetResource ( MAPL, MAXSALINITY, Label="MAX_SALINITY:" , DEFAULT=40.0 ,   _RC)
    call MAPL_GetResource ( MAPL, MINSALINITY, Label="MIN_SALINITY:" , DEFAULT=5.0 ,    _RC)
    call MAPL_GetResource ( MAPL, DO_POND,     Label="CICE_DO_POND:" , DEFAULT=0,       _RC)

    call MAPL_GetResource ( MAPL, loadBalance    , Label="CICE_LOAD_BALANCE:", &
        DEFAULT=.TRUE., _RC)

   call MAPL_GetPointer(EXPORT,EMISS  , 'EMIS' , alloc=.true., _RC)
   call MAPL_GetPointer(EXPORT,ALBVF  , 'ALBVF', alloc=.true., _RC)
   call MAPL_GetPointer(EXPORT,ALBVR  , 'ALBVR', alloc=.true., _RC)
   call MAPL_GetPointer(EXPORT,ALBNF  , 'ALBNF', alloc=.true., _RC)
   call MAPL_GetPointer(EXPORT,ALBNR  , 'ALBNR', alloc=.true., _RC)

   call ESMF_VMGetCurrent(VM, _RC)
   call ESMF_VMGet(VM, mpiCommunicator=COMM, localPet=pet, _RC)
   call ESMF_VMBarrier(VM, _RC)
   call MAPL_TimerOn(MAPL,    "-In_ReDist")
!load balance setup
   if(loadBalance) then

      TILE_WITH_ICE = .true.
      call MAPL_BalanceCreate(OrgLen=NT_ORIGINAL, Comm=COMM, Handle=CICECOREBalanceHandle, BalLen=NT, BufLen=NUMMAX, _RC)
     HorzDims = NT_ORIGINAL   ! Slice size for buffer packing

!****IMPORTANT****!!! Adjust the relevant buffer(s) and pointer assigments BufferPacking.h and BufferUnpacking.h if import/internal/export fields are added/deleted
#include "BufferPacking.h"

   else  ! no load_balance

#include "GetPtr.h"
      NT = NT_ORIGINAL
      LATS => LATS_ORIGINAL
      LONS => LONS_ORIGINAL

   end if
   call MAPL_TimerOff(MAPL,    "-In_ReDist")

! Copy friendly internals into tile-tile local variables
!-------------------------------------------------------

    TS => TI
    allocate( FSWABS (NT), _STAT)
    allocate( ALBVRI (NT), _STAT)
    allocate( ALBVFI (NT), _STAT)
    allocate(  ALBNRI (NT), _STAT)
    allocate( ALBNFI (NT), _STAT)
    allocate(SHF        (NT),                                   _STAT)
    allocate(EVP        (NT),                                   _STAT)
    allocate(SHD        (NT),                                   _STAT)
    allocate(EVD        (NT),                                   _STAT)
    allocate(CFQ        (NT),                                   _STAT)
    allocate(CFT        (NT),                                   _STAT)
    allocate(TXI        (NT),                                   _STAT)
    allocate(TYI        (NT),                                   _STAT)
    allocate(DQS        (NT),                                   _STAT)
    allocate(DTS        (NT),                                   _STAT)
    allocate(DTX        (NT),                                   _STAT)
    allocate(DTY        (NT),                                   _STAT)
    allocate(SWN        (NT),                                   _STAT)
    allocate(PEN        (NT),                                   _STAT)
    allocate(LHF        (NT),                                   _STAT)
    allocate(ZTH        (NT),                                   _STAT)
    allocate(SLR        (NT),                                   _STAT)
    allocate(VSUVR        (NT),                                   _STAT)
    allocate(VSUVF        (NT),                                   _STAT)

! Initialize PAR and UVR beam fluxes
!-----------------------------------

    VSUVR = DRPAR + DRUVR
    VSUVF = DFPAR + DFUVR

    if(associated(SWDNSRFe)) SWDNSRFe = VSUVR+VSUVF+DRNIR+DFNIR
    if(associated(LWDNSRFe)) LWDNSRFe = LWDNSRF

!     allocate arrays for CICE Thermodynamics
    allocate(TRCRTYPE  (NUM_3D_ICE_TRACERS),                   _STAT)
    allocate(TRACERS   (NUM_3D_ICE_TRACERS,NUM_ICE_CATEGORIES),_STAT)
    allocate(TF        (NT),                                   _STAT)
    allocate(MELTLN    (NT),                                   _STAT)
    allocate(FRAZLN    (NT),                                   _STAT)
    allocate(FRESHN    (NT),                                   _STAT)
    allocate(FRESHL    (NT),                                   _STAT)
    allocate(FSALTN    (NT),                                   _STAT)
    allocate(FSALTL    (NT),                                   _STAT)
    allocate(FHOCNN    (NT),                                   _STAT)
    allocate(FHOCNL    (NT),                                   _STAT)
    allocate(RSIDE     (NT),                                   _STAT)
    allocate(FSWTHRU   (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(FCOND     (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(FCONDBOT  (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(TBOT      (NT),                                   _STAT)
    allocate(FBOT      (NT),                                   _STAT)
    allocate(ALBVRN    (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(ALBNRN    (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(ALBVFN    (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(ALBNFN    (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(FSWSFC    (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(FSWINT    (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(ISWABS    (NT,NUM_ICE_LAYERS,NUM_ICE_CATEGORIES), _STAT)
    allocate(SSWABS    (NT,NUM_SNOW_LAYERS,NUM_ICE_CATEGORIES),_STAT)
    allocate(MELTT     (NT),                                   _STAT)
    allocate(MELTS     (NT),                                   _STAT)
    allocate(MELTB     (NT),                                   _STAT)
    allocate(CONGEL    (NT),                                   _STAT)
    allocate(SNOICE    (NT),                                   _STAT)
    allocate(TS_OLD    (NT,NUM_SUBTILES),                      _STAT)
    allocate(ALBIN     (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(ALBSN     (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(ALBPND    (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(DRUVRTHRU (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(DFUVRTHRU (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(DRPARTHRU (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(DFPARTHRU (NT,NUM_ICE_CATEGORIES),                _STAT)
    allocate(TOTALFLUX (NT),                                   _STAT)
    allocate(NEWICEERG (NT),                                   _STAT)
    allocate(SBLX      (NT),                                   _STAT)
    allocate(FSURF     (NT),                                   _STAT)
    allocate(AICENINIT (NT, NUM_ICE_CATEGORIES),               _STAT)
    allocate(VICENINIT (NT, NUM_ICE_CATEGORIES),               _STAT)
    allocate(FR8TMP    (NT, NUM_SUBTILES),                     _STAT)
    allocate(FRCICE    (NT),                                   _STAT)
    allocate(FR_OLD    (NT),                                   _STAT)
    allocate(VOLSNO_OLD(NT),                                   _STAT)
    allocate(VOLICE_OLD(NT),                                   _STAT)
    allocate(VOLICE_DELTA(NT,NUM_ICE_CATEGORIES),              _STAT)

    call MAPL_GetResource ( MAPL, LATSO, Label="LATSO:", DEFAULT=70.0, _RC)
    call MAPL_GetResource ( MAPL, LONSO, Label="LONSO:", DEFAULT=70.0, _RC)

    ! Tracking changes in VOLICE
    TEMPVOLICE = VOLICE

!     initialize arrays for CICE Thermodynamics
    call CICE_PREP_THERMO(TF,TRCRTYPE,TRACERS,MELTLN,FRAZLN,FRESHN,FRESHL,FSALTN,FSALTL,FHOCNN,FHOCNL,RSIDE,  &
                            FSWTHRU,FCOND,FCONDBOT,TBOT,FBOT,ALBIN,ALBSN,ALBPND,ALBVRN,ALBNRN,ALBVFN,ALBNFN,FSWSFC,FSWINT,     &
                            ISWABS,SSWABS,FSWABS,MELTT,MELTS,MELTB,CONGEL,SNOICE,UW,VW,SLMASK,LATS,LONS,LATSO,LONSO,   &
                            FR8,FRCICE,SW,TAUAGE,ICE,NT,VOLPOND,DT,VOLICE,VOLSNO,ERGICE,ERGSNO,TS,VOLICE_DELTA,  &
                            NEWICEERG, SBLX, _RC)

    ! Tracking changes in VOLICE
    if (associated(DELTAVOL1)) DELTAVOL1 = VOLICE - TEMPVOLICE


    FR_OLD     = FRCICE   ! FRCICE is initialized by above subroutine CICE_PREP_THERMO
    TS_OLD     = TS
    VOLICE_OLD = sum(VOLICE,dim=2)
    VOLSNO_OLD = sum(VOLSNO,dim=2)

    AICENINIT  = FR8(:,ICE:)
    VICENINIT  = VOLICE

#if 0
    allocate(vice0(NT),                                   _STAT)
    call ESMF_GridCompGet( GC, VM=VMG, _RC )

    vice0 =  VOLICE_OLD
    TOTALAREAN = sum(vice0*AREA*(MAPL_RADIUS**2), mask=SLMASK<0.5 .and. LATS>0.0)
    TOTALAREAS = sum(vice0*AREA*(MAPL_RADIUS**2), mask=SLMASK<0.5 .and. LATS<0.0)
    deallocate(vice0)

    call ESMF_VMBarrier(VMG, _RC)
    call MAPL_CommsAllReduceSum(VMG, TOTALAREAN, ALLTOTALAREAN, 1, _RC)
    call MAPL_CommsAllReduceSum(VMG, TOTALAREAS, ALLTOTALAREAS, 1, _RC)

    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run2 North ice0  = ', &
                                 ALLTOTALAREAN
    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run2 South ice0  = ', &
                                 ALLTOTALAREAS
#endif

    call MAPL_GetResource ( MAPL, FRZMLT_MAX, Label="CICE_FRZMLT_MAX:" , DEFAULT=1000., _RC)
    DTDB = REAL(DT, kind=MAPL_R8)     ! Convert DT precision: Real4 to Real8 for usage in CICE
    do k=1, NT
          call CICE_INQUIRE_TILE(LATS(K), LONS(K), LATSO, LONSO, OBSERVE, LATSD, LONSD)
          FRZMLTDB  = REAL(FRZMLT(K),            kind=MAPL_R8)
          FRZMLTDB  = min(max(FRZMLTDB,-FRZMLT_MAX),FRZMLT_MAX)
          if(FRZMLT(K)<0.0) then ! heat the already existing ice from below
             TSCDB     = REAL(TW(K)-TFfresh,        kind=MAPL_R8)
             TFDB      = REAL(TF(K),                kind=MAPL_R8)
             TAUXBOTDB = REAL(TAUXBOT(K),           kind=MAPL_R8)
             TAUYBOTDB = REAL(TAUYBOT(K),           kind=MAPL_R8)
             ERGICE_TMP(:,:) = ERGICE(K,:,:)
             ERGSNO_TMP(:,:) = ERGSNO(K,:,:)

             call frzmlt_bottom_lateral (1,1,1,1,1,1,DTDB, &
                  OBSERVE,                                 &
                  FRCICE(K),        FRZMLTDB,       & ! in
                  ERGICE_TMP,       ERGSNO_TMP,     & ! in
                  TSCDB,            TFDB,           & ! in
                  TAUXBOTDB,        TAUYBOTDB,      & ! in
                  TBOTDB, FBOTDB,   RSIDEDB      )    ! out

                  TBOT(K)  =  TBOTDB(1)
                  FBOT(K)  =  FBOTDB(1)
                  RSIDE(K) =  RSIDEDB(1)
          else
                  TBOT(K)  =  TF(K)
                  FBOT(K)  =  0.0
                  RSIDE(K) =  0.0
          endif
    enddo !k

    if(associated(FRZMLTe))   FRZMLTe = FRZMLT

!     Output additional CICE diagnostics?
    call MAPL_GetResource ( MAPL, DIAG_ICE_BUDGET, Label="DIAG_ICE_BUDGET:" , DEFAULT=0    , _RC)

    if (DIAG_ICE_BUDGET /= 0) then
        call ESMF_GridCompGet( GC, VM=VMG, _RC )
        call ESMF_VMGet      (VMG, localpet=MYPE,  _RC)
    endif

    if(associated(RSDSSI_C5))   RSDSSI_C5  = FRCICE * (VSUVR + VSUVF + DRNIR + DFNIR)

    call MAPL_TimerOn(MAPL,    "-Albedo")

    debugzth = .false.

    call ESMF_VMGetCurrent ( VM, _RC )

        ! --------------------------------------------------------------------------
        ! Get the current time.
        ! --------------------------------------------------------------------------

    call ESMF_ClockGet( CLOCK, currTime=CURRENT_TIME, startTime=MODELSTART, TIMESTEP=DELT,  _RC )
    if (MAPL_AM_I_Root(VM).and.debugzth) then
      print *,' start time of clock '
      CALL ESMF_TimePrint ( MODELSTART, OPTIONS="string", _RC )
    endif

        ! --------------------------------------------------------------------------
        ! retrieve the zenith angle
        ! --------------------------------------------------------------------------

!! The next sequence is to make sure that the albedo here and in solar are in sync
!!
! Need to know when Solar was called last, so first get the solar alarm
        call ESMF_ClockGetAlarm ( CLOCK, alarmname="SOLAR_Alarm", ALARM=SOLALARM, _RC )
! Get the interval of the solar alarm - first get it in seconds
        call ESMF_ConfigGetAttribute ( CF, DT_SOLAR, Label="SOLAR_DT:", DEFAULT=DT, _RC )
! Now make an ESMF interval from the increment in seconds
        CALL ESMF_TimeIntervalSet ( TINT, S=NINT(DT_SOLAR), _RC )
! Now print out the solar alarm interval
        if (MAPL_AM_I_Root(VM).and.debugzth) CALL ESMF_TimeIntervalPrint ( TINT, OPTIONS="string", _RC )
! Now find out if it is ringing now: if so, set "BEFORE" to last time it rang before now
         solalarmison = ESMF_AlarmIsRinging(SOLALARM,_RC)
         if (MAPL_AM_I_Root(VM).and.debugzth)print *,' logical for solar alarm ',solalarmison
!     if so, set "BEFORE" to last time it rang before now
        if(solalarmison) then
         if (MAPL_AM_I_Root(VM).and.debugzth)print *,' In catch, solar alarm is ringing '
         NOW = CURRENT_TIME
         BEFORE = NOW - TINT
! Now print out the last time solar alarm rang
         if (MAPL_AM_I_Root(VM).and.debugzth)CALL ESMF_TimePrint ( BEFORE, OPTIONS="string", _RC )
!     If alarm is not ringing now, find out when it rang last
        else
         if (MAPL_AM_I_Root(VM).and.debugzth)print *,' In catch, solar alarm is not ringing '
         call ESMF_AlarmGet ( SOLALARM, prevRingTime=BEFORE, _RC )
! PrevRingTime can lie: if alarm never went off yet it gives next alarm time, not prev.
         if(BEFORE > CURRENT_TIME) then
          BEFORE = BEFORE-TINT
          if (MAPL_AM_I_Root(VM).and.debugzth)print *,' In catch, solar alarm not ringing, prev time lied '
          if (MAPL_AM_I_Root(VM).and.debugzth)CALL ESMF_TimePrint ( BEFORE, OPTIONS="string", _RC )
         else
          if (MAPL_AM_I_Root(VM).and.debugzth)print *,' In catch, solar alarm not ringing, prev time okay '
          if (MAPL_AM_I_Root(VM).and.debugzth)CALL ESMF_TimePrint ( BEFORE, OPTIONS="string", _RC )
         endif
! Now print out the last time solar alarm rang
        endif

! Get the zenith angle at the center of the time between the last solar call and the next one
        call MAPL_SunGetInsolation(LONS, LATS,      &
            ORBIT, ZTH, SLR, &
            INTV = TINT,     &
            currTime=BEFORE+DELT,  &
            _RC )

    ZTH = max(0.0,ZTH)


! Albedo over Sea-Ice. With LANL CICE, it is based on current ice states,
! also compute shortwave radiation passing thru bottom of ice and skin layer bottom (DxxxTHRU; xx=RUVR, FUVR, ...)
!------------------------------------------------------------------------------------------------------------------

    call CICE_ALBSEAICE (ICE,NUM_ICE_CATEGORIES,NUM_ICE_LAYERS,NUM_SNOW_LAYERS,NT,DO_POND,LATSO,LONSO,LATS,LONS,ZTH,FR8,TS,&
                           DRPAR,DFPAR,DRNIR,DFNIR,DRUVR,DFUVR,VSUVR,VSUVF,VOLICE,VOLSNO,APONDN,HPONDN,                      &
                           ISWABS,FSWSFC, FSWINT,FSWTHRU,SSWABS,ALBIN,ALBSN,ALBPND,ALBVRN,ALBVFN,ALBNRN,ALBNFN,              &
                           DRUVRTHRU,DFUVRTHRU,DRPARTHRU,DFPARTHRU,_RC)
!!! Make this call so that during the predictor and corrector we use these albedos to send to radiation
     if(dual_ocean) then
      call ALBSEAICEM2 (ALBVRI,ALBVFI,ALBNRI,ALBNFI,ZTH,LATS,CURRENT_TIME)  ! GEOS albedo over sea ice
      do N=1, NUM_ICE_CATEGORIES
       ALBVRN(:,N) = ALBVRI(:)
       ALBVFN(:,N) = ALBVFI(:)
       ALBNRN(:,N) = ALBNRI(:)
       ALBNFN(:,N) = ALBNFI(:)
      enddo
     else
     endif


    if(associated(PENUVR)) PENUVR = 0.0
    if(associated(PENUVF)) PENUVF = 0.0
    if(associated(PENPAR)) PENPAR = 0.0
    if(associated(PENPAF)) PENPAF = 0.0

    do N=1, NUM_ICE_CATEGORIES
       if(associated(PENUVR)) PENUVR  = PENUVR + FR8(:,N) * DRUVRTHRU(:,N)
       if(associated(PENUVF)) PENUVF  = PENUVF + FR8(:,N) * DFUVRTHRU(:,N)
       if(associated(PENPAR)) PENPAR  = PENPAR + FR8(:,N) * DRPARTHRU(:,N)
       if(associated(PENPAF)) PENPAF  = PENPAF + FR8(:,N) * DFPARTHRU(:,N)
    enddo
    if(associated(PENUVR  )) call Normalize(PENUVR,  FRCICE)
    if(associated(PENUVF  )) call Normalize(PENUVF,  FRCICE)
    if(associated(PENPAR  )) call Normalize(PENPAR,  FRCICE)
    if(associated(PENPAF  )) call Normalize(PENPAF,  FRCICE)

    if(associated(ALBINe  )) then
        ALBINe = ALBIN
        where(FR8 < puny)
           ALBINe = MAPL_UNDEF
        endwhere
        do N=1, NUM_ICE_CATEGORIES
           where(ZTH < puny)
              ALBINe(:,N) = MAPL_UNDEF
           endwhere
        enddo
    endif

    if(associated(ALBSNe  )) then
        ALBSNe = ALBSN
        do N=1, NUM_ICE_CATEGORIES
           do K=1,NT
              if(FR8(K,N) <= puny) then
                 ALBSNe(K,N) = MAPL_UNDEF
              elseif(VOLSNO(K,N)/FR8(K,N) <= puny) then
                 ALBSNe(K,N) = MAPL_UNDEF
              endif
           enddo
        enddo
        do N=1, NUM_ICE_CATEGORIES
           where(ZTH < puny)
              ALBSNe(:,N) = MAPL_UNDEF
           endwhere
        enddo
    endif

    call MAPL_TimerOff(MAPL,    "-Albedo")

! Cycle through sub-tiles doing water and energy budget
!------------------------------------------------------

    if(associated(EVAPOUT)) EVAPOUT = 0.0
    if(associated(SHOUT  )) SHOUT   = 0.0
    if(associated(HLATN  )) HLATN   = 0.0
    if(associated(DELTS  )) DELTS   = 0.0
    if(associated(DELQS  )) DELQS   = 0.0
    if(associated(TST    )) TST     = 0.0
    if(associated(QST    )) QST     = 0.0
    if(associated(HLWUP  )) HLWUP   = 0.0
    if(associated(HLWUPe )) HLWUPe  = 0.0
    if(associated(LWNDSRF)) LWNDSRF = 0.0

    if(associated(SUBLIM )) SUBLIM  = 0.0
    if(associated(HLATICE)) HLATICE = 0.0
    if(associated(FSURFe )) FSURFe  = 0.0
    if(associated(FSURFICE)) FSURFICE = 0.0
    if(associated(SHICE  )) SHICE   = 0.0
    if(associated(MELTTL )) MELTTL  = 0.0
    if(associated(MELTBL )) MELTBL  = 0.0
    if(associated(MELTSL )) MELTSL  = 0.0
    if(associated(CONGELO)) CONGELO = 0.0
    if(associated(SNOICEO)) SNOICEO = 0.0
    if(associated(FBOTL  )) FBOTL   = 0.0
    if(associated(SBLXOUT)) SBLXOUT = 0.0

    if(associated(EVAP_C5))        EVAP_C5         =  0.0
    if(associated(GRCONGEL_C5))    GRCONGEL_C5     =  0.0
    if(associated(GRLATERAL_C5))   GRLATERAL_C5    =  0.0
    if(associated(SNOMELT_C5))     SNOMELT_C5      =  0.0
    if(associated(TMELT_C5))       TMELT_C5        =  0.0
    if(associated(BMELT_C5))       BMELT_C5        =  0.0
    if(associated(IALB_C5))        IALB_C5         =  0.0
    if(associated(RSUSSI_C5))      RSUSSI_C5       =  0.0
    if(associated(FSITHERM_CMIP5)) FSITHERM_CMIP5  =  0.0

    if(associated(HICEUNT )) HICEUNT    = sum(VOLICE, dim=2)
    if(associated(USTARI))   USTARI     = sqrt(sqrt(TAUXBOT**2+TAUYBOT**2)/MAPL_RHO_SEAWATER)

    if(associated(PR_C5)) then
          PR_C5 = FRCICE * (PLS + PCU)
          where(FRCICE == 0.0)
             PR_C5 = 0.0
          endwhere
    endif
    if(associated(PRSN_C5)) then
          PRSN_C5 = FRCICE * SNO
          where(FRCICE == 0.0)
             PRSN_C5 = 0.0
          endwhere
    endif

! Atmospheric surface stresses
!-----------------------------
    !*** already computed in surface
    !UUA = (TAUX/CMATM + UHATM)
    !VVA = (TAUY/CMATM + VHATM)


! Stress over ice
!----------------

    where(FRCICE > puny)
        TXI = sum(CM(:,ICE:)*FR8(:,ICE:), dim=2)/FRCICE*(UUA - UI)
        TYI = sum(CM(:,ICE:)*FR8(:,ICE:), dim=2)/FRCICE*(VVA - VI)
    elsewhere
        TXI = 0.0
        TYI = 0.0
    endwhere

    if(associated(TAUXI)) TAUXI = TXI
    if(associated(TAUYI)) TAUYI = TYI




!xxxxxxxxxxxxxxxxxxxxxxxxxxLANL CICE: 2 step update procedure-- STARTS xxxxxxxxxxxxxxxxxxxxxxxxx
    call MAPL_TimerOn(MAPL,   "-Thermo1")

! 1st Step of LANL CICE Thermodynamics
! ------------------------------------


     FR8TMP = FR8
     categories_th1_: do N=ICE, NUM_SUBTILES   ! Loop over ice catgories.

          NSUB = N

          CFT = (CH(:,N)/CTATM)
          CFQ = (CQ(:,N)/CQATM)
          EVP = CFQ*(EVAP + DEV*(QS(:,N)-QHATM))
          SHF = CFT*(SH   + DSH*(TS(:,N)-THATM))
          SHD = CFT*DSH
          EVD = CFQ*DEV*GEOS_DQSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.)
          LHF = EVP * MAPL_ALHS

          if(associated(SHICEN)) then
               SHICEN(:,N)  = SHF
               where(FR8(:,N) <= puny)
                  SHICEN(:,N) = MAPL_UNDEF
               endwhere
          endif

          call CICE_THERMO1(N,NSUB,NT,ICE,LATS,LONS,LATSO,LONSO,DT,TF,FR8TMP,TS,         &
                           ERGICE,ERGSNO,TAUXBOT,TAUYBOT,TBOT,ISWABS,SSWABS,             &
                           DO_POND,FBOT,RSIDE,PCU,PLS,FSURF,                             &
                           FSWTHRU,FCOND,FCONDBOT,EVP,FRESHN,FSALTN,FHOCNN,              &
                           MELTT,MELTS,MELTB,CONGEL,SNOICE,VOLICE,VOLSNO,SHF,LHF,        &
                           VOLPOND,APONDN,HPONDN,TAUAGE,TRACERS,ALW,BLW,    &
                           FSWSFC,FSWINT,FSWABS,LWDNSRF,EVD,SHD,SNO,SBLX,_RC)

!         Some aggregation of fluxes to the Ocean has to be done now, before using in step2

          FRESHL   = FRESHL   + FRESHN *FR8(:,N)
          FSALTL   = FSALTL   + FSALTN *FR8(:,N)
          FHOCNL   = FHOCNL   + FHOCNN *FR8(:,N)

          NEWICEERG = NEWICEERG + (FCONDBOT(:,NSUB) - FBOT) * FR8(:,N)

!         Update surface temperature and moisture
!         ----------------------------------------

          if(associated(SHICE)  ) SHICE   = SHICE   + SHF        *FR8(:,N)     ! aggregate ice surface fluxes into atm
          if(associated(HLATICE)) HLATICE = HLATICE + LHF        *FR8(:,N)
          if(associated(FSURFe )) FSURFe  = FSURFe  + FSURF      *FR8(:,N)
          if(associated(FSURFICE)) FSURFICE  = FSURFICE  + FSURF *FR8(:,N)
          if(associated(TSURFN )) then
               TSURFN(:,N) = TS(:,N)
               where(FR8(:,N) <= puny)
                 TSURFN(:,N) = MAPL_UNDEF
               endwhere
          endif
          if(associated(FSURFN )) then
               FSURFN(:,N) = FSURF
               where(FR8(:,N) <= puny)
                 FSURFN(:,N) = MAPL_UNDEF
               endwhere
          endif
          if(associated(HLWUPN)) then
               HLWUPN(:,N)  = ALW+BLW*TS(:,N)
               where(FR8(:,N) <= puny)
                  HLWUPN(:,N) = MAPL_UNDEF
               endwhere
          endif
          if(associated(LWNDSRFN)) then
               LWNDSRFN(:,N)  = LWDNSRF-(ALW+BLW*TS(:,N))
               where(FR8(:,N) <= puny)
                  LWNDSRFN(:,N) = MAPL_UNDEF
               endwhere
          endif
          if(associated(EVAPOUT)) EVAPOUT = EVAPOUT + EVP        *FR8(:,N)
          if(associated(SUBLIM )) SUBLIM  = SUBLIM  + EVP        *FR8(:,N)
          if(associated(SHOUT  )) SHOUT   = SHOUT   + SHF        *FR8(:,N)
          if(associated(HLATN  )) HLATN   = HLATN   + LHF        *FR8(:,N)
          if(associated(LWNDSRF)) LWNDSRF = LWNDSRF + (LWDNSRF-ALW-BLW*TS(:,N))*FR8(:,N)
          if(associated(HLWUP  )) HLWUP   = HLWUP   +         (ALW+BLW*TS(:,N))*FR8(:,N)
          if(associated(HLWUPe )) HLWUPe  = HLWUPe  +         (ALW+BLW*TS(:,N))*FR8(:,N)

          if(associated(MELTTL )) MELTTL   = MELTTL   + MELTT   *FR8(:,N) / DT ! m per step -> m s-1
          if(associated(MELTBL )) MELTBL   = MELTBL   + MELTB   *FR8(:,N) / DT ! m per step -> m s-1
          if(associated(MELTSL )) MELTSL   = MELTSL   + MELTS   *FR8(:,N) / DT ! m per step -> m s-1
          if(associated(CONGELO)) CONGELO  = CONGELO  + CONGEL  *FR8(:,N) / DT ! m per step -> m s-1

          if(associated(SBLXOUT)) SBLXOUT  = SBLXOUT  + SBLX    *FR8(:,N) / DT

          if(associated(EVAP_C5))     EVAP_C5      = EVAP_C5     +                      EVP * FR8(:,N)
          if(associated(GRCONGEL_C5)) GRCONGEL_C5  = GRCONGEL_C5 + MAPL_RHO_SEAICE * CONGEL * FR8(:,N) / DT ! kg m-2 s-1
          if(associated(SNOMELT_C5))  SNOMELT_C5   = SNOMELT_C5  + MAPL_RHO_SNOW   * MELTS  * FR8(:,N) / DT ! kg m-2 s-1
          if(associated(TMELT_C5))    TMELT_C5     = TMELT_C5    + MAPL_RHO_SEAICE * MELTT  * FR8(:,N) / DT ! kg m-2 s-1
          if(associated(BMELT_C5))    BMELT_C5     = BMELT_C5    + MAPL_RHO_SEAICE * MELTB  * FR8(:,N) / DT ! kg m-2 s-1

!         Aggregate ts and qs change over ice categories
          DTS     = TS(:,N) - TS_OLD(:,N)
          DQS     = GEOS_QSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.) - QS(:,N)
          QS(:,N) = QS(:,N) + DQS

          if(associated(DELTS  )) DELTS   = DELTS   + DTS*CFT*FR8(:,N)
          if(associated(DELQS  )) DELQS   = DELQS   + DQS*CFQ*FR8(:,N)
     end do categories_th1_

     if(associated(FSWSFCN  )) then
         FSWSFCN  = FSWSFC
         do N=1, NUM_ICE_CATEGORIES
            where(FR8(:,N) <= puny)
               FSWSFCN(:,N) = MAPL_UNDEF
            endwhere
         enddo
     endif
     if(associated(FCONDBOTN))  then
         FCONDBOTN      = FCONDBOT
         do N=1, NUM_ICE_CATEGORIES
            where(FR8(:,N) <= puny)
               FCONDBOTN(:,N) = MAPL_UNDEF
            endwhere
         enddo
     endif
     if(associated(FCONDTOPN))  then
         FCONDTOPN      = FCOND
         do N=1, NUM_ICE_CATEGORIES
            where(FR8(:,N) <= puny)
               FCONDTOPN(:,N) = MAPL_UNDEF
            endwhere
         enddo
     endif

     if(associated(DELTS  )) call Normalize(DELTS,  FRCICE)
     if(associated(DELQS  )) call Normalize(DELQS,  FRCICE)
     if(associated(EVAPOUT)) call Normalize(EVAPOUT,FRCICE)
     if(associated(SHOUT  )) call Normalize(SHOUT,  FRCICE)
     if(associated(HLATN  )) call Normalize(HLATN,  FRCICE)
     if(associated(HLWUP  )) call Normalize(HLWUP,  FRCICE)
     if(associated(LWNDSRF)) call Normalize(LWNDSRF,FRCICE)
     if(associated(FSURFe )) call Normalize(FSURFe, FRCICE)

     if(associated(TST    )) then
          TST = sum(TS*FR8,dim=2)
          call Normalize(TST,    FRCICE)
     endif
     if(associated(QST    )) then
          QST = sum(QS*FR8,dim=2)
          call Normalize(QST,    FRCICE)
     endif

     if(associated(FSWTRUO)) FSWTRUO     = sum(FR8(:,ICE:)*FSWTHRU,dim=2)
     if(associated(FBOTL  )) FBOTL       = FBOT
     if(associated(SNOONICE )) SNOONICE  = FRCICE*SNO

     if(associated(LWNDICE)) then
          where( FRCICE>puny )
             LWNDICE = LWDNSRF - ALW - BLW*(sum(TS(:,ICE:)*FR8(:,ICE:), dim=2) / FRCICE)
          elsewhere
             LWNDICE = MAPL_UNDEF
          end where
     endif

     if(associated(SHICE)) then
          where( FRCICE>puny )
             SHICE = SHICE / FRCICE
          elsewhere
             SHICE = MAPL_UNDEF
          endwhere
     endif

     if(associated(HLATICE)) then
          where( FRCICE>puny )
             HLATICE = HLATICE / FRCICE
          elsewhere
             HLATICE = MAPL_UNDEF
          endwhere
     endif

     if(associated(HLWUPe)) then
          where( FRCICE>puny )
             HLWUPe = HLWUPe / FRCICE
          elsewhere
             HLWUPe = MAPL_UNDEF
          endwhere
     endif

     if(associated(GHTSKIN)) then
          where(FRCICE > puny)
             !*** multiply by -1 to be +ve upward
             GHTSKIN = -1.0*sum(FCOND*FR8(:,ICE:),dim=2)/FRCICE
          elsewhere
             GHTSKIN = MAPL_UNDEF
          endwhere
     endif

     if(associated(FCONDTOP)) then
          where(FRCICE > puny)
             FCONDTOP = sum(FCOND*FR8(:,ICE:),dim=2)/FRCICE
          elsewhere
             FCONDTOP = MAPL_UNDEF
          endwhere
     endif

     if(associated(FSURFICE)) then
          where(FRCICE > puny)
             FSURFICE = FSURFICE / FRCICE
          elsewhere
             FSURFICE = MAPL_UNDEF
          endwhere
     endif

     if(associated(FCONDB)) then
          where(FRCICE > puny)
             FCONDB = sum(FCONDBOT*FR8(:,ICE:),dim=2)/FRCICE
          elsewhere
             FCONDB = MAPL_UNDEF
          endwhere
     end if

     if(associated(SNOMELT_C5)) then
        where(FRCICE == 0.0)
             SNOMELT_C5 = 0.0
        endwhere
     endif

     if(associated(TMELT_C5)) then
          where(FRCICE == 0.0)
             TMELT_C5 = 0.0
          endwhere
     endif

     if(associated(BMELT_C5)) then
          where(FRCICE == 0.0)
             BMELT_C5 = 0.0
          endwhere
     endif

     if(associated(HFSIFRAZIL_C5)) then
          HFSIFRAZIL_C5 = FRZMLT
          where(FRZMLT < 0.0)
             HFSIFRAZIL_C5 = 0.0
          endwhere
     endif

     ALBVR = sum(ALBVRN(:,:)*FR8(:,ICE:),dim=2)
     ALBVF = sum(ALBVFN(:,:)*FR8(:,ICE:),dim=2)
     ALBNR = sum(ALBNRN(:,:)*FR8(:,ICE:),dim=2)
     ALBNF = sum(ALBNFN(:,:)*FR8(:,ICE:),dim=2)

     EMISS = EMSICE*FRCICE

     where( FRCICE>puny )
          EMISS = EMISS / FRCICE
          ALBVR = ALBVR / FRCICE
          ALBVF = ALBVF / FRCICE
          ALBNR = ALBNR / FRCICE
          ALBNF = ALBNF / FRCICE
     endwhere

     if(associated(SWNDICE)) then
        where( FRCICE>0.0 )
             SWNDICE = (1.-ALBVR)*VSUVR + (1.-ALBVF)*VSUVF + &
                       (1.-ALBNR)*DRNIR + (1.-ALBNF)*DFNIR
        elsewhere
             SWNDICE = MAPL_UNDEF
        end where
     end if

    if(associated(SWNDSRF)) then
       SWNDSRF = &
           (1.-ALBVR)*VSUVR + (1.-ALBVF)*VSUVF + &
           (1.-ALBNR)*DRNIR + (1.-ALBNF)*DFNIR
    endif

    if(associated(SIALB)) then
         where(FRCICE > puny)
             SIALB = awtvdr*ALBVR+awtvdf*ALBVF+awtidr*ALBNR+awtidf*ALBNF
         endwhere
         where(FRCICE <= puny)
              SIALB = MAPL_UNDEF
         endwhere
         where(ZTH < 1.e-6)
            SIALB = MAPL_UNDEF
         endwhere
    endif

    ! update internal FR8 here with potentially changed fraction
    FR8 = FR8TMP
    FRCICE = sum(FR8(:,ICE:), dim=2)

    call MAPL_TimerOff(MAPL,  "-Thermo1")

! 2nd Step of LANL CICE Thermodynamics:
! loops over ice categories within the  subroutines,
! redistributing ice and water mass due to freezing and melting
! -------------------------------------------------------------
#if 0
    do k=1,NT
       if((abs(LATS(k)*rad_to_deg-LATSO) < 1.e-3) .and. (abs(LONS(k)*rad_to_deg-LONSO) < 1.e-3)) then
          print*, 'after THERMO1 before THERMO2_STEP1'
          do i=1,NUM_ICE_CATEGORIES
             print*, i, FR8(K,i), VOLICE(K,i)
             print*, i, VOLSNO(K,i), ERGSNO(K,1,i)
             if(VOLSNO(K,i) > 0.0_8) then
                print*, i, ERGSNO(K,1,i)/VOLSNO(K,i), (Lfresh+&
                  ERGSNO(K,1,i)/VOLSNO(K,i)/rhos)/cp_ice
             endif
          enddo
       endif
    enddo
#endif

    call MAPL_TimerOn(MAPL,    "-Thermo2")

    call CICE_THERMO2_STEP1 (NT,ICE,LATS,LONS,LATSO,LONSO,DT,TF,FR8,TS,    &
                             VOLICE,VOLSNO,VOLPOND,ERGICE,ERGSNO,                    &
                             AICENINIT,VICENINIT,TRCRTYPE,FRCICE,FRZMLT,FRAZLN,      &
                             FRESHL,FSALTL,FHOCNL,RSIDE,MELTLN,VOLICE_DELTA,         &
                             TRACERS,TAUAGE,SNOICE,SW,_RC)

    FRCICE       = sum(FR8(:,ICE:), dim=2)

    if(associated(FSITHERM_CMIP5)) FSITHERM_CMIP5 = FRESHL

#if 0
    do k=1,NT
       if((abs(LATS(k)*rad_to_deg-LATSO) < 1.e-3) .and. (abs(LONS(k)*rad_to_deg-LONSO) < 1.e-3)) then
          print*, 'after THERMO2_STEP1 before THERMO2_STEP2'
          do i=1,NUM_ICE_CATEGORIES
             print*, i, FR8(K,i), VOLICE(K,i)
             print*, i, VOLSNO(K,i), ERGSNO(K,1,i)
             if(VOLSNO(K,i) > 0.0_8) then
                print*, i, ERGSNO(K,1,i)/VOLSNO(K,i), (Lfresh+&
                  ERGSNO(K,1,i)/VOLSNO(K,i)/rhos)/cp_ice
             endif
          enddo
       endif
    enddo
#endif

    !*** artificially do a lateral melt step over those frozen lake tiles if the ice gets too thick
    call CICE_THERMO2_STEP2 (NT,ICE,LATS,LONS,LATSO,LONSO,DT,FR8,TS,          &
                             VOLICE,VOLSNO,VOLPOND,ERGICE,ERGSNO,                     &
                             TRCRTYPE,FRCICE,SLMASK,TRACERS,TAUAGE,_RC)

    ! aggregate ice concentration after step2
    ! These are the final area fractions that are in the internal state
    FRCICE       = sum(FR8(:,ICE:), dim=2)

    do k=1,NT
        do N=1,NUM_ICE_CATEGORIES
           if(FR8(K,N) < puny) TI(K,N) = MAPL_TICE+Tocnfrz
        enddo
    enddo

    if(associated(FRACINEW)) FRACINEW = REAL(FRCICE, KIND=MAPL_R4)

    if(associated(GRLATERAL_C5)) then
          GRLATERAL_C5 = MAPL_RHO_SEAICE*sum(VOLICE_DELTA, dim=2)/DT
    endif

    if(associated(SNOTOICE_C5)) then
          SNOTOICE_C5  =  MAPL_RHO_SEAICE * SNOICE / DT ! kg m-2 s-1
          where( sum(VOLSNO, dim=2) == 0.0)
             SNOTOICE_C5 = 0.0
          end where
    endif

    if(associated(SFDSI_C5)) then
          SFDSI_C5  =  FSALTL ! kg m-2 s-1
          where(FRCICE == 0.0)
             SFDSI_C5 = 0.0
          end where
    endif

    if(associated(DAIDTT)) then
          DAIDTT = (FRCICE - FR_OLD) / DT * 8640000
    endif

    if(associated(TINZ   )) then
        TINZ = MAPL_UNDEF
        do K=1, NT
           call diagnose_internal_ice_temp(VOLICE(K,:), ERGICE(K,:,:), TINZ(K,:))
        enddo
    endif


    if(associated(NIERG)) NIERG = NEWICEERG

    if(associated(SNOICEO)) SNOICEO = SNOICE / DT ! m per step -> m s-1
    if(associated(FRESH  )) FRESH   = FRESHL
    if(associated(FSALT  )) FSALT   = FSALTL
    if(associated(FHOCN  )) FHOCN   = FHOCNL
    !if(associated(PICE   )) PICE    = MAPL_GRAV*(sum(VOLICE,dim=2)*MAPL_RHO_SEAICE + sum(VOLSNO,dim=2)*MAPL_RHO_SNOW)
    !fow now, pass zero ice pressure loading
    !fully coupled ice-ocean dynamics not ready yet!!
    if(associated(PICE   )) PICE    = 0.0

    if(associated(TSKINICE)) then
          ! to be consisten with CICE (unit in degC)
          TSKINICE = sum((TS(:,ICE:)-TFfresh)*FR8(:,ICE:),dim=2)
          where(FRCICE > puny)
             TSKINICE = TSKINICE / FRCICE + MAPL_TICE
          elsewhere
             TSKINICE = MAPL_UNDEF
          end where
    endif

    if(associated(IAGE)) then
          ! here ice age is treated as an ice area tracer
          IAGE = sum(TAUAGE(:,ICE:)*FR8(:,ICE:),dim=2) * iage_converter
          where(FRCICE > puny)
             IAGE = IAGE / FRCICE
          elsewhere
             IAGE = MAPL_UNDEF
          end where
    endif

    ! the mean ice/snow thickness is computed as: sum_n_over_ice_categories(FR(n)*H(n)) which is simply
    ! sum_n_over_ice_categories(VOL(n))

    if(associated(HICE  )) HICE    =  sum(VOLICE(:,:),dim=2)
    if(associated(HSNO  )) HSNO    =  sum(VOLSNO(:,:),dim=2)
    if(associated(MELTL )) MELTL   =  MELTLN / DT                   ! m per step -> m s-1
    if(associated(FRAZIL)) FRAZIL  =  FRAZLN / DT                   ! m per step -> m s-1


#if 0
    do k=1,NT
       if((abs(LATS(k)*rad_to_deg-LATSO) < 1.e-3) .and. (abs(LONS(k)*rad_to_deg-LONSO) < 1.e-3)) then
          print*, 'end of cice thermo'
          do i=1,NUM_ICE_CATEGORIES
             print*, i, FR8(K,i), VOLICE(K,i)
             print*, i, VOLSNO(K,i), ERGSNO(K,1,i)
             if(VOLSNO(K,i) > 0.0_8) then
                print*, i, ERGSNO(K,1,i)/VOLSNO(K,i), (Lfresh+&
                  ERGSNO(K,1,i)/VOLSNO(K,i)/rhos)/cp_ice
             endif
          enddo
       endif
    enddo

    allocate(vice1(NT),                                   _STAT)
    vice1 =  sum(VOLICE,dim=2)
    TOTALAREAN1 = sum(vice1*AREA*(MAPL_RADIUS**2),mask=SLMASK<0.5 .and. LATS>0.0)
    TOTALAREAS1 = sum(vice1*AREA*(MAPL_RADIUS**2),mask=SLMASK<0.5 .and. LATS<0.0)
    deallocate(vice1)

    call ESMF_VMBarrier(VMG, _RC)
    call MAPL_CommsAllReduceSum(VMG, TOTALAREAN1, ALLTOTALAREAN1, 1, _RC)
    call MAPL_CommsAllReduceSum(VMG, TOTALAREAS1, ALLTOTALAREAS1, 1, _RC)

    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run2 North ice1  = ', &
                                 ALLTOTALAREAN1
    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run2 South ice1  = ', &
                                 ALLTOTALAREAS1
    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run2 North ice1-ice0  = ', &
                                 ALLTOTALAREAN1-ALLTOTALAREAN
    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run2 South ice1-ice0  = ', &
                                 ALLTOTALAREAS1-ALLTOTALAREAS
    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run2 North ice1-ice0 (%)  = ', &
                                 100*abs(ALLTOTALAREAN1-ALLTOTALAREAN)/ALLTOTALAREAN
    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run2 South ice1-ice0 (%)  = ', &
                                 100*abs(ALLTOTALAREAS1-ALLTOTALAREAS)/ALLTOTALAREAS
    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run2 total ice1-ice0  = ', &
                                 ALLTOTALAREAN1-ALLTOTALAREAN+ALLTOTALAREAS1-ALLTOTALAREAS
    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' Run2 total ice1-ice0 (%)  = ', &
                                 100*abs(ALLTOTALAREAN1-ALLTOTALAREAN+ALLTOTALAREAS1-ALLTOTALAREAS)/ &
                                     (ALLTOTALAREAN+ALLTOTALAREAS)
#endif

    if(associated(DVIDTT))  then
          DVIDTT  = (sum(VOLICE,dim=2) - VOLICE_OLD) / DT * 8640000
       end if

       if(associated(GRFRAZIL_C5)) then
          GRFRAZIL_C5  =  MAPL_RHO_SEAICE * FRAZLN / DT ! kg m-2 s-1
          where(FRCICE == 0.0)
             GRFRAZIL_C5 = 0.0
          end where
       end if

    call MAPL_TimerOff(MAPL,   "-Thermo2")

!xxxxxxxxxxxxxxxxxxxxxxxxxxLANL CICE: 2 step update procedure-- ENDS xxxxxxxxxxxxxxxxxxxxxxxxxxx

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    call MAPL_TimerOn(MAPL,    "-Albedo")

    if(solalarmison) then
       call MAPL_SunGetInsolation(LONS, LATS,      &
            ORBIT, ZTH, SLR,                       &
            INTV = TINT,                           &
            currTime=CURRENT_TIME+DELT,            &
            _RC )

       ZTH = max(0.0,ZTH)

       call CICE_ALBSEAICE (ICE,NUM_ICE_CATEGORIES,NUM_ICE_LAYERS,NUM_SNOW_LAYERS,NT,DO_POND,LATSO,LONSO,LATS,LONS,ZTH,FR8,TS,&
                            DRPAR,DFPAR,DRNIR,DFNIR,DRUVR,DFUVR,VSUVR,VSUVF,VOLICE,VOLSNO,APONDN,HPONDN,                      &
                            ISWABS,FSWSFC, FSWINT,FSWTHRU,SSWABS,ALBIN,ALBSN,ALBPND,ALBVRN,ALBVFN,ALBNRN,ALBNFN,              &
                            DRUVRTHRU,DFUVRTHRU,DRPARTHRU,DFPARTHRU,_RC)

       do N=1,NUM_ICE_CATEGORIES
             do K=1,NT
                if(FR8(K,N) > puny) then
                   if(associated(IALB_C5))   IALB_C5(K)   = IALB_C5(K) + FR8(K,N) * ALBIN(K,N)
                end if
             end do
       end do

          ! report "missing" if there is no sunlight or free of ice
       if(associated(IALB_C5)) then
             where(FRCICE > puny)
                IALB_C5 =  IALB_C5 / FRCICE
             endwhere
             where(FRCICE <= puny)
                IALB_C5 = MAPL_UNDEF
             endwhere
             where(ZTH < 1.e-6)
                IALB_C5 = MAPL_UNDEF
             endwhere
       endif

       ALBVR = sum(ALBVRN(:,:)*FR8(:,ICE:),dim=2)
       ALBVF = sum(ALBVFN(:,:)*FR8(:,ICE:),dim=2)
       ALBNR = sum(ALBNRN(:,:)*FR8(:,ICE:),dim=2)
       ALBNF = sum(ALBNFN(:,:)*FR8(:,ICE:),dim=2)

       where(FRCICE > puny)
         ALBVR = ALBVR / FRCICE
         ALBVF = ALBVF / FRCICE
         ALBNR = ALBNR / FRCICE
         ALBNF = ALBNF / FRCICE
       endwhere

       if(associated(RSUSSI_C5)) then
             RSUSSI_C5 =  ALBVR*VSUVR + ALBVF*VSUVF  + ALBNR*DRNIR + ALBNF*DFNIR
       endif

    endif

    call MAPL_TimerOff(MAPL,    "-Albedo")

    call ESMF_VMBarrier(VM, _RC)
    call MAPL_TimerOn(MAPL,    "-Out_ReDist")
    if(loadBalance) then
#include "BufferUnpacking.h"
       deallocate(BUFIMP,BUFINT,BUFINT8,BUFEXP,_STAT)

       call MAPL_BalanceDestroy(Handle=CICECOREBalanceHandle, _RC)
    endif
    call MAPL_TimerOff(MAPL,    "-Out_ReDist")

    deallocate(FSWABS)
    deallocate(ALBVRI)
    deallocate(ALBVFI)
    deallocate( ALBNRI)
    deallocate(ALBNFI)
    deallocate(SHF)
    deallocate(EVP)
    deallocate(SHD)
    deallocate(EVD)
    deallocate(CFQ)
    deallocate(CFT)
    deallocate(TXI)
    deallocate(TYI)
    deallocate(DQS)
    deallocate(DTS)
    deallocate(DTX)
    deallocate(DTY)
    deallocate(SWN)
    deallocate(PEN)
    deallocate(LHF)
    deallocate(ZTH)
    deallocate(SLR)
    deallocate(VSUVR)
    deallocate(VSUVF)
    deallocate(TRCRTYPE)
    deallocate(TRACERS)
    deallocate(TF)
    deallocate(MELTLN)
    deallocate(FRAZLN)
    deallocate(FRESHN)
    deallocate(FRESHL)
    deallocate(FSALTN)
    deallocate(FSALTL)
    deallocate(FHOCNN)
    deallocate(FHOCNL)
    deallocate(RSIDE)
    deallocate(FSWTHRU)
    deallocate(FCOND)
    deallocate(FCONDBOT)
    deallocate(TBOT)
    deallocate(FBOT)
    deallocate(ALBVRN)
    deallocate(ALBNRN)
    deallocate(ALBVFN)
    deallocate(ALBNFN)
    deallocate(FSWSFC)
    deallocate(FSWINT)
    deallocate(ISWABS)
    deallocate(SSWABS)
    deallocate(MELTT)
    deallocate(MELTS)
    deallocate(MELTB)
    deallocate(CONGEL)
    deallocate(SNOICE)
    deallocate(TS_OLD)
    deallocate(ALBIN)
    deallocate(ALBSN)
    deallocate(ALBPND)
    deallocate(DRUVRTHRU)
    deallocate(DFUVRTHRU)
    deallocate(DRPARTHRU)
    deallocate(DFPARTHRU)
    deallocate(TOTALFLUX)
    deallocate(NEWICEERG)
    deallocate(SBLX)
    deallocate(FSURF)
    deallocate(AICENINIT)
    deallocate(VICENINIT)
    deallocate(FR8TMP)
    deallocate(FRCICE)
    deallocate(FR_OLD)
    deallocate(VOLSNO_OLD)
    deallocate(VOLICE_OLD)
    deallocate(VOLICE_DELTA)
    deallocate(TEMPVOLICE)

!  All done
!-----------

    RETURN_(ESMF_SUCCESS)

  end subroutine CICECORE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: CICE_PREP_THERMO - Initializes for CICE Thermodynamics

! !INTERFACE:

  subroutine CICE_PREP_THERMO(TF, TRCRTYPE,TRACERS,MELTLN,FRAZLN,FRESHN,FRESHL,FSALTN,FSALTL,FHOCNN,FHOCNL,RSIDE,  &
                              FSWTHRU,FCOND,FCONDBOT,TBOT,FBOT,ALBIN,ALBSN,ALBPND,ALBVRN,ALBNRN,ALBVFN,ALBNFN,FSWSFC,FSWINT,     &
                              ISWABS,SSWABS,FSWABS,MELTT,MELTS,MELTB,CONGEL,SNOICE,UW,VW,SLMASK,LATS,LONS,LATSO,LONSO,   &
                              FR ,FRCICE,SW,TAUAGE,ICE,NT,VOLPOND,DT,VOLICE,VOLSNO,ERGICE,ERGSNO,TS,VOLICE_DELTA,  &
                              NEWICEERG, SBLX, RC)
! not passing TFfresh


! !ARGUMENTS:
    integer, optional, intent(OUT) :: RC

    integer, intent(IN)   :: NT                 ! number of tiles
    integer, intent(IN)   :: ICE                ! integer number of 1st ice subtile(s)
    real,    intent(IN)   :: UW         (:)     ! u-current
    real,    intent(IN)   :: VW         (:)     ! v-current
    real,    intent(INOUT):: SLMASK     (:)     ! mask for cice- relating open water and sea-ice
    real,    intent(IN)   :: LATS       (:)     ! lat
    real,    intent(IN)   :: LONS       (:)     ! lon
    real,    intent(IN)   :: LATSO              ! trace CICE computations at this latitude
    real,    intent(IN)   :: LONSO              ! trace CICE computations at this longitude
    real,    intent(IN)   :: SW         (:)     ! Sea Water salinity

    real,    intent(OUT)  :: TF         (:)     !
    integer, intent(OUT)  :: TRCRTYPE   (:)     ! CICE ice tracer type
    real,    intent(OUT)  :: TRACERS    (:,:)   ! CICE ice tracers
    real,    intent(OUT)  :: MELTLN     (:)     ! ?
    real,    intent(OUT)  :: FRAZLN     (:)     ! ?
    real,    intent(OUT)  :: FRESHN     (:)     ! ?
    real,    intent(OUT)  :: FRESHL     (:)     ! ?
    real,    intent(OUT)  :: FSALTN     (:)     ! ?
    real,    intent(OUT)  :: FSALTL     (:)     ! ?
    real,    intent(OUT)  :: FHOCNN     (:)     ! ?
    real,    intent(OUT)  :: FHOCNL     (:)     ! ?
    real,    intent(OUT)  :: RSIDE      (:)     ! ?
    real,    intent(OUT)  :: FSWTHRU    (:,:)   ! SW_flux_thru_ice_to_ocean
    real,    intent(OUT)  :: FCOND      (:,:)   ! ?
    real,    intent(OUT)  :: FCONDBOT   (:,:)   ! ?
    real,    intent(OUT)  :: TBOT       (:)     ! ?
    real,    intent(OUT)  :: FBOT       (:)     ! ?
    real,    intent(OUT)  :: ALBIN      (:,:)   !
    real,    intent(OUT)  :: ALBSN      (:,:)   !
    real,    intent(OUT)  :: ALBPND     (:,:)   !
    real,    intent(OUT)  :: ALBVRN     (:,:)   ! Albedos:
    real,    intent(OUT)  :: ALBNRN     (:,:)   ! Albedos: near IR
    real,    intent(OUT)  :: ALBVFN     (:,:)   ! Albedos:
    real,    intent(OUT)  :: ALBNFN     (:,:)   ! Albedos: near IR
    real,    intent(OUT)  :: FSWSFC     (:,:)   ! ?
    real,    intent(OUT)  :: FSWINT     (:,:)   ! ?
    real,    intent(OUT)  :: ISWABS     (:,:,:) ! ?
    real,    intent(OUT)  :: SSWABS     (:,:,:) ! ?
    real,    intent(OUT)  :: FSWABS     (:)     ! of dimension(1)
    real,    intent(OUT)  :: MELTT      (:)     ! ?
    real,    intent(OUT)  :: MELTS      (:)     ! ?
    real,    intent(OUT)  :: MELTB      (:)     ! ?
    real,    intent(OUT)  :: CONGEL     (:)     ! ?
    real,    intent(OUT)  :: SNOICE     (:)     ! ?
    real,    intent(OUT)  :: NEWICEERG  (:)     ! ?
    real,    intent(OUT)  :: SBLX       (:)     ! ?
    real,    intent(INOUT):: TS         (:,:)   ! skin temperature

    real,    intent(IN)   :: DT
    real(kind=MAPL_R8)    :: DTDB               ! DT (time step) in R8 for CICE
    real(kind=MAPL_R8), intent(INOUT)  :: FRCICE       (:)
    real(kind=MAPL_R8), intent(INOUT)  :: FR           (:,:)
    real(kind=MAPL_R8), intent(INOUT)  :: VOLICE       (:,:)
    real(kind=MAPL_R8), intent(INOUT)  :: VOLSNO       (:,:)
    real(kind=MAPL_R8), intent(INOUT)  :: VOLPOND      (:,:)
    real(kind=MAPL_R8), intent(INOUT)  :: ERGICE       (:,:,:)
    real(kind=MAPL_R8), intent(INOUT)  :: ERGSNO       (:,:,:)
    real(kind=MAPL_R8), intent(INOUT)  :: VOLICE_DELTA (:,:)
    real,               intent(INOUT)  :: TAUAGE       (:,:)

!  !LOCAL VARIABLES

    character(len=ESMF_MAXSTR)            :: IAm
    integer                               :: STATUS

    integer                               :: PRES_ICE, K
    logical                               :: L_STOP
    integer                               :: IDUM, JDUM
    logical,            dimension(1)      :: OBSERVE
    real,               dimension(1)      :: LATSD, LONSD
    real(kind=MAPL_R8), dimension(1)      :: FRWATERDB, FHOCNLDB, FRESHLDB, FSALTLDB, FRCICEDB
    real                                  :: MIN_FREEZE_SALINITY

    real(kind=MAPL_R8), dimension(NUM_3D_ICE_TRACERS, NUM_ICE_CATEGORIES) :: TRACERSDB2

! !DESCRIPTION:

    IAm =  trim(COMP_NAME) // "CICECORE" // "CICE_PREP_THERMO"

    DTDB = REAL(DT, kind=MAPL_R8)       ! Convert DT precision: Real4 to Real8 for usage in CICE

!   PRESCRIBED ICE. 1:AMIP mode, 0: coupled mode
    call MAPL_GetResource ( MAPL, PRES_ICE, Label="PRESCRIBED_ICE:" , DEFAULT=1,    _RC)

    call MAPL_GetResource ( MAPL, MIN_FREEZE_SALINITY, Label="MIN_FREEZE_SALINITY:" , DEFAULT=0.0,    _RC)

    call  FreezingTemperature(TF, SW, MIN_FREEZE_SALINITY, PRES_ICE==1, kelvin=.false.)

    TRCRTYPE(nt_tsfc)  = 0  ! ice/snow surface temperature
    TRCRTYPE(nt_iage)  = 1  ! volume-weighted ice age
    TRCRTYPE(nt_volpn) = 0  ! melt pond volume

    TRACERS            = 0.0

    MELTLN             = 0.0
    FRAZLN             = 0.0
    FRESHN             = 0.0
    FRESHL             = 0.0
    FSALTN             = 0.0
    FSALTL             = 0.0
    FHOCNN             = 0.0
    FHOCNL             = 0.0
    RSIDE              = 0.0
    FSWTHRU            = 0.0
    FCOND              = 0.0
    FCONDBOT           = 0.0

    TBOT               = 0.0
    FBOT               = 0.0
    ALBIN              = 0.0
    ALBSN              = 0.0
    ALBPND             = 0.0
    ALBVRN             = 0.0
    ALBNRN             = 0.0
    ALBVFN             = 0.0
    ALBNFN             = 0.0
    FSWSFC             = 0.0
    FSWINT             = 0.0
    ISWABS             = 0.0
    SSWABS             = 0.0
    FSWABS             = 0.0
    MELTT              = 0.0
    MELTS              = 0.0
    MELTB              = 0.0
    CONGEL             = 0.0
    SNOICE             = 0.0
    NEWICEERG          = 0.0
    SBLX               = 0.0
    VOLICE_DELTA       = 0.0

! determine those tiles where there is no open ocean connection. See note for Atanas.
    where(abs(UW) >  0.0 .or. abs(VW) > 0.0)
        SLMASK = 0.0
    elsewhere
        SLMASK = 1.0
    endwhere

! do a cleanup here, in case transformation from tripolar to tile induces round-off errors
#if 0
    FRCICE = sum(FR (:,ICE:), dim=2)
    do k=1, NT

       call CICE_INQUIRE_TILE(LATS(K), LONS(K), LATSO, LONSO, OBSERVE, LATSD, LONSD)

       TRACERSDB2(nt_tsfc,:) = REAL(TS(K,ICE:)  - TFfresh, kind=MAPL_R8)
       TRACERSDB2(nt_iage,:) = REAL(TAUAGE(K,:),           kind=MAPL_R8)

       TRACERSDB2(nt_volpn,:)= VOLPOND(K,:)
       FRWATERDB             = 1.0_8 - FRCICE(k)

       FHOCNLDB              = REAL(FHOCNL(K),             kind=MAPL_R8)
       FRESHLDB              = REAL(FRESHL(K),             kind=MAPL_R8)
       FSALTLDB              = REAL(FSALTL(K),             kind=MAPL_R8)
       FRCICEDB              = REAL(FRCICE(K),             kind=MAPL_R8)

       call cleanup_itd (1,1,1,1,1,1,DTDB, &
            FR (K,ICE:),   TRACERSDB2,     &
            VOLICE(K,:),   VOLSNO(K,:),    &
            ERGICE(K,:,:), ERGSNO(K,:,:),  &
            FRWATERDB,     FRCICEDB,       &
            TRCRTYPE,      FRESHLDB,       &
            FSALTLDB,      FHOCNLDB,       &
            .true.,        L_STOP,         &
            IDUM,            JDUM,         &
            limit_aice_in=.true.)

       if(L_STOP) then
          print*, 'CICE_PREP_THERMO: Failing at LAT = ', LATSD, 'LON = ', LONSD
       endif

       _ASSERT(.not.L_STOP,'needs informative message')

       FRESHL(K)    = REAL(FRESHLDB(1),                     kind=MAPL_R4)
       FSALTL(K)    = REAL(FSALTLDB(1),                     kind=MAPL_R4)
       FHOCNL(K)    = REAL(FHOCNLDB(1),                     kind=MAPL_R4)

       TS(K,ICE:)   = REAL(TRACERSDB2(nt_tsfc,:) + TFfresh, kind=MAPL_R4)
       TAUAGE(K,:)  = REAL(TRACERSDB2(nt_iage,:),           kind=MAPL_R4)
       VOLPOND(K,:) = REAL(TRACERSDB2(nt_volpn,:),          kind=MAPL_R4)
    enddo
#endif
! freshwater, salt and heat flux accumulated previously is not counted
! because they are not **physical**
    FRESHL = 0.0
    FHOCNL = 0.0
    FSALTL = 0.0

! these lines for fr are efectively same in cmip & amip modes.
!*** FR(:,ICE:) returned from CICEDyna or Data Sea Ice
!*** update FRWATER accordingly
    FRCICE      = sum(FR (:,ICE:), dim=2)

    RETURN_(ESMF_SUCCESS)
  end subroutine CICE_PREP_THERMO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: CICE_THERMO1 - Computes 1st step of CICE Thermodynamics

! !INTERFACE:

  subroutine CICE_THERMO1 (N,NSUB,NT,ICE,LATS,LONS,LATSO,LONSO,DT,TF,FR,TS,              &
                           ERGICE,ERGSNO,TAUXBOT,TAUYBOT,TBOT,ISWABS,SSWABS,             &
                           DO_POND,FBOT,RSIDE,PCU,PLS,FSUR,                              &
                           FSWTHRU,FCOND,FCONDBOT,EVP,FRESHN,FSALTN,FHOCNN,              &
                           MELTT,MELTS,MELTB,CONGEL,SNOICE,VOLICE,VOLSNO,SHF,LHF,        &
                           VOLPOND,APONDN,HPONDN,TAUAGE,TRACERS,ALW,BLW,    &
                           FSWSFC,FSWINT,FSWABS,LWDNSRF,EVD,SHD,SNO,SBLX,RC)
! not passing TFfresh,saltwatercap,nt_tsfc,nt_iage,nt_volpn

! !ARGUMENTS:
    integer, optional, intent(OUT) :: RC

    integer, intent(IN)  ::  N                 ! number of subtile type (or ice category)
    integer, intent(IN)  ::  NSUB              ! number of tiles
    integer, intent(IN)  ::  NT                ! number of tiles
    integer, intent(IN)  ::  ICE               ! subtiles number assigned to surface type: "ICE"
    integer, intent(IN)  ::  DO_POND           ! doing computations for melt ponds explicitly?

    real,    intent(IN)  :: LATS       (:)     ! lat
    real,    intent(IN)  :: LONS       (:)     ! lon
    real,    intent(IN)  :: LATSO              ! trace CICE computations at this latitude
    real,    intent(IN)  :: LONSO              ! trace CICE computations at this longitude

    real,    intent(IN)  :: TF         (:)     ! sea Water freezing temperature in degrees C
    real,    intent(IN)  :: TAUXBOT    (:)     ! zonal      stress at base of sea ice
    real,    intent(IN)  :: TAUYBOT    (:)     ! meridional stress at base of sea ice
    real,    intent(IN)  :: ISWABS     (:,:,:) ! ?
    real,    intent(IN)  :: SSWABS     (:,:,:) ! ?
    real,    intent(IN)  :: PCU        (:)     ! liquid water convective scale
    real,    intent(IN)  :: PLS        (:)     ! liquid water large      scale
    real,    intent(IN)  :: ALW        (:)     ! linearization of \sigma T^4
    real,    intent(IN)  :: BLW        (:)     ! linearization of \sigma T^4
    real,    intent(IN)  :: FSWINT     (:,:)   ! ?
    real,    intent(IN)  :: FSWABS     (:)
    real,    intent(IN)  :: LWDNSRF    (:)     ! longwave at surface
    real,    intent(IN)  :: EVD        (:)     ! related to evap
    real,    intent(IN)  :: SHD        (:)     ! related to sensible heat
    real,    intent(IN)  :: SNO        (:)     ! ?

    real,    intent(INOUT)  :: FSWSFC  (:,:)   ! ?
    real,    intent(INOUT)  :: EVP     (:)     ! evaporation
    real,    intent(INOUT)  :: FBOT    (:)     ! ?
    real,    intent(INOUT)  :: RSIDE   (:)     ! ?
    real,    intent(INOUT)  :: FRESHN  (:)     ! ?
    real,    intent(INOUT)  :: FSALTN  (:)     ! ?
    real,    intent(INOUT)  :: FHOCNN  (:)     ! ?
    real,    intent(INOUT)  :: MELTT   (:)     ! ?
    real,    intent(INOUT)  :: MELTS   (:)     ! ?
    real,    intent(INOUT)  :: MELTB   (:)     ! ?
    real,    intent(INOUT)  :: CONGEL  (:)     ! ?
    real,    intent(INOUT)  :: SNOICE  (:)     ! ?
    real,    intent(INOUT)  :: SHF     (:)     ! sensible heat flux
    real,    intent(INOUT)  :: LHF     (:)     ! latent   heat flux
    real,    intent(INOUT)  :: TBOT    (:)     ! ?
    real,    intent(INOUT)  :: SBLX    (:)     ! ?
    real,    intent(INOUT)  :: TAUAGE  (:,:)   ! ?
    real,    intent(INOUT)  :: TRACERS (:,:)   ! ?
    real,    intent(INOUT)  :: TS      (:,:)   ! skin temperature
    real,    intent(OUT)    :: FSUR    (:)     ! ?

    real(kind=MAPL_R8),    intent(INOUT)  :: FR      (:,:)   ! fractions of water, ice types
    real(kind=MAPL_R8),    intent(INOUT)  :: VOLICE  (:,:)   ! volume of ice
    real(kind=MAPL_R8),    intent(INOUT)  :: VOLSNO  (:,:)   ! volume of snow
    real(kind=MAPL_R8),    intent(INOUT)  :: VOLPOND (:,:)   ! ?
    real(kind=MAPL_R8),    intent(INOUT)  :: APONDN  (:,:)   ! ?
    real(kind=MAPL_R8),    intent(INOUT)  :: HPONDN  (:,:)   ! ?
    real(kind=MAPL_R8),    intent(INOUT)  :: ERGICE  (:,:,:) ! ?
    real(kind=MAPL_R8),    intent(INOUT)  :: ERGSNO  (:,:,:) ! ?

    real,    intent(INOUT)  :: FSWTHRU (:,:)   ! shortwave thru water, how can it be modified?
    real,    intent(INOUT)  :: FCOND   (:,:)   ! ?
    real,    intent(INOUT)  :: FCONDBOT(:,:)   ! ?

    real,    intent(IN)     :: DT
    real(kind=MAPL_R8)      :: DTDB               ! DT (time step) in R8 for CICE

! !LOCAL VARIABLES
    character(len=ESMF_MAXSTR)     :: IAm
    integer                        :: STATUS

    integer   :: K, L
    INTEGER   :: one(1) = (/1/)
    logical   :: OBSERVE(1,1)                    ! could be (1,1) to match cice input
    integer   :: IDUM, JDUM
    logical   :: L_STOP
    real      :: LWUPSRF, DLHDT, DHLAT, DFSDT, LWUP0

    real,                  dimension(1,1)  :: LATSD, LONSD, FSURF, SHF0, LHF0
    character(len=ESMF_MAXSTR)             :: SHORTWAVE

    real(kind=MAPL_R8)                :: YDAYDB
    real(kind=MAPL_R8), dimension(1,1):: FRZMLTDB, TSCDB, TFDB, TAUXBOTDB, TAUYBOTDB, &
                                         TBOTDB, FBOTDB, RSIDEDB, FSURFDB, DLHDTDB,   &
                                         DFSDTDB, SHF0DB, LHF0DB, LWUP0DB, FSWSFCDB,  &
                                         FSWINTDB, LWDNSRFDB, SNODB, FSWABSDB,        &
                                         FSWTHRUDB, FCONDDB, FCONDBOTDB, EVPDB,       &
                                         FRESHNDB, FSALTNDB, FHOCNNDB,                &
                                         MELTTDB, MELTSDB, MELTBDB, CONGELDB,SNOICEDB,&
                                         DSHDB, BLWDB, LATSDB, LONSDB, MLT_ONSETDB,   &
                                         FRZ_ONSETDB, FRDB, VOLICEDB, VOLSNODB,       &
                                         SBLXDB,                                      &
                                         APONDNDB, HPONDNDB, RDUMDB, FRAINDB

    real(kind=MAPL_R8), dimension(NT)                  :: FRCICE
    real(kind=MAPL_R8), dimension(NUM_3D_ICE_TRACERS)  :: TRACERSDB
    real(kind=MAPL_R8), dimension(NUM_ICE_LAYERS)      :: ISWABSDB
    real(kind=MAPL_R8), dimension(NUM_SNOW_LAYERS)     :: SSWABSDB
    real(kind=MAPL_R4), dimension(NUM_ICE_LAYERS*NUM_ICE_CATEGORIES)   :: tint
    real(kind=MAPL_R8), dimension(NUM_ICE_LAYERS)      :: ERGICE_TMP
    real(kind=MAPL_R8), dimension(NUM_SNOW_LAYERS)     :: ERGSNO_TMP


!  !DESCRIPTION:
!        Compute update to TS, FR, SHF, LHF, ...
!          based on CICE Thermodynamics

    IAm =  trim(COMP_NAME) // "CICECORE" // "CICE_THERMO1"

    DTDB = REAL(DT, kind=MAPL_R8)       ! Convert DT precision: Real4 to Real8 for usage in CICE

    call MAPL_GetResource ( MAPL, SHORTWAVE,  Label="CICE_SHORTWAVE:" ,  DEFAULT="shortwave_ccsm" , _RC)

    FSUR = 0.0

! Loop over all tiles
!-----------------------

    TILES: do k=1, NT ! loop over all tiles

       !call CICE_INQUIRE_TILE(LATS(K), LONS(K), LATSO, LONSO, OBSERVE, LATSD, LONSD)
       LATSD(1,1) = LATS(K) * rad_to_deg
       LONSD(1,1) = LONS(K) * rad_to_deg
       OBSERVE(1,1) = (abs(LATSD(1,1)-LATSO) < 1.e-3) .and. (abs(LONSD(1,1)-LONSO) < 1.e-3)

       HAVE_ICE: if(FR(K,N) > puny) then

          ERGICE_TMP(:) = ERGICE(K,:,NSUB)
          ERGSNO_TMP(:) = ERGSNO(K,:,NSUB)

          TAUAGE(K,NSUB) = TAUAGE(K,NSUB) + DT

          TRACERS(nt_tsfc, NSUB) = TS(K,N) - TFfresh
          TRACERS(nt_iage, NSUB) = TAUAGE(K, NSUB)
          TRACERS(nt_volpn,NSUB) = VOLPOND(K,NSUB)

          LWUPSRF        = ALW(K) + BLW(K)*TS(K,N) ! use TS (in Kelvin) here
          FSURF(1,1)     = FSWSFC(K,NSUB) - SHF(K) - LHF(K) + LWDNSRF(K) - LWUPSRF
          DLHDT          = EVD(K) * MAPL_ALHS
          DHLAT          = EVD(K) * MAPL_ALHS
          DFSDT          = (SHD(K) + DHLAT + BLW(K))

          SHF0           = -SHF(K)
          LHF0           = -LHF(K)
          LWUP0          = -LWUPSRF
          FSWSFCDB       =  FSWSFC(K,NSUB)
          FSWINTDB       =  FSWINT(K,NSUB)
          FSURFDB(1,1)   =  FSWSFCDB(1,1) - SHF(K) - LHF(K) + LWDNSRF(K) - LWUPSRF
          DFSDTDB        =  DFSDT
          SHF0DB         = -SHF(K)
          LHF0DB         = -LHF(K)
          LWUP0DB        = -LWUPSRF

          TRACERSDB      =  TRACERS(:,NSUB)
          LWDNSRFDB      =  LWDNSRF(K)
          SNODB          =  SNO(K)
          TBOTDB         =  TBOT(K)
          FBOTDB         =  FBOT(K)
          FSWABSDB       =  FSWABS(K)
          !FSWTHRUDB      =  FSWTHRU(K,N-1)
          FSWTHRUDB      =  FSWTHRU(K,N)
          ISWABSDB       =  ISWABS(K,:,NSUB)
          SSWABSDB       =  SSWABS(K,:,NSUB)
          FCONDDB        =  FCOND(K,NSUB)
          FCONDBOTDB     =  FCONDBOT(K,NSUB)
          EVPDB          =  EVP(K)
          FRESHNDB       =  FRESHN(K)
          FSALTNDB       =  FSALTN(K)
          FHOCNNDB       =  FHOCNN(K)
          MELTTDB        =  MELTT(K)
          MELTSDB        =  MELTS(K)
          MELTBDB        =  MELTB(K)
          CONGELDB       =  CONGEL(K)
          SNOICEDB       =  SNOICE(K)
          DSHDB          =  -SHD(K)
          DLHDTDB        =  -DLHDT
          BLWDB          =  -BLW(K)
          LATSDB         =  LATSD(1,1)
          LONSDB         =  LONSD(1,1)
          MLT_ONSETDB    =  0.0
          FRZ_ONSETDB    =  0.0
          YDAYDB         =  0.0
          SBLXDB         =  0.0
          FRDB           =  FR(K,N)
          VOLICEDB       = VOLICE(K,NSUB)
          VOLSNODB       = VOLSNO(K,NSUB)

          call thermo_vertical(                &
               1,1,DTDB,1,one,one,             &
               FRDB,                           &
               TRACERSDB,                      &
               VOLICEDB,     VOLSNODB,         &
               ERGICE_TMP,   ERGSNO_TMP,       &
               LWDNSRFDB,     RDUMDB,          &
               RDUMDB,        RDUMDB,          &
               SNODB,                          &
               FBOTDB,        TBOTDB,          &
               RDUMDB,        RDUMDB,          &
               FSWSFCDB,      FSWINTDB,        &
               FSWTHRUDB,                      &
               SSWABSDB,                       &
               ISWABSDB,                       &

               FSURFDB,       FCONDDB,         &
               SHF0DB,        LHF0DB,          &
               FSWABSDB,      LWUP0DB,         &
               EVPDB,         FRESHNDB,        &
               FSALTNDB,      FHOCNNDB,        &
               MELTTDB,       MELTSDB,         &
               MELTBDB,                        &
               CONGELDB,      SNOICEDB,        &

               DFSDTDB,DSHDB,DLHDTDB,BLWDB,    &
               LATSDB, LONSDB, OBSERVE,        &
               FCONDBOTDB,    SBLXDB,          &

               MLT_ONSETDB,   FRZ_ONSETDB,     &
               YDAYDB,          L_STOP,        &
               IDUM,          JDUM             )

          if(L_STOP) then
             print*, 'CICE_THERMO1: Failing at LAT = ', LATSD, 'LON = ', LONSD
          endif

          _ASSERT(.not.L_STOP,'needs informative message')

          ERGICE(K,:,NSUB)   = ERGICE_TMP(:)
          ERGSNO(K,:,NSUB)   = ERGSNO_TMP(:)
          SHF0               =  SHF0DB
          LHF0               =  LHF0DB
          TRACERS(:, NSUB)   =  TRACERSDB
          FSWTHRU(K,N)       =  FSWTHRUDB(1,1)
          FCOND(K,NSUB)      =  FCONDDB(1,1)
          FCONDBOT(K,NSUB)   =  FCONDBOTDB(1,1)
          FSURF              =  FSURFDB(1,1)
          FSUR(K)            =  FSURF(1,1)
          LWUPSRF            = -LWUP0DB(1,1)
          !*** EVP computed by CICE has an opposite sign:
          !*** condensation > 0, water vapor goes down
          !*** sublimation  < 0, water vapor goes up
          EVP(K)             =  -EVPDB(1,1)
          SBLX(K)            =  SBLXDB(1,1)
          FRESHN(K)          =  FRESHNDB(1,1)
          FSALTN(K)          =  FSALTNDB(1,1)
          FHOCNN(K)          =  FHOCNNDB(1,1)
          MELTT(K)           =  MELTTDB(1,1)
          MELTS(K)           =  MELTSDB(1,1)
          MELTB(K)           =  MELTBDB(1,1)
          CONGEL(K)          =  CONGELDB(1,1)
          SNOICE(K)          =  SNOICEDB(1,1)
          FR(K,N)            =  FRDB(1,1)
          VOLICE(K,N)        =  VOLICEDB(1,1)
          VOLSNO(K,N)        =  VOLSNODB(1,1)
          ! need to update these for aggregation later
          SHF(K)             = -SHF0(1,1)
          LHF(K)             = -LHF0(1,1)
          FSWSFC(K,N)        =  FSWSFCDB(1,1)

          TS(K,N)            =  TRACERS(nt_tsfc,NSUB) + TFfresh
          TAUAGE(K,NSUB)     =  TRACERS(nt_iage,NSUB)

          if ( (DO_POND==1) .and. trim(SHORTWAVE) == 'dEdd') then
             MELTTDB         =  MELTT(K)
             MELTSDB         =  MELTS(K)
             FRDB            =  FR(K,N)
             VOLICEDB        =  VOLICE(K,NSUB)
             VOLSNODB        =  VOLSNO(K,NSUB)
             APONDNDB        =  APONDN(K,NSUB)
             HPONDNDB        =  HPONDN(K,NSUB)
             FRAINDB         =  PCU(K) + PLS(K)
             call compute_ponds(1, 1,                      &
                           1, 1, 1, 1,                     &
                           MELTTDB, MELTSDB, FRAINDB,      &
                           FRDB, VOLICEDB,                 &
                           VOLSNODB, TRACERSDB,            &
                           APONDNDB, HPONDNDB)
             TRACERS(:,NSUB) = TRACERSDB
             VOLPOND(K,NSUB) = TRACERS(nt_volpn,NSUB)
             APONDN (K,NSUB) = APONDNDB(1,1)
             HPONDN (K,NSUB) = HPONDNDB(1,1)
          endif

       end if HAVE_ICE

    !  if(OBSERVE(1)) then
     !   print*, N, FR(K,N)
      !  print*, FBOT(K), FCONDBOT(K,NSUB)
     !   print*, VOLICE(K,NSUB), VOLSNO(K,NSUB)
      !  print*, FHOCNN(K), FSALTN(K)
     !   print*, TS(K,N), real(FR(K,N)*MELTT(K)/DT*8640000, kind=MAPL_R4)
     !   print*, FSURF(1), FCOND(K,NSUB)
     !   print*, FSWSFC(K,NSUB)
     !   print*, LWDNSRF(K), LWUPSRF
     !   print*, SHF(K), LHF(K)

     !   tint = MAPL_UNDEF
     !   call diagnose_internal_ice_temp(VOLICE(K,:), ERGICE(K,:,:), tint)
     !   do l=1,NUM_ICE_LAYERS
     !       print*, l,  tint(ilyr1(N)+l-1)
     !   enddo
     ! endif

    end do TILES ! K loop

    RETURN_(ESMF_SUCCESS)
  end subroutine CICE_THERMO1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: CICE_THERMO2_STEP1 - Computes part 1 of-- 2nd update of LANL CICE Thermodynamics

! !INTERFACE:

  subroutine CICE_THERMO2_STEP1 (NT,ICE,LATS,LONS,LATSO,LONSO,DT,TF,FR,TS,     &
                                 VOLICE,VOLSNO,VOLPOND,ERGICE,ERGSNO,                    &
                                 AICENINIT,VICENINIT,TRCRTYPE,FRCICE,FRZMLT,FRAZLN,      &
                                 FRESHL,FSALTL,FHOCNL,RSIDE,MELTLN,VOLICE_DELTA,         &
                                 TRACERS,TAUAGE,SNOICE,SW,RC)
! not passing TFfresh,saltwatercap,nt_tsfc,nt_iage,nt_volpn,tauage

! !ARGUMENTS:

    integer, optional, intent(OUT) :: RC

    integer, intent(IN)     :: NT               ! number of tiles
    integer, intent(IN)     :: ICE              ! subtiles number assigned to surface type: "ICE"

    real,    intent(IN)     :: LATS     (:)     ! lat
    real,    intent(IN)     :: LONS     (:)     ! lon
    real,    intent(IN)     :: LATSO            ! trace LANL CICE computations at this latitude
    real,    intent(IN)     :: LONSO            ! trace LANL CICE computations at this longitude
    real,    intent(IN)     :: TF       (:)     ! sea Water freezing temperature in degrees C

    real,    intent(IN)     :: SW       (:)     ! ?

    real,    intent(INOUT)  :: TS       (:,:)   ! skin temperature
    real,    intent(INOUT)  :: TRACERS  (:,:)   ! ?
    real,    intent(INOUT)  :: TAUAGE   (:,:)    ! ?

    integer, intent(INOUT)  :: TRCRTYPE (:)     ! ?
    real,    intent(INOUT)  :: FRZMLT   (:)     ! ?
    real,    intent(INOUT)  :: FRAZLN   (:)     ! ?
    real,    intent(INOUT)  :: FRESHL   (:)     ! ?
    real,    intent(INOUT)  :: FSALTL   (:)     ! ?
    real,    intent(INOUT)  :: FHOCNL   (:)     ! ?
    real,    intent(INOUT)  :: RSIDE    (:)     ! ?
    real,    intent(INOUT)  :: MELTLN   (:)     ! ?
    real,    intent(INOUT)  :: SNOICE   (:)     ! ?

    real,    intent(IN)     :: DT
    real(kind=MAPL_R8)      :: DTDB             ! DT (time step) in R8 for CICE
    real(kind=MAPL_R8),    intent(INOUT)  :: FR       (:,:)   ! fractions of water, ice types
    real(kind=MAPL_R8),    intent(INOUT)  :: VOLICE   (:,:)   ! volume of ice
    real(kind=MAPL_R8),    intent(INOUT)  :: VOLSNO   (:,:)   ! volume of snow
    real(kind=MAPL_R8),    intent(INOUT)  :: VOLPOND  (:,:)   ! ?
    real(kind=MAPL_R8),    intent(INOUT)  :: ERGICE   (:,:,:) ! ?
    real(kind=MAPL_R8),    intent(INOUT)  :: ERGSNO   (:,:,:) ! ?
    real(kind=MAPL_R8),    intent(INOUT)  :: AICENINIT(:,:)   ! initial (after cice_init_thermo) ice concentration
    real(kind=MAPL_R8),    intent(INOUT)  :: VICENINIT(:,:)   ! initial (after cice_init_thermo) volume of ice
    real(kind=MAPL_R8),    intent(INOUT)  :: FRCICE   (:)     ! fraction of ice, surface type
    real(kind=MAPL_R8),    intent(INOUT)  :: VOLICE_DELTA  (:,:)! change in volume of ice

! !LOCAL VARIABLES

    character(len=ESMF_MAXSTR)            :: IAm
    integer                               :: STATUS

    integer                               :: K
    logical                               :: OBSERVE(1,1)        ! could be (1,1) to match cice input
    real,               dimension(1,1)    :: LATSD, LONSD
    logical                               :: L_STOP
    integer                               :: IDUM, JDUM
    real                                  :: MINSWFRESH
    real                                  :: FRZMLT_MAX

    real(kind=MAPL_R8)                    :: YDAYDB
    real(kind=MAPL_R8), dimension(1)      :: FRWATERDB, FRZMLTDB, FRAZLNDB, FRESHLDB, FSALTLDB, TFDB, &
                                             RDUMDB, MELTLNDB, FHOCNLDB, RSIDEDB, SNOICEDB, FRCICEDB

    real(kind=MAPL_R8), dimension(NUM_3D_ICE_TRACERS, NUM_ICE_CATEGORIES) :: TRACERSDB2
    real(kind=MAPL_R8)             :: FR_TMP(NUM_ICE_CATEGORIES)
    real(kind=MAPL_R8)             :: VOLICE_TMP(NUM_ICE_CATEGORIES)
    real(kind=MAPL_R8)             :: VOLSNO_TMP(NUM_ICE_CATEGORIES)
    real(kind=MAPL_R8)             :: AICEN_TMP(NUM_ICE_CATEGORIES)
    real(kind=MAPL_R8)             :: VICEN_TMP(NUM_ICE_CATEGORIES)
    real(kind=MAPL_R8)             :: ERGICE_TMP(NUM_ICE_LAYERS,NUM_ICE_CATEGORIES)
    real(kind=MAPL_R8)             :: ERGSNO_TMP(NUM_SNOW_LAYERS,NUM_ICE_CATEGORIES)
    INTEGER                        :: one(1) = (/1/)
    logical                        :: true(1) = (/.true./)


!  !DESCRIPTION:
!        Compute ...??
!          based on CICE

    IAm =  trim(COMP_NAME) // "CICECORE" // "CICE_THERMO2_STEP1"

    DTDB = REAL(DT, kind=MAPL_R8)       ! Convert DT precision: Real4 to Real8 for usage in CICE
    call MAPL_GetResource ( MAPL, MINSWFRESH, Label="FRESH_NEW_ICE_MIN_SALINITY:" , DEFAULT=5.0,    _RC)
    call MAPL_GetResource ( MAPL, FRZMLT_MAX, Label="CICE_FRZMLT_MAX:" , DEFAULT=1000., _RC)

! Loop over all tiles
!-----------------------

    TILES_1: do k=1, NT ! loop over all tiles

       call CICE_INQUIRE_TILE(LATS(K), LONS(K), LATSO, LONSO, OBSERVE, LATSD, LONSD)


       TRACERS(nt_tsfc,:) = TS(K,ICE:) - TFfresh
       TRACERS(nt_iage,:) = TAUAGE(K,:)
       TRACERS(nt_volpn,:)= VOLPOND(K,:)
       TRACERSDB2         = TRACERS
       FRCICEDB           = REAL(FRCICE(K),             kind=MAPL_R8)

       FR_TMP(:)       = FR(K,ICE:)
       VOLICE_TMP(:)   = VOLICE(K,:)
       VOLSNO_TMP(:)   = VOLSNO(K,:)
       ERGICE_TMP(:,:) = ERGICE(K,:,:)
       ERGSNO_TMP(:,:) = ERGSNO(K,:,:)

       if(FRCICE(K) > 0.0) then
          FRWATERDB  =  1.0_8 - sum(FR(K,ICE:))
          AICEN_TMP(:) = AICENINIT(K,:)
          VICEN_TMP(:) = VICENINIT(K,:)
          call linear_itd (1,1,1,one,one, &
                           TRCRTYPE,          &
                           AICEN_TMP,         &
                           VICEN_TMP,         &
                           FR_TMP,            &
                           TRACERSDB2,        &
                           VOLICE_TMP,        &
                           VOLSNO_TMP,        &
                           ERGICE_TMP,        &
                           ERGSNO_TMP,        &
                           FRCICEDB,          &
                           FRWATERDB,         &
                           LATSD, LONSD,      &
                           L_STOP,            &
                           IDUM, JDUM )

          if(L_STOP) then
             print*, 'CICE_THERMO2_STEP1: after linear_itd. Failing at LAT = ', LATSD, 'LON = ', LONSD
          endif

          _ASSERT(.not.L_STOP,'needs informative message')

       endif

       FRZMLTDB       =  FRZMLT(K)
       FRZMLTDB       =  min(max(FRZMLTDB,-FRZMLT_MAX),FRZMLT_MAX)
       FRAZLNDB       =  FRAZLN(K)
       FRESHLDB       =  FRESHL(K)
       FSALTLDB       =  FSALTL(K)
       TFDB           =  TF(K)
       RDUMDB         =  0.0
       YDAYDB         =  0.0
       FRWATERDB      =  1.0_8 - sum(FR_TMP)

       call add_new_ice (1,1,1,one,one,true, DTDB,           &
                         FR_TMP,                             &
                         TRACERSDB2,                         &
                         VOLICE_TMP,                         &
                         ERGICE_TMP,                         &
                         FRWATERDB,                          &
                         FRCICEDB ,                          &
                         FRZMLTDB,                           &
                         FRAZLNDB,                           &
                         SW(K) < MINSWFRESH,                 &
                         RDUMDB, YDAYDB,                     &
                         FRESHLDB,                           &
                         FSALTLDB,                           &
                         TFDB, L_STOP,                       &
                         IDUM, JDUM)
       if(L_STOP) then
          print*, 'CICE_THERMO2_STEP1: after add_new_ice. Failing at LAT = ', LATSD, 'LON = ', LONSD
       endif

       _ASSERT(.not.L_STOP,'needs informative message')

       FRAZLN(K)   =  FRAZLNDB (1)

       !if(FRZMLT(K) > 0.) NEWICEERG(K) = NEWICEERG(K) - FRZMLT(K)

       VOLICE_DELTA(k,:)  =  VOLICE_TMP
       FHOCNLDB      =  FHOCNL(K)
       RSIDEDB       =  RSIDE(K)
       MELTLNDB      =  MELTLN(K)

       call lateral_melt (1,1,1,1,1,1,DTDB, &
                          FRESHLDB,         &
                          FSALTLDB,         &
                          FHOCNLDB,         &
                          RSIDEDB,          &
                          MELTLNDB,         &
                          FR_TMP,           &
                          VOLICE_TMP,       &
                          VOLSNO_TMP,       &
                          ERGICE_TMP,       &
                          ERGSNO_TMP )

       VOLICE_DELTA(k,:)  = VOLICE_DELTA(k,:) - VOLICE_TMP
       SNOICEDB      = 0.0

       call freeboard_ccsm (1,1,1,1,1,1, DTDB,         &
                            FR_TMP,                    &
                            VOLICE_TMP,  VOLSNO_TMP,   &
                            ERGICE_TMP,                &
                            ERGSNO_TMP,                &
                            SNOICEDB,                  &
                            FSALTLDB)

       SNOICE(K) = SNOICEDB(1)
       FRWATERDB  =  1.0_8 - sum(FR_TMP)

       call cleanup_itd (1,1,1,1,1,1,DTDB, &
                         FR_TMP,           &
                         TRACERSDB2,       &
                         VOLICE_TMP,       &
                         VOLSNO_TMP,       &
                         ERGICE_TMP,       &
                         ERGSNO_TMP,       &
                         FRWATERDB,        &
                         FRCICEDB ,        &
                         TRCRTYPE,         &
                         FRESHLDB,         &
                         FSALTLDB,         &
                         FHOCNLDB,         &
                         .true., L_STOP,   &
                         IDUM, JDUM,       &
                         limit_aice_in=.true. &
                         ,punynum=puny)

       _ASSERT(.not.L_STOP,'needs informative message')

       FR(K,ICE:)    = FR_TMP(:)
       VOLICE(K,:)   = VOLICE_TMP(:)
       VOLSNO(K,:)   = VOLSNO_TMP(:)
       ERGICE(K,:,:) = ERGICE_TMP(:,:)
       ERGSNO(K,:,:) = ERGSNO_TMP(:,:)

       TRACERS       = TRACERSDB2
       FRESHL(K)     = FRESHLDB (1)
       FSALTL(K)     = FSALTLDB (1)
       FHOCNL(K)     = FHOCNLDB (1)
       MELTLN(K)     = MELTLNDB (1)

       TS(K,ICE:)    = TRACERS(nt_tsfc,:) + TFfresh
       TAUAGE(K,:)   = TRACERS(nt_iage,:)
       VOLPOND(K,:)  = TRACERS(nt_volpn,:)

       !if(OBSERVE(1)) then
       !  print*, 'after therm2: ',FRESHL(K)
       !  print*, 'after therm2: ',FSALTL(K)
       !  print*, 'after therm2: ',FHOCNL(K)
       !endif

    end do TILES_1 ! K loop


    RETURN_(ESMF_SUCCESS)
  end subroutine CICE_THERMO2_STEP1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: CICE_THERMO2_STEP2 - Computes 2st step of CICE Thermodynamics

! !INTERFACE:

  subroutine CICE_THERMO2_STEP2 (NT,ICE,LATS,LONS,LATSO,LONSO,DT,FR,TS,                &
                                 VOLICE,VOLSNO,VOLPOND,ERGICE,ERGSNO,                          &
                                 TRCRTYPE,FRCICE,SLMASK,TRACERS,TAUAGE,RC)
! not passing TFfresh,puny,tauage,tracers

! !ARGUMENTS:
    integer, optional, intent(OUT) :: RC

    integer, intent(IN)     :: NT               ! number of tiles
    integer, intent(IN)     :: ICE              ! subtiles number assigned to surface type: "ICE"

    real,    intent(IN)     :: LATS     (:)     ! lat
    real,    intent(IN)     :: LONS     (:)     ! lon
    real,    intent(IN)     :: LATSO            ! trace CICE computations at this latitude
    real,    intent(IN)     :: LONSO            ! trace CICE computations at this longitude
    real,    intent(IN)     :: SLMASK   (:)     ! "salt water lake mask"

    real,    intent(INOUT)  :: TS       (:,:)   ! skin temperature
    real,    intent(INOUT)  :: TRACERS  (:,:)   ! ?
    real,    intent(INOUT)  :: TAUAGE   (:,:)   ! ?

    integer, intent(INOUT)  :: TRCRTYPE (:)     ! ?

    real,    intent(IN)     :: DT
    real(kind=MAPL_R8)      :: DTDB             ! DT (time step) in R8 for CICE
    real(kind=MAPL_R8),    intent(INOUT)  :: VOLICE   (:,:)   ! volume of ice
    real(kind=MAPL_R8),    intent(INOUT)  :: VOLSNO   (:,:)   ! volume of snow
    real(kind=MAPL_R8),    intent(INOUT)  :: FR       (:,:)   ! fractions of water, ice types
    real(kind=MAPL_R8),    intent(INOUT)  :: FRCICE   (:)     ! fraction of ice, surface type
    real(kind=MAPL_R8),    intent(INOUT)  :: VOLPOND  (:,:)   ! ?
    real(kind=MAPL_R8),    intent(INOUT)  :: ERGICE   (:,:,:) ! ?
    real(kind=MAPL_R8),    intent(INOUT)  :: ERGSNO   (:,:,:) ! ?

! !LOCAL VARIABLES
    character(len=ESMF_MAXSTR)     :: IAm
    integer                        :: STATUS


    logical                :: OBSERVE(1,1)        ! could be (1,1) to match cice input
    real,  dimension(1,1)  :: LATSD, LONSD
    logical                :: L_STOP
    integer                :: IDUM, JDUM
    integer                :: k,l,n
    real                   :: ICE_THICKNESS_THRESH
    real                   :: ICE_ARTIFICIAL_MELT

    real(kind=MAPL_R8)     :: hid, hi, hs, hsn

    real(kind=MAPL_R8)     :: qin_save   (NUM_ICE_CATEGORIES, NUM_ICE_LAYERS)
    real(kind=MAPL_R8)     :: qsn_save   (NUM_ICE_CATEGORIES, NUM_SNOW_LAYERS)
    real(kind=MAPL_R8)     :: TRACERSDB2 (NUM_3D_ICE_TRACERS, NUM_ICE_CATEGORIES)

    real(kind=MAPL_R8), dimension(1)  :: FRWATERDB, FRESHLDB, FSALTLDB, FHOCNLDB, FRCICEDB
    real(kind=MAPL_R8)                :: FR_TMP(NUM_ICE_CATEGORIES)
    real(kind=MAPL_R8)                :: VOLICE_TMP(NUM_ICE_CATEGORIES)
    real(kind=MAPL_R8)                :: VOLSNO_TMP(NUM_ICE_CATEGORIES)
    real(kind=MAPL_R8)                :: ERGICE_TMP(NUM_ICE_LAYERS,NUM_ICE_CATEGORIES)
    real(kind=MAPL_R8)                :: ERGSNO_TMP(NUM_SNOW_LAYERS,NUM_ICE_CATEGORIES)


!  !DESCRIPTION:
!        Compute ...??
!          based on CICE

    IAm =  trim(COMP_NAME) // "CICECORE" // "CICE_THERMO2_STEP2"

    DTDB = REAL(DT, kind=MAPL_R8)       ! Convert DT precision: Real4 to Real8 for usage in CICE

    call MAPL_GetResource ( MAPL, ICE_THICKNESS_THRESH, Label="CICE_ICE_THICKNESS_THRESH:", DEFAULT=1.5, _RC)
    call MAPL_GetResource ( MAPL, ICE_ARTIFICIAL_MELT,  Label="CICE_ICE_ARTIFICIAL_MELT:" , DEFAULT=0.1, _RC)
    ! the units of ICE_ARTIFICIAL_MEL are cm/day and it is converted to m/time step
    hid = real(ICE_ARTIFICIAL_MELT*1.e-2*DT/86400.0, kind=8)

! Loop over all tiles
!-----------------------

    TILES_2: do k=1, NT ! loop over all tiles

       FRCICEDB = REAL(FRCICE(k),             kind=MAPL_R8)

       if( (SLMASK(K) > 0.5) .and. (FRCICE(K) > 0.0) .and. &
           (sum(VOLICE(K,:)) > ICE_THICKNESS_THRESH)) then

          call CICE_INQUIRE_TILE(LATS(K), LONS(K), LATSO, LONSO, OBSERVE, LATSD, LONSD)

          qin_save(:,:) = 0.0
          qsn_save(:,:) = 0.0

          loop_over_ice_cat: do n=1, NUM_ICE_CATEGORIES
             hi = 0.0
             hs = 0.0
             !if(VOLICE(k,n) > puny) then
             if(FR(k,n) > puny) then

                hi = VOLICE(k,n) / FR(k,n)
                do l=1,NUM_ICE_LAYERS
                   qin_save(n,l) = ERGICE(k,l,n)*REAL(NUM_ICE_LAYERS,kind=MAPL_R8)/VOLICE(k,n)
                enddo

                if(VOLSNO(k,n) > c0) then
                   hs = VOLSNO(k,n) / FR(k,n)
                   do l=1,NUM_SNOW_LAYERS
                      qsn_save(n,l) = ERGSNO(k,l,n)*REAL(NUM_SNOW_LAYERS,kind=MAPL_R8)/VOLSNO(k,n)
                   enddo
                endif

                if(hi > hid) hi = hi - hid
                hsn = (MAPL_RHO_SEAWATER-MAPL_RHO_SEAICE)/MAPL_RHO_SNOW *hi
                if(hs > hsn) hs = hsn

                VOLICE(k,n) = hi * FR(k,n)
                VOLSNO(k,n) = hs * FR(k,n)

                do l=1,NUM_ICE_LAYERS
                   ERGICE(k,l,n) = qin_save(n,l)*VOLICE(k,n)/real(NUM_ICE_LAYERS,kind=MAPL_R8)
                enddo

                if(VOLSNO(k,n) > puny) then
                   do l=1,NUM_SNOW_LAYERS
                         ERGSNO(k,l,n) = qsn_save(n,l)*VOLSNO(k,n)/real(NUM_SNOW_LAYERS,kind=MAPL_R8)
                   enddo
                endif
             endif
          enddo loop_over_ice_cat

          TRACERS(nt_tsfc,:) = TS(K,ICE:)-TFfresh
          TRACERS(nt_iage,:) = TAUAGE(K,:)
          TRACERS(nt_volpn,:)= VOLPOND(K,:)

          FRWATERDB          = 1.0 - sum(FR(K,:))
          TRACERSDB2         = TRACERS

          FR_TMP(:)          = FR(K,ICE:)
          VOLICE_TMP(:)      = VOLICE(K,:)
          VOLSNO_TMP(:)      = VOLSNO(K,:)
          ERGICE_TMP(:,:)    = ERGICE(K,:,:)
          ERGSNO_TMP(:,:)    = ERGSNO(K,:,:)
          FRESHLDB(:)        = 0.0_8
          FSALTLDB(:)        = 0.0_8
          FHOCNLDB(:)        = 0.0_8

          call cleanup_itd (1,1,1,1,1,1,DTDB, &
                            FR_TMP,           &
                            TRACERSDB2,       &
                            VOLICE_TMP,       &
                            VOLSNO_TMP,       &
                            ERGICE_TMP,       &
                            ERGSNO_TMP,       &
                            FRWATERDB,        &
                            FRCICEDB,         &
                            TRCRTYPE,         &
                            FRESHLDB,         &
                            FSALTLDB,         &
                            FHOCNLDB,         &
                            .true., L_STOP,   &
                            IDUM, JDUM,       &
                            limit_aice_in=.true.)

          if(L_STOP) then
             print*, 'CICE_THERMO2_STEP2: Failing at LAT = ', LATSD, 'LON = ', LONSD
          endif

          _ASSERT(.not.L_STOP,'needs informative message')

          FR(K,ICE:)    = FR_TMP(:)
          VOLICE(K,:)   = VOLICE_TMP(:)
          VOLSNO(K,:)   = VOLSNO_TMP(:)
          ERGICE(K,:,:) = ERGICE_TMP(:,:)
          ERGSNO(K,:,:) = ERGSNO_TMP(:,:)

          TRACERS      = TRACERSDB2

          TS(K,ICE:)   = TRACERS(nt_tsfc,:) + TFfresh
          TAUAGE(K,:)  = TRACERS(nt_iage,:)
          VOLPOND(K,:) = TRACERS(nt_volpn,:)
       endif
    end do TILES_2 ! K loop


    RETURN_(ESMF_SUCCESS)
  end subroutine CICE_THERMO2_STEP2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: CICE_ALBSEAICE - Computes albedos using LANL CICE

! !INTERFACE:

  subroutine CICE_ALBSEAICE (ICE,NUM_ICE_CATEGORIES,NUM_ICE_LAYERS,NUM_SNOW_LAYERS,NT,DO_POND,LATSO,LONSO,LATS,LONS,ZTH,FR,TS, &
                             DRPAR,DFPAR,DRNIR,DFNIR,DRUVR,DFUVR,VSUVR,VSUVF,VOLICE,VOLSNO,APONDN,HPONDN,                      &
                             ISWABS,FSWSFC, FSWINT,FSWTHRU,SSWABS,ALBIN,ALBSN,ALBPND,ALBVRN,ALBVFN,ALBNRN,ALBNFN,              &
                             DRUVRTHRU,DFUVRTHRU,DRPARTHRU,DFPARTHRU,RC)

  use ice_shortwave,      only: shortwave_ccsm3, shortwave_dEdd_set_snow, &
                                shortwave_dEdd_set_pond,                  &
                                shortwave_dEdd

! !ARGUMENTS:
    integer, optional, intent(OUT) :: RC

    integer, intent(IN)  :: ICE                  ! starting index of ICE tiles
    integer, intent(IN)  :: NUM_ICE_CATEGORIES   ! # of ice categories
    integer, intent(IN)  :: NUM_ICE_LAYERS       ! # of ice  layers
    integer, intent(IN)  :: NUM_SNOW_LAYERS      ! # of snow layers
    integer, intent(IN)  :: NT                   ! # tiles in each type
    integer, intent(IN)  :: DO_POND              ! cice_do_pond resource parameter

    real                 :: LATSO                ! trace cice computations at Lat: LATSO
    real                 :: LONSO                ! trace cice computations at Lon: LONSO

    real,    intent(IN)  :: LATS(:)              ! latitudes
    real,    intent(IN)  :: LONS(:)              ! longitudes
    real,    intent(IN)  :: ZTH (:)              ! cosine of solar zenith angle
    real,    intent(IN)  :: TS (:,:)             ! Skin temperature
    real,    intent(IN)  :: DRPAR (:)
    real,    intent(IN)  :: DFPAR (:)
    real,    intent(IN)  :: DRNIR (:)
    real,    intent(IN)  :: DFNIR (:)
    real,    intent(IN)  :: DRUVR (:)
    real,    intent(IN)  :: DFUVR (:)
    real,    intent(IN)  :: VSUVR (:)
    real,    intent(IN)  :: VSUVF (:)

    real,    intent(OUT)  :: ISWABS (:,:,:)
    real,    intent(INOUT):: SSWABS (:,:,:)
    real,    intent(OUT)  :: FSWSFC (:,:)
    real,    intent(OUT)  :: FSWINT (:,:)
    real,    intent(OUT)  :: FSWTHRU(:,:)

    real,    intent(INOUT)  :: ALBIN (:,:)
    real,    intent(INOUT)  :: ALBSN (:,:)
    real,    intent(INOUT)  :: ALBPND(:,:)

    real,    intent(OUT)  :: ALBVRN (:,:)       ! visible direct  albedo
    real,    intent(OUT)  :: ALBVFN (:,:)       ! visible diffuse albedo
    real,    intent(OUT)  :: ALBNRN (:,:)       ! nearIr  direct  albedo
    real,    intent(OUT)  :: ALBNFN (:,:)       ! nearIr  diffuse albedo

    real,    intent(OUT)  :: DRUVRTHRU (:,:)    ! direct  UV  ??
    real,    intent(OUT)  :: DFUVRTHRU (:,:)    ! diffuse UV  ??
    real,    intent(OUT)  :: DRPARTHRU (:,:)    ! direct  PAR ??
    real,    intent(OUT)  :: DFPARTHRU (:,:)    ! diffuse PAR ??

    real(kind=MAPL_R8),    intent(IN)  :: FR     (:,:)             ! Fraction of ice in each category
    real(kind=MAPL_R8),    intent(IN)  :: VOLICE (:,:)
    real(kind=MAPL_R8),    intent(IN)  :: VOLSNO (:,:)
    real(kind=MAPL_R8),    intent(IN)  :: APONDN (:,:)
    real(kind=MAPL_R8),    intent(IN)  :: HPONDN (:,:)

!  !LOCAL VARIABLES
    character(len=ESMF_MAXSTR)     :: IAm
    integer                        :: STATUS

    integer               :: N, NSUB, K
    INTEGER               :: one(1) = (/1/)

    real(kind=MAPL_R8), dimension(NUM_ICE_LAYERS)  :: ISWABSDB
    real(kind=MAPL_R8), dimension(NUM_SNOW_LAYERS) :: SSWABSDB
    real(kind=MAPL_R8), dimension(NUM_SNOW_LAYERS) :: RHOSNWN
    real(kind=MAPL_R8), dimension(NUM_SNOW_LAYERS) :: RSNWN

    real(kind=MAPL_R8), dimension(1) ::                                 &
                                           TSCDB,                       &
                                           DRPARDB,       DFPARDB,      &
                                           DRNIRDB,       DFNIRDB,      &
                                           DRUVRDB,       DFUVRDB,      &
                                           VSUVRDB,       VSUVFDB,      &
                                           ALBVRNDB,      ALBNRNDB,     &
                                           ALBVFNDB,      ALBNFNDB,     &
                                           FSWSFCDB,      FSWINTDB,     &
                                           FSWTHRUDB,                   &
                                           ALBINDB,       ALBSNDB,      &
                                           ALBPNDDB,                    &
                                           FRDB,                        &
                                           VOLICEDB,      VOLSNODB,     &
                                           DRUVRTHRUDB,   DFUVRTHRUDB,  &
                                           DRPARTHRUDB,   DFPARTHRUDB,  &
                                           COSZTH,        FSN,          &
                                           FPN,           HPN

    real,               dimension(1,1)    :: LATSD, LONSD
    logical,            dimension(1,1)    :: OBSERVE           ! could be (1,1) to match cice input
    logical                               :: TR_POND
    character(len=ESMF_MAXSTR)            :: SHORTWAVE

!  !DESCRIPTION:
!        Compute albedo over sea-ice using: Delta-Eddington or CCSM3 (default)
!          based on CICE

    IAm =  trim(COMP_NAME) // "CICECORE" // "CICE_ALBSEAICE"

    call MAPL_GetResource ( MAPL, SHORTWAVE, Label="CICE_SHORTWAVE:" , DEFAULT="shortwave_ccsm" , _RC)

    if (DO_POND == 1) then
       TR_POND = .true.
    else
       TR_POND = .false.
    endif

!  Initialize output
!-------------------

    ALBVRN     = 0.0
    ALBNRN     = 0.0
    ALBVFN     = 0.0
    ALBNFN     = 0.0

    DRUVRTHRU  = 0.0
    DFUVRTHRU  = 0.0
    DRPARTHRU  = 0.0
    DFPARTHRU  = 0.0

! Loop over all subtiles
!-----------------------

    EVERY_SUBTILE: do N=1, NUM_ICE_CATEGORIES
       EVERY_TILE: do K=1, NT

         call CICE_INQUIRE_TILE(LATS(K), LONS(K), LATSO, LONSO, OBSERVE, LATSD, LONSD)

         if(FR(K,N) > puny) then

             TSCDB      =  REAL(TS(K,N) - TFfresh,    kind=MAPL_R8)  ! Convert Inputs precision: Real4 to Real8
             DRPARDB    =  REAL(DRPAR(K),             kind=MAPL_R8)
             DFPARDB    =  REAL(DFPAR(K),             kind=MAPL_R8)
             DRNIRDB    =  REAL(DRNIR(K),             kind=MAPL_R8)
             DFNIRDB    =  REAL(DFNIR(K),             kind=MAPL_R8)
             DRUVRDB    =  REAL(DRUVR(K),             kind=MAPL_R8)
             DFUVRDB    =  REAL(DFUVR(K),             kind=MAPL_R8)
             VSUVRDB    =  REAL(VSUVR(K),             kind=MAPL_R8)
             VSUVFDB    =  REAL(VSUVF(K),             kind=MAPL_R8)

             ALBVRNDB   =  REAL(ALBVRN(K,N),          kind=MAPL_R8)        ! Initialize R8 (Double Precision) outputs
             ALBNRNDB   =  REAL(ALBNRN(K,N),          kind=MAPL_R8)
             ALBVFNDB   =  REAL(ALBVFN(K,N),          kind=MAPL_R8)
             ALBNFNDB   =  REAL(ALBNFN(K,N),          kind=MAPL_R8)

             FSWSFCDB   =  REAL(FSWSFC(K,N),          kind=MAPL_R8)
             FSWINTDB   =  REAL(FSWINT(K,N),          kind=MAPL_R8)
             FSWTHRUDB  =  REAL(FSWTHRU(K,N),         kind=MAPL_R8)
             ISWABSDB   =  REAL(ISWABS(K,:,N),        kind=MAPL_R8)

             ALBINDB    =  REAL(ALBIN(K,N),           kind=MAPL_R8)
             ALBSNDB    =  REAL(ALBSN(K,N),           kind=MAPL_R8)
             ALBPNDDB   =  REAL(ALBPND(K,N),          kind=MAPL_R8)

             FRDB       =  REAL(FR(K,N),              kind=MAPL_R8)

             VOLICEDB   =  REAL(VOLICE(K,N),          kind=MAPL_R8)
             VOLSNODB   =  REAL(VOLSNO(K,N),          kind=MAPL_R8)

             if (trim(SHORTWAVE) == 'dEdd') then
                SSWABSDB   =  SSWABS(K,:,N)
                COSZTH     =  ZTH(K)                         !*** ZTH is actually cos() of solar zenith angle
                call shortwave_dEdd_set_snow(1, 1,        &  ! set snow properties
                              1, one, one,                &
                              FRDB,     VOLSNODB,         &
                              TSCDB,    FSN,              &
                              RHOSNWN,  RSNWN)
                if (.not. TR_POND) then
                   call shortwave_dEdd_set_pond(1, 1,     &  ! set pond properties
                                 1,  one, one,            &
                                 FRDB,   TSCDB,           &
                                 FSN,    FPN,             &
                                 HPN)
                else
                   FPN = REAL(APONDN(K, N), kind=MAPL_R8)
                   HPN = REAL(HPONDN(K, N), kind=MAPL_R8)
                endif
                call shortwave_dEdd(1,        1,          &
                                1, one, one,              &
                                COSZTH,                   &  ! Inputs
                                FRDB,      VOLICEDB,      &
                                VOLSNODB,  FSN,           &
                                RHOSNWN,   RSNWN,         &
                                FPN,       HPN,           &
                                OBSERVE,                  &
                                DRUVRDB,   DFUVRDB,       &
                                DRPARDB,   DFPARDB,       &
                                VSUVRDB,   VSUVFDB,       &
                                DRNIRDB,   DFNIRDB,       &
                                ALBVRNDB,  ALBVFNDB,      &  ! Outputs: order of the following 4 parms is different from shortwave_ccsm3
                                ALBNRNDB,  ALBNFNDB,      &
                                FSWSFCDB,  FSWINTDB,      &
                                FSWTHRUDB, SSWABSDB,      &
                                           ISWABSDB,      &
                                DRUVRTHRUDB,   DFUVRTHRUDB, &
                                DRPARTHRUDB,   DFPARTHRUDB, &
                                ALBINDB,  ALBSNDB,   ALBPNDDB  )

                  SSWABS(K,:,N)     =  REAL(SSWABSDB, kind=MAPL_R4)
                  ALBPND(K,N)       =  REAL(ALBPNDDB(1), kind=MAPL_R4)
             else
                call shortwave_ccsm3 (1,      1,          &
                                1, one, one,              &
                                OBSERVE,                  &  ! Inputs
                                DRUVRDB,       DFUVRDB,   &
                                DRPARDB,       DFPARDB,   &
                                FRDB,          VOLICEDB,  &
                                VOLSNODB,      TSCDB,     &
                                VSUVRDB,       VSUVFDB,   &
                                DRNIRDB,       DFNIRDB,   &
                                ALBVRNDB,      ALBNRNDB,  &  ! Outputs
                                ALBVFNDB,      ALBNFNDB,  &
                                FSWSFCDB,      FSWINTDB,  &
                                FSWTHRUDB,     ISWABSDB,  &
                                DRUVRTHRUDB,   DFUVRTHRUDB, &
                                DRPARTHRUDB,   DFPARTHRUDB, &
                                ALBINDB,       ALBSNDB         )
             endif

             ALBVRN(K,N)       = REAL(ALBVRNDB(1),    kind=MAPL_R4)
             ALBNRN(K,N)       = REAL(ALBNRNDB(1),    kind=MAPL_R4)
             ALBVFN(K,N)       = REAL(ALBVFNDB(1),    kind=MAPL_R4)
             ALBNFN(K,N)       = REAL(ALBNFNDB(1),    kind=MAPL_R4)

             DRUVRTHRU(K,N)    = REAL(DRUVRTHRUDB(1), kind=MAPL_R4)
             DFUVRTHRU(K,N)    = REAL(DFUVRTHRUDB(1), kind=MAPL_R4)
             DRPARTHRU(K,N)    = REAL(DRPARTHRUDB(1), kind=MAPL_R4)
             DFPARTHRU(K,N)    = REAL(DFPARTHRUDB(1), kind=MAPL_R4)

             FSWSFC (K,N)      = REAL(FSWSFCDB(1),    kind=MAPL_R4)
             FSWINT (K,N)      = REAL(FSWINTDB(1),    kind=MAPL_R4)
             FSWTHRU(K,N)      = REAL(FSWTHRUDB(1),   kind=MAPL_R4)
             ISWABS (K,:,N)    = REAL(ISWABSDB,       kind=MAPL_R4)

             ALBIN(K,N)        = REAL(ALBINDB(1),     kind=MAPL_R4)
             ALBSN(K,N)        = REAL(ALBSNDB(1),     kind=MAPL_R4)
         end if

       end do EVERY_TILE
    end do EVERY_SUBTILE


    RETURN_(ESMF_SUCCESS)
  end subroutine CICE_ALBSEAICE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: CICE_ICE_BUDGET - Computes and prints: flux of total mass and fresh water

! !INTERFACE:

  subroutine CICE_ICE_BUDGET (DT, NT, thermo_, VMG, TILEAREA, TOTALFLUX,                      &
                              SLMASK, VOLICE, VOLICE_OLD, VOLSNO, VOLSNO_OLD)

! !ARGUMENTS:

    integer, intent(IN)     :: thermo_
    real,    intent(IN)     :: DT               ! time-step
    integer, intent(IN)     :: NT               ! number of tiles
    type(ESMF_VM), intent(IN) :: VMG
    real,    intent(IN)     :: TILEAREA  (:)
    real,    intent(IN)     :: TOTALFLUX (:)
    real,    intent(IN)     :: SLMASK   (:)     ! "salt water lake mask"

    real(kind=MAPL_R8), intent(IN)  :: VOLICE       (:,:)
    real(kind=MAPL_R8), intent(IN)  :: VOLSNO       (:,:)
    real(kind=MAPL_R8), intent(IN)  :: VOLICE_OLD   (:)
    real(kind=MAPL_R8), intent(IN)  :: VOLSNO_OLD   (:)

!  !LOCAL VARIABLES
    real                    :: TOTALAREA, ALLTOTALAREA

    real(kind=MAPL_R8), dimension(NT)  :: CICEDMASS

!  !DESCRIPTION:
!        Compute total mass flux and fresh water flux based on CICE Thermodynamics

    CICEDMASS = MAPL_RHO_SEAICE*(sum(VOLICE,dim=2)-VOLICE_OLD) + MAPL_RHO_SNOW*(sum(VOLSNO,dim=2)-VOLSNO_OLD)
    where(SLMASK > 0.5)
       CICEDMASS = 0.0
    endwhere
    TOTALAREA = sum(CICEDMASS*TILEAREA)

    call ESMF_VMBarrier(VMG, _RC)
    call MAPL_CommsAllReduceSum(VMG, TOTALAREA, ALLTOTALAREA, 1, _RC)

    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' After Thermo ', thermo_, '******************* '
    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' total ice+sno mass change = ', &
                                 ALLTOTALAREA

    TOTALAREA = sum(TOTALFLUX * DT*TILEAREA)

    call ESMF_VMBarrier(VMG, _RC)
    call MAPL_CommsAllReduceSum(VMG, TOTALAREA, ALLTOTALAREA, 1, _RC)

    if(MAPL_AM_I_ROOT()) print*, trim(Iam), ' total freshwaterflux * dt = ', &
                                 ALLTOTALAREA

    RETURN_(ESMF_SUCCESS)
  end subroutine CICE_ICE_BUDGET

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: CICE_INQUIRE_TILE -
!  Computes OBSERVE which is useful to trace thru LANL CICE computations at a single lat/lon location
!  that location is specified by latso, lonso (resource paramater controlled)

! !INTERFACE:

  subroutine CICE_INQUIRE_TILE(LAT, LON, LATSO, LONSO, OBSERVE, LATSD, LONSD)

  use ice_constants,      only: rad_to_deg

! !ARGUMENTS:

    real,    intent(IN)  :: LAT
    real,    intent(IN)  :: LON
    real,    intent(IN)  :: LATSO
    real,    intent(IN)  :: LONSO

    logical, intent(OUT) :: OBSERVE(1,1)
    real,    intent(OUT) :: LATSD(1,1), LONSD(1,1)

!  !LOCAL VARIABLES

!  !DESCRIPTION:

    LATSD = LAT *  rad_to_deg
    LONSD = LON *  rad_to_deg
    OBSERVE = (abs(LATSD-LATSO) < 1.e-3) .and. (abs(LONSD-LONSO) < 1.e-3)

   RETURN_(ESMF_SUCCESS)
  end subroutine CICE_INQUIRE_TILE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine RUN2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: Finalize        -- Finalize method for CICEThermo wrapper

! !INTERFACE:

  subroutine Finalize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component
  type(ESMF_State),    intent(INOUT) :: import ! Import state
  type(ESMF_State),    intent(INOUT) :: export ! Export state
  type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
  integer, optional,   intent(  OUT) :: rc     ! Error code:

!EOP

    type (MAPL_MetaComp), pointer:: MAPL

! ErrLog Variables

    character(len=ESMF_MAXSTR)       :: IAm
    integer                          :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Finalize"
    call ESMF_GridCompGet( gc, NAME=comp_name, _RC )
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL"   )
    call MAPL_TimerOn(MAPL,"FINALIZE")

    call dealloc_column_physics( MAPL_AM_I_Root(), Iam )

    call MAPL_TimerOff(MAPL,"FINALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL"   )

! Generic Finalize
! ------------------

    call MAPL_GenericFinalize( GC, IMPORT, EXPORT, CLOCK, _RC )

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Finalize


  subroutine Normalize(ptr, frac)

     real,              dimension(:),  intent(inout)  :: ptr
     real(KIND=MAPL_R8),dimension(:),     intent(in)  :: frac

     where(frac > puny)
         ptr = ptr / real(frac, kind=MAPL_R4)
     endwhere

     return
  end subroutine Normalize

  subroutine FreezingTemperature(t, s, smin, prescribed, kelvin)

     real,       dimension(:),  intent(out)  :: t
     real,       dimension(:),  intent(in)   :: s
     real,                      intent(in)   :: smin
     logical,                   intent(in)   :: prescribed
     logical,    optional,      intent(in)   :: kelvin


    if (prescribed) then
       t = -1.8           ! constant freezing temp (C)
    else
       t = -depressT * s  ! CICE default mode: TF = -depressT * SSS. depressT is CICE constant.
       where(s < smin)
          t = -depressT * smin  ! limit TF where salinity too low
       endwhere
    endif

    if(present(kelvin)) then
        if(kelvin) t = t + TFfresh
    endif

    return
  end subroutine FreezingTemperature

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: ALBSEAICEM2 - Computes albedos as a function of  $cos(\zeta)$ over sea-ice surfaces

! !INTERFACE:

  subroutine ALBSEAICEM2 (ALBVR,ALBVF,ALBNR,ALBNF,ZTH,LATS,currTime)

! !ARGUMENTS:

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS,RC
    character(len=ESMF_MAXSTR)              :: COMP_NAME

    type(ESMF_Time), intent(INOUT)    :: currTime
    real,    intent(IN)  :: LATS(:)
    real,    intent(IN)  :: ZTH  (:)
    real,    intent(OUT) :: ALBVR(:) ! visible beam    albedo
    real,    intent(OUT) :: ALBVF(:) ! visible diffuse albedo
    real,    intent(OUT) :: ALBNR(:) ! nearIr  beam    albedo
    real,    intent(OUT) :: ALBNF(:) ! nearIr  diffuse albedo

!  !DESCRIPTION:
!        Compute albedo for ocean points
!          based on Heracles GEOS-AGCM

      !real, parameter :: SEAICEALBVR  = .60
      !real, parameter :: SEAICEALBVF  = .60
      !real, parameter :: SEAICEALBNR  = .60
      !real, parameter :: SEAICEALBNF  = .60
      real                  :: SEAICEALBVRN
      real                  :: SEAICEALBVFN
      real                  :: SEAICEALBNRN
      real                  :: SEAICEALBNFN

      real                  :: SEAICEALBVRS
      real                  :: SEAICEALBVFS
      real                  :: SEAICEALBNRS
      real                  :: SEAICEALBNFS

      real, dimension(0:13) :: shebavis
      real, dimension(0:13) :: shebanir
      real, dimension(0:13) :: nday
      real                  :: afracv,aslopev,abasev
      real                  :: afraci,aslopei,abasei

      character(len=ESMF_MAXSTR)     :: string
      integer               :: YEAR,MONTH,DAY


      shebavis = (/0.820,0.820,0.820,0.820,0.820,0.820,0.751, &
       0.467,0.663,0.820,0.820,0.820,0.820,0.820/)
      shebanir = (/0.820,0.820,0.820,0.820,0.820,0.820,0.751, &
       0.467,0.663,0.820,0.820,0.820,0.820,0.820/)

      !attempt to use spec albedoes had poor results
      !shebavis = (/0.826,0.826,0.826,0.826,0.826,0.762,0.499, &
      ! 0.681,0.719,0.826,0.826,0.826,0.826,0.826/)
      !shebanir = (/0.809,0.809,0.809,0.809,0.809,0.731,0.411, &
      ! 0.632,0.678,0.809,0.809,0.809,0.809,0.809/)
      !Rotate albedoes by 6 months for S Hemis. -not a good idea
      !shem = (/0.751,0.467,0.663,0.820,0.820,0.820,0.820, &
      ! 0.820,0.820,0.820,0.820,0.820,0.751,0.467/)

      nday = (/31.,31.,29.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.,31./)

      call ESMF_TimeGet  ( currTime, TimeString=string  ,_RC )
      read(string( 1: 4),'(i4.4)') YEAR
      read(string( 6: 7),'(i2.2)') MONTH
      read(string( 9:10),'(i2.2)') DAY

      if(mod(YEAR,4).eq.0) then
        nday(2)=29.
      else
        nday(2)=28.
      endif



      if(DAY.ge.15) then
        afracv=(float(DAY)-15.)/nday(MONTH)
        aslopev=(shebavis(MONTH+1)-shebavis(MONTH))/nday(MONTH)
        abasev=shebavis(MONTH)
        afraci=(float(DAY)-15.)/nday(MONTH)
        aslopei=(shebanir(MONTH+1)-shebanir(MONTH))/nday(MONTH)
        abasei=shebanir(MONTH)
      else
        afracv=(float(DAY)+nday(MONTH-1)-15.)/nday(MONTH-1)
        aslopev=(shebavis(MONTH)-shebavis(MONTH-1))/nday(MONTH-1)
        abasev=shebavis(MONTH-1)
        afraci=(float(DAY)+nday(MONTH-1)-15.)/nday(MONTH-1)
        aslopei=(shebanir(MONTH)-shebanir(MONTH-1))/nday(MONTH-1)
        abasei=shebanir(MONTH-1)
      endif


      SEAICEALBVRN=abasev+aslopev*afracv
      SEAICEALBVFN=abasev+aslopev*afracv
      SEAICEALBNRN=abasei+aslopei*afraci
      SEAICEALBNFN=abasei+aslopei*afraci

      SEAICEALBVRS=0.6
      SEAICEALBVFS=0.6
      SEAICEALBNRS=0.6
      SEAICEALBNFS=0.6

      where(LATS.ge.0.)
! Beam albedos
!-------------
        ALBVR = SEAICEALBVRN
        ALBNR = SEAICEALBNRN

! Diffuse albedos
!----------------
        ALBVF = SEAICEALBVFN
        ALBNF = SEAICEALBNFN
      elsewhere
! Beam albedos
!-------------
        ALBVR = SEAICEALBVRS
        ALBNR = SEAICEALBNRS

! Diffuse albedos
!----------------
        ALBVF = SEAICEALBVFS
        ALBNF = SEAICEALBNFS
      end where


   RETURN_(ESMF_SUCCESS)
  end subroutine ALBSEAICEM2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_CICE4ColumnPhysGridComp

subroutine CICEReOrder(Packed, UnPacked, MSK, Pdim, Udim, LM, DIR)
  integer, intent(IN   ) :: Pdim, Udim, LM, DIR
  real,    intent(INOUT) ::   Packed(Pdim,*)
  real,    intent(INOUT) :: UnPacked(Udim,*)
  logical, intent(IN   ) :: MSK(Udim)

  integer :: I, J, L, M

  do L = 1,LM
     M = 1
     do I = 1,Udim
        if (MSK(I)) then
           if(DIR==PACKIT) then
              Packed(M,L) = UnPacked(I,L)
           else
              Unpacked(I,L) = Packed(M,L)
           end if
           M = M+1
        else
           if(DIR/=PACKIT) then
              UnPacked(I,L) = 0
           end if
        end if
     end do
  end do
end subroutine CICEReOrder

subroutine CICEReOrder8(Packed, UnPacked, MSK, Pdim, Udim, LM, DIR)
  use MAPL, only : MAPL_R8
  integer, intent(IN   ) :: Pdim, Udim, LM, DIR
  real(kind=MAPL_R8),    intent(INOUT) ::   Packed(Pdim,*)
  real(kind=MAPL_R8),    intent(INOUT) :: UnPacked(Udim,*)
  logical, intent(IN   ) :: MSK(Udim)

  integer :: I, J, L, M

  do L = 1,LM
     M = 1
     do I = 1,Udim
        if (MSK(I)) then
           if(DIR==PACKIT) then
              Packed(M,L) = UnPacked(I,L)
           else
              Unpacked(I,L) = Packed(M,L)
           end if
           M = M+1
        else
           if(DIR/=PACKIT) then
              UnPacked(I,L) = 0
           end if
        end if
     end do
  end do
end subroutine CICEReOrder8
