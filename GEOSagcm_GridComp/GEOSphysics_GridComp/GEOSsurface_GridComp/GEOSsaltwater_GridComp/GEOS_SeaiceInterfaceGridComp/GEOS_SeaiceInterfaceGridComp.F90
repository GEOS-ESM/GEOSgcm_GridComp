
!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
module GEOS_SeaiceInterfaceGridComp

!BOP
! !MODULE: GEOS_SeaiceInterfaceGridComp -- Implements a simple sea ice interface on ice tiles.

! !DESCRIPTION:
! 
!   {\tt GEOS\_SeaiceInterface} is a light-weight gridded component that updates
!      the skin sub-tiles at saltwater points, be they ocean, estuary, or salt
!      lake. Currently each tile can have multiple subtiles representing ice categories,
!      and includes implementation of a simple sea ice interface layer  module
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
  

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

  integer, parameter :: ICE = 1         
  real,    parameter :: puny = 1.e-11  ! should come from CICE       

  logical ::      DUAL_OCEAN

  type seaice_interface_state
       integer                             :: CHOOSEMOSFC
       logical                             :: retrievedRootGC = .false.
       type (MAPL_LocStreamXform), pointer :: XFORM_A2O => NULL()
       type (MAPL_LocStreamXform), pointer :: XFORM_O2A => NULL()
       type (MAPL_LocStream), pointer      :: locStreamO => NULL()
  end type seaice_interface_state

  type seaice_interface_state_wrap
      type(seaice_interface_state), pointer :: ptr
  end type seaice_interface_state_wrap

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
    integer                                 :: NUM_ICE_CATEGORIES  ! set via resource parameter
    integer ::      iDUAL_OCEAN

    type(seaice_interface_state_wrap) :: wrap
    type(seaice_interface_state), pointer :: mystate
    character(len=ESMF_MAXSTR)     :: SURFRC
    type(ESMF_Config)              :: SCF

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)


    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run1, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,  Finalize, RC=STATUS )
    VERIFY_(STATUS)

! Get constants from CF
! ---------------------

    call ESMF_ConfigGetAttribute(CF, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" , RC=STATUS)
    VERIFY_(STATUS)
    NUM_SUBTILES  = NUM_ICE_CATEGORIES 

    call MAPL_GetResource(MAPL, iDUAL_OCEAN, 'DUAL_OCEAN:', default=0, RC=STATUS )
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
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_visible_beam',   &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVR',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_visible_diffuse',&
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVF',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_near_infrared_beam', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNR',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_near_infrared_diffuse', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNF',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
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
        LONG_NAME          = 'sublimation'               ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'SUBLIM'                    ,&
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
        LONG_NAME          = 'sea_ice_upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHICE'                     ,&
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

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'sea_ice_outgoing_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLWUPICE'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'sea_ice_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWNDICE'                   ,&
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
        LONG_NAME          = 'sea_ice_net_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWNDICE'                   ,&
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
        LONG_NAME          = 'sea_ice_latent_energy_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLATICE'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'total_surface_heat_flux_over_the_whole_tile' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'FSURF'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'total_surface_heat_flux_over_the_ice_tile' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'FSURFICE'                  ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TST',                               &
        LONG_NAME          = 'surface_skin_temperature',          &
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
        LONG_NAME          = 'turbulence_surface_temperature',    &
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
        SHORT_NAME         = 'UH',                                &
        LONG_NAME          = 'turbulence_surface_zonal_velocity', &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'VH',                                &
        LONG_NAME          = 'turbulence_surface_meridional_velocity', &
        UNITS              = 'm s-1',                             &
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
        SHORT_NAME         = 'RET',                               &
        LONG_NAME          = 'surface_reynolds_number',           &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRACI',                             &
        LONG_NAME          = 'ice_covered_fraction_of_tile',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'FRACINEW',                          &
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

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'eastward_stress_over_ice',  &
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUXI'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'northward_stress_over_ice',  &
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUYI'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PENUVR',                             &
        LONG_NAME          = 'penetrative_uvr_direct_flux_through_sea_ice',           &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PENUVF',                             &
        LONG_NAME          = 'penetrative_uvr_diffuse_flux_through_sea_ice',           &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PENPAR',                             &
        LONG_NAME          = 'penetrative_par_direct_flux_through_sea_ice',           &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PENPAF',                             &
        LONG_NAME          = 'penetrative_par_diffuse_flux_through_sea_ice',           &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'TFREEZE',                        &
        LONG_NAME          = 'freezing_temperature_for_interface_layer',&
        UNITS              = 'K',                               &
        DIMS               = MAPL_DimsTileOnly,                 &
        VLOCATION          = MAPL_VLocationNone,                &
        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'surface_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWDNSRF'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'surface_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWDNSRF'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                          &
        SHORT_NAME         = 'GHTSKIN',                   &
        LONG_NAME          = 'Ground_heating_for_skin_temp',&
        UNITS              = 'W m-2',                     &
        DIMS               = MAPL_DimsTileOnly,           &
        VLOCATION          = MAPL_VLocationNone,          &
                                                     _RC  )

!  !INTERNAL STATE:

     call MAPL_AddInternalSpec(GC,                                &
         SHORT_NAME         = 'TSKINI',                            &
         LONG_NAME          = 'ice_skin_temperature',              &
         UNITS              = 'K',                                 &
         UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &    
         DIMS               = MAPL_DimsTileOnly,                   &   
         VLOCATION          = MAPL_VLocationNone,                  &
         FRIENDLYTO         = 'SEAICE',                            &
         DEFAULT            = MAPL_TICE-1.8,                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'SSKINI',                            &
        LONG_NAME          = 'ice_skin_salinity',                 &
        UNITS              = 'psu',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = 30.0,                                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'QS',                                &
        LONG_NAME          = 'surface_specific_humidity',         &
        UNITS              = 'kg kg-1',                           &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.01,                                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CH',                                &
        LONG_NAME          = 'surface_heat_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
        DEFAULT            = 1.0e-4,                              &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CM',                                &
        LONG_NAME          = 'surface_momentum_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
        DEFAULT            = 1.0e-4,                              &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CQ',                                &
        LONG_NAME          = 'surface_moisture_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
        DEFAULT            = 1.0e-4,                              &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'Z0',                                &
        LONG_NAME          = 'aerodynamic_roughness',             &
        UNITS              = 'm',                                 &
        DEFAULT            = 0.00005,                             &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'WW',                                &
        LONG_NAME          = 'vertical_velocity_scale_squared',   &
        UNITS              = 'm+2 s-2',                           &
        DEFAULT            = 0.0,                                 &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

!  !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'ALW',                               &
        LONG_NAME          = 'linearization_of_surface_upwelling_longwave_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'BLW',                               &
        LONG_NAME          = 'linearization_of_surface_upwelling_longwave_flux', &
        UNITS              = 'W m-2 K-1',                         &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'LWDNSRF',                           &
        LONG_NAME          = 'surface_downwelling_longwave_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                             ,&
        LONG_NAME          = 'surface_downwelling_par_beam_flux' ,&
        UNITS              = 'W m-2'                             ,&
        SHORT_NAME         = 'DRPAR'                             ,&
        DIMS               = MAPL_DimsTileOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
        RESTART            = MAPL_RestartSkip                    ,&
                                                       RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_par_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DFPAR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
         RESTART            = MAPL_RestartSkip              ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_nir_beam_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DRNIR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
         RESTART            = MAPL_RestartSkip              ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_nir_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DFNIR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
         RESTART            = MAPL_RestartSkip              ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_uvr_beam_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DRUVR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
         RESTART            = MAPL_RestartSkip              ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_uvr_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DFUVR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
         RESTART            = MAPL_RestartSkip              ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'evaporation',                       &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'EVAP ',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'upward_sensible_heat_flux',         &
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'SH',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'eastward_surface_stress',           &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUX',                              &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'northward_surface_stress',          &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUY',                              &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_evaporation',         &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'DEVAP',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_upward_sensible_heat_flux', &
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'DSH',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'snowfall',                          &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'SNO',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

! Surface air quantities

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_air_temperature',           &
        UNITS              = 'K',                                 &
        SHORT_NAME         = 'TA',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_air_specific_humidity',     &
        UNITS              = 'kg kg-1',                           &
        SHORT_NAME         = 'QA',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_wind_speed',                &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'UU',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'levellm_uwind',                     &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'UWINDLMTILE',                       &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'levellm_vwind',                     &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'VWINDLMTILE',                       &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_layer_height',              &
        UNITS              = 'm',                                 &
        SHORT_NAME         = 'DZ',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_pressure',                  &
        UNITS              = 'Pa',                                &
        SHORT_NAME         = 'PS',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'liquid_water_convective_precipitation',&
        UNITS              = 'kg m-2 s-1'                        ,&
        SHORT_NAME         = 'PCU'                               ,&
        DIMS               = MAPL_DimsTileOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
        RESTART            = MAPL_RestartSkip                    ,&
                                                       RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                            ,&
        LONG_NAME          = 'liquid_water_large_scale_precipitation',&
        UNITS              = 'kg m-2 s-1'                       ,&
        SHORT_NAME         = 'PLS'                              ,&
        DIMS               = MAPL_DimsTileOnly                  ,&
        VLOCATION          = MAPL_VLocationNone                 ,&
        RESTART            = MAPL_RestartSkip                   ,&
                                                      RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'THATM',                             &
        LONG_NAME          = 'effective_surface_skin_temperature',&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'QHATM',                             &
        LONG_NAME          = 'effective_surface_specific_humidity',&
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'UHATM',                             &
        LONG_NAME          = 'effective_surface_zonal_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VHATM',                             &
        LONG_NAME          = 'effective_surface_meridional_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CTATM',                             &
        LONG_NAME          = 'surface_exchange_coefficient_for_heat', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CQATM',                             &
        LONG_NAME          = 'surface_exchange_coefficient_for_moisture', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CMATM',                             &
        LONG_NAME          = 'surface_exchange_coefficient_for_momentum', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'FRACICE',                           &
        LONG_NAME          = 'ice_covered_fraction_of_tile',      &
        UNITS              = '1',                                 &
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'UI',                                &
        LONG_NAME          = 'zonal_velocity_of_surface_ice',     &
        UNITS              = 'm s-1 ',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VI',                                &
        LONG_NAME          = 'meridional_velocity_of_surface_ice',     &
        UNITS              = 'm s-1 ',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


   call MAPL_AddImportSpec(GC,                                    &
        SHORT_NAME         = 'UUA',                               &
        LONG_NAME          = 'interpolated_effective_surface_zonal_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                           __RC__ )

   call MAPL_AddImportSpec(GC,                                    &
        SHORT_NAME         = 'VVA',                               &
        LONG_NAME          = 'interpolated_effective_surface_meridional_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                           __RC__ )

   call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'SS_FOUND',                          &
        LONG_NAME          = 'foundation_salinity_for_interface_layer',               &
        UNITS              = 'psu',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 33.3333,                             &  !corresponding to -1.8 deg c

                                                       _RC  )


   !*CALLBACK*
   ! an ESMF state to pass information b.w. GCs using callback
   call MAPL_AddImportSpec(GC                                    ,&
          SHORT_NAME         = 'SURFSTATE'                       ,&
          LONG_NAME          = 'surface_state_for_seaice_thermo_coupling',  &
          UNITS              = 'W m-2'                           ,&
          !UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),       &
          !DIMS               = MAPL_DimsTileOnly,           &
          !VLOCATION          = MAPL_VLocationNone,          &
          DATATYPE           = MAPL_StateItem                    ,&
          RESTART            = MAPL_RestartSkip                  ,&
                                                           __RC__ )



  call MAPL_AddExportSpec(GC,                    &
          SHORT_NAME         = 'FCONDTOP'                  ,             &
          LONG_NAME          = 'conductive_heat_flux_at_ice_top_surface',&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )

  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC                               ,&
         SHORT_NAME         = 'TSKINICE'                   ,&
         LONG_NAME          = 'snow_or_ice_surface_temperature',&
         UNITS              = 'K'                          ,&
         DIMS               = MAPL_DimsTileOnly            ,&
         VLOCATION          = MAPL_VLocationNone           ,&
                                                     __RC__ )


  call MAPL_AddExportSpec(GC                               ,&
         SHORT_NAME          = 'SUBLIMFLX'                 ,&
          LONG_NAME          = 'heat_flux_associated_with_sublimation_of_snow_ice',&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
                                                     __RC__ )

!  Category dimensional exports


   call MAPL_AddExportSpec(GC,                    &                  
         SHORT_NAME         = 'FCONDTOPN'                 ,                            &
          LONG_NAME          = 'conductive_heat_flux_at_ice_snow_surface_over_ice_categories',&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'FSURFN'                    ,&
        LONG_NAME          = 'net_heat_flux_at_ice_snow_surface_over_ice_categories' ,&
        UNITS              = 'W m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'SHICEN'                    ,&
        LONG_NAME          = 'sea_ice_upward_sensible_heat_flux_over_ice_categories' ,&
        UNITS              = 'W m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'HLWUPN'                     ,&
        LONG_NAME          = 'outgoing_longwave_flux_at_ice_snow_surface_over_ice_categories',&
        UNITS              = 'W m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'LWNDSRFN'                     ,&
        LONG_NAME          = 'net_downward_longwave_flux_at_ice_snow_surface_over_ice_categories',&
        UNITS              = 'W m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'FSWSFCN'                    ,&
        LONG_NAME          = 'SW_absorbed_at_ice_snow_surface_over_ice_categories' ,&
        UNITS              = 'W m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'TSURFN'                    ,&
        LONG_NAME          = 'ice_snow_surface_temperature_over_ice_categories' ,&
        UNITS              = 'K'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'ALBIN'                    ,&
        LONG_NAME          = 'ice_surface_albedo_over_ice_categories' ,&
        UNITS              = '1'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                     &
        SHORT_NAME         = 'ALBSN'                    ,&
        LONG_NAME          = 'snow_surface_albedo_over_ice_categories' ,&
        UNITS              = '1'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/)      ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
   VERIFY_(STATUS)


!EOS

    allocate(mystate,stat=status)
    VERIFY_(status)
    call MAPL_GetResource (MAPL, SURFRC, label = 'SURFRC:', default = 'GEOS_SurfaceGridComp.rc', RC=STATUS) ; VERIFY_(STATUS)
    SCF = ESMF_ConfigCreate(rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigLoadFile     (SCF,SURFRC,rc=status) ; VERIFY_(STATUS)
    call MAPL_GetResource (SCF, mystate%CHOOSEMOSFC, label='CHOOSEMOSFC:', DEFAULT=1, __RC__ )
    call ESMF_ConfigDestroy      (SCF, __RC__)
    wrap%ptr => mystate
    call ESMF_UserCompSetInternalState(gc, 'seaice_interface_private', wrap,status)
    VERIFY_(status)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="INITIALIZE"   ,         RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerAdd(GC,    name="RUN1"   ,               RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN2"  ,                RC=STATUS)
    VERIFY_(STATUS)
  
    call MAPL_TimerAdd(GC,    name="-Thermo1"    ,          RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerAdd(GC,    name="FINALIZE"   ,         RC=STATUS)
    VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC,  RC=STATUS )
    VERIFY_(STATUS)
 
! Set the Run entry point
! -----------------------

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Initialize -- Initialize method for the GEOS Seaice Interface Grid Comp

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm 
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

    integer                             :: NUM_SUBTILES        ! = NUM_ICE_CATEGORIES 
    integer                             :: NUM_ICE_CATEGORIES  ! set via resource parameter
    
! Local derived type aliases

    type (MAPL_MetaComp    ), pointer   :: MAPL => null()

    real                                :: DTI

    integer                             :: PRES_ICE

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" ,     RC=STATUS)
    VERIFY_(STATUS)
    NUM_SUBTILES  = NUM_ICE_CATEGORIES 

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"INITIALIZE")

    call MAPL_Get(MAPL, HEARTBEAT = DTI, RC=STATUS)
    VERIFY_(STATUS)

! Call Initialize for every Child
!--------------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)
 
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
   real, pointer, dimension(:  )  :: FRACIN => null()


! pointers to internal

   real, pointer, dimension(:,:)  :: TI  => null()      ! ice skin temperature  
   real, pointer, dimension(:,:)  :: QS  => null()
   real, pointer, dimension(:,:)  :: CH  => null()
   real, pointer, dimension(:,:)  :: CM  => null()
   real, pointer, dimension(:,:)  :: CQ  => null()
   real, pointer, dimension(:,:)  :: WW  => null()
   real, pointer, dimension(:,:)  :: Z0  => null()

! pointers to import

   real, pointer, dimension(:,:)  :: FR  => null()  
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
   real, pointer, dimension(:)    :: LATS => null()
   real, pointer, dimension(:)    :: LONS => null()

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
   real, allocatable              :: FRI(:)


   integer                        :: K
   logical                        :: L_STOP
   integer                        :: IDUM, JDUM
   real                           :: DT


   real            :: ICEZ0
   real, parameter :: HPBL = 1000.

   integer         :: PRES_ICE
   integer         :: CHOOSEMOSFC
   integer         :: CHOOSEZ0
   type(seaice_interface_state_wrap) :: wrap
   type(seaice_interface_state), pointer :: mystate

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run1"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" ,     RC=STATUS)
    VERIFY_(STATUS)
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
         TILELATS = LATS,         &
         TILELONS = LONS,         &
         RC=STATUS )
    VERIFY_(STATUS)

! Get parameters (0:Louis, 1:Monin-Obukhov)
! -----------------------------------------
    call ESMF_UserCompGetInternalState(gc,'seaice_interface_private',wrap,status)
    VERIFY_(status)
    mystate => wrap%ptr
    CHOOSEMOSFC = mystate%CHOOSEMOSFC

    call MAPL_GetResource ( MAPL, CHOOSEZ0,    Label="CHOOSEZ0:",    DEFAULT=3, RC=STATUS)
    VERIFY_(STATUS)

! Get roughness parameters with and without CICE Thermodynamics
! -------------------------------------------------------------
! icez0 value is based on literature (Ask Bin). It could be revisited, later.
    call MAPL_GetResource ( MAPL, ICEZ0,       Label="ICEZ0:" ,          DEFAULT=5.0e-4, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, PRES_ICE,    Label="PRESCRIBED_ICE:" , DEFAULT=1,      RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_Get(MAPL, HEARTBEAT = DT, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, DT, Label="DT:", DEFAULT=DT, RC=STATUS)
    VERIFY_(STATUS)

! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,UU     , 'UU'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,UWINDLMTILE     , 'UWINDLMTILE'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VWINDLMTILE     , 'VWINDLMTILE'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,UI     , 'UI'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VI     , 'VI'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DZ     , 'DZ'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TA     , 'TA'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,QA     , 'QA'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PS     , 'PS'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PCU    , 'PCU'    ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,FR   , 'FRACICE'  ,    RC=STATUS)
   VERIFY_(STATUS) 
   call MAPL_GetPointer(IMPORT,SW   , 'SS_FOUND' ,    RC=STATUS)
   VERIFY_(STATUS) 


! Pointers to internals
!----------------------

   !call MAPL_GetPointer(INTERNAL,FR   , 'FR'     ,    RC=STATUS)
   !VERIFY_(STATUS) 
   call MAPL_GetPointer(INTERNAL,TI   , 'TSKINI' ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,QS   , 'QS'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CH   , 'CH'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CM   , 'CM'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CQ   , 'CQ'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,Z0   , 'Z0'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,WW   , 'WW'     ,    RC=STATUS)
   VERIFY_(STATUS)

! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,QH    , 'QH'      ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TH    , 'TH'      ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,UH    , 'UH'      ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,VH    , 'VH'      ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,QST   , 'QST'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TST   , 'TST'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CHT   , 'CHT'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CMT   , 'CMT'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CQT   , 'CQT'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CNT   , 'CNT'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RIT   , 'RIT'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RET   , 'RET'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Z0O   , 'Z0'      ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Z0H   , 'Z0H'     ,    RC=STATUS)
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
   call MAPL_GetPointer(EXPORT,GST   , 'GUST'    ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,VNT   , 'VENT'    ,    RC=STATUS)
   VERIFY_(STATUS)

  ! export to openwater
   call MAPL_GetPointer(EXPORT,FRACI , 'FRACI'   ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TF    , 'TFREEZE' ,    _RC)

   NT = size(TA)
  ! if(NT == 0) then
  !    call MAPL_TimerOff(MAPL,"RUN1" )
  !    call MAPL_TimerOff(MAPL,"TOTAL")
  !    RETURN_(ESMF_SUCCESS)
  ! end if

   allocate(RE (NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CN (NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(ZT (NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(T2M (NT)  ,  STAT=STATUS)
   VERIFY_(STATUS)
   allocate(Q2M (NT)  ,  STAT=STATUS)
   VERIFY_(STATUS)
   allocate(U2M (NT)  ,  STAT=STATUS)
   VERIFY_(STATUS)
   allocate(V2M (NT)  ,  STAT=STATUS)
   VERIFY_(STATUS)
   allocate(T10M (NT)  , STAT=STATUS)
   VERIFY_(STATUS)
   allocate(Q10M (NT)  , STAT=STATUS)
   VERIFY_(STATUS)
   allocate(U10M (NT)  , STAT=STATUS)
   VERIFY_(STATUS)
   allocate(V10M (NT)  , STAT=STATUS)
   VERIFY_(STATUS)
   allocate(U50M (NT)  , STAT=STATUS)
   VERIFY_(STATUS)
   allocate(V50M (NT)  , STAT=STATUS)
   VERIFY_(STATUS)
   allocate(ZQ (NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(UUU(NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(RHO(NT) ,    STAT=STATUS)
   VERIFY_(STATUS)
   allocate(PSMB(NT) ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(PSL(NT) ,    STAT=STATUS)
   VERIFY_(STATUS)
   allocate(VKH(NT) ,    STAT=STATUS)
   VERIFY_(STATUS)
   allocate(fakelai(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(VKM(NT) ,    STAT=STATUS)
   VERIFY_(STATUS)
   allocate(USTAR(NT) ,  STAT=STATUS)
   VERIFY_(STATUS)
   allocate(XX(NT)   ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(YY(NT)   ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CU(NT)   ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CT(NT)   ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(RIB(NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(ZETA(NT) ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(WS(NT)   ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(IWATER(NT),  STAT=STATUS)
   VERIFY_(STATUS)
   allocate(LAI(NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CHB(NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CQB(NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CMB(NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(UCN(NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(US (NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(VS (NT)  ,   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(TS (NT,NUM_SUBTILES),   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(FRI (NT) ,   STAT=STATUS)
   VERIFY_(STATUS)


   TS(:,ICE:) = TI           ! TI(in K): returned from Sea Ice Model Update
   US = UI
   VS = VI

   ! refresh QS based on the updated TS due to sea ice dynamic transport
   do N=1,NUM_SUBTILES
     QS(:,N) = GEOS_QSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.) 
   enddo

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

   if(associated(FRACI )) FRACI = min(FRI, 1.0)
   if(associated(TF )) TF = -0.054*SW+MAPL_TICE


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
  integer                             :: NUM_ICE_CATEGORIES  ! set via resource parameter

  real, pointer, dimension(:)         :: LATS => null()
  real, pointer, dimension(:)         :: LONS => null()

  real, pointer, dimension(:)         :: AREA => null()     ! needed to calculate TILEAREA in SaltWaterCore

  type(seaice_interface_state_wrap)     :: wrap
  type(seaice_interface_state), pointer :: mystate


!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run2"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    call ESMF_UserCompGetInternalState(GC,'seaice_interface_private',wrap,status)
    VERIFY_(status)
    mystate => wrap%ptr


! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" ,     RC=STATUS)
    VERIFY_(STATUS)
    NUM_SUBTILES  = NUM_ICE_CATEGORIES 

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN2" )

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,             &
         TILELATS  = LATS ,                      &
         TILELONS  = LONS ,                      &
         TILEAREA  = AREA ,                      &
         ORBIT     = ORBIT,                      &
         INTERNAL_ESMF_STATE = INTERNAL,         &
         CF = CF,                                &
                                       RC=STATUS )
    VERIFY_(STATUS)

! Update the skin variables each step
!------------------------------------

    call CICECORE(NT=size(LONS), RC=STATUS )
    VERIFY_(STATUS)

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN2" )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine CICECORE(NT,RC)

   integer,           intent(IN ) :: NT
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
   real, pointer, dimension(:  )  :: LWDNSRFe=> null()
   real, pointer, dimension(:  )  :: SWDNSRFe=> null()

   real, pointer, dimension(:  )  :: TSKINICE=> null()

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
   real, pointer, dimension(:  )  :: FRACIN => null()
   real, pointer, dimension(:  )  :: GHTSKIN => null()


! pointers to internal

   real, pointer, dimension(:,:)  :: TI    => null()
   real, pointer, dimension(:  )  :: SI    => null()
   real, pointer, dimension(:,:)  :: QS    => null()
   real, pointer, dimension(:,:)  :: CH    => null()
   real, pointer, dimension(:,:)  :: CQ    => null()
   real, pointer, dimension(:,:)  :: CM    => null()


! pointers to import
   real, pointer, dimension(:,:)  :: FR => null()

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
   real, pointer, dimension(:)    :: UI => null()
   real, pointer, dimension(:)    :: VI => null()

   real, pointer, dimension(:,:)  :: TS        => null()

   !*CALLBACK*
   real, pointer, dimension(:,:)  :: AS_PTR_2D => null()
   integer                        :: AS_STATUS
   type(ESMF_State)               :: SURFST
   type(MAPL_LocStream)           :: locStreamO
   type(MAPL_LocStreamXform)      :: xform_A2O
   type(MAPL_LocStreamXform)      :: xform_O2A


   real,    dimension(NT)              :: SHF
   real,    dimension(NT)              :: EVP
   real,    dimension(NT)              :: SHD
   real,    dimension(NT)              :: EVD
   real,    dimension(NT)              :: CFQ
   real,    dimension(NT)              :: CFT
   !real,    dimension(NT)              :: UUA
   !real,    dimension(NT)              :: VVA
   real,    dimension(NT)              :: TXI
   real,    dimension(NT)              :: TYI
   real,    dimension(NT)              :: DQS
   real,    dimension(NT)              :: DTS
   real,    dimension(NT)              :: SWN
   real,    dimension(NT)              :: PEN
   real,    dimension(NT)              :: LHF
   real,    dimension(NT)              :: ZTH
   real,    dimension(NT)              :: SLR
   real,    dimension(NT)              :: VSUVR
   real,    dimension(NT)              :: VSUVF
   real,    dimension(NT)              :: FRCICE

   integer                             :: N
   real                                :: DT
   real                                :: MAXSALINITY
   real                                :: MINSALINITY


   integer                             :: NSUB, I, K, L

   real                                :: LATSO, LONSO
   real,    dimension(1)               :: LATSD, LONSD
   logical, dimension(1)               :: OBSERVE


   real                                :: YDAY 
   real,    dimension(NT)              :: ALBVRI
   real,    dimension(NT)              :: ALBVFI
   real,    dimension(NT)              :: ALBNRI
   real,    dimension(NT)              :: ALBNFI

   real,               allocatable    :: ALBVRN        (:,:)
   real,               allocatable    :: ALBNRN        (:,:)
   real,               allocatable    :: ALBVFN        (:,:)
   real,               allocatable    :: ALBNFN        (:,:)

   real,               allocatable    :: TS_OLD        (:,:)

   real,               allocatable    :: FSURF         (:,:) ! non-solar part 
   real,               allocatable    :: DFSURFDTS     (:,:) ! 
   real,               allocatable    :: DLHFDTS       (:,:) ! 
   real,               allocatable    :: EVAPN         (:,:) ! 
   real,               allocatable    :: LHFN          (:,:) ! 
   real,               allocatable    :: RAIN          (:)   !



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
   logical                             :: debugzth

   real                                :: EMSICE

!  Begin...
!----------

   IAm =  trim(COMP_NAME) // "CICECORE"


! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,ALW    , 'ALW'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,BLW    , 'BLW'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,LWDNSRF, 'LWDNSRF',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DRPAR  , 'DRPAR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DFPAR  , 'DFPAR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DRNIR  , 'DRNIR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DFNIR  , 'DFNIR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DRUVR  , 'DRUVR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DFUVR  , 'DFUVR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,EVAP   , 'EVAP'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SH     , 'SH'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TAUX   , 'TAUX'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TAUY   , 'TAUY'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DEV    , 'DEVAP'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DSH    , 'DSH'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SNO    , 'SNO'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PLS    , 'PLS'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PCU    , 'PCU'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PS     , 'PS'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,UU     , 'UU'     ,    RC=STATUS); VERIFY_(STATUS)

   call MAPL_GetPointer(IMPORT,UI     , 'UI'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VI     , 'VI'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,THATM  , 'THATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,QHATM  , 'QHATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,UHATM  , 'UHATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VHATM  , 'VHATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,UUA    , 'UUA'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VVA    , 'VVA'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,CTATM  , 'CTATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,CQATM  , 'CQATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,CMATM  , 'CMATM'  ,    RC=STATUS); VERIFY_(STATUS)

   call MAPL_GetPointer(IMPORT,FR     , 'FRACICE' ,    RC=STATUS); VERIFY_(STATUS)

! Pointers to internals
!----------------------

   call MAPL_GetPointer(INTERNAL,TI     ,'TSKINI',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,SI     ,'SSKINI',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,QS     , 'QS'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CH     , 'CH'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CQ     , 'CQ'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CM     , 'CM'   ,    RC=STATUS); VERIFY_(STATUS)


! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,EMISS  , 'EMIS' , alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBVF  , 'ALBVF', alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBVR  , 'ALBVR', alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBNF  , 'ALBNF', alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBNR  , 'ALBNR', alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,QST    , 'QST'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TST    , 'TST'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,DELTS  , 'DELTS'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,DELQS  , 'DELQS'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TAUXI  , 'TAUXI'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TAUYI  , 'TAUYI'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,PENUVR , 'PENUVR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,PENUVF , 'PENUVF'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,PENPAR , 'PENPAR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,PENPAF , 'PENPAF'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,EVAPOUT, 'EVAPOUT' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SUBLIM,  'SUBLIM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SHOUT  , 'SHOUT'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SHICE  , 'SHICE'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLATN  , 'HLATN'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLATICE, 'HLATICE' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FSURFe , 'FSURF'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FSURFICE,'FSURFICE',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLWUP  , 'HLWUP'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLWUPe , 'HLWUPICE',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,LWNDSRF, 'LWNDSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SWNDSRF, 'SWNDSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,LWNDICE, 'LWNDICE' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SWNDICE, 'SWNDICE' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FRACI  , 'FRACI'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FRACIN , 'FRACINEW',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,LWDNSRFe,'LWDNSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SWDNSRFe,'SWDNSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TSKINICE,'TSKINICE',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,GHTSKIN, 'GHTSKIN' ,    RC=STATUS); VERIFY_(STATUS)


! Get the time step
! -----------------

    call MAPL_Get(MAPL, HEARTBEAT = DT, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, DT, Label="DT:", DEFAULT=DT, RC=STATUS)
    VERIFY_(STATUS)

! Get parameters
! --------------

    call MAPL_GetResource ( MAPL, EMSICE,      Label="CICE_EMSICE:",   DEFAULT=0.99999, RC=STATUS)
    VERIFY_(STATUS)

! Copy friendly internals into tile-tile local variables
!-------------------------------------------------------

    TS => TI



    allocate(    FSURF(size(TS,1),   size(TS,2)), __STAT__) 
    allocate(DFSURFDTS(size(TS,1),   size(TS,2)), __STAT__) 
    allocate(  DLHFDTS(size(TS,1),   size(TS,2)), __STAT__) 
    allocate(    EVAPN(size(TS,1),   size(TS,2)), __STAT__) 
    allocate(     LHFN(size(TS,1),   size(TS,2)), __STAT__) 
    allocate(   TS_OLD(size(TS,1),   size(TS,2)), __STAT__) 

    allocate(    RAIN(size(TS,1)),  __STAT__) 


! Aggregate imports if required
!-----------------------------------

    RAIN = PLS + PCU

! Initialize PAR and UVR beam fluxes
!-----------------------------------

    VSUVR = DRPAR + DRUVR
    VSUVF = DFPAR + DFUVR

    if(associated(SWDNSRFe)) SWDNSRFe = VSUVR+VSUVF+DRNIR+DFNIR 
    if(associated(LWDNSRFe)) LWDNSRFe = LWDNSRF 


    call MAPL_GetResource ( MAPL, LATSO, Label="LATSO:", DEFAULT=70.0, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, LONSO, Label="LONSO:", DEFAULT=70.0, RC=STATUS)
    VERIFY_(STATUS)


    TS_OLD     = TS


    call MAPL_TimerOn(MAPL,    "-Albedo")

    debugzth = .false.

    call ESMF_VMGetCurrent ( VM, RC=STATUS )

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


! Albedo over Sea-Ice. With LANL CICE, it is based on current ice states, 
! also compute shortwave radiation passing thru bottom of ice and skin layer bottom (DxxxTHRU; xx=RUVR, FUVR, ...)
!------------------------------------------------------------------------------------------------------------------
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
    if(associated(LWNDICE)) LWNDICE = 0.0 


    FRCICE = sum(FR(:,ICE:), dim=2) 

! Atmospheric surface stresses
!-----------------------------
    !*** already computed in surface
    !UUA = (TAUX/CMATM + UHATM)
    !VVA = (TAUY/CMATM + VHATM)


! Stress over ice
!----------------

    where(FRCICE > puny)
        TXI = sum(CM(:,ICE:)*FR(:,ICE:), dim=2)/FRCICE*(UUA - UI)
        TYI = sum(CM(:,ICE:)*FR(:,ICE:), dim=2)/FRCICE*(VVA - VI)
    elsewhere
        TXI = 0.0
        TYI = 0.0
    endwhere

    if(associated(TAUXI)) TAUXI = TXI
    if(associated(TAUYI)) TAUYI = TYI


    call MAPL_TimerOn(MAPL,   "-Thermo1")

! ------------------------------------

    categories_th1_: do N=ICE, NUM_SUBTILES   ! Loop over ice catgories. 

          CFT = (CH(:,N)/CTATM)
          CFQ = (CQ(:,N)/CQATM)
          EVP = CFQ*(EVAP + DEV*(QS(:,N)-QHATM))
          SHF = CFT*(SH   + DSH*(TS(:,N)-THATM))
          SHD = CFT*DSH
          EVD = CFQ*DEV*GEOS_DQSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.)
          LHF = EVP * MAPL_ALHS

          FSURF(:,N)      = LWDNSRF - SHF - LHF - (ALW + BLW*TS(:,N))
          DFSURFDTS(:,N)  = -(SHD + EVD * MAPL_ALHS + BLW) !!! check the sign
          DLHFDTS(:,N)    = -EVD * MAPL_ALHS  !!! check the sign
          EVAPN(:,N)      = -EVP !!! check the sign
          LHFN(:,N)       = -LHF !!! check the sign

    enddo categories_th1_


    !*CALLBACK*
    !==============================================================================================

    ! retrieve the regridding from the "private" wrapped state
    call RetrieveTransforms(mapl, mystate, XFORM_A2O, XFORM_O2A, locStreamO, __RC__) 

    call ESMF_StateGet(IMPORT, 'SURFSTATE', SURFST, __RC__)

    !call RegridA2O_2d(       TS, SURFST,  'TSKINICE', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_2d(    FSURF, SURFST,     'FSURF', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_2d(DFSURFDTS, SURFST, 'DFSURFDTS', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_2d(  DLHFDTS, SURFST,   'DLHFDTS', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_2d(     LHFN, SURFST,       'LHF', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_2d(    EVAPN, SURFST,      'EVAP', XFORM_A2O, locstreamO, __RC__)

    call RegridA2O_1d(     RAIN, SURFST,      'RAIN', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_1d(      SNO, SURFST,      'SNOW', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_1d(      ZTH, SURFST,      'COSZ', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_1d(    DRPAR, SURFST,     'DRPAR', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_1d(    DFPAR, SURFST,     'DFPAR', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_1d(    DRNIR, SURFST,     'DRNIR', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_1d(    DFNIR, SURFST,     'DFNIR', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_1d(    DRUVR, SURFST,     'DRUVR', XFORM_A2O, locstreamO, __RC__)
    call RegridA2O_1d(    DFUVR, SURFST,     'DFUVR', XFORM_A2O, locstreamO, __RC__)

    ! ** execute the sea ice thermo coupling method
    call ESMF_MethodExecute(SURFST, label="thermo_coupling", userRC=AS_STATUS, RC=STATUS)
    VERIFY_(AS_STATUS)
    VERIFY_(STATUS)

    allocate(AS_PTR_2D(size(TS,1),size(TS,2)), __STAT__)


    call RegridO2A_2d(AS_PTR_2D, SURFST, 'DTS', &
         XFORM_O2A, locstreamO, __RC__)

    update_surf: do N=ICE, NUM_SUBTILES   ! Loop over ice catgories. 
          CFT     = (CH(:,N)/CTATM)
          CFQ     = (CQ(:,N)/CQATM)
          EVP     = CFQ*(EVAP + DEV*(QS(:,N)-QHATM))
          SHF     = CFT*(SH   + DSH*(TS(:,N)-THATM))
          SHD     = CFT*DSH
          EVD     = CFQ*DEV*GEOS_DQSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.)
          LHF     = EVP * MAPL_ALHS
!         Aggregate ts and qs change over ice categories
          DTS     = AS_PTR_2D(:,N)           
          TS(:,N) = TS(:,N) + DTS
          DQS     = GEOS_QSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.) - QS(:,N)
          QS(:,N) = QS(:,N) + DQS

          LHF     = LHF + EVD * MAPL_ALHS * DTS
          SHF     = SHF + SHD * DTS

          if(associated(DELTS  )) DELTS   = DELTS   + DTS*CFT*FR(:,N)
          if(associated(DELQS  )) DELQS   = DELQS   + DQS*CFQ*FR(:,N)
          if(associated(TST    )) TST     = TST     + TS(:,N)*FR(:,N)
          if(associated(QST    )) QST     = QST     + QS(:,N)*FR(:,N)
          if(associated(HLATICE)) HLATICE = HLATICE + LHF    *FR(:,N)
          if(associated(SHICE  )) SHICE   = SHICE   + SHF    *FR(:,N)
          if(associated(LWNDICE)) LWNDICE = LWNDICE + (LWDNSRF - ALW - BLW*TS(:,N))*FR(:,N)
    end do update_surf

    EMISS = EMSICE

    if(associated(DELTS  )) call Normalize(DELTS,  FRCICE) 
    if(associated(DELQS  )) call Normalize(DELQS,  FRCICE) 
    if(associated(TST    )) call Normalize(TST,    FRCICE) 
    if(associated(QST    )) call Normalize(QST,    FRCICE) 
    if(associated(HLATICE)) call Normalize(HLATICE,FRCICE)
    if(associated(SHICE  )) call Normalize(SHICE,  FRCICE)

    if(associated(LWNDICE)) call Normalize(LWNDICE,  FRCICE, set_undef=.True.)
          

    call RegridO2A_1d(ALBVR, SURFST, 'ALBVR', &
         XFORM_O2A, locstreamO, __RC__)
    call RegridO2A_1d(ALBVF, SURFST, 'ALBVF', &
         XFORM_O2A, locstreamO, __RC__)
    call RegridO2A_1d(ALBNR, SURFST, 'ALBNR', &
         XFORM_O2A, locstreamO, __RC__)
    call RegridO2A_1d(ALBNF, SURFST, 'ALBNF', &
         XFORM_O2A, locstreamO, __RC__)

    if(associated(PENUVR )) then
        call RegridO2A_1d(PENUVR, SURFST, 'PENUVR', &
             XFORM_O2A, locstreamO, __RC__)
    endif
    if(associated(PENUVF )) then
        call RegridO2A_1d(PENUVF, SURFST, 'PENUVF', &
             XFORM_O2A, locstreamO, __RC__)
    endif
    if(associated(PENPAR )) then
        call RegridO2A_1d(PENPAR, SURFST, 'PENPAR', &
             XFORM_O2A, locstreamO, __RC__)
    endif
    if(associated(PENPAF )) then
        call RegridO2A_1d(PENPAF, SURFST, 'PENPAF', &
             XFORM_O2A, locstreamO, __RC__)
    endif
    !!!*** need to check the sign of GHTSKIN
    if(associated(GHTSKIN )) then
        call RegridO2A_1d(GHTSKIN, SURFST, 'GHTSKIN', &
             XFORM_O2A, locstreamO, __RC__)
    endif


    !************************************************************************************************
    !
    !      surface flux and temperature updated by sea ice
    !      continue updating relevant fields and pass them to surf
    !    .....
    !
    deallocate(AS_PTR_2D)



    call MAPL_TimerOff(MAPL,  "-Thermo1")


    if(associated(FRACI))   FRACI = min(FRCICE, 1.0)
    if(associated(FRACIN)) FRACIN = min(FRCICE, 1.0)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    call MAPL_TimerOn(MAPL,    "-Albedo")

    if(solalarmison) then
       call MAPL_SunGetInsolation(LONS, LATS,      &
            ORBIT, ZTH, SLR,                       &
            INTV = TINT,                           &
            currTime=CURRENT_TIME+DELT,            &
            RC=STATUS )
       VERIFY_(STATUS)

       ZTH = max(0.0,ZTH)

       call RegridA2O_1d(      ZTH, SURFST,      'COSZ', XFORM_A2O, locstreamO, __RC__)

       call ESMF_MethodExecute(SURFST, label="prep_albedo", userRC=AS_STATUS, RC=STATUS)
       VERIFY_(AS_STATUS)
       VERIFY_(STATUS)

       call RegridO2A_1d(ALBVR, SURFST, 'ALBVR', &
             XFORM_O2A, locstreamO, __RC__)
       call RegridO2A_1d(ALBVF, SURFST, 'ALBVF', &
            XFORM_O2A, locstreamO, __RC__)
       call RegridO2A_1d(ALBNR, SURFST, 'ALBNR', &
            XFORM_O2A, locstreamO, __RC__)
       call RegridO2A_1d(ALBNF, SURFST, 'ALBNF', &
            XFORM_O2A, locstreamO, __RC__)
          
    endif

    call MAPL_TimerOff(MAPL,    "-Albedo")

    deallocate(     TS_OLD,      __STAT__)
    deallocate(      FSURF,      __STAT__)
    deallocate(  DFSURFDTS,      __STAT__)
    deallocate(    DLHFDTS,      __STAT__)
    deallocate(      EVAPN,      __STAT__)
    deallocate(       LHFN,      __STAT__)
    deallocate(       RAIN,      __STAT__)

    !deallocate(ALBIN)
    !deallocate(ALBSN)
    !deallocate(ALBPND)

!  All done
!-----------

    RETURN_(ESMF_SUCCESS)
             
  end subroutine CICECORE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !=========CALLBACK==============
  subroutine RegridA2O_1d(ptr_1d, state, name, xform, locstream, rc)
    real :: ptr_1d(:)
    type (ESMF_State) :: state
    character(len=*) :: name
    type (MAPL_LocStreamXform), intent(IN) :: xform
    type (MAPL_LocStream), intent(IN) :: locstream
    integer, optional, intent(OUT) :: rc

    real, pointer, dimension(:,:)   :: PTR_2d  => null()
    real, pointer, dimension(:)     :: ptrO    => null()
    integer :: status
    integer :: k, nc, nt

    call MAPL_LocStreamGet(LocStream, nt_local=nt, __RC__)
!@@    if(NT == 0) then
!@@        RETURN_(ESMF_SUCCESS)
!@@    endif
    call MAPL_GetPointer(state, ptr_2d, name, __RC__)

    allocate(ptrO(NT), __STAT__)
    !
    !   !perform locstreamTrans_T2T+T2G(ts) => tsg (tile)
    !
    call MAPL_LocStreamTransform( ptrO, XFORM, PTR_1D, __RC__)
    call MAPL_LocStreamTransform( locStream, PTR_2D(:,:), ptrO, __RC__)
    deallocate (ptrO)

    RETURN_(ESMF_SUCCESS)
  end subroutine RegridA2O_1d

  subroutine RegridA2O_2d(ptr_2d, state, name, xform, locstream, rc)
    real :: ptr_2d(:,:)
    type (ESMF_State) :: state
    character(len=*) :: name
    type (MAPL_LocStreamXform), intent(IN) :: xform
    type (MAPL_LocStream), intent(IN) :: locstream
    integer, optional, intent(OUT) :: rc

    real, pointer, dimension(:,:,:) :: PTR_3D  => null()
    real, pointer, dimension(:)     :: ptrO    => null()
    integer :: status
    integer :: k, nc, nt

    call MAPL_LocStreamGet(LocStream, nt_local=nt, __RC__)
!@@    if(NT == 0) then
!@@        RETURN_(ESMF_SUCCESS)
!@@    endif
    NC = size(ptr_2d, 2)
    call MAPL_GetPointer(state, ptr_3d, name, __RC__)
    ASSERT_(NC == size(ptr_3d,3)) ! make sure the ungridded dims match

    allocate(ptrO(NT), __STAT__)
    do k = 1, nc
    !
    !   !perform locstreamTrans_T2T+T2G(ts) => tsg (tile)
    !
       call MAPL_LocStreamTransform( ptrO, XFORM, PTR_2D(:,k), __RC__ )
       call MAPL_LocStreamTransform( locStream, PTR_3D(:,:,k), ptrO, __RC__)
    end do
    deallocate (ptrO)

    RETURN_(ESMF_SUCCESS)
  end subroutine RegridA2O_2d

  subroutine RegridO2A_2d(ptr_2d, state, name, xform, locstream, rc)
    real :: ptr_2d(:,:)
    type (ESMF_State) :: state
    character(len=*) :: name
    type (MAPL_LocStreamXform), intent(IN) :: xform
    type (MAPL_LocStream), intent(IN) :: locstream
    integer, optional, intent(OUT) :: rc

    real, pointer, dimension(:,:,:) :: PTR_3D  => null()
    real, pointer, dimension(:)    :: ptrO  => null()
    integer :: status
    integer :: k, nc, nt

    call MAPL_LocStreamGet(LocStream, nt_local=nt, __RC__)
!@@    if(NT == 0) then
!@@       RETURN_(ESMF_SUCCESS)
!@@    endif
    NC = size(ptr_2d, 2)
    call MAPL_GetPointer(state, ptr_3d, name, __RC__)
    ASSERT_(NC == size(ptr_3d,3)) ! make sure the ungridded dims match

    allocate(ptrO(NT), __STAT__)
    do k = 1, nc

       !
       !   !perform locstreamTrans_G2T+T2T(ts) => tsg (tile)
       !
       call MAPL_LocStreamTransform( locStream, ptrO, ptr_3d(:,:,k), __RC__)
       call MAPL_LocStreamTransform( ptr_2d(:,k), XFORM,  ptrO, __RC__ )
    end do
    deallocate (ptrO)

    RETURN_(ESMF_SUCCESS)
  end subroutine RegridO2A_2d

  subroutine RegridO2A_1d(ptr_1d, state, name, xform, locstream, rc)
    real :: ptr_1d(:)
    type (ESMF_State) :: state
    character(len=*) :: name
    type (MAPL_LocStreamXform), intent(IN) :: xform
    type (MAPL_LocStream), intent(IN) :: locstream
    integer, optional, intent(OUT) :: rc

    real, pointer, dimension(:,:) :: PTR_2D  => null()
    real, pointer, dimension(:)   :: ptrO    => null()
    integer :: status
    integer :: nt

    call MAPL_LocStreamGet(LocStream, nt_local=nt, __RC__)
!@@    if(NT == 0) then
!@@       RETURN_(ESMF_SUCCESS)
!@@    endif
    call MAPL_GetPointer(state, ptr_2d, name, __RC__)

    allocate(ptrO(NT), __STAT__)

    !
    !   !perform locstreamTrans_G2T+T2T(ts) => tsg (tile)
    !
    call MAPL_LocStreamTransform( locStream, ptrO, ptr_2d, __RC__)
    call MAPL_LocStreamTransform( ptr_1d, XFORM,  ptrO, __RC__ )
    deallocate (ptrO)

    RETURN_(ESMF_SUCCESS)
  end subroutine RegridO2A_1d


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
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(STATUS)
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
    VERIFY_(STATUS)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL"   )
    call MAPL_TimerOn(MAPL,"FINALIZE")

    call MAPL_TimerOff(MAPL,"FINALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL"   )

! Generic Finalize
! ------------------
    
    call MAPL_GenericFinalize( GC, IMPORT, EXPORT, CLOCK, RC=status )
    VERIFY_(STATUS)

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Finalize

  subroutine Normalize(ptr, frac, set_undef)

     real,              dimension(:),  intent(inout)  :: ptr 
     real,              dimension(:),     intent(in)  :: frac 
     logical,             intent(in),      optional   :: set_undef 


     logical ::  l_set 

     if (present(set_undef)) then
        l_set = set_undef
     else
        l_set = .false.
     endif  
  
     if(l_set) then 
        where(frac > puny)  
           ptr = ptr / frac
        elsewhere
           ptr = MAPL_UNDEF
        endwhere  
     else
        where(frac > puny)  
           ptr = ptr / frac
        endwhere  
     endif 

     return  
  end subroutine Normalize 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine RetrieveTransforms(mapl, mystate, XFORM_A2O, XFORM_O2A, locStreamO, rc) 

    type (MAPL_MetaComp),         pointer, intent(in    ) :: MAPL
    type(seaice_interface_state), pointer, intent(inout ) :: mystate
    type(MAPL_LocStream),                  intent(out   ) :: locStreamO
    type(MAPL_LocStreamXform),             intent(out   ) :: xform_A2O
    type(MAPL_LocStreamXform),             intent(out   ) :: xform_O2A
    integer, optional,                     intent(out   ) :: rc

    type WRAP_X
            type (MAPL_LocStreamXform), pointer :: PTR => null()
    end type WRAP_X
    type WRAP_L
            type (MAPL_LocStream), pointer :: PTR => null()
    end type WRAP_L

    type(ESMF_GridComp) :: rootGC
    type(WRAP_X)        :: wrap_xform
    type(WRAP_L)        :: wrap_locstream
    integer             :: status

    character(len=ESMF_MAXSTR)              :: IAm

    Iam = "RetrieveTransforms"

    if (mystate%retrievedRootGC) then 
       XFORM_A2O = mystate%XFORM_A2O
       XFORM_O2A = mystate%XFORM_O2A
       locStreamO = mystate%locStreamO
       RETURN_(ESMF_SUCCESS)
    endif


    rootGC = MAPL_RootGcRetrieve(MAPL)
!??         ASSERT_(associated(rootGC%ptr))

    call ESMF_UserCompGetInternalState ( rootGC, 'GCM_XFORM_A2O', &
              wrap_xform,status )
    VERIFY_(STATUS)
    mystate%XFORM_A2O => wrap_xform%ptr
    call ESMF_UserCompGetInternalState ( rootGC, 'GCM_XFORM_O2A', &
              wrap_xform, status )
    VERIFY_(STATUS)
    mystate%XFORM_O2A => wrap_xform%ptr
    call ESMF_UserCompGetInternalState ( rootGC, 'GCM_LOCSTREAM_OCEAN', &
              wrap_locstream, status )
    VERIFY_(STATUS)
    mystate%locStreamO => wrap_locstream%ptr

    mystate%retrievedRootGC = .true.
    XFORM_A2O = mystate%XFORM_A2O
    XFORM_O2A = mystate%XFORM_O2A
    locStreamO = mystate%locStreamO

    RETURN_(ESMF_SUCCESS)

  end subroutine RetrieveTransforms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_SeaiceInterfaceGridComp

subroutine SetServices(gc, rc)
   use ESMF
   use GEOS_SeaiceInterfaceGridComp, only : mySetservices=>SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc,rc=rc)
end subroutine
