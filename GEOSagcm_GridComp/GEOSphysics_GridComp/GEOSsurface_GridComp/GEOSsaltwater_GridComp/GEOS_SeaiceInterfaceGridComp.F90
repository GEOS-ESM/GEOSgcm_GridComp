
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

!  !INTERNAL STATE:

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'HSKINI',                            &
        LONG_NAME          = 'ice_skin_layer_mass',               &
        UNITS              = 'kg m-2',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = 0.5*MAPL_RHOWTR,                     &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

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
                                                       RC=STATUS  )
     VERIFY_(STATUS)

!  !IMPORT STATE:

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

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'LWDNSRF',                           &
        LONG_NAME          = 'surface_downwelling_longwave_flux', &
        UNITS              = 'W m-2',                             &
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
        LONG_NAME          = 'eastward_surface_stress',           &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUX',                              &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'northward_surface_stress',          &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUY',                              &
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
        SHORT_NAME         = 'UHATM',                             &
        LONG_NAME          = 'effective_surface_zonal_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VHATM',                             &
        LONG_NAME          = 'effective_surface_meridional_velocity',&
        UNITS              = 'm s-1',                             &
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

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CMATM',                             &
        LONG_NAME          = 'surface_exchange_coefficient_for_momentum', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'FRACICE',                           &
        LONG_NAME          = 'ice_covered_fraction_of_tile',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
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
   real(kind=MAPL_R8), pointer, dimension(:,:)   :: FR      => null()  

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
   real(kind=MAPL_R8), allocatable:: FRI(:)


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
    call MAPL_GetResource ( MAPL, NUM_ICE_LAYERS    , Label="CICE_N_ICE_LAYERS:"     ,     RC=STATUS)
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
    call ESMF_UserCompGetInternalState(gc,'cice_private',wrap,status)
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
    DTDB = REAL(DT, kind=MAPL_R8)

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
   call MAPL_GetPointer(IMPORT,SW     , 'SS_FOUND' ,  RC=STATUS)
   VERIFY_(STATUS)
   ! the call below may be needed in dual-ocean mode 
   !   call MAPL_GetPointer(IMPORT,FI     , 'FRACICE',    RC=STATUS)
   !   VERIFY_(STATUS)


! Pointers to internals
!----------------------

   call MAPL_GetPointer(INTERNAL,FR   , 'FR'     ,    RC=STATUS)
   VERIFY_(STATUS) 
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
   call MAPL_GetPointer(EXPORT,QSAT1 , 'QSAT1'   ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,QSAT2 , 'QSAT2'   ,    RC=STATUS)
   VERIFY_(STATUS)

  ! export to openwater
   call MAPL_GetPointer(EXPORT,FRACI , 'FRACI'   ,    RC=STATUS)
   VERIFY_(STATUS)

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


! pointers to internal

   real, pointer, dimension(:,:)  :: TI    => null()
   real, pointer, dimension(:  )  :: HI    => null()
   real, pointer, dimension(:  )  :: SI    => null()
   real, pointer, dimension(:,:)  :: QS    => null()
   real, pointer, dimension(:,:)  :: CH    => null()
   real, pointer, dimension(:,:)  :: CQ    => null()
   real, pointer, dimension(:,:)  :: CM    => null()


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

   real, pointer, dimension(:)    :: TW        => null()
   real, pointer, dimension(:)    :: SW        => null()
   real, pointer, dimension(:,:)  :: TS        => null()

   !*CALLBACK*
   real, pointer, dimension(:,:)  :: AS_PTR_2D => null()
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
   real,    dimension(NT)              :: DTX
   real,    dimension(NT)              :: DTY
   real,    dimension(NT)              :: SWN
   real,    dimension(NT)              :: PEN
   real,    dimension(NT)              :: LHF
   real,    dimension(NT)              :: ZTH
   real,    dimension(NT)              :: SLR
   real,    dimension(NT)              :: VSUVR
   real,    dimension(NT)              :: VSUVF

   integer                             :: N
   real                                :: DT
   real                                :: MAXSALINITY
   real                                :: MINSALINITY


   integer                             :: NSUB, I, K, L

   real                                :: LATSO, LONSO
   real,    dimension(1)               :: LATSD, LONSD
   logical, dimension(1)               :: OBSERVE


   real,    dimension(NT)              :: FSWABS
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

   real,               allocatable    :: DRUVRTHRU     (:,:)
   real,               allocatable    :: DFUVRTHRU     (:,:)
   real,               allocatable    :: DRPARTHRU     (:,:)
   real,               allocatable    :: DFPARTHRU     (:,:)

   real,               allocatable    :: FSURF         (:) ! 



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


   real                                :: EMSICE

   real, parameter                     :: SALTWATERCAP    = MAPL_CAPWTR
   real, parameter                     :: SALTWATERICECAP = MAPL_CAPICE

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

   call MAPL_GetPointer(IMPORT,UW     , 'UW'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VW     , 'VW'     ,    RC=STATUS); VERIFY_(STATUS)
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

   call MAPL_GetPointer(IMPORT,TAUXBOT, 'TAUXBOT',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TAUYBOT, 'TAUYBOT',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TW     , 'TS_FOUND',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SW     , 'SS_FOUND',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,FRZMLT , 'FRZMLT' ,    RC=STATUS); VERIFY_(STATUS)

! Pointers to internals
!----------------------

   call MAPL_GetPointer(INTERNAL,TI     ,'TSKINI',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,HI     ,'HSKINI',    RC=STATUS); VERIFY_(STATUS)
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
   call MAPL_GetPointer(EXPORT,PENUVR , 'PENUVR'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,PENUVF , 'PENUVF'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,PENPAR , 'PENPAR'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,PENPAF , 'PENPAF'  , RC=STATUS); VERIFY_(STATUS)
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
   call MAPL_GetPointer(EXPORT,LWDNSRFe,'LWDNSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SWDNSRFe,'SWDNSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TSKINICE,'TSKINICE',    RC=STATUS); VERIFY_(STATUS)


   ! category dimensional exports
   call MAPL_GetPointer(EXPORT,SHICEN   ,  'SHICEN'    ,  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLWUPN   ,  'HLWUPN'    ,  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,LWNDSRFN ,  'LWNDSRFN'  ,  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FSURFN   ,  'FSURFN'    ,  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TSURFN   ,  'TSURFN'    ,  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FSWSFCN  ,  'FSWSFCN'   ,  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBINe   ,  'ALBIN'     ,  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBSNe   ,  'ALBSN'     ,  RC=STATUS); VERIFY_(STATUS)

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


    FR_OLD     = FRCICE   ! FRCICE is initialized by above subroutine CICE_PREP_THERMO
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

    call CICE_ALBSEAICE (ICE,NUM_ICE_CATEGORIES,NUM_ICE_LAYERS,NUM_SNOW_LAYERS,NT,DO_POND,LATSO,LONSO,LATS,LONS,ZTH,FR8,TS,&
                           DRPAR,DFPAR,DRNIR,DFNIR,DRUVR,DFUVR,VSUVR,VSUVF,VOLICE,VOLSNO,APONDN,HPONDN,                      &
                           ISWABS,FSWSFC, FSWINT,FSWTHRU,SSWABS,ALBIN,ALBSN,ALBPND,ALBVRN,ALBVFN,ALBNRN,ALBNFN,              &
                           DRUVRTHRU,DFUVRTHRU,DRPARTHRU,DFPARTHRU,RC=STATUS)
      VERIFY_(STATUS)                     
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


    call MAPL_TimerOn(MAPL,   "-Thermo1")

! 1st Step of LANL CICE Thermodynamics
! ------------------------------------

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
                           FSWSFC,FSWINT,FSWABS,LWDNSRF,EVD,SHD,SNO,SBLX,RC=STATUS)
          VERIFY_(STATUS)                 

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


    !*CALLBACK*
    !==============================================================================================

! retrieve the regridding from the "private" wrapped state
! if block only executes during the first time step
! check with Atanas about moving this block to Initialize()?

    if (.not. mystate%retrievedRootGC) then
       block
         type WRAP_X
            type (MAPL_LocStreamXform), pointer :: PTR => null()
         end type WRAP_X
         type WRAP_L
            type (MAPL_LocStream), pointer :: PTR => null()
         end type WRAP_L

         type(ESMF_GridComp) :: rootGC
         type(WRAP_X) :: wrap_xform
         type(WRAP_L) :: wrap_locstream

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
       end block
       mystate%retrievedRootGC = .true.
    end if
    XFORM_A2O = mystate%XFORM_A2O
    XFORM_O2A = mystate%XFORM_O2A
    locStreamO = mystate%locStreamO

    call ESMF_StateGet(IMPORT, 'SURFSTATE', SURFST, __RC__)

    call RegridA2O_2d(TS, SURFST, 'surface_ice_temperature',     &
                      XFORM_A2O, locstreamO, __RC__)

    ! ** execute the sea ice thermo coupling method
    call ESMF_MethodExecute(SURFST, label="thermo_coupling", userRC=AS_STATUS, RC=STATUS)
    VERIFY_(AS_STATUS)
    VERIFY_(STATUS)

    allocate(AS_PTR_2D(size(TS,1),size(TS,2)), __STAT__)
    call RegridO2A_2d(AS_PTR_2D, SURFST, 'surface_ice_temperature', &
         XFORM_O2A, locstreamO, __RC__)

    !************************************************************************************************
    !
    !      surface flux and temperature updated by sea ice
    !      continue updating relevant fields and pass them to surf
    !    DTS = AS_PTR_2D - TS
    !    .....
    !
    deallocate(AS_PTR_2D)


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

     if(associated(TSKINICE)) then
          TSKINICE = sum(TS(:,ICE:)*FR8(:,ICE:),dim=2)
          where(FRCICE > puny)
             TSKINICE = TSKINICE / FRCICE
          elsewhere
             TSKINICE = MAPL_UNDEF
          end where
    endif

    if(associated(FSURFICE)) then
          where(FRCICE > puny)
             FSURFICE = FSURFICE / FRCICE
          elsewhere
             FSURFICE = MAPL_UNDEF
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
          
       call CICE_ALBSEAICE (ICE,NUM_ICE_CATEGORIES,NUM_ICE_LAYERS,NUM_SNOW_LAYERS,NT,DO_POND,LATSO,LONSO,LATS,LONS,ZTH,FR8,TS,&
                            DRPAR,DFPAR,DRNIR,DFNIR,DRUVR,DFUVR,VSUVR,VSUVF,VOLICE,VOLSNO,APONDN,HPONDN,                      &
                            ISWABS,FSWSFC, FSWINT,FSWTHRU,SSWABS,ALBIN,ALBSN,ALBPND,ALBVRN,ALBVFN,ALBNRN,ALBNFN,              &
                            DRUVRTHRU,DFUVRTHRU,DRPARTHRU,DFPARTHRU,RC=STATUS)
       VERIFY_(STATUS)                     

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

!  All done
!-----------

    RETURN_(ESMF_SUCCESS)
             
  end subroutine CICECORE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !=========CALLBACK==============
  subroutine RegridA2O_2d(ptr_2d, state, name, xform, locstream, rc)
    real :: ptr_2d(:,:)
    type (ESMF_State) :: state
    character(len=*) :: name
    type (MAPL_LocStreamXform), intent(IN) :: xform
    type (MAPL_LocStream), intent(IN) :: locstream
    integer, optional, intent(OUT) :: rc

    real, pointer, dimension(:,:,:) :: PTR_3D  => null()
    real, pointer, dimension(:) :: ptrO  => null()
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

    call MAPL_GetResource ( MAPL, SHORTWAVE,  Label="CICE_SHORTWAVE:" ,  DEFAULT="shortwave_ccsm" , RC=STATUS)
    VERIFY_(STATUS)

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

  
    end do TILES ! K loop

    RETURN_(ESMF_SUCCESS)
  end subroutine CICE_THERMO1


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

    call MAPL_GetResource ( MAPL, SHORTWAVE, Label="CICE_SHORTWAVE:" , DEFAULT="shortwave_ccsm" , RC=STATUS)
    VERIFY_(STATUS)

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

      call ESMF_TimeGet  ( currTime, TimeString=string  ,rc=STATUS ) ; VERIFY_(STATUS)
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

end module GEOS_SeaiceInterfaceGridComp


