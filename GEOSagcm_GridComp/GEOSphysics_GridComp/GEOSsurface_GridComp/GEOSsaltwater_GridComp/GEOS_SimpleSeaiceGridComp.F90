!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
module GEOS_SimpleSeaiceGridCompMod

!BOP
! !MODULE: GEOS_SimpleSeaiceGridCompMod -- Implements slab saltwater tiles.

! !DESCRIPTION:
! 
!   {\tt GEOS\_SimpleSeaice} is a light-weight gridded component that updates
!      the skin sub-tiles at saltwater points, be they ocean, estuary, or salt
!      lake. Currently each tile can have only two subtiles, open-water and ice.
!      But the code is easily extensible to multiple ice types,
!      and includes implementation of LANL CICE thermodynamics on salt water tiles.
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
  use atmOcnIntlayer,     only: water_RHO

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

  integer, parameter    :: ICE = 1  
  integer, parameter    :: NUM_SUBTILES = 1      

  type ssi_state
       integer:: CHOOSEMOSFC
  end type ssi_state

  type ssi_state_wrap
      type(ssi_state), pointer :: ptr
  end type 

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

    type(ssi_state_wrap) :: wrap
    type(ssi_state), pointer :: mystate
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

! Sea-Ice Thermodynamics computation: using CICE or not?
!-------------------------------------------------------


    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run1, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2, RC=STATUS )
    VERIFY_(STATUS)

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

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRACINEW',                             &
        LONG_NAME          = 'ice_covered_fraction_of_tile_after_update',      &
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

     call MAPL_AddExportSpec(GC                    ,&
        SHORT_NAME         = 'PICE',                      &
        LONG_NAME          = 'sea_ice_pressure_loading'  ,&
        UNITS              = 'Pa'                        ,&
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

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'TSKINI',                            &
        LONG_NAME          = 'ice_skin_temperature',              &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = 'SEAICE',                            &
        DEFAULT            = MAPL_TICE-1.8,                               &
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


    call MAPL_AddImportSpec (GC,                                   &
         SHORT_NAME = 'DTSDT',                                     &
         LONG_NAME  = 'skin_temperature_analysis_tendency',        &
         UNITS      = 'K s-1',                                     &
         RESTART    = MAPL_RestartSkip,                            &
         DIMS       = MAPL_DimsTileOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'UUA',                             &
        LONG_NAME          = 'interpolated_effective_surface_zonal_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VVA',                             &
        LONG_NAME          = 'interpolated_effective_surface_meridional_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
        RC=STATUS  )
     VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TAUXBOT',                           &
        LONG_NAME          = 'eastward_stress_at_base_of_ice',    &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TAUYBOT',                           &
        LONG_NAME          = 'northward_stress_at_base_of_ice',   &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                                    &
        SHORT_NAME         = 'SSKINW',                            &
        LONG_NAME          = 'water_skin_salinity',               &
        UNITS              = 'psu',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 30.0,                                &

                                                       RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                                    &
        SHORT_NAME         = 'TSKINW',                            &
        LONG_NAME          = 'water_skin_temperature',            &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 280.0,                               &

                                                       RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'SS_FOUND',                          &
        LONG_NAME          = 'foundation_salinity_for_interface_layer',               &
        UNITS              = 'psu',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 30.0,                                &
        RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TS_FOUND',                          &
        LONG_NAME          = 'foundation_temperature_for_interface_layer',            &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 280.0,                               &
        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TFREEZE',                        &
        LONG_NAME          = 'freezing_temperature_for_interface_layer',&
        UNITS              = 'K',                               &
        DIMS               = MAPL_DimsTileOnly,                 &
        VLOCATION          = MAPL_VLocationNone,                &
        DEFAULT            = MAPL_TICE-1.8,                     &
        RC=STATUS  )
     VERIFY_(STATUS)

    
!-------------------Exports---------------------------------------------------------------


  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'FRESH',                     &
    LONG_NAME          = 'fresh_water_flux_to_ocean' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'FSALT',                     &
    LONG_NAME          = 'salt_flux_to_ocean'        ,&
    UNITS              = 'kg m-2 s-1'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'FHOCN',                     &
    LONG_NAME          = 'actual_ocean_ice_flux'     ,&
    UNITS              = 'W m-2'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHTSKIN',                   &
    LONG_NAME          = 'Ground_heating_for_skin_temp',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
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
    call ESMF_UserCompSetInternalState(gc, 'ssi_private', wrap,status)
    VERIFY_(status)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="RUN1"   ,               RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN2"  ,                RC=STATUS)
    VERIFY_(STATUS)
  
    call MAPL_TimerAdd(GC,    name="-Albedo"     ,          RC=STATUS)
    VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC,  RC=STATUS )
    VERIFY_(STATUS)
 
! Set the Run entry point
! -----------------------

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!BOP
! !IROUTINE: RUN1 -- First Run stage for the SimpleSeaice component
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
   real, pointer, dimension(:  )  :: FRACINEW=> null()

! pointers to internal

   real, pointer, dimension(:  )  :: TI  => null()
   real, pointer, dimension(:,:)  :: QS  => null()
   real, pointer, dimension(:,:)  :: CH  => null()
   real, pointer, dimension(:,:)  :: CM  => null()
   real, pointer, dimension(:,:)  :: CQ  => null()
   real, pointer, dimension(:,:)  :: WW  => null()
   real, pointer, dimension(:,:)  :: Z0  => null()

! pointers to import

   real, pointer, dimension(:)    :: UU  => null()
   real, pointer, dimension(:)    :: UWINDLMTILE => null()
   real, pointer, dimension(:)    :: VWINDLMTILE => null()
   real, pointer, dimension(:)    :: UI  => null()
   real, pointer, dimension(:)    :: VI  => null()
   real, pointer, dimension(:)    :: DZ  => null()     
   real, pointer, dimension(:)    :: TA  => null()
   real, pointer, dimension(:)    :: QA  => null()     
   real, pointer, dimension(:)    :: PS  => null()
   real, pointer, dimension(:)    :: PCU => null()
   real, pointer, dimension(:)    :: FI  => null()

   integer                        :: N
   integer                        :: NT
   integer                        :: NC
   integer                        :: niter


   real, allocatable              :: TS (:,:)
   real, allocatable              :: US (:,:)
   real, allocatable              :: VS (:,:)
   real, allocatable              :: FR (:,:)
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
   integer, allocatable           :: IWATER(:)
   real, allocatable              :: PSMB(:)
   real, allocatable              :: PSL(:)

   real            :: OCEANICEZ0
   real, parameter :: HPBL = 1000.

   integer         :: CHOOSEMOSFC
   integer         :: CHOOSEZ0
   type(ssi_state_wrap) :: wrap
   type(ssi_state), pointer :: mystate

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

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN1" )

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,                          &
         INTERNAL_ESMF_STATE = INTERNAL,         &
                                       RC=STATUS )
    VERIFY_(STATUS)

! Get parameters (0:Louis, 1:Monin-Obukhov)
! -----------------------------------------
    call ESMF_UserCompGetInternalState(gc,'ssi_private',wrap,status)
    VERIFY_(status)
    mystate => wrap%ptr
    CHOOSEMOSFC = mystate%CHOOSEMOSFC

    call MAPL_GetResource ( MAPL, CHOOSEZ0,    Label="CHOOSEZ0:",    DEFAULT=3, RC=STATUS)
    VERIFY_(STATUS)

! Get roughness parameters with and without CICE Thermodynamics
! -------------------------------------------------------------
    call MAPL_GetResource ( MAPL, OCEANICEZ0,  Label="OCEANICEZ0:" , DEFAULT=1.0e-3, RC=STATUS) 
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
   call MAPL_GetPointer(IMPORT,FI     , 'FRACICE',    RC=STATUS)
   VERIFY_(STATUS)

! Pointers to internals
!----------------------

   call MAPL_GetPointer(INTERNAL,TI   , 'TSKINI' , RC=STATUS)
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
   call MAPL_GetPointer(EXPORT,TF    , 'TFREEZE' ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FRACI , 'FRACI'   ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FRACINEW , 'FRACINEW'   ,    RC=STATUS)
   VERIFY_(STATUS)

   NT = size(TA)
   if(NT == 0) then
      call MAPL_TimerOff(MAPL,"RUN1" )
      call MAPL_TimerOff(MAPL,"TOTAL")
      RETURN_(ESMF_SUCCESS)
   end if

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
   allocate(US (NT,NUM_SUBTILES),   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(VS (NT,NUM_SUBTILES),   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(TS (NT,NUM_SUBTILES),   STAT=STATUS)
   VERIFY_(STATUS)
   allocate(FR (NT,NUM_SUBTILES),   STAT=STATUS)
   VERIFY_(STATUS)


   TS(:,ICE  ) = TI
   FR(:,ICE  ) = 1.0

   US(:,ICE  ) = UI
   VS(:,ICE  ) = VI

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

   N = ICE

      sfc_layer: if(CHOOSEMOSFC.eq.0) then
         call louissurface(1,N,UU,WW,PS,TA,TS,QA,QS,PCU,LAI,Z0,DZ,CM,CN,RIB,ZT,ZQ,CH,CQ,UUU,UCN,RE)

      elseif (CHOOSEMOSFC.eq.1) then

         niter = 6   ! number of internal iterations in the helfand MO surface layer routine
         IWATER= 5
         Z0(:,N)=OCEANICEZ0

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
         if(associated(MOU50M))MOU50M = U50M(:)*FR(:,N)
         if(associated(MOV50M))MOV50M = V50M(:)*FR(:,N)
         if(associated(MOT10M))MOT10M = T10M(:)*FR(:,N)
         if(associated(MOQ10M))MOQ10M = Q10M(:)*FR(:,N)
         if(associated(MOU10M))MOU10M = U10M(:)*FR(:,N)
         if(associated(MOV10M))MOV10M = V10M(:)*FR(:,N)
         if(associated(MOT2M ))MOT2M  = T2M (:)*FR(:,N)
         if(associated(MOQ2M ))MOQ2M  = Q2M (:)*FR(:,N)
         if(associated(MOU2M ))MOU2M  = U2M (:)*FR(:,N)
         if(associated(MOV2M ))MOV2M  = V2M (:)*FR(:,N)
      endif sfc_layer

      !  Aggregate to tiles
      !--------------------

                             CHB     = CH(:,N)*FR(:,N)
                             CQB     = CQ(:,N)*FR(:,N)
                             CMB     = CM(:,N)*FR(:,N)
         if(associated(TST)) TST     = TS(:,N)*FR(:,N)
         if(associated(QST)) QST     = QS(:,N)*FR(:,N)
         if(associated(CNT)) CNT     = CN(:  )*FR(:,N)
         if(associated(RIT)) RIT     = RIB(: )*FR(:,N)
         if(associated(RET)) RET     = RE(:  )*FR(:,N)
         if(associated(Z0O)) Z0O     = Z0(:,N)*FR(:,N)
         if(associated(Z0H)) Z0H     = ZT(:  )*FR(:,N)
         if(associated(VNT)) VNT     = UUU    *FR(:,N)

      !  Aggregate effective, CD-weighted, surface values of T and Q
      !-------------------------------------------------------------

         if(associated(TH)) TH      = CH(:,N)*TS(:,N)*FR(:,N)
         if(associated(QH)) QH      = CQ(:,N)*QS(:,N)*FR(:,N)
         if(associated(UH)) UH      = CM(:,N)*US(:,N)*FR(:,N)
         if(associated(VH)) VH      = CM(:,N)*VS(:,N)*FR(:,N)


      WW(:,N) = max(CH(:,N)*(TS(:,N)-TA-(MAPL_GRAV/MAPL_CP)*DZ)/TA + MAPL_VIREPS*CQ(:,N)*(QS(:,N)-QA),0.0)
      WW(:,N) = (HPBL*MAPL_GRAV*WW(:,N))**(2./3.)
      if(associated(GST)) GST     = WW(:,N)*FR(:,N)

   if(associated(CHT)) CHT = CHB
   if(associated(CQT)) CQT = CQB
   if(associated(CMT)) CMT = CMB

   if(associated(TF     )) TF      = MAPL_TICE - 5.0 ! set to much lower to prevent Openwater from
                                                     ! adjusting its TW in AMIP mode (zero-diff)
   if(associated(FRACI  )) FRACI   = FI
   if(associated(FRACINEW)) FRACINEW = FI

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
   deallocate(FR )
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
! !IROUTINE: RUN2 -- Second Run stage for the SimpleSeaice component

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

  real, pointer, dimension(:)         :: LATS => null()
  real, pointer, dimension(:)         :: LONS => null()

  real, pointer, dimension(:)         :: AREA => null()     ! needed to calculate TILEAREA in SaltWaterCore

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run2"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

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

    call SEAICECORE(NT=size(LONS), RC=STATUS )
    VERIFY_(STATUS)

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN2" )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine SEAICECORE(NT,RC)
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
   real, pointer, dimension(:  )  :: SHICE   => null()
   real, pointer, dimension(:  )  :: SHOUT   => null()
   real, pointer, dimension(:  )  :: HLATN   => null()
   real, pointer, dimension(:  )  :: HLATICE => null()
   real, pointer, dimension(:  )  :: HLWUP   => null()
   real, pointer, dimension(:  )  :: LWNDSRF => null()
   real, pointer, dimension(:  )  :: SWNDSRF => null()
   real, pointer, dimension(:  )  :: SWNDICE => null()
   real, pointer, dimension(:  )  :: LWNDICE => null()
   real, pointer, dimension(:  )  :: FSURF   => null()

   real, pointer, dimension(:  )  :: DELTS  => null()
   real, pointer, dimension(:  )  :: DELQS  => null()
   real, pointer, dimension(:  )  :: TST    => null()
   real, pointer, dimension(:  )  :: QST    => null()
   real, pointer, dimension(:  )  :: TAUXI  => null()
   real, pointer, dimension(:  )  :: TAUYI  => null()
   real, pointer, dimension(:  )  :: FRI    => null()

   real, pointer, dimension(:  )  :: DRUVRTHRU  => null()
   real, pointer, dimension(:  )  :: DFUVRTHRU  => null()
   real, pointer, dimension(:  )  :: DRPARTHRU  => null()
   real, pointer, dimension(:  )  :: DFPARTHRU  => null()
   real, pointer, dimension(:  )  :: FRESH      => null()
   real, pointer, dimension(:  )  :: FSALT      => null()
   real, pointer, dimension(:  )  :: FHOCN      => null()
   real, pointer, dimension(:  )  :: PICE       => null()
   real, pointer, dimension(:  )  :: GHTSKIN    => null()

! pointers to internal

   real, pointer, dimension(:  )  :: TI    => null()
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
   real, pointer, dimension(:)    :: UI => null()
   real, pointer, dimension(:)    :: VI => null()

   real, pointer, dimension(:)    :: TAUXBOT   => null()                  ! CICE related
   real, pointer, dimension(:)    :: TAUYBOT   => null()

   real, allocatable                   :: TS (:,:)   ! Following 4 Variables: TS to FR need to be 
   real, allocatable                   :: HH (:,:)   ! allocatable because NUM_SUBTILES is NOT a parameter
   real, allocatable                   :: SS (:,:)
   real, allocatable                   :: FR (:,:)
   real,    dimension(NT)              :: SHF
   real,    dimension(NT)              :: EVP
   real,    dimension(NT)              :: SHD
   real,    dimension(NT)              :: EVD
   real,    dimension(NT)              :: CFQ
   real,    dimension(NT)              :: CFT
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
   real,    dimension(NT)              :: ALBVRI
   real,    dimension(NT)              :: ALBVFI
   real,    dimension(NT)              :: ALBNRI
   real,    dimension(NT)              :: ALBNFI
   real,    dimension(NT)              :: VSUVR
   real,    dimension(NT)              :: VSUVF


   integer                             :: N
   real                                :: DT
   real                                :: MAXICEDEPTH
   real                                :: MINICEDEPTH

! following are related  to CICE

   integer                             :: NSUB, I, K, L


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

   real, parameter :: SALTWATERICECAP = MAPL_CAPICE

!  Begin...
!----------

   IAm =  trim(COMP_NAME) // "SEAICECORE"

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

   call MAPL_GetPointer(IMPORT,FI     , 'FRACICE',    RC=STATUS); VERIFY_(STATUS)

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
   call MAPL_GetPointer(EXPORT,TAUXI  , 'TAUXI', alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TAUYI  , 'TAUYI', alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,EVAPOUT, 'EVAPOUT' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SUBLIM,  'SUBLIM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SHOUT  , 'SHOUT'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SHICE  , 'SHICE'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLATN  , 'HLATN'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLATICE, 'HLATICE' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLWUP  , 'HLWUP'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,LWNDSRF, 'LWNDSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SWNDSRF, 'SWNDSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,LWNDICE, 'LWNDICE' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SWNDICE, 'SWNDICE' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FRI    , 'FRACI'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FSURF  , 'FSURF'   ,    RC=STATUS); VERIFY_(STATUS)

   call MAPL_GetPointer(EXPORT,DRUVRTHRU  , 'PENUVR'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,DFUVRTHRU  , 'PENUVF'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,DRPARTHRU  , 'PENPAR'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,DFPARTHRU  , 'PENPAF'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FRESH      , 'FRESH'      ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FSALT      , 'FSALT'      ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FHOCN      , 'FHOCN'      ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,PICE       , 'PICE'       ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,GHTSKIN    , 'GHTSKIN'    ,    RC=STATUS); VERIFY_(STATUS)

   allocate(TS (NT,NUM_SUBTILES),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(HH (NT,NUM_SUBTILES),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(SS (NT,NUM_SUBTILES),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(FR (NT,NUM_SUBTILES),STAT=STATUS)
   VERIFY_(STATUS)

! Get the time step
! -----------------

    call MAPL_Get(MAPL, HEARTBEAT = DT, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, DT, Label="DT:", DEFAULT=DT, RC=STATUS)
    VERIFY_(STATUS)

! Get parameters
! --------------

    call MAPL_GetResource ( MAPL, MAXICEDEPTH  , Label="MAX_SEAICE_DEPTH:", DEFAULT=2.0  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, MINICEDEPTH  , Label="MIN_SEAICE_DEPTH:", DEFAULT=1.E-6, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, EMSICE,        Label="CICE_EMSICE:",      DEFAULT=0.99999, RC=STATUS)
    VERIFY_(STATUS)

    MAXICEDEPTH     = MAXICEDEPTH  * water_RHO('fresh_water')
    MINICEDEPTH     = MINICEDEPTH  * water_RHO('fresh_water')


! Copy friendly internals into tile-tile local variables
!-------------------------------------------------------

    HH(:,ICE  ) = HI
    SS(:,ICE  ) = SI*HI
    TS(:,ICE  ) = TI
    FR(:,ICE  ) = 1.0

! Initialize PAR and UVR beam fluxes
!-----------------------------------

    VSUVR = DRPAR + DRUVR
    VSUVF = DFPAR + DFUVR

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

    ! Albedo over Sea-Ice.
    call ALBSEAICE (ALBVRI,ALBVFI,ALBNRI,ALBNFI,ZTH,LATS,CURRENT_TIME)  ! GEOS albedo over sea ice

    call MAPL_TimerOff(MAPL,    "-Albedo")

! Cycle through sub-tiles doing water and energy budget
!------------------------------------------------------

    if(associated(FRI    )) FRI     = FI
    if(associated(EVAPOUT)) EVAPOUT = 0.0
    if(associated(SHOUT  )) SHOUT   = 0.0
    if(associated(HLATN  )) HLATN   = 0.0
    if(associated(DELTS  )) DELTS   = 0.0
    if(associated(DELQS  )) DELQS   = 0.0

    ! exports to openwater 
    if(associated(FRESH    )) FRESH       = 0.0
    if(associated(FSALT    )) FSALT       = 0.0
    if(associated(FHOCN    )) FHOCN       = 0.0
    if(associated(PICE     )) PICE        = 0.0
    if(associated(DRUVRTHRU)) DRUVRTHRU   = 0.0
    if(associated(DFUVRTHRU)) DFUVRTHRU   = 0.0
    if(associated(DRPARTHRU)) DRPARTHRU   = 0.0
    if(associated(DFPARTHRU)) DFPARTHRU   = 0.0


! Without LANL CICE, i.e., with GEOS CICE, open water and sea-ice tiles were handled together. 
! But to keep the functionality of both (GEOS and LANL) sea-ice, accessible, they are now split.
! Generally speaking, following accounts for change in fluxes and state variables due to 
! fluxes at the top of interface (with LANL CICE, they are in "Thermo1", below). Changes due to 
! bottom of interface layer fluxes is handled later, in an implicit fashion 
! (open-water: SKIN_SST; sea-ice: with LANL CICE: "Thermo2"). 
! Note: 1. with GEOS CICE, there is no account of bottom flux (Max, please verify)
!       2. The sequence of computations: first over ice tiles and then over water is important to be 
!          able to reproduce Heracles-4_0 results (zero diff). 
! ---------------------------------------------------------------------------------------------------

! If using GEOS CICE Thermodynamics. If using LANL CICE, this is accomplished in Thermo1
!---------------------------------------------------------------------------------------

       N   = ICE
       CFT = (CH(:,N)/CTATM)
       CFQ = (CQ(:,N)/CQATM)
       EVP = CFQ*(EVAP + DEV*(QS(:,N)-QHATM))
       SHF = CFT*(SH   + DSH*(TS(:,N)-THATM))
       SHD = CFT*DSH
       EVD = CFQ*DEV*GEOS_DQSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.)
       DTS = LWDNSRF - (ALW + BLW*TS(:,N)) - SHF

       DTX = DT  / ( SALTWATERICECAP*HH(:,N))
       SWN = (1.-ALBVRI)*VSUVR + (1.-ALBVFI)*VSUVF + &
             (1.-ALBNRI)*DRNIR + (1.-ALBNFI)*DFNIR
       DTS = DTX * ( DTS + SWN - EVP*MAPL_ALHS )
       DTS = DTS   / ( 1.0 + DTX*(BLW + SHD + EVD*MAPL_ALHS) )
       DTS = DTS - max((TS(:,N) + DTS)-MAPL_TICE, 0.)
       EVP = EVP + EVD * DTS
       SHF = SHF + SHD * DTS
       LHF = EVP * MAPL_ALHS

       if(associated(HLATICE)) then
          WHERE( FI>0.0 )
             HLATICE = LHF
          ELSEWHERE
             HLATICE = MAPL_UNDEF
          ENDWHERE
       endif
       if(associated(  SHICE)) then
          WHERE( FI>0.0 )
             SHICE   = SHF
          ELSEWHERE
             SHICE   = MAPL_UNDEF
          ENDWHERE
       endif

! Update SEA-ICE surface temperature and moisture
!----------------------------------------

       TS(:,N) = TS(:,N) + DTS
       DQS     = GEOS_QSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.) - QS(:,N)
       QS(:,N) = QS(:,N) + DQS

       HH(:,N) = HH(:,N) + DT*(SNO - EVP)
       HH(:,N) = max(min(HH(:,N),  MAXICEDEPTH),  MINICEDEPTH)

       if(associated(SUBLIM  ))SUBLIM  =           EVP    *FR(:,N)
       if(associated(EVAPOUT)) EVAPOUT = EVAPOUT + EVP    *FR(:,N)
       if(associated(SHOUT  )) SHOUT   = SHOUT   + SHF    *FR(:,N)
       if(associated(HLATN  )) HLATN   = HLATN   + LHF    *FR(:,N)
       if(associated(DELTS  )) DELTS   = DELTS   + DTS*CFT*FR(:,N)
       if(associated(DELQS  )) DELQS   = DELQS   + DQS*CFQ*FR(:,N)

! Copy back to friendly internal variables
!-----------------------------------------

       HI = HH(:,ICE  )
       SI = SS(:,ICE  )/HI
       TI = TS(:,ICE  )

! Stress over ice
!----------------

       where( FI>0.0 )
          TXI = CM(:,ICE)*(UUA - UI)
          TYI = CM(:,ICE)*(VVA - VI)
       elsewhere
          TXI = MAPL_UNDEF
          TYI = MAPL_UNDEF
       end where


    if(associated(TAUXI)) TAUXI = TXI
    if(associated(TAUYI)) TAUYI = TYI
    if(associated(GHTSKIN)) GHTSKIN  = MAPL_UNDEF

! Copies for export
!------------------

    if(associated(TST    )) TST     = 0.0
    if(associated(QST    )) QST     = 0.0
    if(associated(HLWUP  )) HLWUP   = ALW 
    if(associated(LWNDSRF)) LWNDSRF = LWDNSRF - ALW


    if(associated(LWNDICE)) then 
          where( FI>0.0 ) 
             LWNDICE = LWDNSRF - ALW - BLW*TS(:,  ICE) 
          elsewhere
             LWNDICE = MAPL_UNDEF
          end where
    endif

    if(associated(TST    )) TST     = TST     +     TS(:,N)*FR(:,N)
    if(associated(QST    )) QST     = QST     +     QS(:,N)*FR(:,N)
    if(associated(LWNDSRF)) LWNDSRF = LWNDSRF - BLW*TS(:,N)*FR(:,N)
    if(associated(HLWUP  )) HLWUP   = HLWUP   + BLW*TS(:,N)*FR(:,N)
    if(associated(FSURF  )) FSURF   = SWN+LWDNSRF-(ALW+BLW*TS(:,ICE))-SHF-LHF

    EMISS = EMSICE*FR(:,ICE)
    ALBVR = ALBVRI*FR(:,ICE)
    ALBVF = ALBVFI*FR(:,ICE)
    ALBNR = ALBNRI*FR(:,ICE)
    ALBNF = ALBNFI*FR(:,ICE)

    if(associated(SWNDICE)) then 
          where( FI>0.0 ) 
             SWNDICE = (1.-ALBVRI)*VSUVR + (1.-ALBVFI)*VSUVF + &
                       (1.-ALBNRI)*DRNIR + (1.-ALBNFI)*DFNIR
          elsewhere
             SWNDICE = MAPL_UNDEF
          end where
    end if

    if(associated(SWNDSRF)) then 
       SWNDSRF = &
           (1.-ALBVR)*VSUVR + (1.-ALBVF)*VSUVF + &
           (1.-ALBNR)*DRNIR + (1.-ALBNF)*DFNIR
    endif

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
          
       call ALBSEAICE (ALBVRI,ALBVFI,ALBNRI,ALBNFI,ZTH,LATS,CURRENT_TIME)   ! GEOS CICE

       ALBVR = ALBVRI*FR(:,ICE)
       ALBVF = ALBVFI*FR(:,ICE)
       ALBNR = ALBNRI*FR(:,ICE)
       ALBNF = ALBNFI*FR(:,ICE)

    endif

    call MAPL_TimerOff(MAPL,    "-Albedo")
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    deallocate(TS)
    deallocate(HH)
    deallocate(SS)
    deallocate(FR)

!  All done
!-----------

    RETURN_(ESMF_SUCCESS)
             
  end subroutine SEAICECORE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: ALBSEAICE - Computes albedos as a function of  $cos(\zeta)$ over sea-ice surfaces

! !INTERFACE:

  subroutine ALBSEAICE (ALBVR,ALBVF,ALBNR,ALBNF,ZTH,LATS,currTime)

! !ARGUMENTS:

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
  end subroutine ALBSEAICE

end subroutine RUN2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_SimpleSeaiceGridCompMod

