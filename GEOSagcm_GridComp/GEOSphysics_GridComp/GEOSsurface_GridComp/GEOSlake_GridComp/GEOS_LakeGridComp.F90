!  $Id$

#include "MAPL_Generic.h"


!=============================================================================
module GEOS_LakeGridCompMod

!BOP
! !MODULE: GEOS_LakeGridCompMod -- Implements slab lake tiles.
!
! !USES:

  use sfclayer  ! using module that contains sfc layer code
  use ESMF
  use MAPL_Mod
  use GEOS_UtilsMod
  use DragCoefficientsMod
  
  implicit none
  private

  integer, parameter :: ICE   = 1
  integer, parameter :: WATER = 2
  integer, parameter :: NUM_SUBTILES = 2

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !DESCRIPTION:
! 
!   {\tt GEOS\_Lake} is a light-weight gridded component that updates
!      the lake tiles
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

!DESCRIPTION: 
! This version uses the MAPL\_GenericSetServices, which sets
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


!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run1, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

!
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
          LONG_NAME          = 'runoff_flux',&
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'RUNOFF'                    ,&
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
        LONG_NAME          = 'surface_pressure',                  &
        UNITS              = 'Pa',                                &
        SHORT_NAME         = 'PS',                                &
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
        LONG_NAME          = 'ocean_ustar_cubed',         &
        UNITS              = 'm+3 s-3'                   ,&
        SHORT_NAME         = 'OUSTAR3'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'PENUVF',                    &
    LONG_NAME          = 'downwelling_uvr_diffuse_flux_at_skin_base',&
    UNITS              = 'W m-2'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'PENUVR',                    &
    LONG_NAME          = 'downwelling_uvr_direct_flux_at_skin_base',&
    UNITS              = 'W m-2'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'PENPAF',                    &
    LONG_NAME          = 'downwelling_par_diffuse_flux_at_skin_base',&
    UNITS              = 'W m-2'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    SHORT_NAME         = 'PENPAR',                    &
    LONG_NAME          = 'downwelling_par_direct_flux_at_skin_base',&
    UNITS              = 'W m-2'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

!  !Internal state:

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'TS',                                &
        LONG_NAME          = 'surface_skin_temperature',          &
        UNITS              = 'K',                                 &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 280.0,                               &
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
        SHORT_NAME         = 'FR',                                &
        LONG_NAME          = 'ice_fraction',                      &
        UNITS              = '1',                                 &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 1.0/NUM_SUBTILES,                    &
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

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'liquid_water_convective_precipitation', &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'PCU',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'liquid_water_large_scale_precipitation', &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'PLS',                               &
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

! Set the Run entry point
! -----------------------

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!BOP
! !IROUTINE: RUN1 -- First Run stage for the Lake component
!INTERFACE:
subroutine RUN1 ( GC, IMPORT, EXPORT, CLOCK, RC )

  !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

!DESCRIPTION: 
! Periodically refreshes the ozone mixing ratios.

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

! pointers to export

   real, pointer, dimension(:  )  :: TH
   real, pointer, dimension(:  )  :: QH
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
   real, pointer, dimension(:)    :: PCU

! locals

   integer                        :: N
   integer                        :: NT
   integer                        :: niter


   real, allocatable              :: TVA(:)
   real, allocatable              :: TVS(:)
   real, allocatable              :: URA(:)
   real, allocatable              :: UUU(:)
   real, allocatable              :: LAI(:)
   real, allocatable              :: CHH(:)
   real, allocatable              :: CQQ(:)
   real, allocatable              :: CN (:)
   real, allocatable              :: RE (:)
   real, allocatable              :: ZT (:)
   real, allocatable              :: ZQ (:)
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

   real, parameter :: LAKEZ0     = 1.0E-5
   real, parameter :: LAKEICEZ0  = 0.01

   integer                        :: CHOOSEMOSFC
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

! Get parameters (0:Louis, 1:Monin-Obukhov)
! -----------------------------------------
    call MAPL_GetResource ( MAPL, CHOOSEMOSFC, Label="CHOOSEMOSFC:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, CHOOSEZ0, Label="CHOOSEZ0:", DEFAULT=3, RC=STATUS)
    VERIFY_(STATUS)

! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,UU     , 'UU'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DZ     , 'DZ'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TA     , 'TA'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,QA     , 'QA'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PS     , 'PS'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PCU    , 'PCU' , RC=STATUS); VERIFY_(STATUS)

   call MAPL_GetPointer(IMPORT,UWINDLMTILE, 'UWINDLMTILE',    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VWINDLMTILE, 'VWINDLMTILE',    RC=STATUS)
   VERIFY_(STATUS)

! Pointers to internals
!----------------------

   call MAPL_GetPointer(INTERNAL,TS   , 'TS'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,QS   , 'QS'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,FR   , 'FR'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CH   , 'CH'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CM   , 'CM'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CQ   , 'CQ'  , RC=STATUS); VERIFY_(STATUS)

! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,QH    , 'QH'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TH    , 'TH'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,QST   , 'QST'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TST   , 'TST'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CHT   , 'CHT'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CMT   , 'CMT'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CQT   , 'CQT'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CNT   , 'CNT'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RIT   , 'RIT'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,VNT   , 'VENT' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,GST   , 'GUST' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Z0EXP , 'Z0'   , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Z0HEXP, 'Z0H'  , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOT2M, 'MOT2M' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOQ2M, 'MOQ2M' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOU2M, 'MOU2M' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOV2M, 'MOV2M' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOT10M, 'MOT10M',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOQ10M, 'MOQ10M',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOU10M, 'MOU10M',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOV10M, 'MOV10M',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOU50M, 'MOU50M',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOV50M, 'MOV50M',RC=STATUS); VERIFY_(STATUS)

   NT = size(TA)

   allocate(TVA(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(TVS(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(URA(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(UUU(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(LAI(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CHH(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CQQ(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(RE (NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CN (NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(ZT (NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(ZQ (NT),STAT=STATUS)
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

!  Compute drag corfficient at tiles
!-----------------------------------

   CHH = 0.0
   CQQ = 0.0

   if(associated(CMT   )) CMT    = 0.0
   if(associated(CNT   )) CNT    = 0.0
   if(associated(RIT   )) RIT    = 0.0
   if(associated(Z0EXP )) Z0EXP  = 0.0
   if(associated(Z0HEXP)) Z0HEXP = 0.0
   if(associated(TH    )) TH     = 0.0
   if(associated(QH    )) QH     = 0.0
   if(associated(TST   )) TST    = 0.0
   if(associated(QST   )) QST    = 0.0
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

   do N=1,NUM_SUBTILES

! Choose sfc layer: if CHOOSEMOSFC is 1, choose helfand MO,
!                   if CHOOSEMOSFC is 0 (default), choose louis

   if(CHOOSEMOSFC.eq.0) then

         LAI = 0.
    call louissurface(2,N,UU,WW,PS,TA,TS,QA,QS,PCU,LAI,Z0,DZ,CM,CN,RIB,ZT,ZQ,CH,CQ,UUU,UCN,RE)

   elseif (CHOOSEMOSFC.eq.1)then

      niter = 6   ! number of internal iterations in the helfand MO surface layer routine
      IWATER=2
      if(N==ICE) then
        Z0(:,N) = LAKEICEZ0
      else
        Z0(:,N) = LAKEZ0
      end if

    PSMB = PS * 0.01            ! convert to MB
         LAI  = 1.e-4

         ! Approximate pressure at top of surface layer: 
         !  hydrostatic, eqn of state using avg temp and press

    PSL = PSMB * (1. - (DZ*MAPL_GRAV)/(MAPL_RGAS*(TA+TS(:,N)) ) ) /   &
               (1. + (DZ*MAPL_GRAV)/(MAPL_RGAS*(TA+TS(:,N)) ) )

         call helfsurface( UWINDLMTILE,VWINDLMTILE,TA,TS(:,N),QA,QS(:,N),PSL,PSMB,Z0(:,N),LAI,  &
                      IWATER,DZ,niter,nt,RHO,VKH,VKM,USTAR,XX,YY,CU,CT,RIB,ZETA,WS,   &
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

                             CHH    = CHH + CH(:,N)*FR(:,N)
                             CQQ    = CQQ + CQ(:,N)*FR(:,N)
      if(associated(CMT   )) CMT    = CMT + CM(:,N)*FR(:,N)
      if(associated(CNT   )) CNT    = CNT + CN(:  )*FR(:,N)
      if(associated(RIT   )) RIT    = RIT + RIB(:  )*FR(:,N)
      if(associated(Z0EXP )) Z0EXP  = Z0EXP  + Z0(:,N)*FR(:,N)   
      if(associated(Z0HEXP)) Z0HEXP = Z0HEXP + ZT*FR(:,N)

      if(associated(TH    )) TH      = TH  + CH(:,N)*TS(:,N)*FR(:,N)
      if(associated(QH    )) QH      = QH  + CQ(:,N)*QS(:,N)*FR(:,N)
      if(associated(TST   )) TST     = TST + TS(:,N)*FR(:,N)
      if(associated(QST   )) QST     = QST + QS(:,N)*FR(:,N)

   end do

   if(associated(TH )) TH  = TH /CHH
   if(associated(QH )) QH  = QH /CQQ
   if(associated(VNT)) VNT = UUU
   if(associated(GST)) GST = 0.0
   if(associated(CHT)) CHT = CHH
   if(associated(CQT)) CQT = CQQ

   deallocate(TVA)
   deallocate(TVS)
   deallocate(URA)
   deallocate(UUU)
   deallocate(LAI)
   deallocate(CHH)
   deallocate(CQQ)
   deallocate(RE )
   deallocate(CN )
   deallocate(Z0 )
   deallocate(ZT )
   deallocate(ZQ )
   deallocate(WW )
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
! !IROUTINE: RUN2 -- Second Run stage for the Lake component

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
  type (MAPL_SunOrbit)                :: ORBIT
  real                                :: DT_SOLAR
  type (ESMF_TimeInterval)            :: TINT
  type (ESMF_TimeInterval)            :: DELT
  type (ESMF_Time)                    :: CURRENT_TIME
  type (ESMF_Time)                    :: NOW
  type (ESMF_Time)                    :: BEFORE
  type (ESMF_Time)                    :: MODELSTART
  type (ESMF_Alarm)                   :: SOLALARM
  logical                             :: solalarmison
  type(ESMF_VM)                       :: VM
  logical                             :: debugzth

  real, pointer                       :: LONS(:)
  real, pointer                       :: LATS(:)


!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Run"

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
         CF        = CF,                         &
         INTERNAL_ESMF_STATE = INTERNAL,         &
         TILELONS  = LONS,                       &
         TILELATS  = LATS,                       &
         ORBIT     = ORBIT,                      &
         RUNALARM  = ALARM,                      &
                                       RC=STATUS )
    VERIFY_(STATUS)

! Do the calculations
!--------------------

    if ( ESMF_AlarmIsRinging(ALARM, RC=STATUS) ) then
       VERIFY_(STATUS)
       call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
       VERIFY_(STATUS)
       call LAKECORE(NT=size(LATS), RC=STATUS )
       VERIFY_(STATUS)
    end if
    VERIFY_(STATUS)

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN2" )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine LAKECORE(NT,RC)
   integer,           intent( IN) :: NT
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
   real, pointer, dimension(:  )  :: EVAPOUT
   real, pointer, dimension(:  )  :: RUNOFF
   real, pointer, dimension(:  )  :: SUBLIM
   real, pointer, dimension(:  )  :: SHOUT
   real, pointer, dimension(:  )  :: HLATN
   real, pointer, dimension(:  )  :: HLWUP
   real, pointer, dimension(:  )  :: SWNDSRF
   real, pointer, dimension(:  )  :: LWNDSRF

! pointers to internal

   real, pointer, dimension(:,:)  :: TS
   real, pointer, dimension(:,:)  :: QS
   real, pointer, dimension(:,:)  :: FR
   real, pointer, dimension(:,:)  :: CH
   real, pointer, dimension(:,:)  :: CM
   real, pointer, dimension(:,:)  :: CQ

! pointers to import

   real, pointer, dimension(:)    :: ALW
   real, pointer, dimension(:)    :: BLW
   real, pointer, dimension(:)    :: LWDNSRF
   real, pointer, dimension(:)    :: EVAP 
   real, pointer, dimension(:)    :: SH 
   real, pointer, dimension(:)    :: DEV   
   real, pointer, dimension(:)    :: DSH 
   real, pointer, dimension(:)    :: SNO    
   real, pointer, dimension(:)    :: PCU
   real, pointer, dimension(:)    :: PLS    
   real, pointer, dimension(:)    :: PS     
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

   real,          dimension(NT)   :: SHF
   real,          dimension(NT)   :: LHF
   real,          dimension(NT)   :: EVP
   real,          dimension(NT)   :: RNF
   real,          dimension(NT)   :: EVD
   real,          dimension(NT)   :: SHD
   real,          dimension(NT)   :: SWN

   real,          dimension(NT)   :: CFQ
   real,          dimension(NT)   :: CFT
   real,          dimension(NT)   :: DQS
   real,          dimension(NT)   :: DTS
   real,          dimension(NT)   :: DTX
   real,          dimension(NT)   :: ZTH
   real,          dimension(NT)   :: SLR
   real,          dimension(NT)   :: ALBVRO
   real,          dimension(NT)   :: ALBVFO
   real,          dimension(NT)   :: ALBNRO
   real,          dimension(NT)   :: ALBNFO
   real,          dimension(NT)   :: ALBVRI
   real,          dimension(NT)   :: ALBVFI
   real,          dimension(NT)   :: ALBNRI
   real,          dimension(NT)   :: ALBNFI
   real,          dimension(NT)   :: VSUVR
   real,          dimension(NT)   :: VSUVF

   real                           :: DT
   integer                        :: N, I

   real, parameter :: LAKEEMISS    = 0.99070 
   real, parameter :: LAKEICEEMISS = 0.99999 
   real, parameter :: LAKEDEPTH    = 2.00
   real, parameter :: LAKEICEDEPTH = 0.02
   real, parameter :: LAKECAP      = (MAPL_RHOWTR*MAPL_CAPWTR*LAKEDEPTH    )
   real, parameter :: LAKEICECAP   = (MAPL_RHOWTR*MAPL_CAPWTR*LAKEICEDEPTH )


!  Begin...
!----------

   IAm =  trim(COMP_NAME) // "LAKECORE"

! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,ALW    , 'ALW'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,BLW    , 'BLW'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,LWDNSRF, 'LWDNSRF',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,EVAP   , 'EVAP'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SH     , 'SH'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DEV    , 'DEVAP'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DSH    , 'DSH'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SNO    , 'SNO'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PCU    , 'PCU'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PLS    , 'PLS'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PS     , 'PS'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,THATM  , 'THATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,QHATM  , 'QHATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,CTATM  , 'CTATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,CQATM  , 'CQATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DRPAR  , 'DRPAR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DFPAR  , 'DFPAR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DRNIR  , 'DRNIR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DFNIR  , 'DFNIR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DRUVR  , 'DRUVR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DFUVR  , 'DFUVR'  ,    RC=STATUS); VERIFY_(STATUS)

! Pointers to internals
!----------------------

   call MAPL_GetPointer(INTERNAL,TS     , 'TS'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,QS     , 'QS'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,FR     , 'FR'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CH     , 'CH'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CM     , 'CM'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CQ     , 'CQ'   ,    RC=STATUS); VERIFY_(STATUS)

! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,QST   , 'QST'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TST   , 'TST'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,DELTS , 'DELTS'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,DELQS , 'DELQS'   ,    RC=STATUS); VERIFY_(STATUS)

   call MAPL_GetPointer(EXPORT,EVAPOUT, 'EVAPOUT' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RUNOFF , 'RUNOFF'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SHOUT  , 'SHOUT'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLATN  , 'HLATN'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLWUP  , 'HLWUP'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,LWNDSRF, 'LWNDSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SWNDSRF, 'SWNDSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SUBLIM , 'SUBLIM'  ,    RC=STATUS); VERIFY_(STATUS)

   call MAPL_GetPointer(EXPORT,EMISS , 'EMIS' , alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBVF , 'ALBVF', alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBVR , 'ALBVR', alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBNF , 'ALBNF', alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ALBNR , 'ALBNR', alloc=.true., RC=STATUS); VERIFY_(STATUS)


! Get the time step
! -----------------

    call MAPL_Get(MAPL, HEARTBEAT = DT, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, DT, Label="DT:", DEFAULT=DT, RC=STATUS)
    VERIFY_(STATUS)
    
    call ESMF_VMGetCurrent ( VM, RC=STATUS )

        ! --------------------------------------------------------------------------
        ! Get the current time. 
        ! --------------------------------------------------------------------------

    debugzth = .false.

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

    call ALBLAKE   (ALBVRO,ALBVFO,ALBNRO,ALBNFO,ZTH)
    call ALBLAKEICE(ALBVRI,ALBVFI,ALBNRI,ALBNFI,ZTH)

    VSUVR = DRPAR + DRUVR
    VSUVF = DFPAR + DFUVR

    if(associated(EVAPOUT)) EVAPOUT = 0.0
    if(associated(RUNOFF )) RUNOFF  = 0.0
    if(associated(SHOUT  )) SHOUT   = 0.0
    if(associated(HLATN  )) HLATN   = 0.0
    if(associated(DELTS  )) DELTS   = 0.0
    if(associated(DELQS  )) DELQS   = 0.0
    if(associated(SWNDSRF)) SWNDSRF = 0.0
    if(associated(TST    )) TST     = 0.0
    if(associated(QST    )) QST     = 0.0
    if(associated(HLWUP  )) HLWUP   = ALW 
    if(associated(LWNDSRF)) LWNDSRF = LWDNSRF - ALW

    do N=1,NUM_SUBTILES
       CFT = (CH(:,N)/CTATM)
       CFQ = (CQ(:,N)/CQATM)
       EVP = CFQ*(EVAP + DEV*(QS(:,N)-QHATM))
       SHF = CFT*(SH   + DSH*(TS(:,N)-THATM))
       SHD = CFT*DSH
       EVD = CFQ*DEV*GEOS_DQSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.)
       DTS = LWDNSRF - (ALW + BLW*TS(:,N)) - SHF

       if (N==WATER) then
          DTX = (DT/LAKECAP)*FR(:,N) ! FR accounts for skin under ice
          SWN =   (1.-ALBVRO)*VSUVR + (1.-ALBVFO)*VSUVF + &
                  (1.-ALBNRO)*DRNIR + (1.-ALBNFO)*DFNIR
          DTS = DTX * ( DTS + SWN - EVP*MAPL_ALHL - MAPL_ALHF*SNO )
          DTS = DTS   / ( 1.0 + DTX*(BLW + SHD + EVD*MAPL_ALHL) )
          EVP = EVP + EVD * DTS
          SHF = SHF + SHD * DTS
          LHF = EVP * MAPL_ALHL
       else
          DTX = (DT / LAKEICECAP)
          SWN = (1.-ALBVRI)*VSUVR + (1.-ALBVFI)*VSUVF + &
                (1.-ALBNRI)*DRNIR + (1.-ALBNFI)*DFNIR 
          DTS = DTX * ( DTS + SWN - EVP*MAPL_ALHS )
          DTS = DTS   / ( 1.0 + DTX*(BLW + SHD + EVD*MAPL_ALHS) )
          EVP = EVP + EVD * DTS
          SHF = SHF + SHD * DTS
          LHF = EVP * MAPL_ALHS
          if(associated(SUBLIM)) SUBLIM = EVP*FR(:,N)
       end if

       RNF = PCU + PLS + SNO - EVP

! Update surface temperature and moisture
!----------------------------------------

       TS(:,N) = TS(:,N) + DTS
       DQS     = GEOS_QSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.) - QS(:,N)
       QS(:,N) = QS(:,N) + DQS  

! Exports
!--------

       if(associated(EVAPOUT)) EVAPOUT = EVAPOUT + EVP    *FR(:,N)
       if(associated(RUNOFF )) RUNOFF  = RUNOFF  + RNF    *FR(:,N)
       if(associated(SHOUT  )) SHOUT   = SHOUT   + SHF    *FR(:,N)
       if(associated(HLATN  )) HLATN   = HLATN   + LHF    *FR(:,N)
       if(associated(DELTS  )) DELTS   = DELTS   + DTS*CFT*FR(:,N)
       if(associated(DELQS  )) DELQS   = DELQS   + DQS*CFQ*FR(:,N)
       if(associated(SWNDSRF)) SWNDSRF = SWNDSRF + SWN    *FR(:,N)
       if(associated(TST    )) TST     = TST     + TS(:,N)*FR(:,N)
       if(associated(QST    )) QST     = QST     + QS(:,N)*FR(:,N)
       if(associated(LWNDSRF)) LWNDSRF = LWNDSRF - TS(:,N)*FR(:,N)*BLW
       if(associated(HLWUP  )) HLWUP   = HLWUP   + TS(:,N)*FR(:,N)*BLW

    end do

! Update Ice fraction
!--------------------
    do I=1,NT
       if    (TS(I,ICE)>MAPL_TICE .and. FR(I,ICE)>0.0) then
          ! MELT
          FR(I,WATER) = 1.0
          FR(I,ICE  ) = 0.0
          TS(I,WATER) = MAPL_TICE
          TS(I,ICE  ) = MAPL_TICE
       elseif(TS(I,WATER)<MAPL_TICE .and. FR(I,ICE)<1.0) then
          !FREEZE
          FR(I,WATER) = 0.0
          FR(I,ICE  ) = 1.0
          TS(I,WATER) = MAPL_TICE
          TS(I,ICE  ) = MAPL_TICE
       end if
    end do

! Update lake and ice albedos to anticipate
!   next radiation calculation
!-----------------------------------------------
        if(solalarmison) then
          call MAPL_SunGetInsolation(LONS, LATS,      &
            ORBIT, ZTH, SLR, &
            INTV = TINT,     &
            currTime=CURRENT_TIME+DELT,  &
            RC=STATUS )
          VERIFY_(STATUS)

          ZTH = max(0.0,ZTH)
          call ALBLAKE   (ALBVRO,ALBVFO,ALBNRO,ALBNFO,ZTH)
          call ALBLAKEICE(ALBVRI,ALBVFI,ALBNRI,ALBNFI,ZTH)
        endif

! Update albedos and emissivities to account for
!   ice cover changes
!-----------------------------------------------

    EMISS = LAKEEMISS*FR(:,WATER) + LAKEICEEMISS*FR(:,ICE)

    ALBVR = ALBVRO*FR(:,WATER) + ALBVRI*FR(:,ICE)
    ALBVF = ALBVFO*FR(:,WATER) + ALBVFI*FR(:,ICE)
    ALBNR = ALBNRO*FR(:,WATER) + ALBNRI*FR(:,ICE)
    ALBNF = ALBNFO*FR(:,WATER) + ALBNFI*FR(:,ICE)

!  All done
!-----------

    RETURN_(ESMF_SUCCESS)

  end subroutine LAKECORE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: ALBLAKE - Computes albedos as a function of $cos(\zeta)$ over lake surfaces

! !INTERFACE:

  subroutine ALBLAKE (ALBVR,ALBVF,ALBNR,ALBNF,ZTH)

! !ARGUMENTS:

    real,    intent(IN)  :: ZTH  (:)
    real,    intent(OUT) :: ALBVR(:) ! visible beam    albedo
    real,    intent(OUT) :: ALBVF(:) ! visible diffuse albedo
    real,    intent(OUT) :: ALBNR(:) ! nearIr  beam    albedo
    real,    intent(OUT) :: ALBNF(:) ! nearIr  diffuse albedo

!  !DESDRIPTION:
!        Compute albedo for ocean points
!          based on ceres

!  CERES ocean albedo at zth=.5 is 0.052. Our formulation gives .077
!    thus the scaling. The diffuse albedo is given by computing
!    the zth weighted average of the albedo over the hemisphere and
!    then applying the same scaling to match CERES.


!mjs: go back to old formulation 1-5-05
!    real, parameter :: CERESFAC   = 0.052/0.077
    real, parameter :: CERESFAC   = 1.0

    real, parameter :: OCNALBVF   = .08*CERESFAC
    real, parameter :: OCNALBNF   = .08*CERESFAC

    real, parameter :: A0         = 0.40670980*CERESFAC
    real, parameter :: A1         =-1.23236340*CERESFAC
    real, parameter :: A2         = 1.42240510*CERESFAC
    real, parameter :: A3         =-0.55573341*CERESFAC

! Beam albedos
!-------------

    ALBVR = A0+(A1+(A2+A3*ZTH)*ZTH)*ZTH
    ALBNR = ALBVR

! Diffuse albedos
!----------------

    ALBVF = OCNALBVF
    ALBNF = OCNALBNF

   RETURN_(ESMF_SUCCESS)
  end subroutine ALBLAKE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: ALBLAKEICE - Computes albedos as a function of  $cos(\zeta)$ over lake-ice surfaces

! !INTERFACE:

  subroutine ALBLAKEICE (ALBVR,ALBVF,ALBNR,ALBNF,ZTH)

! !ARGUMENTS:

    real,    intent(IN)  :: ZTH  (:)
    real,    intent(OUT) :: ALBVR(:) ! visible beam    albedo
    real,    intent(OUT) :: ALBVF(:) ! visible diffuse albedo
    real,    intent(OUT) :: ALBNR(:) ! nearIr  beam    albedo
    real,    intent(OUT) :: ALBNF(:) ! nearIr  diffuse albedo

!  !DESDRIPTION:
!        Compute albedo for ocean points
!          based on xxx. 

      real, parameter :: LAKEICEALBVR  = .60
      real, parameter :: LAKEICEALBVF  = .60
      real, parameter :: LAKEICEALBNR  = .60
      real, parameter :: LAKEICEALBNF  = .60

! Beam albedos
!-------------

      ALBVR = LAKEICEALBVR
      ALBNR = LAKEICEALBNR

! Diffuse albedos
!----------------

      ALBVF = LAKEICEALBVF
      ALBNF = LAKEICEALBNF

   RETURN_(ESMF_SUCCESS)
  end subroutine ALBLAKEICE


end subroutine RUN2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_LakeGridCompMod

