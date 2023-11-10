
!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_SaltwaterGridCompMod -- Implements slab saltwater tiles.

! !INTERFACE:

module GEOS_SaltwaterGridCompMod

! !DESCRIPTION:
! 
!   {\tt GEOS\_Saltwater} is a light-weight gridded component that updates
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

  use ESMF
  use MAPL
  use GEOS_UtilsMod
  use atmOcnIntlayer,     only: water_RHO

  use GEOS_OpenwaterGridCompMod,            only : OpenWaterSetServices       => SetServices
  use GEOS_SimpleSeaiceGridCompMod,         only : SimpleSeaiceSetServices    => SetServices
  use GEOS_CICE4ColumnPhysGridComp,         only : CICE4ColumnPhysSetServices => SetServices
  use GEOS_SeaiceInterfaceGridComp,         only : SeaiceInterfaceSetServices => SetServices
  use GEOS_ObioGridCompMod,                 only : ObioSetServices            => SetServices
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

  integer, parameter :: NUM_SUBTILES  = 2  ! number of subtiles for each tile (see above prologue)
  integer, parameter :: ICE           = 1  ! index(id) of two children fixed here 
  integer, parameter :: WATER         = 2  ! AddChild needs to adhere to the specification
  integer, parameter :: OBIO          = 3   

    type bandptr
      real, pointer, dimension(:) :: b
    end type bandptr

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
    integer                                 :: I, k
    integer                                 :: DO_OBIO         ! default (=0) is to run saltwater, with no ocean bio and chem
    integer                                 :: DO_CICE_THERMO  ! default (=0) is to run saltwater, with no LANL CICE Thermodynamics
    logical                                 :: seaIceT_extData ! default (=.FALSE.) is to NOT use data sea ice thickness from ExtData

    character(len = 2) :: suffix

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, _RC)
    Iam = trim(COMP_NAME) // 'SetServices'

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

! Sea-Ice Thermodynamics computation: using CICE or not?
!-------------------------------------------------------

    call MAPL_GetResource ( MAPL, DO_CICE_THERMO,     Label="USE_CICE_Thermo:" ,    DEFAULT=0, _RC)

! Ocean biology and chemistry: using OBIO or not?
!------------------------------------------------

    call MAPL_GetResource ( MAPL, DO_OBIO,            Label="USE_OCEANOBIOGEOCHEM:", DEFAULT=0, _RC)

! Data sea ice thickness from ExtData or not?
!--------------------------------------------

    call MAPL_GetResource ( MAPL,  seaIceT_extData, Label="SEAICE_THICKNESS_EXT_DATA:",  DEFAULT=.FALSE., _RC )

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, _RC)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run1, _RC)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2, _RC)

    ! order is important !!!
    ! sea-ice first and openwater second 
    ! changing order requires also changing indices of ICE and WATER (sub-tiles at the top)
    if(DO_CICE_THERMO == 0) then
       I = MAPL_AddChild(GC, NAME='SEAICETHERMO', SS=SimpleSeaiceSetServices,    RC=STATUS)
       VERIFY_(STATUS)
    elseif(DO_CICE_THERMO == 1) then 
       I = MAPL_AddChild(GC, NAME='SEAICETHERMO', SS=CICE4ColumnPhysSetServices, RC=STATUS)
       VERIFY_(STATUS)
    else
       I = MAPL_AddChild(GC, NAME='SEAICETHERMO', SS=SeaiceInterfaceSetServices, RC=STATUS)
       VERIFY_(STATUS)
    endif  

    I = MAPL_AddChild(GC,    NAME='OPENWATER', SS=OpenWaterSetServices,    _RC)

    if(DO_OBIO /= 0) then
       I = MAPL_AddChild(GC, NAME='OBIO', SS=ObioSetServices, _RC)
    endif  

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
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_visible_beam',   &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVR',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_visible_diffuse',&
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVF',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_near_infrared_beam', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNR',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_near_infrared_diffuse', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNF',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'evaporation'               ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'EVAPOUT'                   ,&
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
        LONG_NAME          = 'surface_emitted_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLWUP'                     ,&
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

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TST',                               &
        LONG_NAME          = 'surface_skin_temperature',          &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'QST',                               &
        LONG_NAME          = 'surface_specific_humidity',         &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TH',                                &
        LONG_NAME          = 'turbulence_surface_temperature',    &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'QH',                                &
        LONG_NAME          = 'turbulence_surface_specific_humidity', &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'UH',                                &
        LONG_NAME          = 'turbulence_surface_zonal_velocity', &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'VH',                                &
        LONG_NAME          = 'turbulence_surface_meridional_velocity', &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELTS',                             &
        LONG_NAME          = 'change_of_surface_skin_temperature',&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELQS',                             &
        LONG_NAME          = 'change_of_surface_specific_humidity',&
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CHT',                               &
        LONG_NAME          = 'surface_heat_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CMT',                               &
        LONG_NAME          = 'surface_momentum_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CQT',                               &
        LONG_NAME          = 'surface_moisture_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CNT',                               &
        LONG_NAME          = 'neutral_drag_coefficient',          &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'RIT',                               &
        LONG_NAME          = 'surface_bulk_richardson_number',    &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'RET',                               &
        LONG_NAME          = 'surface_reynolds_number',           &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)
 
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRACI',                             &
        LONG_NAME          = 'ice_covered_fraction_of_tile',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'FRACW',                             &
        LONG_NAME          = 'water_covered_fraction_of_tile',    &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC                     ,&
        SHORT_NAME         = 'SWFLX'                    ,&
        LONG_NAME          = 'surface_net_downward_shortwave_flux_at_ocean_surface',&
        UNITS              = 'W m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'GUST',                      &
        LONG_NAME          = 'gustiness',                 &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly,           &
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'VENT',                      &
        LONG_NAME          = 'surface_ventilation_velocity',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly,           &
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'Z0'                        ,&
        LONG_NAME          = 'surface_roughness'         ,&
        UNITS              = 'm'                         ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'Z0H'                       ,&
        LONG_NAME          = 'surface_roughness_for_heat',&
        UNITS              = 'm'                         ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOT2M',                     &
        LONG_NAME          = 'temperature 2m wind from MO sfc', &
        UNITS              = 'K',                         &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOQ2M',                     &
        LONG_NAME          = 'humidity 2m wind from MO sfc',    &
        UNITS              = 'kg kg-1',                   &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU2M',                    &
        LONG_NAME          = 'zonal 2m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV2M',                    &
        LONG_NAME          = 'meridional 2m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOT10M',                     &
        LONG_NAME          = 'temperature 10m wind from MO sfc', &
        UNITS              = 'K',                         &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOQ10M',                     &
        LONG_NAME          = 'humidity 10m wind from MO sfc',    &
        UNITS              = 'kg kg-1',                   &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU10M',                    &
        LONG_NAME          = 'zonal 10m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV10M',                    &
        LONG_NAME          = 'meridional 10m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU50M',                    &
        LONG_NAME          = 'zonal 50m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV50M',                    &
        LONG_NAME          = 'meridional 50m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               _RC)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'TAUXO'                     ,&
        LONG_NAME          = 'eastward_stress_on_ocean'  ,&
        UNITS              = 'N m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'TAUYO'                     ,&
        LONG_NAME          = 'northward_stress_on_ocean', &
        UNITS              = 'N m-2'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'OUSTAR3'                   ,&
        LONG_NAME          = 'ocean_ustar_cubed',         &
        UNITS              = 'm+3 s-3'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'UU'                        ,&
        LONG_NAME          = 'surface_wind_speed',        &
        UNITS              = 'm s-1'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

! Following export of DISCHARGE is here in saltwater only for the sake of 
! "passing thru" from atmosphere to ocean, no computation is otherwise done with (on) them.
     call MAPL_AddExportSpec(GC,                            &      
          SHORT_NAME         = 'DISCHARGE'                 ,&
          LONG_NAME          = 'river_discharge_at_ocean_points',& 
          UNITS              = 'kg m-2 s-1'                ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PS',                                &
        LONG_NAME          = 'surface_pressure',                  &
        UNITS              = 'Pa',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

! Following 3 exports (HFLUX, WATERFLUX, SALTFLUX) are for ocean model - need to be filled up.
   call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'HFLUX',                            &
        LONG_NAME          = 'heat_flux_below_saltwater_ocean',  &
        UNITS              = 'W m-2',                            &
        DIMS               = MAPL_DimsTileOnly,                  &
        VLOCATION          = MAPL_VLocationNone,                 &
                                                       _RC  )

   call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'WATERFLUX',                        &
        LONG_NAME          = 'water_flux_below_saltwater_ocean', &
        UNITS              = 'kg m-2 s-1',                       &
        DIMS               = MAPL_DimsTileOnly,                  &
        VLOCATION          = MAPL_VLocationNone,                 &
                                                       _RC  )

   call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'SALTFLUX',                         &
        LONG_NAME          = 'salt_flux_below_saltwater_ocean',  &
        UNITS              = 'kg m-2 s-1',                       &
        DIMS               = MAPL_DimsTileOnly,                  &
        VLOCATION          = MAPL_VLocationNone,                 &
                                                       _RC  )

   call MAPL_AddExportSpec(GC                    ,&
          SHORT_NAME         = 'PENUVF',                    &
          LONG_NAME          = 'downwelling_uvr_diffuse_flux_at_skin_base',&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  ) 

   call MAPL_AddExportSpec(GC                    ,&
          SHORT_NAME         = 'PENUVR',                    &
          LONG_NAME          = 'downwelling_uvr_direct_flux_at_skin_base',&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  ) 

   call MAPL_AddExportSpec(GC                    ,&
          SHORT_NAME         = 'PENPAF',                    &
          LONG_NAME          = 'downwelling_par_diffuse_flux_at_skin_base',&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  ) 

   call MAPL_AddExportSpec(GC                    ,&
          SHORT_NAME         = 'PENPAR',                    &
          LONG_NAME          = 'downwelling_par_direct_flux_at_skin_base',&
          UNITS              = 'W m-2'                     ,&
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

!  !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'UU',                                &
        LONG_NAME          = 'surface_wind_speed',                &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PS',                                &
        LONG_NAME          = 'surface_pressure',                  &
        UNITS              = 'Pa',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)
    
  ! Following import is here in saltwater only for the sake of 
  ! "passing thru" from atmosphere to ocean, no computation is otherwise done with (on) them.
    call MAPL_AddImportSpec(GC,                    &
          SHORT_NAME         = 'DISCHARGE'                 ,&
          LONG_NAME          = 'river_discharge_at_ocean_points',&
          UNITS              = 'kg m-2 s-1'                ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartSkip            ,&
          _RC  ) 

  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSKINW'    , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'HSKINW'    , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'SSKINW'    , CHILD_ID = WATER, _RC)  

  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSKINI'    , CHILD_ID = ICE,   _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'SSKINI'    , CHILD_ID = ICE,   _RC)
  if (DO_CICE_THERMO <= 1) then 
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'HSKINI'    , CHILD_ID = ICE,   _RC)
  endif
     
  if(DO_CICE_THERMO == 1) then ! additional exports from CICE4 sea ice thermodynamics
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'FR'     , CHILD_ID =   ICE, _RC)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'VOLICE' , CHILD_ID =   ICE, _RC)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'VOLSNO' , CHILD_ID =   ICE, _RC)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ERGICE' , CHILD_ID =   ICE, _RC)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ERGSNO' , CHILD_ID =   ICE, _RC)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'VOLPOND', CHILD_ID =   ICE, _RC)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TAUAGE' , CHILD_ID =   ICE, _RC)
  else
    if (seaIceT_extData) then
      call MAPL_AddExportSpec(GC, SHORT_NAME = 'SEAICETHICKNESS', CHILD_ID =   ICE, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME = 'HSNO',            CHILD_ID =   ICE, _RC)
    endif
  endif 

  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TAUXW'     , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TAUYW'     , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'HLATWTR'   , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'SWNDWTR'   , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'LWNDWTR'   , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'SHWTR'     , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'SNOWOCN'   , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'RAINOCN'   , CHILD_ID = WATER, _RC)

! Atmosphere-Ocean Interface Layer (AOIL) specific variables
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'DCOOL'     , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'DWARM'     , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TDROP'     , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'QCOOL'     , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'QWARM'     , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'SWCOOL'    , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'SWWARM'    , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'PHIW'      , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'LANGM'     , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'USTARW'    , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TBAR'      , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'LCOOL'     , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'BCOOL'     , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TDEL'      , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TS_FOUND'  , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'SS_FOUND'  , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TAUTW'     , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'ZETA_W'    , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TWMTF'     , CHILD_ID = WATER, _RC)

  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TAUXI'     , CHILD_ID =   ICE, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TAUYI'     , CHILD_ID =   ICE, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'SUBLIM'    , CHILD_ID =   ICE, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'HLATICE'   , CHILD_ID =   ICE, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'SWNDICE'   , CHILD_ID =   ICE, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'LWNDICE'   , CHILD_ID =   ICE, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'SHICE'     , CHILD_ID =   ICE, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'GHTSKIN'   , CHILD_ID =   ICE, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSKINICE'  , CHILD_ID =   ICE, _RC)

  if (DO_CICE_THERMO <= 1) then 
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'FRESH'     , CHILD_ID =   ICE, _RC)
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'FSALT'     , CHILD_ID =   ICE, _RC)
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'FHOCN'     , CHILD_ID =   ICE, _RC)
  endif

! Atmosphere-Ocean Fluxes
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'AO_SHFLX'  , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'AO_QFLUX'  , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'AO_LWFLX'  , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'AO_SNOW'   , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'AO_RAIN'   , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'AO_DRNIR'  , CHILD_ID = WATER, _RC)
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'AO_DFNIR'  , CHILD_ID = WATER, _RC)

  if(DO_OBIO /= 0) then
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'CO2SC'     , CHILD_ID = OBIO, __RC__)
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'DUDP'      , CHILD_ID = OBIO, __RC__)
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'DUWT'      , CHILD_ID = OBIO, __RC__)
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'DUSD'      , CHILD_ID = OBIO, __RC__)
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'BCDP'      , CHILD_ID = OBIO, __RC__)
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'BCWT'      , CHILD_ID = OBIO, __RC__)
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'OCDP'      , CHILD_ID = OBIO, __RC__)
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'OCWT'      , CHILD_ID = OBIO, __RC__)
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'FSWBAND'   , CHILD_ID = OBIO, __RC__)  ! Delete? It's not imported by OBIO, EMS
     call MAPL_AddExportSpec(GC, SHORT_NAME = 'FSWBANDNA' , CHILD_ID = OBIO, __RC__)  ! Delete? It's not imported by OBIO, EMS

!    ! We are setting these variable values for OBIORAD temporarily. These variables should be read from
!    ! file or computed.
     do k=1, 33
        write(unit = suffix, fmt = '(i2.2)') k
        call MAPL_AddExportSpec(GC,                               &
             SHORT_NAME = 'TAUA_'//suffix,                        &
             LONG_NAME  = 'aerosol optical thickness',            &
             UNITS      = '',                                     &
             DIMS       = MAPL_DimsTileOnly,                      &
             VLOCATION  = MAPL_VLocationNone,                     &
             default    = 1.0, &
             __RC__)

        call MAPL_AddExportSpec(GC,                               &
             SHORT_NAME = 'ASYMP_'//suffix,                       &
             LONG_NAME  = 'asymmetry parameter',                  &
             UNITS      = '',                                     &
             DIMS       = MAPL_DimsTileOnly,                      &
             VLOCATION  = MAPL_VLocationNone,                     &
             default    = 0.0, &
             __RC__)

        call MAPL_AddExportSpec(GC,                               &
             SHORT_NAME = 'SSALB_'//suffix,                       &
             LONG_NAME  = 'single scattering albedo',             &
             UNITS      = '',                                     &
             DIMS       = MAPL_DimsTileOnly,                      &
             VLOCATION  = MAPL_VLocationNone,                     &
             default    = 0.95, &
             __RC__)
     enddo

     call MAPL_AddExportSpec(GC,                               &
          SHORT_NAME = 'RH',                                   &
          LONG_NAME  = 'relative humidity',             &
          UNITS      = 'percent',                                     &
          DIMS       = MAPL_DimsTileOnly,                      &
          VLOCATION  = MAPL_VLocationNone,                     &
          default    = 1.0, &
          __RC__)

     call MAPL_AddExportSpec(GC,                               &
          SHORT_NAME = 'CCOVM',                                &
          LONG_NAME  = 'cloud cover',                          &
          UNITS      = 'fraction (dimensionless)',             &
          DIMS       = MAPL_DimsTileOnly,                      &
          VLOCATION  = MAPL_VLocationNone,                     &
          default    = 1.0, &
          __RC__)

     call MAPL_AddExportSpec(GC,                               &
          SHORT_NAME = 'CDREM',                                   &
          LONG_NAME  = 'cloud droplet effective radius',             &
          UNITS      = '',                                     &
          DIMS       = MAPL_DimsTileOnly,                      &
          VLOCATION  = MAPL_VLocationNone,                     &
          default    = 1.0, &
          __RC__)
     
     call MAPL_AddExportSpec(GC,                               &
          SHORT_NAME = 'RLWPM',                                &
          LONG_NAME  = 'cloud liquid water path',              &
          UNITS      = '',                                     &
          DIMS       = MAPL_DimsTileOnly,                      &
          VLOCATION  = MAPL_VLocationNone,                     &
          default    = 1.0, &
          __RC__)

     call MAPL_AddExportSpec(GC,                               &
          SHORT_NAME = 'CLDTCM',                               &
          LONG_NAME  = 'cloud optical thickness',              &
          UNITS      = '',                                     &
          DIMS       = MAPL_DimsTileOnly,                      &
          VLOCATION  = MAPL_VLocationNone,                     &
          default    = 1.0, &
          __RC__)

     call MAPL_AddExportSpec(GC,                               &
          SHORT_NAME = 'OZ',                                   &
          LONG_NAME  = 'ozone thickness',                      &
          UNITS      = 'Dobson units',                         &
          DIMS       = MAPL_DimsTileOnly,                      &
          VLOCATION  = MAPL_VLocationNone,                     &
          default    = 1.0, &
          __RC__)

     call MAPL_AddExportSpec(GC,                               &
          SHORT_NAME = 'WV',                                   &
          LONG_NAME  = 'water vapor',                          &
          UNITS      = 'cm',                                   &
          DIMS       = MAPL_DimsTileOnly,                      &
          VLOCATION  = MAPL_VLocationNone,                     &
          default    = 1.0, &
          __RC__)

  end if


! and that penetrated below ocean model first layer
  call MAPL_AddExportSpec(GC, SHORT_NAME = 'PEN_OCN' ,   CHILD_ID = WATER, _RC)

!EOS

  call MAPL_AddConnectivity ( GC,   &
       SHORT_NAME  = (/'TSKINW','SSKINW'/),  &
       DST_ID = ICE,                &
       SRC_ID = WATER,              &
       _RC  )

  call MAPL_AddConnectivity ( GC,   &
       SHORT_NAME  = [character(len=8) :: 'FRACI', 'FRACINEW','TFREEZE'],     & 
       DST_ID = WATER,              &
       SRC_ID = ICE,                &
       _RC  )

! Set the Profiling timers
! ------------------------
    call MAPL_TimerAdd(GC,    name="INITIALIZE"     ,       _RC)
    call MAPL_TimerAdd(GC,    name="RUN1"   ,               _RC)
    call MAPL_TimerAdd(GC,    name="RUN2"  ,                _RC)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC,  _RC)
 
! Set the Run entry point
! -----------------------

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! !IROUTINE: Initialize -- Initialize method for the composite Surface Gridded Component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Initialize method of the Land Composite Gridded Component.
 

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)           :: IAm 
    integer                              :: STATUS
    character(len=ESMF_MAXSTR)           :: COMP_NAME
    
! Local derived type aliases

    type (MAPL_MetaComp    ), pointer    :: MAPL
    type (MAPL_MetaComp    ), pointer    :: CHILD_MAPL 
    type (MAPL_LocStream   )             :: LOCSTREAM
    type (ESMF_DELayout    )             :: LAYOUT
    type (ESMF_Config      )             :: CF
    type (ESMF_GridComp    ), pointer    :: GCS(:)
  
    integer                              :: I
!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, _RC)
    Iam = trim(COMP_NAME) // "Initialize"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

    call MAPL_TimerOn(MAPL,"INITIALIZE", _RC )
    call MAPL_TimerOn(MAPL,"TOTAL",      _RC )

! Get the ocean tilegrid and the child components
!------------------------------------------------

    call MAPL_Get (MAPL, LOCSTREAM=LOCSTREAM, GCS=GCS, _RC)

! Place the tilegrid in the generic state of each child component
!----------------------------------------------------------------

    do I = 1, SIZE(GCS)
       call MAPL_GetObjectFromGC( GCS(I), CHILD_MAPL, _RC )
       call MAPL_Set (CHILD_MAPL, LOCSTREAM=LOCSTREAM, _RC )
    end do

    call MAPL_TimerOff(MAPL,"TOTAL", _RC )

! Call Initialize for every Child
!--------------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  _RC)

    call MAPL_TimerOff(MAPL,"INITIALIZE", _RC )

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

  type (MAPL_MetaComp), pointer       :: MAPL => null()
  type (ESMF_GridComp), pointer       :: GCS(:)
  type (ESMF_State),    pointer       :: GIM(:)
  type (ESMF_State),    pointer       :: GEX(:)
  character(len=ESMF_MAXSTR),pointer  :: GCnames(:)
  type (ESMF_Config)                  :: CF

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

! pointers to childrens' export

   real, pointer, dimension(:  )  :: FI    => null()
   real, pointer, dimension(:  )  :: QS    => null()
   real, pointer, dimension(:  )  :: CH    => null()
   real, pointer, dimension(:  )  :: CM    => null()
   real, pointer, dimension(:  )  :: CQ    => null()
   real, pointer, dimension(:  )  :: Z0    => null()
   real, pointer, dimension(:  )  :: TS    => null() 
   real, pointer, dimension(:  )  :: US    => null() 
   real, pointer, dimension(:  )  :: VS    => null()
   real, pointer, dimension(:  )  :: CN    => null()
   real, pointer, dimension(:  )  :: RE    => null()
   real, pointer, dimension(:  )  :: ZT    => null()
   real, pointer, dimension(:  )  :: ZQ    => null()
   real, pointer, dimension(:  )  :: U50M  => null()
   real, pointer, dimension(:  )  :: V50M  => null()
   real, pointer, dimension(:  )  :: T10M  => null()
   real, pointer, dimension(:  )  :: Q10M  => null()
   real, pointer, dimension(:  )  :: U10M  => null()
   real, pointer, dimension(:  )  :: V10M  => null()
   real, pointer, dimension(:  )  :: T2M   => null()
   real, pointer, dimension(:  )  :: Q2M   => null()
   real, pointer, dimension(:  )  :: U2M   => null()
   real, pointer, dimension(:  )  :: V2M   => null()
   real, pointer, dimension(:  )  :: THO   => null()
   real, pointer, dimension(:  )  :: QHO   => null()
   real, pointer, dimension(:  )  :: RI    => null()
   real, pointer, dimension(:  )  :: VN    => null()
   real, pointer, dimension(:  )  :: GS    => null()
   real, pointer, dimension(:  )  :: LONS  => null()
   real, pointer, dimension(:  )  :: dummy => null()

   integer                        :: I, N, NT

   real, allocatable              :: FR (:,:)        ! note: rank-2
   real, allocatable              :: CHB(:)
   real, allocatable              :: CQB(:)
   real, allocatable              :: CMB(:)
   real, allocatable              :: UCN(:)

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run1"
    call ESMF_GridCompGet( GC, name=COMP_NAME, _RC)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

    call MAPL_Get(MAPL,                &
                  TILELONS  = LONS ,   &
                  _RC)

    call MAPL_Get (MAPL, GCS=GCS, GIM=GIM, GEX=GEX, GCnames=GCnames,_RC)

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN1" )

! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,QH    , 'QH'      ,    _RC)
   call MAPL_GetPointer(EXPORT,TH    , 'TH'      ,    _RC)
   call MAPL_GetPointer(EXPORT,UH    , 'UH'      ,    _RC)
   call MAPL_GetPointer(EXPORT,VH    , 'VH'      ,    _RC)
   call MAPL_GetPointer(EXPORT,QST   , 'QST'     ,    _RC)
   call MAPL_GetPointer(EXPORT,TST   , 'TST'     ,    _RC)
   call MAPL_GetPointer(EXPORT,CHT   , 'CHT'     ,    _RC)
   call MAPL_GetPointer(EXPORT,CMT   , 'CMT'     ,    _RC)
   call MAPL_GetPointer(EXPORT,CQT   , 'CQT'     ,    _RC)
   call MAPL_GetPointer(EXPORT,CNT   , 'CNT'     ,    _RC)
   call MAPL_GetPointer(EXPORT,RIT   , 'RIT'     ,    _RC)
   call MAPL_GetPointer(EXPORT,RET   , 'RET'     ,    _RC)
   call MAPL_GetPointer(EXPORT,Z0O   , 'Z0'      ,    _RC)
   call MAPL_GetPointer(EXPORT,Z0H   , 'Z0H'     ,    _RC)
   call MAPL_GetPointer(EXPORT,MOT2M, 'MOT2M'    ,    _RC)
   call MAPL_GetPointer(EXPORT,MOQ2M, 'MOQ2M'    ,    _RC)
   call MAPL_GetPointer(EXPORT,MOU2M, 'MOU2M'    ,    _RC)
   call MAPL_GetPointer(EXPORT,MOV2M, 'MOV2M'    ,    _RC)
   call MAPL_GetPointer(EXPORT,MOT10M, 'MOT10M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOQ10M, 'MOQ10M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOU10M, 'MOU10M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOV10M, 'MOV10M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOU50M, 'MOU50M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOV50M, 'MOV50M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,GST   , 'GUST'    ,    _RC)
   call MAPL_GetPointer(EXPORT,VNT   , 'VENT'    ,    _RC)

!  Retrieve pointers to exports from all children except OBIO
   do I = 1, size(GCS)
      if (trim(GCnames(I)) .ne. 'OBIO') then

         if(associated(MOT2M)) then
            call MAPL_GetPointer(GEX(I), dummy, 'MOT2M'  , alloc=.true., _RC)
         endif 
         if(associated(MOQ2M)) then
            call MAPL_GetPointer(GEX(I), dummy, 'MOQ2M'  , alloc=.true., _RC)
         endif 
         if(associated(MOU2M)) then
            call MAPL_GetPointer(GEX(I), dummy, 'MOU2M'  , alloc=.true., _RC)
         endif 
         if(associated(MOV2M)) then
            call MAPL_GetPointer(GEX(I), dummy, 'MOV2M'  , alloc=.true., _RC)
         endif 
         if(associated(MOT10M)) then
            call MAPL_GetPointer(GEX(I), dummy, 'MOT10M' , alloc=.true., _RC)
         endif 
         if(associated(MOQ10M)) then
            call MAPL_GetPointer(GEX(I), dummy, 'MOQ10M' , alloc=.true., _RC)
         endif 
         if(associated(MOU10M)) then
            call MAPL_GetPointer(GEX(I), dummy, 'MOU10M' , alloc=.true., _RC)
         endif 
         if(associated(MOV10M)) then
            call MAPL_GetPointer(GEX(I), dummy, 'MOV10M' , alloc=.true., _RC)
         endif 
         if(associated(MOU50M)) then
            call MAPL_GetPointer(GEX(I), dummy, 'MOU50M' , alloc=.true., _RC)
         endif 
         if(associated(MOV50M)) then
            call MAPL_GetPointer(GEX(I), dummy, 'MOV50M' , alloc=.true., _RC)
         endif 
         if(associated(TH)) then
            call MAPL_GetPointer(GEX(I), dummy, 'TH'     , alloc=.true., _RC)
         endif 
         if(associated(QH)) then
            call MAPL_GetPointer(GEX(I), dummy, 'QH'     , alloc=.true., _RC)
         endif 
         if(associated(UH)) then
            call MAPL_GetPointer(GEX(I), dummy, 'UH'     , alloc=.true., _RC)
         endif 
         if(associated(VH)) then
            call MAPL_GetPointer(GEX(I), dummy, 'VH'     , alloc=.true., _RC)
         endif 
         if(associated(CHT)) then
            call MAPL_GetPointer(GEX(I), dummy, 'CHT'    , alloc=.true., _RC)
         endif 
         if(associated(CQT)) then
            call MAPL_GetPointer(GEX(I), dummy, 'CQT'    , alloc=.true., _RC)
         endif 
         if(associated(CMT)) then
            call MAPL_GetPointer(GEX(I), dummy, 'CMT'    , alloc=.true., _RC)
         endif 
         if(associated(TST)) then
            call MAPL_GetPointer(GEX(I), dummy, 'TST'    , alloc=.true., _RC)
         endif 
         if(associated(QST)) then
            call MAPL_GetPointer(GEX(I), dummy, 'QST'    , alloc=.true., _RC)
         endif 
         if(associated(CNT)) then
            call MAPL_GetPointer(GEX(I), dummy, 'CNT'    , alloc=.true., _RC)
         endif 
         if(associated(CNT)) then
            call MAPL_GetPointer(GEX(I), dummy, 'RIT'    , alloc=.true., _RC)
         endif 
         if(associated(CNT)) then
            call MAPL_GetPointer(GEX(I), dummy, 'RET'    , alloc=.true., _RC)
         endif 
         if(associated(Z0O)) then
            call MAPL_GetPointer(GEX(I), dummy, 'Z0'     , alloc=.true., _RC)
         endif 
         if(associated(Z0H)) then
            call MAPL_GetPointer(GEX(I), dummy, 'Z0H'    , alloc=.true., _RC)
         endif 
         if(associated(GST)) then
            call MAPL_GetPointer(GEX(I), dummy, 'GUST'   , alloc=.true., _RC)
         endif 
         if(associated(VNT)) then
            call MAPL_GetPointer(GEX(I), dummy, 'VENT'   , alloc=.true., _RC)
        endif 
      endif 
   enddo

   NT = size(LONS)

   allocate(CHB(NT)  ,   _STAT)
   allocate(CQB(NT)  ,   _STAT)
   allocate(CMB(NT)  ,   _STAT)
   allocate(UCN(NT)  ,   _STAT)
   allocate(FR (NT,NUM_SUBTILES),_STAT) ! note: rank-2

!  Clear the output tile accumulators
!------------------------------------

                       CHB = 0.0
                       CQB = 0.0
                       CMB = 0.0
                       UCN = 0.0
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

! Call the childrens' RUN1
!-------------------------

   DO I = 1, size(GCS)
       call MAPL_TimerOn(MAPL,trim(GCnames(i)), _RC )

       call ESMF_GridCompRun(GCS(I), importState=GIM(I), exportState=GEX(I), &
                             CLOCK=CLOCK, PHASE=1, userRC=STATUS)
       VERIFY_(STATUS)

       call MAPL_TimerOff(MAPL,trim(GCnames(i)), _RC )
   ENDDO

   call MAPL_GetPointer(GEX(ICE), FI, 'FRACI'  , _RC)

   ! make sure following fractions are bounded in [0.0, 1.0]
   FR(:,WATER) = max(1.0-FI, 0.0)
   FR(:,ICE  ) = min(FI,     1.0)

   ! aggregate over ice and water subtiles
   SUB_TILES: do N=1,NUM_SUBTILES

         call MAPL_GetPointer(GEX(N), U50M, 'MOU50M'  , _RC)
         call MAPL_GetPointer(GEX(N), V50M, 'MOV50M'  , _RC)
         call MAPL_GetPointer(GEX(N), T10M, 'MOT10M'  , _RC)
         call MAPL_GetPointer(GEX(N), Q10M, 'MOQ10M'  , _RC)
         call MAPL_GetPointer(GEX(N), U10M, 'MOU10M'  , _RC)
         call MAPL_GetPointer(GEX(N), V10M, 'MOV10M'  , _RC)
         call MAPL_GetPointer(GEX(N), T2M , 'MOT2M'   , _RC)
         call MAPL_GetPointer(GEX(N), Q2M , 'MOQ2M'   , _RC)
         call MAPL_GetPointer(GEX(N), U2M , 'MOU2M'   , _RC)
         call MAPL_GetPointer(GEX(N), V2M , 'MOV2M'   , _RC)
         call MAPL_GetPointer(GEX(N), CH  , 'CHT'     , _RC)
         call MAPL_GetPointer(GEX(N), CQ  , 'CQT'     , _RC)
         call MAPL_GetPointer(GEX(N), CM  , 'CMT'     , _RC)
         call MAPL_GetPointer(GEX(N), THO , 'TH'      , _RC)
         call MAPL_GetPointer(GEX(N), QHO , 'QH'      , _RC)
         call MAPL_GetPointer(GEX(N), US  , 'UH'      , _RC)
         call MAPL_GetPointer(GEX(N), VS  , 'VH'      , _RC)
         call MAPL_GetPointer(GEX(N), TS  , 'TST'     , _RC)
         call MAPL_GetPointer(GEX(N), QS  , 'QST'     , _RC)
         call MAPL_GetPointer(GEX(N), CN  , 'CNT'     , _RC)
         call MAPL_GetPointer(GEX(N), RI  , 'RIT'     , _RC)
         call MAPL_GetPointer(GEX(N), RE  , 'RET'     , _RC)
         call MAPL_GetPointer(GEX(N), Z0  , 'Z0'      , _RC)
         call MAPL_GetPointer(GEX(N), ZT  , 'Z0H'     , _RC)
         call MAPL_GetPointer(GEX(N), GS  , 'GUST'    , _RC)
         call MAPL_GetPointer(GEX(N), VN  , 'VENT'    , _RC)

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

      !  Aggregate to tiles
      !--------------------

                             CHB     = CHB + CH(:  )*FR(:,N)
                             CQB     = CQB + CQ(:  )*FR(:,N)
                             CMB     = CMB + CM(:  )*FR(:,N)
         if(associated(TST)) TST     = TST + TS(:  )*FR(:,N)
         if(associated(QST)) QST     = QST + QS(:  )*FR(:,N)
         if(associated(CNT)) CNT     = CNT + CN(:  )*FR(:,N)
         if(associated(RIT)) RIT     = RIT + RI(:  )*FR(:,N)
         if(associated(RET)) RET     = RET + RE(:  )*FR(:,N)
         if(associated(Z0O)) Z0O     = Z0O + Z0(:  )*FR(:,N)
         if(associated(Z0H)) Z0H     = Z0H + ZT(:  )*FR(:,N)
         if(associated(GST)) GST     = GST + GS(:  )*FR(:,N)
         if(associated(VNT)) VNT     = VNT + VN(:  )*FR(:,N)

      !  Aggregate effective, CD-weighted, surface values of T and Q
      !-------------------------------------------------------------

         if(associated(TH)) TH      = TH  + THO(:)*FR(:,N)
         if(associated(QH)) QH      = QH  + QHO(:)*FR(:,N)
         if(associated(UH)) UH      = UH  +  US(:)*FR(:,N)
         if(associated(VH)) VH      = VH  +  VS(:)*FR(:,N)

   end do SUB_TILES

   if(associated(TH )) TH  = TH /CHB
   if(associated(QH )) QH  = QH /CQB
   if(associated(UH )) UH  = UH /CMB
   if(associated(VH )) VH  = VH /CMB
   if(associated(CHT)) CHT = CHB
   if(associated(CQT)) CQT = CQB
   if(associated(CMT)) CMT = CMB
   if(associated(GST)) GST = sqrt(max(GST+UCN,0.0))

   deallocate(CHB)
   deallocate(CQB)
   deallocate(CMB)
   deallocate(UCN)
   deallocate(FR )

!  All done with RUN1
!--------------------

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
  type (ESMF_Config      )            :: CF
  type (ESMF_GridComp), pointer       :: GCS(:)
  type (ESMF_State),    pointer       :: GIM(:)
  type (ESMF_State),    pointer       :: GEX(:)
  character(len=ESMF_MAXSTR),pointer  :: GCnames(:)
  type (ESMF_State   )                :: CHILD_INTERNAL
  type (MAPL_MetaComp), pointer       :: CHILD_MAPL => null()

  real, pointer, dimension(:)         :: LATS => null()
  real, pointer, dimension(:)         :: LONS => null()
  real, pointer, dimension(:)         :: AREA => null()     ! needed to calculate TILEAREA in SaltWaterCore

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run2"
    call ESMF_GridCompGet( GC, name=COMP_NAME, _RC)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN2" )

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,             &
         TILELATS  = LATS ,         &
         TILELONS  = LONS ,         &
         CF = CF,                   &
    _RC)

! Update the skin variables at each step
!---------------------------------------

    call SALTWATERCORE(NT=size(LONS), _RC)

!  All done with RUN2
!--------------------

   call MAPL_TimerOff(MAPL,"RUN2" )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine SALTWATERCORE(NT,RC)
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
   real, pointer, dimension(:  )  :: SHOUT   => null()
   real, pointer, dimension(:  )  :: HLATN   => null()
   real, pointer, dimension(:  )  :: HLWUP   => null()
   real, pointer, dimension(:  )  :: SWNDSRF => null()
   real, pointer, dimension(:  )  :: LWNDSRF => null()
   real, pointer, dimension(:  )  :: PENUVR  => null()
   real, pointer, dimension(:  )  :: PENUVF  => null()
   real, pointer, dimension(:  )  :: PENPAR  => null()
   real, pointer, dimension(:  )  :: PENPAF  => null()
   real, pointer, dimension(:  )  :: FSURF   => null()
   real, pointer, dimension(:  )  :: SWFLX   => null()

   real, pointer, dimension(:  )  :: DELTS  => null()
   real, pointer, dimension(:  )  :: DELQS  => null()
   real, pointer, dimension(:  )  :: TST    => null()
   real, pointer, dimension(:  )  :: QST    => null()
   real, pointer, dimension(:  )  :: TAUXO  => null()
   real, pointer, dimension(:  )  :: TAUYO  => null()
   real, pointer, dimension(:  )  :: USTR3  => null()
   real, pointer, dimension(:  )  :: UUEX   => null()
   real, pointer, dimension(:  )  :: PSEX   => null()
   real, pointer, dimension(:  )  :: TSKINI => null()
   real, pointer, dimension(:  )  :: FRI    => null()
   real, pointer, dimension(:  )  :: FRW    => null()

   real, pointer, dimension(:)    :: DISCHARGE   => null()

! pointers to import

   real, pointer, dimension(:)    :: PS => null()
   real, pointer, dimension(:)    :: UU => null()
   real, pointer, dimension(:)    :: FI => null()
   real, pointer, dimension(:)    :: DISCHARGE_IM => null()

! pointers to the childrens' export
   real, pointer, dimension(:)    :: TS        => null() 
   real, pointer, dimension(:)    :: QS        => null()
   real, pointer, dimension(:)    :: EMS       => null()
   real, pointer, dimension(:)    :: AVR       => null()
   real, pointer, dimension(:)    :: AVF       => null()
   real, pointer, dimension(:)    :: ANR       => null()
   real, pointer, dimension(:)    :: ANF       => null()
   real, pointer, dimension(:)    :: SHF       => null()
   real, pointer, dimension(:)    :: EVP       => null()
   real, pointer, dimension(:)    :: LHF       => null()
   real, pointer, dimension(:)    :: DTS       => null()
   real, pointer, dimension(:)    :: DQS       => null()
   real, pointer, dimension(:)    :: TXW       => null()
   real, pointer, dimension(:)    :: TYW       => null()
   real, pointer, dimension(:)    :: TXI       => null()
   real, pointer, dimension(:)    :: TYI       => null() 
   real, pointer, dimension(:)    :: HLW       => null()
   real, pointer, dimension(:)    :: LWND      => null()
   real, pointer, dimension(:)    :: SWND      => null()
   real, pointer, dimension(:)    :: PUR       => null()
   real, pointer, dimension(:)    :: PUF       => null()
   real, pointer, dimension(:)    :: PAR       => null()
   real, pointer, dimension(:)    :: PAF       => null()
   real, pointer, dimension(:)    :: FSUR      => null()
   real, pointer, dimension(:)    :: dummy     => null()

   real,    dimension(NT,NUM_SUBTILES) :: FR     ! note: rank-2
   real,    dimension(NT,NUM_SUBTILES) :: FRNEW  ! note: rank-2

   real,    dimension(NT)              :: TXO
   real,    dimension(NT)              :: TYO

   integer                             :: N
   integer                             :: NSUB, I, K, L
   integer                             :: DO_CICE_THERMO

!  -------------------------------------------------------------------

!  Begin...
!----------

   IAm =  trim(COMP_NAME) // "SALTWATERCORE"

   call MAPL_GetResource ( MAPL, DO_CICE_THERMO,     Label="USE_CICE_Thermo:" ,    DEFAULT=0, _RC)

! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,PS     , 'PS'     ,    _RC)
   call MAPL_GetPointer(IMPORT,UU     , 'UU'     ,    _RC)
!!$   if (DO_GUEST /= 0) then    
      call MAPL_GetPointer(IMPORT, DISCHARGE_IM, 'DISCHARGE', _RC)
!!$   endif

! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,EMISS  , 'EMIS' , alloc=.true., _RC)
   call MAPL_GetPointer(EXPORT,ALBVF  , 'ALBVF', alloc=.true., _RC)
   call MAPL_GetPointer(EXPORT,ALBVR  , 'ALBVR', alloc=.true., _RC)
   call MAPL_GetPointer(EXPORT,ALBNF  , 'ALBNF', alloc=.true., _RC)
   call MAPL_GetPointer(EXPORT,ALBNR  , 'ALBNR', alloc=.true., _RC)

   call MAPL_GetPointer(EXPORT,QST    , 'QST'     ,    _RC)
   call MAPL_GetPointer(EXPORT,TST    , 'TST'     ,    _RC)
   call MAPL_GetPointer(EXPORT,DELTS  , 'DELTS'   ,    _RC)
   call MAPL_GetPointer(EXPORT,DELQS  , 'DELQS'   ,    _RC)
   call MAPL_GetPointer(EXPORT,TAUXO  , 'TAUXO'   ,    _RC)
   call MAPL_GetPointer(EXPORT,TAUYO  , 'TAUYO'   ,    _RC)
   call MAPL_GetPointer(EXPORT,USTR3  , 'OUSTAR3' ,    _RC)
   call MAPL_GetPointer(EXPORT,UUEX   , 'UU'      ,    _RC)
   call MAPL_GetPointer(EXPORT,PSEX   , 'PS'      ,    _RC)
   call MAPL_GetPointer(EXPORT,EVAPOUT, 'EVAPOUT' ,    _RC)
   call MAPL_GetPointer(EXPORT,SHOUT  , 'SHOUT'   ,    _RC)
   call MAPL_GetPointer(EXPORT,HLATN  , 'HLATN'   ,    _RC)
   call MAPL_GetPointer(EXPORT,HLWUP  , 'HLWUP'   ,    _RC)
   call MAPL_GetPointer(EXPORT,LWNDSRF, 'LWNDSRF' ,    _RC)
   call MAPL_GetPointer(EXPORT,SWNDSRF, 'SWNDSRF' ,    _RC)
   call MAPL_GetPointer(EXPORT,FRI    , 'FRACI'   ,    _RC)
   call MAPL_GetPointer(EXPORT,FRW    , 'FRACW'   ,    _RC)
   call MAPL_GetPointer(EXPORT,PENUVR , 'PENUVR'  ,    _RC)
   call MAPL_GetPointer(EXPORT,PENUVF , 'PENUVF'  ,    _RC)
   call MAPL_GetPointer(EXPORT,PENPAR , 'PENPAR'  ,    _RC)
   call MAPL_GetPointer(EXPORT,PENPAF , 'PENPAF'  ,    _RC)
   call MAPL_GetPointer(EXPORT,FSURF  , 'FSURF'   ,    _RC)
   call MAPL_GetPointer(EXPORT,SWFLX  , 'SWFLX'   ,    _RC)

     call MAPL_GetPointer(EXPORT, DISCHARGE, 'DISCHARGE', _RC)
     if(associated(DISCHARGE)) DISCHARGE = DISCHARGE_IM

   call MAPL_Get (MAPL, GCS=GCS, GIM=GIM, GEX=GEX, GCnames=GCnames,_RC)

   do I = 1, size(GCS)
      if (trim(GCnames(I)) .ne. 'OBIO') then
         if(associated(TST)) then
            call MAPL_GetPointer(GEX(I), dummy, 'TST'   , alloc=.true., _RC)
         endif 
         if(associated(QST)) then
            call MAPL_GetPointer(GEX(I), dummy, 'QST'   , alloc=.true., _RC)
         endif 
         if(associated(DELTS)) then
            call MAPL_GetPointer(GEX(I), dummy, 'DELTS' , alloc=.true., _RC)
         endif 
         if(associated(DELQS)) then
            call MAPL_GetPointer(GEX(I), dummy, 'DELQS' , alloc=.true., _RC)
         endif 
         if(associated(EVAPOUT)) then
            call MAPL_GetPointer(GEX(I), dummy, 'EVAPOUT', alloc=.true., _RC)
         endif 
         if(associated(SHOUT)) then
            call MAPL_GetPointer(GEX(I), dummy, 'SHOUT' , alloc=.true., _RC)
         endif 
         if(associated(HLATN)) then
            call MAPL_GetPointer(GEX(I), dummy, 'HLATN' , alloc=.true., _RC)
         endif 
         if(associated(HLWUP)) then
            call MAPL_GetPointer(GEX(I), dummy, 'HLWUP' , alloc=.true., _RC)
         endif 
         if(associated(LWNDSRF)) then
            call MAPL_GetPointer(GEX(I), dummy, 'LWNDSRF', alloc=.true., _RC)
         endif  
         if(associated(SWNDSRF)) then
            call MAPL_GetPointer(GEX(I), dummy, 'SWNDSRF', alloc=.true., _RC)
         endif  
         if(associated(PENUVR)) then
            call MAPL_GetPointer(GEX(I), dummy, 'PENUVR', alloc=.true., _RC)
         endif  
         if(associated(PENUVF)) then
            call MAPL_GetPointer(GEX(I), dummy, 'PENUVF', alloc=.true., _RC)
         endif  
         if(associated(PENPAR)) then
            call MAPL_GetPointer(GEX(I), dummy, 'PENPAR', alloc=.true., _RC)
         endif  
         if(associated(PENPAF)) then
            call MAPL_GetPointer(GEX(I), dummy, 'PENPAF', alloc=.true., _RC)
         endif  
         if(associated(FSURF)) then
            call MAPL_GetPointer(GEX(I), dummy, 'FSURF' , alloc=.true., _RC)
         endif
      endif
   enddo

! Call the childrens' RUN2
!-------------------------

    DO I = 1, size(GCS)
       call MAPL_TimerOn(MAPL,trim(GCnames(i)), _RC )
       call ESMF_GridCompRun(GCS(I), importState=GIM(I), exportState=GEX(I), &
                             CLOCK=CLOCK, PHASE=2, userRC=STATUS)
       VERIFY_(STATUS)
       call MAPL_TimerOff(MAPL,trim(GCnames(i)), _RC )
    ENDDO

    call MAPL_GetPointer(GEX(ICE), FI, 'FRACI'    , _RC)

    ! make sure following fractions are bounded in [0.0, 1.0]
    FR(:,WATER) = max(1.0-FI, 0.0)
    FR(:,ICE  ) = min(FI,     1.0)

    ! FRACINEW is the updated ice fraction by sea ice model in RUN2
    ! it is available here to aggregate some fields
    ! SA: Bin, we should probably add a diagnostic export here to keep track of FRNEW-FR
    call MAPL_GetPointer(GEX(ICE), FI, 'FRACINEW' , _RC)

    FRNEW(:,WATER) = max(1.0-FI, 0.0)
    FRNEW(:,ICE  ) = min(FI,     1.0)

    if(associated(FRI)) FRI = FRNEW(:,  ICE)
    if(associated(FRW)) FRW = FRNEW(:,WATER)


                            EMISS   = 0.0
                            ALBVR   = 0.0
                            ALBVF   = 0.0
                            ALBNR   = 0.0
                            ALBNF   = 0.0
    if(associated(EVAPOUT)) EVAPOUT = 0.0
    if(associated(SHOUT  )) SHOUT   = 0.0
    if(associated(HLATN  )) HLATN   = 0.0
    if(associated(DELTS  )) DELTS   = 0.0
    if(associated(DELQS  )) DELQS   = 0.0
    if(associated(TST    )) TST     = 0.0
    if(associated(QST    )) QST     = 0.0
    if(associated(HLWUP  )) HLWUP   = 0.0
    if(associated(LWNDSRF)) LWNDSRF = 0.0
    if(associated(SWNDSRF)) SWNDSRF = 0.0
    if(associated(PENUVR )) PENUVR  = 0.0
    if(associated(PENUVF )) PENUVF  = 0.0
    if(associated(PENPAR )) PENPAR  = 0.0
    if(associated(PENPAF )) PENPAF  = 0.0
    if(associated(FSURF  )) FSURF   = 0.0
    if(associated(SWFLX  )) SWFLX   = 0.0

! Cycle through sub-tiles aggregating fluxes
!-------------------------------------------
    do N=1,NUM_SUBTILES

       call MAPL_GetPointer(GEX(N), EVP  , 'EVAPOUT'  , _RC)
       call MAPL_GetPointer(GEX(N), SHF  , 'SHOUT'    , _RC)
       call MAPL_GetPointer(GEX(N), LHF  , 'HLATN'    , _RC)
       call MAPL_GetPointer(GEX(N), DTS  , 'DELTS'    , _RC)
       call MAPL_GetPointer(GEX(N), DQS  , 'DELQS'    , _RC)
       call MAPL_GetPointer(GEX(N), TS   , 'TST'      , _RC)
       call MAPL_GetPointer(GEX(N), QS   , 'QST'      , _RC)
       call MAPL_GetPointer(GEX(N), EMS  , 'EMIS'     , _RC)
       call MAPL_GetPointer(GEX(N), AVR  , 'ALBVR'    , _RC)
       call MAPL_GetPointer(GEX(N), AVF  , 'ALBVF'    , _RC)
       call MAPL_GetPointer(GEX(N), ANR  , 'ALBNR'    , _RC)
       call MAPL_GetPointer(GEX(N), ANF  , 'ALBNF'    , _RC)
       call MAPL_GetPointer(GEX(N), HLW  , 'HLWUP'    , _RC)
       call MAPL_GetPointer(GEX(N), LWND , 'LWNDSRF'  , _RC)
       call MAPL_GetPointer(GEX(N), SWND , 'SWNDSRF'  , _RC)
       call MAPL_GetPointer(GEX(N), PUR  , 'PENUVR'   , _RC)
       call MAPL_GetPointer(GEX(N), PUF  , 'PENUVF'   , _RC)
       call MAPL_GetPointer(GEX(N), PAR  , 'PENPAR'   , _RC)
       call MAPL_GetPointer(GEX(N), PAF  , 'PENPAF'   , _RC)
       call MAPL_GetPointer(GEX(N), FSUR , 'FSURF'    , _RC)

                               EMISS   = EMISS   + EMS    *FR(:,N)
                               ! special treatment for albedo 
                               ! use updated ice fraction to aggregate albedos 
                               ! only relevant for coupled; AMIP sees the same fraction
                               ! albedo computed here only correct when solar alarm 
                               ! is ringing
                               ALBVR   = ALBVR   + AVR    *FRNEW(:,N)
                               ALBVF   = ALBVF   + AVF    *FRNEW(:,N)
                               ALBNR   = ALBNR   + ANR    *FRNEW(:,N)
                               ALBNF   = ALBNF   + ANF    *FRNEW(:,N)
       ! fields below need also to be revaluated as to which fraction
       ! should be used for aggregation
       if(associated(EVAPOUT)) EVAPOUT = EVAPOUT + EVP    *FR(:,N)
       if(associated(SHOUT  )) SHOUT   = SHOUT   + SHF    *FR(:,N)
       if(associated(HLATN  )) HLATN   = HLATN   + LHF    *FR(:,N)
       if(associated(DELTS  )) DELTS   = DELTS   + DTS    *FR(:,N)
       if(associated(DELQS  )) DELQS   = DELQS   + DQS    *FR(:,N)
       if(associated(TST    )) TST     = TST     + TS     *FR(:,N)
       if(associated(QST    )) QST     = QST     + QS     *FR(:,N)
       if(associated(HLWUP  )) HLWUP   = HLWUP   + HLW    *FR(:,N)
       if(associated(LWNDSRF)) LWNDSRF = LWNDSRF + LWND   *FR(:,N)
       if(associated(SWNDSRF)) SWNDSRF = SWNDSRF + SWND   *FR(:,N)
       if(associated(PENUVR )) PENUVR  = PENUVR  + PUR    *FR(:,N)
       if(associated(PENUVF )) PENUVF  = PENUVF  + PUF    *FR(:,N)
       if(associated(PENPAR )) PENPAR  = PENPAR  + PAR    *FR(:,N)
       if(associated(PENPAF )) PENPAF  = PENPAF  + PAF    *FR(:,N)
       if(associated(FSURF  )) FSURF   = FSURF   + FSUR   *FR(:,N)
       if(associated(SWFLX  )) SWFLX   = SWFLX   + (PUR+PUF+PAR+PAF)*FR(:,N)
    enddo

    if(associated(SWFLX  )) then
       call MAPL_GetPointer(GEX(WATER), dummy, 'AO_DRNIR' ,  _RC)
       SWFLX  = SWFLX + dummy
       call MAPL_GetPointer(GEX(WATER), dummy, 'AO_DFNIR' ,  _RC)
       SWFLX  = SWFLX + dummy
    endif

! Average ocean stress
!---------------------
    call MAPL_GetPointer(GEX(WATER), TXW  , 'TAUXW'     , _RC)
    call MAPL_GetPointer(GEX(WATER), TYW  , 'TAUYW'     , _RC)
    call MAPL_GetPointer(GEX(ICE)  , TXI  , 'TAUXI'     , _RC)
    call MAPL_GetPointer(GEX(ICE)  , TYI  , 'TAUYI'     , _RC)

    TXO = TXI*FR(:,ICE) + TXW*FR(:,WATER)
    TYO = TYI*FR(:,ICE) + TYW*FR(:,WATER)
    if(associated(TAUXO)) TAUXO = TXO
    if(associated(TAUYO)) TAUYO = TYO

    if(associated(PSEX )) PSEX  = PS 
    if(associated(USTR3)) USTR3 = sqrt(sqrt(TXO*TXO+TYO*TYO)/water_RHO('salt_water'))**3
    if(associated(UUEX))  UUEX  = UU

!  All done with SALTWATERCORE
!-----------------------------

    RETURN_(ESMF_SUCCESS)
             
  end subroutine SALTWATERCORE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine RUN2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_SaltwaterGridCompMod

