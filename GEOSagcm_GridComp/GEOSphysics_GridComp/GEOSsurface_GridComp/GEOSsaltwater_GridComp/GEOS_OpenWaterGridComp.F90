
!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
module GEOS_OpenwaterGridCompMod

!BOP
! !MODULE: GEOS_OpenwaterGridCompMod -- Implements slab saltwater/Water tiles.

! !DESCRIPTION:
! 
!   {\tt GEOS\_Openwater} is a light-weight gridded component that updates
!      the skin sub-tiles at saltwater/Water points, be they ocean, estuary, or salt
!      lake. Currently each tile can have only two subtiles, open-water and ice.
!
!      The component is written with a two stage run method for use with
!      semi-implicit turbulence components. The first run stage computes
!      exchange coefficients for heat, moisture and momentum at each sub-tile
!      and combines these to tile space, accounting for sub tile variability
!      by passing back an effective surface value of the exchanged quantity.
!      ----------------------------------------------------------------------------
!
!      --------- rewrite ---------
!      \noindent The "OPENWATERCORE" implements following AOIL.
!      \noindent If the AOIL is OFF, near-surface temperature variations are neglected, 
!      $T_s = T_o$ from the ocean model, or $=T_f$ from boundary conditions.\\~\\
!
!      \noindent Following are the variables and fluxes used to update $T_s$
!
!      \begin{itemize}
!       \item \textbf{Variables}:
!       \begin{enumerate}
!              \item \textbf{Coupled with Ocean model}:\\
!               $T_o$ and $d, D, \epsilon_d = d/D,$ 
!               where $d$ and $D$ denote the thickness of AOIL and ocean model top level respectively
!               \item \textbf{Uncoupled, i.e., boundary conditions}:\\
!               $T_f$ and $d$
!       \end{enumerate}
!       \item \textbf{Fluxes}:
!       \begin{enumerate}
!              \item \textbf{Coupled with Ocean model}:\\
!               $$Q_{\sigma} = Q_w - \left( \frac{\epsilon_d}{1-\epsilon_d} \right) Q_f$$
!               where $$Q_w = SW_{top} - SW_d + Q^{\downarrow}$$
!               $$Q_f = SW_d - SW_D$$
!               $SW$ and $Q^{\downarrow}$ denote solar and non-solar fluxes respectively\\
!               \textbf{Ocean Model must always receive}  $Q_w+Q_f =  SW_{top} - SW_D + Q^{\downarrow} ;$ it should be 
!               agnostic to the presence/absence of the AOIL. It {\it senses} AOIL via changes to the fluxes brought 
!               about by changes to $T_s$ due to the action of AOIL; recall $T_s = T_o$ if AOIL is OFF.  
!               \item \textbf{Uncoupled, i.e., boundary conditions}:\\
!               $$Q_{\sigma} = Q_w = SW_{top} - SW_d + Q^{\downarrow}$$
!       \end{enumerate}
!      \end{itemize}
!      \noindent The AOIL calculates $T_s = T_{\delta} - \Delta T_c,$ see below for details of $T_{\delta}$.
!      \begin{itemize}
!               \item \textbf{Coupled}:\\
!               $$T_{\delta}= T_o + \left( \frac{1}{\mu} + (1-\epsilon_d) \right) \sigma_T$$
!               \item \textbf{Uncoupled}:\\
!               $$T_{\delta}= T_f + \left( \frac{1+\mu}{\mu}\right) \sigma_T$$
!      \end{itemize} 
!
!      \noindent For details, please see Akella and Suarez, 2018, 
!      "The Atmosphere-Ocean Interface Layer of the NASA Goddard Earth Observing System Model and Data Assimilation System." GMAO Tech Memo, Vol 51. 
!      ----------------------------------------------------------------------------
!

! !USES:

  use sfclayer  ! using module that contains sfc layer code
  use ESMF
  use MAPL
  use GEOS_UtilsMod
  use DragCoefficientsMod
  use atmOcnIntlayer,     only: ALBSEA,            &
                                AOIL_sfcLayer_T,   &
                                water_RHO,         &
                                AOIL_Shortwave_abs,&
                                AOIL_v0

  implicit none
  private

  integer            :: DO_DATASEA               ! =1:uncoupled (AGCM); =0:coupled (AOGCM)
  integer            :: DO_SKIN_LAYER            ! =1:active AOIL (ON); =0:inactive AOIL (OFF)

  public SetServices

!EOP

  integer, parameter            :: WATER        = 1      
  integer, parameter            :: NUM_SUBTILES = 1             ! number of sub-tiles
  real,    parameter            :: KUVR         = 0.09

  type openwater_state
       integer:: CHOOSEMOSFC
  end type openwater_state

  type openwater_state_wrap
      type(openwater_state), pointer :: ptr
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

    integer                                 :: DO_WAVES

    type(openwater_state_wrap) :: wrap
    type(openwater_state), pointer :: mystate
    character(len=ESMF_MAXSTR)     :: SURFRC
    type(ESMF_Config)              :: SCF

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, _RC)
    Iam = trim(COMP_NAME) // 'SetServices'

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run1, _RC)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2, _RC)

! Get constants from CF
! ---------------------

    call MAPL_GetResource ( MAPL, DO_DATASEA,    Label="USE_DATASEA:"     , DEFAULT=1, _RC)

    call MAPL_GetResource ( MAPL, DO_SKIN_LAYER, Label="USE_SKIN_LAYER:"  , DEFAULT=0, _RC)

    call MAPL_GetResource ( MAPL, DO_WAVES,      Label="USE_WAVES:"       , DEFAULT=0, RC=STATUS)
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
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_reflectivity_for_visible_beam',   &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVR',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_reflectivity_for_visible_diffuse',&
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVF',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_reflectivity_for_near_infrared_beam', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNR',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_reflectivity_for_near_infrared_diffuse', &
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
        LONG_NAME          = 'ocean_snowfall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'SNOWOCN'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'ocean_rainfall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'RAINOCN'                   ,&
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
        LONG_NAME          = 'open_water_upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHWTR'                     ,&
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
        LONG_NAME          = 'open_water_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWNDWTR'                   ,&
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
        LONG_NAME          = 'open_water_net_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWNDWTR'                   ,&
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
                                               _RC)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'open_water_latent_energy_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLATWTR'                   ,&
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
        SHORT_NAME         = 'FRACW',                             &
        LONG_NAME          = 'water_covered_fraction_of_tile',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

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
        LONG_NAME          = 'eastward_stress_over_water',&
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUXW'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'northward_stress_over_water',&
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUYW'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
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

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'DCOOL',                     &
          LONG_NAME          = 'depth_of_cool_layer'       ,&
          UNITS              = 'm'                         ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC                          ,&
          SHORT_NAME         = 'DWARM',                      &
          LONG_NAME          = 'depth_at_base_of_warm_layer',&
          UNITS              = 'm'                          ,&
          DIMS               = MAPL_DimsTileOnly            ,&
          VLOCATION          = MAPL_VLocationNone           ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC                                 ,&
          SHORT_NAME         = 'TDROP',                             &
          LONG_NAME          = 'temperature_drop_across_cool_layer',&
          UNITS              = 'K'                                 ,&
          DIMS               = MAPL_DimsTileOnly                   ,&
          VLOCATION          = MAPL_VLocationNone                  ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'QCOOL',                     &
          LONG_NAME          = 'net_heating_in_cool_layer' ,&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'QWARM',                     &
          LONG_NAME          = 'net_heating_in_warm_layer' ,&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC                          ,&
          SHORT_NAME         = 'SWCOOL',                     &
          LONG_NAME          = 'solar_heating_in_cool_layer',&
          UNITS              = 'W m-2'                      ,&
          DIMS               = MAPL_DimsTileOnly            ,&
          VLOCATION          = MAPL_VLocationNone           ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC                          ,&
          SHORT_NAME         = 'SWWARM',                     &
          LONG_NAME          = 'solar_heating_in_warm_layer',&
          UNITS              = 'W m-2'                      ,&
          DIMS               = MAPL_DimsTileOnly            ,&
          VLOCATION          = MAPL_VLocationNone           ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC                           ,&
          SHORT_NAME         = 'PHIW',                        &
          LONG_NAME          = 'Similarity_function_in_warm_layer',&
          UNITS              = '1'                           ,&
          DIMS               = MAPL_DimsTileOnly             ,&
          VLOCATION          = MAPL_VLocationNone            ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'LANGM'                     ,&
          LONG_NAME          = 'Langmuir_number'           ,&
          UNITS              = '1'                         ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC                        ,&
          SHORT_NAME         = 'USTARW',                    &
          LONG_NAME          = 'ustar_over_water'          ,&
          UNITS              = 'm s-1'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC                                  ,&
          SHORT_NAME         = 'TBAR',                               &
          LONG_NAME          = 'mean_temperature_of_interface_layer',&
          UNITS              = 'K'                                  ,&
          DIMS               = MAPL_DimsTileOnly                    ,&
          VLOCATION          = MAPL_VLocationNone                   ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'LCOOL',                     &
          LONG_NAME          = 'Saunders_parameter'        ,&
          UNITS              = '1'                         ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'BCOOL',                     &
          LONG_NAME          = 'bouyancy_generation_in_cool_layer',&
          UNITS              = 'm+2 s-3'                   ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC  ) 

     call MAPL_AddExportSpec(GC,                        &
          SHORT_NAME         = 'TDEL',                     &
          LONG_NAME          = 'temperature_at_base_of_cool_layer', &
          UNITS              = 'K',                        &
          DIMS               = MAPL_DimsTileOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
          _RC)

     call MAPL_AddExportSpec(GC,                        &
          SHORT_NAME         = 'TS_FOUND',                 &
          LONG_NAME          = 'foundation_temperature_for_interface_layer',        &
          UNITS              = 'K',                        &
          DIMS               = MAPL_DimsTileOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
          _RC)

     call MAPL_AddExportSpec(GC,                        &
          SHORT_NAME         = 'TAUTW',                    &
          LONG_NAME          = 'relaxation_time_of_TW_to_TS_FOUND', &
          UNITS              = 's',                        &
          DIMS               = MAPL_DimsTileOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
          _RC)

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'ZETA_W'                    ,&
          LONG_NAME          = 'Stability_parameter_in_Warm_Layer',                 &
          UNITS              = '1'                         ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          _RC)

     call MAPL_AddExportSpec(GC,                           &
          SHORT_NAME         = 'SS_FOUND',                 &
          LONG_NAME          = 'foundation_salinity_for_interface_layer',        &
          UNITS              = 'PSU',                      &
          DIMS               = MAPL_DimsTileOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
          _RC)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'total_surface_heat_flux_over_the_whole_tile' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'FSURF'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     ! Atmosphere-ocean fluxes
     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'atmosphere_ocean_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'AO_SHFLX'                  ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'atmosphere_ocean_evaporation' ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'AO_QFLUX'                  ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'atmosphere_ocean_net_longwave_radiation' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'AO_LWFLX'                  ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'atmosphere_ocean_snowfall' ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'AO_SNOW'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'atmosphere_ocean_rainfall' ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'AO_RAIN'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               _RC  ) 

    call MAPL_AddExportSpec(GC                         ,&
         LONG_NAME          = 'net_surface_downwelling_nir_beam_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'AO_DRNIR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  _RC  ) 

    call MAPL_AddExportSpec(GC                         ,&
         LONG_NAME          = 'net_surface_downwelling_nir_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'AO_DFNIR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  _RC  ) 

    call MAPL_AddExportSpec(GC,                                    &
         LONG_NAME          = 'departure_of_mean_interface_temperature_from_foundation_temperature',   &
         UNITS              = 'K',                                 &
         SHORT_NAME         = 'TWMTF',                             &
         DIMS               = MAPL_DimsTileOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

    call MAPL_AddExportSpec(GC,                                    &
         LONG_NAME          = 'penetrated_shortwave_flux_at_the_bottom_of_first_ocean_model_layer',   &
         UNITS              = 'W m-2',                             &
         SHORT_NAME         = 'PEN_OCN',                           &
         DIMS               = MAPL_DimsTileOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC  )

!  !INTERNAL STATE:

    call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'HSKINW',                            &
        LONG_NAME          = 'water_skin_layer_mass',             &
        UNITS              = 'kg m-2',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = trim(COMP_NAME),                     & ! friendly are a pain to get rid of!
        DEFAULT            = 5.0*MAPL_RHO_SEAWATER,               &
                                                       _RC  )

    call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'TSKINW',                            &
        LONG_NAME          = 'water_skin_temperature',            &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = trim(COMP_NAME),                     & ! friendly are a pain to get rid of!
        DEFAULT            = 280.0,                               &
                                                       _RC  )

    call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'SSKINW',                            &
        LONG_NAME          = 'water_skin_salinity',               &
        UNITS              = 'psu',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        FRIENDLYTO         = trim(COMP_NAME),                     & ! friendly are a pain to get rid of!
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
                                                       _RC)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CH',                                &
        LONG_NAME          = 'surface_heat_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 1.0e-4,                              &
                                                       _RC)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CM',                                &
        LONG_NAME          = 'surface_momentum_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 1.0e-4,                              &
                                                       _RC)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CQ',                                &
        LONG_NAME          = 'surface_moisture_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 1.0e-4,                              &
                                                       _RC)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'Z0',                                &
        LONG_NAME          = 'aerodynamic_roughness',             &
        UNITS              = 'm',                                 &
        DEFAULT            = 0.00005,                             &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'WW',                                &
        LONG_NAME          = 'vertical_velocity_scale_squared',   &
        UNITS              = 'm+2 s-2',                           &
        DEFAULT            = 0.0,                                 &
        NUM_SUBTILES       = NUM_SUBTILES,                        &
        DIMS               = MAPL_DimsTileTile,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddInternalSpec(GC,                                &
        SHORT_NAME         = 'TWMTS',                             &
        LONG_NAME          = 'departure_of_skin_temperature_from_mean_interface_temperature',   &
        UNITS              = 'K',                                 &
        DEFAULT            = 0.0,                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddInternalSpec(GC,                                &
        SHORT_NAME         = 'TWMTF',                             &
        LONG_NAME          = 'departure_of_mean_interface_temperature_from_foundation_temperature',   &
        UNITS              = 'K',                                 &
        DEFAULT            = 0.0,                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'DELTC',                             &
        LONG_NAME          = 'temperature_drop_across_cool_layer',&
        UNITS              = 'K',                                 &
        DEFAULT            = 0.0,                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

!  !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'ALW',                               &
        LONG_NAME          = 'linearization_of_surface_emitted_longwave_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'BLW',                               &
        LONG_NAME          = 'linearization_of_surface_emitted_longwave_flux', &
        UNITS              = 'W m-2 K-1',                         &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'LWDNSRF',                           &
        LONG_NAME          = 'surface_absorbed_longwave_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

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
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'upward_sensible_heat_flux',         &
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'SH',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'eastward_surface_stress',           &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUX',                              &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'northward_surface_stress',          &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUY',                              &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_evaporation',         &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'DEVAP',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_upward_sensible_heat_flux', &
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'DSH',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'snowfall',                          &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'SNO',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

! Surface air quantities

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_air_temperature',           &
        UNITS              = 'K',                                 &
        SHORT_NAME         = 'TA',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_air_specific_humidity',     &
        UNITS              = 'kg kg-1',                           &
        SHORT_NAME         = 'QA',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_wind_speed',                &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'UU',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'levellm_uwind',                     &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'UWINDLMTILE',                       &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'levellm_vwind',                     &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'VWINDLMTILE',                       &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_layer_height',              &
        UNITS              = 'm',                                 &
        SHORT_NAME         = 'DZ',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_pressure',                  &
        UNITS              = 'Pa',                                &
        SHORT_NAME         = 'PS',                                &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

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
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'QHATM',                             &
        LONG_NAME          = 'effective_surface_specific_humidity',&
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'UHATM',                             &
        LONG_NAME          = 'effective_surface_zonal_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VHATM',                             &
        LONG_NAME          = 'effective_surface_meridional_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CTATM',                             &
        LONG_NAME          = 'surface_exchange_coefficient_for_heat', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CQATM',                             &
        LONG_NAME          = 'surface_exchange_coefficient_for_moisture', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CMATM',                             &
        LONG_NAME          = 'surface_exchange_coefficient_for_momentum', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'FRACI',                             &
        LONG_NAME          = 'ice_covered_fraction_of_tile',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'FRACINEW',                             &
        LONG_NAME          = 'ice_covered_fraction_of_tile_after_update',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'UW',                                &
        LONG_NAME          = 'zonal_velocity_of_surface_water',   &
        UNITS              = 'm s-1 ',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VW',                                &
        LONG_NAME          = 'meridional_velocity_of_surface_water',   &
        UNITS              = 'm s-1 ',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'KPAR',                              &
        LONG_NAME          = 'PAR_extinction_coefficient',        &
        UNITS              = 'm-1',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       _RC)

     call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'TS_FOUND',                        &
          LONG_NAME          = 'foundation_temperature_for_interface_layer',&
          UNITS              = 'K',                               &
          DIMS               = MAPL_DimsTileOnly,                 &
          VLOCATION          = MAPL_VLocationNone,                &
          DEFAULT            = 280.0,                             &
          _RC)

     call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'SS_FOUND',                        &
          LONG_NAME          = 'foundation_salinity_for_interface_layer',&
          UNITS              = 'PSU',                             &
          DIMS               = MAPL_DimsTileOnly,                 &
          VLOCATION          = MAPL_VLocationNone,                &
          DEFAULT            = 30.0,                              &
          _RC)

     call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'TFREEZE',                        &
          LONG_NAME          = 'freezing_temperature_for_interface_layer',&
          UNITS              = 'K',                               &
          DIMS               = MAPL_DimsTileOnly,                 &
          VLOCATION          = MAPL_VLocationNone,                &
          DEFAULT            = MAPL_TICE-1.8,                     &
          _RC)

    call MAPL_AddImportSpec (GC,                                   &
         SHORT_NAME = 'DTSDT',                                     &
         LONG_NAME  = 'skin_temperature_analysis_tendency',        &
         UNITS      = 'K s-1',                                     &
         RESTART    = MAPL_RestartSkip,                            &
         DIMS       = MAPL_DimsTileOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               _RC  )
    
    call MAPL_AddImportSpec(GC,                    &
          LONG_NAME          = 'river_discharge_at_ocean_points',&
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'DISCHARGE'                 ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartSkip            ,&
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
        _RC)

     if (DO_WAVES /= 0) then
       call MAPL_AddImportSpec(GC,                                &
          SHORT_NAME       = 'CHARNOCK',                          &
          LONG_NAME        = 'charnock_coefficient',              &
          UNITS            = '1',                                 &
          RESTART          = MAPL_RestartSkip,                    &
          DIMS             = MAPL_DimsTileOnly,                   &
          VLOCATION        = MAPL_VLocationNone,                  &
          RC=STATUS  )
       VERIFY_(STATUS)
     end if

!EOS

    allocate(mystate,_STAT)
    call MAPL_GetResource (MAPL, SURFRC, label = 'SURFRC:', default = 'GEOS_SurfaceGridComp.rc', _RC)
    SCF = ESMF_ConfigCreate(_RC)
    call ESMF_ConfigLoadFile     (SCF,SURFRC,_RC)
    call MAPL_GetResource (SCF, mystate%CHOOSEMOSFC, label='CHOOSEMOSFC:', DEFAULT=1, _RC )
    call ESMF_ConfigDestroy      (SCF, _RC)
    wrap%ptr => mystate
    call ESMF_UserCompSetInternalState(gc, 'openwater_private', wrap,_RC)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="RUN1"  ,                _RC)
    call MAPL_TimerAdd(GC,    name="RUN2"  ,                _RC)
    call MAPL_TimerAdd(GC,    name="-OpenWater",            _RC) 
    call MAPL_TimerAdd(GC,    name="-Albedo"  ,             _RC)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC,  _RC)
 
! Set the Run entry point
! -----------------------

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP
! !IROUTINE: RUN1 -- First Run stage for the Openwater component
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

! pointers to internal

   real, pointer, dimension(:  )  :: TW  => null()
   real, pointer, dimension(:  )  :: HW  => null()
   real, pointer, dimension(:  )  :: TWMTS => null()
   real, pointer, dimension(:  )  :: TWMTF => null()
   real, pointer, dimension(:  )  :: DELTC => null()

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
   real, pointer, dimension(:)    :: UW  => null()
   real, pointer, dimension(:)    :: VW  => null()
   real, pointer, dimension(:)    :: DZ  => null()     
   real, pointer, dimension(:)    :: TA  => null()
   real, pointer, dimension(:)    :: QA  => null()     
   real, pointer, dimension(:)    :: PS  => null()
   real, pointer, dimension(:)    :: PCU => null()
   real, pointer, dimension(:)    :: FI  => null()
   real, pointer, dimension(:)    :: TF  => null()
   real, pointer, dimension(:)    :: TS_FOUNDi => null()
   real, pointer, dimension(:)    :: CHARNOCK  => null()


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

   real            :: OCEANZ0
   real            :: DT
   real, parameter :: HPBL = 1000.

   integer         :: CHOOSEMOSFC
   integer         :: CHOOSEZ0
   real            :: MUSKIN                    ! exponent in T(z) profile in warm layer, based on Zeng & Beljaars, 2005
   real            :: AOIL_depth                ! thickness of atmosphere-ocean interface layer (AOIL) denoted by d in AS2018
   real            :: epsilon_d                 ! ratio: (thickness of AOIL)/(thickness of OGCM top level) = epsilon_d in Akella & Suarez, 2018
   real            :: MaxWaterDepth
   real            :: MinWaterDepth
   real            :: OGCM_top_thickness        ! thickness of OGCM top layer (D) in AS2018

   type(openwater_state_wrap) :: wrap
   type(openwater_state), pointer :: mystate
   character(len=100)             :: WHICH_T_TO_SFCLAYER    ! what temperature does the sfclayer get from AOIL?
   real                           :: DEPTH_T_TO_SFCLAYER    ! temperature (at what depth) does the sfclayer get from AOIL?

   integer                        :: DO_WAVES
   real, allocatable              :: CHARNOCK_(:)

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

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN1" )

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,                          &
         INTERNAL_ESMF_STATE = INTERNAL,         &
                                       _RC)

    call MAPL_Get(MAPL, HEARTBEAT = DT, _RC)
    call MAPL_GetResource ( MAPL, DT, Label="DT:", DEFAULT=DT, _RC)

! Get parameters (0:Louis, 1:Monin-Obukhov)
! -----------------------------------------
! -----------------------------------------
    call ESMF_UserCompGetInternalState(gc,'openwater_private',wrap,_RC)
    mystate => wrap%ptr
    CHOOSEMOSFC = mystate%CHOOSEMOSFC

    call MAPL_GetResource ( MAPL, CHOOSEZ0,    Label="CHOOSEZ0:",    DEFAULT=3, _RC)

! Get roughness parameters 
! -------------------------------------------------------------
    call MAPL_GetResource ( MAPL, OCEANZ0,     Label="OCEANZ0:" ,    DEFAULT=1.0e-3, _RC) 

! Get Thickness of AOIL (m)
! -------------------------
    if (DO_SKIN_LAYER==0) then
       call MAPL_GetResource ( MAPL, MaxWaterDepth, Label="MAX_WATER_DEPTH:" , DEFAULT=1000., _RC)
       call MAPL_GetResource ( MAPL, MinWaterDepth, Label="MIN_WATER_DEPTH:" , DEFAULT=1000., _RC)
    else 
       call MAPL_GetResource ( MAPL, MaxWaterDepth, Label="MAX_WATER_DEPTH:" , DEFAULT=2.,   _RC)
       call MAPL_GetResource ( MAPL, MinWaterDepth, Label="MIN_WATER_DEPTH:" , DEFAULT=2.,   _RC)

!      Exponent in the near-surface temperature profile T(z) within the AOIL
!      ---------------------------------------------------------------------
       call MAPL_GetResource ( MAPL, MUSKIN,        Label="MU_SKIN:"         , DEFAULT=0.2  ,   _RC)
    end if

    AOIL_depth = MAX(MaxWaterDepth, MinWaterDepth)

    if (DO_DATASEA==0) then
!      Thickness of OGCM top level (m)
!      ------------------------------
       call MAPL_GetResource ( MAPL, OGCM_top_thickness, Label="OGCM_TOP_LAYER:" , DEFAULT=10.,   _RC) ! SA: could be an export from GUEST GC
       epsilon_d  = AOIL_depth/OGCM_top_thickness ! < 1. If that is NOT true, AOIL formulation would need revisit; see AS2018
    else
       OGCM_top_thickness = MAPL_UNDEF
       epsilon_d          = 0.0
    end if

! Is the wave model enabled
! -------------------------
    call MAPL_GetResource ( MAPL, DO_WAVES,         Label="USE_WAVES:",           DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)
 
! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,UU     , 'UU'     ,    _RC)
   call MAPL_GetPointer(IMPORT,UWINDLMTILE     , 'UWINDLMTILE'     ,    _RC)
   call MAPL_GetPointer(IMPORT,VWINDLMTILE     , 'VWINDLMTILE'     ,    _RC)
   call MAPL_GetPointer(IMPORT,UW     , 'UW'     ,    _RC)
   call MAPL_GetPointer(IMPORT,VW     , 'VW'     ,    _RC)
   call MAPL_GetPointer(IMPORT,DZ     , 'DZ'     ,    _RC)
   call MAPL_GetPointer(IMPORT,TA     , 'TA'     ,    _RC)
   call MAPL_GetPointer(IMPORT,QA     , 'QA'     ,    _RC)
   call MAPL_GetPointer(IMPORT,PS     , 'PS'     ,    _RC)
   call MAPL_GetPointer(IMPORT,PCU    , 'PCU'    ,    _RC)
   call MAPL_GetPointer(IMPORT,FI     , 'FRACI'  ,    _RC)
   call MAPL_GetPointer(IMPORT,TF     , 'TFREEZE',    _RC)
   call MAPL_GetPointer(IMPORT,TS_FOUNDi, 'TS_FOUND', _RC)

   NT = size(TA)
   allocate(CHARNOCK_(NT),    STAT=STATUS)
   VERIFY_(STATUS)
   
   if (DO_WAVES /= 0) then
     call MAPL_GetPointer(IMPORT,CHARNOCK, 'CHARNOCK', RC=STATUS)
     VERIFY_(STATUS)

     where (CHARNOCK > 0 .and. CHARNOCK < 1.0)
       CHARNOCK_ = CHARNOCK
     elsewhere
       CHARNOCK_ = 0.0185
     end where
   else
     CHARNOCK_ = 0.0185
   end if

! Pointers to internals
!----------------------

   call MAPL_GetPointer(INTERNAL,QS   , 'QS'     ,    _RC)
   call MAPL_GetPointer(INTERNAL,CH   , 'CH'     ,    _RC)
   call MAPL_GetPointer(INTERNAL,CM   , 'CM'     ,    _RC)
   call MAPL_GetPointer(INTERNAL,CQ   , 'CQ'     ,    _RC)
   call MAPL_GetPointer(INTERNAL,Z0   , 'Z0'     ,    _RC)
   call MAPL_GetPointer(INTERNAL,WW   , 'WW'     ,    _RC)
   call MAPL_GetPointer(INTERNAL,TW   , 'TSKINW' ,    _RC)
   call MAPL_GetPointer(INTERNAL,HW   , 'HSKINW' ,    _RC)
   call MAPL_GetPointer(INTERNAL,TWMTF, 'TWMTF'  ,    _RC)
   call MAPL_GetPointer(INTERNAL,DELTC, 'DELTC'  ,    _RC)

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
   call MAPL_GetPointer(EXPORT,MOT2M, 'MOT2M'   ,    _RC)
   call MAPL_GetPointer(EXPORT,MOQ2M, 'MOQ2M'   ,    _RC)
   call MAPL_GetPointer(EXPORT,MOU2M, 'MOU2M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOV2M, 'MOV2M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOT10M, 'MOT10M'   ,    _RC)
   call MAPL_GetPointer(EXPORT,MOQ10M, 'MOQ10M'   ,    _RC)
   call MAPL_GetPointer(EXPORT,MOU10M, 'MOU10M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOV10M, 'MOV10M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOU50M, 'MOU50M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOV50M, 'MOV50M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,GST   , 'GUST'    ,    _RC)
   call MAPL_GetPointer(EXPORT,VNT   , 'VENT'    ,    _RC)

   NT = size(TA)
   if(NT == 0) then
      call MAPL_TimerOff(MAPL,"RUN1" )
      call MAPL_TimerOff(MAPL,"TOTAL")
      RETURN_(ESMF_SUCCESS)
   end if

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
   allocate(US (NT,NUM_SUBTILES),   _STAT)
   allocate(VS (NT,NUM_SUBTILES),   _STAT)
   allocate(TS (NT,NUM_SUBTILES),   _STAT)
   allocate(FR (NT,NUM_SUBTILES),   _STAT)

   call MAPL_GetResource( MAPL, WHICH_T_TO_SFCLAYER, Label="T_from_AOIL_to_SFCLAYER:", DEFAULT="TW_from_internal", _RC)

   call MAPL_GetResource( MAPL, DEPTH_T_TO_SFCLAYER, Label="DEPTH_T_AOIL_to_SFCLAYER:", DEFAULT=0.,                _RC)

   if ( DEPTH_T_TO_SFCLAYER > AOIL_depth) then
     print *, " DEPTH_T_AOIL_to_SFCLAYER must not be greater than the depth of the AOIL, which is currently set =", AOIL_depth, "Exiting!"
     ASSERT_(.false.)
   endif

   call AOIL_sfcLayer_T( WHICH_T_TO_SFCLAYER, DEPTH_T_TO_SFCLAYER, DO_DATASEA, MUSKIN, epsilon_d, &
                          AOIL_depth, TW, TS_FOUNDi, TWMTF, DELTC, TS(:,WATER))

   FR(:,WATER) = 1.0 ! parent(saltwater) will aggregate based on water/ice fraction 

   US(:,WATER) = UW
   VS(:,WATER) = VW

   where (TS(:,WATER) < TF)
       !*** reset TS to freezing point  
       TS(:,WATER) = TF
   endwhere

   TW = TS(:,WATER)

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

   N = WATER  

   sfc_layer: if(CHOOSEMOSFC.eq.0) then
         call louissurface(1,N,UU,WW,PS,TA,TS,QA,QS,PCU,LAI,Z0,DZ,CM,CN,RIB,ZT,ZQ,CH,CQ,UUU,UCN,RE)

   elseif (CHOOSEMOSFC.eq.1) then

         niter = 6   ! number of internal iterations in the helfand MO surface layer routine
         IWATER=1
         Z0(:,N)=OCEANZ0

         PSMB = PS * 0.01            ! convert to MB
         fakelai  = 1.e-4
         ! Approximate pressure at top of surface layer: hydrostatic, eqn of state using avg temp and press
         PSL = PSMB * (1. - (DZ*MAPL_GRAV)/(MAPL_RGAS*(TA+TS(:,N)) ) ) /   &
                      (1. + (DZ*MAPL_GRAV)/(MAPL_RGAS*(TA+TS(:,N)) ) ) 

         call helfsurface( UWINDLMTILE,VWINDLMTILE,TA,TS(:,N),QA,QS(:,N),PSL,PSMB,Z0(:,N),        &
                           fakelai,IWATER,DZ,niter,nt,RHO,VKH,VKM,USTAR,XX,YY,CU,CT,RIB,ZETA,WS,  &
                           t2m,q2m,u2m,v2m,t10m,q10m,u10m,v10m,u50m,v50m,CHOOSEZ0,CHARNOCK_)

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

   deallocate(CHARNOCK_)

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN1" )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)

 end subroutine RUN1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP
! !IROUTINE: RUN2 -- Second Run stage for the Openwater component

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

  real, pointer, dimension(:)         :: AREA => null()     ! needed to calculate TILEAREA in OpenWaterCore

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

    call MAPL_Get(MAPL,                          &
         TILELATS  = LATS ,                      &
         TILELONS  = LONS ,                      &
         TILEAREA  = AREA ,                      &
         ORBIT     = ORBIT,                      &
         INTERNAL_ESMF_STATE = INTERNAL,         &
         CF = CF,                                &
                                       _RC)

! Update the skin variables each step
!------------------------------------

    call OPENWATERCORE(NT=size(LONS), _RC)

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN2" )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine OPENWATERCORE(NT,RC)
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
   real, pointer, dimension(:  )  :: SHWTR   => null()
   real, pointer, dimension(:  )  :: SHOUT   => null()
   real, pointer, dimension(:  )  :: HLATN   => null()
   real, pointer, dimension(:  )  :: HLATWTR => null()
   real, pointer, dimension(:  )  :: HLWUP   => null()
   real, pointer, dimension(:  )  :: SWNDSRF => null()
   real, pointer, dimension(:  )  :: LWNDSRF => null()
   real, pointer, dimension(:  )  :: SWNDWTR => null()
   real, pointer, dimension(:  )  :: LWNDWTR => null()
   real, pointer, dimension(:  )  :: AOSHFLX => null()
   real, pointer, dimension(:  )  :: AOQFLUX => null()
   real, pointer, dimension(:  )  :: AOLWFLX => null()
   real, pointer, dimension(:  )  :: AOSNOW  => null()
   real, pointer, dimension(:  )  :: AORAIN  => null()
   real, pointer, dimension(:  )  :: AODRNIR => null()
   real, pointer, dimension(:  )  :: AODFNIR => null()
   real, pointer, dimension(:  )  :: FSURF   => null()
   real, pointer, dimension(:  )  :: PENOCNe => null()

   real, pointer, dimension(:  )  :: DELTS  => null()
   real, pointer, dimension(:  )  :: DELQS  => null()
   real, pointer, dimension(:  )  :: TST    => null()
   real, pointer, dimension(:  )  :: QST    => null()
   real, pointer, dimension(:  )  :: TAUXW  => null()
   real, pointer, dimension(:  )  :: TAUYW  => null()
   real, pointer, dimension(:  )  :: PENUVR => null()
   real, pointer, dimension(:  )  :: PENUVF => null()
   real, pointer, dimension(:  )  :: PENPAR => null()
   real, pointer, dimension(:  )  :: PENPAF => null()
   real, pointer, dimension(:  )  :: FRI    => null()
   real, pointer, dimension(:  )  :: FRW    => null()

   real, pointer, dimension(:  )  :: Dcool     => null()           ! cool-skin, diurnal warming related exports
   real, pointer, dimension(:  )  :: Dwarm     => null()
   real, pointer, dimension(:  )  :: Tdrop     => null()
   real, pointer, dimension(:  )  :: Tbar      => null()
   real, pointer, dimension(:  )  :: Qcool     => null()
   real, pointer, dimension(:  )  :: SWcool    => null()
   real, pointer, dimension(:  )  :: Qwarm     => null()
   real, pointer, dimension(:  )  :: SWwarm    => null()
   real, pointer, dimension(:  )  :: Phiw      => null()
   real, pointer, dimension(:  )  :: Langm     => null()
   real, pointer, dimension(:  )  :: Ustarw    => null()
   real, pointer, dimension(:  )  :: Lcool     => null()
   real, pointer, dimension(:  )  :: Bcool     => null()
   real, pointer, dimension(:  )  :: Tdel      => null()
   real, pointer, dimension(:  )  :: TAUTW     => null()
   real, pointer, dimension(:  )  :: TS_FOUNDe => null()
   real, pointer, dimension(:  )  :: SS_FOUNDe => null()   
   real, pointer, dimension(:  )  :: ZETA_W    => null()
   real, pointer, dimension(:  )  :: TWMTFe    => null()

! pointers to internal

   real, pointer, dimension(:  )  :: TW    => null()  ! AOIL compatibility
   real, pointer, dimension(:  )  :: HW    => null()
   real, pointer, dimension(:  )  :: SW    => null()

   real, pointer, dimension(:,:)  :: QS    => null()
   real, pointer, dimension(:,:)  :: CH    => null()
   real, pointer, dimension(:,:)  :: CQ    => null()
   real, pointer, dimension(:,:)  :: CM    => null()
   real, pointer, dimension(:  )  :: TWMTS => null()  ! AOIL compatibility

   real, pointer, dimension(:  )  :: TWMTF => null()
   real, pointer, dimension(:  )  :: DELTC => null()

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
   real, pointer, dimension(:)    :: VW => null()
   real, pointer, dimension(:)    :: KPAR => null()
   real, pointer, dimension(:)    :: TS_FOUNDi => null()
   real, pointer, dimension(:)    :: SS_FOUNDi => null()  
   real, pointer, dimension(:)    :: DTSDT => null()
   real, pointer, dimension(:)    :: TF => null()
   real, pointer, dimension(:)    :: DISCHARGE_IM => null()

   real, allocatable                   :: TS (:,:)                  ! Following 4 Variables: TS to FR need to be 
   real, allocatable                   :: HH (:,:)                  ! allocatable because NUM_SUBTILES is NOT a parameter
   real, allocatable                   :: SS (:,:)
   real, allocatable                   :: FR (:,:)
   real,    dimension(NT)              :: FRWATER
   real,    dimension(NT)              :: SHF
   real,    dimension(NT)              :: EVP
   real,    dimension(NT)              :: CFQ
   real,    dimension(NT)              :: CFT
   real,    dimension(NT)              :: TXW
   real,    dimension(NT)              :: TYW
   real,    dimension(NT)              :: DQS
   real,    dimension(NT)              :: DTS
   real,    dimension(NT)              :: SWN
   real,    dimension(NT)              :: SWN_surf
   real,    dimension(NT)              :: PEN
   real,    dimension(NT)              :: PEN_ocean
   real,    dimension(NT)              :: LHF
   real,    dimension(NT)              :: ZTH
   real,    dimension(NT)              :: SLR
   real,    dimension(NT)              :: ALBVRO
   real,    dimension(NT)              :: ALBVFO
   real,    dimension(NT)              :: ALBNRO
   real,    dimension(NT)              :: ALBNFO
   real,    dimension(NT)              :: VSUVR
   real,    dimension(NT)              :: VSUVF

   real,    dimension(NT)              :: LCOOL_
   real,    dimension(NT)              :: BCOOL_
   real,    dimension(NT)              :: USTARW_
   real,    dimension(NT)              :: QCOOL_
   real,    dimension(NT)              :: QWARM_
   real,    dimension(NT)              :: SWCOOL_
   real,    dimension(NT)              :: SWWARM_
   real,    dimension(NT)              :: PHIW_
   real,    dimension(NT)              :: DWARM_
   real,    dimension(NT)              :: TBAR_
   real,    dimension(NT)              :: DCOOL_
   real,    dimension(NT)              :: TDROP_
   real,    dimension(NT)              :: TDEL_
   real,    dimension(NT)              :: LANGM_
   real,    dimension(NT)              :: TAUTW_
   real,    dimension(NT)              :: ZETA_W_
   real,    dimension(NT)              :: uStokes_                  ! Stokes velocity should be an import from Wave Watch

   real                                :: DT
   real                                :: MAXWATERDEPTH
   real                                :: MINWATERDEPTH
   integer                             :: n_iter_cool               ! number of iterations to compute cool-skin layer
   real                                :: fr_ice_thresh             ! threshold on ice fraction, used in diurnal warming and cool-skin layer
   real                                :: MUSKIN                    ! exponent in T(z) profile in warm layer, based on Zeng & Beljaars, 2005
   real                                :: STOKES_SPEED              ! Stokes velocity in m/s. A global constant number until Wave Watch is ready

   real                                :: MAXSALINITY
   real                                :: MINSALINITY
   real                                :: AOIL_depth                ! thickness of atmosphere-ocean interface layer (AOIL) denoted by d in Akella & Suarez, 2018
   real                                :: OGCM_top_thickness        ! thickness of OGCM top layer (D) in AS2018
   real                                :: epsilon_d                 ! ratio: (thickness of AOIL)/(thickness of OGCM top level) = epsilon_d in AS2018
   character(len=3)                    :: DO_GRAD_DECAY_warmLayer   ! simulate gradual decay of warm layer: yes or no. Follows Zeng and Beljaars, 2005.
   character(len=3)                    :: DO_UPDATE_FLUXES_AOIL_SECOND_STEP  ! update fluxes and other variables following second update step in AOIL: yes or no.

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
   logical                             :: debugzth

   real, parameter :: EMSH2O          = 0.99070

!  Begin...
!----------

   IAm =  trim(COMP_NAME) // "OPENWATERCORE"

! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,ALW    , 'ALW'    ,    _RC)
   call MAPL_GetPointer(IMPORT,BLW    , 'BLW'    ,    _RC)
   call MAPL_GetPointer(IMPORT,LWDNSRF, 'LWDNSRF',    _RC)
   call MAPL_GetPointer(IMPORT,DRPAR  , 'DRPAR'  ,    _RC)
   call MAPL_GetPointer(IMPORT,DFPAR  , 'DFPAR'  ,    _RC)
   call MAPL_GetPointer(IMPORT,DRNIR  , 'DRNIR'  ,    _RC)
   call MAPL_GetPointer(IMPORT,DFNIR  , 'DFNIR'  ,    _RC)
   call MAPL_GetPointer(IMPORT,DRUVR  , 'DRUVR'  ,    _RC)
   call MAPL_GetPointer(IMPORT,DFUVR  , 'DFUVR'  ,    _RC)
   call MAPL_GetPointer(IMPORT,EVAP   , 'EVAP'   ,    _RC)
   call MAPL_GetPointer(IMPORT,SH     , 'SH'     ,    _RC)
   call MAPL_GetPointer(IMPORT,TAUX   , 'TAUX'   ,    _RC)
   call MAPL_GetPointer(IMPORT,TAUY   , 'TAUY'   ,    _RC)
   call MAPL_GetPointer(IMPORT,DEV    , 'DEVAP'  ,    _RC)
   call MAPL_GetPointer(IMPORT,DSH    , 'DSH'    ,    _RC)
   call MAPL_GetPointer(IMPORT,SNO    , 'SNO'    ,    _RC)
   call MAPL_GetPointer(IMPORT,PLS    , 'PLS'    ,    _RC)
   call MAPL_GetPointer(IMPORT,PCU    , 'PCU'    ,    _RC)
   call MAPL_GetPointer(IMPORT,PS     , 'PS'     ,    _RC)
   call MAPL_GetPointer(IMPORT,UU     , 'UU'     ,    _RC)

   call MAPL_GetPointer(IMPORT,FI     , 'FRACI'  ,    _RC)
   call MAPL_GetPointer(IMPORT,TF     , 'TFREEZE',    _RC)

   call MAPL_GetPointer(IMPORT,UW     , 'UW'     ,    _RC)
   call MAPL_GetPointer(IMPORT,VW     , 'VW'     ,    _RC)
   call MAPL_GetPointer(IMPORT,THATM  , 'THATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,QHATM  , 'QHATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,UHATM  , 'UHATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,VHATM  , 'VHATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,UUA    , 'UUA'    ,    _RC)
   call MAPL_GetPointer(IMPORT,VVA    , 'VVA'    ,    _RC)
   call MAPL_GetPointer(IMPORT,CTATM  , 'CTATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,CQATM  , 'CQATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,CMATM  , 'CMATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,KPAR   , 'KPAR'   ,    _RC)
   call MAPL_GetPointer(IMPORT,TS_FOUNDi,'TS_FOUND',  _RC)
   call MAPL_GetPointer(IMPORT,SS_FOUNDi,'SS_FOUND',  _RC)
   call MAPL_GetPointer(IMPORT,DTSDT  , 'DTSDT'  ,    _RC)
   call MAPL_GetPointer(IMPORT, DISCHARGE_IM, 'DISCHARGE', _RC)

! Pointers to internals
!----------------------

   call MAPL_GetPointer(INTERNAL,TW     ,'TSKINW',    _RC)
   call MAPL_GetPointer(INTERNAL,HW     ,'HSKINW',    _RC)
   call MAPL_GetPointer(INTERNAL,SW     ,'SSKINW',    _RC)
   call MAPL_GetPointer(INTERNAL,TWMTS  ,'TWMTS',     _RC)
   call MAPL_GetPointer(INTERNAL,TWMTF  ,'TWMTF',     _RC)
   call MAPL_GetPointer(INTERNAL,DELTC  ,'DELTC',     _RC)

   call MAPL_GetPointer(INTERNAL,QS     , 'QS'   ,    _RC)
   call MAPL_GetPointer(INTERNAL,CH     , 'CH'   ,    _RC)
   call MAPL_GetPointer(INTERNAL,CQ     , 'CQ'   ,    _RC)
   call MAPL_GetPointer(INTERNAL,CM     , 'CM'   ,    _RC)

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
   call MAPL_GetPointer(EXPORT,TAUXW  , 'TAUXW', alloc=.true., _RC)
   call MAPL_GetPointer(EXPORT,TAUYW  , 'TAUYW', alloc=.true., _RC)
   call MAPL_GetPointer(EXPORT,PENUVR , 'PENUVR'  ,    _RC)
   call MAPL_GetPointer(EXPORT,PENUVF , 'PENUVF'  ,    _RC)
   call MAPL_GetPointer(EXPORT,PENPAR , 'PENPAR'  ,    _RC)
   call MAPL_GetPointer(EXPORT,PENPAF , 'PENPAF'  ,    _RC)
   call MAPL_GetPointer(EXPORT,EVAPOUT, 'EVAPOUT' ,    _RC)
   call MAPL_GetPointer(EXPORT,SNOWOCN, 'SNOWOCN' ,    _RC)
   call MAPL_GetPointer(EXPORT,RAINOCN, 'RAINOCN' ,    _RC)
   call MAPL_GetPointer(EXPORT,SHOUT  , 'SHOUT'   ,    _RC)
   call MAPL_GetPointer(EXPORT,SHWTR  , 'SHWTR'   ,    _RC)
   call MAPL_GetPointer(EXPORT,HLATN  , 'HLATN'   ,    _RC)
   call MAPL_GetPointer(EXPORT,HLATWTR, 'HLATWTR' ,    _RC)
   call MAPL_GetPointer(EXPORT,HLWUP  , 'HLWUP'   ,    _RC)
   call MAPL_GetPointer(EXPORT,LWNDSRF, 'LWNDSRF' ,    _RC)
   call MAPL_GetPointer(EXPORT,SWNDSRF, 'SWNDSRF' ,    _RC)
   call MAPL_GetPointer(EXPORT,LWNDWTR, 'LWNDWTR' ,    _RC)
   call MAPL_GetPointer(EXPORT,SWNDWTR, 'SWNDWTR' ,    _RC)
   call MAPL_GetPointer(EXPORT,FRW    , 'FRACW'   ,    _RC)
   call MAPL_GetPointer(EXPORT,AOSHFLX, 'AO_SHFLX',    _RC)
   call MAPL_GetPointer(EXPORT,AOQFLUX, 'AO_QFLUX',    _RC)
   call MAPL_GetPointer(EXPORT,AOLWFLX, 'AO_LWFLX',    _RC)
   call MAPL_GetPointer(EXPORT,AOSNOW , 'AO_SNOW' ,    _RC)
   call MAPL_GetPointer(EXPORT,AORAIN , 'AO_RAIN' ,    _RC)
   call MAPL_GetPointer(EXPORT,AODRNIR, 'AO_DRNIR',    _RC)
   call MAPL_GetPointer(EXPORT,AODFNIR, 'AO_DFNIR',    _RC)
   call MAPL_GetPointer(EXPORT,FSURF  , 'FSURF'   ,    _RC)
   call MAPL_GetPointer(EXPORT,PENOCNe, 'PEN_OCN' ,    _RC)

   call MAPL_GetPointer(EXPORT,Dwarm  , 'DWARM'   ,    _RC)
   call MAPL_GetPointer(EXPORT,Dcool  , 'DCOOL'   ,    _RC)
   call MAPL_GetPointer(EXPORT,Tbar   , 'TBAR'    ,    _RC)
   call MAPL_GetPointer(EXPORT,Tdrop  , 'TDROP'   ,    _RC)
   call MAPL_GetPointer(EXPORT,Qcool  , 'QCOOL'   ,    _RC)
   call MAPL_GetPointer(EXPORT,USTARW , 'USTARW'  ,    _RC)
   call MAPL_GetPointer(EXPORT,Lcool  , 'LCOOL'   ,    _RC)
   call MAPL_GetPointer(EXPORT,SWcool , 'SWCOOL'  ,    _RC)
   call MAPL_GetPointer(EXPORT,SWwarm , 'SWWARM'  ,    _RC)
   call MAPL_GetPointer(EXPORT,Qwarm  , 'QWARM'   ,    _RC)
   call MAPL_GetPointer(EXPORT,Phiw   , 'PHIW'    ,    _RC)
   call MAPL_GetPointer(EXPORT,Langm  , 'LANGM'   ,    _RC)
   call MAPL_GetPointer(EXPORT,Bcool  , 'BCOOL'   ,    _RC)
   call MAPL_GetPointer(EXPORT,Tdel   , 'TDEL'    ,    _RC)
   call MAPL_GetPointer(EXPORT,TS_FOUNDe, 'TS_FOUND' , _RC)
   call MAPL_GetPointer(EXPORT,SS_FOUNDe, 'SS_FOUND' , _RC)
   call MAPL_GetPointer(EXPORT,TauTW  , 'TAUTW'   ,    _RC)
   call MAPL_GetPointer(EXPORT,ZETA_W , 'ZETA_W'  ,    _RC)
   call MAPL_GetPointer(EXPORT,TWMTFe , 'TWMTF'   ,    _RC)

   allocate(TS (NT,NUM_SUBTILES),_STAT)
   allocate(HH (NT,NUM_SUBTILES),_STAT)
   allocate(SS (NT,NUM_SUBTILES),_STAT)
   allocate(FR (NT,NUM_SUBTILES),_STAT)

! Get the time step
! -----------------

    call MAPL_Get(MAPL, HEARTBEAT = DT, _RC)
    call MAPL_GetResource ( MAPL, DT, Label="DT:", DEFAULT=DT, _RC)

! Get parameters
! --------------

    if (DO_SKIN_LAYER==0) then
       call MAPL_GetResource ( MAPL, MAXWATERDEPTH, Label="MAX_WATER_DEPTH:" , DEFAULT=1000., _RC)
       call MAPL_GetResource ( MAPL, MINWATERDEPTH, Label="MIN_WATER_DEPTH:" , DEFAULT=1000., _RC)
    else 
       call MAPL_GetResource ( MAPL, MAXWATERDEPTH, Label="MAX_WATER_DEPTH:" , DEFAULT=2.,   _RC)
       call MAPL_GetResource ( MAPL, MINWATERDEPTH, Label="MIN_WATER_DEPTH:" , DEFAULT=2.,   _RC)
    end if

! Exponent in the near-surface temperature profile T(z) within the AOIL
! ---------------------------------------------------------------------
    call MAPL_GetResource ( MAPL, MUSKIN,        Label="MU_SKIN:"         , DEFAULT=0.2  ,   _RC)

! How many cool-skin iterations to do?
! -------------------------------------
    call MAPL_GetResource ( MAPL, n_iter_cool, Label="COOL_SKIN_LAYER_ITERATIONS:" , DEFAULT=3,    _RC)

    AOIL_depth = MAX(MaxWaterDepth, MinWaterDepth)

    if (DO_DATASEA==0) then
!      Thickness of OGCM top level (m)
!      ------------------------------
       call MAPL_GetResource ( MAPL, OGCM_top_thickness, Label="OGCM_TOP_LAYER:" , DEFAULT=10.,   _RC) ! SA: could be an export from GUEST GC
       epsilon_d  = AOIL_depth/OGCM_top_thickness ! < 1. If that is NOT true, AOIL formulation would need revisit; see AS2018
    else
       OGCM_top_thickness = MAPL_UNDEF
       epsilon_d          = 0.0
    end if

!   --------------------------------------------------------------------------------------------------------
!   Treatment of Marginal Ice Zone (MIZ), i.e., threshold on fraction of ice (fraci), to model the SST variations. 
!   One can imagine at least following three possibilities:
!   (i)  SST is NOT allowed to vary within ice extent, 
!        i.e., if fraci    < fr_ice_thresh (= 1.e-11 default), turn AOIL off, set skin SST = TS_FOUND
!   (ii) SST IS allowed to vary with the ice extent, 
!                 FRwater  > fr_ice_thresh (= 0.0 default),    turn AOIL on, only over water, skin SST .ne. TS_FOUND (as it was <= Jason-2_0)
!   (iii) SST is NOT allowed to vary when SST < SST_cold, say, 15C. Turn AOIL off when the water temperature falls below some threshold.
!
!   As already noted, option (ii) was used in versions before and up to Jason-2_0, and (i) was tried in Jason-3_0, which 
!   showed detriment in forecast skill (self verification tests most prominent, and a bit, with respect to ECMWF operations).
!   Hence reverting to option (ii). The final option (iii) has not been tested, just proposed for the sake of completeness.
!   In any case, probably, (iii) will also degrade forecast skill, just as (i) did, because (what SA thinks) -1.7C, set for water temperature is TOO COLD!
!   Unless we understand and model all the processes, we may have to just let diurnal variability (cool-skin+diurnal warming) pick up the tab!
!
!   In Marginal Ice Zone, threshold on fraction: if no LANL CICE, SST IS ALLOWED TO VARY WITHIN ICE EXTENT.
!   
!   ** Revisit when coupled to ocean+sea-ice ** July, 2019.
!   --------------------------------------------------------------------------------------------------------
    call MAPL_GetResource ( MAPL, fr_ice_thresh, Label="THRESHOLD_ICE_FR_SST:" , DEFAULT=0.0,       _RC)   ! i.e., above option (ii)
!   --------------------------------------------------------------------------------------------------------

    call MAPL_GetResource ( MAPL, STOKES_SPEED,  Label="STOKES_VELOCITY:" , DEFAULT=1.E-2,   _RC)

    call MAPL_GetResource ( MAPL, MAXSALINITY,   Label="MAX_SALINITY:" ,    DEFAULT=40.0 ,   _RC)
    call MAPL_GetResource ( MAPL, MINSALINITY,   Label="MIN_SALINITY:" ,    DEFAULT=5.0  ,   _RC)

    call MAPL_GetResource ( MAPL, DO_GRAD_DECAY_warmLayer,   Label="WARM_LAYER_GRAD_DECAY:" ,  DEFAULT="no"  ,   _RC)

    call MAPL_GetResource ( MAPL, DO_UPDATE_FLUXES_AOIL_SECOND_STEP, Label="UPDATE_FLUXES_AOIL_SECOND_STEP:" ,  DEFAULT="no"  ,   _RC)

!   Copy internals into local variables
!   ------------------------------------
    if(DO_SKIN_LAYER==0) then                ! inactive AOIL (OFF)
      HH(:,WATER) = AOIL_depth * water_RHO('salt_water')
      SS(:,WATER) = SW*HW
      TS(:,WATER) = TW

      TWMTS = 0.
      TWMTF = 0.
      DELTC = 0.
    else
      HH(:,WATER) = HW
      SS(:,WATER) = SW*HW
      TS(:,WATER) = TW - TWMTS
    endif

    FR(:,WATER) = 1.0
    FRWATER     = max(1.0 - FI, 0.0)

! Add skin SST analysis increment. This is zero if ANA_TS is false.
!------------------------------------------------------------------
    TS(:,WATER) = TS(:,WATER) + DT*DTSDT

! Update albedo
!--------------
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

    call ALBSEA    (ALBVRO,ALBVFO,ALBNRO,ALBNFO,ZTH)

    call MAPL_TimerOff(MAPL,    "-Albedo")

    call MAPL_TimerOn(MAPL,    "-OpenWater")

! Initialize PAR and UVR beam fluxes
!-----------------------------------
    VSUVR = DRPAR + DRUVR
    VSUVF = DFPAR + DFUVR

! Regardless of interface layer, penetrative solar to ocean always gets surface values 
! ------------------------------------------------------------------------------------
    if(associated(PENUVR)) PENUVR  = (1.-ALBVRO)*DRUVR
    if(associated(PENUVF)) PENUVF  = (1.-ALBVFO)*DFUVR
    if(associated(PENPAR)) PENPAR  = (1.-ALBVRO)*DRPAR
    if(associated(PENPAF)) PENPAF  = (1.-ALBVFO)*DFPAR

! Cycle through sub-tiles doing water and energy budget
!------------------------------------------------------
    if(associated(EVAPOUT)) EVAPOUT = 0.0
    if(associated(SHOUT  )) SHOUT   = 0.0
    if(associated(HLATN  )) HLATN   = 0.0
    if(associated(DELTS  )) DELTS   = 0.0
    if(associated(DELQS  )) DELQS   = 0.0

    CFT = (CH(:,WATER)/CTATM)
    CFQ = (CQ(:,WATER)/CQATM)

! Net Solar insolation (including UV & IR, direct & diffuse) in interface layer
!------------------------------------------------------------------------------
    SWN = (1.-ALBVRO)*VSUVR + (1.-ALBVFO)*VSUVF + &
          (1.-ALBNRO)*DRNIR + (1.-ALBNFO)*DFNIR
    SWN_surf = SWN

    call AOIL_Shortwave_abs (NT, DO_SKIN_LAYER, DO_DATASEA, &
                             AOIL_depth, OGCM_top_thickness, HW, KUVR, KPAR, &
                             ALBVRO, ALBVFO, DRUVR, DFUVR, DRPAR, DFPAR,     &
                             PEN, PEN_ocean)

    SWN   = SWN - PEN
    if (DO_DATASEA == 0) then               ! in coupled mode, first term on RHS of eqn 10 (AS2018)
      SWN   = SWN - (epsilon_d/(1.-epsilon_d))* (PEN-PEN_ocean)
    endif

    call AOIL_v0 (NT, DO_SKIN_LAYER, DO_DATASEA, n_iter_cool, fr_ice_thresh, trim(DO_GRAD_DECAY_warmLayer), &
                  DT, MUSKIN, epsilon_d, MaxWaterDepth, MinWaterDepth, MaxSalinity, MinSalinity,            &
                  STOKES_SPEED, CM(:,WATER), CFT, CFQ, SH, EVAP, DSH, DEV, THATM, QHATM, PS, SNO, PCU+PLS,  &
                  UUA, VVA, UW, VW, FRWATER, SWN, SWN_surf, PEN, PEN_ocean, LWDNSRF, ALW, BLW,              &
                  HH(:,WATER), TS(:,WATER), SS(:,WATER), QS(:,WATER), TS_FOUNDi,                            &
                  DWARM_, TBAR_, USTARW_, DCOOL_, TDROP_, SWCOOL_, QCOOL_, BCOOL_, LCOOL_,                  &
                  TDEL_, SWWARM_, QWARM_, ZETA_W_, PHIW_, LANGM_, TAUTW_, uStokes_, TXW, TYW,               &
                  SHF, LHF, EVP, DTS, DQS, DELTC, HW, TW, SW, TWMTS, TWMTF,                                 &
                  trim(do_update_fluxes_AOIL_second_step))

! Copies for export
!------------------
    if(associated(HLATWTR)) then
          WHERE( FRWATER>0.0 )
             HLATWTR = LHF
          ELSEWHERE
             HLATWTR = MAPL_UNDEF
          ENDWHERE
    endif
    if(associated(  SHWTR)) then
          WHERE( FRWATER>0.0 )
             SHWTR   = SHF
          ELSEWHERE
             SHWTR   = MAPL_UNDEF
          ENDWHERE
    endif

    if(associated(EVAPOUT)) EVAPOUT = EVP    *FR(:,WATER)
    if(associated(SHOUT  )) SHOUT   = SHF    *FR(:,WATER)
    if(associated(HLATN  )) HLATN   = LHF    *FR(:,WATER)
    if(associated(DELTS  )) DELTS   = DTS*CFT*FR(:,WATER)
    if(associated(DELQS  )) DELQS   = DQS*CFQ*FR(:,WATER)

    if(associated(SWcool)) SWcool = SWCOOL_
    if(associated(SWWARM)) SWWARM = SWWARM_
    if(associated(Qcool )) Qcool  = QCOOL_
    if(associated(QWARM )) QWARM  = QWARM_
    if(associated(PHIW  )) PHIW   = PHIW_
    if(associated(LANGM )) LANGM  = LANGM_

    if(associated(Dcool )) Dcool  = DCOOL_
    if(associated(Dwarm )) Dwarm  = AOIL_depth  ! exports from warm-layer
    if(associated(Tdrop )) Tdrop  = TDROP_
    if(associated(Tbar  )) Tbar   = TBAR_
    if(associated(Ustarw)) Ustarw = USTARW_
    if(associated(Lcool )) Lcool  = LCOOL_
    if(associated(Bcool )) Bcool  = BCOOL_
    if(associated(Tdel  )) Tdel   = TDEL_
    if(associated(TauTW )) TauTW  = TAUTW_
    if(associated(ZETA_W)) ZETA_W  = ZETA_W_
    if(associated(TS_FOUNDe)) TS_FOUNDe = TS_FOUNDi
    if(associated(SS_FOUNDe)) SS_FOUNDe = SS_FOUNDi

    if(associated(TAUXW)) TAUXW = TXW
    if(associated(TAUYW)) TAUYW = TYW

    if(associated(AOSHFLX)) AOSHFLX = SHF    *FRWATER
    if(associated(AOQFLUX)) AOQFLUX = EVP    *FRWATER
    if(associated(AOLWFLX)) AOLWFLX = (LWDNSRF-ALW-BLW*TS(:,WATER))*FRWATER
    if(associated(AORAIN )) AORAIN  = PCU + PLS
    if(associated(AOSNOW )) AOSNOW  = SNO    *FRWATER
    if(associated(AODRNIR)) AODRNIR = (1.-ALBNRO)*DRNIR*FRWATER
    if(associated(AODFNIR)) AODFNIR = (1.-ALBNFO)*DFNIR*FRWATER
    if(associated(FSURF  )) FSURF   = SWN+LWDNSRF-(ALW+BLW*TS(:,WATER))-SHF-LHF
    if(associated(PENOCNe)) PENOCNe = PEN_ocean * FRWATER

    if(associated(SNOWOCN)) SNOWOCN = SNO*FR(:,WATER)
    if(associated(RAINOCN)) RAINOCN = PCU + PLS
    if(associated(HLWUP  )) HLWUP   = ALW*FR(:,WATER) 
    if(associated(LWNDSRF)) LWNDSRF = (LWDNSRF - ALW)*FR(:,WATER)

    if(associated(LWNDWTR)) then
          where( FRWATER>0.0 )
             LWNDWTR = LWDNSRF - ALW - BLW*TS(:,WATER)
          elsewhere
             LWNDWTR = MAPL_UNDEF
          end where
    endif

    if(associated(TST    )) TST     = TS(:,WATER)*FR(:,WATER)
    if(associated(QST    )) QST     = QS(:,WATER)*FR(:,WATER)
    if(associated(LWNDSRF)) LWNDSRF = LWNDSRF - BLW*TS(:,WATER)*FR(:,WATER)
    if(associated(HLWUP  )) HLWUP   = HLWUP   + BLW*TS(:,WATER)*FR(:,WATER)

    EMISS = EMSH2O*FR(:,WATER)
    ALBVR = ALBVRO*FR(:,WATER)
    ALBVF = ALBVFO*FR(:,WATER)
    ALBNR = ALBNRO*FR(:,WATER)
    ALBNF = ALBNFO*FR(:,WATER)

    if(associated(SWNDWTR)) then 
          where( FRWATER>0.0 ) 
             SWNDWTR = (1.-ALBVRO)*VSUVR + (1.-ALBVFO)*VSUVF + &
                       (1.-ALBNRO)*DRNIR + (1.-ALBNFO)*DFNIR 
          elsewhere
             SWNDWTR = MAPL_UNDEF
          end where
    end if

    if(associated(SWNDSRF)) then 
       SWNDSRF = &
           (1.-ALBVR)*VSUVR + (1.-ALBVF)*VSUVF + &
           (1.-ALBNR)*DRNIR + (1.-ALBNF)*DFNIR
    endif

    call MAPL_GetPointer(IMPORT,FI     , 'FRACINEW'  ,    _RC)
    if(associated(FRW)) then
       FRW = max(1.0 - FI, 0.0)
    endif

    call MAPL_TimerOff(MAPL,   "-OpenWater")
!------------------------------------------------------

    call MAPL_TimerOn(MAPL,    "-Albedo")

    if(solalarmison) then
       call MAPL_SunGetInsolation(LONS, LATS,      &
            ORBIT, ZTH, SLR,                       &
            INTV = TINT,                           &
            currTime=CURRENT_TIME+DELT,            &
            _RC )

       ZTH = max(0.0,ZTH)

       call ALBSEA    (ALBVRO,ALBVFO,ALBNRO,ALBNFO,ZTH)

       ALBVR = ALBVRO*FR(:,WATER)
       ALBVF = ALBVFO*FR(:,WATER)
       ALBNR = ALBNRO*FR(:,WATER)
       ALBNF = ALBNFO*FR(:,WATER)
          
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
             
  end subroutine OPENWATERCORE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine RUN2

end module GEOS_OpenwaterGridCompMod

