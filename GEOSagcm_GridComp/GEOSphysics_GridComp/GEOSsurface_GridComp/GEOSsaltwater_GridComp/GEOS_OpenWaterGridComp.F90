
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
!      \noindent $\sigma_T$ evolves according to
!      $$ \frac{\partial \sigma_T}{\partial t}= \frac{Q_{\sigma}}{d\,\rho_w\,c_w} - \frac{1}{\tau_{\sigma}} \sigma_T.$$
!      \noindent For complete details, please see Akella and Suarez, 2018, 
!      "The Atmosphere-Ocean Interface Layer of the NASA Goddard Earth Observing System Model and Data Assimilation System." GMAO Tech Memo, Vol 51. 
!
!      \noindent COMPATIBILITY: 
!
!      **************************************************************************
!      Set Services:
!      enables ability to run old and/or new interface(s) by setting compatibility: ON
!      -- this feature adds (back) "old" stuff when set ON
!
!      Whereas Run(1,2) will execute either 
!      -----------      -------------------
!      interface 
!       version        compatibility
!     -----------      -------------------
!        old:               ON
!        new:               OFF
!      **************************************************************************
!
!      ----------------------------------------------------------------------------
!

! !USES:

  use sfclayer  ! using module that contains sfc layer code
  use ESMF
  use MAPL
  use GEOS_UtilsMod
  use DragCoefficientsMod
  

  implicit none
  private

  integer            :: DO_DATASEA   ! between DO_GUEST and DO_DATASEA, we pick one that
                                     ! works in both amip and coupled mode
  integer            :: DO_SKIN_LAYER

  public SetServices

!EOP

  integer, parameter            :: WATER        = 1      
  integer, parameter            :: NUM_SUBTILES = 1             ! number of sub-tiles
  real,    parameter            :: KUVR         = 0.09
  real,    parameter            :: SALTWATERCAP  = MAPL_CAPWTR

  character(len=7)   :: AOIL_COMP_SWITCH  ! Atmosphere-Ocean Interface Layer, compatibility: on/off
                                          ! defualt: OFF, so AOIL is incompatible with "old" interface
                                          ! when it is ON, set servives provides what is needed for both versions
                                          ! whereas RUN(1,2) will use one or other (OFF: AOIL; ON: old interface)
                                                               

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

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run1, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2, RC=STATUS )
    VERIFY_(STATUS)

! Get constants from CF
! ---------------------

    call MAPL_GetResource ( MAPL, DO_DATASEA,    Label="USE_DATASEA:"     , DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, DO_SKIN_LAYER, Label="USE_SKIN_LAYER:"  , DEFAULT=0    , RC=STATUS)
    VERIFY_(STATUS)

! Atmosphere-Ocean Interface Layer compatibility: on/off?
!-------------------------------------------------------

    call MAPL_GetResource( MAPL,  AOIL_COMP_SWITCH,        Label="AOIL_COMP_SWITCH:",     DEFAULT="ON", RC=STATUS)
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
        LONG_NAME          = 'ocean_snowfall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'SNOWOCN'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'ocean_rainfall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'RAINOCN'                   ,&
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
        LONG_NAME          = 'open_water_upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHWTR'                     ,&
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
        LONG_NAME          = 'open_water_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWNDWTR'                   ,&
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
        LONG_NAME          = 'open_water_net_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWNDWTR'                   ,&
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
        LONG_NAME          = 'open_water_latent_energy_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLATWTR'                   ,&
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
        SHORT_NAME         = 'FRACW',                             &
        LONG_NAME          = 'water_covered_fraction_of_tile',      &
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
        LONG_NAME          = 'eastward_stress_over_water',&
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUXW'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'northward_stress_over_water',&
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUYW'                     ,&
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

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'DCOOL',                     &
          LONG_NAME          = 'depth_of_cool_layer'       ,&
          UNITS              = 'm'                         ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                          ,&
          SHORT_NAME         = 'DWARM',                      &
          LONG_NAME          = 'depth_at_base_of_warm_layer',&
          UNITS              = 'm'                          ,&
          DIMS               = MAPL_DimsTileOnly            ,&
          VLOCATION          = MAPL_VLocationNone           ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                                 ,&
          SHORT_NAME         = 'TDROP',                             &
          LONG_NAME          = 'temperature_drop_across_cool_layer',&
          UNITS              = 'K'                                 ,&
          DIMS               = MAPL_DimsTileOnly                   ,&
          VLOCATION          = MAPL_VLocationNone                  ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'QCOOL',                     &
          LONG_NAME          = 'net_heating_in_cool_layer' ,&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'QWARM',                     &
          LONG_NAME          = 'net_heating_in_warm_layer' ,&
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                          ,&
          SHORT_NAME         = 'SWCOOL',                     &
          LONG_NAME          = 'solar_heating_in_cool_layer',&
          UNITS              = 'W m-2'                      ,&
          DIMS               = MAPL_DimsTileOnly            ,&
          VLOCATION          = MAPL_VLocationNone           ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                          ,&
          SHORT_NAME         = 'SWWARM',                     &
          LONG_NAME          = 'solar_heating_in_warm_layer',&
          UNITS              = 'W m-2'                      ,&
          DIMS               = MAPL_DimsTileOnly            ,&
          VLOCATION          = MAPL_VLocationNone           ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                           ,&
          SHORT_NAME         = 'PHIW',                        &
          LONG_NAME          = 'Similarity_function_in_warm_layer',&
          UNITS              = '1'                           ,&
          DIMS               = MAPL_DimsTileOnly             ,&
          VLOCATION          = MAPL_VLocationNone            ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'LANGM'                     ,&
          LONG_NAME          = 'Langmuir_number'           ,&
          UNITS              = '1'                         ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                        ,&
          SHORT_NAME         = 'USTARW',                    &
          LONG_NAME          = 'ustar_over_water'          ,&
          UNITS              = 'm s-1'                     ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                                  ,&
          SHORT_NAME         = 'TBAR',                               &
          LONG_NAME          = 'mean_temperature_of_interface_layer',&
          UNITS              = 'K'                                  ,&
          DIMS               = MAPL_DimsTileOnly                    ,&
          VLOCATION          = MAPL_VLocationNone                   ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'LCOOL',                     &
          LONG_NAME          = 'Saunders_parameter'        ,&
          UNITS              = '1'                         ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'BCOOL',                     &
          LONG_NAME          = 'bouyancy_generation_in_cool_layer',&
          UNITS              = 'm+2 s-3'                   ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                        &
          SHORT_NAME         = 'TDEL',                     &
          LONG_NAME          = 'temperature_at_base_of_cool_layer', &
          UNITS              = 'K',                        &
          DIMS               = MAPL_DimsTileOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                        &
          SHORT_NAME         = 'TS_FOUND',                 &
          LONG_NAME          = 'foundation_temperature_for_interface_layer',        &
          UNITS              = 'K',                        &
          DIMS               = MAPL_DimsTileOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                        &
          SHORT_NAME         = 'TAUTW',                    &
          LONG_NAME          = 'relaxation_time_of_TW_to_TS_FOUND', &
          UNITS              = 's',                        &
          DIMS               = MAPL_DimsTileOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          SHORT_NAME         = 'ZETA_W'                    ,&
          LONG_NAME          = 'Stability_parameter_in_Warm_Layer',                 &
          UNITS              = '1'                         ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

! this is gone because FRZMLT will be from ocean/guest
!     call MAPL_AddExportSpec(GC,                     &
!          SHORT_NAME         = 'FRZMLT'                    ,&
!          LONG_NAME          = 'Freeze_melt_potential',     &
!          UNITS              = 'W m-2'                     ,&
!          DIMS               = MAPL_DimsTileOnly           ,&
!          VLOCATION          = MAPL_VLocationNone          ,&
!                                               RC=STATUS  ) 
!    VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                           &
          SHORT_NAME         = 'SS_FOUND',                 &
          LONG_NAME          = 'foundation_salinity_for_interface_layer',        &
          UNITS              = 'PSU',                      &
          DIMS               = MAPL_DimsTileOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
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

     if( trim(AOIL_COMP_SWITCH) == "OFF") then ! as close as possible to "x0039", while keeping everything as in "x0040"
        call MAPL_AddExportSpec(GC,                    &
             LONG_NAME          = 'saturation_specific_humidity_using_geos_formula',&
             UNITS              = 'kg kg-1'                   ,&
             SHORT_NAME         = 'QSAT1'                     ,&
             DIMS               = MAPL_DimsTileOnly           ,&
             VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
        VERIFY_(STATUS)
 
        call MAPL_AddExportSpec(GC,                    &
             LONG_NAME          = 'saturation_specific_humidity_using_bulk_formula',&
             UNITS              = 'kg kg-1'                   ,&
             SHORT_NAME         = 'QSAT2'                     ,&
             DIMS               = MAPL_DimsTileOnly           ,&
             VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
        VERIFY_(STATUS)
     endif

     ! atmosphere-ocean fluxes
     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'atmosphere_ocean_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'AO_SHFLX'                  ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'atmosphere_ocean_evaporation' ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'AO_QFLUX'                  ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'atmosphere_ocean_net_longwave_radiation' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'AO_LWFLX'                  ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'atmosphere_ocean_snowfall' ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'AO_SNOW'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'atmosphere_ocean_rainfall' ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'AO_RAIN'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC                         ,&
         LONG_NAME          = 'net_surface_downwelling_nir_beam_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'AO_DRNIR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC                         ,&
         LONG_NAME          = 'net_surface_downwelling_nir_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'AO_DFNIR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         LONG_NAME          = 'departure_of_mean_interface_temperature_from_foundation_temperature',   &
         UNITS              = 'K',                                 &
         SHORT_NAME         = 'TWMTF',                             &
         DIMS               = MAPL_DimsTileOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

!  !INTERNAL STATE:

     if( trim(AOIL_COMP_SWITCH) == "ON") then ! as close as possible to "x0039", while keeping everything as in "x0040"
     
        call MAPL_AddInternalSpec(GC,                           &
             SHORT_NAME         = 'HSKINW',                            &
             LONG_NAME          = 'water_skin_layer_mass',             &
             UNITS              = 'kg m-2',                            &
             DIMS               = MAPL_DimsTileOnly,                   &
             VLOCATION          = MAPL_VLocationNone,                  &
             FRIENDLYTO         = 'OCEAN:SEAICE',                      &
             DEFAULT            = 5.0*MAPL_RHO_SEAWATER,               &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddInternalSpec(GC,                           &
             SHORT_NAME         = 'TSKINW',                            &
             LONG_NAME          = 'water_skin_temperature',            &
             UNITS              = 'K',                                 &
             DIMS               = MAPL_DimsTileOnly,                   &
             VLOCATION          = MAPL_VLocationNone,                  &
             FRIENDLYTO         = 'OCEAN:SEAICE',                      &
             DEFAULT            = 280.0,                               &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddInternalSpec(GC,                           &
             SHORT_NAME         = 'SSKINW',                            &
             LONG_NAME          = 'water_skin_salinity',               &
             UNITS              = 'psu',                               &
             DIMS               = MAPL_DimsTileOnly,                   &
             VLOCATION          = MAPL_VLocationNone,                  &
             FRIENDLYTO         = 'OCEAN:SEAICE',                      &
             DEFAULT            = 30.0,                                &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

     endif

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

     if( trim(AOIL_COMP_SWITCH) == "ON") then ! as close as possible to "x0039", while keeping everything as in "x0040"
        call MAPL_AddInternalSpec(GC,                                &
             SHORT_NAME         = 'TWMTS',                             &
             LONG_NAME          = 'departure_of_skin_temperature_from_mean_interface_temperature',   &
             UNITS              = 'K',                                 &
             DEFAULT            = 0.0,                                 &
             DIMS               = MAPL_DimsTileOnly,                   &
             VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)
     endif

     call MAPL_AddInternalSpec(GC,                                &
        SHORT_NAME         = 'TWMTF',                             &
        LONG_NAME          = 'departure_of_mean_interface_temperature_from_foundation_temperature',   &
        UNITS              = 'K',                                 &
        DEFAULT            = 0.0,                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'DELTC',                             &
        LONG_NAME          = 'temperature_drop_across_cool_layer',&
        UNITS              = 'K',                                 &
        DEFAULT            = 0.0,                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
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
        SHORT_NAME         = 'FRACI',                             &
        LONG_NAME          = 'ice_covered_fraction_of_tile',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'FRACINEW',                             &
        LONG_NAME          = 'ice_covered_fraction_of_tile_after_update',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'UW',                                &
        LONG_NAME          = 'zonal_velocity_of_surface_water',   &
        UNITS              = 'm s-1 ',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VW',                                &
        LONG_NAME          = 'meridional_velocity_of_surface_water',   &
        UNITS              = 'm s-1 ',                            &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'KPAR',                              &
        LONG_NAME          = 'PAR_extinction_coefficient',        &
        UNITS              = 'm-1',                               &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'TS_FOUND',                        &
          LONG_NAME          = 'foundation_temperature_for_interface_layer',&
          UNITS              = 'K',                               &
          DIMS               = MAPL_DimsTileOnly,                 &
          VLOCATION          = MAPL_VLocationNone,                &
          DEFAULT            = 280.0,                             &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'SS_FOUND',                        &
          LONG_NAME          = 'foundation_salinity_for_interface_layer',&
          UNITS              = 'PSU',                             &
          DIMS               = MAPL_DimsTileOnly,                 &
          VLOCATION          = MAPL_VLocationNone,                &
          DEFAULT            = 30.0,                              &
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

    call MAPL_AddImportSpec (GC,                                   &
         SHORT_NAME = 'DTSDT',                                     &
         LONG_NAME  = 'skin_temperature_analysis_tendency',        &
         UNITS      = 'K s-1',                                     &
         RESTART    = MAPL_RestartSkip,                            &
         DIMS       = MAPL_DimsTileOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                    &
          LONG_NAME          = 'river_discharge_at_ocean_points',&
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'DISCHARGE'                 ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartSkip            ,&
          RC=STATUS  ) 
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


!    if( trim(AOIL_COMP_SWITCH) == "ON") then ! as close as possible to "x0039", while keeping everything as in "x0040"

! these exports were filled with 0. by SimpleSeaiceGridComp
!       call MAPL_AddImportSpec(GC,                             &
!            SHORT_NAME         = 'DRUVRTHRU',                             &
!            LONG_NAME          = 'penetrative_uvr_beam_flux_through_sea_ice',           &
!            UNITS              = 'W m-2',                             &
!            DIMS               = MAPL_DimsTileOnly,                   &
!            VLOCATION          = MAPL_VLocationNone,                  &
!            RESTART            = MAPL_RestartSkip,                    &
!            RC=STATUS  )
!       VERIFY_(STATUS)

!       call MAPL_AddImportSpec(GC,                             &
!            SHORT_NAME         = 'DFUVRTHRU',                             &
!            LONG_NAME          = 'penetrative_uvr_diffuse_flux_through_sea_ice',           &
!            UNITS              = 'W m-2',                             &
!            DIMS               = MAPL_DimsTileOnly,                   &
!            VLOCATION          = MAPL_VLocationNone,                  &
!            RESTART            = MAPL_RestartSkip,                    &
!            RC=STATUS  )
!       VERIFY_(STATUS)

!       call MAPL_AddImportSpec(GC,                             &
!            SHORT_NAME         = 'DRPARTHRU',                             &
!            LONG_NAME          = 'penetrative_par_beam_flux_through_sea_ice',           &
!            UNITS              = 'W m-2',                             &
!            DIMS               = MAPL_DimsTileOnly,                   &
!            VLOCATION          = MAPL_VLocationNone,                  &
!            RESTART            = MAPL_RestartSkip,                    &
!            RC=STATUS  )
!       VERIFY_(STATUS)

!       call MAPL_AddImportSpec(GC,                             &
!            SHORT_NAME         = 'DFPARTHRU',                             &
!            LONG_NAME          = 'penetrative_par_diffuse_flux_through_sea_ice',           &
!            UNITS              = 'W m-2',                             &
!            DIMS               = MAPL_DimsTileOnly,                   &
!            VLOCATION          = MAPL_VLocationNone,                  &
!            RESTART            = MAPL_RestartSkip,                    &
!            RC=STATUS  )
!       VERIFY_(STATUS)

!       call MAPL_AddImportSpec(GC,                             &
!            SHORT_NAME         = 'FRESH',                             &
!            LONG_NAME          = 'fresh_water_flux_from_sea_ice',     &
!            UNITS              = 'kg m-2 s-1',                        &
!            DIMS               = MAPL_DimsTileOnly,                   &
!            VLOCATION          = MAPL_VLocationNone,                  &
!            RESTART            = MAPL_RestartSkip,                    &
!            RC=STATUS  )
!       VERIFY_(STATUS)

!       call MAPL_AddImportSpec(GC,                             &
!            SHORT_NAME         = 'FSALT',                             &
!            LONG_NAME          = 'salt_flux_from_sea_ice',            &
!            UNITS              = 'kg m-2 s-1',                        &
!            DIMS               = MAPL_DimsTileOnly,                   &
!            VLOCATION          = MAPL_VLocationNone,                  &
!            RESTART            = MAPL_RestartSkip,                    &
!            RC=STATUS  )
!       VERIFY_(STATUS)

!       call MAPL_AddImportSpec(GC,                             &
!            SHORT_NAME         = 'FHOCN',                             &
!            LONG_NAME          = 'heat_flux_from_sea_ice',            &
!            UNITS              = 'W m-2',                             &
!            DIMS               = MAPL_DimsTileOnly,                   &
!            VLOCATION          = MAPL_VLocationNone,                  &
!            RESTART            = MAPL_RestartSkip,                    &
!            RC=STATUS  )
!       VERIFY_(STATUS)

!    endif

     ! added here in case skin layer needs it
     !call MAPL_AddImportSpec(GC,                     &
     !   SHORT_NAME         = 'FRZMLT'                    ,&
     !   LONG_NAME          = 'freeze_melt_potential',     &
     !   UNITS              = 'W m-2'                     ,&
     !   DIMS               = MAPL_DimsTileOnly           ,&
     !   VLOCATION          = MAPL_VLocationNone          ,&
     !   DEFAULT            = 0.0                         ,&
     !   RC=STATUS  ) 
     !VERIFY_(STATUS)

!EOS


! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="RUN1"  ,                RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN2"  ,                RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-OpenWater",            RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_TimerAdd(GC,    name="-Albedo"  ,             RC=STATUS)
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
   real, pointer, dimension(:  )  :: QSAT1  => null()
   real, pointer, dimension(:  )  :: QSAT2  => null()

! pointers to internal

   real, pointer, dimension(:  )  :: TW  => null() ! AOIL compatibility
   real, pointer, dimension(:  )  :: HW  => null()

   real, pointer, dimension(:  )  :: TWMTF => null()
   real, pointer, dimension(:  )  :: DELTC => null()
   real, pointer, dimension(:,:)  :: QS  => null()
   real, pointer, dimension(:,:)  :: CH  => null()
   real, pointer, dimension(:,:)  :: CM  => null()
   real, pointer, dimension(:,:)  :: CQ  => null()
   real, pointer, dimension(:,:)  :: WW  => null()
   real, pointer, dimension(:,:)  :: Z0  => null()

   real, pointer, dimension(:  )  :: TWMTS => null()  ! add here to be able to fix the "bug" with old interface

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
!  real, pointer, dimension(:)    :: FRZMLT    => null()

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
   real            :: QSAT_SCL

   integer         :: iFIX_BUG1                 ! whether to fix the "bug" with old interface? 0: no (default), 1: yes
   integer         :: iMAK_BUG                  ! whether to put a   "bug" in   new interface? 0: no (default), 1: yes
   character(len=ESMF_MAXSTR)     :: SURFRC
   type(ESMF_Config)              :: SCF 

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

    call MAPL_Get(MAPL, HEARTBEAT = DT, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, DT, Label="DT:", DEFAULT=DT, RC=STATUS)
    VERIFY_(STATUS)

! Get parameters (0:Louis, 1:Monin-Obukhov)
! -----------------------------------------
    call MAPL_GetResource (MAPL, SURFRC, label = 'SURFRC:', default = 'GEOS_SurfaceGridComp.rc', RC=STATUS) ; VERIFY_(STATUS)
    SCF = ESMF_ConfigCreate(rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigLoadFile     (SCF,SURFRC,rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute (SCF, label='CHOOSEMOSFC:', value=CHOOSEMOSFC, DEFAULT=1, __RC__ ) 
    call ESMF_ConfigDestroy      (SCF, __RC__)

    call MAPL_GetResource ( MAPL, CHOOSEZ0,    Label="CHOOSEZ0:",    DEFAULT=3, RC=STATUS)
    VERIFY_(STATUS)

! Get roughness parameters 
! -------------------------------------------------------------
    call MAPL_GetResource ( MAPL, OCEANZ0,     Label="OCEANZ0:" ,    DEFAULT=1.0e-3, RC=STATUS) 
    VERIFY_(STATUS)

! Get Thickness of AOIL (m)
! -------------------------
    if (DO_SKIN_LAYER==0) then
       call MAPL_GetResource ( MAPL, MaxWaterDepth, Label="MAX_WATER_DEPTH:" , DEFAULT=1000., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetResource ( MAPL, MinWaterDepth, Label="MIN_WATER_DEPTH:" , DEFAULT=1000., RC=STATUS)
       VERIFY_(STATUS)
    else 
       call MAPL_GetResource ( MAPL, MaxWaterDepth, Label="MAX_WATER_DEPTH:" , DEFAULT=2.,   RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetResource ( MAPL, MinWaterDepth, Label="MIN_WATER_DEPTH:" , DEFAULT=2.,   RC=STATUS)
       VERIFY_(STATUS)

!      Exponent in the near-surface temperature profile T(z) within the AOIL
!      ---------------------------------------------------------------------
       call MAPL_GetResource ( MAPL, MUSKIN,        Label="MU_SKIN:"         , DEFAULT=0.2  ,   RC=STATUS)
       VERIFY_(STATUS)
    end if

    AOIL_depth = MAX(MaxWaterDepth, MinWaterDepth)

    if (DO_DATASEA==0) then
!      Thickness of OGCM top level (m)
!      ------------------------------
       call MAPL_GetResource ( MAPL, OGCM_top_thickness, Label="OGCM_TOP_LAYER:" , DEFAULT=10.,   RC=STATUS) ! SA: could be an export from GUEST GC
       VERIFY_(STATUS)

       epsilon_d  = AOIL_depth/OGCM_top_thickness ! < 1. If that is NOT true, AOIL formulation would need revisit; see AS2018
    end if

! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,UU     , 'UU'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,UWINDLMTILE     , 'UWINDLMTILE'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VWINDLMTILE     , 'VWINDLMTILE'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,UW     , 'UW'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VW     , 'VW'     ,    RC=STATUS)
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
   call MAPL_GetPointer(IMPORT,FI     , 'FRACI'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TF     , 'TFREEZE',    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TS_FOUNDi, 'TS_FOUND', RC=STATUS)
   VERIFY_(STATUS)
!  call MAPL_GetPointer(IMPORT,FRZMLT, 'FRZMLT'  ,    RC=STATUS)
!  VERIFY_(STATUS)

! Pointers to internals
!----------------------

   if( trim(AOIL_COMP_SWITCH) == "ON") then ! as close as possible to "x0039", while keeping everything as in "x0040"
     call MAPL_GetPointer(INTERNAL,TW   , 'TSKINW' ,    RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(INTERNAL,HW   , 'HSKINW' ,    RC=STATUS)
     VERIFY_(STATUS)
   endif

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

   if( trim(AOIL_COMP_SWITCH) == "OFF") then ! as close as possible to "x0039", while keeping everything as in "x0040"
     call MAPL_GetPointer(EXPORT,QSAT1 , 'QSAT1'   ,    RC=STATUS)
     VERIFY_(STATUS)

     call MAPL_GetPointer(EXPORT,QSAT2 , 'QSAT2'   ,    RC=STATUS)
     VERIFY_(STATUS)
   endif

   if( trim(AOIL_COMP_SWITCH) == "ON") then ! as close as possible to "x0039", while keeping everything as in "x0040"
      call MAPL_GetResource ( MAPL, iFIX_BUG1,  Label="FIX1_OLD_INTERFACE:"  , DEFAULT=0,   RC=STATUS)
      VERIFY_(STATUS)
   else
      call MAPL_GetResource ( MAPL, iMAK_BUG,   Label="BREAK_NEW_INTERFACE:" , DEFAULT=0,   RC=STATUS)
      VERIFY_(STATUS)
   endif

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

   if( trim(AOIL_COMP_SWITCH) == "ON") then ! as close as possible to "x0039", while keeping everything as in "x0040"
     TS(:,WATER) = TW  ! TW (in K): returned from MOM/dataocean. Max, should it be TS(:,WATER) = TW - TWMTS?

     if( iFIX_BUG1 == 1) then
       call MAPL_GetPointer(INTERNAL,TWMTS  , 'TWMTS',    RC=STATUS); VERIFY_(STATUS)
       TS(:,WATER)= TS(:,WATER) - TWMTS ! so helfsurface gets SKIN SST (i.e., "TS") as input.
     endif

   else

     if (DO_SKIN_LAYER==0) then
       TS(:,WATER) = TS_FOUNDi
     else
       call MAPL_GetPointer(INTERNAL,TWMTF, 'TWMTF',  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL,DELTC, 'DELTC',  RC=STATUS); VERIFY_(STATUS)

       if (DO_DATASEA == 1) then                                           ! Ocean is from "data"
          TS(:,WATER)= TS_FOUNDi + ((1.+MUSKIN)/MUSKIN) * TWMTF            ! Eqn.(14) of AS2018
       else                                                                ! Ocean is from a model
          TS(:,WATER)= TS_FOUNDi + (1./MUSKIN + (1.-epsilon_d)) * TWMTF    ! RHS is from Eqn.(15) of AS2018; (here) abuse of notation: T_o is from OGCM.
       endif
       TS(:,WATER) = TS(:,WATER) - DELTC                                   ! Eqn.(16) of AS2018

       if( iMAK_BUG == 1) then
         TS(:,WATER) = TS_FOUNDi + TWMTF                                   ! deliberately input helfsurface "TW"
       endif
     endif
   endif


   FR(:,WATER) = 1.0 ! parent(saltwater) will aggregate based on water/ice fraction 

   US(:,WATER) = UW
   VS(:,WATER) = VW

   where (TS(:,WATER) < TF)
       !*** reset TS to freezing point  
       TS(:,WATER) = TF
   endwhere

   if( trim(AOIL_COMP_SWITCH) == "ON") then ! as close as possible to "x0039", while keeping everything as in "x0040"
     TW = TS(:,WATER)
   endif

   if( trim(AOIL_COMP_SWITCH) == "OFF") then ! as close as possible to "x0039", while keeping everything as in "x0040"
   ! it is a minor bug to not having the line below
   ! TS has been updated on the ocean side, so QS should too
   ! commented out for now so Santha can have zero-diff

     !QS(:,WATER) = GEOS_QSAT(TS(:,WATER), PS, RAMP=0.0, PASCALS=.TRUE.) 
     call MAPL_GetResource ( MAPL, QSAT_SCL, Label="QSAT_SALTWATER_SCALING:" , DEFAULT=1.0, RC=STATUS)
     VERIFY_(STATUS)
     QS(:,WATER) = QSAT_SCL*GEOS_QSAT(TS(:,WATER), PS, RAMP=2.0, PASCALS=.TRUE.) 

     if(associated(QSAT1)) QSAT1 = QS(:,WATER)
     !if(associated(QSAT2)) QSAT2 = 1.0/1.22*0.98*640380.0*exp(-5107.4/TS(:,WATER))
   endif

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
                           t2m,q2m,u2m,v2m,t10m,q10m,u10m,v10m,u50m,v50m,CHOOSEZ0)

         CM(:,N)  = VKM
         CH(:,N)  = VKH
         CQ(:,N)  = VKH

         CN = (MAPL_KARMAN/ALOG(DZ/Z0(:,N) + 1.0)) * (MAPL_KARMAN/ALOG(DZ/Z0(:,N) + 1.0))
         ZT = Z0(:,N)
         ZQ = Z0(:,N)
         RE = 0.
         UUU = UU

         if( trim(AOIL_COMP_SWITCH) == "ON") then ! as close as possible to "x0039", while keeping everything as in "x0040"
           UCN = 0.
         endif


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

         if( trim(AOIL_COMP_SWITCH) == "OFF") then ! as close as possible to "x0039", while keeping everything as in "x0040"
           if(associated(QSAT2)) QSAT2 = 1.0/RHO*0.98*640380.0*exp(-5107.4/TS(:,WATER))
         endif
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
         if(associated(GST)) GST     = WW(:,N)*FR(:,N)
         if(associated(VNT)) VNT     = UUU    *FR(:,N)


      !  Aggregate effective, CD-weighted, surface values of T and Q
      !-------------------------------------------------------------

         if(associated(TH)) TH      = CH(:,N)*TS(:,N)*FR(:,N)
         if(associated(QH)) QH      = CQ(:,N)*QS(:,N)*FR(:,N)
         if(associated(UH)) UH      = CM(:,N)*US(:,N)*FR(:,N)
         if(associated(VH)) VH      = CM(:,N)*VS(:,N)*FR(:,N)


         WW(:,N) = max(CH(:,N)*(TS(:,N)-TA-(MAPL_GRAV/MAPL_CP)*DZ)/TA + MAPL_VIREPS*CQ(:,N)*(QS(:,N)-QA),0.0)
         WW(:,N) = (HPBL*MAPL_GRAV*WW(:,N))**(2./3.)

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

  real, parameter                     :: puny = 1.e-11

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

    call OPENWATERCORE(NT=size(LONS), RC=STATUS )
    VERIFY_(STATUS)

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
   real, allocatable                   :: tmp2 (:,:)
   real,    dimension(NT)              :: FRWATER
   real,    dimension(NT)              :: SHF
   real,    dimension(NT)              :: EVP
   real,    dimension(NT)              :: SHD
   real,    dimension(NT)              :: EVD
   real,    dimension(NT)              :: CFQ
   real,    dimension(NT)              :: CFT
   !real,    dimension(NT)              :: UUA
   !real,    dimension(NT)              :: VVA
   real,    dimension(NT)              :: TXW
   real,    dimension(NT)              :: TYW
   real,    dimension(NT)              :: DQS
   real,    dimension(NT)              :: DTS
   real,    dimension(NT)              :: DTX
   real,    dimension(NT)              :: DTY
   real,    dimension(NT)              :: SWN
   real,    dimension(NT)              :: PEN
   real,    dimension(NT)              :: PEN_ocean
   real,    dimension(NT)              :: PUR
   real,    dimension(NT)              :: PUF
   real,    dimension(NT)              :: PPR
   real,    dimension(NT)              :: PPF
   real,    dimension(NT)              :: PENICE
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
   real,    dimension(NT)              :: FRESHATM                  ! fresh water flux from atmosphere

   real,    dimension(NT)              :: ALPH

   integer                             :: N
   real                                :: DT
   real                                :: MAXWATERDEPTH
   real                                :: MINWATERDEPTH
   integer                             :: n_iter_cool               ! number of iterations to compute cool-skin layer
   real                                :: fr_ice_thresh             ! threshold on ice fraction, used in diurnal warming and cool-skin layer
   real                                :: MUSKIN                    ! exponent in T(z) profile in warm layer, based on Zeng & Beljaars, 2005
   real                                :: STOKES_SPEED              ! Stokes velocity in m/s. A global constant number until Wave Watch is ready
   integer                             :: USE_KPAR

   real                                :: MAXSALINITY
   real                                :: MINSALINITY
   real                                :: AOIL_depth                ! thickness of atmosphere-ocean interface layer (AOIL) denoted by d in Akella & Suarez, 2018
   real                                :: OGCM_top_thickness        ! thickness of OGCM top layer (D) in AS2018
   real                                :: epsilon_d                 ! ratio: (thickness of AOIL)/(thickness of OGCM top level) = epsilon_d in AS2018
   real                                :: F_PHI                     ! tunable parameter, used to calculate stability function in AS2018
   real                                :: QSAT_SCL 
   real                                :: fLA
   character(len=10)                   :: diurnal_warming_scheme    ! which formulation of diurnal warming model? AS2018 or AT2017 (DOI:10.1002/qj.2988)

! following are related  to CICE

   integer                             :: NSUB, I, K, L

   real                                :: FSALT, FRESH!, FHOCN

   integer                             :: iFIX_BUG2                 ! whether to fix the "bug" with old interface? 0: no (default), 1: yes


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

   call MAPL_GetPointer(IMPORT,FI     , 'FRACI'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TF     , 'TFREEZE',    RC=STATUS); VERIFY_(STATUS)

   call MAPL_GetPointer(IMPORT,UW     , 'UW'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VW     , 'VW'     ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,THATM  , 'THATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,QHATM  , 'QHATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,UHATM  , 'UHATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VHATM  , 'VHATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,UUA    , 'UUA'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VVA    , 'VVA'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,CTATM  , 'CTATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,CQATM  , 'CQATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,CMATM  , 'CMATM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,KPAR   , 'KPAR'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TS_FOUNDi,'TS_FOUND',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,SS_FOUNDi,'SS_FOUND',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DTSDT  , 'DTSDT'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, DISCHARGE_IM, 'DISCHARGE', RC=STATUS); VERIFY_(STATUS)

! Pointers to internals
!----------------------

   if( trim(AOIL_COMP_SWITCH) == "ON") then ! as close as possible to "x0039", while keeping everything as in "x0040"
     call MAPL_GetPointer(INTERNAL,TW     ,'TSKINW',    RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(INTERNAL,HW     ,'HSKINW',    RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(INTERNAL,SW     ,'SSKINW',    RC=STATUS); VERIFY_(STATUS)
   endif

   call MAPL_GetPointer(INTERNAL,QS     , 'QS'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CH     , 'CH'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CQ     , 'CQ'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CM     , 'CM'   ,    RC=STATUS); VERIFY_(STATUS)

   if( trim(AOIL_COMP_SWITCH) == "ON") then ! as close as possible to "x0039", while keeping everything as in "x0040"
     call MAPL_GetPointer(INTERNAL,TWMTS  , 'TWMTS',    RC=STATUS); VERIFY_(STATUS)
   endif


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
   call MAPL_GetPointer(EXPORT,TAUXW  , 'TAUXW', alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TAUYW  , 'TAUYW', alloc=.true., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,PENUVR , 'PENUVR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,PENUVF , 'PENUVF'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,PENPAR , 'PENPAR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,PENPAF , 'PENPAF'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,EVAPOUT, 'EVAPOUT' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SNOWOCN, 'SNOWOCN' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RAINOCN, 'RAINOCN' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SHOUT  , 'SHOUT'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SHWTR  , 'SHWTR'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLATN  , 'HLATN'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLATWTR, 'HLATWTR' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,HLWUP  , 'HLWUP'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,LWNDSRF, 'LWNDSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SWNDSRF, 'SWNDSRF' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,LWNDWTR, 'LWNDWTR' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SWNDWTR, 'SWNDWTR' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FRW    , 'FRACW'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,AOSHFLX, 'AO_SHFLX',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,AOQFLUX, 'AO_QFLUX',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,AOLWFLX, 'AO_LWFLX',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,AOSNOW , 'AO_SNOW' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,AORAIN , 'AO_RAIN' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,AODRNIR, 'AO_DRNIR',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,AODFNIR, 'AO_DFNIR',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,FSURF  , 'FSURF'   ,    RC=STATUS); VERIFY_(STATUS)


   call MAPL_GetPointer(EXPORT,Dwarm  , 'DWARM'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Dcool  , 'DCOOL'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Tbar   , 'TBAR'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Tdrop  , 'TDROP'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Qcool  , 'QCOOL'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,USTARW , 'USTARW'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Lcool  , 'LCOOL'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SWcool , 'SWCOOL'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SWwarm , 'SWWARM'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Qwarm  , 'QWARM'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Phiw   , 'PHIW'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Langm  , 'LANGM'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Bcool  , 'BCOOL'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Tdel   , 'TDEL'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TS_FOUNDe, 'TS_FOUND' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,SS_FOUNDe, 'SS_FOUND' , RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TauTW  , 'TAUTW'   ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ZETA_W , 'ZETA_W'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TWMTFe , 'TWMTF'   ,    RC=STATUS); VERIFY_(STATUS)


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

    if (DO_SKIN_LAYER==0) then
       call MAPL_GetResource ( MAPL, MAXWATERDEPTH, Label="MAX_WATER_DEPTH:" , DEFAULT=1000., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetResource ( MAPL, MINWATERDEPTH, Label="MIN_WATER_DEPTH:" , DEFAULT=1000., RC=STATUS)
       VERIFY_(STATUS)
    else 
       call MAPL_GetResource ( MAPL, MAXWATERDEPTH, Label="MAX_WATER_DEPTH:" , DEFAULT=2.,   RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetResource ( MAPL, MINWATERDEPTH, Label="MIN_WATER_DEPTH:" , DEFAULT=2.,   RC=STATUS)
       VERIFY_(STATUS)
    end if

! Exponent in the near-surface temperature profile T(z) within the AOIL
! ---------------------------------------------------------------------
    call MAPL_GetResource ( MAPL, MUSKIN,        Label="MU_SKIN:"         , DEFAULT=0.2  ,   RC=STATUS)
    VERIFY_(STATUS)

! How many cool-skin iterations to do?
! -------------------------------------
    call MAPL_GetResource ( MAPL, n_iter_cool, Label="COOL_SKIN_LAYER_ITERATIONS" , DEFAULT=3,    RC=STATUS)
    VERIFY_(STATUS)

    AOIL_depth = MAX(MaxWaterDepth, MinWaterDepth)

    if (DO_DATASEA==0) then
!      Thickness of OGCM top level (m)
!      ------------------------------
       call MAPL_GetResource ( MAPL, OGCM_top_thickness, Label="OGCM_TOP_LAYER:" , DEFAULT=10.,   RC=STATUS) ! SA: could be an export from GUEST GC
       VERIFY_(STATUS)

       epsilon_d  = AOIL_depth/OGCM_top_thickness ! < 1. If that is NOT true, AOIL formulation would need revisit; see AS2018
    end if

!   --------------------------------------------------------------------------------------------------------
!   Treatment of Marginal Ice Zone (MIZ), i.e., threshold on fraction of ice (fraci), to model the SST variations. 
!   One can imagine at least following three possibilities:
!   (i)  SST is NOT allowed to vary within ice extent, 
!        i.e., if fraci    < fr_ice_thresh (1.e-11 default), turn AOIL off, set skin SST = TS_FOUND
!   (ii) SST IS allowed to vary with the ice extent, 
!                 FRwater  > fr_ice_thresh (0.0 default),    turn AOIL on, only over water, skin SST .ne. TS_FOUND (as it was <= Jason-2_0)
!   (iii) SST is NOT allowed to vary when SST < SST_cold, say, 15C. Turn AOIL off when the water temperature falls below some threshold.
!
!   As already noted, option (ii) was used in versions before and up to Jason-2_0, and (i) was tried in Jason-3_0, which 
!   showed detriment in forecast skill (self verification tests most prominent, and a bit, with respect to ECMWF operations).
!   Hence reverting to option (ii). The final option (iii) has not been tested, just proposed for the sake of completeness.
!   In any case, probably, (iii) will also degrade forecast skill, just as (i) did, because (what SA thinks) -1.7C, set for water temperature is TOO COLD!
!   Unless we understand and model all the processes, we may have to just let diurnal variability (cool-skin+diurnal warming) pick up the tab!
!   
!   ** Revisit when coupled to ocean+sea-ice ** July, 2019.
!   --------------------------------------------------------------------------------------------------------
!   call MAPL_GetResource ( MAPL, fr_ice_thresh, Label="THRESHOLD_ICE_FR_SST:" , DEFAULT=1.e-11,    RC=STATUS)   ! above option (i)
    call MAPL_GetResource ( MAPL, fr_ice_thresh, Label="THRESHOLD_ICE_FR_SST:" , DEFAULT=0.0,       RC=STATUS)   ! above option (ii)
    VERIFY_(STATUS)
!   --------------------------------------------------------------------------------------------------------

    call MAPL_GetResource ( MAPL, STOKES_SPEED,  Label="STOKES_VELOCITY:" , DEFAULT=1.E-2,   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, USE_KPAR,      Label="USE_KPAR:",         DEFAULT=1,       RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, MAXSALINITY,   Label="MAX_SALINITY:" ,    DEFAULT=40.0 ,   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, MINSALINITY,   Label="MIN_SALINITY:" ,    DEFAULT=5.0  ,   RC=STATUS)
    VERIFY_(STATUS)

    if( trim(AOIL_COMP_SWITCH) == "ON") then ! as close as possible to "x0039", while keeping everything as in "x0040"
      MaxWaterDepth   = MaxWaterDepth*MAPL_RHOWTR
      MinWaterDepth   = MinWaterDepth*MAPL_RHOWTR
      if(DO_SKIN_LAYER==0) TWMTS = 0.
    else
      MaxWaterDepth   = MaxWaterDepth*MAPL_RHO_SEAWATER
      MinWaterDepth   = MinWaterDepth*MAPL_RHO_SEAWATER
    endif

    if( trim(AOIL_COMP_SWITCH) == "ON") then ! as close as possible to "x0039", while keeping everything as in "x0040"

!   Copy friendly internals into tile-tile local variables
!   -------------------------------------------------------

      HH(:,WATER) = HW
      SS(:,WATER) = SW*HW
      TS(:,WATER) = TW - TWMTS

    else

!   Get TS from internal state exactly as in Run1
!   ---------------------------------------------
       if (DO_SKIN_LAYER==0) then 
         HH(:,WATER) = 1.e+15    ! infinite heat capacity with TS = SST (from data)
         TS(:,WATER) = TS_FOUNDi
         SS(:,WATER) = SS_FOUNDi
         TWMTF       = 0.
         DELTC       = 0.
       else 
         HH(:,WATER) = AOIL_depth*MAPL_RHO_SEAWATER
         SS(:,WATER) = SS_FOUNDi*HH(:,WATER)

         call MAPL_GetPointer(INTERNAL,TWMTF, 'TWMTF',  RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(INTERNAL,DELTC, 'DELTC',  RC=STATUS); VERIFY_(STATUS)

         if (DO_DATASEA == 1) then                                          ! Ocean is from "data"
           TS(:,WATER)= TS_FOUNDi + ((1.+MUSKIN)/MUSKIN) * TWMTF            ! Eqn.(14) of AS2018
         else                                                               ! Ocean is from a model
           TS(:,WATER)= TS_FOUNDi + (1./MUSKIN + (1.-epsilon_d)) * TWMTF    ! RHS is from Eqn.(15) of AS2018; (here) abuse of notation: T_o is from OGCM.
         endif
         TS(:,WATER) = TS(:,WATER) - DELTC                                  ! Eqn.(16) of AS2018
       endif
    endif ! if( trim(AOIL_COMP_SWITCH) == "ON")

    FR(:,WATER) = 1.0
    FRWATER     = max(1.0 - FI, 0.0)

! Initialize PAR and UVR beam fluxes
!-----------------------------------

    VSUVR = DRPAR + DRUVR
    VSUVF = DFPAR + DFUVR
    PUR   = 0.0
    PUF   = 0.0
    PPR   = 0.0
    PPF   = 0.0

! Add analysis increment. This is zero if ANA_TS is false.
!---------------------------------------------------------
    TS(:,WATER) = TS(:,WATER) + DT*DTSDT

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

    call ALBSEA    (ALBVRO,ALBVFO,ALBNRO,ALBNFO,ZTH)

    call MAPL_TimerOff(MAPL,    "-Albedo")

    call MAPL_TimerOn(MAPL,    "-OpenWater")
! Cycle through sub-tiles doing water and energy budget
!------------------------------------------------------

    if(associated(EVAPOUT)) EVAPOUT = 0.0
    if(associated(SHOUT  )) SHOUT   = 0.0
    if(associated(HLATN  )) HLATN   = 0.0
    if(associated(DELTS  )) DELTS   = 0.0
    if(associated(DELQS  )) DELQS   = 0.0

    if(associated(PENUVR))  PENUVR  = 0.0 ! Fill up exports (following 4) related to shortwave absorption
    if(associated(PENUVF))  PENUVF  = 0.0 
    if(associated(PENPAR))  PENPAR  = 0.0 
    if(associated(PENPAF))  PENPAF  = 0.0 

    N   = WATER
    CFT = (CH(:,N)/CTATM)
    CFQ = (CQ(:,N)/CQATM)
    EVP = CFQ*(EVAP + DEV*(QS(:,N)-QHATM))
    SHF = CFT*(SH   + DSH*(TS(:,N)-THATM))
    SHD = CFT*DSH
    EVD = CFQ*DEV*GEOS_DQSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.)
    DTS = LWDNSRF - (ALW + BLW*TS(:,N)) - SHF

    ! FR accounts for skin under ice
    DTX = DT*FRWATER / (SALTWATERCAP*HH(:,N))

    if (DO_DATASEA == 1) then            ! in uncoupled mode
      if (DO_SKIN_LAYER /= 0) DTX = DTX*((MUSKIN+1.)/MUSKIN)
    else
      if (DO_SKIN_LAYER /= 0) DTX = DTX*((MUSKIN+1.-MUSKIN*epsilon_d)/MUSKIN)
    endif

!   Shortwave absorption in water
!   -----------------------------
    if (DO_SKIN_LAYER==0) then 
      PEN = 1.0
    else 
      if( trim(AOIL_COMP_SWITCH) == "ON") then
        PEN = exp(-(KUVR/MAPL_RHOWTR)*HW)           ! replace MAPL_RHOWTR with MAPL_RHO_SEAWATER
      else
        PEN = exp(-KUVR*AOIL_depth)                 ! UV penetration into water
      endif
    endif
    PUR = (1.-ALBVRO)*DRUVR*PEN
    PUF = (1.-ALBVFO)*DFUVR*PEN

    if (DO_SKIN_LAYER==0) then 
      PEN = 1.0
    else 
      if( trim(AOIL_COMP_SWITCH) == "ON") then
        PEN = exp(-(KPAR/MAPL_RHOWTR)*HW)           ! replace MAPL_RHOWTR with MAPL_RHO_SEAWATER
      else
        PEN = exp(-KPAR*AOIL_depth)                 ! near-IR ("blue light" 490nm?)
      endif
    endif
    PPR = (1.-ALBVRO)*DRPAR*PEN
    PPF = (1.-ALBVFO)*DFPAR*PEN

    PEN = PUR + PUF + PPR + PPF          ! total absorbed into water up to AOIL_depth
    
    if (DO_DATASEA == 0) then            ! in coupled mode
      PEN_ocean =  exp(-KUVR*OGCM_top_thickness)
      PUR       =  (1.-ALBVRO)*DRUVR*PEN_ocean
      PUF       =  (1.-ALBVFO)*DFUVR*PEN_ocean
      PEN_ocean =  exp(-KPAR*OGCM_top_thickness)
      PPR       =  (1.-ALBVRO)*DRPAR*PEN_ocean
      PPF       =  (1.-ALBVFO)*DFPAR*PEN_ocean
      PEN_ocean =  PUR + PUF + PPR + PPF    ! total absorbed into water up to OGCM_top_thickness
    else
      PEN_ocean = 0.0
    endif

    SWN = (1.-ALBVRO)*VSUVR + (1.-ALBVFO)*VSUVF + &
          (1.-ALBNRO)*DRNIR + (1.-ALBNFO)*DFNIR


    if( trim(AOIL_COMP_SWITCH) == "ON") then
      ikpar: if ( USE_KPAR /= 1) then    ! (if NOT default) compute penetrated shortwave using below...
        call SIMPLE_SW_ABS(NT, USE_KPAR, (HW/MAPL_RHOWTR), ZTH, SWN, PEN)  ! replace MAPL_RHOWTR with MAPL_RHO_SEAWATER
      end if ikpar
    endif
    SWN   = SWN - PEN

    if (DO_DATASEA == 0) then            ! in coupled mode
      SWN   = SWN - (epsilon_d/(1.-epsilon_d))* (PEN-PEN_ocean)
    endif

    ! regardless of what interface layer is used, penetrative solar to ocean
    ! always gets the surface values 
    if(associated(PENUVR)) PENUVR  = (1.-ALBVRO)*DRUVR
    if(associated(PENUVF)) PENUVF  = (1.-ALBVFO)*DFUVR
    if(associated(PENPAR)) PENPAR  = (1.-ALBVRO)*DRPAR
    if(associated(PENPAF)) PENPAF  = (1.-ALBVFO)*DFPAR

    ! DTY accounts for ice on top of water. Part of Shortwave is absorbed by ice and rest goes to warm water.
    ! Skin layer only absorbs the portion of SW radiation passing thru the bottom of ice MINUS
    ! the portion passing thru the skin layer    

    ! penetrated shortwave from sea ice bottom + associated ocean/ice heat flux
!    DTY = DT / (SALTWATERCAP*HW) * (PENICE * FI + FHOCN)
     DTY = 0. ! SA: revisit above with CICE6 [Nov, 2019]

    DTS = DTX * ( DTS + SWN - EVP*MAPL_ALHL - MAPL_ALHF*SNO ) + DTY
    DTS = DTS   / ( 1.0 + DTX*(BLW + SHD + EVD*MAPL_ALHL) )
    EVP = EVP + EVD * DTS
    SHF = SHF + SHD * DTS
    LHF = EVP * MAPL_ALHL

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

! Update WATER surface temperature and moisture
!----------------------------------------

    TS(:,N) = TS(:,N) + DTS
    DQS     = GEOS_QSAT(TS(:,N), PS, RAMP=0.0, PASCALS=.TRUE.) - QS(:,N)
    QS(:,N) = QS(:,N) + DQS

    if( trim(AOIL_COMP_SWITCH) == "ON") then
      if(DO_SKIN_LAYER/=0) then
        TWMTS = TWMTS - (1.0/(MUSKIN+1.0))*DTS
      endif
    else
      if(DO_SKIN_LAYER/=0) then
       if (DO_DATASEA == 0) then            ! in coupled mode
         TWMTF = TWMTF + (MUSKIN/(1.+MUSKIN-MUSKIN*epsilon_d))*DTS
       else
         TWMTF = TWMTF + (MUSKIN/(1.+MUSKIN))*DTS
       end if
      end if
    end if

! Layer thickness; liquid precip goes right thru ice.
! FRESHATM is useful for mass flux balance.
! freshwater flux from atmosphere needs to be added to HW
! here since it carries zero enthalpy 
!---------------------------------------------------

    FRESHATM    = FRWATER*(SNO - EVP) + PCU + PLS
    FRESH       = 0.

    if( trim(AOIL_COMP_SWITCH) == "ON") then
      HH(:,N) = HH(:,N) + DT*(FRESHATM + FRESH)
      HH(:,N) = max(min(HH(:,N),MaxWaterDepth),MinWaterDepth)
    endif

    if(associated(EVAPOUT)) EVAPOUT = EVP    *FR(:,N)
    if(associated(SHOUT  )) SHOUT   = SHF    *FR(:,N)
    if(associated(HLATN  )) HLATN   = LHF    *FR(:,N)
    if(associated(DELTS  )) DELTS   = DTS*CFT*FR(:,N)
    if(associated(DELQS  )) DELQS   = DQS*CFQ*FR(:,N)

    if( trim(AOIL_COMP_SWITCH) == "ON") then
!    Copy back to friendly internal variables
!    -----------------------------------------

     TW = TS(:,WATER) + TWMTS       ! SA: I don't fully understand how LANL CICE should interact w/Skin Layer. Jul 2015.
     HW = HH(:,WATER)               ! So for now, have the same skin layer interaction as in GEOS CICE
     ! multiply by 1000 to account for g->kg conversion
     FSALT = 0.  ! SA: has to be done this way for compatibility (when it is "ON")
     SW = (SS(:,WATER)+DT*1.e3*FSALT)/HW
     where (.not. (abs(UW) > 0.0 .or. abs(VW) > 0.0))
        SW = max(min(SW,MAXSALINITY),MINSALINITY)
     endwhere
    endif

! Atmospheric surface stresses
!-----------------------------
    !*** already computed in surface
    !UUA = (TAUX/CMATM + UHATM)
    !VVA = (TAUY/CMATM + VHATM)

! Net Solar insolation (including UV & IR, direct & diffuse) in interface layer 
!------------------------------------------------------------------------------

    SWN = (1.-ALBVRO)*VSUVR + (1.-ALBVFO)*VSUVF + &
          (1.-ALBNRO)*DRNIR + (1.-ALBNFO)*DFNIR


! how many cool-skin iterations to do?
!-------------------------------------

    call MAPL_GetResource ( MAPL, n_iter_cool, Label="COOL_SKIN_LAYER_ITERATIONS" , DEFAULT=3,    RC=STATUS)
    VERIFY_(STATUS)

    if( trim(AOIL_COMP_SWITCH) == "ON") then
!     Marginal Ice Zone- threshold on fraction: if no LANL CICE, SST IS ALLOWED TO VARY WITHIN ICE EXTENT.
!     -------------------------------------------------------------------------------------------------

      ! Bin: the default for coupled (with cice) needs to be revisited 
      call MAPL_GetResource ( MAPL, fr_ice_thresh, Label="THRESHOLD_ICE_FR_SST:" , DEFAULT=0.0,    RC=STATUS)
      VERIFY_(STATUS)

!     Cool-skin and diurnal warm layer. It changes TS, TWMTS, TW if DO_SKIN_LAYER = 1
!     --------------------------------------------------------------------------------

      allocate(tmp2(NT, NUM_SUBTILES), STAT=STATUS)
      VERIFY_(STATUS)
      tmp2(:,WATER) = FRWATER

      call SKIN_SST (DO_SKIN_LAYER,NT,CM,UUA,VVA,UW,VW,HW,SWN,LHF,SHF,LWDNSRF,                   &
                     ALW,BLW,PEN,STOKES_SPEED,DT,MUSKIN,TS_FOUNDi,DWARM_,TBAR_,TXW,TYW,USTARW_,  &
                     DCOOL_,TDROP_,SWCOOL_,QCOOL_,BCOOL_,LCOOL_,TDEL_,SWWARM_,QWARM_,ZETA_W_,    &
                     PHIW_,LANGM_,TAUTW_,uStokes_,TS,TWMTS,TW,WATER,tmp2,n_iter_cool,&
                     fr_ice_thresh)

      deallocate(tmp2)

      call MAPL_GetResource ( MAPL, iFIX_BUG2,  Label="FIX2_OLD_INTERFACE:"  , DEFAULT=0,   RC=STATUS)
      VERIFY_(STATUS)
      if( iFIX_BUG2 == 1) then
        DQS         = GEOS_QSAT(TS(:,WATER), PS, RAMP=0.0, PASCALS=.TRUE.) - QS(:,WATER)
        QS(:,WATER) = QS(:,WATER) + DQS
      endif
    endif

    if( trim(AOIL_COMP_SWITCH) == "OFF") then

!     Diagnose cool-skin
!     ---------------------

      call COOL_SKIN (NT,CM,UUA,VVA,UW,VW,SWN,LHF,SHF,LWDNSRF,    &
                      ALW,BLW,TXW,TYW,USTARW_,                    &
                      DCOOL_,TDROP_,SWCOOL_,QCOOL_,BCOOL_,LCOOL_, &
                      TS,WATER,FR,n_iter_cool,fr_ice_thresh)

!     AOIL temperature update due to bottom turbulent flux
!     ----------------------------------------------------

      if (DO_SKIN_LAYER==0) then
        TBAR_   = TS(:,WATER)
        TDEL_   = TS(:,WATER) + TDROP_                   ! for analysis not to die if do_skin_layer is off
        TWMTF   = 0.
        DELTC   = 0.

        SWWARM_ = MAPL_UNDEF
        QWARM_  = MAPL_UNDEF
        ZETA_W_ = MAPL_UNDEF
        PHIW_   = MAPL_UNDEF
        LANGM_  = MAPL_UNDEF
        TAUTW_  = MAPL_UNDEF
      else

!       Formulation of diurnal warming scheme
!       -------------------------------------
        call MAPL_GetResource ( MAPL, diurnal_warming_scheme,        Label="DIURNAL_WARMING_SCHEME:" , DEFAULT='AS2018',   RC=STATUS)
        VERIFY_(STATUS)

!       Parameter in stability function, Eq.(29) of AS2018
!       --------------------------------------------------
        call MAPL_GetResource ( MAPL, F_PHI,        Label="F_PHI:" , DEFAULT=3.0  ,   RC=STATUS)
        VERIFY_(STATUS)

        do N = 1, NT
!         if( FI(N)        < fr_ice_thresh ) then   ! see above note on threshold of MIZ to model SST variations
          if( FR(N, WATER) > fr_ice_thresh ) then
            ALPH(N)   = (0.6 + 0.0935*(TS(N,WATER)-MAPL_TICE))*1.E-4
            SWWARM_(N)= SWN(N) - PEN(N)
            QWARM_ (N)= SWWARM_(N) - (LHF(N) + SHF(N) - (LWDNSRF(N) - ALW(N) - BLW(N)*TS(N,WATER)))
            if (DO_DATASEA == 0) then
              QWARM_(N)= QWARM_(N) - (epsilon_d/(1.-epsilon_d))* (PEN(N)-PEN_ocean(N))
            endif
            ZETA_W_(N)= (AOIL_depth*MAPL_KARMAN/USTARW_(N)**3.)*MAPL_GRAV*ALPH(N)*QWARM_(N)/(MAPL_RHO_SEAWATER*MAPL_CAPWTR)

            if ( trim(diurnal_warming_scheme) == 'ATS2017') then ! See Eq(6) of ATS 2017, DOI:10.1002/qj.2988
              if ( ZETA_W_(N) >= 0.0) then   ! Takaya: Eq(5) or Eq(6) of ATS2017
                PHIW_(N) = 1. + (5*ZETA_W_(N) + 4.*ZETA_W_(N)**2)/(1+3.*ZETA_W_(N)+0.25*ZETA_W_(N)**2)
              else
                PHIW_(N) = 1.0/sqrt(1.-16.*ZETA_W_(N))
              end if
            else ! Following implements AS2018
              PHIW_(N)  = 1.+SQRT(1.+4.*MAPL_KARMAN**2.*(1.+MUSKIN)*F_PHI*AOIL_depth*MAPL_GRAV*ALPH(N)*MAX(TWMTF(N),0.0)/USTARW_(N)**2.)
              PHIW_(N)  = 0.5*PHIW_(N)
            endif

            LANGM_(N) = SQRT( USTARW_(N)/STOKES_SPEED)

            if ( trim(diurnal_warming_scheme) == 'ATS2017') then   ! See Eq(6) of ATS 2017, DOI:10.1002/qj.2988
              fLA           = LANGM_(N)**(-0.66667)        ! Takaya: Eqn(6)
              if (fLA       <= 1.0) fLA = 1.0              ! Limit range of fLa to be >=1
              if (ZETA_W_(N)<= 0.0) fLA = 1.0              ! Apply fLa to stable conditions only 
              TAUTW_(N) = (AOIL_depth*PHIW_(N))/(MAPL_KARMAN*USTARW_(N)*(1.+MUSKIN)*fLA)
            else ! AS2018
              TAUTW_(N) = (AOIL_depth*PHIW_(N))/(MAPL_KARMAN*USTARW_(N)*(1.+MUSKIN))
            endif

            if (DO_DATASEA == 0) then
              TAUTW_(N)= (1.- epsilon_d) * TAUTW_(N)   ! compare \tau_{\sigma} in Eq.(22) and that in section 2.2 of AS2018
            endif
            TAUTW_(N) = MAX(DT, TAUTW_(N)) ! for this time-scale, avoid 0.
            TWMTF(N) = TWMTF(N)/(1.+DT/TAUTW_(N))
!           SA: positivity should be imposed like this to be able to compute above PHIW_(N), that needs TWMTF>0 to arg to SQRT() is <0.
!           if (trim(diurnal_warming_scheme) == 'AS2018') then
!             TWMTF(N) = max( TWMTF(N), 0.)
!           endif

            DELTC(N) = TDROP_(N)

            if (DO_DATASEA == 1) then
              TDEL_(N) = TS_FOUNDi(N) + ((1.+MUSKIN)/MUSKIN) * MAX(TWMTF(N), 0.0)
              TBAR_(N) = TS_FOUNDi(N) +                        MAX(TWMTF(N), 0.0)
            else
              TDEL_(N) = TS_FOUNDi(N) + (1./MUSKIN + (1.-epsilon_d)) * MAX(TWMTF(N), 0.0)
              TBAR_(N) = TS_FOUNDi(N) +              (1.-epsilon_d)  * MAX(TWMTF(N), 0.0)
            endif
            !DTS(N)     = ( TDEL_(N) - TDROP_(N)) - TS(N,WATER)
            TS(N,WATER) = TDEL_(N) - TDROP_(N)                  ! updated skin temperature
          else ! FR(N, WATER) <= fr_ice_thresh
            SWWARM_(N) = MAPL_UNDEF
            QWARM_ (N) = MAPL_UNDEF
            ZETA_W_(N) = MAPL_UNDEF
            PHIW_  (N) = MAPL_UNDEF
            LANGM_ (N) = MAPL_UNDEF
            TAUTW_ (N) = MAPL_UNDEF
            TWMTF  (N) = 0.
            DELTC  (N) = 0.

            TBAR_  (N) = MAPL_UNDEF
            TDEL_  (N) = MAPL_UNDEF
          endif ! if( FR(N, WATER) > fr_ice_thresh )
        end do
        ! associated change in QS
        call MAPL_GetResource ( MAPL, iFIX_BUG2,  Label="FIX2_OLD_INTERFACE:"  , DEFAULT=0,   RC=STATUS)
        VERIFY_(STATUS)
        if( iFIX_BUG2 == 1) then
          DQS         = GEOS_QSAT( TS(:,WATER), PS, RAMP=0.0, PASCALS=.TRUE.) - QS(:,WATER)
          QS(:,WATER) = QS(:,WATER) + DQS
        endif
        !if(associated(DELTS  )) DELTS   = DTS*CFT*FR(:,WATER)
        !if(associated(DELQS  )) DELQS   = DQS*CFQ*FR(:,WATER)
      endif ! if (DO_SKIN_LAYER==0) 
      if(associated(TWMTFe)) TWMTFe = TWMTF
    endif   ! if( trim(AOIL_COMP_SWITCH) == "OFF")

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

! Copies for export
!------------------

    if(associated(SNOWOCN)) SNOWOCN = SNO*FR(:,WATER)

    if(associated(AOSHFLX)) AOSHFLX = SHF    *FRWATER
    if(associated(AOQFLUX)) AOQFLUX = EVP    *FRWATER
    if(associated(AOLWFLX)) AOLWFLX = (LWDNSRF-ALW-BLW*TS(:,WATER))*FRWATER
    if(associated(AORAIN )) AORAIN  = PCU + PLS
    if(associated(AOSNOW )) AOSNOW  = SNO    *FRWATER
    if(associated(AODRNIR)) AODRNIR = (1.-ALBNRO)*DRNIR*FRWATER
    if(associated(AODFNIR)) AODFNIR = (1.-ALBNFO)*DFNIR*FRWATER
    if(associated(FSURF  )) FSURF   = SWN+LWDNSRF-(ALW+BLW*TS(:,WATER))-SHF-LHF

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

    call MAPL_GetPointer(IMPORT,FI     , 'FRACINEW'  ,    RC=STATUS); VERIFY_(STATUS)
    if(associated(FRW)) then
       FRW = max(1.0 - FI, 0.0)
    endif

    call MAPL_TimerOff(MAPL,   "-OpenWater")

    call MAPL_TimerOn(MAPL,    "-Albedo")

    if(solalarmison) then
       call MAPL_SunGetInsolation(LONS, LATS,      &
            ORBIT, ZTH, SLR,                       &
            INTV = TINT,                           &
            currTime=CURRENT_TIME+DELT,            &
            RC=STATUS )
       VERIFY_(STATUS)

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

! !IROUTINE: ALBSEA - Computes albedos as a function of $cos(\zeta)$ over ocean surfaces

! !INTERFACE:

  subroutine ALBSEA (ALBVR,ALBVF,ALBNR,ALBNF,ZTH)

! !ARGUMENTS:

    real,    intent(IN)  :: ZTH  (:)
    real,    intent(OUT) :: ALBVR(:) ! visible beam    albedo
    real,    intent(OUT) :: ALBVF(:) ! visible diffuse albedo
    real,    intent(OUT) :: ALBNR(:) ! nearIr  beam    albedo
    real,    intent(OUT) :: ALBNF(:) ! nearIr  diffuse albedo

!  !DESCRIPTION:
!        Compute albedo for ocean points
!          based on ceres

!  CERES ocean albedo at zth=.5 is 0.052. Our formulation gives .077
!    thus the scaling. The diffuse albedo is given by computing
!    the zth weighted average of the albedo over the hemisphere and
!    then applying the same scaling to match CERES.


!LLT: CERESFAC = 1           reduces to old formulation 1-5-05
!     CERESFAC = 0.052/0.077 is the Original CERES Factor
!     CERESFAC = 0.068/0.077 is the EROS Tuned Value

    real, parameter :: CERESFAC   = 0.068/0.077
!   real, parameter :: CERESFAC   = 1.0

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
  end subroutine ALBSEA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: SIMPLE_SW_ABS - 
!  Implements two simple ways for the absorption of shortwave radiation,
!  as an alternative (for "TESTING" purposes) to using KPAR & KUVR based on PEN,
!  into an interface layer of typical depth = 2m, that is also called "depth of AOIL"

! !INTERFACE:

  subroutine SIMPLE_SW_ABS(NT, USE_KPAR, depth, ZTH, SWN, PEN)

! !ARGUMENTS:

    integer, intent(IN)    :: NT        ! dimension of array
    integer, intent(IN)    :: USE_KPAR  ! absorption profile option
    real,    intent(IN)    :: ZTH(:)    ! cosine of solar zenith angle
    real,    intent(IN)    :: depth(:)  ! depth up to which shortwave needs to be absorbed
    real,    intent(IN)    :: SWN(:)    ! net shortwave at surface of ocean, or @ top of air/sea interface
    real,    intent(OUT)   :: PEN(:)    ! shortwave penetrated below the depth    

!  local variables
    real, dimension(NT)  :: fW

    fW  = 0.0
    PEN = 0.0                ! initialize to zero

    if (USE_KPAR == -1) then
       ! Soloviev, 1982 shortwave absorption profile
       ! --------------------------------------------
       fW = 0.28*exp(-71.5*depth) + 0.27*exp(-2.8*depth) + 0.45*exp(-0.07*depth)

    else if (USE_KPAR == -2) then
       ! Paulson & Simpson, 1981- Taken from Gentemann et al, 2009
       ! ----------------------------------------------------------
       fW    = 0.237*exp(-(depth*ZTH)/34.84)  +  0.36*exp(-(depth*ZTH)/2.266)   + &
               0.179*exp(-(depth*ZTH)/0.0315) + 0.087*exp(-(depth*ZTH)/0.0055)  + &
                0.08*exp(-(depth*ZTH)/8.32e-4)+ 0.025*exp(-(depth*ZTH)/1.26e-4) + &
               0.025*exp(-(depth*ZTH)/3.13e-4)+ 0.007*exp(-(depth*ZTH)/7.82e-4) + &
              0.0004*exp(-(depth*ZTH)/1.44e-5)
    else
       if(MAPL_AM_I_ROOT()) print *, 'ERROR! Unknown use_kpar option: ', USE_KPAR
    end if

    PEN   = SWN * fW

   RETURN_(ESMF_SUCCESS)
  end subroutine SIMPLE_SW_ABS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: COOL_SKIN - Computes variables related to the Cool Skin Layer

! !INTERFACE:

  subroutine COOL_SKIN (NT,CM,UUA,VVA,UW,VW,SWN,LHF,SHF,LWDNSRF,    &
                        ALW,BLW,TXW,TYW,USTARW_,                    &
                        DCOOL_,TDROP_,SWCOOL_,QCOOL_,BCOOL_,LCOOL_, &
                        TS,WATER,FR,n_iter_cool,fr_ice_thresh)

! !ARGUMENTS:

    integer, intent(IN)    :: NT             ! number of tiles
!   real,    intent(IN)    :: FR     (:)     ! fraction of sea ice
    real,    intent(IN)    :: FR     (:,:)   ! fraction of surface (water/ice)
    integer, intent(IN)    :: WATER          ! subtile  number assigned to surface type: "WATER" 
    real,    intent(IN)    :: CM     (:,:)   ! transfer coefficient for wind
    real,    intent(IN)    :: UUA    (:)     ! zonal       wind
    real,    intent(IN)    :: VVA    (:)     ! meridional  wind
    real,    intent(IN)    :: UW     (:)     ! u-current
    real,    intent(IN)    :: VW     (:)     ! v-current
    real,    intent(IN)    :: SWN    (:)     ! net shortwave radiation incident at surface
    real,    intent(IN)    :: LHF    (:)     ! latent   heat flux
    real,    intent(IN)    :: SHF    (:)     ! sensible heat flux
    real,    intent(IN)    :: LWDNSRF(:)     ! downward longwave at surface
    real,    intent(IN)    :: ALW    (:)     ! for linearized \sigma T^4
    real,    intent(IN)    :: BLW    (:)     ! for linearized \sigma T^4
    integer, intent(IN)    :: n_iter_cool    ! number of iterations to compute cool-skin layer 
    real,    intent(IN)    :: fr_ice_thresh  ! threshold on ice fraction, sort of defines Marginal Ice Zone
    real,    intent(IN)    :: TS     (:,:)   ! skin temperature

    real,    intent(OUT)   :: USTARW_(:)     ! u_{*,w} 
    real,    intent(OUT)   :: DCOOL_ (:)     ! depth of cool-skin layer
    real,    intent(OUT)   :: TDROP_ (:)     ! temperature drop across cool-skin
    real,    intent(OUT)   :: SWCOOL_(:)     ! shortwave radiation absorbed in cool-skin 
    real,    intent(OUT)   :: QCOOL_ (:)     ! net heat flux in cool layer
    real,    intent(OUT)   :: BCOOL_ (:)     ! bouyancy in cool layer
    real,    intent(OUT)   :: LCOOL_ (:)     ! Saunder's parameter in cool layer

    real,    intent(INOUT) :: TXW    (:)     ! zonal      stress
    real,    intent(INOUT) :: TYW    (:)     ! meridional stress

!  !LOCAL VARIABLES

    integer         :: N, iter_cool

    real            :: ALPH, Qb, fC

    real, parameter :: NU_WATER        = 1.0E-6  ! kinematic viscosity of water  [m^2/s]
    real, parameter :: TherCond_WATER  = 0.563   ! Thermal conductivity of water [W/m/ K]
    real, parameter :: bigC            = &
          (16.0 * (MAPL_CAPWTR*MAPL_RHO_SEAWATER)**2 * NU_WATER**3) / TherCond_WATER**2

!  !DESCRIPTION:
!        Based on Fairall et al, 1996

    do N = 1, NT  ! N is now looping over all tiles (NOT sub-tiles).

!      Stress over "open" water (or Marginal Ice Zone) depends on ocean currents
!      --------------------------------------------------------------------------
       TXW(N) = CM(N,WATER)*(UUA(N) - UW(N))
       TYW(N) = CM(N,WATER)*(VVA(N) - VW(N))

!      if( FR(N)       < fr_ice_thresh ) then 
       if( FR(N,WATER) > fr_ice_thresh ) then 

!        Ustar in water has a floor of 2 \mu m/s
!        ----------------------------------------
         USTARW_(N) = max( 2.e-6, sqrt(sqrt(TXW(N)*TXW(N)+TYW(N)*TYW(N))/MAPL_RHO_SEAWATER) )

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        ! Cool skin layer- heat loss and temperature drop  @ top of interface layer !
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          DCOOL_(N)  = 1.e-3           ! initial guess for cool-skin layer thickness
          TDROP_(N)  = 0.2             ! guess for cool-skin tdrop. FINAL TDROP IS SENSITIVE TO INITIAL CHOICE. 3 ITER ENOUGH?

          cool_iter: do iter_cool = 1, n_iter_cool

!         Short wave absorbed in the cool layer. This is a modified version of Zeng and Beljaars, 2005
!         ----------------------------------------------------------------------------------------------

             fC  = 0.0685 + 11.0*DCOOL_(N) - (3.3e-5/DCOOL_(N))*(1.0-exp(-DCOOL_(N)/8.0E-4))
             fC  = max( fC, 0.01)        ! absorb at least 1% of shortwave in cool layer
             SWCOOL_(N) = SWN(N)*fC

!            Heat loss at top of skin (cool) layer
!            --------------------------------------

             QCOOL_(N)  = &
                          LHF(N) + SHF(N) - ( LWDNSRF(N) -(ALW(N) + BLW(N)*( TS(N,WATER)-TDROP_(N)))) - &
                          SWCOOL_(N)

!            Bouyancy production in cool layer depends on surface cooling
!            and evap-salinity effect from surface. It does not depend on solar
!            heating, which is assumed to be uniform in cool layer. This last assumption
!            could be improved by including some NIR. For this calculation, we include
!            temperature dependence of the thermal expansion coefficient.
!            -------------------------------------------------------------------------------

             ALPH   = (0.6 + 0.0935*(TS(N,WATER)-MAPL_TICE))*1.E-4
             Qb     = QCOOL_(N) + ( (0.026*MAPL_CAPWTR)/(ALPH*MAPL_ALHL) )*LHF(N)
             BCOOL_(N) = (ALPH*MAPL_GRAV*Qb) / (MAPL_RHO_SEAWATER*MAPL_CAPWTR)

!            Saunders parameter
!            BigC = (16.0 * (MAPL_CAPWTR*MAPL_RHO_SEAWATER)**2 * NU_WATER**3) / TherCond_WATER**2  
!            -------------------------------------------------------------------------------

             if ( BCOOL_(N) > 0.0) then  ! Eqn(14) of F96
                LCOOL_(N)  = 6.0/( 1.0 + ( BCOOL_(N)*bigC / USTARW_(N)**4 )**0.75 )**(1./3.)
                DCOOL_(N)  = LCOOL_(N)*NU_WATER/USTARW_(N)
             else 
                LCOOL_(N)  = 6.0
                DCOOL_(N)  = min( LCOOL_(N)*NU_WATER/USTARW_(N), 1.e-2)    ! Prevent very thick cool layer depth
             end if

             TDROP_(N)    = max( 0.0, DCOOL_(N)*QCOOL_(N)/TherCond_WATER ) ! Eqn(4) & (13) of F96

          end do cool_iter

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        ! Done with Cool skin layer.  Now turbluent heat flux at base of interface layer !
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       else            ! FR(N, WATER) <= fr_ice_thresh
          USTARW_(N)     = MAPL_UNDEF
          DCOOL_ (N)     = MAPL_UNDEF
          TDROP_ (N)     = 0.0 !MAPL_UNDEF
          SWCOOL_(N)     = MAPL_UNDEF
          QCOOL_ (N)     = MAPL_UNDEF
          BCOOL_ (N)     = MAPL_UNDEF
          LCOOL_ (N)     = MAPL_UNDEF
       end if
    end do

   RETURN_(ESMF_SUCCESS)
  end subroutine COOL_SKIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: AOIL_SST - Computes skin SST using AOIL (single column)

! !INTERFACE:

  subroutine AOIL_SST  ( DO_DATASEA, DT, epsilon_d, F_PHI, STOKES_SPEED, &
             AOIL_depth, MUSKIN, SWN, PEN, PEN_ocean, LHF, SHF, LWDNSRF, &
             ALW, BLW, USTARW, TDROP, TS_FOUNDi, SWWARM, QWARM, ZETA_W,  &
             PHIW, LANGM, TAUTW, DELTC, TBAR, TDEL, DTS, TS, TWMTF)

! !ARGUMENTS:

    integer, intent(IN) :: DO_DATASEA  ! =0: coupled with ocean model (CGCM)
                                       ! =1: uncoupled                (AGCM)

    real, intent(IN)    :: DT          ! model time step
    real, intent(IN)    :: epsilon_d   ! ratio: (depth of AOIL)/(ocean model top level)

    real, intent(IN)    :: F_PHI       ! a scaler (tunable parameter) used to compute stability function
    real, intent(IN)    :: STOKES_SPEED! should be input from wave model, dummy for now

    real, intent(IN)    :: AOIL_depth  ! depth of the AOIL
    real, intent(IN)    :: MUSKIN      ! exponent in the prescribed temperature profile within AOIL

    real, intent(IN)    :: SWN         ! net shortwave radiation at surface
    real, intent(IN)    :: PEN         ! net shortwave radiation below the AOIL
    real, intent(IN)    :: PEN_ocean   ! net shortwave radiation below the AOIL with ocean model
    real, intent(IN)    :: LHF         ! latent   heat flux
    real, intent(IN)    :: SHF         ! sensible heat flux
    real, intent(IN)    :: LWDNSRF     ! downward longwave radiation
    real, intent(IN)    :: ALW         ! upward   longwave = ALW + BLW * TS
    real, intent(IN)    :: BLW         ! upward   longwave = ALW + BLW * TS

    real, intent(IN)    :: USTARW      ! friction velocity over water
    real, intent(IN)    :: TDROP       ! temperature drop due to cool skin layer
    real, intent(IN)    :: TS_FOUNDi   ! temperature at base of the AOIL

    real, intent(OUT)   :: SWWARM      ! net shortwave radiation absorbed within the AOIL
    real, intent(OUT)   :: QWARM       ! net heat flux within the AOIL
    real, intent(OUT)   :: ZETA_W      ! similarity parameter = AOIL_depth/(MO length scale)
    real, intent(OUT)   :: PHIW        ! stability function
    real, intent(OUT)   :: LANGM       ! Langmuir `number', dummy for now
    real, intent(OUT)   :: TAUTW       ! time scale at which diurnal warming -> 0, or TWMTF -> 0.
    real, intent(OUT)   :: DELTC       ! temperature drop due to cool skin layer = (above) TDROP
    real, intent(OUT)   :: TBAR        ! depth averaged mean AOIL temperature
    real, intent(OUT)   :: TDEL        ! temperature at top of warm layer (within the AOIL)
    real, intent(OUT)   :: DTS         ! temperature change: next time step - previous

    real, intent(INOUT) :: TS          ! skin SST
    real, intent(INOUT) :: TWMTF       ! AOIL state variable

!  !LOCAL VARIABLES

    real         :: ALPH

!  !DESCRIPTION:
!        Based on Akella and Suarez, 2018
!        "The Atmosphere-Ocean Interface Layer of the NASA
!         Goddard Earth Observing System Model and Data Assimilation System."
!         GMAO Tech Memo, Vol 51.

         ALPH   = (0.6 + 0.0935*(TS-MAPL_TICE))*1.E-4

         SWWARM = SWN - PEN
         QWARM  = SWWARM - (LHF + SHF - (LWDNSRF - ALW - BLW*TS))

         ZETA_W = (AOIL_depth*MAPL_KARMAN/USTARW**3.)*MAPL_GRAV*ALPH*QWARM/(MAPL_RHO_SEAWATER*MAPL_CAPWTR)

         PHIW  = 1.+SQRT(1.+4.*MAPL_KARMAN**2.*(1.+MUSKIN)*F_PHI*AOIL_depth*MAPL_GRAV*ALPH*TWMTF/USTARW**2.)
         PHIW  = 0.5*PHIW

         LANGM = SQRT( USTARW/STOKES_SPEED)

         TAUTW = (AOIL_depth*PHIW)/(MAPL_KARMAN*USTARW*(1.+MUSKIN))

         if (DO_DATASEA == 0) then ! with an OGCM, as in coupled GCM
           QWARM = QWARM - (epsilon_d/(1.-epsilon_d))* (PEN-PEN_ocean)
           TAUTW = (1.- epsilon_d) * TAUTW   ! compare \tau_{\sigma} in Eq.(22) and that in section 2.2 of AS2018
         endif

         TWMTF = (TWMTF+DT*(QWARM/(AOIL_depth*MAPL_RHO_SEAWATER*MAPL_CAPWTR)))/(1.+DT/TAUTW)
         TWMTF = max( TWMTF, 0.)

         DELTC = TDROP

         TBAR  = TWMTF + TS_FOUNDi

         if (DO_DATASEA == 1) then ! Atmospheric GCM, ocean surface is from "data"
          TDEL = TS_FOUNDi + ((1.+MUSKIN)/MUSKIN) * TWMTF
         else                      ! Coupled GCM
          TDEL = TS_FOUNDi + (1./MUSKIN + (1.-epsilon_d)) * TWMTF
         endif

         DTS = (TDEL- TDROP) - TS
         TS  = TDEL - TDROP       ! updated skin temperature

   RETURN_(ESMF_SUCCESS)
  end subroutine AOIL_SST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: SKIN_SST - Computes changes to SST in interface layer due to Cool Skin & Diurnal Warming 

! !INTERFACE:

  subroutine SKIN_SST (DO_SKIN_LAYER,NT,CM,UUA,VVA,UW,VW,HW,SWN,LHF,SHF,LWDNSRF,                   &
                       ALW,BLW,PEN,STOKES_SPEED,DT,MUSKIN,TS_FOUNDi,DWARM_,TBAR_,TXW,TYW,USTARW_,  &
                       DCOOL_,TDROP_,SWCOOL_,QCOOL_,BCOOL_,LCOOL_,TDEL_,SWWARM_,QWARM_,ZETA_W_,    &
                       PHIW_,LANGM_,TAUTW_,uStokes_,TS,TWMTS,TW,WATER,FR,n_iter_cool,fr_ice_thresh)

! !ARGUMENTS:

    integer, intent(IN)    :: DO_SKIN_LAYER  ! 0: No interface layer,     1: active, and accounts for change in SST
    integer, intent(IN)    :: NT             ! number of tiles
    real,    intent(IN)    :: FR     (:,:)   ! fraction of surface (water/ice)
    integer, intent(IN)    :: WATER          ! subtile  number assigned to surface type: "WATER" 
    real,    intent(IN)    :: CM     (:,:)   ! transfer coefficient for wind
    real,    intent(IN)    :: UUA    (:)     ! zonal       wind
    real,    intent(IN)    :: VVA    (:)     ! meridional  wind
    real,    intent(IN)    :: UW     (:)     ! u-current
    real,    intent(IN)    :: VW     (:)     ! v-current
    real,    intent(IN)    :: HW     (:)     ! mass  of skin layer
    real,    intent(IN)    :: SWN    (:)     ! net shortwave radiation incident at surface
    real,    intent(IN)    :: LHF    (:)     ! latent   heat flux
    real,    intent(IN)    :: SHF    (:)     ! sensible heat flux
    real,    intent(IN)    :: LWDNSRF(:)     ! longwave at surface
    real,    intent(IN)    :: ALW    (:)     ! for linearized \sigma T^4
    real,    intent(IN)    :: BLW    (:)     ! for linearized \sigma T^4
    real,    intent(IN)    :: PEN    (:)     ! shortwave radiation that penetrates below interface layer
    real,    intent(IN)    :: STOKES_SPEED   ! scalar value set for Stokes speed- place holder for output from Wave model
    real,    intent(IN)    :: DT             ! time-step
    real,    intent(IN)    :: MUSKIN         ! exponent of temperature: T(z) profile in warm layer
    real,    intent(IN)    :: TS_FOUNDi(:)   ! bulk SST (temperature at base of warm layer)
    integer, intent(IN)    :: n_iter_cool    ! number of iterations to compute cool-skin layer 
    real,    intent(IN)    :: fr_ice_thresh  ! threshold on ice fraction, sort of defines Marginal Ice Zone

    real,    intent(OUT)   :: DWARM_ (:)     ! depth of skin layer
    real,    intent(OUT)   :: TBAR_  (:)     ! copy of TW (also internal state) to export out
    real,    intent(OUT)   :: USTARW_(:)     ! u_{*,w} 
    real,    intent(OUT)   :: DCOOL_ (:)     ! depth of cool-skin layer
    real,    intent(OUT)   :: TDROP_ (:)     ! temperature drop across cool-skin
    real,    intent(OUT)   :: SWCOOL_(:)     ! shortwave radiation absorbed in cool-skin 
    real,    intent(OUT)   :: QCOOL_ (:)     ! net heat flux in cool layer
    real,    intent(OUT)   :: BCOOL_ (:)     ! bouyancy in cool layer
    real,    intent(OUT)   :: LCOOL_ (:)     ! Saunder's parameter in cool layer

    real,    intent(OUT)   :: TDEL_  (:)     ! temperature at top of warm layer
    real,    intent(OUT)   :: SWWARM_(:)     ! shortwave radiation absorbed in warm layer
    real,    intent(OUT)   :: QWARM_ (:)     ! net heat flux in warm layer
    real,    intent(OUT)   :: ZETA_W_(:)     ! stability parameter = dwarm/(Obukhov length)
    real,    intent(OUT)   :: PHIW_  (:)     ! similarity function
    real,    intent(OUT)   :: LANGM_ (:)     ! Langmuir number
    real,    intent(OUT)   :: TAUTW_ (:)     ! time-scale of relaxation to bulk SST (i.e., TS_FOUND)
    real,    intent(OUT)   :: uStokes_(:)    ! Stokes speed

    real,    intent(INOUT) :: TXW    (:)     ! zonal      stress
    real,    intent(INOUT) :: TYW    (:)     ! meridional stress
    real,    intent(INOUT) :: TWMTS  (:)     ! "internal state" variable that has: TW - TS
    real,    intent(INOUT) :: TW     (:)     ! "internal state" variable that has: TW
    real,    intent(INOUT) :: TS     (:,:)   ! skin temperature

!  !LOCAL VARIABLES

    integer         :: N, iter_cool
    real            :: ALPH, Qb, fC, fLA, X1, X2

    real, parameter :: RHO_SEAWATER    = 1022.0  ! sea water density             [kg/m^3]    ! Replace Usage of RHO_SEAWATER with MAPL_RHO_SEAWATER
    real, parameter :: NU_WATER        = 1.0E-6  ! kinematic viscosity of water  [m^2/s]
    real, parameter :: TherCond_WATER  = 0.563   ! Thermal conductivity of water [W/m/ K]
    real, parameter :: bigC            = &
          (16.0 * (MAPL_CAPWTR*MAPL_RHOWTR)**2 * NU_WATER**3) / TherCond_WATER**2

!  !DESCRIPTION:
!        Based on Fairall et al, 1996 for Cool Skin Layer and Takaya et al, 2010 for Warm Layer

! Open water conditions, including computation of skin layer parameters
!----------------------------------------------------------------------

    do N = 1, NT  ! N is now looping over all tiles (NOT sub-tiles).

! Stress over "open" water (or Marginal Ice Zone) depends on ocean currents
!--------------------------------------------------------------------------

       TXW(N) = CM(N,WATER)*(UUA(N) - UW(N))
       TYW(N) = CM(N,WATER)*(VVA(N) - VW(N))

       if( FR(N, WATER) > fr_ice_thresh) then 

! Depth and mean temperature of interface layer
!----------------------------------------------

          DWARM_(N) = HW(N)/MAPL_RHOWTR                                                   ! replace MAPL_RHOWTR with MAPL_RHO_SEAWATER
          TBAR_(N)  = TS(N,WATER) + TWMTS(N)

! Ustar in water has a floor of 2 \mu m/s
!----------------------------------------

          USTARW_(N) = max( 2.e-6, sqrt(sqrt(TXW(N)*TXW(N)+TYW(N)*TYW(N))/MAPL_RHOWTR) )  ! replace MAPL_RHOWTR with MAPL_RHO_SEAWATER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cool skin layer- heat loss and temperature drop  @ top of interface layer !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          DCOOL_(N)  = 1.e-3           ! initial guess for cool-skin layer thickness
          TDROP_(N)  = 0.2             ! guess for cool-skin tdrop. FINAL TDROP IS SENSITIVE TO INITIAL CHOICE. 3 ITER ENOUGH?

          COOL_SKIN: do iter_cool = 1, n_iter_cool

! Short wave absorbed in the cool layer. This is a modified version of Zhang and Beljaars, 2005
!----------------------------------------------------------------------------------------------

             fC  = 0.0685 + 11.0*DCOOL_(N) - (3.3e-5/DCOOL_(N))*(1.0-exp(-DCOOL_(N)/8.0E-4))
             fC  = max( fC, 0.01)        ! absorb at least 1% of shortwave in cool layer
             SWCOOL_(N) = SWN(N)*fC

! Heat loss at top of skin (cool) layer
!--------------------------------------

             X1 = LHF(N) + SHF(N) - ( LWDNSRF(N) -(ALW(N) + BLW(N)*( TS(N,WATER)-TDROP_(N))))
             QCOOL_(N)  = X1 - SWCOOL_(N)

! Bouyancy production in cool layer depends on surface cooling
! and evap-salinity effect from surface. It does not depend on solar
! heating, which is assumed to be uniform in cool layer. This last assumption
! could be improved by including some NIR. For this calculation, we include
! temperature dependence of the thermal expansion coefficient.
!-------------------------------------------------------------------------------

             ALPH   = (0.6 + 0.0935*(TBAR_(N)-MAPL_TICE))*1.E-4
             Qb     = QCOOL_(N) + ( (0.026*MAPL_CAPWTR)/(ALPH*MAPL_ALHL) )*LHF(N)
             BCOOL_(N) = (ALPH*MAPL_GRAV*Qb) / (RHO_SEAWATER*MAPL_CAPWTR)                 ! replace RHO_SEAWATER with MAPL_RHO_SEAWATER

! Saunders parameter
! BigC = (16.0 * (MAPL_CAPWTR*MAPL_RHO_SEAWATER)**2 * NU_WATER**3) / TherCond_WATER**2  
!-------------------------------------------------------------------------------

             if ( BCOOL_(N) > 0.0) then  ! Eqn(14) of F96
                LCOOL_(N)  = 6.0/( 1.0 + ( BCOOL_(N)*bigC / USTARW_(N)**4 )**0.75 )**(1./3.)
                DCOOL_(N)  = LCOOL_(N)*NU_WATER/USTARW_(N)
             else 
                LCOOL_(N)  = 6.0
                DCOOL_(N)  = min( LCOOL_(N)*NU_WATER/USTARW_(N), 1.e-2)  ! Prevent very thick cool layer depth
             end if

             TDROP_(N)    = max( 0.0, DCOOL_(N)*QCOOL_(N)/TherCond_WATER ) ! Eqn(4) & (13) of F96

          end do COOL_SKIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Done with Cool skin layer.  Now turbluent heat flux at base of interface layer  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          WARM_LAYER: if(DO_SKIN_LAYER==0) then   ! Warm layer temperature increase calculated based on definition of mean interface temperature.

             TDEL_(N)    = TS(N,WATER) + TDROP_(N)
             TWMTS(N)    = TBAR_(N)    - TS(N,WATER)

!            fill up with mapl_undef - so that LocStreamMod does NOT die while exporting
             SWWARM_(N)  = MAPL_UNDEF
             QWARM_ (N)  = MAPL_UNDEF
             ZETA_W_(N)  = MAPL_UNDEF
             PHIW_(N)    = MAPL_UNDEF
             LANGM_(N)   = MAPL_UNDEF
             TAUTW_(N)   = MAPL_UNDEF

          else  ! use Takaya et al 2012

! Compute warm layer temperature increase based on Takaya et al, 2012
!--------------------------------------------------------------------

             ALPH        = (0.6 + 0.0935*(TBAR_(N)-MAPL_TICE))*1.E-4

! Short wave absorbed in the warm layer.
!--------------------------------------

             X1         = LHF(N) + SHF(N) - ( LWDNSRF(N) -(ALW(N) + BLW(N)*TS(N,WATER)))
             SWWARM_(N) = SWN(N) - PEN(N)
             QWARM_(N)  = SWWARM_(N) - X1

! Stability parameter & Similarity function
!------------------------------------------

             ZETA_W_(N)  = (DWARM_(N)*MAPL_KARMAN*MAPL_GRAV*ALPH*QWARM_(N)) / &           ! zeta_w = dwarm/obukhov length
                           (RHO_SEAWATER*MAPL_CAPWTR*USTARW_(N)**3)                       ! replace RHO_SEAWATER with MAPL_RHO_SEAWATER

             if ( ZETA_W_(N) >= 0.0) then   ! Takaya: Eqn(5)
                PHIW_(N) = 1. + (5*ZETA_W_(N) + 4.*ZETA_W_(N)**2)/(1+3.*ZETA_W_(N)+0.25*ZETA_W_(N)**2)
             else
                PHIW_(N) = 1.0/sqrt(1.-16.*ZETA_W_(N))
             end if

! Langmuir number- need imports from Wave Model
!----------------------------------------------

             uStokes_(N) = STOKES_SPEED
             LANGM_(N)   = sqrt(USTARW_(N)/uStokes_(N))
             fLA         = LANGM_(N)**(-0.66667)           ! Takaya: Eqn(6)

             IF (fLA       <= 1.0) fLA = 1.0               ! Limit range of fLa to be >=1
             IF (ZETA_W_(N)<= 0.0) fLA = 1.0               ! Apply fLa to stable conditions only

             TAUTW_(N)   = &
                  (DWARM_(N)*PHIW_(N))/(MAPL_KARMAN*USTARW_(N)*fLA*(MUSKIN+1.))

             X2          = DT * &
                  ( (MAPL_KARMAN*USTARW_(N)*fLA*(MUSKIN+1.))/(DWARM_(N)*PHIW_(N)) )

! We DO NOT include cool-skin tdrop in TW, therefore, we now save TW

             TW(N)       = TS_FOUNDi(N) + ( 1.0/(1.+X2))        *    (TBAR_(N) - TS_FOUNDi(N))
             TS(N,WATER) = TS(N,WATER)  + ((1.0+MUSKIN)/MUSKIN) *    (TW(N)    - TBAR_(N))

             TDEL_(N)    = TS_FOUNDi(N) + ((1.0+MUSKIN)/MUSKIN) * MAX(TW(N)    - TS_FOUNDi(N), 0.0)
             TBAR_(N)    = TW(N)

             TS(N,WATER) = TDEL_(N) - TDROP_(N)
             TWMTS(N)    = TW(N) - TS(N,WATER)
          end if WARM_LAYER

       else            ! FR(N, WATER) <= fr_ice_thresh
          DCOOL_ (N)     = MAPL_UNDEF
          LCOOL_ (N)     = MAPL_UNDEF
          DWARM_ (N)     = MAPL_UNDEF
          TBAR_  (N)     = MAPL_UNDEF
          TDROP_ (N)     = MAPL_UNDEF
          QCOOL_ (N)     = MAPL_UNDEF
          USTARW_(N)     = MAPL_UNDEF
          SWCOOL_(N)     = MAPL_UNDEF
          BCOOL_ (N)     = MAPL_UNDEF
          TDEL_  (N)     = MAPL_UNDEF
          TWMTS  (N)     = 0.0
          QWARM_ (N)     = MAPL_UNDEF
          SWWARM_(N)     = MAPL_UNDEF
          PHIW_  (N)     = MAPL_UNDEF
          LANGM_ (N)     = MAPL_UNDEF
          TAUTW_ (N)     = MAPL_UNDEF
          ZETA_W_(N)     = MAPL_UNDEF
       end if
    end do


   RETURN_(ESMF_SUCCESS)
  end subroutine SKIN_SST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine RUN2


end module GEOS_OpenwaterGridCompMod


