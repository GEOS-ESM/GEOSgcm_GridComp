!   $Id$

#include "MAPL_Generic.h"

!=============================================================================

module GEOS_TurbulenceGridCompMod

!BOP

!  !MODULE: GEOS_Turbulence --- An GEOS generic atmospheric turbulence component

! !USES:

  use ESMF
  use GEOS_Mod
  use MAPL
  use LockEntrain
  use shoc

#ifdef _CUDA
  use cudafor
#endif

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !DESCRIPTION:
! 
!   {\tt GEOS\_TurbulenceGridComp} computes atmospheric tendencies due to turbulence.
!   Its physics is a combination of the first-order scheme of Louis---for stable PBLs
!   and free atmospheric turbulence---with a modified version of the non-local-K
!   scheme proposed by Lock for unstable and cloud-topped boundary layers.
!   In addition to diffusive tendencies, it adds the effects orographic form drag
!   for features with horizontal scales of 2 to 20 km following Beljaars et al. (2003,
!   ECMWF Tech. Memo. 427).  
!
!\vspace{12 pt}
!\noindent
!{\bf Grid Considerations}
!
!   Like all GEOS\_Generic-based components, it works on an inherited 
!   3-dimensional ESMF grid. It assumes that the first two (inner) dimensions span the
!   horizontal and the third (outer) dimension is the vertical. In the horizontal,
!   one or both dimensions can be degenerate, effectively supporting 
!   single-columns (1-D), and slices (2-D). No horizontal dimension needs to be
!   aligned with a particular coordinate. In the vertical, the only assumption
!   is that columns are indexed from top to bottom.
!
!\vspace{12 pt}
!\noindent
!{\bf Methods}
!
!   {\tt GEOS\_TurbulenceGridComp} uses the default Initialize and Finalize methods
!   of GEOS\_Generic. It has a 2-stage Run method that can be used in conjunction with
!   two-stage surface calculations to implement semi-implicit time differencing.
!
!\vspace{12 pt}
!\noindent
!{\bf Time Behavior}
!
!   {\tt GEOS\_TurbulenceGridComp} assumes both run stages will be invoked every 
!   RUN\_DT seconds, where RUN\_DT is required in the configuration. On this interval
!   both run stages will perform diffusion updates using diffusivities found in the
!   internal state.  The diffusivities in the internal state may be refreshed intermitently
!   by specifying MY\_STEP and ACCUMINT in the configuration. Accumulated imports used
!   in the intermittent refreshing are valid only on MY\_STEP intervals. Currently the
!   origin of these intervals is the beginning of the run. Accumulation of these imports
!   is done for a period ACCUMINT prior to the valid time. Both ACCUMINT and MY\_STEP are
!   in seconds.
!
!\vspace{12 pt}
!\noindent
!{\bf Working with Bundles and Friendlies}
!
!   {\tt GEOS\_TurbulenceGridComp} works on bundles of quantities to be diffused
!   and with corresponding bundles of their tendencies, surface values, etc.
!   These bundles may contain an arbitrary number of conservative quantities and
!   no requirements or restrictions are placed on what quantities they contain.
!   Quantities required for the calculation, such as pressures, stability, etc
!   are passed separately from the diffused quantities. Little distinction is made
!   of what is in the bundle, except that needed to decide what diffusivity applies
!   to the quantity and in what form its effects are implemented.
!
!   Quantities to be diffused can be marked as "Friendly-for-diffusion". In that case,
!   {\tt GEOS\_TurbulenceGridComp} directly updates the quantity; otherwise it 
!   merely computes its tendency, placing it in the appropriate bundle and treating
!   the quantity itself as read-only.
!
!   In working with bundled quantities, corresponding fields must appear in the 
!   same order in all bundles. Some of these fields, however, 
!   may be ``empty'' in the sense that the data pointer has not been allocated.
!   
!   {\tt GEOS\_TurbulenceGridComp} works with six bundles; three in the import
!   state and three in the export state. The import bundles are:
! \begin{itemize}
!   \item[]
!   \makebox[1in][l]{\bf TR} 
!   \parbox[t]{4in}{The quantity being diffused.}
!   \item[]
!   \makebox[1in][l]{\bf TRG} 
!   \parbox[t]{4in}{The surface (ground) value of the quantity being diffused.
!                   (Used only by Run2)}
!   \item[]
!   \makebox[1in][l]{\bf DTG} 
!   \parbox[t]{4in}{The change of TRG during the time step. (Used only by Run2)}
! \end{itemize}
!
!   The export bundles are:
! \begin{itemize}
!   \item[]
!   \makebox[1in][l]{\bf TRI} 
!   \parbox[t]{4in}{The tendency of the quantity being diffused.
!                   (Produced by Run1, updated by Run2.)  }
!   \item[]
!   \makebox[1in][l]{\bf FSTAR} 
!   \parbox[t]{4in}{After Run1, the ``preliminary'' (i.e., at the original surface
!    value) surface flux of the diffused quantity; after Run2, its final value.
!    (Produced by Run1, updated by Run2)}
!   \item[]
!   \makebox[1in][l]{\bf DFSTAR} 
!   \parbox[t]{4in}{The change of preliminary FSTAR per unit change in the 
!                   surface value. (Produced by Run1)}
! \end{itemize}
!
!   All fields in the export bundles are checked for associated pointers before being
!   updated.
!
!   Fields in the TR bundle can have four attributes:
! \begin{itemize}
! \item FriendlyTo[{\it Component Name}]: default=false --- If true, TR field is updated.
! \item WeightedTendency: default=true --- If true, tendencies (TRI) are pressure-weighted
! \item DiffuseLike: ('S','Q','M') default='S' --- Use mixing coefficients for either
!          heat, moisture or momentum.
! \end{itemize}
!  
!   Only fields in the TR bundle are checked for friendly status. Non-friendly
!   fields in TR and all other bundles are treated with the usual Import/Export
!   rules.
!
!\vspace{12 pt}
!\noindent
!{\bf Other imports and exports}
!
!   In addition to the updates of these bundles, {\tt GEOS\_TurbulenceGridComp} produces
!   a number of diagnostic exports, as well as frictional heating contributions. The latter 
!   are NOT added by {\tt GEOS\_TurbulenceGridComp}, but merely exported to be added
!   elsewhere in the GCM.
!
!\vspace{12 pt}
!\noindent
!{\bf Two-Stage Interactions with the Surface}
!
!   The two-stage scheme for interacting with the surface module is as follows:
! \begin{itemize}
!   \item  The first run stage takes the surface values of the diffused quantities
!      and the surface exchange coefficients as input. These are, of course, on the 
!      grid turbulence is working on.
!   \item  It then does the full diffusion calculation assuming the surface values are
!      fixed, i.e., the explicit surface case. In addition, it also computes derivatives of the
!      tendencies wrt surface values. These are to be used in the second stage.
!   \item The second run stage takes the increments of the surface values as inputs
!      and produces the final results, adding the implicit surface contributions. 
!   \item It also computes the frictional heating due to both implicit and explicit
!       surface contributions.
! \end{itemize}
!
!\vspace{12 pt}
!\noindent
!{\bf GEOS-5 Specific Aspects}
!
!   In GEOS-5, {\tt GEOS\_TurbulenceGridComp} works on the atmosphere's lat-lon grid,
!   while surface quantities are computed during the first run stage of the each of
!   the tiled surface components.  The tiled quantities are properly aggregated to
!   the GEOS-5 lat-lon grid by the first stage of {\tt GEOS\_SurfaceGridComp}, which
!   is called immediately before the first run stage of {\tt GEOS\_TurbulenceGridComp}.
!
!EOP

    logical                             :: dflt_false = .false.
    character(len=ESMF_MAXSTR)          :: dflt_q     = 'Q'
contains

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !DESCRIPTION: This version uses the {\tt GEOS\_GenericSetServices}, which sets
!               the Initialize and Finalize services to generic versions. It also
!   allocates our instance of a generic state and puts it in the 
!   gridded component (GC). Here we only set the two-stage run method and
!   declare the data services.
! \newline
! !REVISION HISTORY: 
!   ??Jul2006 E.Novak./Todling - Added output defining TLM/ADM trajectory

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code
!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Set the Run entry points
! ------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run1, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

!BOS

! !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'PLE',                                       &
        LONG_NAME  = 'air_pressure',                              &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                          &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

      call MAPL_AddImportSpec(GC,                                 &
         SHORT_NAME = 'T',                                        &
         LONG_NAME  = 'air_temperature',                          &
         UNITS      = 'K',                                        &
         DIMS       = MAPL_DimsHorzVert,                          &
         VLOCATION  = MAPL_VLocationCenter,                       &
         RESTART    = MAPL_RestartSkip,                           &
                                                        RC=STATUS  )
      VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'TH',                                        &
        LONG_NAME  = 'potential_temperature',                     &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'QV',                                        &
        LONG_NAME  = 'specific_humidity',                         &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'QLLS',                                      &
        LONG_NAME  = 'liquid_condensate_mixing_ratio',            &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'QILS',                                      &
        LONG_NAME  = 'frozen_condensate_mixing_ratio',            &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'CLLS',                                      &
        LONG_NAME  = 'cloud_fraction',                            &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'QLCN',                                      &
        LONG_NAME  = 'liquid_condensate_mixing_ratio',            &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'QICN',                                      &
        LONG_NAME  = 'frozen_condensate_mixing_ratio',            &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'CLCN',                                      &
        LONG_NAME  = 'cloud_fraction',                            &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'U',                                         &
        LONG_NAME  = 'eastward_wind',                             &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME = 'V',                                         &
        LONG_NAME  = 'northward_wind',                            &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'CT',                                &
        LONG_NAME          = 'surface_heat_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'CQ',                                &
        LONG_NAME          = 'surface_moisture_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'CM',                                &
        LONG_NAME          = 'surface_momentum_exchange_coefficient', &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'BSTAR',                             &
        LONG_NAME          = 'surface_bouyancy_scale',            &
        UNITS              = 'm s-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'USTAR',                             &
        LONG_NAME          = 'surface_velocity_scale',            &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'FRLAND',                            &
        LONG_NAME          = 'land_fraction',                     &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'RADLW',                             &
        LONG_NAME          = 'air_temperature_tendency_due_to_longwave',&
        UNITS              = 'K s-1',                             &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'RADLWC',                            &
        LONG_NAME          = 'clearsky_air_temperature_tendency_lw',&
        UNITS              = 'K s-1',                             &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'PREF',                              &
        LONG_NAME          = 'reference_air_pressure',            &
        UNITS              = 'Pa',                                &
        DIMS               = MAPL_DimsVertOnly,                   &
        VLOCATION          = MAPL_VLocationEdge,                  &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'VARFLT',                            &
        LONG_NAME          = 'variance_of_filtered_topography',   &
        UNITS              = 'm+2',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TR',                                &
        LONG_NAME          = 'diffused_quantities',               &
        UNITS              = 'X',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        DATATYPE           = MAPL_BundleItem,                     &
        RESTART    = MAPL_RestartSkip,                            &

                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TRG',                               &
        LONG_NAME          = 'surface_values_of_diffused_quantity',&
        UNITS              = 'X',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DATATYPE           = MAPL_BundleItem,                     &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'DTG',                               &
        LONG_NAME          = 'change_of_surface_values_of_diffused_quantity',&
        UNITS              = 'X',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        DATATYPE           = MAPL_BundleItem,                     &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'vertical_pressure_velocity',                                        &
       UNITS      = 'Pa s-1',                                                 &
       SHORT_NAME = 'OMEGA',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'EVAP',                            &
        LONG_NAME          = 'surface_evaporation',   &
        UNITS              = 'kg m-2 s-1',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'SH',                            &
        LONG_NAME          = 'surface_sensible_heat_flux',   &
        UNITS              = 'W m-2',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)



! !EXPORT STATE:

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME = 'TRI',                                       &
        LONG_NAME  = 'diffusion_tendencies',                      &
        UNITS      = 'X kg m-2 s-1',                              &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        DATATYPE   = MAPL_BundleItem,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME = 'FSTAR',                                     &
        LONG_NAME  = 'surface_fluxes',                            &
        UNITS      = 'X kg m-2 s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
        DATATYPE   = MAPL_BundleItem,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME = 'DFSTAR',                                    &
        LONG_NAME  = 'change_of_surface_fluxes_for_unit_change_of_surface_value',&
        UNITS      = 'kg m-2 s-1',                                &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
        DATATYPE   = MAPL_BundleItem,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'air_temperature',                                       &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'T',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'eastward_wind',                                         &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'U',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'northward_wind',                                        &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'V',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'specific_humidity',                                     &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'QV',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'total_momentum_diffusivity',                            &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KM',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'total_scalar_diffusivity',                              &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KH',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'Richardson_number_from_Louis',                          &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'RI',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'bulk_shear_from_Louis',                                 &
       UNITS      = 's-1',                                                   &
       SHORT_NAME = 'DU',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'scalar_diffusivity_from_Louis',                         &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KHLS',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'momentum_diffusivity_from_Louis',                       &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KMLS',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_driven_scalar_diffusivity_from_Lock_scheme',    &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KHSFC',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'radiation_driven_scalar_diffusivity_from_Lock_scheme',  &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KHRAD',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'cloudy_LW_radiation_tendency_used_by_Lock_scheme',      &
       UNITS      = 'K s-1',                                                 &
       SHORT_NAME = 'LWCRT',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'entrainment_heat_diffusivity_from_Lock',                &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'EKH',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'entrainment_momentum_diffusivity_from_Lock',            &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'EKM',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'Blackadar_length_scale_for_scalars',                    &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'ALH',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'p-weighted_frictional_heating_rate_from_diffusion',     &
       UNITS      = 'K s-1 Pa',                                              &
       SHORT_NAME = 'INTDIS',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'p-weighted_frictional_heating_rate_from_orographic_drag',&
       UNITS      = 'K s-1 Pa',                                              &
       SHORT_NAME = 'TOPDIS',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME='DPDTTRB',                                      & 
         LONG_NAME ='layer_pressure_thickness_tendency_from_turbulence', &
         UNITS     ='Pa s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                             &
         VLOCATION = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'p-weighted_frictional_heating_rate_from_surface_drag',  &
       UNITS      = 'K s-1 Pa',                                              &
       SHORT_NAME = 'SRFDIS',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'HGTLM5',                                    &
         LONG_NAME  = 'height_at_LM5',&
         UNITS      = 'm',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'LM50M',                                    &
         LONG_NAME  = 'LM_at_50_meters',&
         UNITS      = '1',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'KETRB',                                    &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_across_turbulence',&
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'KESRF',                                    &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_dissipation_due_to_surface_friction',&
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'KEINT',                                    &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_dissipation_due_to_diffusion',&
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'KETOP',                                    &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_dissipation_due_to_topographic_friction',&
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'entrainment_velocity_from_surface_plume',               &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'WESFC',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'entrainment_velocity_from_radiation',                   &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'WERAD',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'entrainment_velocity_from_buoy_rev',                    &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'WEBRV',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'Buoyancy_jump_across_inversion',                        &
       UNITS      = 'm s-2',                                                 &
       SHORT_NAME = 'DBUOY',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'turbulent_velocity_scale_for_sfc',                      &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'VSCSFC',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'turbulent_velocity_scale_for_cooling',                  &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'VSCRAD',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'turbulent_velocity_scale_for_buoy_rev',                 &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'VSCBRV',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'turbulent_entrainment_diff_from_cooling',               &
       UNITS      = 'm+2 s-1',                                               &
       SHORT_NAME = 'KERAD',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'cloud_top_radiative_forcing',                           &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'CLDRF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'pbltop_pressure',                                       &
       UNITS      = 'Pa',                                                    &
       SHORT_NAME = 'PPBL',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'pbltop_height_for_sfc_plume_LOCK',                      &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'ZSML',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'depth_for_rad/brv_plume_LOCK',                          &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'ZRADML',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'hght_of_base_for_rad/brv_plume_LOCK',                   &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'ZRADBS',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'pbltop_cloud_depth_LOCK',                               &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'ZCLD',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'pbltop_cloud_top_height_LOCK',                          &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'ZCLDTOP',                                               &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'optimal_mixture_fraction_for_BRV',                      &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'CHIS',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 's_of_optimal_mixture_for_BRV',                          &
       UNITS      = 'J kg-1',                                                &
       SHORT_NAME = 'SMIXT',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'Scaled_Del_s_at_Cloud_top',                             &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'DELSINV',                                               &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'Siems_buoy_rev_parameter',                              &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DSIEMS',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'Return_codes_for_Lock_top_driven_plume',                &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'RADRCODE',                                              &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'matrix_diagonal_ak_for_scalars_over_dt',                &
       SHORT_NAME = 'AKSODT',                                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'matrix_diagonal_ck_for_scalars_over_dt',                &
       SHORT_NAME = 'CKSODT',                                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'matrix_diagonal_ak_for_moisture_over_dt',               &
       SHORT_NAME = 'AKQODT',                                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'matrix_diagonal_ck_for_moisture_over_dt',               &
       SHORT_NAME = 'CKQODT',                                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'matrix_diagonal_ak_for_winds_over_dt',                  &
       SHORT_NAME = 'AKVODT',                                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'matrix_diagonal_ck_for_winds_over_dt',                  &
       SHORT_NAME = 'CKVODT',                                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'transcom_planetary_boundary_layer_height',              &
       SHORT_NAME = 'TCZPBL',                                                &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'planetary_boundary_layer_height_threshold_2',           &
       SHORT_NAME = 'ZPBL2',                                                 &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'planetary_boundary_layer_height_threshold_10p',         &
       SHORT_NAME = 'ZPBL10p',                                               &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'planetary_boundary_layer_height_horiz_tke',             &
       SHORT_NAME = 'ZPBLHTKE',                                              &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'turbulent_kinetic_energy',                              &
       SHORT_NAME = 'TKE',                                                   &
       UNITS      = 'm+2 s-2',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'planetary_boundary_layer_height_rich_0',                &
       SHORT_NAME = 'ZPBLRI',                                                &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'planetary_boundary_layer_height_rich_02',               &
       SHORT_NAME = 'ZPBLRI2',                                               &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'planetary_boundary_layer_height_thetav',                &
       SHORT_NAME = 'ZPBLTHV',                                               &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'pbltop_level',                                          &
       SHORT_NAME = 'KPBL',                                                  &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'pbltop_level_for_shallow',                              &
       SHORT_NAME = 'KPBL_SC',                                               &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'dry_static_energy_after_diffuse',                       &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'SAFDIFFUSE',                                            &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'dry_static_energy_after_update',                        &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'SAFUPDATE',                                             &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)



! !INTERNAL STATE:

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_ahat_for_scalars',                      &
       SHORT_NAME = 'AKS',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_bhat_for_scalars',                      &
       SHORT_NAME = 'BKS',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_c_for_scalars',                         &
       SHORT_NAME = 'CKS',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'sensitivity_of_tendency_to_surface_value_for_scalars',  &
       SHORT_NAME = 'DKS',                                                   &
       UNITS      = 's-1',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_ahat_for_moisture',                     &
       SHORT_NAME = 'AKQ',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_bhat_for_moisture',                     &
       SHORT_NAME = 'BKQ',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_c_for_moisture',                        &
       SHORT_NAME = 'CKQ',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'sensitivity_of_tendency_to_surface_value_for_moisture', &
       SHORT_NAME = 'DKQ',                                                   &
       UNITS      = 's-1',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_ahat_for_winds',                        &
       SHORT_NAME = 'AKV',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_bhat_for_winds',                        &
       SHORT_NAME = 'BKV',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'matrix_diagonal_c_for_winds',                           &
       SHORT_NAME = 'CKV',                                                   &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'sensitivity_of_tendency_to_surface_value_for_winds',    &
       SHORT_NAME = 'DKV',                                                   &
       UNITS      = 's-1',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'momentum_mixing_factor',                                &
       SHORT_NAME = 'EKV',                                                   &
       UNITS      = 'Pa s-1',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'topographic_roughness_factor',                          &
       SHORT_NAME = 'FKV',                                                   &
       UNITS      = 'Pa s-1',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'geopotential_height_above_surface',                     &
       SHORT_NAME = 'ZLE',                                                   &
       UNITS      = 'm',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'turbulence_tendency_for_dry_static_energy',             &
       SHORT_NAME = 'SINC',                                                  &
       UNITS      = 'm+2 s-3',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       RESTART    = MAPL_RestartSkip,                            &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'ZPBL',                                       &
       LONG_NAME  = 'planetary_boundary_layer_height',            &
       UNITS      = 'm',                                          &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzOnly,                           &
       VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'TKESHOC',                                   &
       LONG_NAME  = 'turbulent_kinetic_energy_from_SHOC',        &
       UNITS      = 'm+2 s-2',                                   &
       DEFAULT    = 0.0,                                           &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'TKEDISS',                                   &
       LONG_NAME  = 'tke_dissipation_from_SHOC',        &
       UNITS      = 'm+2 s-3',                                   &
       DEFAULT    = 0.0,                                           &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'TKEBUOY',                                   &
       LONG_NAME  = 'tke_buoyancy_production_from_SHOC',        &
       UNITS      = 'm+2 s-3',                                   &
       DEFAULT    = 0.0,                                           &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'TKESHEAR',                                   &
       LONG_NAME  = 'tke_shear_production_from_SHOC',        &
       UNITS      = 'm+2 s-3',                                   &
       DEFAULT    = 0.0,                                           &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'TKETRANS',                                  &
       LONG_NAME  = 'tke_transport_from_SHOC',                   &
       UNITS      = 'm+2 s-3',                                   &
       DEFAULT    = 0.0,                                         &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'TKH',                                       &
       LONG_NAME  = 'turbulent_heat_diffusivity_from_SHOC',      &
       UNITS      = 'm+2 s-1',                                   &
       DEFAULT    = 0.0,                                          &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'LSHOC',                                     &
       LONG_NAME  = 'eddy_dissipation_length_from_SHOC',        &
       UNITS      = 'm',                                         &
       DEFAULT    = 0.0,                                         &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'BRUNTSHOC',                                 &
       LONG_NAME  = 'Brunt_Vaisala_frequency_from_SHOC',         &
       UNITS      = 's-1',                                       &
       DEFAULT    = 0.0,                                         &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'SHEARSHOC',                                 &
       LONG_NAME  = 'Shear_from_SHOC',                           &
       UNITS      = 's-1',                                       &
       DEFAULT    = 0.0,                                         &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'CLD',                                       &
       LONG_NAME  = 'cloud_fraction_diagnosed_by_SHOC',          &
       UNITS      = '1',                                         &
       DEFAULT    = 0.0,                                           &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'WTHV2',                               &
       LONG_NAME  = 'Buoyancy_flux_diagnosed_by_SHOC',           &
       UNITS      = '1',                                         &
       DEFAULT    = 0.0,                                           &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'DTDT_SHC',                               &
       LONG_NAME  = 'Temperature_tendency_diagnosed_by_SHOC',           &
       UNITS      = 'K s-1',                                         &
       DEFAULT    = 0.0,                                           &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'DQDT_SHC',                               &
       LONG_NAME  = 'Water_vapor_tendency_diagnosed_by_SHOC',           &
       UNITS      = 'kg kg-1 s-1',                                         &
       DEFAULT    = 0.0,                                           &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'DQLDT_SHC',                               &
       LONG_NAME  = 'Liquid_water_tendency_diagnosed_by_SHOC',           &
       UNITS      = 'kg kg-1 s-1',                                         &
       DEFAULT    = 0.0,                                           &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                &
       SHORT_NAME = 'DQIDT_SHC',                               &
       LONG_NAME  = 'Ice_water_tendency_diagnosed_by_SHOC',           &
       UNITS      = 'kg kg-1 s-1',                                         &
       DEFAULT    = 0.0,                                           &
       FRIENDLYTO = trim(COMP_NAME),                             &
       DIMS       = MAPL_DimsHorzVert,                           &
       VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)
!EOS

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,   name="-RUN1"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="--DIFFUSE"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="--REFRESHKS" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---PRELIMS"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---SHOC"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---LOUIS"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---LOCK"     ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="----LOCK_RUN",RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="----LOCK_DATA",RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="----LOCK_ALLOC",RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="----LOCK_DEALLOC",RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---POSTLOCK" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---BELJAARS" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="---DECOMP"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="-RUN2"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="--UPDATE"    ,RC=STATUS)
    VERIFY_(STATUS)
    
! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================


!BOP

! !IROUTINE: RUN1 -- First run stage for the {\tt MAPL_TurbulenceGridComp} component

! !INTERFACE:

  subroutine RUN1 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

! !DESCRIPTION: The first run stage of {\tt GEOS\_TurbulenceGridComp} computes the diffusivities,
!   sets-up the matrix for a backward-implicit computation of the surface fluxes,
!   and solves this system for a fixed surface value of the diffused quantity. Run1
!   takes as inputs the surface exchange coefficients (i.e., $\rho |U| C_{m,h,q}$) for
!   momentun, heat, and moisture, as well as the pressure, temperature, moisture, 
!   and winds for the sounding. These are used only for computing the diffusivities
!   and, as explained above,  are not the temperatures, moistures, etc. being diffused.
!
!   The computation of turbulence fluxes for fixed surface values is done at every
!   time step in the contained subroutine {\tt DIFFUSE}; but the computation of 
!   diffusivities and orographic drag coefficients, as well as the set-up of the
!   vertical difference matrix and its LU decomposition
!   can be done intermittently for economy in the contained subroutine  {\tt REFRESH}.
!   The results of this calculation are stored in an internal state. 
!   Run1 also computes the sensitivity of the 
!   atmospheric tendencies and the surface flux to changes in the surface value.
!
!   The diffusivities are computed by calls to {\tt LOUIS\_KS} and {\tt ENTRAIN}, which
!   compute the Louis et al. (1983) and Lock (2000) diffusivities. The Louis 
!   diffusivities are computed for all conditions, and {\tt ENTRAIN} overrides them 
!   where appropriate. Lock can be turned off from the resource file.


!

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp), pointer   :: MAPL
    type (ESMF_Config      )            :: CF
    type (ESMF_State       )            :: INTERNAL 
    type (ESMF_Alarm       )            :: ALARM   

! Local variables

    real, dimension(:,:,:), pointer     :: AKS, BKS, CKS, DKS
    real, dimension(:,:,:), pointer     :: AKQ, BKQ, CKQ, DKQ
    real, dimension(:,:,:), pointer     :: AKV, BKV, CKV, DKV, EKV, FKV
    real, dimension(:,:,:), pointer     :: PLE, ZLE, SINC
    real, dimension(:,:  ), pointer     :: CU, CT, CQ, ZPBL
    integer                             :: IM, JM, LM
    real                                :: DT
 
! SHOC-related variables
    real, dimension(:,:,:), pointer     :: TKESHOC,TKH,LSHOC,BRUNTSHOC, & 
                                           SHEARSHOC,CLD, WTHV2,        &
                                           TKEBUOY,TKESHEAR,TKEDISS,TKETRANS,&
                                           DTDT_SHC,DQDT_SHC,DQLDT_SHC, &
                                           DQIDT_SHC
    real, dimension(:,:), pointer       :: EVAP, SH

! Begin... 
!---------

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'Run1'

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"-RUN1")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,        &
         IM=IM, JM=JM, LM=LM,               &
         RUNALARM=ALARM,                    &
         INTERNAL_ESMF_STATE=INTERNAL,      &
                                  RC=STATUS )
    VERIFY_(STATUS)

! Get configuration from component
!---------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

! Get all pointers that are needed by both REFRESH and DIFFUSE
!-------------------------------------------------------------

! Get pressure structure; this is instantaneous.
!-----------------------------------------------

     call MAPL_GetPointer(IMPORT,  PLE,   'PLE',     RC=STATUS)
     VERIFY_(STATUS)

! Get surface exchange coefficients
!----------------------------------

     call MAPL_GetPointer(IMPORT,  CU,     'CM',     RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  CT,     'CT',     RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  CQ,     'CQ',     RC=STATUS)
     VERIFY_(STATUS)

!----- SHOC-related variables -----
    call MAPL_GetPointer(IMPORT, SH,   'SH',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, EVAP, 'EVAP',    RC=STATUS)
    VERIFY_(STATUS)

! Get pointers from internal state
!---------------------------------

    call MAPL_GetPointer(INTERNAL, AKS,   'AKS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, BKS,   'BKS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CKS,   'CKS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DKS,   'DKS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, AKQ,   'AKQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, BKQ,   'BKQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CKQ,   'CKQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DKQ,   'DKQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, AKV,   'AKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, BKV,   'BKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CKV,   'CKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DKV,   'DKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, EKV,   'EKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, FKV,   'FKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, ZLE,   'ZLE',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, SINC,  'SINC',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, ZPBL,  'ZPBL',    RC=STATUS)
    VERIFY_(STATUS)

!----- SHOC-related variables -----
    call MAPL_GetPointer(INTERNAL, TKESHOC,'TKESHOC',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, TKEDISS,'TKEDISS',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, TKEBUOY,'TKEBUOY',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, TKESHEAR,'TKESHEAR',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, TKETRANS,'TKETRANS',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, TKH,   'TKH',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, LSHOC,   'LSHOC',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, BRUNTSHOC,'BRUNTSHOC',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, SHEARSHOC,'SHEARSHOC',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLD,   'CLD',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, WTHV2, 'WTHV2',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DTDT_SHC, 'DTDT_SHC',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DQDT_SHC, 'DQDT_SHC',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DQLDT_SHC, 'DQLDT_SHC',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DQIDT_SHC, 'DQIDT_SHC',    RC=STATUS)
    VERIFY_(STATUS)


! Get application's timestep from configuration
!----------------------------------------------

    call ESMF_ConfigGetAttribute(CF, DT, Label="RUN_DT:" , RC=STATUS)
    VERIFY_(STATUS)

! If its time, do the refresh
! ---------------------------

    if ( ESMF_AlarmIsRinging(ALARM, rc=status) ) then
       VERIFY_(STATUS)
       call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
       VERIFY_(STATUS)

       call MAPL_TimerOn (MAPL,"--REFRESHKS")
        call REFRESH(IM,JM,LM,RC=STATUS)
        VERIFY_(STATUS)
       call MAPL_TimerOff(MAPL,"--REFRESHKS")
    endif

! Solve the free atmosphere problem
! ---------------------------------

    call MAPL_TimerOn (MAPL,"--DIFFUSE")
     call DIFFUSE(IM,JM,LM,RC=STATUS)
     VERIFY_(STATUS)
    call MAPL_TimerOff(MAPL,"--DIFFUSE")

!  All done with RUN1
!--------------------

    call MAPL_TimerOff(MAPL,"-RUN1")
    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  contains

!=============================================================================
!=============================================================================

!BOP

! !CROUTINE: REFRESH -- Refreshes diffusivities.

! !INTERFACE:

   subroutine REFRESH(IM,JM,LM,RC)

! !ARGUMENTS:

     integer,           intent(IN)       :: IM,JM,LM
     integer, optional, intent(OUT)      :: RC

! !DESCRIPTION: 
!   {\tt REFRESH} can be called intermittently to compute new values of the 
!   diffusivities. In addition it does all possible calculations that depend
!   only on these. In particular, it sets up the semi-implicit tridiagonal
!   solver in the vertical and does the LU decomposition. It also includes the
!   local effects of orographic drag, so that it to is done implicitly.
!
!   Diffusivities are first computed with the Louis scheme ({\tt LOUIS\_KS}),
!   and then, where appropriate,
!   they are overridden by the Lock values ({\tt ENTRAIN}).
!   Once diffusivities are computed, {\tt REFRESH} sets-up the tridiagonal
!   matrices for the semi-implicit vertical diffusion calculation and performs
!   their $LU$ decomposition. 
!
!   {\tt REFRESH} requires surface exchange coefficients for heat, moisture, and
!   momentum,  The calculations in the interior are also
!   done for momentum, heat, and water diffusion. Heat and water mixing
!   coefficients differ only at the surface, but these affect the entire $LU$
!   decomposition, and so all three decompositions are saved in the internal state. 
!
!   For a conservatively diffused quantity $q$, we have
!   $$
!   \frac{\partial q}{\partial t} = -g \frac{\partial }{\partial p} 
!       \left(\rho K_q \frac{\partial q}{\partial z} \right)
!   $$
!   In finite difference form, using backward time differencing, this becomes
!   $$
!   \begin{array}{rcl}
!   {q^{n+1}_l - q^{n}_l} & = & - \frac{g}{\delta_l p}^*
!     \delta_l \left[
!      \left( \frac{\Delta t \rho K_q}{\delta_l z} \right)^* (\delta_l q)^{n+1} \right]   \\
!   &&\\
!                         & = & - \alpha_l ( \beta_{l+\frac{1}{2}}(q_{l+1}-q_l)^{n+1} - 
!                                            \beta_{l-\frac{1}{2}}(q_l-q_{l-1})^{n+1} ) \\
!   &&\\
!   \alpha_l & = & \frac{g \Delta t}{(p_{l+\frac{1}{2}}-p_{l-\frac{1}{2}})^*} \\
!   &&\\
!   \beta_{l+\frac{1}{2}} & = & \left( \frac{ (\rho K_q)^*_{l+\frac{1}{2}}}{(z_{l+1}-z_{l})^*} \right) \\
!   \end{array}
!   $$
!   where the subscripts denote levels, superscripts denote times, and the $*$ superscript
!   denotes evaluation at the refresh time.
!   The following tridiagonal set is then solved for $q^{n+1}_l$:
!   $$
!   a_l q_{l-1} + b_l q_l + c_l q_{l+1} = q_l
!   $$
!   where
!   $$
!   \begin{array}{rcl}
!   a_l & = & \alpha_l \beta_{l-\frac{1}{2}} \\
!   c_l & = & \alpha_l \beta_{l+\frac{1}{2}} \\
!   b_l & = & 1 - a_l - c_l.
!   \end{array}
!   $$
!   At the top boundary, we assume $K_q=0$, so  $ \beta_{\frac{1}{2}}=0$ and $a_1=0$.
!   At the surface, $ \beta_{L+\frac{1}{2}}= \rho_s |U|_s C_{m,h,q}$, the surface exchange coefficient.
!   

!EOP
 
     character(len=ESMF_MAXSTR)          :: IAm='Refresh'
     integer                             :: STATUS

     real, dimension(:,:,:), pointer     :: TH, U, V, OMEGA, Q, T, RI, DU, RADLW, RADLWC, LWCRT
     real, dimension(:,:  ), pointer     :: VARFLT
     real, dimension(:,:,:), pointer     :: KH, KM, QLLS, QILS, CLLS, QLCN, QICN, CLCN
     real, dimension(:,:,:), pointer     :: ALH
     real, dimension(:    ), pointer     :: PREF

     real, dimension(IM,JM,1:LM-1)       :: TVE, RDZ
     real, dimension(IM,JM,LM)           :: THV, TV, Z, DMI, PLO, QL, QI, QA
     real, dimension(IM,JM,0:LM)         :: PKE

     real, dimension(:,:,:), pointer     :: EKH, EKM, KHLS, KMLS, KHRAD, KHSFC
     real, dimension(:,:  ), pointer     :: BSTAR, USTAR, PPBL, WERAD, WESFC,VSCRAD,KERAD,DBUOY,ZSML,ZCLD,ZRADML,FRLAND
     real, dimension(:,:  ), pointer     :: TCZPBL => null()
     real, dimension(:,:  ), pointer     :: ZPBL2 => null()
     real, dimension(:,:  ), pointer     :: ZPBL10P => null()
     real, dimension(:,:  ), pointer     :: ZPBLHTKE => null()
     real, dimension(:,:,:), pointer     :: TKE => null()
     real, dimension(:,:  ), pointer     :: ZPBLRI => null()
     real, dimension(:,:  ), pointer     :: ZPBLRI2 => null()
     real, dimension(:,:  ), pointer     :: ZPBLTHV => null()
     real, dimension(:,:  ), pointer     :: KPBL => null()
     real, dimension(:,:  ), pointer     :: KPBL_SC => null()
     real, dimension(:,:  ), pointer     :: WEBRV,VSCBRV,DSIEMS,CHIS,ZCLDTOP,DELSINV,SMIXT,ZRADBS,CLDRF,VSCSFC,RADRCODE

     real, dimension(:,:,:), pointer     :: AKSODT, CKSODT
     real, dimension(:,:,:), pointer     :: AKQODT, CKQODT
     real, dimension(:,:,:), pointer     :: AKVODT, CKVODT


     logical, dimension(IM,JM     )      :: CONVECT
     logical                             :: ALLOC_TCZPBL, CALC_TCZPBL
     logical                             :: ALLOC_ZPBL2, CALC_ZPBL2
     logical                             :: ALLOC_ZPBL10p, CALC_ZPBL10p

     real                                :: LOUIS, ALHFAC, ALMFAC
     real                                :: LAMBDAM, LAMBDAM2
     real                                :: LAMBDAH, LAMBDAH2
     real                                :: ZKMENV, ZKHENV 
     real                                :: MINTHICK
     real                                :: MINSHEAR
     real                                :: AKHMMAX
     real                                :: C_B, LAMBDA_B, ZMAX_B,LOUIS_MEMORY
     real                                :: PRANDTLSFC,PRANDTLRAD,BETA_RAD,BETA_SURF,KHRADFAC,TPFAC_SURF,ENTRATE_SURF
     real                                :: PCEFF_SURF, KHSFCFAC_LND, KHSFCFAC_OCN, ZCHOKE

     integer                             :: I,J,L,LOCK_ON
     integer                             :: KPBLMIN,PBLHT_OPTION

     real,               dimension(LM+1) :: temparray, htke
     real,               dimension(IM,JM,LM  ) :: tcrib !TransCom bulk Ri
     real,               dimension(LM+1) :: thetav

     ! variables associated with SHOC
!     real, dimension( IM, JM )           :: HFLX,EVAP
     real, dimension( IM, JM, LM )       :: PRANDTLSHOC,QPL,QPI,TDUM,QDUM, &
                                            QIDUM,QLDUM
     integer                             :: DO_SHOC,SHC_DO_COND,SHC_DO_TRANS
     real                                :: SHC_LAMBDA,SHC_TSCALE,SHC_VONK,SHC_CK, &
                                            SHC_CEFAC,SHC_CESFAC,SHC_THL2TUNE,    &
                                            SHC_QW2TUNE,SHC_QWTHL2TUNE
     real    :: LAMBDADISS

     integer :: locmax
     real    :: maxkh,minlval
     real, dimension(IM,JM) :: thetavs,thetavh,uv2h,kpbltc,kpbl2,kpbl10p
     real    :: maxdthvdz,dthvdz

     ! PBL-top diagnostic
     ! -----------------------------------------

     real, parameter :: tcri_crit = 0.25
     real, parameter :: ri_crit = 0.00
     real, parameter :: ri_crit2 = 0.20

     real(kind=MAPL_R8), dimension(IM,JM,LM) :: AKX, BKX

#ifdef _CUDA
     type(dim3) :: Grid, Block
     integer :: blocksize_x, blocksize_y
#endif

! Get Sounding from the import state
!-----------------------------------

     call MAPL_GetPointer(IMPORT,     T,       'T', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,     Q,      'QV', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,    TH,      'TH', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,     U,       'U', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,     V,       'V', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT, OMEGA,   'OMEGA', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,VARFLT,  'VARFLT', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  PREF,    'PREF', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT, RADLW,   'RADLW', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,RADLWC,  'RADLWC', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  QLLS,    'QLLS', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  QILS,    'QILS', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  QLCN,    'QLCN', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  QICN,    'QICN', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  CLLS,    'CLLS', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,  CLCN,    'CLCN', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT, BSTAR,   'BSTAR', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT, USTAR,   'USTAR', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT,FRLAND,  'FRLAND', RC=STATUS); VERIFY_(STATUS)

! Get turbulence parameters from configuration
!---------------------------------------------

     if (LM .eq. 72) then
       call MAPL_GetResource (MAPL, LOUIS,        trim(COMP_NAME)//"_LOUIS:",        default=5.0,          RC=STATUS)
       call MAPL_GetResource (MAPL, ALHFAC,       trim(COMP_NAME)//"_ALHFAC:",       default=1.2,          RC=STATUS)
       call MAPL_GetResource (MAPL, ALMFAC,       trim(COMP_NAME)//"_ALMFAC:",       default=1.2,          RC=STATUS)
     else
       call MAPL_GetResource (MAPL, LOUIS,        trim(COMP_NAME)//"_LOUIS:",        default=1.0,          RC=STATUS)
       call MAPL_GetResource (MAPL, ALHFAC,       trim(COMP_NAME)//"_ALHFAC:",       default=1.0,          RC=STATUS)
       call MAPL_GetResource (MAPL, ALMFAC,       trim(COMP_NAME)//"_ALMFAC:",       default=1.0,          RC=STATUS)
     endif
     call MAPL_GetResource (MAPL, LAMBDAM,      trim(COMP_NAME)//"_LAMBDAM:",      default=160.0,        RC=STATUS)
     call MAPL_GetResource (MAPL, LAMBDAM2,     trim(COMP_NAME)//"_LAMBDAM2:",     default=1.0,          RC=STATUS)
     call MAPL_GetResource (MAPL, LAMBDAH,      trim(COMP_NAME)//"_LAMBDAH:",      default=160.0,        RC=STATUS)
     call MAPL_GetResource (MAPL, LAMBDAH2,     trim(COMP_NAME)//"_LAMBDAH2:",     default=1.0,          RC=STATUS)
     call MAPL_GetResource (MAPL, LAMBDADISS,   trim(COMP_NAME)//"_LAMBDADISS:",   default=50.0,         RC=STATUS)
     call MAPL_GetResource (MAPL, ZKMENV,       trim(COMP_NAME)//"_ZKMENV:",       default=3000.,        RC=STATUS)
     call MAPL_GetResource (MAPL, ZKHENV,       trim(COMP_NAME)//"_ZKHENV:",       default=3000.,        RC=STATUS)
     call MAPL_GetResource (MAPL, MINTHICK,     trim(COMP_NAME)//"_MINTHICK:",     default=0.1,          RC=STATUS)
     call MAPL_GetResource (MAPL, MINSHEAR,     trim(COMP_NAME)//"_MINSHEAR:",     default=0.0030,       RC=STATUS)
!    call MAPL_GetResource (MAPL, C_B,          trim(COMP_NAME)//"_C_B:",          default=2.5101471e-8, RC=STATUS)  ! Pre Ganymed-4_1 value
!    call MAPL_GetResource (MAPL, C_B,          trim(COMP_NAME)//"_C_B:",          default=1.02e-7,      RC=STATUS)  ! Value used for Ganymed-4_1 through Heracles-5_4_p3
     call MAPL_GetResource (MAPL, C_B,          trim(COMP_NAME)//"_C_B:",          default=6.00e-7,      RC=STATUS)  ! Value used with updated GMTED TOPO Data
     call MAPL_GetResource (MAPL, LAMBDA_B,     trim(COMP_NAME)//"_LAMBDA_B:",     default=1500.,        RC=STATUS)
     call MAPL_GetResource (MAPL, AKHMMAX,      trim(COMP_NAME)//"_AKHMMAX:",      default=500.,         RC=STATUS)
     call MAPL_GetResource (MAPL, LOCK_ON,      trim(COMP_NAME)//"_LOCK_ON:",      default=1,            RC=STATUS)
     call MAPL_GetResource (MAPL, PRANDTLSFC,   trim(COMP_NAME)//"_PRANDTLSFC:",   default=1.0,          RC=STATUS)
     call MAPL_GetResource (MAPL, PRANDTLRAD,   trim(COMP_NAME)//"_PRANDTLRAD:",   default=0.75,         RC=STATUS)
     call MAPL_GetResource (MAPL, BETA_RAD,     trim(COMP_NAME)//"_BETA_RAD:",     default=0.20,         RC=STATUS)
     call MAPL_GetResource (MAPL, BETA_SURF,    trim(COMP_NAME)//"_BETA_SURF:",    default=0.25,         RC=STATUS)
     call MAPL_GetResource (MAPL, KHRADFAC,     trim(COMP_NAME)//"_KHRADFAC:",     default=0.85,         RC=STATUS)
     call MAPL_GetResource (MAPL, KHSFCFAC_LND, trim(COMP_NAME)//"_KHSFCFAC_LND:", default=0.60,         RC=STATUS)
     call MAPL_GetResource (MAPL, KHSFCFAC_OCN, trim(COMP_NAME)//"_KHSFCFAC_OCN:", default=0.30,         RC=STATUS)
     call MAPL_GetResource (MAPL, TPFAC_SURF,   trim(COMP_NAME)//"_TPFAC_SURF:",   default=20.0,         RC=STATUS)
     call MAPL_GetResource (MAPL, ENTRATE_SURF, trim(COMP_NAME)//"_ENTRATE_SURF:", default=1.5e-3,       RC=STATUS)
     call MAPL_GetResource (MAPL, PCEFF_SURF,   trim(COMP_NAME)//"_PCEFF_SURF:",   default=0.5,          RC=STATUS)
     call MAPL_GetResource (MAPL, LOUIS_MEMORY, trim(COMP_NAME)//"_LOUIS_MEMORY:", default=-999.,        RC=STATUS)

     if (LM .eq. 72) then
     call MAPL_GetResource (MAPL, PBLHT_OPTION, trim(COMP_NAME)//"_PBLHT_OPTION:", default=4,            RC=STATUS)
     else
     call MAPL_GetResource (MAPL, PBLHT_OPTION, trim(COMP_NAME)//"_PBLHT_OPTION:", default=3,            RC=STATUS)
     endif

     call MAPL_GetResource (MAPL, DO_SHOC,      trim(COMP_NAME)//"_DO_SHOC:",      default=0,            RC=STATUS)
     if (DO_SHOC /= 0) then
       call MAPL_GetResource (MAPL, SHC_LAMBDA,   trim(COMP_NAME)//"_SHC_LAMBDA:",   default=0.04,       RC=STATUS)
       call MAPL_GetResource (MAPL, SHC_TSCALE,   trim(COMP_NAME)//"_SHC_TSCALE:",   default=400.,       RC=STATUS)
       call MAPL_GetResource (MAPL, SHC_VONK,     trim(COMP_NAME)//"_SHC_VONK:",     default=0.4,        RC=STATUS)
       call MAPL_GetResource (MAPL, SHC_CK,       trim(COMP_NAME)//"_SHC_CK:",       default=0.1,        RC=STATUS)
       call MAPL_GetResource (MAPL, SHC_CEFAC,    trim(COMP_NAME)//"_SHC_CEFAC:",    default=1.,         RC=STATUS)
       call MAPL_GetResource (MAPL, SHC_CESFAC,   trim(COMP_NAME)//"_SHC_CESFAC:",   default=4.29,        RC=STATUS)
       call MAPL_GetResource (MAPL, SHC_THL2TUNE, trim(COMP_NAME)//"_SHC_THL2TUNE:", default=1.,         RC=STATUS)
       call MAPL_GetResource (MAPL, SHC_QW2TUNE,  trim(COMP_NAME)//"_SHC_QW2TUNE:",  default=1.,         RC=STATUS)
       call MAPL_GetResource (MAPL, SHC_QWTHL2TUNE, trim(COMP_NAME)//"_SHC_QWTHL2TUNE:", default=1.,     RC=STATUS)
       call MAPL_GetResource (MAPL, SHC_DO_COND, trim(COMP_NAME)//"_SHC_DO_COND:", default=0,     RC=STATUS)
       call MAPL_GetResource (MAPL, SHC_DO_TRANS, trim(COMP_NAME)//"_SHC_DO_TRANS:", default=1,     RC=STATUS)
     end if

! Get pointers from export state...
!-----------------------------------

     call MAPL_GetPointer(EXPORT,      KH,      'KH', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,      KM,      'KM', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,      RI,      'RI', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,      DU,      'DU', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,     EKH,     'EKH', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,     EKM,     'EKM', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    KHLS,    'KHLS', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    KMLS,    'KMLS',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   KHSFC,   'KHSFC', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   KHRAD,   'KHRAD', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    PPBL,    'PPBL',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    KPBL,    'KPBL',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    KPBL_SC, 'KPBL_SC',            RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    TCZPBL,  'TCZPBL',             RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZPBL2,  'ZPBL2',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZPBL10p,  'ZPBL10p',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZPBLHTKE,  'ZPBLHTKE',         RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    TKE,  'TKE',         RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZPBLRI,  'ZPBLRI',             RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZPBLRI2,  'ZPBLRI2',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZPBLTHV,  'ZPBLTHV',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   LWCRT,   'LWCRT', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   WERAD,   'WERAD',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   WESFC,   'WESFC',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   DBUOY,   'DBUOY',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  VSCRAD,  'VSCRAD',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  VSCsfc,  'VSCSFC',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   KERAD,   'KERAD',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  VSCBRV,  'VSCBRV',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   WEBRV,   'WEBRV',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    CHIS,    'CHIS',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  DSIEMS,  'DSIEMS',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZCLD,    'ZCLD', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,    ZSML,    'ZSML', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  ZRADML,  'ZRADML', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  ZRADBS,  'ZRADBS', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT, ZCLDTOP, 'ZCLDTOP',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT, DELSINV, 'DELSINV',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,RADRCODE,'RADRCODE',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   SMIXT,   'SMIXT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,   CLDRF,   'CLDRF',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,     ALH,     'ALH',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  AKSODT,  'AKSODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  CKSODT,  'CKSODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  AKQODT,  'AKQODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  CKQODT,  'CKQODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  AKVODT,  'AKVODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT,  CKVODT,  'CKVODT',               RC=STATUS)
     VERIFY_(STATUS)

! Initialize some arrays

      LWCRT = RADLW - RADLWC

      KH    = 0.0
      KM    = 0.0
      RI    = 0.0
      DU    = 0.0
      EKH   = 0.0
      EKM   = 0.0
      KHSFC = 0.0
      KHRAD = 0.0
      if(associated( ALH))  ALH = 0.0
      if(associated(KHLS)) KHLS = 0.0
      if(associated(KMLS)) KMLS = 0.0

      ALLOC_ZPBL2 = .FALSE.
      CALC_ZPBL2 = .FALSE.
      if(associated(ZPBL2).OR.PBLHT_OPTION==1) CALC_ZPBL2 = .TRUE.
      if(.not.associated(ZPBL2 )) then
         allocate(ZPBL2(IM,JM))
         ALLOC_ZPBL2 = .TRUE.
      endif

      ALLOC_ZPBL10p = .FALSE.
      CALC_ZPBL10p = .FALSE.
      if(associated(ZPBL10p).OR.PBLHT_OPTION==2.OR.PBLHT_OPTION==4) CALC_ZPBL10p = .TRUE.
      if(.not.associated(ZPBL10p )) then
         allocate(ZPBL10p(IM,JM))
         ALLOC_ZPBL10p = .TRUE.
      endif

      if(.not.associated(KPBL_SC )) then
         allocate(KPBL_SC(IM,JM))
      end if

      ALLOC_TCZPBL = .FALSE.
      CALC_TCZPBL = .FALSE.
      if(associated(TCZPBL).OR.PBLHT_OPTION==3.OR.PBLHT_OPTION==4) CALC_TCZPBL = .TRUE.
      if(.not.associated(TCZPBL)) then
                allocate(TCZPBL(IM,JM))
                   ALLOC_TCZPBL = .TRUE.
      endif

      call MAPL_TimerOn(MAPL,"---PRELIMS")

      ! Compute the edge heights using Arakawa-Suarez hydrostatic equation
      !---------------------------------------------------------------------------

      PKE = (PLE/MAPL_P00)**MAPL_KAPPA
      ZLE(:,:,LM) = 0.0
      do L = LM, 1, -1
         ZLE(:,:,L-1) = ZLE(:,:,L) + (MAPL_CP/MAPL_GRAV)*TH(:,:,L)*(PKE(:,:,L)-PKE(:,:,L-1))
      end do

      ! Layer height, pressure, and virtual temperatures
      !-------------------------------------------------

      QL  = QLCN + QLLS
      QI  = QICN + QILS
      QA  = CLCN + CLLS ! Currently not used in REFRESHKS
      Z   = 0.5*(ZLE(:,:,0:LM-1)+ZLE(:,:,1:LM))
      PLO = 0.5*(PLE(:,:,0:LM-1)+PLE(:,:,1:LM))

      TV  = T *( 1.0 + MAPL_VIREPS * Q - QL - QI ) 
      THV = TV*(TH/T)

      TVE = (TV(:,:,1:LM-1) + TV(:,:,2:LM))*0.5

      ! Miscellaneous factors
      !----------------------

      RDZ = PLE(:,:,1:LM-1) / ( MAPL_RGAS * TVE )
      RDZ = RDZ(:,:,1:LM-1) / (Z(:,:,1:LM-1)-Z(:,:,2:LM))
      DMI = (MAPL_GRAV*DT)/(PLE(:,:,1:LM)-PLE(:,:,0:LM-1))

      !===> Running 1-2-1 smooth of bottom 5 levels of Virtual Pot. Temp.
      if (LM .eq. 72) then
      THV(:,:,LM  ) = THV(:,:,LM-1)*0.25 + THV(:,:,LM  )*0.75
      THV(:,:,LM-1) = THV(:,:,LM-2)*0.25 + THV(:,:,LM-1)*0.50 + THV(:,:,LM  )*0.25 
      THV(:,:,LM-2) = THV(:,:,LM-3)*0.25 + THV(:,:,LM-2)*0.50 + THV(:,:,LM-1)*0.25 
      THV(:,:,LM-3) = THV(:,:,LM-4)*0.25 + THV(:,:,LM-3)*0.50 + THV(:,:,LM-2)*0.25 
      THV(:,:,LM-4) = THV(:,:,LM-5)*0.25 + THV(:,:,LM-4)*0.50 + THV(:,:,LM-3)*0.25 
      THV(:,:,LM-5) = THV(:,:,LM-6)*0.25 + THV(:,:,LM-5)*0.50 + THV(:,:,LM-4)*0.25 
      endif

      call MAPL_TimerOff(MAPL,"---PRELIMS")


!!!=================================================================
!!!===========================  SHOC  ==============================
!!!=================================================================
!  Description
!
!
!
!!!=================================================================

      if (DO_SHOC /= 0) then

        LOCK_ON = 0

        call MAPL_TimerOn (MAPL,name="---SHOC" ,RC=STATUS)
        VERIFY_(STATUS)

        ! for now just use fixed values
        QPI = 0.
        QPL = 0.
        PRANDTLSHOC = 0.9

        TDUM = T
        QDUM = Q
        QIDUM = QI
        QLDUM = QL

!        print *,'QT before=',Q+QI+QL
!        print *,'HL before=',MAPL_CP*T+MAPL_GRAV*Z-MAPL_ALHL*QL-MAPL_ALHS*QI

        call RUN_SHOC( IM, JM, LM, LM+1, DT, &
                   ! Inputs
                   PLO(:,:,1:LM),         &
                   ZLE(:,:,0:LM),         &
                   Z(:,:,1:LM),           &
                   U(:,:,1:LM),           &
                   V(:,:,1:LM),           &
                   OMEGA(:,:,1:LM),       &
                   SH(:,:),               &
                   EVAP(:,:),             &
                   ! Input-Outputs
                   TDUM(:,:,1:LM),        &
                   QDUM(:,:,1:LM),        &
                   QIDUM(:,:,1:LM),       &
                   QLDUM(:,:,1:LM),       &
                   QPI(:,:,1:LM),         &
                   QPL(:,:,1:LM),         &
                   CLD(:,:,1:LM),         &
                   TKESHOC(:,:,1:LM),     &
                   PRANDTLSHOC(:,:,1:LM), &
                   TKH(:,:,1:LM),         &
                   ! Diagnostics
                   WTHV2(:,:,1:LM),       &
                   TKEDISS(:,:,1:LM),     &
                   TKEBUOY(:,:,1:LM),     &
                   TKESHEAR(:,:,1:LM),    &
                   TKETRANS(:,:,1:LM),    &
                   LSHOC(:,:,1:LM),       &
                   BRUNTSHOC(:,:,1:LM),   &
                   SHEARSHOC(:,:,1:LM),   &
                   ! Tuning params
                   SHC_LAMBDA,            &
                   SHC_TSCALE,            &
                   SHC_VONK,              &
                   SHC_CK,                &
                   SHC_CEFAC,             &
                   SHC_CESFAC,            &
                   SHC_THL2TUNE,          &
                   SHC_QW2TUNE,           &
                   SHC_QWTHL2TUNE,        &
                   SHC_DO_TRANS )

        TKH = max(0.,TKH)

        DTDT_SHC = (TDUM-T)/DT
        DQDT_SHC = (QDUM-Q)/DT
        DQLDT_SHC = (QLDUM-QL)/DT
        DQIDT_SHC = (QIDUM-QI)/DT

        if (SHC_DO_COND /= 0) then
           print *,'SHOC: Applying condensatation tendencies'
           T(:,:,:)  = TDUM
           Q(:,:,:)  = QDUM
           QL(:,:,:) = QLDUM
           QI(:,:,:) = QIDUM
        end if

        KH(:,:,1:LM) = TKH(:,:,1:LM)
        KM(:,:,1:LM) = TKH(:,:,1:LM)*PRANDTLSHOC(:,:,1:LM)

!        print *,'QT after=',Q+QI+QL
!        print *,'HL after=',MAPL_CP*T+MAPL_GRAV*Z-MAPL_ALHL*QL-MAPL_ALHS*QI

        call MAPL_TimerOff (MAPL,name="---SHOC" ,RC=STATUS)
        VERIFY_(STATUS)


      end if

!   Refresh diffusivities: First compute Louis...
!   ---------------------------------------------

      call MAPL_TimerOn (MAPL,name="---LOUIS" ,RC=STATUS)
      VERIFY_(STATUS)

      if (DO_SHOC == 0) then
        call LOUIS_KS(                      &
            Z,ZLE(:,:,1:LM-1),THV,U,V,ZPBL, &    
            KH(:,:,1:LM-1),KM(:,:,1:LM-1),  &
            RI(:,:,1:LM-1),DU(:,:,1:LM-1),  &    
            LOUIS, MINSHEAR, MINTHICK,      &
            LAMBDAM, LAMBDAM2,              & 
            LAMBDAH, LAMBDAH2,              & 
            ALHFAC, ALMFAC,                 &
            ZKMENV, ZKHENV, AKHMMAX,        &
            ALH, KMLS, KHLS                 )
      end if

      call MAPL_TimerOff(MAPL,name="---LOUIS" ,RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (MAPL,name="---LOCK" ,RC=STATUS)
      VERIFY_(STATUS)

      !   ...then add Lock.
      !--------------------

      DO_ENTRAIN: if (LOCK_ON==1) then

#ifdef _CUDA

         _ASSERT(LM <= GPU_MAXLEVS,'needs informative message') !If this is tripped, GNUmakefile
                                    !must be changed

         call MAPL_GetResource(MAPL,BLOCKSIZE_X,'BLOCKSIZE_X:',DEFAULT=16,__RC__)
         call MAPL_GetResource(MAPL,BLOCKSIZE_Y,'BLOCKSIZE_Y:',DEFAULT=8,__RC__)

         Block = dim3(blocksize_x,blocksize_y,1)
         Grid = dim3(ceiling(real(IM)/real(blocksize_x)),ceiling(real(JM)/real(blocksize_y)),1)

         call MAPL_TimerOn (MAPL,name="----LOCK_ALLOC" ,__RC__)

         ! ----------------------
         ! Allocate device arrays
         ! ----------------------

         ! Inputs - Lock
         ! -------------
      
         ALLOCATE(TDTLW_IN_dev(IM,JM,LM), __STAT__)
         ALLOCATE(U_STAR_dev(IM,JM), __STAT__)
         ALLOCATE(B_STAR_dev(IM,JM), __STAT__)
         ALLOCATE(FRLAND_dev(IM,JM), __STAT__)
         ALLOCATE(T_dev(IM,JM,LM), __STAT__)
         ALLOCATE(QV_dev(IM,JM,LM), __STAT__)
         ALLOCATE(QLLS_dev(IM,JM,LM), __STAT__)
         ALLOCATE(QILS_dev(IM,JM,LM), __STAT__)
         ALLOCATE(U_dev(IM,JM,LM), __STAT__)
         ALLOCATE(V_dev(IM,JM,LM), __STAT__)
         ALLOCATE(ZFULL_dev(IM,JM,LM), __STAT__)
         ALLOCATE(PFULL_dev(IM,JM,LM), __STAT__)
         ALLOCATE(ZHALF_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(PHALF_dev(IM,JM,LM+1), __STAT__)

         ! Inoutputs - Lock
         ! ----------------
      
         ALLOCATE(DIFF_M_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(DIFF_T_dev(IM,JM,LM+1), __STAT__)

         ! Outputs - Lock
         ! --------------

         ALLOCATE(K_M_ENTR_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(K_T_ENTR_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(K_SFC_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(K_RAD_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(ZCLOUD_dev(IM,JM), __STAT__)
         ALLOCATE(ZRADML_dev(IM,JM), __STAT__)
         ALLOCATE(ZRADBASE_dev(IM,JM), __STAT__)
         ALLOCATE(ZSML_dev(IM,JM), __STAT__)

         ! Diagnostics - Lock
         ! ------------------

         ! MAT: Using device pointers on CUDA is a bit convoluted. First, we
         ! only allocate the actual working arrays on the device if the 
         ! EXPORT pointer is associated.

         IF (ASSOCIATED(ZCLDTOP))  ALLOCATE(ZCLDTOP_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(WESFC))    ALLOCATE(WENTR_SFC_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(WERAD))    ALLOCATE(WENTR_RAD_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(DBUOY))    ALLOCATE(DEL_BUOY_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(VSCSFC))   ALLOCATE(VSFC_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(VSCRAD))   ALLOCATE(VRAD_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(KERAD))    ALLOCATE(KENTRAD_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(VSCBRV))   ALLOCATE(VBRV_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(WEBRV))    ALLOCATE(WENTR_BRV_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(DSIEMS))   ALLOCATE(DSIEMS_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(CHIS))     ALLOCATE(CHIS_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(DELSINV))  ALLOCATE(DELSINV_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(SMIXT))    ALLOCATE(SLMIXTURE_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(CLDRF))    ALLOCATE(CLDRADF_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(RADRCODE)) ALLOCATE(RADRCODE_DIAG_dev(IM,JM), __STAT__)

         ! Then we associate the CUDA device pointer to the associated device
         ! array. That way CUDA knows what memory that pointer belongs to.
         ! We then pass in the pointer to the subroutine.

         IF (ASSOCIATED(ZCLDTOP))  ZCLDTOP_DIAG_dev_ptr   => ZCLDTOP_DIAG_dev
         IF (ASSOCIATED(WESFC))    WENTR_SFC_DIAG_dev_ptr => WENTR_SFC_DIAG_dev
         IF (ASSOCIATED(WERAD))    WENTR_RAD_DIAG_dev_ptr => WENTR_RAD_DIAG_dev
         IF (ASSOCIATED(DBUOY))    DEL_BUOY_DIAG_dev_ptr  => DEL_BUOY_DIAG_dev
         IF (ASSOCIATED(VSCSFC))   VSFC_DIAG_dev_ptr      => VSFC_DIAG_dev
         IF (ASSOCIATED(VSCRAD))   VRAD_DIAG_dev_ptr      => VRAD_DIAG_dev
         IF (ASSOCIATED(KERAD))    KENTRAD_DIAG_dev_ptr   => KENTRAD_DIAG_dev
         IF (ASSOCIATED(VSCBRV))   VBRV_DIAG_dev_ptr      => VBRV_DIAG_dev
         IF (ASSOCIATED(WEBRV))    WENTR_BRV_DIAG_dev_ptr => WENTR_BRV_DIAG_dev
         IF (ASSOCIATED(DSIEMS))   DSIEMS_DIAG_dev_ptr    => DSIEMS_DIAG_dev
         IF (ASSOCIATED(CHIS))     CHIS_DIAG_dev_ptr      => CHIS_DIAG_dev
         IF (ASSOCIATED(DELSINV))  DELSINV_DIAG_dev_ptr   => DELSINV_DIAG_dev
         IF (ASSOCIATED(SMIXT))    SLMIXTURE_DIAG_dev_ptr => SLMIXTURE_DIAG_dev
         IF (ASSOCIATED(CLDRF))    CLDRADF_DIAG_dev_ptr   => CLDRADF_DIAG_dev
         IF (ASSOCIATED(RADRCODE)) RADRCODE_DIAG_dev_ptr  => RADRCODE_DIAG_dev

         call MAPL_TimerOff(MAPL,name="----LOCK_ALLOC" ,__RC__)

         call MAPL_TimerOn (MAPL,name="----LOCK_DATA" ,__RC__)

         ! ---------------------
         ! Copy inputs to device
         ! ---------------------

         ! Inputs
         ! ------

                  TDTLW_IN_dev = RADLW
                  U_STAR_dev   = USTAR
                  B_STAR_dev   = BSTAR
                  FRLAND_dev   = FRLAND
                  T_dev        = T
                  QV_dev       = Q
                  QLLS_dev     = QLLS
                  QILS_dev     = QILS
                  U_dev        = U
                  V_dev        = V
                  ZFULL_dev    = Z
                  PFULL_dev    = PLO
         ZHALF_dev(:,:,1:LM+1) = ZLE(:,:,0:LM)
         PHALF_dev(:,:,1:LM+1) = PLE(:,:,0:LM)

         ! Inoutputs - Lock
         ! ----------------
      
         DIFF_M_dev(:,:,1:LM+1) = KM(:,:,0:LM)
         DIFF_T_dev(:,:,1:LM+1) = KH(:,:,0:LM)

         call MAPL_TimerOff(MAPL,name="----LOCK_DATA" ,__RC__)

         call MAPL_TimerOn (MAPL,name="----LOCK_RUN" ,__RC__)

         call ENTRAIN<<<Grid,Block>>>(IM, JM, LM,   &
                                      ! Inputs
                                      TDTLW_IN_dev,   &
                                      U_STAR_dev,     &
                                      B_STAR_dev,     &
                                      FRLAND_dev,     &
                                      T_dev,          &
                                      QV_dev,         &
                                      QLLS_dev,       &
                                      QILS_dev,       &
                                      U_dev,          &
                                      V_dev,          &
                                      ZFULL_dev,      &
                                      PFULL_dev,      &
                                      ZHALF_dev,      &
                                      PHALF_dev,      &
                                      ! Inoutputs
                                      DIFF_M_dev,     &
                                      DIFF_T_dev,     &
                                      ! Outputs
                                      K_M_ENTR_dev,   &
                                      K_T_ENTR_dev,   &
                                      K_SFC_dev,      &
                                      K_RAD_dev,      &
                                      ZCLOUD_dev,     &
                                      ZRADML_dev,     &
                                      ZRADBASE_dev,   &
                                      ZSML_dev,       &
                                      ! Diagnostics
                                      ZCLDTOP_DIAG_dev_ptr, &
                                      WENTR_SFC_DIAG_dev_ptr, &
                                      WENTR_RAD_DIAG_dev_ptr, &
                                      DEL_BUOY_DIAG_dev_ptr, &
                                      VSFC_DIAG_dev_ptr, &
                                      VRAD_DIAG_dev_ptr, &
                                      KENTRAD_DIAG_dev_ptr, &
                                      VBRV_DIAG_dev_ptr, &
                                      WENTR_BRV_DIAG_dev_ptr, &
                                      DSIEMS_DIAG_dev_ptr, &
                                      CHIS_DIAG_dev_ptr, &
                                      DELSINV_DIAG_dev_ptr, &
                                      SLMIXTURE_DIAG_dev_ptr, &
                                      CLDRADF_DIAG_dev_ptr, &
                                      RADRCODE_DIAG_dev_ptr, &
                                      ! Input parameter constants
                                      PRANDTLSFC, PRANDTLRAD,   &
                                      BETA_SURF, BETA_RAD,      &
                                      TPFAC_SURF, ENTRATE_SURF, &
                                      PCEFF_SURF, KHRADFAC, KHSFCFAC_LND, KHSFCFAC_OCN )


         STATUS = cudaGetLastError()
         if (STATUS /= 0) then 
            write (*,*) "Error code from ENTRAIN kernel call: ", STATUS
            write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
            _ASSERT(.FALSE.,'needs informative message')
         end if

         ! --------------
         ! Kernel is done
         ! --------------

         call MAPL_TimerOff(MAPL,name="----LOCK_RUN" ,__RC__)

         call MAPL_TimerOn (MAPL,name="----LOCK_DATA" ,__RC__)

         ! ------------------------
         ! Copy outputs to the host
         ! ------------------------

         ! Inoutputs - Lock
         ! ----------------
      
            KM(:,:,0:LM) = DIFF_M_dev(:,:,1:LM+1)
            KH(:,:,0:LM) = DIFF_T_dev(:,:,1:LM+1)

         ! Outputs - Lock
         ! --------------
      
           EKM(:,:,0:LM) = K_M_ENTR_dev(:,:,1:LM+1)
           EKH(:,:,0:LM) = K_T_ENTR_dev(:,:,1:LM+1)
         KHSFC(:,:,0:LM) = K_SFC_dev(:,:,1:LM+1)
         KHRAD(:,:,0:LM) = K_RAD_dev(:,:,1:LM+1)
                  ZCLD   = ZCLOUD_dev
                  ZRADML = ZRADML_dev
                  ZRADBS = ZRADBASE_dev
                  ZSML   = ZSML_dev
      
         ! Diagnostics - Lock
         ! ------------------
      
         IF (ASSOCIATED(ZCLDTOP))  ZCLDTOP  = ZCLDTOP_DIAG_dev
         IF (ASSOCIATED(WESFC))    WESFC    = WENTR_SFC_DIAG_dev
         IF (ASSOCIATED(WERAD))    WERAD    = WENTR_RAD_DIAG_dev
         IF (ASSOCIATED(DBUOY))    DBUOY    = DEL_BUOY_DIAG_dev
         IF (ASSOCIATED(VSCSFC))   VSCSFC   = VSFC_DIAG_dev
         IF (ASSOCIATED(VSCRAD))   VSCRAD   = VRAD_DIAG_dev
         IF (ASSOCIATED(KERAD))    KERAD    = KENTRAD_DIAG_dev
         IF (ASSOCIATED(VSCBRV))   VSCBRV   = VBRV_DIAG_dev
         IF (ASSOCIATED(WEBRV))    WEBRV    = WENTR_BRV_DIAG_dev
         IF (ASSOCIATED(DSIEMS))   DSIEMS   = DSIEMS_DIAG_dev
         IF (ASSOCIATED(CHIS))     CHIS     = CHIS_DIAG_dev
         IF (ASSOCIATED(DELSINV))  DELSINV  = DELSINV_DIAG_dev
         IF (ASSOCIATED(SMIXT))    SMIXT    = SLMIXTURE_DIAG_dev
         IF (ASSOCIATED(CLDRF))    CLDRF    = CLDRADF_DIAG_dev
         IF (ASSOCIATED(RADRCODE)) RADRCODE = RADRCODE_DIAG_dev

         call MAPL_TimerOff(MAPL,name="----LOCK_DATA" ,__RC__)

         call MAPL_TimerOn (MAPL,name="----LOCK_DEALLOC" ,__RC__)

         ! ------------------------
         ! Deallocate device arrays
         ! ------------------------
   
         ! Inputs - Lock
         ! -------------
   
         DEALLOCATE(TDTLW_IN_dev)
         DEALLOCATE(U_STAR_dev)
         DEALLOCATE(B_STAR_dev)
         DEALLOCATE(FRLAND_dev)
         DEALLOCATE(T_dev)
         DEALLOCATE(QV_dev)
         DEALLOCATE(QLLS_dev)
         DEALLOCATE(QILS_dev)
         DEALLOCATE(U_dev)
         DEALLOCATE(V_dev)
         DEALLOCATE(ZFULL_dev)
         DEALLOCATE(PFULL_dev)
         DEALLOCATE(ZHALF_dev)
         DEALLOCATE(PHALF_dev)
      
         ! Inoutputs - Lock
         ! ----------------
   
         DEALLOCATE(DIFF_M_dev)
         DEALLOCATE(DIFF_T_dev)
   
         ! Outputs - Lock
         ! --------------
      
         DEALLOCATE(K_M_ENTR_dev)
         DEALLOCATE(K_T_ENTR_dev)
         DEALLOCATE(K_SFC_dev)
         DEALLOCATE(K_RAD_dev)
         DEALLOCATE(ZCLOUD_dev)
         DEALLOCATE(ZRADML_dev)
         DEALLOCATE(ZRADBASE_dev)
         DEALLOCATE(ZSML_dev)
      
         ! Diagnostics - Lock
         ! ------------------

         ! MAT Again, we only deallocate a device array if the diagnostic
         ! was asked for.
      
         IF (ASSOCIATED(ZCLDTOP))  DEALLOCATE(ZCLDTOP_DIAG_dev)
         IF (ASSOCIATED(WESFC))    DEALLOCATE(WENTR_SFC_DIAG_dev)
         IF (ASSOCIATED(WERAD))    DEALLOCATE(WENTR_RAD_DIAG_dev)
         IF (ASSOCIATED(DBUOY))    DEALLOCATE(DEL_BUOY_DIAG_dev)
         IF (ASSOCIATED(VSCSFC))   DEALLOCATE(VSFC_DIAG_dev)
         IF (ASSOCIATED(VSCRAD))   DEALLOCATE(VRAD_DIAG_dev)
         IF (ASSOCIATED(KERAD))    DEALLOCATE(KENTRAD_DIAG_dev)
         IF (ASSOCIATED(VSCBRV))   DEALLOCATE(VBRV_DIAG_dev)
         IF (ASSOCIATED(WEBRV))    DEALLOCATE(WENTR_BRV_DIAG_dev)
         IF (ASSOCIATED(DSIEMS))   DEALLOCATE(DSIEMS_DIAG_dev)
         IF (ASSOCIATED(CHIS))     DEALLOCATE(CHIS_DIAG_dev)
         IF (ASSOCIATED(DELSINV))  DEALLOCATE(DELSINV_DIAG_dev)
         IF (ASSOCIATED(SMIXT))    DEALLOCATE(SLMIXTURE_DIAG_dev)
         IF (ASSOCIATED(CLDRF))    DEALLOCATE(CLDRADF_DIAG_dev)
         IF (ASSOCIATED(RADRCODE)) DEALLOCATE(RADRCODE_DIAG_dev)

         ! This step is probably unnecessary, but better safe than sorry
         ! as the lifetime of a device pointer is not really specified
         ! by NVIDIA
   
         IF (ASSOCIATED(ZCLDTOP))  NULLIFY(ZCLDTOP_DIAG_dev_ptr)
         IF (ASSOCIATED(WESFC))    NULLIFY(WENTR_SFC_DIAG_dev_ptr)
         IF (ASSOCIATED(WERAD))    NULLIFY(WENTR_RAD_DIAG_dev_ptr)
         IF (ASSOCIATED(DBUOY))    NULLIFY(DEL_BUOY_DIAG_dev_ptr)
         IF (ASSOCIATED(VSCSFC))   NULLIFY(VSFC_DIAG_dev_ptr)
         IF (ASSOCIATED(VSCRAD))   NULLIFY(VRAD_DIAG_dev_ptr)
         IF (ASSOCIATED(KERAD))    NULLIFY(KENTRAD_DIAG_dev_ptr)
         IF (ASSOCIATED(VSCBRV))   NULLIFY(VBRV_DIAG_dev_ptr)
         IF (ASSOCIATED(WEBRV))    NULLIFY(WENTR_BRV_DIAG_dev_ptr)
         IF (ASSOCIATED(DSIEMS))   NULLIFY(DSIEMS_DIAG_dev_ptr)
         IF (ASSOCIATED(CHIS))     NULLIFY(CHIS_DIAG_dev_ptr)
         IF (ASSOCIATED(DELSINV))  NULLIFY(DELSINV_DIAG_dev_ptr)
         IF (ASSOCIATED(SMIXT))    NULLIFY(SLMIXTURE_DIAG_dev_ptr)
         IF (ASSOCIATED(CLDRF))    NULLIFY(CLDRADF_DIAG_dev_ptr)
         IF (ASSOCIATED(RADRCODE)) NULLIFY(RADRCODE_DIAG_dev_ptr)

         call MAPL_TimerOff(MAPL,name="----LOCK_DEALLOC" ,__RC__)

#else

!   ...then add Lock.
!--------------------

         CALL ENTRAIN(IM,JM,LM,                 &
                      ! Inputs
                      RADLW,                    &
                      USTAR,                    &
                      BSTAR,                    &
                      FRLAND,                   &
                      T,                        &
                      Q,                        &
                      QLLS,                     &
                      QILS,                     &
                      U,                        &
                      V,                        &
                      Z,                        &
                      PLO,                      &
                      ZLE,                      &
                      PLE,                      &
                      ! Inoutputs
                      KM,                       &
                      KH,                       &
                      ! Outputs
                      EKM,                      &
                      EKH,                      &
                      KHSFC,                    &
                      KHRAD,                    &
                      ZCLD,                     &
                      ZRADML,                   &
                      ZRADBS,                   &
                      ZSML,                     &
                      ! Diagnostics
                      ZCLDTOP,                  &
                      WESFC,                    &
                      WERAD,                    &
                      DBUOY,                    &
                      VSCSFC,                   &
                      VSCRAD,                   &
                      KERAD,                    &
                      VSCBRV,                   &
                      WEBRV,                    &
                      DSIEMS,                   &
                      CHIS,                     &
                      DELSINV,                  &
                      SMIXT,                    &
                      CLDRF,                    &
                      RADRCODE,                 &
                      ! Input parameter constants
                      PRANDTLSFC, PRANDTLRAD,   &
                      BETA_SURF, BETA_RAD,      &
                      TPFAC_SURF, ENTRATE_SURF, &
                      PCEFF_SURF, KHRADFAC, KHSFCFAC_LND, KHSFCFAC_OCN )

#endif

      else ! Not running ENTRAIN...
         EKM   = 0.0
         EKH   = 0.0
         KHSFC = 0.0
         KHRAD = 0.0
      end if DO_ENTRAIN

      call MAPL_TimerOff(MAPL,name="---LOCK" ,RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn (MAPL,"---POSTLOCK")

      KPBLMIN  = count(PREF < 50000.)

                            ZPBL = MAPL_UNDEF
      if (associated(PPBL)) PPBL = MAPL_UNDEF

      if (CALC_TCZPBL) then
         TCZPBL = MAPL_UNDEF

       if (LM .eq. 72) then
         thetavs = T(:,:,LM)*(1.0+MAPL_VIREPS*Q(:,:,LM)/(1.0-Q(:,:,LM)))*(TH(:,:,LM)/T(:,:,LM))
         tcrib(:,:,LM) = 0.0
         do I = 1, IM
            do J = 1, JM
               do L=LM-1,1,-1
                  thetavh(I,J) = T(I,J,L)*(1.0+MAPL_VIREPS*Q(I,J,L)/(1.0-Q(I,J,L)))*(TH(I,J,L)/T(I,J,L))
                  uv2h(I,J) = max(U(I,J,L)**2+V(I,J,L)**2,1.0E-8)
                  tcrib(I,J,L) = MAPL_GRAV*(thetavh(I,J)-thetavs(I,J))*Z(I,J,L)/(thetavs(I,J)*uv2h(I,J))
                  if (tcrib(I,J,L) >= tcri_crit) then
                     TCZPBL(I,J) = Z(I,J,L+1)+(tcri_crit-tcrib(I,J,L+1))/(tcrib(I,J,L)-tcrib(I,J,L+1))*(Z(I,J,L)-Z(I,J,L+1))
                     KPBLTC(I,J) = float(L)
                     exit
                  end if
               end do
            end do
         end do
       else
         tcrib(:,:,LM) = 0.0
         do I = 1, IM
            do J = 1, JM
               do L=LM-1,1,-1
                  uv2h(I,J) = max(U(I,J,L)**2+V(I,J,L)**2,1.0E-8)
                  tcrib(I,J,L) = MAPL_GRAV*(THV(I,J,L)-THV(I,J,LM))*Z(I,J,L)/(THV(I,J,LM)*uv2h(I,J))
                  if (tcrib(I,J,L) >= tcri_crit) then
                     TCZPBL(I,J) = Z(I,J,L+1)+(tcri_crit-tcrib(I,J,L+1))/(tcrib(I,J,L)-tcrib(I,J,L+1))*(Z(I,J,L)-Z(I,J,L+1))
                     KPBLTC(I,J) = float(L)
                     exit
                  end if
               end do
            end do
         end do
       endif

         where (TCZPBL<0.)
            TCZPBL = Z(:,:,LM)
            KPBLTC = float(LM)
         end where
      end if ! CALC_TCZPBL

      if (CALC_ZPBL2) then
         ZPBL2 = MAPL_UNDEF

         do I = 1, IM
            do J = 1, JM
               do L=LM,2,-1
                  if ((KH(I,J,L-1) < 2.).and.(KH(I,J,L) >= 2.).and.(ZPBL2(I,J)==MAPL_UNDEF)) then
                     ZPBL2(I,J) = Z(I,J,L)
                     KPBL2(I,J) = float(L)
                  end if
               end do
            end do
         end do

         where ( ZPBL2 .eq. MAPL_UNDEF )
            ZPBL2 = Z(:,:,LM)
            KPBL2 = float(LM)
         end where
         ZPBL2 = MIN(ZPBL2,Z(:,:,KPBLMIN))
      end if ! CALC_ZPBL2


     ! Calc KPBL using surface turbulence, for use in shallow scheme
      KPBL_SC = MAPL_UNDEF

      do I = 1, IM
         do J = 1, JM
            temparray(1:LM+1) = KHSFC(I,J,0:LM)

            maxkh = maxval(temparray)
            do L=LM-1,2,-1
               if ( (temparray(L) < 0.1*maxkh) .and. (temparray(L+1) >= 0.1*maxkh)  &
               .and. (KPBL_SC(I,J) == MAPL_UNDEF ) ) then
                  KPBL_SC(I,J) = float(L)
               end if
            end do
            if (  KPBL_SC(I,J) .eq. MAPL_UNDEF .or. (maxkh.lt.1.)) then
               KPBL_SC(I,J) = float(LM)
            endif
         end do
      end do


      if (CALC_ZPBL10p) then
         ZPBL10p = MAPL_UNDEF

         do I = 1, IM
            do J = 1, JM
               temparray(1:LM+1) = KH(I,J,0:LM)
               do L = LM,2,-1
                  locmax = maxloc(temparray,1)
                  minlval = max(0.001,0.0001*maxval(temparray))
                  if(temparray(locmax-1)<minlval.and.temparray(locmax+1)<minlval) temparray(locmax) = minlval
               enddo
               maxkh = temparray(LM)
               do L = LM-1,2,-1
                  if(temparray(L)>maxkh) maxkh = temparray(L)
                  if(temparray(L-1)<minlval) exit
               end do
               do L=LM-1,2,-1
                  if ( (temparray(L) < 0.1*maxkh) .and. (temparray(L+1) >= 0.1*maxkh)  &
                  .and. (ZPBL10p(I,J) == MAPL_UNDEF ) ) then
                     ZPBL10p(I,J) = ZLE(I,J,L)+ &
                  ((ZLE(I,J,L-1)-ZLE(I,J,L))/(temparray(L)-temparray(L+1))) * (0.1*maxkh-temparray(L+1))
                     KPBL10p(I,J) = float(L)
                  end if
               end do
               if (  ZPBL10p(I,J) .eq. MAPL_UNDEF .or. (maxkh.lt.1.)) then
                  ZPBL10p(I,J) = Z(I,J,LM)
                  KPBL10p(I,J) = float(LM)
               endif
            end do
         end do

         ZPBL10p = MIN(ZPBL10p,Z(:,:,KPBLMIN))
      end if ! CALC_ZPBL10p

      ! HTKE pbl height
      if (associated(ZPBLHTKE)) then
         ZPBLHTKE = MAPL_UNDEF
      end if ! ZPBLHTKE

      ! TKE 
      if (associated(TKE)) then         ! Reminder: TKE is on model edges
        if (DO_SHOC /= 0) then
          TKE = TKESHOC
        else
          TKE = MAPL_UNDEF
          do L = 1,LM-1
            TKE(:,:,L) = ( LAMBDADISS * &
            ( -1.*(KH(:,:,L)*MAPL_GRAV/((THV(:,:,L) + THV(:,:,L+1))*0.5) *  ((THV(:,:,L) - THV(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1)))) +  &
            (KM(:,:,L)*((U(:,:,L) - U(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1)))*((U(:,:,L) - U(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1))))  +  &
            (KM(:,:,L)*((V(:,:,L) - V(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1)))*((V(:,:,L) - V(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1)))) )) ** 2
            TKE(:,:,L) = TKE(:,:,L) ** (1./3.)
          enddo
        end if
      end if ! TKE

      ! RI local diagnostic for pbl height thresh 0.
      if (associated(ZPBLRI)) then
         ZPBLRI = MAPL_UNDEF
         where (RI(:,:,LM-1)>ri_crit) ZPBLRI = Z(:,:,LM)

         do I = 1, IM
            do J = 1, JM
               do L=LM-1,1,-1
                  if( (RI(I,J,L-1)>ri_crit) .and. (ZPBLRI(I,J) == MAPL_UNDEF) ) then
                     ZPBLRI(I,J) = Z(I,J,L+1)+(ri_crit-RI(I,J,L))/(RI(I,J,L-1)-RI(I,J,L))*(Z(I,J,L)-Z(I,J,L+1))
                  end if
               end do
            end do 
         end do 

         where ( ZPBLRI .eq. MAPL_UNDEF ) ZPBLRI = Z(:,:,LM)
         ZPBLRI = MIN(ZPBLRI,Z(:,:,KPBLMIN))
         where ( ZPBLRI < 0.0 ) ZPBLRI = Z(:,:,LM)
      end if ! ZPBLRI

      ! RI local diagnostic for pbl height thresh 0.2
      if (associated(ZPBLRI2)) then
         ZPBLRI2 = MAPL_UNDEF
         where (RI(:,:,LM-1) > ri_crit2) ZPBLRI2 = Z(:,:,LM)

         do I = 1, IM
            do J = 1, JM
               do L=LM-1,1,-1
                  if( (RI(I,J,L-1)>ri_crit2) .and. (ZPBLRI2(I,J) == MAPL_UNDEF) ) then
                     ZPBLRI2(I,J) = Z(I,J,L+1)+(ri_crit2-RI(I,J,L))/(RI(I,J,L-1)-RI(I,J,L))*(Z(I,J,L)-Z(I,J,L+1))
                  end if
               end do
            end do
         end do

         where ( ZPBLRI2 .eq. MAPL_UNDEF ) ZPBLRI2 = Z(:,:,LM)
         ZPBLRI2 = MIN(ZPBLRI2,Z(:,:,KPBLMIN))
         where ( ZPBLRI2 < 0.0 ) ZPBLRI2 = Z(:,:,LM)
      end if ! ZPBLRI2

      ! thetav gradient based pbl height diagnostic
      if (associated(ZPBLTHV)) then
         ZPBLTHV = MAPL_UNDEF

         do I = 1, IM
            do J = 1, JM

               do L=LM,1,-1
                  thetav(L) = T(I,J,L)*(1.0*MAPL_VIREPS*Q(I,J,L)/(1.0-Q(I,J,L)))*(TH(I,J,L)/T(I,J,L))
               end do

               maxdthvdz = 0

               do L=LM-1,1,-1
                  if(Z(I,J,L)<=Z(I,J,KPBLMIN)) then
                     dthvdz = (thetav(L+1)-thetav(L))/(Z(I,J,L+1)-Z(I,J,L))
                     if(dthvdz>maxdthvdz) then
                        maxdthvdz = dthvdz
                        ZPBLTHV(I,J) = 0.5*(Z(I,J,L+1)+Z(I,J,L))
                     end if
                  end if
               end do

            end do 
         end do 
      end if ! ZPBLTHV

      SELECT CASE(PBLHT_OPTION)

      CASE( 1 )
         ZPBL = ZPBL2
         KPBL = KPBL2

      CASE( 2 )
         ZPBL = ZPBL10p
         KPBL = KPBL10P

      CASE( 3 )
         ZPBL = TCZPBL
         KPBL = KPBLTC

      CASE( 4 )
         WHERE (FRLAND(:,:)>0)
            ZPBL = TCZPBL
            KPBL = KPBLTC

         ELSEWHERE
            ZPBL = ZPBL10p
            KPBL = KPBL10P

         END WHERE

      END SELECT

      ZPBL = MIN(ZPBL,Z(:,:,KPBLMIN))
      KPBL = MAX(KPBL,float(KPBLMIN))

      if (associated(PPBL)) then
         do I = 1, IM
            do J = 1, JM
               PPBL(I,J) = PLO(I,J,nint(KPBL(I,J)))
            end do
         end do
         PPBL = MAX(PPBL,PLO(:,:,KPBLMIN))
      end if

      ! Second difference coefficients for scalars; RDZ is RHO/DZ, DMI is (G DT)/DP
      ! ---------------------------------------------------------------------------

      CKS(:,:,1:LM-1) = -KH(:,:,1:LM-1) * RDZ(:,:,1:LM-1)
      AKS(:,:,1     ) = 0.0
      AKS(:,:,2:LM  ) = CKS(:,:,1:LM-1) * DMI(:,:,2:LM  )
      CKS(:,:,1:LM-1) = CKS(:,:,1:LM-1) * DMI(:,:,1:LM-1)
      CKS(:,:,  LM  ) = -CT             * DMI(:,:,  LM  )

      ! Fill KH at level LM+1 with CT * RDZ for diagnostic output
      ! ---------------------------------------------------------

      KH(:,:,LM) = CT * (PLE(:,:,LM)/(MAPL_RGAS * TV(:,:,LM))) / Z(:,:,LM)

      ! Water vapor can differ at the surface
      !--------------------------------------

      AKQ         = AKS
      CKQ         = CKS
      CKQ(:,:,LM) = -CQ * DMI(:,:,LM)

      ! Second difference coefficients for winds
      ! EKV is saved to use in the frictional heating calc.
      ! ---------------------------------------------------
      
      EKV(:,:,1:LM-1) = -KM(:,:,1:LM-1) * RDZ(:,:,1:LM-1)
      AKV(:,:,1     ) = 0.0
      AKV(:,:,2:LM  ) = EKV(:,:,1:LM-1) * DMI(:,:,2:LM  )
      CKV(:,:,1:LM-1) = EKV(:,:,1:LM-1) * DMI(:,:,1:LM-1)
      EKV(:,:,1:LM-1) = -MAPL_GRAV      * EKV(:,:,1:LM-1)

      CKV(:,:,  LM  ) = -  CU           * DMI(:,:,  LM  )
      EKV(:,:,  LM  ) =  MAPL_GRAV      * CU

      ! Fill KM at level LM with CU * RDZ for diagnostic output
      ! -------------------------------------------------------

      KM(:,:,LM) = CU * (PLE(:,:,LM)/(MAPL_RGAS * TV(:,:,LM))) / Z(:,:,LM)

      ! Setup the tridiagonal matrix
      ! ----------------------------

      BKS = 1.00 - (AKS+CKS)
      BKQ = 1.00 - (AKQ+CKQ)
      BKV = 1.00 - (AKV+CKV)

      ! Add the topographic roughness term
      ! ----------------------------------

      if (associated(AKSODT)) then
         AKSODT = -AKS/DT
         AKSODT(:,:,1) = 0.0
      end if

      if (associated(CKSODT)) then
         CKSODT = -CKS/DT
         CKSODT(:,:,LM) = 0.0
      end if

      if (associated(AKQODT)) then
         AKQODT = -AKQ/DT
         AKQODT(:,:,1) = 0.0
      end if

      if (associated(CKQODT)) then
         CKQODT = -CKQ/DT
         CKQODT(:,:,LM) = 0.0
      end if

      if (associated(AKVODT)) AKVODT = -AKV/DT
      if (associated(CKVODT)) CKVODT = -CKV/DT

      call MAPL_TimerOff(MAPL,"---POSTLOCK")

!BOP
!
!   Orograpghic drag follows  Beljaars (2003):
!   $$
!   \frac{\partial}{\partial z}\frac{\tau}{\rho} = \frac{C_B}{\lambda_B} |U(z)| U(z) 
!          e^{-\tilde{z}^\frac{3}{2}}\tilde{z}^{-1.2},
!   $$
!   where $z$ is the height above the surface in meters, 
!   $\tilde{z}=\frac{z}{\lambda_B}$, $\tau$ is the orographic stress at $z$,
!   $\rho$ is the air density, $U(z)$ is the wind velocity, and $\lambda_B$ is a vertical length scale.
!   Beljaars uses $\lambda_B = 1500$m, for which the non-dimensional parameter $C_B = 2.5101471 \times 10^{-8}$.
!   These are the default values, but both can be modified from the configuration. To avoid underflow.
!   the tendency is set to zero once $\tilde{z}$ exceeds 4 (i.e., 6 km from the surface for default values). 
!
!EOP

      call MAPL_TimerOn(MAPL,"---BELJAARS")

      call BELJAARS(IM, JM, LM, DT, &
                    LAMBDA_B, C_B,  &
                    U, V, Z,        &
                    VARFLT, PLE,    &
                    BKV, FKV        )

      call MAPL_TimerOff(MAPL,"---BELJAARS")

      call MAPL_TimerOn(MAPL,"---DECOMP")

! Do LU decomposition; C is not modified.
! On exit, B is the main diagonals of the LU
! decomposition, and A is the r.h.s multiplier.
!----------------------------------------------

      AKX = AKS
      BKX = BKS
      call VTRILU(AKX,BKX,CKS)
      AKS = AKX
      BKS = BKX

      AKX = AKQ
      BKX = BKQ
      call VTRILU(AKX,BKX,CKQ)
      AKQ = AKX
      BKQ = BKX

      AKX = AKV
      BKX = BKV
      call VTRILU(AKX,BKX,CKV)
      AKV = AKX
      BKV = BKX

! Get the sensitivity of solution to a unit
! change in the surface value. B and C are
! not modified.
!------------------------------------------

      call VTRISOLVESURF(BKS,CKS,DKS)
      call VTRISOLVESURF(BKQ,CKQ,DKQ)
      call VTRISOLVESURF(BKV,CKV,DKV)

      call MAPL_TimerOff(MAPL,"---DECOMP")

      if(ALLOC_TCZPBL) deallocate(TCZPBL)
      if(ALLOC_ZPBL2) deallocate(ZPBL2)
      if(ALLOC_ZPBL10p) deallocate(ZPBL10p)

      RETURN_(ESMF_SUCCESS)
     end subroutine REFRESH

!=============================================================================
!=============================================================================

!BOP

! !CROUTINE: DIFFUSE -- Solves for semi-implicit diffusive tendencies assuming fixed surface conditions.  

! !INTERFACE:

  subroutine DIFFUSE(IM,JM,LM,RC)

! !ARGUMENTS:

    integer,           intent(IN)       :: IM,JM,LM
    integer, optional, intent(OUT)      :: RC

! !DESCRIPTION: {\tt DIFFUSE} computes semi-implicit tendencies of all fields in
!  the TR bundle. Each field is examined for three attributes: {\tt DiffuseLike},
!  {\tt FriendlyToTURBULENCE}, and {\tt WeightedTendency}. These determine the behavior of
!  {\tt DIFFUSE} for that field. {\tt DiffuseLike} can be either 'U', 'Q', or 'S'; the default is 'Q'.
!  {\tt FriendlyToTURBULENCE}, and {\tt WeightedTendency} are ESMF logicals.
!  If {\tt FriendlyToTURBULENCE} is true, the field in TR is updated directly; otherwise
!  it is left untouched. In either case, If the corresponding pointer TRI bundle is associated, the
!  tendencies are returned there. If {\tt WeightedTendency} is true, the tendency in TRI, if any,
!  is pressure weighted.

!EOP

    character(len=ESMF_MAXSTR)          :: IAm='Diffuse'
    integer                             :: STATUS

    character(len=ESMF_MAXSTR)          :: TYPE
    character(len=ESMF_MAXSTR)          :: NAME
    type (ESMF_Field)                   :: FIELD
    type (ESMF_Array)                   :: ARRAY
    type (ESMF_FieldBundle)             :: TR
    type (ESMF_FieldBundle)             :: TRI
    type (ESMF_FieldBundle)             :: TRG
    type (ESMF_FieldBundle)             :: FSTAR
    type (ESMF_FieldBundle)             :: DFSTAR
    real, dimension(:,:,:), pointer     :: S, SOI, SOD
    real, dimension(:,:),   pointer     :: SG, SF, SDF, CX, SRG
    real, dimension(:,:,:), pointer     :: DX
    real, dimension(:,:,:), pointer     :: AK, BK, CK

    integer                             :: KM, K,L
    logical                             :: FRIENDLY
    logical                             :: WEIGHTED

    real,               dimension(IM,JM,LM) :: DP
    real(kind=MAPL_R8), dimension(IM,JM,LM) :: SX

! AMM pointer to export of S after diffuse
    real, dimension(:,:,:), pointer     :: SAFDIFFUSE

! Get the bundles containing the quantities to be diffused, 
!     their tendencies, their surface values, their surface
!     fluxes, and the derivatives of their surface fluxes
!     wrt the surface values. 
!----------------------------------------------------------

    call ESMF_StateGet(IMPORT, 'TR' ,    TR,     RC=STATUS); VERIFY_(STATUS)
    call ESMF_StateGet(IMPORT, 'TRG',    TRG,    RC=STATUS); VERIFY_(STATUS)

    call ESMF_StateGet(EXPORT, 'TRI',    TRI,    RC=STATUS); VERIFY_(STATUS)
    call ESMF_StateGet(EXPORT, 'FSTAR',  FSTAR,  RC=STATUS); VERIFY_(STATUS)
    call ESMF_StateGet(EXPORT, 'DFSTAR', DFSTAR, RC=STATUS); VERIFY_(STATUS)

! AMM pointer to export of S that diffuse sees
    call MAPL_GetPointer(EXPORT, SAFDIFFUSE ,  'SAFDIFFUSE' , RC=STATUS); VERIFY_(STATUS)

! Count the firlds in TR...
!--------------------------

    call ESMF_FieldBundleGet(TR, fieldCOUNT=KM, RC=STATUS)
    VERIFY_(STATUS)

! ...and make sure the other bundles are the same.
!-------------------------------------------------

    call ESMF_FieldBundleGet(TRI,    FieldCount=K , RC=STATUS)
    VERIFY_(STATUS)
    _ASSERT(KM==K,'needs informative message')
    call ESMF_FieldBundleGet(TRG,    FieldCount=K , RC=STATUS)
    VERIFY_(STATUS)
    _ASSERT(KM==K,'needs informative message')
    call ESMF_FieldBundleGet(FSTAR,  FieldCount=K , RC=STATUS)
    VERIFY_(STATUS)
    _ASSERT(KM==K,'needs informative message')
    call ESMF_FieldBundleGet(DFSTAR, FieldCount=K , RC=STATUS)
    VERIFY_(STATUS)
    _ASSERT(KM==K,'needs informative message')

! Pressure thickness of layers
!-----------------------------

    DP = PLE(:,:,1:LM)-PLE(:,:,0:LM-1)

! Loop over all quantities to be diffused.
!----------------------------------------

    do K=1,KM

! Get the Kth Field and its name from tracer bundle
!--------------------------------------------------

       call ESMF_FieldBundleGet(TR, K, FIELD, RC=STATUS)
       VERIFY_(STATUS)

       call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
       VERIFY_(STATUS)

! Get item's diffusion type (U, S or Q; default is Q)
!----------------------------------------------------

       call ESMF_AttributeGet(FIELD, NAME="DiffuseLike",         &
            VALUE=TYPE, DEFAULTVALUE=dflt_q,    RC=STATUS)
       VERIFY_(STATUS)

! Get item's friendly status (default is not friendly)
!-----------------------------------------------------

       call ESMF_AttributeGet(FIELD, NAME="FriendlyToTURBULENCE", &
            VALUE=FRIENDLY, DEFAULTVALUE=dflt_false,    RC=STATUS)
       VERIFY_(STATUS)

! Get item's weighting (default is unweighted tendencies)
!--------------------------------------------------------

       call ESMF_AttributeGet(FIELD, NAME="WeightedTendency",   &
            VALUE=WEIGHTED, DEFAULTVALUE=dflt_false,    RC=STATUS)
       VERIFY_(STATUS)

! Get pointer to the quantity, its tendency, its surface value,
!   the surface flux, and the sensitivity of the surface flux.
! -------------------------------------------------------------

       call ESMFL_BundleGetPointerToData(TR    ,      NAME,         S  , RC=STATUS)
       VERIFY_(STATUS)
       call ESMFL_BundleGetPointerToData(TRI   , trim(NAME)//'IT' , SOI, RC=STATUS)
       VERIFY_(STATUS)
       call ESMFL_BundleGetPointerToData(TRG   , trim(NAME)//'HAT', SRG, RC=STATUS)
       VERIFY_(STATUS)
       call ESMFL_BundleGetPointerToData(FSTAR , trim(NAME)//'FLX', SF , RC=STATUS)
       VERIFY_(STATUS)
       call ESMFL_BundleGetPointerToData(DFSTAR, trim(NAME)//'DFL', SDF, RC=STATUS)
       VERIFY_(STATUS)

! The quantity must exist; others are optional.
!----------------------------------------------

       _ASSERT(associated(S ),'needs informative message')

! If the surface values does not exists, we assume zero flux.
!------------------------------------------------------------
       
       if(associated(SRG)) then
          SG => SRG
       else
          allocate (SG(0,0), stat=STATUS)
          VERIFY_(STATUS)
       end if

! Pick the right exchange coefficients
!-------------------------------------

       if     ( TYPE=='U' ) then ! Momentum
          CX => CU
          DX => DKV
          AK => AKV; BK => BKV; CK => CKV
       else if( TYPE=='Q' ) then ! Water Vapor or other tracers
          CX => CQ
          DX => DKQ
          AK => AKQ; BK => BKQ; CK => CKQ
       else if( TYPE=='S' ) then ! Heat
          CX => CT
          DX => DKS
          AK => AKS; BK => BKS; CK => CKS
       else
          RETURN_(ESMF_FAILURE)
       endif

! Copy diffused quantity to temp buffer
! ------------------------------------------
       
       SX = S

! Solve for semi-implicit changes. This modifies SX
! -------------------------------------------------

       call VTRISOLVE(AK,BK,CK,SX,SG)

! Compute the surface fluxes
!---------------------------

       if(associated(SF)) then
          if(size(SG)>0) then
             SF = CX*(SG - SX(:,:,LM))
          else
             SF = 0.0
          end if
       end if

! Create tendencies
!------------------

       if(associated(SOI)) then
          if( WEIGHTED ) then
             SOI = ( (SX - S)/DT )*DP
          else
             SOI = ( (SX - S)/DT )
          endif
       end if

       if( NAME=='S' ) then
          SINC = ( (SX - S)/DT )
       end if

! Update friendlies
!------------------

       if(FRIENDLY) then
          S = SX
       end if
!AMM sync up T and Q - so write on S after diffusion
       if( TYPE=='S' ) then 
          if(associated(SAFDIFFUSE)) SAFDIFFUSE = SX
       endif

! Compute the derivative of the surface flux wrt the surface value
!-----------------------------------------------------------------

       if(associated(SDF)) then
          SDF = CX * (1.0-DX(:,:,LM))
       endif

       if(.not.associated(SRG)) then
          deallocate (SG)
       end if

    enddo ! End loop over all quantities to be diffused
! -----------------------------------------------------

    RETURN_(ESMF_SUCCESS)
  end subroutine DIFFUSE

end subroutine RUN1


!*********************************************************************
!*********************************************************************
!*********************************************************************


!BOP

! !IROUTINE: RUN2 -- The second run stage for the TURBULENCE component

! !INTERFACE:

  subroutine RUN2 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Second run stage of {\tt GEOS\_TurbulenceGridComp} performs
!    the updates due to changes in surface quantities. Its input are the changes in 
!    surface quantities during the time step. It can also compute the frictional
!    dissipation terms as exports, but these are not added to the temperatures.


!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp), pointer       :: MAPL
    type (ESMF_Config      )            :: CF
    type (ESMF_State       )            :: INTERNAL 

! Local variables

    integer                             :: IM, JM, LM
    real                                :: DT

    real, pointer, dimension(:,:)       :: VARFLT
    real, pointer, dimension(:,:)       :: LATS

! Begin... 
!---------

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'Run2'

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"-RUN2")

! Get parameters from generic state.
!-----------------------------------

          call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,          &
                               LATS = LATS,                  &
                               INTERNAL_ESMF_STATE=INTERNAL, &
                                                   RC=STATUS )
      VERIFY_(STATUS)

! Get configuration from component
!---------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

! Get application's timestep from configuration
!----------------------------------------------

    call ESMF_ConfigGetAttribute( CF, DT, Label="RUN_DT:" , RC=STATUS)
    VERIFY_(STATUS)


    call MAPL_GetPointer(IMPORT,VARFLT,  'VARFLT', RC=STATUS)
    VERIFY_(STATUS)

! Solve the free atmosphere problem
! ---------------------------------

    call MAPL_TimerOn (MAPL,"--UPDATE")
      call UPDATE(IM,JM,LM,LATS,RC=STATUS)
      VERIFY_(STATUS)
    call MAPL_TimerOff(MAPL,"--UPDATE")

!  All done with RUN
!-------------------

    call MAPL_TimerOff(MAPL,"-RUN2")
    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  contains

!BOP

! !CROUTINE: UPDATE -- Updates diffusive effects for changes at surface.

! !INTERFACE:

    subroutine UPDATE(IM,JM,LM,LATS,RC)

! !ARGUMENTS:

      integer,           intent(IN)       :: IM,JM,LM
      integer, optional, intent(OUT)      :: RC

! !DESCRIPTION: 
!    Some description

!EOP
  
 
      character(len=ESMF_MAXSTR)          :: IAm='Update'
      integer                             :: STATUS

      character(len=ESMF_MAXSTR)          :: TYPE
      character(len=ESMF_MAXSTR)          :: NAME
      type (ESMF_Field)                   :: FIELD
      type (ESMF_FieldBundle)             :: TR
      type (ESMF_FieldBundle)             :: TRI
      type (ESMF_FieldBundle)             :: DTG
      type (ESMF_FieldBundle)             :: FSTAR
      type (ESMF_FieldBundle)             :: DFSTAR
      real, dimension(:,:,:), pointer     :: PLE
      real, dimension(:,:,:), pointer     :: ZLE
      real, dimension(:,:,:), pointer     :: S, SOI, SINC, INTDIS, TOPDIS
      real, dimension(:,:  ), pointer     :: DSG, SF, SDF, SRFDIS
      real, dimension(:,:  ), pointer     :: HGTLM5, LM50M
      real, dimension(:,:  ), pointer     :: KETRB, KESRF, KETOP, KEINT
      real, dimension(:,:,:), pointer     :: DKS, DKV, DKQ, DKX, EKV, FKV
      real, dimension(:,:,:), pointer     :: DPDTTRB

      integer                             :: KM, K, L, I, J
      logical                             :: FRIENDLY
      logical                             :: WEIGHTED
      real, dimension(IM,JM,LM)           :: DP, SX
      real, dimension(IM,JM,LM-1)         :: DF
      integer, allocatable                :: KK(:)
! AMM pointer to export of S that update sees
      real, dimension(:,:,:), pointer     :: SAFUPDATE

! The following variables are for SHVC parameterization

      real, dimension(IM,JM,LM)           :: SOIOFS, XINC
      real,    dimension(IM,JM)           :: z500, z1500, z7000, STDV
      integer, dimension(IM,JM)           :: L500, L1500, L7000, L200
      integer, dimension(IM,JM)           :: LTOPS,LBOT,LTOPQ
      logical, dimension(IM,JM)           :: DidSHVC
      real                                :: REDUFAC, SUMSOI
      real                                :: SHVC_CRIT
      real                                :: SHVC_1500, SHVC_ZDEPTH
      real                                :: lat_in_degrees, lat_effect
      real,  dimension(IM,JM)             :: LATS
      real                                :: SHVC_ALPHA, SHVC_EFFECT, scaling
      logical                             :: DO_SHVC
      integer                             :: KS

      character(len=ESMF_MAXSTR) :: GRIDNAME
      character(len=4)           :: imchar
      character(len=2)           :: dateline
      integer                    :: imsize,nn

! Pressure-weighted dissipation heating rates
!--------------------------------------------

      call MAPL_GetPointer(EXPORT, KETRB ,  'KETRB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KESRF ,  'KESRF' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KETOP ,  'KETOP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KEINT ,  'KEINT' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, DPDTTRB,  'DPDTTRB', RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, SRFDIS,  'SRFDIS',                  &
                       alloc=associated(KETRB) .or. associated(KESRF), &
                                                              RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, INTDIS,  'INTDIS',                  &
                       alloc=associated(KETRB) .or. associated(KEINT), &
                                                              RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TOPDIS,  'TOPDIS',                  &
                       alloc=associated(KETRB) .or. associated(KETOP), &
                                                              RC=STATUS)
      VERIFY_(STATUS)

! SHVC Resource parameters. SHVC_EFFECT can be set to zero to turn-off SHVC.
! SHVC_EFFECT = 1. is the tuned value for  2 degree horizontal resolution.
! It should be set to a lower number at higher resolution.

      call MAPL_GetResource( MAPL, SHVC_EFFECT, 'SHVC_EFFECT:', default=0.   , RC=STATUS )
      VERIFY_(STATUS)

      DO_SHVC = SHVC_EFFECT > 0.0

      if(DO_SHVC) then
         call MAPL_GetResource( MAPL, SHVC_CRIT,   'SHVC_CRIT:'  , default=300. , RC=STATUS )
         VERIFY_(STATUS)
         call MAPL_GetResource( MAPL, SHVC_ALPHA,  'SHVC_ALPHA:' , default=1.   , RC=STATUS )
         VERIFY_(STATUS)
         call MAPL_GetResource( MAPL, SHVC_1500,   'SHVC_1500:'  , default=2100., RC=STATUS )
         VERIFY_(STATUS)
         call MAPL_GetResource( MAPL, SHVC_ZDEPTH, 'SHVC_ZDEPTH:', default=3500., RC=STATUS )
         VERIFY_(STATUS)

         call MAPL_GetResource( MAPL, GRIDNAME, 'AGCM_GRIDNAME:', RC=STATUS )
         VERIFY_(STATUS)
         GRIDNAME =  AdjustL(GRIDNAME)
               nn = len_trim(GRIDNAME)
         dateline = GRIDNAME(nn-1:nn)
           imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
           read(imchar,*) imsize
         if(dateline.eq.'CF') imsize  = imsize*4

         if( imsize.le.200       ) scaling = 1.0  !              Resolution >= 2.000-deg
         if( imsize.gt.200 .and. &
             imsize.le.400       ) scaling = 1.0  !  2.000-deg > Resolution >= 1.000-deg
         if( imsize.gt.400 .and. &
             imsize.le.800       ) scaling = 7.0  !  1.000-deg > Resolution >= 0.500-deg
         if( imsize.gt.800 .and. &
             imsize.le.1600      ) scaling = 7.0  !  0.500-deg > Resolution >= 0.250-deg
         if( imsize.gt.1600      ) scaling = 7.0  !  0.250-deg > Resolution 
      end if

! Get imports
!------------

      call MAPL_GetPointer(IMPORT,    PLE,     'PLE', RC=STATUS); VERIFY_(STATUS)

! Get the tendecy sensitivities computed in RUN1
!-----------------------------------------------

      call MAPL_GetPointer(INTERNAL, DKS,   'DKS',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, DKV,   'DKV',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, DKQ,   'DKQ',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, EKV,   'EKV',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, FKV,   'FKV',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, ZLE,   'ZLE',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, SINC,  'SINC',    RC=STATUS)
      VERIFY_(STATUS)

! Get the bundles containing the quantities to be diffused, 
!     their tendencies, their surface values, their surface
!     fluxes, and the derivatives of their surface fluxes
!     wrt the surface values. 
!----------------------------------------------------------

      call ESMF_StateGet(IMPORT, 'TR' ,    TR,     RC=STATUS); VERIFY_(STATUS)
      call ESMF_StateGet(IMPORT, 'DTG',    DTG,    RC=STATUS); VERIFY_(STATUS)

      call ESMF_StateGet(EXPORT, 'TRI',    TRI,    RC=STATUS); VERIFY_(STATUS)
      call ESMF_StateGet(EXPORT, 'FSTAR' , FSTAR,  RC=STATUS); VERIFY_(STATUS)
      call ESMF_StateGet(EXPORT, 'DFSTAR', DFSTAR, RC=STATUS); VERIFY_(STATUS)

! Count them...
!--------------

      call ESMF_FieldBundleGet(TR , FieldCount=KM, RC=STATUS)
      VERIFY_(STATUS)

! and make sure the other bundles are the same.
!----------------------------------------------

      call ESMF_FieldBundleGet(DTG, FieldCount=K , RC=STATUS)
      VERIFY_(STATUS)

      _ASSERT(KM==K,'needs informative message')

! KK gives the order in which quantities will be process.
!--------------------------------------------------------

      allocate(KK(KM), stat=STATUS)
      VERIFY_(STATUS)

      do K = 1,KM
         KK(K) = K
      end do

! Clear the accumulators for the dissipation.
!--------------------------------------------

      if(associated(SRFDIS)) SRFDIS = 0.0
      if(associated(INTDIS)) INTDIS = 0.0
      if(associated(TOPDIS)) TOPDIS = 0.0
      if(associated(KETRB )) KETRB  = 0.0
      if(associated(KESRF )) KESRF  = 0.0
      if(associated(KETOP )) KETOP  = 0.0
      if(associated(KEINT )) KEINT  = 0.0

! Pressure thickness of layers
!-----------------------------

      DP = PLE(:,:,1:LM)-PLE(:,:,0:LM-1)

! Diagnostics
      call MAPL_GetPointer(EXPORT, HGTLM5 ,  'HGTLM5' , RC=STATUS); VERIFY_(STATUS)
      if(associated(HGTLM5)) then
         HGTLM5 = ZLE(:,:,LM-5)
      end if
      call MAPL_GetPointer(EXPORT, LM50M ,  'LM50M' , RC=STATUS); VERIFY_(STATUS)
      if(associated(LM50M)) then
         LM50M = LM
         do L=LM,2,-1
            where (ZLE(:,:,L) <= 150. .and. ZLE(:,:,L-1) > 150.)
               LM50M=L-1
            endwhere
         enddo
      end if

      L200=LM
      do L=LM,2,-1
         where (ZLE(:,:,L) <= 200. .and. ZLE(:,:,L-1) > 200.)
            L200=L-1
         endwhere
      enddo

! Section 1 of 2. SHVC parameterization (W. Chao, J. Atmos. Sci., May 2012, P.1547) 
!  Defining the top and bottom levels of the heat and moisture redistribution layer
!----------------------------------------------------------------------------------

      SHVC_INIT: if(DO_SHVC) then

! Ensure that S is processed first. This only matters for SHVC
!-------------------------------------------------------------

         KS = 0

         do K = 1,KM
            call ESMF_FieldBundleGet(TR, K, FIELD, RC=STATUS)
            VERIFY_(STATUS)

            call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
            VERIFY_(STATUS)

            if    (NAME == 'S') then
               KS=KK(1); KK(1)=K; KK(K)=KS
            end if
         end do

         _ASSERT(KS /= 0 ,'needs informative message')

! SHVC super-layers
!------------------

         z500  =  500.
         z1500 = 1500.
         z7000 = 7000.

         STDV = sqrt(varflt*scaling)   ! Scaling VARFLT based on resolution

         where (STDV >=700.)
            z1500 = SHVC_1500                   
         endwhere

         where ( (STDV >300.) .and. (STDV <700.) )
            z1500 = 1500.+ (SHVC_1500-1500.)* (STDV - 300.)/400.
         endwhere  

         z7000 = z1500 + SHVC_ZDEPTH



         L500=1.
         do L=LM,2,-1
            where (ZLE(:,:,L) <= z500 .and. ZLE(:,:,L-1) > z500)     
               L500=L-1    
            endwhere
         enddo

         L1500=1.
         do L=LM,2,-1
            where (ZLE(:,:,L) <= z1500 .and. ZLE(:,:,L-1) > z1500)    
               L1500=L-1
            endwhere
         enddo

         L7000=1.
         do L=LM,2,-1
            where (ZLE(:,:,L) <= z7000 .and. ZLE(:,:,L-1) > z7000)    
               L7000=L-1
            endwhere
         enddo

         LBOT  = L1500-1         
         LTOPS = L7000
         LTOPQ = L1500-(LM-L500)*2

         SOIOFS = 0.0

      end if SHVC_INIT

! AMM pointer to export of S that diffuse sees
    call MAPL_GetPointer(EXPORT, SAFUPDATE ,  'SAFUPDATE' , RC=STATUS); VERIFY_(STATUS)

! Loop over all quantities to be diffused.
!-----------------------------------------

       TRACERS: do KS=1,KM

         K = KK(KS)

! Get Kth field from bundle
!--------------------------
          
         call ESMF_FieldBundleGet(TR, K, FIELD, RC=STATUS)
         VERIFY_(STATUS)

         call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
         VERIFY_(STATUS)

! Get item's diffusion type (U, S or Q; default is Q)
!----------------------------------------------------

         call ESMF_AttributeGet(FIELD, NAME="DiffuseLike",         &
              VALUE=TYPE, DEFAULTVALUE=dflt_Q,    RC=STATUS)
         VERIFY_(STATUS)

! Get item's friendly status (default is not friendly)
!-----------------------------------------------------

         call ESMF_AttributeGet(FIELD, NAME="FriendlyToTURBULENCE", &
              VALUE=Friendly, DEFAULTVALUE=dflt_false,    RC=STATUS)
         VERIFY_(STATUS)

! Get item's weighting (default is unweighted tendencies)
!--------------------------------------------------------

         call ESMF_AttributeGet(FIELD, NAME="WeightedTendency",    &
              VALUE=WEIGHTED, DEFAULTVALUE=dflt_false,    RC=STATUS)
         VERIFY_(STATUS)

! Get pointers to the quantity, its tendency, its surface increment,
!   the preliminary surface flux, and the sensitivity of the surface
!   flux to the surface value.
! ------------------------------------------------------------------

         call ESMFL_BundleGetPointerToData(TR    ,      NAME,         S  , RC=STATUS)
         VERIFY_(STATUS)
         call ESMFL_BundleGetPointerToData(TRI   , trim(NAME)//'IT' , SOI, RC=STATUS)
         VERIFY_(STATUS)
         call ESMFL_BundleGetPointerToData(DTG   , trim(NAME)//'DEL', DSG, RC=STATUS)
         VERIFY_(STATUS)
         call ESMFL_BundleGetPointerToData(FSTAR , trim(NAME)//'FLX', SF , RC=STATUS)
         VERIFY_(STATUS)
         call ESMFL_BundleGetPointerToData(DFSTAR, trim(NAME)//'DFL', SDF, RC=STATUS)
         VERIFY_(STATUS)

!  Point to the appropriate sensitivity
!--------------------------------------

         if      ( TYPE=='U' ) then
            DKX => DKV
         else if ( TYPE=='Q' ) then
            DKX => DKQ
         else if ( TYPE=='S' ) then
            DKX => DKS
         else
            RETURN_(ESMF_FAILURE)
         end if

! Update diffused quantity
!-------------------------

         SX = S

         if(associated(DSG)) then
            do L=1,LM 
               SX(:,:,L) = SX(:,:,L) + DKX(:,:,L)*DSG 
            end do
         end if

! Increment the dissipation
!-------------------------- 

         if( TYPE=='U' ) then
            if(associated(KETRB )) KETRB  = 0.0
            if(associated(KESRF )) KESRF  = 0.0
            if(associated(KETOP )) KETOP  = 0.0
            if(associated(KEINT )) KEINT  = 0.0
            if(associated(INTDIS)) then

               DF             = (0.5/(MAPL_CP))*EKV(:,:,1:LM-1)*(SX(:,:,1:LM-1)-SX(:,:,2:LM))**2
               INTDIS(:,:,1:LM-1) = INTDIS(:,:,1:LM-1) + DF
               INTDIS(:,:,2:LM  ) = INTDIS(:,:,2:LM  ) + DF

               !DF(:,:,1) = sum(DP(:,:,LM-10:LM),3)
               !DF(:,:,1) = ((1.0/(MAPL_CP))*EKV(:,:,LM)*SX(:,:,LM)**2)/DF(:,:,1)
               !do L=LM-10,LM
               !   INTDIS(:,:,L) = INTDIS(:,:,L) + DF(:,:,1)*DP(:,:,L)
               !end do

               ! Add surface dissipation to lower 200m
               do J=1,JM
                  do I=1,IM
                     DF(I,J,1) = sum(DP(I,J,L200(I,J):LM))
                     DF(I,J,1) = ((1.0/(MAPL_CP))*EKV(I,J,LM)*SX(I,J,LM)**2)/DF(I,J,1)
                  end do
               end do
               do J=1,JM
                  do I=1,IM
                     do L=L200(I,J),LM
                        INTDIS(I,J,L) = INTDIS(I,J,L) + DF(I,J,1)*DP(I,J,L)
                     end do
                  end do
               end do

               if(associated(KETRB)) then
                  do L=1,LM
                     KETRB = KETRB - INTDIS(:,:,L)* (MAPL_CP/MAPL_GRAV)
                  end do
               end if
               if(associated(KEINT)) then
                  do L=1,LM
                     KEINT = KEINT - INTDIS(:,:,L)* (MAPL_CP/MAPL_GRAV)
                  end do
               end if
            endif
            if(associated(TOPDIS)) then
               TOPDIS = TOPDIS + (1.0/(MAPL_CP))*FKV*SX**2
               if(associated(KETRB)) then
                  do L=1,LM
                     KETRB = KETRB - TOPDIS(:,:,L)* (MAPL_CP/MAPL_GRAV)
                  end do
               end if
               if(associated(KETOP)) then
                  do L=1,LM
                     KETOP = KETOP - TOPDIS(:,:,L)* (MAPL_CP/MAPL_GRAV)
                  end do
               end if
             endif
            if(associated(SRFDIS)) then
               SRFDIS = SRFDIS + (1.0/(MAPL_CP))*EKV(:,:,LM)*SX(:,:,LM)**2
               if(associated(KETRB)) KETRB = KETRB - SRFDIS* (MAPL_CP/MAPL_GRAV)
               if(associated(KESRF)) KESRF = KESRF - SRFDIS* (MAPL_CP/MAPL_GRAV)
            endif
         end if

! Update tendencies
! -----------------

         if(associated(SOI).and.associated(DSG)) then
            if( WEIGHTED ) then
               do L=1,LM
                  SOI(:,:,L) = SOI(:,:,L) +  (DKX(:,:,L)*DSG/DT)*DP(:,:,L)
               end do
            else
               do L=1,LM
                  SOI(:,:,L) = SOI(:,:,L) +  (DKX(:,:,L)*DSG/DT)
               end do
            endif
         end if

! Section 2 of 2. SHVC parameterization   (W. Chao,  J. Atmos. Sci., 2012, p1547)   
!  To use SHVC set SHVC_EFFECT in AGCM.rc to > 0.0.
!--------------------------------------------------------------------------------

         RUN_SHVC: if (DO_SHVC) then

            XINC = 0.0

            S_or_Q: if (NAME=='S') then

               if(associated(DSG)) then
                  do L=1,LM
                     SINC(:,:,L) = SINC(:,:,L) +  (DKX(:,:,L)*DSG/DT)
                  end do
               end if

               do I=1,IM
                  do J=1,JM
                lat_effect = 1.
                lat_in_degrees= ABS(LATS(I,J)/(3.14159/2.)*90.)
                if (lat_in_degrees >=42.) lat_effect=0.
                if (lat_in_degrees >37. .and. lat_in_degrees < 42.)   &
                lat_effect = 1.0 - (lat_in_degrees-37.)/(42.-37.)
                     if (STDV(I,J) > SHVC_CRIT) then

                        SUMSOI = sum(SINC(I,J,L500(I,J):LM)*DP(I,J,L500(I,J):LM))
                        DidSHVC(I,J) = SUMSOI >= 0.0

                        if (DidSHVC(I,J)) then
                           if (STDV(I,J) >= 800.) then
                              REDUFAC = 1.0
                           elseif (STDV(i,j) >700. .and. STDV(I,J) <800.) then
                              REDUFAC = 0.95 + 0.05*(STDV(I,J)-700.)/100.
                           else
                              REDUFAC = max(min((STDV(I,J)-SHVC_CRIT)/100.,0.95),0.0)
                           end if

                           REDUFAC = REDUFAC * SHVC_EFFECT  *lat_effect       

                           SUMSOI = 0.
                           do L=L500(i,j),LM
                              SUMSOI        = SUMSOI + SINC(I,J,L)*REDUFAC*DP(I,J,L)
                              XINC  (I,J,L) = -SINC(I,J,L) * REDUFAC
                              SOIOFS(I,J,L) =  XINC(I,J,L) / SX(I,J,L)
                           enddo   !do L

                           XINC(I,J,LTOPS(I,J):LBOT(I,J)) = SUMSOI/SUM(DP(I,J,LTOPS(I,J):LBOT(I,J)))
                        endif
                     else
                        DidSHVC(I,J) = .false.
                     endif   ! end of if (STDV>SHVC_CRIT)
                  enddo   !do J
               enddo   !do I

            elseif (NAME == 'Q')  then

! SHVC_ALPHA below is the alpha factor mentioned on page 1552 of Chao (2012, cited above)
!----------------------------------------------------------------------------------------

               do J=1,JM
                  do I=1,IM
                     if (DidSHVC(I,J)) then
                        SUMSOI = 0.
                        do L=L500(I,J),LM
                           XINC(I,J,L) = SHVC_ALPHA*(SOIOFS(I,J,L)*SX(I,J,L))
                           SUMSOI      = SUMSOI +  XINC(I,J,L)*DP(I,J,L)
                        enddo

                        XINC(I,J,LTOPQ(I,J):LBOT(I,J)) = - SUMSOI/SUM(DP(I,J,LTOPQ(I,J):LBOT(I,J)))
                     endif
                  enddo
               enddo

            endif S_or_Q

            if (name == 'S' .or. name == 'Q') then
               SX  = SX + XINC * DT

               if(associated(SOI)) then
                  if(WEIGHTED) then
                     SOI = SOI + XINC*DP
                  else
                     SOI = SOI + XINC
                  end if
               end if
            end if


         end if RUN_SHVC

! Replace friendly
!-----------------

         if(FRIENDLY) then
            S = SX
         end if
!AMM sync up T and Q - so write on S after diffusion
       if( name=='S' ) then 
          if(associated(SAFUPDATE)) SAFUPDATE = SX
       endif

! Update surface fluxes
! ---------------------

         if(associated(SF).and.associated(DSG)) then
            SF = SF + DSG*SDF
         end if

         if(associated(DPDTTRB)) then
            if( name=='Q' ) then
               DPDTTRB(:,:,1:LM-1) = 0.0
               DPDTTRB(:,:,LM)     = MAPL_GRAV*SF
            end if
         end if

      enddo TRACERS

! End loop over all quantities to be diffused
!--------------------------------------------

      deallocate(KK)

      RETURN_(ESMF_SUCCESS)
    end subroutine UPDATE

  end subroutine RUN2


!*********************************************************************
!*********************************************************************
!*********************************************************************

!*********************************************************************

!*********************************************************************

!BOP

! !IROUTINE:  LOUIS_KS -- Computes atmospheric diffusivities at interior levels

! !INTERFACE:

   subroutine LOUIS_KS(              &
         ZZ,ZE,PV,UU,VV,ZPBL,        &
         KH,KM,RI,DU,                & 
         LOUIS, MINSHEAR, MINTHICK,  &
         LAMBDAM, LAMBDAM2,          &
         LAMBDAH, LAMBDAH2,          &
         ALHFAC, ALMFAC,             &
         ZKMENV, ZKHENV, AKHMMAX,    &
         ALH_DIAG,KMLS_DIAG,KHLS_DIAG)

! !ARGUMENTS:

      ! Inputs
      real,    intent(IN   ) ::   ZZ(:,:,:) ! Height of layer center above the surface (m).
      real,    intent(IN   ) ::   PV(:,:,:) ! Virtual potential temperature at layer center (K).
      real,    intent(IN   ) ::   UU(:,:,:) ! Eastward velocity at layer center (m s-1).
      real,    intent(IN   ) ::   VV(:,:,:) ! Northward velocity at layer center (m s-1).
      real,    intent(IN   ) ::   ZE(:,:,:) ! Height of layer base above the surface (m).
      real,    intent(IN   ) :: ZPBL(:,:  ) ! PBL Depth (m)

      ! Outputs
      real,    intent(  OUT) ::   KM(:,:,:) ! Momentum diffusivity at base of each layer (m+2 s-1).
      real,    intent(  OUT) ::   KH(:,:,:) ! Heat diffusivity at base of each layer  (m+2 s-1).
      real,    intent(  OUT) ::   RI(:,:,:) ! Richardson number
      real,    intent(  OUT) ::   DU(:,:,:) ! Magnitude of wind shear (s-1).
   
      ! Diagnostic outputs
      real,    pointer       ::  ALH_DIAG(:,:,:) ! Blackadar Length Scale diagnostic (m) [Optional] 
      real,    pointer       :: KMLS_DIAG(:,:,:) ! Momentum diffusivity at base of each layer (m+2 s-1).
      real,    pointer       :: KHLS_DIAG(:,:,:) ! Heat diffusivity at base of each layer  (m+2 s-1).

      ! These are constants
      real,    intent(IN   ) :: LOUIS        ! Louis scheme parameters (usually 5).
      real,    intent(IN   ) :: MINSHEAR     ! Min shear allowed in Ri calculation (s-1).
      real,    intent(IN   ) :: MINTHICK     ! Min layer thickness (m).
      real,    intent(IN   ) :: LAMBDAM      ! Blackadar(1962) length scale parameter for momentum (m).
      real,    intent(IN   ) :: LAMBDAM2     ! Second Blackadar parameter for momentum (m).
      real,    intent(IN   ) :: LAMBDAH      ! Blackadar(1962) length scale parameter for heat (m).
      real,    intent(IN   ) :: LAMBDAH2     ! Second Blackadar parameter for heat (m).
      real,    intent(IN   ) :: ALHFAC
      real,    intent(IN   ) :: ALMFAC
      real,    intent(IN   ) :: ZKMENV       ! Transition height for Blackadar param for momentum (m)
      real,    intent(IN   ) :: ZKHENV       ! Transition height for Blackadar param for heat     (m)
      real,    intent(IN   ) :: AKHMMAX      ! Maximum allowe diffusivity (m+2 s-1).

! !DESCRIPTION: Computes Louis et al.(1979) Richardson-number-based diffusivites,
!                as well as an additional ``entrainment'' diffusivity.
!                The Louis diffusivities for momentum, $K_m$, and for heat
!   and moisture, $K_h$, are defined at the interior layer edges. For LM layers,
!   we define diffusivities at the base of the top LM-1 layers. All indexing
!   is from top to bottom of the atmosphere. 
!
!
!  The Richardson number, Ri, is defined at the same edges as the diffusivities. 
!  $$
!  {\rm Ri}_l = \frac{ \frac{g}{\left(\overline{\theta_v}\right)_l}\left(\frac{\delta \theta_v}{\delta z}\right)_l }
!                    { \left(\frac{\delta {\bf |V|}}{\delta z}\right)^2_l             }, \, \,  l=1,LM-1
!  $$
!  where $\theta_v=\theta(1+\epsilon q)$ is the virtual potential temperature,
!  $\epsilon=\frac{M_a}{M_w}-1$, $M_a$ and $M_w$ are the molecular weights of
!  dry air and water, and $q$ is the specific humidity.
!  $\delta \theta_v$ is the difference of $\theta_v$ in the layers above and below the edge 
!  at which Ri$_l$ is defined; $\overline{\theta_v}$ is their average.
!
!  The diffusivities at the layer edges have the form:
!  $$
!  K^m_l = (\ell^2_m)_l \left(\frac{\delta {\bf |V|}}{\delta z}\right)_l f_m({\rm Ri}_l)
!  $$
!  and
!  $$
!  K^h_l = (\ell^2_h)_l \left(\frac{\delta {\bf |V|}}{\delta z}\right)_l f_h({\rm Ri}_l),
!  $$
!  where $k$ is the Von Karman constant, and $\ell$ is the 
!  Blackdar(1962) length scale, also defined at the layer edges.
!
!  Different turbulent length scales can be used for heat and momentum. 
!  in both cases, we use the  traditional formulation:
!  $$
!  (\ell_{(m,h)})_l  = \frac{kz_l}{1 + \frac{kz_l}{\lambda_{(m,h)}}},
!  $$
!  where, near the surface, the scale is proportional to $z_l$, the height above 
!  the surface of edge level $l$, and far from the surface it approaches $\lambda$.
!  The length scale $\lambda$ is usually taken to be a constant (order 150 m), assuming
!  the same scale for the outre boundary layer and the free atmosphere. We make it
!  a function of height, reducing its value in the free atmosphere. The momentum
!  length scale written as:
!  $$
!  \lambda_m = \max(\lambda_1 e^{\left(\frac{z_l}{z_T}\right)^2}, \lambda_2)
!  $$
!  where $\lambda_2 \le \lambda_1$ and $z_T$ is the top of the boundary layer.
!  The length scale for heat and other scalers is taken as: $\lambda_h =  \sqrt\frac{3d}{2} \lambda_m$,
!  following the scheme used at ECMWF.
!
!  The two universal functions of the Richardson number,  $f_m$ and $f_h$,
!  are taken from Louis et al (1982). For unstable conditions (Ri$\le 0$),
!  they are:
!  $$
!  f_m = (1 - 2b \psi)
!  $$
!  and
!  $$
!  f_h = (1 - 3b \psi),
!  $$
!  where
!  $$
!  \psi = \frac{ {\rm Ri} }{ 1 + 3bC(z)\sqrt{-{\rm Ri}} },
!  $$
!  and
!  $$
!  C(z)=
!  $$

!  For stable condition (Ri$\ge 0$), they are
!  $$
!  f_m = \frac{1}{1.0 + \frac{2b{\rm Ri}}{\psi}}
!  $$
!  and
!  $$
!  f_h = \frac{1}{1.0 + 3b{\rm Ri}\psi},
!  $$
!  where
!  $$
!  \psi =  \sqrt{1+d{\rm Ri}}.
!  $$
!  As in Louis et al (1982), the parameters appearing in these are taken  
!  as $b = c = d = 5$. 


!EOP

! Locals

      real, dimension(size(KM,1),size(KM,2),size(KM,3)) :: ALH, ALM, DZ, DT, TM, PS, LAMBDAM_X, LAMBDAH_X
      real, dimension(size(KM,1),size(KM,2)           ) :: pbllocal

      integer :: L, LM
      !real    :: Zchoke 

! Begin...
!===>   Number of layers; edge levels will be one less (LM-1).

      LM = size(ZZ,3)

!===>   Initialize output arrays

      KH = 0.0
      KM = 0.0
      DU = 0.0
      RI = 0.0

!===>   Initialize pbllocal

      pbllocal = ZPBL
      where ( pbllocal .LE. ZZ(:,:,LM) ) pbllocal = ZZ(:,:,LM)

!===>   Quantities needed for Richardson number

      DZ(:,:,:) = (ZZ(:,:,1:LM-1) - ZZ(:,:,2:LM))
      TM(:,:,:) = (PV(:,:,1:LM-1) + PV(:,:,2:LM))*0.5
      DT(:,:,:) = (PV(:,:,1:LM-1) - PV(:,:,2:LM))
      DU(:,:,:) = (UU(:,:,1:LM-1) - UU(:,:,2:LM))**2 + &
                  (VV(:,:,1:LM-1) - VV(:,:,2:LM))**2

!===>   Limits on distance between layer centers and vertical shear at edges.

      DZ =  max(DZ, MINTHICK)
      DU = sqrt(DU)/DZ

!===>   Richardson number  ( RI = G*(DTheta_v/DZ) / (Theta_v*|DV/DZ|^2) )

      RI = MAPL_GRAV*(DT/DZ)/(TM*( max(DU, MINSHEAR)**2))

!===>   Blackadar(1962) length scale: $1/l = 1/(kz) + 1/\lambda$

!!!   LAMBDAM_X = MAX( LAMBDAM * EXP( -(ZE / ZKMENV )**2 ) , LAMBDAM2 )
!!!   LAMBDAH_X = MAX( LAMBDAH * EXP( -(ZE / ZKHENV )**2 ) , LAMBDAH2 )

      do L = 1, LM-1
         LAMBDAM_X(:,:,L) = MAX( 0.1 * pbllocal(:,:) * EXP( -(ZE(:,:,L) / ZKMENV )**2 ) , LAMBDAM2 )
         LAMBDAH_X(:,:,L) = MAX( 0.1 * pbllocal(:,:) * EXP( -(ZE(:,:,L) / ZKHENV )**2 ) , LAMBDAH2 )
      end do

      ALM = ALMFAC * ( MAPL_KARMAN*ZE/( 1.0 + MAPL_KARMAN*(ZE/LAMBDAM_X) ) )**2
      ALH = ALHFAC * ( MAPL_KARMAN*ZE/( 1.0 + MAPL_KARMAN*(ZE/LAMBDAH_X) ) )**2

      if (associated(ALH_DIAG)) ALH_DIAG(:,:,1:LM-1) = SQRT( ALH )

      where ( RI < 0.0 )
         PS = ( (ZZ(:,:,1:LM-1)/ZZ(:,:,2:LM))**(1./3.) - 1.0 ) ** 3
         PS = ALH*sqrt( PS/(ZE*(DZ**3)) )
         PS = RI/(1.0 + (3.0*LOUIS*LOUIS)*PS*sqrt(-RI))

         KH = 1.0 - (LOUIS*3.0)*PS
         KM = 1.0 - (LOUIS*2.0)*PS
      end where

!===>   Unstable case: Uses (3.14, 3.18, 3.27) in Louis-scheme
!                      should approach (3.13) for small -RI.

!===>   Choke off unstable KH below Zchoke (m). JTB 2/2/06    
!!!   Zchoke = 500.
!!!   where( (RI < 0.0) .and. (ZE < Zchoke ) )  
!!!      KH = KH * (( ZE / Zchoke )**3)            
!!!   endwhere

!===>   Stable case

      where ( RI >= 0.0 )
         PS = sqrt  (1.0 +  LOUIS     *RI   )

         KH = 1.0 / (1.0 + (LOUIS*3.0)*RI*PS)
         KM = PS  / (PS  + (LOUIS*2.0)*RI   )
      end where

!===>   DIMENSIONALIZE Kz and  LIMIT DIFFUSIVITY

      ALM = DU*ALM
      ALH = DU*ALH

      KM  = min(KM*ALM, AKHMMAX)
      KH  = min(KH*ALH, AKHMMAX)

      if (associated(KMLS_DIAG)) KMLS_DIAG(:,:,1:LM-1) = KM(:,:,1:LM-1)
      if (associated(KHLS_DIAG)) KHLS_DIAG(:,:,1:LM-1) = KH(:,:,1:LM-1)

  end subroutine LOUIS_KS

   subroutine BELJAARS(IM, JM, LM, DT, &
                       LAMBDA_B, C_B,  &
                       U, V, Z,        &
                       VARFLT, PLE,    &
                       BKV, FKV        )

!BOP
!
!   Orographic drag follows  Beljaars (2003):
!   $$
!   \frac{\partial}{\partial z}\frac{\tau}{\rho} = \frac{C_B}{\lambda_B} |U(z)| U(z) 
!          e^{-\tilde{z}^\frac{3}{2}}\tilde{z}^{-1.2},
!   $$
!   where $z$ is the height above the surface in meters, 
!   $\tilde{z}=\frac{z}{\lambda_B}$, $\tau$ is the orographic stress at $z$,
!   $\rho$ is the air density, $U(z)$ is the wind velocity, and $\lambda_B$ is a vertical length scale.
!   Beljaars uses $\lambda_B = 1500$m, for which the non-dimensional parameter $C_B = 2.5101471 \times 10^{-8}$.
!   These are the default values, but both can be modified from the configuration. To avoid underflow.
!   the tendency is set to zero once $\tilde{z}$ exceeds 4 (i.e., 6 km from the surface for default values). 
!
!EOP

      integer, intent(IN   )                    :: IM,JM,LM
      real,    intent(IN   )                    :: DT
      real,    intent(IN   )                    :: LAMBDA_B
      real,    intent(IN   )                    :: C_B

      real,    intent(IN   ), dimension(:,:,: ) :: U
      real,    intent(IN   ), dimension(:,:,: ) :: V
      real,    intent(IN   ), dimension(:,:,: ) :: Z
      real,    intent(IN   ), dimension(:,:   ) :: VARFLT
      real,    intent(IN   ), dimension(:,:,0:) :: PLE

      real,    intent(INOUT), dimension(:,:,: ) :: BKV

      real,    intent(  OUT), dimension(:,:,: ) :: FKV

      integer :: I,J,L
      real    :: FKV_temp

      do I = 1, IM
         do J = 1, JM
            do L = LM, 1, -1
               FKV(I,J,L) = 0.0

               if (Z(I,J,L) < 4.0*LAMBDA_B) then
                  FKV_temp = Z(I,J,L)*(1.0/LAMBDA_B)
                  FKV_temp = VARFLT(I,J) * exp(-FKV_temp*sqrt(FKV_temp))*(FKV_temp**(-1.2))
                  FKV_temp = (C_B/LAMBDA_B)*min( sqrt(U(I,J,L)**2+V(I,J,L)**2),5.0 )*FKV_temp

                  BKV(I,J,L) = BKV(I,J,L) + DT*FKV_temp
                  FKV(I,J,L) = FKV_temp * (PLE(I,J,L)-PLE(I,J,L-1))
               end if
            end do
         end do 
      end do 

   end subroutine BELJAARS

!*********************************************************************

!BOP

! !IROUTINE:  VTRILU --  Does LU decomposition of tridiagonal matrix.

! !INTERFACE:

   subroutine VTRILU(A,B,C)

! !ARGUMENTS:

      real,               dimension(:,:,:), intent(IN   ) :: C
      real(kind=MAPL_R8), dimension(:,:,:), intent(INOUT) :: A, B

! !DESCRIPTION: {\tt VTRILU} performs an $LU$ decomposition on
! a tridiagonal matrix $M=LU$.
!
! $$
! M = \left( \begin{array}{ccccccc}
!      b_1 & c_1 & & & & & \\
!      a_2 & b_2 & c_2 & & & &  \\
!      &  \cdot& \cdot & \cdot & & &  \\
!      & & \cdot& \cdot & \cdot & &  \\
!      &&  & \cdot& \cdot & \cdot &  \\
!      &&&& a_{K-1} & b_{K-1} & c_{K-1}   \\
!      &&&&& a_{K} & b_{K}
!    \end{array} \right)
! $$
!
!
! $$
! \begin{array}{lr}
! L = \left( \begin{array}{ccccccc}
!      1 &&&&&& \\
!      \hat{a}_2 & 1 & &&&&  \\
!      &  \cdot& \cdot &  & & &  \\
!      & & \cdot& \cdot &  &&  \\
!      &&  & \cdot& \cdot &  &  \\
!      &&&& \hat{a}_{K-1} & 1 &   \\
!      &&&&& \hat{a}_{K} & 1
!    \end{array} \right)
! &
! U = \left( \begin{array}{ccccccc}
!      \hat{b}_1 & c_1 &&&&& \\
!       & \hat{b}_2 & c_2 &&&&  \\
!      &  & \cdot & \cdot & & &  \\
!      & & & \cdot & \cdot &&  \\
!      &&  & & \cdot & \cdot &  \\
!      &&&&  & \hat{b}_{K-1} & c_{K-1}   \\
!      &&&&&  & \hat{b}_{K}
!    \end{array} \right)
! \end{array}
! $$
!
!
! On input, A, B, and C contain, $a_k$, $b_k$, and $c_k$
! the lower, main, and upper diagonals of the matrix, respectively.
! On output, B contains $1/\hat{b}_k$, the inverse of the main diagonal of $U$,
! and A contains $\hat{a}_k$,
! the lower diagonal of $L$. C contains the upper diagonal of the original matrix and of $U$.
!
! The new diagonals $\hat{a}_k$ and $\hat{b}_k$ are:
! $$
! \begin{array}{rcl}
! \hat{b}_1 & = & b_1, \\
! \hat{a}_k & = & \makebox[2 in][l]{$a_k / \hat{b}_{k-1}$,}  k=2, K, \\
! \hat{b}_k & = & \makebox[2 in][l]{$b_k - c_{k-1} \hat{a}_k$,} k=2, K. 
! \end{array}
! $$
!EOP

      integer :: LM, L

      LM = size(C,3)

      B(:,:,1) = 1. / B(:,:,1)

      do L = 2,LM
         A(:,:,L) = A(:,:,L) * B(:,:,L-1)
         B(:,:,L) = 1. / ( B(:,:,L) - C(:,:,L-1) * A(:,:,L) )
      end do

   end subroutine VTRILU

!*********************************************************************

!BOP

! !IROUTINE:  VTRISOLVESURF -- Solves for sensitivity to surface value


! !INTERFACE:

   subroutine VTRISOLVESURF(B,C,Y)

! !ARGUMENTS:

      real,    dimension(:,:,:), intent(IN   ) :: B, C
      real,    dimension(:,:,:), intent(  OUT) :: Y

! !DESCRIPTION: Solves tridiagonal system that has been LU decomposed
!   for the special case
!   where the surface Y (YG) is 1 and the rest of the input Ys are 0.
!   Everything else is as in {\tt VTRISOLVE}. This gives the sensitivity of the
!   solution to a unit change in the surface values.

!EOP

      integer :: LM, L

      LM = size(B,3)

      Y(:,:,LM) = -C(:,:,LM) * B(:,:,LM)

      do L = LM-1,1,-1
         Y(:,:,L) = -C(:,:,L) * Y(:,:,L+1) * B(:,:,L)
      end do

   end subroutine VTRISOLVESURF

!BOP

! !IROUTINE:  VTRISOLVE -- Solves for tridiagonal system that has been decomposed by VTRILU


! !INTERFACE:

  subroutine VTRISOLVE ( A,B,C,Y,YG )

! !ARGUMENTS:

    real,               dimension(:,:,:),  intent(IN   ) ::  A, B, C
    real(kind=MAPL_R8), dimension(:,:,:),  intent(INOUT) ::  Y
    real,               dimension(:,:),    intent(IN)    ::  YG

! !DESCRIPTION: Solves tridiagonal system that has been LU decomposed
!   $LU x = f$. This is done by first solving $L g = f$ for $g$, and 
!   then solving $U x = g$ for $x$. The solutions are:
! $$
! \begin{array}{rcl}
! g_1 & = & f_1, \\
! g_k & = & \makebox[2 in][l]{$f_k - g_{k-1} \hat{a}_{k}$,}  k=2, K, \\
! \end{array}
! $$
! and  
! $$
! \begin{array}{rcl}
! x_K & = & g_K /\hat{b}_K, \\
! x_k & = & \makebox[2 in][l]{($g_k - c_k g_{k+1}) / \hat{b}_{k}$,}  k=K-1, 1 \\
! \end{array}
! $$
!  
!  On input A contains the $\hat{a}_k$, the lower diagonal of $L$,
!   B contains the $1/\hat{b}_k$, inverse of the  main diagonal of $U$,
!   C contains the $c_k$, the upper diagonal of $U$. The forcing, $f_k$ is
!   
!   It returns the
!   solution in the r.h.s input vector, Y. A has the multiplier from the
!   decomposition, B the 
!   matrix (U), and C the upper diagonal of the original matrix and of U.
!   YG is the LM+1 (Ground) value of Y.

!EOP

    integer :: LM, L

    LM = size(Y,3)

! Sweep down, modifying rhs with multiplier A

    do L = 2,LM
       Y(:,:,L) = Y(:,:,L) - Y(:,:,L-1) * A(:,:,L)
    enddo

! Sweep up, solving for updated value. Note B has the inverse of the main diagonal

    if(size(YG)>0) then
       Y(:,:,LM)   = (Y(:,:,LM) - C(:,:,LM) * YG        )*B(:,:,LM)
    else
       Y(:,:,LM)   =  Y(:,:,LM)*B(:,:,LM-1)/(B(:,:,LM-1) - A(:,:,LM)*(1.0+C(:,:,LM-1)*B(:,:,LM-1) ))
    !  Y(:,:,LM)   =  Y(:,:,LM)*B(:,:,LM)/( 1.0+C(:,:,LM)*B(:,:,LM) )     ! Alternate formulation
    endif

    do L = LM-1,1,-1
       Y(:,:,L) = (Y(:,:,L ) - C(:,:,L ) * Y(:,:,L+1))*B(:,:,L )
    enddo

    return
  end subroutine VTRISOLVE

end module GEOS_TurbulenceGridCompMod

