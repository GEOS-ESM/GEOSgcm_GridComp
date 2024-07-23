
!  $Id$

#include "MAPL_Generic.h"

module GEOS_GwdGridCompMod

!BOP

! !MODULE: GEOS_Gwd -- A Module to compute the forcing due to parameterized gravity wave drag

! !DESCRIPTION:
! 
!   {\tt GWD} is a light-weight gridded component to compute the forcing
! due to gravity wave drags. It operates on the ESMF grid that appears in the
! gridded component passed to its {\tt Initialize} method. Unlike
! heavier gridded components, it does not enforce its own grid.
! The only restrictions are that it be a 3-dimensional grid
! in which one dimension is aligned with the vertical coordinate and
! only the horizontal dimensions are decomposed.
!
! The gravity wave drag scheme is based on NCAR WACCM1b gw\_drag routine.
! The scheme includes parameterizations for orographic (stationary) gravity
! waves (Kiehl et al. 1996), and for a spectrum of traveling gravity waves 
!(Sassi et al. 2003; http://acd.ucar.edu/models/WACCM). Both parameteriz-
! ations are based on Lindzen's [1981] formulation. The interested reader 
! is referred to those publications for details of the mathematical
! derivations.
!

! !USES:

  use ESMF
  use MAPL

  use gw_oro, only : gw_oro_init
  use gw_convect, only : gw_beres_init, BeresSourceDesc
  use gw_common, only: GWBand, gw_common_init
  use gw_drag_ncar, only: gw_intr_ncar

  use gw_drag, only: gw_intr
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP
  type(GWBand)          :: beres_band
  type(BeresSourceDesc) :: beres_desc
  type(GWBand)          :: oro_band

contains

!BOP
! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:
  subroutine SetServices ( GC, RC )

! !ARGUMENTS:
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.

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

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_INITIALIZE,  Initialize,  &
                                      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,  Run,  &
                                      RC=STATUS)
    VERIFY_(STATUS)


! Set the state variable specs.
! -----------------------------

!BOS
! !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'PLE',                                       &
        LONG_NAME  = 'air_pressure',                              &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                          &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'T',                                         &
        LONG_NAME  = 'air_temperature',                           &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'Q',                                         &
        LONG_NAME  = 'specific_humidity',                         &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'U',                                         &
        LONG_NAME  = 'eastward_wind',                             &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'V',                                         &
        LONG_NAME  = 'northward_wind',                            &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'SGH',                                       &
        LONG_NAME  = 'standard_deviation_of_topography',          &
        UNITS      = 'm',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'PREF',                                      &
        LONG_NAME  = 'reference_air_pressure',                    &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsVertOnly,                           &
        VLOCATION  = MAPL_VLocationEdge,                          &
        RESTART    = MAPL_RestartSkip,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


! from moist
        call MAPL_AddImportSpec(GC,                              &
             SHORT_NAME='DTDTCN',                                & 
             LONG_NAME ='T tendency due to convection',          &
             UNITS     ='K s-1',                                 &
             DIMS      = MAPL_DimsHorzVert,                      &
             VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
        VERIFY_(STATUS)  
!WMP: Updated this to be the T tendency due to convection...
!JTB: This was moved (3/25/2020) from imports for NCEP GWD, because 
!     new NCAR code will use it for testing of Beres scheme. Not 
!     sure this is what Beres scheme should actually be using, but OK
!     for now until tuning begins. 
!     from this we can compute QMAX (column maximum value)
!     and KTOP, KBOT near the location of QMAX

! !EXPORT STATE:
  
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'PLE',                                       &
        LONG_NAME  = 'air_pressure',                              &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'T',                                         &
        LONG_NAME  = 'air_temperature',                           &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'Q',                                         &
        LONG_NAME  = 'specific_humidity',                         &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'U',                                         &
        LONG_NAME  = 'eastward_wind',                             &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'V',                                         &
        LONG_NAME  = 'northward_wind',                            &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'SGH',                                       &
        LONG_NAME  = 'standard_deviation_of_topography',          &
        UNITS      = 'm',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'PREF',                                      &
        LONG_NAME  = 'reference_air_pressure',                    &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsVertOnly,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DTDT',                                      &
        LONG_NAME  = 'mass_weighted_air_temperature_tendency_due_to_GWD',    &
        UNITS      = 'Pa K s-1',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TTMGW',                                     &
        LONG_NAME  = 'tendency_of_air_temperature_due_to_GWD',       &
        UNITS      = 'K s-1',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DUDT',                                      &
        LONG_NAME  = 'tendency_of_eastward_wind_due_to_GWD',                 &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DVDT',                                      &
        LONG_NAME  = 'tendency_of_northward_wind_due_to_GWD',                &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DTDT_ORO',                                  &
        LONG_NAME  = 'air_temperature_tendency_due_to_orographic_GWD', &
        UNITS      = 'K s-1',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DUDT_ORO',                                  &
        LONG_NAME  = 'tendency_of_eastward_wind_due_to_orographic_GWD',               &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DVDT_ORO',                                  &
        LONG_NAME  = 'tendency_of_northward_wind_due_to_orographic_GWD',              &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DTDT_BKG',                                  &
        LONG_NAME  = 'air_temperature_tendency_due_to_background_GWD', &
        UNITS      = 'K s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DUDT_BKG',                                  &
        LONG_NAME  = 'tendency_of_eastward_wind_due_to_background_GWD',               &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DVDT_BKG',                                  &
        LONG_NAME  = 'tendency_of_northward_wind_due_to_background_GWD',              &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DTDT_RAY',                                  &
        LONG_NAME  = 'air_temperature_tendency_due_to_Rayleigh_friction',        &
        UNITS      = 'K s-1',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DUDT_RAY',                                  &
        LONG_NAME  = 'tendency_of_eastward_wind_due_to_Rayleigh_friction',       &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DVDT_RAY',                                  &
        LONG_NAME  = 'tendency_of_northward_wind_due_to_Rayleigh_friction',      &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUGWX',                                    &
        LONG_NAME  = 'surface_eastward_gravity_wave_stress',      &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUGWY',                                    &
        LONG_NAME  = 'surface_northward_gravity_wave_stress',     &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUOROX',                                   &
        LONG_NAME  = 'surface_eastward_orographic_gravity_wave_stress',      &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUOROY',                                   &
        LONG_NAME  = 'surface_northward_orographic_gravity_wave_stress',     &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUBKGX',                                   &
        LONG_NAME  = 'surface_eastward_background_gravity_wave_stress',      &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUBKGY',                                   &
        LONG_NAME  = 'surface_northward_background_gravity_wave_stress',     &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUMSTX',                                   &
        LONG_NAME  = 'surface_eastward_gravity_wave_stress_due_to_Moist_Processes',      &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUMSTY',                                   &
        LONG_NAME  = 'surface_northward_gravity_wave_stress_due_to_Moist_Processes',     &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'CLDSTD',                                    &
        LONG_NAME  = 'gravity_wave_drag_standard_deviation_due_to_clouds',     &
        UNITS      = 'm',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'UBASE',                                     &
        LONG_NAME  = 'eastward_component_of_base_level_wind',     &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'VBASE',                                     &
        LONG_NAME  = 'northward_component_of_base_level_wind',    &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'UBAR',                                      &
        LONG_NAME  = 'eastward_component_of_mean_level_wind',     &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'VBAR',                                      &
        LONG_NAME  = 'northward_component_of_mean_level_wind',    &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEGWD',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_across_gwd',         &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEORO',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_orographic_gravity_waves',  &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEBKG',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_gravity_wave_background',   &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PERAY',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_Rayleigh_friction',         &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KEGWD',                                                              &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_across_gwd',           &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                                         &
         SHORT_NAME = 'KEORO',                                                                            &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_dissipation_due_to_orographic_gravity_waves', &
         UNITS      = 'W m-2',                                                                            &
         DIMS       = MAPL_DimsHorzOnly,                                                                  &
         VLOCATION  = MAPL_VLocationNone,                                                      RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                                  &
         SHORT_NAME = 'KERAY',                                                                     &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_dissipation_due_to_Rayleigh_friction', &
         UNITS      = 'W m-2',                                                                     &
         DIMS       = MAPL_DimsHorzOnly,                                                           &
         VLOCATION  = MAPL_VLocationNone,                                               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                                        &
         SHORT_NAME = 'KEBKG',                                                                           &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_dissipation_due_to_gravity_wave_background', &
         UNITS      = 'W m-2',                                                                           &
         DIMS       = MAPL_DimsHorzOnly,                                                                 &
         VLOCATION  = MAPL_VLocationNone,                                                     RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                                        &
         SHORT_NAME = 'KERES',                                                                           &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_residual_for_total_energy_conservation',     &
         UNITS      = 'W m-2',                                                                           &
         DIMS       = MAPL_DimsHorzOnly,                                                                 &
         VLOCATION  = MAPL_VLocationNone,                                                     RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                                        &
         SHORT_NAME = 'BKGERR',                                                                          &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_residual_for_BKG_energy_conservation',       &
         UNITS      = 'W m-2',                                                                           &
         DIMS       = MAPL_DimsHorzOnly,                                                                 &
         VLOCATION  = MAPL_VLocationNone,                                                     RC=STATUS  )
     VERIFY_(STATUS)

!EOS

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="DRIVER"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-DRIVER_RUN"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-INTR"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-DRIVER_DATA"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--DRIVER_DATA_DEVICE"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--DRIVER_DATA_CONST"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-DRIVER_ALLOC"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-DRIVER_DEALLOC"   ,RC=STATUS)
    VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( gc, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !BOP

  ! !IROUTINE: Initialize -- Initialize method for the composite Moist Gridded Component

  ! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

    ! !DESCRIPTION: The Initialize method of the GWD Physics Gridded Component first 
    !   calls the Initialize method of the children.  Then, if using the NCAR GWD
    !   scheme, calls the initialization routines.

    !EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),      pointer  :: MAPL

! NCAR GWD variables

    character(len=ESMF_MAXPATHLEN) :: BERES_FILE_NAME
    character(len=ESMF_MAXSTR)     :: ERRstring
    logical                        :: USE_NCAR_GWD

!=============================================================================

   ! Begin...

   ! Get my name and set-up traceback handle
   ! ---------------------------------------

      Iam = 'Initialize'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // Iam

      ! Get my internal MAPL_Generic state
      !-----------------------------------

      call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
      VERIFY_(STATUS)

      ! Call Generic Initialize for GWD GC
      !-----------------------------------

      call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK, RC=STATUS )
      VERIFY_(STATUS)

      ! Check to see if we are using NCAR GWD
      !--------------------------------------

      call MAPL_GetResource( MAPL, USE_NCAR_GWD, Label="USE_NCAR_GWD:",  default=.false., RC=STATUS)
      VERIFY_(STATUS)

      !++jtb 03/2020
      !-----------------------------------
      if (USE_NCAR_GWD) then
         call gw_common_init( .FALSE. , 1 , & 
                              1.0_MAPL_R8 * MAPL_GRAV , &
                              1.0_MAPL_R8 * MAPL_RGAS , &
                              1.0_MAPL_R8 * MAPL_CP , &
                              0.50_MAPL_R8 , 0.25_MAPL_R8, ERRstring )

         ! Beres Scheme File
         call MAPL_GetResource( MAPL, BERES_FILE_NAME, Label="BERES_FILE_NAME:", RC=STATUS)
         VERIFY_(STATUS)

         call gw_beres_init( BERES_FILE_NAME , beres_band, beres_desc )

         call gw_oro_init ( oro_band )
      end if

      ! All done
      !---------

      RETURN_(ESMF_SUCCESS)
   end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: RUN -- Run method for the GWD component

! !INTERFACE:
subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:
  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),     pointer   :: MAPL
  type (ESMF_Alarm       )            :: ALARM

  integer                             :: IM, JM, LM
  integer                             :: pgwv
  real                                :: effgworo, effgwbkg
  real                                :: CDMBGWD1, CDMBGWD2
  real                                :: bgstressmax
  real, pointer, dimension(:,:)       :: LATS

! Rayleigh friction parameters

  REAL                                :: H0, HH, Z1, TAU1

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

   Iam = "Run"
   call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
!----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)

! Local aliases to the state, grid, and configuration
! ---------------------------------------------------

   call MAPL_TimerOn(MAPL,"TOTAL")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL, &
         IM=IM, JM=JM, LM=LM,        &
         RUNALARM=ALARM, LATS=LATS,  &
                           RC=STATUS )
    VERIFY_(STATUS)

! Gravity wave drag
! -----------------

    call MAPL_GetResource( MAPL, effgworo, Label="EFFGWORO:", default=0.250, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, effgwbkg, Label="EFFGWBKG:", default=0.125, RC=STATUS)
    VERIFY_(STATUS)

    if( LM .eq. 72 ) then
        call MAPL_GetResource( MAPL, pgwv,        Label="PGWV:",        default=4,    RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetResource( MAPL, bgstressmax, Label="BGSTRESSMAX:", default=0.9,  RC=STATUS)
        VERIFY_(STATUS)
     else
        call MAPL_GetResource( MAPL, pgwv,        Label="PGWV:",        default=NINT(4*LM/72.0),    RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetResource( MAPL, bgstressmax, Label="BGSTRESSMAX:", default=0.9, RC=STATUS)
        VERIFY_(STATUS)
     endif

! Rayleigh friction
! -----------------
    CALL MAPL_GetResource( MAPL, Z1,   Label="RAYLEIGH_Z1:",   default=75000.,  RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource( MAPL, TAU1, Label="RAYLEIGH_TAU1:", default=172800., RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource( MAPL, H0,   Label="RAYLEIGH_H0:",   default=7000.,	RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_GetResource( MAPL, HH,   Label="RAYLEIGH_HH:",   default=7500.,	RC=STATUS)
    VERIFY_(STATUS)

! If its time, recalculate the GWD tendency
! -----------------------------------------

   if ( ESMF_AlarmIsRinging( ALARM ) ) then
      call ESMF_AlarmRingerOff(ALARM, RC=STATUS); VERIFY_(STATUS)
      call MAPL_TimerOn (MAPL,"DRIVER")
      call Gwd_Driver(RC=STATUS); VERIFY_(STATUS)
      call MAPL_TimerOff(MAPL,"DRIVER")
   endif

   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)

   contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Gwd_Driver(RC)
      integer, optional, intent(OUT) :: RC

!  Locals

      character(len=ESMF_MAXSTR)      :: IAm
      integer                         :: STATUS

      type (ESMF_TimeInterval)        :: TINT

!  Pointers from Import state

      real, pointer, dimension(:)      :: PREF
      real, pointer, dimension(:,:)    :: SGH
      real, pointer, dimension(:,:,:)  :: PLE, T, Q, U, V
      !++jtb Array for moist/deepconv heating
      real, pointer, dimension(:,:,:)  :: HT_dpc

!  Pointers to Export state

      real, pointer, dimension(:)      :: PREF_EXP
      real, pointer, dimension(:,:)    :: SGH_EXP
      real, pointer, dimension(:,:,:)  :: PLE_EXP, T_EXP, Q_EXP, U_EXP, V_EXP

      real, pointer, dimension(:,:)    :: CLDSTD
      real, pointer, dimension(:,:)    :: UBAR,    VBAR
      real, pointer, dimension(:,:)    :: UBASE,   VBASE
      real, pointer, dimension(:,:)    :: TAUGWX,  TAUGWY
      real, pointer, dimension(:,:)    :: TAUOROX, TAUOROY
      real, pointer, dimension(:,:)    :: TAUBKGX, TAUBKGY
      real, pointer, dimension(:,:,:)  :: TAUOROXZ,TAUOROYZ,FEOROZ,FEPOROZ
      real, pointer, dimension(:,:,:)  :: TAUBKGXZ,TAUBKGYZ,FEBKGZ,FEPBKGZ
      real, pointer, dimension(:,:)    :: TAUOROXT,TAUOROYT,FEOROT,FEPOROT
      real, pointer, dimension(:,:)    :: TAUOROXS,TAUOROYS,FEOROS,FEPOROS
      real, pointer, dimension(:,:)    :: TAUBKGXT,TAUBKGYT,FEBKGT,FEPBKGT
      real, pointer, dimension(:,:)    :: TAUBKGXS,TAUBKGYS,FEBKGS,FEPBKGS
      real, pointer, dimension(:,:)    :: TAUMSTX, TAUMSTY
      real, pointer, dimension(:,:)    :: KEGWD, KEORO,  KERAY,  KEBKG, KERES
      real, pointer, dimension(:,:)    :: PEGWD, PEORO,  PERAY,  PEBKG, BKGERR

      real, pointer, dimension(:,:,:)  :: DTDT, DUDT, DVDT, TTMGW
      real, pointer, dimension(:,:,:)  :: DTDT_ORO, DUDT_ORO, DVDT_ORO
      real, pointer, dimension(:,:,:)  :: DTDT_BKG, DUDT_BKG, DVDT_BKG
      real, pointer, dimension(:,:,:)  :: DTDT_RAY, DUDT_RAY, DVDT_RAY
      real, pointer, dimension(:,:,:)  :: DTGENBKG, DUGENBKG, DVGENBKG
      
! local variables

      real,              dimension(IM,JM,LM  ) :: ZM, PMID, PDEL, RPDEL, PMLN
      real,              dimension(IM,JM,LM  ) :: DUDT_ORG, DVDT_ORG, DTDT_ORG
      real,              dimension(IM,JM,LM  ) :: DUDT_GWD, DVDT_GWD, DTDT_GWD
      real,              dimension(IM,JM,LM  ) :: DUDT_RAH, DVDT_RAH, DTDT_RAH
      real,              dimension(IM,JM,LM  ) :: DUDT_TOT, DVDT_TOT, DTDT_TOT
      real,              dimension(IM,JM,LM+1) :: PILN,   ZI
      real,              dimension(      LM  ) :: ZREF, KRAY
      real,              dimension(IM,JM     ) :: TAUXO_TMP, TAUYO_TMP
      real,              dimension(IM,JM     ) :: TAUXB_TMP, TAUYB_TMP
      real,              dimension(IM,JM,LM+1) :: TAUXO_3D , TAUYO_3D , FEO_3D, FEPO_3D
      real,              dimension(IM,JM,LM+1) :: TAUXB_3D , TAUYB_3D , FEB_3D, FEPB_3D
      real,              dimension(IM,JM,LM  ) :: DUBKGSRC , DVBKGSRC , DTBKGSRC
      real,              dimension(IM,JM)      :: KEGWD_X, KEORO_X,  KERAY_X,  KEBKG_X, KERES_X
      real,              dimension(IM,JM)      :: PEGWD_X, PEORO_X,  PERAY_X,  PEBKG_X, BKGERR_X

      integer                                  :: J, K, L
      real(ESMF_KIND_R8)                       :: DT_R8
      real                                     :: DT     ! time interval in sec

! NCAR GWD vars

      logical :: USE_NCAR_GWD

!  Begin...
!----------

      IAm = "Gwd_Driver"

! Get time step
!-------------------------------------------------

      call ESMF_AlarmGet( ALARM, ringInterval=TINT,    RC=STATUS); VERIFY_(STATUS)
      call ESMF_TimeIntervalGet(TINT, S_R8=DT_R8,      RC=STATUS); VERIFY_(STATUS)

      DT = DT_R8

! Pointers to inputs
!---------------------

      call MAPL_GetPointer( IMPORT, PLE,    'PLE',     RC=STATUS ); VERIFY_(STATUS)
      call MAPL_GetPointer( IMPORT, T,      'T',       RC=STATUS ); VERIFY_(STATUS)
      call MAPL_GetPointer( IMPORT, Q,      'Q',       RC=STATUS ); VERIFY_(STATUS)
      call MAPL_GetPointer( IMPORT, U,      'U',       RC=STATUS ); VERIFY_(STATUS)
      call MAPL_GetPointer( IMPORT, V,      'V',       RC=STATUS ); VERIFY_(STATUS)
      call MAPL_GetPointer( IMPORT, SGH,    'SGH',     RC=STATUS ); VERIFY_(STATUS)
      call MAPL_GetPointer( IMPORT, PREF,   'PREF',    RC=STATUS ); VERIFY_(STATUS)
!++jtb
      call MAPL_GetPointer( IMPORT, HT_dpc, 'DTDTCN',  RC=STATUS ); VERIFY_(STATUS)

! Allocate/refer to the outputs
!------------------------------

      call MAPL_GetPointer(EXPORT,  PLE_EXP, 'PLE'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    T_EXP, 'T'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    Q_EXP, 'Q'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    U_EXP, 'U'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    V_EXP, 'V'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,  SGH_EXP, 'SGH'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PREF_EXP, 'PREF'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    TTMGW, 'TTMGW'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DTDT_ORO, 'DTDT_ORO', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DUDT_ORO, 'DUDT_ORO', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DVDT_ORO, 'DVDT_ORO', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DTDT_BKG, 'DTDT_BKG', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DUDT_BKG, 'DUDT_BKG', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DVDT_BKG, 'DVDT_BKG', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DTDT_RAY, 'DTDT_RAY', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DUDT_RAY, 'DUDT_RAY', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DVDT_RAY, 'DVDT_RAY', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,   TAUGWX, 'TAUGWX'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,   TAUGWY, 'TAUGWY'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,  TAUOROX, 'TAUOROX' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,  TAUOROY, 'TAUOROY' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,  TAUBKGX, 'TAUBKGX' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,  TAUBKGY, 'TAUBKGY' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,  TAUMSTX, 'TAUMSTX' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,  TAUMSTY, 'TAUMSTY' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    UBASE, 'UBASE'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    VBASE, 'VBASE'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,     UBAR, 'UBAR'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,     VBAR, 'VBAR'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,   CLDSTD, 'CLDSTD'  , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT,     DTDT, 'DTDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,     DUDT, 'DUDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,     DVDT, 'DVDT'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT,    PEGWD, 'PEGWD'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    PEORO, 'PEORO'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    PERAY, 'PERAY'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    PEBKG, 'PEBKG'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT,    KEGWD, 'KEGWD'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    KEORO, 'KEORO'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    KERAY, 'KERAY'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    KEBKG, 'KEBKG'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,    KERES, 'KERES'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,   BKGERR, 'BKGERR'  , RC=STATUS); VERIFY_(STATUS)


      CALL PREGEO(IM*JM,   LM,   &
                    PLE, LATS,   PMID,  PDEL, RPDEL,     PILN,     PMLN)

! Compute ZM
!-------------

      call GEOPOTENTIAL( IM*JM, LM,                  &
           PILN,  PMLN,  PLE,   PMID, PDEL, RPDEL,   &
           T,     Q,     ZI,    ZM                   )

! Do gravity wave drag calculations on a list of soundings
!---------------------------------------------------------

    call MAPL_GetResource(MAPL,USE_NCAR_GWD,'USE_NCAR_GWD:', default=.false., RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"-INTR")

    if (USE_NCAR_GWD) then
       ! Use Julio new code
    call gw_intr_ncar(IM*JM,    LM,         DT,                  &
         PGWV,      beres_desc, beres_band, oro_band,            &
         PLE,       T,          U,          V,      HT_dpc,      &
         SGH,       PREF,                                        &
         PMID,      PDEL,       RPDEL,      PILN,   ZM,    LATS, &
         DUDT_GWD,  DVDT_GWD,   DTDT_GWD,                        &
         DUDT_ORG,  DVDT_ORG,   DTDT_ORG,                        &
         TAUXO_TMP, TAUYO_TMP,  TAUXO_3D,   TAUYO_3D,  FEO_3D,   &
         TAUXB_TMP, TAUYB_TMP,  TAUXB_3D,   TAUYB_3D,  FEB_3D,   &
         FEPO_3D,   FEPB_3D,    DUBKGSRC,   DVBKGSRC,  DTBKGSRC, &
         BGSTRESSMAX, effgworo, effgwbkg,   RC=STATUS            )
    VERIFY_(STATUS)
    else
       ! Use GEOS GWD    
       call gw_intr   (IM*JM,      LM,         DT,                  &
            PGWV,                                                   &
            PLE,       T,          U,          V,      SGH,   PREF, &
            PMID,      PDEL,       RPDEL,      PILN,   ZM,    LATS, &
            DUDT_GWD,  DVDT_GWD,   DTDT_GWD,                        &
            DUDT_ORG,  DVDT_ORG,   DTDT_ORG,                        &
            TAUXO_TMP, TAUYO_TMP,  TAUXO_3D,   TAUYO_3D,  FEO_3D,   &
            TAUXB_TMP, TAUYB_TMP,  TAUXB_3D,   TAUYB_3D,  FEB_3D,   &
            FEPO_3D,   FEPB_3D,    DUBKGSRC,   DVBKGSRC,  DTBKGSRC, &
            BGSTRESSMAX, effgworo, effgwbkg,   RC=STATUS            )
      VERIFY_(STATUS)
    end if

    call MAPL_TimerOff(MAPL,"-INTR")

    CALL POSTINTR(IM*JM, LM, DT, H0, HH, Z1, TAU1, &
          PREF,     &
          PDEL,     &
          U,        &
          V,        &
          DUDT_GWD, &
          DVDT_GWD, &
          DTDT_GWD, &
          DUDT_ORG, &
          DVDT_ORG, &
          DTDT_ORG, &

          DUDT_TOT, &
          DVDT_TOT, &
          DTDT_TOT, &
          DUDT_RAH, &
          DVDT_RAH, &
          DTDT_RAH, &
          PEGWD_X,  &
          PEORO_X,  &
          PERAY_X,  &
          PEBKG_X,  &
          KEGWD_X,  &
          KEORO_X,  &
          KERAY_X,  &
          KEBKG_X,  &
          KERES_X,  &
          BKGERR_X  )

!! Tendency diagnostics
!!---------------------

    if(associated(DUDT    )) DUDT     = DUDT_TOT
    if(associated(DVDT    )) DVDT     = DVDT_TOT
    if(associated(DTDT    )) DTDT     = DTDT_TOT

    if(associated(DUDT_RAY)) DUDT_RAY = DUDT_RAH
    if(associated(DVDT_RAY)) DVDT_RAY = DVDT_RAH
    if(associated(DTDT_RAY)) DTDT_RAY = DTDT_RAH

!! KE dIagnostics
!!----------------

    if(associated(PEGWD   )) PEGWD  = PEGWD_X
    if(associated(PEORO   )) PEORO  = PEORO_X
    if(associated(PERAY   )) PERAY  = PERAY_X
    if(associated(PEBKG   )) PEBKG  = PEBKG_X
    if(associated(KEGWD   )) KEGWD  = KEGWD_X
    if(associated(KEORO   )) KEORO  = KEORO_X
    if(associated(KERAY   )) KERAY  = KERAY_X
    if(associated(KEBKG   )) KEBKG  = KEBKG_X
    if(associated(KERES   )) KERES  = KERES_X
    if(associated(BKGERR  )) BKGERR = BKGERR_X

!! Tendency diagnostics
!!---------------------

    if(associated(DUDT_ORO)) DUDT_ORO = DUDT_ORG
    if(associated(DVDT_ORO)) DVDT_ORO = DVDT_ORG
    if(associated(DTDT_ORO)) DTDT_ORO = DTDT_ORG

    if(associated(DUDT_BKG)) DUDT_BKG = DUDT_GWD - DUDT_ORG
    if(associated(DVDT_BKG)) DVDT_BKG = DVDT_GWD - DVDT_ORG
    if(associated(DTDT_BKG)) DTDT_BKG = DTDT_GWD - DTDT_ORG

! Orographic stress
!------------------

    if(associated(TAUGWX  )) TAUGWX  = TAUXO_TMP + TAUXB_TMP
    if(associated(TAUGWY  )) TAUGWY  = TAUYO_TMP + TAUYB_TMP
    if(associated(TAUOROX )) TAUOROX = TAUXO_TMP
    if(associated(TAUOROY )) TAUOROY = TAUYO_TMP
    if(associated(TAUBKGX )) TAUBKGX = TAUXB_TMP
    if(associated(TAUBKGY )) TAUBKGY = TAUYB_TMP

! Export unweighted T Tendency
!-----------------------------

    if(associated(TTMGW)) then
       if(associated(DTDT )) then
          TTMGW = DTDT
       else
          TTMGW = 0.0
       end if
    end if

! AMM modify T_EXP to be the T AFTER GWD, ie., add the tendency*dt
! (need to do this before DTDT is pressure weighted for the dynamics)
    if(associated(T_EXP   )) T_EXP    = T + DTDT*DT

! DTDT has to be pressure weighted and is all due to frictional heating.
!-----------------------------------------------------------------------

    if(associated(DTDT    )) then
       DTDT = DTDT*PDEL 
    end if

    if(associated(PREF_EXP)) PREF_EXP = PREF
    if(associated(SGH_EXP )) SGH_EXP  = SGH
    if(associated(PLE_EXP )) PLE_EXP  = PLE
    if(associated(Q_EXP   )) Q_EXP    = Q
    if(associated(U_EXP   )) U_EXP    = U
    if(associated(V_EXP   )) V_EXP    = V


! All done
!-----------

    RETURN_(ESMF_SUCCESS)
   end subroutine GWD_DRIVER

  end subroutine RUN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine geopotential(pcols  , pver   ,                   &
         piln   , pmln   , pint  , pmid   , pdel   , rpdel  , &
         t      , q      , zi     , zm     )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the geopotential height (above the surface) at the midpoints and 
! interfaces using the input temperatures and pressures.
! Author: B.Boville, Feb 2001 from earlier code by Boville and S.J. Lin
!
!-----------------------------------------------------------------------

    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: pcols                ! Number of longitudes
    integer, intent(in) :: pver                 ! Number of vertical layers

    real,    intent(in) :: piln (pcols,pver+1)  ! Log interface pressures
    real,    intent(in) :: pmln (pcols,pver)    ! Log midpoint pressures
    real,    intent(in) :: pint (pcols,pver+1)  ! Interface pressures
    real,    intent(in) :: pmid (pcols,pver)    ! Midpoint pressures
    real,    intent(in) :: pdel (pcols,pver)    ! layer thickness
    real,    intent(in) :: rpdel(pcols,pver)    ! inverse of layer thickness
    real,    intent(in) :: t    (pcols,pver)    ! temperature
    real,    intent(in) :: q    (pcols,pver)    ! specific humidity

! Output arguments

    real,    intent(out) :: zi(pcols,pver+1)    ! Height above surface at interfaces
    real,    intent(out) :: zm(pcols,pver)      ! Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
!
    logical  :: fvdyn              ! finite volume dynamics
    integer  :: i,k                ! Lon, level indices
    real     :: hkk                ! diagonal element of hydrostatic matrix
    real     :: hkl                ! off-diagonal element
    real     :: tv                 ! virtual temperature
    real     :: tvfac              ! Tv/T

    real, parameter :: ROG     = MAPL_RGAS/MAPL_GRAV
!
!-----------------------------------------------------------------------
!

! Set dynamics flag

    fvdyn = .true.

! The surface height is zero by definition.

    I_LOOP: do i = 1, pcols

       zi(i,pver+1) = 0.0

! Compute zi, zm from bottom up. 
! Note, zi(i,k) is the interface above zm(i,k)

       do k = pver, 1, -1

! First set hydrostatic elements consistent with dynamics

          if (fvdyn) then
             hkl = piln(i,k+1) - piln(i,k)
             hkk = piln(i,k+1) - pmln(i,k)
          else
             hkl = pdel(i,k) / pmid(i,k)
             hkk = 0.5 * hkl
          end if

! Now compute tv, zm, zi

          tvfac   = 1. + MAPL_VIREPS * q(i,k)
          tv      = t(i,k) * tvfac

          zm(i,k) = zi(i,k+1) + ROG * tv * hkk
          zi(i,k) = zi(i,k+1) + ROG * tv * hkl
       end do
    end do I_LOOP

    return
  end subroutine geopotential

!----------------------------------------------------------------------- 

  subroutine pregeo(pcols,pver,&
    ple,lats,pmid,pdel,rpdel,piln,pmln)

    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!

    integer, intent(in) :: pcols    ! Number of longitudes
    integer, intent(in) :: pver     ! Number of vertical layers

    real,    intent(in) :: ple (pcols,pver+1)    ! Interface pressures
    real,    intent(in) :: lats(pcols)           ! latitude in radian

! Output arguments
    
    real,    intent(out) :: pmid  (pcols,pver)   ! Midpoint pressures
    real,    intent(out) :: pdel  (pcols,pver)   ! layer thickness
    real,    intent(out) :: rpdel (pcols,pver)   ! inverse of layer thickness
    real,    intent(out) :: piln  (pcols,pver+1) ! Log interface pressures
    real,    intent(out) :: pmln  (pcols,pver)   ! Log midpoint pressures

!
!---------------------------Local variables-----------------------------
!
    integer :: i,k

    real    :: hvsd  ! Efficiency factor

    real, parameter :: PI_GWD  = 4.0*atan(1.0)  ! This is *not* MAPL_PI

!
!-----------------------------------------------------------------------
!

! Form pressure factors
!----------------------

    I_LOOP: DO I = 1, PCOLS

       DO K = 1, PVER
           PMID(I,K) = 0.5*(  PLE(I,K  ) + PLE(I,K+1) )
           PDEL(I,K) =        PLE(I,K+1) - PLE(I,K  )
          RPDEL(I,K) = 1.0 / PDEL(I,K)
           PILN(I,K) = log(   PLE(I,K) )
           PMLN(I,K) = log(  PMID(I,K) ) !
       END DO
       PILN(I,PVER+1)  = log( PLE(I,PVER+1)  )
    END DO I_LOOP

  end subroutine pregeo

  subroutine postintr(pcols,pver,dt, h0, hh, z1, tau1, &
        pref, &
        pdel, &
        u, &
        v, &
        dudt_gwd, &
        dvdt_gwd, &
        dtdt_gwd, &
        dudt_org, &
        dvdt_org, &
        dtdt_org, &

        ! Outputs
        dudt_tot, &
        dvdt_tot, &
        dtdt_tot, &
        dudt_rah, &
        dvdt_rah, &
        dtdt_rah, &
        pegwd, &
        peoro, &
        peray, &
        pebkg, &
        kegwd, &
        keoro, &
        keray, &
        kebkg, &
        keres, &
        bkgerr )
    
    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!

    integer, intent(in) :: PCOLS ! Number of longitudes
    integer, intent(in) :: PVER  ! Number of vertical layers
    real,    intent(in) :: DT    ! Time step
    real,    intent(in) :: H0, HH, Z1, TAU1 ! Rayleigh friction parameters

    real,    intent(in) :: PREF(PVER+1)
    real,    intent(in) :: PDEL(PCOLS,PVER)
    real,    intent(in) :: U(PCOLS,PVER)
    real,    intent(in) :: V(PCOLS,PVER)

    real,    intent(in) :: DUDT_GWD(PCOLS,PVER)
    real,    intent(in) :: DVDT_GWD(PCOLS,PVER)
    real,    intent(in) :: DTDT_GWD(PCOLS,PVER)
    real,    intent(in) :: DUDT_ORG(PCOLS,PVER)
    real,    intent(in) :: DVDT_ORG(PCOLS,PVER)
    real,    intent(in) :: DTDT_ORG(PCOLS,PVER)

    real,    intent(out) :: DUDT_TOT(PCOLS,PVER)
    real,    intent(out) :: DVDT_TOT(PCOLS,PVER)
    real,    intent(out) :: DTDT_TOT(PCOLS,PVER)
    real,    intent(out) :: DUDT_RAH(PCOLS,PVER)
    real,    intent(out) :: DVDT_RAH(PCOLS,PVER)
    real,    intent(out) :: DTDT_RAH(PCOLS,PVER)
    real,    intent(out) :: PEGWD(PCOLS)
    real,    intent(out) :: PEORO(PCOLS)
    real,    intent(out) :: PERAY(PCOLS)
    real,    intent(out) :: PEBKG(PCOLS)
    real,    intent(out) :: KEGWD(PCOLS)
    real,    intent(out) :: KEORO(PCOLS)
    real,    intent(out) :: KERAY(PCOLS)
    real,    intent(out) :: KEBKG(PCOLS)
    real,    intent(out) :: KERES(PCOLS)
    real,    intent(out) :: BKGERR(PCOLS)

!
!---------------------------Local variables-----------------------------
!
    integer :: i,k
    real :: zref, kray
!
!-----------------------------------------------------------------------
!

    I_LOOP: DO I = 1, PCOLS

       PEGWD(I)  = 0.0
       PEORO(I)  = 0.0
       PERAY(I)  = 0.0
       PEBKG(I)  = 0.0
       KEGWD(I)  = 0.0
       KEORO(I)  = 0.0
       KERAY(I)  = 0.0
       KEBKG(I)  = 0.0
       KERES(I)  = 0.0
       BKGERR(I) = 0.0

       DO K = 1, PVER 

! Rayleigh friction
!------------------

          ZREF     = H0 * LOG(MAPL_P00/(0.5*(PREF(K)+PREF(K+1))))
          KRAY     = (1.0/TAU1)*( 1.0 - TANH( (Z1-ZREF)/HH ) )
          KRAY     = KRAY/(1+DT*KRAY)

          DUDT_RAH(I,K) = -U(I,K)*KRAY
          DVDT_RAH(I,K) = -V(I,K)*KRAY

          DTDT_RAH(I,K) = - ((U(I,K) + (0.5*DT)*DUDT_RAH(I,K))*DUDT_RAH(I,K) + &
                             (V(I,K) + (0.5*DT)*DVDT_RAH(I,K))*DVDT_RAH(I,K)   ) * (1.0/MAPL_CP)

          DUDT_TOT(I,K) = DUDT_RAH(I,K) + DUDT_GWD(I,K)
          DVDT_TOT(I,K) = DVDT_RAH(I,K) + DVDT_GWD(I,K)
          DTDT_TOT(I,K) = DTDT_RAH(I,K) + DTDT_GWD(I,K)

! KE dIagnostics
!----------------

          PEGWD(I) = PEGWD(I) +  DTDT_TOT(I,K)               *PDEL(I,K)*(MAPL_CP/MAPL_GRAV)
          PEORO(I) = PEORO(I) +  DTDT_ORG(I,K)               *PDEL(I,K)*(MAPL_CP/MAPL_GRAV)
          PERAY(I) = PERAY(I) +  DTDT_RAH(I,K)               *PDEL(I,K)*(MAPL_CP/MAPL_GRAV)
          PEBKG(I) = PEBKG(I) + (DTDT_GWD(I,K)-DTDT_ORG(I,K))*PDEL(I,K)*(MAPL_CP/MAPL_GRAV)

          KEGWD(I) = KEGWD(I) + ((U(I,K)+(0.5*DT)*DUDT_TOT(I,K))*DUDT_TOT(I,K) +   &
                                 (V(I,K)+(0.5*DT)*DVDT_TOT(I,K))*DVDT_TOT(I,K) ) * PDEL(I,K)*(1.0/MAPL_GRAV)

          KEORO(I) = KEORO(I) + ((U(I,K)+(0.5*DT)*DUDT_ORG(I,K))*DUDT_ORG(I,K) +   &
                                 (V(I,K)+(0.5*DT)*DVDT_ORG(I,K))*DVDT_ORG(I,K) ) * PDEL(I,K)*(1.0/MAPL_GRAV)

          KERAY(I) = KERAY(I) + ((U(I,K)+(0.5*DT)*DUDT_RAH(I,K))*DUDT_RAH(I,K) +   &
                                 (V(I,K)+(0.5*DT)*DVDT_RAH(I,K))*DVDT_RAH(I,K) ) * PDEL(I,K)*(1.0/MAPL_GRAV)

          KEBKG(I) = KEBKG(I) + ((U(I,K)+(0.5*DT)*(DUDT_GWD(I,K) - DUDT_ORG(I,K)))*(DUDT_GWD(I,K) - DUDT_ORG(I,K)) +     &
                                 (V(I,K)+(0.5*DT)*(DVDT_GWD(I,K) - DVDT_ORG(I,K)))*(DVDT_GWD(I,K) - DVDT_ORG(I,K))   ) * &
                                  PDEL(I,K)*(1.0/MAPL_GRAV)
       END DO

       BKGERR(I) = -( PEBKG(I) + KEBKG(I) )
       KERES(I)  =    PEGWD(I) + KEGWD(I) + BKGERR(I)

    END DO I_LOOP

  end subroutine postintr

end module GEOS_GwdGridCompMod
