
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

#ifdef _CUDA
  use gw_drag, only: &
        ! Subroutines
        GW_INTR, &
        ! Working Arrays
        ZM_DEV, LNPINT_DEV, PMLN_DEV, PMID_DEV, &
        RPDEL_DEV, &
        ! Inputs - PREGEO
        PINT_DEV, RLAT_DEV, &
        ! Outputs - PREGEO
        PDEL_DEV, &
        ! Inputs - GEOPOT
        T_DEV, Q_DEV, &
        ! Outputs - GEOPOT
        ZI_DEV, &
        ! Inputs - INTR
        U_DEV, V_DEV, SGH_DEV, PREF_DEV, &
        ! Outputs - INTR
        DUDT_GWD_DEV, DVDT_GWD_DEV, DTDT_GWD_DEV, &
        DUDT_ORG_DEV, DVDT_ORG_DEV, DTDT_ORG_DEV, &
        TAUGWDX_DEV, TAUGWDY_DEV, TAUOX_DEV, TAUOY_DEV, &
        FEO_DEV, FEPO_DEV,  TAUBKGX_DEV, TAUBKGY_DEV, &
        TAUBX_DEV, TAUBY_DEV,  FEB_DEV, FEPB_DEV, &
        UTBSRC_DEV, VTBSRC_DEV, TTBSRC_DEV, &
        ! Outputs - POSTINTR
        DUDT_TOT_DEV, DVDT_TOT_DEV, DTDT_TOT_DEV, &
        DUDT_RAH_DEV, DVDT_RAH_DEV, DTDT_RAH_DEV, &
        PEGWD_DEV, PEORO_DEV, PERAY_DEV, PEBKG_DEV, &
        KEGWD_DEV, KEORO_DEV, KERAY_DEV, KEBKG_DEV, &
        KERES_DEV, BKGERR_DEV
  use cudafor
#else
  use gw_drag, only: gw_intr
#endif
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP
  logical, parameter :: USE_NCEP_GWD = .false.

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
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'T',                                         &
        LONG_NAME  = 'air_temperature',                           &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'Q',                                         &
        LONG_NAME  = 'specific_humidity',                         &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'U',                                         &
        LONG_NAME  = 'eastward_wind',                             &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'V',                                         &
        LONG_NAME  = 'northward_wind',                            &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'SGH',                                       &
        LONG_NAME  = 'standard_deviation_of_topography',          &
        UNITS      = 'm',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'PREF',                                      &
        LONG_NAME  = 'reference_air_pressure',                    &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsVertOnly,                           &
        VLOCATION  = MAPL_VLocationEdge,                          &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

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
        LONG_NAME  = 'air_temperature_tendency_due_to_GWD',       &
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

     if (USE_NCEP_GWD) then
!ALT: Reminder for myself: we need connections in Physics
! We need some new imports
! from turbulance
        call MAPL_AddImportSpec(GC,                             &
             SHORT_NAME = 'KPBL',                               &
             LONG_NAME  = 'planetary_boundary_layer_level',     &
             UNITS      = '1',                                  &
             DIMS       = MAPL_DimsHorzOnly,                    &
             VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
        VERIFY_(STATUS)      

! from moist
        call MAPL_AddImportSpec(GC,                              &
             SHORT_NAME='DTDT_moist',                            & 
             LONG_NAME ='T tendency due to moist',               &
             UNITS     ='K s-1',                                 &
             DIMS      = MAPL_DimsHorzVert,                      &
             VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
        VERIFY_(STATUS)  
!ALT: from this we can compute QMAX (column maximum value)
!     and KTOP, KBOT near the location of QMAX

        call MAPL_AddImportSpec(GC,                              &
             SHORT_NAME='CNV_FRC',                               &
             LONG_NAME ='convective_fraction',                   &
             UNITS     ='1',                                     &
             DIMS      = MAPL_DimsHorzOnly,                      &
             VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
        VERIFY_(STATUS)

! from dycore
        call MAPL_AddImportSpec ( gc,                            &
             SHORT_NAME = 'DXC',                                 &
             LONG_NAME  = 'cgrid_delta_x',                       &
             UNITS      = 'm'  ,                                 &
             DIMS       = MAPL_DimsHorzOnly,                     &
             VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
        VERIFY_(STATUS)

!        call MAPL_AddImportSpec ( gc,                            &
!             SHORT_NAME = 'DYC',                                 &
!             LONG_NAME  = 'cgrid_delta_y',                       &
!             UNITS      = 'm'  ,                                 &
!             DIMS       = MAPL_DimsHorzOnly,                     &
!             VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
!        VERIFY_(STATUS)


! New internal state for boundary data
        call MAPL_AddInternalSpec ( gc,                       &
             SHORT_NAME = 'STDDEV',                           &
             LONG_NAME  = 'orographic standard deviation',    &
             UNITS      = 'm',                                &
             DIMS       = MAPL_DimsHorzOnly,                  &
             VLOCATION  = MAPL_VLocationNone,                 &
             RC=STATUS  )
        VERIFY_(STATUS)
        
        call MAPL_AddInternalSpec ( gc,                       &
             SHORT_NAME = 'CONVEXITY',                        &
             LONG_NAME  = 'orographic convexity',             &
             UNITS      = '1',                                &
             DIMS       = MAPL_DimsHorzOnly,                  &
             VLOCATION  = MAPL_VLocationNone,                 &
             RC=STATUS  )
        VERIFY_(STATUS)
        
        call MAPL_AddInternalSpec ( gc,                       &
             SHORT_NAME = 'OA4',                              &
             LONG_NAME  = 'orographic assymetry',             &
             UNGRIDDED_DIMS     = (/4/),                      &
             UNITS      = '1',                                &
             DIMS       = MAPL_DimsHorzOnly,                  &
             VLOCATION  = MAPL_VLocationNone,                 &
             RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddInternalSpec ( gc,                       &
             SHORT_NAME = 'CLX4',                             &
             LONG_NAME  = 'fractional area',                  &
             UNGRIDDED_DIMS     = (/4/),                      &
             UNITS      = '1',                                &
             DIMS       = MAPL_DimsHorzOnly,                  &
             VLOCATION  = MAPL_VLocationNone,                 &
             RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddInternalSpec ( gc,                       &
             SHORT_NAME = 'THETA',                            &
             LONG_NAME  = 'angle of mnt with east (x)',       &
             UNITS      = '1',                                &
             DIMS       = MAPL_DimsHorzOnly,                  &
             VLOCATION  = MAPL_VLocationNone,                 &
             RC=STATUS  )
        VERIFY_(STATUS)
        
        call MAPL_AddInternalSpec ( gc,                       &
             SHORT_NAME = 'SIGMA',                            &
             LONG_NAME  = 'orographic slope',                 &
             UNITS      = '1',                                &
             DIMS       = MAPL_DimsHorzOnly,                  &
             VLOCATION  = MAPL_VLocationNone,                 &
             RC=STATUS  )
        VERIFY_(STATUS)
        
        call MAPL_AddInternalSpec ( gc,                       &
             SHORT_NAME = 'GAMMA',                            &
             LONG_NAME  = 'orographic anisotropy',            &
             UNITS      = '1',                                &
             DIMS       = MAPL_DimsHorzOnly,                  &
             VLOCATION  = MAPL_VLocationNone,                 &
             RC=STATUS  )
        VERIFY_(STATUS)
        
        call MAPL_AddInternalSpec ( gc,                       &
             SHORT_NAME = 'ELVMAX',                           &
             LONG_NAME  = 'orographic maximum',               &
             UNITS      = '1',                                &
             DIMS       = MAPL_DimsHorzOnly,                  &
             VLOCATION  = MAPL_VLocationNone,                 &
             RC=STATUS  )
        VERIFY_(STATUS)

     end if ! NCEP
!_end_of_if_ncep
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

  character(len=ESMF_MAXSTR) :: GRIDNAME
  character(len=4)           :: imchar
  character(len=2)           :: dateline
  integer                    :: imsize,nn

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

! Get grid name to determine IMSIZE
    call MAPL_GetResource(MAPL,GRIDNAME,'AGCM_GRIDNAME:', RC=STATUS)
    VERIFY_(STATUS)
    GRIDNAME =  AdjustL(GRIDNAME)
    nn = len_trim(GRIDNAME)
    dateline = GRIDNAME(nn-1:nn)
    imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
    read(imchar,*) imsize
    if(dateline.eq.'CF') imsize = imsize*4

! Gravity wave drag
! -----------------

    if (USE_NCEP_GWD) then
! NCEP  FV3 values
! C768:  cdmbgwd="3.500,0.25"
! C384:  cdmbgwd="1.000,1.20"
! C192:  cdmbgwd="0.200,2.50"
! C96 :  cdmbgwd="0.125,3.00"
! C48 :  cdmbgwd="0.062,3.50"
    if ( imsize.lt.270 ) then
      call MAPL_GetResource( MAPL, CDMBGWD1, Label="CDMBGWD1:", default=0.062, RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetResource( MAPL, CDMBGWD2, Label="CDMBGWD2:", default=3.500, RC=STATUS)
      VERIFY_(STATUS)
    endif
    if ( imsize.ge.270 .and. imsize.lt.540 ) then
      call MAPL_GetResource( MAPL, CDMBGWD1, Label="CDMBGWD1:", default=0.125, RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetResource( MAPL, CDMBGWD2, Label="CDMBGWD2:", default=3.000, RC=STATUS)
      VERIFY_(STATUS)
    endif
    if ( imsize.ge.540 .and. imsize.lt.1080 ) then
      call MAPL_GetResource( MAPL, CDMBGWD1, Label="CDMBGWD1:", default=0.200, RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetResource( MAPL, CDMBGWD2, Label="CDMBGWD2:", default=2.500, RC=STATUS)
      VERIFY_(STATUS)
    endif
    if ( imsize.ge.1080 .and. imsize.lt.3240 ) then
      call MAPL_GetResource( MAPL, CDMBGWD1, Label="CDMBGWD1:", default=1.000, RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetResource( MAPL, CDMBGWD2, Label="CDMBGWD2:", default=1.200, RC=STATUS)
      VERIFY_(STATUS)
    endif
    if ( imsize.ge.3240 ) then
      call MAPL_GetResource( MAPL, CDMBGWD1, Label="CDMBGWD1:", default=3.500, RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetResource( MAPL, CDMBGWD2, Label="CDMBGWD2:", default=0.250, RC=STATUS)
      VERIFY_(STATUS)
    endif
    endif

    call MAPL_GetResource( MAPL, effgworo, Label="EFFGWORO:", default=0.250, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, effgwbkg, Label="EFFGWBKG:", default=0.125, RC=STATUS)
    VERIFY_(STATUS)

    if( LM .le. 72 ) then
        call MAPL_GetResource( MAPL, pgwv,        Label="PGWV:",        default=4,    RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetResource( MAPL, bgstressmax, Label="BGSTRESSMAX:", default=0.9,  RC=STATUS)
        VERIFY_(STATUS)
     else
        call MAPL_GetResource( MAPL, pgwv,        Label="PGWV:",        default=8,    RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetResource( MAPL, bgstressmax, Label="BGSTRESSMAX:", default=2.25, RC=STATUS)
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

#ifdef _CUDA
      type(dim3) :: Grid, Block
      integer :: blocksize
#endif

! NCEP gwd related vars
      real, pointer :: TRATE(:,:,:)=>NULL()
      real          :: CDMBGWD(2)
      logical       :: LPRNT
      logical, allocatable :: KCNV(:,:)
      integer       :: IMX
      integer       :: IPR, ME, LAT, KDT
      integer       :: NMTVR
      integer       :: I,IRUN
      integer       :: IX, IY
      real          :: FV, FHOUR
      real          :: A, QM
      real, allocatable :: PK(:,:,:)
      integer, allocatable :: KPBL(:,:)
      integer, allocatable :: KBOT(:,:)
      integer, allocatable :: KTOP(:,:)
      real, allocatable :: QMAX(:,:)
      real, pointer     :: fPBL(:,:) => NULL()
      real, pointer     :: CLDF(:,:) => NULL()
      real, pointer     :: HPRIME(:,:) => NULL()
      real, pointer     :: OC(:,:) => NULL()
      real, pointer     :: SIGMA(:,:) => NULL()
      real, pointer     :: GAMMA(:,:) => NULL()
      real, pointer     :: THETA(:,:) => NULL()
      real, pointer     :: DLENGTH(:,:) => NULL()
      real, pointer     :: ELVMAX(:,:) => NULL()
      real, pointer     :: OA4(:,:,:) => NULL()
      real, pointer     :: CLX4(:,:,:) => NULL()
      type (ESMF_State) :: INTERNAL
      type (ESMF_Grid)  :: esmfgrid
      integer           :: COUNTS(3)

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

#ifdef _CUDA

      _ASSERT(  LM <= GPU_MAXLEVS,'needs informative message') ! If these are tripped GNUmakefile
      _ASSERT(PGWV <= MAXPGWV,'needs informative message')     ! must be modified.

      call MAPL_GetResource(MAPL,BLOCKSIZE,'BLOCKSIZE:',DEFAULT=128,RC=STATUS)
      VERIFY_(STATUS)

      Block = dim3(blocksize,1,1)
      Grid  = dim3(ceiling(real(IM*JM)/real(blocksize)),1,1)

      call MAPL_TimerOn(MAPL,"-DRIVER_ALLOC",RC=STATUS)
      VERIFY_(STATUS)

      ! ----------------------
      ! Allocate device arrays
      ! ----------------------
  
      ! Working Arrays
      ! --------------
  
      ALLOCATE(ZM_DEV(IM*JM,LM),STAT=STATUS)       ! height above surface at layers
      ALLOCATE(LNPINT_DEV(IM*JM,LM+1),STAT=STATUS) ! log(pint)
      ALLOCATE(PMLN_DEV(IM*JM,LM),STAT=STATUS)     ! log midpoint pressures
      ALLOCATE(PMID_DEV(IM*JM,LM),STAT=STATUS)     ! pressure at the layers
      ALLOCATE(RPDEL_DEV(IM*JM,LM),STAT=STATUS)    ! 1.0 / pdel
  
      ! Inputs - PREGEO
      ! ---------------
  
      ALLOCATE(PINT_DEV(IM*JM,LM+1),STAT=STATUS)   ! pressure at the layer edges
      ALLOCATE(RLAT_DEV(IM*JM),STAT=STATUS)        ! latitude in radian
  
      ! Outputs - PREGEO
      ! ----------------
  
      ALLOCATE(PDEL_DEV(IM*JM,LM),STAT=STATUS)     ! pressure thickness at the layers
  
      ! Inputs - GEOPOT
      ! ---------------
  
      ALLOCATE(T_DEV(IM*JM,LM),STAT=STATUS)        ! temperature at layers
      ALLOCATE(Q_DEV(IM*JM,LM),STAT=STATUS)        ! specific humidity
  
      ! Outputs - GEOPOT
      ! ----------------
  
      ALLOCATE(ZI_DEV(IM*JM,LM+1),STAT=STATUS)     ! Height above surface at interfaces
  
      ! Inputs - INTR
      ! -------------
  
      ALLOCATE(U_DEV(IM*JM,LM),STAT=STATUS)        ! zonal wind at layers
      ALLOCATE(V_DEV(IM*JM,LM),STAT=STATUS)        ! meridional wind at layers
      ALLOCATE(SGH_DEV(IM*JM),STAT=STATUS)         ! standard deviation of orography
      ALLOCATE(PREF_DEV(LM+1),STAT=STATUS)         ! reference pressure at the layeredges
  
      ! Outputs - INTR
      ! --------------
    
      ALLOCATE(DUDT_GWD_DEV(IM*JM,LM),STAT=STATUS) ! zonal wind tendency at layer 
      ALLOCATE(DVDT_GWD_DEV(IM*JM,LM),STAT=STATUS) ! meridional wind tendency at layer 
      ALLOCATE(DTDT_GWD_DEV(IM*JM,LM),STAT=STATUS) ! temperature tendency at layer
      ALLOCATE(DUDT_ORG_DEV(IM*JM,LM),STAT=STATUS) ! zonal wind tendency at layer due to orography GWD
      ALLOCATE(DVDT_ORG_DEV(IM*JM,LM),STAT=STATUS) ! meridional wind tendency at layer  due to orography GWD
      ALLOCATE(DTDT_ORG_DEV(IM*JM,LM),STAT=STATUS) ! temperature tendency at layer  due to orography GWD
      ALLOCATE(TAUGWDX_DEV(IM*JM),STAT=STATUS)     ! zonal      gravity wave surface    stress
      ALLOCATE(TAUGWDY_DEV(IM*JM),STAT=STATUS)     ! meridional gravity wave surface    stress
      ALLOCATE(TAUOX_DEV(IM*JM,LM+1),STAT=STATUS)  ! zonal      orographic gravity wave stress
      ALLOCATE(TAUOY_DEV(IM*JM,LM+1),STAT=STATUS)  ! meridional orographic gravity wave stress
      ALLOCATE(FEO_DEV(IM*JM,LM+1),STAT=STATUS)    ! energy flux of orographic gravity waves
      ALLOCATE(FEPO_DEV(IM*JM,LM+1),STAT=STATUS)   ! pseudoenergy flux of orographic gravity waves
      ALLOCATE(TAUBKGX_DEV(IM*JM),STAT=STATUS)     ! zonal      gravity wave background stress
      ALLOCATE(TAUBKGY_DEV(IM*JM),STAT=STATUS)     ! meridional gravity wave background stress
      ALLOCATE(TAUBX_DEV(IM*JM,LM+1),STAT=STATUS)  ! zonal      background gravity wave stress
      ALLOCATE(TAUBY_DEV(IM*JM,LM+1),STAT=STATUS)  ! meridional background gravity wave stress
      ALLOCATE(FEB_DEV(IM*JM,LM+1),STAT=STATUS)    ! energy flux of background gravity waves
      ALLOCATE(FEPB_DEV(IM*JM,LM+1),STAT=STATUS)   ! pseudoenergy flux of background gravity waves
      ALLOCATE(UTBSRC_DEV(IM*JM,LM),STAT=STATUS)   ! dU/dt below background launch level
      ALLOCATE(VTBSRC_DEV(IM*JM,LM),STAT=STATUS)   ! dV/dt below background launch level
      ALLOCATE(TTBSRC_DEV(IM*JM,LM),STAT=STATUS)   ! dT/dt below background launch level
  
      ! Outputs - POSTINTR
      ! ------------------
  
      ALLOCATE(DUDT_TOT_DEV(IM*JM,LM),STAT=STATUS) ! Tendency of eastward wind due to GWD
      ALLOCATE(DVDT_TOT_DEV(IM*JM,LM),STAT=STATUS) ! Tendency of northward wind due to GWD
      ALLOCATE(DTDT_TOT_DEV(IM*JM,LM),STAT=STATUS) ! Tendency of air temperature due to GWD
      ALLOCATE(DUDT_RAH_DEV(IM*JM,LM),STAT=STATUS) ! Tendency of eastward wind due to Rayleigh friction
      ALLOCATE(DVDT_RAH_DEV(IM*JM,LM),STAT=STATUS) ! Tendency of northward wind due to Rayleigh friction
      ALLOCATE(DTDT_RAH_DEV(IM*JM,LM),STAT=STATUS) ! Tendency of air temperature due to Rayleigh friction
      ALLOCATE(PEGWD_DEV(IM*JM),STAT=STATUS)       ! Potential energy tendency across GWD
      ALLOCATE(PEORO_DEV(IM*JM),STAT=STATUS)       ! Potential energy tendency due to orographic gravity
      ALLOCATE(PERAY_DEV(IM*JM),STAT=STATUS)       ! Potential energy tendency due to Rayleigh friction
      ALLOCATE(PEBKG_DEV(IM*JM),STAT=STATUS)       ! Potential energy tendency due to gw background
      ALLOCATE(KEGWD_DEV(IM*JM),STAT=STATUS)       ! Kinetic energy tendency across GWD
      ALLOCATE(KEORO_DEV(IM*JM),STAT=STATUS)       ! Kinetic energy tendency due to orographic gravity
      ALLOCATE(KERAY_DEV(IM*JM),STAT=STATUS)       ! Kinetic energy tendency due to Rayleigh friction
      ALLOCATE(KEBKG_DEV(IM*JM),STAT=STATUS)       ! Kinetic energy tendency due to gw background
      ALLOCATE(KERES_DEV(IM*JM),STAT=STATUS)       ! Kinetic energy residual for total energy conservation
      ALLOCATE(BKGERR_DEV(IM*JM),STAT=STATUS)      ! Kinetic energy residual for BKG energy conservation

      call MAPL_TimerOff(MAPL,"-DRIVER_ALLOC",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"-DRIVER_DATA",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"--DRIVER_DATA_DEVICE",RC=STATUS)
      VERIFY_(STATUS)

      ! ---------------------
      ! Copy inputs to device
      ! ---------------------
  
      ! Inputs
      ! ------
  
      STATUS = cudaMemcpy(PINT_DEV,PLE,IM*JM*(LM+1))
      STATUS = cudaMemcpy(RLAT_DEV,LATS,IM*JM)
      STATUS = cudaMemcpy(T_DEV,T,IM*JM*LM)
      STATUS = cudaMemcpy(Q_DEV,Q,IM*JM*LM)
      STATUS = cudaMemcpy(U_DEV,U,IM*JM*LM)
      STATUS = cudaMemcpy(V_DEV,V,IM*JM*LM)
      STATUS = cudaMemcpy(SGH_DEV,SGH,IM*JM)
      STATUS = cudaMemcpy(PREF_DEV,PREF,LM+1)

      call MAPL_TimerOff(MAPL,"--DRIVER_DATA_DEVICE",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOff(MAPL,"-DRIVER_DATA",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"-DRIVER_RUN",RC=STATUS)
      VERIFY_(STATUS)

      CALL PREGEO<<<Grid,Block>>>(IM*JM, LM, &
            PINT_DEV, RLAT_DEV, PMID_DEV, PDEL_DEV, RPDEL_DEV, LNPINT_DEV, PMLN_DEV)

      STATUS = cudaGetLastError()
      if (STATUS /= 0) then 
         write (*,*) "Error code from PREGEO kernel call: ", STATUS
         write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
         _ASSERT(.FALSE.,'needs informative message')
      end if

      ! Compute ZM
      !-------------

      call GEOPOTENTIAL<<<Grid,Block>>>(IM*JM,  LM,                        &
            LNPINT_DEV, PMLN_DEV, PINT_DEV, PMID_DEV, PDEL_DEV, RPDEL_DEV, &
            T_DEV,      Q_DEV,    ZI_DEV,   ZM_DEV                         )

      STATUS = cudaGetLastError()
      if (STATUS /= 0) then 
         write (*,*) "Error code from GEOPOTENTIAL kernel call: ", STATUS
         write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
         _ASSERT(.FALSE.,'needs informative message')
      end if

      call MAPL_TimerOn(MAPL,"-INTR",RC=STATUS)
      VERIFY_(STATUS)
  
      !Do gravity wave drag calculations on a list of soundings
      !---------------------------------------------------------
  
      call GW_INTR<<<Grid,Block>>>(IM*JM, LM, DT, PGWV, BGSTRESSMAX, effgworo, effgwbkg)
  
      STATUS = cudaGetLastError()
      if (STATUS /= 0) then 
         write (*,*) "Error code from GW_INTR kernel call: ", STATUS
         write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
         _ASSERT(.FALSE.,'needs informative message')
      end if
  
      call MAPL_TimerOff(MAPL,"-INTR",RC=STATUS)
      VERIFY_(STATUS)

      CALL POSTINTR<<<Grid,Block>>>(IM*JM, LM, DT, H0, HH, Z1, TAU1, &
            ! Inputs
            PREF_DEV, &
            PDEL_DEV, &
            U_DEV, &
            V_DEV, &
            DUDT_GWD_DEV, &
            DVDT_GWD_DEV, &
            DTDT_GWD_DEV, &
            DUDT_ORG_DEV, &
            DVDT_ORG_DEV, &
            DTDT_ORG_DEV, &

            ! Outputs
            DUDT_TOT_DEV, &
            DVDT_TOT_DEV, &
            DTDT_TOT_DEV, &
            DUDT_RAH_DEV, &
            DVDT_RAH_DEV, &
            DTDT_RAH_DEV, &
            PEGWD_DEV,    &
            PEORO_DEV,    &
            PERAY_DEV,    &
            PEBKG_DEV,    &
            KEGWD_DEV,    &
            KEORO_DEV,    &
            KERAY_DEV,    &
            KEBKG_DEV,    &
            KERES_DEV,    &
            BKGERR_DEV    )

      STATUS = cudaGetLastError()
      if (STATUS /= 0) then 
         write (*,*) "Error code from POSTINTR kernel call: ", STATUS
         write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
         _ASSERT(.FALSE.,'needs informative message')
      end if
  
      call MAPL_TimerOff(MAPL,"-DRIVER_RUN",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"-DRIVER_DATA",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"--DRIVER_DATA_DEVICE",RC=STATUS)
      VERIFY_(STATUS)

      ! ------------------------
      ! Copy outputs from device
      ! ------------------------
  
      ! Outputs - PREGEO
      ! ----------------

      ! MAT PDEL might be needed below to weight DTDT
      if(associated(DTDT)) STATUS = cudaMemcpy(PDEL,PDEL_DEV,IM*JM*LM)

      ! Outputs - GEOPOT
      ! ----------------

      !GPU This array is not used anywhere past this point. Uncomment
      !GPU to copy it to the host if needed in future.
      !GPU  STATUS = cudaMemcpy(ZI,ZI_DEV,IM*JM*(LM+1)) !MAT Not used after this

      ! Outputs - INTR
      ! --------------
    
      STATUS = cudaMemcpy(DUDT_GWD,DUDT_GWD_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(DVDT_GWD,DVDT_GWD_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(DTDT_GWD,DTDT_GWD_DEV,IM*JM*LM)

      STATUS = cudaMemcpy(DUDT_ORG,DUDT_ORG_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(DVDT_ORG,DVDT_ORG_DEV,IM*JM*LM)
      STATUS = cudaMemcpy(DTDT_ORG,DTDT_ORG_DEV,IM*JM*LM)

      STATUS = cudaMemcpy(TAUXO_TMP,TAUGWDX_DEV,IM*JM)
      STATUS = cudaMemcpy(TAUYO_TMP,TAUGWDY_DEV,IM*JM)
      STATUS = cudaMemcpy(TAUXB_TMP,TAUBKGX_DEV,IM*JM)
      STATUS = cudaMemcpy(TAUYB_TMP,TAUBKGY_DEV,IM*JM)

      !GPU These arrays are not used anywhere past this point. Uncomment
      !GPU to copy them to the host if needed in future.
      !GPU  STATUS = cudaMemcpy(TAUXO_3D,TAUOX_DEV,IM*JM*(LM+1))
      !GPU  STATUS = cudaMemcpy(TAUYO_3D,TAUOY_DEV,IM*JM*(LM+1))
      !GPU  STATUS = cudaMemcpy(FEO_3D,FEO_DEV,IM*JM*(LM+1))
      !GPU  STATUS = cudaMemcpy(FEPO_3D,FEPO_DEV,IM*JM*(LM+1))
      !GPU  STATUS = cudaMemcpy(TAUXB_3D,TAUBX_DEV,IM*JM*(LM+1))
      !GPU  STATUS = cudaMemcpy(TAUYB_3D,TAUBY_DEV,IM*JM*(LM+1))
      !GPU  STATUS = cudaMemcpy(FEB_3D,FEB_DEV,IM*JM*(LM+1))
      !GPU  STATUS = cudaMemcpy(FEPB_3D,FEPB_DEV,IM*JM*(LM+1))
      !GPU  STATUS = cudaMemcpy(DUBKGSRC,UTBSRC_DEV,IM*JM*LM)
      !GPU  STATUS = cudaMemcpy(DVBKGSRC,VTBSRC_DEV,IM*JM*LM)
      !GPU  STATUS = cudaMemcpy(DTBKGSRC,TTBSRC_DEV,IM*JM*LM)

      ! Outputs - POSTINTR
      ! ------------------
    
      if(associated(DUDT    )) STATUS = cudaMemcpy(DUDT,    DUDT_TOT_DEV,IM*JM*LM)
      if(associated(DVDT    )) STATUS = cudaMemcpy(DVDT,    DVDT_TOT_DEV,IM*JM*LM)
      if(associated(DTDT    )) STATUS = cudaMemcpy(DTDT,    DTDT_TOT_DEV,IM*JM*LM)

      if(associated(DUDT_RAY)) STATUS = cudaMemcpy(DUDT_RAY,DUDT_RAH_DEV,IM*JM*LM)
      if(associated(DVDT_RAY)) STATUS = cudaMemcpy(DVDT_RAY,DVDT_RAH_DEV,IM*JM*LM)
      if(associated(DTDT_RAY)) STATUS = cudaMemcpy(DTDT_RAY,DTDT_RAH_DEV,IM*JM*LM)

      ! KE dIagnostics
      !----------------

      if(associated(PEGWD )) STATUS = cudaMemcpy(PEGWD, PEGWD_DEV,IM*JM)
      if(associated(PEORO )) STATUS = cudaMemcpy(PEORO, PEORO_DEV,IM*JM)
      if(associated(PERAY )) STATUS = cudaMemcpy(PERAY, PERAY_DEV,IM*JM)
      if(associated(PEBKG )) STATUS = cudaMemcpy(PEBKG, PEBKG_DEV,IM*JM)
      if(associated(KEGWD )) STATUS = cudaMemcpy(KEGWD, KEGWD_DEV,IM*JM)
      if(associated(KEORO )) STATUS = cudaMemcpy(KEORO, KEORO_DEV,IM*JM)
      if(associated(KERAY )) STATUS = cudaMemcpy(KERAY, KERAY_DEV,IM*JM)
      if(associated(KEBKG )) STATUS = cudaMemcpy(KEBKG, KEBKG_DEV,IM*JM)
      if(associated(KERES )) STATUS = cudaMemcpy(KERES, KERES_DEV,IM*JM)
      if(associated(BKGERR)) STATUS = cudaMemcpy(BKGERR,BKGERR_DEV,IM*JM)
  
      call MAPL_TimerOff(MAPL,"--DRIVER_DATA_DEVICE",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOff(MAPL,"-DRIVER_DATA",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"-DRIVER_DEALLOC",RC=STATUS)
      VERIFY_(STATUS)

      ! ------------------------
      ! Deallocate device arrays
      ! ------------------------
  
      ! Working Arrays
      ! --------------
  
      DEALLOCATE(ZM_DEV)
      DEALLOCATE(LNPINT_DEV)
      DEALLOCATE(PMLN_DEV)
      DEALLOCATE(PMID_DEV)
      DEALLOCATE(RPDEL_DEV)
  
      ! Inputs - PREGEO
      ! ---------------
  
      DEALLOCATE(PINT_DEV)
      DEALLOCATE(RLAT_DEV)
  
      ! Outputs - PREGEO
      ! ----------------
  
      DEALLOCATE(PDEL_DEV)
  
      ! Inputs - GEOPOT
      ! ---------------
  
      DEALLOCATE(T_DEV)
      DEALLOCATE(Q_DEV)
  
      ! Outputs - GEOPOT
      ! ----------------
  
      DEALLOCATE(ZI_DEV)
  
      ! Inputs - INTR
      ! -------------
  
      DEALLOCATE(U_DEV)
      DEALLOCATE(V_DEV)
      DEALLOCATE(SGH_DEV)
      DEALLOCATE(PREF_DEV)
  
      ! Outputs - INTR
      ! --------------
    
      DEALLOCATE(DUDT_GWD_DEV)
      DEALLOCATE(DVDT_GWD_DEV)
      DEALLOCATE(DTDT_GWD_DEV)
      DEALLOCATE(DUDT_ORG_DEV)
      DEALLOCATE(DVDT_ORG_DEV)
      DEALLOCATE(DTDT_ORG_DEV)
      DEALLOCATE(TAUGWDX_DEV)
      DEALLOCATE(TAUGWDY_DEV)
      DEALLOCATE(TAUOX_DEV)
      DEALLOCATE(TAUOY_DEV)
      DEALLOCATE(FEO_DEV)
      DEALLOCATE(FEPO_DEV)
      DEALLOCATE(TAUBKGX_DEV)
      DEALLOCATE(TAUBKGY_DEV)
      DEALLOCATE(TAUBX_DEV)
      DEALLOCATE(TAUBY_DEV)
      DEALLOCATE(FEB_DEV)
      DEALLOCATE(FEPB_DEV)
      DEALLOCATE(UTBSRC_DEV)
      DEALLOCATE(VTBSRC_DEV)
      DEALLOCATE(TTBSRC_DEV)
  
      ! Outputs - POSTINTR
      ! ------------------
  
      DEALLOCATE(DUDT_TOT_DEV)
      DEALLOCATE(DVDT_TOT_DEV)
      DEALLOCATE(DTDT_TOT_DEV)
      DEALLOCATE(DUDT_RAH_DEV)
      DEALLOCATE(DVDT_RAH_DEV)
      DEALLOCATE(DTDT_RAH_DEV)
      DEALLOCATE(PEGWD_DEV)
      DEALLOCATE(PEORO_DEV)
      DEALLOCATE(PERAY_DEV)
      DEALLOCATE(PEBKG_DEV)
      DEALLOCATE(KEGWD_DEV)
      DEALLOCATE(KEORO_DEV)
      DEALLOCATE(KERAY_DEV)
      DEALLOCATE(KEBKG_DEV)
      DEALLOCATE(KERES_DEV)
      DEALLOCATE(BKGERR_DEV)

      call MAPL_TimerOff(MAPL,"-DRIVER_DEALLOC",RC=STATUS)
      VERIFY_(STATUS)

#else

      CALL PREGEO(IM*JM,   LM,   &
                    PLE, LATS,   PMID,  PDEL, RPDEL,     PILN,     PMLN)

! Compute ZM
!-------------

      call GEOPOTENTIAL( IM*JM, LM,                  &
           PILN,  PMLN,  PLE,   PMID, PDEL, RPDEL,   &
           T,     Q,     ZI,    ZM                   )

! Do gravity wave drag calculations on a list of soundings
!---------------------------------------------------------

    call MAPL_TimerOn(MAPL,"-INTR")
    if (.not. USE_NCEP_GWD) then
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
    else
       ! get pointers from INTERNAL:HPRIME,OC,OA4,CLX4,THETA,SIGMA,GAMMA,ELVMAX
       call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer( INTERNAL, HPRIME, 'STDDEV', RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetPointer( INTERNAL, OC, 'CONVEXITY', RC=STATUS )
       VERIFY_(STATUS)
       ! next 2 have ungridded dim = 4
       call MAPL_GetPointer( INTERNAL, OA4, 'OA4', RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetPointer( INTERNAL, CLX4, 'CLX4', RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetPointer( INTERNAL, THETA, 'THETA', RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetPointer( INTERNAL, SIGMA, 'SIGMA', RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetPointer( INTERNAL, GAMMA, 'GAMMA', RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetPointer( INTERNAL, ELVMAX, 'ELVMAX', RC=STATUS )
       VERIFY_(STATUS)

       KDT = 1 !??? current number of time step (not used)
       IRUN = IM*JM
       IX = IRUN
       IY = IRUN
       DVDT_ORG = 0.0
       DUDT_ORG = 0.0
       DTDT_ORG = 0.0
       ! Make sure we have: U, V, T, Q, KPBL,
       call MAPL_GetPointer( IMPORT, FPBL, 'KPBL', RC=STATUS )
       VERIFY_(STATUS)
       allocate(KPBL(IM,JM),stat=status)
       VERIFY_(STATUS)
       KPBL = nint(fPBL)
!ALT: we need some protection      
       KPBL = MIN(KPBL, LM-1) ! when LM=132, we might need to cap it to LM-3

       ! PLE, PL,DELP
       allocate(PK(IM,JM,0:LM),stat=status)
       VERIFY_(STATUS)

       PK = ( PLE/MAPL_P00 )**MAPL_KAPPA ! Question should P00 be there

       IMX = imsize     ! number of longitudes (along the equator = 4xCnumber)
       nmtvr = 14       ! number of topographic variables such as variances
       CDMBGWD = [CDMBGWD1,CDMBGWD2]  ! multiplication factors for cdmb and gwd
       ! me, lprint, ipr are currently not being used
       me = 0 ! we could get it from VM
       lprnt = .false.
       ipr = 0

       call GWDPS(IRUN,IX,IY,LM,DVDT_ORG,DUDT_ORG,DTDT_ORG,&
            U,V,T,Q,KPBL,              &
            PLE,PDEL,PMID,PK,&
            ZI,ZM,DT,KDT,         &
            HPRIME,OC,OA4,CLX4,THETA,SIGMA,GAMMA,ELVMAX,       &
            TAUXO_TMP,TAUYO_TMP, &
            MAPL_GRAV, MAPL_CP, MAPL_RDRY, MAPL_RVAP, IMX,  &
            nmtvr, cdmbgwd, me, lprnt, ipr)

       ! Get IMPORT DTDT_moist
!ALT: BE careful. The NCEP routines expect 2d packed into 1d 
       ! units
       ! min/max search potentially could fail
       ! scaling coefficient values???
       ! ============================
       ! kind_phys!!! check with Bill
       ! for now 4, otherwise this will create interface problems
       ! ============================
       
       LAT = 1 ! latitude index for debugging prints, not used
       FHOUR = 0.0 ! forecast hour, not used
       FV = MAPL_RDRY/MAPL_RVAP - 1.0     ! con_fvirt = con_rv/con_rd-1

       ! GWDC needs a grid spacing argument, Bill suggesed we use DXC
       call MAPL_GetPointer( IMPORT, DLENGTH, 'DXC', RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_GetPointer( IMPORT, CLDF, 'CNV_FRC', RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetPointer( IMPORT, TRATE, 'DTDT_moist', RC=STATUS )
       VERIFY_(STATUS)
       ! for every i,j search loop over levels to find QMAX, KBOT and KTOP
       ! lat is not used

       allocate(QMAX(IM,JM), KBOT(IM,JM), KTOP(IM,JM), KCNV(IM,JM), stat=status)
       VERIFY_(STATUS)
       KCNV = .false.
       
       DO J=1,JM
          DO I=1,IM
             QM = MAXVAL(TRATE(I,J,:))
             K = MAXLOC(TRATE(I,J,:),DIM=1)
             QMAX(I,J) = QM
             A = QM
             DO L = K+1, LM
                IF(A > TRATE(I,J,L)) then
                   A = TRATE(I,J,L)
                ELSE
                   EXIT
                END IF
             END DO
             KTOP(I,J) = L-1

             A = QM
             DO L = K-1, 1, -1
                IF(A < TRATE(I,J,L)) then
                   A = TRATE(I,J,L)
                ELSE
                   EXIT
                END IF
             END DO
             KBOT(I,J) = L+1

             IF(CLDF(I,J) > 0.15) THEN
                KCNV(I,J) = .TRUE.
             END IF
          END DO
       END DO
       call gwdc(IRUN,IX,IY,LM,LAT,U,V,T,Q,DT, &
            PMID,PLE,PDEL,QMAX,KTOP,KBOT,KCNV,CLDF, &
            MAPL_GRAV,MAPL_CP,MAPL_RDRY,FV,MAPL_PI,&
            DLENGTH,LPRNT,IPR,FHOUR, &
            DUDT_TOT,DVDT_TOT,TAUXB_TMP,TAUYB_TMP)

       ! Adjust to prepare vars for POSTINTR
       DTDT_GWD = DTDT_ORG 
       DUDT_GWD = DUDT_ORG + DUDT_TOT  
       DVDT_GWD = DVDT_ORG + DVDT_TOT  
! reset total; they were used only as temporary variables
       DUDT_TOT = 0.0
       DVDT_TOT = 0.0

       deallocate(KCNV, KTOP, KBOT, QMAX)
       deallocate(PK)
       deallocate(KPBL)


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

#endif

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


#ifdef _CUDA
  attributes(global) &
#endif
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
#ifdef _CUDA
    integer, value :: pcols                ! Number of longitudes
    integer, value :: pver                 ! Number of vertical layers
#else
    integer, intent(in) :: pcols                ! Number of longitudes
    integer, intent(in) :: pver                 ! Number of vertical layers
#endif

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

#ifdef _CUDA
    i = (blockidx%x - 1) * blockdim%x + threadidx%x

    I_LOOP: if ( i <= pcols ) then
#else
    I_LOOP: do i = 1, pcols
#endif
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
#ifndef _CUDA
    end do I_LOOP
#else
    end if I_LOOP
#endif 

    return
  end subroutine geopotential

!----------------------------------------------------------------------- 

#ifdef _CUDA
  attributes(global) &
#endif
  subroutine pregeo(pcols,pver,&
    ple,lats,pmid,pdel,rpdel,piln,pmln)

    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!

#ifdef _CUDA
    integer, value :: pcols         ! Number of longitudes
    integer, value :: pver          ! Number of vertical layers
#else
    integer, intent(in) :: pcols    ! Number of longitudes
    integer, intent(in) :: pver     ! Number of vertical layers
#endif

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

#ifdef _CUDA
    I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

    I_LOOP: IF ( I <= PCOLS ) THEN
#else
    I_LOOP: DO I = 1, PCOLS
#endif
       DO K = 1, PVER
           PMID(I,K) = 0.5*(  PLE(I,K  ) + PLE(I,K+1) )
           PDEL(I,K) =        PLE(I,K+1) - PLE(I,K  )
          RPDEL(I,K) = 1.0 / PDEL(I,K)
           PILN(I,K) = log(   PLE(I,K) )
           PMLN(I,K) = log(  PMID(I,K) ) !
       END DO
       PILN(I,PVER+1)  = log( PLE(I,PVER+1)  )
#ifndef _CUDA
    END DO I_LOOP
#else
    END IF I_LOOP
#endif 

  end subroutine pregeo

#ifdef _CUDA
  attributes(global) &
#endif
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

#ifdef _CUDA
    integer, value :: PCOLS ! Number of longitudes
    integer, value :: PVER  ! Number of vertical layers
    real,    value :: DT    ! Time step
    real,    value :: H0, HH, Z1, TAU1 ! Rayleigh friction parameters
#else    
    integer, intent(in) :: PCOLS ! Number of longitudes
    integer, intent(in) :: PVER  ! Number of vertical layers
    real,    intent(in) :: DT    ! Time step
    real,    intent(in) :: H0, HH, Z1, TAU1 ! Rayleigh friction parameters
#endif

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

#ifdef _CUDA
    I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

    I_LOOP: IF ( I <= PCOLS ) THEN
#else
    I_LOOP: DO I = 1, PCOLS
#endif
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

#ifndef _CUDA
    END DO I_LOOP
#else
    END IF I_LOOP
#endif 

  end subroutine postintr

end module GEOS_GwdGridCompMod
