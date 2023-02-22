
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

  use gw_rdg, only : gw_rdg_init
  use gw_oro, only : gw_oro_init
  use gw_convect, only : gw_beres_init, BeresSourceDesc
  use gw_common, only: GWBand, gw_common_init, gw_newtonian_set
  use gw_drag_ncar, only: gw_intr_ncar

  use gw_drag, only: gw_intr

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP
! config params
  type :: ThreadWorkspace
     type(GWBand)          :: beres_band
     type(BeresSourceDesc) :: beres_dc_desc, beres_sc_desc
     type(GWBand)          :: oro_band
     type(GWBand)          :: rdg_band
  end type ThreadWorkspace

  type       :: GEOS_GwdGridComp
     real :: GEOS_BGSTRESS
     real :: GEOS_EFFGWBKG
     real :: GEOS_EFFGWORO
     integer :: GEOS_PGWV
     real :: NCAR_EFFGWBKG
     real :: NCAR_EFFGWORO
     integer :: NCAR_NRDG
     real :: Z1
     real :: TAU1
     real :: H0
     real :: HH
     real :: HGT_SURFACE
     real :: effbeljaars, limbeljaars
     real, allocatable :: alpha(:) 
     type(ThreadWorkspace), allocatable :: workspaces(:)
  end type GEOS_GwdGridComp

  type wrap_
     type (GEOS_GwdGridComp), pointer     :: PTR
  end type wrap_

  !logical, save      :: FIRST_RUN = .true.

! Beljaars parameters
   real, parameter ::      &
      dxmin_ss =  3000.0, &        ! minimum grid length for Beljaars
      dxmax_ss = 12000.0           ! maximum grid length for Beljaars
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
    type (MAPL_MetaComp),     pointer   :: MAPL
!=============================================================================
    logical :: use_threads
    type (ESMF_Config)                            :: myCF

    type (wrap_)                                :: wrap
    type (GEOS_GwdGridComp), pointer               :: self
    integer :: num_threads

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC )
    Iam = trim(COMP_NAME) // Iam

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, _STAT)
    wrap%ptr => self

    num_threads = MAPL_get_num_threads()
    allocate(self%workspaces(0:num_threads-1), _STAT)

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_INITIALIZE,  Initialize,  _RC)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,  Run,  _RC)
    
    call MAPL_GetObjectFromGC ( GC, MAPL, _RC )
    
     myCF = ESMF_ConfigCreate (_RC)
     call ESMF_ConfigLoadFile (myCF, 'GWD_GridComp.rc', _RC)
     call ESMF_ConfigGetAttribute (myCF, use_threads, label='use_threads:', default=.FALSE., _RC)
!   set use_threads
    call MAPL%set_use_threads(use_threads)
    call ESMF_ConfigDestroy(myCF, _RC)

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
                                                       _RC  )
     
     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'T',                                         &
        LONG_NAME  = 'air_temperature',                           &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       _RC  )
     
     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'Q',                                         &
        LONG_NAME  = 'specific_humidity',                         &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       _RC  )
     
     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'U',                                         &
        LONG_NAME  = 'eastward_wind',                             &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       _RC  )
     
     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'V',                                         &
        LONG_NAME  = 'northward_wind',                            &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &
        RESTART    = MAPL_RestartSkip,                            &
                                                       _RC  )
     
     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PHIS',                              &
        LONG_NAME          = 'surface geopotential height',       &
        UNITS              = 'm+2 s-2',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART    = MAPL_RestartSkip,                            &
                                                        _RC  )
     
     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'SGH',                                       &
        LONG_NAME  = 'standard_deviation_of_topography',          &
        UNITS      = 'm',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
        RESTART    = MAPL_RestartSkip,                            &
                                                       _RC  )
     
     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'VARFLT',                                    &
        LONG_NAME  = 'variance_of_the_filtered_topography',       &
        UNITS      = 'm+2',                                       &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
        RESTART    = MAPL_RestartSkip,                            &
                                                       _RC  )
     
     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'PREF',                                      &
        LONG_NAME  = 'reference_air_pressure',                    &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsVertOnly,                           &
        VLOCATION  = MAPL_VLocationEdge,                          &
        RESTART    = MAPL_RestartSkip,                            &
                                                       _RC  )
     
     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'AREA',                                      &
        LONG_NAME  = 'grid_box_area',                             &
        UNITS      = 'm^2',                                       &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
        RESTART    = MAPL_RestartSkip,                            &
                                                       _RC  )
    
! from moist
     call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME='DTDT_DC',                               &
         LONG_NAME ='T tendency due to deep convection',     &
         UNITS     ='K s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                      &
         VLOCATION = MAPL_VLocationCenter,              _RC  )
     call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME='DTDT_SC',                               &
         LONG_NAME ='T tendency due to shallow convection',  &
         UNITS     ='K s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                      &
         VLOCATION = MAPL_VLocationCenter,              _RC  )
     call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME = 'DQLDT',                                   &
         LONG_NAME = 'total_liq_water_tendency_due_to_moist',       &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         _RC  )
     call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME= 'DQIDT',                                   &
         LONG_NAME = 'total_ice_water_tendency_due_to_moist',       &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         _RC  )
     
! !EXPORT STATE:

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'PLE',                                       &
        LONG_NAME  = 'air_pressure',                              &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'T',                                         &
        LONG_NAME  = 'air_temperature',                           &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'Q',                                         &
        LONG_NAME  = 'specific_humidity',                         &
        UNITS      = 'kg kg-1',                                   &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'U',                                         &
        LONG_NAME  = 'eastward_wind',                             &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'V',                                         &
        LONG_NAME  = 'northward_wind',                            &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'RDG1_MXDIS',                                &
        LONG_NAME  = 'ridge1_mxdis',                              &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'RDG1_HWDTH',                                &
        LONG_NAME  = 'ridge1_hwdth',                              &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'RDG1_CLNGT',                                &
        LONG_NAME  = 'ridge1_clngt',                              &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'RDG1_ANGLL',                                &
        LONG_NAME  = 'ridge1_angll',                              &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'RDG1_ANIXY',                                &
        LONG_NAME  = 'ridge1_anixy',                              &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'RDG1_GBXAR',                                &
        LONG_NAME  = 'ridge1_gridbox_area',                       &
        UNITS      = 'km^2',                                      &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'SGH',                                       &
        LONG_NAME  = 'standard_deviation_of_topography',          &
        UNITS      = 'm',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'PREF',                                      &
        LONG_NAME  = 'reference_air_pressure',                    &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsVertOnly,                           &
        VLOCATION  = MAPL_VLocationEdge,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DTDT',                                      &
        LONG_NAME  = 'mass_weighted_air_temperature_tendency_due_to_GWD',    &
        UNITS      = 'Pa K s-1',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TTMGW',                                     &
        LONG_NAME  = 'air_temperature_tendency_due_to_GWD',       &
        UNITS      = 'K s-1',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DUDT',                                      &
        LONG_NAME  = 'tendency_of_eastward_wind_due_to_GWD',                 &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DVDT',                                      &
        LONG_NAME  = 'tendency_of_northward_wind_due_to_GWD',                &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DUDT_TFD',                                  &
        LONG_NAME  = 'tendency_of_eastward_wind_due_to_topographic_form_drag',               &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DVDT_TFD',                                  &
        LONG_NAME  = 'tendency_of_northward_wind_due_to_topographic_form_dra',              &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DTDT_ORO',                                  &
        LONG_NAME  = 'air_temperature_tendency_due_to_orographic_GWD', &
        UNITS      = 'K s-1',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DUDT_ORO',                                  &
        LONG_NAME  = 'tendency_of_eastward_wind_due_to_orographic_GWD',               &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DVDT_ORO',                                  &
        LONG_NAME  = 'tendency_of_northward_wind_due_to_orographic_GWD',              &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DTDT_BKG',                                  &
        LONG_NAME  = 'air_temperature_tendency_due_to_background_GWD', &
        UNITS      = 'K s-1',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DUDT_BKG',                                  &
        LONG_NAME  = 'tendency_of_eastward_wind_due_to_background_GWD',               &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DVDT_BKG',                                  &
        LONG_NAME  = 'tendency_of_northward_wind_due_to_background_GWD',              &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DTDT_RAY',                                  &
        LONG_NAME  = 'air_temperature_tendency_due_to_Rayleigh_friction',        &
        UNITS      = 'K s-1',                                  &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DUDT_RAY',                                  &
        LONG_NAME  = 'tendency_of_eastward_wind_due_to_Rayleigh_friction',       &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'DVDT_RAY',                                  &
        LONG_NAME  = 'tendency_of_northward_wind_due_to_Rayleigh_friction',      &
        UNITS      = 'm s-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUGWX',                                    &
        LONG_NAME  = 'surface_eastward_gravity_wave_stress',      &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUGWY',                                    &
        LONG_NAME  = 'surface_northward_gravity_wave_stress',     &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUOROX',                                   &
        LONG_NAME  = 'surface_eastward_orographic_gravity_wave_stress',      &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUOROY',                                   &
        LONG_NAME  = 'surface_northward_orographic_gravity_wave_stress',     &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUBKGX',                                   &
        LONG_NAME  = 'surface_eastward_background_gravity_wave_stress',      &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUBKGY',                                   &
        LONG_NAME  = 'surface_northward_background_gravity_wave_stress',     &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUMSTX',                                   &
        LONG_NAME  = 'surface_eastward_gravity_wave_stress_due_to_Moist_Processes',      &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'TAUMSTY',                                   &
        LONG_NAME  = 'surface_northward_gravity_wave_stress_due_to_Moist_Processes',     &
        UNITS      = 'N m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'CLDSTD',                                    &
        LONG_NAME  = 'gravity_wave_drag_standard_deviation_due_to_clouds',     &
        UNITS      = 'm',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'UBASE',                                     &
        LONG_NAME  = 'eastward_component_of_base_level_wind',     &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'VBASE',                                     &
        LONG_NAME  = 'northward_component_of_base_level_wind',    &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'UBAR',                                      &
        LONG_NAME  = 'eastward_component_of_mean_level_wind',     &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME = 'VBAR',                                      &
        LONG_NAME  = 'northward_component_of_mean_level_wind',    &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               _RC  )
     
    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEGWD',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_across_gwd',         &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        _RC  )
     
    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEORO',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_orographic_gravity_waves',  &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        _RC  )
     
    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEBKG',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_gravity_wave_background',   &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        _RC  )
     
    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PERAY',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_Rayleigh_friction',         &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        _RC  )
     
    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KEGWD',                                                              &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_across_gwd',           &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        _RC  )
     
    call MAPL_AddExportSpec ( gc,                                                                         &
         SHORT_NAME = 'KEORO',                                                                            &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_dissipation_due_to_orographic_gravity_waves', &
         UNITS      = 'W m-2',                                                                            &
         DIMS       = MAPL_DimsHorzOnly,                                                                  &
         VLOCATION  = MAPL_VLocationNone,                                                      _RC  )
     
    call MAPL_AddExportSpec ( gc,                                                                  &
         SHORT_NAME = 'KERAY',                                                                     &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_dissipation_due_to_Rayleigh_friction', &
         UNITS      = 'W m-2',                                                                     &
         DIMS       = MAPL_DimsHorzOnly,                                                           &
         VLOCATION  = MAPL_VLocationNone,                                               _RC  )
     
    call MAPL_AddExportSpec ( gc,                                                                        &
         SHORT_NAME = 'KEBKG',                                                                           &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_dissipation_due_to_gravity_wave_background', &
         UNITS      = 'W m-2',                                                                           &
         DIMS       = MAPL_DimsHorzOnly,                                                                 &
         VLOCATION  = MAPL_VLocationNone,                                                     _RC  )
     
    call MAPL_AddExportSpec ( gc,                                                                        &
         SHORT_NAME = 'KERES',                                                                           &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_residual_for_total_energy_conservation',     &
         UNITS      = 'W m-2',                                                                           &
         DIMS       = MAPL_DimsHorzOnly,                                                                 &
         VLOCATION  = MAPL_VLocationNone,                                                     _RC  )
     
    call MAPL_AddExportSpec ( gc,                                                                        &
         SHORT_NAME = 'BKGERR',                                                                          &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_residual_for_BKG_energy_conservation',       &
         UNITS      = 'W m-2',                                                                           &
         DIMS       = MAPL_DimsHorzOnly,                                                                 &
         VLOCATION  = MAPL_VLocationNone,                                                     _RC  )
     
    call MAPL_AddInternalSpec(GC, &
             SHORT_NAME = 'SGH30', &
             LONG_NAME  = 'standard deviation of 30s elevation from 3km cube', &
             UNITS      = 'm', &
             DIMS       = MAPL_DimsHorzOnly,                    &
             VLOCATION  = MAPL_VLocationNone,              _RC  )
    call MAPL_AddInternalSpec(GC, &
             SHORT_NAME = 'KWVRDG', &
             LONG_NAME  = 'horizonal wwavenumber of mountain ridges', &
             UNITS      = 'km', &
             UNGRIDDED_DIMS     = (/16/),                      &
             DIMS       = MAPL_DimsHorzOnly,                    &
             VLOCATION  = MAPL_VLocationNone,              _RC  )
    call MAPL_AddInternalSpec(GC, &
             SHORT_NAME = 'EFFRDG', &
             LONG_NAME  = 'efficiency of mountain ridge scheme', &
             UNITS      = 'km', &
             UNGRIDDED_DIMS     = (/16/),                      &
             DIMS       = MAPL_DimsHorzOnly,                    &
             VLOCATION  = MAPL_VLocationNone,              _RC  )
    call MAPL_AddInternalSpec(GC, &
             SHORT_NAME = 'GBXAR', &
             LONG_NAME  = 'grid box area', &
             UNITS      = 'NA', &
             DIMS       = MAPL_DimsHorzOnly,                    &
             VLOCATION  = MAPL_VLocationNone,              _RC  )
    call MAPL_AddInternalSpec(GC, &
             SHORT_NAME = 'HWDTH', &
             LONG_NAME  = 'width of mountain ridges', &
             UNITS      = 'km', &
             UNGRIDDED_DIMS     = (/16/),                      &
             DIMS       = MAPL_DimsHorzOnly,                    &
             VLOCATION  = MAPL_VLocationNone,              _RC  )
    call MAPL_AddInternalSpec(GC, &
             SHORT_NAME = 'CLNGT', &
             LONG_NAME  = 'width of mountain ridges', &
             UNITS      = 'km', &
             UNGRIDDED_DIMS     = (/16/),                      &
             DIMS       = MAPL_DimsHorzOnly,                    &
             VLOCATION  = MAPL_VLocationNone,              _RC  )
    call MAPL_AddInternalSpec(GC, &
             SHORT_NAME = 'MXDIS', &
             LONG_NAME  = 'NA', &
             UNITS      = 'NA', &
             UNGRIDDED_DIMS     = (/16/),                      &
             DIMS       = MAPL_DimsHorzOnly,                    &
             VLOCATION  = MAPL_VLocationNone,              _RC  )
    call MAPL_AddInternalSpec(GC, &
             SHORT_NAME = 'ANGLL', &
             LONG_NAME  = 'NA', &
             UNITS      = 'NA', &
             UNGRIDDED_DIMS     = (/16/),                      &
             DIMS       = MAPL_DimsHorzOnly,                    &
             VLOCATION  = MAPL_VLocationNone,              _RC  )
    call MAPL_AddInternalSpec(GC, &
             SHORT_NAME = 'ANIXY', &
             LONG_NAME  = 'NA', &
             UNITS      = 'NA', &
             UNGRIDDED_DIMS     = (/16/),                      &
             DIMS       = MAPL_DimsHorzOnly,                    &
             VLOCATION  = MAPL_VLocationNone,              _RC  )
        
!EOS

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="DRIVER"  ,_RC)
    call MAPL_TimerAdd(GC,    name="-DRIVER_RUN"   ,_RC)
    call MAPL_TimerAdd(GC,    name="-INTR"   ,_RC)
    call MAPL_TimerAdd(GC,    name="-INTR_NCAR"   ,_RC)
    call MAPL_TimerAdd(GC,    name="-INTR_GEOS"   ,_RC)
    call MAPL_TimerAdd(GC,    name="-BELJAARS_TOFD"   ,_RC)
    call MAPL_TimerAdd(GC,    name="-DRIVER_DATA"   ,_RC)
    call MAPL_TimerAdd(GC,    name="--DRIVER_DATA_DEVICE"   ,_RC)
    call MAPL_TimerAdd(GC,    name="--DRIVER_DATA_CONST"   ,_RC)
    call MAPL_TimerAdd(GC,    name="-DRIVER_ALLOC"   ,_RC)
    call MAPL_TimerAdd(GC,    name="-DRIVER_DEALLOC"   ,_RC)
    
!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'GEOS_GwdGridComp', wrap, _RC )

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( gc, _RC)
    
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

    integer                             :: IM, JM
    real, pointer, dimension(:,:)       :: LATS

    character(len=ESMF_MAXSTR) :: GRIDNAME
    character(len=4)           :: imchar
    character(len=2)           :: dateline
    integer                    :: imsize,nn
    integer                    :: LM
    real, pointer, dimension(:)      :: PREF

! NCAR GWD variables

    character(len=ESMF_MAXPATHLEN) :: BERES_FILE_NAME
    character(len=ESMF_MAXSTR)     :: ERRstring

    logical :: NCAR_TAU_TOP_ZERO
    real    :: NCAR_PRNDL
    real    :: NCAR_QBO_HDEPTH_SCALING
    integer :: NCAR_ORO_PGWV, NCAR_BKG_PGWV
    real    :: NCAR_ORO_GW_DC, NCAR_BKG_GW_DC
    real    :: NCAR_ORO_FCRIT2, NCAR_BKG_FCRIT2
    real    :: NCAR_ORO_WAVELENGTH, NCAR_BKG_WAVELENGTH
    real    :: NCAR_ORO_SOUTH_FAC
    real    :: NCAR_ORO_TNDMAX
    real    :: NCAR_BKG_TNDMAX
    real    :: NCAR_HR_CF      ! Grid cell convective conversion factor
    real    :: NCAR_ET_TAUBGND ! Extratropical background frontal forcing
    logical :: NCAR_DC_BERES
    logical :: NCAR_SC_BERES
    integer :: GEOS_PGWV
    real :: NCAR_EFFGWBKG
    real :: NCAR_DC_BERES_SRC_LEVEL, NCAR_SC_BERES_SRC_LEVEL

    type (wrap_) :: wrap
    type (GEOS_GwdGridComp), pointer        :: self
    integer :: num_threads, thread

    type(MAPL_Interval), allocatable :: bounds(:)
    integer :: JM_thread

!=============================================================================

   ! Begin...

   ! Get my name and set-up traceback handle
   ! ---------------------------------------

      Iam = 'Initialize'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC )
      Iam = trim(COMP_NAME) // Iam

      ! Get my internal MAPL_Generic state
      !-----------------------------------

      call MAPL_GetObjectFromGC ( GC, MAPL, _RC )
      
!   Get my internal private state
!   -----------------------------
      call ESMF_UserCompGetInternalState(GC, 'GEOS_GwdGridComp', wrap, _RC)
      self => wrap%ptr

      ! Call Generic Initialize for GWD GC
      !-----------------------------------

      call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK, _RC )
      
      call MAPL_Get(MAPL, IM=IM, JM=JM, LM=LM, LATS=LATS, _RC)
      
     ! Get grid name to determine IMSIZE
      call MAPL_GetResource(MAPL,GRIDNAME,'AGCM_GRIDNAME:', _RC)
            GRIDNAME =  AdjustL(GRIDNAME)
      nn = len_trim(GRIDNAME)
      dateline = GRIDNAME(nn-1:nn)
      imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
      read(imchar,*) imsize
      if(dateline.eq.'CF') imsize = imsize*4

! Gravity wave drag
! -----------------

    if (LM .eq. 72) then
                                         GEOS_PGWV = 4
       call MAPL_GetResource( MAPL, self%GEOS_PGWV,     Label="GEOS_PGWV:",     default=GEOS_PGWV, _RC)
       call MAPL_GetResource( MAPL, self%GEOS_BGSTRESS, Label="GEOS_BGSTRESS:", default=0.900, _RC)
       call MAPL_GetResource( MAPL, self%GEOS_EFFGWBKG, Label="GEOS_EFFGWBKG:", default=0.125, _RC)
       call MAPL_GetResource( MAPL, self%GEOS_EFFGWORO, Label="GEOS_EFFGWORO:", default=0.250, _RC)
       call MAPL_GetResource( MAPL, self%NCAR_EFFGWBKG, Label="NCAR_EFFGWBKG:", default=0.000, _RC)
       call MAPL_GetResource( MAPL, self%NCAR_EFFGWORO, Label="NCAR_EFFGWORO:", default=0.000, _RC)
       call MAPL_GetResource( MAPL, self%NCAR_NRDG,     Label="NCAR_NRDG:",     default=0, _RC)
       call MAPL_GetResource( MAPL, self%HGT_SURFACE,   Label="HGT_SURFACE:",   default=0.0, _RC)
       call MAPL_GetResource( MAPL, self%TAU1,          Label="RAYLEIGH_TAU1:", default=172800., _RC)
    else
                                         GEOS_PGWV = NINT(32*LM/181.0)
       call MAPL_GetResource( MAPL, self%GEOS_PGWV,     Label="GEOS_PGWV:",     default=GEOS_PGWV, _RC)
       call MAPL_GetResource( MAPL, self%GEOS_BGSTRESS, Label="GEOS_BGSTRESS:", default=0.000, _RC)
       call MAPL_GetResource( MAPL, self%GEOS_EFFGWBKG, Label="GEOS_EFFGWBKG:", default=0.125, _RC)
       call MAPL_GetResource( MAPL, self%GEOS_EFFGWORO, Label="GEOS_EFFGWORO:", default=0.000, _RC)
       call MAPL_GetResource( MAPL, self%NCAR_EFFGWBKG, Label="NCAR_EFFGWBKG:", default=1.000, _RC)
       call MAPL_GetResource( MAPL, self%NCAR_EFFGWORO, Label="NCAR_EFFGWORO:", default=1.000, _RC)
       call MAPL_GetResource( MAPL, self%NCAR_NRDG,     Label="NCAR_NRDG:",     default=16, _RC)
       call MAPL_GetResource( MAPL, self%HGT_SURFACE,   Label="HGT_SURFACE:",   default=50.0, _RC)
       call MAPL_GetResource( MAPL, self%TAU1,          Label="RAYLEIGH_TAU1:", default=0.00, _RC)
    endif

! Topographic Form Drag [Beljaars et al (2004)]
! ---------------------------------------------
      call MAPL_GetResource( MAPL, self%effbeljaars, Label="BELJAARS_EFF_FACTOR:",  default=0.0, _RC)
      call MAPL_GetResource( MAPL, self%limbeljaars, Label="BELJAARS_LIMITER:",  default=400.0, _RC)
                                   self%limbeljaars = self%limbeljaars/86400.0

! Rayleigh friction
! -----------------
      call MAPL_GetResource( MAPL, self%Z1,   Label="RAYLEIGH_Z1:",   default=75000.,  _RC)
      call MAPL_GetResource( MAPL, self%H0,   Label="RAYLEIGH_H0:",   default=7000.,   _RC)
      call MAPL_GetResource( MAPL, self%HH,   Label="RAYLEIGH_HH:",   default=7500.,   _RC)

! NCAR GWD settings
! -----------------
      call MAPL_GetResource( MAPL, NCAR_TAU_TOP_ZERO, Label="NCAR_TAU_TOP_ZERO:", default=.true., _RC)
      call MAPL_GetResource( MAPL, NCAR_PRNDL, Label="NCAR_PRNDL:", default=0.50, _RC)
                                   NCAR_QBO_HDEPTH_SCALING = min( imsize/1440.0 , 1.0 )
      call MAPL_GetResource( MAPL, NCAR_QBO_HDEPTH_SCALING, Label="NCAR_QBO_HDEPTH_SCALING:", default=NCAR_QBO_HDEPTH_SCALING, _RC)
                                   NCAR_HR_CF = max( 20.0*720.0/imsize , 1.0 )
      call MAPL_GetResource( MAPL, NCAR_HR_CF, Label="NCAR_HR_CF:", default=NCAR_HR_CF, _RC)
         
      call gw_common_init( NCAR_TAU_TOP_ZERO , 1 , &
                           MAPL_GRAV , &
                           MAPL_RGAS , &
                           MAPL_CP , &
                           NCAR_PRNDL, NCAR_QBO_HDEPTH_SCALING, NCAR_HR_CF, ERRstring )

      ! Beres Scheme File
      call MAPL_GetResource( MAPL, BERES_FILE_NAME, Label="BERES_FILE_NAME:", &
            default='ExtData/g5gcm/gwd/newmfspectra40_dc25.nc', _RC)
      call MAPL_GetResource( MAPL, NCAR_BKG_PGWV,       Label="NCAR_BKG_PGWV:",       default=32,    _RC)
      call MAPL_GetResource( MAPL, NCAR_BKG_GW_DC,      Label="NCAR_BKG_GW_DC:",      default=2.5,   _RC)
      call MAPL_GetResource( MAPL, NCAR_BKG_FCRIT2,     Label="NCAR_BKG_FCRIT2:",     default=1.0,   _RC)
      call MAPL_GetResource( MAPL, NCAR_BKG_WAVELENGTH, Label="NCAR_BKG_WAVELENGTH:", default=1.e5,  _RC)
      call MAPL_GetResource( MAPL, NCAR_ET_TAUBGND,     Label="NCAR_ET_TAUBGND:",     default=50.0,  _RC)
      call MAPL_GetResource( MAPL, NCAR_BKG_TNDMAX,     Label="NCAR_BKG_TNDMAX:",     default=500.0, _RC)
      NCAR_BKG_TNDMAX = NCAR_BKG_TNDMAX/86400.0
                 ! Beres DeepCu
      call MAPL_GetResource( MAPL, NCAR_DC_BERES_SRC_LEVEL, "NCAR_DC_BERES_SRC_LEVEL:", DEFAULT=70000.0, _RC)
      call MAPL_GetResource( MAPL, NCAR_DC_BERES, "NCAR_DC_BERES:", DEFAULT=.TRUE., _RC)
      num_threads = MAPL_get_num_threads()
      bounds = MAPL_find_bounds(JM, num_threads)
      do thread = 0, num_threads-1
            JM_thread = bounds(thread+1)%max - bounds(thread+1)%min + 1
            call gw_beres_init( BERES_FILE_NAME ,  &
                                self%workspaces(thread)%beres_band, &
                                self%workspaces(thread)%beres_dc_desc, &
                                NCAR_BKG_PGWV, NCAR_BKG_GW_DC, NCAR_BKG_FCRIT2, &
                                NCAR_BKG_WAVELENGTH, NCAR_DC_BERES_SRC_LEVEL, &
                                1000.0, .TRUE., NCAR_ET_TAUBGND, NCAR_BKG_TNDMAX, NCAR_DC_BERES, &
                                IM*JM_thread, LATS(:,bounds(thread+1)%min:bounds(thread+1)%max))
      end do
      ! Beres ShallowCu
      call MAPL_GetResource( MAPL, NCAR_SC_BERES_SRC_LEVEL, "NCAR_SC_BERES_SRC_LEVEL:", DEFAULT=90000.0, _RC)
      call MAPL_GetResource( MAPL, NCAR_SC_BERES, "NCAR_SC_BERES:", DEFAULT=.FALSE., _RC)
      do thread = 0, num_threads-1
            JM_thread = bounds(thread+1)%max - bounds(thread+1)%min + 1
            call gw_beres_init( BERES_FILE_NAME ,  &
                                self%workspaces(thread)%beres_band,  &
                                self%workspaces(thread)%beres_sc_desc,  &
                                NCAR_BKG_PGWV, NCAR_BKG_GW_DC, NCAR_BKG_FCRIT2,  &
                                NCAR_BKG_WAVELENGTH, NCAR_SC_BERES_SRC_LEVEL, &
                                0.0, .FALSE., NCAR_ET_TAUBGND, NCAR_BKG_TNDMAX, NCAR_SC_BERES, &
                                IM*JM_thread, LATS(:,bounds(thread+1)%min:bounds(thread+1)%max))
      end do

      ! Orographic Scheme
      call MAPL_GetResource( MAPL, NCAR_ORO_PGWV,       Label="NCAR_ORO_PGWV:",       default=0,    _RC)
      call MAPL_GetResource( MAPL, NCAR_ORO_GW_DC,      Label="NCAR_ORO_GW_DC:",      default=2.5,  _RC)
      call MAPL_GetResource( MAPL, NCAR_ORO_FCRIT2,     Label="NCAR_ORO_FCRIT2:",     default=1.0,  _RC)
      call MAPL_GetResource( MAPL, NCAR_ORO_WAVELENGTH, Label="NCAR_ORO_WAVELENGTH:", default=1.e5, _RC)
      call MAPL_GetResource( MAPL, NCAR_ORO_SOUTH_FAC,  Label="NCAR_ORO_SOUTH_FAC:",  default=2.0,  _RC)
      do thread = 0, num_threads-1
            call gw_oro_init ( self%workspaces(thread)%oro_band, NCAR_ORO_GW_DC, &
                               NCAR_ORO_FCRIT2, NCAR_ORO_WAVELENGTH, NCAR_ORO_PGWV, &
                               NCAR_ORO_SOUTH_FAC )
      end do
      ! Ridge Scheme
      if (self%NCAR_NRDG > 0) then
          call MAPL_GetResource( MAPL, NCAR_ORO_TNDMAX,   Label="NCAR_ORO_TNDMAX:",  default=200.0, _RC)
          NCAR_ORO_TNDMAX = NCAR_ORO_TNDMAX/86400.0
          do thread = 0, num_threads-1
             call gw_rdg_init ( self%workspaces(thread)%rdg_band, NCAR_ORO_GW_DC, NCAR_ORO_FCRIT2, NCAR_ORO_WAVELENGTH, NCAR_ORO_TNDMAX, NCAR_ORO_PGWV )
          end do
      endif

      allocate(self%alpha(LM+1), _STAT)
      call MAPL_GetPointer( IMPORT, PREF,     'PREF',    _RC )
      call gw_newtonian_set(LM, PREF, self%alpha)

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
  type (ESMF_Grid        )            :: ESMFGRID

  integer                             :: IM, JM, LM
  !integer                             :: pgwv
  !real                                :: HGT_SURFACE
  !real                                :: effbeljaars, limbeljaars, tcrib
  real                                :: tcrib
  !real                                :: effgworo, effgwbkg
  !real                                :: CDMBGWD1, CDMBGWD2
  !real                                :: bgstressmax
  real, pointer, dimension(:,:)       :: LATS

! Rayleigh friction parameters

  REAL                                :: H0, HH, Z1, TAU1

  type (wrap_) :: wrap
  type (GEOS_GwdGridComp), pointer        :: self
  type(ThreadWorkspace), pointer :: workspace
  integer :: thread

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

   Iam = "Run"
   !call ESMF_GridCompGet( GC, name=COMP_NAME, grid=ESMFGRID, _RC )
   call ESMF_GridCompGet( GC, name=COMP_NAME, _RC )
   Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
!----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, _RC)
   
!   Get my internal private state
!   -----------------------------
    call ESMF_UserCompGetInternalState(GC, 'GEOS_GwdGridComp', wrap, _RC)
    self => wrap%ptr

    H0 = self%H0
    HH = self%HH
    Z1 = self%Z1
    TAU1 = self%TAU1

! Local aliases to the state, grid, and configuration
! ---------------------------------------------------

   !call MAPL_TimerOn(MAPL,"TOTAL")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL, &
         IM=IM, JM=JM, LM=LM,        &
         RUNALARM=ALARM, LATS=LATS,  &
                           _RC )
    
! If its time, recalculate the GWD tendency
! -----------------------------------------

   if ( ESMF_AlarmIsRinging( ALARM ) ) then
      !call ESMF_AlarmRingerOff(ALARM, _RC)
      !call MAPL_TimerOn (MAPL,"DRIVER")
      call Gwd_Driver(_RC)
      !call MAPL_TimerOff(MAPL,"DRIVER")
   endif

   !call MAPL_TimerOff(MAPL,"TOTAL")

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
      real, pointer, dimension(:,:)    :: AREA, SGH, VARFLT, PHIS
      real, pointer, dimension(:,:,:)  :: PLE, T, Q, U, V
      !++jtb Array for moist deep & shallow conv heating
      real, pointer, dimension(:,:,:)  :: HT_dc, HT_sc
      ! Arrays for QL and QI condensate tendencies from Moist
      real, pointer, dimension(:,:,:)  :: QLDT_mst, QIDT_mst
      !++jtb pointers for NCAR Orographic GWP
      !     (in Internal State)
      real, pointer, dimension(:,:,:)  :: MXDIS
      real, pointer, dimension(:,:,:)  :: CLNGT
      real, pointer, dimension(:,:,:)  :: HWDTH
      real, pointer, dimension(:,:,:)  :: ANGLL
      real, pointer, dimension(:,:,:)  :: ANIXY
      real, pointer, dimension(:,:)    :: GBXAR
      real, pointer, dimension(:,:,:)  :: KWVRDG
      real, pointer, dimension(:,:,:)  :: EFFRDG

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
      real, pointer, dimension(:,:,:)  ::           DUDT_TFD, DVDT_TFD
      real, pointer, dimension(:,:,:)  :: DTDT_ORO, DUDT_ORO, DVDT_ORO
      real, pointer, dimension(:,:,:)  :: DTDT_BKG, DUDT_BKG, DVDT_BKG
      real, pointer, dimension(:,:,:)  :: DTDT_RAY, DUDT_RAY, DVDT_RAY
      real, pointer, dimension(:,:,:)  :: DTGENBKG, DUGENBKG, DVGENBKG

      real, pointer, dimension(:,:,:)  :: TMP3D
      real, pointer, dimension(:,:)    :: TMP2D

! local variables

      real,              dimension(IM,JM,LM  ) :: ZM, PMID, PDEL, RPDEL, PMLN
      real,              dimension(IM,JM     ) :: a2, Hefold
      real,              dimension(IM,JM,LM  ) :: DUDT_TOFD, DVDT_TOFD
      real,              dimension(IM,JM,LM  ) :: DUDT_ORG, DVDT_ORG, DTDT_ORG
      real,              dimension(IM,JM,LM  ) :: DUDT_GWD, DVDT_GWD, DTDT_GWD
      real,              dimension(IM,JM,LM  ) :: DUDT_RAH, DVDT_RAH, DTDT_RAH
      real,              dimension(IM,JM,LM  ) :: DUDT_TOT, DVDT_TOT, DTDT_TOT
      real,              dimension(IM,JM,LM+1) :: PILN,   ZI
      real,              dimension(      LM  ) :: ZREF, KRAY
      real,              dimension(IM,JM     ) :: GBXAR_TMP
      real,              dimension(IM,JM     ) :: TAUXO_TMP, TAUYO_TMP
      real,              dimension(IM,JM     ) :: TAUXB_TMP, TAUYB_TMP
      real,              dimension(IM,JM,LM+1) :: TAUXO_3D , TAUYO_3D , FEO_3D, FEPO_3D
      real,              dimension(IM,JM,LM+1) :: TAUXB_3D , TAUYB_3D , FEB_3D, FEPB_3D
      real,              dimension(IM,JM,LM  ) :: DUBKGSRC , DVBKGSRC , DTBKGSRC
      real,              dimension(IM,JM)      :: KEGWD_X, KEORO_X,  KERAY_X,  KEBKG_X, KERES_X
      real,              dimension(IM,JM)      :: PEGWD_X, PEORO_X,  PERAY_X,  PEBKG_X, BKGERR_X

      real,              dimension(IM,JM,LM  ) :: DUDT_GWD_GEOS , DVDT_GWD_GEOS , DTDT_GWD_GEOS
      real,              dimension(IM,JM,LM  ) :: DUDT_ORG_GEOS , DVDT_ORG_GEOS , DTDT_ORG_GEOS
      real,              dimension(IM,JM     ) :: TAUXB_TMP_GEOS, TAUYB_TMP_GEOS
      real,              dimension(IM,JM     ) :: TAUXO_TMP_GEOS, TAUYO_TMP_GEOS

      real,              dimension(IM,JM,LM  ) :: DUDT_GWD_NCAR , DVDT_GWD_NCAR , DTDT_GWD_NCAR
      real,              dimension(IM,JM,LM  ) :: DUDT_ORG_NCAR , DVDT_ORG_NCAR , DTDT_ORG_NCAR
      real,              dimension(IM,JM     ) :: TAUXB_TMP_NCAR, TAUYB_TMP_NCAR
      real,              dimension(IM,JM     ) :: TAUXO_TMP_NCAR, TAUYO_TMP_NCAR

      real,              dimension(IM,JM     ) :: DC_SRC_L, SC_SRC_L

      integer                                  :: J, K, L, nrdg, ikpbl
      real(ESMF_KIND_R8)                       :: DT_R8
      real                                     :: DT     ! time interval in sec
      real                                     :: a1, wsp, var_temp
      !real, allocatable :: THV(:,:,:)
      real :: THV(IM,JM,LM)

      integer           :: I,IRUN
      type (ESMF_State) :: INTERNAL

!  Begin...
!----------

      IAm = "Gwd_Driver"

! Get time step
!-------------------------------------------------

      call ESMF_AlarmGet( ALARM, ringInterval=TINT,    _RC)
      call ESMF_TimeIntervalGet(TINT, S_R8=DT_R8,      _RC)
       
      DT = DT_R8

! Pointers to inputs
!---------------------

      call MAPL_GetPointer( IMPORT, PLE,      'PLE',     _RC )
      call MAPL_GetPointer( IMPORT, T,        'T',       _RC )
      call MAPL_GetPointer( IMPORT, Q,        'Q',       _RC )
      call MAPL_GetPointer( IMPORT, U,        'U',       _RC )
      call MAPL_GetPointer( IMPORT, V,        'V',       _RC )
      call MAPL_GetPointer( IMPORT, PHIS,     'PHIS',    _RC )
      call MAPL_GetPointer( IMPORT, SGH,      'SGH',     _RC )
      call MAPL_GetPointer( IMPORT, PREF,     'PREF',    _RC )
      call MAPL_GetPointer( IMPORT, AREA,     'AREA',    _RC )
      call MAPL_GetPointer( IMPORT, VARFLT,   'VARFLT',  _RC )
      call MAPL_GetPointer( IMPORT, HT_dc,    'DTDT_DC', _RC )
      call MAPL_GetPointer( IMPORT, HT_sc,    'DTDT_SC', _RC )
      call MAPL_GetPointer( IMPORT, QLDT_mst, 'DQLDT'  , _RC )
      call MAPL_GetPointer( IMPORT, QIDT_mst, 'DQIDT'  , _RC )
 
! Allocate/refer to the outputs
!------------------------------

      call MAPL_GetPointer(EXPORT,  PLE_EXP, 'PLE'     , _RC)
      call MAPL_GetPointer(EXPORT,    T_EXP, 'T'       , _RC)
      call MAPL_GetPointer(EXPORT,    Q_EXP, 'Q'       , _RC)
      call MAPL_GetPointer(EXPORT,    U_EXP, 'U'       , _RC)
      call MAPL_GetPointer(EXPORT,    V_EXP, 'V'       , _RC)
      call MAPL_GetPointer(EXPORT,  SGH_EXP, 'SGH'     , _RC)
      call MAPL_GetPointer(EXPORT, PREF_EXP, 'PREF'    , _RC)
      call MAPL_GetPointer(EXPORT,    TTMGW, 'TTMGW'   , _RC)
      call MAPL_GetPointer(EXPORT, DUDT_TFD, 'DUDT_TFD', _RC)
      call MAPL_GetPointer(EXPORT, DVDT_TFD, 'DVDT_TFD', _RC)
      call MAPL_GetPointer(EXPORT, DTDT_ORO, 'DTDT_ORO', _RC)
      call MAPL_GetPointer(EXPORT, DUDT_ORO, 'DUDT_ORO', _RC)
      call MAPL_GetPointer(EXPORT, DVDT_ORO, 'DVDT_ORO', _RC)
      call MAPL_GetPointer(EXPORT, DTDT_BKG, 'DTDT_BKG', _RC)
      call MAPL_GetPointer(EXPORT, DUDT_BKG, 'DUDT_BKG', _RC)
      call MAPL_GetPointer(EXPORT, DVDT_BKG, 'DVDT_BKG', _RC)
      call MAPL_GetPointer(EXPORT, DTDT_RAY, 'DTDT_RAY', _RC)
      call MAPL_GetPointer(EXPORT, DUDT_RAY, 'DUDT_RAY', _RC)
      call MAPL_GetPointer(EXPORT, DVDT_RAY, 'DVDT_RAY', _RC)
      call MAPL_GetPointer(EXPORT,   TAUGWX, 'TAUGWX'  , _RC)
      call MAPL_GetPointer(EXPORT,   TAUGWY, 'TAUGWY'  , _RC)
      call MAPL_GetPointer(EXPORT,  TAUOROX, 'TAUOROX' , _RC)
      call MAPL_GetPointer(EXPORT,  TAUOROY, 'TAUOROY' , _RC)
      call MAPL_GetPointer(EXPORT,  TAUBKGX, 'TAUBKGX' , _RC)
      call MAPL_GetPointer(EXPORT,  TAUBKGY, 'TAUBKGY' , _RC)
      call MAPL_GetPointer(EXPORT,  TAUMSTX, 'TAUMSTX' , _RC)
      call MAPL_GetPointer(EXPORT,  TAUMSTY, 'TAUMSTY' , _RC)
      call MAPL_GetPointer(EXPORT,    UBASE, 'UBASE'   , _RC)
      call MAPL_GetPointer(EXPORT,    VBASE, 'VBASE'   , _RC)
      call MAPL_GetPointer(EXPORT,     UBAR, 'UBAR'    , _RC)
      call MAPL_GetPointer(EXPORT,     VBAR, 'VBAR'    , _RC)
      call MAPL_GetPointer(EXPORT,   CLDSTD, 'CLDSTD'  , _RC)
       
      call MAPL_GetPointer(EXPORT,     DTDT, 'DTDT'    , _RC)
      call MAPL_GetPointer(EXPORT,     DUDT, 'DUDT'    , _RC)
      call MAPL_GetPointer(EXPORT,     DVDT, 'DVDT'    , _RC)
       
      call MAPL_GetPointer(EXPORT,    PEGWD, 'PEGWD'   , _RC)
      call MAPL_GetPointer(EXPORT,    PEORO, 'PEORO'   , _RC)
      call MAPL_GetPointer(EXPORT,    PERAY, 'PERAY'   , _RC)
      call MAPL_GetPointer(EXPORT,    PEBKG, 'PEBKG'   , _RC)
       
      call MAPL_GetPointer(EXPORT,    KEGWD, 'KEGWD'   , _RC)
      call MAPL_GetPointer(EXPORT,    KEORO, 'KEORO'   , _RC)
      call MAPL_GetPointer(EXPORT,    KERAY, 'KERAY'   , _RC)
      call MAPL_GetPointer(EXPORT,    KEBKG, 'KEBKG'   , _RC)
      call MAPL_GetPointer(EXPORT,    KERES, 'KERES'   , _RC)
      call MAPL_GetPointer(EXPORT,   BKGERR, 'BKGERR'  , _RC)
       

      CALL PREGEO(IM*JM,   LM,   &
             PLE, LATS,   PMID,  PDEL, RPDEL,     PILN,     PMLN)

! Compute ZM
!-------------

      call GEOPOTENTIAL( IM*JM, LM,                  &
           PILN,  PMLN,  PLE,   PMID, PDEL, RPDEL,   &
           T,     Q,     ZI,    ZM                   )

! Do gravity wave drag calculations on a list of soundings
!---------------------------------------------------------

    !call MAPL_TimerOn(MAPL,"-INTR")

         ! get pointers from INTERNAL:MXDIS
         call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, _RC)
         call MAPL_GetPointer( INTERNAL, MXDIS, 'MXDIS', _RC )
         call MAPL_GetPointer( INTERNAL, HWDTH, 'HWDTH', _RC )
         call MAPL_GetPointer( INTERNAL, CLNGT, 'CLNGT', _RC )
         call MAPL_GetPointer( INTERNAL, ANGLL, 'ANGLL', _RC )
         call MAPL_GetPointer( INTERNAL, ANIXY, 'ANIXY', _RC )
         call MAPL_GetPointer( INTERNAL, GBXAR, 'GBXAR', _RC )
         call MAPL_GetPointer( INTERNAL, KWVRDG, 'KWVRDG', _RC )
         call MAPL_GetPointer( INTERNAL, EFFRDG, 'EFFRDG', _RC )
         
         GBXAR_TMP = GBXAR * (MAPL_RADIUS/1000.)**2 ! transform to km^2
         WHERE (ANGLL < -180)
           ANGLL = 0.0
         END WHERE

         do nrdg = 1, self%NCAR_NRDG
           KWVRDG(:,:,nrdg) = 0.001/(HWDTH(:,:,nrdg)+0.001)
           EFFRDG(:,:,nrdg) = self%NCAR_EFFGWORO*(HWDTH(:,:,nrdg)*CLNGT(:,:,nrdg))/GBXAR_TMP
         enddo

!         if (FIRST_RUN) then
!           FIRST_RUN = .false.
!           call gw_newtonian_set(LM, PREF)
!!#ifdef DEBUG_GWD
!           if (self%NCAR_NRDG > 0) then
!            IF (MAPL_AM_I_ROOT()) write(*,*) 'GWD internal state: '
!            call Write_Profile(GBXAR_TMP,         AREA, ESMFGRID, 'GBXAR')
!            do nrdg = 1, self%NCAR_NRDG
!             IF (MAPL_AM_I_ROOT()) write(*,*) 'NRDG: ', nrdg
!             call Write_Profile(MXDIS(:,:,nrdg),  AREA, ESMFGRID, 'MXDIS')
!             call Write_Profile(ANGLL(:,:,nrdg),  AREA, ESMFGRID, 'ANGLL')
!             call Write_Profile(ANIXY(:,:,nrdg),  AREA, ESMFGRID, 'ANIXY')
!             call Write_Profile(CLNGT(:,:,nrdg),  AREA, ESMFGRID, 'CLNGT')
!             call Write_Profile(HWDTH(:,:,nrdg),  AREA, ESMFGRID, 'HWDTH')
!             call Write_Profile(KWVRDG(:,:,nrdg), AREA, ESMFGRID, 'KWVRDG')
!             call Write_Profile(EFFRDG(:,:,nrdg), AREA, ESMFGRID, 'EFFRDG')
!            enddo
!          endif
!!#endif
!         endif

         call MAPL_GetPointer(EXPORT, TMP2D, 'RDG1_MXDIS', _RC)
         if(associated(TMP2D)) TMP2D = MXDIS(:,:,1)
         call MAPL_GetPointer(EXPORT, TMP2D, 'RDG1_HWDTH', _RC)
         if(associated(TMP2D)) TMP2D = HWDTH(:,:,1)
         call MAPL_GetPointer(EXPORT, TMP2D, 'RDG1_CLNGT', _RC)
         if(associated(TMP2D)) TMP2D = CLNGT(:,:,1)
         call MAPL_GetPointer(EXPORT, TMP2D, 'RDG1_ANGLL', _RC)
         if(associated(TMP2D)) TMP2D = ANGLL(:,:,1)
         call MAPL_GetPointer(EXPORT, TMP2D, 'RDG1_ANIXY', _RC)
         if(associated(TMP2D)) TMP2D = ANIXY(:,:,1)
         call MAPL_GetPointer(EXPORT, TMP2D, 'RDG1_GBXAR', _RC)
         if(associated(TMP2D)) TMP2D = GBXAR_TMP

         ! Use new NCAR code convective+oro (excludes extratropical bkg sources)
         DUDT_GWD_NCAR = 0.0
         DVDT_GWD_NCAR = 0.0
         DTDT_GWD_NCAR = 0.0
         TAUXB_TMP_NCAR = 0.0
         TAUYB_TMP_NCAR = 0.0
         DUDT_ORG_NCAR = 0.0
         DVDT_ORG_NCAR = 0.0
         DTDT_ORG_NCAR = 0.0
         TAUXO_TMP_NCAR = 0.0
         TAUYO_TMP_NCAR = 0.0
         !call MAPL_TimerOn(MAPL,"-INTR_NCAR")
         if ( (self%NCAR_EFFGWORO /= 0.0) .OR. (self%NCAR_EFFGWBKG /= 0.0) ) then
            thread = MAPL_get_current_thread()
            workspace => self%workspaces(thread)
            call gw_intr_ncar(IM*JM,    LM,         DT,     self%NCAR_NRDG,   &
                 workspace%beres_dc_desc, workspace%beres_sc_desc, &
                 workspace%beres_band, workspace%oro_band, workspace%rdg_band, &
                 PLE,       T,          U,          V,                   &
                 HT_dc,     HT_sc,      QLDT_mst+QIDT_mst,               &
                 SGH,       MXDIS,      HWDTH,      CLNGT,  ANGLL,       &
                 ANIXY,     GBXAR_TMP,  KWVRDG,     EFFRDG, PREF,        &
                 PMID,      PDEL,       RPDEL,      PILN,   ZM,    LATS, &
                 PHIS,                                                   &
                 DUDT_GWD_NCAR,  DVDT_GWD_NCAR,   DTDT_GWD_NCAR,         &
                 DUDT_ORG_NCAR,  DVDT_ORG_NCAR,   DTDT_ORG_NCAR,         &
                 TAUXO_TMP_NCAR, TAUYO_TMP_NCAR,  &
                 TAUXB_TMP_NCAR, TAUYB_TMP_NCAR,  &
                 self%NCAR_EFFGWORO, &
                 self%NCAR_EFFGWBKG, self%alpha, &
                 _RC)
         endif
         !call MAPL_TimerOff(MAPL,"-INTR_NCAR")

         ! Use GEOS GWD only for Extratropical background sources...
         DUDT_GWD_GEOS = 0.0
         DVDT_GWD_GEOS = 0.0
         DTDT_GWD_GEOS = 0.0
         TAUXB_TMP_GEOS = 0.0
         TAUYB_TMP_GEOS = 0.0
         DUDT_ORG_GEOS = 0.0
         DVDT_ORG_GEOS = 0.0
         DTDT_ORG_GEOS = 0.0
         TAUXO_TMP_GEOS = 0.0
         TAUYO_TMP_GEOS = 0.0
         !call MAPL_TimerOn(MAPL,"-INTR_GEOS")
         if ( (self%GEOS_EFFGWORO /= 0.0) .OR. (self%GEOS_EFFGWBKG /= 0.0) ) then
          call gw_intr   (IM*JM,      LM,         DT,                  &
               self%GEOS_PGWV,                                              &
               PLE,       T,          U,          V,      SGH,   PREF, &
               PMID,      PDEL,       RPDEL,      PILN,   ZM,    LATS, &
               DUDT_GWD_GEOS,  DVDT_GWD_GEOS,   DTDT_GWD_GEOS,         &
               DUDT_ORG_GEOS,  DVDT_ORG_GEOS,   DTDT_ORG_GEOS,         &
               TAUXO_TMP_GEOS, TAUYO_TMP_GEOS,  TAUXO_3D,   TAUYO_3D,  FEO_3D,   &
               TAUXB_TMP_GEOS, TAUYB_TMP_GEOS,  TAUXB_3D,   TAUYB_3D,  FEB_3D,   &
               FEPO_3D,   FEPB_3D,    DUBKGSRC,   DVBKGSRC,  DTBKGSRC, &
               self%GEOS_BGSTRESS, &
               self%GEOS_EFFGWORO, &
               self%GEOS_EFFGWBKG, &
               _RC)
         endif
         !call MAPL_TimerOff(MAPL,"-INTR_GEOS")

         ! Total
         DUDT_GWD=DUDT_GWD_GEOS+DUDT_GWD_NCAR
         DVDT_GWD=DVDT_GWD_GEOS+DVDT_GWD_NCAR
         DTDT_GWD=DTDT_GWD_GEOS+DTDT_GWD_NCAR
         ! Background
         TAUXB_TMP=TAUXB_TMP_GEOS+TAUXB_TMP_NCAR
         TAUYB_TMP=TAUYB_TMP_GEOS+TAUYB_TMP_NCAR
         ! Orographic
         DUDT_ORG=DUDT_ORG_GEOS+DUDT_ORG_NCAR
         DVDT_ORG=DVDT_ORG_GEOS+DVDT_ORG_NCAR
         DTDT_ORG=DTDT_ORG_GEOS+DTDT_ORG_NCAR
         TAUXO_TMP=TAUXO_TMP_GEOS+TAUXO_TMP_NCAR
         TAUYO_TMP=TAUYO_TMP_GEOS+TAUYO_TMP_NCAR
    !call MAPL_TimerOff(MAPL,"-INTR")

    !call MAPL_TimerOn(MAPL,"-BELJAARS_TOFD")
    if (self%effbeljaars > 0.0) then
        THV = T * (1.0 + MAPL_VIREPS * Q) / ( (PMID/MAPL_P00)**MAPL_KAPPA )
    DO J=1,JM
       DO I=1,IM
! Find the PBL height
             ikpbl = LM
             do L=LM-1,1,-1
                tcrib = MAPL_GRAV*(THV(I,J,L)-THV(I,J,LM))*ZM(I,J,L)/ &
                        (THV(I,J,LM)*MAX(U(I,J,L)**2+V(I,J,L)**2,1.0E-8))
                if (tcrib >= 0.25) then
                   ikpbl = L
                   exit
                end if
             end do
! determine the efolding height
             a2(i,j)=self%effbeljaars * 1.08371722e-7 * VARFLT(i,j) * &
                     MAX(0.0,MIN(1.0,dxmax_ss*(1.-dxmin_ss/SQRT(AREA(i,j))/(dxmax_ss-dxmin_ss))))
           ! Revise e-folding height based on PBL height and topographic std. dev.
             Hefold(i,j) = 1500.0 !MIN(MAX(2*SQRT(VARFLT(i,j)),ZM(i,j,ikpbl)),1500.)
       END DO
    END DO
    DO L=1, LM
       DO J=1,JM
          DO I=1,IM
               var_temp = 0.0
               if (a2(i,j) > 0.0 .AND. ZM(I,J,L) < 4.0*Hefold(i,j)) then
                  wsp      = SQRT(U(i,j,l)**2 + V(i,j,l)**2)
                  wsp      = SQRT(MIN(wsp/25.0,1.0))*MAX(25.0,wsp) ! enhance winds below 25 m/s
                  var_temp = ZM(I,J,L)/Hefold(i,j)
                  var_temp = exp(-var_temp*sqrt(var_temp))*(var_temp**(-1.2))
                  var_temp = wsp*a2(i,j)*(var_temp/Hefold(i,j))
                 !  Note:  This is a semi-implicit treatment of the time differencing
                 !  per Beljaars et al. (2004, QJRMS) doi: 10.1256/qj.03.73
                  DUDT_TOFD(i,j,l) = - var_temp*U(i,j,l)/(1. + var_temp*DT)
                  DVDT_TOFD(i,j,l) = - var_temp*V(i,j,l)/(1. + var_temp*DT)
                 ! Apply Tendency Limiter
                  if (abs(DUDT_TOFD(i,j,l)) > self%limbeljaars) then
                    DUDT_TOFD(i,j,l) = (self%limbeljaars/abs(DUDT_TOFD(i,j,l))) * DUDT_TOFD(i,j,l)
                  end if
                  if (abs(DVDT_TOFD(i,j,l)) > self%limbeljaars) then
                    DVDT_TOFD(i,j,l) = (self%limbeljaars/abs(DVDT_TOFD(i,j,l))) * DVDT_TOFD(i,j,l)
                  end if
               else
                  DUDT_TOFD(i,j,l) = 0.0
                  DVDT_TOFD(i,j,l) = 0.0
               end if
          END DO
       END DO
    END DO
    DUDT_GWD=DUDT_GWD+DUDT_TOFD
    DVDT_GWD=DVDT_GWD+DVDT_TOFD
    !deallocate( THV )
    else
    DUDT_TOFD=0.0
    DVDT_TOFD=0.0
    endif
    !call MAPL_TimerOff(MAPL,"-BELJAARS_TOFD")

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
    if(associated(DTDT    )) DTDT     = DTDT_TOT*PDEL ! DTDT has to be pressure weighted for dynamics

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

    if(associated(DUDT_TFD)) DUDT_TFD = DUDT_TOFD
    if(associated(DVDT_TFD)) DVDT_TFD = DVDT_TOFD

    if(associated(DUDT_ORO)) DUDT_ORO = DUDT_ORG
    if(associated(DVDT_ORO)) DVDT_ORO = DVDT_ORG
    if(associated(DTDT_ORO)) DTDT_ORO = DTDT_ORG

    if(associated(DUDT_BKG)) DUDT_BKG = DUDT_GWD - DUDT_ORG - DUDT_TOFD
    if(associated(DVDT_BKG)) DVDT_BKG = DVDT_GWD - DVDT_ORG - DVDT_TOFD
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
    if(associated(TTMGW   )) TTMGW    = DTDT_TOT

! Fille additional exports
!-------------------------
    if(associated(    Q_EXP ))    Q_EXP = Q
    if(associated(    U_EXP ))    U_EXP = U + DUDT_TOT*DT
    if(associated(    V_EXP ))    V_EXP = V + DVDT_TOT*DT
    if(associated(    T_EXP ))    T_EXP = T + DTDT_TOT*DT
    if(associated( PREF_EXP )) PREF_EXP = PREF
    if(associated(  SGH_EXP ))  SGH_EXP = SGH
    if(associated(  PLE_EXP ))  PLE_EXP = PLE

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
        if (TAU1 > 0.0) then
          ZREF     = H0 * LOG(MAPL_P00/(0.5*(PREF(K)+PREF(K+1))))
          KRAY     = (1.0/TAU1)*( 1.0 - TANH( (Z1-ZREF)/HH ) )
          KRAY     = KRAY/(1+DT*KRAY)
          DUDT_RAH(I,K) = -U(I,K)*KRAY
          DVDT_RAH(I,K) = -V(I,K)*KRAY
          DTDT_RAH(I,K) = - ((U(I,K) + (0.5*DT)*DUDT_RAH(I,K))*DUDT_RAH(I,K) + &
                             (V(I,K) + (0.5*DT)*DVDT_RAH(I,K))*DVDT_RAH(I,K)   ) * (1.0/MAPL_CP)
        else
          DUDT_RAH(I,K) = 0.0
          DVDT_RAH(I,K) = 0.0
          DTDT_RAH(I,K) = 0.0
        endif

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

  Subroutine Write_Profile(avar, area, grid, name)
    type(ESMF_Grid),  intent(IN) :: grid
    real,             intent(IN) :: avar(:,:)
    real,             intent(IN) :: area(:,:)
    character(len=*), intent(IN) :: name

    real(kind=ESMF_KIND_R8), allocatable :: locArr(:,:)
    real(kind=ESMF_KIND_R8), allocatable :: glbArr(:,:)
    real, allocatable :: area_global(:,:)
    real, allocatable :: avar_global(:,:)
    real :: rng(3)
    integer :: DIMS(3), STATUS, rc

    call MAPL_GridGet(GRID, localCellCountPerDim=DIMS, _RC)
    allocate (      locArr(DIMS(1),DIMS(2)) )

    call MAPL_GridGet(GRID, globalCellCountPerDim=DIMS, _RC)
    allocate (      glbArr(DIMS(1),DIMS(2)) )
    allocate ( area_global(DIMS(1),DIMS(2)) )
    allocate ( avar_global(DIMS(1),DIMS(2)) )

#if 1
    locArr = avar
    call ArrayGather(locArr, glbArr, grid)
    avar_global = glbArr

    locArr = area
    call ArrayGather(locArr, glbArr, grid)
    area_global = glbArr

    IF (MAPL_AM_I_ROOT()) Then
       rng(1) = MINVAL(MINVAL(avar_global,DIM=1),DIM=1)
       rng(2) = MAXVAL(MAXVAL(avar_global,DIM=1),DIM=1)
       rng(3) = SUM(SUM(avar_global*area_global,DIM=1),DIM=1) / &
                SUM(SUM(            area_global,DIM=1),DIM=1)
       Write(*,'(A," ",3(f21.9,1x))'),trim(name),rng(:)
    End IF
#else
    rng(1) = MINVAL(MINVAL(avar,DIM=1),DIM=1)
    rng(2) = MAXVAL(MAXVAL(avar,DIM=1),DIM=1)
    rng(3) = SUM(SUM(avar*area,DIM=1),DIM=1) / &
             SUM(SUM(     area,DIM=1),DIM=1)
    Write(*,'(A," ",3(f21.9,1x))'),trim(name),rng(:)
#endif

    deallocate ( locArr )
    deallocate ( glbArr )
    deallocate ( area_global )
    deallocate ( avar_global )

  End Subroutine Write_Profile

end module GEOS_GwdGridCompMod
