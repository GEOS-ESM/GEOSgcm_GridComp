! $Id: GEOS_AgcmGridComp.F90,v 1.85.12.25.2.1.18.1.2.1.6.1.10.10.2.1 2019/11/18 21:20:23 ltakacs Exp $

#include "MAPL_Generic.h"

!#define PRINT_STATES
#define FULLPHYSICS
#define GCM
#define debug 0

! Held-Suarez is not in module GEOSGCM?????

#if defined HS
#define  GEOS_physicsGridCompMod GEOS_hsGridCompMod
#undef   FULLPHYSICS
#endif


#if defined SCM
#define GEOS_superdynGridCompMod GEOS_singcolGridCompMod
#undef  GCM
#endif


!=============================================================================
!BOP

! !MODULE: GEOS_AgcmGridCompMod -- A Module to combine Supedynamics and Physics Gridded Components

! !INTERFACE:

module GEOS_AgcmGridCompMod

! !USES:

  use ESMF
  use MAPL
  use GEOS_TopoGetMod

  use GEOS_superdynGridCompMod,  only:  SDYN_SetServices => SetServices
  use GEOS_physicsGridCompMod,   only:  PHYS_SetServices => SetServices
  use MAPL_OrbGridCompMod,       only:  ORB_SetServices => SetServices

  use GEOS_RemapMod, only: myremap => remap

  use Chem_GroupMod
  use Bundle_IncrementMod

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION: This gridded component (GC) combines the Superdynamics GC,
!   and Physics GC into a new composite Agcm GC.

!\begin{verbatim}
!       DUDT .... Mass-Weighted U-Wind      Tendency (Pa m /s)
!       DVDT .... Mass-Weighted V-Wind      Tendency (Pa m /s)
!       DPEDT ... Edge-Pressure             Tendency (Pa   /s)
!       DTDT .... Mass-Weighted Temperature Tendency (Pa K /s)
!       TRACER .. Friendly Tracers                   (unknown)
!     If Non-Hydrostatic Dynamics
!       DWDT .... Mass-Weighted W-Wind      Tendency (Pa m /s)
!\end{verbatim}

!EOP

  integer :: SDYN
  integer :: PHYS
  integer :: ORB

  type CONNECT_IAUcoeffs
     real, pointer              :: dfi(:) => NULL()
     integer                    :: istep
  end type CONNECT_IAUcoeffs

! Wrapper for extracting internal state
! -------------------------------------
  type IAU_coeffs
     type (CONNECT_IAUcoeffs), pointer :: PTR
  end type IAU_coeffs

  type CONNECT_ANAnBKG
     private
     class (AbstractRegridder), pointer :: ANA2BKG_regridder => null()
     class (AbstractRegridder), pointer :: BKG2ANA_regridder => null()
     type (ESMF_Grid)           :: GRIDana
     type (ESMF_Alarm)          :: AnaTendStartAlarm
     character(len=ESMF_MAXSTR) :: gridAnaName
     real, pointer              :: phis_bkg(:,:)
     integer                    :: IM
     integer                    :: JM
     integer                    :: LM
     logical                    :: do_transforms=.false.
     logical                    :: initialized=.false.
  end type CONNECT_ANAnBKG

! Wrapper for extracting internal state
! -------------------------------------
  type ANAnBKG
     type (CONNECT_ANAnBKG), pointer :: PTR
  end type ANAnBKG

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION:  The SetServices for the Physics GC needs to register its
!   Initialize and Run.  It uses the MAPL\_Generic construct for defining
!   state specs and couplings among its children.  In addition, it creates the
!   children GCs (SURF, CHEM, RADIATION, MOIST, TURBULENCE) and runs their
!   respective SetServices.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Locals

    integer                       :: I
    integer                       :: RST, SCM_SL
    logical                       :: ANA_TS
    type (ESMF_Config)            :: CF
    type (MAPL_MetaComp), pointer :: MAPL
    character(len=ESMF_MAXSTR)    :: ReplayMode

    type (CONNECT_ANAnBKG), pointer :: upd_internal_state
    type (ANAnBKG)                  :: wrap
    type (CONNECT_IAUcoeffs), pointer :: iau_coeffs_internal_state
    type (IAU_coeffs)                 :: wrap_iau_coeffs

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run       , RC=STATUS )
    VERIFY_(STATUS)

! Get the configuration from the component
!-----------------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, I, Label="ANALYZE_TS:", default=0, RC=STATUS)
    VERIFY_(STATUS)
    ANA_TS = I /= 0

    if (ANA_TS) then
        RST = MAPL_RestartRequired
    else
        RST = MAPL_RestartSkip
    end if

    call MAPL_GetResource(MAPL, SCM_SL, Label="SCM_SL:", default=0, _RC)


    call MAPL_GetResource(MAPL, ReplayMode, Label='REPLAY_MODE:', default="NoReplay", RC=STATUS )
    VERIFY_(STATUS)
    if(ANA_TS .and. ( adjustl(ReplayMode) /= "Exact"      .and. &
                      adjustl(ReplayMode) /= "Regular" ) ) then
             _ASSERT( adjustl(ReplayMode) == "NoReplay"  ,'needs informative message')
    endif

!BOS

! !IMPORT STATE:

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_analysis_increment',          &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_analysis_increment',         &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'temperature_analysis_increment',            &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'edge_pressure_analysis_increment',          &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DQVDT',                                     &
         LONG_NAME  = 'specific_humidity_analysis_increment',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DO3DT',                                     &
         LONG_NAME  = 'ozone_analysis_increment',                  &
         UNITS      = 'mol mol-1',                                 &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DTSDT',                                     &
         LONG_NAME  = 'skin_temperature_increment',                &
         UNITS      = 'K',                                         &
         RESTART    = RST,                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

! !INTERNAL STATE:

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_bias_tendency',               &
         UNITS      = 'm s-2',                                     &
         FRIENDLYTO = trim(COMP_NAME),                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_bias_tendency',              &
         UNITS      = 'm s-2',                                     &
         FRIENDLYTO = trim(COMP_NAME),                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'temperature_bias_tendency',                 &
         UNITS      = 'K s-1',                                     &
         FRIENDLYTO = trim(COMP_NAME),                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'edge_pressure_bias_tendency',               &
         UNITS      = 'Pa s-1',                                    &
         FRIENDLYTO = trim(COMP_NAME),                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DQVDT',                                     &
         LONG_NAME  = 'specific_humidity_bias_tendency',           &
         UNITS      = 'kg kg-1 s-1',                               &
         FRIENDLYTO = trim(COMP_NAME),                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DO3DT',                                     &
         LONG_NAME  = 'ozone_bias_tendency',                       &
         UNITS      = 'mol mol-1 s-1',                             &
         FRIENDLYTO = trim(COMP_NAME),                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DTSDT',                                     &
         LONG_NAME  = 'skin_temperature_tendency',                 &
         UNITS      = 'K s-1',                                     &
         FRIENDLYTO = trim(COMP_NAME),                             &
         RESTART    = RST,                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

! !EXPORT STATE:

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'QVBKG',                                     &
         LONG_NAME  = 'specific_humidity_background',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'QVANA',                                     &
         LONG_NAME  = 'specific_humidity_after_analysis',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'QVCON',                                     &
         LONG_NAME  = 'specific_humidity_after_constraint',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DPBKG',                                     &
         LONG_NAME  = 'delta_pressure_background',      &
         UNITS      = 'Pa',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DPANA',                                     &
         LONG_NAME  = 'delta_pressure_after_analysis',      &
         UNITS      = 'Pa',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DPSDT_CON',                                     &
         LONG_NAME  = 'surface_pressure_adjustment_due_to_constraint', &
         UNITS      = 'Pa s-1',                                        &
         DIMS       = MAPL_DimsHorzOnly,                               &
         VLOCATION  = MAPL_VLocationNone,                   RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DUDT_ANA',                                  &
         LONG_NAME  = 'total_eastward_wind_analysis_tendency',     &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DVDT_ANA',                                  &
         LONG_NAME  = 'total_northward_wind_analysis_tendency',    &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DTDT_ANA',                                  &
         LONG_NAME  = 'total_temperature_analysis_tendency',       &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DPEDT_ANA',                                 &
         LONG_NAME  = 'total_edge_pressure_analysis_tendency',    &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'DQVDT_ANA',                                        &
         LONG_NAME  = 'total_specific_humidity_vapor_analysis_tendency',  &
         UNITS      = 'kg kg-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzVert,                                  &
         VLOCATION  = MAPL_VLocationCenter,                    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'DQVDT_CON',                                        &
         LONG_NAME  = 'total_specific_humidity_vapor_analysis_tendency_due_to_Constraint',  &
         UNITS      = 'kg kg-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzVert,                                  &
         VLOCATION  = MAPL_VLocationCenter,                    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'DQLDT_ANA',                                        &
         LONG_NAME  = 'total_specific_humidity_liquid_analysis_tendency', &
         UNITS      = 'kg kg-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzVert,                                  &
         VLOCATION  = MAPL_VLocationCenter,                    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'DQIDT_ANA',                                        &
         LONG_NAME  = 'total_specific_humidity_ice_analysis_tendency',    &
         UNITS      = 'kg kg-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzVert,                                  &
         VLOCATION  = MAPL_VLocationCenter,                    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'DQRDT_ANA',                                        &
         LONG_NAME  = 'total_suspended_rain_analysis_tendency',    &
         UNITS      = 'kg kg-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzVert,                                  &
         VLOCATION  = MAPL_VLocationCenter,                    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'DQSDT_ANA',                                        &
         LONG_NAME  = 'total_suspended_snow_analysis_tendency',    &
         UNITS      = 'kg kg-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzVert,                                  &
         VLOCATION  = MAPL_VLocationCenter,                    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'DQGDT_ANA',                                        &
         LONG_NAME  = 'total_suspended_graupe_analysis_tendency',    &
         UNITS      = 'kg kg-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzVert,                                  &
         VLOCATION  = MAPL_VLocationCenter,                    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DO3DT_ANA',                                 &
         LONG_NAME  = 'total_ozone_analysis_tendency',             &
         UNITS      = 'mol mol-1 s-1',                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DTSDT_ANA',                                 &
         LONG_NAME  = 'total_skin_temperature_tendency',           &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'DTHVDTFILINT',                                     &
         LONG_NAME  = 'vertically_integrated_thv_adjustment_from_filling',&
         UNITS      = 'K kg m-2 s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'PERES',                                            &
         LONG_NAME  = 'vertically_integrated_cpt_tendency_residual',      &
         UNITS      = 'W m-2',                                            &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'PEFILL',                                           &
         LONG_NAME  = 'vertically_integrated_cpt_adjustment_from_filling',&
         UNITS      = 'W m-2',                                            &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'QTFILL',                                           &
         LONG_NAME  = 'vertically_integrated_total_water_adjustment_from_filling', &
         UNITS      = 'kg m-2 s-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'QVFILL',                                           &
         LONG_NAME  = 'vertically_integrated_qv_adjustment_from_filling', &
         UNITS      = 'kg m-2 s-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'QLFILL',                                           &
         LONG_NAME  = 'vertically_integrated_ql_adjustment_from_filling', &
         UNITS      = 'kg m-2 s-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'QIFILL',                                           &
         LONG_NAME  = 'vertically_integrated_qi_adjustment_from_filling', &
         UNITS      = 'kg m-2 s-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME = 'OXFILL',                                           &
         LONG_NAME  = 'vertically_integrated_ox_adjustment_from_filling', &
         UNITS      = 'kg m-2 s-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                                  &
         VLOCATION  = MAPL_VLocationNone,                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPP_EPV',                                                 &
       LONG_NAME          = 'tropopause_pressure_based_on_EPV_estimate',                 &
       UNITS              = 'Pa',                                                        &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPP_THERMAL',                                             &
       LONG_NAME          = 'tropopause_pressure_based_on_thermal_estimate',             &
       UNITS              = 'Pa',                                                        &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPP_BLENDED',                                             &
       LONG_NAME          = 'tropopause_pressure_based_on_blended_estimate',             &
       UNITS              = 'Pa',                                                        &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPT',                                                     &
       LONG_NAME          = 'tropopause_temperature_using_blended_TROPP_estimate',       &
       UNITS              = 'K',                                                         &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPQ',                                                     &
       LONG_NAME          = 'tropopause_specific_humidity_using_blended_TROPP_estimate', &
       UNITS              = 'kg kg-1',                                                   &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'TQV',                                        &
         LONG_NAME        = 'total_precipitable_water_vapor',             &
         UNITS            = 'kg m-2'  ,                                   &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'TQI',                                        &
         LONG_NAME        = 'total_precipitable_ice_water',               &
         UNITS            = 'kg m-2'  ,                                   &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'TQL',                                        &
         LONG_NAME        = 'total_precipitable_liquid_water',            &
         UNITS            = 'kg m-2'  ,                                   &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'TOX',                                        &
         LONG_NAME        = 'total_column_odd_oxygen',                    &
         UNITS            = 'kg m-2'  ,                                   &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'MASS',                                       &
         LONG_NAME        = 'atmospheric_mass',                           &
         UNITS            = 'kg m-2'  ,                                   &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'KE',                                         &
         LONG_NAME        = 'vertically_integrated_kinetic_energy',       &
         UNITS            = 'J m-2'  ,                                    &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'CPT',                                        &
         LONG_NAME        = 'vertically_integrated_enthalpy',             &
         UNITS            = 'J m-2'  ,                                    &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'THV',                                        &
         LONG_NAME        = 'vertically_integrated_virtual_potential_temperature',             &
         UNITS            = 'K'  ,                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'QLTOT',                                      &
         LONG_NAME        = 'mass_fraction_of_cloud_liquid_water',        &
         UNITS            = 'kg kg-1',                                    &
         DIMS             = MAPL_DimsHorzVert,                            &
         VLOCATION        = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'QITOT',                                      &
         LONG_NAME        = 'mass_fraction_of_cloud_ice_water',           &
         UNITS            = 'kg kg-1',                                    &
         DIMS             = MAPL_DimsHorzVert,                            &
         VLOCATION        = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'PHIS',                                       &
         LONG_NAME        = 'surface geopotential height',                &
         UNITS            = 'm+2 s-2',                                    &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'SGH',                                        &
         LONG_NAME        = 'isotropic stdv of GWD topography',           &
         UNITS            = 'm',                                          &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'GWDVARX',                                    &
         LONG_NAME        = 'east-west variance of GWD topography',       &
         UNITS            = 'm+2',                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'GWDVARY',                                    &
         LONG_NAME        = 'north-south variance of GWD topography',     &
         UNITS            = 'm+2',                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'GWDVARXY',                                   &
         LONG_NAME        = 'SW-NE variance of GWD topography',           &
         UNITS            = 'm+2',                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'GWDVARYX',                                   &
         LONG_NAME        = 'NW-SE variance of GWD topography',           &
         UNITS            = 'm+2',                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'TRBVAR',                                     &
         LONG_NAME        = 'isotropic variance of TRB topography',       &
         UNITS            = 'm+2',                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'VARFLT',                                     &
         LONG_NAME        = 'isotropic variance of filtered topography',  &
         UNITS            = 'm+2',                                        &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                         &
         SHORT_NAME       = 'TRADVI',                                     &
         LONG_NAME        = 'advected_quantities_tendencies',             &
         units            = 'dX s-1',                                     &
         DIMS             = MAPL_DimsHorzVert,                            &
         DATATYPE         = MAPL_BundleItem,                              &
                                                               RC=STATUS  )
    VERIFY_(STATUS)


! Create childrens gridded components and invoke their SetServices
! ----------------------------------------------------------------
#ifdef SCM
    SDYN = MAPL_AddChild(GC, NAME='SCMDYNAMICS', SS=SDYN_SetServices, RC=STATUS)
    VERIFY_(STATUS)
#else
    SDYN = MAPL_AddChild(GC, NAME='SUPERDYNAMICS', SS=SDYN_SetServices, RC=STATUS)
    VERIFY_(STATUS)
#endif
    PHYS = MAPL_AddChild(GC, NAME='PHYSICS', SS=PHYS_SetServices, RC=STATUS)
    VERIFY_(STATUS)

    ORB  = MAPL_AddChild(GC, NAME='ORBIT', SS=ORB_SetServices, RC=STATUS)
    VERIFY_(STATUS)

! Export for IAU or Analysis purposes
! -----------------------------------
!   call MAPL_AddExportSpec ( GC, &
!        SHORT_NAME = 'PHIS', &
!        CHILD_ID = SDYN, &
!        RC=STATUS )
!   VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'AREA', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'AK', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'BK', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'PLE', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

#ifdef HAS_GIGATRAJ
    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'OMEGA', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'PL', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'TH', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)
#endif

    call MAPL_AddExportSpec( GC, &
         SHORT_NAME = 'PS', &
         CHILD_ID = SDYN, &
         RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'DELP', &
         CHILD_ID = SDYN, &
         RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'PE', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'PT', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'TV', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'T', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'U', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'V', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'W', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'U_DGRID', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'V_DGRID', &
         CHILD_ID = SDYN, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'PPBL', &
         CHILD_ID   = PHYS,  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'O3PPMV', &
         CHILD_ID   = PHYS,  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'OX',  &
         CHILD_ID   = PHYS,  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'Q', &
         CHILD_ID = PHYS, &
         RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'QCTOT', &
         CHILD_ID = PHYS, &
         RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'U10N', &
         CHILD_ID = PHYS, &
         RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'V10N', &
         CHILD_ID = PHYS, &
         RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'SNOMAS', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'WET1',  &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'TSOIL1', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'LWI', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'Z0', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'TS', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, &
         SHORT_NAME = 'TRANA', &
         CHILD_ID = PHYS, &
         RC = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'FRLAND', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'FRLANDICE', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'FRLAKE', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'FROCEAN', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
         SHORT_NAME = 'FRACI', &
         CHILD_ID = PHYS, &
         RC=STATUS  )
    VERIFY_(STATUS)
!EOS

! Set internal connections between the childrens IMPORTS and EXPORTS
! ------------------------------------------------------------------

    call MAPL_AddConnectivity ( GC,                                                        &
         SRC_NAME  = (/'U            ','V            ','TH           ','T            ',    &
                       'ZLE          ','PS           ','TA           ','QA           ',    &
                       'US           ','VS           ',                                    &
                       'SPEED        ','DZ           ','PLE          ','W            ',    &
                       'PREF         ','TROPP_BLENDED','S            ','PLK          ',    &
                       'PV           ','TROPK_BLENDED','OMEGA        ','PKE          '/),  &
         DST_NAME  = (/'U     ','V     ','TH    ','T     ',                                &
                       'ZLE   ','PS    ','TA    ','QA    ',                                &
                       'UA    ','VA    ',                                                  &
                       'SPEED ','DZ    ','PLE   ','W     ',                                &
                       'PREF  ','TROPP ','S     ','PLK   ',                                &
                       'PV    ','TROPK ','OMEGA ','PKE   '/),                              &
         DST_ID = PHYS,                                                                    &
         SRC_ID = SDYN,                                                                    &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddConnectivity ( GC,                                                        &
         SRC_NAME  = (/'T            ','PLE          ','ZLE          ','TROPP_BLENDED'/),  &
         DST_NAME  = (/'T_avg24      ','PLE_avg24    ','ZLE_avg24    ','TROPP_avg24  '/),  &
         DST_ID = PHYS,                                                                    &
         SRC_ID = SDYN,                                                                    &
         RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddConnectivity ( GC,              &
         SRC_NAME    = 'PLE',                    &
         DST_NAME    = 'PLEINST',                &
         SRC_ID      = SDYN,                     &
         DST_ID      = PHYS,                     &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddConnectivity ( GC, SRC_NAME = 'AREA', DST_NAME = 'AREA', &
         SRC_ID = SDYN, DST_ID = PHYS, RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddConnectivity ( GC,                                         &
          SRC_NAME  = (/'DTDTDYN      ','DQVDTDYN     ','PLE_DYN_IN   ',    &
                        'U_DYN_IN     ','V_DYN_IN     ','QV_DYN_IN    ',    &
                        'T_DYN_IN     '                                     &
                      /),                                                   &
          DST_NAME  = (/'DTDTDYN      ','DQVDTDYN     ','PLE_DYN_IN   ',    &
                        'U_DYN_IN     ','V_DYN_IN     ','QV_DYN_IN    ',    &
                        'T_DYN_IN     '                                     &
                      /),                                                   &
          SRC_ID = SDYN,                                                    &
          DST_ID = PHYS,                                                    &
          RC=STATUS  )
    VERIFY_(STATUS)

! Bundle of quantities to be advected
!------------------------------------

     call MAPL_AddConnectivity ( GC,                               &
         SRC_NAME  = 'TRADV',                                      &
         DST_NAME  = 'TRADV',                                      &
         SRC_ID      = PHYS,                                       &
         DST_ID      = SDYN,                                       &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

! Orbital Component Bundle
     call MAPL_AddConnectivity( GC,                                &
         SRC_NAME='SATORB',                                        &
         DST_NAME='SATORB',                                        &
         SRC_ID = ORB,                                             &
         DST_ID = PHYS,                                            &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

    if (SCM_SL /= 0 ) then
    print *,'AgcmGC: adding connectivity for LHOBS and SHOBS'
    call MAPL_AddConnectivity ( GC,    &
!         SHORT_NAME  = (/'TSKINOBS','QSKINOBS','LHOBS   ','SHOBS   '/), &
         SHORT_NAME  = (/'LHOBS   ','SHOBS   '/), &
         DST_ID = PHYS,         &
         SRC_ID = SDYN,         &
         RC=STATUS  )
     VERIFY_(STATUS)
     end if

!ALT: we need this if we run with NCEP gwd
    call MAPL_AddConnectivity ( GC,    &
         SHORT_NAME  = (/'DXC'/),      &
         DST_ID = PHYS,                &
         SRC_ID = SDYN,                &
         RC=STATUS  )
     VERIFY_(STATUS)

! We Terminate these IMPORTS which are manually filled
!-----------------------------------------------------

     call MAPL_TerminateImport    ( GC,                                                     &
          SHORT_NAME = (/'DUDT  ','DVDT  ','DWDT  ','DTDT  ','DPEDT ','DQVANA','DQLANA',    &
                         'DQIANA','DQRANA','DQSANA','DQGANA','DOXANA','PHIS  ','VARFLT'/),  &
          CHILD      = SDYN,                                                                &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_TerminateImport    ( GC,                          &
          SHORT_NAME = (/'VARFLT','PHIS  ','SGH   ', 'DTSDT '/), &
          CHILD      = PHYS,                                     &
          RC=STATUS  )
     VERIFY_(STATUS)

! Allocate this instance of the internal state and put it in wrapper
! ------------------------------------------------------------------
    allocate( upd_internal_state, stat=status )
    VERIFY_(STATUS)
    wrap%ptr => upd_internal_state
    allocate( iau_coeffs_internal_state, stat=status )
    VERIFY_(STATUS)
    wrap_iau_coeffs%ptr => iau_coeffs_internal_state

! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------
    call ESMF_UserCompSetInternalState ( GC, 'UPD_STATE', wrap, status )
    VERIFY_(STATUS)
    call ESMF_UserCompSetInternalState ( GC, 'IAU_COEFFS', wrap_iau_coeffs, status )
    VERIFY_(STATUS)


    call MAPL_GenericSetServices    ( GC, RC=STATUS )
    VERIFY_(STATUS)

! Clocks
!-------

    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN"           ,RC=STATUS)
    VERIFY_(STATUS)

! All done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Initialize -- Initialize method for the composite Agcm Gridded Component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION:


!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)           :: IAm
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases

   type (MAPL_MetaComp),  pointer  :: STATE
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_State),         pointer  :: GEX(:)
   type (ESMF_Field)                   :: FIELD
   type (ESMF_Time)                    :: CurrTime, RingTime
   type (ESMF_TimeInterval)            :: TIMEINT
   type (ESMF_Alarm)                   :: ALARM
   type (ESMF_Alarm)                   :: ALARM4D
   type (ESMF_Config)                  :: cf
   integer                             :: I, NQ
   real                                :: POFFSET, DT
   real, pointer, dimension(:,:)       :: PHIS,SGH,VARFLT,PTR
   real, pointer, dimension(:,:,:)     :: TEND!
   character(len=ESMF_MAXSTR)          :: replayMode
   real                                :: RPL_INTERVAL
   real                                :: RPL_SHUTOFF
   real                                :: IAU4dFREQ
   integer                             :: PREDICTOR_DURATION
   integer                             :: MKIAU_FREQUENCY
   character(len=ESMF_MAXSTR), parameter :: INITIALIZED_EXPORTS(3) = &
        (/'PHIS  ', 'SGH   ', 'VARFLT' /)

   logical                             :: DasMode
   character(len=ESMF_MAXSTR)          :: STRING
   character(len=ESMF_MAXSTR)          :: rplMode

! =============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, config=cf, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Initialize"

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)


    call MAPL_TimerOn(STATE,"INITIALIZE")

! Call Initialize for every Child

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(STATE,"TOTAL")

! Get children and their im/ex states from my generic state.
!----------------------------------------------------------

    call MAPL_Get ( STATE, GIM=GIM, GEX=GEX, RC=STATUS )
    VERIFY_(STATUS)

! Initialize the advection increments bundle (TRADVI)
! with tracer increment names
!-----------------------------------------------------

    call Initialize_IncBundle_init(GC, GEX(PHYS), EXPORT, DYNinc, __RC__)

! Make sure that the physics tendencies are allocated
!----------------------------------------------------

    call MAPL_GetPointer(GEX(PHYS), TEND, 'DUDT' , ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GEX(PHYS), TEND, 'DVDT' , ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GEX(PHYS), TEND, 'DWDT' , ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GEX(PHYS), TEND, 'DTDT' , ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GEX(PHYS), TEND, 'DPEDT', ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)

! Fill Childrens TOPO variables and Diagnostics
!----------------------------------------------
    call MAPL_GetPointer(EXPORT, PHIS,   'PHIS',   ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SGH,    'SGH',    ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VARFLT, 'VARFLT', ALLOC=.true., rc=STATUS)
    VERIFY_(STATUS)

! PHIS (topography)...
!---------
    call ESMF_StateGet( GIM(SDYN), 'PHIS', FIELD, rc=STATUS )
    VERIFY_(STATUS)
    Call GEOS_TopoGet ( cf, MEAN=FIELD, rc=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldGet (FIELD, localDE=0, farrayPtr=PTR, rc = status)
    VERIFY_(STATUS)
    PHIS = PTR

! Pass PHIS into PHYS
!---------
    call ESMF_StateGet( GIM(PHYS), 'PHIS', FIELD, rc=STATUS )
    VERIFY_(STATUS)
    Call GEOS_TopoGet ( cf, MEAN=FIELD, rc=STATUS )
    VERIFY_(STATUS)

! GWDVAR (standard deviation)...
!-----------
    call ESMF_StateGet( GIM(PHYS), 'SGH', FIELD, rc=STATUS )
    VERIFY_(STATUS)
    Call GEOS_TopoGet ( cf, GWDVAR=FIELD, rc=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldGet (FIELD, localDE=0, farrayPtr=PTR, rc = status)
    SGH = PTR

! TRBVAR (variance)...
!-----------
    call ESMF_StateGet( GIM(PHYS), 'VARFLT', FIELD, rc=STATUS )
    VERIFY_(STATUS)
    Call GEOS_TopoGet ( cf, TRBVAR=FIELD, rc=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldGet (FIELD, localDE=0, farrayPtr=PTR, rc = status)
    VERIFY_(STATUS)
    VARFLT = PTR

! Pass variance into SDYN
!-----------
    call ESMF_StateGet( GIM(SDYN), 'VARFLT', FIELD, rc=STATUS )
    VERIFY_(STATUS)
    Call GEOS_TopoGet ( cf, TRBVAR=FIELD, rc=STATUS )
    VERIFY_(STATUS)

! ======================================================================
!ALT: the next section addresses the problem when export variables have been
!     assigned values during Initialize. To prevent "connected" exports
!     being overwritten by DEFAULT in the Import spec in the other component
!     we label them as being "initailized by restart". A better solution
!     would be to move the computation to phase 2 of Initialize and
!     eliminate this section alltogether
! ======================================================================
    DO I = 1, size(INITIALIZED_EXPORTS)
       call ESMF_StateGet(EXPORT,INITIALIZED_EXPORTS(I), FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", &
                              VALUE=MAPL_InitialRestart, RC=STATUS)
       VERIFY_(STATUS)
    END DO

! Initialize Predictor Alarm
!---------------------------

   call ESMF_ClockGet(clock, currTime=currTime, rc=status)
   VERIFY_(STATUS)

   call MAPL_GetResource( STATE, DT, Label="RUN_DT:", RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetResource( STATE, POFFSET, Label="PREDICTOR_OFFSET:", default=21600. , RC=STATUS)
   VERIFY_(STATUS)

   call ESMF_TimeIntervalSet(TIMEINT,  S=nint (POFFSET), RC=STATUS)
   VERIFY_(STATUS)

   ringTime = currTime+TIMEINT

   call ESMF_TimeIntervalSet(TIMEINT,  S=nint(DT) , RC=STATUS)
   VERIFY_(STATUS)

   ALARM = ESMF_AlarmCreate( name='PredictorAlarm',    &
                             CLOCK = CLOCK,            &
                             RingTime     = ringTime,  &
                             RingInterval = TIMEINT,   &
                             RC           = STATUS     )
   VERIFY_(STATUS)
   if(ringTime == currTime) then
      call ESMF_AlarmRingerOn(Alarm, rc=status)
      VERIFY_(STATUS)
   end if

 ! if(MAPL_AM_I_ROOT() ) then
 !    PRINT *
 !    PRINT *,TRIM(Iam)//": PredictorAlarm settings"
 !    call ESMF_TimeGet( currTIME, timeString=String, RC=STATUS)
 !    VERIFY_(STATUS)
 !    PRINT *,TRIM(Iam)//": CurrTime: ",trim(string)
 !    call ESMF_TimeGet( RingTIME, timeString=String, RC=STATUS)
 !    VERIFY_(STATUS)
 !    PRINT *,TRIM(Iam)//": RingTime: ",trim(string)
 !    PRINT *,TRIM(Iam)//": Is Ringing: ", ESMF_AlarmIsRinging(ALARM)
 !    PRINT *
 ! endif


   ! Detect if running DasMode (Checking for AGCM_IMPORT)
   ! ----------------------------------------------------
   call MAPL_GetResource( STATE, STRING, LABEL="IMPORT_RESTART_FILE:", RC=STATUS)
   IF (STATUS == ESMF_SUCCESS) THEN
      DasMode = .true.
   ELSE
      DasMode = .false.
   END IF

   ! Detect if running REPLAY
   ! ------------------------
   call MAPL_GetResource( STATE, ReplayMode, 'REPLAY_MODE:', default="NoReplay", RC=STATUS )
   VERIFY_(STATUS)

       rplMode = adjustl(ReplayMode)
   if( rplMode=="Regular" .or. rplMode == "Exact" ) then
       DasMode = .true.
   end if

   ! Disable the predictor alarm if not dasmode
   ! ------------------------------------------
   if (.not. DasMode) then
      call ESMF_AlarmDisable(ALARM, rc=status)
      VERIFY_(STATUS)
   end if

   call MAPL_StateAlarmAdd(STATE,ALARM,RC=status)
   VERIFY_(STATUS)

  ! Note: PREDICTOR_DURATION and MKIAU_FREQUENCY are Initialized in GCM_GridComp
  ! ----------------------------------------------------------------------------
    call MAPL_GetResource( STATE, PREDICTOR_DURATION, Label="PREDICTOR_DURATION:", RC=STATUS ) ; VERIFY_(STATUS)
    call MAPL_GetResource( STATE, MKIAU_FREQUENCY,    Label="MKIAU_FREQUENCY:",    RC=STATUS ) ; VERIFY_(STATUS)

   if(   (adjustl(ReplayMode)=="Exact"  ) .or.   &
       ( (adjustl(ReplayMode)=="Regular") .and. (PREDICTOR_DURATION.gt.MKIAU_FREQUENCY/2) )  ) then

      call MAPL_GetResource(STATE, RPL_INTERVAL, 'REPLAY_INTERVAL:', default=21600., RC=STATUS )
      VERIFY_(STATUS)
      call ESMF_TimeIntervalSet(TIMEINT, S=nint(RPL_INTERVAL), RC=STATUS)
      VERIFY_(STATUS)

      ALARM = ESMF_AlarmCreate ( name='ExactReplay', clock=CLOCK, RingTime=currTime, ringInterval=TIMEINT, sticky=.false., RC=STATUS )
      VERIFY_(STATUS)
      call ESMF_AlarmRingerOn(ALARM, rc=status)
      VERIFY_(STATUS)
      _ASSERT(POFFSET == RPL_INTERVAL,'needs informative message')
   end if

!  Create 4dIAU alarm
!  ------------------
   call ESMF_ClockGet(clock, currTime=currTime, rc=status)
   VERIFY_(STATUS)

   call MAPL_GetResource( STATE, IAU4dFREQ, Label="4DIAU_FREQUENCY:", default=3600. , RC=STATUS)
   VERIFY_(STATUS)

   call ESMF_TimeIntervalSet(TIMEINT,  S=nint (IAU4dFREQ), RC=STATUS)
   VERIFY_(STATUS)

   ALARM4D = ESMF_AlarmCreate( name='4DIAUalarm', &
                               CLOCK = CLOCK, &
                               RingInterval = TIMEINT  ,  &
                               RingTime     = currTime,   &
                               STICKY       = .FALSE.,    &
                               RC           = STATUS      )
   VERIFY_(STATUS)
   call ESMF_AlarmRingerOn(Alarm4D, rc=status)
   VERIFY_(STATUS)

    call MAPL_TimerOff(STATE,"TOTAL")
    call MAPL_TimerOff(STATE,"INITIALIZE")

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": IMPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( IMPORT, rc=STATUS )
    call WRITE_PARALLEL ( trim(Iam)//": EXPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( EXPORT, rc=STATUS )
#endif


    RETURN_(ESMF_SUCCESS)
 end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Run -- Run method for the composite Agcm Gridded Component

! !INTERFACE:

  subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION:


!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)           :: IAm
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases

   type (ESMF_VM)                      :: VMG
   type (MAPL_MetaComp),      pointer  :: STATE
   type (ESMF_GridComp),      pointer  :: GCS(:)
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_State),         pointer  :: GEX(:)
   type (ESMF_State)                   :: INTERNAL
   type (ESMF_Alarm)                   :: ALARM
   type (ESMF_Alarm)                   :: ALARM4D
   type (ESMF_TimeInterval)            :: TINT
   type (ESMF_FieldBundle)             :: Bundle
   type (ESMF_FieldBundle)             :: Advect_Bundle
   type (ESMF_Grid)                    :: grid

   real, pointer, dimension(:)         :: PREF   => null()
   real, pointer, dimension(:,:,:)     :: U      => null()
   real, pointer, dimension(:,:,:)     :: V      => null()
   real, pointer, dimension(:,:,:)     :: W      => null()
   real, pointer, dimension(:,:,:)     :: T      => null()
   real, pointer, dimension(:,:,:)     :: Q      => null()
   real, pointer, dimension(:,:,:)     :: QLLS   => null()
   real, pointer, dimension(:,:,:)     :: QLCN   => null()
   real, pointer, dimension(:,:,:)     :: QILS   => null()
   real, pointer, dimension(:,:,:)     :: QICN   => null()
   real, pointer, dimension(:,:,:)     :: QRAIN   => null()
   real, pointer, dimension(:,:,:)     :: QSNOW   => null()
   real, pointer, dimension(:,:,:)     :: QGRAUPEL   => null()
   real, pointer, dimension(:,:,:)     :: PLE    => null()
   real, pointer, dimension(:,:,:)     :: EPV    => null()
   real, pointer, dimension(:,:,:)     :: QLTOT  => null()
   real, pointer, dimension(:,:,:)     :: QITOT  => null()
   real, pointer, dimension(:,:,:)     :: DPEDT  => null()
   real, pointer, dimension(:,:,:)     :: DTDT   => null()
   real, pointer, dimension(:,:,:)     :: TENDAN => null()

!! real,   allocatable, dimension(:,:)   :: ALPHA2D
   real,   allocatable, dimension(:,:)   :: QFILL
   real,   allocatable, dimension(:,:)   :: QINT
   real,   allocatable, dimension(:,:)   :: ALF_BKS_INT
   real,   allocatable, dimension(:,:)   :: DRY_BKG_INT
   real,   allocatable, dimension(:,:)   :: DRY_ANA_INT
   real,   allocatable, dimension(:,:)   :: QDP_BKG_INT
   real,   allocatable, dimension(:,:)   :: QDP_ANA_INT
   real*8, allocatable, dimension(:,:)   :: SUMKE
   real*8, allocatable, dimension(:,:)   :: SUMCPT1, SUMCPT2
   real*8, allocatable, dimension(:,:)   :: SUMTHV1, SUMTHV2
   real*8, allocatable, dimension(:,:,:) :: PKE
   real*8, allocatable, dimension(:,:,:) :: PKZ

   real, pointer, dimension(:,:)       :: AREA   => null()
   real, pointer, dimension(:,:)       :: OXFILL => null()
   real, pointer, dimension(:,:)       :: QTFILL => null()
   real, pointer, dimension(:,:)       :: QVFILL => null()
   real, pointer, dimension(:,:)       :: QIFILL => null()
   real, pointer, dimension(:,:)       :: QLFILL => null()
   real, pointer, dimension(:,:)       :: TQV    => null()
   real, pointer, dimension(:,:)       :: TQI    => null()
   real, pointer, dimension(:,:)       :: TQL    => null()
   real, pointer, dimension(:,:)       :: TOX    => null()
   real, pointer, dimension(:,:)       :: TROPP1 => null()
   real, pointer, dimension(:,:)       :: TROPP2 => null()
   real, pointer, dimension(:,:)       :: TROPP3 => null()
   real, pointer, dimension(:,:)       :: TROPT  => null()
   real, pointer, dimension(:,:)       :: TROPQ  => null()
   real, pointer, dimension(:,:)       :: MASS   => null()
   real, pointer, dimension(:,:)       :: KE     => null()
   real, pointer, dimension(:,:)       :: CPT    => null()
   real, pointer, dimension(:,:)       :: THV    => null()

   real, pointer, dimension(:,:)       :: PERES        => null()
   real, pointer, dimension(:,:)       :: PEFILL       => null()
   real, pointer, dimension(:,:)       :: PEPHY_SDYN   => null()  ! D(CpT)DT from ADD_INCS (SuperDYNamics)
   real, pointer, dimension(:,:)       :: PEPHY_PHYS   => null()  ! D(CpT)DT from PHYSics

   real, pointer, dimension(:,:)       :: DTHVDTFILINT => null()
   real, pointer, dimension(:,:)       :: DTHVDTPHYINT => null()
   real, pointer, dimension(:,:)       :: DQVDTPHYINT  => null()
   real, pointer, dimension(:,:)       :: DQLDTPHYINT  => null()
   real, pointer, dimension(:,:)       :: DQIDTPHYINT  => null()
   real, pointer, dimension(:,:)       :: DOXDTPHYINT  => null()

   real, pointer, dimension(:,:,:)     :: DP
   real, pointer, dimension(:,:,:)     :: PL
   real, pointer, dimension(:,:,:)     :: FC
   real, pointer, dimension(:,:,:)     :: TROP
   real, pointer, dimension(:,:,:)     :: XXINC
   real, pointer, dimension(:,:,:)     :: DQVANA
   real, pointer, dimension(:,:,:)     :: DQLANA
   real, pointer, dimension(:,:,:)     :: DQIANA
   real, pointer, dimension(:,:,:)     :: DQRANA
   real, pointer, dimension(:,:,:)     :: DQSANA
   real, pointer, dimension(:,:,:)     :: DQGANA
   real, pointer, dimension(:,:,:)     :: DOXANA
   real, pointer, dimension(:,:,:)     :: DQIMPORT => null()

   real, pointer, dimension(:,:)       :: ptr2d
   real, pointer, dimension(:,:,:)     :: ptr3d

   real, pointer, dimension(:)         :: AK
   real, pointer, dimension(:)         :: BK

   real*8, allocatable, dimension(:,:)   :: sumq
   real*8, allocatable, dimension(:,:)   :: sum_qdp_bkg
   real*8, allocatable, dimension(:,:)   :: sum_qdp_ana
   real*8, allocatable, dimension(:,:)   :: sum_dry_ana
   real*8, allocatable, dimension(:,:)   :: sum_dry_bkg
   real,   allocatable, dimension(:,:,:) ::     qdp_bkg
   real,   allocatable, dimension(:,:,:) ::     qdp_ana
   real,   allocatable, dimension(:,:,:) ::     dry_ana
   real,   allocatable, dimension(:,:,:) ::     dry_bkg
   real,   allocatable, dimension(:,:,:) ::  ple_ana
   real,   allocatable, dimension(:,:,:) ::   dp_ana
   real,   allocatable, dimension(:,:,:) :: tdpold
   real,   allocatable, dimension(:,:,:) :: tdpnew
   real,   allocatable, dimension(:,:,:) :: DQVCON

   real*8                                ::   gamma
   real*8                                ::  alf_bks_ave
   real*8                                ::  dry_bkg_ave
   real*8                                ::  dry_ana_ave
   real*8                                ::  qdp_bkg_ave
   real*8                                ::  qdp_ana_ave
   real*8                                :: qint_ana_ave
   real*8                                :: qint_bkg_ave

   real                                :: DT
   integer                             :: IM, JM, LM, L, TYPE, ISFCST
   integer                             :: NumFriendly
   integer                             :: K
   integer                             :: I
   integer                             :: PREDICTOR_DURATION
   integer                             :: MKIAU_FREQUENCY
   logical                             :: DasMode
   logical                             :: DO_PREDICTOR
   logical                             :: Begin_REPLAY_Cycle
   logical                             :: LAST_CORRECTOR
   integer                             :: CONSTRAIN_DAS
   real                                :: ALPHA, BETA, TAUANL, DTX, IAUcoeff
   real                                :: ALPHAQ, BETAQ
   real                                :: ALPHAO, BETAO
   real                                :: ALF, BET
   real                                :: EPS

   character(len=ESMF_MAXSTR)          :: ANA_IS_WEIGHTED
   logical                             ::     IS_WEIGHTED

   character(len=ESMF_MAXSTR), pointer :: Names(:)
   character(len=ESMF_MAXSTR)          :: STRING

   integer, parameter                  :: FREERUN   = 0
   integer, parameter                  :: PREDICTOR = 1
   integer, parameter                  :: CORRECTOR = 2
   integer, parameter                  :: FORECAST  = 3

   integer                             :: unit
   logical                             :: is_ringing
   logical                             ::   is_ExactReplay09_ringing
   logical                             :: is_RegularReplay09_ringing
   logical                             :: is_shutoff
   character(len=ESMF_MAXSTR)          :: FILENAME
   character(len=ESMF_MAXSTR)          :: FILETYPE
   character(len=ESMF_MAXSTR)          :: FileTmpl
   character(len=ESMF_MAXSTR)          :: FileTmpl09
   character(len=ESMF_MAXSTR)          :: replayFile
   character(len=ESMF_MAXSTR)          :: replayFile09
   character(len=ESMF_MAXSTR)          :: replayMode
   character(len=ESMF_MAXSTR)          :: rplMode
   type(ESMF_Time)                     :: currTime

   CHARACTER(LEN=ESMF_MAXSTR)          :: fieldName
   logical :: DO_4DIAU

   type (ESMF_Time)                    :: REPLAY_TIME
   type (ESMF_Time), save              :: REPLAY_TIME0
   logical,save                        :: first=.true.

! For DYN:PHY ratio estimates
   integer       :: START_TIME, END_TIME, DYN_TIME, PHY_TIME
   integer       :: COUNT_MAX, COUNT_RATE
   real(kind=8)  :: CRI

!ALT: for memory leak testing
   logical :: isPresent
   real, allocatable, target :: zero(:,:,:)
!   real, pointer :: zero(:,:,:) => null()

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet ( GC, VM=VMG, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (STATE,"TOTAL")
    call MAPL_TimerOn (STATE,"RUN"  )

! Get children and their im/ex states from my generic state.
!----------------------------------------------------------

    call MAPL_Get ( STATE, GCS=GCS, GIM=GIM, GEX=GEX,  &
                    INTERNAL_ESMF_STATE=INTERNAL,      &
                    IM=IM, JM=JM, LM=LM,               &
                    RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GC, grid=grid, rc=status)
    VERIFY_(STATUS)

! Get the 4DIAU alarm
!--------------------
    call ESMF_ClockGetAlarm(clock, alarmname='4DIAUalarm', alarm=Alarm4D, rc=status)
    VERIFY_(STATUS)

! Set the various time scales
!----------------------------

    call MAPL_GetResource( STATE,     DT,          Label="RUN_DT:",                        RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE,  ALPHA,          Label="ALPHA:",         default=0.0,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE,   BETA,          Label="BETA:",          default=1.0,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE, ALPHAQ,          Label="ALPHAQ:",        default=ALPHA,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE,  BETAQ,          Label="BETAQ:",         default=BETA,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE, ALPHAO,          Label="ALPHAO:",        default=ALPHA,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE,  BETAO,          Label="BETAO:",         default=BETA,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE, TAUANL,          Label="TAUANL:",        default=21600., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE, ISFCST,          Label="IS_FCST:",       default=0,      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( STATE, CONSTRAIN_DAS,   Label="CONSTRAIN_DAS:", default=1,      RC=STATUS); VERIFY_(STATUS)

  ! Note: PREDICTOR_DURATION and MKIAU_FREQUENCY are Initialized in GCM_GridComp
  ! ----------------------------------------------------------------------------
    call MAPL_GetResource( STATE, PREDICTOR_DURATION, Label="PREDICTOR_DURATION:", RC=STATUS ) ; VERIFY_(STATUS)
    call MAPL_GetResource( STATE, MKIAU_FREQUENCY,    Label="MKIAU_FREQUENCY:",    RC=STATUS ) ; VERIFY_(STATUS)

    call MAPL_GetResource( STATE, ANA_IS_WEIGHTED, Label="ANA_IS_WEIGHTED:", default='NO', RC=STATUS)
    VERIFY_(STATUS)
         ANA_IS_WEIGHTED = ESMF_UtilStringUpperCase(ANA_IS_WEIGHTED)
             IS_WEIGHTED =   adjustl(ANA_IS_WEIGHTED)=="YES" .or. adjustl(ANA_IS_WEIGHTED)=="NO"
    _ASSERT( IS_WEIGHTED ,'needs informative message')
             IS_WEIGHTED =   adjustl(ANA_IS_WEIGHTED)=="YES"


   ! Detect if running DasMode (Checking for AGCM_IMPORT)
   ! ----------------------------------------------------
    call MAPL_GetResource( STATE, STRING, LABEL="IMPORT_RESTART_FILE:", RC=STATUS)
    IF (STATUS == ESMF_SUCCESS) THEN
       DasMode = .true.
    ELSE
       DasMode = .false.
    END IF

   ! Detect if running REPLAY
   ! ------------------------
    call MAPL_GetResource( STATE, ReplayMode, 'REPLAY_MODE:', default="NoReplay", RC=STATUS )
    VERIFY_(STATUS)
       rplMode = adjustl(ReplayMode)
    if(rplMode=="Regular" .or. rplMode == "Exact") then
       DasMode = .true.
    end if


! Set type of update
!-------------------

    if     (ISFCST/=0  ) then
       TYPE = FORECAST
    else if(.not. DasMode) then
       TYPE = FREERUN
    else

    ! Get the specific IAU alarm
    !---------------------------
       call MAPL_StateAlarmGet(STATE, ALARM, NAME='PredictorAlarm', RC=STATUS)
       VERIFY_(STATUS)

       DO_PREDICTOR   = ESMF_AlarmIsRinging   ( ALARM, rc=status)
       VERIFY_(STATUS)
       LAST_CORRECTOR = ESMF_AlarmWillRingNext( ALARM, rc=status)
       VERIFY_(STATUS)

REPLAYING: if ( DO_PREDICTOR .and. (rplMode == "Regular") ) then
!-----------------------------------------------------------------------
               call ESMF_ClockGetAlarm(clock, 'startReplay', alarm, rc=status)
               VERIFY_(STATUS)
               LAST_CORRECTOR = ESMF_AlarmWillRingNext( ALARM, rc=status)
               VERIFY_(STATUS)

           else if(  (rplMode=="Exact")   .or.  &
                   ( (rplMode=="Regular") .and. (PREDICTOR_DURATION.gt.MKIAU_FREQUENCY/2) )  ) then

               ! Set Active PREDICTOR_STEP Alarm to OFF
               ! --------------------------------------
               call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
               VERIFY_(STATUS)
               DO_PREDICTOR = .FALSE.

               ! Get file template for READING Exact REPLAY Increment Files
               ! ----------------------------------------------------------
               if ( rplMode=="Exact" ) then
                    call MAPL_GetResource ( STATE, FileTmpl,  'REPLAY_FILE:',                   RC=STATUS )
                    VERIFY_(STATUS)
                    call MAPL_GetResource ( STATE, FileTmpl09,'REPLAY_FILE09:', DEFAULT='NULL', RC=STATUS )
                    VERIFY_(STATUS)
               endif

               ! Get file template for READING Regular REPLAY Increment Files produced during PREDICTOR Step
               ! (This should be consistent with any USER-Supplied MKIAU_CHECKPOINT file)
               ! -------------------------------------------------------------------------------------------
               if ( rplMode=="Regular" .and. (PREDICTOR_DURATION.gt.MKIAU_FREQUENCY/2) ) then
                    call MAPL_GetResource( STATE, FileName, "MKIAU_CHECKPOINT_FILE:", rc=status)
                    VERIFY_(STATUS)
                    call MAPL_GetResource( STATE, FileType, "MKIAU_CHECKPOINT_TYPE:", rc=status)
                    VERIFY_(STATUS)
                    if( FileType == 'binary' ) FileTmpl = trim(FileName) // '.%y4%m2%d2_%h2%n2z.' // 'bin'
                    if( FileType == 'pnc4'   ) FileTmpl = trim(Filename) // '.%y4%m2%d2_%h2%n2z.' // 'nc4'
                    FileTmpl09 = 'NULL'
               endif

! If replay alarm is ringing, we need to reset state
!---------------------------------------------------
               call ESMF_ClockGetAlarm(Clock,'ReplayShutOff',Alarm,rc=Status)
               VERIFY_(status)
               is_shutoff = ESMF_AlarmIsRinging( Alarm,rc=Status)
               VERIFY_(status)

               if (is_shutoff) then ! once this alarm rings, is_shutoff will remain true for the rest of the run
               !  if ( MAPL_am_I_root() ) print *, 'Zeroing AGCM_IMPORT'
                  call MAPL_GetPointer(IMPORT,ptr3d,'DUDT' ,RC=STATUS) ; ptr3d=0.0
                  call MAPL_GetPointer(IMPORT,ptr3d,'DVDT' ,RC=STATUS) ; ptr3d=0.0
                  call MAPL_GetPointer(IMPORT,ptr3d,'DTDT' ,RC=STATUS) ; ptr3d=0.0
                  call MAPL_GetPointer(IMPORT,ptr3d,'DPEDT',RC=STATUS) ; ptr3d=0.0
                  call MAPL_GetPointer(IMPORT,ptr3d,'DQVDT',RC=STATUS) ; ptr3d=0.0
                  call MAPL_GetPointer(IMPORT,ptr3d,'DO3DT',RC=STATUS) ; ptr3d=0.0
                  call MAPL_GetPointer(IMPORT,ptr2d,'DTSDT',RC=STATUS) ; if(associated(ptr2d)) ptr2d=0.0
               else
                  call ESMF_ClockGetAlarm(Clock,'ReplayShutOff',Alarm,rc=Status)
                  VERIFY_(status)
                  is_shutoff = ESMF_AlarmWillRingNext( Alarm,rc=status )
                  VERIFY_(status)
               endif

             ! Check for Beginning of REPLAY cycle
             ! -----------------------------------
               call ESMF_ClockGetAlarm(Clock,'replayCycle',Alarm,rc=Status)
               VERIFY_(status)
               Begin_REPLAY_Cycle = ESMF_AlarmIsRinging( Alarm,rc=status )
               VERIFY_(status)

             ! Check Alarm for Beginning of EXACT_REPLAY09 cycle
             ! -------------------------------------------------
               call ESMF_ClockGetAlarm(Clock,'ExactReplay09',Alarm,rc=Status)
               VERIFY_(status)
               is_ExactReplay09_ringing = ESMF_AlarmIsRinging( Alarm,rc=status )
               VERIFY_(status)

             ! Check Alarm for REGULAR_REPLAY09 cycle
             ! --------------------------------------
               call ESMF_ClockGetAlarm(Clock,'RegularReplay09',Alarm,rc=Status)
               VERIFY_(status)
               is_RegularReplay09_ringing = ESMF_AlarmIsRinging( Alarm,rc=status )
               VERIFY_(status)

             ! Check for Last Corrector
             ! -------------------------------------
               call ESMF_ClockGetAlarm(Clock,'ExactReplay',Alarm,rc=Status)
               VERIFY_(status)
               LAST_CORRECTOR = ESMF_AlarmWillRingNext( ALARM, rc=status)
               VERIFY_(STATUS)

! Force IS_RINGING to be TRUE at Start-Up
! ---------------------------------------
               if( first ) then
                   call ESMF_ClockGet(Clock, CurrTime=currTime, rc=Status)
                   VERIFY_(status)

                   call ESMF_TimeIntervalSet( TINT, S=INT(DT), rc=STATUS )
                   VERIFY_(STATUS)
                   REPLAY_TIME0 = currTime - TINT

                   first = .FALSE.
               endif
               call GET_REPLAY_TIME ( STATE, CLOCK, REPLAY_TIME, Begin_REPLAY_Cycle, RC )
               is_ringing = REPLAY_TIME /= REPLAY_TIME0

               is_ringing = ( is_ringing .or. Begin_REPLAY_Cycle ) .and. (.not. is_shutoff)

               if(is_ringing) then
               ! -----------------
                   ! Read REPLAY file
                   ! ----------------

                   REPLAY_TIME0 = REPLAY_TIME

                   if( rplMode=="Exact" ) then
                       if( is_ExactReplay09_ringing ) then
                           if( filetmpl09.ne.'NULL' ) then
                               call MAPL_GetCurrentFile(FILETMPL=filetmpl09, TIME=REPLAY_TIME, FILENAME=ReplayFile, RC=STATUS)
                               VERIFY_(status)
                           else
                               call MAPL_GetCurrentFile(FILETMPL=filetmpl,   TIME=REPLAY_TIME, FILENAME=ReplayFile, RC=STATUS)
                               VERIFY_(status)
                           endif
                       else
                               call MAPL_GetCurrentFile(FILETMPL=filetmpl,   TIME=REPLAY_TIME, FILENAME=ReplayFile, RC=STATUS)
                               VERIFY_(status)
                       endif
                   endif

                   if( (rplMode=="Regular") .and. (PREDICTOR_DURATION.gt.MKIAU_FREQUENCY/2) ) then
                       if( is_RegularReplay09_ringing ) then
                           if( filetmpl09.ne.'NULL' ) then
                               call MAPL_GetCurrentFile(FILETMPL=filetmpl09, TIME=REPLAY_TIME, FILENAME=ReplayFile, RC=STATUS)
                               VERIFY_(status)
                           else
                               call MAPL_GetCurrentFile(FILETMPL=filetmpl,   TIME=REPLAY_TIME, FILENAME=ReplayFile, RC=STATUS)
                               VERIFY_(status)
                           endif
                       else
                               call MAPL_GetCurrentFile(FILETMPL=filetmpl,   TIME=REPLAY_TIME, FILENAME=ReplayFile, RC=STATUS)
                               VERIFY_(status)
                       endif
                   endif

                   call MAPL_ESMFStateReadFromFile(STATE=IMPORT, CLOCK=CLOCK, FILENAME=ReplayFile, MPL=STATE, HDR=.FALSE., RC=STATUS)
                   VERIFY_(STATUS)
               endif

           end if REPLAYING

       if(DO_PREDICTOR) then
          TYPE = PREDICTOR
       else
          TYPE = CORRECTOR
       end if

    end if

! Get Names Associated with Friendly Analysis Bundle
!---------------------------------------------------

    call ESMF_StateGet(GEX(PHYS), 'TRANA', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet(BUNDLE,FieldCount=NumFriendly,   RC=STATUS)
    VERIFY_(STATUS)

    _ASSERT(NumFriendly==2,'needs informative message')

    allocate(Names(NumFriendly), stat=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet(BUNDLE, fieldNameList=Names, RC=STATUS)
    VERIFY_(STATUS)

! Prepare for update
!-------------------

    call MAPL_GetPointer( GEX(SDYN), PREF,'PREF',rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), PLE, 'PLE', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), U,   'U',   rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), V,   'V',   rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), W,   'W',   rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), T,   'T',   rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), AK,  'AK',  rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), BK,  'BK',  rc=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetPointer( GEX(SDYN), PEPHY_SDYN, 'PEPHY', alloc=.true., rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(PHYS), PEPHY_PHYS, 'PEPHY', alloc=.true., rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), DTHVDTPHYINT, 'DTHVDTPHYINT', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(PHYS), DQVDTPHYINT,  'DQVDTPHYINT',  rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(PHYS), DQLDTPHYINT,  'DQLDTPHYINT',  rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(PHYS), DQIDTPHYINT,  'DQIDTPHYINT',  rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(PHYS), DOXDTPHYINT,  'DOXDTPHYINT',  rc=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetPointer( GEX(SDYN), AREA, 'AREA', rc=STATUS )
    VERIFY_(STATUS)

    allocate( PL(IM,JM,LM),STAT=STATUS )
    VERIFY_(STATUS)
    allocate( DP(IM,JM,LM),STAT=STATUS )
    VERIFY_(STATUS)

    PL  = 0.5*(PLE(:,:,1:LM)+PLE(:,:,0:LM-1))
    DP  =      PLE(:,:,1:LM)-PLE(:,:,0:LM-1)

! --------------------------------------------------------
! ALPHA and BETA BIAS Correction Coefficients
! --------------------------------------------------------
! AGCM_IMPORT   =>  IAU(n)     ALPHA = 0 (Default)
! AGCM_INTERNAL => BIAS(n)     BETA  = 1 (Default)
! BIAS Update:     BIAS(n+1) = ALPHA*IAU(n) + BETA*BIAS(n)
! NOTE:  After LAST_CORRECTOR, the Model (in PREDICTOR mode) is
!        forced by AGCM_INTERNAL (BIAS).
!        Under DEFAULT conditions, the BIAS = 0
! --------------------------------------------------------
    ALF = ALPHA
    BET = BETA

! Get IAU Scaling Coefficient
! ---------------------------
    if( TYPE /= CORRECTOR ) then
        IAUcoeff = 1.0    ! Do NOT modify Forecast/Predictor forcing term
    else
        call get_iau_coeff( IAUcoeff,CLOCK )

      ! If 4DIAU, overwrite increments from analysis by recreating them on the fly
      ! --------------------------------------------------------------------------
        DO_4DIAU = ESMF_AlarmIsRinging( ALARM4D, rc=status)
        VERIFY_(STATUS)
        if(DO_4DIAU .and. (TYPE == CORRECTOR) ) then
           call ESMF_AlarmRingerOff(ALARM4D, RC=STATUS)
           VERIFY_(STATUS)
           call update_ainc_(RC=STATUS)
           VERIFY_(STATUS)
        endif

    endif

! Load Analysis Increments into Imports for Physics and Dynamics State Variables
!-------------------------------------------------------------------------------

    call DO_UPDATE_ANA2D ('DTSDT', PHYS)

! Note: Mass-Weighting is done for DTDT, No Mass-Weighting for DUDT,DVDT,DPEDT
! ----------------------------------------------------------------------------

    call DO_UPDATE_ANA3D ('DUDT' , SDYN, PREF)
    call DO_UPDATE_ANA3D ('DVDT' , SDYN, PREF)
    call DO_UPDATE_ANA3D ('DPEDT', SDYN, PREF, CONSTRAIN_DAS = CONSTRAIN_DAS)

    call MAPL_GetPointer(GIM(SDYN), DPEDT, 'DPEDT', rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GIM(SDYN), DTDT,  'DTDT',  rc=STATUS)
    VERIFY_(STATUS)

    call DO_UPDATE_ANA3D ('DTDT' , SDYN, PREF)

    if( is_weighted ) then
        allocate(  tdpold( IM,JM,  LM),STAT=STATUS ) ; VERIFY_(STATUS)
        allocate(  tdpnew( IM,JM,  LM),STAT=STATUS ) ; VERIFY_(STATUS)
        allocate(  dp_ana( IM,JM,  LM),STAT=STATUS ) ; VERIFY_(STATUS)
        allocate( ple_ana( IM,JM,0:LM),STAT=STATUS ) ; VERIFY_(STATUS)

        ! Create Proxies for Updated Pressure and Temperature due to Analysis
        !--------------------------------------------------------------------
        ple_ana = ple + dt*dpedt
         dp_ana = ple_ana(:,:,1:LM)-ple_ana(:,:,0:LM-1)

        tdpold = T*DP
        tdpnew = ( T + dt*dtdt )*dp_ana
          dtdt = ( tdpnew - tdpold )/dt

        deallocate(  tdpold )
        deallocate(  tdpnew )
        deallocate(  dp_ana )
        deallocate( ple_ana )
    endif

! Add Analysis Increment Directly to Friendlies
!----------------------------------------------

    if(TYPE /= FREERUN) then

       allocate(zero(IM,JM,LM),stat=status)
       VERIFY_(status)
       zero = 0.0

       allocate(FC (IM,JM,LM),STAT=STATUS)
       VERIFY_(STATUS)

       do K=1,NumFriendly

          FC = 1.0

          NULLIFY(Q)  !ALT: ESMF requires that the data pointer is not associated
          call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
          VERIFY_(STATUS)

          STRING = TRIM(Names(K))
	  fieldName = MAPL_RmQualifier(STRING)

          if(TRIM(fieldName) == 'OX') then ! PCHEM OX or GOCART::OX, for example.

             ALF = ALPHAO
             BET = BETAO
             DTX = DT*1.0E-6

! Uncomment damping if problem with ozone
! ---------------------------------------
!            do L=1,LM
!               where(PL(:,:,L) < 100.0 .and. PL(:,:,L) > 0.0 )
!                     FC(:,:,L) = exp(-1.5*(log10(PL(:,:,L))-2.0)**2)
!               end where
!            end do

             call MAPL_GetPointer(GIM(SDYN), DOXANA, 'DOXANA', rc=STATUS)
             VERIFY_(STATUS)
             DOXANA = Q
             call DO_Friendly (Q,'DO3DT',PREF)
             DOXANA = Q - DOXANA

          else

             ALF = ALPHAQ
             BET = BETAQ
             DTX = DT

             if(NAMES(K)=='Q') then ! Q

                ! Initialize DQVANA Diagnostics with Background QV
                !-------------------------------------------------
                call MAPL_GetPointer(GIM(SDYN), DQVANA, 'DQVANA', rc=STATUS)
                VERIFY_(STATUS)
                call MAPL_GetPointer(GIM(SDYN), DQLANA, 'DQLANA', rc=STATUS)
                VERIFY_(STATUS)
                call MAPL_GetPointer(GIM(SDYN), DQIANA, 'DQIANA', rc=STATUS)
                VERIFY_(STATUS)
                call MAPL_GetPointer(GIM(SDYN), DQRANA, 'DQRANA', rc=STATUS)
                VERIFY_(STATUS)
                call MAPL_GetPointer(GIM(SDYN), DQSANA, 'DQSANA', rc=STATUS)
                VERIFY_(STATUS)
                call MAPL_GetPointer(GIM(SDYN), DQGANA, 'DQGANA', rc=STATUS)
                VERIFY_(STATUS)

                ! Get Pointers to QL & QI from Friendly Advection Bundle
                !-------------------------------------------------------
                call ESMF_StateGet(GEX(PHYS), 'TRADV', Advect_BUNDLE, RC=STATUS )
                VERIFY_(STATUS)

                call ESMFL_BundleGetPointerToData( Advect_BUNDLE, 'QLLS', QLLS, RC=STATUS )
                VERIFY_(STATUS)
                call ESMFL_BundleGetPointerToData( Advect_BUNDLE, 'QLCN', QLCN, RC=STATUS )
                VERIFY_(STATUS)
                call ESMFL_BundleGetPointerToData( Advect_BUNDLE, 'QILS', QILS, RC=STATUS )
                VERIFY_(STATUS)
                call ESMFL_BundleGetPointerToData( Advect_BUNDLE, 'QICN', QICN, RC=STATUS )
                VERIFY_(STATUS)
                call ESMF_FieldBundleGet (Advect_BUNDLE, fieldName='QRAIN', isPresent=isPresent, RC=STATUS)
                VERIFY_(STATUS)
                if (isPresent) then
                   call ESMFL_BundleGetPointerToData( Advect_BUNDLE, 'QRAIN', QRAIN, RC=STATUS )
                   VERIFY_(STATUS)
                else
                   QRAIN => zero
                end if
                call ESMF_FieldBundleGet (Advect_BUNDLE, fieldName='QSNOW', isPresent=isPresent, RC=STATUS)
                VERIFY_(STATUS)
                if (isPresent) then
                   call ESMFL_BundleGetPointerToData( Advect_BUNDLE, 'QSNOW', QSNOW, RC=STATUS )
                   VERIFY_(STATUS)
                else
                   QSNOW => zero
                end if
                call ESMF_FieldBundleGet (Advect_BUNDLE, fieldName='QGRAUPEL', isPresent=isPresent, RC=STATUS)
                VERIFY_(STATUS)
                if (isPresent) then
                   call ESMFL_BundleGetPointerToData( Advect_BUNDLE, 'QGRAUPEL', QGRAUPEL, RC=STATUS )
                   VERIFY_(STATUS)
                else
                   QGRAUPEL => zero
                end if

                call MAPL_GetPointer ( EXPORT, ptr3d, 'QVBKG', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(ptr3d)) ptr3d = Q
                call MAPL_GetPointer ( EXPORT, ptr3d, 'DPBKG', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(ptr3d)) ptr3d = DP

                DQVANA = Q
                DQLANA = QLLS + QLCN
                DQIANA = QILS + QICN
                DQRANA = QRAIN
                DQSANA = QSNOW
                DQGANA = QGRAUPEL

                if(TYPE  == CORRECTOR) then

                    ! ---------------------------------------------------
                    ! Create BKG Water Variables
                    ! --------------------------
                    IF( CONSTRAIN_DAS == 1 ) then
                       allocate( qdp_bkg( IM,JM,LM ),STAT=STATUS ) ; VERIFY_(STATUS)
                       qdp_bkg = (         q+qlls+qlcn+qils+qicn+qrain+qsnow+qgraupel   ) * dp
#if debug
                       allocate( dry_bkg( IM,JM,LM ),STAT=STATUS ) ; VERIFY_(STATUS)
                       dry_bkg = ( 1.0 - ( q+qlls+qlcn+qils+qicn+qrain+qsnow+qgraupel ) ) * dp
#endif
                    ENDIF
                    IF( CONSTRAIN_DAS == 2 ) then
                       allocate( dry_bkg( IM,JM,LM ),STAT=STATUS ) ; VERIFY_(STATUS)
                       dry_bkg = ( 1.0 - ( q+qlls+qlcn+qils+qicn+qrain+qsnow+qgraupel ) ) * dp
                    ENDIF

                ENDIF ! End CORRECTOR Test

                ! -----------------------------------------------------------------------------------------
                ! -----------------------------------------------------------------------------------------

                call DO_Friendly (Q,'DQVDT',PREF)

                ! -----------------------------------------------------------------------------------------
                ! -----------------------------------------------------------------------------------------

                ! Create Proxies for Updated Pressure due to Analysis
                !----------------------------------------------------
                allocate(  dp_ana(IM,JM,  LM),STAT=STATUS ) ; VERIFY_(STATUS)
                allocate( ple_ana(IM,JM,0:LM),STAT=STATUS ) ; VERIFY_(STATUS)

                allocate( DQVCON(IM,JM,LM),STAT=STATUS ) ; VERIFY_(STATUS)
                DQVCON = Q                               ! Initialize Constraint Tendency

                ! Proxies for Pressure Changes due to Analysis
                ! --------------------------------------------
                ple_ana = ple + dt*dpedt
                 dp_ana  = ple_ana(:,:,1:LM) - ple_ana(:,:,0:LM-1)

                call MAPL_GetPointer ( EXPORT, ptr3d, 'QVANA', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(ptr3d)) ptr3d = Q
                call MAPL_GetPointer ( EXPORT, ptr3d, 'DPANA', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(ptr3d)) ptr3d = DP_ANA

                if(TYPE  == CORRECTOR) then

                    IF( CONSTRAIN_DAS == 1 ) then
                    ! ---------------------------

                    ! Create ANA Water Mass
                    ! ---------------------
                       allocate(  qdp_ana(IM,JM,LM),STAT=STATUS ); VERIFY_(STATUS)
                       qdp_ana = (         q+qlls+qlcn+qils+qicn+qrain+qsnow+qgraupel   ) * dp_ana

                    ! Vertically Integrate ANA & BKG Water Mass where they Differ
                    ! -----------------------------------------------------------
                       allocate( sum_qdp_bkg( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)
                       allocate( sum_qdp_ana( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)
                       sum_qdp_bkg = 0.0_8
                       sum_qdp_ana = 0.0_8
                       do L=1,lm
                       where( ABS(qdp_ana(:,:,L)-qdp_bkg(:,:,L)) > tiny(1.0_4) )
                              sum_qdp_bkg = sum_qdp_bkg + qdp_bkg(:,:,L)
                              sum_qdp_ana = sum_qdp_ana + qdp_ana(:,:,L)
                       end where
                       enddo

                    ! Compute Area-Mean Vertically Integrated BKG Water Mass
                    ! ------------------------------------------------------
                       allocate( qdp_bkg_int( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)
                       where( sum_qdp_bkg.ne.0.0_8 )
                           qdp_bkg_int = sum_qdp_bkg
                       elsewhere
                           qdp_bkg_int = MAPL_UNDEF
                       end where
                       call MAPL_AreaMean( qdp_bkg_ave, qdp_bkg_int, area, grid, rc=STATUS )
                       VERIFY_(STATUS)

                    ! Compute Area-Mean Vertically Integrated ANA Water Mass
                    ! ------------------------------------------------------
                       allocate( qdp_ana_int( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)
                       where( sum_qdp_ana.ne.0.0_8 )
                           qdp_ana_int = sum_qdp_ana
                       elsewhere
                           qdp_ana_int = MAPL_UNDEF
                       end where
                       call MAPL_AreaMean( qdp_ana_ave, qdp_ana_int, area, grid, rc=STATUS )
                       VERIFY_(STATUS)

                    ! Compute Dry-Mass Scaling Parameter
                    ! ----------------------------------
                       if( qdp_bkg_ave.ne.MAPL_UNDEF .and. &
                           qdp_ana_ave.ne.MAPL_UNDEF       ) then
                           gamma = qdp_bkg_ave / qdp_ana_ave                  ! Prefered Method
                         ! gamma = real( qdp_bkg_ave / qdp_ana_ave, kind=4 )  ! Method for Zero-diff Backward Compatibility
                       else
                           gamma = 1.0_8
                       endif

                    ! Scale ANA Water Variables for Dry-Mass Conservation
                    ! ---------------------------------------------------
                       do L=1,lm
                       where( ABS(qdp_ana(:,:,L)-qdp_bkg(:,:,L)) > tiny(1.0_4) )
                                 q   (:,:,L) =     q   (:,:,L) * gamma
                                 qlls(:,:,L) =     qlls(:,:,L) * gamma
                                 qlcn(:,:,L) =     qlcn(:,:,L) * gamma
                                 qils(:,:,L) =     qils(:,:,L) * gamma
                                 qicn(:,:,L) =     qicn(:,:,L) * gamma
                                qrain(:,:,L) =    qrain(:,:,L) * gamma
                                qsnow(:,:,L) =    qsnow(:,:,L) * gamma
                             qgraupel(:,:,L) = qgraupel(:,:,L) * gamma
                       end where
                       enddo
#if debug
                       i=1
                       call Dry_Mass_Check (im,jm,lm,i)
                       deallocate( dry_bkg )
#endif
                       deallocate( qdp_bkg )
                       deallocate( qdp_ana )
                       deallocate( qdp_bkg_int )
                       deallocate( qdp_ana_int )
                       deallocate( sum_qdp_bkg )
                       deallocate( sum_qdp_ana )

                    ENDIF  ! End CONSTRAIN_DAS ==1 Test

                    IF( CONSTRAIN_DAS == 2 ) then  ! Constrain Dry_Mass Conservation using Least-Squares of P
                    ! ---------------------------------------------------------------------------------------

                       call MAPL_GetPointer ( EXPORT, ptr2d, 'DPSDT_CON', rc=STATUS )
                       VERIFY_(STATUS)
                       if(associated(ptr2d)) ptr2d = ple_ana(:,:,LM)   ! Initialize Constraint Tendency

                       do i=1,10

                    ! Create ANA Dry Mass
                    ! -------------------
                       allocate( dry_ana(IM,JM,LM),STAT=STATUS ); VERIFY_(STATUS)
                       dry_ana = ( 1.0 - ( q+qlls+qlcn+qils+qicn+qrain+qsnow+qgraupel ) ) * dp_ana

                    ! allocate(  ALPHA2D(IM,JM),STAT=STATUS ) ; VERIFY_(STATUS)
                    !            ALPHA2D = 1.0

                    ! Vertically Integrate ANA & BKG Dry Mass
                    ! ---------------------------------------
                       allocate( sum_dry_bkg( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)
                       allocate( sum_dry_ana( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)

                       sum_dry_bkg = 0.0_8
                       sum_dry_ana = 0.0_8
                       do L=1,lm
                              sum_dry_bkg = sum_dry_bkg + dry_bkg(:,:,L)
                              sum_dry_ana = sum_dry_ana + dry_ana(:,:,L)
                       enddo

                    ! Compute Area-Mean Vertically Integrated BKG Water Mass
                    ! ------------------------------------------------------
                       allocate( alf_bks_int( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)
                       allocate( dry_bkg_int( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)
                       allocate( dry_ana_int( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)

                            dry_bkg_int =   sum_dry_bkg
                            dry_ana_int =   sum_dry_ana
                            alf_bks_int = ( sum_dry_ana/ple_ana(:,:,LM) )**2 ! * ALPHA2D    Multiply for generic ALPHA2D NE 1.0

                       call MAPL_AreaMean( alf_bks_ave, alf_bks_int, area, grid, rc=STATUS )
                       VERIFY_(STATUS)
                       call MAPL_AreaMean( dry_bkg_ave, dry_bkg_int, area, grid, rc=STATUS )
                       VERIFY_(STATUS)
                       call MAPL_AreaMean( dry_ana_ave, dry_ana_int, area, grid, rc=STATUS )
                       VERIFY_(STATUS)

                    ! Compute Dry-Mass Constraint Parameter
                    ! -------------------------------------
                                   gamma = ( dry_ana_ave - dry_bkg_ave ) / alf_bks_ave
                         ple_ana(:,:,LM) = ple_ana(:,:,LM) - gamma * sum_dry_ana/ple_ana(:,:,LM) ! * ALPHA2D    Multiply for generic ALPHA2D NE 1.0

                         do L=0,LM-1
                         ple_ana(:,:,L) = AK(L) + BK(L)*ple_ana(:,:,LM)
                         enddo
                         dpedt = ( ple_ana-ple )/dt

#if debug
                       call Dry_Mass_Check (im,jm,lm,i)
#endif

                    !  deallocate( ALPHA2D )
                       deallocate( dry_ana )
                       deallocate( alf_bks_int )
                       deallocate( dry_bkg_int )
                       deallocate( dry_ana_int )
                       deallocate( sum_dry_bkg )
                       deallocate( sum_dry_ana )

                       enddo

                       call MAPL_GetPointer ( EXPORT, ptr2d, 'DPSDT_CON', rc=STATUS )
                       VERIFY_(STATUS)
                       if(associated(ptr2d)) ptr2d = ( ple_ana(:,:,LM) - ptr2d )/DT

                       deallocate( dry_bkg )

                    ENDIF  ! End CONSTRAIN_DAS == 2 Test

                ENDIF  ! End CORRECTOR Test

                call MAPL_GetPointer ( EXPORT, ptr3d, 'QVCON', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(ptr3d)) ptr3d = Q

                DQVCON = Q           - DQVCON
                DQVANA = Q           - DQVANA
                DQLANA = QLLS + QLCN - DQLANA
                DQIANA = QILS + QICN - DQIANA
                DQRANA = QRAIN       - DQRANA
                DQSANA = QSNOW       - DQSANA
                DQGANA = QGRAUPEL    - DQGANA

                ! Update Tendency Diagnostic due to CONSTRAINTS
                ! ---------------------------------------------
                call MAPL_GetPointer ( EXPORT, TENDAN, 'DQVDT_CON', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(TENDAN)) TENDAN = DQVCON/DT
                deallocate( DQVCON )

                call MAPL_GetPointer ( EXPORT, TENDAN, 'DQVDT_ANA', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(TENDAN)) TENDAN = DQVANA/DT

                call MAPL_GetPointer ( EXPORT, TENDAN, 'DQLDT_ANA', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(TENDAN)) TENDAN = DQLANA/DT

                call MAPL_GetPointer ( EXPORT, TENDAN, 'DQIDT_ANA', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(TENDAN)) TENDAN = DQIANA/DT

                call MAPL_GetPointer ( EXPORT, TENDAN, 'DQRDT_ANA', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(TENDAN)) TENDAN = DQRANA/DT

                call MAPL_GetPointer ( EXPORT, TENDAN, 'DQSDT_ANA', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(TENDAN)) TENDAN = DQSANA/DT

                call MAPL_GetPointer ( EXPORT, TENDAN, 'DQGDT_ANA', rc=STATUS )
                VERIFY_(STATUS)
                if(associated(TENDAN)) TENDAN = DQGANA/DT


                deallocate(  dp_ana )
                deallocate( ple_ana )

             else

                call DO_Friendly (Q,'D'//trim(Names(K))//'DT',PREF)

             end if ! End Test for Q Friendly

          end if  ! End Test for OX Friendly

       end do  ! End Friendly Loop List
       deallocate(FC)
       deallocate (zero)

    end if ! not free-running

! Update Total Pressure Tendency due to Analysis + Constraint
! -----------------------------------------------------------
      call MAPL_GetPointer ( EXPORT, ptr3d, 'DPEDT_ANA', rc=STATUS )
      VERIFY_(STATUS)
      if(associated(ptr3d)) ptr3d = dpedt

! Make Sure EPV is Allocated for TROPOPAUSE Diagnostics
!------------------------------------------------------
    call MAPL_GetPointer ( EXPORT, TROPP1, 'TROPP_THERMAL', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TROPP2, 'TROPP_EPV'    , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TROPP3, 'TROPP_BLENDED', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TROPT, 'TROPT', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TROPQ, 'TROPQ', rc=STATUS )
    VERIFY_(STATUS)

    if( associated(TROPP1) .or. &
        associated(TROPP2) .or. &
        associated(TROPP3) .or. &
        associated(TROPT)  .or. &
        associated(TROPQ)       ) then
        call MAPL_GetPointer( GEX(SDYN),EPV,'EPV',ALLOC=.true.,rc=STATUS )
        VERIFY_(STATUS)
     endif

! Initialize TRADVI bundle with TRADV bundle
!--------------------------------------------
    call Initialize_IncBundle_run(GEX(PHYS), EXPORT, DYNinc, __RC__)

! Call basic run phase for both Child
!-------------------------------------

! Call run for satellite orbits
!------------------------------
    call MAPL_TimerOn (STATE,"ORBIT"  )
    call ESMF_GridCompRun(GCS(ORB), importState=GIM(ORB), exportState=GEX(ORB), clock=CLOCK, userRC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,"ORBIT"  )

    call Pack_Chem_Groups( GEX(PHYS) )  ! Prepare to transport chemical families

! ! Call system clock to estmiate Dyn:Phy ratio
!   call SYSTEM_CLOCK(COUNT_MAX=COUNT_MAX)

!   call SYSTEM_CLOCK(START_TIME)
    call MAPL_TimerOn (STATE,"SUPERDYNAMICS"  )
    call ESMF_GridCompRun(GCS(SDYN), importState=GIM(SDYN), exportState=GEX(SDYN), clock=CLOCK, PHASE=1, userRC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerOFF (STATE,"SUPERDYNAMICS"  )
!   call SYSTEM_CLOCK(END_TIME)
!   DYN_TIME = END_TIME-START_TIME
!   if(DYN_TIME<0) then
!      DYN_TIME = DYN_TIME + COUNT_MAX
!   endif

! Compute Tracer Advection increments
!-------------------------------------
    call Compute_IncBundle(GEX(PHYS), EXPORT, DYNinc, STATE, __RC__)

    call Unpack_Chem_Groups( GEX(PHYS), PLE, AREA )  ! Finish transporting chemical families

!   call SYSTEM_CLOCK(START_TIME)
    call MAPL_TimerOn (STATE,"PHYSICS"  )
    call ESMF_GridCompRun(GCS(PHYS), importState=GIM(PHYS), exportState=GEX(PHYS), clock=CLOCK, PHASE=1, userRC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,"PHYSICS"  )
!   call SYSTEM_CLOCK(END_TIME)
!   PHY_TIME = END_TIME-START_TIME
!   if(PHY_TIME<0) then
!      PHY_TIME = PHY_TIME + COUNT_MAX
!   endif

!   if( MAPL_am_I_root() ) write(6,1000) REAL(DYN_TIME,kind=8)/REAL(PHY_TIME,kind=8)
!   1000 format(1x,'DYN:PHY Ratio: ',f7.2)

! Load Physics Tendencies into Imports for RUN2 of Dynamics (ADD_INCS)
!---------------------------------------------------------------------

    call DO_UPDATE_PHY ('DUDT' )
    call DO_UPDATE_PHY ('DVDT' )
    call DO_UPDATE_PHY ('DWDT' )
    call DO_UPDATE_PHY ('DPEDT')
    call DO_UPDATE_PHY ('DTDT' )

    call MAPL_TimerOn (STATE,"AGCM_BARRIER"  )
    call ESMF_VMBarrier(VMG, rc=status); VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,"AGCM_BARRIER"  )

! Run RUN2 of SuperDynamics (ADD_INCS) to add Physics Diabatic Tendencies
!------------------------------------------------------------------------

    call MAPL_TimerOn (STATE,"SUPERDYNAMICS"  )

    call ESMF_GridCompRun(GCS(SDYN), importState=GIM(SDYN), exportState=GEX(SDYN), clock=CLOCK, PHASE=2, userRC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GenericRunCouplers( STATE, SDYN, CLOCK, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOff(STATE,"SUPERDYNAMICS"  )

! Get Names Associated with Friendly Advection Bundle for Final Check for Negative Tracers
!-----------------------------------------------------------------------------------------

    call ESMF_StateGet(GEX(PHYS), 'TRADV', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet(BUNDLE,FieldCount=NumFriendly,   RC=STATUS)
    VERIFY_(STATUS)

    deallocate(Names)
      allocate(Names(NumFriendly), stat=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet(BUNDLE, fieldNameList=Names, RC=STATUS)
    VERIFY_(STATUS)

! Get Pointers to Exports
!------------------------
    call MAPL_GetPointer ( EXPORT, QTFILL, 'QTFILL', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, QVFILL, 'QVFILL', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, QIFILL, 'QIFILL', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, QLFILL, 'QLFILL', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, OXFILL, 'OXFILL', rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TOX   , 'TOX'   , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TQV   , 'TQV'   , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TQI   , 'TQI'   , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, TQL   , 'TQL'   , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, QLTOT , 'QLTOT' , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, QITOT , 'QITOT' , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, PERES       , 'PERES'        , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, PEFILL      , 'PEFILL'       , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, DTHVDTFILINT, 'DTHVDTFILINT' , rc=STATUS )
    VERIFY_(STATUS)

    if(associated(QTFILL)) QTFILL = 0.0
    if(associated(QIFILL)) QIFILL = 0.0
    if(associated(QLFILL)) QLFILL = 0.0
    if(associated(TQI)   ) TQI    = 0.0
    if(associated(TQL)   ) TQL    = 0.0
    if(associated(QLTOT) ) QLTOT  = 0.0
    if(associated(QITOT) ) QITOT  = 0.0

    allocate(QFILL(IM,JM)    ,STAT=STATUS )
    VERIFY_(STATUS)
    allocate(QINT (IM,JM)    ,STAT=STATUS )
    VERIFY_(STATUS)
    allocate( PKE(IM,JM,0:LM),STAT=STATUS )
    VERIFY_(STATUS)
    allocate( PKZ(IM,JM,1:LM),STAT=STATUS )
    VERIFY_(STATUS)

    PL  = 0.5*(PLE(:,:,1:LM)+PLE(:,:,0:LM-1))  ! Recompute Updated Pressure
    DP  =      PLE(:,:,1:LM)-PLE(:,:,0:LM-1)   ! Recompute Updated Pressure Thickness

    PKE = PLE**MAPL_KAPPA
    do L=1,LM
    PKZ(:,:,L) = ( PKE(:,:,L)-PKE(:,:,L-1) ) / ( MAPL_KAPPA*( log(PLE(:,:,L))-log(PLE(:,:,L-1)) ) )
    enddo

! Initialize Vertically Integrated Values of CPT and THV (before QFILL updates)
! -----------------------------------------------------------------------------
    if( associated(PEPHY_SDYN) .or. associated(PEFILL) ) then
        EPS = MAPL_RVAP/MAPL_RGAS-1.0
        do K=1,NumFriendly
           NULLIFY(Q)
           call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
           VERIFY_(STATUS)
           if(NAMES(K)=='Q') then
              allocate( SUMCPT1(IM,JM),STAT=STATUS )
              VERIFY_(STATUS)
              SUMCPT1 = 0.0
              do L=1,LM
              SUMCPT1 = SUMCPT1 + MAPL_CP*T(:,:,L)*(1.0+EPS*Q(:,:,L))*DP(:,:,L)
              enddo
              exit
           endif
        enddo
    endif

    if( associated(DTHVDTPHYINT) .or. associated(DTHVDTFILINT) ) then
        EPS = MAPL_RVAP/MAPL_RGAS-1.0
        do K=1,NumFriendly
           NULLIFY(Q)
           call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
           VERIFY_(STATUS)
           if(NAMES(K)=='Q') then
              allocate( SUMTHV1(IM,JM),STAT=STATUS )
              VERIFY_(STATUS)
              SUMTHV1 = 0.0
              do L=1,LM
              SUMTHV1 = SUMTHV1 + T(:,:,L)/PKZ(:,:,L)*(1.0+EPS*Q(:,:,L))*DP(:,:,L)
              enddo
              exit
           endif
        enddo
    endif

! Perform Final Check for Negative Friendlies
! -------------------------------------------

    do K=1,NumFriendly
       NULLIFY(Q)
       call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
       VERIFY_(STATUS)

       STRING = TRIM(Names(K))
       fieldName = MAPL_RmQualifier(STRING)

! Water Vapor
! -----------
       if(NAMES(K)=='Q') then
          call FILL_Friendly   ( Q,DP,QFILL,QINT )
          if(associated(QVFILL))           QVFILL =               QFILL
          if(associated(QTFILL))           QTFILL = QTFILL      + QFILL
          if(associated(DQVDTPHYINT)) DQVDTPHYINT = DQVDTPHYINT + QFILL
          if(associated(TQV))                 TQV = QINT
       endif

! Ice Water
! ---------
       if(NAMES(K)=='QICN') then
          call FILL_Friendly   ( Q,DP,QFILL,QINT )
          if(associated(QIFILL))           QIFILL = QIFILL      + QFILL
          if(associated(QTFILL))           QTFILL = QTFILL      + QFILL
          if(associated(DQIDTPHYINT)) DQIDTPHYINT = DQIDTPHYINT + QFILL
          if(associated(TQI))                 TQI = TQI         + QINT
          if(associated(QITOT))             QITOT = QITOT       + Q
       endif

       if(NAMES(K)=='QILS') then
          call FILL_Friendly   ( Q,DP,QFILL,QINT )
          if(associated(QIFILL))           QIFILL = QIFILL      + QFILL
          if(associated(QTFILL))           QTFILL = QTFILL      + QFILL
          if(associated(DQIDTPHYINT)) DQIDTPHYINT = DQIDTPHYINT + QFILL
          if(associated(TQI))                 TQI = TQI         + QINT
          if(associated(QITOT))             QITOT = QITOT       + Q
       endif

! Liquid Water
! ------------
       if(NAMES(K)=='QLCN') then
          call FILL_Friendly   ( Q,DP,QFILL,QINT )
          if(associated(QLFILL))           QLFILL = QLFILL      + QFILL
          if(associated(QTFILL))           QTFILL = QTFILL      + QFILL
          if(associated(DQLDTPHYINT)) DQLDTPHYINT = DQLDTPHYINT + QFILL
          if(associated(TQL))                 TQL = TQL         + QINT
          if(associated(QLTOT))             QLTOT = QLTOT       + Q
       endif

       if(NAMES(K)=='QLLS') then
          call FILL_Friendly   ( Q,DP,QFILL,QINT )
          if(associated(QLFILL))           QLFILL = QLFILL      + QFILL
          if(associated(QTFILL))           QTFILL = QTFILL      + QFILL
          if(associated(DQLDTPHYINT)) DQLDTPHYINT = DQLDTPHYINT + QFILL
          if(associated(TQL))                 TQL = TQL         + QINT
          if(associated(QLTOT))             QLTOT = QLTOT       + Q
       endif

! Total Odd-Oxygen
! ----------------
       if(TRIM(fieldName) == 'OX') then
          call FILL_Friendly   ( Q,DP,QFILL,QINT )
          if(associated(OXFILL))           OXFILL = QFILL              *(MAPL_O3MW/MAPL_AIRMW)
          if(associated(DOXDTPHYINT)) DOXDTPHYINT = DOXDTPHYINT + QFILL*(MAPL_O3MW/MAPL_AIRMW)
          if(associated(TOX))                 TOX = QINT               *(MAPL_O3MW/MAPL_AIRMW)
       endif

    enddo
    deallocate(QFILL)
    deallocate(QINT )


! Compute Additional Diagnostics
!-------------------------------

    call MAPL_GetPointer ( EXPORT, MASS , 'MASS' , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, KE   , 'KE'   , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, CPT  , 'CPT'  , rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, THV  , 'THV'  , rc=STATUS )
    VERIFY_(STATUS)

    if( associated(MASS) ) MASS = (PLE(:,:,LM)-PLE(:,:,0)) * (1.0/MAPL_GRAV)

    if( associated(KE) ) then
        allocate( SUMKE(IM,JM),STAT=STATUS )
        VERIFY_(STATUS)
        SUMKE = 0.0
        do L=1,LM
        SUMKE = SUMKE + 0.5*( U(:,:,L)**2 + V(:,:,L)**2 )*DP(:,:,L)
        enddo
           KE = SUMKE * (1.0/MAPL_GRAV)
        deallocate(SUMKE)
    endif

! Instantaneous Values of CPT and THV are done here to include possible QFILL updates
! -----------------------------------------------------------------------------------
    if( associated(CPT) .or. associated(PEPHY_SDYN) .or. associated(PEFILL) ) then
        EPS = MAPL_RVAP/MAPL_RGAS-1.0
        do K=1,NumFriendly
           NULLIFY(Q)
           call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
           VERIFY_(STATUS)
           if(NAMES(K)=='Q') then
              allocate( SUMCPT2(IM,JM),STAT=STATUS )
              VERIFY_(STATUS)
              SUMCPT2 = 0.0
              do L=1,LM
              SUMCPT2 = SUMCPT2 + MAPL_CP*T(:,:,L)*(1.0+EPS*Q(:,:,L))*DP(:,:,L)
              enddo
              if( associated(CPT)         )       CPT = SUMCPT2 * (1.0/MAPL_GRAV)
              if( associated(PEFILL) .or. &
                  associated(PEPHY_SDYN)  )   then
                                              SUMCPT1 = ( SUMCPT2-SUMCPT1 ) / (DT*MAPL_GRAV)
              if( associated(PEFILL) )         PEFILL = SUMCPT1
              if( associated(PEPHY_SDYN) ) PEPHY_SDYN = PEPHY_SDYN + SUMCPT1
              deallocate(SUMCPT1)
              endif
              deallocate(SUMCPT2)
              exit
           endif
        enddo
    endif

    if( associated(PERES)       .and.  &
        associated(PEPHY_SDYN)  .and.  &
        associated(PEPHY_PHYS) ) PERES = PEPHY_SDYN - PEPHY_PHYS

    if( associated(THV) .or. associated(DTHVDTPHYINT) .or. associated(DTHVDTFILINT) ) then
        EPS = MAPL_RVAP/MAPL_RGAS-1.0
        do K=1,NumFriendly
           NULLIFY(Q)
           call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
           VERIFY_(STATUS)
           if(NAMES(K)=='Q') then
              allocate( SUMTHV2(IM,JM),STAT=STATUS )
              VERIFY_(STATUS)
              SUMTHV2 = 0.0
              do L=1,LM
              SUMTHV2 = SUMTHV2 + T(:,:,L)/PKZ(:,:,L)*(1.0+EPS*Q(:,:,L))*DP(:,:,L)
              enddo
              if( associated(THV)          )         THV  = SUMTHV2 * (MAPL_P00**MAPL_KAPPA) * (1.0/MAPL_GRAV)
              if( associated(DTHVDTFILINT)  .or. &
                  associated(DTHVDTPHYINT) ) then
                                                  SUMTHV1 = ( SUMTHV2-SUMTHV1 ) * (MAPL_P00**MAPL_KAPPA) / (DT*MAPL_GRAV)
              if( associated(DTHVDTFILINT) ) DTHVDTFILINT = SUMTHV1
              if( associated(DTHVDTPHYINT) ) DTHVDTPHYINT = DTHVDTPHYINT + SUMTHV1
              deallocate(SUMTHV1)
              endif
              deallocate(SUMTHV2)
              exit
           endif
        enddo
    endif


! Compute Tropopause Diagnostics
! ------------------------------
    if( associated(TROPP1) .or. &
        associated(TROPP2) .or. &
        associated(TROPP3) .or. &
        associated(TROPT)  .or. &
        associated(TROPQ)       ) then
        do K=1,NumFriendly
           NULLIFY(Q)
           call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), Q, RC=STATUS )
           VERIFY_(STATUS)
           if(NAMES(K)=='Q') then
              allocate( TROP(IM,JM,5),STAT=STATUS )
              VERIFY_(STATUS)
              call tropovars ( IM,JM,LM,PLE,PL,T,Q,EPV,TROP(:,:,1),TROP(:,:,2),TROP(:,:,3),TROP(:,:,4),TROP(:,:,5) )
               if( associated(TROPP1) )  TROPP1(:,:) = TROP(:,:,1)
               if( associated(TROPP2) )  TROPP2(:,:) = TROP(:,:,2)
               if( associated(TROPP3) )  TROPP3(:,:) = TROP(:,:,3)
               if( associated(TROPT ) )  TROPT (:,:) = TROP(:,:,4)
               if( associated(TROPQ ) )  TROPQ (:,:) = TROP(:,:,5)
                   deallocate(TROP )
               exit
           endif
        enddo
    endif

! Done
!-----

    deallocate(Names)
    deallocate(PL)
    deallocate(DP)
    deallocate(PKE)
    deallocate(PKZ)

    call MAPL_TimerOff(STATE,"RUN"  )
    call MAPL_TimerOff(STATE,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  contains

    subroutine Dry_Mass_Check (im,jm,lm,n)
    implicit none
    integer im,jm,lm
    integer L,n,status

    real,   allocatable :: dry_ana(:,:,:)
    real,   allocatable ::    qint(:,:)
    real*8, allocatable ::    qsum(:,:)
    real*8  dry_bkg_ave
    real*8  dry_ana_ave
    real*8  dry_diff

                       allocate(    qsum(IM,JM)     ,STAT=STATUS ) ; VERIFY_(STATUS)
                       allocate(    qint(IM,JM)     ,STAT=STATUS ) ; VERIFY_(STATUS)
                       allocate( dry_ana(IM,JM,  LM),STAT=STATUS ) ; VERIFY_(STATUS)

                       ple_ana = real(ple,kind=8) + real(dt,kind=8)*real(dpedt,kind=8)
                       dp_ana  = ple_ana(:,:,1:LM)-ple_ana(:,:,0:LM-1)

                       do L=1,lm
                          dry_ana(:,:,L) = ( 1.0 - (     q(:,:,L)+ qlls(:,:,L)+qlcn(:,:,L) + qils(:,:,L)+qicn(:,:,L)  &
                                                   + qrain(:,:,L)+qsnow(:,:,L)+qgraupel(:,:,L)) )*dp_ana(:,:,L)
                       enddo

                       qsum = 0.0_8
                       do L=1,lm
                       qsum = qsum + dry_bkg(:,:,L)
                       enddo
                       qint = qsum
                       call MAPL_AreaMean( dry_bkg_ave, qint, area, grid, rc=STATUS )
                       VERIFY_(STATUS)

                       qsum = 0.0_8
                       do L=1,lm
                       qsum = qsum + dry_ana(:,:,L)
                       enddo
                       qint = qsum
                       call MAPL_AreaMean( dry_ana_ave, qint, area, grid, rc=STATUS )
                       VERIFY_(STATUS)

                       dry_ana_ave = (dry_ana_ave+1.0)/100
                       dry_bkg_ave = (dry_bkg_ave+1.0)/100
                       dry_diff    =  dry_ana_ave - dry_bkg_ave
#if debug
                       if(MAPL_AM_I_ROOT() ) then
                          write(6,1000) n,dry_ana_ave,dry_bkg_ave,dry_diff
  1000                    format(5x,'n: ',i2,3x,'DRY_ANA: ',g,'  DRY_BKG: ',g,'  DIFF: ',g)
                       endif
#endif
                       deallocate(   qsum)
                       deallocate(   qint)
                       deallocate(dry_ana)

                       return

    end subroutine Dry_Mass_Check

    subroutine DO_UPDATE_ANA3D(NAME, COMP, PREF, CONSTRAIN_DAS)
      character*(*),     intent(IN) :: NAME
      integer,           intent(IN) :: COMP
      real,    pointer,  intent(IN) :: PREF(:)
      integer, optional, intent(IN) :: CONSTRAIN_DAS

      real,   pointer,     dimension(:,:)   :: DPSDT_CON => null()
      real,   pointer,     dimension(:,:,:) :: TENDSD    => null()
      real,   pointer,     dimension(:,:,:) :: TENDPH    => null()
      real,   pointer,     dimension(:,:,:) :: TENDBS    => null()
      real,   pointer,     dimension(:,:,:) :: TENDAN    => null()
      real,   pointer,     dimension(:,:,:) :: ANAINC    => null()
      real,   allocatable, dimension(:,:,:) :: TENDANAL
      real,   allocatable, dimension(:,:)   :: dummy
      real,   allocatable, dimension(:)     :: ALFZ, BETZ
      real*8                                :: qave1
      real*8                                :: qave2
      real*8                                :: qave3
      integer                               :: L,LL,LU

      call MAPL_GetPointer(GIM(COMP), TENDSD, trim(NAME)        , rc=STATUS)
      VERIFY_(STATUS)

      select case (TYPE)
         case (FREERUN)

            TENDSD = 0.0

      case (PREDICTOR)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)
         TENDSD = TENDBS

      case (FORECAST)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)

         TENDBS = BET * TENDBS
         TENDSD =       TENDBS

      case (CORRECTOR)

         LL = lbound(TENDSD,3)
         LU = ubound(TENDSD,3)

         allocate(TENDANAL(size(TENDSD,1),size(TENDSD,2),LL:LU), stat=STATUS)
         VERIFY_(STATUS)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME),  rc=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT   , ANAINC, trim(NAME),  rc=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT, DPSDT_CON, 'DPSDT_CON', rc=STATUS)
         VERIFY_(STATUS)

         TENDANAL = ANAINC*IAUcoeff  ! No Constraints

         IF( present(CONSTRAIN_DAS) ) then
                  IF(CONSTRAIN_DAS == 1 ) then
                  ! --------------------------
                     allocate(dummy(size(TENDSD,1),size(TENDSD,2)), stat=STATUS)
                     VERIFY_(STATUS)

                     ! Method to Simply Re-Scale
                     ! -------------------------
                     where( ANAINC(:,:,LU).ne.0.0 )
                            dummy = ANAINC(:,:,LU) ! ANAINC = PS_ANA - PS_BKG
                     elsewhere
                            dummy = MAPL_UNDEF
                     endwhere
                     call MAPL_AreaMean( qave1, dummy, area, grid, rc=STATUS )  ! qave1 = AreaMean( ANAINC )
                     VERIFY_(STATUS)

                     where( ANAINC(:,:,LU).ne.0.0 )
                            dummy = PLE(:,:,LU)    ! P_n
                     elsewhere
                            dummy = MAPL_UNDEF
                     endwhere
                     call MAPL_AreaMean( qave2, dummy, area, grid, rc=STATUS )  ! qave2 = AreaMean( P_n )

                     if(      qave1        .ne.MAPL_UNDEF  .and. &
                              qave2        .ne.MAPL_UNDEF ) then
                         qave3 = qave2 + qave1*dt*IAUcoeff ! qave3 = AreaMean( P_n+1 = P_n + ANAINC*dt/tau )
                     else
                         qave3 = MAPL_UNDEF
                     endif

                     if(      qave3        .ne.MAPL_UNDEF ) then
                         where( ANAINC(:,:,LU).ne.0.0 )
                                dummy = ANAINC(:,:,LU)*(qave2/qave3) - dummy*(qave1/qave3)   ! Preferred Method
                              ! dummy = ANAINC(:,:,LU)*real(qave2/qave3,kind=4) - dummy * real(qave1/qave3,kind=4)
                         elsewhere
                                dummy = ANAINC(:,:,LU)
                         endwhere
                     else
                                dummy = ANAINC(:,:,LU)
                     endif

                     if(associated(DPSDT_CON)) then
                                   DPSDT_CON = ( dummy - ANAINC(:,:,LU) )*IAUcoeff
                     endif

                     DO L=LL,LU
                     TENDANAL(:,:,L) = dummy(:,:)*BK(L)
                     ENDDO
                     TENDANAL = TENDANAL*IAUcoeff
                     deallocate(dummy)
                  ELSE
                     if(associated(DPSDT_CON)) then
                                   DPSDT_CON = 0.0
                     endif
                  ENDIF
         ENDIF

         TENDSD = (TENDBS + TENDANAL)

         if (LAST_CORRECTOR) then
            allocate( ALFZ(LL:LU), stat=STATUS)
            VERIFY_(STATUS)
            allocate( BETZ(LL:LU), stat=STATUS)
            VERIFY_(STATUS)

            DO L=LL,LU
            if( PREF(L).GT.10000.0 ) then
                ALFZ(L) = ALF      ! PREF > 100 mb
                BETZ(L) = BET      ! PREF > 100 mb
            elseif( PREF(L).LT.1000.0 ) then
                ALFZ(L) = 0.0      ! PREF <  10 mb
                BETZ(L) = 0.0      ! PREF <  10 mb
            else
                ALFZ(L) = ALF*( PREF(L)-1000.0 )/9000.0
                BETZ(L) = BET*( PREF(L)-1000.0 )/9000.0
            endif
            ENDDO

            DO L=LL,LU
            TENDBS(:,:,L) = BETZ(L)*TENDBS(:,:,L) + ALFZ(L)*TENDANAL(:,:,L)
            ENDDO

            deallocate(ALFZ)
            deallocate(BETZ)
         end if

         deallocate(TENDANAL)

      case default

         _ASSERT(.false.,'needs informative message')

      end select

! Fill Total Increment Tendency (Current + Bias) Diagnostic
! ---------------------------------------------------------

      call MAPL_GetPointer ( EXPORT, TENDAN, trim(NAME)//'_ANA', rc=STATUS )
      VERIFY_(STATUS)

      if(associated(TENDAN)) then
                    TENDAN = TENDSD
      endif

    end subroutine DO_UPDATE_ANA3D

    subroutine DO_UPDATE_ANA2D(NAME, COMP)
      character*(*), intent(IN) :: NAME
      integer,       intent(IN) :: COMP

      real, pointer,     dimension(:,:)     :: TENDSD   => null()
      real, pointer,     dimension(:,:)     :: TENDPH   => null()
      real, pointer,     dimension(:,:)     :: TENDBS   => null()
      real, pointer,     dimension(:,:)     :: TENDAN   => null()
      real, pointer,     dimension(:,:)     :: ANAINC   => null()
      real, allocatable, dimension(:,:)     :: TENDANAL

      call MAPL_GetPointer(GIM(COMP), TENDSD, trim(NAME)        , rc=STATUS)
      VERIFY_(STATUS)

      select case (TYPE)
         case (FREERUN)

            TENDSD = 0.0

      case (PREDICTOR)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)
         TENDSD = TENDBS

      case (FORECAST)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)

         TENDBS = BET * TENDBS
         TENDSD = TENDBS

      case (CORRECTOR)

         allocate( TENDANAL(size(TENDSD,1),size(TENDSD,2)), stat=STATUS )
         VERIFY_(STATUS)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT   , ANAINC, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)

         TENDANAL = ANAINC*IAUcoeff
         TENDSD   = (TENDBS + TENDANAL)

         if (LAST_CORRECTOR) then
            TENDBS = BET*TENDBS + ALF*TENDANAL
         end if

         deallocate(TENDANAL)

      case default

         _ASSERT(.false.,'needs informative message')

      end select

! Fill Total Increment Tendency (Current + Bias) Diagnostic
! ---------------------------------------------------------

      call MAPL_GetPointer ( EXPORT, TENDAN, trim(NAME)//'_ANA', rc=STATUS )
      VERIFY_(STATUS)

      if(associated(TENDAN)) then
                    TENDAN = TENDSD
      endif

    end subroutine DO_UPDATE_ANA2D


    subroutine DO_UPDATE_PHY (NAME)
      character*(*), intent(IN) :: NAME

      real, pointer,     dimension(:,:,:)     :: TENDSD   => null()
      real, pointer,     dimension(:,:,:)     :: TENDPH   => null()

      call MAPL_GetPointer(GIM(SDYN), TENDSD, trim(NAME)        , rc=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(GEX(PHYS), TENDPH, trim(NAME)        , rc=STATUS)
      VERIFY_(STATUS)

      TENDSD = TENDPH

    end subroutine DO_UPDATE_PHY


    subroutine DO_Friendly(Q, NAME, PREF)
      character*(*), intent(IN   ) :: NAME
      real, pointer, intent(IN)    :: PREF(:)
      real,          intent(INOUT) :: Q(:,:,:)

      real, pointer,     dimension(:,:,:)     :: TENDBS   => null()
      real, pointer,     dimension(:,:,:)     :: TENDAN   => null()
      real, pointer,     dimension(:,:,:)     :: ANAINC   => null()
      real, allocatable, dimension(:,:,:)     :: TENDANAL
      real, allocatable, dimension(:,:,:)     :: QOLD
      real, allocatable, dimension(:)         :: ALFZ, BETZ

      allocate( QOLD(IM,JM,LM), stat=STATUS)
      VERIFY_(STATUS)

      QOLD = Q   ! Initialize Old Value for Total Tendency Diagnostic

      select case (TYPE)
      case (PREDICTOR)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)

         Q = Q + max( DTX*TENDBS*FC, -Q )  ! Prevent Negative Q

      case (FORECAST)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)

         TENDBS = BET * TENDBS

         Q = Q + max( DTX*TENDBS*FC, -Q )  ! Prevent Negative Q

      case (CORRECTOR)

         allocate(TENDANAL(IM,JM,LM), stat=STATUS)
         VERIFY_(STATUS)

         call MAPL_GetPointer(INTERNAL , TENDBS, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT   , ANAINC, trim(NAME), rc=STATUS)
         VERIFY_(STATUS)

         TENDANAL = ANAINC*IAUcoeff

         Q = Q + max( DTX*(TENDBS + TENDANAL)*FC, -Q )  ! Prevent Negative Q

         if (LAST_CORRECTOR) then
            allocate( ALFZ(LM), stat=STATUS)
            VERIFY_(STATUS)
            allocate( BETZ(LM), stat=STATUS)
            VERIFY_(STATUS)

            DO L=1,LM
            if( PREF(L).GT.10000.0 ) then
                ALFZ(L) = ALF      ! PREF > 100 mb
                BETZ(L) = BET      ! PREF > 100 mb
            elseif( PREF(L).LT.1000.0 ) then
                ALFZ(L) = 0.0      ! PREF <  10 mb
                BETZ(L) = 0.0      ! PREF <  10 mb
            else
                ALFZ(L) = ALF*( PREF(L)-1000.0 )/9000.0
                BETZ(L) = BET*( PREF(L)-1000.0 )/9000.0
            endif
            ENDDO

            DO L=1,LM
            TENDBS(:,:,L) = BETZ(L)*TENDBS(:,:,L) + ALFZ(L)*TENDANAL(:,:,L)
            ENDDO

            deallocate(ALFZ)
            deallocate(BETZ)
         end if

         deallocate(TENDANAL)

      case default

         _ASSERT(.false.,'needs informative message')

      end select

! Fill Total Increment Tendency (Current + Bias) Diagnostic
! ---------------------------------------------------------
      call MAPL_GetPointer ( EXPORT, TENDAN, trim(NAME)//'_ANA', rc=STATUS )
      VERIFY_(STATUS)
      if(associated(TENDAN)) TENDAN = (Q-QOLD)/DT

      deallocate(QOLD)

    end subroutine DO_FRIENDLY

    subroutine FILL_Friendly ( Q,DP,QFILL,QINT )
      real, intent(INOUT) ::   Q(:,:,:)
      real, intent(IN   ) ::  DP(:,:,:)
      real, intent(OUT  ) :: QFILL(:,:)
      real, intent(OUT  ) :: QINT (:,:)

      real*8, allocatable, dimension(:,:) :: QTEMP1
      real*8, allocatable, dimension(:,:) :: QTEMP2

      allocate(QTEMP1(IM,JM), stat=STATUS)
      VERIFY_(STATUS)
      allocate(QTEMP2(IM,JM), stat=STATUS)
      VERIFY_(STATUS)

      QTEMP1 = 0.0
      do L=1,LM
      QTEMP1(:,:) = QTEMP1(:,:) + Q(:,:,L)*DP(:,:,L)
      enddo

      where( Q < 0.0 ) Q = 0.0

      QTEMP2 = 0.0
      do L=1,LM
      QTEMP2(:,:) = QTEMP2(:,:) + Q(:,:,L)*DP(:,:,L)
      enddo

      where( qtemp2.ne.0.0_8 )
             qtemp2 = max( qtemp1/qtemp2, 0.0_8 )
      end where

      do L=1,LM
      Q(:,:,L) = Q(:,:,L)*qtemp2(:,:)
      enddo

      QTEMP2 = 0.0
      do L=1,LM
      QTEMP2(:,:) = QTEMP2(:,:) + Q(:,:,L)*DP(:,:,L)
      enddo

      WHERE( QTEMP1 >= 0.0 )
              QFILL  = 0.0
      ELSEWHERE
              QFILL = -QTEMP1 / (DT*MAPL_GRAV)
      END WHERE
      QINT  =  QTEMP2         /     MAPL_GRAV

      deallocate(QTEMP1)
      deallocate(QTEMP2)

    end subroutine FILL_FRIENDLY

    subroutine get_iau_coeff( TNDCoeff,CLOCK )
    implicit none

    real, intent(OUT) :: TNDCoeff
    type(ESMF_Clock), intent(inout) :: CLOCK
    type(ESMF_Time)                 :: currtime
    type(ESMF_Time)                 :: MKIAU_RefTime
    type(ESMF_Calendar)             :: cal
    type(ESMF_TimeInterval)         :: MKIAU_HALF_FREQUENCY
    type(ESMF_TimeInterval)         :: TIME_Offset

    real*8  :: TIME_Fraction
    integer :: nsteps
    integer :: kstep, kshift
    integer :: MKIAU_RingDate
    integer :: MKIAU_RingTime
    integer :: rep_YY, rep_MM, rep_DD
    integer :: rep_H,  rep_M,  rep_S
    integer :: MKIAU_FREQUENCY
    logical :: IAU_DIGITAL_FILTER

    character(len=ESMF_MAXSTR) :: STRING
    character(len=ESMF_MAXSTR) :: REPLAY_MODE

    type (CONNECT_IAUcoeffs), pointer :: myCoeffs => NULL()
    type (IAU_coeffs)                 :: wrap
    real, allocatable                 :: shifted_dfi(:)

    call ESMF_ClockGet( CLOCK, currTime=currTime, calendar=cal, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(STATE, REPLAY_MODE, Label='REPLAY_MODE:', default="NoReplay", RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetResource(STATE, STRING, LABEL="IAU_DIGITAL_FILTER:", default="YES", RC=STATUS)
    VERIFY_(STATUS)
    STRING = ESMF_UtilStringUpperCase(STRING)
    IAU_DIGITAL_FILTER = trim(STRING)=="YES"

!   Standard Constant IAU Scaling (1/TAU)
!   -------------------------------------
    if ( .not. IAU_DIGITAL_FILTER ) then
          TNDcoeff = 1.0/TAUANL
          return
    else

!   Digital Filter IAU Scaling
!   --------------------------
       call ESMF_UserCompGetInternalState(GC, 'IAU_COEFFS', wrap, status)
       VERIFY_(STATUS)
       myCoeffs => wrap%ptr

     ! Note: MKIAU_FREQUENCY is Initialized in GCM_GridComp
     ! ----------------------------------------------------------
       call MAPL_GetResource( STATE,MKIAU_FREQUENCY, Label="MKIAU_FREQUENCY:", rc=STATUS )
       VERIFY_(STATUS)

       nsteps = nint( MKIAU_FREQUENCY/DT ) + 1

       if (.not.associated(myCoeffs%dfi)) then
            allocate(myCoeffs%dfi(nsteps))
            myCoeffs%istep=0

            call dfi_coeffs (DT,MKIAU_FREQUENCY,TAUANL,nsteps,myCoeffs%dfi)

        ! Shift DFI Coefficients if Necessary
        ! -----------------------------------
           call MAPL_GetResource( STATE, MKIAU_RingDate, Label="MKIAU_RingDate:", RC=STATUS )
           VERIFY_(STATUS)
           call MAPL_GetResource( STATE, MKIAU_RingTime, Label="MKIAU_RingTime:", RC=STATUS )
           VERIFY_(STATUS)

           call ESMF_TimeIntervalSet( MKIAU_HALF_FREQUENCY, S=MKIAU_FREQUENCY/2, rc=STATUS )
           VERIFY_(STATUS)

         ! REPACK MKIAU_RingDate and MKIAU_RingTime
         ! ----------------------------------------
           rep_YY =     MKIAU_RingDate /10000
           rep_MM = mod(MKIAU_RingDate ,10000)/100
           rep_DD = mod(MKIAU_RingDate ,100)
           rep_H  =     MKIAU_RingTime /10000
           rep_M  = mod(MKIAU_RingTime ,10000)/100
           rep_S  = mod(MKIAU_RingTime ,100)

           call ESMF_TimeSet( MKIAU_RefTime, YY = rep_YY, &
                                             MM = rep_MM, &
                                             DD = rep_DD, &
                                              H = rep_H , &
                                              M = rep_M , &
                                              S = rep_S , &
                              calendar=cal,  rc = STATUS  )
           VERIFY_(STATUS)

           TIME_Offset   = MKIAU_RefTime - CurrTime
           TIME_Fraction = Time_Offset   / MKIAU_HALF_FREQUENCY

           kshift = abs( 1.0 - Time_Fraction )*(nsteps-1)/2

           allocate( shifted_dfi(nsteps) )
           do i=1,nsteps-1
                                     kstep =     i +  kshift
             if( kstep.gt.nsteps-1 ) kstep = kstep - (nsteps-1)
             shifted_dfi(i) = myCoeffs%dfi(kstep)
           enddo
           shifted_dfi(nsteps) = shifted_dfi(1)
           myCoeffs%dfi        = shifted_dfi
           deallocate( shifted_dfi )

           if (MAPL_am_I_root()) then
              print*, 'DFI initialized for',nsteps,' steps'
              do i=1,nsteps-1
              print '(1x,a,i4,5x,a,g13.6,5x,a,i6)', 'i: ',i,'DFI Coeff: ',TAUANL*myCoeffs%dfi(i)
              enddo
           endif
       endif

                                myCoeffs%istep = myCoeffs%istep + 1
       TNDcoeff = myCoeffs%dfi( myCoeffs%istep )

     ! if( MAPL_am_I_root() ) print '(1x,a,i4,5x,a,g13.6)', 'i: ',myCoeffs%istep,'DFI Coeff: ',TAUANL*TNDcoeff
       if( myCoeffs%istep == nsteps-1 ) myCoeffs%istep = 0

    endif

    end subroutine get_iau_coeff

    subroutine update_ainc_(RC)

    use ESMF_CFIOFileMod
    use GEOS_UtilsMod
    use GEOS_RemapMod, only: myremap => remap
    implicit none

    integer,optional, intent(OUT) :: RC

! Brief description: This routine converts either analysis increments or
! analysis fields to model increments corresponding to what IAU
! tendencies require.
! In case of handling analysis full fields:
!    - convert background fields corresponding to analysis fields to
!      analysis resolution.
!    - calculate (IAU) increment on analysis grid
!    - if needed, interpolate increments to model grid
! In case of handling analysis increment:
!    - if needed, interpolate analysis increment to model grid
!    - convert analysis increment to IAU increment (e.g., tv to td)


! The following name declarations will be easily unwired ...
    character(len=*), parameter :: incnames(7) = (/ 'sphu ', &
                                                    'u    ', &
                                                    'v    ', &
                                                    'tv   ', &
                                                    'ozone', &
                                                    'ps   ', &
                                                    'ts   ' /)

    character(len=*), parameter :: ananames(9) = (/ 'sphu ', &
                                                    'u    ', &
                                                    'v    ', &
                                                    'tv   ', &
                                                    'ozone', &
                                                    'delp ', &
                                                    'phis ', &
                                                    'ps   ', &
                                                    'ts   ' /)
    integer rank,ni,nt,nvars,natts
    integer IMana_World,JMana_World
    integer IMbkg_World,JMbkg_World
    integer IMana,JMana,LMana
    integer IMbkg,JMbkg,LMbkg
    integer NX,NY
    integer nymd,nhms
    integer vm_comm
    integer :: DIMS(ESMF_MAXGRIDDIM)
    integer :: K,NQ,FID
    integer :: idum
    integer :: method
    integer :: AINC_TIME(6)
    integer :: NPHIS, NPHIS_MAX
    logical do_transforms
    logical l_cube
    logical l_nudge
    logical l_remap
    logical l_store_transforms
    logical l_use_ana_delp
    logical fromANA2BKG
    logical IuseTS
    logical l_windfix
    logical done_remap
    type(ESMF_Config)                   :: CF
    type(ESMF_Field)                    :: Field
    type(ESMF_FieldBundle)              :: RBUNDLE
    type(ESMF_Time)                     :: AincTime
    type(ESMF_VM)                       :: VM
    real,allocatable, dimension(:,:)    :: dps
    real,allocatable, dimension(:,:,:)  :: dp_inc
    real,allocatable, dimension(:,:,:)  :: uptr
    real,allocatable, dimension(:,:,:)  :: vptr
    real sclinc
!
    real, pointer, dimension(:,:,:) :: qdum1
    real, pointer, dimension(:,:,:) :: qdum2
!
    real, pointer, dimension(:,:,:)     ::   u_bkg => NULL()
    real, pointer, dimension(:,:,:)     ::   u_ana => NULL()
    real, pointer, dimension(:,:,:)     ::  du_inc => NULL()
    real, pointer, dimension(:,:,:)     ::   v_bkg => NULL()
    real, pointer, dimension(:,:,:)     ::   v_ana => NULL()
    real, pointer, dimension(:,:,:)     ::  dv_inc => NULL()
    real, pointer, dimension(:,:,:)     ::  tv_bkg => NULL()
    real, pointer, dimension(:,:,:)     ::  tv_ana => NULL()
    real, pointer, dimension(:,:,:)     ::  dt_inc => NULL()
    real, pointer, dimension(:,:,:)     ::   q_bkg => NULL()
    real, pointer, dimension(:,:,:)     ::   q_ana => NULL()
    real, pointer, dimension(:,:,:)     ::  dq_inc => NULL()
    real, pointer, dimension(:,:,:)     ::  o3_bkg => NULL()
    real, pointer, dimension(:,:,:)     ::  o3_ana => NULL()
    real, pointer, dimension(:,:,:)     :: do3_inc => NULL()
    real, pointer, dimension(:,:,:)     :: ple_bkg => NULL()
    real, pointer, dimension(:,:,:)     ::delp_bkg => NULL()
    real, pointer, dimension(:,:,:)     :: ple_ana => NULL()
    real, pointer, dimension(:,:,:)     ::dple_inc => NULL()
    real, pointer, dimension(:,:,:)     ::  pk_ana => NULL()
    real, pointer, dimension(:,:,:)     :: pke_ana => NULL()
    real, pointer, dimension(:,:,:)     ::  pk_bkg => NULL()
    real, pointer, dimension(:,:,:)     :: pke_bkg => NULL()
    real, pointer, dimension(:,:,:)     ::     dpk => NULL()
    real, pointer, dimension(:,:,:)     :: thv_ana => NULL()
    real, pointer, dimension(:,:)       ::  ts_bkg => NULL()
    real, pointer, dimension(:,:)       ::  ts_ana => NULL()
    real, pointer, dimension(:,:)       :: dts_inc => NULL()
    real, pointer, dimension(:,:)       ::  ps_bkg => NULL()
    real, pointer, dimension(:,:)       ::  ps_ana => NULL()
    real, pointer, dimension(:,:)       ::phis_bkg
    real, pointer, dimension(:,:)       ::phis_ana => NULL()
!
    real, pointer, dimension(:,:,:)     :: dpkz => NULL()
!
    real, pointer, dimension(:,:)   ::  dts_aux => NULL()
    real, pointer, dimension(:,:)   ::  dps_aux => NULL()
    real, pointer, dimension(:,:,:) ::   du_aux => NULL()
    real, pointer, dimension(:,:,:) ::   dv_aux => NULL()
    real, pointer, dimension(:,:,:) ::   dt_aux => NULL()
    real, pointer, dimension(:,:,:) ::   dp_aux => NULL()
    real, pointer, dimension(:,:,:) ::   dq_aux => NULL()
    real, pointer, dimension(:,:,:) ::  do3_aux => NULL()
    real, pointer, dimension(:,:,:) :: dple_aux => NULL()
!
    real,pointer :: aptr2d(:,:)  => NULL()  ! analysis increment pointers
    real,pointer :: aptr3d(:,:,:)=> NULL()  ! analysis increment pointers
    real,pointer :: gptr2d(:,:)  => NULL()  ! gcm background pointers
    real,pointer :: gptr3d(:,:,:)=> NULL()  ! gcm background pointers
!
    real, allocatable, dimension(:,:)   ::  vintdiva
    real, allocatable, dimension(:,:)   ::  vintdivb
    real, allocatable, dimension(:,:)   ::  vintdivc
!
    character(len=ESMF_MAXSTR) :: name
    character(len=ESMF_MAXSTR) :: imstr, jmstr
    character(len=ESMF_MAXSTR) :: FILETMPL,AINCFILE
    character(len=ESMF_MAXSTR) :: STRING
    character(len=ESMF_MAXSTR) :: DATE
    character(len=ESMF_MAXSTR) :: NUDGE
    character(len=ESMF_MAXSTR) :: NUDGE_REMAP
    character(len=ESMF_MAXSTR) :: NUDGE_WINDFIX
    character(len=ESMF_MAXSTR) :: NUDGE_STORE_TRANSFORMS
    character(len=ESMF_MAXSTR) :: USE_ANA_DELP
    character(len=ESMF_MAXSTR) :: SKIP_FIRST

    type(ESMF_TimeInterval)    :: Frequency

    class (AbstractRegridder), pointer :: ANA2BKG => null()
    class (AbstractRegridder), pointer :: BKG2ANA => null()

    type (CONNECT_ANAnBKG), pointer :: myANA => NULL()
    type (ANAnBKG)                  :: wrap

    type(ESMF_Grid), EXTERNAL :: AppGridCreateF

    integer                             :: ifirst
    data                                   ifirst /0/

    INTEGER DATE_TIME(8)
    CHARACTER (LEN = 12) REAL_CLOCK(3)

!   When assimilation period is over, do not even bother ...
!   --------------------------------------------------------
    call MAPL_GetResource(STATE, FILETMPL,    LABEL="AINC_FILE:", default="NULL", RC=STATUS)
    VERIFY_(STATUS)

!   If analysis (or increment) file not specified, move on ...
!   ----------------------------------------------------------
    if( trim(FILETMPL)=="NULL") return

!   Get pointers to analysis tendencies
!   -----------------------------------
    call MAPL_GetPointer(import,   du_inc, 'DUDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,   dv_inc, 'DVDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,   dt_inc, 'DTDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,   dq_inc, 'DQVDT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,  do3_inc, 'DO3DT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import, dple_inc, 'DPEDT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,  dts_inc, 'DTSDT', RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource( STATE, sclinc, label ='SCLINC:', default=-1.0, rc=status )

    call MAPL_GetResource(STATE, SKIP_FIRST,  LABEL="SKIP_FIRST:", default="NO", RC=STATUS)
    VERIFY_(STATUS)
    if(trim(SKIP_FIRST)=="YES") then  ! when trying to "replay" to trajectory of free model
                                      ! either have agcm_import filled w/ zeros
                                      ! or bypass zero tendencies out in first
                                      ! pass
      if(ifirst<1)then
        du_inc  =0.0
        dv_inc  =0.0
        dt_inc  =0.0
        dq_inc  =0.0
        do3_inc =0.0
        dple_inc=0.0
        dts_inc =0.0
        ifirst =ifirst+1
        return
      endif

    else                              ! typicall, will expect tendency at initial time
                                      ! to be meaningful so, use what is in agcm_import
                                      ! in the first pass, from then on analysis should
                                      ! be read in at desired frequency
      if(ifirst<1)then
        ifirst =ifirst+1
        if (ifirst==1.and.sclinc>0.0) then  ! if neeed, only once, IAU tendency from RST
            du_inc   = sclinc * du_inc
            dv_inc   = sclinc * dv_inc
            dq_inc   = sclinc * dq_inc
            dple_inc = sclinc * dple_inc
            dt_inc   = sclinc * dt_inc
            do3_inc  = sclinc * do3_inc
            dts_inc  = sclinc * dts_inc
            if(MAPL_AM_I_ROOT()) then
               print *
               print *, 'Scaled initial IAU tendency by: ', sclinc
            endif
        endif
        return
      endif

    endif

!  Get my internal private state. This contains the transforms
!  between the background grid and the ANA grid, as well as ANA grid
!  itself.
!-------------------------------------------------------------
    call ESMF_UserCompGetInternalState(GC, 'UPD_STATE', wrap, status)
    VERIFY_(STATUS)
    myANA => wrap%ptr

    call MAPL_GetResource(STATE, NUDGE,    LABEL="NUDGE_STATE:", default="NO", RC=STATUS)
    VERIFY_(STATUS)
    NUDGE = ESMF_UtilStringUpperCase(NUDGE)
    l_nudge=trim(NUDGE)=="YES"

    call MAPL_GetResource(STATE, NUDGE_REMAP, LABEL="NUDGE_REMAP:", default="YES", RC=STATUS)
    VERIFY_(STATUS)
    NUDGE_REMAP = ESMF_UtilStringUpperCase(NUDGE_REMAP)
    l_remap=trim(NUDGE_REMAP)=="YES"

    call MAPL_GetResource(STATE, NUDGE_WINDFIX, LABEL="NUDGE_WINDFIX:", default="YES", RC=STATUS)
    VERIFY_(STATUS)
    NUDGE_WINDFIX = ESMF_UtilStringUpperCase(NUDGE_WINDFIX)
    l_windfix=trim(NUDGE_WINDFIX)=="YES"

    call MAPL_GetResource(STATE, USE_ANA_DELP, LABEL="USE_ANA_DELP:", default="NO", RC=STATUS)
    VERIFY_(STATUS)
    USE_ANA_DELP = ESMF_UtilStringUpperCase(USE_ANA_DELP)
    l_use_ana_delp=trim(USE_ANA_DELP)=="YES"

    call MAPL_GetResource(STATE, NUDGE_STORE_TRANSFORMS, LABEL="NUDGE_STORE_TRANSFORMS:", default="YES", RC=STATUS)
    VERIFY_(STATUS)
    NUDGE_STORE_TRANSFORMS = ESMF_UtilStringUpperCase(NUDGE_STORE_TRANSFORMS)
    l_store_transforms=trim(NUDGE_STORE_TRANSFORMS)=="YES"

    IMbkg=IM
    JMbkg=JM
    LMbkg=LM

!   Validate grid
!   -------------
    call ESMF_GridValidate(GRID,RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridGet(GRID, globalCellCountPerDim=DIMS, RC=STATUS)
    VERIFY_(STATUS)
    IMbkg_World=DIMS(1)
    JMbkg_World=DIMS(2)
    L_CUBE = JMbkg_World==6*IMbkg_World

    call ESMF_ClockGet(clock, currTime=currTime, rc=status)
    VERIFY_(STATUS)
    call ESMF_TimeGet(currTIME, timeString=DATE, RC=STATUS)
    VERIFY_(STATUS)
    call strToInt(DATE, nymd, nhms)

    call ESMF_CFIOstrTemplate ( AINCFILE, FILETMPL, 'GRADS', nymd=nymd, nhms=nhms, stat=STATUS )
    VERIFY_(STATUS)
    if(MAPL_AM_I_ROOT()) then
       print *
       print *, 'Overwriting IAU-inc with ANA-Inc: ', trim(AINCFILE)
    endif

    call MAPL_GetResource( STATE, NX,  Label="NX:", RC=status )
    VERIFY_(STATUS)
    call MAPL_GetResource( STATE, NY,  Label="NY:", RC=status )
    VERIFY_(STATUS)
    call MAPL_GetResource( STATE, idum, 'ANALYZE_TS:', default=0, RC=STATUS )
    VERIFY_(STATUS)
    IuseTS=idum/=0

    if (.not.myANA%initialized) then

       call CFIO_Open       ( AINCFILE, 1, fid, STATUS )
       VERIFY_(STATUS)
       call CFIO_DimInquire ( fid, IMana_World, JMana_world, LM, nt, nvars, natts, rc=STATUS )
       VERIFY_(STATUS)
       call CFIO_Close      ( fid, STATUS )
       VERIFY_(STATUS)

       call make_ana_grid(myANA, IMana_World, JMana_World, NX, NY, LM, rc)
!      Validate grid
!      -------------
       call ESMF_GridValidate(myAna%GRIDana,RC=STATUS)
       VERIFY_(STATUS)

       call MAPL_GridGet(myANA%GRIDana, localCellCountPerDim=DIMS, RC=STATUS)
       VERIFY_(STATUS)
       myANA%IM = DIMS(1)
       myANA%JM = DIMS(2)
       myANA%LM = DIMS(3)

!      Create transforms handles
!      -------------------------
       myANA%do_transforms = ( IMbkg_World /= IMana_World ) .or. &
                             ( JMbkg_World /= JMana_World )

      if (myANA%do_transforms) then
         if(MAPL_AM_I_ROOT()) then
            print *
            print *, 'Ana res: ', IMana_World,JMana_World,' Bkg res: ',IMbkg_World,JMbkg_World
            print *
         endif
         myANA%ANA2BKG_regridder => regridder_manager%make_regridder(myANA%GRIDana, GRID, REGRID_METHOD_BILINEAR, rc=status)
         VERIFY_(status)

         myANA%BKG2ANA_regridder => regridder_manager%make_regridder(GRID, myANA%GRIDana, REGRID_METHOD_BILINEAR, rc=status)
          VERIFY_(status)
         call WRITE_PARALLEL("Created transforms Ana2Bkg/Bkg2Ana ...")
      endif

!     When remap required, need to get PHIS from model state onto analysis grid
!     -------------------------------------------------------------------------
      if (l_remap) then

!        Get PHIS from background
!        ------------------------
         call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
         VERIFY_(STATUS)
         call ESMF_StateGet( GIM(SDYN), 'PHIS', FIELD, rc=STATUS )
         VERIFY_(STATUS)
         Call GEOS_TopoGet ( CF, MEAN=FIELD, rc=STATUS )
         VERIFY_(STATUS)

         call ESMF_StateGet( GIM(PHYS), 'PHIS', FIELD, rc=STATUS )
         VERIFY_(STATUS)
         Call GEOS_TopoGet ( CF, MEAN=FIELD, rc=STATUS )
         VERIFY_(STATUS)
         call ESMF_FieldGet (FIELD, localDE=0, farrayPtr=phis_bkg, rc = status)

!        Derive PHIS as consequence of going from BKG to ANA back to BKG grid
!        --------------------------------------------------------------------
         allocate(myANA%phis_bkg(myANA%IM,myANA%JM),stat=STATUS);VERIFY_(STATUS)
         if (myANA%do_transforms) then
             allocate(qdum1(myANA%IM,myANA%JM,1))
             allocate(qdum2(IMbkg,JMbkg,1))
             qdum2(:,:,1)=phis_bkg
             call myAna%bkg2ana_regridder%regrid(qdum2, qdum1, rc=status); VERIFY_(STATUS)
             myANA%phis_bkg=qdum1(:,:,1)
             deallocate(qdum2)
             deallocate(qdum1)
         else
             myANA%phis_bkg=phis_bkg
         endif
      endif

      myAna%Initialized=.true.
    endif

    IMana=myANA%IM
    JMana=myANA%JM
    LMana=myANA%LM
    do_transforms = myANA%do_transforms
    if (do_transforms) then
        ANA2BKG => myANA%ANA2BKG_regridder
        BKG2ANA => myANA%BKG2ANA_regridder
    endif

    AINC_TIME(1) =     nymd/10000
    AINC_TIME(2) = mod(nymd,10000)/100
    AINC_TIME(3) = mod(nymd,100)
    AINC_TIME(4) =     nhms/10000
    AINC_TIME(5) = mod(nhms,10000)/100
    AINC_TIME(6) = mod(nhms,100)

!   Set analysis increment time
!   ---------------------------
    call ESMF_TimeSet(  AincTime, YY =  AINC_TIME(1), &
                                  MM =  AINC_TIME(2), &
                                  DD =  AINC_TIME(3), &
                                  H  =  AINC_TIME(4), &
                                  M  =  AINC_TIME(5), &
                                  S  =  AINC_TIME(6), rc=status ); VERIFY_(STATUS)

!   Get MPI communicator
!   --------------------
    call ESMF_VMGetCurrent(vm, rc=status)
    VERIFY_(STATUS)
    call ESMF_VmGet(VM, mpicommunicator=vm_comm, rc=status)
    VERIFY_(STATUS)

! *****************************************************************************
! ****   READ Internal STATE (ie. ANA.ETA) from REPLAY File into BUNDLE    ****
! *****************************************************************************

    RBundle = ESMF_FieldBundleCreate( RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldBundleSet(RBundle, grid=myAna%GRIDana, rc=status)
    VERIFY_(STATUS)
    call MAPL_CFIORead ( AINCFILE, AincTime, Rbundle , RC=status)
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet ( RBUNDLE, fieldCount=NQ, RC=STATUS )
    VERIFY_(STATUS)

!   Get pointer to background fields (current GCM state)
!   ----------------------------------------------------
    call MAPL_GetPointer( GEX(SDYN),        AK,'AK', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN),        BK,'BK', RC=STATUS)
    VERIFY_(STATUS)
    if (l_use_ana_delp) then
       call MAPL_GetPointer( GEX(SDYN), delp_bkg,'DELP',RC=STATUS)
       VERIFY_(STATUS)
    else
       call MAPL_GetPointer( GEX(SDYN),  ple_bkg,'PLE', RC=STATUS)
       VERIFY_(STATUS)
    endif
    call MAPL_GetPointer( GEX(SDYN),      u_bkg,'U', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN),      v_bkg,'V', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN),    tv_bkg,'TV', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(PHYS),o3_bkg,'O3PPMV', RC=STATUS)
    VERIFY_(STATUS)
    do K=1,NumFriendly
       NULLIFY(q_bkg)
       call ESMFL_BundleGetPointerToData(BUNDLE, Names(K), q_bkg, RC=STATUS )
       VERIFY_(STATUS)
       if(Names(K)=='Q') exit
    enddo
    call MAPL_GetPointer( GEX(PHYS), ts_bkg, 'TS', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( GEX(SDYN), ps_bkg, 'PS', RC=STATUS)
    VERIFY_(STATUS)

!   Loop over GSI increment fields
!   ------------------------------
    if (l_nudge) then ! Bundle has analysis, leave it on its own grid
       allocate(gptr3d(IMana,JMana,LMana), stat=status );VERIFY_(STATUS)
       allocate(gptr2d(IMana,JMana),       stat=status );VERIFY_(STATUS)
       allocate(dp_aux  (IMana,JMana,LMana),stat=status );VERIFY_(STATUS)
       allocate(dts_aux (IMana,JMana),      stat=status );VERIFY_(STATUS)
       allocate(du_aux  (IMana,JMana,LMana),stat=status );VERIFY_(STATUS)
       allocate(dv_aux  (IMana,JMana,LMana),stat=status );VERIFY_(STATUS)
       allocate(dt_aux  (IMana,JMana,LMana),stat=status );VERIFY_(STATUS)
       allocate(dq_aux  (IMana,JMana,LMana),stat=status );VERIFY_(STATUS)
       allocate(do3_aux (IMana,JMana,LMana),stat=status );VERIFY_(STATUS)
       allocate(dple_aux(IMana,JMana,0:LMana),stat=status );VERIFY_(STATUS)
       if(l_use_ana_delp) then
          allocate(ple_bkg (IMbkg,JMbkg,0:LMbkg),stat=status );VERIFY_(STATUS)
       else
          allocate(dps_aux (IMana,JMana),        stat=status );VERIFY_(STATUS)
       endif
    else              ! Bundle has increment, therefore bring it to GCM grid
       allocate(gptr3d(IMbkg,JMbkg,LMbkg),stat=STATUS); VERIFY_(STATUS)
       allocate(gptr2d(IMbkg,JMbkg),  stat=STATUS); VERIFY_(STATUS)
       allocate(qdum1 (IMbkg,JMbkg,1),stat=STATUS); VERIFY_(STATUS)
       allocate(qdum2 (IMana,JMana,1),stat=STATUS); VERIFY_(STATUS)
    endif
    allocate(phis_ana(IMana,JMana),stat=STATUS)
    VERIFY_(STATUS)
    allocate(dp_inc(size(gptr3d,1),size(gptr3d,2),size(gptr3d,3)),stat=STATUS)
    VERIFY_(STATUS)
    allocate(dps(size(gptr2d,1),size(gptr2d,2)),stat=STATUS)
    VERIFY_(STATUS)
    fromANA2BKG=(.not.l_nudge) .and. do_transforms
    do ni = 1, nq
       call ESMF_FieldBundleGet(RBUNDLE, ni, Field, __RC__ )
       call ESMF_FieldGet(Field, NAME=NAME, dimCount = rank, __RC__ )
       if (l_nudge) then
          if (.not.check_list_(NAME,ananames)) cycle
       else
          if (.not.check_list_(NAME,incnames)) cycle
       endif
       if (rank==2) then
           call ESMF_FieldGet(Field, farrayPtr=aptr2d, __RC__ )
           if (fromANA2BKG) then
               qdum2(:,:,1)=aptr2d
               call ana2bkg%regrid(qdum2, qdum1, rc=status); VERIFY_(STATUS)
               gptr2d=qdum1(:,:,1)
           else
               gptr2d=aptr2d
           endif
           if (fromANA2BKG) then
              if(trim(NAME)=='ts') dts_inc = gptr2d
           else
              if(trim(NAME)=='ts') dts_aux = gptr2d
           endif
           if(trim(NAME)=='ps') then
              dps=gptr2d
           endif
           if(trim(NAME)=='phis') then
              phis_ana=gptr2d
           endif
       else
           call ESMF_FieldGet(Field, farrayPtr=aptr3d, __RC__ )
           if(trim(NAME)=='u'.or.trim(NAME)=='v') then
              if(trim(NAME)=='u') then
                 allocate(uptr(IMana,JMana,LMana))
                 uptr=aptr3d
                 cycle
              endif
              if(trim(NAME)=='v') then
                 allocate(vptr(IMana,JMana,LMana))
                 vptr=aptr3d
                 cycle
              endif
           endif
           if (fromANA2BKG) then
             call ANA2BKG%regrid(aptr3d, gptr3d, rc=status); VERIFY_(STATUS)
           else
               gptr3d=aptr3d
           endif
           if (fromANA2BKG) then
              if(trim(NAME)=='ozone') do3_inc=gptr3d
              if(trim(NAME)=='sphu' ) dq_inc =gptr3d
              if(trim(NAME)=='tv'   ) dt_inc =gptr3d
              if(trim(NAME)=='delp'.and.l_use_ana_delp ) dp_inc =gptr3d
           else
              if(trim(NAME)=='ozone') do3_aux=gptr3d
              if(trim(NAME)=='sphu' ) dq_aux =gptr3d
              if(trim(NAME)=='tv'   ) dt_aux =gptr3d
              if(trim(NAME)=='delp'.and.l_use_ana_delp ) dp_aux =gptr3d
           endif
           if(.not.l_use_ana_delp) then
              if(trim(NAME)=='delp') then
                 dp_inc =gptr3d
              endif
           endif

       endif
    enddo
    if (fromANA2BKG) then
       deallocate(qdum2)
       deallocate(qdum1)
    endif

!   U and V
!   -------
    ! could apply wind fix here ... but conversion to cubed
    ! TBD
    if (fromANA2BKG) then
       if( L_CUBE ) then
          call ANA2BKG%regrid(uptr, vptr, du_aux, dv_aux, rc=status); VERIFY_(STATUS)
       else
          call ANA2BKG%regrid(uptr, du_aux, rc=status); VERIFY_(STATUS)
          call ANA2BKG%regrid(vptr, dv_aux, rc=status); VERIFY_(STATUS)
!          call POLEFIX ( du_aux,dv_aux,VM,GRID )
       endif
    else
        du_aux=uptr
        dv_aux=vptr
    endif

!   Calculate 3d-pressure change
!   -----------------------------
    if (l_use_ana_delp) then
       if (.not.l_nudge) then ! increment
           dple_inc(:,:,0)=0.0
           do L=1,LM
              dple_inc(:,:,L)=dp_inc(:,:,L)
           enddo
       endif
    else
       if (l_nudge) then ! full fields
           dps_aux = dps
       else              ! increments
           do L=0,LMbkg
              dple_inc(:,:,L)=bk(L)*dps(:,:)
           enddo
       endif
    endif

!   Convert virtual temperature increment into dry temperature increment
!   -------------------------------------------------------------------
    EPS = MAPL_RVAP/MAPL_RGAS-1.0
    if (.not.l_nudge) then ! in this case, using the background fields is not
                           ! quite legitimate since these refer to the really
                           ! current trajectory and not quite the original
                           ! background used by the analysis
       dt_inc = dt_inc /(1.0+eps*q_bkg) - eps*dq_inc*tv_bkg/((1.0+eps*q_bkg)*(1.0+eps*q_bkg)) ! now dt is inc in dry temperature
    endif

!   When nugding, RBundle carries full analysis state.
!   In this case the follow will take place:
!     1. Convert present state of GCM to analysis grid
!     2. Calculate difference of analysis and present GCM state on analysis grid
!     3. Convert difference field back to GCM grid and overwrite analysis "tendencies"
!   -----------------------------------------------------------------------------------
    if (l_nudge) then

       allocate(u_ana  (IMana,JMana,  LMana),stat=STATUS); VERIFY_(STATUS)
       allocate(v_ana  (IMana,JMana,  LMana),stat=STATUS); VERIFY_(STATUS)
       allocate(tv_ana (IMana,JMana,  LMana),stat=STATUS); VERIFY_(STATUS)
       allocate(q_ana  (IMana,JMana,  LMana),stat=STATUS); VERIFY_(STATUS)
       allocate(o3_ana (IMana,JMana,  LMana),stat=STATUS); VERIFY_(STATUS)
       allocate(ple_ana(IMana,JMana,0:LMana),stat=STATUS); VERIFY_(STATUS)

       u_ana = du_aux
       v_ana = dv_aux
       q_ana = dq_aux
       tv_ana= dt_aux
       o3_ana= do3_aux

       if (l_use_ana_delp) then
!         Analyzed pressure edges
!         -----------------------
          ple_ana(:,:,0) = ak(0)
          do L=1,LM
             ple_ana(:,:,L) = ple_ana(:,:,L-1) + dp_aux(:,:,L)
          enddo

!         Background pressure edges
!         -------------------------
          ple_bkg(:,:,0) = ak(0)
          do L=1,LM
             ple_bkg(:,:,L) = ple_bkg(:,:,L-1) + delp_bkg(:,:,L)
          enddo
       else
!         Analyzed surface pressure
!         -------------------------
          allocate( ps_ana(IMana,JMana),stat=STATUS); VERIFY_(STATUS)
          ps_ana= dps_aux
       endif

       if (do_transforms) then
! Winds:
          if( L_CUBE ) then
             call BKG2ANA%regrid(u_bkg, v_bkg, du_aux, dv_aux, rc=status); VERIFY_(STATUS)
          else
             call BKG2ANA%regrid(u_bkg, du_aux, rc=status); VERIFY_(STATUS)
             call BKG2ANA%regrid(v_bkg, dv_aux, rc=status); VERIFY_(STATUS)
!             call POLEFIX ( du_aux,dv_aux,VM,myANA%GRIDana )
          endif

! Specific humidity:
          call BKG2ANA%regrid(q_bkg, dq_aux, rc=status); VERIFY_(STATUS)

! Pressure edges:
          call BKG2ANA%regrid(ple_bkg, dple_aux, rc=status); VERIFY_(STATUS)
          if (.not.l_use_ana_delp) then
!              Surface pressure:
               allocate(qdum2(IMana,JMana,1))
               allocate(qdum1(IMbkg,JMbkg,1))
               qdum1(:,:,1)=ps_bkg
               call BKG2ANA%regrid(qdum1, qdum2, RC=STATUS ); VERIFY_(STATUS)
               dps_aux = qdum2(:,:,1)
               deallocate(qdum1)
               deallocate(qdum2)
!              Analyzed pressure edges:
               do L=0,LM
                  ple_ana(:,:,L) = dple_aux(:,:,L) + bk(L)*( ps_ana(:,:)-dps_aux(:,:) )
               enddo
          endif

! Virtutal Temperature:
          call BKG2ANA%regrid(tv_bkg, dt_aux, rc=status); VERIFY_(STATUS)
! Ozone:
          call BKG2ANA%regrid(o3_bkg, do3_aux, rc=status); VERIFY_(STATUS)

          done_remap=.false.
          if (l_remap) then

              NPHIS = count( phis_ana.ne.myANA%phis_bkg )
              call MAPL_CommsAllReduceMax(vm,sendbuf=NPHIS,recvbuf=NPHIS_MAX,cnt=1,rc=status)
              VERIFY_(STATUS)
              if( NPHIS_MAX > 0 ) then

                 if(MAPL_AM_I_ROOT()) then
                   print *, 'Remapping ANA to Internal State Topography'
                   print *
                 endif

                 allocate(pk_ana (IMana,JMana,  LMana),stat=STATUS);VERIFY_(STATUS)
                 allocate(pke_ana(IMana,JMana,0:LMana),stat=STATUS);VERIFY_(STATUS)
                 pke_ana(:,:,:)  = ple_ana(:,:,:)**MAPL_KAPPA
                 do L=1,lm
                    pk_ana(:,:,L)  = ( pke_ana(:,:,L)-pke_ana(:,:,L-1) ) &
                                   / ( MAPL_KAPPA*log(ple_ana(:,:,L)/ple_ana(:,:,L-1)) )
                 enddo

                 allocate(thv_ana(IMana,JMana,LMana),stat=STATUS);VERIFY_(STATUS)
                 thv_ana = tv_ana/pk_ana

                 call myremap ( ple_ana, &
                                  u_ana, &
                                  v_ana, &
                                thv_ana, &
                                  q_ana, &
                                 o3_ana, &
                               phis_ana,myANA%phis_bkg,ak,bk,IMana,JMana,LMana )

                 ! Re-create ANA Dry Temperature
                 ! -----------------------------
                 pke_ana(:,:,:) = ple_ana(:,:,:)**MAPL_KAPPA
                 do L=1,LMana
                    pk_ana(:,:,L) = ( pke_ana(:,:,L)-pke_ana(:,:,L-1) ) &
                                  / ( MAPL_KAPPA*log(ple_ana(:,:,L)/ple_ana(:,:,L-1)) )
                 enddo
                 tv_ana= thv_ana*pk_ana

                 deallocate(thv_ana,stat=STATUS);VERIFY_(STATUS)
                 deallocate(pke_ana,stat=STATUS);VERIFY_(STATUS)
                 deallocate(pk_ana ,stat=STATUS);VERIFY_(STATUS)

                 done_remap=.true.
              else
                 if(MAPL_AM_I_ROOT()) then
                    print *
                    print *, 'Vertical Remapping not necessary since ANA and BKG Topographies are identical.'
                    print *
                 endif
              endif ! <phis_comp>

          endif ! <l_remap>

          ! Modify Vertically Integrated Mass-Divergence Increment
          ! ------------------------------------------------------
          if ( l_windfix )  then
             if(MAPL_AM_I_ROOT()) then
                CALL DATE_AND_TIME (REAL_CLOCK(1), REAL_CLOCK(2), &
                                    REAL_CLOCK(3), DATE_TIME)
                print *, 'Applying Mass Divergence Fix ...'
                print *, '     YYYYMMDD: ', REAL_CLOCK(1)
                print *, '     HHMMSSSS: ', REAL_CLOCK(2)
                print *
             endif
             method = 1
             allocate ( vintdiva(IMana,JMana) )
             allocate ( vintdivb(IMana,JMana) )
             allocate ( vintdivc(IMana,JMana) )
             call windfix ( u_ana,v_ana,ple_ana,                                              &
                            du_aux,dv_aux,dple_aux,IMana,JMana,LMana,VM,myANA%GRIDana,method, &
                            vintdiva,vintdivb,vintdivc                                        )
             deallocate ( vintdivc )
             deallocate ( vintdivb )
             deallocate ( vintdiva )
          endif

          if(MAPL_AM_I_ROOT()) then
             CALL DATE_AND_TIME (REAL_CLOCK(1), REAL_CLOCK(2), &
                                 REAL_CLOCK(3), DATE_TIME)
             print *, 'Calculate wind increments ...'
             print *, '     YYYYMMDD: ', REAL_CLOCK(1)
             print *, '     HHMMSSSS: ', REAL_CLOCK(2)
             print *
          endif

!         Calculate wind increment on analysis grid
          du_aux = u_ana - du_aux
          dv_aux = v_ana - dv_aux

!         Bring wind increments from analysis grid to model grid
          if( L_CUBE ) then
             call ANA2BKG%regrid(du_aux, dv_aux, du_inc, dv_inc, rc=status); VERIFY_(STATUS)
          else
             call ANA2BKG%regrid(du_aux, du_inc, rc=status); VERIFY_(STATUS)
             call ANA2BKG%regrid(dv_aux, dv_inc, rc=status); VERIFY_(STATUS)
!             call POLEFIX ( du_inc,dv_inc,VM,myANA%GRIDana )
          endif

          if(MAPL_AM_I_ROOT()) then
             CALL DATE_AND_TIME (REAL_CLOCK(1), REAL_CLOCK(2), &
                                 REAL_CLOCK(3), DATE_TIME)
             print *, 'Calculate sphum increments ...'
             print *, '     YYYYMMDD: ', REAL_CLOCK(1)
             print *, '     HHMMSSSS: ', REAL_CLOCK(2)
             print *
          endif

!         Calculate specific humdity increment on analysis grid
          dq_aux = q_ana - dq_aux
!         Bring specific humidity increment from analysis grid to model grid
          call ANA2BKG%regrid(dq_aux, dq_inc, rc=status); VERIFY_(STATUS)

          if(MAPL_AM_I_ROOT()) then
             CALL DATE_AND_TIME (REAL_CLOCK(1), REAL_CLOCK(2), &
                                 REAL_CLOCK(3), DATE_TIME)
             print *, 'Calculate pressure increments ...'
             print *, '     YYYYMMDD: ', REAL_CLOCK(1)
             print *, '     HHMMSSSS: ', REAL_CLOCK(2)
             print *
          endif

          if (l_use_ana_delp) then
!            Calculate increment to pressure thickness
             dple_aux = ple_ana - dple_aux
          else
!            Bring pressure increment from analysis grid to model grid
             if (done_remap) then
                dple_aux = ple_ana - dple_aux
             else
                do L=0,LMana
                   dple_aux(:,:,L) = bk(L)*( ps_ana(:,:)-dps_aux(:,:) )
                enddo
             endif
          endif
          call ANA2BKG%regrid(dple_aux, dple_inc, rc=status); VERIFY_(STATUS)

          if(MAPL_AM_I_ROOT()) then
             CALL DATE_AND_TIME (REAL_CLOCK(1), REAL_CLOCK(2), &
                                 REAL_CLOCK(3), DATE_TIME)
             print *, 'Calculate temperature increments ...'
             print *, '     YYYYMMDD: ', REAL_CLOCK(1)
             print *, '     HHMMSSSS: ', REAL_CLOCK(2)
             print *
          endif

!         Calculate virtual temperature increment
          dt_aux = (tv_ana - dt_aux)                                                    ! virtual temperature increment
          call ANA2BKG%regrid(dt_aux, dt_inc, rc=status); VERIFY_(STATUS)
!         Convert virtual temperature increment into dry temperature increment
          dt_inc = dt_inc/(1.0+eps*q_bkg) - eps*dq_inc*tv_bkg/((1.0+eps*q_bkg)*(1.0+eps*q_bkg)) ! dt_inc now has inc on dry temperature

!         Calculate ozone increment on analysis grid
          do3_aux = o3_ana - do3_aux
!         Bring specific humidity increment from analysis grid to model grid
          call ANA2BKG%regrid(do3_aux, do3_inc, rc=status); VERIFY_(STATUS)

          if(MAPL_AM_I_ROOT()) then
             CALL DATE_AND_TIME (REAL_CLOCK(1), REAL_CLOCK(2), &
                                 REAL_CLOCK(3), DATE_TIME)
             print *, 'Calculate tskin increments ...'
             print *, '     YYYYMMDD: ', REAL_CLOCK(1)
             print *, '     HHMMSSSS: ', REAL_CLOCK(2)
             print *
          endif

!         Calcuate T-skin increment
          if (IuseTS) then
             allocate(qdum2(IMana,JMana,1))
             allocate(qdum1(IMbkg,JMbkg,1))
             qdum1(:,:,1)=ts_bkg
             call BKG2ANA%regrid(qdum1, qdum2, rc=status); VERIFY_(STATUS)
             qdum2(:,:,1) = dts_aux - qdum2(:,:,1)
             call ANA2BKG%regrid(qdum2, qdum1, rc=status); VERIFY_(STATUS)
             VERIFY_(STATUS)
             dts_inc = qdum1(:,:,1)
             deallocate(qdum1)
             deallocate(qdum2)
          else
             dts_inc = 0.0
          endif

          if(MAPL_AM_I_ROOT()) then
             CALL DATE_AND_TIME (REAL_CLOCK(1), REAL_CLOCK(2), &
                                 REAL_CLOCK(3), DATE_TIME)
             print *, 'Finished with increments ...'
             print *, '     YYYYMMDD: ', REAL_CLOCK(1)
             print *, '     HHMMSSSS: ', REAL_CLOCK(2)
             print *
          endif

       else ! both analysis and GCM are on the same grid (likely lat/lon)
          if(MAPL_AM_I_ROOT()) then
             print *
             print *, 'Handling same resolution case'
             print *
          endif
!         Wind increment
          du_inc = du_aux - u_bkg
          dv_inc = dv_aux - v_bkg
!         Specific humidity increment
          dq_inc = dq_aux - q_bkg
!         Pressure increments
          if (l_use_ana_delp) then
             dple_inc = ple_ana - ple_bkg
          else
             do L=0,LMbkg
                dple_inc(:,:,L) = bk(L)*( dps_aux(:,:)-ps_bkg(:,:) )
             enddo
          endif
!         Virtual Temperature
          dt_inc = (tv_ana - tv_bkg)                                                            ! virtual temperature increment
          dt_inc = dt_inc/(1.0+eps*q_bkg) - eps*dq_inc*tv_bkg/((1.0+eps*q_bkg)*(1.0+eps*q_bkg)) ! dry temperature increment
!         Ozone
          do3_inc = do3_aux - o3_bkg
!         Skin temperature
          if (IuseTS) then
             dts_inc = dts_aux - ts_bkg
          else
             dts_inc = 0.0
          endif
       endif
       if (sclinc>0.0) then
           du_inc   = sclinc * du_inc
           dv_inc   = sclinc * dv_inc
           dq_inc   = sclinc * dq_inc
           dple_inc = sclinc * dple_inc
           dt_inc   = sclinc * dt_inc
           do3_inc  = sclinc * do3_inc
           dts_inc  = sclinc * dts_inc
            if(MAPL_AM_I_ROOT()) then
               print *
               print *, 'Scaled refreshed IAU tendency by: ', sclinc
            endif
       endif

       if(.not.l_use_ana_delp) then
          deallocate(ps_ana ,stat=STATUS); VERIFY_(STATUS)
       endif
       deallocate(ple_ana,stat=STATUS); VERIFY_(STATUS)
       deallocate(o3_ana, stat=STATUS); VERIFY_(STATUS)
       deallocate(q_ana,  stat=STATUS); VERIFY_(STATUS)
       deallocate(tv_ana, stat=STATUS); VERIFY_(STATUS)
       deallocate(v_ana,  stat=STATUS); VERIFY_(STATUS)
       deallocate(u_ana,  stat=STATUS); VERIFY_(STATUS)

     endif ! <nudge>

!   Clean up
!   --------
    call MAPL_FieldBundleDestroy ( RBundle, RC=STATUS)
    VERIFY_(STATUS)

    if(l_nudge) then
       if (l_use_ana_delp) then
          deallocate(ple_bkg, stat=STATUS); VERIFY_(STATUS)
          deallocate(dp_aux  ,stat=STATUS); VERIFY_(STATUS)
       else
          deallocate(dps_aux ,stat=STATUS); VERIFY_(STATUS)
       endif
       deallocate(dts_aux ,stat=STATUS); VERIFY_(STATUS)
       deallocate(du_aux  ,stat=STATUS); VERIFY_(STATUS)
       deallocate(dv_aux  ,stat=STATUS); VERIFY_(STATUS)
       deallocate(dt_aux  ,stat=STATUS); VERIFY_(STATUS)
       deallocate(dq_aux  ,stat=STATUS); VERIFY_(STATUS)
       deallocate(do3_aux ,stat=STATUS); VERIFY_(STATUS)
       deallocate(dple_aux,stat=STATUS); VERIFY_(STATUS)
    endif
    if(associated(phis_ana))deallocate(phis_ana)
    if(associated(gptr2d))deallocate(gptr2d)
    if(associated(gptr3d))deallocate(gptr3d)
    if(allocated(dps))   deallocate(dps)
    if(allocated(dp_inc))deallocate(dp_inc)
    if(allocated(uptr))  deallocate(uptr)
    if(allocated(vptr))  deallocate(vptr)


    if ( (.not.l_store_transforms) ) then
      if (associated(myANA%phis_bkg)) then
          deallocate(myANA%phis_bkg,stat=STATUS);VERIFY_(STATUS)
          nullify(myANA%phis_bkg)
      endif
      if( do_transforms ) then
        ! Currently no 'destroy' option in new regridders.
!$$         call MAPL_HorzTransformDestroy(BKG2ANA, RC=STATUS);VERIFY_(STATUS)
!$$         call MAPL_HorzTransformDestroy(ANA2BKG, RC=STATUS);VERIFY_(STATUS)
         call WRITE_PARALLEL("Destroyed transforms ANA2BKG/BKG2ANA")
      endif
      call ESMF_GridDestroy(myAna%GRIDana, rc=status);VERIFY_(STATUS)
      call WRITE_PARALLEL("Destroyed myAna%GRIDana")
      myANA%do_transforms=.false.
      myANA%initialized=.false.
    endif

    if( LAST_CORRECTOR .and. myANA%initialized ) then
      if (associated(myANA%phis_bkg)) then
          deallocate(myANA%phis_bkg,stat=STATUS);VERIFY_(STATUS)
          nullify(myANA%phis_bkg)
      endif
      if( do_transforms ) then
        ! Currently no 'destroy' option in new regridders.
!$$         call MAPL_HorzTransformDestroy(BKG2ANA, RC=STATUS);VERIFY_(STATUS)
!$$         call MAPL_HorzTransformDestroy(ANA2BKG, RC=STATUS);VERIFY_(STATUS)
         call WRITE_PARALLEL("Destroyed transforms ANA2BKG/BKG2ANA")
      endif
      call ESMF_GridDestroy(myAna%GRIDana, rc=status);VERIFY_(STATUS)
      call WRITE_PARALLEL("Destroyed myAna%GRIDana")
      myANA%do_transforms=.false.
      myANA%initialized=.false.
    endif

    RETURN_(ESMF_SUCCESS)
    end subroutine update_ainc_


    subroutine make_ana_grid(myANA, IM_world, JM_world, NX, NY, LM, rc)
       type (CONNECT_AnanBKG), intent(inout) :: myANA
       integer, intent(in) :: IM_world
       integer, intent(in) :: JM_world
       integer, intent(in) :: NX
       integer, intent(in) :: NY
       integer, intent(in) :: LM
       integer, optional, intent(out) :: rc

       type (ESMF_Config) :: tmp_config
       integer :: status
       character(len=ESMF_MAXSTR) :: imstr, jmstr

       write(imstr,*) IM_World
       write(jmstr,*) JM_World
       if(JM_World==6*IM_World) then
          myANA%gridAnaName='PE'//trim(adjustl(imstr))//'x'//trim(adjustl(jmstr))//'-CF'
          tmp_config =  MAPL_ConfigCreate()
          call MAPL_ConfigSetAttribute(tmp_config, myANA%gridAnaName, 'GRID_NAME:')
          call MAPL_ConfigSetAttribute(tmp_config, 'Cubed-Sphere','GRID_TYPE:')
          call MAPL_ConfigSetAttribute(tmp_config, NX, 'NX:')
          call MAPL_ConfigSetAttribute(tmp_config, NY, 'NY:')
          call MAPL_ConfigSetAttribute(tmp_config, IM_World, 'IM_WORLD:')
          call MAPL_ConfigSetAttribute(tmp_config, JM_World,  'JM_WORLD:')
          myANA%GridAna = grid_manager%make_grid(tmp_config)
          call ESMF_ConfigDestroy(tmp_config)
          call WRITE_PARALLEL("Created cube myAna%GRIDana...")
       else
          myANA%gridAnaName='PC'//trim(adjustl(imstr))//'x'//trim(adjustl(jmstr))//'-DC'
          myANA%GRIDana = grid_manager%make_grid( &
               & LatLonGridFactory(im_world=IM_World, jm_world=JM_World, lm=LM, &
               & nx=NX, ny=NY, pole='PC', dateline= 'DC', rc=status))
          VERIFY_(STATUS)
          call WRITE_PARALLEL("Created lat/lon myAna%GRIDana...")
       endif

    end subroutine make_ana_grid

    logical function check_list_(name,vars)
    implicit none
    character(len=*) :: name
    character(len=*) :: vars(:)
    integer ii
    check_list_=.false.
    do ii = 1,size(vars)
       if(trim(name)==trim(vars(ii))) then
          check_list_=.true.
          exit
       endif
    enddo
    end function check_list_

    subroutine dfi_coeffs (DT,FILE_FREQUENCY,TAUIAU,nsteps,dfi)
    implicit none

    real,   intent(in)  :: DT     ! model time step
    real,   intent(in)  :: TAUIAU ! IAU time scale
    integer,intent(in)  :: nsteps         ! number of steps:  FILE_FREQUENCY/DT+1
    integer,intent(in)  :: FILE_FREQUENCY
    real,   intent(out) :: dfi(nsteps)

    real pi,arg,wc
    integer n,i,k,nhlf,np1

!   Calculate DFI coefficients
!   --------------------------
    pi = 4.0*atan(1.)
    nhlf = (nsteps+1)/2
    do k = 1, nhlf-1
       n   = k-nhlf
       arg = n*pi/nhlf
       wc  = sin(arg)/arg ! Lanczos window
       dfi(k) = wc*sin(n*2.0*pi*DT/FILE_FREQUENCY)/(n*pi)
    end do
    dfi(nhlf) = 2*DT/FILE_FREQUENCY
    do i = nhlf+1, nsteps
       dfi(i) = dfi(nsteps-i+1)
    end do

!   Normalize coefficients
!   ----------------------
    dfi = dfi/sum(dfi)
    dfi = dfi*(FILE_FREQUENCY/TAUIAU)/DT ! remember: dynamics multiplies by DT

   end subroutine dfi_coeffs

  end subroutine Run

  function  my_nearest_time(TimeLeft, TimeRight, TimeNow, RC) result(TimeNearest)
  type(ESMF_Time),   intent(IN)  :: TimeLeft
  type(ESMF_Time),   intent(IN)  :: TimeRight
  type(ESMF_Time),   intent(IN)  :: TimeNow
  integer, optional, intent(OUT) :: RC
  type(ESMF_Time) :: TimeNearest

  if (TimeLeft==TimeRight) then
      TimeNearest=TimeLeft
      if(present(RC)) then
         RC = 0
      endif
      return
  endif
  if (TimeNow<=TimeLeft) then
      TimeNearest=TimeLeft
      if(present(RC)) then
         RC = 0
      endif
      return
  endif
  if (TimeNow>=TimeRight) then
      TimeNearest=TimeRight
      if(present(RC)) then
         RC = 0
      endif
      return
  endif
  if ((TimeNow-TimeLeft) < (TimeNow-TimeRight)) then
      TimeNearest=TimeLeft
  else
      TimeNearest=TimeRight
  endif
  if(present(RC)) then
     RC = 0
  endif
  end function my_nearest_time

  subroutine GET_REPLAY_TIME ( MAPL, CLOCK, REPLAY_TIME, Begin_REPLAY_Cycle, RC )


    type(ESMF_Clock),    intent(inout) :: CLOCK
    type(ESMF_Time),     intent(  out) :: REPLAY_TIME
    logical                               Begin_REPLAY_Cycle
    integer, optional,   intent(  out) :: RC

! Locals

    type(MAPL_MetaComp),      pointer   :: MAPL
    type(ESMF_Calendar)                 :: cal
    type(ESMF_Time)                     :: currtime
    type(ESMF_Time)                     :: REPLAY_TIMEP0
    type(ESMF_Time)                     :: REPLAY_TIMEM1

    character(len=ESMF_MAXSTR)          :: IAm
    character(len=ESMF_MAXSTR)          :: TimeString
    type(ESMF_TimeInterval)             :: FileFreq
    integer                             :: FileFreq_SEC
    integer                             :: FileReft_HMS
    integer                             :: FileReft_SEC
    integer                             ::    TOTAL_SEC

    integer                             :: PREDICTOR_DURATION
    real*8                              :: facp0, facm1
    integer                             :: CUR_YY,CUR_MM,CUR_DD,CUR_H,CUR_M,CUR_S
    integer                             :: nymd,  nhms
    integer                             :: rymd,  rhms
    integer                             :: Pymd,  Phms
    integer                             :: Mymd,  Mhms
    integer                             :: STATUS
    integer nsecf
            nsecf(nhms) = nhms/10000*3600 + mod(nhms,10000)/100*60 + mod(nhms,100)

!=============================================================================

   Iam = 'REPLAY_Time'

! Note: MKIAU_FREQUENCY and MKIAU_REFERENCE_TIME are initialized within GEOS_GcmGridComp
! --------------------------------------------------------------------------------------
   call MAPL_GetResource( MAPL,FileFreq_SEC, Label="MKIAU_FREQUENCY:",      rc=STATUS )
   VERIFY_(STATUS)
   call MAPL_GetResource( MAPL,FileReft_HMS, Label="MKIAU_REFERENCE_TIME:", rc=STATUS )
   VERIFY_(STATUS)

 ! Note: PREDICTOR Duration is Initialized in GCM_GridComp
 ! -------------------------------------------------------
   call MAPL_GetResource( MAPL, PREDICTOR_DURATION, Label="PREDICTOR_DURATION:", RC=STATUS )
   VERIFY_(STATUS)

   call ESMF_ClockGet( CLOCK, currTime=currTime, calendar=cal, RC=STATUS)
   VERIFY_(STATUS)

        FileReft_SEC = NSECF( FileReft_HMS )
        FileReft_SEC =   mod( FileReft_SEC,FileFreq_SEC )

        call ESMF_TimeIntervalSet( FileFreq, S=FileFreq_SEC, rc=STATUS )
        VERIFY_(STATUS)
        REPLAY_TIME = currTime + ( FileFreq / 2 )

        call ESMF_TimeGet( REPLAY_TIME, YY = CUR_YY, &
                                        MM = CUR_MM, &
                                        DD = CUR_DD, &
                                         H = CUR_H , &
                                         M = CUR_M , &
                                         S = CUR_S , &
                                        rc = STATUS  )
        VERIFY_(STATUS)
        TOTAL_SEC = CUR_H * 3600 + CUR_M * 60 + CUR_S
        TOTAL_SEC = FileFreq_SEC * ( TOTAL_SEC/FileFreq_SEC ) + FileReft_SEC

        CUR_H =       TOTAL_SEC/3600
        CUR_M =  mod( TOTAL_SEC,3600 )/60
        CUR_S =  mod( TOTAL_SEC, 60  )

        call ESMF_TimeSet( REPLAY_TIME, YY = CUR_YY, &
                                        MM = CUR_MM, &
                                        DD = CUR_DD, &
                                         H = CUR_H , &
                                         M = CUR_M , &
                                         S = CUR_S , &
                          calendar=cal, rc = STATUS  )
        VERIFY_(STATUS)

! --------------------------------------------------------------------------------------------------------

    if( currTime /= REPLAY_TIME .or. Begin_REPLAY_Cycle ) then

        if( currTime < REPLAY_TIME ) then
            REPLAY_TIMEP0 = REPLAY_TIME
            REPLAY_TIMEM1 = REPLAY_TIME - FileFreq
        else
            REPLAY_TIMEP0 = REPLAY_TIME + FileFreq
            REPLAY_TIMEM1 = REPLAY_TIME
        endif

        facp0 = (currTime-REPLAY_TIMEM1) / (REPLAY_TIMEP0-REPLAY_TIMEM1)
        facm1 = (REPLAY_TIMEP0-currTime) / (REPLAY_TIMEP0-REPLAY_TIMEM1)

      ! Backward Time
      ! -------------
      ! if( PREDICTOR_DURATION == 0 ) then
      !     if( facm1 > 0.0 ) then
      !         REPLAY_TIME = REPLAY_TIMEM1
      !     else
      !         REPLAY_TIME = REPLAY_TIMEP0
      !     endif
      ! endif

      ! Forward Time
      ! ------------
      ! if( PREDICTOR_DURATION == FileFreq_SEC ) then
      !     if( facp0 > 0.0 ) then
      !         REPLAY_TIME = REPLAY_TIMEP0
      !     else
      !         REPLAY_TIME = REPLAY_TIMEM1
      !     endif
      ! endif

      ! Nearest Time
      ! ------------
      ! if( PREDICTOR_DURATION == FileFreq_SEC/2 ) then
            if( facm1 > 0.5 ) then
                REPLAY_TIME = REPLAY_TIMEM1
            else
                REPLAY_TIME = REPLAY_TIMEP0
            endif
      ! endif

    endif

 ! if(MAPL_AM_I_ROOT() ) then
 !   call ESMF_TimeGet( currTIME, timeString=TimeString, RC=STATUS)
 !   VERIFY_(STATUS)
 !   call strToInt(TimeString, nymd, nhms)

 !   call ESMF_TimeGet(REPLAY_TIME, timeString=TimeString, RC=STATUS)
 !   VERIFY_(STATUS)
 !   call strToInt(TimeString, rymd, rhms)

 !   call ESMF_TimeGet(REPLAY_TIMEP0, timeString=TimeString, RC=STATUS)
 !   VERIFY_(STATUS)
 !   call strToInt(TimeString, Pymd, Phms)

 !   call ESMF_TimeGet(REPLAY_TIMEM1, timeString=TimeString, RC=STATUS)
 !   VERIFY_(STATUS)
 !   call strToInt(TimeString, Mymd, Mhms)

 !   write(6,'(1x,a,i8.8,a,i6.6,a,f5.3,a,a,i8.8,a,i6.6,a,i8.8,a,i6.6,a,i8.8,a,i6.6,a,f5.3,a,l)') &
 !                                            ' Current_Time: ',nymd,' ',nhms,' (',facm1,') ', &
 !                                            ' -Replay_Time: ',Mymd,' ',Mhms, &
 !                                            '  Replay_Time: ',rymd,' ',rhms, &
 !                                            ' +Replay_Time: ',Pymd,' ',Phms,' (',facp0,') Begin_REPLAY_Cycle: ',Begin_REPLAY_Cycle
 ! endif

! --------------------------------------------------------------------------------------------------------

  ! RETURN_(ESMF_SUCCESS)
    return
  end subroutine GET_REPLAY_TIME
end module GEOS_AgcmGridCompMod
