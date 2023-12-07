! $Id: GEOS_PhysicsGridComp.F90,v 1.119.4.3.2.7.2.1.8.1.2.13.6.1.20.1.2.1.2.5.2.2.2.1.12.1 2019/07/23 15:30:41 mmanyin Exp $

#include "MAPL_Generic.h"

#define debug 0

!=============================================================================
!BOP

! !MODULE: GEOS_PhysicsGridCompMod -- A Module to combine Short-Wave, Long-Wave Radiation Moist-Physics and Turbulence Gridded Components

! !INTERFACE:

module GEOS_PhysicsGridCompMod

! !USES:

  use ESMF
  use MAPL
  use stoch_module

  use GEOS_SurfaceGridCompMod,    only : SurfSetServices      => SetServices
  use GEOS_MoistGridCompMod,      only : MoistSetServices     => SetServices
  use GEOS_TurbulenceGridCompMod, only : TurblSetServices     => SetServices
  use GEOS_RadiationGridCompMod,  only : RadiationSetServices => SetServices
  use GEOS_ChemGridCompMod,       only : AChemSetServices     => SetServices
  use GEOS_GwdGridCompMod,        only : GwdSetServices       => SetServices

  use GEOS_UtilsMod, only: GEOS_Qsat
  use Bundle_IncrementMod
  use MBundle_IncrementMod

! PGI Module that contains the initialization
! routines for the GPUs
#ifdef _CUDA
  use cudafor
#endif

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION: This gridded component (GC) combines the Radiation (Short-Wave and Long-Wave),
!   Moist-Physics, Chem, Surface and Turbulence GCs into a new composite Physics GC.
!   The Export Couplings of the Physics GC are the union of the Export
!   Couplings of the individual child GCs, plus the combined tendencies needed by
!   the dynamics. These last are the pressure-weighted tendencies of the atmospheric
!   state variables U,V,T (due to external diabatic forcing), the tendency of the
!   edge pressures, and a collection of "Friendly" tracers for advection. In the current
!   version, the only friendly tracers are variables from Moist-Physics and Chem.
!
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

  integer ::        GWD
  integer ::        SURF
  integer ::        CHEM
  integer ::        MOIST
  integer ::        TURBL
  integer ::        RAD

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer,             intent(  OUT) :: RC  ! return code

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

    type (MAPL_MetaComp),  pointer          :: MAPL
    CHARACTER(LEN=ESMF_MAXSTR)              :: RATsProviderName
    integer                                 :: I
    type (ESMF_Config)                      :: CF

    integer                                 :: DO_OBIO, DO_CO2CNNEE, ATM_CO2, nCols, NQ

    real                                    :: SYNCTQ
    character(len=ESMF_MAXSTR), allocatable :: NAMES(:)
    character(len=ESMF_MAXSTR)              :: TendUnits
    character(len=ESMF_MAXSTR)              :: SURFRC
    type(ESMF_Config)                       :: SCF

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "::" // Iam

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,        RC=STATUS )
    VERIFY_(STATUS)

! Create children`s gridded components and invoke their SetServices
! -----------------------------------------------------------------

! Note chemistry must be added before surface so that chemistry has
! a change to put the fields in the AERO_DP bundle. Otherwise when
! surface reads the import restart AERO_DP will be empty and it will
! not be properly restarted
    GWD = MAPL_AddChild(GC, NAME='GWD', SS=GwdSetServices, RC=STATUS)
    VERIFY_(STATUS)
    MOIST = MAPL_AddChild(GC, NAME='MOIST', SS=MoistSetServices, RC=STATUS)
    VERIFY_(STATUS)
    TURBL = MAPL_AddChild(GC, NAME='TURBULENCE', SS=TurblSetServices, RC=STATUS)
    VERIFY_(STATUS)
    CHEM = MAPL_AddChild(GC, NAME='CHEMISTRY', SS=AChemSetServices, RC=STATUS)
    VERIFY_(STATUS)
    SURF = MAPL_AddChild(GC, NAME='SURFACE', SS=SurfSetServices, RC=STATUS)
    VERIFY_(STATUS)
    RAD = MAPL_AddChild(GC, NAME='RADIATION', SS=RadiationSetServices, RC=STATUS)
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, DO_OBIO, Label="USE_OCEANOBIOGEOCHEM:",DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource (MAPL, SURFRC, label = 'SURFRC:', default = 'GEOS_SurfaceGridComp.rc', RC=STATUS) ; VERIFY_(STATUS)
    SCF = ESMF_ConfigCreate(rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigLoadFile(SCF,SURFRC,rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute (SCF, label='ATM_CO2:',   value=ATM_CO2,   DEFAULT=0, __RC__ )
    call ESMF_ConfigDestroy      (SCF, __RC__)

    SCF = ESMF_ConfigCreate(rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigLoadFile(SCF,'CO2_GridComp.rc',rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute (SCF, label='USE_CNNEE:', value=DO_CO2CNNEE,   DEFAULT=0, __RC__ )
    call ESMF_ConfigDestroy      (SCF, __RC__)

! Get SYNCTQ flag from config to know whether to terminate some imports
! ---------------------------------------------------------------------------
    call MAPL_GetResource ( MAPL, SYNCTQ, Label="SYNCTQ:", DEFAULT= 1.0, RC=STATUS)
    VERIFY_(STATUS)
!BOS

! !INTERNAL STATE
    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QW',                                        &
         LONG_NAME  = 'mass_fraction_of_wet_air',                  &
         UNITS      = 'kg kg-1',                                   &
         RESTART    = MAPL_RestartSkip,                            &
         FRIENDLYTO = 'TURBULENCE:MOIST',                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)


! !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'W',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'DZ',                                        &
         LONG_NAME  = 'surface_layer_height',                      &
         UNITS      = 'm',                                         &
         DIMS       =  MAPL_DimsHorzOnly,                          &
         VLOCATION  =  MAPL_VLocationNone,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'T',                                         &
         LONG_NAME  = 'air_temperature',                           &
         UNITS      = 'K',                                         &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'S',                                         &
         LONG_NAME  = 'dry_static_energy',                         &
         UNITS      = 'm+2 s-2',                                   &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'ZLE',                                       &
         LONG_NAME  = 'geopotential_height',                       &
         UNITS      = 'm',                                         &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationEdge,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'PLE',                                       &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationEdge,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'DTDTDYN',                                     &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_dynamics', &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                      &
         SHORT_NAME = 'DQVDTDYN',                                      &
         LONG_NAME  = 'tendency_of_specific_humidity_due_to_dynamics', &
         UNITS      = 'kg/kg/s',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'QV_DYN_IN',                                 &
         LONG_NAME  = 'spec_humidity_at_begin_of_time_step',       &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'T_DYN_IN',                                 &
         LONG_NAME  = 'temperature_at_begin_of_time_step',       &
         UNITS      = 'K',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'U_DYN_IN',                                 &
         LONG_NAME  = 'u_wind_at_begin_of_time_step',       &
         UNITS      = 'm s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'V_DYN_IN',                                 &
         LONG_NAME  = 'v_wind_at_begin_of_time_step',       &
         UNITS      = 'm s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'PLE_DYN_IN',                                 &
         LONG_NAME  = 'edge_pressure_at_begin_of_time_step',       &
         UNITS      = 'Pa',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'AREA',                                      &
         LONG_NAME  = 'Grid-Cell Area',                            &
         UNITS      = 'm+2',                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

! !EXPORT STATE:

!   Add export states for turbulence increments
!-------------------------------------------------
    call ESMF_ConfigGetDim (cf, NQ, nCols, label=('TRI_increments::'), rc=STATUS)
    if (NQ > 0) then
      call ESMF_ConfigFindLabel (cf, ('TRI_increments::'), rc=STATUS)
      VERIFY_(STATUS)
      allocate (NAMES(NQ), stat=STATUS)
      VERIFY_(STATUS)
      do i = 1, NQ
        call ESMF_ConfigNextLine(cf, rc=STATUS)
        VERIFY_(STATUS)
        call ESMF_ConfigGetAttribute(cf, NAMES(i), rc=STATUS)
        VERIFY_(STATUS)
      enddo
      do i = 1, NQ
        if (NAMES(i) == 'AOADAYS') then
          TendUnits = 'days s-1'
        else
          TendUnits = 'UNITS'
        end if
        call MAPL_AddExportSpec(GC,                                           &
          SHORT_NAME =  trim(NAMES(i))//'IT',                                 &
          LONG_NAME  = 'tendency_of_'//trim(NAMES(i))//'_due_to_turbulence',  &
          UNITS      =  TendUnits,                                            &
          DIMS       =  MAPL_DimsHorzVert,                                    &
          VLOCATION  =  MAPL_VLocationCenter,                                 &
          RC=STATUS  )
        VERIFY_(STATUS)
      end do
      deallocate(NAMES)
    end if !NQ > 0

!-----------------------------------------------------------
    call MAPL_AddExportSpec(GC,                                                       &
         SHORT_NAME = 'DTDT',                                                         &
         LONG_NAME  = 'pressure_weighted_tendency_of_air_temperature_due_to_physics', &
         UNITS      = 'Pa K s-1',                                                     &
         DIMS       =  MAPL_DimsHorzVert,                                             &
         VLOCATION  =  MAPL_VLocationCenter,                                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'DTDTTOT',                                    &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_physics', &
         UNITS      = 'K s-1',                                      &
         DIMS       =  MAPL_DimsHorzVert,                           &
         VLOCATION  =  MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DTDTRAD',                                      &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_radiation', &
         UNITS      = 'K s-1',                                        &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'tendency_of_eastward_wind_due_to_physics',  &
         UNITS      = 'm s-2',                                     &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'tendency_of_northward_wind_due_to_physics', &
         UNITS      = 'm s-2',                                     &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'DWDT',                                      &
         LONG_NAME  = 'tendency_of_vertical_velocity_due_to_physics',  &
         UNITS      = 'm s-2',                                     &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'DPDTPHY',                                   &
         LONG_NAME  = 'tendency_of_pressure_at_bottom_edges_levels_due_to_physics',&
         UNITS      = 'Pa s-1',                                    &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'tendency_of_pressure_at_layer_edges_due_to_physics',&
         UNITS      = 'Pa s-1',                                    &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'DMDT',                                      &
         LONG_NAME  = 'vertically_integrated_mass_tendency_due_to_physics', &
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       =  MAPL_DimsHorzOnly,                          &
         VLOCATION  =  MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'TIM',                                       &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_moist_processes',&
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'TIMFRIC',                                   &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_moist_processes_friction',&
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'SIT',                                       &
         LONG_NAME  = 'pressure_weighted_tendency_of_dry_static_energy_due_to_turbulence',&
         UNITS      = 'Pa m+2 s-3',                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                        &
         SHORT_NAME = 'TIT',                                           &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_turbulence', &
         UNITS      = 'K s-1',                                         &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                      &
         SHORT_NAME = 'UIT',                                         &
         LONG_NAME  = 'tendency_of_eastward_wind_due_to_turbulence', &
         UNITS      = 'm s-2',                                       &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'VIT',                                          &
         LONG_NAME  = 'tendency_of_northward_wind_due_to_turbulence', &
         UNITS      = 'm s-2',                                        &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,                           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME = 'QVIT',                                            &
         LONG_NAME  = 'tendency_of_specific_humidity_due_to_turbulence', &
         UNITS      = 'kg kg-1 s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME = 'QLLSIT',                                          &
         LONG_NAME  = 'tendency_of_liquid_condensate_due_to_turbulence', &
         UNITS      = 'kg kg-1 s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME = 'QILSIT',                                          &
         LONG_NAME  = 'tendency_of_frozen_condensate_due_to_turbulence', &
         UNITS      = 'kg kg-1 s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                      &
         SHORT_NAME = 'TIF',                                         &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_friction', &
         UNITS      = 'K s-1',                                       &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'TRADV ',                                    &
         LONG_NAME  = 'advected_quantities',                       &
         UNITS      = 'X',                                         &
         DATATYPE   = MAPL_BundleItem,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'H2ORTRI',                                   &
         LONG_NAME  = 'H2O_rescale_increments',                    &
         UNITS      = 'UNITS s-1',                                 &
         DATATYPE   = MAPL_BundleItem,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'MTRI',                                      &
         LONG_NAME  = 'moist_quantities',                          &
         UNITS      = 'UNITS s-1',                                 &
         DATATYPE   = MAPL_BundleItem,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'MCHEMTRI',                                      &
         LONG_NAME  = 'moist_quantities',                          &
         UNITS      = 'UNITS s-1',                                 &
         DATATYPE   = MAPL_BundleItem,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         LONG_NAME  = 'upward_net_turbulence_heat_flux',           &
         UNITS      = 'W m-2',                                     &
         SHORT_NAME = 'FTB',                                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         LONG_NAME  = 'upward_net_turbulence_eastward_momentum_flux', &
         UNITS      = 'm+2 s-2',                                      &
         SHORT_NAME = 'FTU',                                          &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationEdge,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                        &
         LONG_NAME  = 'upward_net_turbulence_northward_momentum_flux', &
         UNITS      = 'm+2 s-2',                                       &
         SHORT_NAME = 'FTV',                                           &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationEdge,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'TRANA',                                     &
         LONG_NAME  = 'analyzed_quantities',                       &
         UNITS      = 'X',                                         &
         DATATYPE   = MAPL_BundleItem,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KEPHY',                                                              &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_across_physics',       &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEPHY',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_across_physics',     &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PERAD',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_across_radiation',   &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PETRB',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_across_turbulence',  &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEMST',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_across_moist',       &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEFRI',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_friction',    &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEGWD',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_across_gwd',         &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PECUF',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_cumulus_friction', &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                            &
         SHORT_NAME = 'DQVDTTRBINT',                                                       &
         LONG_NAME  = 'vertically_integrated_water_vapor_tendency_due_to_turbulence',      &
         UNITS      = 'kg m-2 s-1',                                                        &
         DIMS       = MAPL_DimsHorzOnly,                                                   &
         VLOCATION  = MAPL_VLocationNone,                                                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                            &
         SHORT_NAME = 'DQVDTMSTINT',                                                       &
         LONG_NAME  = 'vertically_integrated_water_vapor_tendency_due_to_moist_processes', &
         UNITS      = 'kg m-2 s-1',                                                        &
         DIMS       = MAPL_DimsHorzOnly,                                                   &
         VLOCATION  = MAPL_VLocationNone,                                                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                            &
         SHORT_NAME = 'DQVDTCHMINT',                                                       &
         LONG_NAME  = 'vertically_integrated_water_vapor_tendency_due_to_chemistry',       &
         UNITS      = 'kg m-2 s-1',                                                        &
         DIMS       = MAPL_DimsHorzOnly,                                                   &
         VLOCATION  = MAPL_VLocationNone,                                                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                            &
         SHORT_NAME = 'DQLDTMSTINT',                                                       &
         LONG_NAME  = 'vertically_integrated_liquid_water_tendency_due_to_moist_processes',&
         UNITS      = 'kg m-2 s-1',                                                        &
         DIMS       = MAPL_DimsHorzOnly,                                                   &
         VLOCATION  = MAPL_VLocationNone,                                                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                            &
         SHORT_NAME = 'DQIDTMSTINT',                                                       &
         LONG_NAME  = 'vertically_integrated_ice_tendency_due_to_moist_processes',         &
         UNITS      = 'kg m-2 s-1',                                                        &
         DIMS       = MAPL_DimsHorzOnly,                                                   &
         VLOCATION  = MAPL_VLocationNone,                                                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                            &
         SHORT_NAME = 'DOXDTCHMINT',                                                       &
         LONG_NAME  = 'vertically_integrated_odd_oxygen_tendency_due_to_chemistry',        &
         UNITS      = 'kg m-2 s-1',                                                        &
         DIMS       = MAPL_DimsHorzOnly,                                                   &
         VLOCATION  = MAPL_VLocationNone,                                                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                            &
         SHORT_NAME = 'DQVDTPHYINT',                                                       &
         LONG_NAME  = 'vertically_integrated_water_vapor_tendency_due_to_physics',         &
         UNITS      = 'kg m-2 s-1',                                                        &
         DIMS       = MAPL_DimsHorzOnly,                                                   &
         VLOCATION  = MAPL_VLocationNone,                                                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                            &
         SHORT_NAME = 'DQLDTPHYINT',                                                       &
         LONG_NAME  = 'vertically_integrated_liquid_water_tendency_due_to_physics',        &
         UNITS      = 'kg m-2 s-1',                                                        &
         DIMS       = MAPL_DimsHorzOnly,                                                   &
         VLOCATION  = MAPL_VLocationNone,                                                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                            &
         SHORT_NAME = 'DQIDTPHYINT',                                                       &
         LONG_NAME  = 'vertically_integrated_ice_tendency_due_to_physics',                 &
         UNITS      = 'kg m-2 s-1',                                                        &
         DIMS       = MAPL_DimsHorzOnly,                                                   &
         VLOCATION  = MAPL_VLocationNone,                                                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                            &
         SHORT_NAME = 'DOXDTPHYINT',                                                       &
         LONG_NAME  = 'vertically_integrated_odd_oxygen_tendency_due_to_physics',          &
         UNITS      = 'kg m-2 s-1',                                                        &
         DIMS       = MAPL_DimsHorzOnly,                                                   &
         VLOCATION  = MAPL_VLocationNone,                                                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DQVDTSCL',                                     &
         LONG_NAME  = 'tendency_of_water_vapor_due_to_mass_scaling',  &
         UNITS      = 'kg m-2 s-1',                                   &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DQLDTSCL',                                     &
         LONG_NAME  = 'tendency_of_cloud_water_due_to_mass_scaling',  &
         UNITS      = 'kg m-2 s-1',                                   &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DQIDTSCL',                                     &
         LONG_NAME  = 'tendency_of_cloud_ice_due_to_mass_scaling',    &
         UNITS      = 'kg m-2 s-1',                                   &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DTDTUNPERT',                                   &
         LONG_NAME  = 'unperturbed_air_temperature_tendency',         &
         UNITS      = 'K s-1',                                        &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DUDTUNPERT',                                   &
         LONG_NAME  = 'unperturtbed_tendency_of_eastward_wind',       &
         UNITS      = 'm s-2 s-1',                                    &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DVDTUNPERT',                                   &
         LONG_NAME  = 'unperturtbed_tendency_of_northward_wind',      &
         UNITS      = 'm s-2 s-1',                                    &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'QVUNPERT',                                     &
         LONG_NAME  = 'unperturbed_water_vapor',                      &
         UNITS      = 'kg kg-1',                                      &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DQVDTUNPERT',                                  &
         LONG_NAME  = 'unperturbed_tendency_of_water_vapor',          &
         UNITS      = 'kg kg-1 s-1',                                  &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'QVPERT',                                       &
         LONG_NAME  = 'stochastically_pert_water_vapor',              &
         UNITS      = 'kg kg-1',                                      &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DQVDTPERT',                                    &
         LONG_NAME  = 'stochastically_pert_tendency_of_water_vapor',  &
         UNITS      = 'kg kg-1 s-1',                                   &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DUDTSTOCH',                                    &
         LONG_NAME  = 'eastward_wind_tendency_due_to_stochastic_physics', &
         UNITS      = 'm s-2 s-1',                                    &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DVDTSTOCH',                                    &
         LONG_NAME  = 'northward_wind_tendency_due_to_stochastic_physics', &
         UNITS      = 'm s-2 s-1',                                    &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DTDTSTOCH',                                    &
         LONG_NAME  = 'air_temperature_tendency_due_to_stochastic_physics', &
         UNITS      = 'K s-1',                                        &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'DQVDTSTOCH',                                   &
         LONG_NAME  = 'water_vapor_tendency_due_to_stochastic_physics', &
         UNITS      = 'kg kg-1 s-1',                                  &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'RNDPTR',                                       &
         LONG_NAME  = 'sppt_stochastic_pattern',                      &
         UNITS      = 'n/a',                                          &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'SKEBU',                                        &
         LONG_NAME  = 'skeb_perturbation_for_eastward_wind_tendency', &
         UNITS      = 'm s-2',                                        &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME = 'SKEBV',                                        &
         LONG_NAME  = 'skeb_perturbation_for_northward_wind_tendency',&
         UNITS      = 'm s-2',                                        &
         DIMS       =  MAPL_DimsHorzVert,                             &
         VLOCATION  =  MAPL_VLocationCenter,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

! Ozone (ppmv) and Odd Oxygen (mol/mol)
!   Note: GMI currently provides just O3 as Odd Oxygen
! ----------------------------------------------------------
    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'O3PPMV',                                    &
         CHILD_ID   = CHEM,                                        &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'OX',                                        &
         CHILD_ID   = CHEM,                                        &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

! The following are exported up for Atmos Ana purposes
! ----------------------------------------------------
    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PPBL',                                      &
         CHILD_ID   = TURBL,                                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'Q',                                         &
         CHILD_ID = MOIST,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'QCTOT',                                     &
         CHILD_ID = MOIST,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'U10M',                                      &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'V10M',                                      &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'U10N',                                      &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'V10N',                                      &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'SNOMAS',                                    &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'WET1',                                      &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'TSOIL1',                                    &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'LWI',                                       &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'TS',                                        &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'FRLAND',                                    &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'FRLANDICE',                                 &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'FRLAKE',                                    &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'FROCEAN',                                   &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'FRACI',                                     &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'Z0',                                        &
         CHILD_ID = SURF,                                          &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

!EOS


! Set internal connections between the children`s IMPORTS and EXPORTS
! -------------------------------------------------------------------


! !CONNECTIONS:

! Turbulence imports
!-------------------

    call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/ 'RADLW ', 'RADLWC' /),                   &
         DST_ID      =  TURBL,                                     &
         SRC_ID      =  RAD,                                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'QV    ','QLTOT ','QITOT ','QCTOT ',      &
                         'WTHV2 ','WQT_DC'                   /),   &
         DST_ID      = TURBL,                                      &
         SRC_ID      = MOIST,                                      &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddConnectivity ( GC,                               &
         SHORT_NAME  = (/'CT   ','CM   ','CQ   ',                  &
                         'BSTAR','USTAR'              /),          &
         DST_ID      = TURBL,                                      &
         SRC_ID      = SURF,                                       &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'FRLAND','EVAP  ','SH    '/),             &
         DST_ID      = TURBL,                                      &
         SRC_ID      = SURF,                                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

! Radiation Imports
!-------------------

     call MAPL_AddConnectivity ( GC,                            &
         SHORT_NAME  = [character(len=4) ::                     &
                           'QV','QL','QI','QR','QS','QG',       &
                         'FCLD','RL','RI','RR','RS','RG'],      &
         DST_ID      = RAD,                                     &
         SRC_ID      = MOIST,                                   &
                                                     RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddConnectivity ( GC,                               &
         SHORT_NAME  = (/'ALBVR  ','ALBVF  ','ALBNR  ','ALBNF  ',  &
                         'EMIS   ','TS     ',                      &
                         'FRLAND ','FROCEAN'                  /),  &
         DST_ID      = RAD,                                        &
         SRC_ID      = SURF,                                       &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

! -----------------------------------------------------------------
! Radiatively active species and required units
!  Specie  Units   Contents
!  ------ -------  -------------------
!  OX     mol/mol  Odd oxygen or ozone volume mixing ratio
!  O3      kg/kg   Ozone mass fraction
!  CH4    mol/mol  Methane
!  N2O    mol/mol  Nitrous oxide
!  CFC11  mol/mol  CFC-11 (CCl3F)
!  CFC12  mol/mol  CFC-12 (CCl2F2)
!  HCFC22 mol/mol  HCFC-22 (CHClF2)
! -----------------------------------------------------------------
      CALL MAPL_AddConnectivity( GC, &
            SHORT_NAME  = (/'OX    ','O3    ','CH4   ','N2O   ', &
                            'CFC11 ','CFC12 ','HCFC22'       /), &
                           DST_ID=RAD, SRC_ID=CHEM, RC=STATUS    )
      VERIFY_(STATUS)
! -----------------------------------------------------------------

     call MAPL_AddConnectivity ( GC,                               &
         SHORT_NAME  = (/'AERO'/),                                 &
         DST_ID      =  RAD,                                       &
         SRC_ID      =  CHEM,                                      &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddConnectivity ( GC,                               &
         SRC_NAME    = 'TS',                                       &
         DST_NAME    = 'TSINST',                                   &
         SRC_ID      = SURF,                                       &
         DST_ID      = RAD,                                        &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

! Surface Imports
!----------------

    call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'PCU    ', 'PLS    ', 'SNO    ',          &
                         'ICE    ', 'FRZR   ', 'TPREC  ',          &
                         'CN_PRCP' /),                             &
         DST_ID      = SURF,                                       &
         SRC_ID      = MOIST,                                      &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddConnectivity ( GC,                               &
         SHORT_NAME  = (/'ALW   ','BLW   ',                        &
                         'DRPARN','DFPARN','DRNIRN',               &
                         'DFNIRN','DRUVRN','DFUVRN'    /),         &
         DST_ID      = SURF,                                       &
         SRC_ID      = RAD,                                        &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddConnectivity ( GC,                               &
         SRC_NAME    = 'LWS0',                                     &
         DST_NAME    = 'LWDNSRF',                                  &
         SRC_ID      = RAD,                                        &
         DST_ID      = SURF,                                       &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

     IF((DO_OBIO /= 0) .OR. (ATM_CO2 == 4)) THEN
        call MAPL_AddConnectivity ( GC,                               &
             SRC_NAME    = 'CO2SC001',                                &
             DST_NAME    = 'CO2SC',                                   &
             SRC_ID      = CHEM,                                      &
             DST_ID      = SURF,                                      &
             RC=STATUS  )
        VERIFY_(STATUS)
     ENDIF

     call MAPL_AddConnectivity ( GC,                               &
         SHORT_NAME  = (/'AERO_DP'/),                              &
         SRC_ID      = CHEM,                                       &
         DST_ID      = SURF,                                       &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddConnectivity ( GC,                               &
         SHORT_NAME  = (/'FSWBAND  ', 'FSWBANDNA'/),               &
         SRC_ID      = RAD,                                        &
         DST_ID      = SURF,                                       &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

! Imports for GWD
!----------------
    call MAPL_AddConnectivity ( GC,                                    &
         SHORT_NAME  = [character(len=7):: 'Q', 'DTDT_DC', 'DTDT_SC'], &
         DST_ID      = GWD,                                            &
         SRC_ID      = MOIST,                                          &
                                                        RC=STATUS      )
    VERIFY_(STATUS)
    call MAPL_AddConnectivity ( GC,                                      &
         SRC_NAME    = 'DQIDT_micro',                                    &
         DST_NAME    = 'DQIDT',                                          &
         DST_ID      = GWD,                                              &
         SRC_ID      = MOIST,                                            &
                                                       RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddConnectivity ( GC,                                      &
         SRC_NAME    = 'DQLDT_micro',                                    &
         DST_NAME    = 'DQLDT',                                          &
         DST_ID      = GWD,                                              &
         SRC_ID      = MOIST,                                            &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

! Chemistry Imports
! -----------------

     call MAPL_AddConnectivity ( GC,                                      &
        SHORT_NAME  = (/ 'RL      ',  'QL      ', 'QLTOT   ', 'DQLDT   ', &
                         'RI      ',  'QI      ', 'QITOT   ', 'DQIDT   ', &
                         'QLCN    ',  'PFL_CN  ', 'PFL_LSAN',             &
                         'QICN    ',  'PFI_CN  ', 'PFI_LSAN',             &
                         'FCLD    ',  'QCTOT   ', 'CNV_QC  ',             &
                         'REV_LS  ',  'REV_AN  ', 'REV_CN  ', 'TPREC   ', &
                         'Q       ',  'DQDT    ', 'DQRL    ', 'DQRC    ', &
                         'CNV_MFC ',  'CNV_MFD ', 'CNV_CVW ', 'CNV_FRC ', &
                         'LFR_GCC ',  'RH2     ', 'CN_PRCP ',             &
                         'BYNCY   ',  'CAPE    ', 'INHB    ' /),          &
        DST_ID      = CHEM,                                               &
        SRC_ID      = MOIST,                                              &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddConnectivity ( GC,                              &
        SRC_NAME  = (/ 'Q         ','FCLD      ' /),              &
        DST_NAME  = (/ 'Q_avg24   ','FCLD_avg24' /),              &
        DST_ID      = CHEM,                                       &
        SRC_ID      = MOIST,                                      &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddConnectivity ( GC,                              &
         SHORT_NAME  = (/'ZPBL','PPBL'/),                         &
         DST_ID      = CHEM,                                      &
         SRC_ID      = TURBL,                                     &
                                                        RC=STATUS )
     VERIFY_(STATUS)

     call MAPL_AddConnectivity ( GC,                              &
        SHORT_NAME  = (/ 'LWI      ', 'FRLAND   ', 'FRLANDICE',   &
                         'FROCEAN  ', 'FRLAKE   ', 'WET1     ',   &
                         'GRN      ', 'USTAR    ', 'U10M     ',   &
                         'V10M     ', 'SH       ', 'Z0H      ',   &
                         'LAI      ', 'TSOIL1   ', 'FRACI    ',   &
                         'TA       ', 'T2M      ', 'SWNDSRF  ',   &
                         'ALBVF    ', 'ASNOW    ', 'U10N     ',   &
                         'V10N     ', 'TS       ', 'CM       ',   &
                         'CN       ', 'RHOS     ', 'WET2     ',   &
                         'SNOMAS   ', 'SNOWDP   ', 'ITY      ',   &
                         'LHFX     ', 'Q2M      ', 'Q10M     ',   &
                         'T10M     ', 'WCSF     ', 'CN_PRCP  ',   &
                         'PRECTOT  '                          /), &
        DST_ID      = CHEM,                                       &
        SRC_ID      = SURF,                                       &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     if (DO_CO2CNNEE == 1) then
        call MAPL_AddConnectivity ( GC,                           &
             SHORT_NAME  = (/'CNNEE'/),                           &
             DST_ID      = CHEM,                                  &
             SRC_ID      = SURF,                                  &
             RC=STATUS )
        VERIFY_(STATUS)
     endif

     call MAPL_AddConnectivity ( GC,                              &
        SHORT_NAME  = (/ 'TAUCLI', 'TAUCLW', 'CLDTT ',            &
                         'DFPAR ', 'DRPAR '                   /), &
        DST_ID      =  CHEM,                                      &
        SRC_ID      =  RAD,                                       &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddConnectivity ( GC,                              &
        SRC_NAME  = (/ 'TAUCLI      ', 'TAUCLW      '/),          &
        DST_NAME  = (/ 'TAUCLI_avg24', 'TAUCLW_avg24'/),          &
        DST_ID      =  CHEM,                                      &
        SRC_ID      =  RAD,                                       &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

! Moist Imports
!--------------

    call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'U','V','T'/),                            &
         DST_ID      = MOIST,                                      &
         SRC_ID      = GWD,                                        &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddConnectivity ( GC,                                          &
         SHORT_NAME  = (/'KH           ', 'KPBL         ', 'KPBL_SC      ',     &
                         'TKE          ', 'TKESHOC      ', 'PDF_A        ',     &
                         'SL2          ', 'SL3          ', 'W2           ',     &
                         'W3           ', 'SLQT         ', 'WQT          ',     &
                         'WSL          ', 'QT2          ', 'QT3          '/),    &
         DST_ID      = MOIST,                                      &
         SRC_ID      = TURBL,                                      &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'TS' /),                                  &
         DST_ID      = MOIST,                                      &
         SRC_ID      = SURF,                                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'SNOMAS   ','FRLAND   ','FROCEAN  ',      &
                         'FRLANDICE','FRACI    '/),                &
         DST_ID      = MOIST,                                      &
         SRC_ID      = SURF,                                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'T2M ', 'Q2M ', 'TA  ', 'QA  ', 'SH  ',   &
                         'EVAP'                                    &
                       /),                                         &
         DST_ID      = MOIST,                                      &
         SRC_ID      = SURF,                                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)


    !-------------- DONIF Additional Moist Imports


    call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'VSCSFC'/),                         &
         DST_ID      = MOIST,                                      &
         SRC_ID      = TURBL,                                      &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    !Aerosol
    call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'AERO'/),                           &
         DST_ID      =  MOIST,                                     &
         SRC_ID      =  CHEM,                                      &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

   call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'TAUOROX'/),                                 &
         DST_ID      =  MOIST,                                     &
         SRC_ID      =  GWD,                                      &
                                                        RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'TAUOROY'/),                                 &
         DST_ID      =  MOIST,                                     &
         SRC_ID      =  GWD,                                      &
                                                        RC=STATUS  )
   VERIFY_(STATUS)


   call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'RADLW'/),                                 &
         DST_ID      =  MOIST,                                     &
         SRC_ID      =  RAD,                                      &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'RADSW'/),                                 &
         DST_ID      =  MOIST,                                     &
         SRC_ID      =  RAD,                                      &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'ALH'/),                                 &
         DST_ID      =  MOIST,                                     &
         SRC_ID      =  TURBL,                                      &
                                                        RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/'TAUX', 'TAUY'/),                                 &
         DST_ID      =  MOIST,                                     &
         SRC_ID      =  SURF,                                      &
                                                        RC=STATUS  )

   VERIFY_(STATUS)


!EOP

! Disable connectivities of Surface imports that are filled manually from
!  turbulence bundles.
!------------------------------------------------------------------------

     call MAPL_TerminateImport    ( GC,   &
          SHORT_NAME = (/'SH   ','TAUX ','TAUY ','EVAP ','DEWL ','FRSL ',     &
                         'DSH  ','DFU  ','DFV  ','DEVAP','DDEWL','DFRSL'/),   &
          CHILD      = SURF,           &
          RC=STATUS  )
     VERIFY_(STATUS)

! terminate imports to SURF for SYNCTQ
     if ( SYNCTQ.ge.1.) then
       call MAPL_TerminateImport    ( GC,  &
          SHORT_NAME = [character(len=5) :: 'UA','VA','TA','QA','SPEED' ], &
          CHILD      = SURF,               &
          RC=STATUS  )
       VERIFY_(STATUS)
     endif

     call MAPL_TerminateImport    ( GC,        &
          SHORT_NAME = (/'TR ','TRG','DTG' /), &
          CHILD      = TURBL,                  &
          RC=STATUS  )
     VERIFY_(STATUS)

! terminate imports to turb for SYNCTQ
     if ( SYNCTQ.ge.1.) then
       call MAPL_TerminateImport    ( GC,  &
          SHORT_NAME = (/'U ','V ','T ','TH' /),     &
          CHILD      = TURBL,              &
          RC=STATUS  )
       VERIFY_(STATUS)
     endif

     call MAPL_TerminateImport    ( GC,    &
          SHORT_NAME = (/'MTR'/),          &
          CHILD = MOIST,                   &
          RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_TerminateImport    ( GC,    &
          SHORT_NAME = (/'DQDT_BL','DTDT_BL'/),          &
          CHILD = MOIST,                   &
          RC=STATUS)
     VERIFY_(STATUS)

! terminate imports to chem for SYNCTQ
     if ( SYNCTQ.eq.1.) then
       call MAPL_TerminateImport  ( GC,    &
          SHORT_NAME = (/'T ','TH'/),      &
          CHILD = CHEM,                    &
          RC=STATUS)
       VERIFY_(STATUS)
     endif

! terminate imports to RAD for SYNCTQ
     if ( SYNCTQ.ge.1.) then
       call MAPL_TerminateImport  ( GC,    &
          SHORT_NAME = (/'T'/),            &
          CHILD = RAD,                     &
          RC=STATUS)
       VERIFY_(STATUS)
     endif

    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)

! Set a profiling timer for GPU initialization
! --------------------------------------------
#ifdef _CUDA
    call MAPL_TimerAdd(GC, name="-GPUINIT"      ,RC=STATUS)
    VERIFY_(STATUS)
#endif

    call MAPL_TimerAdd(GC, name="RUN"           ,RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GenericSetServices ( GC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Initialize -- Initialize method for the composite Physics Gridded Component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Initialize method of the Physics Composite Gridded Component.
!  It acts as a driver for the initializtion of the five children: Radiation,
!  Turbulence, Moist, Chem, and Surface. It also sets up the frieldly connections
!  between the children and their sibling Turbulence, as well as with their
!  ``uncles'' Advection and Analysis.
!
!   For the turbulence tracer bundle, U and V come from
!   the import state, S is computed here from T and Z and kept
!   in the export state, the rest are friendlies from MOIST and CHEM.
!
!   The turbulence default behavior is a friendly with a zero flux
!   lower boundary condition and not producing a tendency.
!   Default tracers are put at the end of the bundles with a single
!   call; all others have to be done manually.
!
!   Any of the children`s exports that are friendly to advection or analysis
!   are put in the respective bundles by a single MAPL_Generic call. Remember
!   that friendly exports are were automatically allocated by the children
!   during the initialization sequence of the entire tree below Physics, which
!   is the first thing done here.

!   The increment tracer bundles for Moist and Turbulence are created with empty fields
!   except for those tracers which have explicit tendency Exports.
!
!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)           :: IAm
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases

   type (MAPL_MetaComp),       pointer :: STATE
   type (ESMF_GridComp),       pointer :: GCS(:)
   type (ESMF_State),          pointer :: GIM(:)
   type (ESMF_State),          pointer :: GEX(:)
   type (ESMF_FieldBundle)             :: BUNDLE, iBUNDLE
   type (ESMF_Field)                   :: FIELD, TempField
   type (ESMF_Grid)                    :: GRID

   integer                             :: NUM_TRACERS
   integer                             :: I
   integer                             :: NA
   character(len=ESMF_MAXSTR), pointer :: NAMES(:)
   character(len=ESMF_MAXSTR)          :: myNAME
   character(len=ESMF_MAXSTR)          ::  iNAME
   character(len=ESMF_MAXSTR)          :: fieldname

! Variables needed for GPU initialization

#ifdef _CUDA
   type (ESMF_VM)                      :: VM
   integer                             :: MYID  ! MPI Rank
   integer                             :: NCPUS ! MPI Size
   integer                             :: num_devices
   integer                             :: devicenum, inum
#endif

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

!define PRINT_STATES

    Iam = "Initialize"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, GRID=GRID, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "::" // Iam

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

#ifdef _CUDA

! Initialize the GPUs
! -------------------

!! Step One: Get our rank and size
!! -------------------------------

    call ESMF_VmGetCurrent(VM, rc=status)
    VERIFY_(STATUS)

    call ESMF_VmGet(VM, localPet=MYID, petCount=NCPUS, rc=STATUS)
    VERIFY_(STATUS)

!! Step Two: Initialize the GPUs
!! -----------------------------

    call MAPL_TimerOn(STATE,"TOTAL")
    call MAPL_TimerOn(STATE,"INITIALIZE")

    call MAPL_TimerOn(STATE,"-GPUINIT")

    STATUS = cudaGetDeviceCount(num_devices)
    if (STATUS /= 0) then
       write (*,*) "cudaGetDeviceCount failed: ", cudaGetErrorString(STATUS)
       _ASSERT(.FALSE.,'needs informative message')
    end if

    devicenum = mod(MYID, num_devices)

    STATUS = cudaSetDevice(devicenum)
    if (STATUS /= 0) then
       write (*,*) "cudaSetDevice failed: ", cudaGetErrorString(STATUS)
       _ASSERT(.FALSE.,'needs informative message')
    end if

    STATUS = cudaDeviceSetCacheConfig(cudaFuncCachePreferL1)
    if (STATUS /= 0) then
       write (*,*) "cudaDeviceSetCacheConfig failed: ", cudaGetErrorString(STATUS)
       _ASSERT(.FALSE.,'needs informative message')
    end if

    call MAPL_TimerOff(STATE,"-GPUINIT")

    call MAPL_TimerOff(STATE,"INITIALIZE")
    call MAPL_TimerOff(STATE,"TOTAL")

#endif

! Call Initialize for every Child
!--------------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(STATE,"TOTAL")
    call MAPL_TimerOn(STATE,"INITIALIZE")

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": IMPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( IMPORT, rc=STATUS )
    call WRITE_PARALLEL ( trim(Iam)//": EXPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( EXPORT, rc=STATUS )
#endif


! Get children and their im/ex states from my generic state.
!----------------------------------------------------------

    call MAPL_Get ( STATE, GCS=GCS, GIM=GIM, GEX=GEX, RC=STATUS )
    VERIFY_(STATUS)


!   Fill the turbulence tracer bundle: S, U, and V come from
!   the import state, the rest are friendlies from MOIST and CHEM.
!   For now, only S, U, V, and QV are non-default. Default tracers go last
!   in the bundle. This will have to be done better later by using
!   a default attribute.
!
!   The turbulence default behavior is a friendly with a zero flux
!   boundary condition, Default tracers do not expect a surface values
!   and do not produce produce a tendency or other products.
!
!   Default tracers are put at the end of the bundles with a single
!   call; all others have to be done manually.
! -----------------------------------------------------------------

    call ESMF_StateGet   (GIM(TURBL),  'TR' , BUNDLE,                     RC=STATUS )
    VERIFY_(STATUS)

! Add Non-Friendlies from Dynamics

    call ESMF_StateGet    (IMPORT,     'S'   , FIELD,                      RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_AttributeSet(FIELD, NAME="DiffuseLike"     ,VALUE="S",       RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_AttributeSet(FIELD, NAME="WeightedTendency",VALUE=.TRUE., RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_FieldBundleAdd   (BUNDLE,   FIELD,                                RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_StateGet    (IMPORT,     'U'   , FIELD,                      RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_AttributeSet(FIELD, NAME="DiffuseLike"     ,VALUE="U",       RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_FieldBundleAdd   (BUNDLE,   FIELD,                                RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_StateGet    (IMPORT,     'V'   , FIELD,                      RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_AttributeSet(FIELD, NAME="DiffuseLike"     ,VALUE="U",       RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_FieldBundleAdd   (BUNDLE,   FIELD,                                RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_StateGet    (GEX(TURBL),  'TKESHOC'   , FIELD,    RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_AttributeSet(FIELD, NAME="DiffuseLike"     ,VALUE="Q",       RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_FieldBundleAdd   (BUNDLE,   FIELD,                       RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_StateGet    (GEX(TURBL),  'QT2'   , FIELD,    RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_AttributeSet(FIELD, NAME="DiffuseLike"     ,VALUE="Q",       RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_FieldBundleAdd   (BUNDLE,   FIELD,                       RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_StateGet    (GEX(TURBL),  'QT3'   , FIELD,    RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_AttributeSet(FIELD, NAME="DiffuseLike"     ,VALUE="Q",       RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_FieldBundleAdd   (BUNDLE,   FIELD,                       RC=STATUS )
    VERIFY_(STATUS)

! Add Friendlies from Physics
    call MAPL_GridCompGetFriendlies(GC, "TURBULENCE", BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
! Add Friendlies from Moist (We assume QV is among these, all others are treated as default)
    call MAPL_GridCompGetFriendlies(GCS(MOIST) , "TURBULENCE", BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
! Add Friendlies from Chem (These are default tracers--zero surface flux)
    call MAPL_GridCompGetFriendlies(GCS(CHEM), "TURBULENCE", BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
! Print what is in the BUNDLE
    call ESMF_FieldBundleGet ( BUNDLE, fieldCount=NUM_TRACERS, RC=STATUS )
    VERIFY_(STATUS)
    allocate( NAMES(NUM_TRACERS),STAT=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet ( BUNDLE, fieldNameList=NAMES, rc=STATUS )
    VERIFY_(STATUS)
    do I=1,NUM_TRACERS
       call WRITE_PARALLEL ( trim(NAMES(I))//": in the TURBULENCE Bundle" )
    end do
    deallocate ( NAMES )

! Count tracers
!--------------

    call ESMF_FieldBundleGet(BUNDLE,FieldCount=NUM_TRACERS, RC=STATUS)
    VERIFY_(STATUS)

! Get the names of all tracers to fill other turbulence bundles.
!---------------------------------------------------------------

    allocate(NAMES(NUM_TRACERS),STAT=STATUS)
    VERIFY_(STATUS)

    call ESMF_FieldBundleGet(BUNDLE, fieldNameList=NAMES,  RC=STATUS)
    VERIFY_(STATUS)

! Fill the increments bundle that turbulence will export.
!  These fields come from the physics EXPORT, ones not there
!  will have an empty field.
!------------------------------------------------------------

    call ESMF_StateGet ( GEX(TURBL), 'TRI', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)

    do I=1,NUM_TRACERS
       select case (trim(NAMES(I)))
       case ('Q')
          iNAME = 'QVIT'
       case default
          iNAME = trim(NAMES(I)) // 'IT'
       end select
       call ESMFL_StateGetField  (EXPORT, (/iNAME/), &
            BUNDLE, (/trim(NAMES(I))//'IT'/), RC=STATUS )
       VERIFY_(STATUS)
    end do

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": Turbulence Increment Bundle" )
    if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, rc=STATUS )
#endif

! Fill the turbulence bundle of surface skin values
!--------------------------------------------------

    call ESMF_StateGet ( GIM(TURBL), 'TRG', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)

    do I=1,NUM_TRACERS
       iNAME = trim(NAMES(I)) // 'HAT'
       call ESMFL_StateGetField  (GEX(SURF), (/iNAME/), BUNDLE, RC=STATUS )
       VERIFY_(STATUS)
    end do

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": Turbulence Surface Bundle" )
    if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, rc=STATUS )
#endif

! Fill the turbulence bundle of changes of surface skin values
!-------------------------------------------------------------

    call ESMF_StateGet ( GIM(TURBL), 'DTG', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)

    do I=1,NUM_TRACERS
       select case (trim(NAMES(I)))
       case('S')
          iNAME = 'DELSS'
       case('U')
          iNAME = 'DELUS'
       case('V')
          iNAME = 'DELVS'
       case('Q')
          iNAME = 'DELQS'
       case default
          iNAME = NAMES(I)
       end select

       call ESMFL_StateGetField  (GEX(SURF), (/iNAME/), &
            BUNDLE, (/trim(NAMES(I))//'DEL'/), RC=STATUS )
       VERIFY_(STATUS)
    end do

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": Turbulence DTG Bundle" )
    if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, rc=STATUS )
#endif

! Fill the turbulence FSTAR bundle (surface fluxes)
!-----------------------------------

    call ESMF_StateGet ( GEX(TURBL), 'FSTAR', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)

    do I=1,NUM_TRACERS
       select case (trim(NAMES(I)))
       case('S')
          iNAME = 'SH'
       case('U')
          iNAME = 'TAUX'
       case('V')
          iNAME = 'TAUY'
       case('Q')
          iNAME = 'EVAP'
       case default
          iNAME = NAMES(I)
       end select

       call ESMFL_StateGetField  (GIM(SURF), (/iNAME/), &
            BUNDLE, (/trim(NAMES(I))//'FLX'/), RC=STATUS )
       VERIFY_(STATUS)
    end do

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": Surface Flux Bundle" )
    if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, rc=STATUS )
#endif

! Fill the turbulence DFSTAR bundle d(surface fluxes)/d(TG)
!----------------------------------------------------------

    call ESMF_StateGet ( GEX(TURBL), 'DFSTAR', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)

    do I=1,NUM_TRACERS
       select case (trim(NAMES(I)))
       case('S')
          iNAME = 'DSH'
       case('U')
          iNAME = 'DFU'
       case('V')
          iNAME = 'DFV'
       case('Q')
          iNAME = 'DEVAP'
       case default
          iNAME = NAMES(I)
       end select

       call ESMFL_StateGetField  (GIM(SURF), (/iNAME/), &
            BUNDLE, (/trim(NAMES(I))//'DFL'/), RC=STATUS )
       VERIFY_(STATUS)
    end do

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": Surface Flux Derivative Bundle" )
    if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( BUNDLE, rc=STATUS )
#endif

    deallocate(NAMES)

! Fill export bundle of child quantities to be advected
!------------------------------------------------------

    call ESMF_StateGet(EXPORT, 'TRADV', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
! Add Friendlies to DYNAMICS from all children of physics
    call MAPL_GridCompGetFriendlies(GC,  "DYNAMICS", BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
! Add Friendlies to DYNAMICS from all children of physics
    call MAPL_GridCompGetFriendlies(GCS, "DYNAMICS", BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
! Print what is in the BUNDLE
    call ESMF_FieldBundleGet ( BUNDLE, fieldCount=NUM_TRACERS, RC=STATUS )
    VERIFY_(STATUS)
    allocate( NAMES(NUM_TRACERS),STAT=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet ( BUNDLE, fieldNameList=NAMES, rc=STATUS )
    VERIFY_(STATUS)
    do I=1,NUM_TRACERS
       call WRITE_PARALLEL ( trim(NAMES(I))//": in the DYNAMICS Advection Bundle" )
    end do
    deallocate ( NAMES )

! Initialize Water rescale tendency bundle
!--------------------------------------------
    call Initialize_IncBundle_init(GC, EXPORT, EXPORT, H2Oinc, __RC__)

! Fill export bundle of child quantities to be analyzed
!------------------------------------------------------

    call ESMF_StateGet(EXPORT, 'TRANA', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
! Add Friendlies to MOIST from Physics
    call MAPL_GridCompGetFriendlies(GC, "ANALYSIS", BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
! Add Friendlies ta ANALYSIS from all children of physics
    call MAPL_GridCompGetFriendlies(GCS, "ANALYSIS", BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
! Print what is in the BUNDLE
    call ESMF_FieldBundleGet ( BUNDLE, fieldCount=NUM_TRACERS, RC=STATUS )
    VERIFY_(STATUS)
    allocate( NAMES(NUM_TRACERS),STAT=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet ( BUNDLE, fieldNameList=NAMES, rc=STATUS )
    VERIFY_(STATUS)
    do I=1,NUM_TRACERS
       call WRITE_PARALLEL ( trim(NAMES(I))//": in the ANALYSIS Bundle" )
    end do
    deallocate ( NAMES )

! Fill export bundle of child quantities to go thru CONVECTIVE transport
!  No need for tendencies at this point; scavenging may be controled by
!  field attributes (TBD)
!-----------------------------------------------------------------------

    call ESMF_StateGet       (GIM(MOIST), 'MTR', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
! Add Friendlies to MOIST from Physics
    call MAPL_GridCompGetFriendlies(GC, "MOIST", BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
! Add Friendlies to MOIST from all children of physics
    call MAPL_GridCompGetFriendlies(GCS, "MOIST", BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
! Print what is in the BUNDLE
    call ESMF_FieldBundleGet ( BUNDLE, fieldCount=NUM_TRACERS, RC=STATUS )
    VERIFY_(STATUS)
    allocate( NAMES(NUM_TRACERS),STAT=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet ( BUNDLE, fieldNameList=NAMES, rc=STATUS )
    VERIFY_(STATUS)
    do I=1,NUM_TRACERS
       call WRITE_PARALLEL ( trim(NAMES(I))//": in the MOIST Convective transport Bundle" )
    end do
    deallocate ( NAMES )

! Fill the moist increments bundle
!---------------------------------

    call Initialize_IncMBundle_init(GC, GIM(MOIST), EXPORT, __RC__)

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": Convective Transport Tendency Bundle" )
    if ( MAPL_am_I_root() ) call ESMF_FieldBundlePrint ( iBUNDLE, rc=STATUS )
#endif

    call MAPL_TimerOff(STATE,"INITIALIZE")
    call MAPL_TimerOff(STATE,"TOTAL")


! All Done
!---------

    RETURN_(ESMF_SUCCESS)
 end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Run -- Run method for the composite Physics Gridded Component

! !INTERFACE:

   subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The run method for the physics calls the children`s
!   run methods. It also prepares inputs and couplings amongst them.
!   Its main outputs are the combined tendencies needed by the dynamics.

!EOP

! ErrLog Variables

   character(len=ESMF_MAXSTR)          :: IAm
   integer                             :: STATUS
   character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

   type (MAPL_MetaComp),      pointer  :: STATE
   type (MAPL_MetaComp),      pointer  :: CMETA
   type (ESMF_GridComp),      pointer  :: GCS(:)
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_State),         pointer  :: GEX(:)
   type (ESMF_State)                   :: INTERNAL
   type (ESMF_Grid)                    :: grid
   type( ESMF_VM )                     :: VMG

   type (ESMF_Config)                  :: CF
   type (ESMF_FieldBundle)             :: BUNDLE
   character(len=ESMF_MAXSTR),pointer  :: GCNames(:)
   character(len=ESMF_MAXSTR)          :: DUMMY
   integer                             :: I, J, L, K, N
   integer                             :: IM, JM, LM, NQ
   integer                             :: ISPPT,ISKEB
   logical                             :: DO_SPPT,DO_SKEB
   logical                             :: NEED_TOT
   logical                             :: NEED_FRI
   logical                             :: NEED_STN
   logical                             :: DPEDT_PHYS
   real                                :: DT
   real                                :: SYNCTQ, DOPHYSICS
   real                                :: HGT_SURFACE

   real, pointer, dimension(:,:,:)     :: S, T, ZLE, PLE, PK, U, V, W
   real, pointer, dimension(:,:,:)     :: DM, DPI, TOT, FRI, TTN, STN,TMP
   real, pointer, dimension(:,:,:)     :: QW, QV, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL
   real, pointer, dimension(:,:,:)     :: ptr3d

   real, pointer, dimension(:,:,:)     :: DUDT
   real, pointer, dimension(:,:,:)     :: DVDT
   real, pointer, dimension(:,:,:)     :: DWDT
   real, pointer, dimension(:,:,:)     :: DTDT
   real, pointer, dimension(:,:,:)     :: DTDTTOT
   real, pointer, dimension(:,:,:)     :: DTDTRAD
   real, pointer, dimension(:,:,:)     :: DPDTPHY
   real, pointer, dimension(:,:,:)     :: DPDT
   real, pointer, dimension(:,:)       :: DMDT
   real, pointer, dimension(:    )     :: PREF

   real, pointer, dimension(:,:,:)     :: DOXDTCHM
   real, pointer, dimension(:,:,:)     :: DQVDTMST, DQVDTTRB, DQVDTCHM
   real, pointer, dimension(:,:,:)     :: DQLDTTRB, DQIDTTRB
   real, pointer, dimension(:,:,:)     :: DQLDTSCL, DQIDTSCL, DQVDTSCL
   real, pointer, dimension(:,:,:)     :: DQLDTMST, DQIDTMST
   real, pointer, dimension(:,:,:)     :: DPDTMST,  DPDTTRB

   real, pointer, dimension(:,:,:)     :: RNDPERT,RNDPTR
   real, pointer, dimension(:,:,:)     :: SKEBU_WT,SKEBV_WT
   real, pointer, dimension(:,:,:)     :: SKEBU,SKEBV
   real, pointer, dimension(:,:,:)     :: DUDTSTOCH, DVDTSTOCH, DTDTSTOCH, DQVDTSTOCH
   real, pointer, dimension(:,:,:)     :: DQVDTPERT
   real, pointer, dimension(:,:,:)     :: QVPERT, QVUNPERT
   real, pointer, dimension(:,:,:)     :: DTDTUNPERT, DUDTUNPERT, DVDTUNPERT, DQVDTUNPERT
   real, pointer, dimension(:,:,:)     :: TIR, TIM, TIMFRIC, TIT, TIF
   real, pointer, dimension(:,:,:)     :: UIM, VIM, WIM
   real, pointer, dimension(:,:,:)     :: UIT, VIT, SIT
   real, pointer, dimension(:,:,:)     :: UIG, VIG, TIG, TICU
   real, pointer, dimension(:,:,:)     :: FTU, FTV
   real, pointer, dimension(:,:,:)     :: INTDIS, TOPDIS
   real, pointer, dimension(:,:  )     :: SRFDIS

   real, pointer, dimension(:,:  )     :: DQVDTPHYINT, DQLDTPHYINT, DQIDTPHYINT, DOXDTPHYINT
   real, pointer, dimension(:,:  )     :: DQVDTTRBINT, DQVDTMSTINT, DQVDTCHMINT
   real, pointer, dimension(:,:  )     :: DQLDTMSTINT, DQIDTMSTINT, DOXDTCHMINT
   real, pointer, dimension(:,:  )     :: PERAD,PETRB,PEMST,PEFRI,PEGWD,PECUF
   real, pointer, dimension(:,:  )     :: PEPHY
   real, pointer, dimension(:,:  )     :: KEPHY
   real, pointer, dimension(:,:  )     :: AREA

   real*8, allocatable, dimension(:,:)   :: sumq
   real*8, allocatable, dimension(:,:,:) :: ple_new

   character(len=ESMF_MAXSTR), allocatable  :: NAMES(:)



! SYNCTQ & UV pointers
   real, pointer, dimension(:,:,:)     :: UAFMOIST, VAFMOIST,  TAFMOIST, QAFMOIST, THAFMOIST, SAFMOIST
   real, pointer, dimension(:,:)       ::  UFORSURF, VFORSURF, TFORSURF, QFORSURF, SPD4SURF
   real, pointer, dimension(:,:,:)     ::  UFORCHEM, VFORCHEM, TFORCHEM, THFORCHEM
   real, pointer, dimension(:,:,:)     ::  UFORTURB, VFORTURB, TFORTURB, THFORTURB, SFORTURB
   real, pointer, dimension(:,:,:)     ::                      TFORRAD
   real, pointer, dimension(:,:,:)     :: UAFDIFFUSE, VAFDIFFUSE, QAFDIFFUSE, SAFDIFFUSE, SAFUPDATE

   real, allocatable, dimension(:,:,:) :: HGT
   real, allocatable, dimension(:,:,:) :: TDPOLD, TDPNEW
   real, allocatable, dimension(:,:,:) :: TFORQS
   real, allocatable, dimension(:,:)   :: qs,pmean

   logical :: isPresent, SCM_NO_RAD
   real, allocatable, target :: zero(:,:,:)

   real(kind=MAPL_R8), allocatable, dimension(:,:) :: sumdq
   real(kind=MAPL_R8), allocatable, dimension(:,:) ::  dpe
   real(kind=MAPL_R8), allocatable, dimension(:,:,:) :: dq

   real, pointer, dimension(:,:,:)     :: DTDT_BL, DQDT_BL

   real*8, allocatable, dimension(:,:)   :: sum_qdp_b4
   real*8, allocatable, dimension(:,:)   :: sum_qdp_af
   real,   allocatable, dimension(:,:,:) ::     qdp_b4
   real,   allocatable, dimension(:,:,:) ::     qdp_af
   real,   allocatable, dimension(:,:,:) ::    ple_mst
   real,   allocatable, dimension(:,:,:) ::     dp_mst
   real,   allocatable, dimension(:,:)   :: qdp_b4_int
   real,   allocatable, dimension(:,:)   :: qdp_af_int

   real*8                                ::       gamma
   real*8                                ::  qdp_b4_ave
   real*8                                ::  qdp_af_ave

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, VM=VMG, config=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "::" // Iam

    call ESMF_GridCompGet ( GC, grid=grid, rc=status )
    VERIFY_(STATUS)

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(STATE,"TOTAL")
    call MAPL_TimerOn(STATE,"RUN")

    call MAPL_GetResource(STATE, SCM_NO_RAD, Label="SCM_NO_RAD:", default=.FALSE., RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(STATE, DUMMY, Label="DPEDT_PHYS:", default='YES', RC=STATUS)
    VERIFY_(STATUS)
         DUMMY = ESMF_UtilStringUpperCase(DUMMY)
    DPEDT_PHYS = TRIM(DUMMY).eq.'YES'

! Get the children`s states from the generic state
!-------------------------------------------------

    call MAPL_Get ( STATE,   &
        GCS=GCS, GIM=GIM, GEX=GEX,       &
        IM = IM, JM = JM, LM = LM,       &
        GCNames = GCNames,               &
        INTERNAL_ESMF_STATE = INTERNAL,  &
                               RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(CF, DT, Label="RUN_DT:" , RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(CF, ISPPT, Label="SPPT:", DEFAULT = 0 , RC=STATUS)
    VERIFY_(STATUS)
    DO_SPPT = (ISPPT/=0)
    call ESMF_ConfigGetAttribute(CF, ISKEB, Label="SKEB:", DEFAULT = 0 , RC=STATUS)
    VERIFY_(STATUS)
    DO_SKEB = (ISKEB/=0)

    call ESMF_StateGet (EXPORT, 'TRADV', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)

    allocate(zero(IM,JM,LM),stat=status)
    VERIFY_(status)
    zero = 0.0

    call ESMFL_BundleGetPointertoData( BUNDLE,'Q'   ,QV  , RC=STATUS)
    VERIFY_(STATUS)
    call ESMFL_BundleGetPointertoData( BUNDLE,'QLLS',QLLS, RC=STATUS)
    VERIFY_(STATUS)
    call ESMFL_BundleGetPointertoData( BUNDLE,'QLCN',QLCN, RC=STATUS)
    VERIFY_(STATUS)
    call ESMFL_BundleGetPointertoData( BUNDLE,'QILS',QILS, RC=STATUS)
    VERIFY_(STATUS)
    call ESMFL_BundleGetPointertoData( BUNDLE,'QICN',QICN, RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet ( BUNDLE, fieldName='QRAIN', isPresent=isPresent, RC=STATUS)
    if (isPresent) then
       call ESMFL_BundleGetPointerToData( BUNDLE, 'QRAIN', QRAIN, RC=STATUS )
       VERIFY_(STATUS)
    else
       QRAIN => zero
    end if
    call ESMF_FieldBundleGet ( BUNDLE, fieldName='QSNOW', isPresent=isPresent, RC=STATUS)
    if (isPresent) then
       call ESMFL_BundleGetPointerToData( BUNDLE, 'QSNOW', QSNOW, RC=STATUS )
       VERIFY_(STATUS)
    else
       QSNOW => zero
    end if
    call ESMF_FieldBundleGet ( BUNDLE, fieldName='QGRAUPEL', isPresent=isPresent, RC=STATUS)
    if (isPresent) then
       call ESMFL_BundleGetPointerToData( BUNDLE, 'QGRAUPEL', QGRAUPEL, RC=STATUS )
       VERIFY_(STATUS)
    else
       QGRAUPEL => zero
    end if

! Initialize Passive Tracer QW
! ----------------------------
    call MAPL_GetPointer(INTERNAL, QW, 'QW', RC=STATUS); VERIFY_(STATUS)
    QW = QV+QLLS+QLCN+QILS+QICN+QRAIN+QSNOW+QGRAUPEL

! Get Global PHYSICS Parameters
! -----------------------------
    call MAPL_GetResource(STATE, SYNCTQ,    'SYNCTQ:',    DEFAULT= 1.0, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(STATE, DOPHYSICS, 'DOPHYSICS:', DEFAULT= 1.0, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(STATE, HGT_SURFACE, Label="HGT_SURFACE:", DEFAULT= 50.0, RC=STATUS)
    VERIFY_(STATUS)


! Pointers to Imports
!--------------------

    call MAPL_GetPointer(IMPORT,  U,       'U'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,  V,       'V'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,  W,       'W'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,  T,       'T'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,  S,       'S'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,  ZLE,     'ZLE'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,  PLE,     'PLE'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,  AREA,    'AREA'   , RC=STATUS); VERIFY_(STATUS)

    allocate( TDPOLD(IM,JM,LM),stat=STATUS )
    VERIFY_(STATUS)

    TDPOLD = T(:,:,1:LM) * (PLE(:,:,1:LM)-PLE(:,:,0:LM-1))

    allocate(DM(IM,JM,LM),stat=STATUS)
    VERIFY_(STATUS)
    DM = (PLE(:,:,1:LM)-PLE(:,:,0:LM-1))*(1.0/MAPL_GRAV)

    allocate(DPI(IM,JM,LM),stat=STATUS)
    VERIFY_(STATUS)
    DPI = 1./(PLE(:,:,1:LM)-PLE(:,:,0:LM-1))

   ! Create Old Dry Mass Variables
   ! -----------------------------
     allocate(   sumq( IM,JM ),    STAT=STATUS ) ; VERIFY_(STATUS)
     allocate( ple_new(IM,JM,0:LM),STAT=STATUS ) ; VERIFY_(STATUS)

! Pointers to Exports
!--------------------

    call MAPL_GetPointer(EXPORT, DUDT,     'DUDT'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DVDT,     'DVDT'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DWDT,     'DWDT'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DTDT,     'DTDT'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DTDTTOT,  'DTDTTOT' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DTDTRAD,  'DTDTRAD' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DPDTPHY,  'DPDTPHY' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DPDT,     'DPEDT'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DMDT,     'DMDT'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TIT,      'TIT'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TIM,      'TIM'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TIMFRIC,  'TIMFRIC' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TIF,      'TIF'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FTU,      'FTU'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FTV,      'FTV'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, KEPHY,    'KEPHY'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PEPHY,    'PEPHY'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PERAD,    'PERAD'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PETRB,    'PETRB'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PEMST,    'PEMST'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PEFRI,    'PEFRI'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PEGWD,    'PEGWD'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PECUF,    'PECUF'   , RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT, DQVDTMSTINT, 'DQVDTMSTINT', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQLDTMSTINT, 'DQLDTMSTINT', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQIDTMSTINT, 'DQIDTMSTINT', RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT, DQVDTTRBINT, 'DQVDTTRBINT', RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT, DQVDTCHMINT, 'DQVDTCHMINT', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DOXDTCHMINT, 'DOXDTCHMINT', RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT, DQVDTPHYINT, 'DQVDTPHYINT', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQLDTPHYINT, 'DQLDTPHYINT', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQIDTPHYINT, 'DQIDTPHYINT', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DOXDTPHYINT, 'DOXDTPHYINT', RC=STATUS); VERIFY_(STATUS)

! Get and allocate pointers to Exports that have been put in turbulence
!   bundle, as well as the required tendencies in the children`s exports.
!------------------------------------------------------------------------

    if(associated(DUDT)) then
       call MAPL_GetPointer(EXPORT     ,   UIT,     'UIT', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GEX(MOIST) ,   UIM,    'DUDT', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GEX(GWD)   ,   UIG,    'DUDT', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(DVDT)) then
       call MAPL_GetPointer(EXPORT     ,   VIT,     'VIT', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GEX(MOIST) ,   VIM,    'DVDT', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GEX(GWD)   ,   VIG,    'DVDT', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(DWDT)) then
       call MAPL_GetPointer(GEX(MOIST) ,   WIM,    'DWDT', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(DTDT) .or. associated(TIM) .or. associated(DTDTTOT)) then
       call MAPL_GetPointer(GEX(MOIST) ,   TTN,   'DTDT', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(DTDT) .or. associated(TIT) .or. associated(DTDTTOT)) then
       call MAPL_GetPointer(EXPORT     ,   SIT,     'SIT', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(DTDT) .or. associated(DTDTRAD) .or. associated(DTDTTOT)) then
       call MAPL_GetPointer(GEX(RAD ) ,    TIR,    'DTDT', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(DTDT) .or. associated(TIF) .or. associated(DTDTTOT)) then
       call MAPL_GetPointer(GEX(TURBL), INTDIS,  'INTDIS', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GEX(TURBL), TOPDIS,  'TOPDIS', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GEX(TURBL), SRFDIS,  'SRFDIS', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GEX(GWD ) ,    TIG,    'DTDT', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GEX(MOIST),   TICU,'DTDTFRIC', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

       call MAPL_GetPointer ( EXPORT,     DQVDTTRB, 'QVIT',     alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer ( EXPORT,     DQLDTTRB, 'QLLSIT',   alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer ( EXPORT,     DQIDTTRB, 'QILSIT',   alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer ( GEX(TURBL), DPDTTRB , 'DPDTTRB',  alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer ( GEX(MOIST), DPDTMST , 'DPDTMST',  alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer ( GEX(MOIST), DQVDTMST, 'DQDT',     alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer ( GEX(CHEM),  DQVDTCHM, 'H2O_TEND', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer ( GEX(MOIST), DQLDTMST, 'DQLDT',    alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer ( GEX(MOIST), DQIDTMST, 'DQIDT',    alloc=.true., RC=STATUS)
       VERIFY_(STATUS)

       call MAPL_GetPointer (EXPORT, DQVDTSCL, 'DQVDTSCL', RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, DQLDTSCL, 'DQLDTSCL', RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, DQIDTSCL, 'DQIDTSCL', RC=STATUS)
       VERIFY_(STATUS)
    if (DO_SPPT) then
       call MAPL_GetPointer (EXPORT, RNDPTR,    'RNDPTR',    RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, DTDTUNPERT, 'DTDTUNPERT',   RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, DUDTUNPERT, 'DUDTUNPERT',  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, DVDTUNPERT, 'DVDTUNPERT',  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, DQVDTUNPERT, 'DQVDTUNPERT',   RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, DQVDTPERT, 'DQVDTPERT',  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, QVUNPERT,  'QVUNPERT',   RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, QVPERT,    'QVPERT',     RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, DUDTSTOCH, 'DUDTSTOCH', RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, DVDTSTOCH, 'DVDTSTOCH', RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, DTDTSTOCH, 'DTDTSTOCH', RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, DQVDTSTOCH, 'DQVDTSTOCH', RC=STATUS)
       VERIFY_(STATUS)
    endif
    if (DO_SKEB) then
       call MAPL_GetPointer (EXPORT, SKEBU,    'SKEBU',    alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer (EXPORT, SKEBV,    'SKEBV',    alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(DOXDTCHMINT)) then
       call MAPL_GetPointer ( GEX(CHEM),  DOXDTCHM,  'OX_TEND', alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

!----------------------

    if ( DOPHYSICS.eq.1. ) then

     if ( SYNCTQ.ge.1. ) then
      !  Will need PK to get from T to TH and back
      allocate(PK(IM,JM,LM),stat=STATUS);VERIFY_(STATUS)
      PK = ((0.5*(PLE(:,:,0:LM-1)+PLE(:,:,1:LM))) / MAPL_P00)**MAPL_KAPPA
      if ( (LM .ne. 72) .and. (HGT_SURFACE .gt. 0.0) ) then
         allocate(HGT(IM,JM,LM+1),stat=STATUS);VERIFY_(STATUS)
         do k = 1,LM+1
           HGT(:,:,k) = (ZLE(:,:,k-1) - ZLE(:,:,LM))
         enddo
      endif
     endif

! Gravity Wave Drag
!  (must be first to use Q from dynamics,
!    as it was when it was in superdyn.)
!----------------------------------------

    I=GWD

    call MAPL_TimerOn (STATE,GCNames(I))
     call ESMF_GridCompRun (GCS(I), importState=GIM(I), exportState=GEX(I), clock=CLOCK, userRC=STATUS ); VERIFY_(STATUS)
     call MAPL_GenericRunCouplers (STATE, I,        CLOCK,    RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,GCNames(I))

! Moist Processes
!----------------

    call Initialize_IncMBundle_run(GIM(MOIST), EXPORT, DM=DM,__RC__)

    I=MOIST

    call MAPL_TimerOn (STATE,GCNames(I))
     call ESMF_GridCompRun (GCS(I), importState=GIM(I), exportState=GEX(I), clock=CLOCK, userRC=STATUS ); VERIFY_(STATUS)
     call MAPL_GenericRunCouplers (STATE, I,        CLOCK,    RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,GCNames(I))

    call MAPL_GetObjectFromGC ( GCS(I), CMETA, _RC)

    call Compute_IncMBundle(GIM(MOIST), EXPORT, CMETA, DM=DM,__RC__)

    call MAPL_GetPointer(GIM(MOIST), DTDT_BL, 'DTDT_BL', alloc = .true. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(GIM(MOIST), DQDT_BL, 'DQDT_BL', alloc = .true. ,RC=STATUS); VERIFY_(STATUS)

!  SYNCTQ - Stage 1 SYNC of T/Q and U/V
!--------------------------------------
    if ( SYNCTQ.ge.1. ) then
    ! From Moist
     call MAPL_GetPointer ( GEX(MOIST),  UAFMOIST,  'UAFMOIST', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GEX(MOIST),  VAFMOIST,  'VAFMOIST', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GEX(MOIST),  TAFMOIST,  'TAFMOIST', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GEX(MOIST), THAFMOIST, 'THAFMOIST', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GEX(MOIST),  SAFMOIST,  'SAFMOIST', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GEX(MOIST),  QAFMOIST,  'QAFMOIST', RC=STATUS); VERIFY_(STATUS)
    ! Boundary Layer Tendencies for GF
     DTDT_BL=TAFMOIST
     DQDT_BL=QV
    ! For SURF
     call MAPL_GetPointer ( GIM(SURF),  UFORSURF,  'UA',    RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(SURF),  VFORSURF,  'VA',    RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(SURF),  TFORSURF,  'TA',    RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(SURF),  QFORSURF,  'QA',    RC=STATUS); VERIFY_(STATUS)
     if ( (LM .ne. 72) .and. (HGT_SURFACE .gt. 0.0) ) then
       call VertInterp(UFORSURF,UAFMOIST,-HGT,-HGT_SURFACE, status); VERIFY_(STATUS)
       call VertInterp(VFORSURF,VAFMOIST,-HGT,-HGT_SURFACE, status); VERIFY_(STATUS)
       call VertInterp(TFORSURF,TAFMOIST,-HGT,-HGT_SURFACE, status); VERIFY_(STATUS)
       call VertInterp(QFORSURF,QAFMOIST,-HGT,-HGT_SURFACE, status); VERIFY_(STATUS)
     else
       UFORSURF = UAFMOIST(:,:,LM)
       VFORSURF = VAFMOIST(:,:,LM)
       TFORSURF = TAFMOIST(:,:,LM)
       QFORSURF = QAFMOIST(:,:,LM)
     endif
     call MAPL_GetPointer ( GIM(SURF),  SPD4SURF,  'SPEED', RC=STATUS); VERIFY_(STATUS)
     SPD4SURF = SQRT( UFORSURF*UFORSURF + VFORSURF*VFORSURF )
    ! For CHEM
     if ( SYNCTQ.eq.1. ) then
       call MAPL_GetPointer ( GIM(CHEM), TFORCHEM,   'T',  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer ( GIM(CHEM), THFORCHEM,  'TH', RC=STATUS); VERIFY_(STATUS)
        TFORCHEM =  TAFMOIST
       THFORCHEM = THAFMOIST
     endif
    ! For TURBL
     call ESMF_StateGet(GIM(TURBL), 'TR', BUNDLE, RC=STATUS ); VERIFY_(STATUS)
     call ESMFL_BundleGetPointerToData(BUNDLE,'S',SFORTURB, RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(TURBL), UFORTURB,   'U', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(TURBL), VFORTURB,   'V', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(TURBL), TFORTURB,   'T', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(TURBL), THFORTURB, 'TH', RC=STATUS); VERIFY_(STATUS)
      UFORTURB =  UAFMOIST
      VFORTURB =  VAFMOIST
      TFORTURB =  TAFMOIST
     THFORTURB = THAFMOIST
      SFORTURB =  SAFMOIST
    endif

! Surface Stage 1
!----------------

    I=SURF

    call MAPL_TimerOn(STATE,GCNames(I))
     call ESMF_GridCompRun (GCS(I), importState=GIM(I), exportState=GEX(I), clock=CLOCK, PHASE=1, userRC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,GCNames(I))

! Aerosol/Chemistry Stage 1
!--------------------------

    I=CHEM

    call MAPL_TimerOn(STATE,GCNames(I))
     call ESMF_GridCompRun (GCS(I), importState=GIM(I), exportState=GEX(I), clock=CLOCK, phase=1, userRC=STATUS ); VERIFY_(STATUS)
     call MAPL_GenericRunCouplers (STATE, I,        CLOCK,    RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,GCNames(I))

! Turbulence Stage 1
!-------------------

    I=TURBL

    call MAPL_TimerOn(STATE,GCNames(I))
     call ESMF_GridCompRun (GCS(I), importState=GIM(I), exportState=GEX(I), clock=CLOCK, PHASE=1, userRC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,GCNames(I))

!  SYNCTQ - Stage 2 SYNC of T/Q and U/V
!--------------------------------------
    if ( SYNCTQ.ge.1. ) then
    ! From TURBL Run 1
     call MAPL_GetPointer ( GEX(TURBL), UAFDIFFUSE, 'UAFDIFFUSE', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GEX(TURBL), VAFDIFFUSE, 'VAFDIFFUSE', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GEX(TURBL), SAFDIFFUSE, 'SAFDIFFUSE', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GEX(TURBL), QAFDIFFUSE, 'QAFDIFFUSE', RC=STATUS); VERIFY_(STATUS)
    ! For TURBL
     call ESMF_StateGet(GIM(TURBL), 'TR', BUNDLE, RC=STATUS ); VERIFY_(STATUS)
     call ESMFL_BundleGetPointerToData(BUNDLE,'S',SFORTURB, RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(TURBL), UFORTURB,   'U', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(TURBL), VFORTURB,   'V', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(TURBL), TFORTURB,   'T', RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(TURBL), THFORTURB, 'TH', RC=STATUS); VERIFY_(STATUS)
      UFORTURB = UAFDIFFUSE
      VFORTURB = VAFDIFFUSE
     ! For Stage 2 - Changes in S from TURBL assumed to be all in T
      TFORTURB = TFORTURB + (SAFDIFFUSE-SFORTURB)/MAPL_CP
     THFORTURB = TFORTURB/PK
      SFORTURB = SAFDIFFUSE
    ! For SURF
     call MAPL_GetPointer ( GIM(SURF),  UFORSURF,  'UA',    RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(SURF),  VFORSURF,  'VA',    RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(SURF),  TFORSURF,  'TA',    RC=STATUS); VERIFY_(STATUS)
     call MAPL_GetPointer ( GIM(SURF),  QFORSURF,  'QA',    RC=STATUS); VERIFY_(STATUS)
     if ( (LM .ne. 72) .and. (HGT_SURFACE .gt. 0.0) ) then
       call VertInterp(TFORSURF,TFORTURB,-HGT,-HGT_SURFACE, status); VERIFY_(STATUS)
       call VertInterp(UFORSURF,UAFDIFFUSE,-HGT,-HGT_SURFACE, status); VERIFY_(STATUS)
       call VertInterp(VFORSURF,VAFDIFFUSE,-HGT,-HGT_SURFACE, status); VERIFY_(STATUS)
       call VertInterp(QFORSURF,QAFDIFFUSE,-HGT,-HGT_SURFACE, status); VERIFY_(STATUS)
     else
       TFORSURF =   TFORTURB(:,:,LM)
       UFORSURF = UAFDIFFUSE(:,:,LM)
       VFORSURF = VAFDIFFUSE(:,:,LM)
       QFORSURF = QAFDIFFUSE(:,:,LM)
     endif
     call MAPL_GetPointer ( GIM(SURF),  SPD4SURF,  'SPEED', RC=STATUS); VERIFY_(STATUS)
     SPD4SURF = SQRT( UFORSURF*UFORSURF + VFORSURF*VFORSURF )
    endif

! Surface Stage 2
!----------------

    I=SURF

    call MAPL_TimerOn(STATE,GCNames(I))
     call ESMF_GridCompRun (GCS(I), importState=GIM(I), exportState=GEX(I), clock=CLOCK, PHASE=2, userRC=STATUS ); VERIFY_(STATUS)
     call MAPL_GenericRunCouplers (STATE, I,        CLOCK,    RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,GCNames(I))

! Turbulence Stage 2
!-------------------

    I=TURBL

    call MAPL_TimerOn(STATE,GCNames(I))
     call ESMF_GridCompRun (GCS(I), importState=GIM(I), exportState=GEX(I), clock=CLOCK, PHASE=2, userRC=STATUS ); VERIFY_(STATUS)
     call MAPL_GenericRunCouplers (STATE, I,        CLOCK,    RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,GCNames(I))

    if ( SYNCTQ.ge.1. ) then
     ! From TURBL Stage 2
     call MAPL_GetPointer ( GEX(TURBL), SAFUPDATE,  'SAFUPDATE', RC=STATUS); VERIFY_(STATUS)
     ! For RAD
     call MAPL_GetPointer ( GIM(RAD), TFORRAD, 'T', RC=STATUS); VERIFY_(STATUS)
     ! For Stage 2 - Changes in S from TURBL assumed to be all in T
      TFORRAD = TFORTURB + (SAFUPDATE-SAFDIFFUSE)/MAPL_CP
     ! For CHEM use the same T as CHEM
     if ( SYNCTQ.eq.1. ) then
       call MAPL_GetPointer ( GIM(CHEM), TFORCHEM,   'T',  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer ( GIM(CHEM), THFORCHEM,  'TH', RC=STATUS); VERIFY_(STATUS)
        TFORCHEM = TFORRAD
       THFORCHEM = TFORRAD/PK
     endif
    endif

! Boundary Layer Tendencies for GF
!--------------------------
    if ( SYNCTQ.ge.1. ) then
       DTDT_BL=(TFORRAD-DTDT_BL)/DT
       DQDT_BL=(QV-DQDT_BL)/DT
    endif

! Aerosol/Chemistry Stage 2
!--------------------------

    I=CHEM

    call MAPL_TimerOn(STATE,GCNames(I))
     call ESMF_GridCompRun (GCS(I), importState=GIM(I), exportState=GEX(I), clock=CLOCK, PHASE=2, userRC=STATUS ); VERIFY_(STATUS)
     call MAPL_GenericRunCouplers (STATE, I,        CLOCK,    RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,GCNames(I))

! Radiation
!----------

    I=RAD

    call MAPL_TimerOn (STATE,GCNames(I))
     call ESMF_GridCompRun (GCS(I), importState=GIM(I), exportState=GEX(I), clock=CLOCK, userRC=STATUS ); VERIFY_(STATUS)
     call MAPL_GenericRunCouplers (STATE, I,        CLOCK,    RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOff(STATE,GCNames(I))

! Clean up SYNTQ things
    if ( SYNCTQ.ge.1. ) then
      deallocate(PK)
      if ( (LM .ne. 72) .and. (HGT_SURFACE .gt. 0.0) ) deallocate(HGT)
    endif

    endif     !   end of if do physics condition

!-------------------------------
!-stochastic-physics
! Update SPPT and/or SKEB pattern and perts

    if(DO_SPPT) then
       allocate(RNDPERT(IM,JM,LM),stat=STATUS)
       VERIFY_(STATUS)
       RNDPERT = 0.
       call MAPL_GetPointer ( GIM(GWD),  PREF,  'PREF', RC=STATUS)
       VERIFY_(STATUS)
       call sppt_pattern(CF,grid,RNDPERT,PREF,IM,JM,LM,DT)
       if( associated(RNDPTR) ) RNDPTR = RNDPERT
    endif

    if(DO_SKEB) then
       allocate(SKEBU_WT(IM,JM,LM),stat=STATUS)
       VERIFY_(STATUS)
       allocate(SKEBV_WT(IM,JM,LM),stat=STATUS)
       VERIFY_(STATUS)

       call MAPL_GetPointer ( GIM(GWD),  PREF,  'PREF', RC=STATUS)
       VERIFY_(STATUS)

       SKEBU_WT = 0. ; SKEBV_WT = 0.
       call skeb_pattern(CF,grid,SKEBU_WT,SKEBV_WT,PREF,IM,JM,LM,DT)
       if(associated (SKEBU) .and. associated(SKEBV)) then
          SKEBU = SKEBU_WT
          SKEBV = SKEBV_WT
       endif
    endif
!-stochastic-physics

! Fill the physics tendencies for the dynamics state variables.
! Q and other tracers are updated by their respective components
! and may be friendly to dynamics.
!---------------------------------------------------------------

!   NEED_TOT = associated(DTDTTOT) .or. associated(DTDT)
    NEED_TOT = .TRUE.
    NEED_FRI = associated(    TIF) .or. NEED_TOT
    NEED_STN = associated(    TIT) .or. NEED_TOT

    if(NEED_FRI) then
       allocate(FRI(IM,JM,LM),stat=STATUS)
       VERIFY_(STATUS)
       FRI         = INTDIS + TOPDIS
    end if

    if(NEED_STN) then
       allocate(STN(IM,JM,LM),stat=STATUS)
       VERIFY_(STATUS)
       STN = SIT*(1./MAPL_CP)
    end if

    if(associated(DUDT   )) DUDT    = UIM + UIT + UIG
    if(associated(DVDT   )) DVDT    = VIM + VIT + VIG
    if(associated(DWDT   )) DWDT    = WIM

!-stochastic-physics
    IF( DO_SPPT ) THEN
       allocate(TMP(IM,JM,LM),stat=STATUS)
       VERIFY_(STATUS)
       if( associated(DUDT) .and. associated(DVDT) ) then
          if( associated(DUDTUNPERT) ) DUDTUNPERT = DUDT
          if( associated(DVDTUNPERT) ) DVDTUNPERT = DVDT
          TMP = 0.
          DO L=1,LM
              TMP(:,:,L) = DUDT(:,:,L)*RNDPERT(:,:,L)
             DUDT(:,:,L) = DUDT(:,:,L) +   TMP(:,:,L)
          ENDDO
          if( associated(DUDTSTOCH) ) DUDTSTOCH = TMP
          TMP = 0.
          DO L=1,LM
              TMP(:,:,L) = DVDT(:,:,L)*RNDPERT(:,:,L)
             DVDT(:,:,L) = DVDT(:,:,L) +   TMP(:,:,L)
          ENDDO
          if( associated(DVDTSTOCH) ) DVDTSTOCH = TMP
       endif
    ENDIF

    IF ( DO_SKEB .and. associated(DUDT) .and. associated(DVDT) ) THEN
       DO L=1,LM
          DUDT(:,:,L)= DUDT(:,:,L) + SKEBU(:,:,L)
          DVDT(:,:,L)= DVDT(:,:,L) + SKEBV(:,:,L)
       ENDDO
    ENDIF
!-stochastic-physics

    if(associated(KEPHY   )) KEPHY = 0.0
    if(associated(PEPHY   )) PEPHY = 0.0
    if(associated(PERAD   )) PERAD = 0.0
    if(associated(PETRB   )) PETRB = 0.0
    if(associated(PEMST   )) PEMST = 0.0
    if(associated(PEFRI   )) PEFRI = 0.0
    if(associated(PEGWD   )) PEGWD = 0.0
    if(associated(PECUF   )) PECUF = 0.0

    if(associated(DUDT) .and. associated(DVDT) .and. associated(KEPHY)) then
       do L=1,LM
          KEPHY  = KEPHY  +  ((U(:,:,L)+(0.5*DT)*DUDT(:,:,L))*DUDT(:,:,L) +     &
                              (V(:,:,L)+(0.5*DT)*DVDT(:,:,L))*DVDT(:,:,L)   ) * &
                              DM(:,:,L)
       end do
    end if

    if(NEED_TOT) then
       allocate(TOT(IM,JM,LM),stat=STATUS)
       VERIFY_(STATUS)

       if ( .not.associated(TIR) .or. .not.associated(STN) .or. &
            .not.associated(TTN) .or. .not.associated(FRI) .or. &
            .not.associated(TIG) .or. .not.associated(TICU) ) then
            status=99
            if( MAPL_am_I_root() ) print*, "GEOS_PhysicsGridComp: missing T-tend pointer, aborting ..."
            VERIFY_(STATUS)
       endif

       if (SCM_NO_RAD) then
          TOT = STN   &  ! Mass-Weighted Temperature Tendency due to Turbulent Mixing
              + TTN   &  ! Mass-Weighted Temperature Tendency due to Moist Processes
              + FRI   &  ! Mass-Weighted Temperature Tendency due to Friction (Turbulence)
              + TIG   &  ! Mass-Weighted Temperature Tendency due to GWD
              + TICU     ! Mass-Weighted Temperature Tendency due to Cumulus Friction
       else
          TOT = TIR   &  ! Mass-Weighted Temperature Tendency due to Radiation
              + STN   &  ! Mass-Weighted Temperature Tendency due to Turbulent Mixing
              + TTN   &  ! Mass-Weighted Temperature Tendency due to Moist Processes
              + FRI   &  ! Mass-Weighted Temperature Tendency due to Friction (Turbulence)
              + TIG   &  ! Mass-Weighted Temperature Tendency due to GWD
              + TICU     ! Mass-Weighted Temperature Tendency due to Cumulus Friction
       end if

       IF(DO_SPPT) THEN
          allocate(TFORQS(IM,JM,LM))
          TFORQS = T + DT*TOT*DPI
          if( associated(DTDTUNPERT) ) DTDTUNPERT = TOT
          DO L=1,LM
            TMP(:,:,L) = (TOT(:,:,L) - TIG(:,:,L) ) *RNDPERT(:,:,L) ! Remove contribution from GWD before rndpert
            TOT(:,:,L) =  TOT(:,:,L) + TMP(:,:,L)
          ENDDO
          if( associated(DTDTSTOCH) ) DTDTSTOCH = TMP * DPI
       ENDIF

       if(associated(PERAD   )) then
          do L=1,LM
             PERAD = PERAD + TIR(:,:,L)*(MAPL_CP/MAPL_GRAV)
          end do
       end if

       if(associated(PETRB   )) then
          do L=1,LM
             PETRB = PETRB + STN(:,:,L)*(MAPL_CP/MAPL_GRAV)
          end do
       end if

       if(associated(PEMST   )) then
          do L=1,LM
             PEMST = PEMST + TTN(:,:,L)*(MAPL_CP/MAPL_GRAV)
          end do
       end if

       if(associated(PEFRI   )) then
          do L=1,LM
             PEFRI = PEFRI + FRI(:,:,L)*(MAPL_CP/MAPL_GRAV)
          end do
       end if

       if(associated(PEGWD   )) then
          do L=1,LM
             PEGWD = PEGWD + TIG(:,:,L)*(MAPL_CP/MAPL_GRAV)
          end do
       end if

       if(associated(PECUF   )) then
          do L=1,LM
             PECUF = PECUF + TICU(:,:,L)*(MAPL_CP/MAPL_GRAV)
          end do
       end if

       if(associated(DTDT    )) then
          DTDT     = TOT
          if(associated(PEPHY   )) then
             do L=1,LM
                PEPHY = PEPHY + DTDT(:,:,L)*(MAPL_CP/MAPL_GRAV)
             end do
          end if
       end if

       if(associated(DTDTTOT)) DTDTTOT = TOT * DPI
       if (DO_SPPT) then
          if(associated(DTDTUNPERT)) DTDTUNPERT = DTDTUNPERT * DPI
       end if
       deallocate(TOT)
    end if


    if(associated(DTDTRAD)) DTDTRAD = TIR * DPI
    if(associated(TIM    )) TIM     = TTN * DPI
    if(associated(TIMFRIC)) TIMFRIC = TICU* DPI
    if(associated(TIT    )) TIT     = STN * DPI
    if(associated(TIF    )) TIF     = FRI * DPI

   !  Compute Total Water Mass Change due to Physics Sources and Sinks
   !  ----------------------------------------------------------------
    allocate( DQ(IM,JM,LM) )
    DQ = QV+QLLS+QLCN+QILS+QICN+QRAIN+QSNOW+QGRAUPEL - QW

    if( DPEDT_PHYS ) then
       allocate(sumdq(IM,JM))
       allocate(  dpe(IM,JM))

       call ESMF_StateGet (EXPORT, 'TRADV', BUNDLE, RC=STATUS )
       VERIFY_(STATUS)
       call ESMF_FieldBundleGet ( BUNDLE, fieldCount=NQ, RC=STATUS )
       VERIFY_(STATUS)

       allocate( NAMES(NQ),STAT=STATUS )
       VERIFY_(STATUS)
       call ESMF_FieldBundleGet ( BUNDLE, fieldNameList=NAMES, rc=STATUS )
       VERIFY_(STATUS)

       !Re-initialize water rescale increment bundle
       !----------------------------------------------
       call Initialize_IncBundle_run(EXPORT, EXPORT, H2Oinc, __RC__)

       ! Add diagnostic for scaling tendency of QI and QL -> fill diagnostic with "before" value
       ! ------------------------------------------------
       if( associated(DQVDTSCL) ) DQVDTSCL =  QV
       if( associated(DQLDTSCL) ) DQLDTSCL = (QLLS+QLCN)
       if( associated(DQIDTSCL) ) DQIDTSCL = (QILS+QICN)

   !  Modify P and Q such that Pdry is conserved
   !  ------------------------------------------
       ple_new = ple*1.0_8
       sumdq = 0.0
       DPDT(:,:,0) = 0.0
       do l=1,lm
                    sumdq = sumdq + dq(:,:,L) * ( ple(:,:,L)-ple(:,:,L-1) ) / DT
              dpdt(:,:,L) = sumdq
           ple_new(:,:,L) = ple_new(:,:,L) + dt*dpdt(:,:,L)
               dpe(:,:)   = (ple(:,:,L)-ple(:,:,L-1)) / (ple_new(:,:,L)-ple_new(:,:,L-1))
             do N=1,NQ
                call ESMFL_BundleGetPointertoData( BUNDLE, trim(NAMES(N)), PTR3D, RC=STATUS)
                VERIFY_(STATUS)
                if( trim(NAMES(N)) /= 'CLCN'       .and. & ! Exclude: Advected Convective and Large-Scale
                    trim(NAMES(N)) /= 'CLLS'     )  then   ! -------- Cloud Fractions
                PTR3D(:,:,L) = PTR3d(:,:,L) * dpe(:,:)
                endif
             end do
       end do

       ! Compute water rescale increments
       !----------------------------------
       call Compute_IncBundle(EXPORT, EXPORT, H2Oinc, STATE, __RC__)

       ! Add diagnostic for scaling tendency of QI and QL -> update diagnostic with "after" value
       ! ------------------------------------------------
       if( associated(DQVDTSCL) ) DQVDTSCL = (    QV     - DQVDTSCL )/DT
       if( associated(DQLDTSCL) ) DQLDTSCL = ( QLLS+QLCN - DQLDTSCL )/DT
       if( associated(DQIDTSCL) ) DQIDTSCL = ( QILS+QICN - DQIDTSCL )/DT

       deallocate( sumdq   )
       deallocate( dpe     )
       deallocate( names   )
       deallocate( sumq    )
       deallocate( ple_new )

    else
                                      DPDT = 0.0
       if( associated(DQVDTSCL) ) DQVDTSCL = 0.0
       if( associated(DQLDTSCL) ) DQLDTSCL = 0.0
       if( associated(DQIDTSCL) ) DQIDTSCL = 0.0
    endif

    deallocate( DQ )

    if(associated(DMDT)) DMDT(:,:) = DPDT(:,:,LM)*(1.0/MAPL_GRAV)

    if( associated(DPDTPHY) ) then
        do L=1,LM
           DPDTPHY(:,:,L) = dpdt(:,:,L)
        enddo
    endif

    allocate( TDPNEW(IM,JM,LM),stat=STATUS )
    VERIFY_(STATUS)
    do L=1,LM
       TDPNEW(:,:,L) = ( T(:,:,L) + DT*DTDT(:,:,L)*DPI(:,:,L) ) * ( PLE(:,:,L)-PLE(:,:,L-1) + DT*(DPDT(:,:,L)-DPDT(:,:,L-1)) )
    enddo
       DTDT = ( TDPNEW - TDPOLD )/DT
    deallocate( TDPNEW )
    deallocate( TDPOLD )

    if(associated(FTU)) then
       FTU(:,:,0) = 0.0
       do L=1,LM
          FTU(:,:,L) = FTU(:,:,L-1) - UIT(:,:,L)*(ZLE(:,:,L)-ZLE(:,:,L-1))
       end do
    end if

    if(associated(FTV)) then
       FTV(:,:,0) = 0.0
       do L=1,LM
          FTV(:,:,L) = FTV(:,:,L-1) - VIT(:,:,L)*(ZLE(:,:,L)-ZLE(:,:,L-1))
       end do
    end if

! QV Tendencies
! -------------
    IF(DO_SPPT) THEN

        if ( .not.associated(DQVDTMST) .or. .not.associated(DQVDTTRB) .or. .not.associated(DQVDTCHM) ) then
           status=99
           if( MAPL_am_I_root() ) print*, "GEOS_PhysicsGridComp: missing Q-tend pointer, aborting ..."
            VERIFY_(STATUS)
        endif

        ! Create Proxy for Updated Pressure due to Moist Physics
        !-------------------------------------------------------
         allocate(  dp_mst(IM,JM,  LM),STAT=STATUS ) ; VERIFY_(STATUS)
         allocate( ple_mst(IM,JM,0:LM),STAT=STATUS ) ; VERIFY_(STATUS)
         ple_mst = ple + dt*dpdt
          dp_mst = ple_mst(:,:,1:LM)-ple_mst(:,:,0:LM-1)

        ! Create Water Mass before Stochastic Perturbation
        ! ------------------------------------------------
          allocate( qdp_b4( IM,JM,LM ),STAT=STATUS ) ; VERIFY_(STATUS)
          do L=1,lm
              qdp_b4(:,:,L) = (   qv(:,:,L)                + &
                                qlls(:,:,L) +  qlcn(:,:,L) + &
                                qils(:,:,L) +  qicn(:,:,L) + &
                               qrain(:,:,L) + qsnow(:,:,L) + qgraupel(:,:,L) )*dp_mst(:,:,L)
          enddo

        ! Compute Stochastic Perturbation
        ! -------------------------------
          if( associated(QVUNPERT)    ) QVUNPERT    = QV
          if( associated(DQVDTPERT)   ) DQVDTPERT   = QV
          if( associated(DQVDTSTOCH)  ) DQVDTSTOCH  = QV
          if( associated(DQVDTUNPERT) ) DQVDTUNPERT = DQVDTMST + DQVDTTRB + DQVDTCHM
          DO L=1,LM
            TMP(:,:,L) = (DQVDTMST(:,:,L)+DQVDTTRB(:,:,L))*RNDPERT(:,:,L)
            if( associated(DQVDTSTOCH) ) DQVDTSTOCH(:,:,L) = TMP(:,:,L)
            TMP(:,:,L) = QV(:,:,L) + max(DT*TMP(:,:,L),-QV(:,:,L))
          ENDDO
          allocate(qs(IM,JM))
          allocate(pmean(IM,JM))
          DO L=1,LM
             pmean = 0.5*(PLE(:,:,L-1)+PLE(:,:,L))
             qs = GEOS_Qsat(TFORQS(:,:,L),pmean,PASCALS=.TRUE.)
             where (TMP(:,:,L) < qs)
                QV(:,:,L) = TMP(:,:,L)
             endwhere
          ENDDO
!         if( associated(DQVDTSTOCH) ) DQVDTSTOCH = (QV - DQVDTSTOCH)/DT

        ! Create Water Mass after Stochastic Perturbation
        ! -----------------------------------------------
          allocate( qdp_af( IM,JM,LM ),STAT=STATUS ) ; VERIFY_(STATUS)
          do L=1,lm
              qdp_af(:,:,L) = (   qv(:,:,L)                + &
                                qlls(:,:,L) +  qlcn(:,:,L) + &
                                qils(:,:,L) +  qicn(:,:,L) + &
                               qrain(:,:,L) + qsnow(:,:,L) + qgraupel(:,:,L) )*dp_mst(:,:,L)
          enddo

        ! Vertically Integrate Water Mass where they Differ
        ! -------------------------------------------------
          allocate( sum_qdp_b4( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)
          allocate( sum_qdp_af( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)
          sum_qdp_b4 = 0.0_8
          sum_qdp_af = 0.0_8
          do L=1,lm
          where( qdp_b4(:,:,L).ne.qdp_af(:,:,L) )
                 sum_qdp_b4 = sum_qdp_b4 + qdp_b4(:,:,L)
                 sum_qdp_af = sum_qdp_af + qdp_af(:,:,L)
          end where
          enddo

        ! Compute Area-Mean Vertically Integrated B4 Water Mass
        ! -----------------------------------------------------
          allocate( qdp_b4_int( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)
          where( sum_qdp_b4.ne.0.0_8 )
                 qdp_b4_int = sum_qdp_b4
          elsewhere
                 qdp_b4_int = MAPL_UNDEF
          end where
          call MAPL_AreaMean( qdp_b4_ave, qdp_b4_int, area, grid, rc=STATUS )
          VERIFY_(STATUS)

        ! Compute Area-Mean Vertically Integrated AF Water Mass
        ! -----------------------------------------------------
          allocate( qdp_af_int( IM,JM ),STAT=STATUS ) ; VERIFY_(STATUS)
          where( sum_qdp_af.ne.0.0_8 )
                 qdp_af_int = sum_qdp_af
          elsewhere
                 qdp_af_int = MAPL_UNDEF
          end where
          call MAPL_AreaMean( qdp_af_ave, qdp_af_int, area, grid, rc=STATUS )
          VERIFY_(STATUS)

        ! Compute Dry-Mass Scaling Parameter
        ! ----------------------------------
          if( real(qdp_b4_ave,kind=4).ne.MAPL_UNDEF .and. &
              real(qdp_af_ave,kind=4).ne.MAPL_UNDEF       ) then
              gamma = real( qdp_b4_ave / qdp_af_ave , kind=4 )
          else
              gamma = 1.0_8
          endif

        ! Scale Water Variables for Dry-Mass Conservation
        ! -----------------------------------------------
          do L=1,lm
              where(  qdp_af(:,:,L).ne. qdp_b4(:,:,L) )
                        qv  (:,:,L) =     qv  (:,:,L) * gamma
                        qlls(:,:,L) =     qlls(:,:,L) * gamma
                        qlcn(:,:,L) =     qlcn(:,:,L) * gamma
                        qils(:,:,L) =     qils(:,:,L) * gamma
                        qicn(:,:,L) =     qicn(:,:,L) * gamma
                       qrain(:,:,L) =    qrain(:,:,L) * gamma
                       qsnow(:,:,L) =    qsnow(:,:,L) * gamma
                    qgraupel(:,:,L) = qgraupel(:,:,L) * gamma
              end where
          enddo

          if( associated(QVPERT)    ) QVPERT = QV
          if( associated(DQVDTPERT) ) DQVDTPERT = (QV-DQVDTPERT)/DT

          deallocate(  dp_mst    )
          deallocate( ple_mst    )
          deallocate( qdp_b4     )
          deallocate( qdp_af     )
          deallocate( qdp_b4_int )
          deallocate( qdp_af_int )
          deallocate( sum_qdp_b4 )
          deallocate( sum_qdp_af )

          deallocate(pmean,qs,TFORQS)
    ENDIF

    if(associated(DQVDTPHYINT)) then
       DQVDTPHYINT = 0.0
       do L=1,LM
       DQVDTPHYINT = DQVDTPHYINT + ( DQVDTMST(:,:,L) &
                                   + DQVDTTRB(:,:,L) &
                                   + DQVDTCHM(:,:,L) ) * DM(:,:,L)
       end do
    end if

    if(associated(DQVDTTRBINT)) then
       DQVDTTRBINT = 0.0
       do L=1,LM
       DQVDTTRBINT = DQVDTTRBINT + DQVDTTRB(:,:,L)*DM(:,:,L)
       end do
    end if

    if(associated(DQVDTMSTINT)) then
       DQVDTMSTINT = 0.0
       do L=1,LM
       DQVDTMSTINT = DQVDTMSTINT + DQVDTMST(:,:,L)*DM(:,:,L)
       end do
    end if

    if(associated(DQVDTCHMINT)) then
       DQVDTCHMINT = 0.0
       do L=1,LM
       DQVDTCHMINT = DQVDTCHMINT + DQVDTCHM(:,:,L)*DM(:,:,L)
       end do
    end if

! QL Tendencies
! -------------
    if(associated(DQLDTPHYINT)) then
       DQLDTPHYINT = 0.0
       do L=1,LM
       DQLDTPHYINT = DQLDTPHYINT + DQLDTMST(:,:,L)*DM(:,:,L)
       end do
    end if

    if(associated(DQLDTMSTINT)) then
       DQLDTMSTINT = 0.0
       do L=1,LM
       DQLDTMSTINT = DQLDTMSTINT + DQLDTMST(:,:,L)*DM(:,:,L)
       end do
    end if

! QI Tendencies
! -------------
    if(associated(DQIDTPHYINT)) then
       DQIDTPHYINT = 0.0
       do L=1,LM
       DQIDTPHYINT = DQIDTPHYINT + DQIDTMST(:,:,L)*DM(:,:,L)
       end do
    end if

    if(associated(DQIDTMSTINT)) then
       DQIDTMSTINT = 0.0
       do L=1,LM
       DQIDTMSTINT = DQIDTMSTINT + DQIDTMST(:,:,L)*DM(:,:,L)
       end do
    end if

! OX Tendencies
! -------------
    if(associated(DOXDTPHYINT)) then
       DOXDTPHYINT = 0.0
       do L=1,LM
       DOXDTPHYINT = DOXDTPHYINT + DOXDTCHM(:,:,L)*DM(:,:,L)
       end do
       DOXDTPHYINT = DOXDTPHYINT * (MAPL_O3MW/MAPL_AIRMW)
    end if

    if(associated(DOXDTCHMINT)) then
       DOXDTCHMINT = 0.0
       do L=1,LM
       DOXDTCHMINT = DOXDTCHMINT + DOXDTCHM(:,:,L)*DM(:,:,L)
       end do
       DOXDTCHMINT = DOXDTCHMINT * (MAPL_O3MW/MAPL_AIRMW)
    end if

    if(SYNCTQ.eq.0.) then
      !- save 'boundary layer' tendencies of Q and T for the convection scheme
      DQDT_BL = DQVDTTRB
      DTDT_BL = 0.
      !- for SCM setup, TIT/TIF are not associated
      if( associated(TIF)) DTDT_BL = DTDT_BL + TIF
      if( associated(TIT)) DTDT_BL = DTDT_BL + TIT
    endif

    if(associated(DM )) deallocate(DM )
    if(associated(DPI)) deallocate(DPI)
    if(associated(FRI)) deallocate(FRI)
    if(associated(STN)) deallocate(STN)

!-stochastic-physics
    if(DO_SPPT) deallocate(TMP,RNDPERT)
    if(DO_SKEB) deallocate(SKEBU_WT,SKEBV_WT)
!-stochastic-physics

    deallocate( zero )

    call MAPL_TimerOff(STATE,"RUN")


    call MAPL_TimerOff(STATE,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  end subroutine Run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine VertInterp(v2,v3,ple,pp,rc)

    real    , intent(OUT) :: v2(:,:)
    real    , intent(IN ) :: v3(:,:,:)
    real    , intent(IN ) :: ple(:,:,:)
    real    , intent(IN ) :: pp
    integer, optional, intent(OUT) :: rc

    real, dimension(size(v2,1),size(v2,2)) :: al,PT,PB
    integer k,km
    logical edge

    character*(10) :: Iam='VertInterp'

    km   = size(ple,3)-1
    edge = size(v3,3)==km+1

    _ASSERT(edge .or. size(v3,3)==km,'needs informative message')

    v2   = MAPL_UNDEF

    if(EDGE) then
       pb   = ple(:,:,km+1)
       do k=km,1,-1
          pt = ple(:,:,k)
          if(all(pb<pp)) exit
          where(pp>pt .and. pp<=pb)
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,k)*al + v3(:,:,k+1)*(1.0-al)
          end where
          pb = pt
       end do
    else
       pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
       do k=km,2,-1
          pt = 0.5*(ple(:,:,k-1)+ple(:,:,k))
          if(all(pb<pp)) exit
          where( (pp>pt.and.pp<=pb) )
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,k-1)*al + v3(:,:,k)*(1.0-al)
          end where
          pb = pt
       end do
       pt = 0.5*(ple(:,:,km)+ple(:,:,km-1))
       pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
          where( (pp>pb.and.pp<=ple(:,:,km+1)) )
             v2 = v3(:,:,km)
          end where
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine VertInterp

end module GEOS_PhysicsGridCompMod
