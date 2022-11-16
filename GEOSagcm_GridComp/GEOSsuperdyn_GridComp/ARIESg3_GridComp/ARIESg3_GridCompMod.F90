!  $Id: ARIESg3_GridCompMod.F90,v 1.19.112.1.2.1 2019/07/23 15:31:37 mmanyin Exp $

#include "MAPL_Generic.h"



!-----------------------------------------------------------------------
!              ESMA - Earth System Modeling Applications
!-----------------------------------------------------------------------
   Module ARIESg3_GridCompMod

!BOP
!
! !MODULE: ARIESg3_GridCompMod --- ARIESr/GEOS3 Dynamical Core Grid Component
!


! !USES:

   use ESMF                ! ESMF base class
   use MAPL                ! GEOS base class
   use dynamics_vars, only : T_TRACERS, T_FVDYCORE_VARS, &
                             T_FVDYCORE_GRID, T_FVDYCORE_STATE

! !PUBLIC MEMBER FUNCTIONS:

  implicit none
  private

  public  SetServices      ! Register component methods

! !DESCRIPTION: This module implements the FVCAM Dynamical Core as
!               an ESMF gridded component.
!
! \paragraph*{Overview}
!
!   This module contains an ESMF wrapper for the Finite-Volume
!   Dynamical Core used in the Community Atmospheric Model
!   (FVCAM). This component will hereafter be referred
!   to as the ``FVdycore'' ESMF gridded component.  FVdycore
!   consists of four sub-components,
!
!   \begin{itemize}
!      \item {\tt cd\_core:}  The C/D-grid dycore component
!      \item {\tt te\_map:}   Vertical remapping algorithm
!      \item {\tt trac2d:}    Tracer advection
!      \item {\tt benergy:}   Energy balance
!   \end{itemize}
!
!   Subsequently the ESMF component design for FV dycore
!   will be described.   
!
! \paragraph*{Internal State}
!
!  FVdycore maintains an internal state consisting of the
!  following fields:  control variables
!
!   \begin{itemize}
!     \item {\tt U}:    U winds on a D-grid (m/s)
!     \item {\tt V}:    V winds on a D-grid (m/s)
!     \item {\tt PT}:   Scaled Virtual Potential Temperature(T$_v$/PKZ)
!     \item {\tt PE}:   Edge pressures
!     \item {\tt Q}:    Tracers
!     \item {\tt PKZ}:  Consistent mean for p$^\kappa$
!   \end{itemize}
!
!  as well as a GRID (to be mentioned later) 
!  and same additional run-specific variables 
!  (dt, iord, jord, nsplit -- to be mentioned later)
!
! Note: {\tt PT} is not updated if the flag {\tt CONVT} is true.
!
! The internal state is updated each time FVdycore is called.
!
! \paragraph*{Import State}
!
! The import state consists of the tendencies of the 
! control variables plus the surface geopotential heights:
!
!   \begin{itemize}
!     \item {\tt DUDT}:    U wind tendency on a A-grid (m/s)
!     \item {\tt DVDT}:    V wind tendency on a A-grid (m/s)
!     \item {\tt DTDT}:    Delta-pressure-weighted temperature tendency
!     \item {\tt DPEDT}:   Edge pressure tendency
!     \item {\tt PHIS}:    Surface Geopotential Heights
!   \end{itemize}
!
! These are by definition on an A-grid and have an XY
! domain decomposition.
!
! \paragraph*{Export State}
!
!   The export state can provide the following variables:
!
!   \begin{itemize}
!     \item {\tt U}:         U winds on a A-grid (m/s)
!     \item {\tt V}:         V winds on a A-grid (m/s)
!     \item {\tt U\_CGRID}:   U winds on a C-grid (m/s)
!     \item {\tt V\_CGRID}:   V winds on a C-grid (m/s)
!     \item {\tt U\_DGRID}:   U winds on a D-grid (m/s)
!     \item {\tt V\_DGRID}:   V winds on a D-grid (m/s)
!     \item {\tt T}:         Temperature (K)
!     \item {\tt Q}:         Tracers
!     \item {\tt TH}:        Potential Temperature (K)
!     \item {\tt ZL}:        Mid-Layer Heights (m)
!     \item {\tt ZLE}:       Edge Heights (m)
!     \item {\tt PLE}:       Edge pressures (Pa)
!     \item {\tt PLK}:       $P^\kappa$ at Mid-Layers
!     \item {\tt OMEGA}:     Vertical pressure velocity (pa/s)
!     \item {\tt PTFX}:      Mass-Weighted PT flux on C-Grid (K Pa m$^2$/s)
!     \item {\tt PTFY}:      Mass-Weighted PT flux on C-Grid (K Pa m$^2$/s)
!     \item {\tt MFX\_UR}:    Mass-Weighted U-Wind on C-Grid (Pa m$^2$/s)
!     \item {\tt MFY\_UR}:    Mass-Weighted V-wind on C-Grid (Pa m$^2$/s)
!     \item {\tt MFX}:       Remapped Mass-Weighted U-Wind on C-Grid (Pa m$^2$/s)
!     \item {\tt MFY}:       Remapped Mass-Weighted V-wind on C-Grid (Pa m$^2$/s)
!     \item {\tt MFZ}:       Remapped Vertical mass flux (kg/(m$^2$*s))
!     \item {\tt MFX\_A}:     Remapped Mass-Weighted U-Wind on A-Grid (Pa m$^2$/s)
!     \item {\tt MFY\_A}:     Remapped Mass-Weighted V-wind on A-Grid (Pa m$^2$/s)
!     \item {\tt PV}:        Ertel's Potential Vorticity (m$^2$ / kg*s)
!     \item {\tt DUDT}:      U-wind Tendency (m/s/s)
!     \item {\tt DVDT}:      V-wind Tendency (m/s/s)
!     \item {\tt DTDT}:      Mass-Weighted Temperature Tendency (Pa K/s)
!     \item {\tt AREA}:      Cell areas on the A-Grid (m$^2$, polar caps at J=1, J=JM) 
!   \end{itemize}
!
!   All variables are on an A-grid with points at the poles, and have an XY decomposition.
!
! \paragraph*{Grids and Decompositions}
!
!   The current version supports only a 1D latitude-based
!   decomposition of the domain (with OMP task-parallelism
!   in the vertical, resulting in reasonable scalability 
!   on large PE configurations).  In the near future it will 
!   support a 2D domain decomposition, in which import and
!   export state are decomposed in longitude and latitude,
!   while the internal state (for the most part) is 
!   decomposed in latitude and level.  When needed, 
!   the data is redistributed (``transposed'') internally.
!
!   There are two fundamental ESMF grids in use;
!   \begin{itemize}
!     \item {GRIDXY}: longitude-latitude ESMF grid (public)
!     \item {GRIDYZ}: A latitude-level cross-sectional
!                     decomposition (private to this module) 
!   \end{itemize}
!
!   PILGRIM will be used for communication until ESMF has 
!   sufficient functionality and performance to take over 
!   the task.  The use of pilgrim requires a call to 
!   {\tt INIT\_SPMD} to set SPMD parameters, decompositions,
!   etc.
!
!   Currently, only a 1D decomposition in latitude is employed.
!   Thus GRIDXY and GRIDYZ actually represent the same 
!   decomposition and no transposes are employed.
!
! \paragraph*{Required Files}
!
!  The following files are needed for a standard restart run:
!
!  \begin{itemize}
!    \item Layout file
!      \begin{itemize}
!        \item {\tt nprxy\_x, nprxy\_y, npryz\_y, npryz\_z}:
!          process dimensions in XY and YZ.
!        \item {\tt imxy, jmxy, jmyz, kmyz}: distributions for XY and YZ
!        \item {\tt iord, jord}: the order of the lon. and lat. algorithms
!        \item {\tt dtime}:  The large (advection) time step
!        \item {\tt nsplit}: the ratio between the large and small time step 
!          (possibly zero for automatic determination), 
!      \end{itemize}
!    \item Restart file
!      \begin{itemize}
!        \item date in standard format yy, mm, dd, hh, mm, ss
!        \item dimensions im, jm, km, nq
!        \item control variables {\tt U, V, PT, PE, Q}
!      \end{itemize}
!    \item Topography file
!
!  \end{itemize}
!
! \paragraph*{Future Additions}
!
!  \begin{itemize}
!    \item  Conservation of energy (CONSV  == .TRUE. ) 
!    \item  2D decomposition (requires transposes in the coupler)
!    \item  Use r8 instead of r4 (currently supported in StopGap)
!  \end{itemize}
!
!EOP

! !REVISION HISTORY:
!
! 11Jul2003  Sawyer    From Trayanov/da Silva EVAC 
! 23Jul2003  Sawyer    First informal tiptoe-through
! 29Jul2003  Sawyer    Modifications based on comments from 23Jul2003
! 28Aug2003  Sawyer    First check-in; Internal state to D-grid
! 15Sep2003  Sawyer    Extensive bug fixes, revisions
! 24Sep2003  Sawyer    Modified names; corrected weighting of T, Q
! 22Oct2003  Sawyer    pmgrid removed (data now in spmd_dyn)
! 25Nov2003  Sawyer    Optimization for 1D decomposition (as in FVCAM)
! 03Dec2003  Sawyer    Switched over to specified decompositions
! 04Dec2003  Sawyer    Moved T_FVDYCORE_GRID to dynamics_vars
! 21Jan2004  Takacs    Modified Import/Export, Added Generic State, Added TOPO utility
! 20Sep2004  Sawyer    Revised cd_core, trac2d interfaces, refactoring
! 06Oct2004  Sawyer    More refactoring, removed spmd_dyn
! 17Feb2005  Sawyer    Added Ertel's potential vorticity to diagnostics
! 20Mar2005  Sawyer    Tracers are now pointers into import state
! 12Apr2005  Sawyer    Extensive changes to minimize tracer memory
! 18May2005  Sawyer    Put FVdycore_wrapper in separate file; CAM/GEOS5 merge
! 16Nov2005  Takacs    Added option for DCADJ, Merge with Daedalus_p5
! 18Jan2006  Putman    Added mass fluxes to export state
! 01Apr2009  Sawyer    Upgraded to PILGRIM from cam3_6_33
!
!----------------------------------------------------------------------

  integer,  parameter :: r8           = 8
  integer,  parameter :: r4           = 4

  real(r8), parameter :: RADIUS       = MAPL_RADIUS
  real(r8), parameter :: CP           = MAPL_CP
  real(r8), parameter :: PI           = MAPL_PI_R8
  real(r8), parameter :: OMEGA        = MAPL_OMEGA
  real(r8), parameter :: KAPPA        = MAPL_KAPPA
  real(r8), parameter :: P00          = MAPL_P00
  real(r8), parameter :: GRAV         = MAPL_GRAV
  real(r8), parameter :: RGAS         = MAPL_RGAS
  real(r8), parameter :: RVAP         = MAPL_RVAP
  real(r8), parameter :: EPS          = RVAP/RGAS-1.0

  integer,  parameter :: TIME_TO_RUN  = 1
  integer,  parameter :: CHECK_MAXMIN = 2

  integer :: I, J, K  !  Default declaration for loops.

! Wrapper for extracting internal state
! -------------------------------------

  type DYN_wrap
     type (T_FVDYCORE_STATE), pointer :: DYN_STATE
  end type DYN_wrap

contains

!----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Set services for FVCAM Dynamical Core
 
! !INTERFACE:

   Subroutine SetServices ( gc, rc )

! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: gc     ! gridded component
   integer, intent(out), optional     :: rc     ! return code
    

! !DESCRIPTION: Set services (register) for the FVCAM Dynamical Core
!               Grid Component.
!         
!EOP         
!----------------------------------------------------------------------
  
    type (T_FVDYCORE_STATE), pointer :: dyn_internal_state 
    type (DYN_wrap)                  :: wrap

    integer                          :: status
    character(len=ESMF_MAXSTR)       :: IAm
    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Begin
!------

    Iam = "SetServices"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam
    
! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( dyn_internal_state, stat=status )
    VERIFY_(STATUS)
    wrap%dyn_state => dyn_internal_state

! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC,'FVstate',wrap,status )
    VERIFY_(STATUS)

!BOS
! !IMPORT STATE:

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_tendency',                    &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_tendency',                   &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'delta-p_weighted_temperature_tendency',     &
         UNITS      = 'Pa K s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQVANA',                                          &
         LONG_NAME  = 'specific_humidity_vapor_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                         &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                   RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQLANA',                                           &
         LONG_NAME  = 'specific_humidity_liquid_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                          &
         DIMS       = MAPL_DimsHorzVert,                                  &
         VLOCATION  = MAPL_VLocationCenter,                    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQIANA',                                        &
         LONG_NAME  = 'specific_humidity_ice_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                       &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DOXANA',                                    &
         LONG_NAME  = 'ozone_increment_from_analysis',             &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'edge_pressure_tendency',                    &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'PHIS',                                      &
         LONG_NAME  = 'surface_geopotential_height',               &
         UNITS      = 'm+2 s-2',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec( gc,                              &
        SHORT_NAME = 'TRADV',                                        &
        LONG_NAME  = 'advected_quantities',                        &
        UNITS      = 'unknown',                                    &
        DATATYPE   = MAPL_BundleItem,               &
        RC=STATUS  )
    VERIFY_(STATUS)

! !EXPORT STATE:

    call MAPL_AddExportSpec ( gc,                                         &
         SHORT_NAME       = 'KE',                                         &
         LONG_NAME        = 'vertically_integrated_kinetic_energy',       &
         UNITS            = 'J m-2'  ,                                    &
         DIMS             = MAPL_DimsHorzOnly,                            &
         VLOCATION        = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'TAVE',                                                               &
         LONG_NAME  = 'vertically_averaged_dry_temperature',                                &
         UNITS      = 'K',                                                                  &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'UAVE',                                                               &
         LONG_NAME  = 'vertically_averaged_zonal_wind',                                     &
         UNITS      = 'm s-1',                                                            &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         FIELD_TYPE = MAPL_VectorField,                                                     &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KEPHY',                                                              &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_physics',       &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                &
       SHORT_NAME   = 'PEPHY',                                                   &
       LONG_NAME    = 'total_potential_energy_tendency_due_to_physics',          &
       UNITS        = 'W m-2',                                                   &
       DIMS         = MAPL_DimsHorzOnly,                                         &
       VLOCATION    = MAPL_VLocationNone,                              RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                &
       SHORT_NAME   = 'TEPHY',                                                   &
       LONG_NAME    = 'mountain_work_tendency_due_to_physics',                   &
       UNITS        = 'W m-2',                                                   &
       DIMS         = MAPL_DimsHorzOnly,                                         &
       VLOCATION    = MAPL_VLocationNone,                              RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                &
       SHORT_NAME         = 'KEANA',                                             &
       LONG_NAME          = 'total_kinetic_energy_tendency_due_to_analysis',     &
       UNITS              = 'W m-2',                                             &
       DIMS               = MAPL_DimsHorzOnly,                                   &
       VLOCATION          = MAPL_VLocationNone,                        RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                &
       SHORT_NAME         = 'PEANA',                                             &
       LONG_NAME          = 'total_potential_energy_tendency_due_to_analysis',   &
       UNITS              = 'W m-2',                                             &
       DIMS               = MAPL_DimsHorzOnly,                                   &
       VLOCATION          = MAPL_VLocationNone,                        RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                &
       SHORT_NAME         = 'TEANA',                                             &
       LONG_NAME          = 'mountain_work_tendency_due_to_analysis',            &
       UNITS              = 'W m-2',                                             &
       DIMS               = MAPL_DimsHorzOnly,                                   &
       VLOCATION          = MAPL_VLocationNone,                        RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                             &
         SHORT_NAME = 'KEHOT',                                                                &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_HOT',             &
         UNITS      = 'W m-2',                                                                &
         DIMS       = MAPL_DimsHorzOnly,                                                      &
         VLOCATION  = MAPL_VLocationNone,                                          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                             &
         SHORT_NAME = 'KEDP',                                                                 &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_pressure_change', &
         UNITS      = 'W m-2',                                                                &
         DIMS       = MAPL_DimsHorzOnly,                                                      &
         VLOCATION  = MAPL_VLocationNone,                                          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                             &
         SHORT_NAME = 'KEADV',                                                                &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_dynamics_advection',      &
         UNITS      = 'W m-2',                                                                &
         DIMS       = MAPL_DimsHorzOnly,                                                      &
         VLOCATION  = MAPL_VLocationNone,                                          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                             &
         SHORT_NAME = 'KEPG',                                                                 &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_pressure_gradient',      &
         UNITS      = 'W m-2',                                                                &
         DIMS       = MAPL_DimsHorzOnly,                                                      &
         VLOCATION  = MAPL_VLocationNone,                                          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KEDYN',                                                              &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_dynamics',      &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEDYN',                                                              &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_dynamics',    &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'TEDYN',                                                              &
         LONG_NAME  = 'mountain_work_tendency_due_to_dynamics',                             &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KECDCOR',                                                            &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_cdcore',        &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PECDCOR',                                                            &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_cdcore',      &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'TECDCOR',                                                            &
         LONG_NAME  = 'mountain_work_tendency_due_to_cdcore',                               &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'QFIXER',                                                             &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_CONSV',       &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KEREMAP',                                                            &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_due_to_remap',         &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'PEREMAP',                                                            &
         LONG_NAME  = 'vertically_integrated_potential_energy_tendency_due_to_remap',       &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'TEREMAP',                                                            &
         LONG_NAME  = 'mountain_work_tendency_due_to_remap',                                &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'KEGEN',                                                              &
         LONG_NAME  = 'vertically_integrated_generation_of_kinetic_energy',                 &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DKERESIN',                                                           &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_residual_from_inertial_terms',  &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DKERESPG',                                                           &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_residual_from_PG_terms',        &
         UNITS      = 'W m-2',                                                              &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DMDTANA',                                                            &
         LONG_NAME  = 'vertically_integrated_mass_tendency_due_to_analysis',                &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DOXDTANAINT',                                                        &
         LONG_NAME  = 'vertically_integrated_ozone_tendency_due_to_analysis',               &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQVDTANAINT',                                                        &
         LONG_NAME  = 'vertically_integrated_water_vapor_tendency_due_to_analysis',         &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQLDTANAINT',                                                        &
         LONG_NAME  = 'vertically_integrated_liquid_water_tendency_due_to_analysis',        &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQIDTANAINT',                                                        &
         LONG_NAME  = 'vertically_integrated_ice_water_tendency_due_to_analysis',           &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DMDTDYN',                                                            &
         LONG_NAME  = 'vertically_integrated_mass_tendency_due_to_dynamics',                &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DOXDTDYNINT',                                                        &
         LONG_NAME  = 'vertically_integrated_ozone_tendency_due_to_dynamics',               &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DTHVDTDYNINT',                                                       &
         LONG_NAME  = 'vertically_integrated_THV_tendency_due_to_dynamics',                 &
         UNITS      = 'K kg m-2 s-1',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DTHVDTREMAP',                                                        &
         LONG_NAME  = 'vertically_integrated_THV_tendency_due_to_vertical_remapping',       &
         UNITS      = 'K kg m-2 s-1',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DTHVDTCONSV',                                                        &
         LONG_NAME  = 'vertically_integrated_THV_tendency_due_to_TE_conservation',          &
         UNITS      = 'K kg m-2 s-1',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DTHVDTPHYINT',                                                       &
         LONG_NAME  = 'vertically_integrated_THV_tendency_due_to_physics',                  &
         UNITS      = 'K kg m-2 s-1',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DTHVDTANAINT',                                                       &
         LONG_NAME  = 'vertically_integrated_THV_tendency_due_to_analysis',                 &
         UNITS      = 'K kg m-2 s-1',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQVDTDYNINT',                                                        &
         LONG_NAME  = 'vertically_integrated_water_vapor_tendency_due_to_dynamics',         &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQLDTDYNINT',                                                        &
         LONG_NAME  = 'vertically_integrated_liquid_water_tendency_due_to_dynamics',        &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQIDTDYNINT',                                                        &
         LONG_NAME  = 'vertically_integrated_ice_water_tendency_due_to_dynamics',           &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                    &
         SHORT_NAME = 'CONVKE',                                                      &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_convergence',            &
         UNITS      = 'W m-2',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                             &
         VLOCATION  = MAPL_VLocationCenter,                               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                    &
         SHORT_NAME = 'CONVTHV',                                                     &
         LONG_NAME  = 'vertically_integrated_thetav_convergence',                    &
         UNITS      = 'W m-2',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                             &
         VLOCATION  = MAPL_VLocationCenter,                               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                    &
         SHORT_NAME = 'CONVCPT',                                                     &
         LONG_NAME  = 'vertically_integrated_enthalpy_convergence',                  &
         UNITS      = 'W m-2',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                             &
         VLOCATION  = MAPL_VLocationCenter,                               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                    &
         SHORT_NAME = 'CONVPHI',                                                     &
         LONG_NAME  = 'vertically_integrated_geopotential_convergence',              &
         UNITS      = 'W m-2',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                             &
         VLOCATION  = MAPL_VLocationCenter,                               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  ) 
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'T',                                         &
         LONG_NAME  = 'air_temperature',                           &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PL',                                        &
         LONG_NAME  = 'mid_level_pressure',                        &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'ZLE',                                       &
         LONG_NAME  = 'edge_heights',                              &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'ZL',                                        &
         LONG_NAME  = 'mid_layer_heights',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'S',                                         &
         LONG_NAME  = 'mid_layer_dry_static_energy',               &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PLE',                                       &
         LONG_NAME  = 'edge_pressure',                             &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'TH',                                        &
         LONG_NAME  = 'potential_temperature',                     &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PLK',                                       &
         LONG_NAME  = 'mid_layer_$p^\kappa$',                      &
         UNITS      = 'Pa$^\kappa$',                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'OMEGA',                                     &
         LONG_NAME  = 'vertical_pressure_velocity',                &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                    &
         SHORT_NAME = 'PTFX',                                        &
         LONG_NAME  = 'pressure_weighted_eastward_potential_temperature_flux_unremapped',  &
         UNITS      = 'K Pa m+2 s-1',                                &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,                          &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                    &
         SHORT_NAME = 'PTFY',                                        &
         LONG_NAME  = 'pressure_weighted_northward_potential_temperature_flux_unremapped', &
         UNITS      = 'K Pa m+2 s-1',                                &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,                          &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                    &
         SHORT_NAME = 'MFX_UR',                                      &
         LONG_NAME  = 'pressure_weighted_eastward_wind_unremapped',  &
         UNITS      = 'Pa m+2 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,                          &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                    &
         SHORT_NAME = 'MFY_UR',                                      &
         LONG_NAME  = 'pressure_weighted_northward_wind_unremapped', &
         UNITS      = 'Pa m+2 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,                          &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'MFX',                                       &
         LONG_NAME  = 'pressure_weighted_eastward_wind',           &
         UNITS      = 'Pa m+2 s-1',                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'MFY',                                       &
         LONG_NAME  = 'pressure_weighted_northward_wind',          &
         UNITS      = 'Pa m+2 s-1',                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'MFZ',                                       &
         LONG_NAME  = 'vertical_mass_flux',                        &
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'MFX_A',                                     &
         LONG_NAME  = 'zonal_mass_flux',                           &
         UNITS      = 'Pa m+2 s-1',                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'MFY_A',                                     &
         LONG_NAME  = 'meridional_mass_flux',                      &
         UNITS      = 'Pa m+2 s-1',                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PV',                                        &
         LONG_NAME  = 'ertels_isentropic_potential_vorticity',     &
         UNITS      = 'm+2 kg-1 s-1',                            &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'EPV',                                       &
         LONG_NAME  = 'ertels_potential_vorticity',                &
         UNITS      = 'K m+2 kg-1 s-1',                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'Q',                                         &
         LONG_NAME  = 'specific_humidity',                         &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DUDTANA',                                        &
         LONG_NAME  = 'tendency_of_eastward_wind_due_to_analysis',      &
         UNITS      = 'm s-2',                                      &
         DIMS       = MAPL_DimsHorzVert,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DVDTANA',                                        &
         LONG_NAME  = 'tendency_of_northward_wind_due_to_analysis',     &
         UNITS      = 'm s-2',                                      &
         DIMS       = MAPL_DimsHorzVert,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DTDTANA',                                        &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_analysis',    &
         UNITS      = 'K s-1',                                        &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DDELPDTANA',                                     &
         LONG_NAME  = 'tendency_of_pressure_thickness_due_to_analysis', &
         UNITS      = 'K s-1',                                        &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DUDTDYN',                                   &
         LONG_NAME  = 'tendency_of_eastward_wind_due_to_dynamics', &
         UNITS      = 'm s-2',                                 &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DVDTDYN',                                   &
         LONG_NAME  = 'tendency_of_northward_wind_due_to_dynamics',&
         UNITS      = 'm s-2',                                 &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                    &
         SHORT_NAME = 'DTDTDYN',                                     &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_dynamics', &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DQVDTDYN',                                      &
         LONG_NAME  = 'tendency_of_specific_humidity_due_to_dynamics', &
         UNITS      = 'kg kg-1 s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DQIDTDYN',                                      &
         LONG_NAME  = 'tendency_of_ice_water_due_to_dynamics',         &
         UNITS      = 'kg kg-1 s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DQLDTDYN',                                      &
         LONG_NAME  = 'tendency_of_liquid_water_due_to_dynamics',      &
         UNITS      = 'kg kg-1 s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DOXDTDYN',                                      &
         LONG_NAME  = 'tendency_of_ozone_due_to_dynamics',             &
         UNITS      = 'mol mol-1 s-1',                                 &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PREF',                                      &
         LONG_NAME  = 'reference_air_pressure',                    &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'AK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_a',                   &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'BK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_b',                   &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'PS',                                  &
       LONG_NAME          = 'surface_pressure',                    &
       UNITS              = 'Pa',                                  &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'TA',                                  &
       LONG_NAME          = 'surface_air_temperature',             &
       UNITS              = 'K',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'QA',                                  &
       LONG_NAME          = 'surface_specific_humidity',           &
       UNITS              = 'kg kg-1',                             &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'US',                                  &
       LONG_NAME          = 'surface_eastward_wind',               &
       UNITS              = 'm s-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       FIELD_TYPE         = MAPL_VectorField,                      &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'VS',                                  &
       LONG_NAME          = 'surface_northward_wind',              &
       UNITS              = 'm s-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       FIELD_TYPE         = MAPL_VectorField,                      &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'SPEED',                               &
       LONG_NAME          = 'surface_wind_speed',                  &
       UNITS              = 'm s-1',                               &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'DZ',                                  &
       LONG_NAME          = 'surface_layer_height',                &
       UNITS              = 'm',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'SLP',                                 &
       LONG_NAME          = 'sea_level_pressure',                  &
       UNITS              = 'Pa',                                  &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
       SHORT_NAME         = 'H1000',                               &
       LONG_NAME          = 'height_at_1000_mb',                   &
       UNITS              = 'm',                                   &
       DIMS               = MAPL_DimsHorzOnly,                     &
       VLOCATION          = MAPL_VLocationNone,          RC=STATUS )
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
       UNITS              = 'kg kg-1',                                               &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DELP',                                      &
         LONG_NAME  = 'pressure_thickness',                        &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U_CGRID',                                   &
         LONG_NAME  = 'eastward_wind_on_C-Grid',                   &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V_CGRID',                                   &
         LONG_NAME  = 'northward_wind_on_C-Grid',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U_DGRID',                                   &
         LONG_NAME  = 'eastward_wind_on_native_D-Grid',            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V_DGRID',                                   &
         LONG_NAME  = 'northward_wind_on_native_D-Grid',           &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'TV',                                        &
         LONG_NAME  = 'air_virtual_temperature',                   &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'THV',                                       &
         LONG_NAME  = 'scaled_virtual_potential_temperature',      &
         UNITS      = 'K/Pa$^\kappa$',                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'DDELPDTDYN',                                     &
         LONG_NAME  = 'tendency_of_pressure_thickness_due_to_dynamics', &
         UNITS      = 'Pa s-1',                                       &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'UKE',                                            &
         LONG_NAME  = 'eastward_flux_of_atmospheric_kinetic_energy',    &
         UNITS      = 'J m-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzOnly,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationNone,                    RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'VKE',                                            &
         LONG_NAME  = 'northward_flux_of_atmospheric_kinetic_energy',   &
         UNITS      = 'J m-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzOnly,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationNone,                    RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'UCPT',                                      &
         LONG_NAME  = 'eastward_flux_of_atmospheric_enthalpy',     &
         UNITS      = 'J m-1 s-1',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'VCPT',                                      &
         LONG_NAME  = 'northward_flux_of_atmospheric_enthalpy',    &
         UNITS      = 'J m-1 s-1',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'UPHI',                                           &
         LONG_NAME  = 'eastward_flux_of_atmospheric_potential_energy',  &
         UNITS      = 'J m-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzOnly,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationNone,                    RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                       &
         SHORT_NAME = 'VPHI',                                           &
         LONG_NAME  = 'northward_flux_of_atmospheric_potential_energy', &
         UNITS      = 'J m-1 s-1',                                      &
         DIMS       = MAPL_DimsHorzOnly,                                &
         FIELD_TYPE = MAPL_VectorField,                                 &
         VLOCATION  = MAPL_VLocationNone,                    RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'UQV',                                       &
         LONG_NAME  = 'eastward_flux_of_atmospheric_water_vapor',  &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'VQV',                                       &
         LONG_NAME  = 'northward_flux_of_atmospheric_water_vapor', &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'UQL',                                       &
         LONG_NAME  = 'eastward_flux_of_atmospheric_liquid_water', &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'VQL',                                       &
         LONG_NAME  = 'northward_flux_of_atmospheric_liquid_water',&
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'UQI',                                       &
         LONG_NAME  = 'eastward_flux_of_atmospheric_ice',          &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'VQI',                                       &
         LONG_NAME  = 'northward_flux_of_atmospheric_ice',         &
         UNITS      = 'kg m-1 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DKE',                                       &
         LONG_NAME  = 'tendency_of_atmosphere_kinetic_energy_content_due_to_dynamics',&
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DCPT',                                      &
         LONG_NAME  = 'tendency_of_atmosphere_dry_energy_content_due_to_dynamics',&
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DPET',                                      &
         LONG_NAME  = 'tendency_of_atmosphere_topographic_potential_energy_due_to_dynamics',&
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'WRKT',                                      &
         LONG_NAME  = 'work_done_by_atmosphere_at_top',            &
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DQV',                                       &
         LONG_NAME  = 'tendency_of_atmosphere_water_vapor_content_due_to_dynamics',&
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DQL',                                       &
         LONG_NAME  = 'tendency_of_atmosphere_liquid_water_content_due_to_dynamics',&
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DQI',                                       &
         LONG_NAME  = 'tendency_of_atmosphere_ice_content_due_to_dynamics',&
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'CNV',                                       &
         LONG_NAME  = 'generation_of_atmosphere_kinetic_energy_content',&
         UNITS      = 'W m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U850',                                      &
         LONG_NAME  = 'eastward_wind_at_850_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U500',                                      &
         LONG_NAME  = 'eastward_wind_at_500_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U250',                                      &
         LONG_NAME  = 'eastward_wind_at_250_hPa',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V850',                                      &
         LONG_NAME  = 'northward_wind_at_850_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V500',                                      &
         LONG_NAME  = 'northward_wind_at_500_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V250',                                      &
         LONG_NAME  = 'northward_wind_at_250_hPa',                 &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'T850',                                      &
         LONG_NAME  = 'air_temperature_at_850_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'T500',                                      &
         LONG_NAME  = 'air_temperature_at_500_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'T250',                                      &
         LONG_NAME  = 'air_temperature_at_250_hPa',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'Q850',                                      &
         LONG_NAME  = 'specific_humidity_at_850_hPa',              &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'Q500',                                      &
         LONG_NAME  = 'specific_humidity_at_500_hPa',              &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'Q250',                                      &
         LONG_NAME  = 'specific_humidity_at_250_hPa',              &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'H850',                                      &
         LONG_NAME  = 'height_at_850_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'H500',                                      &
         LONG_NAME  = 'height_at_500_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'H250',                                      &
         LONG_NAME  = 'height_at_250_hPa',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'OMEGA500',                                  &
         LONG_NAME  = 'omega_at_500_hPa',                          &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U50M',                                      &
         LONG_NAME  = 'eastward_wind_at_50_meters',                &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V50M',                                      &
         LONG_NAME  = 'northward_wind_at_50_meters',               &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'AREA',                                      &
         LONG_NAME  = 'agrid_cell_area',                           &
         UNITS      = 'm+2'  ,                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                &
         SHORT_NAME = 'PT',                                        &
         LONG_NAME  = 'scaled_potential_temperature',              &
         UNITS      = 'K Pa$^{-\kappa}$',                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                &
         SHORT_NAME = 'PE',                                        &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

! !INTERNAL STATE:

!ALT: technically the first 2 records of "old" style FV restart have 
!     6 ints: YYYY MM DD H M S
!     5 ints: I,J,K, KS (num true pressure levels), NQ (num tracers) headers

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'AK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_a',                   &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsVertOnly,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'BK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_b',                   &
         UNITS      = '1',                                         &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsVertOnly,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PT',                                        &
         LONG_NAME  = 'scaled_potential_temperature',              &
         UNITS      = 'K Pa$^{-\kappa}$',                          &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PE',                                        &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PKZ',                                       &
         LONG_NAME  = 'pressure_to_kappa',                         &
         UNITS      = 'Pa$^\kappa$',                               &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         RESTART    = MAPL_RestartRequired,                        &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )

!EOS


! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="INITIALIZE"            , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN1"                  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-WRAPPER"              , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--CDCORE"              , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--OMEGA"               , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--BUDGETS"             , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--EPVD"                , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="---PRE_C_CORE"         , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="----PRE_C_CORE_COMM"   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="----C_DELP_LOOP"       , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="---C_CORE"             , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="---C_GEOP"             , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="----YZ_TO_XY_C_GEOP"   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="----XY_TO_YZ_C_GEOP"   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="---PRE_D_CORE"         , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="----PRE_D_CORE_COMM"   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="----C_U_LOOP"          , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="----C_V_PGRAD"         , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="---D_CORE"             , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="---D_GEOP"             , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="----XY_TO_YZ_D_GEOP"   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="----YZ_TO_XY_D_GEOP"   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="---PRE_D_PGRAD"        , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="----PRE_D_PGRAD_COMM_1", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="----D_DELP_LOOP"       , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="---D_PGRAD_1"          , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="---D_PGRAD_2"          , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="---PRE_D_PGRAD_COMM_2" , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--TRAC2D"              , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="---TRAC2D_COMM"        , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="---TRAC2D_TRACER"      , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="----TRAC2D_TRACER_COMM", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--REMAP"               , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_TimerAdd(GC,    name="--BENERGY"             , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--TRANSPOSE_FWD"       , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN2"                  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="FINALIZE"              , RC=STATUS)
    VERIFY_(STATUS)

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_INITIALIZE,  Initialize, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,   Run1, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,   Run2, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_FINALIZE, Finalize, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_READRESTART, Coldstart, rc=status)
    VERIFY_(STATUS)
 
! Generic SetServices
!--------------------

    call MAPL_GenericSetServices( GC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine Initialize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  use g3_dynamics_state_module
  type ( dynamics_grid_type  ) g3_grid
  include 'mpif.h'

  type(ESMF_GridComp), intent(inout) :: gc       ! composite gridded component 
  type(ESMF_State),    intent(inout) :: import   ! import state
  type(ESMF_State),    intent(inout) :: export   ! export state
  type(ESMF_Clock),    intent(inout) :: clock    ! the clock
  
  integer, intent(out), OPTIONAL     :: rc       ! Error code:
                                                 ! = 0 all is well
                                                 ! otherwise, error
  integer                            :: I,J
  type (ESMF_Grid)                   :: grid
  type (ESMF_Config)                 :: cf
  type (ESMF_Config), pointer        :: config

  type (DYN_wrap)                    :: wrap
  type (T_FVDYCORE_STATE),  pointer  :: STATE
  type (T_FVDYCORE_GRID),   pointer  :: FVGRID

  type (MAPL_MetaComp),      pointer :: mapl 

  character (len=ESMF_MAXSTR)        :: restart_file

  type (ESMF_Field)                  :: field
  type (ESMF_Array)                  :: array
  type (ESMF_VM)                     :: VM
  real, pointer                      :: pref(:), ak4(:), bk4(:)
  real(r8), pointer                  :: ak(:), bk(:)
  real(r8), pointer                  ::  pe(:,:,:)
  real(r4), pointer                  :: ple(:,:,:)
  real(r4), pointer                  :: temp2d(:,:)
  character(len=ESMF_MAXSTR)         :: ReplayMode
  real                               :: DNS_INTERVAL
  type (ESMF_TimeInterval)           :: Intv
  type (ESMF_Alarm)                  :: Alarm

  
  integer                            :: status
  character(len=ESMF_MAXSTR)         :: IAm
  character(len=ESMF_MAXSTR)         :: COMP_NAME

  type (ESMF_State)                  :: INTERNAL
  
  real(r8), allocatable              :: tmp2d(:,:)
  integer                            :: ifirstxy, ilastxy, jfirstxy, jlastxy
  integer                            :: im,jm
  integer                            :: NX,NY
  integer                            :: imglobal,jmglobal,lmglobal


! Begin
!------

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Call Generic Initialize
!------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

! Retrieve the pointer to the state
! ---------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

! Start the timers
!-----------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"INITIALIZE")

! Get the private internal state
!-------------------------------

    call ESMF_UserCompGetInternalState(GC, 'FVstate', wrap, status)
    VERIFY_(STATUS)

     state => wrap%dyn_state
    fvgrid => state%grid   ! direct handle to grid

! Set Private Internal State from ESMF internal state in MAPL object
! ------------------------------------------------------------------

    call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)

    call FV_InitState ( STATE, CLOCK, INTERNAL, GC )

! Create PLE and PREF EXPORT Coupling (Needs to be done only once per run)
! ------------------------------------------------------------------------

    call MAPL_GetPointer(EXPORT,PREF,'PREF',ALLOC=.true.,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,AK4 ,'AK'  ,ALLOC=.true.,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,BK4 ,'BK'  ,ALLOC=.true.,RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetPointer(INTERNAL, AK, 'AK', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, BK, 'BK', RC=STATUS)
    VERIFY_(STATUS)

     AK4 = AK
     BK4 = BK
    PREF = AK + BK * P00

    call MAPL_GetPointer(EXPORT,PLE,'PLE',ALLOC=.true.,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL,PE,'PE',RC=STATUS)
    VERIFY_(STATUS)

    PLE = PE

! **********************************************************************
! ****                      Create G3 Grid                          ****
! **********************************************************************

    call MAPL_GetResource( MAPL, NX, 'NX:', default=0, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, NY, 'NY:', default=0, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, imglobal, 'AGCM_IM:', RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, jmglobal, 'AGCM_JM:', RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, lmglobal, 'AGCM_LM:', RC=STATUS )
    VERIFY_(STATUS)

    call create_dynamics_lattice ( g3_grid%lattice,nx,ny )
    call   init_dynamics_lattice ( g3_grid%lattice,mpi_comm_world,imglobal,jmglobal,lmglobal )
    call create_dynamics_grid    ( g3_grid,imglobal,jmglobal,lmglobal )
    call   init_dynamics_grid    ( g3_grid,imglobal,jmglobal,lmglobal,0,state%grid%ak,state%grid%bk )

! Compute Grid-Cell Area
! ----------------------
    call MAPL_GetPointer(export,temp2d,'AREA', ALLOC=.true., rc=status)
    VERIFY_(STATUS)

    ifirstxy = fvgrid%ifirstxy
    ilastxy  = fvgrid%ilastxy
    jfirstxy = fvgrid%jfirstxy
    jlastxy  = fvgrid%jlastxy

    ALLOCATE( tmp2d(ifirstxy:ilastxy,jfirstxy:jlastxy) )

    do j=MAX(2,jfirstxy),MIN(jlastxy,jmglobal-1)
       tmp2d(:,j) = fvgrid%dl*fvgrid%cosp(j)*RADIUS * fvgrid%dp*RADIUS
    enddo
    if ( jfirstxy == 1  ) then
         j=1
         tmp2d(:,j) = fvgrid%acap*(  fvgrid%dl*RADIUS * fvgrid%dp*RADIUS)/imglobal
    endif
    if ( jlastxy  == jmglobal ) then
         j=jmglobal
         tmp2d(:,j) = fvgrid%acap*(  fvgrid%dl*RADIUS * fvgrid%dp*RADIUS)/imglobal
    endif
    temp2d = tmp2d

    DEALLOCATE( tmp2d )

! ======================================================================
!ALT: the next section addresses the problem when export variables have been
!     assigned values during Initialize. To prevent "connected" exports
!     being overwritten by DEFAULT in the Import spec in the other component
!     we label them as being "initailized by restart". A better solution
!     would be to move the computation to phase 2 of Initialize and
!     eliminate this section alltogether
! ======================================================================
    call ESMF_StateGet(EXPORT, 'PREF', FIELD, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", &
                           VALUE=MAPL_InitialRestart, RC=STATUS)
    VERIFY_(STATUS)      

    call ESMF_StateGet(EXPORT, 'PLE', FIELD, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_AttributeSet(field, NAME="MAPL_InitStatus", &
                           VALUE=MAPL_InitialRestart, RC=STATUS)
    VERIFY_(STATUS)      

!=====Begin intemittent replay=======================

! Set the intermittent replay alarm, if needed.
! Note that it is a non-sticky alarm
! and is set to ringing on first step. So it will
! work whether the clock is backed-up ans ticked
! or not.

    call MAPL_GetResource(MAPL, ReplayMode, 'REPLAY_MODE:', default="NoReplay", RC=STATUS )
    VERIFY_(STATUS)

    if(adjustl(ReplayMode)=="Intermittent") then
       call MAPL_GetResource(MAPL, DNS_INTERVAL,'REPLAY_INTERVAL:', default=21600., RC=STATUS )
       VERIFY_(STATUS)
       call ESMF_TimeIntervalSet(Intv, S=nint(DNS_INTERVAL), RC=STATUS)
       VERIFY_(STATUS)

       ALARM = ESMF_AlarmCreate(name='INTERMITTENT', clock=CLOCK,      &
                                ringInterval=Intv, sticky=.false.,    &
                                                            RC=STATUS )
       VERIFY_(STATUS)
       call ESMF_AlarmRingerOn(ALARM, rc=status)
       VERIFY_(STATUS)
    end if

!========End intermittent replay========================

    call MAPL_TimerOff(MAPL,"INITIALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 subroutine FV_InitState (STATE, CLOCK, INTERNAL, GC)

  use dynamics_vars, only : dynamics_init
  use parutilitiesmodule, only : gsize, gid, parinit

  type (T_FVDYCORE_STATE),pointer              :: STATE

  type (ESMF_Clock), target,     intent(INOUT) :: CLOCK
  type (ESMF_GridComp)         , intent(INout) :: GC
  type (ESMF_State)            , intent(INOUT) :: INTERNAL

! Local variables

  type (ESMF_TimeInterval)     :: Time2Run
  type (ESMF_TimeInterval)     :: CheckMaxMin
  type (ESMF_VM)               :: VM
  type (T_FVDYCORE_GRID) , pointer :: GRID
  integer              :: rc
  integer              :: status
  integer              :: len
  real(r8) :: REAL_PACK(6)

  integer :: NPRXY_X, NPRXY_Y, NPRYZ_Y, NPRYZ_Z, &
             DT, IORD, JORD, KORD, TE_METHOD, NSPLIT
  integer :: force_2d, geopktrans
  integer, allocatable :: IMXY(:), JMXY(:), JMYZ(:), KMYZ(:)

  integer :: nx, ny
  integer :: nstep, nymd, nhms
  integer :: yr, mm, dd, h, m, s, itmp
  integer :: INT_PACK(6)

  type(ESMF_DELayout) :: layoutYZ
  integer             :: I, nDEs
  integer             :: img
  integer             :: jmg
  integer             :: kmg

  integer   :: im, jm, km         !  Global dims
  integer   :: nq                 !  No. advected tracers
  integer   :: ntotq              !  No. total tracers
  integer   :: ks                 !  True # press. levs
  integer   :: ifirstxy, ilastxy  !  Interval
  integer   :: jfirstxy, jlastxy  !  Interval
  integer   :: jfirst, jlast      !  Interval
  integer   :: kfirst, klast      !  Interval
  integer   :: k                  !  Vertical loop index
  integer   :: srcCellCountPerDim(3), srcStartPerDEPerDim(gsize,3)

  character(len=ESMF_MAXSTR)       :: IAm='FV:Init_State'

  real(r8), pointer                   :: AK(:), BK(:)
  real(r8), dimension(:,:,:), pointer :: U, V, PT, PE, PKZ
  type (MAPL_MetaComp),       pointer :: mapl 
  integer                             :: comm

  real  ple,ples,plet,sig,dpl

! Retrieve the pointer to the state
! ---------------------------------

  call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
  VERIFY_(STATUS)

! Save the mapl state for FVperf_module
! -------------------------------------

  STATE%GRID%FVgenstate => MAPL
                   GRID => STATE%GRID     ! For convenience

! Initialize Layout based on 2-D decomposition
! --------------------------------------------

  call MAPL_GetResource( MAPL, NX, 'NX:', default=0, RC=STATUS )
  VERIFY_(STATUS)
  call MAPL_GetResource( MAPL, NY, 'NY:', default=0, RC=STATUS )
  VERIFY_(STATUS)

  NPRXY_X = NX
  NPRXY_Y = NY
  NPRYZ_Y = NY
  NPRYZ_Z = NX

  _ASSERT( NPRXY_X>0 .AND. NPRXY_Y>0          ,'needs informative message')
  _ASSERT( NPRYZ_Y>0 .AND. NPRYZ_Z>0          ,'needs informative message')
  _ASSERT( NPRXY_X*NPRXY_Y == NPRYZ_Y*NPRYZ_Z ,'needs informative message')

  call MAPL_GetResource( MAPL, force_2d, 'force_2d:', default=0, RC=STATUS )
  VERIFY_(STATUS)

! Get the layout and store directly in the GRID data structure
! ------------------------------------------------------------

  grid%twod_decomp = 1

  if (npryz_z .eq. 1 .and. nprxy_x .eq. 1 .and. force_2d .eq. 0) then
      grid%twod_decomp = 0
      call WRITE_PARALLEL('Code operating with 1D decomposition')
  endif

! Pilgrim initialization: pass the 2D decomposition and other parameters for FV optimization
! ------------------------------------------------------------------------------------------

   call ESMF_VMGetCurrent(vm, rc=rc)
   call ESMF_VMGet(vm, mpiCommunicator=comm, rc=rc)
   call parinit( comm=comm, npryzxy = (/ npryz_y, npryz_z, nprxy_x, nprxy_y/), &
                 mod_method  = grid%mod_method,                     &
                 mod_geopk   = grid%mod_geopk,                      &
                 mod_gatscat = grid%mod_gatscat )
           
! Get Global Dimensions
! ---------------------
  call MAPL_GetResource( MAPL, IMG, 'AGCM_IM:', default=0, RC=STATUS )
  VERIFY_(STATUS)
  call MAPL_GetResource( MAPL, JMG, 'AGCM_JM:', default=0, RC=STATUS )
  VERIFY_(STATUS)
  call MAPL_GetResource( MAPL, KMG, 'AGCM_LM:', default=0, RC=STATUS )
  VERIFY_(STATUS)

! Create IMXY, JMXY, JMYZ, KMYZ vectors
! -------------------------------------

  allocate( imxy(0:nprxy_x-1) )  
  allocate( jmxy(0:nprxy_y-1) )  
  allocate( jmyz(0:npryz_y-1) )  
  allocate( kmyz(0:npryz_z-1) )  

  call MAPL_DecomposeDim ( img,imxy,nprxy_x )
  call MAPL_DecomposeDim ( jmg,jmxy,nprxy_y )
  call MAPL_DecomposeDim ( jmg,jmyz,npryz_y )
  call MAPL_DecomposeDim ( kmg,kmyz,npryz_z )

! Get other scalars
! -----------------

  call MAPL_GetResource( MAPL, dt,   'RUN_DT:', default=0, RC=STATUS )
  VERIFY_(STATUS)
  call MAPL_GetResource( MAPL, iord,   'iord:', default=3, RC=STATUS )
  VERIFY_(STATUS)
  call MAPL_GetResource( MAPL, jord,   'jord:', default=3, RC=STATUS )
  VERIFY_(STATUS)
  call MAPL_GetResource( MAPL, kord,   'kord:', default=4, RC=STATUS )
  VERIFY_(STATUS)

! Vertical Remapping Method for Total Energy (default=1 is cubic interpolation)
! -----------------------------------------------------------------------------

  call MAPL_GetResource( MAPL, te_method, 'te_method:', default=1, RC=STATUS )
  VERIFY_(STATUS)

! Ratio of Large/Small Timesteps (default=0 implies automatic calculation)
! ------------------------------------------------------------------------

  call MAPL_GetResource( MAPL, nsplit, 'nsplit:', default=0, RC=STATUS )
  VERIFY_(STATUS)

! Heritage Code for Tracers
! -------------------------
  ntotq = 1        ! Total Number of Tracers
     nq = ntotq    ! Total Number of Advected Tracers

! Other assertions
!
  _ASSERT(maxval(IMXY)>0 .AND. maxval(JMXY)>0,'needs informative message')
  _ASSERT(maxval(JMYZ)>0 .AND. maxval(KMYZ)>0,'needs informative message')
  _ASSERT(DT > 0.0                           ,'needs informative message')

  call WRITE_PARALLEL('Dynamics PE Layout')
  call WRITE_PARALLEL(IMG        ,format='("IM_Global: ",(   I4))')
  call WRITE_PARALLEL(JMG        ,format='("JM_Global: ",(   I4))')
  call WRITE_PARALLEL(KMG        ,format='("LM_Global: ",(   I4))')
  call WRITE_PARALLEL(NPRXY_X    ,format='("NPRXY_X  : ",(   I4))')
  call WRITE_PARALLEL(NPRXY_Y    ,format='("NPRXY_Y  : ",(   I4))')
  call WRITE_PARALLEL(NPRYZ_Y    ,format='("NPRYZ_Y  : ",(   I4))')
  call WRITE_PARALLEL(NPRYZ_Z    ,format='("NPRYZ_Z  : ",(   I4))')
  call WRITE_PARALLEL(IMXY(0:NPRXY_X-1),format='("IMXY : ",(256I3))')
  call WRITE_PARALLEL(JMXY(0:NPRXY_Y-1),format='("JMXY : ",(256I3))')
  call WRITE_PARALLEL(JMYZ(0:NPRYZ_Y-1),format='("JMYZ : ",(256I3))')
  call WRITE_PARALLEL(KMYZ(0:NPRYZ_Z-1),format='("KMYZ : ",(256I3))')

  call WRITE_PARALLEL(iord,format='(/,"IORD: ",(I2))')
  call WRITE_PARALLEL(jord,format='(  "JORD: ",(I2))')
  call WRITE_PARALLEL(kord,format='(  "KORD: ",(I2))')
  call WRITE_PARALLEL(te_method,format='(  "TE_METHOD: ",(I2),/)')

! These are run-specific variables:  
!     DT              Time step
!     IORD            Order (mode) of X interpolation (1,..,6)
!     JORD            Order (mode) of Y interpolation (1,..,6)
!     NSPLIT          Ratio of big to small timestep (set to zero if in doubt)
!

  STATE%DOTIME    = .TRUE.
  STATE%CHECK_DT  = 21600.   ! Check max and min of arrays every 6 hours.
  STATE%DT        = DT
  STATE%IORD      = IORD
  STATE%JORD      = JORD
  STATE%KORD      = KORD
  STATE%TE_METHOD = TE_METHOD

! Calculation of orders for the C grid is fixed by D-grid IORD, JORD
!-------------------------------------------------------------------

  if( iord <= 2 ) then
    STATE%ICD =  1
  else
    STATE%ICD = -2
  endif

  if( jord <= 2 ) then
    STATE%JCD =  1
  else
    STATE%JCD = -2
  endif

  call WRITE_PARALLEL(STATE%DT,format='("Dynamics time step: ",(F10.4))')


! Get the main GRIDXY grid from the application (no longer set in this module)
!-----------------------------------------------------------------------------

  call ESMF_GridCompGet(gc, grid=GRID%GRIDXY, vm=vm, rc=STATUS)

! Get size, grid, and coordinate specifications
!----------------------------------------------

!MJS: we should get these from the MAPL object

  im  = SUM(IMXY)
  jm  = SUM(JMXY)
  km  = SUM(KMYZ)

! Calculate NSPLIT if it was specified as 0
! -----------------------------------------
  if ( NSPLIT == 0 ) then
       STATE%NSPLIT = INIT_NSPLIT(STATE%DT,IM,JM)
  else
       STATE%NSPLIT = NSPLIT
       call WRITE_PARALLEL(STATE%NSPLIT,format='("Dynamics    NSPLIT: ",(I3),/)')
  endif

  call WRITE_PARALLEL((/im,jm,km/)       , &
    format='("Resolution of dynamics restart     =",3I5)'  )

  ks = 0 ! ALT: this was the value when we read "old" style FV_internal restart
         !      if needed, we could compute, ks by count(BK==0.0)
         !      then FV will try to run slightly more efficient code
         !      So far, GEOS-5 has used ks = 0
  _ASSERT(ks <= KM+1,'needs informative message')
  call WRITE_PARALLEL(ks                          , &
     format='("Number of true pressure levels =", I5)'   )

!
! Make sure that IM, JM, KM are the sums of the (exclusive) dist.
!
  _ASSERT(jm == SUM(JMYZ),'needs informative message')

!
!
! Note: it is necessary to create GRIDXY and GRIDYZ now in 
!       order to access the first and last local indices.
!       This makes it difficult to cleanly separate the 
!       grid initialization into init_fvdycore_grid.
!
!

!ALT???
! we need to check if the grid is OK

  GRID%GRIDYZ = ESMF_GridCreate(             &
            name="FVCORE_YZ_grid",           &
            countsPerDEDim1=JMYZ,            &
            countsPerDEDim2=KMYZ,            &
            indexFlag = ESMF_INDEX_GLOBAL, &
            coordDep1 = (/1,2/),             &
            coordDep2 = (/1,2/),             &
            gridEdgeLWidth = (/0,0/),    &
            gridEdgeUWidth = (/0,0/),    &
            rc=status)
       VERIFY_(STATUS)


  call MAPL_GRID_INTERIOR(GRID%GRIDXY, ifirstxy, ilastxy, &
                                       jfirstxy, jlastxy )

  call MAPL_GRID_INTERIOR(GRID%GRIDYZ, jfirst, jlast, &
                                       kfirst, klast )

! Get pointers to internal state vars
!------------------------------------

  call MAPL_GetPointer(internal, ak,  "AK", rc=status)
  VERIFY_(STATUS) 
  call MAPL_GetPointer(internal, bk,  "BK", rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, u,   "U",  rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, v,   "V",  rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, pt,  "PT", rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, pe,  "PE", rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, pkz, "PKZ",rc=status)
  VERIFY_(STATUS)

!
! WS: CREATE_VARS moved here to define STATE%VARS soon enough
!
  call CREATE_VARS ( ifirstxy, ilastxy,      &
                     jfirstxy, jlastxy,      &
                     1, km, km+1,            &
                     U, V, PT, PE, PKZ,      &
                     STATE%VARS  )


! Report
!-------
 
      if( gid.eq.0 ) then
        print *
        write(6,100)
100     format(2x,' k ','      A(k)    ',2x,'   B(k)   ',2x,'  Pref    ',2x,'  DelP   ',2x,'  Sige  ',/, &
               1x,'----',3x,'----------',2x,'----------',2x,'----------',2x,'---------',2x,'--------' )
            k=0
        plet = ak(k )*0.01 + 1000.0*bk(k)
        ple  = ak(k )*0.01 + 1000.0*bk(k)
        ples = ak(km)*0.01 + 1000.0*bk(km)
        write(6,101) k+1,ak(k)*0.01, bk(k), ple
        do k=1,km
        dpl = ple
        ple = ak(k)*0.01 + 1000.0*bk(k)
        dpl =  ple-dpl
        sig = (ple-plet)/(ples-plet)
        write(6,102) k+1,ak(k)*0.01, bk(k), ple, dpl, sig
        enddo
 
        print *
101     format(2x,i3,2x,f10.6,2x,f10.6,2x,f10.4)
102     format(2x,i3,2x,f10.6,2x,f10.6,2x,f10.4,3x,f8.4,2x,f10.6)
      endif

! Initialize the FVDYCORE static variables in the GRID
!-----------------------------------------------------

  call dynamics_init( state%dt, state%jord, im, jm, km,       &
                      PI, RADIUS, OMEGA, nq, ntotq, ks,       &
                      ifirstxy,ilastxy, jfirstxy,jlastxy,     &
                      jfirst,  jlast,   kfirst,  klast,       &
                      nprxy_x, nprxy_y, npryz_y, npryz_z,     &
                      imxy,    jmxy,    jmyz,    kmyz,        &
                      ak,      bk,      0,    grid )

  STATE%CLOCK => CLOCK

  call ESMF_TimeIntervalSet(Time2Run, &
                            S=nint(STATE%DT), rc=status)
  VERIFY_(status)

  STATE%ALARMS(TIME_TO_RUN) = ESMF_AlarmCreate(name="Time2Run", clock=clock, &
                              ringInterval=Time2Run, &
                              Enabled=.TRUE., rc=status)
  VERIFY_(status)

  call ESMF_AlarmEnable  (STATE%ALARMS(TIME_TO_RUN), rc=status); VERIFY_(status)
  call ESMF_AlarmRingerOn(STATE%ALARMS(TIME_TO_RUN), rc=status); VERIFY_(status)

  call ESMF_TimeIntervalSet(CheckMaxMin, S=nint(STATE%CHECK_DT), rc=status)
  VERIFY_(status)

  STATE%ALARMS(CHECK_MAXMIN) = ESMF_AlarmCreate(name="CheckMaxMin", clock=clock, &
                               RingInterval=CheckMaxMin, &
                               Enabled=.TRUE., rc=status)
  VERIFY_(status)

  call WRITE_PARALLEL(' ')

  call WRITE_PARALLEL(STATE%DT, &
    format='("INITIALIZED ALARM: DYN_TIME_TO_RUN EVERY   ",F9.1," secs.")')
  call WRITE_PARALLEL(STATE%CHECK_DT, &
    format='("INITIALIZED ALARM: CHECK MAX AND MIN EVERY ",F9.1," secs.")')

  return

contains

!-----------------------------------------------------------------------
! BOP
! !IROUTINE:  init_nsplit --- find proper value for nsplit if not specified
!
! !INTERFACE:
  integer function INIT_NSPLIT(dtphy,im,jm) 
!
! !USES:
    implicit none

! !INPUT PARAMETERS:
    real (r8), intent(in) :: dtphy      !  Physics Time Step
    integer,   intent(in) :: im, jm     !  Global horizontal resolution

! !DESCRIPTION:
!
!    If nsplit=0 (module variable) then determine a good value 
!    for ns (used in dynpkg) based on resolution and the large-time-step 
!    (pdt). The user may have to set this manually if instability occurs.
!
! !REVISION HISTORY:
!   00.10.19   Lin     Creation
!   01.03.26   Sawyer  ProTeX documentation
!   01.06.10   Sawyer  Modified for dynamics_init framework
!   03.12.04   Sawyer  Moved here from dynamics_vars.  Now a function
!
! EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
    
    integer dtdyn

      dtdyn = nint(18*dtphy/im)

      do while ( mod(dtphy,real(dtdyn,kind=8)).ne.0 )
      dtdyn = dtdyn - 1
      if( dtdyn.lt.1 ) then
          print *
          print *, 'Cannot determine Dynamics Timestep'
          print *
          stop
      endif
      enddo

    init_nsplit = dtphy/dtdyn
    
    call WRITE_PARALLEL ( init_nsplit ,format='("Dynamics    NSPLIT: ",(I3),/)' )
    
    return
  end function INIT_NSPLIT
!---------------------------------------------------------------------
  subroutine CREATE_VARS (I1, IN, J1, JN, K1, KN, KP, &
       U, V, PT, PE, PKZ, VARS )

    integer, intent(IN   ) :: I1, IN, J1, JN, K1, KN, KP
    real(r8), target ::   U(I1:IN,J1:JN,K1:KN  )
    real(r8), target ::   V(I1:IN,J1:JN,K1:KN  )
    real(r8), target ::  PT(I1:IN,J1:JN,K1:KN  )
    real(r8), target ::  PE(I1:IN,J1:JN,K1:KP  )
    real(r8), target :: PKZ(I1:IN,J1:JN,K1:KN  )
    
    type (T_FVDYCORE_VARS), intent(INOUT) :: VARS
    
    VARS%U => U
    VARS%V => V
    VARS%PT => PT
    VARS%PE => PE
    VARS%PKZ => PKZ

    return
  end subroutine CREATE_VARS



end subroutine FV_INITSTATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






  subroutine Run1(gc, import, export, clock, rc)
   use dynamics_vars, only : c2a3d

    use g3_dynamics_state_module
    type ( dynamics_state_type ) dynamics
    save                         dynamics

    type(ESMF_GridComp), intent(inout) :: gc
    type (ESMF_State),   intent(inout) :: import
    type (ESMF_State),   intent(inout) :: export
    type (ESMF_Clock),   intent(INout) :: clock
    integer, intent(out), optional     :: rc 

    include 'mpif.h'
    integer  nx,ny,ierror
    integer  imglobal,jmglobal
    real(r8) pi,dl,dp
    real(r4) alpha
    integer  nsplit
    character*4 scheme

! !Local Variables:
  
    integer                                          :: status
    type (ESMF_FieldBundle)                          :: bundle
    type (ESMF_FieldBundle)                          :: DNS_Bundle
    type (ESMF_Field)                                :: field
    type (ESMF_Field)                                :: dns_field
    type (ESMF_Config)                               :: cf
    type (ESMF_Alarm)                                :: Alarm
    type (ESMF_Grid)                                 :: ESMFGRID
    type (ESMF_Time)                                 :: currentTime

    type (MAPL_MetaComp), pointer :: mapl 

    type (DYN_wrap) :: wrap
    type (T_FVDYCORE_STATE), pointer :: STATE
    type (T_FVDYCORE_GRID),  pointer :: GRID
    type (T_FVDYCORE_VARS),  pointer :: VARS
    
    integer  :: J1, JN, K1, KN, NQ, KQ
    integer  :: IM, JM, KM
    integer  :: NKE, NPHI
    integer  :: NUMVARS
    integer  :: ifirstxy, ilastxy, jfirstxy, jlastxy
    integer  :: I, J, K, L, n, pos
    logical, parameter :: convt = .false. ! Until this is run with full physics
    logical  :: is_ringing
    logical  first
    data     first /.true./

    real(r8),     pointer :: phisxy(:,:)
    real(kind=4), pointer ::   phis(:,:)

    real(r8), allocatable ::    pke(:,:,:) ! pe**kappa
    real(r8), allocatable ::    pk (:,:,:) ! mid-level pressure

    real(r8), allocatable ::   pkxy(:,:,:) ! pe**kappa
    real(r8), allocatable ::     pl(:,:,:) ! mid-level pressure
    real(r8), allocatable :: tempxy(:,:,:) ! mid-level temperature
    real(r8), allocatable ::     ua(:,:,:) ! temporary array
    real(r8), allocatable ::     va(:,:,:) ! temporary array
    real(r8), allocatable ::     qv(:,:,:) ! temporary array
    real(r8), allocatable ::     ql(:,:,:) ! temporary array
    real(r8), allocatable ::     qi(:,:,:) ! temporary array
    real(r8), allocatable ::  qdnew(:,:,:) ! temporary array
    real(r8), allocatable ::  qdold(:,:,:) ! temporary array
    real(r8), allocatable ::  qvold(:,:,:) ! temporary array
    real(r8), allocatable ::  qlold(:,:,:) ! temporary array
    real(r8), allocatable ::  qiold(:,:,:) ! temporary array
    real(r8), allocatable ::     ox(:,:,:) ! temporary array
    real(r8), allocatable ::     zl(:,:,:) ! temporary array
    real(r8), allocatable ::    zle(:,:,:) ! temporary array
    real(r8), allocatable ::   delp(:,:,:) ! temporary array
    real(r8), allocatable ::   dudt(:,:,:) ! temporary array
    real(r8), allocatable ::   dvdt(:,:,:) ! temporary array
    real(r8), allocatable ::   dtdt(:,:,:) ! temporary array
    real(r8), allocatable ::   dqdt(:,:,:) ! temporary array
    real(r8), allocatable ::  dthdt(:,:,:) ! temporary array
    real(r8), allocatable ::  ddpdt(:,:,:) ! temporary array
    real(r8), allocatable ::     dmdt(:,:) ! temporary array
    real(r8), allocatable ::    gze(:,:,:) ! temporary array

    real(r8), allocatable, target :: ke    (:,:,:) ! Kinetic    Energy
    real(r8), allocatable, target :: cpt   (:,:,:) ! Internal   Energy
    real(r8), allocatable, target :: phi   (:,:,:) ! Potential  Energy
    real(r8), allocatable :: qsum1 (:,:)   ! Vertically Integrated Kinetic   Energy Tracer
    real(r8), allocatable :: qsum2 (:,:)   ! Vertically Integrated Internal  Energy Tracer

    real(r8), allocatable :: phi00 (:,:)   ! Vertically Integrated phi
    real(r8), allocatable :: penrg (:,:)   ! Vertically Integrated Cp*T
    real(r8), allocatable :: kenrg (:,:)   ! Vertically Integrated K
    real(r8), allocatable :: tenrg (:,:)   ! PHIS*(Psurf-Ptop)
    real(r8), allocatable :: penrg0(:,:)   ! Vertically Integrated Cp*T
    real(r8), allocatable :: kenrg0(:,:)   ! Vertically Integrated K
    real(r8), allocatable :: tenrg0(:,:)   ! PHIS*(Psurf-Ptop)
    real(r8), allocatable :: penrga(:,:)   ! Vertically Integrated Cp*T
    real(r8), allocatable :: kenrga(:,:)   ! Vertically Integrated K
    real(r8), allocatable :: tenrga(:,:)   ! PHIS*(Psurf-Ptop)
    real(r8), allocatable :: penrgb(:,:)   ! Vertically Integrated Cp*T
    real(r8), allocatable :: kenrgb(:,:)   ! Vertically Integrated K
    real(r8), allocatable :: tenrgb(:,:)   ! PHIS*(Psurf-Ptop)
    real(r8), allocatable :: kehot (:,:)   ! Vertically Integrated K due to higher-order-terms
    real(r8), allocatable :: kedp  (:,:)   ! Vertically Integrated K due to pressure change
    real(r8), allocatable :: keadv (:,:)   ! Vertically Integrated K due to advection
    real(r8), allocatable :: kepg  (:,:)   ! Vertically Integrated K due to pressure gradient

    real(r8), allocatable :: kegen  (:,:)
    real(r8), allocatable :: kedyn  (:,:)
    real(r8), allocatable :: pedyn  (:,:)
    real(r8), allocatable :: tedyn  (:,:)
    real(r8), allocatable :: kecdcor(:,:)
    real(r8), allocatable :: pecdcor(:,:)
    real(r8), allocatable :: tecdcor(:,:)
    real(r8), allocatable :: keremap(:,:)
    real(r8), allocatable :: peremap(:,:)
    real(r8), allocatable :: teremap(:,:)
    real(r8), allocatable :: convke (:,:)
    real(r8), allocatable :: convcpt(:,:)
    real(r8), allocatable :: convphi(:,:)
    real(r8), allocatable :: convthv(:,:)

    real(r8), allocatable :: DNS_phis(:,:)
    real(r8), allocatable :: DNS_thv (:,:,:)

    real(r8),     allocatable :: dthdtremap  (:,:)   ! Vertically Integrated THV tendency due to vertical remapping
    real(r8),     allocatable :: dthdtconsv  (:,:)   ! Vertically Integrated THV tendency due to TE conservation
    real(kind=4), allocatable :: dqvdtanaint1(:,:)
    real(kind=4), allocatable :: dqvdtanaint2(:,:)
    real(kind=4), allocatable :: dqldtanaint1(:,:)
    real(kind=4), allocatable :: dqldtanaint2(:,:)
    real(kind=4), allocatable :: dqidtanaint1(:,:)
    real(kind=4), allocatable :: dqidtanaint2(:,:)
    real(kind=4), allocatable :: doxdtanaint1(:,:)
    real(kind=4), allocatable :: doxdtanaint2(:,:)
    real(kind=4), allocatable :: dthdtanaint1(:,:)
    real(kind=4), allocatable :: dthdtanaint2(:,:)

    real(kind=4), allocatable :: dummy (:,:,:) ! Dummy 3-D  Variable
    real(kind=4), allocatable :: tropp1(:,:)   ! Tropopause Pressure
    real(kind=4), allocatable :: tropp2(:,:)   ! Tropopause Pressure
    real(kind=4), allocatable :: tropp3(:,:)   ! Tropopause Pressure
    real(kind=4), allocatable :: tropt (:,:)   ! Tropopause Temperature
    real(kind=4), allocatable :: tropq (:,:)   ! Tropopause Specific Humidity

    real(r8), allocatable :: pelnxz(:,:,:) ! log pressure (pe) at layer edges
    real(r8), allocatable :: omaxyz(:,:,:) ! vertical pressure velocity (pa/sec)
    real(r8), allocatable :: cptxyz(:,:,:) ! Cp*Tv
    real(r8), allocatable :: thvxyz(:,:,:) ! Thetav
    real(r8), allocatable :: epvxyz(:,:,:) ! ertel's potential vorticity
    real(r8), allocatable :: cxxyz(:,:,:)  ! Accumulated zonal winds
    real(r8), allocatable :: cyxyz(:,:,:)  ! Accumulated meridional winds
    real(r8), allocatable :: ptfxxyz(:,:,:) ! zonal mass-weighted PT flux
    real(r8), allocatable :: ptfyxyz(:,:,:) ! meridional mass-weighted PT flux
    real(r8), allocatable :: mfxxyz_ur(:,:,:) ! zonal mass flux
    real(r8), allocatable :: mfyxyz_ur(:,:,:) ! meridional mass flux
    real(r8), allocatable :: mfxxyz(:,:,:) ! zonal mass flux
    real(r8), allocatable :: mfyxyz(:,:,:) ! meridional mass flux
    real(r8), allocatable :: mfzxyz(:,:,:) ! vertical mass flux
    real(r8), allocatable :: mfxxyz_a(:,:,:) ! zonal mass flux A-Grid
    real(r8), allocatable :: mfyxyz_a(:,:,:) ! meridional mass flux A-Grid
    real(r8)              :: dt            ! Dynamics time step
    real(r8)              :: kinetic       ! local kinetic   energy
    real(r8)              :: potential     ! local potential energy
    real(r8)              :: dtmp          ! Temperature Change due to CONSV=TRUE
    real(r8), allocatable :: tempr8(:,:,:) ! Cp*Tv
    real(r8), allocatable :: trsum1(:)     ! Global Sum of Tracers before Add_Incs
    real(r8), allocatable :: trsum2(:)     ! Global Sum of Tracers after  Add_Incs

    real(kind=4), pointer ::      dudtana(:,:,:)
    real(kind=4), pointer ::      dvdtana(:,:,:)
    real(kind=4), pointer ::      dtdtana(:,:,:)
    real(kind=4), pointer ::     ddpdtana(:,:,:)
    real(kind=4), pointer ::       dqldt (:,:,:)
    real(kind=4), pointer ::       dqidt (:,:,:)
    real(kind=4), pointer ::       doxdt (:,:,:)
    real(kind=4), pointer ::      dqvana (:,:,:)
    real(kind=4), pointer ::      dqlana (:,:,:)
    real(kind=4), pointer ::      dqiana (:,:,:)
    real(kind=4), pointer ::      doxana (:,:,:)
    real(kind=4), pointer ::       temp3d(:,:,:)
    real(kind=4), pointer ::       temp2d(:,:)
    real(kind=4), pointer ::       tempu (:,:)
    real(kind=4), pointer ::       tempv (:,:)

    character(len=ESMF_MAXSTR), ALLOCATABLE       :: NAMES (:)
    character(len=ESMF_MAXSTR), ALLOCATABLE, save :: NAMES0(:)
    character(len=ESMF_MAXSTR) :: IAm
    character(len=ESMF_MAXSTR) :: COMP_NAME
    character(len=ESMF_MAXSTR) :: STRING
    character(len=ESMF_MAXSTR) :: ReplayFile
    character(len=ESMF_MAXSTR) :: ReplayMode

    type(T_TRACERS)            :: qqq       ! Specific Humidity
    type(T_TRACERS)            :: ooo       ! OX
    integer                    :: NXQ       ! Number of Additional Budget Tracers

    type (MAPL_SunOrbit)       :: ORBIT
    real(kind=4), pointer      :: LATS(:,:)
    real(kind=4), pointer      :: LONS(:,:)
    real(kind=4), allocatable  ::  ZTH(:,:)
    real(kind=4), allocatable  ::  SLR(:,:)

    logical LCONSV, LFILL
    integer  CONSV,  FILL

  Iam = "Run1"
  call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, Grid=ESMFGRID, RC=STATUS )
  VERIFY_(STATUS)
  Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the generic state
! -----------------------------------------

  call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
  VERIFY_(STATUS)

  call MAPL_TimerOn(MAPL,"TOTAL")
  call MAPL_TimerOn(MAPL,"RUN1")

! Retrieve the pointer to the internal state
! ------------------------------------------

  call ESMF_UserCompGetInternalState(gc, 'FVstate', wrap, status)
  VERIFY_(STATUS)
  state => wrap%dyn_state

  vars  => state%vars   ! direct handle to control variables
  grid  => state%grid   ! direct handle to grid
  dt    =  state%dt     ! dynamics time step (large)

  ifirstxy = grid%ifirstxy
  ilastxy  = grid%ilastxy
  jfirstxy = grid%jfirstxy
  jlastxy  = grid%jlastxy

! im       = grid%im
! jm       = grid%jm
  km       = grid%km


  is_ringing = ESMF_AlarmIsRinging( STATE%ALARMS(TIME_TO_RUN),rc=status); VERIFY_(status) 
  if (.not. is_ringing) return

! Allocate Arrays
! ---------------
      ALLOCATE( tempr8(ifirstxy:ilastxy,jfirstxy:jlastxy,2 ) )
      ALLOCATE(  dummy(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   delp(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   dudt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   dvdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   dtdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(   dqdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  dthdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  ddpdt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE( tempxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     pl(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     ua(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     va(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     qv(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     ql(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     qi(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qdnew(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qdold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qvold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qlold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(  qiold(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(     ox(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )

      ALLOCATE(     ke(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(    cpt(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(    phi(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      ALLOCATE(    gze(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )

      ALLOCATE(  qsum1(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(  qsum2(ifirstxy:ilastxy,jfirstxy:jlastxy)    )

      ALLOCATE(   dmdt(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(  phi00(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(  kenrg(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(  penrg(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE(  tenrg(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( kenrg0(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( penrg0(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( tenrg0(ifirstxy:ilastxy,jfirstxy:jlastxy)    )

      ALLOCATE( kepg  (ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( keadv (ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( kedp  (ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( kehot (ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( kenrga(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( penrga(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( tenrga(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( kenrgb(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( penrgb(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      ALLOCATE( tenrgb(ifirstxy:ilastxy,jfirstxy:jlastxy)    )

      ALLOCATE( kegen  (ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( kedyn  (ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( pedyn  (ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( tedyn  (ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( kecdcor(ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( pecdcor(ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( tecdcor(ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( keremap(ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( peremap(ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( teremap(ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( convke (ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( convcpt(ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( convphi(ifirstxy:ilastxy,jfirstxy:jlastxy)   )
      ALLOCATE( convthv(ifirstxy:ilastxy,jfirstxy:jlastxy)   )

      ALLOCATE( tropp1      (ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( tropp2      (ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( tropp3      (ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( tropt       (ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( tropq       (ifirstxy:ilastxy,jfirstxy:jlastxy) )

      ALLOCATE( dqvdtanaint1(ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( dqvdtanaint2(ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( dqldtanaint1(ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( dqldtanaint2(ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( dqidtanaint1(ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( dqidtanaint2(ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( doxdtanaint1(ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( doxdtanaint2(ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( dthdtanaint1(ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( dthdtanaint2(ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( dthdtremap  (ifirstxy:ilastxy,jfirstxy:jlastxy) )
      ALLOCATE( dthdtconsv  (ifirstxy:ilastxy,jfirstxy:jlastxy) )

      ALLOCATE( pelnxz   (ifirstxy:ilastxy,km+1,jfirstxy:jlastxy) )

      ALLOCATE( phisxy   (ifirstxy:ilastxy,jfirstxy:jlastxy     ) )
      ALLOCATE(   pkxy   (ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE(     zl   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE(    zle   (ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE( omaxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( cptxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( thvxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( epvxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE(  cxxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE(  cyxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( ptfxxyz  (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( ptfyxyz  (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfxxyz_ur(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfyxyz_ur(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfxxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfyxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfzxyz   (ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      ALLOCATE( mfxxyz_a (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      ALLOCATE( mfyxyz_a (ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )


! Report advected friendlies
!---------------------------

      call ESMF_StateGet ( IMPORT, 'TRADV' , BUNDLE,   RC=STATUS )
      VERIFY_(STATUS)
      call ESMF_FieldBundleGet ( BUNDLE, fieldCount=NQ, RC=STATUS )
      VERIFY_(STATUS)

      allocate( NAMES(NQ),STAT=STATUS )
      VERIFY_(STATUS)
      call ESMF_FieldBundleGet ( BUNDLE, itemorderflag=ESMF_ITEMORDER_ADDORDER, fieldNameList=NAMES, rc=STATUS )
      VERIFY_(STATUS)

      if( .not.allocated( names0 ) ) then
           allocate( NAMES0(NQ),STAT=STATUS )
           VERIFY_(STATUS)
           write(STRING,'(A,I5,A)') "Advecting the following ", nq, " tracers in FV:"
              call WRITE_PARALLEL( trim(STRING)   )
           do k=1,nq
              call WRITE_PARALLEL( trim(NAMES(K)) )
           end do
           NAMES0 = NAMES
      endif

     !if( size(names0).ne.size(names) ) then
     !     deallocate( NAMES0 )
     !       allocate( NAMES0(NQ),STAT=STATUS )
     !     VERIFY_(STATUS)
     !     write(STRING,'(A,I5,A)') "Advecting the following ", nq, " tracers in FV:"
     !        call WRITE_PARALLEL( trim(STRING)   )
     !     do k=1,nq
     !        call WRITE_PARALLEL( trim(NAMES(K)) )
     !     end do
     !     NAMES0 = NAMES
     !endif

! Surface Geopotential from IMPORT state
!---------------------------------------

      call MAPL_GetPointer ( IMPORT, PHIS, 'PHIS', RC=STATUS )
      VERIFY_(STATUS)

      phisxy = real(phis,kind=r8)


! Set Addition Tracers for Exact Budget Diagnostics
!--------------------------------------------------

      call MAPL_GetPointer ( export,temp2D,'KEHOT',rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
         NXQ = 2
      else
         NXQ = 0
      endif


! Get tracers from IMPORT State (Note: Contains Updates from Analysis)
!---------------------------------------------------------------------

      call PULL_Q ( STATE, IMPORT, qqq, NXQ, rc )

      do k=1,size(names)
         pos = index(names(k),'::')
         if(pos > 0) then
           if( (names(k)(pos+2:))=='OX' ) then
            ooo = vars%tracer(k)
           endif
         endif
         if( trim(names(k))=='Q'  ) then
            qqq = vars%tracer(k)
             kq = k
         endif
      enddo

! If requested, do Intermittent Replay
!-------------------------------------

      call MAPL_GetResource(MAPL, ReplayMode, 'REPLAY_MODE:', default="NoReplay", RC=STATUS )
      VERIFY_(STATUS)

      REPLAYING: if(adjustl(ReplayMode)=="Intermittent") then

! It is an error not to specify a replay file at this point.
!-----------------------------------------------------------

         call MAPL_GetResource ( MAPL,ReplayFile,'REPLAY_FILE:', RC=STATUS )
         VERIFY_(status)

! If replay alarm is ringing, we need to reset state
!---------------------------------------------------

         call ESMF_ClockGetAlarm(Clock,'INTERMITTENT',Alarm,rc=Status)
         VERIFY_(status) 

         is_ringing = ESMF_AlarmIsRinging( Alarm,rc=status )
         VERIFY_(status) 

         TIME_TO_REPLAY: if(is_ringing) then

            ALLOCATE( DNS_phis(ifirstxy:ilastxy,jfirstxy:jlastxy    ) )
            ALLOCATE( DNS_thv (ifirstxy:ilastxy,jfirstxy:jlastxy,km ) )

! Read the fields to be reset into a bundle
!------------------------------------------

            DNS_Bundle = ESMF_FieldBundleCreate( RC=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldBundleSet(DNS_bundle, grid=ESMFGRID, RC=STATUS)
            VERIFY_(STATUS)

            call ESMF_ClockGet(CLOCK,      CurrTime=currentTIME,    RC=STATUS)
            VERIFY_(STATUS)
            call MAPL_CFIORead(ReplayFile, currentTime, DNS_Bundle,  &
!                ONLY_VARS='uwnd,vwnd,delp,ozone,sphu,theta,phis',   &
                                                            RC=STATUS)
            VERIFY_(STATUS)

! Fill the state variables from the bundle only if
!  the corresponding fields are there.
!-------------------------------------------------

! U
            call ESMFL_BundleGetPointertoData(DNS_Bundle,'uwnd',TEMP3D, RC=STATUS)
!            VERIFY_(STATUS)
            if(STATUS==ESMF_SUCCESS) then
               if(grid%iam==0) print *, 'Replaying uwnd'
               vars%u = TEMP3D
            endif
! V
            call ESMFL_BundleGetPointertoData(DNS_Bundle,'vwnd',TEMP3D, RC=STATUS)
!            VERIFY_(STATUS)
            if(STATUS==ESMF_SUCCESS) then
               if(grid%iam==0) print *, 'Replaying vwnd'
               vars%v = TEMP3D 
            endif
! PE
            call ESMFL_BundleGetPointertoData(DNS_Bundle,'delp',TEMP3D, RC=STATUS)
!            VERIFY_(STATUS)
            DNS_PRESSURE: if(STATUS==ESMF_SUCCESS) then
               if(grid%iam==0) print *, 'Replaying ple'
               vars%pe(:,:,1) = grid%ak(1)
               do k=2,km+1
                  vars%pe(:,:,k) = vars%pe(:,:,k-1) + temp3d(:,:,k-1)
               enddo
            end if DNS_PRESSURE
! O3
            call ESMFL_BundleGetPointertoData(DNS_Bundle,'ozone',TEMP3D, RC=STATUS)
!            VERIFY_(STATUS)
            DNS_OZONE: if(STATUS==ESMF_SUCCESS) then

!   Ozone needs to be adjusted to OX
!-----------------------------------
               if(grid%iam==0) print *, 'Replaying ozone'
               
               call MAPL_Get(MAPL, LONS=LONS, LATS=LATS, ORBIT=ORBIT, RC=STATUS )
               VERIFY_(STATUS)

               allocate( ZTH( size(LONS,1),size(LONS,2) ) )
               allocate( SLR( size(LONS,1),size(LONS,2) ) )

               call MAPL_SunGetInsolation( LONS,LATS,ORBIT,ZTH,SLR, CLOCK=CLOCK,RC=STATUS  )
               VERIFY_(STATUS)

               pl = ( vars%pe(:,:,2:) + vars%pe(:,:,:km) ) * 0.5

               do L=1,km
                  if( ooo%is_r4 ) then 
                     where(PL(:,:,L) >= 100.0 .or. ZTH <= 0.0) &
                          ooo%content_r4(:,:,L) =  TEMP3D(:,:,L)*1.0E-6
                  else
                     where(PL(:,:,L) >= 100.0 .or. ZTH <= 0.0) &
                          ooo%content   (:,:,L) =  TEMP3D(:,:,L)*1.0E-6
                  end if
               enddo

               deallocate( ZTH, SLR )

            end if DNS_OZONE
! QV
            call ESMFL_BundleGetPointertoData(DNS_Bundle,'sphu',TEMP3D, RC=STATUS)
!            VERIFY_(STATUS)
            DSN_HUMIDITY: if(STATUS==ESMF_SUCCESS) then
               if(grid%iam==0) print *, 'Replaying sphu'
               if( qqq%is_r4 ) then 
                  qqq%content_r4 = TEMP3D
               else
                  qqq%content    = TEMP3D
               endif
            end if DSN_HUMIDITY
! PT
            call ESMFL_BundleGetPointertoData(DNS_Bundle,'theta',TEMP3D, RC=STATUS)
!            VERIFY_(STATUS)
            DSN_THETAV: if(STATUS==ESMF_SUCCESS) then
               if(grid%iam==0) print *, 'Replaying thetav'
               DNS_thv = TEMP3D 
            else
               if( qqq%is_r4 ) then 
                  DNS_thv = vars%pt*(1.0+eps*qqq%content_r4)
               else
                  DNS_thv = vars%pt*(1.0+eps*qqq%content   )
               endif
            end if DSN_THETAV

! If there is a topo in the file, remap fields
!---------------------------------------------

            call ESMFL_BundleGetPointertoData(DNS_Bundle,'phis',TEMP2D, RC=STATUS)
!            VERIFY_(STATUS)
            if(STATUS==ESMF_SUCCESS) then
               if(grid%iam==0) print *, 'Remapping ...'
               DNS_phis = TEMP2D
               call remap ( vars%pe, vars%u, vars%v, DNS_thv, vars%tracer, DNS_phis, phisxy, &
                            grid%ak, grid%bk, size(DNS_thv,1), size(DNS_thv,2), km, nq )
            end if

            if( qqq%is_r4 ) then 
               vars%pt = dns_thv/(1.0+eps*qqq%content_r4)
            else
               vars%pt = dns_thv/(1.0+eps*qqq%content   )
            endif

            pkxy = vars%pe**kappa
            do k=1,km
               vars%pkz(:,:,k) = ( pkxy(:,:,k+1)-pkxy(:,:,k) ) &
                    / ( kappa*( log(vars%pe(:,:,k+1))-log(vars%pe(:,:,k)) ) )
            enddo

! Done with replay; clean-up
!---------------------------

            call ESMF_FieldBundleGet(DNS_Bundle , FieldCount=NUMVARS,      RC=STATUS)
            VERIFY_(STATUS)

            do k=1,NUMVARS
               call ESMF_FieldBundleGet (DNS_Bundle, k, DNS_FIELD,    RC=STATUS)
               VERIFY_(STATUS)
               call MAPL_FieldDestroy   (DNS_Field,                   RC=STATUS)
               VERIFY_(STATUS)
            end do

            call ESMF_FieldBundleDestroy(DNS_Bundle,                       RC=STATUS)
            VERIFY_(STATUS)

            DEALLOCATE( DNS_phis )
            DEALLOCATE( DNS_thv  )

         end if TIME_TO_REPLAY
      end if REPLAYING

! Create Local Copy of QV and OX (Contains Updates from Analysis)
!----------------------------------------------------------------

                  ox = 0.0d0  ! Initialize in case no OX advection
      do k=1,size(names)
         pos = index(names(k),'::')
         if(pos > 0) then
           if( (names(k)(pos+2:))=='OX' ) then
             if ( ooo%is_r4 ) then
                  ox = ooo%content_r4
             else
                  ox = ooo%content
             endif
           endif
         endif
         if( trim(names(k))=='Q'  ) then
             if ( qqq%is_r4 ) then
                  qv = qqq%content_r4
             else
                  qv = qqq%content
             endif
         endif
      enddo

! G3 Initialization
! -----------------
  if( first ) then

      call MAPL_GetResource( MAPL, NX,       'NX:',      RC=STATUS ) ; VERIFY_(STATUS)
      call MAPL_GetResource( MAPL, NY,       'NY:',      RC=STATUS ) ; VERIFY_(STATUS)
      call MAPL_GetResource( MAPL, imglobal, 'AGCM_IM:', RC=STATUS ) ; VERIFY_(STATUS)
      call MAPL_GetResource( MAPL, jmglobal, 'AGCM_JM:', RC=STATUS ) ; VERIFY_(STATUS)

      call create_dynamics_lattice ( dynamics%grid%lattice,nx,ny )
      call   init_dynamics_lattice ( dynamics%grid%lattice,mpi_comm_world,imglobal,jmglobal,km )

      im = dynamics%grid%lattice%im( dynamics%grid%lattice%pei )
      jm = dynamics%grid%lattice%jm( dynamics%grid%lattice%pej )

! Create G3 Dynamics State
! ------------------------
      call create_dynamics ( dynamics,im,jm,km,nq )

! Initialize G3 Dynamics Grid
! ---------------------------
      call init_dynamics_grid ( dynamics%grid,imglobal,jmglobal,km,nq,grid%ak,grid%bk )

      first = .false.
  endif

      im = dynamics%grid%im
      jm = dynamics%grid%jm

! Diagnostics Before Analysis Increments are Added
!-------------------------------------------------

      call MAPL_GetPointer ( IMPORT, dqvana, 'DQVANA', RC=STATUS )   ! Get QV Increment from Analysis
      VERIFY_(STATUS)
      call MAPL_GetPointer ( IMPORT, dqlana, 'DQLANA', RC=STATUS )   ! Get QL Increment from Analysis
      VERIFY_(STATUS)
      call MAPL_GetPointer ( IMPORT, dqiana, 'DQIANA', RC=STATUS )   ! Get QI Increment from Analysis
      VERIFY_(STATUS)
      call MAPL_GetPointer ( IMPORT, doxana, 'DOXANA', RC=STATUS )   ! Get OX Increment from Analysis
      VERIFY_(STATUS)

      QL = 0.0
      QI = 0.0
      do N = 1,size(names)
           if( trim(names(N)).eq.'QLCN' .or. &
               trim(names(N)).eq.'QLLS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then
                     QL = QL + state%vars%tracer(N)%content_r4
                 else
                     QL = QL + state%vars%tracer(N)%content
                 endif
           endif
           if( trim(names(N)).eq.'QICN' .or. &
               trim(names(N)).eq.'QILS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then
                     QI = QI + state%vars%tracer(N)%content_r4
                 else
                     QI = QI + state%vars%tracer(N)%content
                 endif
           endif
      enddo
      QVOLD = QV-DQVANA
      QLOLD = QL-DQLANA
      QIOLD = QI-DQIANA

      QDOLD = 1.0 - (QVOLD+QLOLD+QIOLD)
      QDNEW = 1.0 - (QV   +QL   +QI   )

      call ctoa_winds ( vars%u,vars%v,ua,va,                                                 &
                        dynamics%grid%dlam,dynamics%grid%dphi,im,jm,km,dynamics%grid%lattice )

      delp   = vars%pe(:,:,2:)  -vars%pe(:,:,:km)   ! Pressure Thickness
      dmdt   = vars%pe(:,:,km+1)-vars%pe(:,:,1)     ! Psurf-Ptop
      tempxy = vars%pt * (1.0+eps*(qv-dqvana))      ! Compute THV Before Analysis Update
 
      call Energetics (state,vars%u,vars%v,tempxy,vars%pe,delp,vars%pkz,phisxy,kenrg,penrg,tenrg,dynamics%grid)

! DUDTANA
! -------
      call MAPL_GetPointer ( export, dudtana, 'DUDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(dudtana) ) dudtana = ua

! DVDTANA
! -------
      call MAPL_GetPointer ( export, dvdtana, 'DVDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(dvdtana) ) dvdtana = va

! DTDTANA
! -------
      call MAPL_GetPointer ( export, dtdtana, 'DTDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(dtdtana) ) dtdtana = vars%pt * vars%pkz

! DDELPDTANA
! ----------
      call MAPL_GetPointer ( export, ddpdtana, 'DDELPDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(ddpdtana) ) ddpdtana = delp

! DTHVDTANAINT
! ------------
      call MAPL_GetPointer ( export, temp2D, 'DTHVDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          tempxy       = vars%pt*(1+eps*(qv-dqvana))   ! Set tempxy = TH*QVold (Before Analysis Update)
          dthdtanaint1 = 0.0
          do k=1,km
          dthdtanaint1 = dthdtanaint1 + tempxy(:,:,k)*delp(:,:,k)
          enddo
      endif

! DQVDTANAINT
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DQVDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          tempxy       = qv-dqvana   ! Set tempxy = QVold (Before Analysis Update)
          dqvdtanaint1 = 0.0
          do k=1,km
          dqvdtanaint1 = dqvdtanaint1 + tempxy(:,:,k)*delp(:,:,k)
          enddo
      endif

! DQLDTANAINT
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DQLDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dqldtanaint1 = 0.0
          do N = 1,size(names)
             if( trim(names(N)).eq.'QLCN' .or. &
                 trim(names(N)).eq.'QLLS' ) then
                 do k=1,km
                 if( state%vars%tracer(N)%is_r4 ) then 
                     dqldtanaint1 = dqldtanaint1 + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                 else
                     dqldtanaint1 = dqldtanaint1 + state%vars%tracer(N)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo
             endif
          enddo
          do k=1,km
             dqldtanaint1 = dqldtanaint1 - dqlana(:,:,k)*delp(:,:,k)
          enddo
      endif

! DQIDTANAINT
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DQIDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dqidtanaint1 = 0.0
          do N = 1,size(names)
             if( trim(names(N)).eq.'QICN' .or. &
                 trim(names(N)).eq.'QILS' ) then
                 do k=1,km
                 if( state%vars%tracer(N)%is_r4 ) then 
                     dqidtanaint1 = dqidtanaint1 + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                 else
                     dqidtanaint1 = dqidtanaint1 + state%vars%tracer(N)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo
             endif
          enddo
          do k=1,km
             dqidtanaint1 = dqidtanaint1 - dqiana(:,:,k)*delp(:,:,k)
          enddo
      endif

! DOXDTANAINT
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DOXDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          tempxy       = ox-doxana   ! Set tempxy = OXold (Before Analysis Update)
          doxdtanaint1 = 0.0
          do k=1,km
          doxdtanaint1 = doxdtanaint1 + tempxy(:,:,k)*delp(:,:,k)
          enddo
      endif

! Add Diabatic Forcing from Analysis to State Variables
! -----------------------------------------------------

      allocate( trsum1(nq) )
      allocate( trsum2(nq) )

      ! Compute Global Mass of Aerosol Constituents Before ANA Updates
      ! --------------------------------------------------------------
      call glosum   ( STATE,NQ,TRSUM1 )

      call ADD_INCS ( STATE,IMPORT,DT,dynamics%grid )

! Update Specific Mass of Aerosol Constituents Keeping Mixing_Ratio Constant WRT_Dry_Air After ANA Updates
! --------------------------------------------------------------------------------------------------------
      do n=1,NQ
      if( (trim(names(n)).ne.'Q'   ) .and. &
          (trim(names(n)).ne.'QLLS') .and. &
          (trim(names(n)).ne.'QLCN') .and. &
          (trim(names(n)).ne.'QILS') .and. &
          (trim(names(n)).ne.'QICN') .and. &
          (trim(names(n)).ne.'CLLS') .and. &
          (trim(names(n)).ne.'CLCN')       ) then
           if( STATE%VARS%TRACER(N)%IS_R4 ) then
               state%vars%tracer(n)%content_r4 = state%vars%tracer(n)%content_r4 * ( QDNEW/QDOLD )
           else
               state%vars%tracer(n)%content    = state%vars%tracer(n)%content    * ( QDNEW/QDOLD )
           endif
      endif
      enddo

      ! Compute Global Mass of Aerosol Constituents After ANA Updates
      ! -------------------------------------------------------------
      call glosum   ( STATE,NQ,TRSUM2 )

      ! Ensure Conservation of Global Mass of Aerosol Constituents After ANA Updates
      ! ----------------------------------------------------------------------------
      do n=1,NQ
      if( (trim(names(n)).ne.'Q'   ) .and. &
          (trim(names(n)).ne.'QLLS') .and. &
          (trim(names(n)).ne.'QLCN') .and. &
          (trim(names(n)).ne.'QILS') .and. &
          (trim(names(n)).ne.'QICN') .and. &
          (trim(names(n)).ne.'CLLS') .and. &
          (trim(names(n)).ne.'CLCN')       ) then

           if( trsum2(n).ne.0.0d0 ) then
               trsum2(n) = trsum1(n)/trsum2(n)
           else
               trsum2(n) = 1.0d0
           endif
          !IF (MAPL_AM_I_ROOT()) print *, trim(names(n)),' ratio is: ',trsum2(n)

           if( STATE%VARS%TRACER(N)%IS_R4 ) then
               state%vars%tracer(n)%content_r4 = state%vars%tracer(n)%content_r4 * trsum2(n)
           else
               state%vars%tracer(n)%content    = state%vars%tracer(n)%content    * trsum2(n)
           endif
      endif
      enddo

      deallocate( trsum1 )
      deallocate( trsum2 )

! Update Local Copy of QV and OX to account for Global Sum Adjustment
!--------------------------------------------------------------------

      do k=1,size(names)
         pos = index(names(k),'::')
         if(pos > 0) then
           if( (names(k)(pos+2:))=='OX' ) then
             if ( ooo%is_r4 ) then
                  ox = ooo%content_r4
             else
                  ox = ooo%content
             endif
           endif
         endif
         if( trim(names(k))=='Q'  ) then
             if ( qqq%is_r4 ) then
                  qv = qqq%content_r4
             else
                  qv = qqq%content
             endif
         endif
      enddo

! Diagnostics After Analysis Increments are Added
!------------------------------------------------

      call MAPL_GetPointer ( export, temp2D, 'DMDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) temp2D = ( (vars%pe(:,:,km+1)-vars%pe(:,:,1)) - dmdt )/(grav*dt)

      call ctoa_winds ( vars%u,vars%v,ua,va,                                                 &
                        dynamics%grid%dlam,dynamics%grid%dphi,im,jm,km,dynamics%grid%lattice )

      delp = vars%pe(:,:,2:)  -vars%pe(:,:,:km)   ! Pressure Thickness
      dmdt = vars%pe(:,:,km+1)-vars%pe(:,:,1)     ! Psurf-Ptop

! DUDTANA
! -------
      call MAPL_GetPointer ( export, dudtana, 'DUDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(dudtana) ) then
                     dummy   =  ua
                     dudtana = (dummy-dudtana)/dt
      endif

! DVDTANA
! -------
      call MAPL_GetPointer ( export, dvdtana, 'DVDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(dvdtana) ) then
                     dummy   =  va
                     dvdtana = (dummy-dvdtana)/dt
      endif

! DTDTANA
! -------
      call MAPL_GetPointer ( export, dtdtana, 'DTDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(dtdtana) ) then
                     dummy   =  vars%pt*vars%pkz
                     dtdtana = (dummy-dtdtana)/dt
      endif

! DDELPDTANA
! ----------
      call MAPL_GetPointer ( export, ddpdtana, 'DDELPDTANA', rc=status )
      VERIFY_(STATUS)
      if( associated(ddpdtana) ) then
                     dummy    =  delp
                     ddpdtana = (dummy-ddpdtana)/dt
      endif

! DTHVDTANAINT
! ------------
      call MAPL_GetPointer ( export, temp2D, 'DTHVDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          tempxy       = vars%pt*(1+eps*qv)   ! Set tempxy = TH*QVnew (After Analysis Update)
          dthdtanaint2 = 0.0
          do k=1,km
          dthdtanaint2 = dthdtanaint2 + tempxy(:,:,k)*delp(:,:,k)
          enddo
          temp2D       = (dthdtanaint2-dthdtanaint1) * MAPL_P00**MAPL_KAPPA / (MAPL_GRAV*DT)
      endif

! DQVDTANAINT
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DQVDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          tempxy       = qv         ! Set tempxy = QNEW (After Analysis Update)
          dqvdtanaint2 = 0.0
          do k=1,km
          dqvdtanaint2 = dqvdtanaint2 + tempxy(:,:,k)*delp(:,:,k)
          enddo
          temp2D       = (dqvdtanaint2-dqvdtanaint1) / (MAPL_GRAV*DT)
      endif

! DQLDTANAINT
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DQLDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dqldtanaint2 = 0.0
          do N = 1,size(names)
             if( trim(names(N)).eq.'QLCN' .or. &
                 trim(names(N)).eq.'QLLS' ) then
                 do k=1,km
                 if( state%vars%tracer(N)%is_r4 ) then 
                     dqldtanaint2 = dqldtanaint2 + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                 else
                     dqldtanaint2 = dqldtanaint2 + state%vars%tracer(N)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo
             endif
          enddo
          temp2D = (dqldtanaint2-dqldtanaint1) / (MAPL_GRAV*DT)
      endif

! DQIDTANAINT
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DQIDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dqidtanaint2 = 0.0
          do N = 1,size(names)
             if( trim(names(N)).eq.'QICN' .or. &
                 trim(names(N)).eq.'QILS' ) then
                 do k=1,km
                 if( state%vars%tracer(N)%is_r4 ) then 
                     dqidtanaint2 = dqidtanaint2 + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                 else
                     dqidtanaint2 = dqidtanaint2 + state%vars%tracer(N)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo
             endif
          enddo
          temp2D = (dqidtanaint2-dqidtanaint1) / (MAPL_GRAV*DT)
      endif

! DOXDTANAINT
! -----------
      call MAPL_GetPointer ( export, temp2D, 'DOXDTANAINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          tempxy       = ox         ! Set tempxy = OXnew (After Analysis Update)
          doxdtanaint2 = 0.0
          do k=1,km
          doxdtanaint2 = doxdtanaint2 + tempxy(:,:,k)*delp(:,:,k)
          enddo
          temp2D       = (doxdtanaint2-doxdtanaint1) * (MAPL_O3MW/MAPL_AIRMW) / (MAPL_GRAV*DT)
      endif

! Create FV Thermodynamic Variables
!----------------------------------

      tempxy = vars%pt * vars%pkz      ! Compute Dry Temperature
     vars%pt = vars%pt * (1.0+eps*qv)  ! Compute Virtual Potential Temperature

! Initialize Diagnostic Dynamics Tendencies
! -----------------------------------------

      ddpdt  =   delp       ! Pressure Thickness Tendency
      dudt   =     ua       ! U-Wind on A-Grid   Tendency
      dvdt   =     va       ! V-Wind on A-Grid   Tendency
      dtdt   = tempxy       ! Dry Temperature    Tendency
      dqdt   =     qv       ! Specific Humidity  Tendency

! Initialize 3-D Tracer Dynamics Tendencies
! -----------------------------------------

      call MAPL_GetPointer( export,dqldt,'DQLDTDYN', rc=status )
      VERIFY_(STATUS)
      if( associated(dqldt) ) then
          dqldt = 0.0
          do N = 1,size(names)
             if( trim(names(N)).eq.'QLCN' .or. &
                 trim(names(N)).eq.'QLLS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then 
                     dqldt = dqldt - state%vars%tracer(N)%content_r4
                 else
                     dqldt = dqldt - state%vars%tracer(N)%content
                 endif
             endif
          enddo
      endif

      call MAPL_GetPointer( export,dqidt,'DQIDTDYN', rc=status )
      VERIFY_(STATUS)
      if( associated(dqidt) ) then
          dqidt = 0.0
          do N = 1,size(names)
             if( trim(names(N)).eq.'QICN' .or. &
                 trim(names(N)).eq.'QILS' ) then    
                 if( state%vars%tracer(N)%is_r4 ) then 
                     dqidt = dqidt - state%vars%tracer(N)%content_r4
                 else
                     dqidt = dqidt - state%vars%tracer(N)%content
                 endif
             endif
          enddo
      endif

      call MAPL_GetPointer( export,doxdt,'DOXDTDYN', rc=status )
      VERIFY_(STATUS)
      if( associated(doxdt) ) then
          doxdt = 0.0
          do N = 1,size(names)
             pos = index(names(N),'::')
             if(pos > 0) then
               if( (names(N)(pos+2:))=='OX' ) then
                 if( state%vars%tracer(N)%is_r4 ) then 
                     doxdt = doxdt - state%vars%tracer(N)%content_r4
                 else
                     doxdt = doxdt - state%vars%tracer(N)%content
                 endif
               endif
             endif
          enddo
      endif

! Initialize 2-D Vertically Integrated Tracer Dynamics Tendencies
! ---------------------------------------------------------------

      call MAPL_GetPointer ( export, temp2D, 'DQVDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          tempr8(:,:,1) = 0.0
          do k=1,km
          tempr8(:,:,1) = tempr8(:,:,1) + qv(:,:,k)*delp(:,:,k)
          enddo
      endif

      call MAPL_GetPointer ( export, temp2D, 'DQLDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          temp2d = 0.0
          do N = 1,size(names)
             if( trim(names(N)).eq.'QLCN' .or. &
                 trim(names(N)).eq.'QLLS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then 
                     do k=1,km
                     temp2d = temp2d - state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                 else
                     do k=1,km
                     temp2d = temp2d - state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                 endif
             endif
          enddo
      endif

      call MAPL_GetPointer ( export, temp2D, 'DQIDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          temp2d = 0.0
          do N = 1,size(names)
             if( trim(names(N)).eq.'QICN' .or. &
                 trim(names(N)).eq.'QILS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then 
                     do k=1,km
                     temp2d = temp2d - state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                 else
                     do k=1,km
                     temp2d = temp2d - state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                 endif
             endif
          enddo
      endif

      call MAPL_GetPointer ( export, temp2D, 'DOXDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          temp2d = 0.0
          do N = 1,size(names)
             pos = index(names(N),'::')
             if(pos > 0) then
               if( (names(N)(pos+2:))=='OX' ) then
                 if( state%vars%tracer(N)%is_r4 ) then 
                     do k=1,km
                     temp2d = temp2d - state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                 else
                     do k=1,km
                     temp2d = temp2d - state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                 endif
               endif
             endif
          enddo
      endif


! Compute Energetics After Analysis (and Before Dycore)
! -----------------------------------------------------

    call Energetics (state,vars%u,vars%v,vars%pt,vars%pe,delp,vars%pkz,phisxy, &
                     kenrg0,penrg0,tenrg0,dynamics%grid,ke=ke,cpt=cpt,gze=gze)

    kenrg = (kenrg0-kenrg)/DT
    penrg = (penrg0-penrg)/DT
    tenrg = (tenrg0-tenrg)/DT

    call FILLOUT2 (export, 'KEANA', kenrg, rc=status); VERIFY_(STATUS)
    call FILLOUT2 (export, 'PEANA', penrg, rc=status); VERIFY_(STATUS)
    call FILLOUT2 (export, 'TEANA', tenrg, rc=status); VERIFY_(STATUS)

! Add Passive Tracers for KE and PHI+CpT
! --------------------------------------

      nq = STATE%GRID%NQ

      phi00 = 0.0
      do k=1,km
       phi(:,:,k) = ( gze(:,:,k+1)*vars%pe(:,:,k+1)-gze(:,:,k)*vars%pe(:,:,k) )/delp(:,:,k) + (1+kappa)*cpt(:,:,k)
      phi00 = phi00 + phi(:,:,k)*delp(:,:,k)
      enddo
      phi00 = phi00 / grav

      if( NXQ.eq.2 ) then
          NKE  = nq-1
          NPHI = nq
          state%vars%tracer(NKE )%content => KE
          state%vars%tracer(NPHI)%content => PHI
          state%vars%tracer(NKE )%is_r4 = .false.
          state%vars%tracer(NPHI)%is_r4 = .false.
      else
          NKE  = -999
          NPHI = -999
      endif

      dthdt = vars%pt*delp

! Clear mass fluxes
!------------------

      ptfxxyz  (:,:,:) = 0.
      ptfyxyz  (:,:,:) = 0.
      mfxxyz_ur(:,:,:) = 0.
      mfyxyz_ur(:,:,:) = 0.
      mfxxyz   (:,:,:) = 0.
      mfyxyz   (:,:,:) = 0.
      mfzxyz   (:,:,:) = 0.
      mfxxyz_a (:,:,:) = 0.
      mfyxyz_a (:,:,:) = 0.

! Call Wrapper for FVDycore
! -------------------------
      call MAPL_GetResource( MAPL, CONSV, 'CONSV:', default=1, RC=STATUS )
      VERIFY_(STATUS)
      call MAPL_GetResource( MAPL,  FILL,  'FILL:', default=0, RC=STATUS )
      VERIFY_(STATUS)

      LCONSV = CONSV.eq.1
      LFILL  =  FILL.eq.1

! Load FV Variables into G3 Dynamics State
! ----------------------------------------
      dynamics%vars(1)%p(:,:)   = ( state%vars%pe(:,:,km+1) - dynamics%grid%ptop  )  ! Convert to PI=(PS-PTOP) in Pa
      dynamics%vars(1)%u(:,:,:) =   state%vars%u (:,:,:)
      dynamics%vars(1)%v(:,:,:) =   state%vars%v (:,:,:)

      allocate( pke(im,jm,km+1) )
      allocate( pk (im,jm,km  ) )

      call getpke( dynamics%vars(1)%p,pke   ,dynamics%grid,im,jm )
      call getpk ( dynamics%vars(1)%p,pke,pk,dynamics%grid,im,jm )

      dynamics%vars(1)%t(:,:,:) = tempxy(:,:,:) / pk(:,:,:)

      do n=1,dynamics%grid%ntracer
         if( state%vars%tracer(n)%is_r4 ) then
               dynamics%vars(1)%q(:,:,:,n) = state%vars%tracer(n)%content_r4(:,:,:)
         else
               dynamics%vars(1)%q(:,:,:,n) = state%vars%tracer(n)%content   (:,:,:)
         endif
      enddo

! Initialize TRSUMs for Conservation of Tracers
! ---------------------------------------------

      allocate( trsum1(nq) )
      allocate( trsum2(nq) )

      call glosum ( STATE,NQ,TRSUM1 )

! Call Dynamics Wrapper
! ---------------------
      call MAPL_TimerOn (MAPL,"-WRAPPER")

      call MAPL_GetResource( MAPL, scheme, 'TSCHEME:', DEFAULT='LEAP',RC=STATUS ) ; VERIFY_(STATUS)
      call MAPL_GetResource( MAPL,  alpha,   'ALPHA:', DEFAULT=0.05,  RC=STATUS ) ; VERIFY_(STATUS)

      nsplit = state%nsplit
                       dynamics%vars(1)%t(:,:,:) = dynamics%vars(1)%t(:,:,:)*(1.0+eps*dynamics%vars(1)%q(:,:,:,kq))
      call g3_wrapper (dynamics,phisxy,scheme,dt,nsplit,alpha,omaxyz)
                       dynamics%vars(1)%t(:,:,:) = dynamics%vars(1)%t(:,:,:)/(1.0+eps*dynamics%vars(1)%q(:,:,:,kq))

      vars%pe(:,:,km+1) = ( dynamics%vars(1)%p(:,:) + dynamics%grid%ptop )  ! Construct PS = PI+PTOP in Pa

      pelnxz     = 0
      cptxyz     = 0
      thvxyz     = 0
      epvxyz     = 0
      cxxyz      = 0
      cyxyz      = 0
      ptfxxyz    = 0
      ptfyxyz    = 0
      mfxxyz_ur  = 0
      mfyxyz_ur  = 0
      mfxxyz     = 0
      mfyxyz     = 0
      mfzxyz     = 0
      kenrga     = 0
      penrga     = 0
      tenrga     = 0
      kenrgb     = 0
      penrgb     = 0
      tenrgb     = 0
      keadv      = 0
      kepg       = 0
      kedp       = 0
      kehot      = 0
      dthdtremap = 0
      dthdtconsv = 0
      dtmp       = 0

      call MAPL_TimerOff(MAPL,"-WRAPPER")

! Copy G3 Dynamics State back to FV State
! ---------------------------------------
      do k=1,km
      vars%pe(:,:,k) = dynamics%grid%sige(k)*( vars%pe(:,:,km+1)-dynamics%grid%ptop ) + dynamics%grid%ptop
      enddo

      pkxy = vars%pe**kappa
      do k=1,km
      vars%pkz(:,:,k) = ( pkxy(:,:,k+1)-pkxy(:,:,k) ) &
                      / ( kappa*( log(vars%pe(:,:,k+1))-log(vars%pe(:,:,k)) ) )
!     vars%pkz(:,:,k) = ( pkxy(:,:,k+1)*vars%pe(:,:,k+1)-pkxy(:,:,k)*vars%pe(:,:,k) ) &
!                     / (               vars%pe(:,:,k+1)-            vars%pe(:,:,k) ) / (1+kappa)
      enddo

      vars%u(:,:,:) = dynamics%vars(1)%u(:,:,:)
      vars%v(:,:,:) = dynamics%vars(1)%v(:,:,:)

      call getpke( dynamics%vars(1)%p,pke   ,dynamics%grid,im,jm )
      call getpk ( dynamics%vars(1)%p,pke,pk,dynamics%grid,im,jm )
      tempxy(:,:,:) = dynamics%vars(1)%t(:,:,:)*pk(:,:,:)*(1+eps*dynamics%vars(1)%q(:,:,:,kq))  ! Virtual Temperature

      vars%pt(:,:,:) = tempxy/vars%pkz  ! Virtual Potential Temperature
 
      do n=1,dynamics%grid%ntracer
         if( state%vars%tracer(n)%is_r4 ) then
               state%vars%tracer(n)%content_r4(:,:,:) = dynamics%vars(1)%q(:,:,:,n)
         else
               state%vars%tracer(n)%content   (:,:,:) = dynamics%vars(1)%q(:,:,:,n)
         endif
      enddo

      deallocate (pke)
      deallocate (pk )

! Adjust Tracers for Conservation (due to Shapiro Filter)
! -------------------------------------------------------

      call glosum ( STATE,NQ,TRSUM2 )

      do n=1,NQ
      if( (trim(names(n)).ne.'CLLS') .and. &
          (trim(names(n)).ne.'CLCN')       ) then
         if( trsum2(n).ne.0.0d0 ) then
             trsum2(n) = trsum1(n)/trsum2(n)
         else
             trsum2(n) = 1.0d0
         endif
         if( STATE%VARS%TRACER(N)%IS_R4 ) then
             state%vars%tracer(n)%content_r4 = state%vars%tracer(n)%content_r4 * trsum2(n)
         else
             state%vars%tracer(n)%content    = state%vars%tracer(n)%content    * trsum2(n)
         endif
      endif
      enddo

      deallocate( trsum1 )
      deallocate( trsum2 )

! Vertically Integrated THV Tendency Diagnostic
! ---------------------------------------------
      delp  = ( vars%pe(:,:,2:) - vars%pe(:,:,:km) )
      dthdt = ( vars%pt*delp-dthdt )/dt

      call MAPL_GetPointer(export,temp2d,'DTHVDTDYNINT', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         qsum1 = 0.0
      do k=1,km
         qsum1 = qsum1 + dthdt(:,:,k)
      enddo
         temp2d = qsum1 * (MAPL_P00**MAPL_KAPPA) / grav
      end if

      call MAPL_GetPointer(export,temp2d,'DTHVDTREMAP', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = dthdtremap

      call MAPL_GetPointer(export,temp2d,'DTHVDTCONSV', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = dthdtconsv

! Unify Poles for Tracers
! -----------------------

      call PUSH_Q ( STATE )

! Load Local Variable with Vapor Specific Humidity
! ------------------------------------------------

      if ( qqq%is_r4 ) then
           qv = qqq%content_r4
      else
           qv = qqq%content
      endif

! Compute Dry Theta and T with Unified Poles
! ------------------------------------------

      vars%pt = vars%pt / (1.0+eps*qv )
      tempxy  = vars%pt * vars%pkz

! Compute Mid-Layer Pressure and Pressure Thickness
! -------------------------------------------------

      delp = ( vars%pe(:,:,2:) - vars%pe(:,:,:km) )
      pl   = ( vars%pe(:,:,2:) + vars%pe(:,:,:km) ) * 0.5

! Compute Tropopause Pressure, Temperature, and Moisture
! ------------------------------------------------------

      call tropovars ( ilastxy-ifirstxy+1,jlastxy-jfirstxy+1,km, &
                       real(vars%pe            ,kind=4),         &
                       real(pl                 ,kind=4),         &
                       real(tempxy             ,kind=4),         &
                       real(qv                 ,kind=4),         &
                       real(epvxyz*(p00**kappa),kind=4),         &
                       tropp1,tropp2,tropp3,tropt,tropq          )

      call MAPL_GetPointer(export,temp2D,'TROPP_THERMAL',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2D)) temp2D = tropp1

      call MAPL_GetPointer(export,temp2D,'TROPP_EPV',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2D)) temp2D = tropp2

      call MAPL_GetPointer(export,temp2D,'TROPP_BLENDED',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2D)) temp2D = tropp3

      call MAPL_GetPointer(export,temp2D,'TROPT',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2D)) temp2D = tropt

      call MAPL_GetPointer(export,temp2D,'TROPQ',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2D)) temp2D = tropq

! Compute A-Grid Winds
! --------------------
      call ctoa_winds ( dynamics%vars(1)%u,dynamics%vars(1)%v,ua,va,                         &
                        dynamics%grid%dlam,dynamics%grid%dphi,im,jm,km,dynamics%grid%lattice )

! Compute A-Grid Mass Fluxes
! --------------------------
      call c2a3d( grid, mfxxyz, mfyxyz, mfxxyz_a, mfyxyz_a )

! Compute Diagnostic Dynamics Tendencies
!  (Note: initial values of d(m,u,v,T,q)/dt are progs m,u,v,T,q)
! --------------------------------------------------------------

      dmdt = ( vars%pe(:,:,km+1)-vars%pe(:,:,1) - dmdt )/(grav*dt)

      dudt = (    ua-dudt )/dt
      dvdt = (    va-dvdt )/dt
      dtdt = (tempxy-dtdt )/dt
      dqdt = (    qv-dqdt )/dt

      ddpdt = ( delp - ddpdt )/dt ! Pressure Thickness Tendency

      call FILLOUT3 (export, 'DELP'      ,delp , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'DUDTDYN'   ,dudt , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'DVDTDYN'   ,dvdt , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'DTDTDYN'   ,dtdt , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'DQVDTDYN'  ,dqdt , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'DDELPDTDYN',ddpdt, rc=status); VERIFY_(STATUS)

      call FILLOUT3 (export, 'U_CGRID'   ,cxxyz, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'V_CGRID'   ,cyxyz, rc=status); VERIFY_(STATUS)

      call FILLOUT3 (export, 'PTFX' , ptfxxyz , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'PTFY' , ptfyxyz , rc=status); VERIFY_(STATUS)

      call FILLOUT3 (export, 'MFX_UR' , mfxxyz_ur  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'MFY_UR' , mfyxyz_ur  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'MFX'    , mfxxyz  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'MFY'    , mfyxyz  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'MFZ'    , mfzxyz  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'MFX_A'  , mfxxyz_a, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'MFY_A'  , mfyxyz_a, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'U'      , ua      , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'V'      , va      , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'T'      , tempxy  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'Q'      , qv      , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'PL'     , pl      , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'PLE'    , vars%pe , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'PLK'    , vars%pkz, rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'U_DGRID', vars%u  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'V_DGRID', vars%v  , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'PT'     , vars%pt , rc=status); VERIFY_(STATUS)
      call FILLOUT3 (export, 'PE'     , vars%pe , rc=status); VERIFY_(STATUS)

      call MAPL_GetPointer(export, temp3D, 'EPV', rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = epvxyz*(p00**kappa)
 
      call MAPL_GetPointer(export, temp3D, 'PV', rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = epvxyz/vars%pt
 
      call MAPL_GetPointer(export, temp3D, 'S', rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = tempxy*cp

      call MAPL_GetPointer(export, temp3d, 'TH',rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = vars%pt*(p00**kappa)

      call MAPL_GetPointer(export, temp2d, 'DMDTDYN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = dmdt


! Compute 3-D Tracer Dynamics Tendencies
! --------------------------------------

      if( associated(dqldt) ) then
          do N = 1,size(names)
             if( trim(names(N)).eq.'QLCN' .or. &
                 trim(names(N)).eq.'QLLS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then 
                     dqldt = dqldt + state%vars%tracer(N)%content_r4
                 else
                     dqldt = dqldt + state%vars%tracer(N)%content
                 endif
             endif
          enddo
          dqldt = dqldt/dt
      endif

      if( associated(dqidt) ) then
          do N = 1,size(names)
             if( trim(names(N)).eq.'QICN' .or. &
                 trim(names(N)).eq.'QILS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then 
                     dqidt = dqidt + state%vars%tracer(N)%content_r4
                 else
                     dqidt = dqidt + state%vars%tracer(N)%content
                 endif
             endif
          enddo
          dqidt = dqidt/dt
      endif

      if( associated(doxdt) ) then
          do N = 1,size(names)
             pos = index(names(N),'::')
             if(pos > 0) then
               if( (names(N)(pos+2:))=='OX' ) then
                 if( state%vars%tracer(N)%is_r4 ) then 
                     doxdt = doxdt + state%vars%tracer(N)%content_r4
                 else
                     doxdt = doxdt + state%vars%tracer(N)%content
                 endif
               endif
             endif
          enddo
          doxdt = doxdt/dt
      endif

! Compute 2-D Vertically Integrated Tracer Dynamics Tendencies
! ------------------------------------------------------------

      call MAPL_GetPointer ( export, temp2D, 'DQVDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          tempr8(:,:,2) = 0.0
          do k=1,km
          tempr8(:,:,2) = tempr8(:,:,2) + qv(:,:,k)*delp(:,:,k)
          enddo
          tempr8(:,:,2) = ( tempr8(:,:,2)-tempr8(:,:,1) )/(grav*dt)
          temp2d = tempr8(:,:,2)
      endif

      call MAPL_GetPointer ( export, temp2D, 'DQLDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          do N = 1,size(names)
             if( trim(names(N)).eq.'QLCN' .or. &
                 trim(names(N)).eq.'QLLS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then 
                     do k=1,km
                     temp2d = temp2d + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                 else
                     do k=1,km
                     temp2d = temp2d + state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                 endif
             endif
          enddo
          temp2d = temp2d/(grav*dt)
      endif

      call MAPL_GetPointer ( export, temp2D, 'DQIDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          do N = 1,size(names)
             if( trim(names(N)).eq.'QICN' .or. &
                 trim(names(N)).eq.'QILS' ) then
                 if( state%vars%tracer(N)%is_r4 ) then 
                     do k=1,km
                     temp2d = temp2d + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                 else
                     do k=1,km
                     temp2d = temp2d + state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                 endif
             endif
          enddo
          temp2d = temp2d/(grav*dt)
      endif

      call MAPL_GetPointer ( export, temp2D, 'DOXDTDYNINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          do N = 1,size(names)
             pos = index(names(N),'::')
             if(pos > 0) then
               if( (names(N)(pos+2:))=='OX' ) then
                 if( state%vars%tracer(N)%is_r4 ) then 
                     do k=1,km
                     temp2d = temp2d + state%vars%tracer(N)%content_r4(:,:,k)*delp(:,:,k)
                     enddo
                 else
                     do k=1,km
                     temp2d = temp2d + state%vars%tracer(N)%content(:,:,k)*delp(:,:,k)
                     enddo
                 endif
               endif
             endif
          enddo
          temp2d = temp2d * (MAPL_O3MW/MAPL_AIRMW) / (MAPL_GRAV*DT)
      endif

! Fill Surface and Near-Surface Variables
! ---------------------------------------

      call MAPL_GetPointer(export,temp2d,'PS',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d =  vars%pe(:,:,km+1)    

      call MAPL_GetPointer(export,temp2d,'US',  rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d =       ua(:,:,km)      

      call MAPL_GetPointer(export,temp2d,'VS'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d =       va(:,:,km)      

      call MAPL_GetPointer(export,temp2d,'TA'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d =   tempxy(:,:,km)      

      call MAPL_GetPointer(export,temp2d,'QA'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d =       qv(:,:,km)      

      call MAPL_GetPointer(export,temp2d,'SPEED',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = sqrt( ua(:,:,km)**2 + va(:,:,km)**2 ) 


! Virtual temperature
! -------------------

      tempxy =  tempxy*(1.0+eps*qv)

      call MAPL_GetPointer(export,temp3D,'TV'   ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp3D)) temp3D = tempxy


! Fluxes: UCPT & VCPT
!--------------------
      call MAPL_GetPointer(export,temp2d,'UCPT',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
         do k=1,km
            temp2d = temp2d + ua(:,:,k)*tempxy(:,:,k)*delp(:,:,k)
         enddo
         temp2d = temp2d*(cp/grav)
      end if

      call MAPL_GetPointer(export,temp2d,'VCPT',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
         do k=1,km
            temp2d = temp2d + va(:,:,k)*tempxy(:,:,k)*delp(:,:,k)
         enddo
         temp2d = temp2d*(cp/grav)
      end if


! Compute Energetics After Dycore
! -------------------------------

      tempxy = vars%pt*(1.0+eps*qv)  ! Convert TH to THV

      call MAPL_GetPointer(export,temp3d,'THV',rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = tempxy

      call Energetics (state,vars%u,vars%v,tempxy,vars%pe,delp,vars%pkz,phisxy,kenrg,penrg,tenrg,dynamics%grid)

      kedyn   = (kenrg -kenrg0)/DT
      pedyn   = (penrg -penrg0)/DT
      tedyn   = (tenrg -tenrg0)/DT

      kecdcor = (kenrga-kenrg0)/DT
      pecdcor = (penrga-penrg0)/DT
      tecdcor = (tenrga-tenrg0)/DT

      keremap = (kenrgb-kenrga)/DT
      peremap = (penrgb-penrga)/DT
      teremap = (tenrgb-tenrga)/DT

      call MAPL_GetPointer(export,temp2d,'KEDYN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = kedyn

      call MAPL_GetPointer(export,temp2d,'PEDYN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = pedyn

      call MAPL_GetPointer(export,temp2d,'TEDYN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = tedyn

      call MAPL_GetPointer(export,temp2d,'KECDCOR',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = kecdcor

      call MAPL_GetPointer(export,temp2d,'PECDCOR',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = pecdcor

      call MAPL_GetPointer(export,temp2d,'TECDCOR',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = tecdcor

      call MAPL_GetPointer(export,temp2d,'KEREMAP',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = keremap

      call MAPL_GetPointer(export,temp2d,'PEREMAP',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = peremap

      call MAPL_GetPointer(export,temp2d,'TEREMAP',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = teremap


! Fluxes: UKE & VKE
! -----------------
      call MAPL_GetPointer(export,tempu,'UKE',rc=status); VERIFY_(STATUS)
      call MAPL_GetPointer(export,tempv,'VKE',rc=status); VERIFY_(STATUS)

      if(associated(tempu) .or. associated(tempv)) then
         ke = 0.5*(ua**2 + va**2)
      end if

      if(associated(tempu)) then
         tempu = 0.0
         do k=1,km
            tempu = tempu + ua(:,:,k)*ke(:,:,k)*delp(:,:,k)
         enddo
         tempu = tempu / grav
      end if

      if(associated(tempv)) then
         tempv = 0.0
         do k=1,km
            tempv = tempv + va(:,:,k)*ke(:,:,k)*delp(:,:,k)
         enddo
         tempv = tempv / grav
      end if

! Fluxes: UQV & VQV
! -----------------
      call MAPL_GetPointer(export,temp2d,'UQV',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
         do k=1,km
            temp2d = temp2d + ua(:,:,k)*QV(:,:,k)*delp(:,:,k)
         enddo
         temp2d = temp2d / grav
      end if

      call MAPL_GetPointer(export,temp2d,'VQV',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
         do k=1,km
            temp2d = temp2d + va(:,:,k)*QV(:,:,k)*delp(:,:,k)
         enddo
         temp2d = temp2d / grav
      end if

! Fluxes: UQL & VQL
! -----------------
      call MAPL_GetPointer(export,temp2d,'UQL',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
          do N = 1,size(names)
             if( trim(names(n)).eq.'QLCN' .or. &
                 trim(names(n)).eq.'QLLS' ) then
                 do k=1,km
                 if( state%vars%tracer(n)%is_r4 ) then 
                      temp2d = temp2d + ua(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                 else
                      temp2d = temp2d + ua(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo
             endif
          enddo
         temp2d = temp2d / grav
      end if

      call MAPL_GetPointer(export,temp2d,'VQL',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
          do N = 1,size(names)
             if( trim(names(n)).eq.'QLCN' .or. &
                 trim(names(n)).eq.'QLLS' ) then
                 do k=1,km
                 if( state%vars%tracer(n)%is_r4 ) then 
                      temp2d = temp2d + va(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                 else
                      temp2d = temp2d + va(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo
             endif
          enddo
         temp2d = temp2d / grav
      end if

! Fluxes: UQI & VQI
! -----------------
      call MAPL_GetPointer(export,temp2d,'UQI',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
          do N = 1,size(names)
             if( trim(names(n)).eq.'QICN' .or. &
                 trim(names(n)).eq.'QILS' ) then
                 do k=1,km
                 if( state%vars%tracer(n)%is_r4 ) then 
                      temp2d = temp2d + ua(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                 else
                      temp2d = temp2d + ua(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo
             endif
          enddo
         temp2d = temp2d / grav
      end if

      call MAPL_GetPointer(export,temp2d,'VQI',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         temp2d = 0.0
          do N = 1,size(names)
             if( trim(names(n)).eq.'QICN' .or. &
                 trim(names(n)).eq.'QILS' ) then
                 do k=1,km
                 if( state%vars%tracer(n)%is_r4 ) then 
                      temp2d = temp2d + va(:,:,k)*state%vars%tracer(n)%content_r4(:,:,k)*delp(:,:,k)
                 else
                      temp2d = temp2d + va(:,:,k)*state%vars%tracer(n)%content   (:,:,k)*delp(:,:,k)
                 endif
                 enddo
             endif
          enddo
         temp2d = temp2d / grav
      end if

! Height related diagnostics
! --------------------------
      zle(:,:,km+1) = phisxy(:,:)
      do k=km,1,-1
        zle(:,:,k) = zle(:,:,k+1) + cp*tempxy(:,:,k)*( pkxy(:,:,k+1)-pkxy(:,:,k) )
      enddo
      zle = zle/grav

      call MAPL_GetPointer(export,temp3d,'ZLE',rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = zle

      call MAPL_GetPointer(export,temp2d,'DZ', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = 0.5*( zle(:,:,km)-zle(:,:,km+1) )

      call MAPL_GetPointer(export,temp3d,'ZL' ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = 0.5*( zle(:,:,:km)+zle(:,:,2:) )

      call MAPL_GetPointer(export,temp3d,'S'  ,rc=status)
      VERIFY_(STATUS)
      if(associated(temp3d)) temp3d = temp3d + grav*(0.5*( zle(:,:,:km)+zle(:,:,2:) ))

! Energy Budget Calculations
! --------------------------
      convke  = 0.0
      convthv = 0.0
      convcpt = 0.0
      convphi = 0.0
      kegen   = 0.0
      do k=1,km
      kegen   = kegen   + omaxyz(:,:,k)*grav*( zle(:,:,k+1)-zle(:,:,k) )
      convke  = convke  +     ke(:,:,k)*delp(:,:,k)
      convphi = convphi +    phi(:,:,k)*delp(:,:,k)
      convthv = convthv + thvxyz(:,:,k)
      convcpt = convcpt + cptxyz(:,:,k)
      enddo
      kegen   =  kegen  /grav
      convthv =  convthv/grav * (MAPL_P00**MAPL_KAPPA)
      convcpt =  convcpt/grav
      convke  = (convke /grav-kenrg0)/dt
      convphi = (convphi/grav- phi00)/dt - convcpt

      call MAPL_GetPointer(export,temp2d,'KEHOT',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = kehot

      call MAPL_GetPointer(export,temp2d,'KEDP',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = kedp

      call MAPL_GetPointer(export,temp2d,'KEADV',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = keadv

      call MAPL_GetPointer(export,temp2d,'KEPG',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = kepg

      call MAPL_GetPointer(export,temp2d,'KEGEN', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = kegen

      call MAPL_GetPointer(export,temp2d,'CONVKE',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = convke

      call MAPL_GetPointer(export,temp2d,'CONVTHV', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = convthv

      call MAPL_GetPointer(export,temp2d,'CONVCPT', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = convcpt

      call MAPL_GetPointer(export,temp2d,'CONVPHI',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = convphi

      call MAPL_GetPointer(export,temp2d,'DKERESIN',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = keadv + kedp + kehot - convke

      call MAPL_GetPointer(export,temp2d,'DKERESPG',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = kepg - ( convphi + kegen - tedyn )

      call MAPL_GetPointer(export,temp2d,'QFIXER',rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) temp2d = pedyn - convcpt + kegen

! Fluxes: UPHI & VPHI
! -------------------
      call MAPL_GetPointer(export,tempu,'UPHI',rc=status); VERIFY_(STATUS)
      call MAPL_GetPointer(export,tempv,'VPHI',rc=status); VERIFY_(STATUS)

      if( associated(tempu).or.associated(tempv) ) zl = 0.5*( zle(:,:,:km)+zle(:,:,2:) )

      if(associated(tempu)) then
         tempu = 0.0
         do k=1,km
            tempu = tempu + ua(:,:,k)*zl(:,:,k)*delp(:,:,k)
         enddo
      end if

      if(associated(tempv)) then
         tempv = 0.0
         do k=1,km
            tempv = tempv + va(:,:,k)*zl(:,:,k)*delp(:,:,k)
         enddo
      end if


! Compute Omega
! -------------
      call FILLOUT3 (export,'OMEGA',omaxyz,rc=status)
      VERIFY_(STATUS)

      call MAPL_GetPointer(export,temp2d,'OMEGA500', rc=status)
      VERIFY_(STATUS)
      if(associated(temp2d)) then
         call VertInterp(temp2d,omaxyz,log(pl),log(50000.),log(vars%pe(:,:,km+1)),status)
         VERIFY_(STATUS)
      end if
         
! De-Allocate Arrays
! ------------------

      DEALLOCATE( tempr8       )
      DEALLOCATE( dummy        )
      DEALLOCATE( dqvdtanaint1 )
      DEALLOCATE( dqvdtanaint2 )
      DEALLOCATE( dqldtanaint1 )
      DEALLOCATE( dqldtanaint2 )
      DEALLOCATE( dqidtanaint1 )
      DEALLOCATE( dqidtanaint2 )
      DEALLOCATE( doxdtanaint1 )
      DEALLOCATE( doxdtanaint2 )
      DEALLOCATE( dthdtanaint1 )
      DEALLOCATE( dthdtanaint2 )
      DEALLOCATE( dthdtremap   )
      DEALLOCATE( dthdtconsv   )

      DEALLOCATE( TROPP1 )
      DEALLOCATE( TROPP2 )
      DEALLOCATE( TROPP3 )
      DEALLOCATE( TROPT  )
      DEALLOCATE( TROPQ  )

      DEALLOCATE( KEPG   )
      DEALLOCATE( KEADV  )
      DEALLOCATE( KEDP   )
      DEALLOCATE( KEHOT  )
      DEALLOCATE( KENRG  )
      DEALLOCATE( PENRG  )
      DEALLOCATE( TENRG  )
      DEALLOCATE( KENRG0 )
      DEALLOCATE( PENRG0 )
      DEALLOCATE( TENRG0 )

      DEALLOCATE( KEGEN  )
      DEALLOCATE( KEDYN  )
      DEALLOCATE( PEDYN  )
      DEALLOCATE( TEDYN  )
      DEALLOCATE( KECDCOR)
      DEALLOCATE( PECDCOR)
      DEALLOCATE( TECDCOR)
      DEALLOCATE( KEREMAP)
      DEALLOCATE( PEREMAP)
      DEALLOCATE( TEREMAP)

      DEALLOCATE( CONVKE )
      DEALLOCATE( CONVCPT)
      DEALLOCATE( CONVPHI)
      DEALLOCATE( CONVTHV)

      DEALLOCATE( KENRGA )
      DEALLOCATE( PENRGA )
      DEALLOCATE( TENRGA )
      DEALLOCATE( KENRGB )
      DEALLOCATE( PENRGB )
      DEALLOCATE( TENRGB )

      DEALLOCATE( ke    )
      DEALLOCATE( cpt   )
      DEALLOCATE( phi   )
      DEALLOCATE( gze   )
      DEALLOCATE( qsum1 )
      DEALLOCATE( qsum2 )

      DEALLOCATE( ZL     )
      DEALLOCATE( ZLE    )
      DEALLOCATE( PKXY   )
      DEALLOCATE( pelnxz )
      DEALLOCATE( omaxyz )
      DEALLOCATE( cptxyz )
      DEALLOCATE( thvxyz )
      DEALLOCATE( epvxyz )
      DEALLOCATE(  cxxyz )
      DEALLOCATE(  cyxyz )
      DEALLOCATE( ptfxxyz )
      DEALLOCATE( ptfyxyz )
      DEALLOCATE( mfxxyz_ur )
      DEALLOCATE( mfyxyz_ur )
      DEALLOCATE( mfxxyz )
      DEALLOCATE( mfyxyz )
      DEALLOCATE( mfzxyz )
      DEALLOCATE( mfxxyz_a )
      DEALLOCATE( mfyxyz_a )
      DEALLOCATE( tempxy )
      DEALLOCATE( pl     )
      DEALLOCATE( va     )
      DEALLOCATE( ua     )
      DEALLOCATE( qv     )
      DEALLOCATE( ql     )
      DEALLOCATE( qi     )
      DEALLOCATE( qdnew  )
      DEALLOCATE( qdold  )
      DEALLOCATE( qvold  )
      DEALLOCATE( qlold  )
      DEALLOCATE( qiold  )
      DEALLOCATE( ox     )
      DEALLOCATE( delp   )
      DEALLOCATE( dmdt   )
      DEALLOCATE( dudt   )
      DEALLOCATE( dvdt   )
      DEALLOCATE( dtdt   )
      DEALLOCATE( dqdt   )
      DEALLOCATE( dthdt  )
      DEALLOCATE( ddpdt  )
      DEALLOCATE( phisxy )
      DEALLOCATE( names  )
      DEALLOCATE( phi00  )
    
      do i=1,size(STATE%VARS%tracer)-NXQ
         if ( .not. STATE%VARS%TRACER(I)%IS_R4 ) then
            DEALLOCATE(STATE%VARS%tracer(i)%content, STAT=STATUS)   ! TEMPORARY, till pointers are passed
         end if
      enddo

      DEALLOCATE( STATE%VARS%tracer, STAT=STATUS )   ! Comment out to output tracer to checkpoint file

  call MAPL_TimerOff(MAPL,"RUN1")
  call MAPL_TimerOff(MAPL,"TOTAL")

  RETURN_(ESMF_SUCCESS)

end subroutine RUN1

!-----------------------------------------------------------------------

  subroutine PULL_Q(STATE, IMPORT, QQQ, NXQ, RC)
   
    type (T_FVDYCORE_STATE)        :: STATE
    type (ESMF_State)              :: IMPORT
    type (T_TRACERS)               :: QQQ       ! Specific Humidity
    integer,           intent(IN)  :: NXQ
    integer, optional, intent(OUT) :: RC

    integer                          :: STATUS
    character(len=ESMF_MAXSTR)       :: IAm="Pull_Q"
    character(len=ESMF_MAXSTR)       :: FIELDNAME
    type (ESMF_FieldBundle)          :: BUNDLE
    type (ESMF_Field)                :: field
    type (ESMF_Array)                :: array
    type (ESMF_TypeKind_Flag)        :: kind
    real(r4),              pointer   :: ptr_r4(:,:,:), humidity(:,:,:)
    real(r8),              pointer   :: ptr(:,:,:)
    integer                          :: I,K,N,NQ
    logical                          :: EMPTY
    integer                          :: i1,in,j1,jn,im,jm,km
    real(r8)                         :: sumout

    i1 = state%grid%ifirstxy
    in = state%grid%ilastxy
    j1 = state%grid%jfirstxy
    jn = state%grid%jlastxy
    im = state%grid%im
    jm = state%grid%jm
    km = state%grid%km

    call ESMF_StateGet(IMPORT, 'TRADV' , BUNDLE,   RC=STATUS)
    VERIFY_(STATUS)

! Count the friendlies
!---------------------

    call ESMF_FieldBundleGet(BUNDLE, fieldCount=NQ, RC=STATUS)
    VERIFY_(STATUS)

               NQ = NQ + NXQ
    STATE%GRID%NQ = NQ       ! GRID%NQ is now the "official" NQ

!
! Tracer pointer array
!
    ALLOCATE(STATE%VARS%tracer(nq), STAT=STATUS)
    VERIFY_(STATUS)

    DO n = 1, NQ-NXQ
       call ESMF_FieldBundleGet(bundle, fieldIndex=n, field=field, rc=status)
       VERIFY_(STATUS)
       call ESMF_FieldGet(field, name=fieldname, rc=status)
       VERIFY_(STATUS)
       call ESMF_FieldGet(field, Array=array, rc=status)
       VERIFY_(STATUS)

       call ESMF_ArrayGet(array,typekind=kind,rc=status)
       VERIFY_(STATUS)

       STATE%VARS%TRACER(N)%IS_R4  = (kind == ESMF_TYPEKIND_R4)   ! Is real*4?

       if ( STATE%VARS%TRACER(N)%IS_R4 ) then
          call ESMF_ArrayGet(array, localDE=0, farrayptr=ptr_r4, rc=status)
          VERIFY_(STATUS)
          state%vars%tracer(n)%content_r4 => MAPL_RemapBounds(PTR_R4, i1,in,j1,jn, &
                                                              1, STATE%GRID%KM)
          if (fieldname == "Q") then 
             qqq%is_r4 = .true.
             qqq%content_r4 => state%vars%tracer(n)%content_r4
          end if

! Constrain Poles 
! ----------------
          if ( j1 == 1 ) then
             do k=1,km
                call par_xsum_r4 ( state%grid, state%vars%tracer(n)%content_r4(i1:in,1,k), &
                                   1, sumout )
                sumout = sumout/im
                do i=i1,in
                   state%vars%tracer(n)%content_r4(i,1,k) = sumout
                enddo
              enddo
           endif
           if ( jn == jm ) then
              do k=1,km
                 call par_xsum_r4 ( state%grid, state%vars%tracer(n)%content_r4(i1:in,jm,k), &
                                    1, sumout )
                 sumout = sumout/im
                 do i=i1,in
                   state%vars%tracer(n)%content_r4(i,jm,k) = sumout
                 enddo
               enddo
           endif

        else   ! Tracer is R8

            call ESMF_ArrayGet(array, localDE=0, farrayptr=ptr, rc=status)
            VERIFY_(STATUS)

            state%vars%tracer(n)%content => PTR
            if (fieldname == "Q") then 
               qqq%is_r4   = .false.
               qqq%content => state%vars%tracer(n)%content
            end if

! Constrain Poles 
! ----------------
            if ( j1 == 1 ) then
               do k=1,km
                   call par_xsum ( state%grid, state%vars%tracer(n)%content(i1:in,1,k), &
                                   1, sumout )
                   sumout = sumout/im
                   do i=i1,in
                      state%vars%tracer(n)%content(i,1,k) = sumout
                   enddo
                enddo
             endif
             if ( jn == jm ) then
                do k=1,km
                   call par_xsum ( state%grid, state%vars%tracer(n)%content(i1:in,jm,k), &
                                   1, sumout )
                   sumout = sumout/im
                   do i=i1,in
                      state%vars%tracer(n)%content(i,jm,k) = sumout
                   enddo
                enddo
             endif
          endif
       END DO

  end subroutine PULL_Q

!-----------------------------------------------------------------------

  subroutine PUSH_Q (STATE)
   
    type (T_FVDYCORE_STATE)        :: STATE


    integer                          :: STATUS
    integer                          :: I,K,N
    integer                          :: i1,in,j1,jn,im,jm,km
    real(r8)                         :: sumout

    i1 = state%grid%ifirstxy
    in = state%grid%ilastxy
    j1 = state%grid%jfirstxy
    jn = state%grid%jlastxy
    im = state%grid%im
    jm = state%grid%jm
    km = state%grid%km

! Count the friendlies
!---------------------

    DO N = 1, state%grid%NQ

! Constrain Poles
! ---------------
    if ( state%vars%tracer(n)%is_r4 ) then
      if ( j1 == 1 ) then
         do k=1,km
            call par_xsum_r4 ( state%grid, state%vars%tracer(n)%content_r4(i1:in,1,k), &
                               1, sumout )
            sumout = sumout/im
            do i=i1,in
              state%vars%tracer(n)%content_r4(i,1,k) = sumout
            enddo
         enddo
      endif
      if ( jn == jm ) then
         do k=1,km
            call par_xsum_r4 ( state%grid, state%vars%tracer(n)%content_r4(i1:in,jm,k),   &
                               1, sumout )
            sumout = sumout/im
            do i=i1,in
              state%vars%tracer(n)%content_r4(i,jm,k) = sumout
            enddo
         enddo
      endif

    else                        ! Content is R8

      if ( j1 == 1 ) then
         do k=1,km
            call par_xsum ( state%grid, state%vars%tracer(n)%content(i1:in,1,k), 1, sumout )
            sumout = sumout/im
            do i=i1,in
              state%vars%tracer(n)%content(i,1,k) = sumout
            enddo
         enddo
      endif
      if ( jn == jm ) then
         do k=1,km
            call par_xsum ( state%grid, state%vars%tracer(n)%content(i1:in,jm,k), 1, sumout )
            sumout = sumout/im
            do i=i1,in
              state%vars%tracer(n)%content(i,jm,k) = sumout
            enddo
         enddo
      endif

    endif
    END DO

  end subroutine PUSH_Q
!-----------------------------------------------------------------------

   subroutine FILLOUT3(export, name, V, RC)
     type (ESMF_State),  intent(inout) :: export
     character(len=*),   intent(IN   ) :: name
     real(r8),           intent(IN   ) :: V(:,:,:)
     integer, optional,  intent(  out) :: rc

     real(r4), pointer          :: CPL(:,:,:)
     integer                    :: status
     character(len=ESMF_MAXSTR) :: IAm="Fillout3"

     call MAPL_GetPointer(export, cpl, name, RC=STATUS)
     VERIFY_(STATUS)
     if(associated(cpl)) cpl=v

   end subroutine FILLOUT3

!-----------------------------------------------------------------------

   subroutine FILLOUT2(export, name, V, rc)
     type (ESMF_State),  intent(inout) :: export
     character(len=*),   intent(IN   ) :: name
     real(r8),           intent(IN   ) :: V(:,:)
     integer, optional,  intent(  out) :: rc

     real(kind=4), pointer      :: CPL(:,:)
     integer                    :: status
     character(len=ESMF_MAXSTR) :: IAm="Fillout2"

     call MAPL_GetPointer(export, cpl, name, RC=STATUS)
     VERIFY_(STATUS)
     if(associated(cpl)) cpl=v

     return
   end subroutine FILLOUT2

!-----------------------------------------------------------------------

  subroutine Energetics (state,ud,vd,thv,ple,delp,pk,phiS,keint,peint,teint,grid,ke,cpt,gze)
  use g3_dynamics_state_module
  implicit none
  type ( dynamics_grid_type ) grid

  type (T_FVDYCORE_STATE) :: STATE
  real(8), optional, intent(out) ::   ke(:,:,:)
  real(8), optional, intent(out) ::  cpt(:,:,:)
  real(8), optional, intent(out) ::  gze(:,:,:)
  real(8)   ud(:,:,:)
  real(8)   vd(:,:,:)
  real(8)  thv(:,:,:)
  real(8)  ple(:,:,:)
  real(8) delp(:,:,:)
  real(8)   pk(:,:,:)
  real(8)   keint(:,:)
  real(8)   peint(:,:)
  real(8)   teint(:,:)
  real(8) phiS(:,:)

  real(8) kinetic, potential, sump
  integer i,ifirst,ilast
  integer j,jfirst,jlast
  integer im,jm,km,k

  real(8), allocatable ::   ud2(:,:,:)
  real(8), allocatable ::   vd2(:,:,:)
  real(8), allocatable ::   ua2(:,:,:)
  real(8), allocatable ::   va2(:,:,:)
  real(8), allocatable ::   pke(:,:,:)
  real(8), allocatable ::  phiT(:,:)

  ifirst = lbound( ud,1 )
  ilast  = ubound( ud,1 )
  jfirst = lbound( ud,2 )
  jlast  = ubound( ud,2 )
  km     = ubound( ud,3 )
  im     = state%grid%im 
  jm     = state%grid%jm 

  allocate( ua2  ( ifirst:ilast, jfirst:jlast , 1:km   ) )
  allocate( va2  ( ifirst:ilast, jfirst:jlast , 1:km   ) )
  allocate( ud2  ( ifirst:ilast, jfirst:jlast , 1:km   ) )
  allocate( vd2  ( ifirst:ilast, jfirst:jlast , 1:km   ) )
  allocate( pke  ( ifirst:ilast, jfirst:jlast , 1:km+1 ) )
  allocate( phiT ( ifirst:ilast, jfirst:jlast ) )

! Compute Model Edge Heights
! --------------------------
    pke  = ple**kappa
    phiT = phiS
    if( present(gze) ) gze(:,:,km+1) = phiS
    do k=km,1,-1
                       phiT = phiT + cp*thv(:,:,k)*( pke(:,:,k+1)-pke(:,:,k) )
    if( present(gze) ) gze(:,:,k) = phiT
    enddo

! Compute C-Grid Kinetic Energy (using FV D-Grid Variable Names)
! --------------------------------------------------------------
    ud2 = ud*ud
    vd2 = vd*vd
    call ctoa_winds ( ud2,vd2,ua2,va2,                                    &
                      grid%dlam,grid%dphi,grid%im,grid%jm,km,grid%lattice )

    if( state%grid%jfirstxy.eq.1 ) then
        ua2(:,1,:) = ud2(:,2,:)
        va2(:,1,:) = vd2(:,2,:)
    endif
    if( state%grid%jlastxy.eq.jm ) then
        ua2(:,jlast,:) = ud2(:,jlast-1,:)
        va2(:,jlast,:) = vd2(:,jlast-1,:)
    endif

! Compute Energetics:  Cp*Tv + K + PHI
! ------------------------------------
       keint = 0.0
       peint = 0.0
  do k=1,km
  do j=jfirst,jlast
  do i=ifirst,ilast
       kinetic      = 0.5_r8*( ua2(i,j,k) + va2(i,j,k) )
       potential    =  cp*thv(i,j,k)*pk(i,j,k)
       keint(i,j)   =   keint(i,j) +   kinetic  *delp(i,j,k)
       peint(i,j)   =   peint(i,j) +   potential*delp(i,j,k)
       if( present(ke)  )  ke(i,j,k) = kinetic
       if( present(cpt) ) cpt(i,j,k) = potential
  enddo
  enddo
  enddo
       keint(:,:) =    keint(:,:)/grav
       peint(:,:) =    peint(:,:)/grav
       teint(:,:) = (phiS(:,:)*ple(:,:,km+1)-phiT(:,:)*ple(:,:,1))/grav

  if( state%grid%jfirstxy.eq.1 ) then
      call par_xsum ( state%grid, keint(ifirst:ilast,1), 1, sump )        ! Unify Pole Estimate
!     call par_xsum ( state%grid, keint(ifirst:ilast,2), 1, sump )        ! Average Surrounding Points to Pole Location
      sump = sump/im
      do i=ifirst,ilast
         keint(i,1) = sump
      enddo
  endif
  if( state%grid%jlastxy.eq.jm ) then
      call par_xsum ( state%grid, keint(ifirst:ilast,jlast  ), 1, sump )  ! Unify Pole Estimate
!     call par_xsum ( state%grid, keint(ifirst:ilast,jlast-1), 1, sump )  ! Average Surrounding Points to Pole Location
      sump = sump/im
      do i=ifirst,ilast
         keint(i,jlast) = sump
      enddo
  endif

  deallocate ( ua2  )
  deallocate ( va2  )
  deallocate ( ud2  )
  deallocate ( vd2  )
  deallocate ( pke  )
  deallocate ( phiT )

  return
  end subroutine Energetics

!-----------------------------------------------------------------------

  subroutine Run2(gc, import, export, clock, rc)
    use mod_comm, only: commglobal, mp_swapirr

    use g3_dynamics_state_module
    type ( dynamics_grid_type ) g3_grid
    save                        g3_grid

    include 'mpif.h'
    integer imglobal, jmglobal, nx, ny
    character*2 cnx,cny
    logical first
    data    first /.true./

    type(ESMF_GridComp), intent(inout) :: gc
    type (ESMF_State),   intent(inout) :: import
    type (ESMF_State),   intent(inout) :: export
    type (ESMF_Clock),   intent(in)    :: clock
    integer, intent(out), optional     :: rc 

! !Local Variables:
  
    integer                                          :: status
    character(len=ESMF_MAXSTR) :: IAm

    type (MAPL_MetaComp), pointer :: MAPL 

    type (DYN_wrap) :: wrap
    type (T_FVDYCORE_STATE), pointer :: STATE
    type (T_FVDYCORE_GRID),  pointer :: GRID
    type (T_FVDYCORE_VARS),  pointer :: VARS
    type (T_TRACERS)                 :: qqq     ! Specific Humidity
    
    real(r8), allocatable :: penrg (:,:)   ! Vertically Integrated Cp*T
    real(r8), allocatable :: kenrg (:,:)   ! Vertically Integrated K
    real(r8), allocatable :: tenrg (:,:)   ! PHIS*(Psurf-Ptop)
    real(r8), allocatable :: penrg0(:,:)   ! Vertically Integrated Cp*T
    real(r8), allocatable :: kenrg0(:,:)   ! Vertically Integrated K
    real(r8), allocatable :: tenrg0(:,:)   ! PHIS*(Psurf-Ptop)

    real(r8),     pointer :: phisxy(:,:)
    real(r4),     pointer ::   phis(:,:)
    real(r8), allocatable ::    slp(:,:)
    real(r8), allocatable ::  H1000(:,:)
    real(r8), allocatable ::  H850 (:,:)
    real(r8), allocatable ::  H500 (:,:)
    real(r8), allocatable ::    pke(:,:,:)
    real(r8), allocatable ::     pl(:,:,:)
    real(r8), allocatable ::     ua(:,:,:)
    real(r8), allocatable ::     va(:,:,:)
    real(r8), allocatable ::  va_yz(:,:,:)
    real(r8), allocatable ::  vd_yz(:,:,:)
    real(r8), allocatable ::     qv(:,:,:)
    real(r8), allocatable ::     dp(:,:,:)
    real(r8), allocatable ::    thv(:,:,:)
    real(r8), allocatable ::    zle(:,:,:)
    real(r8), allocatable :: tempxy(:,:,:)

    real(r8), allocatable ::  logpl(:,:,:)
    real(r8), allocatable ::  logpe(:,:,:)
    real(r8), allocatable ::  logps(:,:)

    real(r8)              :: dt
    real(r8)              :: delp          ! delta pressure thickness
    real(r8)              :: kinetic       ! local kinetic   energy
    real(r8)              :: potential     ! local potential energy


    real(r4), pointer     :: QOLD(:,:,:)
    real(r4), pointer     :: temp3d(:,:,:)
    real(r4), pointer     :: temp2d(:,:  )
    real(r4), pointer     :: ztemp1(:,:  )
    real(r4), pointer     :: ztemp2(:,:  )
    real(r4), pointer     :: ztemp3(:,:  )

    real(kind=4), allocatable :: dthdtphyint1(:,:)
    real(kind=4), allocatable :: dthdtphyint2(:,:)

    integer ifirstxy, ilastxy, ifirst, ilast
    integer jfirstxy, jlastxy, jfirst, jlast
    integer kfirst,   klast
    integer im,jm,km, nxq, im1, ng_s, ng_c, ng_d
    integer i,j,k

    character(len=ESMF_MAXSTR) :: COMP_NAME

    Iam = "Run2"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the generic state
! -----------------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RUN2")

! Retrieve the pointer to the internal state
! ------------------------------------------

    call ESMF_UserCompGetInternalState(gc, 'FVstate', wrap, status)
    VERIFY_(STATUS)
    state => wrap%dyn_state

    vars  => state%vars   ! direct handle to control variables
    grid  => state%grid   ! direct handle to grid
    dt    =  state%dt     ! dynamics time step (large)

    ifirstxy = grid%ifirstxy
    ilastxy  = grid%ilastxy
    jfirstxy = grid%jfirstxy
    jlastxy  = grid%jlastxy
    jfirst   = GRID%jfirst
    jlast    = GRID%jlast
    kfirst   = GRID%kfirst
    klast    = GRID%klast
    ng_s     = GRID%ng_s
    ng_c     = GRID%ng_c
    ng_d     = GRID%ng_d

    imglobal = grid%im
    jmglobal = grid%jm
    km       = grid%km
    nxq      = 0

! **********************************************************************
! ****                      Create G3 Grid                          ****
! **********************************************************************

    if( first ) then
        call MAPL_GetResource( MAPL, NX, 'NX:', default=0, RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GetResource( MAPL, NY, 'NY:', default=0, RC=STATUS )
        VERIFY_(STATUS)

        call create_dynamics_lattice ( g3_grid%lattice,nx,ny )
        call   init_dynamics_lattice ( g3_grid%lattice,mpi_comm_world,imglobal,jmglobal,km )
        call create_dynamics_grid    ( g3_grid,imglobal,jmglobal,km )
        call   init_dynamics_grid    ( g3_grid,imglobal,jmglobal,km,0,state%grid%ak,state%grid%bk )

        first = .false.
    endif

    im = g3_grid%im
    jm = g3_grid%jm

    ALLOCATE( dthdtphyint1(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE( dthdtphyint2(ifirstxy:ilastxy,jfirstxy:jlastxy) )

    ALLOCATE(  kenrg(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE(  penrg(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE(  tenrg(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE( kenrg0(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE( penrg0(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE( tenrg0(ifirstxy:ilastxy,jfirstxy:jlastxy) )

    ALLOCATE( phisxy(ifirstxy:ilastxy,jfirstxy:jlastxy) )
    ALLOCATE(  logps(ifirstxy:ilastxy,jfirstxy:jlastxy) )

    ALLOCATE(     ua(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE(     va(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE(     qv(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE(     pl(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE(  logpl(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE(     dp(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE(    thv(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )
    ALLOCATE( tempxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km)   )

    ALLOCATE(    pke(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
    ALLOCATE(  logpe(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
    ALLOCATE(    zle(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )

    call MAPL_GetPointer ( IMPORT, PHIS, 'PHIS', RC=STATUS )
    VERIFY_(STATUS)

    phisxy = real(phis,kind=r8)

! Compute Pressure Thickness
! --------------------------

    dp = ( vars%pe(:,:,2:) - vars%pe (:,:,:km) )

! Create A-Grid Winds
! -------------------
    call ctoa_winds ( vars%u,vars%v,ua,va,                               &
                      g3_grid%dlam,g3_grid%dphi,im,jm,km,g3_grid%lattice )

! Specific humidity before and after physics updates
! --------------------------------------------------

    call MAPL_GetPointer(export,QOLD,'Q',  rc=status)

    call PULL_Q ( STATE, IMPORT, qqq, NXQ, rc )
      
    if ( qqq%is_r4 ) then
         qv = qqq%content_r4
    else
         qv = qqq%content
    endif

! Compute Energetics Before Diabatic Forcing
! ------------------------------------------

    thv = vars%pt*(1.0+eps*QOLD)

    call Energetics (state,vars%u,vars%v,thv,vars%pe,dp,vars%pkz,phisxy,kenrg0,penrg0,tenrg0,g3_grid)

! DTHVDTPHYINT
! ------------
      call MAPL_GetPointer ( export, temp2D, 'DTHVDTPHYINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dthdtphyint1 = 0.0
          do k=1,km
          dthdtphyint1 = dthdtphyint1 + thv(:,:,k)*dp(:,:,k)
          enddo
      endif

! Add Diabatic Forcing to State Variables
! ---------------------------------------

    call ADD_INCS ( STATE,IMPORT,DT,g3_grid )

! Update Mid-Layer Pressure and Pressure Thickness
! ------------------------------------------------

    dp = ( vars%pe(:,:,2:) - vars%pe (:,:,:km) )
    pl = ( vars%pe(:,:,2:) + vars%pe (:,:,:km) )*0.5

    logpl = log(pl)
    logpe = log(vars%pe)
    logps = log(vars%pe(:,:,km+1))
    
! Create A-Grid Winds
! -------------------
    call ctoa_winds ( vars%u,vars%v,ua,va,                               &
                      g3_grid%dlam,g3_grid%dphi,im,jm,km,g3_grid%lattice )

! Compute Energetics After Diabatic Forcing
! -----------------------------------------

    thv = vars%pt*(1.0+eps*qv)

    call Energetics (state,vars%u,vars%v,thv,vars%pe,dp,vars%pkz,phisxy,kenrg,penrg,tenrg,g3_grid)

    call MAPL_GetPointer(export,temp2d,'KE',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) temp2d = kenrg

    kenrg = (kenrg-kenrg0)/DT
    penrg = (penrg-penrg0)/DT
    tenrg = (tenrg-tenrg0)/DT

    call FILLOUT2 (export, 'KEPHY', kenrg, rc=status); VERIFY_(STATUS)
    call FILLOUT2 (export, 'PEPHY', penrg, rc=status); VERIFY_(STATUS)
    call FILLOUT2 (export, 'TEPHY', tenrg, rc=status); VERIFY_(STATUS)

! DTHVDTPHYINT
! ------------
      call MAPL_GetPointer ( export, temp2D, 'DTHVDTPHYINT', rc=status )
      VERIFY_(STATUS)
      if( associated(temp2D) ) then
          dthdtphyint2 = 0.0
          do k=1,km
          dthdtphyint2 = dthdtphyint2 + thv(:,:,k)*dp(:,:,k)
          enddo
          temp2D       = (dthdtphyint2-dthdtphyint1) * MAPL_P00**MAPL_KAPPA / (MAPL_GRAV*DT)
      endif

! Fill V_DGRID with Reasonable Pole Values Averaged from V_AGRID (For use by utilities outside Model)
! ---------------------------------------------------------------------------------------------------

     allocate( va_yz(imglobal,jfirst:jlast,kfirst:klast) )
     allocate( vd_yz(imglobal,jfirst:jlast,kfirst:klast) )

     if( grid%twod_decomp /= 0 ) then
         call mp_swapirr( commglobal, grid%ijk_xy_to_yz%SendDesc, &
                          grid%ijk_xy_to_yz%RecvDesc, va, va_yz,  &
                          a2in=vars%v, a2out=vd_yz )

! Question: why should this be grid%ijk_xy_to_yz and not grid%vxy_to_v ??  Ghosting??
!!!         call mp_sendirr( vars%v, grid%ijk_xy_to_yz%SendDesc, grid%ijk_xy_to_yz%RecvDesc, vd_yz )
!!!         call mp_recvirr(                              vd_yz, grid%ijk_xy_to_yz%RecvDesc        )
     else
         do k=1,km
         do j=jfirst,jlast
         do i=1,imglobal
            va_yz(i,j,k) =     va(i,j,k)
            vd_yz(i,j,k) = vars%v(i,j,k)
         enddo
         enddo
         enddo
     endif

     if ( jfirst == 1 ) then
          do k=kfirst,klast
          im1 =  imglobal
          do i=1,imglobal
          vd_yz(i,jfirst,k) = 0.5_r8*( va_yz(i,jfirst,k)+va_yz(im1,jfirst,k) )
          im1 = i
          enddo
          enddo
     endif
     if ( jlast == jmglobal ) then
          do k=kfirst,klast
          im1 =  imglobal
          do i=1,imglobal
          vd_yz(i,jlast,k) = 0.5_r8*( va_yz(i,jlast,k)+va_yz(im1,jlast,k) )
          im1 = i
          enddo
          enddo
     endif

     if( grid%twod_decomp /= 0 ) then
         call mp_swapirr( commglobal, grid%ijk_yz_to_xy%SendDesc, &
                          grid%ijk_yz_to_xy%RecvDesc, vd_yz, vars%v )
     else
         do k=1,km
         do j=jfirst,jlast
         do i=1,imglobal
            vars%v(i,j,k) = vd_yz(i,j,k)
         enddo
         enddo
         enddo
     endif

     deallocate( va_yz )
     deallocate( vd_yz )

! Fill Diagnostics
! ----------------

    tempxy = vars%pt * vars%pkz   ! Dry Temperature

    call FILLOUT3 (export, 'DELP'   , dp      , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'U'      , ua      , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'V'      , va      , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'T'      , tempxy  , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'Q'      , qv      , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PL'     , pl      , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PLE'    , vars%pe , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PLK'    , vars%pkz, rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'THV'    , thv     , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'U_DGRID', vars%u  , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'V_DGRID', vars%v  , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PT'     , vars%pt , rc=status); VERIFY_(STATUS)
    call FILLOUT3 (export, 'PE'     , vars%pe , rc=status); VERIFY_(STATUS)

    call MAPL_GetPointer(export,temp3d,'TH',rc=status)
    VERIFY_(STATUS)
    if(associated(temp3d)) temp3d = vars%pt*(p00**kappa)

! Compute Edge Heights
! --------------------

    pke           = vars%pe**kappa
    zle(:,:,km+1) = phisxy(:,:)
    do k=km,1,-1
       zle(:,:,k) = zle(:,:,k+1) + cp*thv(:,:,k)*( pke(:,:,k+1)-pke(:,:,k) )
    enddo
       zle(:,:,:) = zle(:,:,:)/grav

    call FILLOUT3 (export, 'ZLE', zle, rc=status); VERIFY_(STATUS)

! Compute Mid-Layer Heights
! -------------------------

    call MAPL_GetPointer(export,temp3d,'ZL',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp3d)) temp3d = 0.5*( zle(:,:,2:) + zle(:,:,:km) )

! Fill Single Level Variables
! ---------------------------

    call MAPL_GetPointer(export,temp2d,'U250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,ua,logpl,log(25000.),logps,status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'U500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,ua,logpl,log(50000.),logps,status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'U850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,ua,logpl,log(85000.),logps,status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,va,logpl,log(25000.),logps,status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,va,logpl,log(50000.),logps,status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,va,logpl,log(85000.),logps,status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,tempxy,logpl,log(25000.),logps,status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,tempxy,logpl,log(50000.),logps,status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'T850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,tempxy,logpl,log(85000.),logps,status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Q250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,qv,logpl,log(25000.),logps,status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Q500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,qv,logpl,log(50000.),logps,status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'Q850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,qv,logpl,log(85000.),logps,status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H250',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,zle,logpe,log(25000.),rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H500',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,zle,logpe,log(50000.),rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H850',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,zle,logpe,log(85000.),rc=status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'H1000',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,zle,logpe,log(100000.),rc=status)
       VERIFY_(STATUS)
    end if

! Compute Mid-Level Heights Above Surface
! ---------------------------------------
    do k=1,km
    zle(:,:,k) = 0.5*( zle(:,:,k)+zle(:,:,k+1) ) - zle(:,:,km+1)
    enddo
    zle(:,:,km+1) = 0.0
    
    call MAPL_GetPointer(export,temp2d,'U50M',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,ua,zle(:,:,1:km),50.,zle(:,:,km+1),status)
       VERIFY_(STATUS)
    end if

    call MAPL_GetPointer(export,temp2d,'V50M',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       call VertInterp(temp2d,va,zle(:,:,1:km),50.,zle(:,:,km+1),status)
       VERIFY_(STATUS)
    end if

! Compute Surface Pressure
! ------------------------

    call MAPL_GetPointer(export,temp2d,'PS',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) temp2d=vars%pe(:,:,km+1)

! Compute Vertically Averaged T,U
! -------------------------------
    call MAPL_GetPointer(export,temp2d,'TAVE',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       temp2d = 0.0
       do k=1,km
       temp2d = temp2d + tempxy(:,:,k)*dp(:,:,k)
       enddo
       temp2d = temp2d / (vars%pe(:,:,km+1)-vars%pe(:,:,1))
    endif

    call MAPL_GetPointer(export,temp2d,'UAVE',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp2d)) then
       temp2d = 0.0
       do k=1,km
       temp2d = temp2d + ua(:,:,k)*dp(:,:,k)
       enddo
       temp2d = temp2d / (vars%pe(:,:,km+1)-vars%pe(:,:,1))
    endif

! Convert T to Tv
! ---------------

    tempxy = tempxy*(1.0+eps*qv)

    call MAPL_GetPointer(export,temp3d,'TV',  rc=status)
    VERIFY_(STATUS)
    if(associated(temp3d)) temp3d=tempxy

! Compute Sea-Level Pressure
! --------------------------
    
    call MAPL_GetPointer(export,temp2d,'SLP'  ,rc=status)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,Ztemp1,'H1000',rc=status)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,Ztemp2,'H850' ,rc=status)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,Ztemp3,'H500' ,rc=status)
    VERIFY_(STATUS)

    if(associated(temp2d) .or. associated(ztemp1) &
                          .or. associated(ztemp2) &
                          .or. associated(ztemp3) ) then
       ALLOCATE(  slp(ifirstxy:ilastxy,jfirstxy:jlastxy) )
       ALLOCATE(H1000(ifirstxy:ilastxy,jfirstxy:jlastxy) )
       ALLOCATE(H850 (ifirstxy:ilastxy,jfirstxy:jlastxy) )
       ALLOCATE(H500 (ifirstxy:ilastxy,jfirstxy:jlastxy) )
       do j=jfirstxy,jlastxy
          do i=ifirstxy,ilastxy
             call get_slp ( km,vars%pe (i,j,  km+1),phisxy(i,j),  slp(i,j), &
                               vars%pe (i,j,1:km+1),                        &
                               vars%pkz(i,j,1:km  ),                        &
                                 tempxy(i,j,1:km  ),                        &
                                  H1000(i,j), H850(i,j), H500(i,j)          )
          enddo
       enddo
       if(associated(temp2d)) temp2d = slp
       if(associated(ztemp1)) where( ztemp1.eq.MAPL_UNDEF ) ztemp1 = H1000
       if(associated(ztemp2)) where( ztemp2.eq.MAPL_UNDEF ) ztemp2 = H850
       if(associated(ztemp3)) where( ztemp3.eq.MAPL_UNDEF ) ztemp3 = H500
       DEALLOCATE(slp,H1000,H850,H500)
    end if

! Deallocate Memory
! -----------------

    DEALLOCATE(  kenrg )
    DEALLOCATE(  penrg )
    DEALLOCATE(  tenrg )
    DEALLOCATE( kenrg0 )
    DEALLOCATE( penrg0 )
    DEALLOCATE( tenrg0 )

    DEALLOCATE( phisxy )

    DEALLOCATE(     ua )
    DEALLOCATE(     va )
    DEALLOCATE(     qv )
    DEALLOCATE(     pl )
    DEALLOCATE(     dp )
    DEALLOCATE( tempxy )

    DEALLOCATE(    thv )
    DEALLOCATE(    pke )
    DEALLOCATE(  logpl )
    DEALLOCATE(  logpe )
    DEALLOCATE(  logps )
    DEALLOCATE(    zle )
    DEALLOCATE( dthdtphyint1 )
    DEALLOCATE( dthdtphyint2 )

    DEALLOCATE( STATE%VARS%tracer, STAT=STATUS )   ! Allocated by call to PULL_Q

    call MAPL_TimerOff(MAPL,"RUN2")
    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)
end subroutine Run2

!-----------------------------------------------------------------------
  subroutine ADD_INCS ( STATE,IMPORT,DT,g3_grid,RC )

    use g3_dynamics_state_module

    include 'mpif.h'
    integer imglobal, jmglobal
!
! !INPUT PARAMETERS:

   type(T_FVDYCORE_STATE), intent(INOUT)  :: STATE
   type(ESMF_State),       intent(INOUT)  :: IMPORT
   real(r8),               intent(IN   )  :: DT
   type ( dynamics_grid_type )            :: g3_grid
   integer,  optional,     intent(OUT  )  :: RC

!
! !DESCRIPTION:  This routine adds the tendencies to the state,
!                weighted appropriately by the time step.  Temperature
!                tendencies are pressure weighted (ie., DELP*DT/Dt).
!                All tendencies are on the A-grid, and have an XY decomposition.
!

    integer                          :: status

    integer               :: I1, IN, J1, JN, K, im, jm, km
    integer               :: KL, KU, NX, NY
    real(r8)              :: SUMOUT
    real(r8), allocatable ::    dum(:,:,:), pkzold(:,:,:)
    real(r8), allocatable ::    pke(:,:,:),  dpinv(:,:,:)
    real(r8), allocatable :: tend_u(:,:,:), tend_v(:,:,:)
    real(kind=4), pointer :: tend(:,:,:)

    character(len=ESMF_MAXSTR)         :: IAm="ADD_INCS"
    character*2 cnx,cny

    i1 = state%grid%ifirstxy
    in = state%grid%ilastxy
    j1 = state%grid%jfirstxy
    jn = state%grid%jlastxy

    imglobal = state%grid%im
    jmglobal = state%grid%jm
    km       = state%grid%km

    im       = g3_grid%im
    jm       = g3_grid%jm

! **********************************************************************
! ****                       Update Winds                           ****
! **********************************************************************

    ALLOCATE( tend_u(i1:in,j1:jn,km) )
    ALLOCATE( tend_v(i1:in,j1:jn,km) )

    call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DUDT',RC=STATUS )
    VERIFY_(STATUS)
    tend_u = tend

    call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DVDT',RC=STATUS )
    VERIFY_(STATUS)
    tend_v = tend

    _ASSERT( im == in-i1+1 ,'needs informative message')
    _ASSERT( jm == jn-j1+1 ,'needs informative message')

! Put the wind tendencies on the C-grid
! -------------------------------------
    call atoc ( tend_u,tend_u,g3_grid%dlam,g3_grid%dphi,im,jm,km,1,g3_grid%lattice )
    call atoc ( tend_v,tend_v,g3_grid%dlam,g3_grid%dphi,im,jm,km,2,g3_grid%lattice )

! Add the wind tendencies to the control variables
! ------------------------------------------------
    STATE%VARS%U = STATE%VARS%U + DT*TEND_U
    STATE%VARS%V = STATE%VARS%V + DT*TEND_V

! Set D-GRID U at the South Pole to UNDEF
! ---------------------------------------
!   if ( j1 == 1 ) STATE%VARS%U(:,1,:) = MAPL_UNDEF

! Set D-GRID V at Both Poles to UNDEF
! -----------------------------------
!   if ( j1 == 1  ) STATE%VARS%V(:, 1,:) = MAPL_UNDEF
!   if ( jn == jm ) STATE%VARS%V(:,jm,:) = MAPL_UNDEF

    DEALLOCATE( tend_u )
    DEALLOCATE( tend_v )

! **********************************************************************
! ****        Compute Pressure Thickness Using Old Pressures        ****
! **********************************************************************

    ALLOCATE( dpinv(i1:in,j1:jn,km) )
    do k=1,km
       dpinv(:,:,k) = 1.0/( state%vars%pe(:,:,k+1)-state%vars%pe(:,:,k) )
    enddo

! **********************************************************************
! ****                     Update Edge Pressures                    ****
! **********************************************************************

    call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DPEDT',RC=STATUS )
    VERIFY_(STATUS)

    KL = lbound( tend,3 )
    KU = ubound( tend,3 )

    allocate( dum(i1:in,j1:jn,KL:KU) )

    DUM = DT*TEND

! Constrain Poles
! ---------------
    if ( j1 == 1 ) then
    do k=KL,KU
         call par_xsum ( state%grid,  DUM(i1:in,1,k), 1, sumout )
         sumout = sumout/imglobal
         do i=i1,in
         DUM(i,1,k) = sumout
         enddo
    enddo
    endif
    if ( jn == jmglobal ) then
    do k=KL,KU
         call par_xsum ( state%grid,  DUM(i1:in,jn,k), 1, sumout )
         sumout = sumout/imglobal
         do i=i1,in
         DUM(i,jn,k) = sumout
         enddo
    enddo
    endif

    STATE%VARS%PE = STATE%VARS%PE + DUM
    DEALLOCATE (DUM)

! **********************************************************************
! ****              Update P*Kappa at Mid-Levels                    ****
! **********************************************************************
  
    ALLOCATE( pke   (i1:in,j1:jn,km+1) )
    ALLOCATE( pkzold(i1:in,j1:jn,1:km) )

    pke    = STATE%VARS%PE**kappa
    pkzold = STATE%VARS%PKZ

    do k=1,km
    STATE%VARS%PKZ(:,:,k) = ( pke(:,:,k+1)-pke(:,:,k) ) &
                          / ( kappa*( log(STATE%VARS%PE(:,:,k+1))-log(STATE%VARS%PE(:,:,k)) ) )
    enddo
 
! *********************************************************************
! ****                Update Dry Potential Temperature             ****
! ****                --------------------------------             ****
! ****  Note: State Variable is Potential Temperature T/P**kappa   ****
! ****             while IMPORT Coupling is (Delta_P)*DTDt         ****
! *********************************************************************

    call ESMFL_StateGetPointerToData ( IMPORT,TEND,'DTDT',RC=STATUS )
    VERIFY_(STATUS)

    KL = lbound( tend,3 )
    KU = ubound( tend,3 )

    allocate( dum(i1:in,j1:jn,KL:KU) )

    DUM = DT*TEND*DPINV/STATE%VARS%PKZ                   &
        +    STATE%VARS%PT*( PKZOLD/STATE%VARS%PKZ - 1.0 )

! Constrain Poles
! ---------------
    if ( j1 == 1 ) then
    do k=KL,KU
         call par_xsum ( state%grid,  DUM(i1:in,1,k), 1, sumout )
         sumout = sumout/imglobal
         do i=i1,in
         DUM(i,1,k) = sumout
         enddo
    enddo
    endif
    if ( jn == jmglobal ) then
    do k=KL,KU
         call par_xsum ( state%grid,  DUM(i1:in,jn,k), 1, sumout )
         sumout = sumout/imglobal
         do i=i1,in
         DUM(i,jn,k) = sumout
         enddo
    enddo
    endif

    STATE%VARS%PT =  STATE%VARS%PT              + DUM
!   STATE%VARS%PT = (STATE%VARS%PT*(1+EPS*QOLD) + DUM)/(1+EPS*QNEW)
    DEALLOCATE (DUM)

    DEALLOCATE( PKE    )
    DEALLOCATE( PKZOLD )
    DEALLOCATE( DPINV  )

   return

 end subroutine ADD_INCS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Finalize

! !DESCRIPTION: Writes restarts and cleans-up through MAPL\_GenericFinalize and
!   deallocates memory from the Private Internal state. 
!
! !INTERFACE:

subroutine Finalize(gc, import, export, clock, rc)
   use dynamics_vars, only : dynamics_clean

! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: gc
    type (ESMF_State),    intent(inout) :: import
    type (ESMF_State),    intent(inout) :: export
    type (ESMF_Clock),    intent(inout) :: clock
    integer, optional,    intent(  out) :: rc
 
!EOP

! Local variables
    type (DYN_wrap) :: wrap
    type (T_FVDYCORE_STATE), pointer  :: STATE
    character (len=ESMF_MAXSTR)       :: restart_file
 
    character(len=ESMF_MAXSTR)        :: IAm
    character(len=ESMF_MAXSTR)        :: COMP_NAME
    integer                           :: status

    type (MAPL_MetaComp),     pointer :: MAPL 
    type (ESMF_Config)                :: cf


! BEGIN

    Iam = "Finalize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, config=cf, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
! ---------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"FINALIZE")

! Retrieve the pointer to the state
!----------------------------------

    call ESMF_UserCompGetInternalState(gc, 'FVstate', wrap, status)
    VERIFY_(STATUS)
  
    state => wrap%dyn_state
 
    call dynamics_clean    (STATE%GRID)
 
! Call Generic Finalize
!----------------------

    call MAPL_TimerOff(MAPL,"FINALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL")

    call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine FINALIZE


!=======================================================================

      subroutine get_slp ( km,ps,phis,slp,pe,pk,tv,H1000,H850,H500)
      implicit   none
      integer  km
      real(r8)   pk(km)    ! layer-mean P**kappa
      real(r8)   tv(km)    ! layer-mean virtual Temperature
      real(r8)   pe(km+1)  ! press at layer edges (Pa)
      real(r8)   ps        ! surface pressure (Pa)
      real(r8) phis        ! surface geopotential
      real(r8)  slp        ! sea-level pressure (hPa)
      real(r8)  H1000      ! 1000mb height
      real(r8)  H850       !  850mb height
      real(r8)  H500       !  500mb height
      real(r8)  tstar                 ! extrapolated temperature (K)
      real(r8) p_bot
      real(r8) tref                   ! Reference virtual temperature (K)
      real(r8) pref                   ! Reference pressure level (Pa)
      real(r8) pkref                  ! Reference pressure level (Pa) ** kappa
      real(r8) dp1, dp2

      real(r8), parameter :: gamma    = 6.5e-3
      real(r8), parameter :: p_offset = 15000.
      real(r8), parameter :: gg       = gamma/MAPL_GRAV

      real(r8), parameter :: factor   = MAPL_grav / ( MAPL_Rgas * gamma ) 
      real(r8), parameter :: yfactor  = MAPL_Rgas * gg

      integer k_bot, k, k1, k2

      p_bot = ps - p_offset
      k_bot = -1

      do k = km, 2, -1
         if ( pe(k+1) .lt. p_bot ) then
              k_bot = k
              exit
         endif
      enddo

      k1    = k_bot - 1
      k2    = k_bot
      dp1   = pe(k_bot)   - pe(k_bot-1)
      dp2   = pe(k_bot+1) - pe(k_bot)
      pkref = ( pk(k1)*dp1 + pk(k2)*dp2 ) / (dp1+dp2)
       tref = ( tv(k1)*dp1 + tv(k2)*dp2 ) / (dp1+dp2)
       pref = 0.5 * ( pe(k_bot+1) + pe(k_bot-1) )
      tstar = tref*( ps/pref )**yfactor

      slp   = ps*( 1.0+gg*phis/tstar )**factor
      H1000 = (phis/MAPL_grav) - (tstar/gamma)*((100000.0/ps)**(1./factor)-1.0)
      H850  = (phis/MAPL_grav) - (tstar/gamma)*(( 85000.0/ps)**(1./factor)-1.0)
      H500  = (phis/MAPL_grav) - (tstar/gamma)*(( 50000.0/ps)**(1./factor)-1.0)
      return
  end subroutine get_slp

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine VertInterp(v2,v3,pl,pp,ps,rc)

    real(r4),           intent(OUT) :: v2(:,:)
    real(r8),           intent(IN ) :: v3(:,:,:)
    real(r8),  target,  intent(IN ) :: pl(:,:,:)
    real,               intent(IN ) :: pp
    real(r8), optional, intent(IN ) :: ps(:,:)
    integer,  optional, intent(OUT) :: rc

    real, dimension(size(v2,1),size(v2,2)) :: al,PT,PB
    integer km, K, msn
    logical flip
    real    ppx
    real(r8), pointer   :: plx(:,:,:)

    integer        :: status
    character*(10) :: Iam='VertInterp'

    km   = size(pl,3)

    flip = pl(1,1,km) < pl(1,1,km-1)

    if(flip) then
       allocate(plx(size(pl,1),size(pl,2),size(pl,3)),stat=status)
       VERIFY_(STATUS)
       plx = -pl
       ppx = -pp
       msn = -1
    else
       plx => pl
       ppx = pp
       msn = 1
    end if

    v2   = MAPL_UNDEF

       pb   = plx(:,:,km)
       do k=km-1,1,-1
          pt = plx(:,:,k)
          if(all(pb<ppx)) exit
          where(pt<ppx .and. pb>=ppx)
             al = (pb-ppx)/(pb-pt)
             v2 = v3(:,:,k)*al + v3(:,:,k+1)*(1.0-al)
          end where
          pb = pt
       end do

! Extend Lowest Level Value to the Surface
! ----------------------------------------
    if( present(ps) ) then
        where( (plx(:,:,km)<ppx.and.ps*msn>=ppx) )
                v2 = v3(:,:,km)
        end where
    end if

    if(flip) then
       deallocate(plx,stat=status)
       VERIFY_(STATUS)
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine VertInterp

!BOP

! !IROUTINE: Coldstart

! !DESCRIPTION:
!   Routine to coldstart from an isothermal state of rest.
!   The temperature can be specified in the config, otherwise
!   it is 300K. The surface pressure is assumed to be 1000 hPa.
!
! !INTERFACE:

subroutine Coldstart(gc, import, export, clock, rc)

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc
    type(ESMF_State),    intent(inout) :: import
    type(ESMF_State),    intent(inout) :: export
    type (ESMF_Clock),   intent(in)    :: clock
    integer, intent(out), optional     :: rc
 
!EOP

    character(len=ESMF_MAXSTR)        :: IAm
    character(len=ESMF_MAXSTR)        :: COMP_NAME
    integer                           :: status

    type (MAPL_MetaComp),     pointer :: MAPL 
    type (ESMF_State)                 :: INTERNAL

    real(r8), pointer                 :: AK(:), BK(:)
    real(r8), pointer                 :: Ptr3(:,:,:)
    real(r8), pointer                 :: PKL (:,:,:)
    real, pointer                     :: LATS (:,:)
    real                              :: T0
    integer                           :: L
    type(ESMF_Config)                 :: CF

! Begin

    Iam = "Coldstart"
    call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the state
! ---------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")

!BOR    
! !RESOURCE_ITEM: K :: Value of isothermal temperature on coldstart
    call MAPL_GetResource ( MAPL, T0, 'T0:', default=300., RC=STATUS )
    VERIFY_(STATUS)
!EOR

    call MAPL_Get ( MAPL,                &
           INTERNAL_ESMF_STATE=INTERNAL, &
           lats = LATS,                  &
                               RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetPointer(Internal,Ptr3,'U'  ,rc=STATUS)
    VERIFY_(STATUS)
    Ptr3 = 0.0

    Ptr3(1,:,ubound(Ptr3,3)) = .001*abs(lats(1,:))

    call MAPL_GetPointer(Internal,Ptr3,'V'  ,rc=STATUS)
    VERIFY_(STATUS)
    Ptr3 = 0.0

    call MAPL_GetPointer(Internal,Ptr3,'PE',rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(Internal,PKL ,'PKZ',rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetPointer(Internal,ak  ,'AK' ,rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(Internal,bk  ,'BK' ,rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigFindLabel( cf, 'AK:', rc = status )
    VERIFY_(STATUS)
    do L = 0, SIZE(AK)-1
       call ESMF_ConfigNextLine  ( CF, rc=STATUS )
       call ESMF_ConfigGetAttribute( cf, AK(L), rc = status )
       VERIFY_(STATUS)
    enddo

    call ESMF_ConfigFindLabel( cf, 'BK:', rc = status )
    VERIFY_(STATUS)
    do L = 0, SIZE(bk)-1
       call ESMF_ConfigNextLine  ( CF, rc=STATUS )
       call ESMF_ConfigGetAttribute( cf, BK(L), rc = status )
       VERIFY_(STATUS)
    enddo

   _ASSERT(ANY(AK /= 0.0) .or. ANY(BK /= 0.0),'needs informative message')
    do L=lbound(Ptr3,3),ubound(Ptr3,3)
       Ptr3(:,:,L) = AK(L) + BK(L)*MAPL_P00
    enddo
 
    PKL = 0.5*(Ptr3(:,:,lbound(Ptr3,3)  :ubound(Ptr3,3)-1) + &
               Ptr3(:,:,lbound(Ptr3,3)+1:ubound(Ptr3,3)  ) )
    PKL = PKL**MAPL_KAPPA

    call MAPL_GetPointer(Internal,Ptr3,'PT',rc=STATUS)
    VERIFY_(STATUS)

    Ptr3 = T0/PKL

    call MAPL_TimerOff(MAPL,"TOTAL")


    RETURN_(ESMF_SUCCESS)
  end subroutine COLDSTART

end module ARIESg3_GridCompMod
