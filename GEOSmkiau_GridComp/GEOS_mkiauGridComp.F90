!  $Id: GEOS_mkiauGridComp.F90,v 1.38.2.21.18.5.2.5 2019/10/22 20:53:09 ltakacs Exp $

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_mkiau -- A Module to compute the IAU forcing

! !INTERFACE:

module GEOS_mkiauGridCompMod

! !USES:

  use ESMF
  use MAPL
  use ESMF_CFIOFileMod
  use GEOS_UtilsMod
! use GEOS_RemapMod, only: myremap => remap
  use m_set_eta, only: set_eta
#ifdef PYMLINC_INTEGRATION
  use pyMLINC_interface_mod
  use ieee_exceptions, only: ieee_get_halting_mode, ieee_set_halting_mode, ieee_all
#endif
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================
  type T_MKIAU_STATE
     private
     class (AbstractRegridder), pointer :: ANA2BKG_regridder => null()
     class (AbstractRegridder), pointer :: BKG2ANA_regridder => null()
     type(ESMF_Grid)            :: GRIDana    ! Analysis    Data using Horizontal:ANA  Vertical:BKG
     type(ESMF_Grid)            :: GRIDrep    ! Replay File Data using Horizontal:ANA  Vertical:ANA
     integer                    :: IM
     integer                    :: JM
     integer                    :: LM
  end type T_MKIAU_STATE

! Wrapper for extracting internal state
! -------------------------------------
  type MKIAU_WRAP
     type (T_MKIAU_STATE), pointer :: PTR
  end type MKIAU_WRAP


!=============================================================================

! !DESCRIPTION:
!
!

!EOP

contains

!BOP

! ! IROUTINE: SetServices -- Sets ESMF services for this component

! ! INTERFACE:

  subroutine SetServices ( GC, RC )

! ! ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! ! DESCRIPTION: This version uses the MAPL_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF_State INTERNAL, which is in the MAPL_MetaComp.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME
    type (MAPL_MetaComp),         pointer   :: MAPL
    type (T_MKIAU_STATE),         pointer   :: mkiau_internal_state
    type (MKIAU_wrap)                       :: wrap
    type (ESMF_Config)                      :: CF

    logical                                 :: BLEND_AT_PBL
#ifdef PYMLINC_INTEGRATION
    ! IEEE trapping see below
    logical                                 :: halting_mode(5)
    ! BOGUS DATA TO SHOW USAGE
    type(a_pod_struct_type) :: options
    real, allocatable, dimension(:,:,:) :: in_buffer
    real, allocatable, dimension(:,:,:) :: out_buffer
#endif
  !=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
!----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, BLEND_AT_PBL,    LABEL="REPLAY_BLEND_AT_PBL:",   default=.FALSE., RC=status)
    VERIFY_(STATUS)

! Set the Run entry points (phase 1 for regular IAU and phase 2 for clearing
! --------------------------------------------------------------------------

    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,  Run,  &
                                      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,  Clear,  &
                                      RC=STATUS)
    VERIFY_(STATUS)


! Set the state variable specs.
! -----------------------------

! !IMPORT STATE:
! --------------
    call MAPL_AddImportSpec ( GC,                                  &
         SHORT_NAME = 'PHIS',                                      &
         LONG_NAME  = 'surface geopotential height',               &
         UNITS      = 'm2 sec-2',                                  &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( GC,                                  &
         SHORT_NAME = 'AK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_a',                   &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( GC,                                  &
         SHORT_NAME = 'BK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_b',                   &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'PS',                                        &
         LONG_NAME  = 'surface_air_pressure',                      &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'DELP',                                      &
         LONG_NAME  = 'air_pressure_thickness',                    &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'TV',                                        &
         LONG_NAME  = 'virtual_air_temperature',                   &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'O3PPMV',                                    &
         LONG_NAME  = 'ozone_volume_mixing_ratio',                 &
         UNITS      = 'ppmv',                                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'TS',                                        &
         LONG_NAME  = 'skin_temperature',                          &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME = 'QV',                                        &
         LONG_NAME  = 'water_vapor_specific_humdity',              &
         UNITS      = 'kg/kg',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    if( BLEND_AT_PBL ) then
    call MAPL_AddImportSpec(GC,                                        &
         SHORT_NAME = 'PPBL',                                          &
         LONG_NAME  = 'pbl_top_pressure',                              &
         UNITS      = 'Pa',                                            &
         DIMS       = MAPL_DimsHorzOnly,                               &
         VLOCATION  = MAPL_VLocationNone,                              &
         RC=STATUS  )
    VERIFY_(STATUS)
    endif


! !EXPORT STATE:
! --------------
    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_analysis_increment',          &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_analysis_increment',         &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'temperature_analysis_increment',            &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'edge_pressure_analysis_increment',          &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DQVDT',                                     &
         LONG_NAME  = 'specific_humidity_analysis_increment',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DO3DT',                                     &
         LONG_NAME  = 'ozone_analysis_increment',                  &
         UNITS      = 'ppmv',                                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DTSDT',                                     &
         LONG_NAME  = 'skin_temparature_increment',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DUWINDFIX',                                 &
         LONG_NAME  = 'eastward_windfix_increment',                &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DVWINDFIX',                                 &
         LONG_NAME  = 'northward_windfix_increment',               &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'PSBKG',                                     &
         LONG_NAME  = 'surface_air_pressure_of_background',        &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'DELPBKG',                                   &
         LONG_NAME  = 'air_pressure_thickness_of_background',      &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'TVBKG',                                     &
         LONG_NAME  = 'virtual_air_temperature_of_background',     &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'UBKG',                                      &
         LONG_NAME  = 'eastward_wind_of_background',               &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'VBKG',                                      &
         LONG_NAME  = 'northward_wind_of_background',              &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'QVBKG',                                     &
         LONG_NAME  = 'water_vapor_specific_humdity_of_background',&
         UNITS      = 'kg/kg',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'O3PPMVBKG',                                 &
         LONG_NAME  = 'ozone_volume_mixing_ratio_of_background',   &
         UNITS      = 'ppmv',                                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'TSBKG',                                     &
         LONG_NAME  = 'skin_temperature_of_background',            &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                                &
         SHORT_NAME = 'VINTDIV_ANA',                                                             &
         LONG_NAME  = 'vertically_integrated_mass_divergence_increment_from_analysis',           &
         UNITS      = 'Pa s-1',                                                                  &
         DIMS       = MAPL_DimsHorzOnly,                                                         &
         VLOCATION  = MAPL_VLocationNone,                                                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                                &
         SHORT_NAME = 'VINTDIV_BKG',                                                             &
         LONG_NAME  = 'vertically_integrated_mass_divergence_increment_from_background',         &
         UNITS      = 'Pa s-1',                                                                  &
         DIMS       = MAPL_DimsHorzOnly,                                                         &
         VLOCATION  = MAPL_VLocationNone,                                                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                                                &
         SHORT_NAME = 'VINTDIV_COR',                                                             &
         LONG_NAME  = 'vertically_integrated_mass_divergence_increment_from_analysis_corrected', &
         UNITS      = 'Pa s-1',                                                                  &
         DIMS       = MAPL_DimsHorzOnly,                                                         &
         VLOCATION  = MAPL_VLocationNone,                                                        &
         RC=STATUS  )
    VERIFY_(STATUS)


! Internal State (None)
! ---------------------


! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC, name="INITIALIZE" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="DRIVER"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-INTR"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-RUN"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-REGRIDSPC"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--WINDFIX"    ,RC=STATUS)
    VERIFY_(STATUS)

! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( mkiau_internal_state, stat=status )
    VERIFY_(STATUS)
    wrap%ptr => mkiau_internal_state

! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'MKIAU_state', wrap, status )
    VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( gc, RC=STATUS)
    VERIFY_(STATUS)

#ifdef PYMLINC_INTEGRATION
    ! Spin the interface - we have to deactivate the ieee error
    ! to be able to load numpy, scipy and other numpy packages
    ! that generate NaN as an init mechanism for numerical solving
    call ieee_get_halting_mode(ieee_all, halting_mode)
    call ieee_set_halting_mode(ieee_all, .false.)
    call pyMLINC_interface_f_setservice()
    call ieee_set_halting_mode(ieee_all, halting_mode)

    ! BOGUS CODE TO SHOW USAGE
    options%npx = 10
    options%npy = 11
    options%npz = 12
    allocate (in_buffer(10,11,12), source = 42.42 )
    allocate (out_buffer(10,11,12), source = 0.0 )
    call pyMLINC_interface_f_run(options, in_buffer, out_buffer)
    write(*,*) "[pyMLINC] From fortran OUT[5,5,5] is ", out_buffer(5,5,5)
#endif

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! ! IROUTINE: RUN -- Run method for the MAKEIAU component

! !INTERFACE:

subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: This version uses the MAPL_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),     pointer   :: MAPL
  type (ESMF_State       )            :: INTERNAL

  integer                             :: IM,    JM,    LM
  integer                             :: IMbkg, JMbkg, LMbkg
  integer                             :: IMana, JMana, LMana

  integer                             :: IMbkg_World, JMbkg_World
  integer                             :: IMana_World, JMana_World

  real, pointer, dimension(:,:,:)     :: uptr3d
  real, pointer, dimension(:,:,:)     :: vptr3d
  real, pointer, dimension(:,:,:)     ::  ptr3d, temp3d
  real, pointer, dimension(:,:)       ::  ptr2d, temp2d

! Background Variables from IMPORT State
! --------------------------------------
  real, pointer, dimension(:,:)       :: vintdiv_ana
  real, pointer, dimension(:,:)       :: vintdiv_bkg
  real, pointer, dimension(:,:)       :: vintdiv_cor
  real, pointer, dimension(:,:,:)     ::   u_bkg, du, duwindfix
  real, pointer, dimension(:,:,:)     ::   v_bkg, dv, dvwindfix
  real, pointer, dimension(:,:,:)     ::  tv_bkg, dt, t_bkg
  real, pointer, dimension(:,:,:)     ::   q_bkg, dq
  real, pointer, dimension(:,:,:)     ::  o3_bkg, do3
  real, pointer, dimension(:,:,:)     :: ple_bkg, dple
  real, pointer, dimension(:,:)       ::  ts_bkg, dts
  real, pointer, dimension(:,:)       ::  ps_bkg
  real, pointer, dimension(:,:)       ::phis_bkg
  real, pointer, dimension(:)         ::  ak,bk

  real, allocatable, dimension(:,:,:) ::  dp_bkg

! Analysis Variables from REPLAY Files
! ------------------------------------
  real, allocatable, dimension(:,:)   ::  ps_rep
  real, allocatable, dimension(:,:,:) ::  dp_rep
  real, allocatable, dimension(:,:,:) ::   u_rep
  real, allocatable, dimension(:,:,:) ::   v_rep
  real, allocatable, dimension(:,:,:) ::   t_rep
  real, allocatable, dimension(:,:,:) ::   q_rep
  real, allocatable, dimension(:,:,:) ::  o3_rep
  real, allocatable, dimension(:,:,:) :: thv_rep
  real, allocatable, dimension(:,:,:) :: ple_rep
  real, allocatable, dimension(:,:,:) ::  pk_rep
  real, allocatable, dimension(:,:,:) :: pke_rep
  real, allocatable, dimension(:)     ::  ak_rep
  real, allocatable, dimension(:)     ::  bk_rep

! Analysis Variables from REPLAY Files REMAPPED to Background Vertical Resolution
! -------------------------------------------------------------------------------
  real, pointer, dimension(:,:)       ::phis_ana
  real, pointer, dimension(:,:)       ::  ts_ana
  real, pointer, dimension(:,:)       ::  ps_ana
  real, pointer, dimension(:,:,:)     ::   u_ana
  real, pointer, dimension(:,:,:)     ::   v_ana
  real, pointer, dimension(:,:,:)     ::   t_ana
  real, pointer, dimension(:,:,:)     :: thv_ana
  real, pointer, dimension(:,:,:)     ::   q_ana
  real, pointer, dimension(:,:,:)     ::  o3_ana

  real, allocatable, dimension(:,:,:) :: ple_ana
  real, allocatable, dimension(:,:,:) ::  pk_ana
  real, allocatable, dimension(:,:,:) :: pke_ana
  real, allocatable, dimension(:,:,:) :: qdum1
  real, allocatable, dimension(:,:,:) :: qdum2

  real,     pointer, dimension(:,:,:) :: pdum1 => null()
  real,     pointer, dimension(:,:,:) :: pdum2 => null()
  real,     pointer, dimension(:,:)   :: blnpp => null()

  real, allocatable, dimension(:,:,:) ::  du_fix
  real, allocatable, dimension(:,:,:) ::  dv_fix

  real  ptopdum
  real  pintdum
  integer ksdum

  real, allocatable, dimension(:,:)   ::  vintdiva
  real, allocatable, dimension(:,:)   ::  vintdivb
  real, allocatable, dimension(:,:)   ::  vintdivc

  real, parameter                     :: EPS = MAPL_RVAP/MAPL_RGAS-1.0

  character(len=ESMF_MAXSTR), save    :: FILEP1
  character(len=ESMF_MAXSTR), save    :: FILEP0
  character(len=ESMF_MAXSTR), save    :: FILEM1
  character(len=ESMF_MAXSTR), save    :: FILEM2
  character(len=ESMF_MAXSTR)          :: REPLAY_FILEP1
  character(len=ESMF_MAXSTR)          :: REPLAY_FILEP0
  character(len=ESMF_MAXSTR)          :: REPLAY_FILEM1
  character(len=ESMF_MAXSTR)          :: REPLAY_FILEM2
  character(len=ESMF_MAXSTR)          :: REPLAY_TIME_INTERP
  character(len=ESMF_MAXSTR)          :: FILETMPL
  character(len=ESMF_MAXSTR)          :: GRIDINC
  character(len=ESMF_MAXSTR)          :: cremap
  character(len=ESMF_MAXSTR)          :: fixwind

  type(ESMF_FieldBundle), save        :: RBUNDLEP1
  type(ESMF_FieldBundle), save        :: RBUNDLEP0
  type(ESMF_FieldBundle), save        :: RBUNDLEM1
  type(ESMF_FieldBundle), save        :: RBUNDLEM2
  type(ESMF_Time),        save        :: FILE_TIMEP1
  type(ESMF_Time),        save        :: FILE_TIMEP0
  type(ESMF_Time),        save        :: FILE_TIMEM1
  type(ESMF_Time),        save        :: FILE_TIMEM2
  type(ESMF_Time)                     :: REPLAY_TIME
  type(ESMF_Time)                     :: REPLAY_TIMEP1
  type(ESMF_Time)                     :: REPLAY_TIMEP0
  type(ESMF_Time)                     :: REPLAY_TIMEM1
  type(ESMF_Time)                     :: REPLAY_TIMEM2

  type(ESMF_Alarm)                    :: ALARM
  logical                             :: is_RegularReplay09_ringing

  integer                                 :: K,NQ,FID
  character(len=ESMF_MAXSTR), ALLOCATABLE :: RNAMES(:)
  character(len=ESMF_MAXSTR)              :: STRING
  character(len=ESMF_MAXSTR)              :: DATE

  type(ESMF_FieldBundle)              :: bundle
  type(ESMF_Grid)                     :: GRIDbkg
  type(ESMF_Grid)                     :: GRIDana
  type(ESMF_Grid)                     :: GRIDrep
  type(ESMF_Time)                     :: currtime
  type(ESMF_Calendar)                 :: cal
  type(ESMF_TimeInterval)             :: FileFreq
  integer                             :: FileFreq_SEC
  integer                             :: FileReft_SEC, FileReft_HMS
  integer                             :: CUR_YY,CUR_MM,CUR_DD,CUR_H,CUR_M,CUR_S
  integer                             :: TOTAL_SEC

  real                                :: FACP1, FACP0, FACM1, FACM2
  real                                :: DAMPBEG, DAMPEND
  logical                             :: BLEND_AT_PBL
  integer                             :: i,j,L,n
  integer                             :: nt,nvars,natts
  integer                             :: nymd, nhms
  integer                             :: nymd1,nhms1
  integer                             :: nymd2,nhms2
  integer                             :: nymdp1,nhmsp1
  integer                             :: nymdp0,nhmsp0
  integer                             :: nymdm1,nhmsm1
  integer                             :: nymdm2,nhmsm2
  integer                             :: NX,NY,IMG,JMG
  integer                             :: method
  integer                             :: DIMS(ESMF_MAXGRIDDIM)
  integer                             :: JCAP,LMP1
  logical                             :: dowindfix
  logical                             :: doremap
  logical                             :: FOUND
  logical                             :: ANALYZE_TS
  character(len=ESMF_MAXSTR)          :: REPLAY_MODE
  character(len=ESMF_MAXSTR)          :: REPLAY_U, REPLAY_V,  REPLAY_T, REPLAY_QV, REPLAY_TS, REPLAY_O3
  character(len=ESMF_MAXSTR)          :: REPLAY_P, REPLAY_PS, REPLAY_DP
  character(len=ESMF_MAXSTR)          :: REPLAY_PHIS
  character(len=ESMF_MAXSTR)          :: REPLAY_T_TYPE
  logical                             :: L_REPLAY_PHIS, L_REPLAY_P, L_REPLAY_U, L_REPLAY_V, L_REPLAY_T, L_REPLAY_QV, L_REPLAY_TS, L_REPLAY_O3
  logical                             :: L_CUBE
  logical                             :: do_transforms
  logical                             :: USE_SPECFILT
  real                                :: REPLAY_P_FACTOR
  real                                :: REPLAY_U_FACTOR
  real                                :: REPLAY_V_FACTOR
  real                                :: REPLAY_T_FACTOR
  real                                :: REPLAY_QV_FACTOR
  real                                :: REPLAY_O3_FACTOR
  real                                :: REPLAY_TS_FACTOR

  class (AbstractRegridder), pointer :: ANA2BKG => null()
  class (AbstractRegridder), pointer :: BKG2ANA => null()
  integer                             :: NPHIS, NPHIS_MAX

  type (ESMF_VM)                      :: VM
  integer                             :: vm_comm
  integer                             :: IHAVEAINC

  type (T_MKIAU_STATE), pointer       :: mkiau_internal_state
  type (MKIAU_wrap)                   :: wrap
  logical                             :: refresh_internal_state
  logical                             :: bkg2anaConsrv
  logical                             :: ana2bkgConsrv
  character(len=ESMF_MAXSTR)          :: imstr, jmstr, gridAnaName

  logical                             :: first
  data                                   first /.true./
  logical                             :: NEED_BUNDLE1
  data                                   NEED_BUNDLE1 /.true./
  logical                             :: NEED_BUNDLE2
  data                                   NEED_BUNDLE2 /.true./
  logical                             :: NEED_BUNDLEP1
  data                                   NEED_BUNDLEP1 /.true./
  logical                             :: NEED_BUNDLEP0
  data                                   NEED_BUNDLEP0 /.true./
  logical                             :: NEED_BUNDLEM1
  data                                   NEED_BUNDLEM1 /.true./
  logical                             :: NEED_BUNDLEM2
  data                                   NEED_BUNDLEM2 /.true./

  integer nsecf
          nsecf(nhms) = nhms/10000*3600 + mod(nhms,10000)/100*60 + mod(nhms,100)

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
   call MAPL_TimerOn(MAPL,"-RUN")

! Get my internal private state. This contains the transforms
!  between the background grid and the ANA grid, as well as ANA grid
!  itself.
!-------------------------------------------------------------

    call ESMF_UserCompGetInternalState(gc, 'MKIAU_state', wrap, status)
    VERIFY_(STATUS)
    mkiau_internal_state => wrap%ptr

! Get parameters from generic background state
!---------------------------------------------

    call MAPL_Get( MAPL, IM=IMbkg, JM=JMbkg, LM=LMbkg, &
                   INTERNAL_ESMF_STATE=INTERNAL,       &
                                       RC=STATUS       )
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GC, grid=GRIDbkg, rc=status)
    VERIFY_(STATUS)

    call MAPL_GridGet(GRIDbkg, globalCellCountPerDim=DIMS, RC=STATUS)
    VERIFY_(STATUS)
    IMbkg_World=DIMS(1)
    JMbkg_World=DIMS(2)
    L_CUBE = JMbkg_World==6*IMbkg_World

! Get Resource Parameters
!------------------------
    call MAPL_GetResource(MAPL, GRIDINC,     LABEL="REPLAY_GRIDINC:", default="ANA",  RC=STATUS)
    VERIFY_(STATUS)
    GRIDINC = ESMF_UtilStringUpperCase(GRIDINC)
    _ASSERT( trim(GRIDINC) == "ANA" .or. trim(GRIDINC) == "BKG" ,'needs informative message')

    call MAPL_GetResource(MAPL, REPLAY_TIME_INTERP, LABEL="REPLAY_TIME_INTERP:", default="LINEAR",  RC=STATUS)
    VERIFY_(STATUS)
    REPLAY_TIME_INTERP = ESMF_UtilStringUpperCase(REPLAY_TIME_INTERP)
    _ASSERT( trim(REPLAY_TIME_INTERP) == "LINEAR" .or. trim(REPLAY_TIME_INTERP) == "CUBIC" ,'needs informative message')

! Check for 09 Files
! ------------------
    call ESMF_ClockGetAlarm(Clock,'RegularReplay09',Alarm,rc=Status)
    if(STATUS==ESMF_SUCCESS) then
       is_RegularReplay09_ringing = ESMF_AlarmIsRinging( Alarm,rc=status )
       VERIFY_(status)
       if( is_RegularReplay09_ringing ) then
           call MAPL_GetResource(MAPL, FILETMPL, LABEL="REPLAY_FILE09:", default="NULL", RC=STATUS)
           VERIFY_(STATUS)
           if( trim(FILETMPL) == "NULL" ) then
           call MAPL_GetResource(MAPL, FILETMPL, LABEL="REPLAY_FILE:", RC=STATUS)
           VERIFY_(STATUS)
           endif
       else
           call MAPL_GetResource(MAPL, FILETMPL, LABEL="REPLAY_FILE:", RC=STATUS)
           VERIFY_(STATUS)
       endif
    else
           call MAPL_GetResource(MAPL, FILETMPL, LABEL="REPLAY_FILE:", RC=STATUS)
           VERIFY_(STATUS)
    endif

    call MAPL_GetResource(MAPL, REPLAY_MODE, LABEL="REPLAY_MODE:",    default="NULL", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CREMAP,      LABEL="REPLAY_REMAP:",   default="yes",  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, FIXWIND,     LABEL="REPLAY_WINDFIX:", default="yes",  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, K,           Label="ANALYZE_TS:",     default=0,      RC=STATUS)
    VERIFY_(STATUS)

        ANALYZE_TS = (K /= 0)
    if( ANALYZE_TS ) then
         REPLAY_TS = 'YES'
    else
         REPLAY_TS = 'NO'
    endif

    call MAPL_GetResource(MAPL, IHAVEAINC,     Label='REPLAY_TO_ANAINC:', default=0,  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_PHIS,   Label="REPLAY_PHIS:",   default='YES',           RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_TS,     Label="REPLAY_TS:",     default=trim(REPLAY_TS), RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_P ,     Label="REPLAY_P:",      default='YES',           RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_PS,     Label="REPLAY_PS:",     default='YES',           RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_DP,     Label="REPLAY_DP:",     default='YES',           RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_U ,     Label="REPLAY_U:",      default='YES',           RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_V ,     Label="REPLAY_V:",      default='YES',           RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_T ,     Label="REPLAY_T:",      default='YES',           RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_QV,     Label="REPLAY_QV:",     default='YES',           RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_O3,     Label="REPLAY_O3:",     default='YES',           RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_T_TYPE, Label="REPLAY_T_TYPE:", default='NULL',          RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, REPLAY_P_FACTOR,  Label="REPLAY_P_FACTOR:",  default=1.0 , RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_U_FACTOR,  Label="REPLAY_U_FACTOR:",  default=1.0 , RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_V_FACTOR,  Label="REPLAY_V_FACTOR:",  default=1.0 , RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_T_FACTOR,  Label="REPLAY_T_FACTOR:",  default=1.0 , RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_QV_FACTOR, Label="REPLAY_QV_FACTOR:", default=1.0 , RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_O3_FACTOR, Label="REPLAY_O3_FACTOR:", default=1.0 , RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, REPLAY_TS_FACTOR, Label="REPLAY_TS_FACTOR:", default=1.0 , RC=STATUS) ; VERIFY_(STATUS)

      REPLAY_PHIS   = ESMF_UtilStringUpperCase(REPLAY_PHIS)
      REPLAY_TS     = ESMF_UtilStringUpperCase(REPLAY_TS)
      REPLAY_P      = ESMF_UtilStringUpperCase(REPLAY_P )
      REPLAY_U      = ESMF_UtilStringUpperCase(REPLAY_U )
      REPLAY_V      = ESMF_UtilStringUpperCase(REPLAY_V )
      REPLAY_T      = ESMF_UtilStringUpperCase(REPLAY_T )
      REPLAY_QV     = ESMF_UtilStringUpperCase(REPLAY_QV)
      REPLAY_O3     = ESMF_UtilStringUpperCase(REPLAY_O3)
      REPLAY_T_TYPE = ESMF_UtilStringUpperCase(REPLAY_T_TYPE)

    L_REPLAY_PHIS   = trim(REPLAY_PHIS)  .ne.'NO' .and. trim(REPLAY_PHIS)  .ne.'NULL'
    L_REPLAY_TS     = trim(REPLAY_TS)    .ne.'NO' .and. trim(REPLAY_TS)    .ne.'NULL'
    L_REPLAY_P      = trim(REPLAY_P)     .ne.'NO' .and. trim(REPLAY_P)     .ne.'NULL'
    L_REPLAY_U      = trim(REPLAY_U)     .ne.'NO' .and. trim(REPLAY_U)     .ne.'NULL'
    L_REPLAY_V      = trim(REPLAY_V)     .ne.'NO' .and. trim(REPLAY_V)     .ne.'NULL'
    L_REPLAY_T      = trim(REPLAY_T)     .ne.'NO' .and. trim(REPLAY_T)     .ne.'NULL'
    L_REPLAY_QV     = trim(REPLAY_QV)    .ne.'NO' .and. trim(REPLAY_QV)    .ne.'NULL'
    L_REPLAY_O3     = trim(REPLAY_O3)    .ne.'NO' .and. trim(REPLAY_O3)    .ne.'NULL'

    call MAPL_GetResource(MAPL, DAMPBEG,  LABEL="REPLAY_DAMPBEG:", default=1.0, RC=status)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, DAMPEND,  LABEL="REPLAY_DAMPEND:", default=1.0, RC=status)
    VERIFY_(STATUS)
    _ASSERT(DAMPBEG.le.DAMPEND   ,'needs informative message')

    call MAPL_GetResource(MAPL, BLEND_AT_PBL,  LABEL="REPLAY_BLEND_AT_PBL:", default=.FALSE., RC=status)
    VERIFY_(STATUS)

       CREMAP = ESMF_UtilStringUpperCase(CREMAP)
      FIXWIND = ESMF_UtilStringUpperCase(FIXWIND)
    DOWINDFIX = trim(FIXWIND)=="YES"

     CALL MAPL_GetResource(MAPL,JCAP,LABEL="MKIAU_JCAP:",default=-1,RC=STATUS)
     VERIFY_(STATUS)
     USE_SPECFILT = (JCAP /= -1)

    if( DOWINDFIX .or. USE_SPECFILT ) then
        _ASSERT( trim(GRIDINC) == "ANA" ,'needs informative message')
    endif

! **********************************************************************
! ****      Check REPLAY_FILE Resolution and Create ANA Grid        ****
! **********************************************************************

    call ESMF_ClockGet(clock, currTime=currTime, calendar=cal, rc=status)
    VERIFY_(STATUS)
    call ESMF_TimeGet(currTIME, timeString=DATE, RC=STATUS)
    VERIFY_(STATUS)
    call strToInt(DATE, nymd, nhms)

    call MAPL_GetResource(MAPL,FileFreq_SEC, Label="REPLAY_FILE_FREQUENCY:",      default=-999,   rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL,FileReft_HMS, Label="REPLAY_FILE_REFERENCE_TIME:", default=000000, rc=STATUS )
    VERIFY_(STATUS)

    if( FileFreq_SEC == -999 ) then
        REPLAY_TIME  = currTime
    else
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
    endif

! --------------------------------------------------------------------------------------------------------
    if( currTime == REPLAY_TIME ) then

        REPLAY_TIMEP0 = REPLAY_TIME

        call ESMF_CFIOstrTemplate ( REPLAY_FILEP0, FILETMPL, 'GRADS', nymd=nymd, nhms=nhms, stat=STATUS )
        VERIFY_(STATUS)

        if(MAPL_AM_I_ROOT() ) then
           print *, 'Current nymd: ',nymd,'  nhms: ',nhms,'  FAC:  1.00000'
        endif

    else

        if( currTime  < REPLAY_TIME ) then
        REPLAY_TIMEP1 = REPLAY_TIME + FileFreq
        REPLAY_TIMEP0 = REPLAY_TIME
        REPLAY_TIMEM1 = REPLAY_TIME - FileFreq
        REPLAY_TIMEM2 = REPLAY_TIME - FileFreq*2

        else
        REPLAY_TIMEP1 = REPLAY_TIME + FileFreq*2
        REPLAY_TIMEP0 = REPLAY_TIME + FileFreq
        REPLAY_TIMEM1 = REPLAY_TIME
        REPLAY_TIMEM2 = REPLAY_TIME - FileFreq
        endif

        call ESMF_TimeGet(REPLAY_TIMEP1, timeString=DATE, RC=STATUS)
        VERIFY_(STATUS)
        call strToInt(DATE, nymdp1, nhmsp1)
        call ESMF_CFIOstrTemplate ( REPLAY_FILEP1, FILETMPL, 'GRADS', nymd=nymdp1, nhms=nhmsp1, stat=STATUS )
        VERIFY_(STATUS)

        call ESMF_TimeGet(REPLAY_TIMEP0, timeString=DATE, RC=STATUS)
        VERIFY_(STATUS)
        call strToInt(DATE, nymdp0, nhmsp0)
        call ESMF_CFIOstrTemplate ( REPLAY_FILEP0, FILETMPL, 'GRADS', nymd=nymdp0, nhms=nhmsp0, stat=STATUS )
        VERIFY_(STATUS)

        call ESMF_TimeGet(REPLAY_TIMEM1, timeString=DATE, RC=STATUS)
        VERIFY_(STATUS)
        call strToInt(DATE, nymdm1, nhmsm1)
        call ESMF_CFIOstrTemplate ( REPLAY_FILEM1, FILETMPL, 'GRADS', nymd=nymdm1, nhms=nhmsm1, stat=STATUS )
        VERIFY_(STATUS)

        call ESMF_TimeGet(REPLAY_TIMEM2, timeString=DATE, RC=STATUS)
        VERIFY_(STATUS)
        call strToInt(DATE, nymdm2, nhmsm2)
        call ESMF_CFIOstrTemplate ( REPLAY_FILEM2, FILETMPL, 'GRADS', nymd=nymdm2, nhms=nhmsm2, stat=STATUS )
        VERIFY_(STATUS)

        if( REPLAY_TIME_INTERP == "CUBIC" ) then
        facp1 = ( (currTime-REPLAY_TIMEP0) / (REPLAY_TIMEP1-REPLAY_TIMEP0) ) &
              * ( (currTime-REPLAY_TIMEM1) / (REPLAY_TIMEP1-REPLAY_TIMEM1) ) &
              * ( (currTime-REPLAY_TIMEM2) / (REPLAY_TIMEP1-REPLAY_TIMEM2) )

        facp0 = ( (REPLAY_TIMEP1-currTime) / (REPLAY_TIMEP1-REPLAY_TIMEP0) ) &
              * ( (currTime-REPLAY_TIMEM1) / (REPLAY_TIMEP0-REPLAY_TIMEM1) ) &
              * ( (currTime-REPLAY_TIMEM2) / (REPLAY_TIMEP0-REPLAY_TIMEM2) )

        facm1 = ( (REPLAY_TIMEP1-currTime) / (REPLAY_TIMEP1-REPLAY_TIMEM1) ) &
              * ( (REPLAY_TIMEP0-currTime) / (REPLAY_TIMEP0-REPLAY_TIMEM1) ) &
              * ( (currTime-REPLAY_TIMEM2) / (REPLAY_TIMEM1-REPLAY_TIMEM2) )

        facm2 = ( (REPLAY_TIMEP1-currTime) / (REPLAY_TIMEP1-REPLAY_TIMEM2) ) &
              * ( (REPLAY_TIMEP0-currTime) / (REPLAY_TIMEP0-REPLAY_TIMEM2) ) &
              * ( (REPLAY_TIMEM1-currTime) / (REPLAY_TIMEM1-REPLAY_TIMEM2) )

        if(MAPL_AM_I_ROOT() ) then
           write(6,'(1x,a,i8.8,a,i6.6,a,e17.10)') 'Previous nymd: ',nymdm2,'  nhms: ',nhmsm2,'  FACM2: ',facm2
           write(6,'(1x,a,i8.8,a,i6.6,a,e17.10)') 'Previous nymd: ',nymdm1,'  nhms: ',nhmsm1,'  FACM1: ',facm1
           write(6,'(1x,a,i8.8,a,i6.6         )') 'Current  nymd: ',nymd , '  nhms: ',nhms
           write(6,'(1x,a,i8.8,a,i6.6,a,e17.10)') 'Next     nymd: ',nymdp0,'  nhms: ',nhmsp0,'  FACP0: ',facp0
           write(6,'(1x,a,i8.8,a,i6.6,a,e17.10)') 'Next     nymd: ',nymdp1,'  nhms: ',nhmsp1,'  FACP1: ',facp1
           print *
        endif
        else ! LINEAR Case
            facp0 = (currTime-REPLAY_TIMEM1) / (REPLAY_TIMEP0-REPLAY_TIMEM1)

            facm1 = (REPLAY_TIMEP0-currTime) / (REPLAY_TIMEP0-REPLAY_TIMEM1)

            if(MAPL_AM_I_ROOT() ) then
               write(6,'(1x,a,i8.8,a,i6.6,a,e17.10)') 'Previous nymd: ',nymdm1,'  nhms: ',nhmsm1,'  FACM1: ',facm1
               write(6,'(1x,a,i8.8,a,i6.6         )') 'Current  nymd: ',nymd , '  nhms: ',nhms
               write(6,'(1x,a,i8.8,a,i6.6,a,e17.10)') 'Next     nymd: ',nymdp0,'  nhms: ',nhmsp0,'  FACP0: ',facp0
               print *
            endif
        endif

    endif
! --------------------------------------------------------------------------------------------------------

    call CFIO_Open       ( REPLAY_FILEP0, 1, fid, STATUS )
    VERIFY_(STATUS)
    call CFIO_DimInquire ( fid, IMana_World, JMana_world, LMana, nt, nvars, natts, rc=STATUS )
    VERIFY_(STATUS)
    call CFIO_Close      ( fid, STATUS )
    VERIFY_(STATUS)

    call MAPL_MakeDecomposition(nx,ny,rc=status)
    VERIFY_(status)

    do_transforms = ( IMbkg_World /= IMana_World ) .or. &
                    ( JMbkg_World /= JMana_World ) .or. &
                    ( LMbkg       /= LMana       )

    refresh_internal_state = .false. ! Default
    if (first) then
       refresh_internal_state = .true. ! first time: do it!
    else
       if ( mkiau_internal_state%IM /= IMana_World .or. &
            mkiau_internal_state%JM /= JMana_World .or. &
            mkiau_internal_state%LM /= LMana ) then
          refresh_internal_state = .true. ! Resolution of Analysis File has changed since last update
       end if
    end if

    if (refresh_internal_state) then
       if (.not. first) then
          call WRITE_PARALLEL("Destroying GRIDana...")
          call ESMF_GridDestroy(mkiau_internal_state%GRIDana, rc=status)
          VERIFY_(STATUS)
          call ESMF_GridDestroy(mkiau_internal_state%GRIDrep, rc=status)
          VERIFY_(STATUS)
       end if

       call WRITE_PARALLEL("Creating GRIDana...")
       write(imstr,*) IMana_World
       write(jmstr,*) JMana_World
       gridAnaName='PC'//trim(adjustl(imstr))//'x'//trim(adjustl(jmstr))//'-DC'

       ! Get grid_dimensions from file.
       call CFIO_Open(REPLAY_FILEP0, 1, fid, rc=status)
       VERIFY_(status)
       call CFIO_DimInquire (fid, IMana_World, JMana_World, LMana, nt, nvars, natts, rc=status)
       VERIFY_(status)
       call CFIO_Close(fid, rc=status)
       VERIFY_(status)

       block
         use MAPL_LatLonGridFactoryMod
         GRIDrep = grid_manager%make_grid(                                                 &
                   LatLonGridFactory(im_world=IMana_World, jm_world=JMana_World, lm=LMana, &
                   nx=NX, ny=NY, pole='PC', dateline= 'DC', rc=status)                     )
         VERIFY_(STATUS)
         GRIDana = grid_manager%make_grid(                                                 &
                   LatLonGridFactory(im_world=IMana_World, jm_world=JMana_World, lm=LMbkg, &
                   nx=NX, ny=NY, pole='PC', dateline= 'DC', rc=status)                     )
         VERIFY_(STATUS)
       end block

       mkiau_internal_state%im      =   IMana_World
       mkiau_internal_state%jm      =   JMana_World
       mkiau_internal_state%lm      =   LMana
       mkiau_internal_state%GRIDana = GRIDana
       mkiau_internal_state%GRIDrep = GRIDrep

       call MAPL_GetResource(MAPL, K, Label="BKG2ANACNSRV:", default=0, RC=STATUS)
       VERIFY_(STATUS)
       BKG2ANAConsrv = (K /= 0)

       call MAPL_GetResource(MAPL, K, Label="ANA2BKGCNSRV:", default=0, RC=STATUS)
       VERIFY_(STATUS)
       ANA2BKGConsrv = (K /= 0)

       if (ana2bkgconsrv) then
          mkiau_internal_state%ana2bkg_regridder => new_regridder_manager%make_regridder(GRIDana, GRIDbkg, REGRID_METHOD_CONSERVE, rc=status)
          VERIFY_(status)
       else
          mkiau_internal_state%ana2bkg_regridder => new_regridder_manager%make_regridder(GRIDana, GRIDbkg, REGRID_METHOD_BILINEAR, rc=status)
          VERIFY_(status)
       end if

       if (bkg2anaConsrv) then
          mkiau_internal_state%bkg2ana_regridder => new_regridder_manager%make_regridder(GRIDbkg, GRIDana, REGRID_METHOD_CONSERVE, rc=status)
          VERIFY_(status)
       else
          mkiau_internal_state%bkg2ana_regridder => new_regridder_manager%make_regridder(GRIDbkg, GRIDana, REGRID_METHOD_BILINEAR, rc=status)
          VERIFY_(status)
       end if

    else
       if(first) call WRITE_PARALLEL("Using stored GRIDana...")
       GRIDana = mkiau_internal_state%GRIDana
       GRIDrep = mkiau_internal_state%GRIDrep
    end if

    !ALT: Get current VM and the mpi communicator
    !--------------------------------------------
    call ESMF_VMGetCurrent(vm, rc=status)
    VERIFY_(STATUS)
    call ESMF_VmGet(VM, mpicommunicator=vm_comm, rc=status)
    VERIFY_(STATUS)

    ANA2BKG => mkiau_internal_state%ANA2BKG_regridder
    BKG2ANA => mkiau_internal_state%BKG2ANA_regridder

!   Set Local Dimensions to GRIDana and GRIDbkg
!   -------------------------------------------
    call MAPL_GridGet(GRIDrep, localCellCountPerDim=DIMS, RC=STATUS)
    VERIFY_(STATUS)
    IMana   = DIMS(1)
    JMana   = DIMS(2)
    LMana   = DIMS(3)

    call MAPL_GridGet(GRIDbkg, localCellCountPerDim=DIMS, RC=STATUS)
    VERIFY_(STATUS)
    IMbkg   = DIMS(1)
    JMbkg   = DIMS(2)
    LMbkg   = DIMS(3)

!   Set Local Dimensions to GRIDINC (i.e., the GRID on which the increments are computed)
!   Note:  In all cases, the vertical resolution is defined by the Background
!   -------------------------------------------------------------------------------------
    if( trim(GRIDINC) == "ANA" ) then
        IM   = IMana
        JM   = JMana
    endif
    if( trim(GRIDINC) == "BKG" ) then
        IM   = IMbkg
        JM   = JMbkg
    endif
        LM   = LMbkg
        LMP1 = LMbkg+1

    if ( IHAVEAINC/=0 ) then
       call handleINC_
    else
       call handleANA_
    endif

    call MAPL_TimerOff(MAPL,"-RUN")
    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

CONTAINS

! *****************************************************************************
!   This interface incorporates pre-calculated Analysis Increments

    subroutine handleINC_

! *****************************************************************************

    real,pointer :: aptr2d(:,:), aptr3d(:,:,:)  ! analysis increment pointers
    real,pointer :: gptr2d(:,:), gptr3d(:,:,:)  ! gcm background pointers

    character(len=*), parameter :: incnames(7) = (/ 'sphu ', &
                                                    'u    ', &
                                                    'v    ', &
                                                    'tv   ', &
                                                    'ozone', &
                                                    'delp ', &
                                                    'ts   ' /)
    integer rank,ni
    type(ESMF_Field)  :: Field
    real,allocatable, dimension(:,:,:)  :: dp
    real,allocatable, dimension(:,:,:)  :: uptr
    real,allocatable, dimension(:,:,:)  :: vptr
    real,allocatable, dimension(:,:,:)  :: pke
    real,allocatable, dimension(:,:,:)  :: pkz
    real,allocatable, dimension(:,:,:)  :: dpkz
    character(len=ESMF_MAXSTR) :: name

! *****************************************************************************
! ****   READ Internal STATE (ie. ANA.ETA) from REPLAY File into BUNDLE    ****
! *****************************************************************************

    RBUNDLEP0 = ESMF_FieldBundleCreate( RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldBundleSet(RBUNDLEP0, grid=GRIDana, rc=status)
    VERIFY_(STATUS)
    call MAPL_CFIORead ( REPLAY_FILEP0, REPLAY_TIMEP0, RBUNDLEP0, RC=status)
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet ( RBUNDLEP0, fieldCount=NQ, RC=STATUS )
    VERIFY_(STATUS)

!   Get pointers to hold IAU increment
!   ----------------------------------
    call MAPL_GetPointer(export,   du, 'DUDT',  alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,   dv, 'DVDT',  alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,   dt, 'DTDT',  alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,   dq, 'DQVDT', alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,  do3, 'DO3DT', alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export, dple, 'DPEDT', alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,  dts, 'DTSDT', alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)

!   Convert increment fields to background grid
!   -------------------------------------------


! Loop over GSI increment fields
! ------------------------------
    allocate(dp(IMbkg,JMbkg,LMbkg))
    allocate(gptr3d(IMbkg,JMbkg,LMbkg))
    allocate(gptr2d(IMbkg,JMbkg))
    allocate(qdum1(IMbkg,JMbkg,1))
    allocate(qdum2(IM,   JM,   1))
    do ni = 1, nq
       call ESMF_FieldBundleGet(RBUNDLEP0, ni, Field, __RC__ )
       call ESMF_FieldGet(Field, NAME=NAME, dimCount = rank, __RC__ )
       if (.not.check_list_(NAME,incnames)) cycle
       if (rank==2) then
           call ESMF_FieldGet(Field, farrayPtr=aptr2d, __RC__ )
           if (do_transforms) then
               qdum2(:,:,1)=aptr2d
               call mkiau_internal_state%ana2bkg_regridder%regrid(qdum2, qdum1, rc=status)
               VERIFY_(status)
               gptr2d=qdum1(:,:,1)
           else
               gptr2d=aptr2d
           endif
           if(trim(NAME)=='ts') then
              if( L_REPLAY_TS ) then
                 dts = gptr2d
               else
                 dts = 0.0
               endif
           endif
       else
           call ESMF_FieldGet(Field, farrayPtr=aptr3d, __RC__ )
           if(trim(NAME)=='u'.or.trim(NAME)=='v') then
              if(trim(NAME)=='u') then
                 allocate(uptr(IM,JM,LM))
                 uptr=aptr3d
              endif
              if(trim(NAME)=='v') then
                 allocate(vptr(IM,JM,LM))
                 vptr=aptr3d
              endif
           else
           if (do_transforms) then
               call mkiau_internal_state%ana2bkg_regridder%regrid(aptr3d, gptr3d, rc=status)
               VERIFY_(status)
           else
               gptr3d=aptr3d
           endif
           if(trim(NAME)=='ozone') do3=gptr3d
           if(trim(NAME)=='sphu' ) dq =gptr3d
           if(trim(NAME)=='tv'   ) dt =gptr3d
           if(trim(NAME)=='delp' ) dp =gptr3d
           endif
       endif
    enddo
    deallocate(qdum2)
    deallocate(qdum1)

!   U and V
!   -------
    ! could apply wind fix here ... but conversion to cubed
    ! TBD
    if (do_transforms) then
       if( L_CUBE ) then
          call mkiau_internal_state%ana2bkg_regridder%regrid(uptr, vptr, du, dv, rotate=.false., rc=status)
          VERIFY_(status)
       else
          call mkiau_internal_state%ana2bkg_regridder%regrid(uptr, du, rc=status)
          VERIFY_(status)
          call mkiau_internal_state%ana2bkg_regridder%regrid(vptr, dv, rc=status)
          VERIFY_(status)
          call POLEFIX ( du,dv,VM,GRIDbkg )
       endif
    else
        du=uptr
        dv=vptr
    endif

!   Calculate 3d-pressure change
!   -----------------------------
    dple(:,:,0) = 0.0
    do L=1,LMbkg
       dple(:,:,L) = dp(:,:,L)
    enddo

    call MAPL_GetPointer( import, tv_bkg, 'TV',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( import, q_bkg,'QV',   RC=STATUS)
    VERIFY_(STATUS)

!   Convert virtual temperature increment into dry temperature increment
!   -------------------------------------------------------------------
    dt = ( dt - eps*dq*(tv_bkg/(1.0+eps*q_bkg)) ) / (1.0+eps*q_bkg) ! now dt is inc on dry temperature

!   Clean up
!   --------
    call MAPL_FieldBundleDestroy ( RBUNDLEP0, RC=STATUS)
    VERIFY_(STATUS)
    deallocate(dp)
    deallocate(gptr2d)
    deallocate(gptr3d)
    if (allocated(uptr)) deallocate(uptr)
    if (allocated(vptr)) deallocate(vptr)

    end subroutine handleINC_

! *****************************************************************************
!   This interface computes Analysis Increments based on BKG and ANA variables

    subroutine handleANA_

! *****************************************************************************

    allocate( phis_bkg(IM,JM)     , source = 0.0 )
    allocate(   ts_bkg(IM,JM)     , source = 0.0 )
    allocate(   ps_bkg(IM,JM)     , source = 0.0 )
    allocate(    u_bkg(IM,JM,1:LM), source = 0.0 )
    allocate(    v_bkg(IM,JM,1:LM), source = 0.0 )
    allocate(    t_bkg(IM,JM,1:LM), source = 0.0 )
    allocate(   tv_bkg(IM,JM,1:LM), source = 0.0 )
    allocate(    q_bkg(IM,JM,1:LM), source = 0.0 )
    allocate(   o3_bkg(IM,JM,1:LM), source = 0.0 )
    allocate(   dp_bkg(IM,JM,1:LM), source = 0.0 )
    allocate(  ple_bkg(IM,JM,0:LM), source = 0.0 )

! **********************************************************************
! ****      Transform Import Data (ie. BKG.ETA) to ANA Grid         ****
! **********************************************************************

    call MAPL_GetPointer( import, uptr3d, 'U',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( export, temp3d, 'UBKG', RC=STATUS)
    VERIFY_(STATUS)
    if( associated(temp3d) ) temp3d = uptr3d

    call MAPL_GetPointer( import, vptr3d, 'V',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( export, temp3d, 'VBKG', RC=STATUS)
    VERIFY_(STATUS)
    if( associated(temp3d) ) temp3d = vptr3d

    if ( trim(GRIDINC)=="ANA" .and. do_transforms) then
       if (L_CUBE) then
          call mkiau_internal_state%bkg2ana_regridder%regrid(uptr3d, vptr3d, u_bkg, v_bkg, rotate=.false., rc=status)
          VERIFY_(status)
       else
          call mkiau_internal_state%bkg2ana_regridder%regrid(uptr3d, u_bkg, rc=status)
          VERIFY_(status)
          call mkiau_internal_state%bkg2ana_regridder%regrid(vptr3d, v_bkg, rc=status)
          VERIFY_(status)
          call POLEFIX ( u_bkg,v_bkg,VM,GRIDana )
      endif
   else
       u_bkg=uptr3d
       v_bkg=vptr3d
   endif

    call MAPL_GetPointer( import, ptr3d,  'TV',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( export, temp3d, 'TVBKG', RC=STATUS)
    VERIFY_(STATUS)
    if( associated(temp3d) ) temp3d = ptr3d
    if (trim(GRIDINC)=="ANA" .and. do_transforms) then
       call mkiau_internal_state%bkg2ana_regridder%regrid(ptr3d, tv_bkg, rc=status)
       VERIFY_(status)
    else
       tv_bkg=ptr3d
    endif

    call MAPL_GetPointer( import, ptr3d,  'DELP',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( export, temp3d, 'DELPBKG', RC=STATUS)
    VERIFY_(STATUS)
    if( associated(temp3d) ) temp3d = ptr3d
    if (trim(GRIDINC)=="ANA" .and. do_transforms) then
       call mkiau_internal_state%bkg2ana_regridder%regrid(ptr3d, dp_bkg, rc=status)
       VERIFY_(status)
    else
       dp_bkg=ptr3d
    endif

    call MAPL_GetPointer( import, ptr3d,  'O3PPMV',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( export, temp3d, 'O3PPMVBKG', RC=STATUS)
    VERIFY_(STATUS)
    if( associated(temp3d) ) temp3d = ptr3d
    if (trim(GRIDINC)=="ANA" .and. do_transforms) then
       call mkiau_internal_state%bkg2ana_regridder%regrid(ptr3d, o3_bkg, rc=status)
       VERIFY_(status)
    else
        o3_bkg=ptr3d
    endif

    call MAPL_GetPointer(import,  ptr3d,   'QV',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( export, temp3d, 'QVBKG', RC=STATUS)
    VERIFY_(STATUS)
    if( associated(temp3d) ) temp3d = ptr3d

    if (trim(GRIDINC)=="ANA" .and. do_transforms) then
       call mkiau_internal_state%bkg2ana_regridder%regrid(ptr3d, q_bkg, rc=status)
       VERIFY_(status)
    else
       q_bkg=ptr3d
    endif

    allocate ( qdum1(IMbkg,JMbkg,3) )
    allocate ( qdum2(IM,   JM,   3) )
    qdum1=0.0

    call MAPL_GetPointer( import, ptr2d, 'PHIS', RC=STATUS)
    VERIFY_(STATUS)
    qdum1(:,:,1) = ptr2d
    call MAPL_GetPointer( import, ptr2d,  'TS',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( export, temp2d, 'TSBKG', RC=STATUS)
    VERIFY_(STATUS)
    if( associated(temp2d) ) temp2d = ptr2d
    qdum1(:,:,2) = ptr2d

    call MAPL_GetPointer( import, ptr2d,  'PS',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer( export, temp2d, 'PSBKG', RC=STATUS)
    VERIFY_(STATUS)
    if( associated(temp2d) ) temp2d = ptr2d
    qdum1(:,:,3) = ptr2d

    if (trim(GRIDINC)=="ANA" .and. do_transforms) then
       call mkiau_internal_state%bkg2ana_regridder%regrid(qdum1, qdum2, rc=status)
       VERIFY_(status)
    else
       qdum2=qdum1
    endif
    phis_bkg = qdum2(:,:,1)
      ts_bkg = qdum2(:,:,2)
      ps_bkg = qdum2(:,:,3)

    deallocate ( qdum1 )
    deallocate ( qdum2 )

    call MAPL_GetPointer(import, ak, 'AK', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import, bk, 'BK', RC=STATUS)
    VERIFY_(STATUS)

! BKG Pressure Variables
! ----------------------
    ple_bkg(:,:,0) = ak(0)
    do L=1,lm
    ple_bkg(:,:,L) = ple_bkg(:,:,L-1) + dp_bkg(:,:,L)
    end do

! BKG Dry Temperature
! -------------------
    t_bkg =  tv_bkg/(1.0+eps*q_bkg)

! *****************************************************************************
! ****   READ Internal STATE (ie. ANA.ETA) from REPLAY File into BUNDLE    ****
! *****************************************************************************

    if( NEED_BUNDLEP0 ) then
        RBUNDLEP0 = ESMF_FieldBundleCreate( RC=STATUS)
        VERIFY_(STATUS)
        if ( trim(GRIDINC)=="ANA" ) call ESMF_FieldBundleSet(RBUNDLEP0, grid=GRIDrep, rc=status)
        if ( trim(GRIDINC)=="BKG" ) call ESMF_FieldBundleSet(RBUNDLEP0, grid=GRIDbkg, rc=status)
        VERIFY_(STATUS)
        call MAPL_CFIORead ( REPLAY_FILEP0, REPLAY_TIMEP0, RBUNDLEP0 , RC=status)
        VERIFY_(STATUS)
             FILEP0 = REPLAY_FILEP0
        FILE_TIMEP0 = REPLAY_TIMEP0
        NEED_BUNDLEP0 = .FALSE.
    else if( (FILE_TIMEP0 .ne. REPLAY_TIMEP0) .or. (FILEP0 .ne. REPLAY_FILEP0) ) then
        call MAPL_CFIORead ( REPLAY_FILEP0, REPLAY_TIMEP0, RBUNDLEP0 , RC=status)
        VERIFY_(STATUS)
             FILEP0 = REPLAY_FILEP0
        FILE_TIMEP0 = REPLAY_TIMEP0
    endif

    if( currTime /= REPLAY_TIMEP0 ) then
        if( NEED_BUNDLEM1 ) then
            RBUNDLEM1 = ESMF_FieldBundleCreate( RC=STATUS)
            VERIFY_(STATUS)
            if ( trim(GRIDINC)=="ANA" ) call ESMF_FieldBundleSet(RBUNDLEM1, grid=GRIDrep, rc=status)
            if ( trim(GRIDINC)=="BKG" ) call ESMF_FieldBundleSet(RBUNDLEM1, grid=GRIDbkg, rc=status)
            VERIFY_(STATUS)
            call MAPL_CFIORead ( REPLAY_FILEM1, REPLAY_TIMEM1, RBUNDLEM1 , RC=status)
            VERIFY_(STATUS)
                 FILEM1 = REPLAY_FILEM1
            FILE_TIMEM1 = REPLAY_TIMEM1
            NEED_BUNDLEM1 = .FALSE.
        else if ( (FILE_TIMEM1 .ne. REPLAY_TIMEM1) .or. (FILEM1 .ne. REPLAY_FILEM1) ) then
            call MAPL_CFIORead ( REPLAY_FILEM1, REPLAY_TIMEM1, RBUNDLEM1 , RC=status)
            VERIFY_(STATUS)
                 FILEM1 = REPLAY_FILEM1
            FILE_TIMEM1 = REPLAY_TIMEM1
        endif

        if( REPLAY_TIME_INTERP == "CUBIC" ) then
            if( NEED_BUNDLEP1 ) then
                RBUNDLEP1 = ESMF_FieldBundleCreate( RC=STATUS)
                VERIFY_(STATUS)
                if ( trim(GRIDINC)=="ANA" ) call ESMF_FieldBundleSet(RBUNDLEP1, grid=GRIDrep, rc=status)
                if ( trim(GRIDINC)=="BKG" ) call ESMF_FieldBundleSet(RBUNDLEP1, grid=GRIDbkg, rc=status)
                VERIFY_(STATUS)
                call MAPL_CFIORead ( REPLAY_FILEP1, REPLAY_TIMEP1, RBUNDLEP1 , RC=status)
                VERIFY_(STATUS)
                     FILEP1 = REPLAY_FILEP1
                FILE_TIMEP1 = REPLAY_TIMEP1
                NEED_BUNDLEP1 = .FALSE.
            else if ( FILE_TIMEP1 .ne. REPLAY_TIMEP1 .or. (FILEP1 .ne. REPLAY_FILEP1) ) then
                call MAPL_CFIORead ( REPLAY_FILEP1, REPLAY_TIMEP1, RBUNDLEP1 , RC=status)
                VERIFY_(STATUS)
                     FILEP1 = REPLAY_FILEP1
                FILE_TIMEP1 = REPLAY_TIMEP1
            endif

            if( NEED_BUNDLEM2 ) then
                RBUNDLEM2 = ESMF_FieldBundleCreate( RC=STATUS)
                VERIFY_(STATUS)
                if ( trim(GRIDINC)=="ANA" ) call ESMF_FieldBundleSet(RBUNDLEM2, grid=GRIDrep, rc=status)
                if ( trim(GRIDINC)=="BKG" ) call ESMF_FieldBundleSet(RBUNDLEM2, grid=GRIDbkg, rc=status)
                VERIFY_(STATUS)
                call MAPL_CFIORead ( REPLAY_FILEM2, REPLAY_TIMEM2, RBUNDLEM2 , RC=status)
                VERIFY_(STATUS)
                     FILEM2 = REPLAY_FILEM2
                FILE_TIMEM2 = REPLAY_TIMEM2
                NEED_BUNDLEM2 = .FALSE.
            else if ( FILE_TIMEM2 .ne. REPLAY_TIMEM2 .or. (FILEM2 .ne. REPLAY_FILEM2) ) then
                call MAPL_CFIORead ( REPLAY_FILEM2, REPLAY_TIMEM2, RBUNDLEM2 , RC=status)
                VERIFY_(STATUS)
                     FILEM2 = REPLAY_FILEM2
                FILE_TIMEM2 = REPLAY_TIMEM2
            endif
        endif
    endif

    call ESMF_FieldBundleGet ( RBUNDLEP0, fieldCount=NQ, RC=STATUS )
    VERIFY_(STATUS)

    if( .not.allocated( rnames ) ) then
         allocate( RNAMES(NQ),STAT=STATUS )
         VERIFY_(STATUS)
         call ESMF_FieldBundleGet ( RBUNDLEP0, fieldNameList=RNAMES, rc=STATUS )
         VERIFY_(STATUS)
         if( first ) then
             if(MAPL_AM_I_ROOT() ) then
             print *
             print *, 'REPLAY File Dimensions: ', IMana_World,JMana_World,LMana
             print *
             print *, 'REPLAY File Variables, NQ: ', nq
             print *, '--------------------------'
             do k=1,nq
             print *, k,')  ',trim(rnames(k))
             enddo
             print *
             print *, 'REPLAY Options:'
             print *, '---------------'
             print *, 'REPLAY_TS ....... ',trim(REPLAY_TS)
             print *, 'REPLAY_P ........ ',trim(REPLAY_P)
             print *, 'REPLAY_U ........ ',trim(REPLAY_U)
             print *, 'REPLAY_V ........ ',trim(REPLAY_V)
             print *, 'REPLAY_T ........ ',trim(REPLAY_T)
             print *, 'REPLAY_QV ....... ',trim(REPLAY_QV)
             print *, 'REPLAY_O3 ....... ',trim(REPLAY_O3)
             print *
             endif
          endif
    endif

! **********************************************************************
! ****     Get Pointers to Internal STATE (ANA.ETA) from BUNDLE     ****
! **********************************************************************

    allocate ( phis_ana(IM,JM)     , source = 0.0 )
    allocate (   ts_ana(IM,JM)     , source = 0.0 )
    allocate (   ps_ana(IM,JM)     , source = 0.0 )
    allocate (   du_fix(IM,JM,  LM), source = 0.0 )
    allocate (   dv_fix(IM,JM,  LM), source = 0.0 )
    allocate (    u_ana(IM,JM,  LM), source = 0.0 )
    allocate (    v_ana(IM,JM,  LM), source = 0.0 )
    allocate (    t_ana(IM,JM,  LM), source = 0.0 )
    allocate (  thv_ana(IM,JM,  LM), source = 0.0 )
    allocate (    q_ana(IM,JM,  LM), source = 0.0 )
    allocate (   o3_ana(IM,JM,  LM), source = 0.0 )
    allocate (   pk_ana(IM,JM,  LM), source = 0.0 )
    allocate (  ple_ana(IM,JM,0:LM), source = 0.0 )
    allocate (  pke_ana(IM,JM,0:LM), source = 0.0 )

    allocate (   dp_rep(IM,JM,  LMana), source = 0.0 )
    allocate (    u_rep(IM,JM,  LMana), source = 0.0 )
    allocate (    v_rep(IM,JM,  LMana), source = 0.0 )
    allocate (    t_rep(IM,JM,  LMana), source = 0.0 )
    allocate (  thv_rep(IM,JM,  LMana), source = 0.0 )
    allocate (    q_rep(IM,JM,  LMana), source = 0.0 )
    allocate (   o3_rep(IM,JM,  LMana), source = 0.0 )
    allocate (   pk_rep(IM,JM,  LMana), source = 0.0 )
    allocate (  ple_rep(IM,JM,0:LMana), source = 0.0 )
    allocate (  pke_rep(IM,JM,0:LMana), source = 0.0 )

    allocate ( ak_rep(0:LMana), source = 0.0 )
    allocate ( bk_rep(0:LMana), source = 0.0 )

    doremap = (trim(cremap).eq.'YES') .or. (LMana.ne.LMbkg)

! Initialize ANA.ETA variables to Transformed BKG IMPORT State (In case REPLAY Variables are turned OFF)
! ------------------------------------------------------------------------------------------------------
    phis_ana = phis_bkg
      ts_ana =   ts_bkg
      ps_ana =   ps_bkg
       u_ana =    u_bkg
       v_ana =    v_bkg
       t_ana =    t_bkg
       q_ana =    q_bkg
      o3_ana =   o3_bkg

! Surface Temperature
!--------------------
        FOUND = .false.
        do k=1,nq
           if( match('ts',REPLAY_TS,rnames(k)) ) then
               call ESMFL_BundleGetPointertoData(RBUNDLEP0,trim(rnames(k)),ptr2d, RC=STATUS)
               if(STATUS==ESMF_SUCCESS) then
                          ts_ana =  ptr2d
                          if( currTime /= REPLAY_TIMEP0 ) then
                              call ESMFL_BundleGetPointertoData(RBUNDLEM1,trim(rnames(k)),ptr2d, RC=STATUS)
                              VERIFY_(STATUS)
                              ts_ana =  facp0*ts_ana + facm1*ptr2d
                              if( REPLAY_TIME_INTERP == "CUBIC" ) then
                                  call ESMFL_BundleGetPointertoData(RBUNDLEP1,trim(rnames(k)),ptr2d, RC=STATUS)
                                  VERIFY_(STATUS)
                                  ts_ana =  ts_ana + facp1*ptr2d
                              call ESMFL_BundleGetPointertoData(RBUNDLEM2,trim(rnames(k)),ptr2d, RC=STATUS)
                              VERIFY_(STATUS)
                              ts_ana =  ts_ana + facm2*ptr2d
                          endif
                          endif
                          if( REPLAY_TS_FACTOR.ne.1.0 ) ts_ana = ts_ana * REPLAY_TS_FACTOR
                          FOUND  = .true.
                          exit
               endif
           endif
        enddo
        if( .not.FOUND .and. L_REPLAY_TS ) then
           write(STRING,'(A)') "ANA Variable: TS Not Found!"
           call WRITE_PARALLEL( trim(STRING)   )
           RETURN_(ESMF_FAILURE)
        endif

! ANA Surface Geopotential Heights
!---------------------------------
        FOUND = .false.
        do k=1,nq
           if( match('phis',REPLAY_PHIS,rnames(k)) ) then
               call ESMFL_BundleGetPointertoData(RBUNDLEP0,trim(rnames(k)),ptr2d, RC=STATUS)
               if(STATUS==ESMF_SUCCESS) then
                          phis_ana =  ptr2d
                          FOUND    = .true.
                          exit
               endif
           endif
        enddo
        if( .not.FOUND .and. L_REPLAY_PHIS ) then
           write(STRING,'(A)') "ANA Variable: PHIS Not Found!"
           call WRITE_PARALLEL( trim(STRING)   )
           RETURN_(ESMF_FAILURE)
        endif

! ANA Surface Pressure
!---------------------
        FOUND = .false.
        do k=1,nq
           if( match('ps',REPLAY_PS,rnames(k)) ) then
               call ESMFL_BundleGetPointertoData(RBUNDLEP0,trim(rnames(k)),ptr2d, RC=STATUS)
               if(STATUS==ESMF_SUCCESS) then
                          ps_ana =  ptr2d
                          if( currTime /= REPLAY_TIMEP0 ) then
                              call ESMFL_BundleGetPointertoData(RBUNDLEM1,trim(rnames(k)),ptr2d, RC=STATUS)
                              VERIFY_(STATUS)
                              ps_ana =  facp0*ps_ana + facm1*ptr2d
                              if( REPLAY_TIME_INTERP == "CUBIC" ) then
                                  call ESMFL_BundleGetPointertoData(RBUNDLEP1,trim(rnames(k)),ptr2d, RC=STATUS)
                                  VERIFY_(STATUS)
                                  ps_ana =  ps_ana + facp1*ptr2d
                              call ESMFL_BundleGetPointertoData(RBUNDLEM2,trim(rnames(k)),ptr2d, RC=STATUS)
                              VERIFY_(STATUS)
                              ps_ana =  ps_ana + facm2*ptr2d
                          endif
                          endif
                          if( REPLAY_P_FACTOR.ne.1.0 ) ps_ana = ps_ana * REPLAY_P_FACTOR
                          FOUND  = .true.
                          exit
               endif
           endif
        enddo
        if( .not.FOUND .and. L_REPLAY_P ) then
           write(STRING,'(A)') "ANA Variable: PS Not Found!"
           call WRITE_PARALLEL( trim(STRING)   )
           RETURN_(ESMF_FAILURE)
        endif

! ANA Pressure Thickness
!-----------------------
        FOUND = .false.
        do k=1,nq
           if( match('dp',REPLAY_DP,rnames(k)) ) then
               call ESMFL_BundleGetPointertoData(RBUNDLEP0,trim(rnames(k)),ptr3d, RC=STATUS)
               if(STATUS==ESMF_SUCCESS) then
                          dp_rep =  ptr3d
                          if( currTime /= REPLAY_TIMEP0 ) then
                              call ESMFL_BundleGetPointertoData(RBUNDLEM1,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              dp_rep =  facp0*dp_rep + facm1*ptr3d
                              if( REPLAY_TIME_INTERP == "CUBIC" ) then
                                  call ESMFL_BundleGetPointertoData(RBUNDLEP1,trim(rnames(k)),ptr3d, RC=STATUS)
                                  VERIFY_(STATUS)
                                  dp_rep =  dp_rep + facp1*ptr3d
                              call ESMFL_BundleGetPointertoData(RBUNDLEM2,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              dp_rep =  dp_rep + facm2*ptr3d
                          endif
                          endif
                          if( REPLAY_P_FACTOR.ne.1.0 ) dp_rep = dp_rep * REPLAY_P_FACTOR
                          FOUND  = .true.
                          exit
               endif
           endif
        enddo
        if( .not.FOUND .and. L_REPLAY_P ) then
           write(STRING,'(A)') "ANA Variable: DP Not Found!"
           call WRITE_PARALLEL( trim(STRING)   )
           RETURN_(ESMF_FAILURE)
        endif

! ANA U-Wind
!-----------
        FOUND = .false.
        do k=1,nq
           if( match('u',REPLAY_U,rnames(k)) ) then
               call ESMFL_BundleGetPointertoData(RBUNDLEP0,trim(rnames(k)),ptr3d, RC=STATUS)
               if(STATUS==ESMF_SUCCESS) then
                          u_rep =  ptr3d
                          if( currTime /= REPLAY_TIMEP0 ) then
                              call ESMFL_BundleGetPointertoData(RBUNDLEM1,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              u_rep =  facp0*u_rep + facm1*ptr3d
                              if( REPLAY_TIME_INTERP == "CUBIC" ) then
                                  call ESMFL_BundleGetPointertoData(RBUNDLEP1,trim(rnames(k)),ptr3d, RC=STATUS)
                                  VERIFY_(STATUS)
                                  u_rep =  u_rep + facp1*ptr3d
                              call ESMFL_BundleGetPointertoData(RBUNDLEM2,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              u_rep =  u_rep + facm2*ptr3d
                          endif
                          endif
                          if( REPLAY_U_FACTOR.ne.1.0 ) u_rep = u_rep * REPLAY_U_FACTOR
                          FOUND = .true.
                          exit
               endif
           endif
        enddo
        if( .not.FOUND .and. L_REPLAY_U ) then
           write(STRING,'(A)') "ANA Variable: U Not Found!"
           call WRITE_PARALLEL( trim(STRING)   )
           RETURN_(ESMF_FAILURE)
        endif

! ANA V-Wind
!-----------
        FOUND = .false.
        do k=1,nq
           if( match('v',REPLAY_V,rnames(k)) ) then
               call ESMFL_BundleGetPointertoData(RBUNDLEP0,trim(rnames(k)),ptr3d, RC=STATUS)
               if(STATUS==ESMF_SUCCESS) then
                          v_rep =  ptr3d
                          if( currTime /= REPLAY_TIMEP0 ) then
                              call ESMFL_BundleGetPointertoData(RBUNDLEM1,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              v_rep =  facp0*v_rep + facm1*ptr3d
                              if( REPLAY_TIME_INTERP == "CUBIC" ) then
                                  call ESMFL_BundleGetPointertoData(RBUNDLEP1,trim(rnames(k)),ptr3d, RC=STATUS)
                                  VERIFY_(STATUS)
                                  v_rep =  v_rep + facp1*ptr3d
                              call ESMFL_BundleGetPointertoData(RBUNDLEM2,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              v_rep =  v_rep + facm2*ptr3d
                          endif
                          endif
                          if( REPLAY_V_FACTOR.ne.1.0 ) v_rep = v_rep * REPLAY_V_FACTOR
                          FOUND = .true.
                          exit
               endif
           endif
        enddo
        if( .not.FOUND .and. L_REPLAY_V ) then
           write(STRING,'(A)') "ANA Variable: V Not Found!"
           call WRITE_PARALLEL( trim(STRING)   )
           RETURN_(ESMF_FAILURE)
        endif

   ! the following assumes GRIDana is lat-lon
   ! ----------------------------------------
   ! call POLEFIX ( u_ana,v_ana,VM,GRIDana )

! ANA Moisture Variable
!----------------------
        FOUND = .false.
        do k=1,nq
           if( match('qv',REPLAY_QV,rnames(k)) ) then
               call ESMFL_BundleGetPointertoData(RBUNDLEP0,trim(rnames(k)),ptr3d, RC=STATUS)
               if(STATUS==ESMF_SUCCESS) then
                          q_rep =  ptr3d
                          if( currTime /= REPLAY_TIMEP0 ) then
                              call ESMFL_BundleGetPointertoData(RBUNDLEM1,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              q_rep =  facp0*q_rep + facm1*ptr3d
                              if( REPLAY_TIME_INTERP == "CUBIC" ) then
                                  call ESMFL_BundleGetPointertoData(RBUNDLEP1,trim(rnames(k)),ptr3d, RC=STATUS)
                                  VERIFY_(STATUS)
                                  q_rep =  q_rep + facp1*ptr3d
                              call ESMFL_BundleGetPointertoData(RBUNDLEM2,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              q_rep =  q_rep + facm2*ptr3d
                          !   q_rep =  max( q_rep, 0.0 )
                          endif
                          endif
                          if( REPLAY_QV_FACTOR.ne.1.0 ) q_rep = q_rep * REPLAY_QV_FACTOR
                          FOUND = .true.
                          exit
               endif
           endif
        enddo
        if( .not.FOUND .and. L_REPLAY_QV ) then
           write(STRING,'(A)') "ANA Variable: QV Not Found!"
           call WRITE_PARALLEL( trim(STRING)   )
           RETURN_(ESMF_FAILURE)
        endif

! ANA Ozone
!----------
        FOUND = .false.
        do k=1,nq
           if( match('o3',REPLAY_O3,rnames(k)) ) then
               call ESMFL_BundleGetPointertoData(RBUNDLEP0,trim(rnames(k)),ptr3d, RC=STATUS)
               if(STATUS==ESMF_SUCCESS) then
                          o3_rep =  ptr3d
                          if( currTime /= REPLAY_TIMEP0 ) then
                              call ESMFL_BundleGetPointertoData(RBUNDLEM1,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              o3_rep =  facp0*o3_rep + facm1*ptr3d
                              if( REPLAY_TIME_INTERP == "CUBIC" ) then
                                  call ESMFL_BundleGetPointertoData(RBUNDLEP1,trim(rnames(k)),ptr3d, RC=STATUS)
                                  VERIFY_(STATUS)
                                  o3_rep =  o3_rep + facp1*ptr3d
                              call ESMFL_BundleGetPointertoData(RBUNDLEM2,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              o3_rep =  o3_rep + facm2*ptr3d
                            ! o3_rep =  max( o3_rep, 0.0 )
                          endif
                          endif
                          if( REPLAY_O3_FACTOR.ne.1.0 ) o3_rep = o3_rep * REPLAY_O3_FACTOR
                          FOUND  = .true.
                          exit
               endif
           endif
        enddo
        if( .not.FOUND .and. L_REPLAY_O3 ) then
           write(STRING,'(A)') "ANA Variable: O3 Not Found!"
           call WRITE_PARALLEL( trim(STRING)   )
           RETURN_(ESMF_FAILURE)
        endif

! ANA Pressure Variables
! ----------------------
        if( LMana.eq.LMbkg ) then
            ak_rep = ak
            bk_rep = bk
        else
            call set_eta ( LMana,ksdum,ptopdum,pintdum,ak_rep,bk_rep )
        endif

        ple_rep(:,:,0) = ak_rep(0)
        do L=1,LMana
        ple_rep(:,:,L) = ple_rep(:,:,L-1) + dp_rep(:,:,L)
        enddo

        pke_rep(:,:,:) = ple_rep(:,:,:)**MAPL_KAPPA
        do L=1,LMana
         pk_rep(:,:,L) = ( pke_rep(:,:,L)-pke_rep(:,:,L-1) ) &
                        / ( MAPL_KAPPA*log(ple_rep(:,:,L)/ple_rep(:,:,L-1)) )
        enddo

! ANA Temperature Variable
! ------------------------
        if( trim(REPLAY_T).ne.'NO' ) then
        if( trim(REPLAY_T).ne.'YES' .and. trim(REPLAY_T_TYPE).eq.'NULL' ) then
            if(MAPL_AM_I_ROOT()) then
               print *
               print *, 'You must specify REPLAY_T_TYPE when setting REPLAY_T name.'
               print *, 'REPLAY_T_TYPE OPTIONS:  T, TV, TH, THV'
               print *
            endif
            RETURN_(ESMF_FAILURE)
        endif
        endif

        FOUND = .false.
        do k=1,nq
           if( match('t',REPLAY_T,rnames(k)) ) then
               call ESMFL_BundleGetPointertoData(RBUNDLEP0,trim(rnames(k)),ptr3d, RC=STATUS)
               if(STATUS==ESMF_SUCCESS) then
                  if( trim(REPLAY_T_TYPE).eq.'NULL' .or. trim(REPLAY_T_TYPE).eq.'T' ) then
                          t_rep =  ptr3d
                          if( currTime /= REPLAY_TIMEP0 ) then
                              call ESMFL_BundleGetPointertoData(RBUNDLEM1,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              t_rep =  facp0*t_rep + facm1*ptr3d
                              if( REPLAY_TIME_INTERP == "CUBIC" ) then
                                  call ESMFL_BundleGetPointertoData(RBUNDLEP1,trim(rnames(k)),ptr3d, RC=STATUS)
                                  VERIFY_(STATUS)
                                  t_rep =  t_rep + facp1*ptr3d
                              call ESMFL_BundleGetPointertoData(RBUNDLEM2,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              t_rep =  t_rep + facm2*ptr3d
                          endif
                          endif
                          if( REPLAY_T_FACTOR.ne.1.0 ) T_rep = T_rep * REPLAY_T_FACTOR
                          FOUND = .true.
                          exit
                  endif
               endif
           endif
           if( match('tv',REPLAY_T,rnames(k)) ) then
               call ESMFL_BundleGetPointertoData(RBUNDLEP0,trim(rnames(k)),ptr3d, RC=STATUS)
               if(STATUS==ESMF_SUCCESS) then
                  if( trim(REPLAY_T_TYPE).eq.'NULL' .or. trim(REPLAY_T_TYPE).eq.'TV' ) then
                          t_rep =  ptr3d
                          if( currTime /= REPLAY_TIMEP0 ) then
                              call ESMFL_BundleGetPointertoData(RBUNDLEM1,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              t_rep =  facp0*t_rep + facm1*ptr3d
                              if( REPLAY_TIME_INTERP == "CUBIC" ) then
                                  call ESMFL_BundleGetPointertoData(RBUNDLEP1,trim(rnames(k)),ptr3d, RC=STATUS)
                                  VERIFY_(STATUS)
                                  t_rep =  t_rep + facp1*ptr3d
                              call ESMFL_BundleGetPointertoData(RBUNDLEM2,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              t_rep =  t_rep + facm2*ptr3d
                          endif
                          endif
                          t_rep =  t_rep/(1.0+eps*q_rep)
                          if( REPLAY_T_FACTOR.ne.1.0 ) T_rep = T_rep * REPLAY_T_FACTOR
                          FOUND = .true.
                          exit
                  endif
               endif
           endif
           if( match('th',REPLAY_T,rnames(k)) ) then
               call ESMFL_BundleGetPointertoData(RBUNDLEP0,trim(rnames(k)),ptr3d, RC=STATUS)
               if(STATUS==ESMF_SUCCESS) then
                  if( trim(REPLAY_T_TYPE).eq.'NULL' .or. trim(REPLAY_T_TYPE).eq.'TH' ) then
                          t_rep =  ptr3d
                          if( currTime /= REPLAY_TIMEP0 ) then
                              call ESMFL_BundleGetPointertoData(RBUNDLEM1,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              t_rep =  facp0*t_rep + facm1*ptr3d
                              if( REPLAY_TIME_INTERP == "CUBIC" ) then
                                  call ESMFL_BundleGetPointertoData(RBUNDLEP1,trim(rnames(k)),ptr3d, RC=STATUS)
                                  VERIFY_(STATUS)
                                  t_rep =  t_rep + facp1*ptr3d
                              call ESMFL_BundleGetPointertoData(RBUNDLEM2,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              t_rep =  t_rep + facm2*ptr3d
                          endif
                          endif
                          t_rep =  t_rep*pk_rep
                          if( REPLAY_T_FACTOR.ne.1.0 ) T_rep = T_rep * REPLAY_T_FACTOR
                          FOUND = .true.
                          exit
                  endif
               endif
           endif
           if( match('thv',REPLAY_T,rnames(k)) ) then
               call ESMFL_BundleGetPointertoData(RBUNDLEP0,trim(rnames(k)),ptr3d, RC=STATUS)
               if(STATUS==ESMF_SUCCESS) then
                  if( trim(REPLAY_T_TYPE).eq.'NULL' .or. trim(REPLAY_T_TYPE).eq.'THV' ) then
                          t_rep =  ptr3d
                          if( currTime /= REPLAY_TIMEP0 ) then
                              call ESMFL_BundleGetPointertoData(RBUNDLEM1,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              t_rep =  facp0*t_rep + facm1*ptr3d
                              if( REPLAY_TIME_INTERP == "CUBIC" ) then
                                  call ESMFL_BundleGetPointertoData(RBUNDLEP1,trim(rnames(k)),ptr3d, RC=STATUS)
                                  VERIFY_(STATUS)
                                  t_rep =  t_rep + facp1*ptr3d
                              call ESMFL_BundleGetPointertoData(RBUNDLEM2,trim(rnames(k)),ptr3d, RC=STATUS)
                              VERIFY_(STATUS)
                              t_rep =  t_rep + facm2*ptr3d
                              endif
                          endif
                          t_rep =  t_rep*pk_rep/(1.0+eps*q_rep)
                          if( REPLAY_T_FACTOR.ne.1.0 ) T_rep = T_rep * REPLAY_T_FACTOR
                          FOUND = .true.
                          exit
                  endif
               endif
           endif
        enddo
        if( .not.FOUND .and. L_REPLAY_T ) then
            write(STRING,'(A)') "ANA Variable: T Not Found!"
            call WRITE_PARALLEL( trim(STRING)   )
            RETURN_(ESMF_FAILURE)
        endif

! Test for Re-Mapping
! -------------------
    if (doremap) then

        if( LMana.eq.LMbkg ) then
            NPHIS = count( phis_ana.ne.phis_bkg )
            call MAPL_CommsAllReduceMax(vm,sendbuf=NPHIS,recvbuf=NPHIS_MAX,cnt=1,rc=status)
            VERIFY_(STATUS)
        else
            NPHIS_MAX = 999 ! Force Vertical Remapping when LMana != LMbkg
        endif

            if( NPHIS_MAX > 0 ) then

                if(first .and. MAPL_AM_I_ROOT()) then
                   print *, 'Vertical Remapping ANA Data to BKG Topography and Levels ...'
                   print *
                endif
                thv_rep = t_rep*(1.0+eps*q_rep)/pk_rep
                call myremap ( ple_rep,  ple_ana, &
                                 u_rep,    u_ana, &
                                 v_rep,    v_ana, &
                               thv_rep,  thv_ana, &
                                 q_rep,    q_ana, &
                                o3_rep,   o3_ana, &
                              phis_ana, phis_bkg, &
                            ak_rep,bk_rep, ak,bk, &
                                im,jm,LMana,LMbkg )

                ! Create ANA Dry Temperature
                ! --------------------------
                   ps_ana(:,:)   = ple_ana(:,:,LMbkg)
                  pke_ana(:,:,:) = ple_ana(:,:,:)**MAPL_KAPPA
                  do L=1,LMbkg
                     pk_ana(:,:,L) = ( pke_ana(:,:,L)-pke_ana(:,:,L-1) ) &
                                   / ( MAPL_KAPPA*log(ple_ana(:,:,L)/ple_ana(:,:,L-1)) )
                  enddo
                      t_ana = thv_ana*pk_ana/(1.0+eps*q_ana)

            else

                if(first .and. MAPL_AM_I_ROOT()) then
                   print *, 'Vertical Remapping not necessary since ANA and BKG Topographies and Levels are identical.'
                   print *
                endif
                ple_ana = ple_rep
                  u_ana =   u_rep
                  v_ana =   v_rep
                  t_ana =   t_rep
                  q_ana =   q_rep
                 o3_ana =  o3_rep

            endif

    else
            if(first .and. MAPL_AM_I_ROOT()) then
               print *
               print *, 'Vertical Remapping ANA Data to BKG Topography and Levels is disabled.'
               print *
            endif
                ple_ana = ple_rep
                  u_ana =   u_rep
                  v_ana =   v_rep
                  t_ana =   t_rep
                  q_ana =   q_rep
                 o3_ana =  o3_rep
    endif

! **********************************************************************
! ****   Blend ANA and BKG Variables between DAMPBEG and DAMPEND,   ****
! ****   with option to blend QV specially, starting at tropopause. ****
! **********************************************************************

      if( DAMPBEG.ne.DAMPEND .or. BLEND_AT_PBL ) then

          if(first .and. MAPL_AM_I_ROOT()) then
             if(DAMPBEG.ne.DAMPEND) then
                print *, 'Blending Between ',DAMPBEG/100,' & ',DAMPEND/100,' mb'
             else
                print *, 'No Upper-Air ANA Blending to BKG will be done'
             endif
             if(BLEND_AT_PBL) then
                print *, 'Blending ANA and BKG based on PBL'
             else
                print *, 'No blending based on PBL'
             endif
             print *
          endif

          if( BLEND_AT_PBL ) then
             allocate ( pdum1(IMbkg,JMbkg,1) )
             allocate ( pdum2(IM,   JM,   1) )
             pdum1=0.0

             call MAPL_GetPointer(import, ptr2d, 'PPBL', RC=STATUS)
             VERIFY_(STATUS)
             pdum1(:,:,1) = ptr2d

             if (trim(GRIDINC)=="ANA" .and. do_transforms) then
                call mkiau_internal_state%bkg2ana_regridder%regrid(pdum1, pdum2, RC=STATUS)
                VERIFY_(STATUS)
             else
                pdum2=pdum1
             endif
             blnpp => pdum2(:,:,1)
          endif

          call blend ( ple_ana,u_ana,v_ana,t_ana,q_ana,o3_ana,     &
                       ple_bkg,u_bkg,v_bkg,t_bkg,q_bkg,o3_bkg,     &
                       im,jm,LMbkg, DAMPBEG,DAMPEND, BLEND_AT_PBL,  &
                       blnpp=blnpp )

          if( BLEND_AT_PBL ) then
             deallocate ( pdum1 )
             deallocate ( pdum2 )
          endif

      endif

    ! Modify Vertically Integrated Mass-Divergence Increment
    ! ------------------------------------------------------
    allocate ( vintdiva(IM,JM) )
    allocate ( vintdivb(IM,JM) )
    allocate ( vintdivc(IM,JM) )

    vintdiva = MAPL_UNDEF
    vintdivb = MAPL_UNDEF
    vintdivc = MAPL_UNDEF

    du_fix = u_ana
    dv_fix = v_ana

    if( trim(GRIDINC)=="ANA" .and. DOWINDFIX )  then
        if(first .and. MAPL_AM_I_ROOT()) then
           print *, 'Applying Mass Divergence Fix ...'
           print *
         endif
         method = 1
         call MAPL_TimerON(MAPL,"--WINDFIX")
         call windfix ( u_ana,v_ana,ple_ana,                            &
                        u_bkg,v_bkg,ple_bkg,im,jm,LMbkg,VM,GRIDana,method, &
                        vintdiva,vintdivb,vintdivc                      )
         call MAPL_TimerOFF(MAPL,"--WINDFIX")
    endif

    du_fix = u_ana-du_fix
    dv_fix = v_ana-dv_fix

! **********************************************************************
! ****                     Create IAU Increment                     ****
! **********************************************************************

    call MAPL_TimerON(MAPL,"-REGRIDSPC")

    call MAPL_GetPointer(export,   du, 'DUDT',  alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,   dv, 'DVDT',  alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,   dt, 'DTDT',  alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,   dq, 'DQVDT', alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,  do3, 'DO3DT', alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export, dple, 'DPEDT', alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,  dts, 'DTSDT', alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetPointer(export,duwindfix, 'DUWINDFIX', alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,dvwindfix, 'DVWINDFIX', alloc=.TRUE., RC=STATUS)
    VERIFY_(STATUS)

    allocate(  ptr3d(IM,JM,LMbkg),stat=STATUS)
    VERIFY_(STATUS)
    allocate( uptr3d(IM,JM,LMbkg),stat=STATUS)
    VERIFY_(STATUS)
    allocate( vptr3d(IM,JM,LMbkg),stat=STATUS)
    VERIFY_(STATUS)

    uptr3d = du_fix
    vptr3d = dv_fix

    if (trim(GRIDINC)=="ANA" .and. USE_SPECFILT .and. (L_REPLAY_U .or. L_REPLAY_V) ) then
        call Spectrans_VectorPar (im,jm,LMbkg,uptr3d,vptr3d,JCAP,GRIDana,RC=STATUS)
        VERIFY_(STATUS)
    endif
    if (trim(GRIDINC)=="ANA" .and. do_transforms) then
       if( L_CUBE ) then
          call mkiau_internal_state%ana2bkg_regridder%regrid(uptr3d, vptr3d, duwindfix, dvwindfix, rotate=.false., rc=status)
          VERIFY_(status)
       else
          call mkiau_internal_state%ana2bkg_regridder%regrid(uptr3d, duwindfix,rc=status)
          VERIFY_(status)
          call mkiau_internal_state%ana2bkg_regridder%regrid(vptr3d, dvwindfix,rc=status)
          VERIFY_(status)
          call POLEFIX ( duwindfix,dvwindfix,VM,GRIDbkg )
       endif
    else
       duwindfix=uptr3d
       dvwindfix=vptr3d
    endif

    if( L_REPLAY_U ) then
       uptr3d = u_ana-u_bkg
    else
       uptr3d = 0.0
    endif

    if( L_REPLAY_V ) then
       vptr3d = v_ana-v_bkg
    else
       vptr3d = 0.0
    endif

    if (trim(GRIDINC)=="ANA" .and. USE_SPECFILT .and. (L_REPLAY_U .or. L_REPLAY_V) ) then
        call Spectrans_VectorPar (im,jm,LMbkg,uptr3d,vptr3d,JCAP,GRIDana,RC=STATUS)
        VERIFY_(STATUS)
    endif
    if (trim(GRIDINC)=="ANA" .and. do_transforms) then
       if( L_CUBE ) then
          call mkiau_internal_state%ana2bkg_regridder%regrid(uptr3d, vptr3d, du, dv, rotate=.false., rc=status)
          VERIFY_(STATUS)
       else
          call mkiau_internal_state%ana2bkg_regridder%regrid(uptr3d, du, rc=status)
          VERIFY_(STATUS)
          call mkiau_internal_state%ana2bkg_regridder%regrid(vptr3d, dv, rc=status)
          VERIFY_(STATUS)
          call POLEFIX ( du,dv,VM,GRIDbkg )
       endif
    else
       du=uptr3d
       dv=vptr3d
    endif
    deallocate( uptr3d, vptr3d )

! Temperature
! -----------
    if( L_REPLAY_T ) then
        ptr3d = t_ana-t_bkg
    else
        ptr3d = 0.0
    endif
    if (trim(GRIDINC)=="ANA" .and. USE_SPECFILT .and. L_REPLAY_T) then
        call Spectrans_ScalarPar (im,jm,LMbkg,ptr3d,JCAP,GRIDana,RC=STATUS)
        VERIFY_(STATUS)
    endif
    if (trim(GRIDINC)=="ANA" .and. do_transforms) then
       call mkiau_internal_state%ana2bkg_regridder%regrid(ptr3d, dt, rc=status)
       VERIFY_(status)
    else
       dt=ptr3d
    endif

! Moisture
! --------
    if( L_REPLAY_QV ) then
        ptr3d = q_ana-q_bkg
    else
        ptr3d = 0.0
    endif
    if (trim(GRIDINC)=="ANA" .and. USE_SPECFILT .and. L_REPLAY_QV) then
        call Spectrans_ScalarPar (im,jm,LMbkg,ptr3d,JCAP,GRIDana,RC=STATUS)
        VERIFY_(STATUS)
    endif
    if (trim(GRIDINC)=="ANA" .and. do_transforms) then
       call mkiau_internal_state%ana2bkg_regridder%regrid(ptr3d, dq, rc=status)
       VERIFY_(status)
    else
       dq=ptr3d
    endif

! Ozone
! -----
    if( L_REPLAY_O3 ) then
        ptr3d = o3_ana-o3_bkg
    else
        ptr3d = 0.0
    endif
    if (trim(GRIDINC)=="ANA" .and. USE_SPECFILT .and. L_REPLAY_O3) then
        call Spectrans_ScalarPar (im,jm,LMbkg,ptr3d,JCAP,GRIDana,RC=STATUS)
        VERIFY_(STATUS)
    endif
    if (trim(GRIDINC)=="ANA" .and. do_transforms) then
       call mkiau_internal_state%ana2bkg_regridder%regrid(ptr3d, do3, rc=status)
       VERIFY_(status)
    else
       do3=ptr3d
    endif

! PLE
! ---
    deallocate( ptr3d )
      allocate( ptr3d(IM,JM,0:LMbkg),stat=STATUS)
    VERIFY_(STATUS)
    if( L_REPLAY_P ) then
        ptr3d = ple_ana-ple_bkg
    else
        ptr3d = 0.0
    endif
    if (trim(GRIDINC)=="ANA" .and. USE_SPECFILT .and. L_REPLAY_P) then
        call Spectrans_ScalarPar (im,jm,lmp1,ptr3d,JCAP,GRIDana,RC=STATUS)
        VERIFY_(STATUS)
    endif
    if (trim(GRIDINC)=="ANA" .and. do_transforms) then
       call mkiau_internal_state%ana2bkg_regridder%regrid(ptr3d, dple, rc=status)
       VERIFY_(status)
    else
       dple=ptr3d
    endif
    deallocate( ptr3d )

! 2-Dimensional EXPORTS need to be computed using 3-D Arrays
! ----------------------------------------------------------
    allocate ( qdum1(IMbkg,JMbkg,1) )
    allocate ( qdum2(IM,   JM,   1) )

    if( L_REPLAY_TS ) then
        qdum2(:,:,1) = ts_ana-ts_bkg
    else
        qdum2(:,:,1) = 0.0
    endif
    if (trim(GRIDINC)=="ANA" .and. USE_SPECFILT .and. L_REPLAY_TS) then
        call Spectrans_ScalarPar (im,jm,1,qdum2(1,1,1),JCAP,GRIDana,RC=STATUS)
        VERIFY_(STATUS)
    endif
    if (trim(GRIDINC)=="ANA" .and. do_transforms) then
       call mkiau_internal_state%ana2bkg_regridder%regrid(qdum2, qdum1, rc=status)
       VERIFY_(status)
    else
       qdum1=qdum2
    endif
    dts = qdum1(:,:,1)

! Vertically Integrated Divergence Increment Diagnostics
! ------------------------------------------------------
    call MAPL_GetPointer(export, vintdiv_ana, 'VINTDIV_ANA', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export, vintdiv_bkg, 'VINTDIV_BKG', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export, vintdiv_cor, 'VINTDIV_COR', RC=STATUS)
    VERIFY_(STATUS)

    if( associated(vintdiv_ana) ) then
        qdum2(:,:,1) = vintdiva
        if (trim(GRIDINC)=="ANA" .and. USE_SPECFILT .and. (L_REPLAY_U .or. L_REPLAY_V) ) then
            call Spectrans_ScalarPar (im,jm,1,qdum2(1,1,1),JCAP,GRIDana,RC=STATUS)
            VERIFY_(STATUS)
        endif
        if (trim(GRIDINC)=="ANA" .and. do_transforms) then
           call mkiau_internal_state%ana2bkg_regridder%regrid(qdum2, qdum1, rc=status)
           VERIFY_(status)
        else
           qdum1 = qdum2
        endif
        vintdiv_ana = qdum1(:,:,1)
    endif

    if( associated(vintdiv_bkg) ) then
        qdum2(:,:,1) = vintdivb
        if (trim(GRIDINC)=="ANA" .and. USE_SPECFILT .and. (L_REPLAY_U .or. L_REPLAY_V) ) then
            call Spectrans_ScalarPar (im,jm,1,qdum2(1,1,1),JCAP,GRIDana,RC=STATUS)
            VERIFY_(STATUS)
        endif
        if (trim(GRIDINC)=="ANA" .and. do_transforms) then
           call mkiau_internal_state%ana2bkg_regridder%regrid(qdum2, qdum1, rc=status)
           VERIFY_(status)
        else
           qdum1 = qdum2
        endif
        vintdiv_bkg = qdum1(:,:,1)
    endif

    if( associated(vintdiv_cor) ) then
        qdum2(:,:,1) = vintdivc
        if (trim(GRIDINC)=="ANA" .and. USE_SPECFILT .and. (L_REPLAY_U .or. L_REPLAY_V) ) then
            call Spectrans_ScalarPar (im,jm,1,qdum2(1,1,1),JCAP,GRIDana,RC=STATUS)
            VERIFY_(STATUS)
        endif
        if (trim(GRIDINC)=="ANA" .and. do_transforms) then
           call mkiau_internal_state%ana2bkg_regridder%regrid(qdum2, qdum1, rc=status)
           VERIFY_(status)
        else
           qdum1 = qdum2
        endif
        vintdiv_cor = qdum1(:,:,1)
    endif

    call MAPL_TimerOff(MAPL,"-REGRIDSPC")

    deallocate ( qdum1 )
    deallocate ( qdum2 )

! Clean-Up
! --------
    deallocate(  phis_bkg )
    deallocate(    ps_bkg )
    deallocate(    ts_bkg )
    deallocate(     u_bkg )
    deallocate(     v_bkg )
    deallocate(     t_bkg )
    deallocate(    tv_bkg )
    deallocate(     q_bkg )
    deallocate(    o3_bkg )
    deallocate (   dp_bkg )
    deallocate(   ple_bkg )

    deallocate (   RNAMES )
    deallocate (  ple_ana )
    deallocate ( phis_ana )
    deallocate (   ts_ana )
    deallocate (    u_ana )
    deallocate (    v_ana )
    deallocate (    t_ana )
    deallocate (    q_ana )
    deallocate (   o3_ana )
    deallocate (   ps_ana )
    deallocate (   pk_ana )
    deallocate (  pke_ana )
    deallocate (  thv_ana )
    deallocate (   du_fix )
    deallocate (   dv_fix )

    deallocate ( vintdiva )
    deallocate ( vintdivb )
    deallocate ( vintdivc )

    deallocate (   dp_rep )
    deallocate (    u_rep )
    deallocate (    v_rep )
    deallocate (    t_rep )
    deallocate (  thv_rep )
    deallocate (    q_rep )
    deallocate (   o3_rep )
    deallocate (   pk_rep )
    deallocate (  ple_rep )
    deallocate (  pke_rep )
    deallocate (   ak_rep )
    deallocate (   bk_rep )

    first = .false.
    end subroutine handleANA_

  end subroutine RUN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! ! IROUTINE: Clear -- Run method for the MAKEIAU component to clear increments

! !INTERFACE:

  subroutine CLEAR ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:


! ! DESCRIPTION: This version uses the MAPL_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP


! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),     pointer   :: MAPL
    type (ESMF_Field)                   :: FIELD

    real, pointer, dimension(:,:)       :: ptr2d
    real, pointer, dimension(:,:,:)     :: ptr3d

    integer                             :: I, N, fieldRank
    type(ESMF_FieldStatus_Flag)         :: fieldStatus
    character (len=ESMF_MAXSTR), allocatable  :: itemNameList(:)

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

   Iam = "Clear"
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

! get all export names; for each of them get pointer and zero it out

    call ESMF_StateGet(EXPORT, ITEMCOUNT=N, RC=STATUS)
    VERIFY_(STATUS)

    allocate(itemNameList(N), STAT=STATUS)
    VERIFY_(STATUS)

    call ESMF_StateGet(EXPORT, ItemNameList=itemNameList, RC=STATUS)
    VERIFY_(STATUS)

    do I = 1, N
       call ESMF_StateGet(Export, itemNameList(i), FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_FieldGet(field, status=fieldStatus, rc=status)
       VERIFY_(STATUS)
       if (fieldStatus /= ESMF_FIELDSTATUS_COMPLETE) cycle
       call ESMF_FieldGet(FIELD, dimCount=fieldRank, RC=status)
       VERIFY_(STATUS)
       if(fieldRank == 2) then
          call ESMF_FieldGet(field, farrayPtr=ptr2d, rc=status)
          VERIFY_(STATUS)
          if(associated(ptr2d)) ptr2d = 0.0
       else if (fieldRank == 3) then
          call ESMF_FieldGet(field, farrayPtr=ptr3d, rc=status)
          VERIFY_(STATUS)
          if(associated(ptr3d)) ptr3d = 0.0
       else
          _ASSERT(.false.,'needs informative message') ! not yet implemented
       endif
    end do

    deallocate(itemNameList)

    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)
  end subroutine CLEAR


! ***************************************************************************
  subroutine blend ( plea,ua,va,ta,qa,oa,     &
                     pleb,ub,vb,tb,qb,ob,     &
                     im,jm,lm, pabove,pbelow, &
                     BLEND_AT_PBL, blnpp    )

! Blends Anaylsis and Background values.
! This routine is called if pabove /= pbelow or BLEND_AT_PBL
! ***************************************************************************

      implicit none
      integer, intent(IN)    :: im,jm,lm
      real,    intent(IN)    :: pabove,pbelow
      logical, intent(IN)    :: BLEND_AT_PBL

      real,    intent(IN)    :: pleb(im,jm,lm+1)
      real,    intent(IN)    ::   ub(im,jm,lm)
      real,    intent(IN)    ::   vb(im,jm,lm)
      real,    intent(IN)    ::   tb(im,jm,lm)
      real,    intent(IN)    ::   qb(im,jm,lm)
      real,    intent(IN)    ::   ob(im,jm,lm)

      real,    intent(INOUT) :: plea(im,jm,lm+1)
      real,    intent(INOUT) ::   ua(im,jm,lm)
      real,    intent(INOUT) ::   va(im,jm,lm)
      real,    intent(INOUT) ::   ta(im,jm,lm)
      real,    intent(INOUT) ::   qa(im,jm,lm)
      real,    intent(INOUT) ::   oa(im,jm,lm)

      real,    intent(IN), optional, pointer :: blnpp(:,:)   ! blending pressure when BLEND_AT_PBL is TRUE

! Locals
! ------
      real pkea(im,jm,lm+1)
      real pkeb(im,jm,lm+1)
      real phia(im,jm,lm+1)
      real phib(im,jm,lm+1)

      real  thva(im,jm,lm)
      real  thvb(im,jm,lm)
      real   pka(im,jm,lm)
      real   pkb(im,jm,lm)

      real pabove_BL,pbelow_BL
      real bl_press

      real alf,eps,p
      integer i,j,L

       eps = MAPL_RVAP/MAPL_RGAS-1.0
      pkea = plea**MAPL_KAPPA
      pkeb = pleb**MAPL_KAPPA

      do L=1,lm
       pka(:,:,L) = ( pkea(:,:,L+1)-pkea(:,:,L) ) / ( MAPL_KAPPA*log(plea(:,:,L+1)/plea(:,:,L)) )
       pkb(:,:,L) = ( pkeb(:,:,L+1)-pkeb(:,:,L) ) / ( MAPL_KAPPA*log(pleb(:,:,L+1)/pleb(:,:,L)) )
      thva(:,:,L) = ta(:,:,L)*(1.0+eps*qa(:,:,L)) / pka(:,:,L)
      thvb(:,:,L) = tb(:,:,L)*(1.0+eps*qb(:,:,L)) / pkb(:,:,L)
      enddo

      phia(:,:,lm+1) = 0.0
      phib(:,:,lm+1) = 0.0
      do L=lm,1,-1
      phia(:,:,L) = phia(:,:,L+1) + MAPL_CP*thva(:,:,L)*( pkea(:,:,L+1)-pkea(:,:,L) )
      phib(:,:,L) = phib(:,:,L+1) + MAPL_CP*thvb(:,:,L)*( pkeb(:,:,L+1)-pkeb(:,:,L) )
      enddo

      if ( pabove /= pbelow ) then

! Blend mid-level u,v,q and o3
! ----------------------------
      do L=1,lm
      do j=1,jm
      do i=1,im
         p = 0.5*( plea(i,j,L)+plea(i,j,L+1) )
         if( p.le.pabove ) then
             alf = 0.0
         else if( p.gt.pabove .and. p.le.pbelow ) then
             alf = (p-pabove)/(pbelow-pabove)
         else
             alf = 1.0
         endif
                                   ua(i,j,L) =   ub(i,j,L) + alf*(   ua(i,j,L)-  ub(i,j,L) )
                                   va(i,j,L) =   vb(i,j,L) + alf*(   va(i,j,L)-  vb(i,j,L) )
                                   oa(i,j,L) =   ob(i,j,L) + alf*(   oa(i,j,L)-  ob(i,j,L) )
                                   qa(i,j,L) =   qb(i,j,L) + alf*(   qa(i,j,L)-  qb(i,j,L) )
      enddo
      enddo
      enddo

! Blend edge-level phi
! --------------------
      do L=1,lm+1
      do j=1,jm
      do i=1,im
         p = plea(i,j,L)
         if( p.le.pabove ) then
             alf = 0.0
         else if( p.gt.pabove .and. p.le.pbelow ) then
             alf = (p-pabove)/(pbelow-pabove)
         else
             alf = 1.0
         endif
         phia(i,j,L) = phib(i,j,L) + alf*( phia(i,j,L)-phib(i,j,L) )
      enddo
      enddo
      enddo

! Compute T based on blended phi
! ------------------------------
      do L=1,lm
        ta(:,:,L) = ( phia(:,:,L)-phia(:,:,L+1) )/( pkea(:,:,L+1)-pkea(:,:,L) ) &
                  / (MAPL_CP*(1.0+eps*qa(:,:,L))) * pka(:,:,L)
      enddo
      endif

! Blend from surface to blnpp
! -------------------------------------
      if ( BLEND_AT_PBL ) then
           do j=1,jm
           do i=1,im

           IF ( blnpp(i,j) == MAPL_UNDEF ) THEN
                pbelow_BL = 500.0 * 100.0   ! 500 hPa
           ELSE
                pbelow_BL = blnpp(i,j)
           ENDIF
           ! blend in the 200hPa above the PBL
           pabove_BL = pbelow_BL - 20000.0

           do L=1,lm
             p = 0.5*( pleb(i,j,L)+pleb(i,j,L+1) )
             if( p.le.pabove_BL ) then
                 alf = 0.0   !  use the analysis value
             else if( p.gt.pabove_BL .and. p.le.pbelow_BL ) then
                 alf = ((LOG(p)        -LOG(pabove_BL))/ &
                        (LOG(pbelow_BL)-LOG(pabove_BL)))**3
             else
                 alf = 1.0   !  use the background value
             endif

           plea(i,j,L) = plea(i,j,L) + alf*( pleb(i,j,L)- plea(i,j,L) )
             ua(i,j,L) =   ua(i,j,L) + alf*(   ub(i,j,L)-   ua(i,j,L) )
             va(i,j,L) =   va(i,j,L) + alf*(   vb(i,j,L)-   va(i,j,L) )
             ta(i,j,L) =   ta(i,j,L) + alf*(   tb(i,j,L)-   ta(i,j,L) )
             qa(i,j,L) =   qa(i,j,L) + alf*(   qb(i,j,L)-   qa(i,j,L) )
             oa(i,j,L) =   oa(i,j,L) + alf*(   ob(i,j,L)-   oa(i,j,L) )

           enddo
           plea(i,j,LM+1) = pleb(i,j,LM+1)

           enddo
           enddo
      endif

      return
  end subroutine blend


  subroutine POLEFIX (  ua,va,VM,GRID )
     implicit none
     real             :: ua(:,:,:)
     real             :: va(:,:,:)
     type (ESMF_VM)   :: VM
     type (ESMF_Grid) :: GRID
     integer          :: DIMS(ESMF_MAXGRIDDIM)
     integer          :: IM,JM,LM,IMG,JMG
     integer          :: vm_comm

     character(len=ESMF_MAXSTR)     :: IAm
     real, allocatable              :: uglo(:,:)
     real, allocatable              :: vglo(:,:)
     real, allocatable              :: sinl(:)
     real, allocatable              :: cosl(:)
     real                           :: LON,DL,UP,VP
     integer                        :: i,J,L,M,N,myid
     integer                        :: RC,STATUS

     IAM = "POLEFIX"

     call ESMF_VMGet  (VM, localpet=myid, RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GridGet(GRID, globalCellCountPerDim=DIMS, RC=STATUS)
     VERIFY_(STATUS)
     IMG = DIMS(1)
     JMG = DIMS(2)
     LM  = DIMS(3)
     DL  = 2*MAPL_PI/IMG

     allocate( uglo(IMG,JMG) )
     allocate( vglo(IMG,JMG) )
     allocate( sinl(IMG)     )
     allocate( cosl(IMG)     )

     do i=1,IMG
            LON = -MAPL_PI + (i-1)*DL
        cosl(i) = cos(LON)
        sinl(i) = sin(LON)
     enddo

     do L=1,LM
        call ArrayGather (local_array= ua(:,:,L),global_array= uglo(:,:), grid=GRID, rc=status)
        VERIFY_(STATUS)
        call ArrayGather (local_array= va(:,:,L),global_array= vglo(:,:), grid=GRID, rc=status)
        VERIFY_(STATUS)
        if( myid.eq.0 ) then
            do m = 1,2
                       N = (-1)**m
            if(m.eq.1) J = 1
            if(m.eq.2) J = JMG
               UP = 0.0
               VP = 0.0
               do i = 1,IMG
               UP = UP -   uglo(i,J-N)*sinl(i) - N*vglo(i,J-N)*cosl(i)
               VP = VP + N*uglo(i,J-N)*cosl(i) -   vglo(i,J-N)*sinl(i)
               enddo
               UP = UP / IMG
               VP = VP / IMG
               do i = 1, IMG
               uglo(i,J) = -   UP*sinl(i) + N*VP*cosl(i)
               vglo(i,J) = - N*UP*cosl(i) -   VP*sinl(i)
               enddo
            enddo
        endif
        call ArrayScatter (local_array=ua(:,:,L), global_array=uglo(:,:), grid=GRID, rc=status)
        VERIFY_(STATUS)
        call ArrayScatter (local_array=va(:,:,L), global_array=vglo(:,:), grid=GRID, rc=status)
        VERIFY_(STATUS)
     enddo

     deallocate( uglo )
     deallocate( vglo )
     deallocate( sinl )
     deallocate( cosl )

  end subroutine POLEFIX

      function match (replay_name,replay_alias,replay_var)
      character(*)               :: replay_name,replay_alias,replay_var
      character(len=ESMF_MAXSTR) :: name,alias,var
      logical         match
                      match = .false.

      name  = ESMF_UtilStringUpperCase( trim(replay_name ) )
      alias = ESMF_UtilStringUpperCase( trim(replay_alias) )
      var   = ESMF_UtilStringUpperCase( trim(replay_var  ) )

      if(     trim(var) == trim(alias) ) match = .true.

      if(     trim(name) == 'U'        ) then
          if( trim(var)  == 'U'        ) match = .true.
          if( trim(var)  == 'UGRD'     ) match = .true.
      endif

      if(     trim(name) == 'V'        ) then
          if( trim(var)  == 'V'        ) match = .true.
          if( trim(var)  == 'VGRD'     ) match = .true.
      endif

      if(     trim(name) == 'T'        ) then
          if( trim(var)  == 'T'        ) match = .true.
          if( trim(var)  == 'TMPU'     ) match = .true.
          if( trim(var)  == 'TMP'      ) match = .true.
      endif

      if(     trim(name) == 'QV'       ) then
          if( trim(var)  == 'Q'        ) match = .true.
          if( trim(var)  == 'QV'       ) match = .true.
          if( trim(var)  == 'SPHU'     ) match = .true.
      endif

      if(     trim(name) == 'TH'       ) then
          if( trim(var)  == 'TH'       ) match = .true.
          if( trim(var)  == 'THETA'    ) match = .true.
      endif

      if(     trim(name) == 'TV'       ) then
          if( trim(var)  == 'TV'       ) match = .true.
      endif

      if(     trim(name) == 'THV'      ) then
          if( trim(var)  == 'THV'      ) match = .true.
          if( trim(var)  == 'THETAV'   ) match = .true.
      endif

      if(     trim(name) == 'PS'       ) then
          if( trim(var)  == 'PS'       ) match = .true.
      endif

      if(     trim(name) == 'TS'       ) then
          if( trim(var)  == 'TS'       ) match = .true.
      endif

      if(     trim(name) == 'O3'       ) then
          if( trim(var)  == 'O3'       ) match = .true.
          if( trim(var)  == 'OZONE'    ) match = .true.
      endif

      if(     trim(name) == 'DP'       ) then
          if( trim(var)  == 'DP'       ) match = .true.
          if( trim(var)  == 'DELP'     ) match = .true.
      endif

      if(     trim(name) == 'PHIS'     ) then
          if( trim(var)  == 'PHIS'     ) match = .true.
          if( trim(var)  == 'GZ'       ) match = .true.
      endif

      return
      end function match

! subroutines for spectra filter
  subroutine spectrans_vector(im,jm,lm,var1,var2,jcap,GRID,RC)
  use mkiau_specmod, only: sptezv_s, sptez_s, init_spec_vars,destroy_spec_vars,isinitialized
  integer, intent(in)                      :: im,jm,lm
  integer, intent(in)                      :: jcap
  real, intent(inout), dimension(im,jm,lm) :: var1
  real, intent(inout), dimension(im,jm,lm) :: var2
  type(ESMF_Grid)                          :: GRID
  integer, optional, intent(out)           :: RC
! local variables
  integer                                  :: im_world, jm_world
  real, allocatable, dimension(:,:)        :: var_world1,var_world2
  real, allocatable, dimension(:)          :: grd1,grd2,spc1,spc2
  integer                                  :: k,nc
  integer                                  :: DIMS(ESMF_MAXGRIDDIM)
  integer                                  :: status
  character(len=ESMF_MAXSTR)               :: IAm

  IAm='spectrans_vector'
! Get global dimensions
  call MAPL_GridGet(GRID, globalCellCountPerDim=DIMS, RC=STATUS)
  VERIFY_(STATUS)
  IM_WORLD = DIMS(1)
  JM_WORLD = DIMS(2)
  nc = (jcap+1)*(jcap+2)
  isinitialized=.false.
  call init_spec_vars(im_world,jm_world,jcap,0)
  allocate(var_world1(im_world,jm_world))
  allocate(var_world2(im_world,jm_world))
  do k=1,lm
! gather array
   call ArrayGather(var1(:,:,k),var_world1,grid,rc=status)
   VERIFY_(STATUS)
   call ArrayGather(var2(:,:,k),var_world2,grid,rc=status)
   VERIFY_(STATUS)
   if (MAPL_AM_I_ROOT()) then
    allocate(grd1(im_world*jm_world))
    allocate(spc1(nc))
    allocate(grd2(im_world*jm_world))
    allocate(spc2(nc))
    grd1 = reshape(var_world1(:,:),(/im_world*jm_world/))
    grd2 = reshape(var_world2(:,:),(/im_world*jm_world/))
    call sptezv_s(spc1,spc2,grd1,grd2,-1)
    grd1 = 0.0
    grd2 = 0.0
    call sptezv_s(spc1,spc2,grd1,grd2,1)
    var_world1(:,:) = reshape(grd1,(/im_world,jm_world/))
    var_world2(:,:) = reshape(grd2,(/im_world,jm_world/))
    deallocate(grd1)
    deallocate(spc1)
    deallocate(grd2)
    deallocate(spc2)
   endif

! scatter array
   call ArrayScatter(var1(:,:,k),var_world1,grid,rc=status)
   VERIFY_(STATUS)
   call ArrayScatter(var2(:,:,k),var_world2,grid,rc=status)
   VERIFY_(STATUS)
  enddo

  call destroy_spec_vars
  deallocate(var_world1)
  deallocate(var_world2)
  RETURN_(ESMF_SUCCESS)
  end subroutine spectrans_vector

  subroutine spectrans_vectorpar(im,jm,lm,var1,var2,jcap,GRID,RC)
  use mkiau_specmod, only: sptezv_s, sptez_s, init_spec_vars,destroy_spec_vars,isinitialized
  integer, intent(in)                      :: im,jm,lm
  integer, intent(in)                      :: jcap
  real, intent(inout), dimension(im,jm,lm) :: var1
  real, intent(inout), dimension(im,jm,lm) :: var2
  type(ESMF_Grid)                          :: GRID
  integer, optional, intent(out)           :: RC
! local variables
  integer                                  :: im_world, jm_world
  integer                                  :: k,nc
  integer                                  :: DIMS(ESMF_MAXGRIDDIM)
  real, pointer                            :: InGlob1 (:,:,:)
  real, allocatable                        :: OutGlob1(:,:,:)
  real, pointer                            :: InGlob2 (:,:,:)
  real, allocatable                        :: OutGlob2(:,:,:)
  integer                                  :: status
  character(len=ESMF_MAXSTR)               :: IAm

  IAm='spectrans_vectorpar'
! Get global dimensions
  call MAPL_GridGet(GRID, globalCellCountPerDim=DIMS, RC=STATUS)
  VERIFY_(STATUS)
  IM_WORLD = DIMS(1)
  JM_WORLD = DIMS(2)
  nc = (jcap+1)*(jcap+2)
  isinitialized=.false.
  call init_spec_vars(im_world,jm_world,jcap,0)

  nullify(InGlob1)
  nullify(InGlob2)
  call MAPL_CollectiveGather3D(grid,var1,InGlob1,rc=status)
  VERIFY_(STATUS)
  call MAPL_CollectiveGather3D(grid,var2,InGlob2,rc=status)
  VERIFY_(STATUS)

  if (size(InGlob1) > 1) then
     allocate(OutGlob1(Im_World,JM_World,size(InGlob1,3)),stat=status)
     VERIFY_(STATUS)
     allocate(OutGlob2(Im_World,JM_World,size(InGlob2,3)),stat=status)
     VERIFY_(STATUS)
     call spectrans_vectorglob(InGlob1,InGlob2,OutGlob1,OutGlob2,jcap,rc=status)
     VERIFY_(STATUS)
  end if

  deallocate(InGlob1)
  deallocate(InGlob2)
  call MAPL_CollectiveScatter3d(Grid,OutGlob1,var1,rc=status)
  VERIFY_(STATUS)
  call MAPL_CollectiveScatter3d(Grid,OutGlob2,var2,rc=status)
  VERIFY_(STATUS)

  call destroy_spec_vars
  RETURN_(ESMF_SUCCESS)
  end subroutine spectrans_vectorpar

  subroutine spectrans_vectorglob(InGlob1,InGlob2,OutGlob1,OutGlob2,jcap,RC)
  use mkiau_specmod, only: sptezv_s, sptez_s, init_spec_vars,destroy_spec_vars,isinitialized
  integer, intent(in)                      :: jcap
  real, intent(inout), dimension(:,:,:)    :: InGlob1
  real, intent(inout), dimension(:,:,:)    :: OutGlob1
  real, intent(inout), dimension(:,:,:)    :: InGlob2
  real, intent(inout), dimension(:,:,:)    :: OutGlob2
  integer, optional, intent(out)           :: RC
! local variables
  integer                                  :: im_world, jm_world
  real, allocatable, dimension(:)          :: grd1,grd2,spc1,spc2
  integer                                  :: k,nc
  integer                                  :: status
  character(len=ESMF_MAXSTR)               :: IAm

  IAm='spectrans_vectorglob'
! Get global dimensions
  IM_WORLD = size(InGlob1,1)
  JM_WORLD = size(InGlob1,2)

  nc = (jcap+1)*(jcap+2)

  allocate(grd1(im_world*jm_world),stat=status)
  VERIFY_(STATUS)
  allocate(spc1(nc),stat=status)
  VERIFY_(STATUS)
  allocate(grd2(im_world*jm_world),stat=status)
  VERIFY_(STATUS)
  allocate(spc2(nc),stat=status)
  VERIFY_(STATUS)
  do k=1,size(InGlob1,3)
    grd1 = reshape(InGlob1(:,:,k),(/im_world*jm_world/))
    grd2 = reshape(InGlob2(:,:,k),(/im_world*jm_world/))
    call sptezv_s(spc1,spc2,grd1,grd2,-1)
    grd1 = 0.0
    grd2 = 0.0
    call sptezv_s(spc1,spc2,grd1,grd2,1)
    OutGlob1(:,:,k) = reshape(grd1,(/im_world,jm_world/))
    OutGlob2(:,:,k) = reshape(grd2,(/im_world,jm_world/))
  enddo
  deallocate(grd1)
  deallocate(spc1)
  deallocate(grd2)
  deallocate(spc2)

  RETURN_(ESMF_SUCCESS)
  end subroutine spectrans_vectorglob

  subroutine spectrans_scalar(im,jm,lm,var,jcap,GRID,RC)
  use mkiau_specmod, only: sptezv_s, sptez_s, init_spec_vars,destroy_spec_vars,isinitialized
  integer, intent(in)                      :: im,jm,lm
  integer, intent(in)                      :: jcap
  real, intent(inout), dimension(im,jm,lm) :: var
  type(ESMF_Grid)                          :: GRID
  integer, optional, intent(out)           :: RC
! local variables
  real, allocatable, dimension(:,:)        :: var_world
  real, allocatable, dimension(:)          :: grd,spc
  integer                                  :: k,nc
  integer                                  :: im_world, jm_world
  integer                                  :: DIMS(ESMF_MAXGRIDDIM)
  integer                                  :: status
  character(len=ESMF_MAXSTR)               :: IAm

  IAm='spectrans_scalar'
! Get global dimensions
  call MAPL_GridGet(GRID, globalCellCountPerDim=DIMS, RC=STATUS)
  VERIFY_(STATUS)
  IM_WORLD = DIMS(1)
  JM_WORLD = DIMS(2)
  nc = (jcap+1)*(jcap+2)
  isinitialized=.false.
  call init_spec_vars(im_world,jm_world,jcap,0)
  allocate(var_world(im_world,jm_world))
  do k=1,lm
   ! gather array
   call ArrayGather(var(:,:,k),var_world,grid,rc=status)
   VERIFY_(STATUS)
   if (MAPL_AM_I_ROOT()) then
    allocate(grd(im_world*jm_world))
    allocate(spc(nc))
    grd = reshape(var_world,(/im_world*jm_world/))
    call sptez_s(spc,grd,-1)
    grd = 0.0
    call sptez_s(spc,grd,1)
    var_world = reshape(grd,(/im_world,jm_world/))
    deallocate(grd)
    deallocate(spc)
   endif

   ! scatter array
   call ArrayScatter(var(:,:,k),var_world,grid,rc=status)
   VERIFY_(STATUS)
  enddo

  call destroy_spec_vars
  deallocate(var_world)
  RETURN_(ESMF_SUCCESS)
  end subroutine spectrans_scalar

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

  subroutine spectrans_scalarpar(im,jm,lm,var,jcap,GRID,RC)
  use mkiau_specmod, only: sptezv_s, sptez_s, init_spec_vars,destroy_spec_vars,isinitialized
  integer, intent(in)                      :: im,jm,lm
  integer, intent(in)                      :: jcap
  real, intent(inout), dimension(im,jm,lm) :: var
  type(ESMF_Grid)                          :: GRID
  integer, optional, intent(out)           :: RC
! local variables
  integer                                  :: im_world, jm_world
  integer                                  :: DIMS(ESMF_MAXGRIDDIM)
  real, pointer                            :: InGlob (:,:,:)
  real, allocatable                        :: OutGlob(:,:,:)
  integer                                  :: status
  character(len=ESMF_MAXSTR)               :: IAm

  IAm='spectrans_scalarpar'
! Get global dimensions
  call MAPL_GridGet(GRID, globalCellCountPerDim=DIMS, RC=STATUS)
  VERIFY_(STATUS)
  IM_WORLD = DIMS(1)
  JM_WORLD = DIMS(2)
  isinitialized=.false.
  call init_spec_vars(im_world,jm_world,jcap,0)

  nullify(InGlob)
  call MAPL_CollectiveGather3D(grid,var,InGlob,rc=status)
  VERIFY_(STATUS)

  if (size(InGlob) > 1) then
     allocate(OutGlob(Im_World,JM_World,size(InGlob,3)),stat=status)
     VERIFY_(STATUS)
     call spectrans_scalarglob(InGlob,OutGlob,jcap,rc=status)
     VERIFY_(STATUS)
  end if

  deallocate(InGlob)
  call MAPL_CollectiveScatter3d(Grid,OutGlob,var,rc=status)
  VERIFY_(STATUS)

  call destroy_spec_vars
  RETURN_(ESMF_SUCCESS)
  end subroutine spectrans_scalarpar

  ! perform scalar spectral transform on global data
  ! assuming calling routine already intitalized the spectral filter
  subroutine spectrans_scalarglob(InGlob,OutGLob,jcap,RC)
  use mkiau_specmod, only: sptezv_s, sptez_s, init_spec_vars,destroy_spec_vars,isinitialized
  integer, intent(in)                      :: jcap
  real, intent(inout), dimension(:,:,:   ) :: InGlob
  real, intent(inout), dimension(:,:,:   ) :: OutGlob
  integer, optional, intent(out)           :: RC
! local variables
  real, allocatable, dimension(:)          :: grd,spc
  integer                                  :: k,nc
  integer                                  :: im_world, jm_world
  integer                                  :: status
  character(len=ESMF_MAXSTR)               :: IAm

  IAm='spectrans_scalarglob'
! Get global dimensions
  nc = (jcap+1)*(jcap+2)
  im_world = size(InGlob,1)
  jm_world = size(InGlob,2)
  allocate(grd(im_world*jm_world),stat=status)
  VERIFY_(STATUS)
  allocate(spc(nc),stat=status)
  VERIFY_(STATUS)
  do k=1,size(InGlob,3)
   ! gather array
    grd = reshape(InGlob(:,:,k),(/im_world*jm_world/))
    call sptez_s(spc,grd,-1)
    grd = 0.0
    call sptez_s(spc,grd,1)
    OutGlob(:,:,k) = reshape(grd,(/im_world,jm_world/))
  enddo
  deallocate(grd)
  deallocate(spc)

  RETURN_(ESMF_SUCCESS)
  end subroutine spectrans_scalarglob

      subroutine myremap ( ple_in,ple_out,    &
                             u_in,  u_out,    &
                             v_in,  v_out,    &
                           thv_in,thv_out,    &
                            qv_in, qv_out,    &
                            o3_in, o3_out,    &
                            phis_in,phis_out, &
                            ak_in,bk_in, ak_out,bk_out,im,jm,LM_in,LM_out )

!***********************************************************************
!
!  Purpose
!     Driver for Remapping Fields to New Topography and Levels
!
!  Argument Description
!
!     ple_in ...... input edge pressure
!     u_in  ....... input zonal      wind
!     v_in  ....... input meridional wind
!     thv_in  ..... input virtual potential  temperature
!     qv_in ....... input specific   humidity
!     o3_in  ...... input ozone

!     ple_out...... output edge pressure
!     u_out ....... output zonal      wind
!     v_out ....... output meridional wind
!     thv_out ..... output virtual potential  temperature
!     qv_out ...... output specific   humidity
!     o3_out ...... output ozone

!     phis_in... input  surface geopotential
!     phis_out.. output surface geopotential
!     ak_in .... input  vertical   dimension
!     bk_in .... input  vertical   dimension
!     ak_out ... output vertical   dimension
!     bk_out ... output vertical   dimension
!
!     im ....... zonal      dimension
!     jm ....... meridional dimension
!     LM_in .... input  vertical dimension
!     LM_out ... output vertical dimension
!
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************

      use GEOS_GmapMod, only: gmap
      implicit none
      integer  im,jm,LM_in,LM_out

! Input variables
! ---------------
      real      ple_in(im,jm,LM_in+1)
      real        u_in(im,jm,LM_in)
      real        v_in(im,jm,LM_in)
      real      thv_in(im,jm,LM_in)
      real       qv_in(im,jm,LM_in)
      real       o3_in(im,jm,LM_in)

      real      ple_out(im,jm,LM_out+1)
      real        u_out(im,jm,LM_out)
      real        v_out(im,jm,LM_out)
      real      thv_out(im,jm,LM_out)
      real       qv_out(im,jm,LM_out)
      real       o3_out(im,jm,LM_out)

      real phis_in (im,jm)
      real phis_out(im,jm)

      real    ak_in (LM_in +1)
      real    bk_in (LM_in +1)
      real    ak_out(LM_out+1)
      real    bk_out(LM_out+1)

! Local variables
! ---------------
      real, allocatable ::  phi_in (:,:,:)
      real, allocatable ::  pke_in (:,:,:)

      real, allocatable ::   ps_out(:,:)
      real, allocatable ::  pke_out(:,:,:)

      real, allocatable ::    q_in (:,:,:,:)
      real, allocatable ::    q_out(:,:,:,:)

      real    kappa,cp,rgas,eps,rvap
      integer i,j,L

      kappa = 2.0/7.0
      rgas  = 8314.3/28.97
      rvap  = 8314.3/18.01
      eps   = rvap/rgas-1.0
      cp    = rgas/kappa

      allocate(  phi_in (im,jm,LM_in +1) )
      allocate(  pke_in (im,jm,LM_in +1) )

      allocate(  ps_out (im,jm)          )
      allocate(  pke_out(im,jm,LM_out+1) )

      allocate(    q_in (im,jm,LM_in ,2) )
      allocate(    q_out(im,jm,LM_out,2) )

! Construct Input Heights
! -----------------------
      pke_in(:,:,:) = ple_in(:,:,:)**kappa

      phi_in(:,:,LM_in+1) = phis_in(:,:)
      do L=LM_in,1,-1
      phi_in(:,:,L) = phi_in(:,:,L+1) + cp*thv_in(:,:,L)*( pke_in(:,:,L+1)-pke_in(:,:,L) )
      enddo

! Compute new surface pressure consistent with output topography
! --------------------------------------------------------------
      do j=1,jm
      do i=1,im
           L = LM_in
           do while ( phi_in(i,j,L).lt.phis_out(i,j) )
           L = L-1
           enddo
           ps_out(i,j) = ple_in(i,j,L+1)*( 1+(phi_in(i,j,L+1)-phis_out(i,j))/(cp*thv_in(i,j,L)*pke_in(i,j,L+1)) )**(1.0/kappa)
      enddo
      enddo

! Construct pressure variables using new surface pressure
! -------------------------------------------------------
      if( LM_in .eq. LM_out ) then
          do L=1,LM_in+1
          do j=1,jm
          do i=1,im
           ple_out(i,j,L) = ple_in(i,j,L) + bk_in(L)*( ps_out(i,j)-ple_in(i,j,LM_in+1) )
          enddo
          enddo
          enddo
      else
          do L=1,LM_out+1
          do j=1,jm
          do i=1,im
           ple_out(i,j,L) = ak_out(L) + bk_out(L)*ps_out(i,j)
          enddo
          enddo
          enddo
      endif

      pke_out(:,:,:) = ple_out(:,:,:)**kappa

! Map original fv state onto new eta grid
! ---------------------------------------
      q_in(:,:,:,1) = qv_in(:,:,:)
      q_in(:,:,:,2) = o3_in(:,:,:)

      call gmap ( im,jm,2 , kappa, &
                  LM_in,  pke_in ,ple_in ,u_in ,v_in ,thv_in ,q_in , &
                  LM_out, pke_out,ple_out,u_out,v_out,thv_out,q_out)

       qv_out(:,:,:) = q_out(:,:,:,1)
       o3_out(:,:,:) = q_out(:,:,:,2)

      deallocate(  phi_in  )
      deallocate(  pke_in  )
      deallocate(   ps_out )
      deallocate(  pke_out )

      deallocate( q_in  )
      deallocate( q_out )

      return
      end subroutine myremap

end module GEOS_mkiauGridCompMod
