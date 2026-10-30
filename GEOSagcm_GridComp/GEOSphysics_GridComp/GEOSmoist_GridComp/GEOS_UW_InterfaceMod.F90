! $Id$

#include "MAPL.h"

!#define UWDIAG 1

!=============================================================================
!BOP

! !MODULE: GEOS_UW_InterfaceMod -- A Module to interface with the
!   UW convection

module GEOS_UW_InterfaceMod

  use ESMF
  use MAPL
  use UWSHCU   ! using module that contains uwshcu code
  use GEOSmoist_Process_Library
  use moist_dsl_workarounds

  implicit none

  integer USE_TRACER_TRANSP_UW      ! transport tracers in UW
  real    :: SCLM_SHALLOW
  logical :: JASON_UW, JASON_MFD_SC
  logical :: REPORT_UW_NEGATIVES
  logical :: USE_EIS, USE_SIG
  logical :: USE_PYMOIST_UW = .false.

  private

  integer                                 :: STATUS

  public :: UW_Setup, UW_Initialize, UW_Run, UW_Finalize
   
contains

subroutine UW_Setup (GC, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

    call MAPL_GridCompAddSpec(GC,                                          &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                     &
         SHORT_NAME       = 'CUSH',                                        &
         STANDARD_NAME    = 'Cumulus_scale_height_from_UW_shlw_convection', &
         UNITS            = 'm',                                           &
         DIMS             = 'xy',                                          &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_NONE,                    &
         FILL_VALUE       = 1000.0,                                        &
         _RC)

end subroutine UW_Setup

subroutine UW_Initialize (GC, CLOCK, IMPORT, EXPORT, RC)
    type (ESMF_GridComp), intent(inout) :: GC     ! gridded component
    type (ESMF_State),    intent(inout) :: IMPORT
    type (ESMF_State),    intent(inout) :: EXPORT
    type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional                   :: RC  ! return code
    integer :: LM

    real                            :: MOIST_DT
    real                            :: UW_DT

    type(ESMF_Calendar)     :: calendar
    type(ESMF_Time)         :: currentTime
    type(ESMF_Alarm)        :: UW_RunAlarm
    type(ESMF_Time)         :: ringTime
    type(ESMF_TimeInterval) :: ringInterval
    integer                 :: year, month, day, hh, mm, ss

    call MAPL_GridCompGet(GC, num_levels=LM, _RC)
    ! Moist's timestep (MAPL2: the RUNALARM ring interval)
    call MAPL_ClockGet(CLOCK, dt=MOIST_DT, _RC)
    call MAPL_GridCompGetResource( GC, 'UW_DT', UW_DT, default=MOIST_DT, _RC)

    call ESMF_ClockGet(CLOCK, calendar=calendar, currTime=currentTime, _RC)
    call ESMF_TimeGet(currentTime, YY=year, MM=month, DD=day, H=hh, M=mm, S=ss, _RC)
    call ESMF_TimeSet(ringTime, YY=year, MM=month, DD=day, H=0, M=0, S=0, _RC)
    call ESMF_TimeIntervalSet(ringInterval, S=nint(UW_DT), calendar=calendar, _RC)

    UW_RunAlarm = ESMF_AlarmCreate(Clock        = CLOCK,        &
                                   Name         = 'UW_RunAlarm', &
                                   RingInterval = ringInterval, &
                                   RingTime     = currentTime,  &
                                   Enabled      = .true.   ,    &
                                   Sticky       = .false.  , _RC)

                JASON_UW = .FALSE.
    if (LM==72) JASON_UW = .TRUE.
    call MAPL_GridCompGetResource( GC, 'JASON_UW', JASON_UW, default=(JASON_UW), _RC)

                JASON_MFD_SC = .FALSE.
    if (LM==72) JASON_MFD_SC = .TRUE.
    call MAPL_GridCompGetResource( GC, 'JASON_MFD_SC', JASON_MFD_SC, default=(JASON_MFD_SC), _RC)

    call MAPL_GridCompGetResource( GC, 'REPORT_UW_NEGATIVES', REPORT_UW_NEGATIVES, default=.FALSE., _RC)

    call MAPL_GridCompGetResource( GC, 'USE_PYMOIST_UW', USE_PYMOIST_UW, default=.FALSE., _RC)

    if (USE_PYMOIST_UW) then
      call MAPL_ConfigSetAttribute(CF, UW_DT, 'DSL__UW_DT:', _RC)
      call MAPL_pybridge_gcinit( "pyMoist.fortran.param_interfaces.convection.UW_interface", MAPL, IMPORT, EXPORT )
    else
    call MAPL_GridCompGetResource( GC, 'USE_TRACER_TRANSP_UW', USE_TRACER_TRANSP_UW,default= 1      , _RC)
    if (JASON_UW) then
      call MAPL_GridCompGetResource( GC, 'WINDSRCAVG', SHLWPARAMS%WINDSRCAVG      ,DEFAULT=0,      _RC)
      call MAPL_GridCompGetResource( GC, 'MIXSCALE', SHLWPARAMS%MIXSCALE        ,DEFAULT=0.0,    _RC)
      call MAPL_GridCompGetResource( GC, 'CRIQC', SHLWPARAMS%CRIQC           ,DEFAULT=1.0e-3, _RC)
      call MAPL_GridCompGetResource( GC, 'THLSRC_FAC', SHLWPARAMS%THLSRC_FAC      ,DEFAULT= 0.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'QTSRC_FAC', SHLWPARAMS%QTSRC_FAC       ,DEFAULT= 0.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'QTSRCHGT', SHLWPARAMS%QTSRCHGT        ,DEFAULT= 0.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'RKFRE', SHLWPARAMS%RKFRE           ,DEFAULT= 1.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'RKM', SHLWPARAMS%RKM             ,DEFAULT= 12.0,  _RC)
      call MAPL_GridCompGetResource( GC, 'RMAXFRAC', SHLWPARAMS%RMAXFRAC        ,DEFAULT= 0.1,   _RC)
      call MAPL_GridCompGetResource( GC, 'FRC_RASN', SHLWPARAMS%FRC_RASN        ,DEFAULT= 0.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'RPEN', SHLWPARAMS%RPEN            ,DEFAULT= 3.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'SCLM_SHALLOW', SCLM_SHALLOW    ,DEFAULT= 1.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'NITER_XC', SHLWPARAMS%NITER_XC        ,DEFAULT=2,      _RC)
      call MAPL_GridCompGetResource( GC, 'UW_USE_EIS', USE_EIS      ,DEFAULT=.FALSE.,_RC)
      call MAPL_GridCompGetResource( GC, 'UW_USE_SIG', USE_SIG      ,DEFAULT=.FALSE.,_RC)
    else
      call MAPL_GridCompGetResource( GC, 'WINDSRCAVG', SHLWPARAMS%WINDSRCAVG      ,DEFAULT=0,      _RC)
      call MAPL_GridCompGetResource( GC, 'MIXSCALE', SHLWPARAMS%MIXSCALE        ,DEFAULT=0.0,    _RC)
      call MAPL_GridCompGetResource( GC, 'CRIQC', SHLWPARAMS%CRIQC           ,DEFAULT=1.0e-3, _RC)
      call MAPL_GridCompGetResource( GC, 'THLSRC_FAC', SHLWPARAMS%THLSRC_FAC      ,DEFAULT= 0.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'QTSRC_FAC', SHLWPARAMS%QTSRC_FAC       ,DEFAULT= 0.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'QTSRCHGT', SHLWPARAMS%QTSRCHGT        ,DEFAULT= 0.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'RKFRE', SHLWPARAMS%RKFRE           ,DEFAULT= 1.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'RKFRE_HR', SHLWPARAMS%RKFRE_HR        ,DEFAULT= 0.75,  _RC)
      call MAPL_GridCompGetResource( GC, 'RKM', SHLWPARAMS%RKM             ,DEFAULT= 12.0,  _RC)
      call MAPL_GridCompGetResource( GC, 'RKM_HR', SHLWPARAMS%RKM_HR          ,DEFAULT= 12.0,  _RC)
      call MAPL_GridCompGetResource( GC, 'RMAXFRAC', SHLWPARAMS%RMAXFRAC        ,DEFAULT= 0.1,   _RC)
      call MAPL_GridCompGetResource( GC, 'RMAXFRAC_HR', SHLWPARAMS%RMAXFRAC_HR     ,DEFAULT= 0.1,   _RC)
      call MAPL_GridCompGetResource( GC, 'FRC_RASN', SHLWPARAMS%FRC_RASN        ,DEFAULT= 0.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'RPEN', SHLWPARAMS%RPEN            ,DEFAULT= 3.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'SCLM_SHALLOW', SCLM_SHALLOW    ,DEFAULT= 1.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'NITER_XC', SHLWPARAMS%NITER_XC        ,DEFAULT=2,      _RC)
      call MAPL_GridCompGetResource( GC, 'UW_USE_EIS', USE_EIS      ,DEFAULT=.FALSE.,_RC)
      call MAPL_GridCompGetResource( GC, 'UW_USE_SIG', USE_SIG      ,DEFAULT=.FALSE.,_RC)
    endif
    call MAPL_GridCompGetResource( GC, 'ITER_CIN', SHLWPARAMS%ITER_CIN        ,DEFAULT=2,      _RC)
    call MAPL_GridCompGetResource( GC, 'USE_CINCIN', SHLWPARAMS%USE_CINCIN      ,DEFAULT=1,      _RC)
    call MAPL_GridCompGetResource( GC, 'CRIDIST_OPT', SHLWPARAMS%CRIDIST_OPT     ,DEFAULT=0,      _RC)
    call MAPL_GridCompGetResource( GC, 'USE_SELF_DETRAIN', SHLWPARAMS%USE_SELF_DETRAIN,DEFAULT=0,      _RC)
    call MAPL_GridCompGetResource( GC, 'USE_MOMENFLX', SHLWPARAMS%USE_MOMENFLX    ,DEFAULT=1,      _RC)
    call MAPL_GridCompGetResource( GC, 'USE_CUMPENENT', SHLWPARAMS%USE_CUMPENENT   ,DEFAULT=1,      _RC)
    call MAPL_GridCompGetResource( GC, 'SCVERBOSE', SHLWPARAMS%SCVERBOSE       ,DEFAULT=0,      _RC)
    call MAPL_GridCompGetResource( GC, 'RLE', SHLWPARAMS%RLE             ,DEFAULT=0.1,    _RC)
    call MAPL_GridCompGetResource( GC, 'MUMIN1', SHLWPARAMS%MUMIN1          ,DEFAULT=0.906,  _RC)
    call MAPL_GridCompGetResource( GC, 'RBUOY', SHLWPARAMS%RBUOY           ,DEFAULT=1.0,    _RC)
    call MAPL_GridCompGetResource( GC, 'RDRAG', SHLWPARAMS%RDRAG           ,DEFAULT=1.0,    _RC)
    call MAPL_GridCompGetResource( GC, 'EPSVARW', SHLWPARAMS%EPSVARW         ,DEFAULT=5.e-4,  _RC)
    call MAPL_GridCompGetResource( GC, 'PGFC', SHLWPARAMS%PGFC            ,DEFAULT=0.7,    _RC)
    call MAPL_GridCompGetResource( GC, 'KEVP', SHLWPARAMS%KEVP            ,DEFAULT=2.e-6,  _RC)
    call MAPL_GridCompGetResource( GC, 'SHLW_RDROP', SHLWPARAMS%RDROP      ,DEFAULT=8.e-6,  _RC)
    call MAPL_GridCompGetResource( GC, 'DETRHGT', SHLWPARAMS%DETRHGT         ,DEFAULT=1800.0, _RC)

    endif ! USE_PYMOIST_UW

end subroutine UW_Initialize

subroutine UW_Run (GC, IMPORT, EXPORT, CLOCK, RC)
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN
    real, pointer, dimension(:,:)   :: CUSH

    ! Imports
    real, pointer, dimension(:,:)   :: FRLAND, SH, EVAP, KPBL_SC
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, TKE
    real, pointer, contiguous, dimension(:,:,:) :: p3d ! scratch for the edge-field bounds remap

    ! allocatable derived quantities
    real,    allocatable, dimension(:,:,:) :: ZLE0, ZL0
    real,    allocatable, dimension(:,:,:) :: PL, PK, PKE, DP
    real,    allocatable, dimension(:,:,:) :: MASS
    real,    allocatable, dimension(:,:,:) :: DQLDT_SC_, DQIDT_SC_
    real,    allocatable, dimension(:,:)   :: RKM2D, RKFRE, MIX2D, RMAXFRAC2D
    real,    allocatable, dimension(:,:,:) :: TMP3D

    ! Re-exported internals (were ALLOC'd exports shared with the driver and the other schemes)
    real, pointer, dimension(:,:)   :: CNPCPRATE
    real, pointer, dimension(:,:)   :: CNV_FRC, SRF_TYPE
    real, pointer, dimension(:,:,:) :: UMF_SC, MFD_SC
    real, pointer, dimension(:,:,:) :: DTDT_SC, DQVDT_SC, DQRDT_SC, DQSDT_SC, &
                                       DUDT_SC,  DVDT_SC, DQIDT_SC, DQLDT_SC, DQADT_SC
    real, pointer, dimension(:,:,:) :: QLDET_SC, QIDET_SC, QLENT_SC, QIENT_SC, &
                                       QLSUB_SC, QISUB_SC
    real, pointer, dimension(:,:)   :: LTS, EIS
    real, pointer, dimension(:,:,:) :: QLTOT, QITOT
    ! Exports computed in local arrays (were ALLOC'd), copied to the export if it is connected
    real, allocatable, dimension(:,:,:) :: CUFRC_SC, DCM_SC
    real, allocatable, dimension(:,:,:) :: QTFLX_SC, SLFLX_SC, UFLX_SC, VFLX_SC
    real, allocatable, dimension(:,:,:) :: ENTR_SC, DETR_SC, SC_NDROP, SC_NICE
    real, allocatable, dimension(:,:)   :: TPERT_SC, QPERT_SC
    real, allocatable, dimension(:,:)   :: CBMF_SC, PLCL_SC, PLFC_SC,     &
                                           PINV_SC, PREL_SC, PBUP_SC,    &
                                           CLDTOP_SC
    ! Their export pointers
    real, pointer, dimension(:,:,:) :: CUFRC_SC_ptr, DCM_SC_ptr
    real, pointer, dimension(:,:,:) :: QTFLX_SC_ptr, SLFLX_SC_ptr, UFLX_SC_ptr, VFLX_SC_ptr
    real, pointer, dimension(:,:,:) :: ENTR_SC_ptr, DETR_SC_ptr, SC_NDROP_ptr, SC_NICE_ptr
    real, pointer, dimension(:,:)   :: TPERT_SC_ptr, QPERT_SC_ptr
    real, pointer, dimension(:,:)   :: CBMF_SC_ptr, PLCL_SC_ptr, PLFC_SC_ptr,     &
                                       PINV_SC_ptr, PREL_SC_ptr, PBUP_SC_ptr,    &
                                       CLDTOP_SC_ptr
    ! Exports
    real, pointer, dimension(:,:)   :: SC_QT, SC_MSE
#ifdef UWDIAG
    real, pointer, dimension(:,:)   :: CIN_SC, CNT_SC, CNB_SC,      &
                                       WLCL_SC, QTSRC_SC, THLSRC_SC, &
                                       THVLSRC_SC, TKEAVG_SC
    real, pointer, dimension(:,:,:) :: SHL_DQCDT, WUP_SC, QTUP_SC,   &
                                       THLUP_SC, THVUP_SC, UUP_SC,   & 
                                       VUP_SC, XC_SC, QCU_SC,        &
                                       QLU_SC, QIU_SC
#endif
    real, pointer, dimension(:,:,:) ::   DQVDT_FILL
    real, pointer, dimension(:,:,:) :: DQLLSDT_FILL
    real, pointer, dimension(:,:,:) :: DQLCNDT_FILL
    real, pointer, dimension(:,:,:) :: DQILSDT_FILL
    real, pointer, dimension(:,:,:) :: DQICNDT_FILL
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:)   :: PTR2D

    real, allocatable, dimension(:,:) :: LONS, LATS

    type (ESMF_Grid    )            :: GRID
    type (ESMF_State   )            :: INTERNAL
    type (ESMF_TimeInterval)        :: TINT
    real(ESMF_KIND_R8)              :: DT_R8
    real                            :: UW_DT, MOIST_DT
    real                            :: DX, SIG, mix2d_phys
    type(ESMF_Alarm)                :: alarm
    logical                         :: alarm_is_ringing
    type( ESMF_VM )                 :: VMG

    ! Local variables

    real :: fac_eis                    ! Estimated enversion strength 0:1 factor
    real :: rkfre_base                 ! Base fractional entrainment rate before EIS modification
    real :: rkm_base                   ! Base momentum entrainment rate before EIS modification  
    real :: rkm_scale_fac
    real :: mix2d_base                 ! Base mixing length scale before EIS modification
    real :: rmaxfrac_base              ! Base maximum updraft area fraction before EIS modification
    real :: eis_rkfre_factor           ! EIS modification factor for RKFRE [0-1]
    real :: eis_rkm_factor             ! EIS modification factor for RKM
    real :: eis_mix2d_factor           ! EIS modification factor for MIX2D [0-1]
    real :: eis_rmaxfrac_factor        ! EIS modification factor for RMAXFRAC [1.0-1.1]

    integer                         :: I, J, L, K
    integer                         :: IM,JM,LM

    call MAPL_GridCompGet(GC, num_levels=LM, grid=GRID, _RC)
    call MAPL_GridGet(GRID, IM=IM, JM=JM, _RC)
    call MAPL_GridGetCoordinates(GRID, longitudes=LONS, latitudes=LATS, _RC)
    call MAPL_GridCompGetInternalState(GC, INTERNAL, _RC)
    ! Moist's timestep (MAPL2: the RUNALARM ring interval)
    call MAPL_ClockGet(CLOCK, dt=MOIST_DT, _RC)

    call ESMF_GridCompGet ( GC, VM=VMG, _RC )

    ! Internals
    call MAPL_StateGetPointer(INTERNAL, Q,      'Q'       , _RC)
    call MAPL_StateGetPointer(INTERNAL, QLLS,   'QLLS'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, QLCN,   'QLCN'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, CLCN,   'CLCN'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, CLLS,   'CLLS'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, QILS,   'QILS'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, QICN,   'QICN'    , _RC)
    ! Imports
    call MAPL_StateGetPointer(IMPORT, T,'T'         , _RC)
    call MAPL_StateGetPointer(IMPORT, U,'U'         , _RC)
    call MAPL_StateGetPointer(IMPORT, V,'V'         , _RC)
    ! Tendency Export (internals: they keep their values between UW_RunAlarm rings)
    call MAPL_StateGetPointer(INTERNAL, DUDT_SC,    'DUDT_SC'   , _RC)
    call MAPL_StateGetPointer(INTERNAL, DVDT_SC,    'DVDT_SC'   , _RC)
    call MAPL_StateGetPointer(INTERNAL, DTDT_SC,    'DTDT_SC'   , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQVDT_SC,   'DQVDT_SC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQIDT_SC,   'DQIDT_SC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQLDT_SC,   'DQLDT_SC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQRDT_SC,   'DQRDT_SC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQSDT_SC,   'DQSDT_SC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQADT_SC,   'DQADT_SC'  , _RC)
    ! Lower tropospheric stability and estimated inversion strength from MoistGC
    call MAPL_StateGetPointer(INTERNAL, LTS,   'LTS'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, EIS,   'EIS'  , _RC)

#ifdef UWDIAG
      call MAPL_StateGetPointer(EXPORT, SHL_DQCDT, 'SHL_DQCDT' , _RC)
      call MAPL_StateGetPointer(EXPORT, CNT_SC,    'CNT_SC'    , _RC)
      call MAPL_StateGetPointer(EXPORT, CNB_SC,    'CNB_SC'    , _RC)
      call MAPL_StateGetPointer(EXPORT, WLCL_SC,   'WLCL_SC'   , _RC)
      call MAPL_StateGetPointer(EXPORT, QTSRC_SC,  'QTSRC_SC'  , _RC)
      call MAPL_StateGetPointer(EXPORT, THLSRC_SC, 'THLSRC_SC' , _RC)
      call MAPL_StateGetPointer(EXPORT, THVLSRC_SC,'THVLSRC_SC', _RC)
      call MAPL_StateGetPointer(EXPORT, TKEAVG_SC, 'TKEAVG_SC' , _RC)
      call MAPL_StateGetPointer(EXPORT, WUP_SC,    'WUP_SC'    , _RC)
      call MAPL_StateGetPointer(EXPORT, QTUP_SC,   'QTUP_SC'   , _RC)
      call MAPL_StateGetPointer(EXPORT, THLUP_SC,  'THLUP_SC'  , _RC)
      call MAPL_StateGetPointer(EXPORT, THVUP_SC,  'THVUP_SC'  , _RC)
      call MAPL_StateGetPointer(EXPORT, UUP_SC,    'UUP_SC'    , _RC)
      call MAPL_StateGetPointer(EXPORT, VUP_SC,    'VUP_SC'    , _RC)
      call MAPL_StateGetPointer(EXPORT, XC_SC,     'XC_SC'     , _RC)
      call MAPL_StateGetPointer(EXPORT, QCU_SC,    'QCU_SC'    , _RC)
      call MAPL_StateGetPointer(EXPORT, QLU_SC,    'QLU_SC'    , _RC)
      call MAPL_StateGetPointer(EXPORT, QIU_SC,    'QIU_SC'    , _RC)
      ! Edge fields: MAPL3 pointers are 1:LM+1; MAPL2's were 0:LM (null unless requested)
      if (associated(WUP_SC))   then; p3d => WUP_SC;   WUP_SC(1:IM,1:JM,0:LM)   => p3d; end if
      if (associated(QTUP_SC))  then; p3d => QTUP_SC;  QTUP_SC(1:IM,1:JM,0:LM)  => p3d; end if
      if (associated(THLUP_SC)) then; p3d => THLUP_SC; THLUP_SC(1:IM,1:JM,0:LM) => p3d; end if
      if (associated(THVUP_SC)) then; p3d => THVUP_SC; THVUP_SC(1:IM,1:JM,0:LM) => p3d; end if
      if (associated(UUP_SC))   then; p3d => UUP_SC;   UUP_SC(1:IM,1:JM,0:LM)   => p3d; end if
      if (associated(VUP_SC))   then; p3d => VUP_SC;   VUP_SC(1:IM,1:JM,0:LM)   => p3d; end if
#endif
    
    call MAPL_GridCompTimerStart(GC, "--UW", _RC)
    
    if (USE_PYMOIST_UW) then
      call CNV_Tracers_To_SOA()
      call MAPL_pybridge_gcrun_with_internal( "pyMoist.fortran.param_interfaces.convection.UW_interface", MAPL, IMPORT, EXPORT, INTERNAL )
      call CNV_Tracers_To_AOS()
    else
      ! Internals
    call MAPL_StateGetPointer(INTERNAL, CUSH,   'CUSH'    , _RC)
   
    ! Imports
    call MAPL_StateGetPointer(IMPORT, FRLAND,'FRLAND'    , _RC)
    call MAPL_StateGetPointer(IMPORT, ZLE,'ZLE'       , _RC)
    call MAPL_StateGetPointer(IMPORT, PLE,'PLE'       , _RC)
    call MAPL_StateGetPointer(IMPORT, SH,'SH'        , _RC)
    call MAPL_StateGetPointer(IMPORT, EVAP,'EVAP'      , _RC)
    call MAPL_StateGetPointer(IMPORT, KPBL_SC,'KPBL_SC'   , _RC)
    call MAPL_StateGetPointer(IMPORT, TKE,'TKE'       , _RC)
    ! Edge fields: MAPL3 pointers are 1:LM+1; MAPL2's were 0:LM, which this code assumes
    p3d => ZLE; ZLE(1:IM,1:JM,0:LM) => p3d
    p3d => PLE; PLE(1:IM,1:JM,0:LM) => p3d
    p3d => TKE; TKE(1:IM,1:JM,0:LM) => p3d

    ! Allocatables
     ! Edge variables 
    ALLOCATE ( ZLE0 (IM,JM,0:LM) )
    ALLOCATE ( PKE  (IM,JM,0:LM) )
     ! Layer variables
    ALLOCATE ( ZL0  (IM,JM,LM  ) )
    ALLOCATE ( PL   (IM,JM,LM  ) )
    ALLOCATE ( PK   (IM,JM,LM  ) )
    ALLOCATE ( DP   (IM,JM,LM  ) )
    ALLOCATE ( MASS (IM,JM,LM  ) )
     ! Temporary UW exports
    ALLOCATE ( DQLDT_SC_(IM,JM,LM  ) )
    ALLOCATE ( DQIDT_SC_(IM,JM,LM  ) )
     ! 2D Variables
    ALLOCATE ( RKFRE  (IM,JM) )
    ALLOCATE ( RKM2D  (IM,JM) )
    ALLOCATE ( MIX2D  (IM,JM) )
    ALLOCATE ( RMAXFRAC2D (IM,JM) )

    ! Derived States
    !--------------------------------------------------------------
      
    ! 1. Compute the top-of-atmosphere edge (k = 0)
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP SHARED(IM, JM, LM, PKE, PLE, ZLE0, ZLE) &
    !$OMP PRIVATE(i, j)
    do j = 1, JM
       !DIR$ IVDEP
       do i = 1, IM
          PKE(i,j,0)  = (PLE(i,j,0) / MAPL_P00)**(MAPL_KAPPA)
          ZLE0(i,j,0) = ZLE(i,j,0) - ZLE(i,j,LM)
       end do
    end do
    !$OMP END PARALLEL DO

    ! 2. Compute the remaining edges and all layer variables (k = 1 to LM)
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP SHARED(IM, JM, LM, PKE, PLE, PL, PK, &
    !$OMP        ZLE0, ZLE, ZL0, DP, MASS) &
    !$OMP PRIVATE(i, j, k)
    do k = 1, LM
       do j = 1, JM
          !DIR$ IVDEP
          do i = 1, IM
               
             ! Edge variables
             PKE(i,j,k)  = (PLE(i,j,k) / MAPL_P00)**(MAPL_KAPPA)
             ZLE0(i,j,k) = ZLE(i,j,k) - ZLE(i,j,LM)
               
             ! Layer variables
             PL(i,j,k)   = 0.5 * (PLE(i,j,k-1) + PLE(i,j,k))
             PK(i,j,k)   = (PL(i,j,k) / MAPL_P00)**(MAPL_KAPPA)
               
             ZL0(i,j,k)  = 0.5 * (ZLE0(i,j,k-1) + ZLE0(i,j,k))
               
             DP(i,j,k)   = PLE(i,j,k) - PLE(i,j,k-1)
             MASS(i,j,k) = DP(i,j,k) / MAPL_GRAV
               
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    call ESMF_ClockGetAlarm(clock, 'UW_RunAlarm', alarm, _RC)
    alarm_is_ringing = ESMF_AlarmIsRinging(alarm, _RC)

    if (alarm_is_ringing) then

!!! call WRITE_PARALLEL('UW is Running')
    call ESMF_AlarmRingerOff(alarm, _RC)
    call ESMF_AlarmGet(alarm, RingInterval=TINT, _RC)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,_RC)
    UW_DT = DT_R8

    ! Required Exports (connectivities to moist siblings): re-exported internals
    call MAPL_StateGetPointer(INTERNAL, MFD_SC,     'MFD_SC'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, QLDET_SC,   'QLDET_SC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, QIDET_SC,   'QIDET_SC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, CNPCPRATE,  'CNPCPRATE' , _RC)
    call MAPL_StateGetPointer(INTERNAL, CNV_FRC ,   'CNV_FRC'   , _RC)
    call MAPL_StateGetPointer(INTERNAL, SRF_TYPE,   'SRF_TYPE'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, UMF_SC,     'UMF_SC'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, QLSUB_SC,   'QLSUB_SC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, QISUB_SC,   'QISUB_SC'  , _RC)
    ! Edge field: MAPL3 pointers are 1:LM+1; MAPL2's were 0:LM, which this code assumes
    p3d => UMF_SC; UMF_SC(1:IM,1:JM,0:LM) => p3d
    ! Exports computed locally, zero-filled as MAPL2's ALLOC did
    allocate(PLCL_SC(IM,JM), PLFC_SC(IM,JM), PINV_SC(IM,JM), PREL_SC(IM,JM), PBUP_SC(IM,JM), source=0.0)
    allocate(CBMF_SC(IM,JM), CLDTOP_SC(IM,JM), TPERT_SC(IM,JM), QPERT_SC(IM,JM), source=0.0)
    allocate(CUFRC_SC(IM,JM,LM), DCM_SC(IM,JM,LM), ENTR_SC(IM,JM,LM), DETR_SC(IM,JM,LM), source=0.0)
    allocate(SC_NDROP(IM,JM,LM), SC_NICE(IM,JM,LM), source=0.0)
    allocate(QTFLX_SC(IM,JM,0:LM), SLFLX_SC(IM,JM,0:LM), UFLX_SC(IM,JM,0:LM), VFLX_SC(IM,JM,0:LM), source=0.0) ! edge

    ! 1. Fetch all pointers first
    !--------------------------------------------------------------
    call MAPL_StateGetPointer(IMPORT, PTR2D, 'AREA', _RC)
    call MAPL_StateGetPointer(INTERNAL, QLTOT, 'QLTOT', _RC)
    call MAPL_StateGetPointer(INTERNAL, QITOT, 'QITOT', _RC)

    ! 2. 2D parameters for UW
    !--------------------------------------------------------------
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP SHARED(IM, JM, SHLWPARAMS, RKFRE, RKM2D, MIX2D, RMAXFRAC2D, &
    !$OMP        USE_SIG, USE_EIS, EIS, srf_type, PTR2D) &
    !$OMP PRIVATE(i, j, fac_eis, DX, SIG, rkfre_base, rkm_base, mix2d_base, &
    !$OMP        rmaxfrac_base, eis_rkfre_factor, eis_rkm_factor, eis_mix2d_factor, &
    !$OMP        eis_rmaxfrac_factor)
    do j = 1, JM
       !DIR$ IVDEP
       do i = 1, IM
         ! Optional sigma dependence
         if (USE_SIG) then 
            DX  = SQRT(PTR2D(i,j))
            SIG = SIGMA(DX)
         else
            SIG = 1.0 
         endif
         ! Optional EIS dependence
         if (USE_EIS) then
            fac_eis = get_fac_eis(EIS(i,j), srf_type(i,j))
         else
            fac_eis = 0.0
         endif
         ! Base resolution-dependent parameters
         ! Support for varying UW parameters by resolution  ! Coarse*SIG -> Fine*(1.0-SIG)
         rkfre_base    = SHLWPARAMS%RKFRE   *SIG  + SHLWPARAMS%RKFRE_HR   *(1.0-SIG)
         rkm_base      = SHLWPARAMS%RKM     *SIG  + SHLWPARAMS%RKM_HR     *(1.0-SIG) 
         mix2d_base    = SHLWPARAMS%MIXSCALE
         rmaxfrac_base = SHLWPARAMS%RMAXFRAC*SIG  + SHLWPARAMS%RMAXFRAC_HR*(1.0-SIG)
         ! EIS-based regime modifications for marine stratocumulus enhancement
         ! Reduce shallow convection activity in high EIS (stable inversion) regions
         eis_rkfre_factor    = 1.0 - 0.8*fac_eis           ! Reduce RKFRE by up to 80% in stable regimes
         eis_rkm_factor      = 1.0 + 0.2*fac_eis           ! Increase RKM by up to 40% in stable regimes
         eis_mix2d_factor    = 1.0 - 0.3*fac_eis           ! Reduce mixing scale by up to 30% in stable regimes
         eis_rmaxfrac_factor = 1.0 + 0.1*fac_eis           ! INCREASE rmaxfrac in stable (high EIS) regimes
         ! Apply EIS modifications
         RKFRE(i,j) = rkfre_base * eis_rkfre_factor
         RKM2D(i,j) = rkm_base   * eis_rkm_factor
         MIX2D(i,j) = mix2d_base * eis_mix2d_factor
         RMAXFRAC2D(i,j) = rmaxfrac_base * eis_rmaxfrac_factor
         ! Optional: Add minimum limits to prevent unrealistic values
         RKFRE(i,j) = max(RKFRE(i,j), 0.1)                 ! Minimum RKFRE threshold
         RKM2D(i,j) = min(RKM2D(i,j), 20.0)                ! Maximum RKM threshold
         MIX2D(i,j) = max(MIX2D(i,j), 0.0)                 ! Minimum mixing scale threshold
         RMAXFRAC2D(i,j) = max(min(RMAXFRAC2D(i,j), 0.9), 0.05) ! Bounds: 5% to 90%
       end do
    end do
    !$OMP END PARALLEL DO

    ! 3. Combine condensates for input (not updated within UW) 
    !--------------------------------------------------------------
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP SHARED(IM, JM, LM, QLTOT, DQLDT_SC, QLLS, QLCN, &
    !$OMP                    QITOT, DQIDT_SC, QILS, QICN) &
    !$OMP PRIVATE(i, j, k)
    do k = 1, LM
       do j = 1, JM
          !DIR$ IVDEP
          do i = 1, IM
             QLTOT(i,j,k) = QLLS(i,j,k) + QLCN(i,j,k)
             QITOT(i,j,k) = QILS(i,j,k) + QICN(i,j,k)
             ! Initialize tendencies
             DQLDT_SC(i,j,k) = QLTOT(i,j,k)
             DQIDT_SC(i,j,k) = QITOT(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

      !  Call UW shallow convection
      !----------------------------------------------------------------
      call compute_uwshcu_inv(IM*JM, LM, UW_DT,           & ! IN
            PL, ZL0, PK, PLE, ZLE0, PKE, DP,              &
            U, V, Q, QLTOT, QITOT, T, TKE, RKFRE, KPBL_SC,&
            SH, EVAP, CNPCPRATE, FRLAND, RKM2D, MIX2D, RMAXFRAC2D, &
            CUSH,                                         & ! INOUT
            UMF_SC, DCM_SC, DQVDT_SC, DQLDT_SC_, DQIDT_SC_, & ! OUT
            DTDT_SC, DUDT_SC, DVDT_SC, DQRDT_SC,          &
            DQSDT_SC, CUFRC_SC, ENTR_SC, DETR_SC,         &
            QLDET_SC, QIDET_SC, QLSUB_SC, QISUB_SC,       &
            SC_NDROP, SC_NICE, TPERT_SC, QPERT_SC,        &
            QTFLX_SC, SLFLX_SC, UFLX_SC, VFLX_SC,         &
            CBMF_SC, PLCL_SC, PLFC_SC, PINV_SC,           & ! DIAG ONLY
            PREL_SC, PBUP_SC, CLDTOP_SC,                  &
#ifdef UWDIAG 
            QCU_SC, QLU_SC,                               & 
            QIU_SC, SHL_DQCDT, CNT_SC, CNB_SC,            &
            CIN_SC, WLCL_SC, QTSRC_SC, THLSRC_SC,         &
            THVLSRC_SC, TKEAVG_SC, CLDTOP_SC, WUP_SC,     &
            QTUP_SC, THLUP_SC, THVUP_SC, UUP_SC, VUP_SC,  &
            XC_SC,                                        &
#endif 
            USE_TRACER_TRANSP_UW)

      ! 1. Fetch ALL pointers at the top
      !--------------------------------------------------------------
      call MAPL_StateGetPointer(INTERNAL, QLENT_SC, 'QLENT_SC', _RC)
      call MAPL_StateGetPointer(INTERNAL, QIENT_SC, 'QIENT_SC', _RC)
      call MAPL_StateGetPointer(EXPORT, SC_QT,   'SC_QT', _RC)
      call MAPL_StateGetPointer(EXPORT, SC_MSE,  'SC_MSE', _RC)

      ! 2. Fused 3D Loop for Detrainment and Conversions
      !--------------------------------------------------------------
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(IM, JM, LM, JASON_MFD_SC, DETR_SC, UMF_SC, DP, MFD_SC, &
      !$OMP        DCM_SC, DQADT_SC, SCLM_SHALLOW, MASS, QLENT_SC, QLDET_SC, QIENT_SC, QIDET_SC) &
      !$OMP PRIVATE(i, j, k)
      do k = 1, LM
         do j = 1, JM
            !DIR$ IVDEP 
            do i = 1, IM
               ! Calculate detrained mass flux
               if (JASON_MFD_SC) then
                  if (DETR_SC(i,j,k) /= MAPL_UNDEF) then
                     MFD_SC(i,j,k) = 0.5 * (UMF_SC(i,j,k) + UMF_SC(i,j,k-1)) * DETR_SC(i,j,k) * DP(i,j,k)
                  else
                     MFD_SC(i,j,k) = 0.0
                  end if
               else
                  MFD_SC(i,j,k) = DCM_SC(i,j,k)
               end if
               
               DQADT_SC(i,j,k) = MFD_SC(i,j,k) * SCLM_SHALLOW / MASS(i,j,k)

               ! Convert detrained water units before passing to cloud
               QLENT_SC(i,j,k) = 0.0
               QIENT_SC(i,j,k) = 0.0
               
               if (QLDET_SC(i,j,k) < 0.0) then
                  QLENT_SC(i,j,k) = QLDET_SC(i,j,k)
                  QLDET_SC(i,j,k) = 0.0
               end if
               
               if (QIDET_SC(i,j,k) < 0.0) then
                  QIENT_SC(i,j,k) = QIDET_SC(i,j,k)
                  QIDET_SC(i,j,k) = 0.0
               end if

               ! Scale the detrained fluxes before exporting
               QLDET_SC(i,j,k) = QLDET_SC(i,j,k) * MASS(i,j,k)
               QIDET_SC(i,j,k) = QIDET_SC(i,j,k) * MASS(i,j,k)
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      ! 3. Whole-array copies for direct 3D variables
      !--------------------------------------------------------------
      call MAPL_StateGetPointer(INTERNAL, PTR3D, 'SHLW_PRC3', _RC)
      if (associated(PTR3D)) PTR3D = DQRDT_SC    ! [kg/kg/s]
      call MAPL_StateGetPointer(INTERNAL, PTR3D, 'SHLW_SNO3', _RC)
      if (associated(PTR3D)) PTR3D = DQSDT_SC    ! [kg/kg/s]
      call MAPL_StateGetPointer(EXPORT, PTR2D,  'CUSH_SC', _RC)
      if (associated(PTR2D)) PTR2D = CUSH

      ! 4. Fused 2D Loop for Column Integrals
      !--------------------------------------------------------------
      ! We parallelize over j,i and accumulate over k internally
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(IM, JM, LM, SC_QT, SC_MSE, DQSDT_SC, DQRDT_SC, DQVDT_SC, &
      !$OMP        QLENT_SC, QLSUB_SC, QIENT_SC, QISUB_SC, MASS, QLDET_SC, QIDET_SC, &
      !$OMP        DTDT_SC, DQIDT_SC) &
      !$OMP PRIVATE(i, j, k)
      do j = 1, JM
         do i = 1, IM
            if (associated(SC_QT)) SC_QT(i,j) = 0.0
            if (associated(SC_MSE)) SC_MSE(i,j) = 0.0
            
            do k = 1, LM
               if (associated(SC_QT)) then
                  SC_QT(i,j) = SC_QT(i,j) + &
                     ( DQSDT_SC(i,j,k) + DQRDT_SC(i,j,k) + DQVDT_SC(i,j,k) + &
                       QLENT_SC(i,j,k) + QLSUB_SC(i,j,k) + QIENT_SC(i,j,k) + &
                       QISUB_SC(i,j,k) ) * MASS(i,j,k) + &
                       QLDET_SC(i,j,k) + QIDET_SC(i,j,k)
               end if
               
               if (associated(SC_MSE)) then
                  SC_MSE(i,j) = SC_MSE(i,j) + &
                     ( MAPL_CP   * DTDT_SC(i,j,k) + &
                       MAPL_ALHL * DQVDT_SC(i,j,k) - &
                       MAPL_ALHF * DQIDT_SC(i,j,k) ) * MASS(i,j,k)
               end if
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      ! 5. Copy the local arrays to their exports, if connected
      !--------------------------------------------------------------
      call MAPL_StateGetPointer(EXPORT, PLCL_SC_ptr,   'PLCL_SC'  , _RC)
      call MAPL_StateGetPointer(EXPORT, PLFC_SC_ptr,   'PLFC_SC'  , _RC)
      call MAPL_StateGetPointer(EXPORT, PINV_SC_ptr,   'PINV_SC'  , _RC)
      call MAPL_StateGetPointer(EXPORT, PREL_SC_ptr,   'PREL_SC'  , _RC)
      call MAPL_StateGetPointer(EXPORT, PBUP_SC_ptr,   'PBUP_SC'  , _RC)
      call MAPL_StateGetPointer(EXPORT, CBMF_SC_ptr,   'CBMF_SC'  , _RC)
      call MAPL_StateGetPointer(EXPORT, CLDTOP_SC_ptr, 'CLDTOP_SC', _RC)
      call MAPL_StateGetPointer(EXPORT, CUFRC_SC_ptr,  'CUFRC_SC' , _RC)
      call MAPL_StateGetPointer(EXPORT, DCM_SC_ptr,    'DCM_SC'   , _RC)
      call MAPL_StateGetPointer(EXPORT, ENTR_SC_ptr,   'ENTR_SC'  , _RC)
      call MAPL_StateGetPointer(EXPORT, DETR_SC_ptr,   'DETR_SC'  , _RC)
      call MAPL_StateGetPointer(EXPORT, SC_NDROP_ptr,  'SC_NDROP' , _RC)
      call MAPL_StateGetPointer(EXPORT, SC_NICE_ptr,   'SC_NICE'  , _RC)
      call MAPL_StateGetPointer(EXPORT, TPERT_SC_ptr,  'TPERT_SC' , _RC)
      call MAPL_StateGetPointer(EXPORT, QPERT_SC_ptr,  'QPERT_SC' , _RC)
      call MAPL_StateGetPointer(EXPORT, QTFLX_SC_ptr,  'QTFLX_SC' , _RC)
      call MAPL_StateGetPointer(EXPORT, SLFLX_SC_ptr,  'SLFLX_SC' , _RC)
      call MAPL_StateGetPointer(EXPORT, UFLX_SC_ptr,   'UFLX_SC'  , _RC)
      call MAPL_StateGetPointer(EXPORT, VFLX_SC_ptr,   'VFLX_SC'  , _RC)
      if (associated(PLCL_SC_ptr))   PLCL_SC_ptr   = PLCL_SC
      if (associated(PLFC_SC_ptr))   PLFC_SC_ptr   = PLFC_SC
      if (associated(PINV_SC_ptr))   PINV_SC_ptr   = PINV_SC
      if (associated(PREL_SC_ptr))   PREL_SC_ptr   = PREL_SC
      if (associated(PBUP_SC_ptr))   PBUP_SC_ptr   = PBUP_SC
      if (associated(CBMF_SC_ptr))   CBMF_SC_ptr   = CBMF_SC
      if (associated(CLDTOP_SC_ptr)) CLDTOP_SC_ptr = CLDTOP_SC
      if (associated(CUFRC_SC_ptr))  CUFRC_SC_ptr  = CUFRC_SC
      if (associated(DCM_SC_ptr))    DCM_SC_ptr    = DCM_SC
      if (associated(ENTR_SC_ptr))   ENTR_SC_ptr   = ENTR_SC
      if (associated(DETR_SC_ptr))   DETR_SC_ptr   = DETR_SC
      if (associated(SC_NDROP_ptr))  SC_NDROP_ptr  = SC_NDROP
      if (associated(SC_NICE_ptr))   SC_NICE_ptr   = SC_NICE
      if (associated(TPERT_SC_ptr))  TPERT_SC_ptr  = TPERT_SC
      if (associated(QPERT_SC_ptr))  QPERT_SC_ptr  = QPERT_SC
      if (associated(QTFLX_SC_ptr))  QTFLX_SC_ptr  = QTFLX_SC
      if (associated(SLFLX_SC_ptr))  SLFLX_SC_ptr  = SLFLX_SC
      if (associated(UFLX_SC_ptr))   UFLX_SC_ptr   = UFLX_SC
      if (associated(VFLX_SC_ptr))   VFLX_SC_ptr   = VFLX_SC

  endif


  ! 1. Fetch all pointers FIRST before doing any math
  !--------------------------------------------------------------
  call MAPL_StateGetPointer(INTERNAL, QLDET_SC, 'QLDET_SC', _RC)
  call MAPL_StateGetPointer(INTERNAL, QIDET_SC, 'QIDET_SC', _RC)
  call MAPL_StateGetPointer(INTERNAL, QLSUB_SC, 'QLSUB_SC', _RC)
  call MAPL_StateGetPointer(INTERNAL, QLENT_SC, 'QLENT_SC', _RC)
  call MAPL_StateGetPointer(INTERNAL, QISUB_SC, 'QISUB_SC', _RC)
  call MAPL_StateGetPointer(INTERNAL, QIENT_SC, 'QIENT_SC', _RC)
  
  ! 2. Apply tendencies in a single fused loop with OpenMP
  !--------------------------------------------------------------
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP SHARED(IM, JM, LM, Q, DQVDT_SC, MOIST_DT, T, DTDT_SC, U, DUDT_SC, V, DVDT_SC, &
  !$OMP        CLCN, DQADT_SC, QLCN, QLDET_SC, DQLDT_SC, MASS, QICN, QIDET_SC, DQIDT_SC, &
  !$OMP        QLLS, QLSUB_SC, QLENT_SC, QILS, QISUB_SC, QIENT_SC) &
  !$OMP PRIVATE(i, j, k)
  do k = 1, LM
     do j = 1, JM
        !DIR$ IVDEP 
        !DIR$ VECTOR ALWAYS
        do i = 1, IM
           ! Apply tendencies
           Q(i,j,k) = Q(i,j,k) + DQVDT_SC(i,j,k) * MOIST_DT
           T(i,j,k) = T(i,j,k) + DTDT_SC(i,j,k)  * MOIST_DT
           U(i,j,k) = U(i,j,k) + DUDT_SC(i,j,k)  * MOIST_DT
           V(i,j,k) = V(i,j,k) + DVDT_SC(i,j,k)  * MOIST_DT
           
           ! Tiedtke-style cloud fraction 
           CLCN(i,j,k) = MAX(0.0, MIN(CLCN(i,j,k) + DQADT_SC(i,j,k)*MOIST_DT, 1.0))
           
           ! Add detrained shallow convective ice/liquid source
           QLCN(i,j,k) = QLCN(i,j,k) + QLDET_SC(i,j,k)*MOIST_DT/MASS(i,j,k)
           QICN(i,j,k) = QICN(i,j,k) + QIDET_SC(i,j,k)*MOIST_DT/MASS(i,j,k)
           
           ! Apply condensate tendency from subsidence, and sink from
           ! condensate entrained into shallow updraft. 
           QLLS(i,j,k) = QLLS(i,j,k) + (QLSUB_SC(i,j,k)+QLENT_SC(i,j,k))*MOIST_DT
           QILS(i,j,k) = QILS(i,j,k) + (QISUB_SC(i,j,k)+QIENT_SC(i,j,k))*MOIST_DT

           ! Get export QL/QI tendencies
           DQLDT_SC(i,j,k) = (QLLS(i,j,k) + QLCN(i,j,k) - DQLDT_SC(i,j,k)) / MOIST_DT
           DQIDT_SC(i,j,k) = (QILS(i,j,k) + QICN(i,j,k) - DQIDT_SC(i,j,k)) / MOIST_DT
        end do
     end do
  end do
  !$OMP END PARALLEL DO

! Cleanup negative water species
! ------------------------------
  call MAPL_StateGetPointer(EXPORT,   DQVDT_FILL,   'DQVDT_FILL_SC', _RC)
  call MAPL_StateGetPointer(EXPORT, DQLLSDT_FILL, 'DQLLSDT_FILL_SC', _RC)
  call MAPL_StateGetPointer(EXPORT, DQLCNDT_FILL, 'DQLCNDT_FILL_SC', _RC)
  call MAPL_StateGetPointer(EXPORT, DQILSDT_FILL, 'DQILSDT_FILL_SC', _RC)
  call MAPL_StateGetPointer(EXPORT, DQICNDT_FILL, 'DQICNDT_FILL_SC', _RC)
  if (REPORT_UW_NEGATIVES) then
    call FILLQ2ZERO( Q       , MASS, DT=MOIST_DT, DQDT=  DQVDT_FILL, WARNING_LABEL="QV   After UW ShallowCu", VM=VMG, _RC)
    call FILLQ2ZERO( QLLS    , MASS, DT=MOIST_DT, DQDT=DQLLSDT_FILL, WARNING_LABEL="QLLS After UW ShallowCu", VM=VMG, _RC)
    call FILLQ2ZERO( QLCN    , MASS, DT=MOIST_DT, DQDT=DQLCNDT_FILL, WARNING_LABEL="QLCN After UW ShallowCu", VM=VMG, _RC)
    call FILLQ2ZERO( QILS    , MASS, DT=MOIST_DT, DQDT=DQILSDT_FILL, WARNING_LABEL="QILS After UW ShallowCu", VM=VMG, _RC)
    call FILLQ2ZERO( QICN    , MASS, DT=MOIST_DT, DQDT=DQICNDT_FILL, WARNING_LABEL="QICN After UW ShallowCu", VM=VMG, _RC)
  else
    call FILLQ2ZERO( Q       , MASS, DT=MOIST_DT, DQDT=  DQVDT_FILL, VM=VMG, _RC)
    call FILLQ2ZERO( QLLS    , MASS, DT=MOIST_DT, DQDT=DQLLSDT_FILL, VM=VMG, _RC)
    call FILLQ2ZERO( QLCN    , MASS, DT=MOIST_DT, DQDT=DQLCNDT_FILL, VM=VMG, _RC)
    call FILLQ2ZERO( QILS    , MASS, DT=MOIST_DT, DQDT=DQILSDT_FILL, VM=VMG, _RC)
    call FILLQ2ZERO( QICN    , MASS, DT=MOIST_DT, DQDT=DQICNDT_FILL, VM=VMG, _RC)
  endif


  if (DEBUG_TQ_ERRORS) then
        do L=1,LM                
          do J=1,JM              
           do I=1,IM             
             if (T(I,J,L) > 333.0) then
                 print *, "Temperature spike detected : ", T(I,J,L)
                 print *, "         UW Temp Increment : ", DTDT_SC(I,J,L)*MOIST_DT
                 print *, "    AFTER UW ShallowCu"
                 print *, "  Latitude       =", LATS(I,J)*180.0/MAPL_PI
                 print *, "  Longitude      =", LONS(I,J)*180.0/MAPL_PI
                 print *, "  Pressure (mb)  =", 0.5*(PLE(I,J,L)+PLE(I,J,L-1))/100.0                  
                 print *, "                       CLLS=",  CLLS(I,J,L),             "CLCN=", CLCN(I,J,L)
                 print *, "  QV=",     Q(I,J,L), "  QL=",  QLLS(I,J,L)+QLCN(I,J,L), "  QI=", QILS(I,J,L)+QICN(I,J,L)
             endif               
           end do ! IM loop      
         end do ! JM loop        
       end do ! LM loop          
  endif

  endif ! endif USE_PYMOIST_UW

  call MAPL_GridCompTimerStop(GC, "--UW", _RC)


end subroutine UW_Run

subroutine UW_Finalize(gc, import, export, rc)

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  integer, optional,   intent(  out) :: RC     ! Error code
  
  type (MAPL_MetaComp), pointer   :: MAPL
  
  ! Get my internal MAPL_Generic state
  !-----------------------------------
  call MAPL_GetObjectFromGC ( GC, MAPL, _RC)


  if (USE_PYMOIST_UW) then
    call MAPL_pybridge_gcfinalize( "pyMoist.fortran.param_interfaces.convection.UW_interface", MAPL, IMPORT, EXPORT )
  endif

end subroutine UW_Finalize

end module GEOS_UW_InterfaceMod
