! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_UW_InterfaceMod -- A Module to interface with the
!   UW convection

module GEOS_UW_InterfaceMod

  use ESMF
  use MAPL

  implicit none

  integer USE_TRACER_TRANSP_UW      ! transport tracers in UW
  type SHLWPARAM_TYPE
     integer  :: niter_xc           ! Number xc iterations
     integer  :: iter_cin           ! Number iterations for implicit CIN
     integer  :: use_CINcin         ! if true, calc CIN thru L..
     integer  :: use_self_detrain   ! 
     integer  :: use_momenflx       ! Perform momentum transport
     integer  :: use_cumpenent      ! Cumulative penetrative entrainment
     integer  :: scverbose          ! activate print statements
     integer  :: windsrcavg         ! Source air uses PBL mean momentum
     real     :: rpen               ! Penentrative entrainment factor
     real     :: rle
     real     :: rkm                ! Factor controlling lateral mixing rate
     real     :: mixscale           ! Controls vertical structure of mixing
     real     :: detrhgt            ! Mixing rate increases above this height
     real     :: rkfre              ! Vertical velocity fraction of tke
     real     :: rmaxfrac           ! Maximum core updraft fraction
     real     :: mumin1             ! 
     real     :: rbuoy              ! Non-hydro pressure effect on updraft
     real     :: rdrag              ! Drag coefficient
     real     :: epsvarw            ! Variance of PBL w by mesoscale
     real     :: PGFc               ! Pressure gradient force
     real     :: criqc              ! Updraft maximum condensate 
     real     :: frc_rasn           ! Precip fraction of expelled condensate
     real     :: kevp               ! Evaporative efficiency
     real     :: rdrop              ! liquid drop radius
     real     :: thlsrc_fac         ! Scaling factor for thlsrc perturbation
     real     :: qtsrc_fac          ! Scaling factor for qtsrc perturbation
     real     :: qtsrchgt           ! Interpolation height for total water source
     integer  :: cridist_opt
  endtype SHLWPARAM_TYPE
  type   (SHLWPARAM_TYPE) :: SHLWPARAMS

  private

  character(len=ESMF_MAXSTR)              :: IAm
  integer                                 :: STATUS

  public :: UW_Setup, UW_Initialize, UW_Run

contains

subroutine UW_Setup (GC, CF, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(ESMF_Config),   intent(inout) :: CF
    integer, optional                  :: RC  ! return code
    character(len=ESMF_MAXSTR)         :: COMP_NAME

    IAm = "GEOS_UW_InterfaceMod"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    call MAPL_AddInternalSpec(GC,                                   &
         SHORT_NAME ='PDF_A',                                       &
          LONG_NAME = 'SHOC_PDF_relative_area_fraction',            &
         UNITS      ='1',                                           &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationCenter,                         &
         DEFAULT= 0.5,                                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                    &
         SHORT_NAME ='CUSH',                                         &
         LONG_NAME  = 'Cumulus_scale_height_from_UW_shlw_convection',&
         UNITS      ='m',                                            &
         DIMS       = MAPL_DimsHorzOnly,                             &
         VLOCATION  = MAPL_VLocationNone,                            &
         DEFAULT= 1000.0,                                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_TimerAdd(GC, name="--UW", RC=STATUS)
    VERIFY_(STATUS)

end subroutine UW_Setup

subroutine UW_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

    call MAPL_GetResource(MAPL, USE_TRACER_TRANSP_UW,        'USE_TRACER_TRANSP_UW:',default= 1, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%NITER_XC,         'NITER_XC:'        ,DEFAULT=2,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%ITER_CIN,         'ITER_CIN:'        ,DEFAULT=2,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%USE_CINCIN,       'USE_CINCIN:'      ,DEFAULT=1,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%CRIDIST_OPT,      'CRIDIST_OPT:'     ,DEFAULT=0,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%USE_SELF_DETRAIN, 'USE_SELF_DETRAIN:',DEFAULT=0,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%USE_MOMENFLX,     'USE_MOMENFLX:'    ,DEFAULT=1,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%USE_CUMPENENT,    'USE_CUMPENENT:'   ,DEFAULT=1,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%SCVERBOSE,        'SCVERBOSE:'       ,DEFAULT=0,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%WINDSRCAVG,       'WINDSRCAVG:'      ,DEFAULT=1,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RPEN,             'RPEN:'            ,DEFAULT=3.0,    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RLE,              'RLE:'             ,DEFAULT=0.1,    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RMAXFRAC,         'RMAXFRAC:'        ,DEFAULT=0.1,    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%MUMIN1,           'MUMIN1:'          ,DEFAULT=0.906,  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RBUOY,            'RBUOY:'           ,DEFAULT=1.0,    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RDRAG,            'RDRAG:'           ,DEFAULT=1.0,    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%EPSVARW,          'EPSVARW:'         ,DEFAULT=5.e-4,  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%PGFC,             'PGFC:'            ,DEFAULT=0.7,    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%KEVP,             'KEVP:'            ,DEFAULT=2.e-6,  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RDROP,            'SHLW_RDROP:'      ,DEFAULT=8.e-6,  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%DETRHGT,          'DETRHGT:'         ,DEFAULT=1800.0, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%MIXSCALE,         'MIXSCALE:'        ,DEFAULT=3000.0, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%CRIQC,            'CRIQC:'           ,DEFAULT=0.9e-3, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%THLSRC_FAC,       'THLSRC_FAC:'      ,DEFAULT= 2.0,   RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%QTSRC_FAC,        'QTSRC_FAC:'       ,DEFAULT= 0.0,   RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%QTSRCHGT,         'QTSRCHGT:'        ,DEFAULT=40.0,   RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%FRC_RASN,         'FRC_RASN:'        ,DEFAULT= 1.0,   RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RKM,              'RKM:'             ,DEFAULT= 8.0,   RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RKFRE,            'RKFRE:'           ,DEFAULT= 1.0,   RC=STATUS) ; VERIFY_(STATUS)

end subroutine UW_Initialize

subroutine UW_Run (GC, IMPORT, EXPORT, CLOCK, RC)
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

    ! Local derived type aliases

    type (MAPL_MetaComp), pointer   :: MAPL
    type (ESMF_Config  )            :: CF
    type (ESMF_State   )            :: INTERNAL
    type (ESMF_Alarm   )            :: ALARM
    type (ESMF_TimeInterval)        :: TINT
    real(ESMF_KIND_R8)              :: DT_R8
    real                            :: DT_MOIST

    ! Local variables

    integer                         :: IM,JM,LM
    real, pointer, dimension(:,:)   :: LONS
    real, pointer, dimension(:,:)   :: LATS

    call ESMF_GridCompGet( GC, CONFIG=CF, RC=STATUS ) 
    VERIFY_(STATUS)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"--UW")

    ! Get parameters from generic state.
    !-----------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,   &
         RUNALARM = ALARM,             &
         CF       = CF,                &
         LONS     = LONS,              &
         LATS     = LATS,              &
         INTERNAL_ESMF_STATE=INTERNAL, &
         RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_AlarmGet(ALARM, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)
    DT_MOIST = DT_R8

    call MAPL_GetPointer(INTERNAL, CUSH, 'CUSH', RC=STATUS); VERIFY_(STATUS) 

      !  Call UW shallow convection
      !----------------------------------------------------------------
      call compute_uwshcu_inv(IDIM, K0, ITRCR, DT_MOIST,  & ! IN
            PLO*100., ZLO, PK, PLE, ZLE, PKE, DP,         &
            U1, V1, Q1, QLLS, QILS, TH1, TKE, KPBL_SC,    &
            SH, EVAP, CNVPRCP, FRLAND,                    &
            CUSH, XHO,                                    & ! INOUT
            UMF_SC, DCM_SC, DQVDT_SC, DQLDT_SC, DQIDT_SC, & ! OUT
            DTHDT_SC, DUDT_SC, DVDT_SC, DQRDT_SC,         &
            DQSDT_SC, CUFRC_SC, ENTR_SC, DETR_SC,         &
            QLDET_SC, QIDET_SC, QLSUB_SC, QISUB_SC,       &
            SC_NDROP, SC_NICE, TPERT_SC, QPERT_SC,        &
            QTFLX_SC, SLFLX_SC, UFLX_SC, VFLX_SC,         &
#ifdef UWDIAG 
            QCU_SC, QLU_SC,                               & ! DIAG ONLY 
            QIU_SC, CBMF_SC, DQCDT_SC, CNT_SC, CNB_SC,    &
            CIN_SC, PLCL_SC, PLFC_SC, PINV_SC, PREL_SC,   &
            PBUP_SC, WLCL_SC, QTSRC_SC, THLSRC_SC,        &
            THVLSRC_SC, TKEAVG_SC, CLDTOP_SC, WUP_SC,     &
            QTUP_SC, THLUP_SC, THVUP_SC, UUP_SC, VUP_SC,  &
            XC_SC,                                        &
#endif 
            USE_TRACER_TRANSP_UW, SHLWPARAMS )
      !  Apply tendencies
      !--------------------------------------------------------------
      Q1  = Q1  + DQVDT_SC * DT_MOIST    ! note this adds to the convective
      TH1 = TH1 + DTHDT_SC * DT_MOIST    !  tendencies calculated below
      U1  = U1  + DUDT_SC * DT_MOIST
      V1  = V1  + DVDT_SC * DT_MOIST
      if (associated(DTDT_SC)) DTDT_SC = DTHDT_SC*PK
      !  Calculate detrained mass flux
      !--------------------------------------------------------------
      where (DETR_SC.ne.MAPL_UNDEF)
        MFD_SC = 0.5*(UMF_SC(:,:,1:LM)+UMF_SC(:,:,0:LM-1))*DETR_SC*DP
      elsewhere
        MFD_SC = 0.0
      end where
      !  Convert detrained water units before passing to cloud
      !---------------------------------------------------------------
        QLENT_SC = 0.
        QIENT_SC = 0.
        WHERE (QLDET_SC.lt.0.)
          QLENT_SC = QLDET_SC
          QLDET_SC = 0.
        END WHERE
        WHERE (QIDET_SC.lt.0.)
          QIENT_SC = QIDET_SC
          QIDET_SC = 0.
        END WHERE
        QLDET_SC = QLDET_SC*MASS
        QIDET_SC = QIDET_SC*MASS
      !  Apply condensate tendency from subsidence, and sink from
      !  condensate entrained into shallow updraft. 
      !  Detrained condensate added in microphysics below.
      !-------------------------------------------------------------
        QLLS = QLLS + (QLSUB_SC+QLENT_SC)*DT_MOIST
        QILS = QILS + (QISUB_SC+QIENT_SC)*DT_MOIST
      !  Calculate updraft core fraction from cumulus fraction.
      !  CUFRC is assumed in compute_uwshcu to be twice updraft frac
      !--------------------------------------------------------------
      UFRC_SC = 0.5 * CUFRC_SC
      !  Number concentrations for 2-moment microphysics
      !--------------------------------------------------------------
      SC_NDROP = SC_NDROP*MASS
      SC_NICE = SC_NICE*MASS
      !  Precipitation
      !--------------------------------------------------------------
      SHLW_PRC3 = DQRDT_SC    ! [kg/kg/s]
      SHLW_SNO3 = DQSDT_SC    ! [kg/kg/s]
      if (associated(SC_QT)) then
        ! column integral of UW total water tendency, for checking conservation
        dum2d = 0.
        DO K = 1,LM
           dum2d = dum2d + (DQSDT_SC(:,:,k)+DQRDT_SC(:,:,k)+DQVDT_SC(:,:,k) &
                         + QLENT_SC(:,:,k)+QLSUB_SC(:,:,k)+QIENT_SC(:,:,k)  &
                         + QISUB_SC(:,:,k))*MASS(:,:,k)+QLDET_SC(:,:,k)     &
                         + QIDET_SC(:,:,k)
        END DO
        SC_QT = dum2d
      end if
      if (associated(SC_MSE)) then
        ! column integral of UW moist static energy tendency
        dum2d = 0.
        DO K = 1,LM
           dum2d = dum2d + (MAPL_CP*DTHDT_SC(:,:,k)*PK(:,:,k) &
                         + MAPL_ALHL*DQVDT_SC(:,:,k)          &
                         - MAPL_ALHF*DQIDT_SC(:,:,k))*MASS(:,:,k)
        END DO
        SC_MSE = dum2d
      end if
      if (associated(CUSH_SC)) CUSH_SC = CUSH

    call MAPL_TimerOff (MAPL,"--UW")

end subroutine UW_Run

end module GEOS_UW_InterfaceMod
