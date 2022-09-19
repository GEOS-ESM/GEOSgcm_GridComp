! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_UW_InterfaceMod -- A Module to interface with the
!   UW convection

module GEOS_UW_InterfaceMod

  use ESMF
  use MAPL
  use UWSHCU   ! using module that contains uwshcu code

  implicit none

  integer USE_TRACER_TRANSP_UW      ! transport tracers in UW
  real    :: SCLM_SHALLOW
  logical :: JASON_UW

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

    call MAPL_GetResource(MAPL, USE_TRACER_TRANSP_UW,        'USE_TRACER_TRANSP_UW:',default= 1      , RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, JASON_UW,                    'JASON_UW:'            ,default= .FALSE., RC=STATUS) ; VERIFY_(STATUS)
    if (JASON_UW) then
      call MAPL_GetResource(MAPL, SHLWPARAMS%WINDSRCAVG,       'WINDSRCAVG:'      ,DEFAULT=0,      RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%MIXSCALE,         'MIXSCALE:'        ,DEFAULT=0.0,    RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%CRIQC,            'CRIQC:'           ,DEFAULT=1.0e-3, RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%THLSRC_FAC,       'THLSRC_FAC:'      ,DEFAULT= 0.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%FRC_RASN,         'FRC_RASN:'        ,DEFAULT= 0.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RKM,              'RKM:'             ,DEFAULT= 12.0,  RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SCLM_SHALLOW,                'SCLM_SHALLOW:'    ,DEFAULT= 1.0,   RC=STATUS) ; VERIFY_(STATUS)
    else
      call MAPL_GetResource(MAPL, SHLWPARAMS%WINDSRCAVG,       'WINDSRCAVG:'      ,DEFAULT=1,      RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%MIXSCALE,         'MIXSCALE:'        ,DEFAULT=3000.0, RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%CRIQC,            'CRIQC:'           ,DEFAULT=0.9e-3, RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%THLSRC_FAC,       'THLSRC_FAC:'      ,DEFAULT= 2.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%FRC_RASN,         'FRC_RASN:'        ,DEFAULT= 1.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RKM,              'RKM:'             ,DEFAULT= 8.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SCLM_SHALLOW,                'SCLM_SHALLOW:'    ,DEFAULT= 2.0,   RC=STATUS) ; VERIFY_(STATUS)
    endif
    call MAPL_GetResource(MAPL, SHLWPARAMS%NITER_XC,         'NITER_XC:'        ,DEFAULT=2,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%ITER_CIN,         'ITER_CIN:'        ,DEFAULT=2,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%USE_CINCIN,       'USE_CINCIN:'      ,DEFAULT=1,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%CRIDIST_OPT,      'CRIDIST_OPT:'     ,DEFAULT=0,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%USE_SELF_DETRAIN, 'USE_SELF_DETRAIN:',DEFAULT=0,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%USE_MOMENFLX,     'USE_MOMENFLX:'    ,DEFAULT=1,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%USE_CUMPENENT,    'USE_CUMPENENT:'   ,DEFAULT=1,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%SCVERBOSE,        'SCVERBOSE:'       ,DEFAULT=0,      RC=STATUS) ; VERIFY_(STATUS)
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
    call MAPL_GetResource(MAPL, SHLWPARAMS%QTSRC_FAC,        'QTSRC_FAC:'       ,DEFAULT= 0.0,   RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%QTSRCHGT,         'QTSRCHGT:'        ,DEFAULT=40.0,   RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RKFRE,            'RKFRE:'           ,DEFAULT= 1.0,   RC=STATUS) ; VERIFY_(STATUS)

end subroutine UW_Initialize

subroutine UW_Run (GC, IMPORT, EXPORT, CLOCK, RC)
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN
    real, pointer, dimension(:,:,:) :: QRAIN, QSNOW
    real, pointer, dimension(:,:)   :: CUSH

    ! Imports
    real, pointer, dimension(:,:)   :: FRLAND, SH, EVAP, KPBL_SC
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, TKE

    ! allocatable derived quantities
    real,    allocatable, dimension(:,:,:) :: ZLE0, ZL0
    real,    allocatable, dimension(:,:,:) :: PL, PK, PKE, DP
    real,    allocatable, dimension(:,:,:) :: MASS
    real,    allocatable, dimension(:,:,:) :: TH
    real,    allocatable, dimension(:,:,:) :: TMP3D
    real,    allocatable, dimension(:,:)   :: TMP2D

    ! Required Exports (connectivities to moist siblings)
    real, pointer, dimension(:,:)   :: CNPCPRATE

    ! Exports
    real, pointer, dimension(:,:,:) :: CUFRC_SC
    real, pointer, dimension(:,:,:) :: UMF_SC, MFD_SC, DCM_SC
    real, pointer, dimension(:,:,:) :: QTFLX_SC, SLFLX_SC, UFLX_SC, VFLX_SC
    real, pointer, dimension(:,:,:) :: DTDT_SC, DTHDT_SC, DQVDT_SC, DQRDT_SC, DQSDT_SC, &
                                       DUDT_SC,  DVDT_SC, DQIDT_SC, DQLDT_SC, DQADT_SC
    real, pointer, dimension(:,:,:) :: ENTR_SC, DETR_SC, QLDET_SC, &
                                       QIDET_SC, QLENT_SC, QIENT_SC, &
                                       QLSUB_SC, QISUB_SC, SC_NDROP, SC_NICE
    real, pointer, dimension(:,:)   :: TPERT_SC, QPERT_SC
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:)   :: PTR2D

    type (MAPL_MetaComp), pointer   :: MAPL
    type (ESMF_Config  )            :: CF
    type (ESMF_State   )            :: INTERNAL
    type (ESMF_Alarm   )            :: ALARM
    type (ESMF_TimeInterval)        :: TINT
    real(ESMF_KIND_R8)              :: DT_R8
    real                            :: DT_MOIST

    ! Local variables

    integer                         :: I, J, L
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

    ! Internals
    call MAPL_GetPointer(INTERNAL, Q,      'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,   'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,   'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLCN,   'CLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLLS,   'CLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,   'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,   'QICN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CUSH,   'CUSH'    , RC=STATUS); VERIFY_(STATUS)

    ! Imports
    call MAPL_GetPointer(IMPORT, FRLAND    ,'FRLAND'    ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, ZLE       ,'ZLE'       ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PLE       ,'PLE'       ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, T         ,'T'         ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, U         ,'U'         ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, V         ,'V'         ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SH        ,'SH'        ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, EVAP      ,'EVAP'      ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, KPBL_SC   ,'KPBL_SC'   ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TKE       ,'TKE'       ,RC=STATUS); VERIFY_(STATUS)

    ! Allocatables
     ! Edge variables 
    ALLOCATE ( ZLE0 (IM,JM,0:LM) )
    ALLOCATE ( PKE  (IM,JM,0:LM) )
     ! Layer variables
    ALLOCATE ( ZL0  (IM,JM,LM  ) )
    ALLOCATE ( PL   (IM,JM,LM  ) )
    ALLOCATE ( PK   (IM,JM,LM  ) )
    ALLOCATE ( TH   (IM,JM,LM  ) )
    ALLOCATE ( DP   (IM,JM,LM  ) )
    ALLOCATE ( MASS (IM,JM,LM  ) )
    ALLOCATE ( TMP3D(IM,JM,LM  ) )
     ! 2D Variables
    ALLOCATE ( TMP2D  (IM,JM) )

    ! Derived States
    PKE      = (PLE/MAPL_P00)**(MAPL_KAPPA)
    PL       = 0.5*(PLE(:,:,0:LM-1) + PLE(:,:,1:LM))
    PK       = (PL/MAPL_P00)**(MAPL_KAPPA)
    DO L=0,LM
       ZLE0(:,:,L)= ZLE(:,:,L) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
    END DO
    ZL0      = 0.5*(ZLE0(:,:,0:LM-1) + ZLE0(:,:,1:LM) ) ! Layer Height (m) above the surface
    TH       = T/PK
    DP       = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )
    MASS     = DP/MAPL_GRAV

    ! Required Exports (connectivities to moist siblings)
    call MAPL_GetPointer(EXPORT, MFD_SC,     'MFD_SC'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QLDET_SC,   'QLDET_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QIDET_SC,   'QIDET_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CUFRC_SC,   'CUFRC_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNPCPRATE,  'CNPCPRATE' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! Tendency Export
    call MAPL_GetPointer(EXPORT, DUDT_SC,    'DUDT_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DVDT_SC,    'DVDT_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DTDT_SC,    'DTDT_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DTHDT_SC,   'DTHDT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQVDT_SC,   'DQVDT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQIDT_SC,   'DQIDT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQLDT_SC,   'DQLDT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQRDT_SC,   'DQRDT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQSDT_SC,   'DQSDT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQADT_SC,   'DQADT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! Exports
    call MAPL_GetPointer(EXPORT, UMF_SC,     'UMF_SC'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DCM_SC,     'DCM_SC'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ENTR_SC,    'ENTR_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DETR_SC,    'DETR_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QLSUB_SC,   'QLSUB_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QISUB_SC,   'QISUB_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SC_NDROP,   'SC_NDROP'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SC_NICE,    'SC_NICE'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TPERT_SC,   'TPERT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QPERT_SC,   'QPERT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QTFLX_SC,   'QTFLX_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SLFLX_SC,   'SLFLX_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UFLX_SC,    'UFLX_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFLX_SC,    'VFLX_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

      !  Call UW shallow convection
      !----------------------------------------------------------------
      call compute_uwshcu_inv(IM*JM, LM,       DT_MOIST,  & ! IN
            PL, ZL0, PK, PLE, ZLE0, PKE, DP,              &
            U, V, Q, QLLS, QILS, TH, TKE, KPBL_SC,        &
            SH, EVAP, CNPCPRATE, FRLAND,                  &
            CUSH,                                         & ! INOUT
            UMF_SC, DCM_SC, DQVDT_SC, DQLDT_SC, DQIDT_SC, & ! OUT
            DTHDT_SC, DUDT_SC, DVDT_SC, DQRDT_SC,         &
            DQSDT_SC, CUFRC_SC, ENTR_SC, DETR_SC,         &
            QLDET_SC, QIDET_SC, QLSUB_SC, QISUB_SC,       &
            SC_NDROP, SC_NICE, TPERT_SC, QPERT_SC,        &
            QTFLX_SC, SLFLX_SC, UFLX_SC, VFLX_SC,         &
#ifdef UWDIAG 
            QCU_SC, QLU_SC,                               & ! DIAG ONLY 
            QIU_SC, CBMF_SC, SHL_DQCDT, CNT_SC, CNB_SC,   &
            CIN_SC, PLCL_SC, PLFC_SC, PINV_SC, PREL_SC,   &
            PBUP_SC, WLCL_SC, QTSRC_SC, THLSRC_SC,        &
            THVLSRC_SC, TKEAVG_SC, CLDTOP_SC, WUP_SC,     &
            QTUP_SC, THLUP_SC, THVUP_SC, UUP_SC, VUP_SC,  &
            XC_SC,                                        &
#endif 
            USE_TRACER_TRANSP_UW)

      !  Apply tendencies
      !--------------------------------------------------------------
        Q  = Q  + DQVDT_SC * DT_MOIST    ! note this adds to the convective
        TH = TH + DTHDT_SC * DT_MOIST    !  tendencies calculated below
        T  = TH*PK
        U  = U  +  DUDT_SC * DT_MOIST
        V  = V  +  DVDT_SC * DT_MOIST
      !  Update the temperature tendency
      !--------------------------------------------------------------
        DTDT_SC = DTHDT_SC*PK
      !  Calculate detrained mass flux
      !--------------------------------------------------------------
        where (DETR_SC.ne.MAPL_UNDEF)
          MFD_SC = 0.5*(UMF_SC(:,:,1:LM)+UMF_SC(:,:,0:LM-1))*DETR_SC*DP
        elsewhere
          MFD_SC = 0.0
        end where
       ! Tiedtke-style cloud fraction !!
        if (JASON_UW) then
          DQADT_SC= MFD_SC*SCLM_SHALLOW/MASS
        else
          DQADT_SC= DCM_SC*SCLM_SHALLOW/MASS
        endif
        CLCN = CLCN + DQADT_SC*DT_MOIST
        CLCN = MIN( CLCN , 1.0 )
      !  Convert detrained water units before passing to cloud
      !---------------------------------------------------------------
        call MAPL_GetPointer(EXPORT, QLENT_SC, 'QLENT_SC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, QIENT_SC, 'QIENT_SC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
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
       ! add detrained shallow convective ice/liquid source
        QLCN = QLCN + QLDET_SC*DT_MOIST
        QICN = QICN + QIDET_SC*DT_MOIST
       ! scale the detrained fluxes before exporting
        QLDET_SC = QLDET_SC*MASS
        QIDET_SC = QIDET_SC*MASS
      !  Apply condensate tendency from subsidence, and sink from
      !  condensate entrained into shallow updraft. 
      !-------------------------------------------------------------
        QLLS = QLLS + (QLSUB_SC+QLENT_SC)*DT_MOIST
        QILS = QILS + (QISUB_SC+QIENT_SC)*DT_MOIST
      !  Number concentrations for 2-moment microphysics
      !--------------------------------------------------------------
        SC_NDROP = SC_NDROP*MASS
        SC_NICE  = SC_NICE *MASS
      !  Precipitation
      !--------------------------------------------------------------
        call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_PRC3', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) PTR3D = DQRDT_SC    ! [kg/kg/s]
        call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_SNO3', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) PTR3D = DQSDT_SC    ! [kg/kg/s]
      ! Other exports
        call MAPL_GetPointer(EXPORT, PTR2D, 'SC_QT', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) then
        ! column integral of UW total water tendency, for checking conservation
        PTR2D = 0.
        DO L = 1,LM
           PTR2D = PTR2D + ( DQSDT_SC(:,:,L)+DQRDT_SC(:,:,L)+DQVDT_SC(:,:,L) &
                         +   QLENT_SC(:,:,L)+QLSUB_SC(:,:,L)+QIENT_SC(:,:,L) &
                         +   QISUB_SC(:,:,L) )*MASS(:,:,L) &
                         +  QLDET_SC(:,:,L)+QIDET_SC(:,:,L)
        END DO
        end if
        call MAPL_GetPointer(EXPORT, PTR2D, 'SC_MSE', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) then
        ! column integral of UW moist static energy tendency
        PTR2D = 0.
        DO L = 1,LM
           PTR2D = PTR2D + (MAPL_CP  *DTHDT_SC(:,:,L)*PK(:,:,L) &
                         +  MAPL_ALHL*DQVDT_SC(:,:,L)          &
                         -  MAPL_ALHF*DQIDT_SC(:,:,L))*MASS(:,:,L)
        END DO
        end if

       !--------------------------------------------------------------
       !  For Now add ShallowCu contribution to total/detraining mass flux exports
       !--------------------------------------------------------------
        call MAPL_GetPointer(EXPORT, PTR3D, 'CNV_MFC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        PTR3D = PTR3D + UMF_SC
        call MAPL_GetPointer(EXPORT, PTR3D, 'CNV_MFD', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        PTR3D = PTR3D + MFD_SC

        call MAPL_GetPointer(EXPORT, PTR2D,  'CUSH_SC', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) PTR2D = CUSH

    call MAPL_TimerOff (MAPL,"--UW")

end subroutine UW_Run

end module GEOS_UW_InterfaceMod
