! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_RAS_InterfaceMod -- A Module to interface with the
!   RAS convection

module GEOS_RAS_InterfaceMod

  use ESMF
  use MAPL

  implicit none

  private

  character(len=ESMF_MAXSTR)              :: IAm
  integer                                 :: STATUS

  type RASPARAM_TYPE
      real               :: CUFRICFAC             ! 1
      real               :: SHR_LAMBDA_FAC        ! 2

      real               :: QC_CRIT_CN            ! 4
      real               :: RASAL1                ! 5
      real               :: RASAL2                ! 6
      real               :: RASNCL                ! 7
      real               :: LAMBDA_FAC            ! 8
      real               :: LAMBMX_FAC            ! 9
      real               :: MIN_DIAMETER          ! 10
      real               :: CUFRICLAMBDA          ! 11
      real               :: RDTLEXPON             ! 12
      real               :: STRAPPING             ! 13
      real               :: SDQV2                 ! 14
      real               :: SDQV3                 ! 15
      real               :: SDQVT1                ! 16
      real               :: ACRITFAC              ! 17
      real               :: HMINTRIGGER           ! 18
      real               :: LLDISAGGXP            ! 19
      real               :: PBLFRAC               ! 20
      real               :: RASAUTORAMPB          ! 21
      real               :: AUTOC_CN_ZDEP         ! 22
      real               :: MAXDALLOWED_S         ! 23
      real               :: MAXDALLOWED_D         ! 24
      real               :: MAXDALLOWED_E         ! 25
      real               :: RAS_RHMIN             ! 26
      real               :: RAS_RHFULL            ! 27
      real               :: CLDMICRO              ! 28
      real               :: FDROP_DUST            ! 29
      real               :: FDROP_SOOT            ! 30
      real               :: FDROP_SEASALT         ! 31
      real               :: RASAL_SLOPE
  endtype RASPARAM_TYPE
  type   (RASPARAM_TYPE) :: RASPARAMS

  public :: RAS_Setup, RAS_Initialize, RAS_Run

contains

subroutine RAS_Setup (GC, CF, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(ESMF_Config),   intent(inout) :: CF
    integer, optional                  :: RC  ! return code
    character(len=ESMF_MAXSTR)         :: COMP_NAME

    IAm = "GEOS_RAS_InterfaceMod"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    call MAPL_TimerAdd(GC, name="--RAS", RC=STATUS)
    VERIFY_(STATUS)

end subroutine RAS_Setup

subroutine RAS_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

#ifdef NODISABLE

      call MAPL_GetResource(MAPL, RAS_NO_NEG,              'RAS_NO_NEG:', default=.FALSE. , RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%CUFRICFAC,     'CUFRICFAC:',      DEFAULT= 1.000, RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%SHR_LAMBDA_FAC,'SHR_LAMBDA_FAC:', DEFAULT= 0.05,  RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%QC_CRIT_CN,    'QC_CRIT_CN:',    DEFAULT= 8.0e-4,  RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RASAL1,        'RASAL1:',        DEFAULT=   1800.0,RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RASAL2,        'RASAL2:',        DEFAULT= -86400.0,RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RASNCL,        'RASNCL:',        DEFAULT= -300.,   RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%LAMBDA_FAC,    'LAMBDA_FAC:',    DEFAULT= 4.0 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%LAMBMX_FAC,    'LAMBMX_FAC:',    DEFAULT= 0.0 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%STRAPPING,     'STRAPPING:',     DEFAULT=-1, RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%MIN_DIAMETER,  'MIN_DIAMETER:',  DEFAULT= 400.,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%CUFRICLAMBDA,  'CUFRICLAMBDA:',  DEFAULT= 7.5e-4,  RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RDTLEXPON,     'RDTLEXPON:',     DEFAULT= 1.0 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%SDQV2,         'SDQV2:',         DEFAULT= 1.3 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%SDQV3,         'SDQV3:',         DEFAULT= 1.3 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%SDQVT1,        'SDQVT1:',        DEFAULT= 263.,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%ACRITFAC,      'ACRITFAC:',      DEFAULT= 0.5 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%HMINTRIGGER,   'HMINTRIGGER:',   DEFAULT= 1.0 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%LLDISAGGXP,    'LLDISAGGXP:',    DEFAULT= 0.0 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%PBLFRAC,       'PBLFRAC:',       DEFAULT= 0.1 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RASAUTORAMPB,  'RASAUTORAMPB:',  DEFAULT= 0.8 ,    RC=STATUS)
      call MAPL_GetResource(MAPL,           AUTOC_CN_OCN,  'AUTOC_CN:',      DEFAULT= 2.5e-3,  RC=STATUS)
      call MAPL_GetResource(MAPL,           AUTOC_CN_LAND, 'AUTOC_CN_LAND:', DEFAULT= AUTOC_CN_OCN, RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%AUTOC_CN_ZDEP, 'AUTOC_CN_ZDEP:', DEFAULT= 1.0          ,RC=STATUS)
      if( imsize.le.200       ) call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT= 4000.0 ,RC=STATUS)
      if( imsize.gt.200 .and. &
          imsize.le.400       ) call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT= 2000.0 ,RC=STATUS)
      if( imsize.gt.400 .and. &
          imsize.le.800       ) call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT=  700.0 ,RC=STATUS)
      if( imsize.gt.800 .and. &
          imsize.le.1600      ) call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT=  450.0 ,RC=STATUS)
      if( imsize.gt.1600      ) call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT=  450.0 ,RC=STATUS)
                                call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_D, 'MAXDALLOWED_D:', DEFAULT= RASPARAMS%MAXDALLOWED_S ,RC=STATUS)
                                call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_E, 'MAXDALLOWED_E:', DEFAULT= -0.500 ,RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RASAL_SLOPE ,  'RASAL_SLOPE:',     DEFAULT= 8000.0  ,  RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RAS_RHMIN ,    'RAS_RHMIN:',       DEFAULT= 0.5  ,  RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RAS_RHFULL,    'RAS_RHFULL:',      DEFAULT= 0.65 ,  RC=STATUS)
      call MAPL_GetResource(MAPL, CBL_QPERT,               'CBL_QPERT:',       DEFAULT= 0.0   , RC=STATUS)
      call MAPL_GetResource(MAPL, CBL_TPERT,               'CBL_TPERT:',       DEFAULT=-1.0   , RC=STATUS)
      call MAPL_GetResource(MAPL, CBL_TPERT_MXOCN,         'CBL_TPERT_MXOCN:', DEFAULT= 2.0   , RC=STATUS)
      call MAPL_GetResource(MAPL, CBL_TPERT_MXLND,         'CBL_TPERT_MXLND:', DEFAULT= 0.0   , RC=STATUS)
#endif

end subroutine RAS_Initialize


subroutine RAS_Run (GC, IMPORT, EXPORT, CLOCK, RC)
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

    call MAPL_TimerOn (MAPL,"--RAS")

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

#ifdef NODISABLE

      IRAS       = nint(LONS*100)
      JRAS       = nint(LATS*100)
      RASPARAMS%CLDMICRO = 0.0
      if(adjustl(CLDMICRO)=="MGB2_2M") then
         RASPARAMS%CLDMICRO = 1.0
         RASPARAMS%FDROP_DUST = FDROP_DUST
         RASPARAMS%FDROP_SOOT = FDROP_SOOT
         RASPARAMS%FDROP_SEASALT = SS_INFAC
      endif
      SEEDINI(:,:,1) = 1000000 * ( 100*TEMP(:,:,LM)   - INT( 100*TEMP(:,:,LM) ) )
      SEEDINI(:,:,2) = 1000000 * ( 100*TEMP(:,:,LM-1) - INT( 100*TEMP(:,:,LM-1) ) )

      if (RASPARAMS%RASAL2 > 0.0) then
         RASAL2_2d(:,:) = RASPARAMS%RASAL2
      else
         ! include CNV dependence
         DO J=1, JM
            DO I=1, IM
            RASAL2_2d(I,J) = CNV_FRACTION(I,J)*ABS(RASPARAMS%RASAL2) + (1-CNV_FRACTION(I,J))*RASPARAMS%RASAL1
            END DO
         END DO
      endif
      ! Cheat by adding a kick to CB temp and q
      ! ---------------------------------------
      TPERT  = ABS(CBL_TPERT) * ( TS - ( TEMP(:,:,LM)+ MAPL_GRAV*ZLO(:,:,LM)/MAPL_CP )  )
      if (CBL_TPERT < 0) then
       ! Make TPERT 0 in areas of deep convection
        TPERT = TPERT*(1.0-CNV_FRACTION)
      endif
      QPERT  = CBL_QPERT * ( QSSFC - Q(:,:,LM) )
      TPERT  = MAX( TPERT , 0.0 )
      QPERT  = MAX( QPERT , 0.0 )
      where (FRLAND<0.1)
         TPERT = MIN( TPERT , CBL_TPERT_MXOCN ) ! ocean
      elsewhere
         TPERT = MIN( TPERT , CBL_TPERT_MXLND ) ! land
      end where
      SIGE = PREF/PREF(LM) ! this should eventually change

      ! MATMAT Export out the inputs before RAS needed for RAStest
      call MAPL_GetPointer(EXPORT, THOI,   'THOI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QHOI,   'QHOI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QSSI,   'QSSI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQSI,   'DQSI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PLEI,   'PLEI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TPERTI, 'TPERTI', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KCBLI,  'KCBLI' , RC=STATUS); VERIFY_(STATUS)
      if(associated(THOI  )) THOI   = TH1
      if(associated(QHOI  )) QHOI   = Q1
      if(associated(QSSI  )) QSSI   = QSS
      if(associated(DQSI  )) DQSI   = DQS
      if(associated(PLEI  )) PLEI   = CNV_PLE
      if(associated(TPERTI)) TPERTI = TPERT
      if(associated(KCBLI )) KCBLI  = KCBL
      ! temp kluge to differentiate ocean,land convective autoc (jtb 6/29/05)
      ! ---------------------------------------------------------------------
      where (FRLAND<0.1)
         CO_AUTO = AUTOC_CN_OCN   ! ocean value
      elsewhere
         CO_AUTO = AUTOC_CN_LAND  ! land value
      end where
      RAS_ALPHA   = MAPL_UNDEF
      RAS_TAU     = MAPL_UNDEF
      RAS_TIME    = MAPL_UNDEF
      RAS_TRG     = MAPL_UNDEF
      RAS_TOKI    = MAPL_UNDEF
      RAS_PBL     = MAPL_UNDEF
      RAS_WFN     = MAPL_UNDEF

      call RASE(                        &
           IDIM                 , &
           IRUN                 , &
           K0                   , &
           ICMIN                , &
           DT_MOIST             , &  !!where? see below.
           MAPL_CP              , &
           MAPL_ALHL            , &
           MAPL_ALHS            , &
           MAPL_TICE            , &
           MAPL_GRAV            , &
           SEEDINI              , &
           IRAS                 , &
           JRAS                 , &
           SIGE                 , &
           ! inputs for CBL
           KCBL                 , &
           WGT0                 , &
           WGT1                 , &
           ZCBLx                , &
           MXDIAMx              , &
           TPERT                , &
           QPERT                , &
           ! inputs
           TH1                  , &
           Q1                   , &
           U1                   , &
           V1                   , &
           QSS                  , &
           DQS                  , &
           CNV_FRACTION         , &
           RASAL2_2d            , &
           ! Pass in CO_AUTO
           CO_AUTO              , &
           ! - new for ras 2
           PK                   , &
           PLO                  , &
           GZLO                 , &
           GZLE                 , &
           QLCN                 , &
           QICN                 , &
           !
           CNV_PLE              , &
           PKE                  , &
           ! outputs
           CNV_DQLDT            , &   ! -> progno_cloud
           CNV_MF0              , &   ! -> diag
           CNV_MFD              , &   ! -> progno_cloud
           CNV_MFC              , &   ! -> diag
           CNV_PRC3             , &   ! -> progno_cloud 
           CNV_UPDF             , &   ! -> progno_cloud
           CNV_CVW              , &   ! 
           CNV_QC               , &   ! -> progno_cloud ???
           ENTLAM               , &   ! -> diag
           CLCN                 , &   ! -> upd if RAS-2 
           HHO                  , &
           HSO                  , &   ! -> diag
           CNVPRCP              , &

           RASPARAMS            , & ! params
           RAS_NO_NEG           , &
           RAS_TIME, RAS_TRG, RAS_TOKI, RAS_PBL, RAS_WFN, &
           RAS_TAU        , &

!!!=======AER_CLOUD=======

         !  AeroProps    , &  !-> Aerosol properties
           CNV_FICE       , & !-> Fraction of ice in detrainment 
           CNV_NICE       , & !-> Detrained ice crystal concentration
           CNV_NDROP      , & !-> Detrained cloud droplet concentration
           RAS_ALPHA      , &


!!!==============              
           ITRCR                , &
           irccode              , &
           XHO = XHO            , &
           TRIEDLEV_DIAG = trdlx, &
           FSCAV  = FSCAV       , &
           DISSKE = KEX           )

      if(associated(QVRAS  )) QVRAS    = Q1
      if(associated(THRAS  )) THRAS    = TH1
      if(associated(URAS   )) URAS     = U1
      if(associated(VRAS   )) VRAS     = V1
      if(associated(MXDIAM )) MXDIAM  = MXDIAMx
      if(associated(RCCODE )) RCCODE  = 1.0*IRCCODE
      if(associated(TRIEDLV)) TRIEDLV = TRDLX
      if(associated(TVEX   )) TVEX     = SUM( (MAPL_CP*TEMP + MAPL_ALHL*Q)*MASS, 3 )
#endif

    call MAPL_TimerOff(MAPL,"--RAS")

end subroutine RAS_Run

end module GEOS_RAS_InterfaceMod
