! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_GFDL_1M_InterfaceMod -- A Module to interface with the
!   GFDL_1M cloud microphysics

module GEOS_GFDL_1M_InterfaceMod

  use ESMF
  use MAPL
  use gfdl2_cloud_microphys_mod

  implicit none

  logical :: LHYDROSTATIC
  logical :: LPHYS_HYDROSTATIC

  private

  character(len=ESMF_MAXSTR)              :: IAm
  integer                                 :: STATUS

  ! specify how to handle friendlies with DYN:TRB:CHM:ANA
  type FRIENDLIES_TYPE
         character(len=ESMF_MAXSTR) :: QV
         character(len=ESMF_MAXSTR) :: QW
         character(len=ESMF_MAXSTR) :: CLLS
         character(len=ESMF_MAXSTR) :: CLCN
         character(len=ESMF_MAXSTR) :: QLLS
         character(len=ESMF_MAXSTR) :: QLCN
         character(len=ESMF_MAXSTR) :: QILS
         character(len=ESMF_MAXSTR) :: QICN
         character(len=ESMF_MAXSTR) :: QRAIN
         character(len=ESMF_MAXSTR) :: QSNOW
         character(len=ESMF_MAXSTR) :: QGRAUPEL
  end type FRIENDLIES_TYPE
  type (FRIENDLIES_TYPE) FRIENDLIES

  character(len=ESMF_MAXSTR)        :: COMP_NAME

  public :: GFDL_1M_Setup, GFDL_1M_Initialize, GFDL_1M_Run
  public :: LHYDROSTATIC, LPHYS_HYDROSTATIC

contains

subroutine GFDL_1M_Setup (GC, CF, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(ESMF_Config),   intent(inout) :: CF
    integer, optional                  :: RC  ! return code
    character(len=ESMF_MAXSTR)         :: COMP_NAME

    IAm = "GEOS_GFDL_1M_InterfaceMod"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    ! !INTERNAL STATE:

      FRIENDLIES%QV       = trim(COMP_NAME)//":DYNAMICS:TURBULENCE:CHEMISTRY:ANALYSIS"
      FRIENDLIES%QW       = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%CLLS     = trim(COMP_NAME)//":DYNAMICS"
      FRIENDLIES%CLCN     = trim(COMP_NAME)//":DYNAMICS"
      FRIENDLIES%QLLS     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QLCN     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QILS     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QICN     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QRAIN    = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QSNOW    = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QGRAUPEL = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"

    !BOS

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'Q',                                         &
         LONG_NAME  = 'specific_humidity',                         &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES%QV),                         &
         default    = 1.0e-6,                                      &
         RESTART    = MAPL_RestartRequired,                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                        &
         SHORT_NAME = 'QLLS',                                            &
         LONG_NAME  = 'mass_fraction_of_large_scale_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                         &
         FRIENDLYTO = trim(FRIENDLIES%QLLS),                             &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                   RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                       &
         SHORT_NAME = 'QLCN',                                           &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                        &
         FRIENDLYTO = trim(FRIENDLIES%QLCN),                            &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'CLLS',                                      &
         LONG_NAME  = 'large_scale_cloud_area_fraction',           &
         UNITS      = '1',                                         &
         FRIENDLYTO = trim(FRIENDLIES%CLLS),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'CLCN',                                      &
         LONG_NAME  = 'convective_cloud_area_fraction',            &
         UNITS      = '1',                                         &
         FRIENDLYTO = trim(FRIENDLIES%CLCN),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                     &
         SHORT_NAME = 'QILS',                                         &
         LONG_NAME  = 'mass_fraction_of_large_scale_cloud_ice_water', &
         UNITS      = 'kg kg-1',                                      &
         FRIENDLYTO = trim(FRIENDLIES%QILS),                          &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                    &
         SHORT_NAME = 'QICN',                                        &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_ice_water', &
         UNITS      = 'kg kg-1',                                     &
         FRIENDLYTO = trim(FRIENDLIES%QICN),                         &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QRAIN',                                     &
         LONG_NAME  = 'mass_fraction_of_rain',                     & 
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES%QRAIN),                      &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QSNOW',                                     &
         LONG_NAME  = 'mass_fraction_of_snow',                     &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES%QSNOW),                      &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QGRAUPEL',                                  &
         LONG_NAME  = 'mass_fraction_of_graupel',                  &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES%QGRAUPEL),                   &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QW',                                        &
         LONG_NAME  = 'mass_fraction_of_wet_air',                  &
         UNITS      = 'kg kg-1',                                   &
         RESTART    = MAPL_RestartSkip,                            &
         FRIENDLYTO = trim(FRIENDLIES%QW),                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                      

    call MAPL_TimerAdd(GC, name="--GFDL_1M", RC=STATUS)
    VERIFY_(STATUS)

end subroutine GFDL_1M_Setup

subroutine GFDL_1M_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

    type (ESMF_Grid )                   :: GRID
    type (ESMF_State),         pointer  :: GIM(:)
    type (ESMF_State),         pointer  :: GEX(:)
    type (ESMF_State)                   :: INTERNAL

    real, pointer, dimension(:,:,:)     :: Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL, QW

    call MAPL_GetResource( MAPL, LHYDROSTATIC, Label="HYDROSTATIC:",  default=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LPHYS_HYDROSTATIC, Label="PHYS_HYDROSTATIC:",  default=.TRUE., RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_Get ( MAPL, LM=LM, GIM=GIM, GEX=GEX, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QRAIN,    'QRAIN'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QSNOW,    'QSNOW'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QW,       'QW'      , RC=STATUS); VERIFY_(STATUS)
    QW = Q+QLLS+QLCN+QILS+QICN+QRAIN+QSNOW+QGRAUPEL

    call gfdl_cloud_microphys_init()
    call WRITE_PARALLEL ("INITIALIZED GFDL_1M microphysics in non-generic GC INIT")

    call MAPL_GetResource(MAPL,USE_AEROSOL_NN,'USE_AEROSOL_NN:',default=.TRUE., RC=STATUS )
    VERIFY_(STATUS)
    call aer_cloud_init()
    call WRITE_PARALLEL ("INITIALIZED aer_cloud_init for GFDL_1M")

    call MAPL_GetResource (MAPL, JASON_TUNING, trim(COMP_NAME)//"_JASON_TUNING:", default=0, RC=STATUS); VERIFY_(STATUS)

end subroutine GFDL_1M_Initialize

subroutine GFDL_1M_Run (GC, IMPORT, EXPORT, CLOCK, RC)
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

    call MAPL_TimerOn (MAPL,"--GFDL_1M",RC=STATUS)

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

    call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QRAIN,    'QRAIN'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QSNOW,    'QSNOW'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLCN,     'CLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLLS,     'CLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)

        call MAPL_TimerOn(MAPL,"---CLDMACRO")
        TEMP = TH1*PK
         DTDT_macro=TEMP
        DQVDT_macro=Q1
        DQLDT_macro=QLCN+QLLS
        DQIDT_macro=QICN+QILS
        DQADT_macro=CLCN+CLLS
        DQRDT_macro=QRAIN
        DQSDT_macro=QSNOW
        DQST3 = GEOS_DQSAT(TEMP, PLO, QSAT=QST3)
        do K=1,LM
          do J=1,JM
           do I=1,IM
       ! add DeepCu & ShallowCu QL/QI/CL to Convective
             call cnvsrc (            &
                  DT_MOIST          , &
                  CLDPARAMS%CNV_ICEPARAM, &
                  CLDPARAMS%SCLM_DEEP, &
                  CLDPARAMS%SCLM_SHALLOW, &
                  MASS(I,J,K)       , &
                  iMASS(I,J,K)      , &
                  PLO(I,J,K)        , &
                  TEMP(I,J,K)       , &
                  Q1(I,J,K)         , &
                  CNV_DQLDT(I,J,K)  , &   ! <- DeepCu
                  CNV_MFD(I,J,K)    , &   ! <- DeepCu
                  QLDET_SC(I,J,K)   , &   ! <- ShlwCu 
                  QIDET_SC(I,J,K)   , &   ! <- ShlwCu   
!                  MFD_SC(I,J,K)     , &   ! <- ShlwCu   
                  DCM_SC(I,J,K)     , &   ! <- ShlwCu   
                  QLCN(I,J,K)       , &
                  QICN(I,J,K)       , &
                  CLLS(I,J,K)       , &
                  CLCN(I,J,K)       , &
                  QST3(I,J,K)       , &
                  CNV_FRACTION(I,J), SNOMAS(I,J), FRLANDICE(I,J), FRLAND(I,J), &
                  CONVPAR_OPTION )
       ! Send the condensates through the pdf after convection
       !  Use Slingo-Ritter (1985) formulation for critical relative humidity
             TURNRHCRIT   = TURNRHCRIT2D(I,J)
             tmpminrhcrit = minrhcrit2D(I,J)
             tmpmaxrhcrit = maxrhcrit2D(I,J)
             ALPHA = tmpmaxrhcrit
           ! lower turn from maxrhcrit
             if (PLO(i,j,k) .le. TURNRHCRIT) then
                ALPHAl = tmpminrhcrit
             else
                if (k.eq.LM) then
                   ALPHAl = tmpmaxrhcrit
                else
                   ALPHAl = tmpminrhcrit + (tmpmaxrhcrit-tmpminrhcrit)/(19.) * &
                           ((atan( (2.*(PLO(i,j,k)-TURNRHCRIT)/min(100., CNV_PLE(i,j,LM)-TURNRHCRIT)-1.) * &
                           tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
                endif
             endif
           ! upper turn back to maxrhcrit
             TURNRHCRIT_UP = TROPP(i,j)/100.0
             IF (TURNRHCRIT_UP == MAPL_UNDEF) TURNRHCRIT_UP = 100.
             if (PLO(i,j,k) .le. TURNRHCRIT_UP) then
                ALPHAu = tmpmaxrhcrit
             else
                ALPHAu = tmpmaxrhcrit - (tmpmaxrhcrit-tmpminrhcrit)/(19.) * &
                        ((atan( (2.*(PLO(i,j,k)-TURNRHCRIT_UP)/( TURNRHCRIT-TURNRHCRIT_UP)-1.) * &
                        tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
             endif
           ! combine and limit
             ALPHA = min( 0.25, 1.0 - min(max(ALPHAl,ALPHAu),1.) ) ! restrict RHcrit to > 75% 
             RHCRIT = 1.0 - ALPHA
       ! evaporation/sublimation for CN/LS
             EVAPC_X(I,J,K) = Q1(I,J,K)
             call evap3 (                 &
                  DT_MOIST              , &
                  CLDPARAMS%CCW_EVAP_EFF, &
                  RHCRIT                , &
                  PLO(I,J,K)            , &
                  TEMP(I,J,K)           , &
                  Q1(I,J,K)             , &
                  QLCN(I,J,K)           , &
                  QICN(I,J,K)           , &
                  CLCN(I,J,K)           , &
                  CLLS(I,J,K)           , &
                  NACTL(I,J,K)          , &
                  NACTI(I,J,K)          , &
                  QST3(I,J,K)           )
             EVAPC_X(I,J,K) = ( Q1(I,J,K) - EVAPC_X(I,J,K) ) / DT_MOIST
             SUBLC_X(I,J,K) = Q1(I,J,K)
             call subl3 (                 &
                  DT_MOIST              , &
                  CLDPARAMS%CCI_EVAP_EFF, &
                  RHCRIT                , &
                  PLO(I,J,K)            , &
                  TEMP(I,J,K)           , &
                  Q1(I,J,K)             , &
                  QLCN(I,J,K)           , &
                  QICN(I,J,K)           , &
                  CLCN(I,J,K)           , &
                  CLLS(I,J,K)           , &
                  NACTL(I,J,K)          , &
                  NACTI(I,J,K)          , &
                  QST3(I,J,K)           )
             SUBLC_X(I,J,K) = ( Q1(I,J,K) - SUBLC_X(I,J,K) ) / DT_MOIST
       ! cleanup clouds
             call fix_up_clouds( Q1(I,J,K), TEMP(I,J,K), QLLS(I,J,K), QILS(I,J,K), CLLS(I,J,K), QLCN(I,J,K), QICN(I,J,K), CLCN(I,J,K) )
            end do ! IM loop
          end do ! JM loop
        end do ! LM loop
      ! add ShallowCu rain/snow tendencies
        QRAIN = QRAIN + SHLW_PRC3*DT_MOIST
        QSNOW = QSNOW + SHLW_SNO3*DT_MOIST
      ! add DeepCu rain/snow
        if (ADJUSTL(CONVPAR_OPTION) /= 'GF') then
            where (TEMP < MAPL_TICE) !SNOW
              QSNOW = QSNOW + CNV_PRC3*iMASS*DT_MOIST
              TEMP  = TEMP  + CNV_PRC3*iMASS*(MAPL_ALHS-MAPL_ALHL) / MAPL_CP
            else where ! RAIN
              QRAIN = QRAIN + CNV_PRC3*iMASS*DT_MOIST
            end where
        end if
        DTDT_macro=(TEMP-DTDT_macro)/DT_MOIST
        DQVDT_macro=(Q1-DQVDT_macro)/DT_MOIST
        DQLDT_macro=((QLCN+QLLS)-DQLDT_macro)/DT_MOIST
        DQIDT_macro=((QICN+QILS)-DQIDT_macro)/DT_MOIST
        DQADT_macro=((CLCN+CLLS)-DQADT_macro)/DT_MOIST
        DQRDT_macro=(QRAIN-DQRDT_macro)/DT_MOIST
        DQSDT_macro=(QSNOW-DQSDT_macro)/DT_MOIST
        TH1 = TEMP/PK
     ! Zero-out 3D CNV/ANV/SHL CLDMACRO Precipitation & Fluxes
        REV_AN_X = 0.0
        RSU_AN_X = 0.0
        PFI_CN_X = 0.0
        PFL_CN_X = 0.0
        PFI_AN_X = 0.0
        PFL_AN_X = 0.0
        PFI_SC_X = 0.0
        PFL_SC_X = 0.0
        AN_PRC2  = 0.0
        CN_PRC2  = 0.0
        SC_PRC2  = 0.0
        AN_SNR   = 0.0
        CN_SNR   = 0.0
        SC_SNR   = 0.0
        call MAPL_TimerOff(MAPL,"---CLDMACRO")

        call MAPL_TimerOn(MAPL,"---CLDMICRO")
        ! Temperature (K)
         TEMP = TH1*PK
        ! Delta-Z layer thickness (gfdl expects this to be negative)
         DZ = -1.0*DZET
        ! Get cloud nuclei particle numbers
         if (USE_AEROSOL_NN) then
           CFX =100.*PLO*r_air/TEMP !density times conversion factor
           NCPL = NACTL/CFX ! kg-1
           NCPI = NACTI/CFX ! kg-1
         else
           NCPL = 0.
           NCPI = 0.
         endif
        ! Zero-out microphysics tendencies
         if (associated(DQVDT_micro)) DQVDT_micro = Q1
         if (associated(DQLDT_micro)) DQLDT_micro = QLLS + QLCN
         if (associated(DQIDT_micro)) DQIDT_micro = QILS + QICN
         if (associated(DQRDT_micro)) DQRDT_micro = QRAIN
         if (associated(DQSDT_micro)) DQSDT_micro = QSNOW
         if (associated(DQGDT_micro)) DQGDT_micro = QGRAUPEL
         if (associated(DQADT_micro)) DQADT_micro = CLLS + CLCN
         if (associated(DUDT_micro) ) DUDT_micro  = U1
         if (associated(DVDT_micro) ) DVDT_micro  = V1
         if (associated(DTDT_micro) ) DTDT_micro  = TH1*PK
        ! Zero-out local microphysics tendencies
         DQVDTmic = 0.
         DQLDTmic = 0.
         DQRDTmic = 0.
         DQIDTmic = 0.
         DQSDTmic = 0.
         DQGDTmic = 0.
         DQADTmic = 0.
          DUDTmic = 0.
          DVDTmic = 0.
          DTDTmic = 0.
       ! Zero-out 3D Precipitation Fluxes 
        ! Ice
         PFI_LS_X = 0.
        ! Liquid
         PFL_LS_X = 0.
     ! Put condensates in touch with the PDF
        do K=1,LM
          do J=1,JM
            do I=1,IM
       ! Send the condensates through the pdf after convection
       !  Use Slingo-Ritter (1985) formulation for critical relative humidity
             TURNRHCRIT   = TURNRHCRIT2D(I,J)
             tmpminrhcrit = minrhcrit2D(I,J)
             tmpmaxrhcrit = maxrhcrit2D(I,J)
             ALPHA = tmpmaxrhcrit
           ! lower turn from maxrhcrit
             if (PLO(i,j,k) .le. TURNRHCRIT) then
                ALPHAl = tmpminrhcrit
             else
                if (k.eq.LM) then
                   ALPHAl = tmpmaxrhcrit
                else
                   ALPHAl = tmpminrhcrit + (tmpmaxrhcrit-tmpminrhcrit)/(19.) * &
                           ((atan( (2.*(PLO(i,j,k)-TURNRHCRIT)/min(100., CNV_PLE(i,j,LM)-TURNRHCRIT)-1.) * &
                           tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
                endif
             endif
           ! upper turn back to maxrhcrit
             TURNRHCRIT_UP = TROPP(i,j)/100.0
             IF (TURNRHCRIT_UP == MAPL_UNDEF) TURNRHCRIT_UP = 100.
             if (PLO(i,j,k) .le. TURNRHCRIT_UP) then
                ALPHAu = tmpmaxrhcrit
             else
                ALPHAu = tmpmaxrhcrit - (tmpmaxrhcrit-tmpminrhcrit)/(19.) * &
                        ((atan( (2.*(PLO(i,j,k)-TURNRHCRIT_UP)/( TURNRHCRIT-TURNRHCRIT_UP)-1.) * &
                        tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
             endif
           ! combine and limit
             ALPHA = min( 0.25, 1.0 - min(max(ALPHAl,ALPHAu),1.) ) ! restrict RHcrit to > 75% 
             RHCRIT = 1.0 - ALPHA
             IF(USE_AEROSOL_NN) THEN
                   call hystpdf_new( &
                      DT_MOIST    , &
                      ALPHA       , &
                      INT(CLDPARAMS%PDFSHAPE), &
                      PLO(I,J,K)  , &
                      ZLO(I,J,K)  , &
                      Q1(I,J,K)   , &
                      QLLS(I,J,K) , &
                      QLCN(I,J,K) , &
                      QILS(I,J,K) , &
                      QICN(I,J,K) , &
                      TEMP(I,J,K) , &
                      CLLS(I,J,K) , &
                      CLCN(I,J,K) , &
                      NACTL(I,J,K), &
                      NACTI(I,J,K), &
                      WHL(I,J,K),   &
                      WQT(I,J,K),   &
                      HL2(I,J,K),   &
                      QT2(I,J,K),   &
                      HLQT(I,J,K),  &
                      W3(I,J,K),    &
                      W2(I,J,K),    &
                      QT3(I,J,K),   &
                      HL3(I,J,K),   &
                      EDMF_FRC(I,J,K),  &
                      PDF_A(I,J,K), &
#ifdef PDFDIAG
                      PDF_SIGW1(I,J,K), PDF_SIGW2(I,J,K), PDF_W1(I,J,K), PDF_W2(I,J,K), &
                      PDF_SIGTH1(I,J,K), PDF_SIGTH2(I,J,K), PDF_TH1(I,J,K), PDF_TH2(I,J,K), &
                      PDF_SIGQT1(I,J,K), PDF_SIGQT2(I,J,K), PDF_QT1(I,J,K), PDF_QT2(I,J,K), &
                      PDF_RQTTH(I,J,K), PDF_RWTH(I,J,K), PDF_RWQT(I,J,K),            &
#endif
                      WTHV2(I,J,K), WQL(I,J,K), &
                      CNV_FRACTION(I,J), SNOMAS(I,J), FRLANDICE(I,J), FRLAND(I,J))
             ELSE
                   call hystpdf( &
                      DT_MOIST    , &
                      ALPHA       , &
                      INT(CLDPARAMS%PDFSHAPE), &
                      PLO(I,J,K)  , &
                      Q1(I,J,K)   , &
                      QLLS(I,J,K) , &
                      QLCN(I,J,K) , &
                      QILS(I,J,K) , &
                      QICN(I,J,K) , &
                      TEMP(I,J,K) , &
                      CLLS(I,J,K) , &
                      CLCN(I,J,K) , &
                      CNV_FRACTION(I,J), SNOMAS(I,J), FRLANDICE(I,J), FRLAND(I,J))
             ENDIF
             call fix_up_clouds( Q1(I,J,K), TEMP(I,J,K), QLLS(I,J,K), QILS(I,J,K), CLLS(I,J,K), QLCN(I,J,K), QICN(I,J,K), CLCN(I,J,K) )
             RHX_X(I,J,K) = Q1(I,J,K)/GEOS_QSAT( TEMP(I,J,K), PLO(I,J,K) )
            end do ! IM loop
          end do ! JM loop
        end do ! LM loop
        ! Cloud
         RAD_CF = MIN(CLCN+CLLS,1.0)
        ! Liquid
         RAD_QL = QLCN+QLLS
        ! Ice
         RAD_QI = QICN+QILS
        ! VAPOR
         RAD_QV = Q1
        ! RAIN
         RAD_QR = QRAIN
        ! SNOW
         RAD_QS = QSNOW
        ! GRAUPEL
         RAD_QG = QGRAUPEL
        ! Vertical velocity
         W1 = W ! directly imported from FV3 in m/s
      !  W1 = -OMEGA*(1.0+MAPL_VIREPS*RAD_QV) * TEMP / PLO * (MAPL_RDRY/MAPL_GRAV)
     ! Preserve CNV/TOT ratios
        ! Cloud
         FQA  =  MIN(1.0,MAX(CLCN/MAX(RAD_CF,1.e-5),0.0))
        ! Liquid
         FQAl =  MIN(1.0,MAX(QLCN/MAX(RAD_QL,1.E-8),0.0))
        ! Ice
         FQAi =  MIN(1.0,MAX(QICN/MAX(RAD_QI,1.E-8),0.0))
        ! Execute GFDL_1M microphysics
         call gfdl_cloud_microphys_driver( &
                             ! Input water/cloud species and liquid+ice CCN [NACTL+NACTI]
                               RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, (NACTL+NACTI)/1.e6, &
                             ! Output tendencies
                               DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
                               DQSDTmic, DQGDTmic, DQADTmic, DTDTmic, &
                             ! Input fields
                               TEMP, W1, U1, V1, DUDTmic, DVDTmic, DZ, DP, &
                             ! constant inputs
                               AREA, DT_MOIST, FRLAND, CNV_FRACTION, &
                               CLDPARAMS%ANV_ICEFALL, CLDPARAMS%LS_ICEFALL, &
                             ! Output rain re-evaporation and sublimation
                               REV_MC_X, RSU_MC_X, & !EVAPC_X, & 
                             ! Output precipitates
                               PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
                             ! Output mass flux during sedimentation (Pa kg/kg)
                               PFL_LS_X(:,:,1:LM), PFI_LS_X(:,:,1:LM), &
                             ! constant grid/time information
                               LHYDROSTATIC, LPHYS_HYDROSTATIC, &
                               1,IM, 1,JM, 1,LM, 1, LM)
     ! Apply tendencies
         TEMP = TEMP + DTDTmic  * DT_MOIST
         U1   = U1   + DUDTmic  * DT_MOIST
         V1   = V1   + DVDTmic  * DT_MOIST
     ! W1 was updated in gfdl_microphsycis, fill WI tendency export
         if (associated(WI)) WI =  (W1 - W)/DT_MOIST
     ! Apply moist/cloud species tendencies
         RAD_QV = RAD_QV + DQVDTmic * DT_MOIST
         RAD_QL = RAD_QL + DQLDTmic * DT_MOIST
         RAD_QR = RAD_QR + DQRDTmic * DT_MOIST
         RAD_QI = RAD_QI + DQIDTmic * DT_MOIST
         RAD_QS = RAD_QS + DQSDTmic * DT_MOIST
         RAD_QG = RAD_QG + DQGDTmic * DT_MOIST
         RAD_CF = RAD_CF + DQADTmic * DT_MOIST
     ! Redistribute CN/LS CF/QL/QI
         call REDISTRIBUTE_CLOUDS(RAD_CF, RAD_QL, RAD_QI, CLCN, CLLS, QLCN, QLLS, QICN, QILS, RAD_QV, TEMP)
     ! Get new CNV/TOT ratios
        ! Cloud
         FQA  =  MIN(1.0,MAX(CLCN/MAX(RAD_CF,1.e-5),0.0))
        ! Liquid
         FQAl =  MIN(1.0,MAX(QLCN/MAX(RAD_QL,1.E-8),0.0))
        ! Ice
         FQAi =  MIN(1.0,MAX(QICN/MAX(RAD_QI,1.E-8),0.0))
     ! Cloud liquid & Ice tendencies (these exports are confusing, for now keep them zeros)
         REV_LS_X = REV_MC_X
         RSU_LS_X = RSU_MC_X
     ! Convert precip diagnostics from mm/day to kg m-2 s-1
         PRCP_RAIN    = PRCP_RAIN    / 86400.
         PRCP_SNOW    = PRCP_SNOW    / 86400.
         PRCP_ICE     = PRCP_ICE     / 86400.
         PRCP_GRAUPEL = PRCP_GRAUPEL / 86400.
         where (PRCP_RAIN .lt. 0.0)
            PRCP_RAIN = 0.0
         end where
         where (PRCP_SNOW .lt. 0.0)
            PRCP_SNOW = 0.0
         end where
         where (PRCP_ICE .lt. 0.0)
            PRCP_ICE = 0.0
         end where
         where (PRCP_GRAUPEL .lt. 0.0)
            PRCP_GRAUPEL = 0.0
         end where
     ! Convert precipitation fluxes from (Pa kg/kg) to (kg m-2 s-1)
         PFL_LS_X = PFL_LS_X/(MAPL_GRAV*DT_MOIST)
         PFI_LS_X = PFI_LS_X/(MAPL_GRAV*DT_MOIST)
     ! Redistribute precipitation fluxes for chemistry
         PFL_AN_X(:,:,1:LM) = PFL_LS_X(:,:,1:LM) * FQAl
         PFL_LS_X(:,:,1:LM) = PFL_LS_X(:,:,1:LM) - PFL_AN_X(:,:,1:LM)
         PFI_AN_X(:,:,1:LM) = PFI_LS_X(:,:,1:LM) * FQAi
         PFI_LS_X(:,:,1:LM) = PFI_LS_X(:,:,1:LM) - PFI_AN_X(:,:,1:LM)
     ! Fill precip diagnostics
         LS_PRC2 = PRCP_RAIN
         LS_SNR  = PRCP_SNOW + PRCP_ICE + PRCP_GRAUPEL
     ! cleanup suspended precipitation condensates
         call FIX_NEGATIVE_PRECIP(RAD_QR, RAD_QS, RAD_QG, RAD_QV, TEMP)
     ! Fill vapor/rain/snow/graupel state
         Q1       = RAD_QV
         QRAIN    = RAD_QR
         QSNOW    = RAD_QS
         QGRAUPEL = RAD_QG
     ! Fill rain/snow exports
         if (associated(QRTOT)) QRTOT = QRAIN
         if (associated(QSTOT)) QSTOT = QSNOW
         if (associated(QGTOT)) QGTOT = QGRAUPEL
     ! Convert back to PT
         TH1 = TEMP/PK
     ! Radiation Coupling
      if (CLDPARAMS%DISABLE_RAD==1) then
         RAD_QL     = 0.
         RAD_QI     = 0.
         RAD_QR     = 0.
         RAD_QS     = 0.
         RAD_QG     = 0.
         RAD_CF     = 0.
         CLDREFFL   = 0.
         CLDREFFI   = 0.
      else
         do K = 1, LM
           do J = 1, JM
             do I = 1, IM
              ! cleanup clouds
               call fix_up_clouds( Q1(I,J,K), TEMP(I,J,K), QLLS(I,J,K), QILS(I,J,K), CLLS(I,J,K), QLCN(I,J,K), QICN(I,J,K), CLCN(I,J,K) )
              ! get radiative properties
               call RADCOUPLE ( TEMP(I,J,K), PLO(I,J,K), CLLS(I,J,K), CLCN(I,J,K), &
                     Q1(I,J,K), QLLS(I,J,K), QILS(I,J,K), QLCN(I,J,K), QICN(I,J,K), QRAIN(I,J,K), QSNOW(I,J,K), QGRAUPEL(I,J,K), NACTL(I,J,K), NACTI(I,J,K), &
                     RAD_QV(I,J,K), RAD_QL(I,J,K), RAD_QI(I,J,K), RAD_QR(I,J,K), RAD_QS(I,J,K), RAD_QG(I,J,K), RAD_CF(I,J,K), &
                     CLDREFFL(I,J,K), CLDREFFI(I,J,K), FRLAND(I,J), CNV_FRACTION(I,J), INT(CLDPARAMS%FR_AN_WAT), & 
                     CLDPARAMS%FAC_RL, CLDPARAMS%MIN_RL, CLDPARAMS%MAX_RL, CLDPARAMS%FAC_RI, CLDPARAMS%MIN_RI, CLDPARAMS%MAX_RI, &
                     CLDPARAMS%CCN_OCEAN, CLDPARAMS%CCN_LAND)
            enddo
          enddo
        enddo
      ! Clean up radiation clouds after microphysics
        CALL FIX_UP_RAD_CLOUDS( RAD_QV, RAD_QL, RAD_QI, RAD_CF, RAD_QR, RAD_QS, RAD_QG)
      endif

         if (USE_AEROSOL_NN) then
           CFX =100.*PLO*r_air/TEMP !density times conversion factor
           NCPL = NACTL/CFX ! kg-1
           NCPI = NACTI/CFX ! kg-1
         else
           NCPL = 0.
           NCPI = 0.
         endif
       ! Exports required
         CFLIQ  = 0.0
         CFICE  = 0.0
         QTOT   = RAD_QL+RAD_QI
         QL_TOT = RAD_QL
         QI_TOT = RAD_QI
         WHERE (QTOT .gt. 1.e-8)
            CFLIQ=RAD_CF*QL_TOT/QTOT
            CFICE=RAD_CF*QI_TOT/QTOT
         END WHERE
         CFLIQ=MAX(MIN(CFLIQ, 1.0), 0.0)
         CFICE=MAX(MIN(CFICE, 1.0), 0.0)
         where (QI_TOT .le. 1.e-8)
            CFICE =0.0
            NCPI  =0.0
            CLDREFFI = MAPL_UNDEF
         end where
         where (QL_TOT .le. 1.e-8)
            CFLIQ =0.0
            NCPL  =0.0
            CLDREFFL = MAPL_UNDEF
         end where
         if (associated(DQVDT_micro)) DQVDT_micro = (Q1 - DQVDT_micro         ) / DT_MOIST
         if (associated(DQLDT_micro)) DQLDT_micro = ((QLLS+QLCN) - DQLDT_micro) / DT_MOIST
         if (associated(DQIDT_micro)) DQIDT_micro = ((QILS+QICN) - DQIDT_micro) / DT_MOIST
         if (associated(DQRDT_micro)) DQRDT_micro = (QRAIN - DQRDT_micro      ) / DT_MOIST
         if (associated(DQSDT_micro)) DQSDT_micro = (QSNOW - DQSDT_micro      ) / DT_MOIST
         if (associated(DQGDT_micro)) DQGDT_micro = (QGRAUPEL - DQGDT_micro   ) / DT_MOIST
         if (associated(DQADT_micro)) DQADT_micro = ((CLLS+CLCN) - DQADT_micro) / DT_MOIST
         if (associated(DUDT_micro) ) DUDT_micro  = (U1 - DUDT_micro - U1     ) / DT_MOIST
         if (associated(DVDT_micro) ) DVDT_micro  = (V1 - DVDT_micro - V1     ) / DT_MOIST
         if (associated(DTDT_micro) ) DTDT_micro  = (TH1*PK - DTDT_micro      ) / DT_MOIST
        call MAPL_TimerOff(MAPL,"---CLDMICRO")

     call MAPL_TimerOff(MAPL,"--GFDL_1M",RC=STATUS)

end subroutine GFDL_1M_Run

end module GEOS_GFDL_1M_InterfaceMod
