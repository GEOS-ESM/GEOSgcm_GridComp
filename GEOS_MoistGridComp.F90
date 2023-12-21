
! VERIFY_ and RETURN_ macros for error handling.

!#define UWDIAG 1
!#define PDFDIAG 1

! #include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_Moist -- A Module to compute moist processes, including convection,
!   large-scale condensation and precipitation and cloud parameters.

! !INTERFACE:

module GEOS_MoistGridCompMod

  ! !USES:

!   use ESMF
!   use MAPL
  use GEOS_GFDL_1M_InterfaceMod
!   use GEOS_BACM_1M_InterfaceMod
!   use GEOS_MGB2_2M_InterfaceMod
!   use GEOS_RAS_InterfaceMod
  use GEOS_GF_InterfaceMod
  use GEOS_UW_InterfaceMod

  use aer_cloud
  use Aer_Actv_Single_Moment
!   use Lightning_mod, only: HEMCO_FlashRate
  use GEOSmoist_Process_Library
  use GEOS_UtilsMod

  implicit none

  character(LEN=ESMF_MAXSTR):: CONVPAR_OPTION  ! GF, RAS, NONE
  character(LEN=ESMF_MAXSTR):: SHALLOW_OPTION  ! UW, NONE
  character(LEN=ESMF_MAXSTR):: CLDMICR_OPTION  ! BACM_1M, GFDL_1M, MGB2_2M

  private

  logical :: DEBUG = .false.
  logical :: LDIAGNOSE_PRECIP_TYPE
  logical :: LUPDATE_PRECIP_TYPE
  logical :: LHYDROSTATIC
  logical :: USE_AERO_BUFFER
  real    :: CCN_OCN
  real    :: CCN_LND

  ! !PUBLIC MEMBER FUNCTIONS:

!   public SetServices

  ! !DESCRIPTION:
  !
  !   {\tt GEOS\_MoistGridCompMod} implements moist processes in GEOS-5. These
  !   include all processes that involve phase changes in the atmosphere, such
  !   as large-scale condensation, convective clouds, and all rain and cloud
  !   formation. Its state consists of water vapor, various types of condensate,
  !   and fractions of various cloud types.
  !   two moment cloud microphysics (Barahona et al., GMD, 2014.) can be run by setting CLDMACRO==MGB2_2M.
  !   When using 2-moment microphysics the number concentration of ice crystals and cloud droplets
  !   are part of the state.
  !

  !EOP

contains

  ! !IROUTINE: RUN -- Run method for the CONVECT component

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

    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    type( ESMF_VM )                 :: VMG

    ! Local derived type aliases

    type (MAPL_MetaComp), pointer   :: MAPL
    type (ESMF_Config  )            :: CF
    type (ESMF_State   )            :: INTERNAL
    type (ESMF_Alarm   )            :: ALARM
    type (ESMF_TimeInterval)        :: TINT
    real(ESMF_KIND_R8)              :: DT_R8
    real                            :: DT_MOIST

    ! Local variables
    real                                :: Tmax
    real, allocatable, dimension(:,:,:) :: PLEmb, PKE, ZLE0, PK, MASS
    real, allocatable, dimension(:,:,:) :: PLmb,  ZL0, DZET
    real, allocatable, dimension(:,:,:) :: QST3, DQST3, MWFA
    real, allocatable, dimension(:,:,:) :: TMP3D
    real, allocatable, dimension(:,:)   :: TMP2D
    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN
    real, pointer, dimension(:,:,:) :: NACTL, NACTI
    ! Imports
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, W
    real, pointer, dimension(:,:)   :: FRLAND, FRLANDICE, FRACI, SNOMAS
    real, pointer, dimension(:,:)   :: SH, EVAP, KPBL
    real, pointer, dimension(:,:,:) :: TKE, OMEGA
    type(ESMF_State)                :: AERO
    type(ESMF_FieldBundle)          :: TR
    ! Exports
    real, pointer, dimension(:,:,:) :: DQDT, DQADT, DQIDT, DQLDT, DQRDT, DQSDT, DQGDT
    real, pointer, dimension(:,:,:) :: DTDT, DUDT,  DVDT,  DWDT
    real, pointer, dimension(:,:,:) :: DPDTMST, PFL_LSAN, PFI_LSAN
    real, pointer, dimension(:,:,:) :: DTDT_ER, DQVDT_ER
    real, pointer, dimension(:,:  ) :: PTYPE, TPREC, CN_PRCP, LS_PRCP, AN_PRCP, SC_PRCP, PLS, PCU
    real, pointer, dimension(:,:  ) :: RAIN, SNOW, ICE, FRZR, PREC_STRAT, PREC_CONV
    real, pointer, dimension(:,:,:) :: BYNCY
    real, pointer, dimension(:,:  ) :: CAPE, INHB, MLCAPE, SBCAPE, MLCIN, MUCAPE, MUCIN, SBCIN, LFC, LNB
    real, pointer, dimension(:,:  ) :: CNV_FRC, SRF_TYPE
    real, pointer, dimension(:,:,:) :: CFICE, CFLIQ
    real, pointer, dimension(:,:,:  ) :: NWFA
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D



    integer :: IM,JM,LM
    integer :: I, J, L

    !=============================================================================

    ! Begin...

    ! Get my name and set-up traceback handle
    ! ---------------------------------------

    Iam = 'Run'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, VM=VMG, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS) ; VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"TOTAL")

    ! If its time, call run methods
    ! --------------------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,   &
         RUNALARM = ALARM,             &
         INTERNAL_ESMF_STATE=INTERNAL, &
         RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_AlarmGet(ALARM, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)
    DT_MOIST = DT_R8

    if ( ESMF_AlarmIsRinging( ALARM, RC=STATUS) ) then

       call ESMF_AlarmRingerOff(ALARM, RC=STATUS) ; VERIFY_(STATUS)

       ! Internal State
       call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, CLCN,     'CLCN'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, CLLS,     'CLLS'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, NACTL,   'NACTL'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, NACTI,   'NACTI'    , RC=STATUS); VERIFY_(STATUS)

       ! Import State
       call MAPL_GetPointer(IMPORT, PLE,     'PLE'     , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, ZLE,     'ZLE'     , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, T,       'T'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, U,       'U'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, V,       'V'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, W,       'W'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, KPBL,    'KPBL'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, SH,      'SH'      , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, EVAP,    'EVAP'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, OMEGA,   'OMEGA'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, TKE,     'TKE'     ,RC=STATUS); VERIFY_(STATUS)
       call   ESMF_StateGet(IMPORT,'AERO',    AERO     , RC=STATUS); VERIFY_(STATUS)
       call   ESMF_StateGet(IMPORT,'MTR',     TR       , RC=STATUS); VERIFY_(STATUS)

       ! Update SRF_TYPE for ice_fraction
       call MAPL_GetPointer(IMPORT, FRLAND,    'FRLAND'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, FRLANDICE, 'FRLANDICE' , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, FRACI,     'FRACI'     , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, SNOMAS,    'SNOMAS'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, SRF_TYPE,  'SRF_TYPE'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       SRF_TYPE = 0.0 ! Ocean
       where (FRLAND > 0.1)
         SRF_TYPE = 1.0 ! Land
       end where
       where ( (SNOMAS > 0.1) .OR. (FRLANDICE > 0.5) .OR. (FRACI > 0.5) )
         SRF_TYPE = 2.0 ! Ice/Snow
       end where

       ! Allocatables
        ! Edge variables
       ALLOCATE ( ZLE0 (IM,JM,0:LM) )
       ALLOCATE ( PLEmb(IM,JM,0:LM) )
       ALLOCATE ( PKE  (IM,JM,0:LM) )
        ! Layer variables
       ALLOCATE ( ZL0  (IM,JM,LM  ) )
       ALLOCATE ( DZET (IM,JM,LM  ) )
       ALLOCATE ( PLmb (IM,JM,LM  ) )
       ALLOCATE ( PK   (IM,JM,LM  ) )
       ALLOCATE ( DQST3(IM,JM,LM  ) )
       ALLOCATE (  QST3(IM,JM,LM  ) )
       ALLOCATE ( MASS (IM,JM,LM  ) )
       ALLOCATE ( TMP3D(IM,JM,LM  ) )
       ALLOCATE ( TMP2D(IM,JM     ) )

       ! Save input winds
       call MAPL_GetPointer(EXPORT, PTR3D, 'UMST0', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       PTR3D = U
       call MAPL_GetPointer(EXPORT, PTR3D, 'VMST0', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       PTR3D = V

       ! Derived States
       MASS     = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )/MAPL_GRAV
       call FILLQ2ZERO(Q, MASS, TMP2D)
       call MAPL_GetPointer(EXPORT, PTR2D, 'FILLNQV_IN', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = TMP2D
       PLEmb    =  PLE*.01
       PKE      = (PLE/MAPL_P00)**(MAPL_KAPPA)
       PLmb     = 0.5*(PLEmb(:,:,0:LM-1) + PLEmb(:,:,1:LM))
       PK       = (100.0*PLmb/MAPL_P00)**(MAPL_KAPPA)
       DO L=0,LM
          ZLE0(:,:,L)= ZLE(:,:,L) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
       END DO
       ZL0      = 0.5*(ZLE0(:,:,0:LM-1) + ZLE0(:,:,1:LM) ) ! Layer Height (m) above the surface
       DZET     =     (ZLE0(:,:,0:LM-1) - ZLE0(:,:,1:LM) ) ! Layer thickness (m)
       DQST3    = GEOS_DQSAT(T, PLmb, QSAT=QST3)

       ! These may be used by children
       call MAPL_GetPointer(EXPORT, NWFA,    'NWFA'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, CNV_FRC, 'CNV_FRC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, BYNCY,   'BYNCY'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, CAPE,    'CAPE'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, INHB,    'INHB'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, SBCAPE,  'SBCAPE' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, SBCIN,   'SBCIN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, MLCAPE,  'MLCAPE' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, MLCIN,   'MLCIN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, MUCAPE,  'MUCAPE' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, MUCIN,   'MUCIN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, LFC,     'ZLFC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, LNB,     'ZLNB'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call BUOYANCY2( IM, JM, LM, T, Q, QST3, DQST3, DZET, ZL0, PLmb, PLEmb(:,:,LM), SBCAPE, MLCAPE, MUCAPE, SBCIN, MLCIN, MUCIN, BYNCY, LFC, LNB )
       call BUOYANCY( T, Q, QST3, DQST3, DZET, ZL0, BYNCY, CAPE, INHB)

       CNV_FRC = 0.0
       if( CNV_FRACTION_MAX > CNV_FRACTION_MIN ) then
         WHERE (CAPE .ne. MAPL_UNDEF)
            CNV_FRC =(MAX(1.e-6,MIN(1.0,(CAPE-CNV_FRACTION_MIN)/(CNV_FRACTION_MAX-CNV_FRACTION_MIN))))
         END WHERE
       endif
       if (CNV_FRACTION_EXP /= 1.0) then
          CNV_FRC = CNV_FRC**CNV_FRACTION_EXP
       endif

       ! Extract convective tracers from the TR bundle
       call MAPL_TimerOn (MAPL,"---CONV_TRACERS")
       call CNV_Tracers_Init(TR, RC)
       call MAPL_TimerOff(MAPL,"---CONV_TRACERS")

       ! Get aerosol activation properties
       call MAPL_TimerOn (MAPL,"---AERO_ACTIVATE")
       if (USE_AEROSOL_NN) then
         allocate ( AeroProps(IM,JM,LM) )
         ! get veritical velocity
         if (LHYDROSTATIC) then
           TMP3D = -OMEGA/(MAPL_GRAV*PLmb*100.0/(MAPL_RGAS*T))
         else
           TMP3D = W
         endif
         ! Pressures in Pa
         call Aer_Activation(IM,JM,LM, Q, T, PLmb*100.0, PLE, ZL0, ZLE0, QLCN, QICN, QLLS, QILS, &
                             SH, EVAP, KPBL, TKE, TMP3D, FRLAND, USE_AERO_BUFFER, &
                             AeroProps, AERO, NACTL, NACTI, NWFA, CCN_LND*1.e6, CCN_OCN*1.e6)
       else
         do L=1,LM
           NACTL(:,:,L) = (CCN_LND*FRLAND + CCN_OCN*(1.0-FRLAND))*1.e6 ! #/m^3
           NACTI(:,:,L) = (CCN_LND*FRLAND + CCN_OCN*(1.0-FRLAND))*1.e6 ! #/m^3
         end do
       endif
       call MAPL_GetPointer(EXPORT, PTR3D, 'NCCN_LIQ', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = NACTL*1.e-6
       call MAPL_GetPointer(EXPORT, PTR3D, 'NCCN_ICE', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = NACTI*1.e-6

       call MAPL_TimerOff(MAPL,"---AERO_ACTIVATE")

       if (adjustl(CONVPAR_OPTION)=="RAS"    ) call     RAS_Run(GC, IMPORT, EXPORT, CLOCK, RC=STATUS) ; VERIFY_(STATUS)
       if (adjustl(CONVPAR_OPTION)=="GF"     ) call      GF_Run(GC, IMPORT, EXPORT, CLOCK, RC=STATUS) ; VERIFY_(STATUS)
       if (adjustl(SHALLOW_OPTION)=="UW"     ) call      UW_Run(GC, IMPORT, EXPORT, CLOCK, RC=STATUS) ; VERIFY_(STATUS)
       if (adjustl(CLDMICR_OPTION)=="BACM_1M") call BACM_1M_Run(GC, IMPORT, EXPORT, CLOCK, RC=STATUS) ; VERIFY_(STATUS)
       if (adjustl(CLDMICR_OPTION)=="GFDL_1M") call GFDL_1M_Run(GC, IMPORT, EXPORT, CLOCK, RC=STATUS) ; VERIFY_(STATUS)
       if (adjustl(CLDMICR_OPTION)=="MGB2_2M") call MGB2_2M_Run(GC, IMPORT, EXPORT, CLOCK, RC=STATUS) ; VERIFY_(STATUS)

       ! Exports
         ! Cloud fraction exports


         call MAPL_GetPointer(EXPORT, CFICE, 'CFICE', ALLOC=.true., RC=STATUS); VERIFY_(STATUS)
         if (associated(CFICE)) then
           CFICE=0.0
           WHERE (QILS+QICN .gt. 1.0e-12)
              CFICE=(CLLS+CLCN)*(QILS+QICN)/(QLLS+QLCN+QILS+QICN)
           END WHERE
           CFICE=MAX(MIN(CFICE, 1.0), 0.0)
         endif
         call MAPL_GetPointer(EXPORT, CFLIQ, 'CFLIQ', RC=STATUS); VERIFY_(STATUS)
         if (associated(CFLIQ)) then
           CFLIQ=0.0
           WHERE (QLLS+QLCN .gt. 1.0e-12)
              CFLIQ=(CLLS+CLCN)*(QLLS+QLCN)/(QLLS+QLCN+QILS+QICN)
           END WHERE
           CFLIQ=MAX(MIN(CFLIQ, 1.0), 0.0)
         endif

         ! Rain-out and Relative Humidity where RH > 110%
         call MAPL_GetPointer(EXPORT,  DTDT_ER,  'DTDT_ER', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT, DQVDT_ER, 'DQVDT_ER', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
          DTDT_ER = T
         DQVDT_ER = Q

         ! some diagnostics to export
         if (.FALSE.) then
          QST3  = GEOS_QsatICE (T, PLmb*100.0, DQ=DQST3)
         else
          DQST3 = GEOS_DQSAT   (T, PLmb, QSAT=QST3)      ! this qsat function expects hPa...
         end if
         call MAPL_GetPointer(EXPORT, PTR3D, 'RHICE', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR3D)) then
           PTR3D = Q/QST3
           where (T>MAPL_TICE)
             PTR3D=0.0
           end where
         endif

         call MAPL_GetPointer(EXPORT, PTR3D, 'SAT_RAT', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR3D)) then
           where (CFICE .lt. 0.99 .and. QST3 .gt. 1.0e-20)
            TMP3D = max((Q - QST3*CFICE), 0.0)/(1.0-CFICE)
            PTR3D = min(TMP3D/QST3, 2.0)
           elsewhere
            PTR3D = 1.0
           end where
         endif

         if (.FALSE.) then
          QST3  = GEOS_QsatLQU (T, PLmb*100.0, DQ=DQST3) !clean up only with respect to liquid water
         else
          DQST3 = GEOS_DQSAT   (T, PLmb, QSAT=QST3)      ! this qsat function expects hPa...
         end if
         call MAPL_GetPointer(EXPORT, PTR3D, 'RHLIQ', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR3D)) PTR3D = Q/QST3

         ! rainout excesive RH
         call MAPL_GetPointer(EXPORT, LS_PRCP, 'LS_PRCP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT, PTR2D,   'ER_PRCP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
         where ( Q > 1.1*QST3 )
            TMP3D = (Q - 1.1*QST3)/( 1.0 + 1.1*DQST3*MAPL_ALHL/MAPL_CP )
         elsewhere
            TMP3D = 0.0
         endwhere
         PTR2D = SUM(TMP3D*MASS,3)/DT_MOIST
         LS_PRCP = LS_PRCP + PTR2D
         Q = Q - TMP3D
         T = T + (MAPL_ALHL/MAPL_CP)*TMP3D
          DTDT_ER = (T -  DTDT_ER)/DT_MOIST
         DQVDT_ER = (Q - DQVDT_ER)/DT_MOIST

         ! cleanup any negative QV/QC/CF
         call FILLQ2ZERO(Q, MASS, TMP2D)
         call MAPL_GetPointer(EXPORT, PTR2D, 'FILLNQV', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR2D)) PTR2D = TMP2D/DT_MOIST

       if (USE_AEROSOL_NN) then
         deallocate ( AeroProps )
       endif


       ! Export Total Moist Tendencies

       call MAPL_GetPointer(EXPORT, DUDT, 'DUDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DUDT)) then
          DUDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DUDT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DUDT = DUDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DUDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DUDT = DUDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DUDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DUDT = DUDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DUDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DUDT = DUDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DVDT, 'DVDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DVDT)) then
          DVDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DVDT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DVDT = DVDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DVDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DVDT = DVDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DVDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DVDT = DVDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DVDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DVDT = DVDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DTDT, 'DTDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DTDT)) then
          DTDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DTDT = DTDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DTDT = DTDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DTDT = DTDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DTDT = DTDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_ER'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DTDT = DTDT + PTR3D
          DTDT = DTDT*(PLE(:,:,1:LM)-PLE(:,:,0:LM-1)) ! Pressure weighted tendency
       endif

       call MAPL_GetPointer(EXPORT, DQDT, 'DQDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DQDT)) then
          DQDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQDT = DQDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQDT = DQDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQDT = DQDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQDT = DQDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_ER'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQDT = DQDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DQLDT, 'DQLDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DQLDT)) then
          DQLDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DQIDT, 'DQIDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DQIDT)) then
          DQIDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DQRDT, 'DQRDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DQRDT)) then
          DQRDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQRDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQRDT = DQRDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQRDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQRDT = DQRDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQRDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQRDT = DQRDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DQSDT, 'DQSDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DQSDT)) then
          DQSDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQSDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQSDT = DQSDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQSDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQSDT = DQSDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQSDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQSDT = DQSDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DQGDT, 'DQGDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DQGDT)) then
          DQGDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQGDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQGDT = DQGDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQGDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQGDT = DQGDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DQADT, 'DQADT'  , RC=STATUS); VERIFY_(STATUS)
       if (associated(DQADT)) then
          DQADT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQADT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQADT = DQADT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQADT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQADT = DQADT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQADT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQADT = DQADT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQADT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQADT = DQADT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DPDTMST, 'DPDTMST'  , RC=STATUS); VERIFY_(STATUS)
       if (associated(DPDTMST)) then
          DPDTMST = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_CN'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_AN'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_LS'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_CN'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_AN'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_LS'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          DPDTMST = MAPL_GRAV * DPDTMST
        endif

       call MAPL_GetPointer(EXPORT, PFL_LSAN, 'PFL_LSAN'  , RC=STATUS); VERIFY_(STATUS)
       if (associated(PFL_LSAN)) then
          PFL_LSAN = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_AN'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) PFL_LSAN = PFL_LSAN + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_LS'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) PFL_LSAN = PFL_LSAN + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, PFI_LSAN, 'PFI_LSAN'  , RC=STATUS); VERIFY_(STATUS)
       if (associated(PFI_LSAN)) then
          PFI_LSAN = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_AN'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) PFI_LSAN = PFI_LSAN + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_LS'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) PFI_LSAN = PFI_LSAN + PTR3D
       endif

       ! Combine Precip Exports

       ! liquid convective precip
       call MAPL_GetPointer(EXPORT, PCU, 'PCU', RC=STATUS); VERIFY_(STATUS)
       if (associated(PCU)) then
          PCU = 0.0
          call MAPL_GetPointer(EXPORT, PTR2D, 'CN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) PCU = PCU + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'SC_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) PCU = PCU + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'CNPCPRATE' , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) PCU = PCU + PTR2D
          PCU = MAX(PCU, 0.0)
       endif

       ! liquid large-scale precip
       call MAPL_GetPointer(EXPORT, PLS, 'PLS', RC=STATUS); VERIFY_(STATUS)
       if (associated(PLS)) then
          PLS = 0.0
          call MAPL_GetPointer(EXPORT, PTR2D, 'LS_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) PLS = PLS + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'AN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) PLS = PLS + PTR2D
          PLS = MAX(PLS, 0.0)
       endif

       ! all liquid precip
       call MAPL_GetPointer(EXPORT, RAIN, 'RAIN', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(RAIN)) then
          RAIN = 0.0
          call MAPL_GetPointer(EXPORT, PTR2D, 'LS_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) RAIN = RAIN + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'AN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) RAIN = RAIN + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'CN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) RAIN = RAIN + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'SC_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) RAIN = RAIN + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'CNPCPRATE' , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) RAIN = RAIN + PTR2D
          RAIN = MAX(RAIN, 0.0)
       endif

       ! all frozen precip (snow at this point)
       call MAPL_GetPointer(EXPORT, SNOW, 'SNO', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(SNOW)) then
          SNOW = 0.0
          call MAPL_GetPointer(EXPORT, PTR2D, 'LS_SNR'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) SNOW = SNOW + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'AN_SNR'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) SNOW = SNOW + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'CN_SNR'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) SNOW = SNOW + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'SC_SNR'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) SNOW = SNOW + PTR2D
          SNOW = MAX(SNOW, 0.0)
       endif

       ! all deep convective precip (rain+snow)
       call MAPL_GetPointer(EXPORT, CN_PRCP, 'CN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
       if (associated(CN_PRCP)) then
          call MAPL_GetPointer(EXPORT, PTR2D, 'CNPCPRATE' , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) CN_PRCP = CN_PRCP + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'CN_SNR'    , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) CN_PRCP = CN_PRCP + PTR2D
       endif

       ! all large-scale precip (rain+snow)
       call MAPL_GetPointer(EXPORT, LS_PRCP, 'LS_PRCP'   , RC=STATUS); VERIFY_(STATUS)
       if (associated(LS_PRCP)) then
          call MAPL_GetPointer(EXPORT, PTR2D, 'LS_SNR'    , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) LS_PRCP = LS_PRCP + PTR2D
       endif

       ! all anvil precip (rain+snow)
       call MAPL_GetPointer(EXPORT, AN_PRCP, 'AN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
       if (associated(AN_PRCP)) then
          call MAPL_GetPointer(EXPORT, PTR2D, 'AN_SNR'    , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) AN_PRCP = AN_PRCP + PTR2D
       endif

       ! all shallow precip (rain+snow)
       call MAPL_GetPointer(EXPORT, SC_PRCP, 'SC_PRCP'   , RC=STATUS); VERIFY_(STATUS)
       if (associated(SC_PRCP)) then
          call MAPL_GetPointer(EXPORT, PTR2D, 'SC_SNR'    , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) SC_PRCP = SC_PRCP + PTR2D
       endif

       ! Total - all precip (rain+snow)
       call MAPL_GetPointer(EXPORT, TPREC, 'TPREC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(TPREC)) then
          TPREC = 0.0
          call MAPL_GetPointer(EXPORT, PTR2D, 'LS_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) TPREC = TPREC + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'AN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) TPREC = TPREC + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'CN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) TPREC = TPREC + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'SC_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) TPREC = TPREC + PTR2D
          TPREC = MAX(TPREC, 0.0)
       endif

       ! diagnosed stratiform precip (rain+snow)
       call MAPL_GetPointer(EXPORT, PREC_STRAT, 'PREC_STRAT', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PREC_STRAT)) then
          if( CNV_FRACTION_MAX > CNV_FRACTION_MIN ) then
             PREC_STRAT = (1.0-CNV_FRC)*TPREC
          else
             PREC_STRAT = 0.0
             call MAPL_GetPointer(EXPORT, PTR2D, 'LS_PRCP'   , RC=STATUS); VERIFY_(STATUS)
             if (associated(PTR2D)) PREC_STRAT = PREC_STRAT + PTR2D
             call MAPL_GetPointer(EXPORT, PTR2D, 'AN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
             if (associated(PTR2D)) PREC_STRAT = PREC_STRAT + PTR2D
             PREC_STRAT = MAX(PREC_STRAT, 0.0)
          endif
       endif

       ! diagnosed convective precip (rain+snow)
       call MAPL_GetPointer(EXPORT, PREC_CONV, 'PREC_CONV', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PREC_CONV)) then
          if( CNV_FRACTION_MAX > CNV_FRACTION_MIN ) then
             PREC_CONV = CNV_FRC*TPREC
          else
             PREC_CONV = 0.0
             call MAPL_GetPointer(EXPORT, PTR2D, 'CN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
             if (associated(PTR2D)) PREC_CONV = PREC_CONV + PTR2D
             call MAPL_GetPointer(EXPORT, PTR2D, 'SC_PRCP'   , RC=STATUS); VERIFY_(STATUS)
             if (associated(PTR2D)) PREC_CONV = PREC_CONV + PTR2D
             PREC_CONV = MAX(PREC_CONV, 0.0)
          endif
       endif

     ! Diagnostic precip types:
       call MAPL_GetPointer(EXPORT, ICE,   'ICE',   ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, FRZR,  'FRZR',  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (LUPDATE_PRECIP_TYPE .OR. LDIAGNOSE_PRECIP_TYPE) then
          call MAPL_GetPointer(EXPORT, PTYPE, 'PTYPE', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
          call DIAGNOSE_PRECIP_TYPE(IM, JM, LM, TPREC, PLS, PCU, RAIN, SNOW, ICE, FRZR, &
                                    PTYPE, PLE, T/PK, PK, PKE, ZL0, LUPDATE_PRECIP_TYPE)
       endif
     ! Get Kuchera snow:rain ratios
       do I = 1,IM
          do J = 1,JM
              Tmax = 0.0
              do L =  LM, 1, -1
                 if (PLmb(I,J,L).gt.500.) then
                    Tmax = MAX(Tmax,T(I,J,L))
                 end if
              end do
              if (Tmax <= 271.16) then
                 TMP2D(I,J) = 12.0 + (271.16 - Tmax)
              else
                 TMP2D(I,J) = 12.0 + 2*(271.16 - Tmax)
              end if
              TMP2D(I,J) = max(0.0,TMP2D(I,J))
          end do
       end do
       call MAPL_GetPointer(EXPORT, PTR2D,'KUCHERA_RATIO', RC=STATUS); VERIFY_(STATUS)
       if(associated(PTR2D)) PTR2D=TMP2D
     ! Accumulated precip totals (mm), apply KUCHERA_RATIO for SNOW
       call MAPL_GetPointer(EXPORT, PTR2D,'SNOWTOTAL', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = TMP2D*DT_MOIST*(SNOW+ICE)
       call MAPL_GetPointer(EXPORT, PTR2D,'PRECTOTAL', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = DT_MOIST*TPREC

       call MAPL_GetPointer(EXPORT, PTR3D, 'QLTOT', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = QLLS+QLCN

       call MAPL_GetPointer(EXPORT, PTR3D, 'QITOT', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = QILS+QICN

       call MAPL_GetPointer(EXPORT, PTR3D, 'QCTOT', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = MIN(CLLS+CLCN,1.0)

       ! Cloud condensate exports
       call MAPL_GetPointer(EXPORT, PTR3D, 'QLLSX1', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = QLLS

       call MAPL_GetPointer(EXPORT, PTR3D, 'QILSX1', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = QILS

       call MAPL_GetPointer(EXPORT, PTR3D, 'QLCNX1', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = QLCN

       call MAPL_GetPointer(EXPORT, PTR3D, 'QICNX1', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = QICN

       ! Fill wind, temperature & RH exports needed for SYNCTQ

       call MAPL_GetPointer(EXPORT, PTR3D, 'UAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = U

       call MAPL_GetPointer(EXPORT, PTR3D, 'VAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = V

       call MAPL_GetPointer(EXPORT, PTR3D, 'TAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = T

       call MAPL_GetPointer(EXPORT, PTR3D, 'QAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = Q

       call MAPL_GetPointer(EXPORT, PTR3D, 'THAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = T/PK

       call MAPL_GetPointer(EXPORT, PTR3D, 'SAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) then
          do L=1,LM
            PTR3D(:,:,L) = MAPL_CP*T(:,:,L) + MAPL_GRAV*(ZL0(:,:,L)+ZLE(:,:,LM))
          enddo
       endif

       call MAPL_GetPointer(EXPORT, PTR3D, 'RH2', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = MAX(MIN( Q/GEOS_QSAT (T, PLmb) , 1.02 ),0.0)
       call MAPL_GetPointer(EXPORT, PTR2D, 'CWP', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = SUM( ( QLCN+QLLS+QICN+QILS )*MASS , 3 )
       call MAPL_GetPointer(EXPORT, PTR2D, 'CLWP', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = SUM( ( QLCN+QLLS ) *MASS , 3 )
       call MAPL_GetPointer(EXPORT, PTR2D, 'LWP', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = SUM( ( QLCN+QLLS) *MASS , 3 )
       call MAPL_GetPointer(EXPORT, PTR2D, 'IWP', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = SUM( ( QICN+QILS ) *MASS , 3 )
       call MAPL_GetPointer(EXPORT, PTR2D, 'TPW', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = SUM( ( Q         ) *MASS , 3 )

       ! Lightning Exports
       call MAPL_GetPointer(EXPORT, PTR2D, 'LFR_GCC', NotFoundOk=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = 0.0

    else

       ! Internal State
       call MAPL_GetPointer(INTERNAL, Q,        'Q'    , RC=STATUS); VERIFY_(STATUS)
       ! Import State
       call MAPL_GetPointer(IMPORT, PLE,     'PLE'     , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, ZLE,     'ZLE'     , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, U,       'U'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, V,       'V'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, T,       'T'       , RC=STATUS); VERIFY_(STATUS)
       ! Allocatables
        ! Edge variables
       ALLOCATE ( ZLE0 (IM,JM,0:LM) )
       ALLOCATE ( PLEmb(IM,JM,0:LM) )
        ! Layer variables
       ALLOCATE ( ZL0  (IM,JM,LM  ) )
       ALLOCATE ( PLmb (IM,JM,LM  ) )
       ALLOCATE ( PK   (IM,JM,LM  ) )
       ALLOCATE ( MASS (IM,JM,LM  ) )
       ALLOCATE ( TMP2D(IM,JM     ) )
       ! dervied states
       MASS     = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )/MAPL_GRAV
       call FILLQ2ZERO(Q, MASS, TMP2D)
       call MAPL_GetPointer(EXPORT, PTR2D, 'FILLNQV_IN', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = TMP2D
       PLEmb    = PLE*.01
       PLmb     = 0.5*(PLEmb(:,:,0:LM-1) + PLEmb(:,:,1:LM))
       PK       = (100.0*PLmb/MAPL_P00)**(MAPL_KAPPA)
       DO L=0,LM
          ZLE0(:,:,L)= ZLE(:,:,L) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
       END DO
       ZL0      = 0.5*(ZLE0(:,:,0:LM-1) + ZLE0(:,:,1:LM) ) ! Layer Height (m) above the surface

       ! Fill Wind, Temperature & RH exports needed for SYNCTQ

       call MAPL_GetPointer(EXPORT, PTR3D, 'UAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = U

       call MAPL_GetPointer(EXPORT, PTR3D, 'VAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = V

       call MAPL_GetPointer(EXPORT, PTR3D, 'TAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = T

       call MAPL_GetPointer(EXPORT, PTR3D, 'QAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = Q

       call MAPL_GetPointer(EXPORT, PTR3D, 'THAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = T/PK

       call MAPL_GetPointer(EXPORT, PTR3D, 'SAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) then
          do L=1,LM
            PTR3D(:,:,L) = MAPL_CP*T(:,:,L) + MAPL_GRAV*(ZL0(:,:,L)+ZLE(:,:,LM))
          enddo
       endif

       call MAPL_GetPointer(EXPORT, PTR3D, 'RH2', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = MAX(MIN( Q/GEOS_QSAT (T, PLmb) , 1.02 ),0.0)

    endif

    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  end subroutine RUN

end module GEOS_MoistGridCompMod

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”
  
! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”
  
! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at
  
! http://www.apache.org/licenses/LICENSE-2.0
  
! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.