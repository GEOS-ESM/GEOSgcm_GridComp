! _VERIFY and _RETURN macros for error handling.

!#define UWDIAG 1
!#define PDFDIAG 1

#include "MAPL.h"

!=============================================================================
!BOP

! !MODULE: GEOS_Moist -- A Module to compute moist processes, including convection,
!   large-scale condensation and precipitation and cloud parameters.

! !INTERFACE:

module GEOS_MoistGridCompMod

   ! !USES:

   use ESMF
   use MAPL
   use gFTL2_StringVector, only: StringVector
   use GEOS_BACM_1M_InterfaceMod
   use GEOS_GFDL_1M_InterfaceMod
   use GEOS_THOM_1M_InterfaceMod
   use GEOS_MGB2_2M_InterfaceMod
   use GEOS_RAS_InterfaceMod
   use GEOS_GF_InterfaceMod
   use GEOS_UW_InterfaceMod

   use Aer_Actv_Single_Moment
   use Lightning_mod, only: HEMCO_FlashRate
   use GEOSmoist_Process_Library
   use GEOS_UtilsMod

   implicit none

   character(len=:), allocatable :: CONVPAR_OPTION ! GF, RAS, NONE
   character(len=:), allocatable :: SHALLOW_OPTION ! UW, NONE
   character(len=:), allocatable :: CLDMICR_OPTION ! BACM_1M, GFDL_1M, THOM_1M, MGB2_2M

   private

   logical :: DEBUG_MST
   logical :: LDIAGNOSE_PRECIP_TYPE
   logical :: LUPDATE_PRECIP_TYPE
   real :: CCN_OCN, CIN_OCN
   real :: CCN_LND, CIN_LND
   logical :: MOIST_USE_NCLOUD_CLIM = .false.
   real :: DETRAIN_INACTIVE_CNV
   real :: TAU_DETRAIN_CNV

   ! COND for the PDF_* exports in Moist_StateSpecs.rc (ACG has no #ifdef)
#ifdef PDFDIAG
   logical, parameter :: PDF_DIAG = .true.
#else
   logical, parameter :: PDF_DIAG = .false.
#endif

   ! !PUBLIC MEMBER FUNCTIONS:

   public SetServices

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

   !BOP

   ! !IROUTINE: SetServices -- Sets ESMF services for this component

   ! !INTERFACE:

   subroutine SetServices(GC, RC)

      ! !ARGUMENTS:

      type(ESMF_GridComp) :: GC ! gridded component
      integer, intent(out) :: RC ! return code

      ! !DESCRIPTION:  {\tt GEOS\_MoistGridCompMod} uses the default Initialize and Finalize
      !                services, but registers its own Run method.

      !EOP

      !=============================================================================
      !
      ! ErrLog Variables

      integer :: STATUS

      ! Local derived type aliases

      type(ESMF_Config) :: CF ! only passed to the unported schemes' *_Setup
      type(StringVector) :: ADVECTED

      logical :: LCONVPAR
      logical :: LSHALLOW
      logical :: LCLDMICR

      !=============================================================================

      ! Begin...

      ! Set the Init entry point
      !-----------------------------------------
      call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_INITIALIZE, Initialize, _RC)

      ! Set the Run entry point
      ! -----------------------
      call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN, Run, phase_name='run', _RC)

      ! Set the Finalize entry point
      ! -----------------------
      call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_FINALIZE, Finalize, _RC)

      ! Inititialize deep convective parameterizations (Options: RAS, GF or NONE)
      !----------------------------------------------------------------------
      call MAPL_GridCompGetResource(GC, 'CONVPAR_OPTION', CONVPAR_OPTION, default="GF", _RC)
      LCONVPAR = adjustl(CONVPAR_OPTION) == "RAS" .or. &
           adjustl(CONVPAR_OPTION) == "GF" .or. &
           adjustl(CONVPAR_OPTION) == "NONE"
      _ASSERT(LCONVPAR, 'Unsupported Deep Convection Option')
      ! default=(X): a copy of the Process_Library default (value and default must not alias)
      call MAPL_GridCompGetResource(GC, 'SIGMA_DX', SIGMA_DX, default=(SIGMA_DX), _RC)
      call MAPL_GridCompGetResource(GC, 'SIGMA_EXP', SIGMA_EXP, default=(SIGMA_EXP), _RC)

      ! Inititialize shallow convective parameterizations (Options: UW or NONE)
      !----------------------------------------------------------------------
      call MAPL_GridCompGetResource(GC, 'SHALLOW_OPTION', SHALLOW_OPTION, default="UW", _RC)
      LSHALLOW = adjustl(SHALLOW_OPTION) == "UW" .or. &
           adjustl(SHALLOW_OPTION) == "NONE"
      _ASSERT(LSHALLOW, 'Unsupported Shallow Convection Option')

      ! Inititialize cloud microphysics (Options: BACM_1M, GFDL_1M, THOM_1M or MGB2_2M)
      !--------------------------------------------------------------
      call MAPL_GridCompGetResource(GC, 'CLDMICR_OPTION', CLDMICR_OPTION, default="BACM_1M", _RC)
      LCLDMICR = adjustl(CLDMICR_OPTION) == "BACM_1M" .or. &
           adjustl(CLDMICR_OPTION) == "GFDL_1M" .or. &
           adjustl(CLDMICR_OPTION) == "THOM_1M" .or. &
           adjustl(CLDMICR_OPTION) == "MGB2_2M"
      _ASSERT(LCLDMICR, 'Unsupported Cloud Microphysics Option')

      call MAPL_GridCompGetResource(GC, 'DEBUG_MST', DEBUG_MST, default=.false., _RC)

      call MAPL_GridCompGetResource(GC, 'DEBUG_TQ_ERRORS', DEBUG_TQ_ERRORS, default=.false., _RC)

      !***********Aerosol-Cloud related

      call MAPL_GridCompGetResource(GC, 'USE_NCLOUD_CLIM', MOIST_USE_NCLOUD_CLIM, default=.false., _RC)
      call MAPL_GridCompGetResource(GC, 'WSUB_OPTION', WSUB_OPTION, default=3, _RC)!0- param 1- Use Wsub climatology 2-USE WNET` 3-Symbolic Wnet

      ! MAT These have to be defined as they are passed into Aer_Activate below and are intent(in)
      !     Note: It's possible these aren't *used* if USE_AEROSOL_NN=.TRUE. but they are still passed
      !           in so they have to be defined
      call MAPL_GridCompGetResource(GC, 'NCCN_OCN', CCN_OCN, default=50., _RC)
      call MAPL_GridCompGetResource(GC, 'NCCN_LND', CCN_LND, default=300., _RC)

      call MAPL_GridCompGetResource(GC, 'NCIN_OCN', CIN_OCN, default=0.005, _RC)
      call MAPL_GridCompGetResource(GC, 'NCIN_LND', CIN_LND, default=0.05, _RC)

      ! NOTE: Binary restarts expect Q to be the first field in the moist_internal_rst. Thus,
      !       the first MAPL_AddInternalSpec call must be from the microphysics

      ! MAPL3: only GFDL_1M, GF and UW are ported so far. The other schemes
      !        fail here, before any of their code runs.

      if (adjustl(CLDMICR_OPTION) == "BACM_1M") then
         _FAIL('CLDMICR_OPTION BACM_1M is not yet ported to MAPL3')
         call BACM_1M_Setup(GC, CF, _RC)
      end if
      if (adjustl(CLDMICR_OPTION) == "MGB2_2M") then
         _FAIL('CLDMICR_OPTION MGB2_2M is not yet ported to MAPL3')
         call MGB2_2M_Setup(GC, CF, _RC)
      end if
      if (adjustl(CLDMICR_OPTION) == "GFDL_1M") call GFDL_1M_Setup(GC, ADVECTED, _RC)
      if (adjustl(CLDMICR_OPTION) == "THOM_1M") then
         _FAIL('CLDMICR_OPTION THOM_1M is not yet ported to MAPL3')
         call THOM_1M_Setup(GC, CF, _RC)
      end if
      if (adjustl(CONVPAR_OPTION) == "RAS") then
         _FAIL('CONVPAR_OPTION RAS is not yet ported to MAPL3')
         call RAS_Setup(GC, CF, _RC)
      end if
      if (adjustl(CONVPAR_OPTION) == "GF") call GF_Setup(GC, ADVECTED, _RC)
      if (adjustl(SHALLOW_OPTION) == "UW") call UW_Setup(GC, _RC)

      ! Advection (replaces MAPL2's FRIENDLYTO=DYNAMICS): one subscription to the TRADV
      ! service for the internals the schemes listed in ADVECTED (items must be internals)
      if (ADVECTED%size() > 0) then
         call MAPL_GridCompAddSpec(GC, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name='TRADV', &
              standard_name='advected_quantities', &
              itemtype=MAPL_STATEITEM_SERVICE, &
              service_items=ADVECTED, _RC)
      end if

      ! !IMPORT STATE:

#include "Moist_Import___.h"

      call MAPL_GridCompAddSpec(GC, &
           state_intent=ESMF_STATEINTENT_IMPORT, &
           short_name='AERO', &
           standard_name='aerosols', &
           itemtype=MAPL_STATEITEM_STATE, _RC)

      ! !EXPORT STATE:

#include "Moist_Export___.h"

      ! Convective transport service (replaces MAPL2's MTR import, which Physics filled
      ! with the fields FRIENDLYTO=MOIST). Subscribers list their internals as service
      ! items; MAPL adds the bundle to both the import and the export state.
      call MAPL_GridCompAddSpec(GC, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name='MTR', &
           standard_name='tracers_for_moist', &
           itemtype=MAPL_STATEITEM_SERVICE, _RC)

      !EOS

      _RETURN(_SUCCESS)

   end subroutine SetServices

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!-!-!-!!!!!!!!!!!!!!!

   !EOP

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !BOP

   ! !IROUTINE: Initialize -- Initialize method for the composite Moist Gridded Component

   ! !INTERFACE:

   subroutine Initialize(GC, IMPORT, EXPORT, CLOCK, RC)

      ! !ARGUMENTS:

      type(ESMF_GridComp), intent(inout) :: GC ! Gridded component
      type(ESMF_State), intent(inout) :: IMPORT ! Import state
      type(ESMF_State), intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock), intent(inout) :: CLOCK ! The clock
      integer, optional, intent(out) :: RC ! Error code

      ! !DESCRIPTION: The Initialize method of the Moist Physics Gridded Component first
      !   calls the Initialize method of the child Dynamics.  The Dynamics Initialize method will
      !   create the ESMF GRID, which will then be used to set the GRID associated with the
      !   SuperDyn Composite Component itself.  It should be noted that the
      !   SuperDyn Initialize method also invokes the GEOS Topo Utility which creates all
      !   topography related quantities.

      !EOP

      ! ErrLog Variables

      integer :: STATUS

      ! Local derived type aliases

      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_Grid) :: GRID
      type(ESMF_State), pointer :: GIM(:)
      type(ESMF_State), pointer :: GEX(:)
      type(ESMF_State) :: INTERNAL

      type(ESMF_Config) :: CF

      type(ESMF_Alarm) :: ALARM
      type(ESMF_TimeInterval) :: TINT
      real(kind=ESMF_KIND_R8) :: DT_R8
      real :: DT_MOIST
      real :: DBZ_DT
      type(ESMF_Calendar) :: calendar
      type(ESMF_Time) :: currTime
      type(ESMF_Alarm) :: DBZ_RunAlarm
      type(ESMF_TimeInterval) :: ringInterval

      !=============================================================================

      ! Begin...

      ! Get the configuration
      ! ---------------------
      call ESMF_GridCompGet(GC, config=CF, _RC)

      ! Get my internal MAPL_Generic state
      !-----------------------------------

      call MAPL_GetObjectFromGC(GC, MAPL, _RC)

      ! Call Generic Initialize for MOIST GC
      !----------------------------------------

      call MAPL_GenericInitialize(GC, IMPORT, EXPORT, CLOCK, _RC)

      ! Get parameters from generic state.
      !-----------------------------------
      call MAPL_GetResource(MAPL, LDIAGNOSE_PRECIP_TYPE, Label="DIAGNOSE_PRECIP_TYPE:", default=.false., _RC)
      call MAPL_GetResource(MAPL, LUPDATE_PRECIP_TYPE, Label="UPDATE_PRECIP_TYPE:", default=.false., _RC)

      DETRAIN_INACTIVE_CNV = 0.0
      if (adjustl(CLDMICR_OPTION) == "GFDL_1M") DETRAIN_INACTIVE_CNV = 5.0e-4
      call MAPL_GetResource(MAPL, DETRAIN_INACTIVE_CNV, Label="DETRAIN_INACTIVE_CNV:", default=DETRAIN_INACTIVE_CNV, &
           _RC)
      call MAPL_GetResource(MAPL, TAU_DETRAIN_CNV, Label="TAU_DETRAIN_CNV:", default=3600.0, _RC)

      if (adjustl(CONVPAR_OPTION) == "RAS") call RAS_Initialize(MAPL, _RC)
      if (adjustl(CONVPAR_OPTION) == "GF") call GF_Initialize(MAPL, CF, CLOCK, IMPORT, EXPORT, _RC)
      if (adjustl(SHALLOW_OPTION) == "UW") call UW_Initialize(MAPL, CF, CLOCK, IMPORT, EXPORT, _RC)
      if (adjustl(CLDMICR_OPTION) == "BACM_1M") call BACM_1M_Initialize(MAPL, _RC)
      if (adjustl(CLDMICR_OPTION) == "GFDL_1M") call GFDL_1M_Initialize(MAPL, CF, CLOCK, IMPORT, EXPORT, _RC)
      if (adjustl(CLDMICR_OPTION) == "THOM_1M") call THOM_1M_Initialize(MAPL, _RC)
      if (adjustl(CLDMICR_OPTION) == "MGB2_2M") call MGB2_2M_Initialize(MAPL, _RC)

      call MAPL_Get(MAPL, RUNALARM=ALARM, _RC)
      call ESMF_AlarmGet(ALARM, ringInterval=TINT, _RC)
      call ESMF_TimeIntervalGet(TINT, S_R8=DT_R8, _RC)
      DT_MOIST = DT_R8

      DBZ_DT = max(DT_MOIST, 900.0)
      call MAPL_GetResource(MAPL, DBZ_DT, 'DBZ_DT:', default=DBZ_DT, _RC)
      ! Get the current time in addition to the calendar
      call ESMF_ClockGet(CLOCK, currTime=currTime, calendar=calendar, _RC)
      call ESMF_TimeIntervalSet(ringInterval, S=nint(DBZ_DT), calendar=calendar, _RC)
      ! Add RingTime = currTime to anchor the alarm
      DBZ_RunAlarm = ESMF_AlarmCreate(CLOCK=CLOCK, &
           Name='DBZ_RunAlarm', &
           RingTime=currTime - TINT, &
           ringInterval=ringInterval, &
           Sticky=.false., _RC)
      call init_refl10cm()
      call MAPL_GetResource(MAPL, refl10cm_allow_wet_graupel, 'refl10cm_allow_wet_graupel:', default=.false., _RC)
      call MAPL_GetResource(MAPL, refl10cm_allow_wet_snow, 'refl10cm_allow_wet_snow:', default=.false., _RC)
      call MAPL_GetResource(MAPL, DBZ_VAR_INTERCP, 'DBZ_VAR_INTERCP:', default=DBZ_VAR_INTERCP, _RC)
      call MAPL_GetResource(MAPL, LIQUID_SKIN_SNOW, 'LIQUID_SKIN_SNOW:', default=.false., _RC)
      call MAPL_GetResource(MAPL, LIQUID_SKIN_GRAUPEL, 'LIQUID_SKIN_GRAUPEL:', default=.false., _RC)
      call MAPL_GetResource(MAPL, LIQUID_SKIN_HAIL, 'LIQUID_SKIN_HAIL:', default=.true., _RC)

      ! All done
      !---------

      _RETURN(_SUCCESS)
   end subroutine Initialize

   !===================================================================================

   !BOP

   ! !IROUTINE: RUN -- Run method for the CONVECT component

   ! !INTERFACE:

   subroutine Run(GC, IMPORT, EXPORT, CLOCK, RC)

      ! !ARGUMENTS:

      type(ESMF_GridComp), intent(inout) :: GC ! Gridded component
      type(ESMF_State), intent(inout) :: IMPORT ! Import state
      type(ESMF_State), intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock), intent(inout) :: CLOCK ! The clock
      integer, optional, intent(out) :: RC ! Error code:

      ! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
      !                the Initialize and Finalize services, as well as allocating

      !EOP

      ! ErrLog Variables

      integer :: STATUS

      type(ESMF_VM) :: VMG

      ! Local derived type aliases

      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_Config) :: CF
      type(ESMF_State) :: INTERNAL
      type(ESMF_Alarm) :: ALARM
      type(ESMF_TimeInterval) :: TINT
      real(kind=ESMF_KIND_R8) :: DT_R8
      real :: DT_MOIST

      ! Local variables
      real :: MFC ! Layer-centered mass flux [kg m-2 s-1]
      real :: inactivity_weight ! Scale from 0.0 to 1.0 based on MFC
      real :: transfer_rate ! Fraction of mass to transfer this timestep
      real :: dq_l ! Liquid mass being transferred
      real :: dq_i ! Ice mass being transferred
      real :: d_cf ! Cloud fraction being transferred
      real :: Tmax, KCBLMIN, PMIN_CBL
      real :: CNV_CAPE_NORM, CNV_CAPE_SCALE
      real, allocatable, dimension(:, :, :) :: PLEmb, PKE, ZLE0, PK, MASS
      real, allocatable, dimension(:, :, :) :: PLmb, ZL0, DZET
      real, allocatable, dimension(:, :, :) :: QST3, DQST3, MWFA
      real, allocatable, dimension(:, :, :) :: TMP3D, TMP3Dp1
      real, allocatable, dimension(:, :) :: TMP2D
      integer, allocatable, dimension(:, :) :: KLCL
      ! Internals
      real, pointer, dimension(:, :, :) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN, QRTOT, QSTOT, QGTOT
      real, pointer, dimension(:, :, :) :: NACTL, NACTI
      ! Imports
      real, pointer, dimension(:, :, :) :: ZLE, PLE, T, U, V, W
      real, pointer, dimension(:, :) :: FRLAND, FRLANDICE, FRACI, SNOMAS
      real, pointer, dimension(:, :) :: SH, TS, EVAP, KPBL
      real, pointer, dimension(:, :, :) :: KH, TKE, OMEGA
      integer :: n_modes
      type(ESMF_State) :: AERO
      type(ESMF_FieldBundle) :: TR
      ! Exports
      real, pointer, dimension(:, :, :) :: DQDT, DQADT, DQIDT, DQLDT, DQRDT, DQSDT, DQGDT
      real, pointer, dimension(:, :, :) :: DTDT, DUDT, DVDT, DWDT
      real, pointer, dimension(:, :, :) :: DPDTMST, PFL_LSAN, PFI_LSAN
      real, pointer, dimension(:, :, :) :: DTDT_ER, DQVDT_ER
      real, pointer, dimension(:, :) :: PTYPE, TPREC, CN_PRCP, LS_PRCP, AN_PRCP, SC_PRCP, PLS, PCU
      real, pointer, dimension(:, :) :: RAIN, SNOW, ICE, FRZR, PREC_STRAT, PREC_CONV
      real, pointer, dimension(:, :, :) :: BYNCY
      real, pointer, dimension(:, :) :: CAPE, INHB, MLCAPE, SBCAPE, MLCIN, MUCAPE, MUCIN, SBCIN, LFC, LNB, LCL_AGL
      real, pointer, dimension(:, :) :: CNV_FRC, SRF_TYPE
      real, pointer, dimension(:, :, :) :: CFICE, CFLIQ
      real, pointer, dimension(:, :, :) :: NWFA, NIFA
      real, pointer, dimension(:, :) :: EIS, LTS
      real, pointer, dimension(:, :, :) :: PTRDC, PTRSC
      real, pointer, dimension(:, :, :) :: PTR3D
      real, pointer, dimension(:, :) :: PTR2D
      real, pointer, dimension(:) :: PTR1D

      integer :: IM, JM, LM
      integer :: I, J, L, n

      !=============================================================================

      ! Begin...

      ! Get my configuration and VM
      ! ---------------------------

      call ESMF_GridCompGet(GC, config=CF, VM=VMG, _RC)

      ! Get my internal MAPL_Generic state
      !-----------------------------------

      call MAPL_GetObjectFromGC(GC, MAPL, _RC)

      call MAPL_TimerOn(MAPL, "TOTAL")

      ! If its time, call run methods
      ! --------------------------------------------

      call MAPL_Get(MAPL, IM=IM, JM=JM, LM=LM, &
           RUNALARM=ALARM, &
           INTERNAL_ESMF_STATE=INTERNAL, &
           _RC)

      call ESMF_AlarmGet(ALARM, ringInterval=TINT, _RC)
      call ESMF_TimeIntervalGet(TINT, S_R8=DT_R8, _RC)
      DT_MOIST = DT_R8

      if (ESMF_AlarmIsRinging(ALARM, RC=STATUS)) then

         call MAPL_TimerOn(MAPL, "---MOIST_PROLOGUE")

         call ESMF_AlarmRingerOff(ALARM, _RC)

         ! Internal State
         call MAPL_GetPointer(INTERNAL, Q, 'Q', _RC)
         call MAPL_GetPointer(INTERNAL, QLLS, 'QLLS', _RC)
         call MAPL_GetPointer(INTERNAL, QLCN, 'QLCN', _RC)
         call MAPL_GetPointer(INTERNAL, CLCN, 'CLCN', _RC)
         call MAPL_GetPointer(INTERNAL, CLLS, 'CLLS', _RC)
         call MAPL_GetPointer(INTERNAL, QILS, 'QILS', _RC)
         call MAPL_GetPointer(INTERNAL, QICN, 'QICN', _RC)
         call MAPL_GetPointer(INTERNAL, NACTL, 'NACTL', _RC)
         call MAPL_GetPointer(INTERNAL, NACTI, 'NACTI', _RC)

         ! Import State
         call MAPL_GetPointer(IMPORT, PLE, 'PLE', _RC)
         call MAPL_GetPointer(IMPORT, ZLE, 'ZLE', _RC)
         call MAPL_GetPointer(IMPORT, TS, 'TS', _RC)
         call MAPL_GetPointer(IMPORT, T, 'T', _RC)
         call MAPL_GetPointer(IMPORT, U, 'U', _RC)
         call MAPL_GetPointer(IMPORT, V, 'V', _RC)
         call MAPL_GetPointer(IMPORT, W, 'W', _RC)
         call MAPL_GetPointer(IMPORT, KPBL, 'KPBL', _RC)
         call MAPL_GetPointer(IMPORT, KH, 'KH', _RC)
         call MAPL_GetPointer(IMPORT, SH, 'SH', _RC)
         call MAPL_GetPointer(IMPORT, EVAP, 'EVAP', _RC)
         call MAPL_GetPointer(IMPORT, OMEGA, 'OMEGA', _RC)
         call MAPL_GetPointer(IMPORT, TKE, 'TKE', _RC)
         call ESMF_StateGet(IMPORT, 'AERO', AERO, _RC)
         call ESMF_StateGet(IMPORT, 'MTR', TR, _RC)

         ! Update SRF_TYPE for ice_fraction
         call MAPL_GetPointer(IMPORT, FRLAND, 'FRLAND', _RC)
         call MAPL_GetPointer(IMPORT, FRLANDICE, 'FRLANDICE', _RC)
         call MAPL_GetPointer(IMPORT, FRACI, 'FRACI', _RC)
         call MAPL_GetPointer(IMPORT, SNOMAS, 'SNOMAS', _RC)
         call MAPL_GetPointer(EXPORT, SRF_TYPE, 'SRF_TYPE', ALLOC=.true., _RC)

         where (FRLANDICE > 0.5)
            SRF_TYPE = SRF_TYPE_LANDICE
            elsewhere (FRACI > 0.5)
            SRF_TYPE = SRF_TYPE_ICE
            elsewhere (SNOMAS > 0.1 .and. SNOMAS /= MAPL_UNDEF)
            ! NOTE: SNOMAS has UNDEFs so we need to make sure we don't
            !       allow that to infect this comparison
            SRF_TYPE = SRF_TYPE_SNOW
            elsewhere (FRLAND > 0.1)
            SRF_TYPE = SRF_TYPE_LAND
            elsewhere
            SRF_TYPE = SRF_TYPE_OCEAN
         end where

         ! Allocatables
         ! Edge variables
         allocate(ZLE0(IM, JM, 0:LM))
         allocate(PLEmb(IM, JM, 0:LM))
         allocate(PKE(IM, JM, 0:LM))
         ! Layer variables
         allocate(ZL0(IM, JM, LM))
         allocate(DZET(IM, JM, LM))
         allocate(PLmb(IM, JM, LM))
         allocate(PK(IM, JM, LM))
         allocate(DQST3(IM, JM, LM))
         allocate(QST3(IM, JM, LM))
         allocate(MASS(IM, JM, LM))
         allocate(TMP3D(IM, JM, LM))
         allocate(TMP2D(IM, JM))
         allocate(KLCL(IM, JM))

         ! Save input winds
         call MAPL_GetPointer(EXPORT, PTR3D, 'UMST0', ALLOC=.true., _RC)
         PTR3D = U
         call MAPL_GetPointer(EXPORT, PTR3D, 'VMST0', ALLOC=.true., _RC)
         PTR3D = V

         ! pre-fill default condensate radii
         call MAPL_GetPointer(EXPORT, PTR3D, 'RL', _RC)
         if (associated(PTR3D)) PTR3D = MAPL_UNDEF
         call MAPL_GetPointer(EXPORT, PTR3D, 'RI', _RC)
         if (associated(PTR3D)) PTR3D = MAPL_UNDEF
         call MAPL_GetPointer(EXPORT, PTR3D, 'RR', _RC)
         if (associated(PTR3D)) PTR3D = 50.e-6
         call MAPL_GetPointer(EXPORT, PTR3D, 'RS', _RC)
         if (associated(PTR3D)) PTR3D = 50.e-6
         call MAPL_GetPointer(EXPORT, PTR3D, 'RG', _RC)
         if (associated(PTR3D)) PTR3D = 50.e-6

         ! Derived States
         MASS = (PLE(:, :, 1:LM) - PLE(:, :, 0:LM - 1)) / MAPL_GRAV
         PLEmb = PLE * .01
         PKE = (PLE / MAPL_P00)**(MAPL_KAPPA)
         PLmb = 0.5 * (PLEmb(:, :, 0:LM - 1) + PLEmb(:, :, 1:LM))
         PK = (100.0 * PLmb / MAPL_P00)**(MAPL_KAPPA)
         do L = 0, LM
            ZLE0(:, :, L) = ZLE(:, :, L) - ZLE(:, :, LM) ! Edge Height (m) above the surface
         end do
         ZL0 = 0.5 * (ZLE0(:, :, 0:LM - 1) + ZLE0(:, :, 1:LM)) ! Layer Height (m) above the surface
         DZET = (ZLE0(:, :, 0:LM - 1) - ZLE0(:, :, 1:LM)) ! Layer thickness (m)
         DQST3 = GEOS_DQSAT(T, PLmb, QSAT=QST3)

         ! Lower tropospheric stability and estimated inversion strength
         call MAPL_GetPointer(EXPORT, LTS, 'LTS', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, EIS, 'EIS', ALLOC=.true., _RC)
         KLCL = FIND_KLCL(T, Q, PLmb, IM, JM, LM)
         call MAPL_GetPointer(EXPORT, PTR2D, 'ZLCL', _RC)
         if (associated(PTR2D)) then
            do J = 1, JM
               do I = 1, IM
                  PTR2D(I, J) = ZL0(I, J, KLCL(I, J))
               end do
            end do
         end if
         TMP3D = (100.0 * PLmb / MAPL_P00)**(MAPL_KAPPA)
         call FIND_EIS(T / TMP3D, QST3, T, ZL0, PLEmb, KLCL, IM, JM, LM, LTS, EIS)

         ! Recording of import/internal vars into export if desired
         !---------------------------------------------------------
         call MAPL_GetPointer(EXPORT, PTR3D, 'QX0', _RC)
         if (associated(PTR3D)) PTR3D = Q
         call MAPL_GetPointer(EXPORT, PTR3D, 'QLLSX0', _RC)
         if (associated(PTR3D)) PTR3D = QLLS
         call MAPL_GetPointer(EXPORT, PTR3D, 'QLCNX0', _RC)
         if (associated(PTR3D)) PTR3D = QLCN
         call MAPL_GetPointer(EXPORT, PTR3D, 'CLCNX0', _RC)
         if (associated(PTR3D)) PTR3D = CLCN
         call MAPL_GetPointer(EXPORT, PTR3D, 'CLLSX0', _RC)
         if (associated(PTR3D)) PTR3D = CLLS
         call MAPL_GetPointer(EXPORT, PTR3D, 'QILSX0', _RC)
         if (associated(PTR3D)) PTR3D = QILS
         call MAPL_GetPointer(EXPORT, PTR3D, 'QICNX0', _RC)
         if (associated(PTR3D)) PTR3D = QICN
         call MAPL_GetPointer(EXPORT, PTR3D, 'QCCNX0', _RC)
         if (associated(PTR3D)) PTR3D = QICN + QLCN
         call MAPL_GetPointer(EXPORT, PTR3D, 'QCLSX0', _RC)
         if (associated(PTR3D)) PTR3D = QILS + QLLS

         ! Trajectory for the moist TLM/ADJ
         !---------------------------------------------------------
         call MAPL_GetPointer(EXPORT, PTR3D, 'TH_moist', _RC)
         if (associated(PTR3D)) PTR3D = T / PK
         call MAPL_GetPointer(EXPORT, PTR3D, 'Q_moist', _RC)
         if (associated(PTR3D)) PTR3D = Q
         call MAPL_GetPointer(EXPORT, PTR2D, 'TS_moist', _RC)
         if (associated(PTR2D)) PTR2D = TS
         call MAPL_GetPointer(EXPORT, PTR2D, 'KHu_moist', _RC)
         if (associated(PTR2D)) then
            PTR2D = -1
            do L = 0, LM
               where ((PTR2D == -1) .and. (KH(:, :, L) > 2.0))
                  PTR2D = L * 1.0
               end where
            end do
         end if
         call MAPL_GetPointer(EXPORT, PTR2D, 'KHl_moist', _RC)
         if (associated(PTR2D)) then
            PTR2D = -1
            do L = LM, 0, -1
               where ((PTR2D == -1) .and. (KH(:, :, L) > 2.0))
                  PTR2D = L * 1.0
               end where
            end do
         end if

         call MAPL_GetPointer(EXPORT, PTR2D, 'KCBL_moist', _RC)
         if (associated(PTR2D)) then
            call MAPL_GetResource(MAPL, PMIN_CBL, 'PMIN_CBL', default=50000.0, _RC)
            call MAPL_GetPointer(IMPORT, PTR1D, 'PREF', _RC)
            KCBLMIN = count(PTR1D < PMIN_CBL)
            do J = 1, JM
               do I = 1, IM
                  if (nint(KPBL(I, J)) /= 0) then
                     PTR2D(I, J) = max(min(nint(KPBL(I, J)) + 1, LM - 1), 1)
                  else
                     PTR2D(I, J) = LM - 1
                  end if
                  PTR2D(I, J) = max(PTR2D(I, J), KCBLMIN)
               end do
            end do
         end if

         ! These may be used by children
         call MAPL_GetPointer(EXPORT, NWFA, 'NWFA', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, NIFA, 'NIFA', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, CNV_FRC, 'CNV_FRC', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, BYNCY, 'BYNCY', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, CAPE, 'CAPE', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, INHB, 'INHB', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, SBCAPE, 'SBCAPE', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, SBCIN, 'SBCIN', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, MLCAPE, 'MLCAPE', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, MLCIN, 'MLCIN', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, MUCAPE, 'MUCAPE', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, MUCIN, 'MUCIN', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, LFC, 'ZLFC', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, LNB, 'ZLNB', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, LCL_AGL, 'LCL_AGL', ALLOC=.true., _RC)

         call MAPL_TimerOn(MAPL, "-----BUOYANCY2")
         call BUOYANCY2(IM, JM, LM, T, Q, QST3, DQST3, DZET, ZL0, PLmb, PLEmb(:, :, LM), &
              SBCAPE, MLCAPE, MUCAPE, SBCIN, MLCIN, MUCIN, BYNCY, LFC, LNB, LCL_AGL)
         call MAPL_TimerOff(MAPL, "-----BUOYANCY2")

         call BUOYANCY(T, Q, QST3, DQST3, DZET, ZL0, BYNCY, CAPE, INHB)

         ! initialize diagnosed convective fraction
         CNV_FRC = 0.0
         if (CNV_FRACTION_MAX > CNV_FRACTION_MIN) then
            if (CNV_FRACTION_EXP /= 1.0) then
               where (CAPE /= MAPL_UNDEF)
                  CNV_FRC = (max(1.e-6, min(1.0, (CAPE - CNV_FRACTION_MIN) / (CNV_FRACTION_MAX - CNV_FRACTION_MIN))))
               end where
               CNV_FRC = CNV_FRC**CNV_FRACTION_EXP
            else
               where (CAPE /= MAPL_UNDEF)
                  CNV_FRC = (1.0 - COS(MAPL_PI * (CAPE - CNV_FRACTION_MIN) / (CNV_FRACTION_MAX - CNV_FRACTION_MIN))) &
                       / 2.0
               end where
               where ((CAPE <= CNV_FRACTION_MIN) .and. (CAPE /= MAPL_UNDEF))
                  CNV_FRC = 0.0
               end where
               where ((CAPE >= CNV_FRACTION_MAX) .and. (CAPE /= MAPL_UNDEF))
                  CNV_FRC = 1.0
               end where
            end if
         end if

         call MAPL_TimerOff(MAPL, "---MOIST_PROLOGUE")

         ! Extract convective tracers from the TR bundle
         call MAPL_TimerOn(MAPL, "---CONV_TRACERS")
         call CNV_Tracers_Init(TR, RC)
         call MAPL_TimerOff(MAPL, "---CONV_TRACERS")

         ! Get aerosol activation properties
         call MAPL_TimerOn(MAPL, "---AERO_ACTIVATE")

         if (MOIST_USE_NCLOUD_CLIM) then !Setup ND/NI climatology from GiOcean
            call MAPL_GetPointer(IMPORT, PTR3D, 'NCPL_CLIM', _RC)
            NACTL = PTR3D
            call MAPL_GetPointer(IMPORT, PTR3D, 'NCPI_CLIM', _RC)
            NACTI = PTR3D
         else
            if (USE_AEROSOL_NN) then
               ! get veritical velocity
               if (all(W == 0.0)) then
                  TMP3D = -OMEGA / (MAPL_GRAV * PLmb * 100.0 / (MAPL_RGAS * T))
               else
                  TMP3D = W
               end if
               ! Call aerosol activation (pressures in Pa)
               call Aer_Activation(MAPL, IM, JM, LM, Q, T, PLmb * 100.0, PLE, TKE, TMP3D, FRLAND, &
                    AERO, NACTL, NACTI, NWFA, NIFA, CCN_LND * 1.e6, CCN_OCN * 1.e6, &
                    .true., _RC)
               ! Apply scaling factors
               NACTL = NACTL * NN_FAC_LIQ
               NACTI = NACTI * NN_FAC_ICE
               ! Apply min/max limits
               NACTL = max(NN_MIN_LIQ, min(NACTL, NN_MAX_LIQ))
               NACTI = max(NN_MIN_ICE, min(NACTI, NN_MAX_ICE))
            else
               do L = 1, LM
                  NACTL(:, :, L) = (CCN_LND * FRLAND + CCN_OCN * (1.0 - FRLAND)) * 1.e6 ! #/m^3
                  NACTI(:, :, L) = (CIN_LND * FRLAND + CIN_OCN * (1.0 - FRLAND)) * 1.e6 ! #/m^3
               end do
            end if
         end if

         call MAPL_GetPointer(EXPORT, PTR3D, 'NCCN_LIQ', _RC)
         if (associated(PTR3D)) PTR3D = NACTL
         call MAPL_GetPointer(EXPORT, PTR3D, 'NCCN_ICE', _RC)
         if (associated(PTR3D)) PTR3D = NACTI

         call MAPL_TimerOff(MAPL, "---AERO_ACTIVATE")

         if (SH_MD_DP) then
            if (adjustl(SHALLOW_OPTION) == "UW") call UW_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
            if (adjustl(SHALLOW_OPTION) == "UW") then
               if (DEBUG_MST) then
                  call MAPL_MaxMin('MST: Q_AF_UW  ', Q)
                  call MAPL_MaxMin('MST: T_AF_UW  ', T)
                  call MAPL_MaxMin('MST: U_AF_UW  ', U)
                  call MAPL_MaxMin('MST: V_AF_UW  ', V)
               end if
            end if
            if (adjustl(CONVPAR_OPTION) == "RAS") call RAS_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
            if (adjustl(CONVPAR_OPTION) == "GF") call GF_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
            if (adjustl(CONVPAR_OPTION) == "GF") then
               if (DEBUG_MST) then
                  call MAPL_MaxMin('MST: Q_AF_GF  ', Q)
                  call MAPL_MaxMin('MST: T_AF_GF  ', T)
                  call MAPL_MaxMin('MST: U_AF_GF  ', U)
                  call MAPL_MaxMin('MST: V_AF_GF  ', V)
               end if
            end if
         else
            if (adjustl(CONVPAR_OPTION) == "RAS") call RAS_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
            if (adjustl(CONVPAR_OPTION) == "GF") call GF_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
            if (adjustl(CONVPAR_OPTION) == "GF") then
               if (DEBUG_MST) then
                  call MAPL_MaxMin('MST: Q_AF_GF  ', Q)
                  call MAPL_MaxMin('MST: T_AF_GF  ', T)
                  call MAPL_MaxMin('MST: U_AF_GF  ', U)
                  call MAPL_MaxMin('MST: V_AF_GF  ', V)
               end if
            end if
            if (adjustl(SHALLOW_OPTION) == "UW") call UW_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
            if (adjustl(SHALLOW_OPTION) == "UW") then
               if (DEBUG_MST) then
                  call MAPL_MaxMin('MST: Q_AF_UW  ', Q)
                  call MAPL_MaxMin('MST: T_AF_UW  ', T)
                  call MAPL_MaxMin('MST: U_AF_UW  ', U)
                  call MAPL_MaxMin('MST: V_AF_UW  ', V)
               end if
            end if
         end if

         call MAPL_TimerOn(MAPL, "---MOIST_EPILOGUE")

         ! Mass fluxes
         allocate(TMP3Dp1(IM, JM, LM + 1))
         ! Edge updraft mass flux over deep and shalow convection
         call MAPL_GetPointer(EXPORT, PTR3D, 'CNV_MFC', _RC)
         call MAPL_GetPointer(EXPORT, PTRDC, 'UMF_DC', _RC)
         call MAPL_GetPointer(EXPORT, PTRSC, 'UMF_SC', _RC)
         TMP3Dp1 = 0.0
         if (associated(PTRDC)) TMP3Dp1 = TMP3Dp1 + PTRDC
         if (associated(PTRSC)) TMP3Dp1 = TMP3Dp1 + PTRSC
         if (associated(PTR3D)) PTR3D = TMP3Dp1
         if (DETRAIN_INACTIVE_CNV > 0.0) then
            do L = 1, LM
               do J = 1, JM
                  do I = 1, IM
                     ! Calculate local mass flux
                     MFC = 0.5 * (TMP3Dp1(I, J, L) + TMP3Dp1(I, J, L + 1))
                     if (MFC < DETRAIN_INACTIVE_CNV) then
                        ! 1. Calculate a smooth inactivity factor (0.0 at threshold, 1.0 when MFC is 0)
                        ! 2. Scale it by the timestep vs relaxation time (DT_MOIST / TAU)
                        inactivity_weight = 1.0 - (MFC / DETRAIN_INACTIVE_CNV)
                        transfer_rate = inactivity_weight * (DT_MOIST / TAU_DETRAIN_CNV)
                        ! Bound the rate safely between 0 and 1
                        transfer_rate = min(1.0, max(0.0, transfer_rate))
                        ! Calculate the exact amounts to transfer this timestep
                        dq_l = QLCN(I, J, L) * transfer_rate
                        dq_i = QICN(I, J, L) * transfer_rate
                        d_cf = CLCN(I, J, L) * transfer_rate
                        ! Move the Liquid
                        QLLS(I, J, L) = QLLS(I, J, L) + dq_l
                        QLCN(I, J, L) = QLCN(I, J, L) - dq_l
                        ! Move the Ice
                        QILS(I, J, L) = QILS(I, J, L) + dq_i
                        QICN(I, J, L) = QICN(I, J, L) - dq_i
                        ! Move the Cloud Fraction using Random Overlap for the transferred piece
                        CLLS(I, J, L) = CLLS(I, J, L) + d_cf - (CLLS(I, J, L) * d_cf)
                        CLCN(I, J, L) = CLCN(I, J, L) - d_cf
                     end if
                  end do
               end do
            end do
         end if
         deallocate(TMP3Dp1)

         call MAPL_GetPointer(EXPORT, PTR3D, 'CNV_MFD', _RC)
         call MAPL_GetPointer(EXPORT, PTRDC, 'MFD_DC', _RC)
         call MAPL_GetPointer(EXPORT, PTRSC, 'MFD_SC', _RC)
         TMP3D = 0.0
         if (associated(PTRDC)) TMP3D = TMP3D + PTRDC
         if (associated(PTRSC)) TMP3D = TMP3D + PTRSC
         if (associated(PTR3D)) PTR3D = TMP3D

         call MAPL_TimerOff(MAPL, "---MOIST_EPILOGUE")

         if (adjustl(CLDMICR_OPTION) == "BACM_1M") call BACM_1M_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
         if (adjustl(CLDMICR_OPTION) == "GFDL_1M") call GFDL_1M_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
         if (adjustl(CLDMICR_OPTION) == "THOM_1M") call THOM_1M_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
         if (adjustl(CLDMICR_OPTION) == "MGB2_2M") call MGB2_2M_Run(GC, IMPORT, EXPORT, CLOCK, _RC)

         call MAPL_TimerOn(MAPL, "---MOIST_EPILOGUE")

         if (DEBUG_MST) then
            call MAPL_MaxMin('MST: Q_AF_MP  ', Q)
            call MAPL_MaxMin('MST: T_AF_MP  ', T)
            call MAPL_MaxMin('MST: U_AF_MP  ', U)
            call MAPL_MaxMin('MST: V_AF_MP  ', V)
         end if

         ! Exports
         ! Cloud fraction exports
         call MAPL_GetPointer(EXPORT, CFICE, 'CFICE', ALLOC=.true., _RC)
         if (associated(CFICE)) then
            CFICE = 0.0
            where (QILS + QICN > 1.0e-12)
               CFICE = (CLLS + CLCN) * (QILS + QICN) / (QLLS + QLCN + QILS + QICN)
            end where
            CFICE = max(min(CFICE, 1.0), 0.0)
         end if
         call MAPL_GetPointer(EXPORT, CFLIQ, 'CFLIQ', _RC)
         if (associated(CFLIQ)) then
            CFLIQ = 0.0
            where (QLLS + QLCN > 1.0e-12)
               CFLIQ = (CLLS + CLCN) * (QLLS + QLCN) / (QLLS + QLCN + QILS + QICN)
            end where
            CFLIQ = max(min(CFLIQ, 1.0), 0.0)
         end if

         ! some diagnostics to export
         QST3 = GEOS_QsatICE(T, PLmb * 100.0, DQ=DQST3)
         call MAPL_GetPointer(EXPORT, PTR3D, 'RHICE', _RC)
         if (associated(PTR3D)) then
            PTR3D = Q / QST3
            where (T > MAPL_TICE)
               PTR3D = 0.0
            end where
         end if
         call MAPL_GetPointer(EXPORT, PTR3D, 'SAT_RAT', _RC)
         if (associated(PTR3D)) then
            where (CFICE < 0.99 .and. QST3 > 1.0e-20)
               TMP3D = max((Q - QST3 * CFICE), 0.0) / (1.0 - CFICE)
               PTR3D = min(TMP3D / QST3, 2.0)
               elsewhere
               PTR3D = 1.0
            end where
         end if

         QST3 = GEOS_QsatLQU(T, PLmb * 100.0, DQ=DQST3)
         call MAPL_GetPointer(EXPORT, PTR3D, 'RHLIQ', _RC)
         if (associated(PTR3D)) PTR3D = Q / QST3

         ! Clean up supersaturation only with respect to liquid water
         call MAPL_GetPointer(EXPORT, DTDT_ER, 'DTDT_ER', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, DQVDT_ER, 'DQVDT_ER', ALLOC=.true., _RC)
         DTDT_ER = T
         DQVDT_ER = Q

         ! Rain-out of Relative Humidity where RH > 110%
         call MAPL_GetPointer(EXPORT, LS_PRCP, 'LS_PRCP', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, PTR2D, 'ER_PRCP', ALLOC=.true., _RC)
         where (Q > 1.1 * QST3)
            TMP3D = (Q - 1.1 * QST3) / (1.0 + 1.1 * DQST3 * MAPL_ALHL / MAPL_CP)
            elsewhere
            TMP3D = 0.0
         end where
         ! protect against producing negatice QV
         where (TMP3D > Q)
            TMP3D = Q
         end where
         PTR2D = SUM(TMP3D * MASS, 3) / DT_MOIST
         LS_PRCP = LS_PRCP + PTR2D
         Q = Q - TMP3D
         T = T + (MAPL_ALHL / MAPL_CP) * TMP3D

         DTDT_ER = (T - DTDT_ER) / DT_MOIST
         DQVDT_ER = (Q - DQVDT_ER) / DT_MOIST

         ! Export Total Moist Tendencies

         call MAPL_GetPointer(EXPORT, DUDT, 'DUDT', _RC)
         if (associated(DUDT)) then
            DUDT = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'DUDT_DC', _RC)
            if (associated(PTR3D)) DUDT = DUDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DUDT_SC', _RC)
            if (associated(PTR3D)) DUDT = DUDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DUDT_macro', _RC)
            if (associated(PTR3D)) DUDT = DUDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DUDT_micro', _RC)
            if (associated(PTR3D)) DUDT = DUDT + PTR3D
         end if

         call MAPL_GetPointer(EXPORT, DVDT, 'DVDT', _RC)
         if (associated(DVDT)) then
            DVDT = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'DVDT_DC', _RC)
            if (associated(PTR3D)) DVDT = DVDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DVDT_SC', _RC)
            if (associated(PTR3D)) DVDT = DVDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DVDT_macro', _RC)
            if (associated(PTR3D)) DVDT = DVDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DVDT_micro', _RC)
            if (associated(PTR3D)) DVDT = DVDT + PTR3D
         end if

         call MAPL_GetPointer(EXPORT, DTDT, 'DTDT', _RC)
         if (associated(DTDT)) then
            DTDT = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_DC', _RC)
            if (associated(PTR3D)) DTDT = DTDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_SC', _RC)
            if (associated(PTR3D)) DTDT = DTDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_macro', _RC)
            if (associated(PTR3D)) DTDT = DTDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_micro', _RC)
            if (associated(PTR3D)) DTDT = DTDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_ER', _RC)
            if (associated(PTR3D)) DTDT = DTDT + PTR3D
            DTDT = DTDT * (PLE(:, :, 1:LM) - PLE(:, :, 0:LM - 1)) ! Pressure weighted tendency
         end if

         call MAPL_GetPointer(EXPORT, DQDT, 'DQDT', _RC)
         if (associated(DQDT)) then
            DQDT = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_DC', _RC)
            if (associated(PTR3D)) DQDT = DQDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_SC', _RC)
            if (associated(PTR3D)) DQDT = DQDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_macro', _RC)
            if (associated(PTR3D)) DQDT = DQDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_micro', _RC)
            if (associated(PTR3D)) DQDT = DQDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_ER', _RC)
            if (associated(PTR3D)) DQDT = DQDT + PTR3D
         end if

         call MAPL_GetPointer(EXPORT, DQLDT, 'DQLDT', _RC)
         if (associated(DQLDT)) then
            DQLDT = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_DC', _RC)
            if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_SC', _RC)
            if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_macro', _RC)
            if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_micro', _RC)
            if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
         end if

         call MAPL_GetPointer(EXPORT, DQIDT, 'DQIDT', _RC)
         if (associated(DQIDT)) then
            DQIDT = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_DC', _RC)
            if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_SC', _RC)
            if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_macro', _RC)
            if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_micro', _RC)
            if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
         end if

         call MAPL_GetPointer(EXPORT, DQRDT, 'DQRDT', _RC)
         if (associated(DQRDT)) then
            DQRDT = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQRDT_SC', _RC)
            if (associated(PTR3D)) DQRDT = DQRDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQRDT_macro', _RC)
            if (associated(PTR3D)) DQRDT = DQRDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQRDT_micro', _RC)
            if (associated(PTR3D)) DQRDT = DQRDT + PTR3D
         end if

         call MAPL_GetPointer(EXPORT, DQSDT, 'DQSDT', _RC)
         if (associated(DQSDT)) then
            DQSDT = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQSDT_SC', _RC)
            if (associated(PTR3D)) DQSDT = DQSDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQSDT_macro', _RC)
            if (associated(PTR3D)) DQSDT = DQSDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQSDT_micro', _RC)
            if (associated(PTR3D)) DQSDT = DQSDT + PTR3D
         end if

         call MAPL_GetPointer(EXPORT, DQGDT, 'DQGDT', _RC)
         if (associated(DQGDT)) then
            DQGDT = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQGDT_macro', _RC)
            if (associated(PTR3D)) DQGDT = DQGDT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQGDT_micro', _RC)
            if (associated(PTR3D)) DQGDT = DQGDT + PTR3D
         end if

         call MAPL_GetPointer(EXPORT, DQADT, 'DQADT', _RC)
         if (associated(DQADT)) then
            DQADT = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQADT_DC', _RC)
            if (associated(PTR3D)) DQADT = DQADT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQADT_SC', _RC)
            if (associated(PTR3D)) DQADT = DQADT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQADT_macro', _RC)
            if (associated(PTR3D)) DQADT = DQADT + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'DQADT_micro', _RC)
            if (associated(PTR3D)) DQADT = DQADT + PTR3D
         end if

         call MAPL_GetPointer(EXPORT, DPDTMST, 'DPDTMST', _RC)
         if (associated(DPDTMST)) then
            DPDTMST = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_CN', _RC)
            if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
            call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_SC', _RC)
            if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
            call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_AN', _RC)
            if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
            call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_LS', _RC)
            if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
            call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_CN', _RC)
            if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
            call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_SC', _RC)
            if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
            call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_AN', _RC)
            if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
            call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_LS', _RC)
            if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
            DPDTMST = MAPL_GRAV * DPDTMST
         end if

         call MAPL_GetPointer(EXPORT, PFL_LSAN, 'PFL_LSAN', _RC)
         if (associated(PFL_LSAN)) then
            PFL_LSAN = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_AN', _RC)
            if (associated(PTR3D)) PFL_LSAN = PFL_LSAN + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_LS', _RC)
            if (associated(PTR3D)) PFL_LSAN = PFL_LSAN + PTR3D
         end if

         call MAPL_GetPointer(EXPORT, PFI_LSAN, 'PFI_LSAN', _RC)
         if (associated(PFI_LSAN)) then
            PFI_LSAN = 0.0
            call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_AN', _RC)
            if (associated(PTR3D)) PFI_LSAN = PFI_LSAN + PTR3D
            call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_LS', _RC)
            if (associated(PTR3D)) PFI_LSAN = PFI_LSAN + PTR3D
         end if

         ! Combine Precip Exports

         ! liquid convective precip
         call MAPL_GetPointer(EXPORT, PCU, 'PCU', _RC)
         if (associated(PCU)) then
            PCU = 0.0
            call MAPL_GetPointer(EXPORT, PTR2D, 'CN_PRCP', _RC)
            if (associated(PTR2D)) PCU = PCU + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'SC_PRCP', _RC)
            if (associated(PTR2D)) PCU = PCU + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'CNPCPRATE', _RC)
            if (associated(PTR2D)) PCU = PCU + PTR2D
            PCU = max(PCU, 0.0)
         end if

         ! liquid large-scale precip
         call MAPL_GetPointer(EXPORT, PLS, 'PLS', _RC)
         if (associated(PLS)) then
            PLS = 0.0
            call MAPL_GetPointer(EXPORT, PTR2D, 'LS_PRCP', _RC)
            if (associated(PTR2D)) PLS = PLS + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'AN_PRCP', _RC)
            if (associated(PTR2D)) PLS = PLS + PTR2D
            PLS = max(PLS, 0.0)
         end if

         ! all liquid precip
         call MAPL_GetPointer(EXPORT, RAIN, 'RAIN', ALLOC=.true., _RC)
         if (associated(RAIN)) then
            RAIN = 0.0
            call MAPL_GetPointer(EXPORT, PTR2D, 'LS_PRCP', _RC)
            if (associated(PTR2D)) RAIN = RAIN + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'AN_PRCP', _RC)
            if (associated(PTR2D)) RAIN = RAIN + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'CN_PRCP', _RC)
            if (associated(PTR2D)) RAIN = RAIN + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'SC_PRCP', _RC)
            if (associated(PTR2D)) RAIN = RAIN + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'CNPCPRATE', _RC)
            if (associated(PTR2D)) RAIN = RAIN + PTR2D
            RAIN = max(RAIN, 0.0)
         end if

         ! all frozen precip (snow at this point)
         call MAPL_GetPointer(EXPORT, SNOW, 'SNO', ALLOC=.true., _RC)
         if (associated(SNOW)) then
            SNOW = 0.0
            call MAPL_GetPointer(EXPORT, PTR2D, 'LS_SNR', _RC)
            if (associated(PTR2D)) SNOW = SNOW + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'AN_SNR', _RC)
            if (associated(PTR2D)) SNOW = SNOW + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'CN_SNR', _RC)
            if (associated(PTR2D)) SNOW = SNOW + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'SC_SNR', _RC)
            if (associated(PTR2D)) SNOW = SNOW + PTR2D
            SNOW = max(SNOW, 0.0)
         end if

         ! all deep convective precip (rain+snow)
         call MAPL_GetPointer(EXPORT, CN_PRCP, 'CN_PRCP', _RC)
         if (associated(CN_PRCP)) then
            call MAPL_GetPointer(EXPORT, PTR2D, 'CNPCPRATE', _RC)
            if (associated(PTR2D)) CN_PRCP = CN_PRCP + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'CN_SNR', _RC)
            if (associated(PTR2D)) CN_PRCP = CN_PRCP + PTR2D
            CN_PRCP = max(CN_PRCP, 0.0)
         end if

         ! all large-scale precip (rain+snow)
         call MAPL_GetPointer(EXPORT, LS_PRCP, 'LS_PRCP', _RC)
         if (associated(LS_PRCP)) then
            call MAPL_GetPointer(EXPORT, PTR2D, 'LS_SNR', _RC)
            if (associated(PTR2D)) LS_PRCP = LS_PRCP + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'ICE', _RC)
            if (associated(PTR2D)) LS_PRCP = LS_PRCP + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'FRZR', _RC)
            if (associated(PTR2D)) LS_PRCP = LS_PRCP + PTR2D
            LS_PRCP = max(LS_PRCP, 0.0)
         end if

         ! all anvil precip (rain+snow)
         call MAPL_GetPointer(EXPORT, AN_PRCP, 'AN_PRCP', _RC)
         if (associated(AN_PRCP)) then
            call MAPL_GetPointer(EXPORT, PTR2D, 'AN_SNR', _RC)
            if (associated(PTR2D)) AN_PRCP = AN_PRCP + PTR2D
         end if

         ! all shallow precip (rain+snow)
         call MAPL_GetPointer(EXPORT, SC_PRCP, 'SC_PRCP', _RC)
         if (associated(SC_PRCP)) then
            call MAPL_GetPointer(EXPORT, PTR2D, 'SC_SNR', _RC)
            if (associated(PTR2D)) SC_PRCP = SC_PRCP + PTR2D
            SC_PRCP = max(SC_PRCP, 0.0)
         end if

         ! Total - all precip (rain+snow)
         call MAPL_GetPointer(EXPORT, TPREC, 'TPREC', ALLOC=.true., _RC)
         if (associated(TPREC)) then
            TPREC = 0.0
            call MAPL_GetPointer(EXPORT, PTR2D, 'LS_PRCP', _RC)
            if (associated(PTR2D)) TPREC = TPREC + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'AN_PRCP', _RC)
            if (associated(PTR2D)) TPREC = TPREC + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'CN_PRCP', _RC)
            if (associated(PTR2D)) TPREC = TPREC + PTR2D
            call MAPL_GetPointer(EXPORT, PTR2D, 'SC_PRCP', _RC)
            if (associated(PTR2D)) TPREC = TPREC + PTR2D
            TPREC = max(TPREC, 0.0)
         end if

         ! diagnosed stratiform precip (rain+snow)
         call MAPL_GetPointer(EXPORT, PREC_STRAT, 'PREC_STRAT', ALLOC=.true., _RC)
         if (associated(PREC_STRAT)) PREC_STRAT = max((1.0 - CNV_FRC) * TPREC, 0.0)

         ! diagnosed convective precip (rain+snow)
         call MAPL_GetPointer(EXPORT, PREC_CONV, 'PREC_CONV', ALLOC=.true., _RC)
         if (associated(PREC_CONV)) PREC_CONV = max(CNV_FRC * TPREC, 0.0)

         ! Diagnostic precip types:
         call MAPL_GetPointer(EXPORT, ICE, 'ICE', ALLOC=.true., _RC)
         call MAPL_GetPointer(EXPORT, FRZR, 'FRZR', ALLOC=.true., _RC)
         if (LUPDATE_PRECIP_TYPE .or. LDIAGNOSE_PRECIP_TYPE) then
            call MAPL_GetPointer(EXPORT, PTYPE, 'PTYPE', ALLOC=.true., _RC)
            call DIAGNOSE_PRECIP_TYPE(IM, JM, LM, TPREC, PLS, PCU, RAIN, SNOW, ICE, FRZR, &
                 PTYPE, PLE, T / PK, PK, PKE, ZL0, LUPDATE_PRECIP_TYPE)
         end if
         ICE = max(ICE, 0.0)
         FRZR = max(FRZR, 0.0)
         ! Get Kuchera snow:rain ratios
         do I = 1, IM
            do J = 1, JM
               Tmax = 0.0
               do L = LM, 1, -1
                  if (PLmb(I, J, L) > 500.) then
                     Tmax = max(Tmax, T(I, J, L))
                  end if
               end do
               if (Tmax <= 271.16) then
                  TMP2D(I, J) = 12.0 + (271.16 - Tmax)
               else
                  TMP2D(I, J) = 12.0 + 2 * (271.16 - Tmax)
               end if
               TMP2D(I, J) = max(0.0, TMP2D(I, J))
            end do
         end do
         call MAPL_GetPointer(EXPORT, PTR2D, 'KUCHERA_RATIO', _RC)
         if (associated(PTR2D)) PTR2D = TMP2D
         ! Accumulated precip totals (mm), apply KUCHERA_RATIO for SNOW
         call MAPL_GetPointer(EXPORT, PTR2D, 'SNOWTOTAL', _RC)
         if (associated(PTR2D)) PTR2D = TMP2D * DT_MOIST * (SNOW + ICE)
         call MAPL_GetPointer(EXPORT, PTR2D, 'PRECTOTAL', _RC)
         if (associated(PTR2D)) PTR2D = DT_MOIST * TPREC

         call MAPL_GetPointer(EXPORT, PTR3D, 'QLTOT', _RC)
         if (associated(PTR3D)) PTR3D = QLLS + QLCN

         call MAPL_GetPointer(EXPORT, PTR3D, 'QITOT', _RC)
         if (associated(PTR3D)) PTR3D = QILS + QICN

         call MAPL_GetPointer(EXPORT, PTR3D, 'QCTOT', _RC)
         if (associated(PTR3D)) then
            PTR3D = QLLS + QLCN + QILS + QICN
            call MAPL_GetPointer(EXPORT, QRTOT, 'QRTOT', _RC)
            if (associated(QRTOT)) PTR3D = PTR3D + QRTOT
            call MAPL_GetPointer(EXPORT, QSTOT, 'QSTOT', _RC)
            if (associated(QSTOT)) PTR3D = PTR3D + QSTOT
            call MAPL_GetPointer(EXPORT, QGTOT, 'QGTOT', _RC)
            if (associated(QGTOT)) PTR3D = PTR3D + QGTOT
         end if

         ! Cloud condensate exports
         call MAPL_GetPointer(EXPORT, PTR3D, 'QLLSX1', _RC)
         if (associated(PTR3D)) PTR3D = QLLS

         call MAPL_GetPointer(EXPORT, PTR3D, 'QILSX1', _RC)
         if (associated(PTR3D)) PTR3D = QILS

         call MAPL_GetPointer(EXPORT, PTR3D, 'QLCNX1', _RC)
         if (associated(PTR3D)) PTR3D = QLCN

         call MAPL_GetPointer(EXPORT, PTR3D, 'QICNX1', _RC)
         if (associated(PTR3D)) PTR3D = QICN

         call MAPL_GetPointer(EXPORT, PTR3D, 'CLLSX1', _RC)
         if (associated(PTR3D)) PTR3D = CLLS

         call MAPL_GetPointer(EXPORT, PTR3D, 'CLCNX1', _RC)
         if (associated(PTR3D)) PTR3D = CLCN

         ! Fill wind, temperature & RH exports needed for SYNCTQ

         call MAPL_GetPointer(EXPORT, PTR3D, 'UAFMOIST', ALLOC=.true., _RC)
         if (associated(PTR3D)) PTR3D = U

         call MAPL_GetPointer(EXPORT, PTR3D, 'VAFMOIST', ALLOC=.true., _RC)
         if (associated(PTR3D)) PTR3D = V

         call MAPL_GetPointer(EXPORT, PTR3D, 'TAFMOIST', ALLOC=.true., _RC)
         if (associated(PTR3D)) PTR3D = T

         call MAPL_GetPointer(EXPORT, PTR3D, 'QAFMOIST', ALLOC=.true., _RC)
         if (associated(PTR3D)) PTR3D = Q

         call MAPL_GetPointer(EXPORT, PTR3D, 'THAFMOIST', ALLOC=.true., _RC)
         if (associated(PTR3D)) PTR3D = T / PK

         call MAPL_GetPointer(EXPORT, PTR3D, 'SAFMOIST', ALLOC=.true., _RC)
         if (associated(PTR3D)) then
            do L = 1, LM
               PTR3D(:, :, L) = MAPL_CP * T(:, :, L) + MAPL_GRAV * (ZL0(:, :, L) + ZLE(:, :, LM))
            end do
         end if

         call MAPL_GetPointer(EXPORT, PTR3D, 'RH2', _RC)
         if (associated(PTR3D)) PTR3D = max(min(Q / GEOS_QSAT(T, PLmb), 1.02), 0.0)

         ! Cloud Water Path
         call MAPL_GetPointer(EXPORT, PTR2D, 'CWP', _RC)
         if (associated(PTR2D)) PTR2D = SUM((QLCN + QLLS + QICN + QILS) * MASS, 3)
         ! Cloud Liquid Water Path
         call MAPL_GetPointer(EXPORT, PTR2D, 'CLWP', _RC)
         if (associated(PTR2D)) PTR2D = SUM((QLCN + QLLS) * MASS, 3)
         ! Cloud Ice Water Path
         call MAPL_GetPointer(EXPORT, PTR2D, 'CIWP', _RC)
         if (associated(PTR2D)) PTR2D = SUM((QICN + QILS) * MASS, 3)
         ! Total LWP and IWP need to be defined inside the microphysics depends on number was water species

         ! Total Precipitable Water
         call MAPL_GetPointer(EXPORT, PTR2D, 'TPW', _RC)
         if (associated(PTR2D)) PTR2D = SUM((Q) * MASS, 3)

         ! Lightning Exports
         call MAPL_GetPointer(EXPORT, PTR2D, 'LFR_GCC', NotFoundOk=.true., _RC)
         if (associated(PTR2D)) PTR2D = 0.0

         call MAPL_TimerOff(MAPL, "---MOIST_EPILOGUE")

      else ! Alarm ringing

         ! Internal State
         call MAPL_GetPointer(INTERNAL, Q, 'Q', _RC)
         ! Import State
         call MAPL_GetPointer(IMPORT, PLE, 'PLE', _RC)
         call MAPL_GetPointer(IMPORT, ZLE, 'ZLE', _RC)
         call MAPL_GetPointer(IMPORT, U, 'U', _RC)
         call MAPL_GetPointer(IMPORT, V, 'V', _RC)
         call MAPL_GetPointer(IMPORT, T, 'T', _RC)
         ! Allocatables
         ! Edge variables
         allocate(ZLE0(IM, JM, 0:LM))
         allocate(PLEmb(IM, JM, 0:LM))
         ! Layer variables
         allocate(ZL0(IM, JM, LM))
         allocate(PLmb(IM, JM, LM))
         allocate(PK(IM, JM, LM))
         allocate(MASS(IM, JM, LM))
         allocate(TMP2D(IM, JM))
         ! dervied states
         MASS = (PLE(:, :, 1:LM) - PLE(:, :, 0:LM - 1)) / MAPL_GRAV
         PLEmb = PLE * .01
         PLmb = 0.5 * (PLEmb(:, :, 0:LM - 1) + PLEmb(:, :, 1:LM))
         PK = (100.0 * PLmb / MAPL_P00)**(MAPL_KAPPA)
         do L = 0, LM
            ZLE0(:, :, L) = ZLE(:, :, L) - ZLE(:, :, LM) ! Edge Height (m) above the surface
         end do
         ZL0 = 0.5 * (ZLE0(:, :, 0:LM - 1) + ZLE0(:, :, 1:LM)) ! Layer Height (m) above the surface

         ! Fill Wind, Temperature & RH exports needed for SYNCTQ

         call MAPL_GetPointer(EXPORT, PTR3D, 'UAFMOIST', ALLOC=.true., _RC)
         if (associated(PTR3D)) PTR3D = U

         call MAPL_GetPointer(EXPORT, PTR3D, 'VAFMOIST', ALLOC=.true., _RC)
         if (associated(PTR3D)) PTR3D = V

         call MAPL_GetPointer(EXPORT, PTR3D, 'TAFMOIST', ALLOC=.true., _RC)
         if (associated(PTR3D)) PTR3D = T

         call MAPL_GetPointer(EXPORT, PTR3D, 'QAFMOIST', ALLOC=.true., _RC)
         if (associated(PTR3D)) PTR3D = Q

         call MAPL_GetPointer(EXPORT, PTR3D, 'THAFMOIST', ALLOC=.true., _RC)
         if (associated(PTR3D)) PTR3D = T / PK

         call MAPL_GetPointer(EXPORT, PTR3D, 'SAFMOIST', ALLOC=.true., _RC)
         if (associated(PTR3D)) then
            do L = 1, LM
               PTR3D(:, :, L) = MAPL_CP * T(:, :, L) + MAPL_GRAV * (ZL0(:, :, L) + ZLE(:, :, LM))
            end do
         end if

         call MAPL_GetPointer(EXPORT, PTR3D, 'RH2', _RC)
         if (associated(PTR3D)) PTR3D = max(min(Q / GEOS_QSAT(T, PLmb), 1.02), 0.0)

      end if ! Alarm

      call MAPL_TimerOff(MAPL, "TOTAL")

      _RETURN(_SUCCESS)

   end subroutine Run

   subroutine Finalize(GC, IMPORT, EXPORT, CLOCK, RC)
      type(ESMF_GridComp), intent(inout) :: GC ! Gridded component
      type(ESMF_State), intent(inout) :: IMPORT ! Import state
      type(ESMF_State), intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock), intent(inout) :: CLOCK ! The clock
      integer, optional, intent(out) :: RC ! Error code

      ! ErrLog variables
      integer :: STATUS
      ! Begin...

      if (adjustl(CLDMICR_OPTION) == "GFDL_1M") call GFDL_1M_FINALIZE(GC, IMPORT, EXPORT, _RC)
      if (adjustl(SHALLOW_OPTION) == "UW") call UW_Finalize(GC, IMPORT, EXPORT, _RC)
      if (adjustl(CONVPAR_OPTION) == "GF") call GF_Finalize(GC, IMPORT, EXPORT, _RC)

      ! Call Finalize for every child
      call MAPL_GenericFinalize(GC, IMPORT, EXPORT, CLOCK, _RC)
      ! End
      _RETURN(_SUCCESS)

   end subroutine Finalize

end module GEOS_MoistGridCompMod
