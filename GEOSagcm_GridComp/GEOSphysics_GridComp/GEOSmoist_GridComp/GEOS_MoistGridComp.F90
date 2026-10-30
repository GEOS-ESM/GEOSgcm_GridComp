!#define UWDIAG 1
!#define PDFDIAG 1

#include "MAPL.h"

!BOP

!MODULE:
! GEOS_Moist -- A Module to compute moist processes, including convection,
! large-scale condensation and precipitation and cloud parameters.

!INTERFACE:

module GEOS_MoistGridCompMod

   !USES:

   use ESMF
   use MAPL
   use pflogger, only: logger_t => logger
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

   ! MAPL2's MAPL_UNDEF (MAPL3's is -huge), for SNOMAS: its undefined points come from
   ! Surface or from inputs converted from MAPL2
   real, parameter :: MAPL2_UNDEF = 1.0e15

   ! COND for the PDF_* exports in Moist_StateSpecs.rc (ACG has no #ifdef)
#ifdef PDFDIAG
   logical, parameter :: PDF_DIAG = .true.
#else
   logical, parameter :: PDF_DIAG = .false.
#endif

   !PUBLIC MEMBER FUNCTIONS:

   public SetServices

   !DESCRIPTION:
   ! {\tt GEOS\_MoistGridCompMod} implements moist processes in GEOS-5. These
   ! include all processes that involve phase changes in the atmosphere, such
   ! as large-scale condensation, convective clouds, and all rain and cloud
   ! formation. Its state consists of water vapor, various types of condensate,
   ! and fractions of various cloud types.
   ! Two moment cloud microphysics (Barahona et al., GMD, 2014.) can be run by
   ! setting CLDMACRO==MGB2_2M. When using 2-moment microphysics the number
   ! concentration of ice crystals and cloud droplets are part of the state.

   !EOP

contains

   !BOP
   !IROUTINE: SetServices -- Sets ESMF services for this component
   !INTERFACE:
   subroutine SetServices(GC, RC)
      !ARGUMENTS:
      type(ESMF_GridComp) :: GC ! gridded component
      integer, intent(out) :: RC ! return code

      !DESCRIPTION:
      ! {\tt GEOS\_MoistGridCompMod} uses the default Initialize and Finalize
      ! services, but registers its own Run method.

      !EOP

      integer :: STATUS

      ! Local derived type aliases

      type(ESMF_Config) :: CF ! only passed to the unported schemes' *_Setup
      type(StringVector) :: ADVECTED

      logical :: LCONVPAR
      logical :: LSHALLOW
      logical :: LCLDMICR

      ! Begin...

      ! Set the Init entry point
      call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_INITIALIZE, Initialize, _RC)

      ! Set the Run entry point
      call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN, Run, phase_name='run', _RC)

      ! Set the Finalize entry point
      call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_FINALIZE, Finalize, _RC)

      ! Inititialize deep convective parameterizations (Options: RAS, GF or NONE)
      call MAPL_GridCompGetResource(GC, 'CONVPAR_OPTION', CONVPAR_OPTION, default="GF", _RC)
      LCONVPAR = adjustl(CONVPAR_OPTION) == "RAS" .or. &
           adjustl(CONVPAR_OPTION) == "GF" .or. &
           adjustl(CONVPAR_OPTION) == "NONE"
      _ASSERT(LCONVPAR, 'Unsupported Deep Convection Option')
      ! default=(X): a copy of the Process_Library default (value and default must not alias)
      call MAPL_GridCompGetResource(GC, 'SIGMA_DX', SIGMA_DX, default=(SIGMA_DX), _RC)
      call MAPL_GridCompGetResource(GC, 'SIGMA_EXP', SIGMA_EXP, default=(SIGMA_EXP), _RC)

      ! Inititialize shallow convective parameterizations (Options: UW or NONE)
      call MAPL_GridCompGetResource(GC, 'SHALLOW_OPTION', SHALLOW_OPTION, default="UW", _RC)
      LSHALLOW = adjustl(SHALLOW_OPTION) == "UW" .or. &
           adjustl(SHALLOW_OPTION) == "NONE"
      _ASSERT(LSHALLOW, 'Unsupported Shallow Convection Option')

      ! Inititialize cloud microphysics (Options: BACM_1M, GFDL_1M, THOM_1M or MGB2_2M)
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
      ! Note: It's possible these aren't *used* if USE_AEROSOL_NN=.TRUE. but they are still passed
      ! in so they have to be defined
      call MAPL_GridCompGetResource(GC, 'NCCN_OCN', CCN_OCN, default=50., _RC)
      call MAPL_GridCompGetResource(GC, 'NCCN_LND', CCN_LND, default=300., _RC)

      call MAPL_GridCompGetResource(GC, 'NCIN_OCN', CIN_OCN, default=0.005, _RC)
      call MAPL_GridCompGetResource(GC, 'NCIN_LND', CIN_LND, default=0.05, _RC)

      ! NOTE: Binary restarts expect Q to be the first field in the moist_internal_rst. Thus,
      ! the first MAPL_AddInternalSpec call must be from the microphysics

      ! MAPL3: only GFDL_1M, GF and UW are ported so far. The other schemes
      ! fail here, before any of their code runs.

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

      !IMPORT STATE:
#include "Moist_Import___.h"

      call MAPL_GridCompAddSpec(GC, &
           state_intent=ESMF_STATEINTENT_IMPORT, &
           short_name='AERO', &
           standard_name='aerosols', &
           itemtype=MAPL_STATEITEM_STATE, _RC)

      !EXPORT STATE:
#include "Moist_Export___.h"

      ! INTERNAL STATE: the ALLOC'd exports that need storage beyond one routine call,
      ! re-exported (ADD2EXPORT). After the scheme Setups, so Q stays the first internal.
#include "Moist_Internal___.h"

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

   !BOP
   !IROUTINE: Initialize -- Initialize method for the composite Moist Gridded Component
   !INTERFACE:
   subroutine Initialize(GC, IMPORT, EXPORT, CLOCK, RC)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: GC ! Gridded component
      type(ESMF_State), intent(inout) :: IMPORT ! Import state
      type(ESMF_State), intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock), intent(inout) :: CLOCK ! The clock
      integer, optional, intent(out) :: RC ! Error code

      !DESCRIPTION: The Initialize method of the Moist Physics Gridded Component first
      ! calls the Initialize method of the child Dynamics.  The Dynamics Initialize method will
      ! create the ESMF GRID, which will then be used to set the GRID associated with the
      ! SuperDyn Composite Component itself.  It should be noted that the
      ! SuperDyn Initialize method also invokes the GEOS Topo Utility which creates all
      ! topography related quantities.

      !EOP

      integer :: STATUS

      ! Local derived type aliases

      type(ESMF_TimeInterval) :: TINT
      real :: DT_MOIST
      real :: DBZ_DT
      type(ESMF_Calendar) :: calendar
      type(ESMF_Time) :: currTime
      type(ESMF_Alarm) :: DBZ_RunAlarm
      type(ESMF_TimeInterval) :: ringInterval

      ! Begin...

      ! Get parameters from the component
      call MAPL_GridCompGetResource(GC, 'DIAGNOSE_PRECIP_TYPE', LDIAGNOSE_PRECIP_TYPE, default=.false., _RC)
      call MAPL_GridCompGetResource(GC, 'UPDATE_PRECIP_TYPE', LUPDATE_PRECIP_TYPE, default=.false., _RC)

      ! default=(X): a copy (value is intent(inout) and must not alias default)
      DETRAIN_INACTIVE_CNV = 0.0
      if (adjustl(CLDMICR_OPTION) == "GFDL_1M") DETRAIN_INACTIVE_CNV = 5.0e-4
      call MAPL_GridCompGetResource(GC, 'DETRAIN_INACTIVE_CNV', DETRAIN_INACTIVE_CNV, &
           default=(DETRAIN_INACTIVE_CNV), _RC)
      call MAPL_GridCompGetResource(GC, 'TAU_DETRAIN_CNV', TAU_DETRAIN_CNV, default=3600.0, _RC)

      ! MAPL3: the unported schemes still take MAPL2's MetaComp; SetServices already
      ! fails for them, so these branches are unreachable until they are ported (step 11)
      if (adjustl(CONVPAR_OPTION) == "RAS") then
         _FAIL('CONVPAR_OPTION RAS is not yet ported to MAPL3') ! was RAS_Initialize(MAPL, _RC)
      end if
      if (adjustl(CONVPAR_OPTION) == "GF") call GF_Initialize(GC, CLOCK, IMPORT, EXPORT, _RC)
      if (adjustl(SHALLOW_OPTION) == "UW") call UW_Initialize(GC, CLOCK, IMPORT, EXPORT, _RC)
      if (adjustl(CLDMICR_OPTION) == "BACM_1M") then
         _FAIL('CLDMICR_OPTION BACM_1M is not yet ported to MAPL3') ! was BACM_1M_Initialize(MAPL, _RC)
      end if
      if (adjustl(CLDMICR_OPTION) == "GFDL_1M") call GFDL_1M_Initialize(GC, CLOCK, IMPORT, EXPORT, _RC)
      if (adjustl(CLDMICR_OPTION) == "THOM_1M") then
         _FAIL('CLDMICR_OPTION THOM_1M is not yet ported to MAPL3') ! was THOM_1M_Initialize(MAPL, _RC)
      end if
      if (adjustl(CLDMICR_OPTION) == "MGB2_2M") then
         _FAIL('CLDMICR_OPTION MGB2_2M is not yet ported to MAPL3') ! was MGB2_2M_Initialize(MAPL, _RC)
      end if

      ! Moist's timestep (MAPL2: the RUNALARM ring interval). MAPL3 gives Moist its own
      ! clock, advanced by its own timestep.
      call ESMF_ClockGet(CLOCK, timeStep=TINT, _RC)
      call MAPL_ClockGet(CLOCK, dt=DT_MOIST, _RC)

      DBZ_DT = max(DT_MOIST, 900.0)
      call MAPL_GridCompGetResource(GC, 'DBZ_DT', DBZ_DT, default=(DBZ_DT), _RC)
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
      call MAPL_GridCompGetResource(GC, 'refl10cm_allow_wet_graupel', refl10cm_allow_wet_graupel, default=.false., _RC)
      call MAPL_GridCompGetResource(GC, 'refl10cm_allow_wet_snow', refl10cm_allow_wet_snow, default=.false., _RC)
      call MAPL_GridCompGetResource(GC, 'DBZ_VAR_INTERCP', DBZ_VAR_INTERCP, default=(DBZ_VAR_INTERCP), _RC)
      call MAPL_GridCompGetResource(GC, 'LIQUID_SKIN_SNOW', LIQUID_SKIN_SNOW, default=.false., _RC)
      call MAPL_GridCompGetResource(GC, 'LIQUID_SKIN_GRAUPEL', LIQUID_SKIN_GRAUPEL, default=.false., _RC)
      call MAPL_GridCompGetResource(GC, 'LIQUID_SKIN_HAIL', LIQUID_SKIN_HAIL, default=.true., _RC)

      ! All done

      _RETURN(_SUCCESS)
   end subroutine Initialize

   !BOP
   !IROUTINE: RUN -- Run method for the CONVECT component
   !INTERFACE:
   subroutine Run(GC, IMPORT, EXPORT, CLOCK, RC)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: GC ! Gridded component
      type(ESMF_State), intent(inout) :: IMPORT ! Import state
      type(ESMF_State), intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock), intent(inout) :: CLOCK ! The clock
      integer, optional, intent(out) :: RC ! Error code:

      !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
      ! the Initialize and Finalize services, as well as allocating

      !EOP

      integer :: STATUS

      type(ESMF_VM) :: VMG

      ! Local derived type aliases

      type(ESMF_Grid) :: GRID
      type(ESMF_State) :: INTERNAL
      class(logger_t), pointer :: logger
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
      type(ESMF_FieldBundle) :: TR, TR_IMPORT
      type(ESMF_FieldBundle) :: TR_CNV ! MTR's fields plus GF's CNV_TR, for CNV_Tracers_Init
      type(ESMF_Field), allocatable :: TR_FIELDS(:)
      type(ESMF_Field) :: CNV_TR_FIELD
      integer :: TR_COUNT
      logical :: same_mtr_data
      ! Re-exported internals (were ALLOC'd exports shared with the schemes)
      real, pointer, dimension(:, :) :: CN_PRCP, LS_PRCP, AN_PRCP, SC_PRCP, ICE, FRZR
      real, pointer, dimension(:, :, :) :: BYNCY
      real, pointer, dimension(:, :) :: CNV_FRC, SRF_TYPE
      real, pointer, dimension(:, :) :: EIS, LTS
      ! Exports computed in local arrays (were ALLOC'd), copied to the export if it is connected
      real, allocatable, dimension(:, :, :) :: DTDT_ER, DQVDT_ER
      real, allocatable, dimension(:, :) :: PTYPE, TPREC, RAIN, SNOW, ER_PRCP, PLS, PCU
      real, allocatable, dimension(:, :) :: CAPE, INHB, SBCAPE, SBCIN, LFC, LNB, LCL_AGL
      ! Pointers, not allocatables: BUOYANCY2 takes them as pointer dummies (Process_Library.F90)
      real, pointer, dimension(:, :) :: MLCAPE, MLCIN, MUCAPE, MUCIN
      real, allocatable, dimension(:, :, :) :: CFICE
      real, allocatable, dimension(:, :, :) :: NWFA, NIFA
      ! Their export pointers
      real, pointer, dimension(:, :, :) :: DTDT_ER_ptr, DQVDT_ER_ptr
      real, pointer, dimension(:, :) :: PTYPE_ptr, TPREC_ptr, RAIN_ptr, SNOW_ptr, ER_PRCP_ptr, PLS_ptr, PCU_ptr
      real, pointer, dimension(:, :) :: CAPE_ptr, INHB_ptr, SBCAPE_ptr, SBCIN_ptr, LFC_ptr, LNB_ptr, LCL_AGL_ptr
      real, pointer, dimension(:, :) :: MLCAPE_ptr, MLCIN_ptr, MUCAPE_ptr, MUCIN_ptr
      real, pointer, dimension(:, :, :) :: CFICE_ptr
      real, pointer, dimension(:, :, :) :: NWFA_ptr, NIFA_ptr
      ! Exports
      real, pointer, dimension(:, :, :) :: DQDT, DQADT, DQIDT, DQLDT, DQRDT, DQSDT, DQGDT
      real, pointer, dimension(:, :, :) :: DTDT, DUDT, DVDT, DWDT
      real, pointer, dimension(:, :, :) :: DPDTMST, PFL_LSAN, PFI_LSAN
      real, pointer, dimension(:, :) :: PREC_STRAT, PREC_CONV
      real, pointer, dimension(:, :, :) :: CFLIQ
      real, pointer, dimension(:, :, :) :: PTRDC, PTRSC
      real, pointer, contiguous, dimension(:, :, :) :: p3d ! scratch for the edge-field bounds remap
      real, pointer, dimension(:, :, :) :: PTR3D
      real, pointer, dimension(:, :) :: PTR2D
      real, pointer, dimension(:) :: PTR1D

      integer :: IM, JM, LM
      integer :: I, J, L, n

      ! Begin...

      ! Get my VM

      call ESMF_GridCompGet(GC, VM=VMG, _RC)

      call MAPL_GridCompTimerStart(GC, "TOTAL", _RC)

      call MAPL_GridCompGet(GC, num_levels=LM, grid=GRID, logger=logger, _RC)
      call MAPL_GridGet(GRID, IM=IM, JM=JM, _RC)
      call MAPL_GridCompGetInternalState(GC, INTERNAL, _RC)

      ! Moist's timestep (MAPL2: the RUNALARM ring interval). MAPL3 calls Run only when
      ! Moist's own run alarm rings (MAPL2 called it every heartbeat and Run checked
      ! RUNALARM; its not-ringing branch only refilled the *AFMOIST and RH2 exports).
      call MAPL_ClockGet(CLOCK, dt=DT_MOIST, _RC)

      call MAPL_GridCompTimerStart(GC, "---MOIST_PROLOGUE", _RC)

      ! Internal State
      call MAPL_StateGetPointer(INTERNAL, Q, 'Q', _RC)
      call MAPL_StateGetPointer(INTERNAL, QLLS, 'QLLS', _RC)
      call MAPL_StateGetPointer(INTERNAL, QLCN, 'QLCN', _RC)
      call MAPL_StateGetPointer(INTERNAL, CLCN, 'CLCN', _RC)
      call MAPL_StateGetPointer(INTERNAL, CLLS, 'CLLS', _RC)
      call MAPL_StateGetPointer(INTERNAL, QILS, 'QILS', _RC)
      call MAPL_StateGetPointer(INTERNAL, QICN, 'QICN', _RC)
      call MAPL_StateGetPointer(INTERNAL, NACTL, 'NACTL', _RC)
      call MAPL_StateGetPointer(INTERNAL, NACTI, 'NACTI', _RC)

      ! Import State
      call MAPL_StateGetPointer(IMPORT, PLE, 'PLE', _RC)
      call MAPL_StateGetPointer(IMPORT, ZLE, 'ZLE', _RC)
      call MAPL_StateGetPointer(IMPORT, TS, 'TS', _RC)
      call MAPL_StateGetPointer(IMPORT, T, 'T', _RC)
      call MAPL_StateGetPointer(IMPORT, U, 'U', _RC)
      call MAPL_StateGetPointer(IMPORT, V, 'V', _RC)
      call MAPL_StateGetPointer(IMPORT, W, 'W', _RC)
      call MAPL_StateGetPointer(IMPORT, KPBL, 'KPBL', _RC)
      call MAPL_StateGetPointer(IMPORT, KH, 'KH', _RC)
      call MAPL_StateGetPointer(IMPORT, SH, 'SH', _RC)
      call MAPL_StateGetPointer(IMPORT, EVAP, 'EVAP', _RC)
      call MAPL_StateGetPointer(IMPORT, OMEGA, 'OMEGA', _RC)
      call MAPL_StateGetPointer(IMPORT, TKE, 'TKE', _RC)
      ! Edge fields: MAPL3 pointers are 1:LM+1; MAPL2's were 0:LM, which this code assumes
      p3d => PLE
      PLE(1:IM, 1:JM, 0:LM) => p3d
      p3d => ZLE
      ZLE(1:IM, 1:JM, 0:LM) => p3d
      p3d => KH
      KH(1:IM, 1:JM, 0:LM) => p3d
      p3d => TKE
      TKE(1:IM, 1:JM, 0:LM) => p3d
      call ESMF_StateGet(IMPORT, 'AERO', AERO, _RC)
      ! MTR: Moist provides this service; MAPL puts the bundle in both states. The fields
      ! are the subscribers' internals, so check the two bundles share the data (as DynCore
      ! does for TRADV) and work on the export one.
      call ESMF_StateGet(EXPORT, 'MTR', TR, _RC)
      call ESMF_StateGet(IMPORT, 'MTR', TR_IMPORT, _RC)
      same_mtr_data = MAPL_FieldBundleSameData(TR_IMPORT, TR, _RC)
      _ASSERT(same_mtr_data, "MTR bundles in import and export do not point to the same data")

      ! Update SRF_TYPE for ice_fraction
      call MAPL_StateGetPointer(IMPORT, FRLAND, 'FRLAND', _RC)
      call MAPL_StateGetPointer(IMPORT, FRLANDICE, 'FRLANDICE', _RC)
      call MAPL_StateGetPointer(IMPORT, FRACI, 'FRACI', _RC)
      call MAPL_StateGetPointer(IMPORT, SNOMAS, 'SNOMAS', _RC)
      call MAPL_StateGetPointer(INTERNAL, SRF_TYPE, 'SRF_TYPE', _RC)

      where (FRLANDICE > 0.5)
         SRF_TYPE = SRF_TYPE_LANDICE
         elsewhere (FRACI > 0.5)
         SRF_TYPE = SRF_TYPE_ICE
         elsewhere (SNOMAS > 0.1 .and. SNOMAS /= MAPL2_UNDEF)
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
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'UMST0', _RC)
      if (associated(PTR3D)) PTR3D = U
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'VMST0', _RC)
      if (associated(PTR3D)) PTR3D = V

      ! pre-fill default condensate radii
      call MAPL_StateGetPointer(INTERNAL, PTR3D, 'RL', _RC)
      if (associated(PTR3D)) PTR3D = MAPL_UNDEF
      call MAPL_StateGetPointer(INTERNAL, PTR3D, 'RI', _RC)
      if (associated(PTR3D)) PTR3D = MAPL_UNDEF
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'RR', _RC)
      if (associated(PTR3D)) PTR3D = 50.e-6
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'RS', _RC)
      if (associated(PTR3D)) PTR3D = 50.e-6
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'RG', _RC)
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
      call MAPL_StateGetPointer(INTERNAL, LTS, 'LTS', _RC)
      call MAPL_StateGetPointer(INTERNAL, EIS, 'EIS', _RC)
      KLCL = FIND_KLCL(T, Q, PLmb, IM, JM, LM)
      call MAPL_StateGetPointer(EXPORT, PTR2D, 'ZLCL', _RC)
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
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QX0', _RC)
      if (associated(PTR3D)) PTR3D = Q
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QLLSX0', _RC)
      if (associated(PTR3D)) PTR3D = QLLS
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QLCNX0', _RC)
      if (associated(PTR3D)) PTR3D = QLCN
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'CLCNX0', _RC)
      if (associated(PTR3D)) PTR3D = CLCN
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'CLLSX0', _RC)
      if (associated(PTR3D)) PTR3D = CLLS
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QILSX0', _RC)
      if (associated(PTR3D)) PTR3D = QILS
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QICNX0', _RC)
      if (associated(PTR3D)) PTR3D = QICN
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QCCNX0', _RC)
      if (associated(PTR3D)) PTR3D = QICN + QLCN
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QCLSX0', _RC)
      if (associated(PTR3D)) PTR3D = QILS + QLLS

      ! Trajectory for the moist TLM/ADJ
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'TH_moist', _RC)
      if (associated(PTR3D)) PTR3D = T / PK
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'Q_moist', _RC)
      if (associated(PTR3D)) PTR3D = Q
      call MAPL_StateGetPointer(EXPORT, PTR2D, 'TS_moist', _RC)
      if (associated(PTR2D)) PTR2D = TS
      call MAPL_StateGetPointer(EXPORT, PTR2D, 'KHu_moist', _RC)
      if (associated(PTR2D)) then
         PTR2D = -1
         do L = 0, LM
            where ((PTR2D == -1) .and. (KH(:, :, L) > 2.0))
               PTR2D = L * 1.0
            end where
         end do
      end if
      call MAPL_StateGetPointer(EXPORT, PTR2D, 'KHl_moist', _RC)
      if (associated(PTR2D)) then
         PTR2D = -1
         do L = LM, 0, -1
            where ((PTR2D == -1) .and. (KH(:, :, L) > 2.0))
               PTR2D = L * 1.0
            end where
         end do
      end if

      call MAPL_StateGetPointer(EXPORT, PTR2D, 'KCBL_moist', _RC)
      if (associated(PTR2D)) then
         call MAPL_GridCompGetResource(GC, 'PMIN_CBL', PMIN_CBL, default=50000.0, _RC)
         call MAPL_StateGetPointer(IMPORT, PTR1D, 'PREF', _RC)
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

      ! These may be used by children: re-exported internals
      call MAPL_StateGetPointer(INTERNAL, CNV_FRC, 'CNV_FRC', _RC)
      call MAPL_StateGetPointer(INTERNAL, BYNCY, 'BYNCY', _RC)
      ! Computed locally, zero-filled as MAPL2's ALLOC did (NWFA, NIFA and LCL_AGL are not
      ! always written)
      allocate(NWFA(IM, JM, LM), NIFA(IM, JM, LM), source=0.0)
      allocate(CAPE(IM, JM), INHB(IM, JM), SBCAPE(IM, JM), SBCIN(IM, JM), source=0.0)
      allocate(MLCAPE(IM, JM), MLCIN(IM, JM), MUCAPE(IM, JM), MUCIN(IM, JM), source=0.0)
      allocate(LFC(IM, JM), LNB(IM, JM), LCL_AGL(IM, JM), source=0.0)

      call MAPL_GridCompTimerStart(GC, "-----BUOYANCY2", _RC)
      call BUOYANCY2(IM, JM, LM, T, Q, QST3, DQST3, DZET, ZL0, PLmb, PLEmb(:, :, LM), &
           SBCAPE, MLCAPE, MUCAPE, SBCIN, MLCIN, MUCIN, BYNCY, LFC, LNB, LCL_AGL)
      call MAPL_GridCompTimerStop(GC, "-----BUOYANCY2", _RC)

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

      call MAPL_GridCompTimerStop(GC, "---MOIST_PROLOGUE", _RC)

      ! Extract convective tracers from the TR bundle
      call MAPL_GridCompTimerStart(GC, "---CONV_TRACERS", _RC)
      ! Moist cannot subscribe to its own service, so GF's CNV_TR internal (MAPL2:
      ! FRIENDLYTO MOIST) is added here, to a local bundle: MAPL's MTR bundle is left as is.
      call ESMF_FieldBundleGet(TR, fieldCount=TR_COUNT, _RC)
      allocate(TR_FIELDS(TR_COUNT))
      if (TR_COUNT > 0) then
         call ESMF_FieldBundleGet(TR, fieldList=TR_FIELDS, itemorderflag=ESMF_ITEMORDER_ADDORDER, _RC)
      end if
      if (adjustl(CONVPAR_OPTION) == "GF") then
         call ESMF_StateGet(INTERNAL, 'CNV_TR', CNV_TR_FIELD, _RC)
         TR_FIELDS = [TR_FIELDS, CNV_TR_FIELD]
      end if
      TR_CNV = ESMF_FieldBundleCreate(fieldList=TR_FIELDS, name='MTR_CNV', _RC)
      call CNV_Tracers_Init(TR_CNV, RC)
      call ESMF_FieldBundleDestroy(TR_CNV, _RC) ! the fields are not destroyed
      call MAPL_GridCompTimerStop(GC, "---CONV_TRACERS", _RC)

      ! Get aerosol activation properties
      call MAPL_GridCompTimerStart(GC, "---AERO_ACTIVATE", _RC)

      if (MOIST_USE_NCLOUD_CLIM) then !Setup ND/NI climatology from GiOcean
         call MAPL_StateGetPointer(IMPORT, PTR3D, 'NCPL_CLIM', _RC)
         NACTL = PTR3D
         call MAPL_StateGetPointer(IMPORT, PTR3D, 'NCPI_CLIM', _RC)
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
            call Aer_Activation(GC, IM, JM, LM, Q, T, PLmb * 100.0, PLE, TKE, TMP3D, FRLAND, &
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

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'NCCN_LIQ', _RC)
      if (associated(PTR3D)) PTR3D = NACTL
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'NCCN_ICE', _RC)
      if (associated(PTR3D)) PTR3D = NACTI

      call MAPL_GridCompTimerStop(GC, "---AERO_ACTIVATE", _RC)

      if (SH_MD_DP) then
         if (adjustl(SHALLOW_OPTION) == "UW") call UW_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
         if (adjustl(SHALLOW_OPTION) == "UW") then
            if (DEBUG_MST) then
               call report_maxmin(logger, VMG, 'MST: Q_AF_UW  ', Q, _RC)
               call report_maxmin(logger, VMG, 'MST: T_AF_UW  ', T, _RC)
               call report_maxmin(logger, VMG, 'MST: U_AF_UW  ', U, _RC)
               call report_maxmin(logger, VMG, 'MST: V_AF_UW  ', V, _RC)
            end if
         end if
         if (adjustl(CONVPAR_OPTION) == "RAS") call RAS_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
         if (adjustl(CONVPAR_OPTION) == "GF") call GF_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
         if (adjustl(CONVPAR_OPTION) == "GF") then
            if (DEBUG_MST) then
               call report_maxmin(logger, VMG, 'MST: Q_AF_GF  ', Q, _RC)
               call report_maxmin(logger, VMG, 'MST: T_AF_GF  ', T, _RC)
               call report_maxmin(logger, VMG, 'MST: U_AF_GF  ', U, _RC)
               call report_maxmin(logger, VMG, 'MST: V_AF_GF  ', V, _RC)
            end if
         end if
      else
         if (adjustl(CONVPAR_OPTION) == "RAS") call RAS_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
         if (adjustl(CONVPAR_OPTION) == "GF") call GF_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
         if (adjustl(CONVPAR_OPTION) == "GF") then
            if (DEBUG_MST) then
               call report_maxmin(logger, VMG, 'MST: Q_AF_GF  ', Q, _RC)
               call report_maxmin(logger, VMG, 'MST: T_AF_GF  ', T, _RC)
               call report_maxmin(logger, VMG, 'MST: U_AF_GF  ', U, _RC)
               call report_maxmin(logger, VMG, 'MST: V_AF_GF  ', V, _RC)
            end if
         end if
         if (adjustl(SHALLOW_OPTION) == "UW") call UW_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
         if (adjustl(SHALLOW_OPTION) == "UW") then
            if (DEBUG_MST) then
               call report_maxmin(logger, VMG, 'MST: Q_AF_UW  ', Q, _RC)
               call report_maxmin(logger, VMG, 'MST: T_AF_UW  ', T, _RC)
               call report_maxmin(logger, VMG, 'MST: U_AF_UW  ', U, _RC)
               call report_maxmin(logger, VMG, 'MST: V_AF_UW  ', V, _RC)
            end if
         end if
      end if

      call MAPL_GridCompTimerStart(GC, "---MOIST_EPILOGUE", _RC)

      ! Mass fluxes
      allocate(TMP3Dp1(IM, JM, LM + 1))
      ! Edge updraft mass flux over deep and shalow convection
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'CNV_MFC', _RC)
      call MAPL_StateGetPointer(INTERNAL, PTRDC, 'UMF_DC', _RC)
      call MAPL_StateGetPointer(INTERNAL, PTRSC, 'UMF_SC', _RC)
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

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'CNV_MFD', _RC)
      call MAPL_StateGetPointer(INTERNAL, PTRDC, 'MFD_DC', _RC)
      call MAPL_StateGetPointer(INTERNAL, PTRSC, 'MFD_SC', _RC)
      TMP3D = 0.0
      if (associated(PTRDC)) TMP3D = TMP3D + PTRDC
      if (associated(PTRSC)) TMP3D = TMP3D + PTRSC
      if (associated(PTR3D)) PTR3D = TMP3D

      call MAPL_GridCompTimerStop(GC, "---MOIST_EPILOGUE", _RC)

      if (adjustl(CLDMICR_OPTION) == "BACM_1M") call BACM_1M_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
      if (adjustl(CLDMICR_OPTION) == "GFDL_1M") call GFDL_1M_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
      if (adjustl(CLDMICR_OPTION) == "THOM_1M") call THOM_1M_Run(GC, IMPORT, EXPORT, CLOCK, _RC)
      if (adjustl(CLDMICR_OPTION) == "MGB2_2M") call MGB2_2M_Run(GC, IMPORT, EXPORT, CLOCK, _RC)

      call MAPL_GridCompTimerStart(GC, "---MOIST_EPILOGUE", _RC)

      if (DEBUG_MST) then
         call report_maxmin(logger, VMG, 'MST: Q_AF_MP  ', Q, _RC)
         call report_maxmin(logger, VMG, 'MST: T_AF_MP  ', T, _RC)
         call report_maxmin(logger, VMG, 'MST: U_AF_MP  ', U, _RC)
         call report_maxmin(logger, VMG, 'MST: V_AF_MP  ', V, _RC)
      end if

      ! Exports
      ! Cloud fraction exports
      allocate(CFICE(IM, JM, LM), source=0.0)
      where (QILS + QICN > 1.0e-12)
         CFICE = (CLLS + CLCN) * (QILS + QICN) / (QLLS + QLCN + QILS + QICN)
      end where
      CFICE = max(min(CFICE, 1.0), 0.0)
      call MAPL_StateGetPointer(EXPORT, CFLIQ, 'CFLIQ', _RC)
      if (associated(CFLIQ)) then
         CFLIQ = 0.0
         where (QLLS + QLCN > 1.0e-12)
            CFLIQ = (CLLS + CLCN) * (QLLS + QLCN) / (QLLS + QLCN + QILS + QICN)
         end where
         CFLIQ = max(min(CFLIQ, 1.0), 0.0)
      end if

      ! some diagnostics to export
      QST3 = GEOS_QsatICE(T, PLmb * 100.0, DQ=DQST3)
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'RHICE', _RC)
      if (associated(PTR3D)) then
         PTR3D = Q / QST3
         where (T > MAPL_TICE)
            PTR3D = 0.0
         end where
      end if
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'SAT_RAT', _RC)
      if (associated(PTR3D)) then
         where (CFICE < 0.99 .and. QST3 > 1.0e-20)
            TMP3D = max((Q - QST3 * CFICE), 0.0) / (1.0 - CFICE)
            PTR3D = min(TMP3D / QST3, 2.0)
            elsewhere
            PTR3D = 1.0
         end where
      end if

      QST3 = GEOS_QsatLQU(T, PLmb * 100.0, DQ=DQST3)
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'RHLIQ', _RC)
      if (associated(PTR3D)) PTR3D = Q / QST3

      ! Clean up supersaturation only with respect to liquid water
      allocate(DTDT_ER(IM, JM, LM), DQVDT_ER(IM, JM, LM), source=0.0)
      DTDT_ER = T
      DQVDT_ER = Q

      ! Rain-out of Relative Humidity where RH > 110%
      call MAPL_StateGetPointer(INTERNAL, LS_PRCP, 'LS_PRCP', _RC)
      allocate(ER_PRCP(IM, JM), source=0.0)
      where (Q > 1.1 * QST3)
         TMP3D = (Q - 1.1 * QST3) / (1.0 + 1.1 * DQST3 * MAPL_ALHL / MAPL_CP)
         elsewhere
         TMP3D = 0.0
      end where
      ! protect against producing negatice QV
      where (TMP3D > Q)
         TMP3D = Q
      end where
      ER_PRCP = SUM(TMP3D * MASS, 3) / DT_MOIST
      LS_PRCP = LS_PRCP + ER_PRCP
      Q = Q - TMP3D
      T = T + (MAPL_ALHL / MAPL_CP) * TMP3D

      DTDT_ER = (T - DTDT_ER) / DT_MOIST
      DQVDT_ER = (Q - DQVDT_ER) / DT_MOIST

      ! Export Total Moist Tendencies

      call MAPL_StateGetPointer(EXPORT, DUDT, 'DUDT', _RC)
      if (associated(DUDT)) then
         DUDT = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DUDT_DC', _RC)
         if (associated(PTR3D)) DUDT = DUDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DUDT_SC', _RC)
         if (associated(PTR3D)) DUDT = DUDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DUDT_macro', _RC)
         if (associated(PTR3D)) DUDT = DUDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DUDT_micro', _RC)
         if (associated(PTR3D)) DUDT = DUDT + PTR3D
      end if

      call MAPL_StateGetPointer(EXPORT, DVDT, 'DVDT', _RC)
      if (associated(DVDT)) then
         DVDT = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DVDT_DC', _RC)
         if (associated(PTR3D)) DVDT = DVDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DVDT_SC', _RC)
         if (associated(PTR3D)) DVDT = DVDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DVDT_macro', _RC)
         if (associated(PTR3D)) DVDT = DVDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DVDT_micro', _RC)
         if (associated(PTR3D)) DVDT = DVDT + PTR3D
      end if

      call MAPL_StateGetPointer(EXPORT, DTDT, 'DTDT', _RC)
      if (associated(DTDT)) then
         DTDT = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DTDT_DC', _RC)
         if (associated(PTR3D)) DTDT = DTDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DTDT_SC', _RC)
         if (associated(PTR3D)) DTDT = DTDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DTDT_macro', _RC)
         if (associated(PTR3D)) DTDT = DTDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DTDT_micro', _RC)
         if (associated(PTR3D)) DTDT = DTDT + PTR3D
         DTDT = DTDT + DTDT_ER
         DTDT = DTDT * (PLE(:, :, 1:LM) - PLE(:, :, 0:LM - 1)) ! Pressure weighted tendency
      end if

      call MAPL_StateGetPointer(EXPORT, DQDT, 'DQDT', _RC)
      if (associated(DQDT)) then
         DQDT = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQVDT_DC', _RC)
         if (associated(PTR3D)) DQDT = DQDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQVDT_SC', _RC)
         if (associated(PTR3D)) DQDT = DQDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQVDT_macro', _RC)
         if (associated(PTR3D)) DQDT = DQDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQVDT_micro', _RC)
         if (associated(PTR3D)) DQDT = DQDT + PTR3D
         DQDT = DQDT + DQVDT_ER
      end if

      call MAPL_StateGetPointer(EXPORT, DQLDT, 'DQLDT', _RC)
      if (associated(DQLDT)) then
         DQLDT = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQLDT_DC', _RC)
         if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQLDT_SC', _RC)
         if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQLDT_macro', _RC)
         if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQLDT_micro', _RC)
         if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
      end if

      call MAPL_StateGetPointer(EXPORT, DQIDT, 'DQIDT', _RC)
      if (associated(DQIDT)) then
         DQIDT = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQIDT_DC', _RC)
         if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQIDT_SC', _RC)
         if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQIDT_macro', _RC)
         if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQIDT_micro', _RC)
         if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
      end if

      call MAPL_StateGetPointer(EXPORT, DQRDT, 'DQRDT', _RC)
      if (associated(DQRDT)) then
         DQRDT = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQRDT_SC', _RC)
         if (associated(PTR3D)) DQRDT = DQRDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQRDT_macro', _RC)
         if (associated(PTR3D)) DQRDT = DQRDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQRDT_micro', _RC)
         if (associated(PTR3D)) DQRDT = DQRDT + PTR3D
      end if

      call MAPL_StateGetPointer(EXPORT, DQSDT, 'DQSDT', _RC)
      if (associated(DQSDT)) then
         DQSDT = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQSDT_SC', _RC)
         if (associated(PTR3D)) DQSDT = DQSDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQSDT_macro', _RC)
         if (associated(PTR3D)) DQSDT = DQSDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQSDT_micro', _RC)
         if (associated(PTR3D)) DQSDT = DQSDT + PTR3D
      end if

      call MAPL_StateGetPointer(EXPORT, DQGDT, 'DQGDT', _RC)
      if (associated(DQGDT)) then
         DQGDT = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQGDT_macro', _RC)
         if (associated(PTR3D)) DQGDT = DQGDT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQGDT_micro', _RC)
         if (associated(PTR3D)) DQGDT = DQGDT + PTR3D
      end if

      call MAPL_StateGetPointer(EXPORT, DQADT, 'DQADT', _RC)
      if (associated(DQADT)) then
         DQADT = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQADT_DC', _RC)
         if (associated(PTR3D)) DQADT = DQADT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQADT_SC', _RC)
         if (associated(PTR3D)) DQADT = DQADT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQADT_macro', _RC)
         if (associated(PTR3D)) DQADT = DQADT + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'DQADT_micro', _RC)
         if (associated(PTR3D)) DQADT = DQADT + PTR3D
      end if

      call MAPL_StateGetPointer(EXPORT, DPDTMST, 'DPDTMST', _RC)
      if (associated(DPDTMST)) then
         DPDTMST = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'PFI_CN', _RC)
         p3d => PTR3D
         PTR3D(1:IM, 1:JM, 0:LM) => p3d ! edge: MAPL2's 0:LM bounds
         if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
         call MAPL_StateGetPointer(EXPORT, PTR3D, 'PFI_SC', _RC)
         if (associated(PTR3D)) then
            p3d => PTR3D
            PTR3D(1:IM, 1:JM, 0:LM) => p3d ! edge: MAPL2's 0:LM bounds
            DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
         end if
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'PFI_AN', _RC)
         p3d => PTR3D
         PTR3D(1:IM, 1:JM, 0:LM) => p3d ! edge: MAPL2's 0:LM bounds
         if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'PFI_LS', _RC)
         p3d => PTR3D
         PTR3D(1:IM, 1:JM, 0:LM) => p3d ! edge: MAPL2's 0:LM bounds
         if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'PFL_CN', _RC)
         p3d => PTR3D
         PTR3D(1:IM, 1:JM, 0:LM) => p3d ! edge: MAPL2's 0:LM bounds
         if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
         call MAPL_StateGetPointer(EXPORT, PTR3D, 'PFL_SC', _RC)
         if (associated(PTR3D)) then
            p3d => PTR3D
            PTR3D(1:IM, 1:JM, 0:LM) => p3d ! edge: MAPL2's 0:LM bounds
            DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
         end if
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'PFL_AN', _RC)
         p3d => PTR3D
         PTR3D(1:IM, 1:JM, 0:LM) => p3d ! edge: MAPL2's 0:LM bounds
         if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'PFL_LS', _RC)
         p3d => PTR3D
         PTR3D(1:IM, 1:JM, 0:LM) => p3d ! edge: MAPL2's 0:LM bounds
         if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:, :, 0:LM - 1) - PTR3D(:, :, 1:LM)
         DPDTMST = MAPL_GRAV * DPDTMST
      end if

      call MAPL_StateGetPointer(EXPORT, PFL_LSAN, 'PFL_LSAN', _RC)
      if (associated(PFL_LSAN)) then
         PFL_LSAN = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'PFL_AN', _RC)
         if (associated(PTR3D)) PFL_LSAN = PFL_LSAN + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'PFL_LS', _RC)
         if (associated(PTR3D)) PFL_LSAN = PFL_LSAN + PTR3D
      end if

      call MAPL_StateGetPointer(EXPORT, PFI_LSAN, 'PFI_LSAN', _RC)
      if (associated(PFI_LSAN)) then
         PFI_LSAN = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'PFI_AN', _RC)
         if (associated(PTR3D)) PFI_LSAN = PFI_LSAN + PTR3D
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'PFI_LS', _RC)
         if (associated(PTR3D)) PFI_LSAN = PFI_LSAN + PTR3D
      end if

      ! Combine Precip Exports

      ! liquid convective precip
      ! PCU and PLS are local arrays: DIAGNOSE_PRECIP_TYPE takes them unguarded (in MAPL2
      ! Surface always connected them)
      allocate(PCU(IM, JM), source=0.0)
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'CN_PRCP', _RC)
      if (associated(PTR2D)) PCU = PCU + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'SC_PRCP', _RC)
      if (associated(PTR2D)) PCU = PCU + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'CNPCPRATE', _RC)
      if (associated(PTR2D)) PCU = PCU + PTR2D
      PCU = max(PCU, 0.0)

      ! liquid large-scale precip
      allocate(PLS(IM, JM), source=0.0)
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'LS_PRCP', _RC)
      if (associated(PTR2D)) PLS = PLS + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'AN_PRCP', _RC)
      if (associated(PTR2D)) PLS = PLS + PTR2D
      PLS = max(PLS, 0.0)

      ! all liquid precip
      allocate(RAIN(IM, JM), source=0.0)
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'LS_PRCP', _RC)
      if (associated(PTR2D)) RAIN = RAIN + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'AN_PRCP', _RC)
      if (associated(PTR2D)) RAIN = RAIN + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'CN_PRCP', _RC)
      if (associated(PTR2D)) RAIN = RAIN + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'SC_PRCP', _RC)
      if (associated(PTR2D)) RAIN = RAIN + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'CNPCPRATE', _RC)
      if (associated(PTR2D)) RAIN = RAIN + PTR2D
      RAIN = max(RAIN, 0.0)

      ! all frozen precip (snow at this point)
      allocate(SNOW(IM, JM), source=0.0)
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'LS_SNR', _RC)
      if (associated(PTR2D)) SNOW = SNOW + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'AN_SNR', _RC)
      if (associated(PTR2D)) SNOW = SNOW + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'CN_SNR', _RC)
      if (associated(PTR2D)) SNOW = SNOW + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'SC_SNR', _RC)
      if (associated(PTR2D)) SNOW = SNOW + PTR2D
      SNOW = max(SNOW, 0.0)

      ! all deep convective precip (rain+snow)
      call MAPL_StateGetPointer(INTERNAL, CN_PRCP, 'CN_PRCP', _RC)
      if (associated(CN_PRCP)) then
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'CNPCPRATE', _RC)
         if (associated(PTR2D)) CN_PRCP = CN_PRCP + PTR2D
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'CN_SNR', _RC)
         if (associated(PTR2D)) CN_PRCP = CN_PRCP + PTR2D
         CN_PRCP = max(CN_PRCP, 0.0)
      end if

      ! all large-scale precip (rain+snow)
      call MAPL_StateGetPointer(INTERNAL, LS_PRCP, 'LS_PRCP', _RC)
      if (associated(LS_PRCP)) then
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'LS_SNR', _RC)
         if (associated(PTR2D)) LS_PRCP = LS_PRCP + PTR2D
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'ICE', _RC)
         if (associated(PTR2D)) LS_PRCP = LS_PRCP + PTR2D
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'FRZR', _RC)
         if (associated(PTR2D)) LS_PRCP = LS_PRCP + PTR2D
         LS_PRCP = max(LS_PRCP, 0.0)
      end if

      ! all anvil precip (rain+snow)
      call MAPL_StateGetPointer(INTERNAL, AN_PRCP, 'AN_PRCP', _RC)
      if (associated(AN_PRCP)) then
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'AN_SNR', _RC)
         if (associated(PTR2D)) AN_PRCP = AN_PRCP + PTR2D
      end if

      ! all shallow precip (rain+snow)
      call MAPL_StateGetPointer(INTERNAL, SC_PRCP, 'SC_PRCP', _RC)
      if (associated(SC_PRCP)) then
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'SC_SNR', _RC)
         if (associated(PTR2D)) SC_PRCP = SC_PRCP + PTR2D
         SC_PRCP = max(SC_PRCP, 0.0)
      end if

      ! Total - all precip (rain+snow)
      allocate(TPREC(IM, JM), source=0.0)
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'LS_PRCP', _RC)
      if (associated(PTR2D)) TPREC = TPREC + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'AN_PRCP', _RC)
      if (associated(PTR2D)) TPREC = TPREC + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'CN_PRCP', _RC)
      if (associated(PTR2D)) TPREC = TPREC + PTR2D
      call MAPL_StateGetPointer(INTERNAL, PTR2D, 'SC_PRCP', _RC)
      if (associated(PTR2D)) TPREC = TPREC + PTR2D
      TPREC = max(TPREC, 0.0)

      ! diagnosed stratiform precip (rain+snow)
      call MAPL_StateGetPointer(EXPORT, PREC_STRAT, 'PREC_STRAT', _RC)
      if (associated(PREC_STRAT)) PREC_STRAT = max((1.0 - CNV_FRC) * TPREC, 0.0)

      ! diagnosed convective precip (rain+snow)
      call MAPL_StateGetPointer(EXPORT, PREC_CONV, 'PREC_CONV', _RC)
      if (associated(PREC_CONV)) PREC_CONV = max(CNV_FRC * TPREC, 0.0)

      ! Diagnostic precip types:
      call MAPL_StateGetPointer(INTERNAL, ICE, 'ICE', _RC)
      call MAPL_StateGetPointer(INTERNAL, FRZR, 'FRZR', _RC)
      if (LUPDATE_PRECIP_TYPE .or. LDIAGNOSE_PRECIP_TYPE) then
         allocate(PTYPE(IM, JM), source=0.0)
         call DIAGNOSE_PRECIP_TYPE(IM, JM, LM, TPREC, PLS, PCU, RAIN, SNOW, ICE, FRZR, &
              PTYPE, PLE, T / PK, PK, PKE, ZL0, LUPDATE_PRECIP_TYPE)
         call MAPL_StateGetPointer(EXPORT, PTYPE_ptr, 'PTYPE', _RC)
         if (associated(PTYPE_ptr)) PTYPE_ptr = PTYPE
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
      call MAPL_StateGetPointer(EXPORT, PTR2D, 'KUCHERA_RATIO', _RC)
      if (associated(PTR2D)) PTR2D = TMP2D
      ! Accumulated precip totals (mm), apply KUCHERA_RATIO for SNOW
      call MAPL_StateGetPointer(EXPORT, PTR2D, 'SNOWTOTAL', _RC)
      if (associated(PTR2D)) PTR2D = TMP2D * DT_MOIST * (SNOW + ICE)
      call MAPL_StateGetPointer(EXPORT, PTR2D, 'PRECTOTAL', _RC)
      if (associated(PTR2D)) PTR2D = DT_MOIST * TPREC

      call MAPL_StateGetPointer(INTERNAL, PTR3D, 'QLTOT', _RC)
      if (associated(PTR3D)) PTR3D = QLLS + QLCN

      call MAPL_StateGetPointer(INTERNAL, PTR3D, 'QITOT', _RC)
      if (associated(PTR3D)) PTR3D = QILS + QICN

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QCTOT', _RC)
      if (associated(PTR3D)) then
         PTR3D = QLLS + QLCN + QILS + QICN
         call MAPL_StateGetPointer(EXPORT, QRTOT, 'QRTOT', _RC)
         if (associated(QRTOT)) PTR3D = PTR3D + QRTOT
         call MAPL_StateGetPointer(EXPORT, QSTOT, 'QSTOT', _RC)
         if (associated(QSTOT)) PTR3D = PTR3D + QSTOT
         call MAPL_StateGetPointer(EXPORT, QGTOT, 'QGTOT', _RC)
         if (associated(QGTOT)) PTR3D = PTR3D + QGTOT
      end if

      ! Cloud condensate exports
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QLLSX1', _RC)
      if (associated(PTR3D)) PTR3D = QLLS

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QILSX1', _RC)
      if (associated(PTR3D)) PTR3D = QILS

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QLCNX1', _RC)
      if (associated(PTR3D)) PTR3D = QLCN

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QICNX1', _RC)
      if (associated(PTR3D)) PTR3D = QICN

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'CLLSX1', _RC)
      if (associated(PTR3D)) PTR3D = CLLS

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'CLCNX1', _RC)
      if (associated(PTR3D)) PTR3D = CLCN

      ! Fill wind, temperature & RH exports needed for SYNCTQ

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'UAFMOIST', _RC)
      if (associated(PTR3D)) PTR3D = U

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'VAFMOIST', _RC)
      if (associated(PTR3D)) PTR3D = V

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'TAFMOIST', _RC)
      if (associated(PTR3D)) PTR3D = T

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'QAFMOIST', _RC)
      if (associated(PTR3D)) PTR3D = Q

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'THAFMOIST', _RC)
      if (associated(PTR3D)) PTR3D = T / PK

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'SAFMOIST', _RC)
      if (associated(PTR3D)) then
         do L = 1, LM
            PTR3D(:, :, L) = MAPL_CP * T(:, :, L) + MAPL_GRAV * (ZL0(:, :, L) + ZLE(:, :, LM))
         end do
      end if

      call MAPL_StateGetPointer(EXPORT, PTR3D, 'RH2', _RC)
      if (associated(PTR3D)) PTR3D = max(min(Q / GEOS_QSAT(T, PLmb), 1.02), 0.0)

      ! Cloud Water Path
      call MAPL_StateGetPointer(EXPORT, PTR2D, 'CWP', _RC)
      if (associated(PTR2D)) PTR2D = SUM((QLCN + QLLS + QICN + QILS) * MASS, 3)
      ! Cloud Liquid Water Path
      call MAPL_StateGetPointer(EXPORT, PTR2D, 'CLWP', _RC)
      if (associated(PTR2D)) PTR2D = SUM((QLCN + QLLS) * MASS, 3)
      ! Cloud Ice Water Path
      call MAPL_StateGetPointer(EXPORT, PTR2D, 'CIWP', _RC)
      if (associated(PTR2D)) PTR2D = SUM((QICN + QILS) * MASS, 3)
      ! Total LWP and IWP need to be defined inside the microphysics depends on number was water species

      ! Total Precipitable Water
      call MAPL_StateGetPointer(EXPORT, PTR2D, 'TPW', _RC)
      if (associated(PTR2D)) PTR2D = SUM((Q) * MASS, 3)

      ! Lightning Exports
      call MAPL_StateGetPointer(EXPORT, PTR2D, 'LFR_GCC', _RC)
      if (associated(PTR2D)) PTR2D = 0.0

      ! Copy the local arrays to their exports, if connected (PTYPE is copied above)
      call MAPL_StateGetPointer(EXPORT, NWFA_ptr, 'NWFA', _RC)
      call MAPL_StateGetPointer(EXPORT, NIFA_ptr, 'NIFA', _RC)
      call MAPL_StateGetPointer(EXPORT, CAPE_ptr, 'CAPE', _RC)
      call MAPL_StateGetPointer(EXPORT, INHB_ptr, 'INHB', _RC)
      call MAPL_StateGetPointer(EXPORT, SBCAPE_ptr, 'SBCAPE', _RC)
      call MAPL_StateGetPointer(EXPORT, SBCIN_ptr, 'SBCIN', _RC)
      call MAPL_StateGetPointer(EXPORT, MLCAPE_ptr, 'MLCAPE', _RC)
      call MAPL_StateGetPointer(EXPORT, MLCIN_ptr, 'MLCIN', _RC)
      call MAPL_StateGetPointer(EXPORT, MUCAPE_ptr, 'MUCAPE', _RC)
      call MAPL_StateGetPointer(EXPORT, MUCIN_ptr, 'MUCIN', _RC)
      call MAPL_StateGetPointer(EXPORT, LFC_ptr, 'ZLFC', _RC)
      call MAPL_StateGetPointer(EXPORT, LNB_ptr, 'ZLNB', _RC)
      call MAPL_StateGetPointer(EXPORT, LCL_AGL_ptr, 'LCL_AGL', _RC)
      call MAPL_StateGetPointer(EXPORT, CFICE_ptr, 'CFICE', _RC)
      call MAPL_StateGetPointer(EXPORT, DTDT_ER_ptr, 'DTDT_ER', _RC)
      call MAPL_StateGetPointer(EXPORT, DQVDT_ER_ptr, 'DQVDT_ER', _RC)
      call MAPL_StateGetPointer(EXPORT, ER_PRCP_ptr, 'ER_PRCP', _RC)
      call MAPL_StateGetPointer(EXPORT, RAIN_ptr, 'RAIN', _RC)
      call MAPL_StateGetPointer(EXPORT, SNOW_ptr, 'SNO', _RC)
      call MAPL_StateGetPointer(EXPORT, TPREC_ptr, 'TPREC', _RC)
      call MAPL_StateGetPointer(EXPORT, PCU_ptr, 'PCU', _RC)
      call MAPL_StateGetPointer(EXPORT, PLS_ptr, 'PLS', _RC)
      if (associated(NWFA_ptr)) NWFA_ptr = NWFA
      if (associated(NIFA_ptr)) NIFA_ptr = NIFA
      if (associated(CAPE_ptr)) CAPE_ptr = CAPE
      if (associated(INHB_ptr)) INHB_ptr = INHB
      if (associated(SBCAPE_ptr)) SBCAPE_ptr = SBCAPE
      if (associated(SBCIN_ptr)) SBCIN_ptr = SBCIN
      if (associated(MLCAPE_ptr)) MLCAPE_ptr = MLCAPE
      if (associated(MLCIN_ptr)) MLCIN_ptr = MLCIN
      if (associated(MUCAPE_ptr)) MUCAPE_ptr = MUCAPE
      if (associated(MUCIN_ptr)) MUCIN_ptr = MUCIN
      if (associated(LFC_ptr)) LFC_ptr = LFC
      if (associated(LNB_ptr)) LNB_ptr = LNB
      if (associated(LCL_AGL_ptr)) LCL_AGL_ptr = LCL_AGL
      if (associated(CFICE_ptr)) CFICE_ptr = CFICE
      if (associated(DTDT_ER_ptr)) DTDT_ER_ptr = DTDT_ER
      if (associated(DQVDT_ER_ptr)) DQVDT_ER_ptr = DQVDT_ER
      if (associated(ER_PRCP_ptr)) ER_PRCP_ptr = ER_PRCP
      if (associated(RAIN_ptr)) RAIN_ptr = RAIN
      if (associated(SNOW_ptr)) SNOW_ptr = SNOW
      if (associated(TPREC_ptr)) TPREC_ptr = TPREC
      if (associated(PCU_ptr)) PCU_ptr = PCU
      if (associated(PLS_ptr)) PLS_ptr = PLS
      ! Allocated as pointers (see their declaration), so not freed on return
      deallocate(MLCAPE, MLCIN, MUCAPE, MUCIN)

      call MAPL_GridCompTimerStop(GC, "---MOIST_EPILOGUE", _RC)

      call MAPL_GridCompTimerStop(GC, "TOTAL", _RC)

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

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(CLOCK)
   end subroutine Finalize

   ! DEBUG_MST print of the global max and min of a field (MAPL2's MAPL_MaxMin printed them;
   ! MAPL3's returns [max, min])
   subroutine report_maxmin(logger, vm, label, x, rc)
      class(logger_t), pointer, intent(in) :: logger
      type(ESMF_VM), intent(in) :: vm
      character(len=*), intent(in) :: label
      real, intent(in) :: x(:, :, :)
      integer, optional, intent(out) :: rc

      integer :: status, comm
      real :: mm(2)

      call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
      mm = MAPL_MaxMin(x, comm, _RC)
      call logger%info(label // ' max, min: %es12.5 %es12.5', mm(1), mm(2))

      _RETURN(_SUCCESS)
   end subroutine report_maxmin

end module GEOS_MoistGridCompMod
