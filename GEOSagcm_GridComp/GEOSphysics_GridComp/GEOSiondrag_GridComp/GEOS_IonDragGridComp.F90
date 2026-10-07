#include "MAPL_Generic.h"

module GEOS_IonDragGridCompMod

   !BOP
   ! !MODULE: GEOS_IonDrag -- Ion drag forcing for GEOS-MLT
   !
   ! !DESCRIPTION:
   !
   ! The component computes horizontal neutral-wind tendencies from
   ! magnetized ion-neutral coupling using the WACCM-X/TIE-GCM conductivity
   ! formulation. IRI supplies major-ion densities and plasma temperatures,
   ! IGRF supplies the local magnetic field, MSIS supplies absolute neutral
   ! O/O2/N2 densities and neutral mass density, and the ML ion-velocity model
   ! supplies geographic eastward/northward ion drift velocities.
   !
   ! Momentum tendencies are returned to GEOS_PhysicsGridComp. The associated
   ! ion-neutral friction/Joule heating is exported only as a diagnostic
   ! temperature tendency because GEOS-MLT already applies MLRADJH.
   !

   use ESMF
   use MAPL
   use MAPL_PythonBridge, only: MAPL_pybridge_gcinit, &
                                MAPL_pybridge_gcrun_with_internal

   use iri_input_module, only: get_iri_state
   use igrf_input_module, only: get_igrf_field
   use ion_drag_module, only: compute_drag_fields, N_MAJOR_ION_SPECIES, &
                              ION_OP, ION_O2P, ION_NOP
   use msis_wrapper, only: msis_wrapper_init, msis_prepare_time, msis_point, &
                           msis_get_current_f107
   use calc_gas_specific_heat_mlt_mod, only: mlt_mixture_thermo_from_number_density
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   implicit none
   private

   public SetServices

   ! Ion drag is a fixed part of the GEOS-MLT configuration. The lower
   ! boundary is 1 Pa = 0.01 hPa; ion drag is negligible below this region.
   real, parameter :: IONDRAG_BOTTOM_PRESSURE_PA = 1.0
   real, parameter :: FALLBACK_ALT_KM = 220.0

   type :: GEOS_IonDragGridComp
      logical :: MLION_PYBRIDGE_INITIALIZED = .false.

      integer :: MLION_LAST_YEAR = -1
      integer :: MLION_LAST_DOY = -1
      integer :: MLION_LAST_HOUR = -1

      integer :: PLASMA_LAST_YEAR = -1
      integer :: PLASMA_LAST_DOY = -1
      integer :: PLASMA_LAST_HOUR = -1


      real, allocatable :: NE_CACHE(:,:,:)
      real, allocatable :: ION_DENSITY_CACHE(:,:,:,:)
      real, allocatable :: TI_CACHE(:,:,:)
      real, allocatable :: TE_CACHE(:,:,:)
      real, allocatable :: BNORTH_CACHE(:,:,:)
      real, allocatable :: BEAST_CACHE(:,:,:)
      real, allocatable :: BDOWN_CACHE(:,:,:)
      real, allocatable :: BMAG_CACHE(:,:,:)
   end type GEOS_IonDragGridComp

   type wrap_
      type (GEOS_IonDragGridComp), pointer :: PTR
   end type wrap_

contains

   subroutine SetServices ( GC, RC )

      type(ESMF_GridComp), intent(INOUT) :: GC
      integer, optional                  :: RC

      character(len=ESMF_MAXSTR)          :: IAm
      integer                             :: STATUS
      character(len=ESMF_MAXSTR)          :: COMP_NAME
      type (MAPL_MetaComp), pointer       :: MAPL
      type (wrap_)                        :: wrap
      type (GEOS_IonDragGridComp), pointer :: self

      Iam = 'SetServices'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC )
      Iam = trim(COMP_NAME) // Iam

      allocate (self, _STAT)
      wrap%ptr => self

      call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN, Run, _RC)

      call MAPL_GetObjectFromGC ( GC, MAPL, _RC )

      call register_state_specs(GC, _RC)

      call MAPL_TimerAdd(GC, name="DRIVER", _RC)
      call MAPL_TimerAdd(GC, name="-MLION", _RC)
      call MAPL_TimerAdd(GC, name="-IRI", _RC)
      call MAPL_TimerAdd(GC, name="-IGRF", _RC)
      call MAPL_TimerAdd(GC, name="-MSIS", _RC)
      call MAPL_TimerAdd(GC, name="-DRAG", _RC)

      call ESMF_UserCompSetInternalState ( GC, 'GEOS_IonDragGridComp', wrap, _RC )

      call MAPL_GenericSetServices ( gc, _RC)

      RETURN_(ESMF_SUCCESS)

   end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine register_state_specs(GC, RC)
      ! Register MAPL imports, exports, and internal fields explicitly.
      ! This replaces IonDrag_StateSpecs.rc and the mapl_acg-generated files.

      type(ESMF_GridComp), intent(inout) :: GC
      integer, optional, intent(out) :: RC

      character(len=ESMF_MAXSTR) :: IAm
      integer :: STATUS

      IAm = 'register_state_specs'

      ! ==================
      ! Internal State
      ! ==================
      call MAPL_AddInternalSpec(GC, SHORT_NAME='MLION_LATS', &
           LONG_NAME='ml_ion_velocity_helper_latitude', UNITS='radians', &
           DIMS=MAPL_DimsHorzOnly, VLOCATION=MAPL_VLocationNone, _RC)
      call MAPL_AddInternalSpec(GC, SHORT_NAME='MLION_LONS', &
           LONG_NAME='ml_ion_velocity_helper_longitude', UNITS='radians', &
           DIMS=MAPL_DimsHorzOnly, VLOCATION=MAPL_VLocationNone, _RC)
      call MAPL_AddInternalSpec(GC, SHORT_NAME='MLION_YY', &
           LONG_NAME='ml_ion_velocity_helper_year', UNITS='1', &
           DIMS=MAPL_DimsHorzOnly, VLOCATION=MAPL_VLocationNone, _RC)
      call MAPL_AddInternalSpec(GC, SHORT_NAME='MLION_DOY', &
           LONG_NAME='ml_ion_velocity_helper_day_of_year', UNITS='1', &
           DIMS=MAPL_DimsHorzOnly, VLOCATION=MAPL_VLocationNone, _RC)
      call MAPL_AddInternalSpec(GC, SHORT_NAME='MLION_HH', &
           LONG_NAME='ml_ion_velocity_helper_hour_utc', UNITS='hour', &
           DIMS=MAPL_DimsHorzOnly, VLOCATION=MAPL_VLocationNone, _RC)

      ! ==================
      ! Import State
      ! ==================
      call MAPL_AddImportSpec(GC, SHORT_NAME='T', LONG_NAME='air_temperature', &
           UNITS='K', DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, RESTART=MAPL_RestartSkip, _RC)
      call MAPL_AddImportSpec(GC, SHORT_NAME='U', LONG_NAME='eastward_wind', &
           UNITS='m s-1', DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, RESTART=MAPL_RestartSkip, _RC)
      call MAPL_AddImportSpec(GC, SHORT_NAME='V', LONG_NAME='northward_wind', &
           UNITS='m s-1', DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, RESTART=MAPL_RestartSkip, _RC)
      call MAPL_AddImportSpec(GC, SHORT_NAME='PLE', LONG_NAME='air_pressure', &
           UNITS='Pa', DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationEdge, RESTART=MAPL_RestartSkip, _RC)
      call MAPL_AddImportSpec(GC, SHORT_NAME='ZLE', &
           LONG_NAME='geopotential_height', UNITS='m', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationEdge, &
           RESTART=MAPL_RestartSkip, _RC)
      call MAPL_AddImportSpec(GC, SHORT_NAME='PREF', &
           LONG_NAME='reference_air_pressure', UNITS='Pa', &
           DIMS=MAPL_DimsVertOnly, VLOCATION=MAPL_VLocationEdge, &
           RESTART=MAPL_RestartSkip, _RC)

      ! ==================
      ! Export State
      ! ==================
      call MAPL_AddExportSpec(GC, SHORT_NAME='UI_IONDRAG', &
           LONG_NAME='eastward_ion_velocity_from_ML', UNITS='m s-1', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='VI_IONDRAG', &
           LONG_NAME='northward_ion_velocity_from_ML', UNITS='m s-1', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='DUDT_IONDRAG', &
           LONG_NAME='eastward_wind_tendency_due_to_ion_drag', UNITS='m s-2', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='DVDT_IONDRAG', &
           LONG_NAME='northward_wind_tendency_due_to_ion_drag', UNITS='m s-2', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='DTDT_IONDRAG', &
           LONG_NAME='diagnostic_joule_heating_temperature_tendency', &
           UNITS='K s-1', DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='NE_IONDRAG', &
           LONG_NAME='electron_number_density_from_IRI', UNITS='m-3', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='OP_IONDRAG', &
           LONG_NAME='Oplus_number_density_from_IRI', UNITS='m-3', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='O2P_IONDRAG', &
           LONG_NAME='O2plus_number_density_from_IRI', UNITS='m-3', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='NOP_IONDRAG', &
           LONG_NAME='NOplus_number_density_from_IRI', UNITS='m-3', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='TI_IONDRAG', &
           LONG_NAME='ion_temperature_from_IRI', UNITS='K', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='TE_IONDRAG', &
           LONG_NAME='electron_temperature_from_IRI', UNITS='K', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='RHO_MSIS_IONDRAG', &
           LONG_NAME='neutral_mass_density_from_MSIS', UNITS='kg m-3', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='BMAG_IONDRAG', &
           LONG_NAME='IGRF_magnetic_field_magnitude', UNITS='T', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='SIGMAPED_IONDRAG', &
           LONG_NAME='Pedersen_conductivity', UNITS='S m-1', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='SIGMAHALL_IONDRAG', &
           LONG_NAME='Hall_conductivity', UNITS='S m-1', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='LXX_IONDRAG', &
           LONG_NAME='ion_drag_tensor_xx', UNITS='s-1', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='LYY_IONDRAG', &
           LONG_NAME='ion_drag_tensor_yy', UNITS='s-1', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='LXY_IONDRAG', &
           LONG_NAME='ion_drag_tensor_xy', UNITS='s-1', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)
      call MAPL_AddExportSpec(GC, SHORT_NAME='LYX_IONDRAG', &
           LONG_NAME='ion_drag_tensor_yx', UNITS='s-1', &
           DIMS=MAPL_DimsHorzVert, VLOCATION=MAPL_VLocationCenter, _RC)

      RETURN_(ESMF_SUCCESS)
   end subroutine register_state_specs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

      type(ESMF_GridComp), intent(inout) :: GC
      type(ESMF_State),    intent(inout) :: IMPORT
      type(ESMF_State),    intent(inout) :: EXPORT
      type(ESMF_Clock),    intent(inout) :: CLOCK
      integer, optional,   intent(out)   :: RC

      character(len=ESMF_MAXSTR)          :: IAm
      integer                             :: STATUS
      character(len=ESMF_MAXSTR)          :: COMP_NAME
      type (MAPL_MetaComp), pointer       :: MAPL
      type (wrap_)                        :: wrap
      type (GEOS_IonDragGridComp), pointer :: self

      Iam = 'Initialize'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC )
      Iam = trim(COMP_NAME) // Iam

      call MAPL_GetObjectFromGC ( GC, MAPL, _RC )
      call ESMF_UserCompGetInternalState(GC, 'GEOS_IonDragGridComp', wrap, _RC)
      self => wrap%ptr

      call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK, _RC )

      self%MLION_PYBRIDGE_INITIALIZED = .false.
      self%MLION_LAST_YEAR = -1
      self%MLION_LAST_DOY = -1
      self%MLION_LAST_HOUR = -1
      self%PLASMA_LAST_YEAR = -1
      self%PLASMA_LAST_DOY = -1
      self%PLASMA_LAST_HOUR = -1

      call msis_wrapper_init()

      RETURN_(ESMF_SUCCESS)
   end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

      type(ESMF_GridComp), intent(inout) :: GC
      type(ESMF_State),    intent(inout) :: IMPORT
      type(ESMF_State),    intent(inout) :: EXPORT
      type(ESMF_Clock),    intent(inout) :: CLOCK
      integer, optional,   intent(out)   :: RC

      character(len=ESMF_MAXSTR)          :: IAm
      integer                             :: STATUS
      character(len=ESMF_MAXSTR)          :: COMP_NAME
      type (MAPL_MetaComp), pointer       :: MAPL
      type (ESMF_Alarm)                   :: ALARM
      integer                             :: IM, JM, LM
      type (wrap_)                        :: wrap
      type (GEOS_IonDragGridComp), pointer :: self

      Iam = "Run"
      call ESMF_GridCompGet( GC, name=COMP_NAME, _RC )
      Iam = trim(COMP_NAME) // Iam

      call MAPL_GetObjectFromGC ( GC, MAPL, _RC)
      call ESMF_UserCompGetInternalState(GC, 'GEOS_IonDragGridComp', wrap, _RC)
      self => wrap%ptr


      call MAPL_Get(MAPL, IM=IM, JM=JM, LM=LM, RUNALARM=ALARM, _RC )

      if ( ESMF_AlarmIsRinging( ALARM ) ) then
         call MAPL_TimerOn(MAPL, "DRIVER")
         call IonDrag_Driver(_RC)
         call MAPL_TimerOff(MAPL, "DRIVER")
      endif

      RETURN_(ESMF_SUCCESS)

   contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine IonDrag_Driver(RC)
         integer, optional, intent(OUT) :: RC

         character(len=ESMF_MAXSTR)      :: IAm
         integer                         :: STATUS

         ! Import pointers used directly by the Fortran driver. PLE remains
         ! registered as an import because PythonBridge reads it directly.
         real, pointer :: T(:,:,:), U(:,:,:), V(:,:,:)
         real, pointer :: ZLE(:,:,:), PREF(:)

         ! Export pointers.
         real, pointer :: UI_IONDRAG(:,:,:), VI_IONDRAG(:,:,:)
         real, pointer :: DUDT_IONDRAG(:,:,:), DVDT_IONDRAG(:,:,:)
         real, pointer :: DTDT_IONDRAG(:,:,:), NE_IONDRAG(:,:,:)
         real, pointer :: OP_IONDRAG(:,:,:), O2P_IONDRAG(:,:,:)
         real, pointer :: NOP_IONDRAG(:,:,:), TI_IONDRAG(:,:,:)
         real, pointer :: TE_IONDRAG(:,:,:), RHO_MSIS_IONDRAG(:,:,:)
         real, pointer :: BMAG_IONDRAG(:,:,:), SIGMAPED_IONDRAG(:,:,:)
         real, pointer :: SIGMAHALL_IONDRAG(:,:,:)
         real, pointer :: LXX_IONDRAG(:,:,:), LYY_IONDRAG(:,:,:)
         real, pointer :: LXY_IONDRAG(:,:,:), LYX_IONDRAG(:,:,:)

         type (ESMF_State)        :: INTERNAL
         type (ESMF_Time)         :: CURRENT_TIME
         type (ESMF_TimeInterval) :: MODEL_TIMESTEP

         real, pointer, dimension(:,:) :: LATS_2D, LONS_2D
         real, pointer, dimension(:,:) :: MLION_LATS_2D, MLION_LONS_2D
         real, pointer, dimension(:,:) :: MLION_YY_2D, MLION_DOY_2D
         real, pointer, dimension(:,:) :: MLION_HH_2D
         real, allocatable :: LATS_DEG(:,:), LONS_DEG(:,:)

         integer :: IYEAR, DOY, HH, MN, SS
         integer :: MLION_HOUR
         real :: UT_HOUR
         real :: F107_DAILY_NOW, F107_81DAY_NOW
         real(kind=8) :: TIMESTEP_SECONDS_R8
         real :: TIMESTEP_SECONDS
         logical :: UPDATE_MLION, UPDATE_PLASMA
         logical :: RESET_PLASMA_CACHE

         real, allocatable :: alt_km(:,:,:)
         real, allocatable :: ui(:,:,:), vi(:,:,:)
         real, allocatable :: n_o_msis(:,:,:), n_o2_msis(:,:,:), n_n2_msis(:,:,:)
         real, allocatable :: rho_msis(:,:,:), cp_msis(:,:,:)
         real, allocatable :: drag_u(:,:,:), drag_v(:,:,:)
         real, allocatable :: joule_heating_wkg(:,:,:)
         real, allocatable :: sigma_pedersen(:,:,:), sigma_hall(:,:,:)
         real, allocatable :: lxx(:,:,:), lyy(:,:,:), lxy(:,:,:), lyx(:,:,:)

         integer :: i, j, k
         integer :: nlev_active
         real :: pref_mid_pa

         IAm = "IonDrag_Driver"

         call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, &
                       LATS=LATS_2D, LONS=LONS_2D, _RC)
         ! Explicit state pointers replace the ACG-generated GetPointer file.
         call MAPL_GetPointer(IMPORT, T,    'T',    _RC)
         call MAPL_GetPointer(IMPORT, U,    'U',    _RC)
         call MAPL_GetPointer(IMPORT, V,    'V',    _RC)
         call MAPL_GetPointer(IMPORT, ZLE,  'ZLE',  _RC)
         call MAPL_GetPointer(IMPORT, PREF, 'PREF', _RC)

         call MAPL_GetPointer(EXPORT, UI_IONDRAG, 'UI_IONDRAG', &
                              alloc=.true., _RC)
         call MAPL_GetPointer(EXPORT, VI_IONDRAG, 'VI_IONDRAG', &
                              alloc=.true., _RC)
         call MAPL_GetPointer(EXPORT, DUDT_IONDRAG, 'DUDT_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, DVDT_IONDRAG, 'DVDT_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, DTDT_IONDRAG, 'DTDT_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, NE_IONDRAG, 'NE_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, OP_IONDRAG, 'OP_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, O2P_IONDRAG, 'O2P_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, NOP_IONDRAG, 'NOP_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, TI_IONDRAG, 'TI_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, TE_IONDRAG, 'TE_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, RHO_MSIS_IONDRAG, &
                              'RHO_MSIS_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, BMAG_IONDRAG, 'BMAG_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, SIGMAPED_IONDRAG, &
                              'SIGMAPED_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, SIGMAHALL_IONDRAG, &
                              'SIGMAHALL_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, LXX_IONDRAG, 'LXX_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, LYY_IONDRAG, 'LYY_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, LXY_IONDRAG, 'LXY_IONDRAG', _RC)
         call MAPL_GetPointer(EXPORT, LYX_IONDRAG, 'LYX_IONDRAG', _RC)


         call ESMF_ClockGet(CLOCK, CurrTime=CURRENT_TIME, &
                            TimeStep=MODEL_TIMESTEP, _RC)
         call ESMF_TimeGet(CURRENT_TIME, YY=IYEAR, DayOfYear=DOY, &
                           H=HH, M=MN, S=SS, _RC)
         call ESMF_TimeIntervalGet(MODEL_TIMESTEP, &
                                   S_R8=TIMESTEP_SECONDS_R8, _RC)

         UT_HOUR = real(HH) + real(MN)/60.0 + real(SS)/3600.0
         TIMESTEP_SECONDS = real(TIMESTEP_SECONDS_R8)

         if (TIMESTEP_SECONDS <= 0.0) then
            error stop 'IONDRAG: model timestep must be positive'
         end if

         nlev_active = 0
         do k = 1, LM
            pref_mid_pa = 0.5 * (PREF(k) + PREF(k+1))
            if (pref_mid_pa <= IONDRAG_BOTTOM_PRESSURE_PA) then
               nlev_active = k
            else
               exit
            end if
         end do

         if (nlev_active < 1) then
            if (MAPL_am_I_root()) then
               print *, 'IONDRAG_ERROR: pressure cutoff selects no model levels:', &
                        IONDRAG_BOTTOM_PRESSURE_PA
            end if
            error stop 'IONDRAG: pressure cutoff selects no model levels'
         end if

         allocate(LATS_DEG(IM,JM), LONS_DEG(IM,JM))
         LATS_DEG = LATS_2D * (180.0/MAPL_PI)
         LONS_DEG = LONS_2D * (180.0/MAPL_PI)

         allocate(alt_km(IM, JM, nlev_active))
         allocate(ui(IM, JM, nlev_active), vi(IM, JM, nlev_active))
         allocate(n_o_msis(IM, JM, nlev_active))
         allocate(n_o2_msis(IM, JM, nlev_active))
         allocate(n_n2_msis(IM, JM, nlev_active))
         allocate(rho_msis(IM, JM, nlev_active))
         allocate(cp_msis(IM, JM, nlev_active))
         allocate(drag_u(IM, JM, nlev_active), drag_v(IM, JM, nlev_active))
         allocate(joule_heating_wkg(IM, JM, nlev_active))
         allocate(sigma_pedersen(IM, JM, nlev_active))
         allocate(sigma_hall(IM, JM, nlev_active))
         allocate(lxx(IM, JM, nlev_active), lyy(IM, JM, nlev_active))
         allocate(lxy(IM, JM, nlev_active), lyx(IM, JM, nlev_active))

         do k = 1, nlev_active
            do j = 1, JM
               do i = 1, IM
                  if (ieee_is_finite(ZLE(i,j,k-1)) .and. &
                      ieee_is_finite(ZLE(i,j,k))) then
                     alt_km(i,j,k) = 0.5 * &
                          (ZLE(i,j,k-1) + ZLE(i,j,k)) / 1000.0
                  else
                     alt_km(i,j,k) = FALLBACK_ALT_KM
                  end if
               end do
            end do
         end do

         MLION_HOUR = max(0, min(23, int(UT_HOUR)))

         UPDATE_MLION = &
              self%MLION_LAST_YEAR /= IYEAR .or. &
              self%MLION_LAST_DOY  /= DOY   .or. &
              self%MLION_LAST_HOUR /= MLION_HOUR

         if (UPDATE_MLION) then
            call MAPL_TimerOn(MAPL, "-MLION")

            call MAPL_GetPointer(INTERNAL, MLION_LATS_2D, 'MLION_LATS', _RC)
            call MAPL_GetPointer(INTERNAL, MLION_LONS_2D, 'MLION_LONS', _RC)
            call MAPL_GetPointer(INTERNAL, MLION_YY_2D,   'MLION_YY',   _RC)
            call MAPL_GetPointer(INTERNAL, MLION_DOY_2D,  'MLION_DOY',  _RC)
            call MAPL_GetPointer(INTERNAL, MLION_HH_2D,   'MLION_HH',   _RC)

            MLION_LATS_2D(:,:) = LATS_2D(:,:)
            MLION_LONS_2D(:,:) = LONS_2D(:,:)
            MLION_YY_2D(:,:)   = real(IYEAR)
            MLION_DOY_2D(:,:)  = real(DOY)
            MLION_HH_2D(:,:)   = real(MLION_HOUR)

            if (.not. self%MLION_PYBRIDGE_INITIALIZED) then
               call MAPL_pybridge_gcinit( &
                    "geos_mlionvel_driver", MAPL, IMPORT, EXPORT)
               self%MLION_PYBRIDGE_INITIALIZED = .true.
            end if

            call MAPL_pybridge_gcrun_with_internal( &
                 "geos_mlionvel_driver", MAPL, IMPORT, EXPORT, INTERNAL)

            self%MLION_LAST_YEAR = IYEAR
            self%MLION_LAST_DOY  = DOY
            self%MLION_LAST_HOUR = MLION_HOUR

            call MAPL_TimerOff(MAPL, "-MLION")
         end if

         if (nlev_active < LM) then
            UI_IONDRAG(:,:,nlev_active+1:LM) = 0.0
            VI_IONDRAG(:,:,nlev_active+1:LM) = 0.0
         end if

         ui = UI_IONDRAG(:,:,1:nlev_active)
         vi = VI_IONDRAG(:,:,1:nlev_active)

         ! Prepare the same hourly space-weather forcing used by MSIS and IRI.
         call msis_prepare_time(IYEAR, DOY, nint(UT_HOUR*3600.0))
         call msis_get_current_f107(F107_DAILY_NOW, F107_81DAY_NOW)

         ! Allocate persistent hourly IRI/IGRF caches on the local GEOS tile.
         RESET_PLASMA_CACHE = .false.
         if (.not. allocated(self%NE_CACHE)) then
            RESET_PLASMA_CACHE = .true.
         else if (size(self%NE_CACHE,1) /= IM .or. &
                  size(self%NE_CACHE,2) /= JM .or. &
                  size(self%NE_CACHE,3) /= nlev_active) then
            RESET_PLASMA_CACHE = .true.
         end if

         if (RESET_PLASMA_CACHE) then
            if (allocated(self%NE_CACHE)) then
               deallocate(self%NE_CACHE, self%ION_DENSITY_CACHE, &
                          self%TI_CACHE, self%TE_CACHE, &
                          self%BNORTH_CACHE, self%BEAST_CACHE, &
                          self%BDOWN_CACHE, self%BMAG_CACHE)
            end if

            allocate(self%NE_CACHE(IM,JM,nlev_active))
            allocate(self%ION_DENSITY_CACHE(N_MAJOR_ION_SPECIES,IM,JM,nlev_active))
            allocate(self%TI_CACHE(IM,JM,nlev_active))
            allocate(self%TE_CACHE(IM,JM,nlev_active))
            allocate(self%BNORTH_CACHE(IM,JM,nlev_active))
            allocate(self%BEAST_CACHE(IM,JM,nlev_active))
            allocate(self%BDOWN_CACHE(IM,JM,nlev_active))
            allocate(self%BMAG_CACHE(IM,JM,nlev_active))

            self%PLASMA_LAST_YEAR = -1
            self%PLASMA_LAST_DOY = -1
            self%PLASMA_LAST_HOUR = -1
         end if

         UPDATE_PLASMA = &
              self%PLASMA_LAST_YEAR /= IYEAR .or. &
              self%PLASMA_LAST_DOY  /= DOY   .or. &
              self%PLASMA_LAST_HOUR /= MLION_HOUR

         if (UPDATE_PLASMA) then
            call MAPL_TimerOn(MAPL, "-IRI")
            call get_iri_state( &
                 LATS_DEG, LONS_DEG, IYEAR, DOY, UT_HOUR, &
                 F107_DAILY_NOW, F107_81DAY_NOW, alt_km, &
                 self%NE_CACHE, self%ION_DENSITY_CACHE, &
                 self%TI_CACHE, self%TE_CACHE)
            call MAPL_TimerOff(MAPL, "-IRI")

            call MAPL_TimerOn(MAPL, "-IGRF")
            call get_igrf_field( &
                 LATS_DEG, LONS_DEG, IYEAR, DOY, alt_km, &
                 self%BNORTH_CACHE, self%BEAST_CACHE, &
                 self%BDOWN_CACHE, self%BMAG_CACHE)
            call MAPL_TimerOff(MAPL, "-IGRF")

            self%PLASMA_LAST_YEAR = IYEAR
            self%PLASMA_LAST_DOY = DOY
            self%PLASMA_LAST_HOUR = MLION_HOUR
         end if

         ! MSIS is evaluated every physics call. Its absolute O/O2/N2 number
         ! densities and resulting neutral mass density are used directly by
         ! the collision/conductivity calculation.
         call MAPL_TimerOn(MAPL, "-MSIS")
         call compute_msis_state( &
              IM, JM, nlev_active, IYEAR, DOY, UT_HOUR, &
              LATS_2D, LONS_2D, alt_km, &
              n_o_msis, n_o2_msis, n_n2_msis, &
              rho_msis, cp_msis, _RC)
         call MAPL_TimerOff(MAPL, "-MSIS")

         call MAPL_TimerOn(MAPL, "-DRAG")
         call compute_drag_fields( &
              ui, vi, U(:,:,1:nlev_active), V(:,:,1:nlev_active), &
              T(:,:,1:nlev_active), self%TI_CACHE, self%TE_CACHE, &
              self%ION_DENSITY_CACHE, &
              n_o_msis, n_o2_msis, n_n2_msis, rho_msis, &
              self%BNORTH_CACHE, self%BEAST_CACHE, self%BDOWN_CACHE, &
              TIMESTEP_SECONDS, drag_u, drag_v, joule_heating_wkg, &
              sigma_pedersen, sigma_hall, lxx, lyy, lxy, lyx)
         call MAPL_TimerOff(MAPL, "-DRAG")

         if (associated(DUDT_IONDRAG)) then
            DUDT_IONDRAG = 0.0
            DUDT_IONDRAG(:,:,1:nlev_active) = drag_u
         end if
         if (associated(DVDT_IONDRAG)) then
            DVDT_IONDRAG = 0.0
            DVDT_IONDRAG(:,:,1:nlev_active) = drag_v
         end if

         ! Diagnostic only. GEOS_PhysicsGridComp must continue to exclude this
         ! field from total heating while MLRADJH is active.
         if (associated(DTDT_IONDRAG)) then
            DTDT_IONDRAG = 0.0
            where (cp_msis > 0.0)
               DTDT_IONDRAG(:,:,1:nlev_active) = joule_heating_wkg / cp_msis
            elsewhere
               DTDT_IONDRAG(:,:,1:nlev_active) = 0.0
            end where
         end if

         if (associated(NE_IONDRAG)) then
            NE_IONDRAG = 0.0
            NE_IONDRAG(:,:,1:nlev_active) = self%NE_CACHE
         end if
         if (associated(OP_IONDRAG)) then
            OP_IONDRAG = 0.0
            OP_IONDRAG(:,:,1:nlev_active) = self%ION_DENSITY_CACHE(ION_OP,:,:,:)
         end if
         if (associated(O2P_IONDRAG)) then
            O2P_IONDRAG = 0.0
            O2P_IONDRAG(:,:,1:nlev_active) = self%ION_DENSITY_CACHE(ION_O2P,:,:,:)
         end if
         if (associated(NOP_IONDRAG)) then
            NOP_IONDRAG = 0.0
            NOP_IONDRAG(:,:,1:nlev_active) = self%ION_DENSITY_CACHE(ION_NOP,:,:,:)
         end if
         if (associated(TI_IONDRAG)) then
            TI_IONDRAG = 0.0
            TI_IONDRAG(:,:,1:nlev_active) = self%TI_CACHE
         end if
         if (associated(TE_IONDRAG)) then
            TE_IONDRAG = 0.0
            TE_IONDRAG(:,:,1:nlev_active) = self%TE_CACHE
         end if
         if (associated(RHO_MSIS_IONDRAG)) then
            RHO_MSIS_IONDRAG = 0.0
            RHO_MSIS_IONDRAG(:,:,1:nlev_active) = rho_msis
         end if
         if (associated(BMAG_IONDRAG)) then
            BMAG_IONDRAG = 0.0
            BMAG_IONDRAG(:,:,1:nlev_active) = self%BMAG_CACHE
         end if
         if (associated(SIGMAPED_IONDRAG)) then
            SIGMAPED_IONDRAG = 0.0
            SIGMAPED_IONDRAG(:,:,1:nlev_active) = sigma_pedersen
         end if
         if (associated(SIGMAHALL_IONDRAG)) then
            SIGMAHALL_IONDRAG = 0.0
            SIGMAHALL_IONDRAG(:,:,1:nlev_active) = sigma_hall
         end if
         if (associated(LXX_IONDRAG)) then
            LXX_IONDRAG = 0.0
            LXX_IONDRAG(:,:,1:nlev_active) = lxx
         end if
         if (associated(LYY_IONDRAG)) then
            LYY_IONDRAG = 0.0
            LYY_IONDRAG(:,:,1:nlev_active) = lyy
         end if
         if (associated(LXY_IONDRAG)) then
            LXY_IONDRAG = 0.0
            LXY_IONDRAG(:,:,1:nlev_active) = lxy
         end if
         if (associated(LYX_IONDRAG)) then
            LYX_IONDRAG = 0.0
            LYX_IONDRAG(:,:,1:nlev_active) = lyx
         end if

         deallocate(LATS_DEG, LONS_DEG, alt_km, ui, vi, &
                    n_o_msis, n_o2_msis, n_n2_msis, &
                    rho_msis, cp_msis, drag_u, drag_v, &
                    joule_heating_wkg, sigma_pedersen, sigma_hall, &
                    lxx, lyy, lxy, lyx)

         RETURN_(ESMF_SUCCESS)

      end subroutine IonDrag_Driver

   end subroutine RUN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine compute_msis_state(IM, JM, NLEV_ACTIVE, IYEAR, DOY, &
                                 UT_HOUR, LATS_2D, LONS_2D, alt_km_in, &
                                 n_o_out, n_o2_out, n_n2_out, &
                                 rho_out, cp_out, RC)
      ! Diagnose the neutral state required by ion drag from MSIS.
      !
      ! n_o_out/n_o2_out/n_n2_out : absolute species number densities [m-3]
      ! rho_out                   : O + O2 + N2 mass density [kg m-3]
      ! cp_out                    : mixture specific heat [J kg-1 K-1]

      integer, intent(in) :: IM, JM, NLEV_ACTIVE
      integer, intent(in) :: IYEAR, DOY
      real, intent(in) :: UT_HOUR
      real, pointer, dimension(:,:), intent(in) :: LATS_2D, LONS_2D
      real, intent(in) :: alt_km_in(:,:,:)
      real, intent(out) :: n_o_out(:,:,:), n_o2_out(:,:,:), n_n2_out(:,:,:)
      real, intent(out) :: rho_out(:,:,:), cp_out(:,:,:)
      integer, optional, intent(OUT) :: RC

      character(len=ESMF_MAXSTR) :: IAm
      integer :: STATUS

      real, parameter :: AMU_KG = 1.66053906660e-27
      real, parameter :: CM3_TO_M3 = 1.0e6
      real, parameter :: MASS_O  = 16.0
      real, parameter :: MASS_N2 = 28.0
      real, parameter :: MASS_O2 = 32.0

      real(4) :: O_out, N2_out, O2_out, T_msis_out
      real(4) :: alt_r4, glat_r4, glong_r4, stl_r4
      real :: r_mix, cp_mix, cv_mix, kappa_mix
      real :: phi_o, phi_n2, phi_o2
      integer :: i, j, k

      IAm = "compute_msis_state"

      n_o_out = 0.0
      n_o2_out = 0.0
      n_n2_out = 0.0
      rho_out = 0.0
      cp_out = 0.0

      do k = 1, NLEV_ACTIVE
         do j = 1, JM
            do i = 1, IM
               alt_r4 = real(alt_km_in(i,j,k), kind=4)
               glat_r4 = real(LATS_2D(i,j) * (180.0/MAPL_PI), kind=4)
               glong_r4 = real(LONS_2D(i,j) * (180.0/MAPL_PI), kind=4)
               stl_r4 = modulo(real(UT_HOUR, kind=4) + glong_r4/15.0_4, 24.0_4)

               call msis_point(IYEAR, DOY, nint(UT_HOUR*3600.0), &
                    alt_r4, glat_r4, glong_r4, stl_r4, &
                    O_out, N2_out, O2_out, T_msis_out)

               if (.not. ieee_is_finite(O_out) .or. &
                   .not. ieee_is_finite(N2_out) .or. &
                   .not. ieee_is_finite(O2_out)) then
                  O_out = 0.0_4
                  N2_out = 0.0_4
                  O2_out = 0.0_4
               end if

               n_o_out(i,j,k) = max(real(O_out), 0.0) * CM3_TO_M3
               n_o2_out(i,j,k) = max(real(O2_out), 0.0) * CM3_TO_M3
               n_n2_out(i,j,k) = max(real(N2_out), 0.0) * CM3_TO_M3

               rho_out(i,j,k) = &
                    (max(real(O_out),0.0)*MASS_O + &
                     max(real(N2_out),0.0)*MASS_N2 + &
                     max(real(O2_out),0.0)*MASS_O2) * AMU_KG * CM3_TO_M3

               call mlt_mixture_thermo_from_number_density( &
                    max(real(O_out),0.0), max(real(N2_out),0.0), &
                    max(real(O2_out),0.0), r_mix, cp_mix, cv_mix, &
                    kappa_mix, phi_o, phi_n2, phi_o2)
               cp_out(i,j,k) = cp_mix
            end do
         end do
      end do

      RETURN_(ESMF_SUCCESS)
   end subroutine compute_msis_state

end module GEOS_IonDragGridCompMod
