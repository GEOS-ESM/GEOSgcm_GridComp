#include "MAPL_Generic.h"

module GEOS_IonDragGridCompMod

   !BOP

   ! !MODULE: GEOS_IonDrag -- A Module to compute ion drag forcing on the
   ! neutral atmosphere in the thermosphere (GEOS-MLT)

   ! !DESCRIPTION:
   !
   ! IonDrag is a light-weight gridded component that computes the
   ! momentum drag and frictional heating tendencies on the neutral wind and
   ! temperature fields due to collisions with ions above a configurable
   ! reference-pressure cutoff. Ion densities and species fractions are
   ! obtained from IRI (evaluated at real per-column geometric altitude,
   ! derived from the GEOS geopotential-height field ZLE); neutral number
   ! density, mass density, and mixture heat capacity come from MSIS
   ! (via the same msis_wrapper module used on the dynamics
   ! side); ion winds are currently a constant placeholder (to be replaced
   ! by an ML model output after initial commits). Like GEOSgwd_GridComp,
   ! this component only exports tendencies (DUDT_IONDRAG, DVDT_IONDRAG,
   ! DTDT_IONDRAG) -- it does not mutate U/V/T directly. Those tendencies
   ! are collected by the parent GEOS_PhysicsGridComp into the combined
   ! physics DUDT/DVDT/DTDT applied by the dynamics.
   !

   ! !USES:

   use ESMF
   use MAPL

   use iri_input_module, only: get_iri_densities, N_ION_SPECIES
   use ion_drag_module,  only: compute_drag_fields
   use msis_wrapper,     only: msis_wrapper_init, msis_prepare_time, msis_point, &
                               msis_get_current_f107
   use calc_gas_specific_heat_mlt_mod, only: mlt_mixture_thermo_from_number_density
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   
   implicit none
   private

   ! !PUBLIC MEMBER FUNCTIONS:

   public SetServices

   !EOP

   type :: GEOS_IonDragGridComp
      logical :: IONDRAG_ON
      real    :: BOTTOM_PRESSURE_PA      ! Maximum midpoint pressure for ion drag [Pa]
      real    :: HBEG, HEND, HSTEP       ! IRI altitude profile range/step (km)
      real    :: TEST_UI_MS, TEST_VI_MS  ! placeholder constant ion winds
   end type GEOS_IonDragGridComp

   type wrap_
      type (GEOS_IonDragGridComp), pointer :: PTR
   end type wrap_

contains

   !BOP
   ! !IROUTINE: SetServices -- Sets ESMF services for this component

   ! !INTERFACE:
   subroutine SetServices ( GC, RC )

      ! !ARGUMENTS:
      type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
      integer, optional                  :: RC  ! return code

      !EOP

      character(len=ESMF_MAXSTR)              :: IAm
      integer                                 :: STATUS
      character(len=ESMF_MAXSTR)              :: COMP_NAME
      type (MAPL_MetaComp),     pointer       :: MAPL

      type (wrap_)                                :: wrap
      type (GEOS_IonDragGridComp), pointer         :: self

      ! Begin...

      Iam = 'SetServices'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC )
      Iam = trim(COMP_NAME) // Iam

      !   Wrap internal state for storing in GC
      !   -------------------------------------
      allocate (self, _STAT)
      wrap%ptr => self

      ! Set the Run entry point
      ! -----------------------

      call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_INITIALIZE,  Initialize,  _RC)
      call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,  Run,  _RC)

      call MAPL_GetObjectFromGC ( GC, MAPL, _RC )

      ! Set the state variable specs (generated from IonDrag_StateSpecs.rc).
      ! ---------------------------------------------------------------------
#include "IonDrag_Import___.h"
#include "IonDrag_Export___.h"
#include "IonDrag_Internal___.h"

      ! Set the Profiling timers
      ! ------------------------

      call MAPL_TimerAdd(GC,    name="DRIVER"     ,_RC)
      call MAPL_TimerAdd(GC,    name="-IRI"       ,_RC)
      call MAPL_TimerAdd(GC,    name="-MSIS"      ,_RC)
      call MAPL_TimerAdd(GC,    name="-DRAG"      ,_RC)

      !   Store internal state in GC
      !   --------------------------
      call ESMF_UserCompSetInternalState ( GC, 'GEOS_IonDragGridComp', wrap, _RC )

      ! Set generic init and final methods
      ! ----------------------------------

      call MAPL_GenericSetServices    ( gc, _RC)

      RETURN_(ESMF_SUCCESS)

   end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !BOP
   ! !IROUTINE: Initialize -- Initialize method for the IonDrag Gridded Component

   ! !INTERFACE:
   subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

      ! !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
      integer, optional,   intent(  out) :: RC     ! Error code

      !EOP

      character(len=ESMF_MAXSTR)              :: IAm
      integer                                 :: STATUS
      character(len=ESMF_MAXSTR)              :: COMP_NAME

      type (MAPL_MetaComp),      pointer  :: MAPL

      type (wrap_) :: wrap
      type (GEOS_IonDragGridComp), pointer :: self

      ! Begin...

      Iam = 'Initialize'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC )
      Iam = trim(COMP_NAME) // Iam

      call MAPL_GetObjectFromGC ( GC, MAPL, _RC )

      call ESMF_UserCompGetInternalState(GC, 'GEOS_IonDragGridComp', wrap, _RC)
      self => wrap%ptr

      call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK, _RC )

      ! Resource config
      ! ---------------
      call MAPL_GetResource( MAPL, self%IONDRAG_ON, Label="IONDRAG_ON:", default=.true., _RC)
      call MAPL_GetResource( MAPL, self%BOTTOM_PRESSURE_PA, &
           Label="IONDRAG_BOTTOM_PRESSURE_PA:", default=1.0, _RC)
      call MAPL_GetResource( MAPL, self%HBEG, Label="IRI_HBEG:", default=70.0, _RC)
      call MAPL_GetResource( MAPL, self%HEND,          Label="IRI_HEND:",      default=250.0,   _RC)
      call MAPL_GetResource( MAPL, self%HSTEP,         Label="IRI_HSTEP:",     default=10.0,    _RC)
      call MAPL_GetResource( MAPL, self%TEST_UI_MS,    Label="TEST_UI_MS:",    default=50.0,    _RC)
      call MAPL_GetResource( MAPL, self%TEST_VI_MS,    Label="TEST_VI_MS:",    default=0.0,     _RC)

      ! Validate the pressure and IRI configuration before Run.
      if (self%BOTTOM_PRESSURE_PA <= 0.0) then
         error stop 'IONDRAG: IONDRAG_BOTTOM_PRESSURE_PA must be positive'
      end if
      if (self%HSTEP <= 0.0) then
         error stop 'IONDRAG: IRI_HSTEP must be positive'
      end if
      if (self%HEND < self%HBEG) then
         error stop 'IONDRAG: IRI_HEND must be >= IRI_HBEG'
      end if

      ! Initialize MSIS only when ion drag is enabled.
      if (self%IONDRAG_ON) then
         call msis_wrapper_init()
      end if

      RETURN_(ESMF_SUCCESS)
   end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !BOP
   ! !IROUTINE: RUN -- Run method for the IonDrag component

   ! !INTERFACE:
   subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

      ! !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
      integer, optional,   intent(  out) :: RC     ! Error code

      !EOP

      character(len=ESMF_MAXSTR)          :: IAm
      integer                             :: STATUS
      character(len=ESMF_MAXSTR)          :: COMP_NAME

      type (MAPL_MetaComp),     pointer   :: MAPL
      type (ESMF_Alarm       )            :: ALARM

      integer                             :: IM, JM, LM

      type (wrap_) :: wrap
      type (GEOS_IonDragGridComp), pointer :: self

      ! Begin...

      Iam = "Run"
      call ESMF_GridCompGet( GC, name=COMP_NAME, _RC )
      Iam = trim(COMP_NAME) // Iam

      call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

      call ESMF_UserCompGetInternalState(GC, 'GEOS_IonDragGridComp', wrap, _RC)
      self => wrap%ptr

      if (.not. self%IONDRAG_ON) then
         RETURN_(ESMF_SUCCESS)
      end if

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

#include "IonDrag_DeclarePointer___.h"

         type (ESMF_State) :: INTERNAL
         type (ESMF_Time)  :: CURRENT_TIME

         real, pointer, dimension(:,:) :: LATS_2D, LONS_2D
         real, allocatable :: LATS_DEG(:,:), LONS_DEG(:,:)

         integer :: IYEAR, DOY, HH, MN, SS
         real    :: UT_HOUR
         real    :: F107_DAILY_NOW, F107_81DAY_NOW

         real, parameter :: FALLBACK_ALT_KM = 220.0

         real, allocatable :: alt_km(:,:,:)
         real, allocatable :: ui(:,:,:), vi(:,:,:)
         real, allocatable :: ne_m3(:,:,:)
         real, allocatable :: species_fraction(:,:,:,:)
         real, allocatable :: n_msis(:,:,:), rho_msis(:,:,:), cp_msis(:,:,:)
         real, allocatable :: drag_u(:,:,:), drag_v(:,:,:)
         real, allocatable :: frictional_heating(:,:,:)

         integer :: i, j, k
         integer :: nlev_active
         real    :: pref_mid_pa

         IAm = "IonDrag_Driver"

         call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, LATS=LATS_2D, LONS=LONS_2D, _RC)
#include "IonDrag_GetPointer___.h"

         ! Current model time -> year, day-of-year, UT hour.
         ! ESMF_TimeGet's DayOfYear argument does this natively -- no custom
         ! calendar helper needed (matches GEOS_SolarGridComp.F90's pattern).
         call ESMF_ClockGet(CLOCK, CurrTime=CURRENT_TIME, _RC)
         call ESMF_TimeGet(CURRENT_TIME, YY=IYEAR, DayOfYear=DOY, H=HH, M=MN, S=SS, _RC)
         UT_HOUR = real(HH) + real(MN)/60.0 + real(SS)/3600.0

         ! Determine the active ion-drag domain from the reference-pressure
         ! grid. Level 1 is the model-top layer and pressure increases downward.
         ! Using PREF makes the cutoff independent of the number of model levels.
         nlev_active = 0
         do k = 1, LM
            pref_mid_pa = 0.5 * (PREF(k) + PREF(k+1))
            if (pref_mid_pa <= self%BOTTOM_PRESSURE_PA) then
               nlev_active = k
            else
               exit
            end if
         end do

         if (nlev_active < 1) then
            if (MAPL_am_I_root()) then
               print *, 'IONDRAG_ERROR: pressure cutoff selects no model levels:', &
                        self%BOTTOM_PRESSURE_PA
            end if
            error stop 'IONDRAG: pressure cutoff selects no model levels'
         end if

         ! Cubed-sphere latitude/longitude are two-dimensional fields and are
         ! not separable latitude and longitude axes. Keep each (i,j) pair.
         allocate(LATS_DEG(IM,JM), LONS_DEG(IM,JM))
         LATS_DEG = LATS_2D * (180.0/MAPL_PI)
         LONS_DEG = LONS_2D * (180.0/MAPL_PI)

         allocate(alt_km(IM, JM, nlev_active))
         allocate(ui(IM, JM, nlev_active), vi(IM, JM, nlev_active))
         allocate(ne_m3(IM, JM, nlev_active))
         allocate(species_fraction(N_ION_SPECIES, IM, JM, nlev_active))
         allocate(n_msis(IM, JM, nlev_active))
         allocate(rho_msis(IM, JM, nlev_active))
         allocate(cp_msis(IM, JM, nlev_active))
         allocate(drag_u(IM, JM, nlev_active), &
                  drag_v(IM, JM, nlev_active))
         allocate(frictional_heating(IM, JM, nlev_active))

         ! ZLE is geopotential height at model interfaces in meters. Use the
         ! layer midpoint as the altitude supplied to IRI and MSIS. Guard
         ! invalid top-edge values before calling either empirical model.
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

         ! Step 1: Ion winds: replace with Andrew's ML model output.
         ui = self%TEST_UI_MS
         vi = self%TEST_VI_MS
         ! MSIS space-weather indices must be prepared once per timestep,
         ! before any msis_point calls and before the IRI call,
         ! since IRI reuses these same real F10.7/F10.7A values rather than
         ! static .rc placeholders.
         call msis_prepare_time(IYEAR, DOY, nint(UT_HOUR*3600.0))
         call msis_get_current_f107(F107_DAILY_NOW, F107_81DAY_NOW)

         ! Step 2: IRI ion densities / species fractions. Latitude and
         ! longitude are paired two-dimensional cubed-sphere coordinates, and
         ! alt_km contains the requested GEOS model-level altitude at each
         ! horizontal grid point.
         call MAPL_TimerOn(MAPL, "-IRI")
         call get_iri_densities( &
              LATS_DEG, LONS_DEG, IYEAR, DOY, UT_HOUR, &
              F107_DAILY_NOW, F107_81DAY_NOW, &
              alt_km, self%HBEG, self%HEND, self%HSTEP, &
              ne_m3, species_fraction)
         call MAPL_TimerOff(MAPL, "-IRI")

         ! Step 3: Diagnose neutral number density, mass density, and Cp from
         ! the same MSIS O/N2/O2 composition used by GEOS-MLT thermodynamics.
         call MAPL_TimerOn(MAPL, "-MSIS")
         call compute_msis_state( &
              IM, JM, nlev_active, IYEAR, DOY, UT_HOUR, &
              LATS_2D, LONS_2D, alt_km, n_msis, rho_msis, cp_msis, _RC)
         call MAPL_TimerOff(MAPL, "-MSIS")

         ! Step 4: Ion drag physics -- returns tendencies only.
         call MAPL_TimerOn(MAPL, "-DRAG")
         call compute_drag_fields(ui, vi, ne_m3, species_fraction, &
                                   U(:,:,1:nlev_active), &
                                   V(:,:,1:nlev_active), &
                                   n_msis, rho_msis, &
                                   T(:,:,1:nlev_active), &
                                   drag_u, drag_v, frictional_heating)
         call MAPL_TimerOff(MAPL, "-DRAG")

         ! Step 5: Populate exports ONLY -- U/V/T are read-only imports here.
         ! Tendencies are collected by GEOS_PhysicsGridComp into the combined
         ! physics DUDT/DVDT/DTDT applied by the dynamics (same pattern as
         ! GEOSgwd_GridComp's DUDT/DVDT/DTDT and GEOS_SolarGridComp's MLRADJH).
         if (associated(DUDT_IONDRAG)) then
            DUDT_IONDRAG = 0.0
            DUDT_IONDRAG(:,:,1:nlev_active) = drag_u
         end if
         if (associated(DVDT_IONDRAG)) then
            DVDT_IONDRAG = 0.0
            DVDT_IONDRAG(:,:,1:nlev_active) = drag_v
         end if
         if (associated(DTDT_IONDRAG)) then
            DTDT_IONDRAG = 0.0
            where (rho_msis > 0.0 .and. cp_msis > 0.0)
               DTDT_IONDRAG(:,:,1:nlev_active) = &
                    frictional_heating / (rho_msis * cp_msis)
            elsewhere
               DTDT_IONDRAG(:,:,1:nlev_active) = 0.0
            end where
         end if
         if (associated(NE_IONDRAG)) then
            NE_IONDRAG = 0.0
            NE_IONDRAG(:,:,1:nlev_active) = ne_m3
         end if
         if (associated(RHO_MSIS_IONDRAG)) then
            RHO_MSIS_IONDRAG = 0.0
            RHO_MSIS_IONDRAG(:,:,1:nlev_active) = rho_msis
         end if

         deallocate(LATS_DEG, LONS_DEG, alt_km, ui, vi, ne_m3, &
                    species_fraction, n_msis, rho_msis, cp_msis, &
                    drag_u, drag_v, frictional_heating)

         RETURN_(ESMF_SUCCESS)

      end subroutine IonDrag_Driver

   end subroutine RUN

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine compute_msis_state(IM, JM, NLEV_ACTIVE, IYEAR, DOY, &
                                 UT_HOUR, LATS_2D, LONS_2D, alt_km_in, &
                                 n_out, rho_out, cp_out, RC)
      ! Diagnose the neutral state required by ion drag from MSIS.
      !
      ! n_out   : total O + N2 + O2 number density [m-3]
      ! rho_out : O + N2 + O2 mass density [kg m-3]
      ! cp_out  : mixture specific heat at constant pressure [J kg-1 K-1]
      integer, intent(in) :: IM, JM, NLEV_ACTIVE
      integer, intent(in) :: IYEAR, DOY
      real,    intent(in) :: UT_HOUR
      real, pointer, dimension(:,:), intent(in) :: LATS_2D, LONS_2D
      real, intent(in)  :: alt_km_in(:,:,:)
      real, intent(out) :: n_out(:,:,:), rho_out(:,:,:), cp_out(:,:,:)
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

      do k = 1, NLEV_ACTIVE
         do j = 1, JM
            do i = 1, IM
               alt_r4   = real(alt_km_in(i,j,k), kind=4)
               glat_r4  = real(LATS_2D(i,j) * (180.0/MAPL_PI), kind=4)
               glong_r4 = real(LONS_2D(i,j) * (180.0/MAPL_PI), kind=4)
               stl_r4 = modulo( &
                    real(UT_HOUR, kind=4) + glong_r4/15.0_4, 24.0_4)

               call msis_point(IYEAR, DOY, nint(UT_HOUR*3600.0), &
                    alt_r4, glat_r4, glong_r4, stl_r4, &
                    O_out, N2_out, O2_out, T_msis_out)

               if (.not. ieee_is_finite(O_out) .or. &
                   .not. ieee_is_finite(N2_out) .or. &
                   .not. ieee_is_finite(O2_out)) then
                  O_out  = 0.0_4
                  N2_out = 0.0_4
                  O2_out = 0.0_4
               end if

               n_out(i,j,k) = &
                    (real(O_out) + real(N2_out) + real(O2_out)) * CM3_TO_M3

               rho_out(i,j,k) = &
                    (real(O_out)*MASS_O + real(N2_out)*MASS_N2 + &
                     real(O2_out)*MASS_O2) * AMU_KG * CM3_TO_M3

               call mlt_mixture_thermo_from_number_density( &
                    real(O_out), real(N2_out), real(O2_out), &
                    r_mix, cp_mix, cv_mix, kappa_mix, &
                    phi_o, phi_n2, phi_o2)
               cp_out(i,j,k) = cp_mix
            end do
         end do
      end do

      RETURN_(ESMF_SUCCESS)
   end subroutine compute_msis_state

end module GEOS_IonDragGridCompMod
