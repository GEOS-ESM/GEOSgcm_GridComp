#include "MAPL_Generic.h"

module GEOS_IonDragGridCompMod

   !BOP

   ! !MODULE: GEOS_IonDrag -- A Module to compute ion drag forcing on the
   ! neutral atmosphere in the thermosphere (GEOS-MLT)

   ! !DESCRIPTION:
   !
   ! This Ion drag implementation is a light-weight gridded component that computes the
   ! momentum drag and frictional heating tendencies on the neutral wind and
   ! temperature fields due to collisions with ions, over the top
   ! NLEV_IONDRAG model levels. Ion densities and species fractions are
   ! obtained from IRI (evaluated at real per-column geometric altitude,
   ! derived from geopotential height GZ the same way mol_mom_diff_mod
   ! derives altitude on the dynamics side); neutral mass density comes
   ! from MSIS (via the same msis_wrapper module used on the dynamics
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
   use msis_wrapper,     only: msis_wrapper_init, msis_prepare_time, msis_point

   implicit none
   private

   ! !PUBLIC MEMBER FUNCTIONS:

   public SetServices

   !EOP

   type :: GEOS_IonDragGridComp
      logical :: IONDRAG_ON
      integer :: NLEV_IONDRAG
      real    :: MEAN_MASS
      real    :: HBEG, HEND, HSTEP       ! IRI altitude profile range/step (km)
      real    :: F107_DAILY, F107_81DAY  ! climatological placeholders for v1
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
      call MAPL_GetResource( MAPL, self%IONDRAG_ON,    Label="IONDRAG_ON:",    default=.true.,  _RC)
      call MAPL_GetResource( MAPL, self%NLEV_IONDRAG,  Label="NLEV_IONDRAG:",  default=10,      _RC)
      call MAPL_GetResource( MAPL, self%MEAN_MASS,     Label="MEAN_MASS:",     default=16.0,    _RC)
      call MAPL_GetResource( MAPL, self%HBEG,          Label="IRI_HBEG:",      default=100.0,   _RC)
      call MAPL_GetResource( MAPL, self%HEND,          Label="IRI_HEND:",      default=700.0,   _RC)
      call MAPL_GetResource( MAPL, self%HSTEP,         Label="IRI_HSTEP:",     default=10.0,    _RC)
      call MAPL_GetResource( MAPL, self%F107_DAILY,    Label="F107_DAILY:",    default=150.0,   _RC)
      call MAPL_GetResource( MAPL, self%F107_81DAY,    Label="F107_81DAY:",    default=150.0,   _RC)
      call MAPL_GetResource( MAPL, self%TEST_UI_MS,    Label="TEST_UI_MS:",    default=50.0,    _RC)
      call MAPL_GetResource( MAPL, self%TEST_VI_MS,    Label="TEST_VI_MS:",    default=0.0,     _RC)

      ! Initialize MSIS (loads msis21.parm and F107_ap_appended.txt).
      ! Must happen once, here, not in Run.
      call msis_wrapper_init()

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
         type (ESMF_TimeInterval) :: TINT
         real(ESMF_KIND_R8) :: DT_R8
         real :: DT_PHYSICS

         real, pointer, dimension(:,:) :: LATS_2D, LONS_2D
         real, allocatable :: LATS_DEG(:), LONS_DEG(:)

         integer :: IYEAR, DOY, HH, MN, SS
         real    :: UT_HOUR

         real, parameter :: GRAV0 = 9.80665

         real, allocatable :: alt_km(:,:,:)   ! (IM, JM, NLEV_IONDRAG) -- per column
         real, allocatable :: ui(:,:,:), vi(:,:,:)
         real, allocatable :: ne_m3(:,:,:)
         real, allocatable :: species_fraction(:,:,:,:)
         real, allocatable :: rho_msis(:,:,:)
         real, allocatable :: drag_u(:,:,:), drag_v(:,:,:)
         real, allocatable :: frictional_heating(:,:,:)

         integer :: i, j, k

         IAm = "IonDrag_Driver"

         call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, LATS=LATS_2D, LONS=LONS_2D, _RC)
#include "IonDrag_GetPointer___.h"

         ! Time step
         ! ---------
         call ESMF_AlarmGet( ALARM, ringInterval=TINT, _RC)
         call ESMF_TimeIntervalGet(TINT, S_R8=DT_R8, _RC)
         DT_PHYSICS = DT_R8

         ! Current model time -> year, day-of-year, UT hour.
         ! calendar helper needed (matches GEOS_SolarGridComp.F90's pattern).
         call ESMF_ClockGet(CLOCK, CurrTime=CURRENT_TIME, _RC)
         call ESMF_TimeGet(CURRENT_TIME, YY=IYEAR, DayOfYear=DOY, H=HH, M=MN, S=SS, _RC)
         UT_HOUR = real(HH) + real(MN)/60.0 + real(SS)/3600.0

         ! Flatten lat/lon (radians -> degrees) for the per-column IRI call.
         allocate(LATS_DEG(IM*JM), LONS_DEG(IM*JM))
         LATS_DEG = reshape(LATS_2D, (/IM*JM/)) * (180.0/MAPL_PI)
         LONS_DEG = reshape(LONS_2D, (/IM*JM/)) * (180.0/MAPL_PI)

         allocate(alt_km(IM, JM, self%NLEV_IONDRAG))
         allocate(ui(IM, JM, self%NLEV_IONDRAG), vi(IM, JM, self%NLEV_IONDRAG))
         allocate(ne_m3(IM, JM, self%NLEV_IONDRAG))
         allocate(species_fraction(N_ION_SPECIES, IM, JM, self%NLEV_IONDRAG))
         allocate(rho_msis(IM, JM, self%NLEV_IONDRAG))
         allocate(drag_u(IM, JM, self%NLEV_IONDRAG), drag_v(IM, JM, self%NLEV_IONDRAG))
         allocate(frictional_heating(IM, JM, self%NLEV_IONDRAG))

         ! Real per-column altitude (km) for the top NLEV_IONDRAG levels,
         ! derived from geopotential height at interfaces (GZ), following
         ! the same approach as mol_mom_diff_mod on the dynamics side. GZ
         ! interfaces (k, k+1) unambiguously bound layer k regardless of
         ! top-down/bottom-up orientation, since we take the midpoint.
         do k = 1, self%NLEV_IONDRAG
            do j = 1, JM
               do i = 1, IM
                  alt_km(i,j,k) = 0.5 * ( GZ(i,j,k) + GZ(i,j,k+1) ) / GRAV0 / 1000.0
               end do
            end do
         end do

         ! Step 1: Ion winds -- placeholder constants for testing.
         ! TODO: replace with Andrew's ML model.
         ui = self%TEST_UI_MS
         vi = self%TEST_VI_MS

         ! Step 2: IRI ion densities / species fractions (real per-column
         ! altitude via alt_km).
         call MAPL_TimerOn(MAPL, "-IRI")
         call get_iri_densities(LATS_DEG, LONS_DEG, IYEAR, DOY, UT_HOUR, &
                                 self%F107_DAILY, self%F107_81DAY, &
                                 alt_km, self%HBEG, self%HEND, self%HSTEP, &
                                 ne_m3, species_fraction)
         call MAPL_TimerOff(MAPL, "-IRI")

         ! Step 3: MSIS neutral mass density (top NLEV_IONDRAG levels),
         ! also evaluated at real per-column altitude via alt_km.
         ! msis_prepare_time must be called once per timestep before any
         ! msis_point calls (msis_point errors out otherwise).
         call MAPL_TimerOn(MAPL, "-MSIS")
         call msis_prepare_time(IYEAR, DOY, nint(UT_HOUR*3600.0))
         call compute_msis_density(IM, JM, self%NLEV_IONDRAG, IYEAR, DOY, UT_HOUR, LATS_2D, LONS_2D, alt_km, rho_msis, _RC)
         call MAPL_TimerOff(MAPL, "-MSIS")

         ! Step 4: Ion drag physics -- returns tendencies only.
         call MAPL_TimerOn(MAPL, "-DRAG")
         call compute_drag_fields(ui, vi, LATS_DEG, ne_m3, species_fraction, &
                                   U(:,:,1:self%NLEV_IONDRAG), &
                                   V(:,:,1:self%NLEV_IONDRAG), &
                                   rho_msis, &
                                   T(:,:,1:self%NLEV_IONDRAG), &
                                   self%MEAN_MASS, drag_u, drag_v, frictional_heating)
         call MAPL_TimerOff(MAPL, "-DRAG")

         ! Step 5: Populate exports ONLY -- U/V/T are read-only imports here.
         ! Tendencies are collected by GEOS_PhysicsGridComp into the combined
         ! physics DUDT/DVDT/DTDT applied by the dynamics (same pattern as
         ! GEOSgwd_GridComp's DUDT/DVDT/DTDT and GEOS_SolarGridComp's MLRADJH).
         if (associated(DUDT_IONDRAG)) DUDT_IONDRAG = drag_u
         if (associated(DVDT_IONDRAG)) DVDT_IONDRAG = drag_v
         if (associated(DTDT_IONDRAG)) then
            where (rho_msis > 0.0 .and. CP_MLT(:,:,1:self%NLEV_IONDRAG) > 0.0)
               DTDT_IONDRAG = frictional_heating / (rho_msis * CP_MLT(:,:,1:self%NLEV_IONDRAG))
            elsewhere
               DTDT_IONDRAG = 0.0
            end where
         end if
         if (associated(NE_IONDRAG))       NE_IONDRAG       = ne_m3
         if (associated(RHO_MSIS_IONDRAG)) RHO_MSIS_IONDRAG = rho_msis

         deallocate(LATS_DEG, LONS_DEG, alt_km, ui, vi, ne_m3, species_fraction, &
                    rho_msis, drag_u, drag_v, frictional_heating)

         RETURN_(ESMF_SUCCESS)

      end subroutine IonDrag_Driver

   end subroutine RUN

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine compute_msis_density(IM, JM, NLEV_IONDRAG, IYEAR, DOY, UT_HOUR, &
                                    LATS_2D, LONS_2D, alt_km_in, rho_out, RC)
      ! Computes neutral mass density (kg/m^3) over the top NLEV_IONDRAG
      ! levels from MSIS number densities (O, N2, O2), following the same
      ! molecular-mass-weighted approach as calc_gas_specific_heat_MLT.
      ! Units of O_out/N2_out/O2_out from msis_point have been confirmed
      ! (cm^-3), so the CM3_TO_M3 conversion below is correct as written.
      integer, intent(in) :: IM, JM, NLEV_IONDRAG
      integer, intent(in) :: IYEAR, DOY
      real,    intent(in) :: UT_HOUR
      real, pointer, dimension(:,:), intent(in) :: LATS_2D, LONS_2D
      real, intent(in)  :: alt_km_in(:,:,:)   ! (IM, JM, NLEV_IONDRAG), per column
      real, intent(out) :: rho_out(:,:,:)
      integer, optional, intent(OUT) :: RC

      character(len=ESMF_MAXSTR) :: IAm
      integer :: STATUS

      real, parameter :: AMU_KG = 1.66053906660e-27
      real, parameter :: CM3_TO_M3 = 1.0e6   ! cm^-3 -> m^-3
      real, parameter :: MASS_O  = 16.0
      real, parameter :: MASS_N2 = 28.0
      real, parameter :: MASS_O2 = 32.0

      real(4) :: O_out, N2_out, O2_out, T_msis_out
      real(4) :: alt_r4, glat_r4, glong_r4, stl_r4
      integer :: i, j, k

      IAm = "compute_msis_density"

      do k = 1, NLEV_IONDRAG
         do j = 1, JM
            do i = 1, IM
               alt_r4  = real(alt_km_in(i,j,k), kind=4)
               glat_r4 = real(LATS_2D(i,j) * (180.0/MAPL_PI), kind=4)
               glong_r4 = real(LONS_2D(i,j) * (180.0/MAPL_PI), kind=4)
               ! Solar local time (hours) = UT hour + longitude/15
               stl_r4  = real(UT_HOUR + glong_r4/15.0, kind=4)

               call msis_point(IYEAR, DOY, nint(UT_HOUR*3600.0), &
                    alt_r4, glat_r4, glong_r4, stl_r4, &
                    O_out, N2_out, O2_out, T_msis_out)

               ! Number density (cm^-3) -> mass density (kg/m^3)
               rho_out(i,j,k) = ( real(O_out,kind=8)*MASS_O   + &
                                   real(N2_out,kind=8)*MASS_N2 + &
                                   real(O2_out,kind=8)*MASS_O2 ) &
                                 * AMU_KG * CM3_TO_M3
            end do
         end do
      end do

      RETURN_(ESMF_SUCCESS)
   end subroutine compute_msis_density

end module GEOS_IonDragGridCompMod
