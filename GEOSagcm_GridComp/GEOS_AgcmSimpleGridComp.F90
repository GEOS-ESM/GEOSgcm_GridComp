#include "MAPL.h"

module GEOS_AgcmSimpleGridCompMod

   !=============================================================================
   !BOP
   ! \renewcommand{\comp}{\tt GEOS\_AgcmSimpleGridCompMod}
   !
   ! !MODULE: GEOS_AgcmSimpleGridCompMod
   !
   ! !DESCRIPTION: This gridded component (GC) combines the the GC that
   !   implements the Finite-Volume (FV) dynamics, with a simple physics
   !   component that implements the Held-Suarez benchmark forcing for
   !   testing dry dynamical cores.

   !USES:
   use ESMF
   use MAPL, only: MAPL_GridCompSetEntryPoint
   use MAPL, only: MAPL_GridCompAddSpec, user_setservices, MAPL_GridCompReexport
   use MAPL, only: MAPL_GridCompAddChild, MAPL_GridCompAddConnection
   use MAPL, only: user_setservices, MAPL_StateGetPointer
   use MAPL, only: MAPL_GridCompGet, MAPL_GridCompGetResource, MAPL_GridCompGetChildName
   use MAPL, only: MAPL_GridCompGetInternalState, MAPL_GridCompRunChild
   use MAPL, only: VERTICAL_STAGGER_CENTER, VERTICAL_STAGGER_EDGE, VERTICAL_STAGGER_NONE
   use MAPL, only: MAPL_RESTART_SKIP, MAPL_STATEITEM_SERVICE
   use MAPL, only: MAPL_Verify, MAPL_Return
   ! use GEOS_TopoGetMod, only: GEOS_TopoGet
   use GEOS_superdynGridCompMod, only: SDYN_SetServices => SetServices
   use GEOS_hsGridCompMod, only: PHYS_SetServices => SetServices

   implicit none
   private

   !PUBLIC MEMBER FUNCTIONS:
   public SetServices
   !EOP

contains

   !BOP
   !IROUTINE: SetServices

   !DESCRIPTION:  This is the only method of this component, since all of its
   !   registered methods can be defaulted to their MAPL generic versions.
   !   SetServices merely creates the children through MAPL and connects
   !   their Import-Export states. FV has some Imports that are
   !   not used in the Held-Suarez benchmark. Since these have proper defaults
   !   in FV, they can simply be terminated here.

   !INTERFACE:
   subroutine SetServices(gc, rc)
      !ARGUMENTS:

      type(ESMF_GridComp) :: gc ! gridded component
      integer, intent(out) :: rc ! return code
      !EOP

      integer :: status

      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, _RC)

      ! NOTE: The tracer advection is provided as a service by DYN, which we originally subscribed
      ! to in this gridcomp. Subsequently, the subscription was moved to the gridcomp DataMoist,
      ! using the configurable gridcomp and data-moist.yaml config file. The original code
      ! implementing the subscriber in this component is commented out and retained, for reference
      !
      ! #include "AgcmSimple_Internal___.h"
      ! Erstwhile AgcmSimple_StateSpecs.rc, for reference
      ! schema_version: 2.0.0
      ! component: AgcmSimple
      !
      ! # FAKE internals for bundling tracers for FV3
      ! category: INTERNAL
      ! #-------------------------------------------------------------------
      ! SHORT_NAME | UNITS   | FILLVAL | DIMS | VLOC | RESTART | LONG_NAME
      ! #-------------------------------------------------------------------
      ! Q          | kg kg-1 | 1.0e-6  | xyz  | C    | SKIP   | specific_humidity
      ! QLLS       | kg kg-1 |         | xyz  | C    | SKIP   |
      ! mass_fraction_of_large_scale_cloud_liquid_water
      ! QLCN       | kg kg-1 |         | xyz  | C    | SKIP   |
      ! mass_fraction_of_convective_cloud_liquid_water
      ! QILS       | kg kg-1 |         | xyz  | C    | SKIP   |
      ! mass_fraction_of_large_scale_cloud_ice_water
      ! QICN       | kg kg-1 |         | xyz  | C    | SKIP   |
      ! mass_fraction_of_convective_cloud_ice_water
      ! CLLS       | 1       |         | xyz  | C    | SKIP   | large_scale_cloud_area_fraction
      ! CLCN       | 1       |         | xyz  | C    | SKIP   | convective_cloud_area_fraction
      ! QRAIN      | kg kg-1 | 0.0     | xyz  | C    | SKIP   | mass_fraction_of_rain
      ! QSNOW      | kg kg-1 | 0.0     | xyz  | C    | SKIP   | mass_fraction_of_snow
      ! QGRAUPEL   | kg kg-1 | 0.0     | xyz  | C    | SKIP   | mass_fraction_of_graupel

      ! "FAKE" specs to provide PHIS and VARFLT to FV3
      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_INTERNAL, &
           short_name="PHIS", &
           standard_name="surface_geopotential_height", &
           units="m+2 s-2", &
           dims="xy", &
           vstagger=VERTICAL_STAGGER_NONE, &
           add_to_export=.true., _RC)
      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_INTERNAL, &
           short_name="VARFLT", &
           standard_name="variance_of_filtered_topography", &
           units="m+2", &
           dims="xy", &
           vstagger=VERTICAL_STAGGER_NONE, &
           add_to_export=.true., _RC)

      ! SUBSCRIBED "FAKE" spec to bundle tracers for FV3
      ! tracer_list = [ &
      !      "Q       ", &
      !      "QLLS    ", "QLCN    ", "QILS    ", "QICN    ", &
      !      "CLLS    ", "CLCN    ", "QRAIN   ", "QSNOW   ", "QGRAUPEL"]
      ! do iter = 1, size(tracer_list)
      !    call service_items%push_back(trim(tracer_list(iter)))
      ! end do
      ! call MAPL_GridCompAddSpec(gc, &
      !      state_intent=ESMF_STATEINTENT_IMPORT, &
      !      short_name="TRADV", &
      !      standard_name='advected_quantities', &
      !      units="unknown", &
      !      dims="xyz", & ! TODO: we shouldn't need dims/vstagger for bundles
      !      vstagger=VERTICAL_STAGGER_NONE, &
      !      itemtype=MAPL_STATEITEM_SERVICE, &
      !      service_items=service_items, _RC)

      call MAPL_GridCompAddChild(gc, "SDYN", user_setservices(SDYN_SetServices), "superdyn.yaml", _RC)
      call MAPL_GridCompAddChild(gc, "PHYS", user_setservices(PHYS_SetServices), "held-suarez.yaml", _RC)

      ! TODO: pchakrab - we don't really need these, do we?
      ! call MAPL_GridCompReexport(gc, src_comp="SDYN", src_name="T", _RC)
      ! call MAPL_GridCompReexport(gc, src_comp="SDYN", src_name="PS", _RC)

      ! Register connections between children
      call MAPL_GridCompAddConnection(gc, &
           src_comp="PHYS", &
           dst_comp="SDYN", &
           src_names="DUDT, DVDT, DTDT", _RC)
      call MAPL_GridCompAddConnection(gc, &
           src_comp="SDYN", &
           dst_comp="PHYS", &
           src_names="U, V, T, PLE", &
           dst_names="U, V, TEMP, PLE", _RC)
      call MAPL_GridCompAddConnection(gc, &
           src_comp="<self>", &
           dst_comp="SDYN", &
           src_names="PHIS, VARFLT", _RC)

      ! CONNECTION between provider and subscriber of the "FAKE" tracer bundle service
      ! call MAPL_GridCompAddConnection(gc, &
      !      src_comp="SDYN", &
      !      dst_comp="<self>", &
      !      src_names="TRADV", _RC)

      _RETURN(_SUCCESS)
   end subroutine SetServices

   subroutine Initialize(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp) :: gc ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_Clock) :: clock ! The clock
      integer, intent(out) :: rc ! Error code

      !DESCRIPTION: The Initialize method of this Gridded Component.

      type(ESMF_State) :: internal
      real, pointer, dimension(:, :) :: phis
      integer :: status

      ! PHIS ... (zeroed out, instead of reading from a zero file)
      call MAPL_GridCompGetInternalState(gc, internal, _RC)
      call MAPL_StateGetPointer(internal, phis, "PHIS", _RC)
      phis = 0.0

      _RETURN(_SUCCESS)
   end subroutine Initialize

   subroutine Run(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp) :: gc ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_Clock) :: clock ! The clock
      integer, intent(out) :: rc ! Error code
      !EOP

      ! character(len=:), allocatable :: child_name
      ! integer :: iter, num_children
      integer :: status

      ! call MAPL_GridCompGet(gc, num_children=num_children, _RC)
      ! do iter = 1, num_children
      !    child_name = MAPL_GridCompGetChildName(gc, iter, _RC)
      !    _HERE, "Child ", iter, ": ", trim(child_name)
      ! end do

      call MAPL_GridCompRunChild(gc, "SDYN", phase_name="run", _RC)
      call MAPL_GridCompRunChild(gc, "PHYS", phase_name="run", _RC)
      ! call MAPL_GridCompRunChild(gc, "DataMoist", phase_name="run", _RC)
      call MAPL_GridCompRunChild(gc, "SDYN", phase_name="run_add_incs", _RC)

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(import)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)
   end subroutine Run

end module GEOS_AgcmSimpleGridCompMod

subroutine SetServices(gc, rc)
   use ESMF
   use GEOS_AgcmSimpleGridCompMod, only : mySetservices => SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine SetServices
