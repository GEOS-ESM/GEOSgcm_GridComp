#include "MAPL.h"

module GEOS_SuperdynGridCompMod

   !BOP
   !MODULE: GEOS_SuperdynGridCompMod -- A Module to combine Dynamics and Gravity-Wave-Drag Gridded Components
   !USES:
   use ESMF
   use MAPL, only: MAPL_GridCompSetEntryPoint
   use MAPL, only: MAPL_GridCompAddChild, user_setservices, MAPL_GridCompRunChild
   use MAPL, only: MAPL_GridCompGetResource
   use MAPL, only: MAPL_GridCompAddSpec, MAPL_GridCompReexport, MAPL_GridCompAddConnection
   use MAPL, only: VERTICAL_STAGGER_CENTER
   use MAPL, only: MAPL_Verify, MAPL_Return, MAPL_Assert

   ! use FVdycore_GridCompMod, only : FV_SetServices => SetServices
   use FVdycoreCubed_GridComp, only : FV3_SetServices => SetServices
   ! use ARIESg3_GridCompMod, only : ARIES_SetServices => SetServices
   ! use GEOS_DatmoDynGridCompMod, only : DATMO_SetServices => SetServices
   ! use AdvCore_GridCompMod, only : ADV_SetServices => SetServices

   implicit none
   private

   !PUBLIC MEMBER FUNCTIONS:
   public SetServices

   !DESCRIPTION: This gridded component (GC) combines the Dynamics GC and
   !   the Gravity Wave Drag GC into a new composite SuperDyn GC.
   !
   !\vspace{5mm}
   !   {\bf Import Couplings}:
   !
   !   The Import Couplings of the SuperDyn GC are the tendencies of the atmospheric
   !   state variables U,V,T,PE (due to external diabatic forcing) in addition to
   !   a collection of "Friendly" tracers for advection.  The Friendlies will be searched
   !   for moisture for use in virtual effects in both the Dynamics and Gravity Wave Drag
   !   parameterization.  If no moisture is found, the SuperDyn will be run dry.
   !
   !\begin{verbatim}
   !       DUDT .... U-Wind                    Tendency (m/s)
   !       DVDT .... V-Wind                    Tendency (m/s)
   !       DPEDT ... Edge-Pressure             Tendency (Pa/s)
   !       DTDT .... Mass-Weighted Temperature Tendency (Pa K/s)
   !       TRACER .. Friendly Tracers                   (unknown)
   !     If Non-Hydrostatic Dynamics
   !       DWDT .... W-Wind                    Tendency (m/s)
   !\end{verbatim}
   !
   !\vspace{5mm}
   !   {\bf Run Method}:
   !
   !   The run method first calls the Gravity Wave Drag parameterization.
   !   The tendencies of the atmospheric state variables created by the GWD are then ADDED to the
   !   SuperDyn Import Couplings (i.e., state variable tendencies due to external diabatic forcing),
   !   which are then used to force the Dynamics GC.
   !
   !\vspace{5mm}
   !   {\bf Export Couplings}:
   !
   !   The Export Couplings of the SuperDyn GC are the union of the Export
   !   Couplings of the individual GCs.  It should be noted that the SuperDyn GC
   !   controls the GEOS Topo Utility and produces Topo Variables based on the Grid
   !   defined by the DYN GC.  The Topo Variables are computed during the SuperDyn Initialize
   !   method, and are part of the SuperDyn Export Couplings.  ESMF utilities may be used to regrid
   !   these Topo Variables to other components with differing Grids.
   !EOP

   integer :: DYN
   integer :: ADV = -1

contains

   !BOP
   !IROUTINE: SetServices -- Sets ESMF services for this component
   !INTERFACE:
   subroutine SetServices(gc, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc ! gridded component
      integer, intent(out) :: rc ! return code

      !DESCRIPTION:  The SetServices for the SuperDyn GC needs to register its
      !   Initialize, Run, and Finalize methods.  In addition, we need to create the
      !   children GCs (DYN and GWD) and run their respective SetServices.
      !EOP

      character(len=:), allocatable :: dycore
      integer :: scm_sl, status

      ! Register services for this component
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, phase_name="Run", _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, RunAddIncs, phase_name="RunAddIncs", _RC)

      ! Create children's gridded components and invoke their SetServices
      call MAPL_GridCompGetResource(gc, "DYCORE", dycore, default="FV3", _RC)
      call MAPL_GridCompGetResource(gc, "SCM_SL", scm_sl, default=0, _RC)

      select case (trim(dycore))
      ! case ("FV")
      !    call MAPL_GridCompAddChild(gc, "DYN", user_setservices(FV_SetServices), "dyn.yaml", _RC)
      case ("FV3")
         call MAPL_GridCompAddChild(gc, "DYN", user_setservices(FV3_SetServices), "dyn.yaml", _RC)
      ! case ("FV3+ADV")
      !    call MAPL_GridCompAddChild(gc, "DYN", user_setservices(FV3_SetServices), "dyn.yaml", _RC)
      !    call MAPL_GridCompAddChild(gc, "ADV", user_setservices(ADV_SetServices), "dyn.yaml", _RC)
      ! case ("ARIES")
      !    call MAPL_GridCompAddChild(gc, "DYN", user_setservices(ARIES_SetServices), "dyn.yaml", _RC)
      ! case ("DATMO")
      !    call MAPL_GridCompAddChild(gc, "DYN", user_setservices(DATMO_SetServices), "dyn.yaml", _RC)
      case default
         _FAIL("Invalid DYCORE specified for SuperDyn: " // trim(dycore))
      end select

      ! Add Exports promoted from child (FV) exports
      if (SCM_SL /= 0) then
         ! print *,'SuperDyn: adding LHOBS and SHOBS exports'
         call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="LHOBS", _RC)
         call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="SHOBS", _RC)
      end if

      ! Re-export, use MAPL_GridCompReexport(gc, src_comp="DYN", src_name="U", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="U", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="V", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="W", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="T", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="S", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="TH", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PLE", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PL", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="ZLE", _RC)
#ifdef HAS_GIGATRAJ
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="ZL", _RC)
#endif
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PREF", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="AK", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="BK", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PLK", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PKE", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PS", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="DELP", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="US", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="VS", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="TA", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="QA", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="SPEED", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="DZ", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="TROPP_BLENDED", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="TROPK_BLENDED", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PV", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="TV", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="OMEGA", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="EPV", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PEANA", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="DTHVDTANAINT", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PEPHY", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="DTHVDTPHYINT", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="DQVDTANAINT", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="DQLDTANAINT", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="DQIDTANAINT", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="DOXDTANAINT", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="AREA", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="U_DGRID", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="V_DGRID", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PT", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PE", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="QV_DYN_IN", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="T_DYN_IN", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="U_DYN_IN", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="V_DYN_IN", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PLE_DYN_IN", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="DTDTDYN", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="DQVDTDYN", _RC)
      ! Service
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="TRADV", _RC)
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="PLE4", _RC)
      ! Needed for NCEP GWD
      call MAPL_GridCompReexport(gc, src_comp="DYN", src_name="DXC", _RC)

      ! Create DYN+ADV connectivities
      if (adjustl(DYCORE) == "FV3+ADV") then
         call MAPL_GridCompAddConnection(gc, &
              src_comp="DYN", &
              src_names="MFX, MFY, CX, CY, PLE0, PLE1", &
              dst_comp="ADV", &
              _RC)
      end if

      _RETURN(_SUCCESS)

   end subroutine SetServices

   !BOP
   !IROUTINE: Initialize -- Initialize method for the composite SuperDyn Gridded Component
   !INTERFACE:
   subroutine Initialize(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION: The Initialize method of the SuperDyn Composite Gridded Component first
      !   calls the Initialize method of the child Dynamics.  The Dynamics Initialize method will
      !   create the ESMF GRID, which will then be used to set the GRID associated with the
      !   SuperDyn Composite Component itself.  It should be noted that the
      !   SuperDyn Initialize method also invokes the GEOS Topo Utility which creates all
      !   topography related quantities.
      !EOP

#ifdef PRINT_STATES
      character(len=*), parameter :: IAm = "GEOS_SuperdynGridCompMod:Initialize"
      call WRITE_PARALLEL(trim(IAm) // ": IMPORT State")
      if (MAPL_am_I_root()) call ESMF_StatePrint(import, _RC)
      call WRITE_PARALLEL(trim(IAm) // ": EXPORT State")
      if (MAPL_am_I_root()) call ESMF_StatePrint(export, _RC)
#endif

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(gc)
      _UNUSED_DUMMY(import)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)
   end subroutine Initialize

   !BOP
   !IROUTINE: Run -- Run method for the composite SuperDyn Gridded Component
   !INTERFACE:
   subroutine Run(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION: The run method first calls the Gravity Wave Drag parameterization.
      !   The tendencies of the atmospheric state variables created by the GWD are then ADDED to the
      !   SuperDyn Import Couplings (i.e., state variable tendencies due to external diabatic forcing),
      !   which are then used to force the Dynamics GC.
      !EOP

      integer :: status

      call MAPL_GridCompRunChild(gc, "DYN", phase_name="Run", _RC)
      if (ADV /= -1) then
         call MAPL_GridCompRunChild(gc, "ADV", phase_name="Run", _RC)
      end if

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(import)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)
   end subroutine Run

   !IROUTINE: RunAddIncs -- Run method to add increments
   !INTERFACE:
   subroutine RunAddIncs(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION: The run method to add increments
      !EOP

      integer :: status

      call MAPL_GridCompRunChild(gc, "DYN", phase_name="RunAddIncs", _RC)

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(import)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)
   end subroutine RunAddIncs

end module GEOS_SuperdynGridCompMod
