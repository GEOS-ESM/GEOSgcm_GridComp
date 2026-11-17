#include "MAPL_Generic.h"

module LIS_GridComp

   !BOP
   !MODULE: LIS_GridComp -- couples MAPL/GEOS directly to LISF.

   !DESCRIPTION:
   ! MAPL/ESMF Gridded Component wrapping the Land Information System
   ! (LISF). Unlike LISF's own NUOPC cap
   ! (@LISF/lis/runmodes/nuopc_cpl_mode/LIS_NUOPC_Cap.F90), this Plug talks to
   ! LISF's core Fortran API directly -- the same routines LISF's own offline
   ! driver (@LISF/lis/offline/lisdrv.F90) and its "retrospective" runmode
   ! (@LISF/lis/runmodes/retrospective/retrospective_runMod.F90) use -- and
   ! hand-maps GEOS Import/Export fields into/out of LIS's tile space using
   ! LIS's own LIS_grid2tile/LIS_tile2grid utilities. No ESMF/NUOPC coupling
   ! layer is involved.
   !
   ! This component assumes a single LIS nest (n=1) running on the same
   ! decomposition/PET layout as its parent GEOS grid.

   !USES:
   use ESMF
   use MAPL

   use LIS_coreMod, only: LIS_rc, LIS_config_init, LIS_core_init, LIS_finalize
   use LIS_domainMod, only: LIS_domain_init
   use LIS_paramsMod, only: LIS_param_init
   use LIS_surfaceModelMod, only: LIS_surfaceModel_init, LIS_surfaceModel_setup
   use LIS_surfaceModelMod, only: LIS_surfaceModel_readrestart, LIS_surfaceModel_writerestart
   use LIS_surfaceModelMod, only: LIS_surfaceModel_f2t, LIS_surfaceModel_run
   use LIS_metforcingMod, only: LIS_metforcing_init, LIS_FORC_State
   use LIS_FORC_AttributesMod, only: LIS_FORC_Tair, LIS_FORC_Qair
   use LIS_FORC_AttributesMod, only: LIS_FORC_SWdown, LIS_FORC_LWdown
   use LIS_FORC_AttributesMod, only: LIS_FORC_Wind_E, LIS_FORC_Wind_N
   use LIS_FORC_AttributesMod, only: LIS_FORC_Psurf, LIS_FORC_Rainf, LIS_FORC_Snowf
   use LIS_historyMod, only: LIS_tile2grid, LIS_grid2tile
   use LIS_histDataMod, only: LIS_histData, LIS_metadataEntry
   use LIS_histDataMod, only: LIS_MOC_QH, LIS_MOC_QLE, LIS_MOC_QS, LIS_MOC_SWE
   use LIS_timeMgrMod, only: LIS_timemgr_set

   use iso_fortran_env, only: real32

   implicit none
   private

   !PUBLIC MEMBER FUNCTIONS:
   public :: SetServices

   ! LISF nest index this Plug drives. GEOS coupling only needs one nest.
   integer, parameter :: NEST = 1

contains

   !BOP
   !IROUTINE: SetServices -- Sets ESMF services for this component
   !INTERFACE:

   subroutine SetServices(gc, rc)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! gridded component
      integer, optional :: rc ! return code

      !DESCRIPTION: Registers Initialize/Run/Finalize/Record and the
      !   Import/Export state specs for the direct-call LISF Plug.
      !EOP

      character(len=ESMF_MAXSTR) :: IAm
      integer :: status
      character(len=ESMF_MAXSTR) :: comp_name

      IAm = 'SetServices'
      call ESMF_GridCompGet(gc, name=comp_name, _RC)
      IAm = trim(comp_name) // IAm

      ! Set the Initialize, Run, Finalize entry points
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_FINALIZE, Finalize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_WRITERESTART, Record, _RC)

      ! Set the state variable specs (generated from LIS_StateSpecs.rc)
#include "LIS_Import___.h"
#include "LIS_Export___.h"
#include "LIS_Internal___.h"

      ! Set the Profiling timers
      call MAPL_TimerAdd(gc, name="INITIALIZE", _RC)
      call MAPL_TimerAdd(gc, name="RUN", _RC)
      call MAPL_TimerAdd(gc, name="FINALIZE", _RC)

      ! Generic SetServices
      call MAPL_GenericSetServices(gc, _RC)

      RETURN_(ESMF_SUCCESS)
   end subroutine SetServices

   !BOP
   !IROUTINE: Initialize -- Initialize method for the direct-call LISF Plug
   !INTERFACE:
   subroutine Initialize(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, optional, intent(out) :: rc ! Error code:

      !EOP

      character(len=ESMF_MAXSTR) :: IAm
      integer :: status
      character(len=ESMF_MAXSTR) :: comp_name

      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_VM) :: vm
      type(ESMF_Grid) :: grid
      type(ESMF_Time) :: current_time
      integer :: comm
      integer :: year, month, day, hour, minute, second
      integer :: local_counts(3)
      character(len=ESMF_MAXSTR) :: lis_config_file

      IAm = "Initialize"
      call ESMF_GridCompGet(gc, name=comp_name, _RC)
      IAm = trim(comp_name) // IAm

      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      call MAPL_TimerOn(MAPL, "INITIALIZE")

      call MAPL_GenericInitialize(gc, import, export, clock, _RC)

      ! LISF configuration file (lis.config-style resource file)
      call MAPL_GetResource(MAPL, lis_config_file, Label="LISF_CONFIG_FILE:", DEFAULT="lis.config", _RC)

      ! Reuse GEOS's VM/communicator so LIS runs on the same PET layout as its parent
      call ESMF_VMGetCurrent(vm, _RC)
      call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)

      ! Reuse grid and decomposition from parent GEOS GridComp
      call ESMF_GridCompGet(gc, grid=grid, _RC)
      call MAPL_GridGet(grid, localCellCountPerDim=local_counts, _RC)

      ! Initialize LISF: configuration, domain, parameters, surface model(s),
      ! and met-forcing structures. This mirrors LIS_init_retrospective in
      ! @LISF/lis/runmodes/retrospective/retrospective_runMod.F90, minus the
      ! pieces (DA, routing, RTM, irrigation, app model) not yet wired into
      ! this Plug. We call these directly instead of the dynamic
      ! lisinit(trim(LIS_rc%runmode)//char(0)) dispatch (available after
      ! LIS_config_init registers the runmode plugins -- see
      ! LIS_runmode_pluginMod.F90) because that would unconditionally pull in
      ! the excluded subsystems above and follow whatever RUNMODE: is set in
      ! the config file rather than this Plug's fixed minimal sequence.
      LIS_rc%lis_config_file = trim(lis_config_file)
      call LIS_config_init(cmd_args=.false., vm=vm, comm=comm)

      call LIS_domain_init()
      call LIS_param_init()
      call LIS_surfaceModel_init()
      call LIS_surfaceModel_setup()
      call LIS_metforcing_init()

      ASSERT_(LIS_rc%lnc(NEST) == local_counts(1))
      ASSERT_(LIS_rc%lnr(NEST) == local_counts(2))

      if (trim(LIS_rc%startcode) == "restart") then
         call LIS_surfaceModel_readrestart()
      end if

      ! Sync LIS's internal clock to the incoming ESMF clock.
      call ESMF_ClockGet(clock, currTime=current_time, _RC)
      call ESMF_TimeGet(current_time, YY=year, MM=month, DD=day, H=hour, M=minute, S=second, _RC)
      call LIS_timemgr_set(LIS_rc, year, month, day, hour, minute, second, 0, LIS_rc%ts)

      call MAPL_TimerOff(MAPL, "INITIALIZE")

      RETURN_(ESMF_SUCCESS)
   end subroutine Initialize

   !BOP
   !IROUTINE: Run -- Run method for the direct-call LISF Plug
   !INTERFACE:
   subroutine Run(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc
      type(ESMF_State), intent(inout) :: import
      type(ESMF_State), intent(inout) :: export
      type(ESMF_Clock), intent(inout) :: clock
      integer, optional, intent(out) :: rc
      !EOP

      character(len=ESMF_MAXSTR) :: IAm
      integer :: status
      character(len=ESMF_MAXSTR) :: comp_name

      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_Time) :: current_time
      integer :: year, month, day, hour, minute, second

      real(kind=real32), pointer, dimension(:, :) :: TAIR_F, QAIR_F, SWDOWN_F, LWDOWN_F
      real(kind=real32), pointer, dimension(:, :) :: EWIND_F, NWIND_F, PSURF_F, RAINF_F, SNOWF_F
      real(kind=real32), pointer, dimension(:, :) :: AVGSURFT, QH, QLE, QS, SWE

      IAm = "Run"
      call ESMF_GridCompGet(gc, name=comp_name, _RC)
      IAm = trim(comp_name) // IAm

      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      call MAPL_TimerOn(MAPL, "RUN")

      ! Get pointers to the GEOS 2D (lon,lat) Import/Export arrays
      call MAPL_GetPointer(import, TAIR_F, 'TAIR_F', _RC)
      call MAPL_GetPointer(import, QAIR_F, 'QAIR_F', _RC)
      call MAPL_GetPointer(import, SWDOWN_F, 'SWDOWN_F', _RC)
      call MAPL_GetPointer(import, LWDOWN_F, 'LWDOWN_F', _RC)
      call MAPL_GetPointer(import, EWIND_F, 'EWIND_F', _RC)
      call MAPL_GetPointer(import, NWIND_F, 'NWIND_F', _RC)
      call MAPL_GetPointer(import, PSURF_F, 'PSURF_F', _RC)
      call MAPL_GetPointer(import, RAINF_F, 'RAINF_F', _RC)
      call MAPL_GetPointer(import, SNOWF_F, 'SNOWF_F', _RC)

      call MAPL_GetPointer(export, AVGSURFT, 'AVGSURFT', ALLOC=.true., _RC)
      call MAPL_GetPointer(export, QH, 'QH', _RC)
      call MAPL_GetPointer(export, QLE, 'QLE', _RC)
      call MAPL_GetPointer(export, QS, 'QS', _RC)
      call MAPL_GetPointer(export, SWE, 'SWE', _RC)

      ! Sync LIS's clock to the current GEOS time, then scatter GEOS's gridded
      ! forcing into LIS's tile space.
      call ESMF_ClockGet(clock, currTime=current_time, _RC)
      call ESMF_TimeGet(current_time, YY=year, MM=month, DD=day, H=hour, M=minute, S=second, _RC)
      call LIS_timemgr_set(LIS_rc, year, month, day, hour, minute, second, 0, LIS_rc%ts)

      call PutForcing(trim(LIS_FORC_Tair%varname(1)), TAIR_F, _RC)
      call PutForcing(trim(LIS_FORC_Qair%varname(1)), QAIR_F, _RC)
      call PutForcing(trim(LIS_FORC_SWdown%varname(1)), SWDOWN_F, _RC)
      call PutForcing(trim(LIS_FORC_LWdown%varname(1)), LWDOWN_F, _RC)
      call PutForcing(trim(LIS_FORC_Wind_E%varname(1)), EWIND_F, _RC)
      call PutForcing(trim(LIS_FORC_Wind_N%varname(1)), NWIND_F, _RC)
      call PutForcing(trim(LIS_FORC_Psurf%varname(1)), PSURF_F, _RC)
      call PutForcing(trim(LIS_FORC_Rainf%varname(1)), RAINF_F, _RC)
      call PutForcing(trim(LIS_FORC_Snowf%varname(1)), SNOWF_F, _RC)

      ! Advance LIS one GEOS heartbeat: forcing-to-tile-space interpolation,
      ! then run the surface model(s) for this nest.
      call LIS_surfaceModel_f2t(NEST)
      call LIS_surfaceModel_run(NEST)

      if (associated(QH)) then
         call GetDiagnostic(LIS_MOC_QH, 'Qh', QH, _RC)
      end if
      if (associated(QLE)) then
         call GetDiagnostic(LIS_MOC_QLE, 'Qle', QLE, _RC)
      end if
      if (associated(QS)) then
         call GetDiagnostic(LIS_MOC_QS, 'Qs', QS, _RC)
      end if
      if (associated(SWE)) then
         call GetDiagnostic(LIS_MOC_SWE, 'SWE', SWE, _RC)
      end if

      ! TODO: fill AVGSURFT and the remaining exports the same way via GetDiagnostic
      call MAPL_TimerOff(MAPL, "RUN")

      RETURN_(ESMF_SUCCESS)
   end subroutine Run

   subroutine PutForcing(varname, gridvar, rc)
      character(len=*), intent(in) :: varname
      real(kind=real32), pointer, intent(in) :: gridvar(:, :)
      integer, optional, intent(out) :: rc

      integer :: status
      type(ESMF_Field) :: field
      real, pointer :: tilevar(:)

      call ESMF_StateGet(LIS_FORC_State(NEST), trim(varname), field, _RC)
      call ESMF_FieldGet(field, localDE=0, farrayPtr=tilevar, _RC)
      call LIS_grid2tile(NEST, real(gridvar), tilevar)

      _RETURN(_SUCCESS)
   end subroutine PutForcing

   subroutine GetDiagnostic(moc_index, moc_name, gridvar, rc)
      integer, intent(in) :: moc_index
      character(len=*), intent(in) :: moc_name
      real(kind=real32), pointer, intent(in) :: gridvar(:, :)
      integer, optional, intent(out) :: rc

      type(LIS_metadataEntry), pointer :: diag
      integer :: slot
      logical :: instantaneous

      _ASSERT(moc_index > 0, moc_name // ' not enabled in MODEL_OUTPUT_LIST')
      diag => LIS_histData(NEST)%ptr_into_lsm_list(moc_index)%dataEntryPtr
      _ASSERT(diag%selectOpt /= 0, moc_name // ' not enabled in MODEL_OUTPUT_LIST')
      ! For timeAvgOpt 1/3, modelOutput holds a running sum, not the current value
      instantaneous = diag%timeAvgOpt == 0 .or. diag%timeAvgOpt == 2
      _ASSERT(instantaneous, moc_name // ' needs time-avg option 0 or 2')
      slot = merge(2, 1, diag%timeAvgOpt == 2)
      call LIS_tile2grid(NEST, gridvar, diag%modelOutput(slot, :, 1))

      _RETURN(_SUCCESS)
   end subroutine GetDiagnostic

   !BOP
   !IROUTINE: Finalize -- Finalize method for the direct-call LISF Plug
   !INTERFACE:
   subroutine Finalize(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc
      type(ESMF_State), intent(inout) :: import
      type(ESMF_State), intent(inout) :: export
      type(ESMF_Clock), intent(inout) :: clock
      integer, optional, intent(out) :: rc
      !EOP

      character(len=ESMF_MAXSTR) :: IAm
      integer :: status
      character(len=ESMF_MAXSTR) :: comp_name
      type(MAPL_MetaComp), pointer :: MAPL

      IAm = "Finalize"
      call ESMF_GridCompGet(gc, name=comp_name, _RC)
      IAm = trim(comp_name) // IAm

      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      call MAPL_TimerOn(MAPL, "FINALIZE")
      call LIS_finalize()
      call MAPL_TimerOff(MAPL, "FINALIZE")

      call MAPL_GenericFinalize(gc, import, export, clock, _RC)

      RETURN_(ESMF_SUCCESS)
   end subroutine Finalize

   !BOP
   !IROUTINE: Record -- write LISF restarts
   !INTERFACE:
   subroutine Record(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc
      type(ESMF_State), intent(inout) :: import
      type(ESMF_State), intent(inout) :: export
      type(ESMF_Clock), intent(inout) :: clock
      integer, optional, intent(out) :: rc
      !EOP

      character(len=ESMF_MAXSTR) :: IAm
      integer :: status
      character(len=ESMF_MAXSTR) :: comp_name

      IAm = "Record"
      call ESMF_GridCompGet(gc, name=comp_name, _RC)
      IAm = trim(comp_name) // IAm

      call LIS_surfaceModel_writerestart(NEST)

      RETURN_(ESMF_SUCCESS)
   end subroutine Record

end module LIS_GridComp
