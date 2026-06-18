! $Id: GEOS_AgcmSimpleGridComp.F90,v 1.1 2007/05/16 15:33:09 trayanov Exp $

#include "MAPL_Generic.h"

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
   !
   !USES:

   use ESMF
   use MAPL, only: MAPL_GridCompSetEntryPoint
   use MAPL, only: MAPL_AddImportSpec, MAPL_AddExportSpec, MAPL_AddInternalSpec
   use MAPL, only: MAPL_DimsHorzVert, MAPL_VLocationCenter
   use MAPL, only: MAPL_GenericSetServices, MAPL_GenericInitialize, MAPL_GenericRunCouplers
   use MAPL, only: MAPL_AddChild, MAPL_AddConnectivity, MAPL_TerminateImport
   use MAPL, only: MAPL_MetaComp
   use MAPL, only: MAPL_GetObjectFromGC, MAPL_Get, MAPL_GridCompGetFriendlies, MAPL_GetResource
   use MAPL, only: MAPL_GridCreate, MAPL_TimerOn, MAPL_TimerOff
   use MAPL, only: MAPL_Verify, MAPL_Return
   use GEOS_TopoGetMod, only: GEOS_TopoGet

   use GEOS_superdynGridCompMod, only: SDYN_SetServices => SetServices
   use GEOS_hsGridCompMod, only: PHS_SetServices => SetServices

   implicit none
   private

   integer :: SDYN
   integer :: PHS

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
      type(ESMF_GridComp), intent(inout) :: gc ! gridded component
      integer, optional, intent(out) :: rc ! return code
      !EOP

      !BOC

      ! ErrLog Variables
      character(len=ESMF_MAXSTR) :: IAm
      integer :: status
      character(len=ESMF_MAXSTR) :: comp_name

      ! Get my name and set-up traceback handle
      IAm = 'SetServices'
      call ESMF_GridCompGet(gc, NAME=comp_name, _RC)
      IAm = trim(comp_name) // trim(IAm)

      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, _RC)

      ! dummy import for testing comcurrent ens
      !    call MAPL_AddImportSpec ( gc, &
      !         SHORT_NAME='DTDT', &
      !         LONG_NAME='temperature increment', &
      !         UNITS='K s-1', &
      !         DIMS=MAPL_DimsHorzVert, &
      !         VLOCATION=MAPL_VLocationCenter, _RC)

      ! These internal specs are "fake" and here only to provide moisture to FV
      ! this component will not touch them
      call MAPL_AddInternalSpec(gc, &
           SHORT_NAME='Q', &
           LONG_NAME='specific_humidity', &
           UNITS='kg kg-1', &
           FRIENDLYTO='DYNAMICS', &
           default=1.0e-6, &
           ! RESTART=MAPL_RestartRequired, &
           DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, _RC)

      call MAPL_AddInternalSpec(gc, &
           SHORT_NAME='QLLS', &
           LONG_NAME='mass_fraction_of_large_scale_cloud_liquid_water', &
           UNITS='kg kg-1', &
           FRIENDLYTO='DYNAMICS', &
           DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, _RC)

      call MAPL_AddInternalSpec(gc, &
           SHORT_NAME='QLCN', &
           LONG_NAME='mass_fraction_of_convective_cloud_liquid_water', &
           UNITS='kg kg-1', &
           FRIENDLYTO='DYNAMICS', &
           DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, _RC)

      call MAPL_AddInternalSpec(gc, &
           SHORT_NAME='QILS', &
           LONG_NAME='mass_fraction_of_large_scale_cloud_ice_water', &
           UNITS='kg kg-1', &
           FRIENDLYTO='DYNAMICS', &
           DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, _RC)

      call MAPL_AddInternalSpec(gc, &
           SHORT_NAME='QICN', &
           LONG_NAME='mass_fraction_of_convective_cloud_ice_water', &
           UNITS='kg kg-1', &
           FRIENDLYTO='DYNAMICS', &
           DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, _RC)

      call MAPL_AddInternalSpec(gc, &
           SHORT_NAME='CLLS', &
           LONG_NAME='large_scale_cloud_area_fraction', &
           UNITS='1', &
           FRIENDLYTO='DYNAMICS', &
           DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, _RC)

      call MAPL_AddInternalSpec(gc, &
           SHORT_NAME='CLCN', &
           LONG_NAME='convective_cloud_area_fraction', &
           UNITS='1', &
           FRIENDLYTO='DYNAMICS', &
           DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, _RC)

      call MAPL_AddInternalSpec(gc, &
           SHORT_NAME='QRAIN', &
           LONG_NAME='mass_fraction_of_rain', &
           UNITS='kg kg-1', &
           FRIENDLYTO='DYNAMICS', &
           default=0.0, &
           DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, _RC)

      call MAPL_AddInternalSpec(gc, &
           SHORT_NAME='QSNOW', &
           LONG_NAME='mass_fraction_of_snow', &
           UNITS='kg kg-1', &
           FRIENDLYTO='DYNAMICS', &
           default=0.0, &
           DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, _RC)

      call MAPL_AddInternalSpec(gc, &
           SHORT_NAME='QGRAUPEL', &
           LONG_NAME='mass_fraction_of_graupel', &
           UNITS='kg kg-1', &
           FRIENDLYTO='DYNAMICS', &
           default=0.0, &
           DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, _RC)

      ! Register children with MAPL and go down their SS hierarchy
      SDYN = MAPL_AddChild(gc, NAME='SUPERDYNAMICS', SS=SDYN_SetServices, _RC)
      PHS = MAPL_AddChild(gc, NAME='HSPHYSICS', SS=PHS_SetServices, _RC)

      call MAPL_AddExportSpec(gc, &
           SHORT_NAME='T', &
           CHILD_ID=SDYN, &
           _RC)

      call MAPL_AddExportSpec(gc, &
           SHORT_NAME='PS', &
           CHILD_ID=SDYN, &
           _RC)

      ! Register connections between children
      call MAPL_AddConnectivity(gc, &
           SHORT_NAME=(/ 'DUDT', 'DVDT', 'DTDT' /), &
           SRC_ID=PHS, &
           DST_ID=SDYN, &
           _RC)

      call MAPL_AddConnectivity(gc, &
           SRC_NAME=(/ 'U    ', 'V    ', 'T    ', 'PLE  ' /), &
           DST_NAME=(/ 'U    ', 'V    ', 'TEMP ', 'PLE  ' /), &
           SRC_ID=SDYN, &
           DST_ID=PHS, &
           _RC)

      ! SetServices clean-up on the way back up through the hierarchy
      call MAPL_TerminateImport(gc, SHORT_NAME=(/'PHIS ', 'DPEDT'/), CHILD=SDYN, _RC)

      call MAPL_GenericSetServices(gc, _RC)

      _RETURN(_SUCCESS)
   end subroutine SetServices

   !EOC

   subroutine Initialize(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, optional, intent(out) :: rc ! Error code

      !DESCRIPTION: The Initialize method of this Gridded Component.

      ! ErrLog Variables
      character(len=ESMF_MAXSTR) :: IAm
      integer :: status
      character(len=ESMF_MAXSTR) :: comp_name

      ! local vars
      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_State), pointer :: GIM(:)
      type(ESMF_Field) :: FIELD
      type(ESMF_FieldBundle) :: BUNDLE
      type(ESMF_Config) :: cf
      type(ESMF_Alarm) :: replay_shutoff_alarm
      type(ESMF_TimeInterval) :: shutoff
      integer :: rplshut

      IAm = "Initialize"
      call ESMF_GridCompGet(gc, NAME=comp_name, config=cf, _RC)
      IAm = trim(comp_name) // IAm

      call MAPL_GridCreate(gc, _RC)

      call MAPL_GenericInitialize(gc, import, export, clock, _RC)

      ! Get my MAPL_Generic state
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      ! Fill Childrens TOPO variables and Diagnostics
      call MAPL_Get(MAPL, GIM=GIM, _RC)

      ! PHIS ...
      call ESMF_StateGet(GIM(SDYN), 'PHIS', FIELD, _RC)
      call GEOS_TopoGet(cf, MEAN=FIELD, _RC)

      ! TRADV ...
      call ESMF_StateGet(GIM(SDYN), 'TRADV', BUNDLE, _RC)
      call MAPL_GridCompGetFriendlies(gc, "DYNAMICS", BUNDLE, _RC)

      ! Initialize alarms
      call MAPL_GetResource(MAPL, rplshut, Label="REPLAY_SHUTOFF:", default=-3600, _RC)
      call ESMF_TimeIntervalSet(shutoff, S=abs(rplshut), _RC)
      replay_shutoff_alarm = ESMF_AlarmCreate( &
           NAME="ReplayShutOff", &
           clock=clock, &
           ringInterval=shutoff, &
           sticky=.true., &
           _RC)

      _RETURN(_SUCCESS)
   end subroutine Initialize

   subroutine Run(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, optional, intent(out) :: rc ! Error code
      !EOP

      ! ErrLog Variables
      character(len=ESMF_MAXSTR) :: IAm
      integer :: status
      character(len=ESMF_MAXSTR) :: comp_name

      ! Local derived type aliases
      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_GridComp), pointer :: GCS(:)
      type(ESMF_State), pointer :: GIM(:)
      type(ESMF_State), pointer :: GEX(:)

      ! Get the target components name and set-up traceback handle.
      IAm = "Run"
      call ESMF_GridCompGet(gc, NAME=comp_name, _RC)
      IAm = trim(comp_name) // trim(IAm)

      ! Get my MAPL_MetaComp
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      ! Start the TOTAL timer
      call MAPL_TimerOn(MAPL, "TOTAL")
      call MAPL_TimerOn(MAPL, "RUN")

      ! Get esmf internal state from generic state.
      call MAPL_Get(MAPL, GCS=GCS, GIM=GIM, GEX=GEX, _RC)

      call ESMF_GridCompRun(GCS(SDYN), importState=GIM(SDYN), exportState=GEX(SDYN), clock=clock, PHASE=1, userRC=status)
      _VERIFY(status)

      call MAPL_GenericRunCouplers(MAPL, CHILD=SDYN, clock=clock, _RC)

      call ESMF_GridCompRun(GCS(PHS), importState=GIM(PHS), exportState=GEX(PHS), clock=clock, userRC=status)
      _VERIFY(status)

      call ESMF_GridCompRun(GCS(SDYN), importState=GIM(SDYN), exportState=GEX(SDYN), clock=clock, PHASE=2, userRC=status)
      _VERIFY(status)

      call MAPL_TimerOff(MAPL, "RUN")
      call MAPL_TimerOff(MAPL, "TOTAL")

      _RETURN(_SUCCESS)
   end subroutine Run

end module GEOS_AgcmSimpleGridCompMod
