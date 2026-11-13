#include "MAPL.h"

! Root of the turbulence-shocmf-sa/run2 regression case (test only, not part of
! the model). TURBULENCE's Run2 continues from what its Run1 left in the
! internal state (RESTART: SKIP fields such as DKS, SINC and the tendencies)
! and in its exports (SFLX, QFLX, SDFL, QDFL, ...), which MAPL3 cannot
! restart, and Surface changes some of its imports in between. So this
! component runs Run1, then sets the fields listed under phase2_fields to
! their values before MAPL2's Run2, then runs Run2.
!
! Its states, children and connections come from its configuration, as for
! the configurable gridcomp. For each name X in phase2_fields it imports X,
! connected to the same provider as TURBULENCE's import X (so it is the same
! field), and X_PHASE2, holding the value captured before MAPL2's Run2.

module TurbulencePhasesDriver_GridCompMod

   use ESMF
   use MAPL

   implicit none
   private

   public :: SetServices

contains

   subroutine SetServices(gc, rc)
      type(ESMF_GridComp) :: gc
      integer, intent(out) :: rc

      integer :: status

      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, phase_name='run', _RC)

      _RETURN(_SUCCESS)
   end subroutine SetServices

   recursive subroutine Run(gc, import, export, clock, rc)
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      integer :: status, i, rank
      type(ESMF_HConfig) :: hconfig
      character(len=ESMF_MAXSTR), allocatable :: names(:)
      type(ESMF_Field) :: field
      real, pointer :: x2(:,:), y2(:,:), x3(:,:,:), y3(:,:,:)

      call MAPL_GridCompGet(gc, hconfig=hconfig, _RC)
      names = ESMF_HConfigAsStringSeq(hconfig, keyString='phase2_fields', stringLen=ESMF_MAXSTR, _RC)

      call MAPL_GridCompRunChild(gc, child_name='TURBULENCE', phase_name='Run1', _RC)

      do i = 1, size(names)
         call ESMF_StateGet(import, trim(names(i)), field, _RC)
         call ESMF_FieldGet(field, rank=rank, _RC)
         select case (rank)
         case (2)
            call MAPL_StateGetPointer(import, x2, trim(names(i)), _RC)
            call MAPL_StateGetPointer(import, y2, trim(names(i))//'_PHASE2', _RC)
            x2 = y2
         case (3)
            call MAPL_StateGetPointer(import, x3, trim(names(i)), _RC)
            call MAPL_StateGetPointer(import, y3, trim(names(i))//'_PHASE2', _RC)
            x3 = y3
         case default
            _FAIL('phase2_fields: rank 2 or 3 only: '//trim(names(i)))
         end select
      end do

      call MAPL_GridCompRunChild(gc, child_name='TURBULENCE', phase_name='Run2', _RC)

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)
   end subroutine Run

end module TurbulencePhasesDriver_GridCompMod

subroutine TurbulencePhasesDriver_SetServices(gc, rc)
   use ESMF
   use TurbulencePhasesDriver_GridCompMod, only: mySetservices => SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetservices(gc, rc=rc)
end subroutine TurbulencePhasesDriver_SetServices
