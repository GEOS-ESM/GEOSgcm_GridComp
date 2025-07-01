module perf_mod

!-----------------------------------------------------------------------
!
! Purpose: This module is responsible for controlling the performance
!          timer logic.
!
! Author:  P. Worley, January 2007
!
! $Id$
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- Uses ----------------------------------------------------------------
!-----------------------------------------------------------------------

   use shr_sys_mod,       only: shr_sys_abort
   use shr_kind_mod,      only: SHR_KIND_CS, SHR_KIND_CM, SHR_KIND_CX, &
                                SHR_KIND_R8, SHR_KIND_I8

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make the default access private
   save

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public t_startf
   public t_stopf

!=======================================================================
contains
!=======================================================================
!========================================================================
!
   subroutine t_startf(event, handle)
!-----------------------------------------------------------------------
! Purpose: Start an event timer
! Author: P. Worley
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
!
   ! performance timer event name
   character(len=*), intent(in) :: event
!
!---------------------------Input/Output arguments----------------------
!
   ! GPTL event handle
   integer,  optional :: handle
!
!---------------------------Local workspace-----------------------------
!
   integer  ierr                          ! GPTL error return
   integer  str_length, i                 ! support for adding
                                          !  detail prefix
   character(len=2) cdetail               ! char variable for detail
!
!-----------------------------------------------------------------------
!

   return
   end subroutine t_startf
!
!========================================================================
!
   subroutine t_stopf(event, handle)
!-----------------------------------------------------------------------
! Purpose: Stop an event timer
! Author: P. Worley
!-----------------------------------------------------------------------
!---------------------------Input arguments-----------------------------
!
   ! performance timer event name
   character(len=*), intent(in) :: event
!
!---------------------------Input/Output arguments----------------------
!
   ! GPTL event handle
   integer, optional :: handle
!
!---------------------------Local workspace-----------------------------
!
   integer  ierr                          ! GPTL error return
   integer  str_length, i                 ! support for adding
                                          !  detail prefix
   character(len=2) cdetail               ! char variable for detail
!
!-----------------------------------------------------------------------
!

   return
   end subroutine t_stopf
!
end module perf_mod
