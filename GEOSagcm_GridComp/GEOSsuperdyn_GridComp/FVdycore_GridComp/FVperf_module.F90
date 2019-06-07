!===============================================================================
! CVS: $Id$
! CVS: $Source$
! CVS: $Name$
!===============================================================================

!-----------------------------------------------------------------------
!              ESMA - Earth System Modeling Applications
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: FVperf_module --- Simple interfaces for performance profiling
!
! !INTERFACE:

MODULE FVperf_module

! !USES:
   use dynamics_vars, only : T_FVDYCORE_GRID
#if defined(MAPL_MODE)
   use MAPL_Mod            ! MAPL base class
#elif defined(CAM_MODE)
   use perf_mod
#else
   use perf_mod
#endif
#if defined( SPMD )
   use mod_comm, only: mp_barrier
#endif

CONTAINS

! !DESCRIPTION: A hack to toggle between GEOS5 and CAM profiling
!
!  The basic problem solved here is to access GENSTATE in GEOS\_MODE
!  without being overly intrusive (e.g. putting GEOS\_MODE in every
!  file in which MAPL_TimerOn is used.  If GEOS\_MODE
!  is defined, the GENSTATE must be initialized outside this file,
!  the various timing markers registered with GenericStateClockAdd,
!  and the genstate in this module manual set.  If CAM\_MODE is 
!  defined, the user may use FVstartclock/FVstopclock exactly like
!  the CAM utilities t\_startf and t\_stopf.
!
!  This module will be removed as soon as there is consensus on a
!  unified profiling utility.
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: FVstartclock --- start the clock
!
! !INTERFACE:
      subroutine FVstartclock(grid,marker)
! !USES:
      implicit none
! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(inout) :: grid
      character(LEN=*)      , intent(in) :: marker
!EOP
!-----------------------------------------------------------------------
!BOC
#if defined(MAPL_MODE)
      call MAPL_TimerOn(grid%FVgenstate,marker)
#elif defined(CAM_MODE)
      call t_startf(marker)    
#else
      call t_startf(marker)
#endif
!EOC
      end subroutine FVstartclock
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: FVstopclock --- stop the clock
!
! !INTERFACE:
      subroutine FVstopclock(grid,marker)
! !USES:
      implicit none
! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(inout) :: grid
      character(LEN=*)      , intent(in) :: marker
!EOP
!-----------------------------------------------------------------------
!BOC
#if defined(MAPL_MODE)
      call MAPL_TimerOff(grid%FVgenstate,marker)
#elif defined(CAM_MODE)
      call t_stopf(marker)    
#else
      call t_stopf(marker)
#endif
!EOC
      end subroutine FVstopclock
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: FVbarrierclock --- instrumented timing barrier
!
! !INTERFACE:
      subroutine FVbarrierclock(grid,marker)
! !USES:
      implicit none
! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(inout) :: grid
      character(LEN=*)      , intent(in) :: marker
!EOP
!-----------------------------------------------------------------------
!BOC
#if defined(MAPL_MODE)
#if ( defined SPMD ) && ( defined TIMING_BARRIERS )
      call MAPL_TimerOn(grid%FVgenstate,marker)
      call mp_barrier()
      call MAPL_TimerOff(grid%FVgenstate,marker)
#endif
#elif defined(CAM_MODE)
#if ( defined SPMD )
      if (t_profile_onf()) then
         if (t_barrier_onf()) then
            call t_startf(marker)    
            call mp_barrier()
            call t_stopf(marker)    
         endif
      endif
#endif
#else
#if ( defined SPMD )
      if (t_profile_onf()) then
         if (t_barrier_onf()) then
            call t_startf(marker)    
            call mp_barrier()
            call t_stopf(marker)    
         endif
      endif
#endif
#endif
!EOC
      end subroutine FVbarrierclock
!-----------------------------------------------------------------------

END MODULE FVperf_module


