C $Header: /u/gcmpack/MITgcm_contrib/PRM/multi_comp_setup/comp_mitgcm/code_basic/main.F,v 1.2 2006/12/13 19:56:59 cnh Exp $
C $Name:  $

CBOI
C
C
CEOI


C--   Get C preprocessor options
#include "CPP_OPTIONS.h"
#include "CPP_EEOPTIONS.h"

CBOP
C !ROUTINE: MAIN

C !INTERFACE:
      PROGRAM MAIN

      USE ESMF

      USE MITGCM_STATE_MOD , ONLY :
     & MITGCM_ISTATE_CONTAINER,
     & GETDP

      USE MITGCM_DRIVER_MOD , ONLY :
     & DRIVER_INIT,
     & DRIVER_RUN

      IMPLICIT NONE

#include "mpif.h"

C !DESCRIPTION:
C     *==========================================================*
C     | PROGRAM MAIN                                            
C     | o Test top-level for external driver to MITgcm kernel.
C     |   Supports storing and restoring of state from upper
C     |   level. Can be adapted to drive ensembles of configurations
C     |   with different parameters. It also illustrates some of the 
C     |   machinery used in multi-scale runs to drive core MITgcm
C     |   dynamics.
C     *==========================================================*

!     Variable to hold model state for each instance
      TYPE(MITGCM_ISTATE_CONTAINER) :: mitgcmIState(2)

!     Loop counters
      INTEGER nTimeSteps
      INTEGER I

!     MPI context variable
      INTEGER myComm

!     Variables for accessing MITgcm internal state
      REAL*8, POINTER :: uVel(:,:,:,:,:)
      REAL*8, POINTER :: eta( :,:,:,:  )
      REAL*8, POINTER :: deltaTClock
      REAL*8, POINTER :: deltaTMom   
      REAL*8, POINTER :: currentTime
      REAL*8, POINTER :: dTTracerLev(:)
!     Variables for temporary timestep reduction
      REAL*8          :: deltaTScal
      REAL*8          :: deltaTClockFrac
      REAL*8          :: deltaTMomFrac
      REAL*8, POINTER :: dTTracerLevFrac(:)
      INTEGER         :: iLo, iHi

!     Variables for booting ESMF
      type(ESMF_VM)       :: vm

!     ESMF clock setting variables
      TYPE(ESMF_Calendar) :: myCalendar
      TYPE(ESMF_TimeInterval) :: timeStep
      TYPE(ESMF_Time)         :: startTime
      TYPE(ESMF_Time)         :: stopTime
      TYPE(ESMF_Clock)        :: clock
      TYPE(ESMF_Time)         :: alarmTime
      TYPE(ESMF_Alarm)        :: alarm
!     ESMF clock reading variables
      TYPE(ESMF_Time)         :: theCurrentTime
      TYPE(ESMF_TimeInterval) :: theCurrentTimeStep
      REAL*8                  :: theSecs
      REAL*8                  :: theDt

!     Variable for error codes
      INTEGER RC

!     Setup ESMF (only do this once)
      call ESMF_Initialize (vm=vm, logKindFlag=ESMF_LOGKIND_NONE, rc=rc)
      CALL ESMF_VMGet(VM, mpiCommunicator=myComm, rc=RC)

!     Create and set the clock and an alarm
      call ESMF_CalendarSetDefault(ESMF_CALKIND_GREGORIAN, RC=RC)
      CALL ESMF_TimeIntervalSet( timeStep, s_r8=3600.0d0, rc=RC )
!     CALL ESMF_TimeIntervalSet( timeStep, d=1      , rc=RC )
      CALL ESMF_TimeSet( startTime, yy=1999, mm=12, dd=31, h=23, rc=RC )
      CALL ESMF_TimeSet(  stopTime, yy=2000, mm= 1, dd= 1, h= 2, rc=RC )
      clock=ESMF_ClockCreate(name="Top clock", timeStep=timeStep, 
     &     startTime=startTime, stopTime=stopTime, rc=RC )
      CALL ESMF_TimeSet( alarmTime, yy=2000, mm=1, dd=1 )
      alarm = ESMF_AlarmCreate( clock=clock, ringTime=alarmTime, rc=RC )

      CALL DRIVER_INIT( mitgcmIState(1)%p, myComm)
      CALL DRIVER_INIT( mitgcmIState(2)%p, myComm)

      nTimeSteps = 1
      DO I=1,10
!      - Get the current time and how long to run for.
       CALL ESMF_ClockGet( clock, 
     &                     timeStep=theCurrentTimeStep, 
     &                     currTime=theCurrentTime,
     &                     rc=rc)
       CALL ESMF_TimeGet( theCurrentTime, s_r8=theSecs )
       CALL ESMF_TimeIntervalGet( theCurrentTimeStep, s_r8=theDt )
       WRITE(0,*) ' theSecs = ', theSecs
       WRITE(0,*) ' theDt   = ', theDt
       CALL DRIVER_RUN( mitgcmIState(1)%p, nTimeSteps )
!      Get data pointers to internal fields.
       CALL GETDP( mitgcmIState(1)%p, 'UVEL',         uVel,        rc )
       CALL GETDP( mitgcmIState(1)%p, 'ETAN',         eta,         rc )
       CALL GETDP( mitgcmIState(1)%p, 'DELTATCLOCK',    deltaTClock, rc )
       CALL GETDP( mitgcmIState(1)%p, 'DELTATMOM',      deltaTMom,   rc )
       CALL GETDP( mitgcmIState(1)%p, 'DTTRACERLEV',    dTTracerLev, rc )
       CALL GETDP( mitgcmIState(1)%p, 'CURRENTTIME',  currentTime, rc )

       WRITE (0,*) 'MITgcm currentTime = ', currentTime
       WRITE (0,*) 'MITgcm deltaTClock = ', deltaTClock
       WRITE (0,*) 'Whole time steps = ', INT(theDt/deltaTClock)
       WRITE (0,*) 'Fractional residual time step = ', theDt/deltaTClock - INT(theDt/deltaTClock) 

       deltaTScal = theDt/deltaTClock - INT(theDt/deltaTClock)
       deltaTClockFrac = deltaTClock*deltaTScal
       deltaTMomFrac   = deltaTMom*deltaTScal
       WRITE (0,*) ' Fractional deltaTClock = ', deltaTClockFrac
       WRITE (0,*) ' Fractional deltaTMom   = ', deltaTMomFrac
       iLo=LBOUND(dtTracerLev,1)
       iHi=UBOUND(dtTracerLev,1)
       ALLOCATE(   dtTracerLevFrac(iLo:iHi) )
       dtTracerLevFrac=dTTracerLev*deltaTScal
       WRITE (0,*) ' Fractional dtTracerLev   = ', dtTracerLevFrac
       DEALLOCATE( dtTracerLevFrac          )
       
       CALL ESMF_ClockAdvance(clock, rc=RC)
      ENDDO
      nTimeSteps = 1
      DO I=1,10
       CALL DRIVER_RUN( mitgcmIState(2)%p, nTimeSteps )
      ENDDO

      CALL ESMF_Finalize(rc=rc)


      END

!     -- Test that we can link in the MAPL library --
      SUBROUTINE FOO

      use ESMF
      use MAPL_Mod

      type(ESMF_GridComp) :: GC  ! gridded component
      INTEGER status
      

      CALL MAPL_AddImportSpec(GC,                            
     &   SHORT_NAME         = 'TAUY',                             
     &   LONG_NAME          = 'Agrid_northward_stress_on_skin',   
     &   UNITS              = 'N m-2',                            
     &   DIMS               = MAPL_DimsHorzOnly,                  
     &   VLOCATION          = MAPL_VLocationNone,                 
     &   RC=STATUS  )

      RETURN
      END
