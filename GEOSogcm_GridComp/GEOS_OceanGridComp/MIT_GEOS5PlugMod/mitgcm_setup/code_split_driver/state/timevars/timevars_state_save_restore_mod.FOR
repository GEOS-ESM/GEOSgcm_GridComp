      MODULE TIMEVARS_STATE_SAVE_RESTORE_MOD
#include "CPP_OPTIONS.h"

      USE TIMEVARS_STATE_TYPES_MOD
      IMPLICIT NONE

      INTERFACE SAVE_STATE
       MODULE PROCEDURE TIMEVARS_SAVE_STATE
      END INTERFACE
      INTERFACE RESTORE_STATE
       MODULE PROCEDURE TIMEVARS_RESTORE_STATE
      END INTERFACE

      CONTAINS

      SUBROUTINE TIMEVARS_SAVE_STATE( timeVarsPtr )          
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
!     -- Deallocate memory for MITgcm timevars variables. --
      TYPE(MITGCM_TIMEVARS), POINTER :: timeVarsPtr
 
!     -- Local variables --
      TYPE(MITGCM_TIMEVARS), POINTER :: p

      p => timeVarsPtr

      p%nIter0        = nIter0
      p%startTime     = startTime
      p%nTimeSteps    = nTimeSteps
      p%deltaT        = deltaT
      p%deltaTMom     = deltaTMom
      p%dTTracerLev   = dTTracerLev
      p%deltaTClock   = deltaTClock 

      RETURN
      END SUBROUTINE

      SUBROUTINE TIMEVARS_RESTORE_STATE( timeVarsPtr )          
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
!     -- Deallocate memory for MITgcm timevars variables. --
      TYPE(MITGCM_TIMEVARS), POINTER :: timeVarsPtr

!     -- Local variables --
      TYPE(MITGCM_TIMEVARS), POINTER :: p

      p => timeVarsPtr

      nIter0        = p%nIter0    
      startTime     = p%startTime 
      nTimeSteps    = p%nTimeSteps
      deltaT        = p%deltaT
      deltaTMom     = p%deltaTMom
      dTTracerLev   = p%dTTracerLev
      deltaTClock   = p%deltaTClock 

      RETURN
      END SUBROUTINE

      END MODULE
