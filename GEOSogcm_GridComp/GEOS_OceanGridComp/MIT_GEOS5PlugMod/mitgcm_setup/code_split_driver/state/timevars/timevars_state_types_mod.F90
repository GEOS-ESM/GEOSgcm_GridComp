      MODULE TIMEVARS_STATE_TYPES_MOD
#include "CPP_OPTIONS.h"

!     Dynamically allocated type for the TIMEVARS information in MITgcm
      IMPLICIT NONE

      TYPE MITGCM_TIMEVARS
       SEQUENCE
       REAL*8          :: currentTime
       INTEGER         :: currentIter
       INTEGER         :: nIter0
       _RL             :: startTime
       INTEGER         :: nTimeSteps
       INTEGER         :: iLoop
       CHARACTER*4     :: pad1
       REAL*8          :: deltaT
       REAL*8          :: deltaTMom
       REAL*8, POINTER :: dTTracerLev(:) => NULL()
       REAL*8          :: deltaTClock 
       REAL*8          :: deltaTFreeSurf
      END TYPE

      END MODULE
