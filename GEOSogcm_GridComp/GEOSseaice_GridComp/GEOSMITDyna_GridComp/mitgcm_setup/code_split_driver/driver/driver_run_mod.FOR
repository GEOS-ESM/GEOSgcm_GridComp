      MODULE DRIVER_RUN_MOD
C $Header: /u/gcmpack/MITgcm_contrib/PRM/multi_comp_setup/comp_mitgcm/code_basic/driver_run.F,v 1.1.1.1 2006/10/10 18:17:28 cnh Exp $

#include "CPP_OPTIONS.h"

      CONTAINS

CBOP
C     !ROUTINE: DRIVER_RUN
C     !INTERFACE:
      SUBROUTINE DRIVER_RUN( mitgcmIState, nTimeS)
C     !DESCRIPTION: \bv
C     *==================================================================
C     | SUBROUTINE driver_run
C     | o External driver control routine for MITgcm forward step
C     |   execution phase. 
C     *==================================================================
C     |
C     | DRIVER routines are used to control the MITgcm code from an external
C     | driver. This routine invokes the forward phase of code execution.
C     | The driver here is written for an ocean configuration and is designed
C     | for use with either an argument based call/return interface or with a
C     | messaging based "event loop" interface.
C     |
C     *==================================================================
C     \ev

C     !USES:
      USE MITGCM_STATE_MOD , ONLY :
     &    MITGCM_ISTATE,
     &    MITGCM_STACKVARS,
     &    RESTORE_STATE,
     &    SAVE_STATE
      USE EXPORT_STATE_MOD , ONLY :
     &     FILL
      USE IMPORT_STATE_MOD , ONLY :
     &     FILL

      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"

C     == Routine arguments ==
      INTEGER nTimeS
      TYPE(MITGCM_ISTATE), POINTER :: mitgcmIState

C     == Local variables ==
      INTEGER iStep
      TYPE(MITGCM_STACKVARS) sv
      INTEGER myThid
      INTEGER myCurrentIter
      _RL     myCurrentTime
      INTEGER iLoop
CEOP

      myThid = 1

      CALL RESTORE_STATE( mitgcmIState, sv )
      myCurrentTime = sv%myTime
      myCurrentIter = sv%myIter
      iLoop         = sv%iLoop

C     Fill the import fields
      CALL FILL( mitgcmIState%import )

      DO iStep = 1, nTimeS
       CALL FORWARD_STEP( iLoop, myCurrentTime, 
     &                    myCurrentIter, myThid )
       iLoop = iLoop+1
      ENDDO

      sv%myTime = myCurrentTime
      sv%myIter = myCurrentIter
      sv%iLoop  = iLoop
      CALL SAVE_STATE( mitgcmIState, sv )

C     Fill the export fields
      CALL FILL( mitgcmIState%export )

C     CALL COMP_OCN_POP_CONTEXT(1)
      CALL POPDIR

      RETURN
      END SUBROUTINE

      END MODULE
