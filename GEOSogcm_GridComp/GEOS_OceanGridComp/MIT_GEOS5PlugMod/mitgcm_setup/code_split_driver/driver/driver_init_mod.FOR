       MODULE DRIVER_INIT_MOD
C $Header: /u/gcmpack/MITgcm_contrib/PRM/multi_comp_setup/comp_mitgcm/code_basic/driver_init.F,v 1.1.1.1 2006/10/10 18:17:28 cnh Exp $
C $Name:  $

#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: DRIVER_INIT
C     !INTERFACE:

      CONTAINS

      SUBROUTINE DRIVER_INIT(  mitgcmIState, myComm)
!     SUBROUTINE DRIVER_INIT(  mitgcmIState, myComm,
!    &                         iLoop,
!    &                         myCurrentTime,
!    &                         myCurrentIter )
C     !DESCRIPTION: \bv
C     *==================================================================
C     | SUBROUTINE driver_init
C     | o External driver control routine for MITgcm init phase.
C     *==================================================================
C     | DRIVER routines are used to control the MITgcm code from an external
C     | driver. This routine invokes the init phase of code execution
C     *==================================================================
C     \ev
C     !USES:
      USE MITGCM_STATE_MOD , ONLY :
     &    MITGCM_ISTATE,
     &    MITGCM_STACKVARS,
     &    CREATE,
     &    DESTROY,
     &    SAVE_STATE
      USE EXPORT_STATE_MOD , ONLY :
     &     FILL

      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
#endif

C     == Routine arguments ==
      TYPE(MITGCM_ISTATE), POINTER :: mitgcmIState
      INTEGER myComm
CEOP

C     == Local variables ====
      INTEGER n3dSize(4)
      INTEGER iLoop
      INTEGER myThid
      INTEGER myCurrentIter
      _RL     myCurrentTime
      TYPE(MITGCM_STACKVARS) sv

C     CALL COMP_OCN_PUSH_CONTEXT(1)
!     CALL SETDIR(myRank, iNumber)

!     ALLOCATE( mitgcmIstate )
!     STOP

!     Pass down a communicator from ESMF for EEBOOT
      CALL EEBOOT( myComm )
      myThid = 1
      n3dSize(1:4) = 1
      n3dSize(1) = Nr
#ifdef ALLOW_SEAICE
      n3dSize(2) = nITD
      n3dSize(3) = nIceLayers
      n3dSize(4) = nSnowLayers
#endif

C     Extract the import public and private fields, filling internal
C     state as appropriate.

      CALL INITIALISE_FIXED(    myThid )
      CALL INITIALISE_VARIA(    myThid )
      myCurrentTime = startTime
      myCurrentIter = nIter0
      iLoop  = 1

C     Allocate a variable to save the internal state
C     and then fill it from the current global heap and stack state.
C     o create
      CALL CREATE( mitgcmIState,
     &      sNx, sNy, OLx, OLy, n3dSize, nSx, nSy  )
C     o save
      sv%myTime = myCurrentTime
      sv%myIter = myCurrentIter
      sv%iLoop  = iLoop
      CALL SAVE_STATE( mitgcmIState, sv )
      CALL FILL( mitgcmIState%export )

!     CLOSE(errorMessageUnit)
!     CLOSE(standardMessageUnit)
C     CALL COMP_OCN_POP_CONTEXT(1)
!     CALL POPDIR

      RETURN
      END SUBROUTINE
      END MODULE
