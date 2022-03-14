      MODULE STACKVARS_STATE_MOD
#include "CPP_OPTIONS.h"

!     -- Dynamically allocated type for stack information in MITgcm --
      IMPLICIT NONE

      TYPE MITGCM_STACKVARS
       SEQUENCE
       _RL     myTime
       INTEGER iloop
       INTEGER myThid
       INTEGER myIter
       CHARACTER*4     :: pad1
      END TYPE

      INTERFACE CREATE
       MODULE PROCEDURE MITGCM_STACKVARS_CREATE
      END INTERFACE
      INTERFACE DESTROY
       MODULE PROCEDURE MITGCM_STACKVARS_DESTROY
      END INTERFACE
      INTERFACE SAVE_STATE
       MODULE PROCEDURE MITGCM_STACKVARS_SAVE_STATE
      END INTERFACE
      INTERFACE RESTORE_STATE
       MODULE PROCEDURE MITGCM_STACKVARS_RESTORE_STATE
      END INTERFACE

      CONTAINS

      SUBROUTINE MITGCM_STACKVARS_CREATE( varsPtr )
!     -- Allocate memory for MITgcm stackvars variables. --
      TYPE(MITGCM_STACKVARS), POINTER :: varsPtr

      ALLOCATE( varsPtr )
     
      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_STACKVARS_DESTROY( varsPtr )          
!     -- Deallocate memory for MITgcm stackvars variables. --
      TYPE(MITGCM_STACKVARS), POINTER :: varsPtr

      IF ( ASSOCIATED(varsPtr) ) THEN
       DEALLOCATE( varsPtr )
      ENDIF
 
      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_STACKVARS_SAVE_STATE( varsPtr, varsVal )
!     -- Save values for MITgcm stackvars variables. --
      TYPE(MITGCM_STACKVARS), POINTER :: varsPtr
      TYPE(MITGCM_STACKVARS)          :: varsVal
 
!     -- Local variables --
      TYPE(MITGCM_STACKVARS), POINTER :: p

      p => varsPtr

      p%myTime     = varsVal%myTime
      p%iLoop      = varsVal%iLoop
      p%myThid     = varsVal%myThid
      p%myIter     = varsVal%myIter

      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_STACKVARS_RESTORE_STATE( varsPtr, varsVal )
!     -- Restore values for MITgcm stackvars variables. --
      TYPE(MITGCM_STACKVARS), POINTER :: varsPtr
      TYPE(MITGCM_STACKVARS)          :: varsVal

!     -- Local variables --
      TYPE(MITGCM_STACKVARS), POINTER :: p

      p => varsPtr

      varsVal%myTime  = p%myTime
      varsVal%iLoop   = p%iLoop
      varsVal%myThid  = p%myThid
      varsVal%myIter  = p%myIter

      RETURN
      END SUBROUTINE

      END MODULE
