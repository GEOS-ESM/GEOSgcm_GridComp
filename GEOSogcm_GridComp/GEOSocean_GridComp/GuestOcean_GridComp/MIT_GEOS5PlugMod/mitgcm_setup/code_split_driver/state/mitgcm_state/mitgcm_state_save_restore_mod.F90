      MODULE MITGCM_STATE_SAVE_RESTORE_MOD

      USE DYNVARS_H_STATE_MOD , ONLY :  &
          restore_state,                &
          save_state

      USE TIMEVARS_STATE_MOD , ONLY :   &
          restore_state,                &
          save_state

      USE STACKVARS_STATE_MOD , ONLY :  &
          MITGCM_STACKVARS,             &
          restore_state,                &
          save_state

      USE MITGCM_STATE_TYPES_MOD

      IMPLICIT NONE

      INTERFACE SAVE_STATE
       MODULE PROCEDURE MITGCM_ISTATE_SAVE    
      END INTERFACE
      INTERFACE RESTORE_STATE
       MODULE PROCEDURE MITGCM_ISTATE_RESTORE
      END INTERFACE

      CONTAINS

     
      SUBROUTINE MITGCM_ISTATE_SAVE( iState, sv )
!     -- Save internal state from global state --
      TYPE(MITGCM_ISTATE), POINTER :: iState
      TYPE(MITGCM_STACKVARS)       :: sv

      CALL SAVE_STATE( iState%dynvars_h     )
      CALL SAVE_STATE( iState%timevars      )
      CALL SAVE_STATE( iState%stackvars, sv )

      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_ISTATE_RESTORE( iState, sv )
!     -- Restore internal state to global state --
      TYPE(MITGCM_ISTATE), POINTER :: iState
      TYPE(MITGCM_STACKVARS)       :: sv

      CALL RESTORE_STATE( iState%dynvars_h)
      CALL RESTORE_STATE( iState%timevars )
      CALL RESTORE_STATE( iState%stackvars, sv )

      RETURN
      END SUBROUTINE

      END MODULE
