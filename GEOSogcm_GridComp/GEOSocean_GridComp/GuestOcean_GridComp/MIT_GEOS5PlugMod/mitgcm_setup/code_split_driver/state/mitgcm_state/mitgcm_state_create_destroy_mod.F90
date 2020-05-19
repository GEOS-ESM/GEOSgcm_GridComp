      MODULE MITGCM_STATE_CREATE_DESTROY_MOD
!--   Module of code that handles creation and destruction of of MITgcm state type variables --
!--   The MITgcm state variables contain child types that are created (destroyed) here by
!--   calls to those child types create (detroy) functions.

!--   DYNVARS state child functions
      USE DYNVARS_H_STATE_MOD , ONLY :  &
          create,                       &
          destroy

!--   TIMEVARS state child functions
      USE TIMEVARS_STATE_MOD , ONLY :    &
          create,                        &
          destroy

!--   STACKVARS state child functions
      USE STACKVARS_STATE_MOD , ONLY :   &
          create,                        &
          destroy

!--   FFIELDS state child functions
      USE FFIELDS_H_STATE_MOD , ONLY :   &
          create,                        &
          destroy

!--   Export state child functions
      USE EXPORT_STATE_MOD , ONLY :      &
          create,                        &
          destroy

!--   Import state child functions
      USE IMPORT_STATE_MOD , ONLY :      &
          create,                        &
          destroy

      USE MITGCM_STATE_TYPES_MOD

      IMPLICIT NONE

      INTERFACE CREATE
       MODULE PROCEDURE MITGCM_MULTI_ISTATE_CREATE
       MODULE PROCEDURE MITGCM_ISTATE_CREATE
      END INTERFACE
      INTERFACE DESTROY
       MODULE PROCEDURE MITGCM_MULTI_ISTATE_DESTROY
       MODULE PROCEDURE MITGCM_ISTATE_DESTROY
      END INTERFACE

      CONTAINS

      SUBROUTINE MITGCM_MULTI_ISTATE_CREATE( statesArr, nStates,  &
                 snx, sny, olx, oly, n3ds, nsx, nsy               &
                 )
!     -- Allocate memory for an array of MITgcm states of a specific size --
      TYPE(MITGCM_MULTI_ISTATE), POINTER :: statesArr
      INTEGER                            :: nStates
      INTEGER                            :: snx, sny, olx, oly
      INTEGER, intent(in), dimension(:)  :: n3ds
      INTEGER                            :: nsx, nsy

!     -- Local variables --
      INTEGER I
      TYPE(MITGCM_ISTATE), POINTER :: p

      ALLOCATE( statesArr                      )
      ALLOCATE( statesArr%istate_list(nStates) )
      DO I=1,nStates
       p => statesArr%istate_list(I)
       CALL CREATE( p,snx, sny, olx, oly, n3ds, nsx, nsy )
      ENDDO

      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_MULTI_ISTATE_DESTROY( statesArr )
!     -- Allocate memory for an array of MITgcm states of a specific size --
      TYPE(MITGCM_MULTI_ISTATE), POINTER :: statesArr

!     -- Local variables --
      INTEGER I
      TYPE(MITGCM_ISTATE), POINTER :: p

      IF ( .NOT. ASSOCIATED(statesArr) ) RETURN

      DO I=LBOUND(statesArr%istate_list,1),UBOUND(statesArr%istate_list,1)
       p => statesArr%istate_list(I)
       CALL DESTROY( p  )
      ENDDO
      DEALLOCATE( statesArr%istate_list )
      DEALLOCATE( statesArr )

      RETURN
      END SUBROUTINE
      SUBROUTINE MITGCM_ISTATE_CREATE( iState,         &
                 snx, sny, olx, oly, n3ds, nsx, nsy    &
                 )
!     -- Allocate memory an MITgcm internal state of a specific size --
      TYPE(MITGCM_ISTATE), POINTER :: iState
      INTEGER                      :: snx, sny, olx, oly
      INTEGER, intent(in), dimension(:)  :: n3ds
      INTEGER                      :: nsx, nsy

      ALLOCATE( iState )
      CALL CREATE( iState%dynvars_h,                                  &
       snx, sny, olx, oly, n3ds, nsx, nsy                             )
      CALL CREATE( iState%timevars, n3ds                              )
      CALL CREATE( iState%stackvars                                   )
      CALL CREATE( iState%export,                                     &
       snx, sny, olx, oly, n3ds, nsx, nsy                             )
      CALL CREATE( iState%import,                                     &
       snx, sny, olx, oly, n3ds, nsx, nsy                             )

      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_ISTATE_DESTROY( iState )
!     -- Free memory being used by MITgcm iState --
      TYPE(MITGCM_ISTATE), POINTER :: iState

      IF ( .NOT. ASSOCIATED(iState) ) RETURN

      CALL DESTROY( iState%dynvars_h )
      CALL DESTROY( iState%timevars  )
      CALL DESTROY( iState%stackvars )
      CALL DESTROY( iState%export    )
      CALL DESTROY( iState%import    )
      DEALLOCATE( iState )

      RETURN
      END SUBROUTINE

      END MODULE
