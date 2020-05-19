#define _HARD_ASSERT( _a ) IF ( _a .NE. 0 ) THEN; WRITE(0,*) code; STOP 'ABNORMAL END'; ENDIF
      MODULE MITGCM_STATE_GETDP_MOD
!--   Routine to get pointer to data for a given entry in an MITgcm state variable --

      USE MITGCM_STATE_TYPES_MOD

      IMPLICIT NONE

      INTERFACE GETDP
       MODULE PROCEDURE MITGCM_STATE_GET_R85D
       MODULE PROCEDURE MITGCM_STATE_GET_R84D
       MODULE PROCEDURE MITGCM_STATE_GET_R81D
       MODULE PROCEDURE MITGCM_STATE_GET_R8
       MODULE PROCEDURE MITGCM_STATE_GET_R45D
       MODULE PROCEDURE MITGCM_STATE_GET_R44D
       MODULE PROCEDURE MITGCM_STATE_GET_R41D
       MODULE PROCEDURE MITGCM_STATE_GET_R4
      END INTERFACE

      CONTAINS

      SUBROUTINE MITGCM_STATE_GET_R85D( miState, code, dPtr, rc )
!--   Get pointer to data for real*8 5-d array
      TYPE(MITGCM_ISTATE), POINTER :: miState
      CHARACTER*(*)                :: code
      REAL*8, POINTER              :: dPtr(:,:,:,:,:)
      INTEGER                      :: rc
    
      RC = -1
      IF ( code .EQ. 'UVEL' ) THEN
       IF ( ASSOCIATED(miState%dynvars_h%dynvars%uVel) ) THEN
        dPtr => miState%dynvars_h%dynvars%uVel
        RC   =  0
       ENDIF
      ENDIF

      _HARD_ASSERT( RC )
      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_STATE_GET_R84D( miState, code, dPtr, rc )
!--   Get pointer to data for real*8 4-d array
      TYPE(MITGCM_ISTATE), POINTER :: miState
      CHARACTER*(*)                :: code
      REAL*8, POINTER              :: dPtr(:,:,:,:)
      INTEGER                      :: rc

      RC = -1
      IF ( code .EQ. 'ETAN' ) THEN
       IF ( ASSOCIATED(miState%dynvars_h%dynvars%etaN) ) THEN
        dPtr => miState%dynvars_h%dynvars%etaN
        RC   =  0
       ENDIF
      ENDIF

      _HARD_ASSERT( RC )
      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_STATE_GET_R81D(   miState, code, dPtr, rc )
!--   Get pointer to data for real*8 1-d array
      TYPE(MITGCM_ISTATE), POINTER :: miState
      CHARACTER*(*)                :: code
      REAL*8, POINTER              :: dPtr(:)
      INTEGER                      :: rc

      RC = -1

      IF ( code .EQ. 'DTTRACERLEV' ) THEN
       IF ( ASSOCIATED(miState%timeVars ) ) THEN
        dPtr => miState%timeVars%dTTracerLev
        RC = 0
       ENDIF
      ENDIF

      _HARD_ASSERT( RC )
      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_STATE_GET_R8(   miState, code, dPtr, rc )
!--   Get pointer to data for real*8 scalar
      TYPE(MITGCM_ISTATE), POINTER :: miState
      CHARACTER*(*)                :: code
      REAL*8, POINTER              :: dPtr
      INTEGER                      :: rc

      RC = -1

      IF ( code .EQ. 'CURRENTTIME' ) THEN
       IF ( ASSOCIATED(miState%stackVars ) ) THEN
        dPtr => miState%stackVars%myTime
        RC = 0
       ENDIF
      ENDIF

      IF ( code .EQ. 'DELTATCLOCK' ) THEN
       IF ( ASSOCIATED(miState%timeVars ) ) THEN
        dPtr => miState%timeVars%deltaTClock
        RC = 0
       ENDIF
      ENDIF

      IF ( code .EQ. 'DELTATMOM' ) THEN
       IF ( ASSOCIATED(miState%timeVars ) ) THEN
        dPtr => miState%timeVars%deltaTMom
        RC = 0
       ENDIF
      ENDIF

      _HARD_ASSERT( RC )
      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_STATE_GET_R45D( miState, code, dPtr, rc )
!--   Get pointer to data for real*4 5-d array
      TYPE(MITGCM_ISTATE), POINTER :: miState
      CHARACTER*(*)                :: code
      REAL*4, POINTER              :: dPtr(:,:,:,:,:)
      INTEGER                      :: rc

      RC = -1

      _HARD_ASSERT( RC )
      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_STATE_GET_R44D( miState, code, dPtr, rc )
!--   Get pointer to data for real*4 4-d array
      TYPE(MITGCM_ISTATE), POINTER :: miState
      CHARACTER*(*)                :: code
      REAL*4, POINTER              :: dPtr(:,:,:,:)
      INTEGER                      :: rc

      RC = -1

      _HARD_ASSERT( RC )
      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_STATE_GET_R41D(   miState, code, dPtr, rc )
!--   Get pointer to data for real*4 1-d array
      TYPE(MITGCM_ISTATE), POINTER :: miState
      CHARACTER*(*)                :: code
      REAL*4, POINTER              :: dPtr(:)
      INTEGER                      :: rc

      RC = -1

      _HARD_ASSERT( RC )
      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_STATE_GET_R4(   miState, code, dPtr, rc )
!--   Get pointer to data for real*8 5-d array
      TYPE(MITGCM_ISTATE), POINTER :: miState
      CHARACTER*(*)                :: code
      REAL*4, POINTER              :: dPtr
      INTEGER                      :: rc

      RC = -1

      _HARD_ASSERT( RC )
      RETURN
      END SUBROUTINE

      END MODULE
