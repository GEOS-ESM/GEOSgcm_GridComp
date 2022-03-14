      MODULE TIMEVARS_STATE_CREATE_DESTROY_MOD
#include "CPP_OPTIONS.h"

      USE TIMEVARS_STATE_TYPES_MOD
      IMPLICIT NONE

      INTERFACE CREATE
       MODULE PROCEDURE TIMEVARS_CREATE
      END INTERFACE
      INTERFACE DESTROY
       MODULE PROCEDURE TIMEVARS_DESTROY
      END INTERFACE

      CONTAINS

      SUBROUTINE TIMEVARS_CREATE( timeVarsPtr, n3ds )
!     -- Allocate memory for MITgcm timevars variables. --
      TYPE(MITGCM_TIMEVARS), POINTER :: timeVarsPtr
      INTEGER, intent(in), dimension(:) :: n3ds

!     -- Local variables --
      INTEGER nr
      TYPE(MITGCM_TIMEVARS),  POINTER :: p

      nr = n3ds(1)
      ALLOCATE( timeVarsPtr )
      ALLOCATE( timeVarsPtr%dtTracerLev(nr) )

      RETURN
      END SUBROUTINE

      SUBROUTINE TIMEVARS_DESTROY( timeVarsPtr )
!     -- Deallocate memory for MITgcm timevars variables. --
      TYPE(MITGCM_TIMEVARS), POINTER :: timeVarsPtr

!     -- Local variables --
      TYPE(MITGCM_TIMEVARS),  POINTER :: p

      IF ( ASSOCIATED(timeVarsPtr) ) THEN
       DEALLOCATE( timeVarsPtr%dtTracerLev )
       DEALLOCATE( timeVarsPtr )
       NULLIFY( timeVarsPtr )
      ENDIF

      RETURN
      END SUBROUTINE

      END MODULE
