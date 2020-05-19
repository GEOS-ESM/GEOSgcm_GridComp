! $Id: $

#include "MEM_MANAGE_MACROS.h"

      MODULE DYNVARS_H_STATE_CREATE_DESTROY_MOD

      USE DYNVARS_H_STATE_TYPES_MOD
      IMPLICIT NONE

      INTERFACE CREATE
       MODULE PROCEDURE MITGCM_DYNVARS_H_CREATE
      END INTERFACE
      INTERFACE DESTROY
       MODULE PROCEDURE MITGCM_DYNVARS_H_DESTROY
      END INTERFACE

      CONTAINS

      SUBROUTINE MITGCM_DYNVARS_H_CREATE( dynVarsPtr,      &
                 snx, sny, olx, oly, n3ds, nsx, nsy        &
                 )
!     -- Allocate memory for MITgcm dynvars_h type of a specific size --
      TYPE(MITGCM_DYNVARS_H), POINTER :: dynVarsPtr
      INTEGER                         :: snx, sny, olx, oly
      INTEGER, intent(in), dimension(:) :: n3ds
      INTEGER                         :: nsx, nsy

!     -- Local variables --
      INTEGER                         :: nr
      TYPE(MITGCM_DYNVARS),    POINTER :: p
      TYPE(MITGCM_DYNVARS_R),  POINTER :: p_dr

      nr = n3ds(1)
      ALLOCATE( dynVarsPtr )
      ALLOCATE( dynVarsPtr%dynvars   )
      ALLOCATE( dynVarsPtr%dynvars_r )
      p => dynVarsPtr%dynvars
      ALLOCATE ( p%uVel(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)  )
      ALLOCATE ( p%vVel(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)  )
      ALLOCATE ( p%wVel(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)  )
      ALLOCATE ( p%theta(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy) )
      ALLOCATE ( p%salt(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)  )
      ALLOCATE ( p%gUNM1(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy) )
      ALLOCATE ( p%gVNM1(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy) )
      ALLOCATE ( p%gTNM1(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy) )
      ALLOCATE ( p%gSNM1(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy) )
      ALLOCATE ( p%gU(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)    )
      ALLOCATE ( p%gV(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)    )
      ALLOCATE ( p%uVelD(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy) )
      ALLOCATE ( p%vVelD(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy) )
      ALLOCATE ( p%uNM1(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)  )
      ALLOCATE ( p%vNM1(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)  )
      ALLOCATE ( p%etaN(1-olx:snx+olx,1-oly:sny+oly,nsx,nsy)     )
      ALLOCATE ( p%etaNM1(1-olx:snx+olx,1-oly:sny+oly,nsx,nsy)   )
      p_dr => dynVarsPtr%dynvars_r
      ALLOCATE( p_dr%etaH(1-olx:snx+olx,1-oly:sny+oly,nsx,nsy)   )

      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_DYNVARS_H_DESTROY( dynVarsPtr )
!     -- Allocate memory for an array of MITgcm states of a specific size --
      TYPE(MITGCM_DYNVARS_H), POINTER :: dynVarsPtr

!     -- Local variables --
      TYPE(MITGCM_DYNVARS),  POINTER :: p
      TYPE(MITGCM_DYNVARS_R),  POINTER :: p_dr

      IF ( ASSOCIATED(dynVarsPtr) ) THEN
       p => dynVarsPtr%dynvars
       _DEALLOC ( p%uVel   )
       _DEALLOC ( p%vVel   )
       _DEALLOC ( p%wVel   )
       _DEALLOC ( p%theta  )
       _DEALLOC ( p%salt   )
       _DEALLOC ( p%gUNM1  )
       _DEALLOC ( p%gVNM1  )
       _DEALLOC ( p%gTNM1  )
       _DEALLOC ( p%gSNM1  )
       _DEALLOC ( p%gU     )
       _DEALLOC ( p%gV     )
       _DEALLOC ( p%uVelD  )
       _DEALLOC ( p%vVelD  )
       _DEALLOC ( p%uNM1   )
       _DEALLOC ( p%vNM1   )
       _DEALLOC ( p%etaN   )
       _DEALLOC ( p%etaNM1 )
       p_dr => dynVarsPtr%dynvars_r
       _DEALLOC ( p_dr%etaH )
       _DEALLOC( dynVarsPtr )
      ENDIF

      RETURN
      END SUBROUTINE

      END MODULE
