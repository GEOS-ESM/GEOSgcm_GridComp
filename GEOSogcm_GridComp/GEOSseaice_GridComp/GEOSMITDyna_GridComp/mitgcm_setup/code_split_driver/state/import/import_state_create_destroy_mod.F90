! $Id: $

#include "MEM_MANAGE_MACROS.h"

      MODULE IMPORT_STATE_CREATE_DESTROY_MOD

      USE IMPORT_STATE_TYPES_MOD
      IMPLICIT NONE

      INTERFACE CREATE
       MODULE PROCEDURE MITGCM_IMPORT_CREATE
      END INTERFACE
      INTERFACE DESTROY
       MODULE PROCEDURE MITGCM_IMPORT_DESTROY
      END INTERFACE

      CONTAINS

      SUBROUTINE MITGCM_IMPORT_CREATE( importPtr,        &
                 snx, sny, olx, oly, n3ds, nsx, nsy      &
                 )
!     -- Allocate memory for MITgcm import state type of a specific size --
      TYPE(MITGCM_IMPORT), POINTER :: importPtr
      INTEGER                      :: snx, sny, olx, oly
      INTEGER, intent(in), dimension(:) :: n3ds
      INTEGER                      :: nsx, nsy

!     -- Local variables --
      INTEGER                      :: nr
      INTEGER                      :: nItd, nIla, nSla
      TYPE(MITGCM_IMPORT), POINTER :: p

      nr = n3ds(1)
      nItd = n3ds(2)
      nIla = n3ds(3)*nItd
      nSla = n3ds(4)*nItd
      ALLOCATE( importPtr )
      p => importPtr
      ALLOCATE ( p%TAUX( 1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%TAUY( 1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%TAUXI(1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%TAUYI(1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%PS(   1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%SWHEAT(1:snx*nSx,1:sny*nSy) )
      ALLOCATE ( p%QFLX( 1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%DISCHARGE( 1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%HFLX( 1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%SFLX( 1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%LAT(  1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%LON(  1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%WGHT( 1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%FRAICE(1:snx*nSx,1:sny*nSy,1:nItd) )
      ALLOCATE ( p%VOLICE(1:snx*nSx,1:sny*nSy,1:nItd) )
      ALLOCATE ( p%VOLSNO(1:snx*nSx,1:sny*nSy,1:nItd) )
      ALLOCATE ( p%ERGICE(1:snx*nSx,1:sny*nSy,1:nIla) )
      ALLOCATE ( p%ERGSNO(1:snx*nSx,1:sny*nSy,1:nSla) )
      ALLOCATE ( p%MPOND( 1:snx*nSx,1:sny*nSy,1:nItd) )
      ALLOCATE ( p%TAUAGE(1:snx*nSx,1:sny*nSy,1:nItd) )
      ALLOCATE ( p%TI(   1:snx*nSx,1:sny*nSy,1:nItd) )
      ALLOCATE ( p%SI(   1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%HI(   1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%TS(   1:snx*nSx,1:sny*nSy)  )

      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_IMPORT_DESTROY( importPtr )
!     -- Deallocate memory for an array of MITgcm imports.
      TYPE(MITGCM_IMPORT), POINTER :: importPtr

!     -- Local variables --
      TYPE(MITGCM_IMPORT), POINTER :: p

      IF ( ASSOCIATED(importPtr) ) THEN
       p => importPtr
       _DEALLOC ( p%TAUX   )
       _DEALLOC ( p%TAUY   )
       _DEALLOC ( p%TAUXI  )
       _DEALLOC ( p%TAUYI  )
       _DEALLOC ( p%PS     )
       _DEALLOC ( p%SWHEAT )
       _DEALLOC ( p%QFLX   )
       _DEALLOC ( p%DISCHARGE)
       _DEALLOC ( p%HFLX   )
       _DEALLOC ( p%SFLX   )
       _DEALLOC ( p%LAT    )
       _DEALLOC ( p%LON    )
       _DEALLOC ( p%WGHT   )
       _DEALLOC ( p%FRAICE )
       _DEALLOC ( p%VOLICE )
       _DEALLOC ( p%VOLSNO )
       _DEALLOC ( p%ERGICE )
       _DEALLOC ( p%ERGSNO )
       _DEALLOC ( p%MPOND  )
       _DEALLOC ( p%TAUAGE )
       _DEALLOC ( p%TI     )
       _DEALLOC ( p%SI     )
       _DEALLOC ( p%HI     )
       _DEALLOC ( p%TS     )
       _DEALLOC( importPtr )
      ENDIF

      RETURN
      END SUBROUTINE

      END MODULE
