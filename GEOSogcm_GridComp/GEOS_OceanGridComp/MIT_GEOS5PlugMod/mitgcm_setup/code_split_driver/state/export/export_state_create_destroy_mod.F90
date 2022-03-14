! $Id: $

#include "MEM_MANAGE_MACROS.h"

      MODULE EXPORT_STATE_CREATE_DESTROY_MOD

      USE EXPORT_STATE_TYPES_MOD
      IMPLICIT NONE

      INTERFACE CREATE
       MODULE PROCEDURE MITGCM_EXPORT_CREATE
      END INTERFACE
      INTERFACE DESTROY
       MODULE PROCEDURE MITGCM_EXPORT_DESTROY
      END INTERFACE

      CONTAINS

      SUBROUTINE MITGCM_EXPORT_CREATE( exportPtr,        &
                 snx, sny, olx, oly, n3ds, nSx, nSy      &
                 )
!     -- Allocate memory for MITgcm dynvars_h type of a specific size --
      TYPE(MITGCM_EXPORT), POINTER :: exportPtr
      INTEGER                      :: snx, sny, olx, oly
      INTEGER, intent(in), dimension(:) :: n3ds
      INTEGER                      :: nSx, nSy

!     -- Local variables --
      INTEGER                      :: Nr
      INTEGER                      :: nItd, nIla, nSla
      TYPE(MITGCM_EXPORT), POINTER :: p

      Nr = n3ds(1)
      nItd = n3ds(2)
      nIla = n3ds(3)*nItd
      nSla = n3ds(4)*nItd
      ALLOCATE( exportPtr )
      p => exportPtr
      ALLOCATE ( p%US(1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%VS(1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%CA(1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%SA(1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%TS(1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%SS(1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%MASK(1:snx*nSx,1:sny*nSy,1:Nr)  )

      ALLOCATE ( p%UI(1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%VI(1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%DELFRAICE(1:snx*nSx,1:sny*nSy,1:nItd) )
      ALLOCATE ( p%DELVOLICE(1:snx*nSx,1:sny*nSy,1:nItd) )
      ALLOCATE ( p%DELVOLSNO(1:snx*nSx,1:sny*nSy,1:nItd) )
      ALLOCATE ( p%DELERGICE(1:snx*nSx,1:sny*nSy,1:nIla) )
      ALLOCATE ( p%DELERGSNO(1:snx*nSx,1:sny*nSy,1:nSla) )
      ALLOCATE ( p%DELMPOND (1:snx*nSx,1:sny*nSy,1:nItd) )
      ALLOCATE ( p%DELTAUAGE(1:snx*nSx,1:sny*nSy,1:nItd) )
      ALLOCATE ( p%DELTI(   1:snx*nSx,1:sny*nSy,1:nItd)  )
      ALLOCATE ( p%DELSI(   1:snx*nSx,1:sny*nSy)  )
      ALLOCATE ( p%DELHI(   1:snx*nSx,1:sny*nSy)  )

      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_EXPORT_DESTROY( exportPtr )
!     -- Deallocate memory for an array of MITgcm exports.
      TYPE(MITGCM_EXPORT), POINTER :: exportPtr

!     -- Local variables --
      TYPE(MITGCM_EXPORT), POINTER :: p

      IF ( ASSOCIATED(exportPtr) ) THEN
       p => exportPtr
       _DEALLOC ( p%US     )
       _DEALLOC ( p%VS     )
       _DEALLOC ( p%CA     )
       _DEALLOC ( p%SA     )
       _DEALLOC ( p%TS     )
       _DEALLOC ( p%SS     )
       _DEALLOC ( p%MASK   )
       _DEALLOC ( p%UI     )
       _DEALLOC ( p%VI     )
       _DEALLOC ( p%DELFRAICE )
       _DEALLOC ( p%DELVOLICE )
       _DEALLOC ( p%DELVOLSNO )
       _DEALLOC ( p%DELERGICE )
       _DEALLOC ( p%DELERGSNO )
       _DEALLOC ( p%DELMPOND  )
       _DEALLOC ( p%DELTAUAGE )
       _DEALLOC ( p%DELTI )
       _DEALLOC ( p%DELSI )
       _DEALLOC ( p%DELHI )
       _DEALLOC( exportPtr )
      ENDIF

      RETURN
      END SUBROUTINE

      END MODULE
