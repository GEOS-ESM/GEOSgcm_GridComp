! $Id: $

#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h"
#endif

      MODULE EXPORT_STATE_FILL_MOD

      USE EXPORT_STATE_TYPES_MOD
      IMPLICIT NONE

      INTERFACE FILL
       MODULE PROCEDURE MITGCM_EXPORT_FILL
      END INTERFACE

      CONTAINS

      SUBROUTINE MITGCM_EXPORT_FILL  ( exportPtr )
!     -- Fill the exports for an MITgcm instance.
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
# include "SEAICE_PARAMS.h"
# include "SEAICE.h"
# include "SEAICE_LAYERS.h"
#endif

!     -- Routine arguments    --
!     -- Export state pointer --
      TYPE(MITGCM_EXPORT), POINTER :: exportPtr

!     -- Local variables --
      TYPE(MITGCM_EXPORT), POINTER :: p
      INTEGER bi, bj, i, j, jSLo, iSLo
      INTEGER k, l, n
      _RL tmpU, tmpV

!     -- Put flow vectors on a-grid and rotate to geographical coordinate system.
      p => exportPtr

!     store surface fields to export
      k=1
      DO bj=1,nSy
       DO bi=1,nSx
        jSLo = (bj-1)*sNy
        iSLo = (bi-1)*sNx
        DO j=1,sNy
         DO i=1,sNx

           p%TS(iSLo+i,jSLo+j) = theta(i,j,k,bi,bj) + celsius2K
           p%SS(iSLo+i,jSLo+j) =  salt(i,j,k,bi,bj)

!     Interpolate ocean surface current from C-grid to A-grid locations.
           tmpU = 0.5*(uVel(i  ,j,k,bi,bj) + uVel(i+1,j,k,bi,bj))
           tmpV = 0.5*(vVel(i,j  ,k,bi,bj) + vVel(i,j+1,k,bi,bj))

!     Rotate ocean surface current from MITgcm grid orientation to East/North.
           p%US(iSLo+i,jSLo+j) = angleCosC(i,j,bi,bj) * tmpU
     &                         - angleSinC(i,j,bi,bj) * tmpV
           p%VS(iSLo+i,jSLo+j) = angleSinC(i,j,bi,bj) * tmpU
     &                         + angleCosC(i,j,bi,bj) * tmpV

         ENDDO
        ENDDO

       ENDDO
      ENDDO

!     store 3-D mask to export
      DO bj=1,nSy
       DO bi=1,nSx
        jSLo = (bj-1)*sNy
        iSLo = (bi-1)*sNx
        DO k=1,Nr
         DO j=1,sNy
          DO i=1,sNx
            p%MASK(iSLo+i,jSLo+j,k) =  maskC(i,j,k,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#ifdef ALLOW_SEAICE
      DO bj=1,nSy
       DO bi=1,nSx
        jSLo = (bj-1)*sNy
        iSLo = (bi-1)*sNx

        IF ( SEAICEupdateOceanStress ) THEN
         DO j=1,sNy
          DO i=1,sNx
!     Interpolate seaice velocity from C-grid to A-grid locations.
            tmpU = 0.5*(UICE(i  ,j,bi,bj) + UICE(i+1,j,bi,bj))
            tmpV = 0.5*(VICE(i,j  ,bi,bj) + VICE(i,j+1,bi,bj))
!     Rotate seaice velocity from MITgcm grid orientation to East/North.
            p%UI(iSLo+i,jSLo+j) = angleCosC(i,j,bi,bj) * tmpU
     &                          - angleSinC(i,j,bi,bj) * tmpV
            p%VI(iSLo+i,jSLo+j) = angleSinC(i,j,bi,bj) * tmpU
     &                          + angleCosC(i,j,bi,bj) * tmpV
          ENDDO
         ENDDO
        ENDIF

!     store advection increment of SEAICE thermodynamical fields
        DO n=1,nITD
         DO j=1,sNy
          DO i=1,sNx
!     increment of seaice fractional cover
            p%DELFRAICE(iSLo+i,jSLo+j,n) = SIadv_Area (i,j,n,bi,bj)
!     increment of seaice and snow volume per grid-cell area [m]
            p%DELVOLICE(iSLo+i,jSLo+j,n) = SIadv_Heff (i,j,n,bi,bj)
            p%DELVOLSNO(iSLo+i,jSLo+j,n) = SIadv_Hsnow(i,j,n,bi,bj)
!     increment of melt-pond volume per grid-cell area [m]
            p%DELMPOND (iSLo+i,jSLo+j,n) = SIadv_meltPd(i,j,n,bi,bj)
!     increment of seaice mean age [s]
            p%DELTAUAGE(iSLo+i,jSLo+j,n) = SIadv_iceAge(i,j,n,bi,bj)
!     increment of seaice skin temperature
            p%DELTI    (iSLo+i,jSLo+j,n) = SIadv_tIces(i,j,n,bi,bj)
          ENDDO
         ENDDO
        ENDDO
!     increment of seaice energy
        DO n=1,nITD
         DO l=1,nIceLayers
          k = l + (n-1)*nIceLayers
          DO j=1,sNy
           DO i=1,sNx
            p%DELERGICE(iSLo+i,jSLo+j,k) = SIadv_qIce(i,j,l,n,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
!     increment of snow energy
        DO n=1,nITD
         DO l=1,nSnowLayers
          k = l + (n-1)*nSnowLayers
          DO j=1,sNy
           DO i=1,sNx
            p%DELERGSNO(iSLo+i,jSLo+j,k) = SIadv_qSnow(i,j,l,n,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        DO j=1,sNy
          DO i=1,sNx
!     increment of seaice skin salinity
            p%DELSI(iSLo+i,jSLo+j) = SIadv_skinS(i,j,bi,bj)
!     increment of seaice skin layer depth [m]
            p%DELHI(iSLo+i,jSLo+j) = SIadv_skinH(i,j,bi,bj)
! Add FRZMLT
            p%FRZMLT(iSLo+i,jSLo+j) = SI_FRZMLT(i,j,bi,bj)
          ENDDO
        ENDDO

       ENDDO
      ENDDO
#endif /* ALLOW_SEAICE */

      RETURN
      END SUBROUTINE

      END MODULE
