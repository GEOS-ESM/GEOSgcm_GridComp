! $Id: $

#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h"
#endif

      MODULE IMPORT_STATE_FILL_MOD

      USE IMPORT_STATE_TYPES_MOD
      IMPLICIT NONE

      INTERFACE FILL
       MODULE PROCEDURE MITGCM_IMPORT_FILL
      END INTERFACE

      CONTAINS

      SUBROUTINE MITGCM_IMPORT_FILL  ( importPtr )
!     -- Fill the imports for an MITgcm instance.
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "DYNVARS.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
# include "SEAICE.h"
# include "SEAICE_LAYERS.h"
#endif
#ifdef ALLOW_MYPACKAGE
# include "MYPACKAGE.h"
#endif

!     -- Routine arguments    --
!     -- Import state pointer --
      TYPE(MITGCM_IMPORT), POINTER :: importPtr

!     -- Local variables --
      TYPE(MITGCM_IMPORT), POINTER :: p
      INTEGER bi, bj, i, j, jSLo, iSLo
      INTEGER k, l, n
      INTEGER iMin, iMax, jMin, jMax
      INTEGER n3d, myThid
      PARAMETER ( iMin = 2-OLx , iMax = sNx+OLx )
      PARAMETER ( jMin = 2-OLy , jMax = sNy+OLy )
      _RL tmpU (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tmpV (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tmpW (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tmpD (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL exclVal
      Real*4 exclSng, advlim

      myThid  = 1
      exclSng = 1.E+15
      exclVal = exclSng
!US   limit for seaice advection in terms of land-ocean fraction
!     below this ocean fraction no seaice advection.
!     some land-ocean mask points are not exactly 1.
!     Therefore I chose a arbitrary number
!     smaller than one
      advlim = 0.99

!     Initialize local arrays:
      DO bj=1,nSy
       DO bi=1,nSx
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          tmpU(i,j,bi,bj) = 0.
          tmpV(i,j,bi,bj) = 0.
          tmpW(i,j,bi,bj) = 0.
          tmpD(i,j,bi,bj) = 0.
         ENDDO
        ENDDO
       ENDDO
      ENDDO

!     Transfer fields from importPtr to MITgcm forcing variables.
      p => importPtr
      DO bj=1,nSy
       DO bi=1,nSx
        jSLo = (bj-1)*sNy
        iSLo = (bi-1)*sNx
        DO j=1,sNy
         DO i=1,sNx

!     Transfer WGHT to tmpW (in case we need to mask some imports)
          tmpW(i,j,bi,bj) = p%WGHT(iSLo+i,jSLo+j)

!     Transfer SWHEAT to Qsw.
!     SWHEAT : Solar heating rate spread over several levels in W/m^2.
!              > 0 for increase in theta (ocean warming)
!     Qsw    : Net upward shortwave radiation in W/m^2.
!              Qsw = - ( downward - ice and snow absorption - reflected )
!              > 0 for decrease in theta (ocean cooling)
!              Typical range: -350 < Qsw < 0
!              Southwest C-grid tracer point (i.e., A-grid)
          Qsw(i,j,bi,bj) = -p%SWHEAT(iSLo+i,jSLo+j)

!     Transfer HFLX + SWHEAT to Qnet.
!     HFLX   : Turbulent heat flux from skin to ocean in W/m^2.
!              Accounts for the effect of sensible, latent and radiative
!              fluxes at the interface between the skin and the atmosphere.
!              > 0 for increase in theta (ocean warming)
!     Qnet   : Net upward surface heat flux (including shortwave) in W/m^2.
!              Qnet = latent + sensible + net longwave + net shortwave
!              > 0 for decrease in theta (ocean cooling)
!              Typical range: -250 < Qnet < 600
!              Southwest C-grid tracer point (i.e., A-grid)
          Qnet(i,j,bi,bj) = p%HFLX(iSLo+i,jSLo+j)

!     Transfer QFLX to EMPMR.
!     QFLX   : Freshwater flux from skin to ocean in kg/m2/s.
!              > 0 for decrease in salt (ocean salinity)
!     EmPmR  : Net upward freshwater flux in kg/m2/s.
!              EmPmR = Evaporation - precipitation - runoff
!              > 0 for increase in salt (ocean salinity)
!              Typical range: -1e-4 < EmPmR < 1e-4
!              Southwest C-grid tracer point (i.e., A-grid)
! Note: runoff (=DISCHARGE) has been added to QFLX later on.
          EmPmR(i,j,bi,bj) = - p%QFLX(iSLo+i,jSLo+j)
          tmpD(i,j,bi,bj)  = - p%DISCHARGE(iSLo+i,jSLo+j)

!     Transfer of SFLX to saltFlux
!     SFLX   : Salt flux from skin to ocean in kg/m^2/s
!   saltFlux : Net upward salt flux in g/m^2/s (=psu.kg/m^2/s) .
!              flux of Salt taken out of the ocean per time unit (second).
!              Note: a) only used when salty sea-ice forms or melts.
!                    b) units: when salinity (unit= psu) is expressed
!                       in g/kg, saltFlux unit becomes g/m^2/s.
!              > 0 for decrease in SSS.
!              Southwest C-grid tracer point (i.e., A-grid)
! I am guessing that when there is evaporation:
!    (i) "salty water" flows into the skin layer from the ocean,
!    (ii) freshwater flows from skin layer to atmosphere, and
!    (iii) salt flows from skin layer to ocean.
! If SFLX is this third component of evaporation, it can be ignored for driving
! MITgcm since MITgcm only exports freshwater from ocean during evaporation.
! But SFLUX has both positive and negative numbers, so above decomposition may not be correct.
! Another detail that we will need to worry about later on is that in the presence of sea-ice
! the MITgcm exports “salty water” to the ice and then rejects some of that salt back to
! the ocean.
          saltFlux(i,j,bi,bj) = -1000.D0 * p%SFLX(iSLo+i,jSLo+j)
!         saltFlux(i,j,bi,bj) = 0.0D0

!     Transfer of PS to pLoad
!     PS     : Surface atmospheric pressure in Pa
!     pLoad  : for the ocean:      atmospheric pressure at z=eta
!              Units are           Pa=N/m^2
! Code crash because of GEOS/MITgcm land mask inconsistency if pLoad=p%PS.
!         pLoad(i,j,bi,bj) = p%PS(iSLo+i,jSLo+j)
!         pLoad(i,j,bi,bj) = 0.0D0
         IF ( tmpW(i,j,bi,bj).EQ.0.D0 ) THEN 
           pLoad(i,j,bi,bj) = 1.D+5
         ELSE
           pLoad(i,j,bi,bj) = p%PS(iSLo+i,jSLo+j)/tmpW(i,j,bi,bj)
         ENDIF

!     Transfer TS to theta (for dual ocean).
!     TS      : SST in K
!     theta   : SST in DegC
!         theta(i,j,1,bi,bj) = p%TS(iSLo+i,jSLo+j)-celsius2K


!     Transfer TAUX and TAUY to fu and fv.
!     TAUX   : A-grid eastward stress on skin in N/m^2.
!     TAUY   : A-grid northward stress on skin in N/m^2.
!     fu     : Zonal surface wind stress in N/m^2.
!              > 0 for increase in uVel
!              Typical range: -0.5 < fu < 0.5
!              Southwest C-grid U point
!     fv     : Meridional surface wind stress in N/m^2.
!              > 0 for increase in vVel
!              Typical range: -0.5 < fv < 0.5
!              Southwest C-grid V point

!     Rotate stress from East/North to MITgcm grid orientation.
          tmpU(i,j,bi,bj) =
     &         angleCosC(i,j,bi,bj) * p%TAUX(iSLo+i,jSLo+j)
     &       + angleSinC(i,j,bi,bj) * p%TAUY(iSLo+i,jSLo+j)
          tmpV(i,j,bi,bj) =
     &       - angleSinC(i,j,bi,bj) * p%TAUX(iSLo+i,jSLo+j)
     &       + angleCosC(i,j,bi,bj) * p%TAUY(iSLo+i,jSLo+j)

         ENDDO
        ENDDO
       ENDDO
      ENDDO

!     Interpolate stress from A-grid to C-grid locations.
      CALL EXCH_UV_AGRID_3D_RL( tmpU, tmpV, .TRUE., 1, myThid )
      DO bj=1,nSy
       DO bi=1,nSx
        DO j = jMin,jMax
         DO i = iMin,iMax
          fu(i,j,bi,bj) = (tmpU(i,j,bi,bj)+tmpU(i-1,j,bi,bj))
     &                  * 0.5D0 * maskW(i,j,1,bi,bj)
          fv(i,j,bi,bj) = (tmpV(i,j,bi,bj)+tmpV(i,j-1,bi,bj))
     &                  * 0.5D0 * maskS(i,j,1,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#ifdef ALLOW_MYPACKAGE
      IF ( useMYPACKAGE ) THEN
       IF ( myPa_index1.GT.0 ) THEN
        DO n=1,myPa_index1
C--
         CALL EXCH_XY_RL( tmpD, myThid )
         DO bj=1,nSy
          DO bi=1,nSx

           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             tmpU(i,j,bi,bj) = tmpD(i,j,bi,bj)*rA(i,j,bi,bj)
            ENDDO
           ENDDO
           DO j=1,sNy
            DO i=1,sNx
             tmpD(i,j,bi,bj) = tmpD(i,j,bi,bj) - myPa_param1
     &        * ( ( tmpU(i,j,bi,bj) - tmpU(i-1,j,bi,bj) )
     &            *maskW(i,j,1,bi,bj)
     &           +( tmpU(i,j,bi,bj) - tmpU(i+1,j,bi,bj) )
     &            *maskW(i+1,j,1,bi,bj)
     &           +( tmpU(i,j,bi,bj) - tmpU(i,j-1,bi,bj) )
     &            *maskS(i,j,1,bi,bj)
     &           +( tmpU(i,j,bi,bj) - tmpU(i,j+1,bi,bj) )
     &            *maskS(i,j+1,1,bi,bj)
     &          )*recip_rA(i,j,bi,bj)
            ENDDO
           ENDDO

          ENDDO
         ENDDO
C--
        ENDDO

       ENDIF
      ENDIF
#endif
      DO bj=1,nSy
       DO bi=1,nSx
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           EmPmR(i,j,bi,bj) = EmPmR(i,j,bi,bj)+tmpD(i,j,bi,bj)
          ENDDO
         ENDDO
       ENDDO
      ENDDO
      CALL EXCH_XY_RS( EmPmR, myThid )
      CALL EXCH_XY_RS( pLoad, myThid )

#ifdef ALLOW_SEAICE
!     Transfer fields from importPtr to MITgcm seaice variables.
      DO bj=1,nSy
       DO bi=1,nSx
        jSLo = (bj-1)*sNy
        iSLo = (bi-1)*sNx

        DO j=1,sNy
          DO i=1,sNx
!     Rotate stress from East/North to MITgcm grid orientation.
            tmpU(i,j,bi,bj) =
     &         angleCosC(i,j,bi,bj) * p%TAUXI(iSLo+i,jSLo+j)
     &       + angleSinC(i,j,bi,bj) * p%TAUYI(iSLo+i,jSLo+j)
            tmpV(i,j,bi,bj) =
     &       - angleSinC(i,j,bi,bj) * p%TAUXI(iSLo+i,jSLo+j)
     &       + angleCosC(i,j,bi,bj) * p%TAUYI(iSLo+i,jSLo+j)
          ENDDO
        ENDDO

        DO n=1,nITD
         DO j=1,sNy
          DO i=1,sNx
!     fractional_cover_of_seaice
            AREAITD (i,j,n,bi,bj) = p%FRAICE(iSLo+i,jSLo+j,n)
!     seaice and snow volume per grid-cell area [m]
            HEFFITD (i,j,n,bi,bj) = p%VOLICE(iSLo+i,jSLo+j,n)
            HSNOWITD(i,j,n,bi,bj) = p%VOLSNO(iSLo+i,jSLo+j,n)
!     melt-pond volume per grid-cell area [m]
            SImeltPd(i,j,n,bi,bj) = p%MPOND(iSLo+i,jSLo+j,n)
!     seaice mean age [s]
            SIiceAge(i,j,n,bi,bj) = p%TAUAGE(iSLo+i,jSLo+j,n)
!     seaice skin temperature
            TICES   (i,j,n,bi,bj) = p%TI(iSLo+i,jSLo+j,n)
          ENDDO
         ENDDO
        ENDDO
!     seaice energy
        DO n=1,nITD
         DO l=1,nIceLayers
          k = l + (n-1)*nIceLayers
          DO j=1,sNy
           DO i=1,sNx
            SIqIce(i,j,l,n,bi,bj) = p%ERGICE(iSLo+i,jSLo+j,k)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
!     snow energy
        DO n=1,nITD
         DO l=1,nSnowLayers
          k = l + (n-1)*nSnowLayers
          DO j=1,sNy
           DO i=1,sNx
            SIqSnow(i,j,l,n,bi,bj) = p%ERGSNO(iSLo+i,jSLo+j,k)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        DO j=1,sNy
          DO i=1,sNx
!     seaice skin salinity
            SIskinS(i,j,bi,bj) = p%SI(iSLo+i,jSLo+j)
!     seaice skin-layer depth
            SIskinH(i,j,bi,bj) = p%HI(iSLo+i,jSLo+j)
!     Transfer WGHT to oceWeight (just for diagnostic purpose)
            oceWeight(i,j,bi,bj) = tmpW(i,j,bi,bj)
          ENDDO
        ENDDO

!     Remove Exclude value from seaice fields:
        DO n=1,nITD
         DO j=1,sNy
          DO i=1,sNx
c          IF ( AREAITD(i,j,n,bi,bj).EQ.exclVal ) THEN
           IF ( tmpW(i,j,bi,bj).LT.advlim ) THEN
            AREAITD (i,j,n,bi,bj) = 0.D0
            HEFFITD (i,j,n,bi,bj) = 0.D0
            HSNOWITD(i,j,n,bi,bj) = 0.D0
            SImeltPd(i,j,n,bi,bj) = 0.D0
            SIiceAge(i,j,n,bi,bj) = 0.D0
            TICES   (i,j,n,bi,bj) = 0.D0
            DO l=1,nIceLayers
             SIqIce(i,j,l,n,bi,bj) = 0.D0
            ENDDO
            DO l=1,nSnowLayers
             SIqSnow(i,j,l,n,bi,bj)= 0.D0
            ENDDO
            SIskinS(i,j,bi,bj) = 0.D0
            SIskinH(i,j,bi,bj) = 0.D0
           ENDIF
          ENDDO
         ENDDO
        ENDDO

!     end bi,bj loops
       ENDDO
      ENDDO

!     Apply EXCH to seaice fields:
      CALL EXCH_UV_AGRID_3D_RL( tmpU, tmpV, .TRUE., 1, myThid )
      CALL EXCH_3D_RL( AREAITD , nITD, myThid )
      CALL EXCH_3D_RL( HEFFITD , nITD, myThid )
      CALL EXCH_3D_RL( HSNOWITD, nITD, myThid )
      CALL EXCH_3D_RL( SImeltPd, nITD, myThid )
      CALL EXCH_3D_RL( SIiceAge, nITD, myThid )
      CALL EXCH_3D_RL( TICES   , nITD, myThid )
      n3d = nIceLayers*nITD
      CALL EXCH_3D_RL( SIqIce , n3d, myThid )
      n3d = nSnowLayers*nITD
      CALL EXCH_3D_RL( SIqSnow, n3d, myThid )
      CALL EXCH_XY_RL( SIskinS, myThid )
      CALL EXCH_XY_RL( SIskinH, myThid )

!     Initialize seaice mask
      DO bj=1,nSy
       DO bi=1,nSx
        DO j=1,sNy
         DO i=1,sNx
          ! US some landocean mask points are not exactly 1
          ! Therefore I chose a arbitrary number
          ! smaller than one
          IF ( tmpW(i,j,bi,bj).LT.advlim ) THEN 
           HEFFM(i,j,bi,bj) = 0.D0
          ELSE
           HEFFM(i,j,bi,bj) = 1.D0
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      CALL EXCH_XY_RL( HEFFM, myThid )

      DO bj=1,nSy
       DO bi=1,nSx
        DO j=1,sNy
         DO i=1,sNx
          SIMaskU(i,j,bi,bj)=min(HEFFM(i-1,j,bi,bj),HEFFM(i,j,bi,bj))
          SIMaskV(i,j,bi,bj)=min(HEFFM(i,j-1,bi,bj),HEFFM(i,j,bi,bj))
          seaiceMaskU(i,j,bi,bj)=SIMaskU(i,j,bi,bj)
          seaiceMaskV(i,j,bi,bj)=SIMaskV(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      CALL EXCH_UV_XY_RL( SIMaskU, SIMaskV, .FALSE., myThid )


      DO bj=1,nSy
       DO bi=1,nSx
!     Interpolate stress from A-grid to C-grid locations.
        DO j = jMin,jMax
         DO i = iMin,iMax
          SIwindTauX(i,j,bi,bj) = (tmpU(i,j,bi,bj)+tmpU(i-1,j,bi,bj))
     &                  * 0.5D0 * SIMaskU(i,j,bi,bj)
          SIwindTauY(i,j,bi,bj) = (tmpV(i,j,bi,bj)+tmpV(i,j-1,bi,bj))
     &                  * 0.5D0 * SIMaskV(i,j,bi,bj)
         ENDDO
        ENDDO
!     Fill-in grid-cell mean from all seaice categories
        CALL SEAICE_ITD_SUM( bi, bj, zeroRL, 0, myThid )
!     end bi,bj loops
       ENDDO
      ENDDO

#endif /* ALLOW_SEAICE */

      RETURN
      END SUBROUTINE

      END MODULE
