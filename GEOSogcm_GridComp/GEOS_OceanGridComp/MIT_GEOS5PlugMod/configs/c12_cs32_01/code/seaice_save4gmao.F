#include "SEAICE_OPTIONS.h"

CBOP
C     !ROUTINE: SEAICE_SAVE4GMAO
C     !INTERFACE:
      SUBROUTINE SEAICE_SAVE4GMAO(
     I                  myTime, seqFlag, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R SEAICE_SAVE4GMAO
C     | o Save Seaice Advective Increment for GMAO Coupling
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
# include "SEAICE_PARAMS.h"
# include "SEAICE.h"
# ifdef HACK_FOR_GMAO_CPL
#  include "SEAICE_LAYERS.h"
# endif
#endif /* ALLOW_SEAICE */

C     !INPUT/OUTPUT PARAMETERS:
C     myTime  :: Current time of simulation ( s )
C     seqFlag :: flag that indicate where this S/R is called from:
C             :: =0 reset inc. to zero (from initialise_varia)
C             :: =1 called from the beginning of SEAICE_MODEL
C             :: =2 called from the end of SEAICE_MODEL
C     myIter  :: Iteration number
C     myThid  :: my Thread Id number
      _RL     myTime
      INTEGER seqFlag
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef HACK_FOR_GMAO_CPL
# ifdef ALLOW_SEAICE
C     !LOCAL VARIABLES:
C     i,j,bi,bj :: Loop counters
      INTEGER i, j, bi, bj
      INTEGER l, n
c     _RL tmpFld (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
c     _RL fld    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( seqFlag.EQ.0 ) THEN
C--   Initialise to zero Advective Increments

        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)

          DO n=1,nITD
C-        start loop on ice-category "n"
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
              SIadv_Area  (i,j,n,bi,bj) = 0. _d 0
              SIadv_Heff  (i,j,n,bi,bj) = 0. _d 0
              SIadv_Hsnow (i,j,n,bi,bj) = 0. _d 0
              SIadv_meltPd(i,j,n,bi,bj) = 0. _d 0
              SIadv_iceAge(i,j,n,bi,bj) = 0. _d 0
              SIadv_tIces (i,j,n,bi,bj) = 0. _d 0
            ENDDO
           ENDDO
           DO l=1,nIceLayers
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              SIadv_qIce  (i,j,l,n,bi,bj) = 0. _d 0
             ENDDO
            ENDDO
           ENDDO
           DO l=1,nSnowLayers
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              SIadv_qSnow (i,j,l,n,bi,bj) = 0. _d 0
             ENDDO
            ENDDO
           ENDDO
C-        end loop on ice-category "n"
          ENDDO

          DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
              SIadv_skinS (i,j,bi,bj) = 0. _d 0
              SIadv_skinH (i,j,bi,bj) = 0. _d 0
            ENDDO
          ENDDO

C-        end bi,bj loops
         ENDDO
        ENDDO

      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( seqFlag.EQ.1 ) THEN
C--   Save seaice fields as entering seaice model

        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)

          DO n=1,nITD
C-        start loop on ice-category "n"
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
              SIadv_Area  (i,j,n,bi,bj) = AREAITD (i,j,n,bi,bj)
              SIadv_Heff  (i,j,n,bi,bj) = HEFFITD (i,j,n,bi,bj)
              SIadv_Hsnow (i,j,n,bi,bj) = HSNOWITD(i,j,n,bi,bj)
              SIadv_meltPd(i,j,n,bi,bj) = SImeltPd(i,j,n,bi,bj)
              SIadv_iceAge(i,j,n,bi,bj) = SIiceAge(i,j,n,bi,bj)
              SIadv_tIces (i,j,n,bi,bj) = TICES   (i,j,n,bi,bj)
            ENDDO
           ENDDO
           DO l=1,nIceLayers
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              SIadv_qIce (i,j,l,n,bi,bj) = SIqIce (i,j,l,n,bi,bj)
             ENDDO
            ENDDO
           ENDDO
           DO l=1,nSnowLayers
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              SIadv_qSnow(i,j,l,n,bi,bj) = SIqSnow(i,j,l,n,bi,bj)
             ENDDO
            ENDDO
           ENDDO
C-        end loop on ice-category "n"
          ENDDO

c         DO j=1-OLy,sNy+OLy
c           DO i=1-OLx,sNx+OLx
c             SIadv_skinS (i,j,bi,bj) = SIskinS(i,j,bi,bj)
c             SIadv_skinH (i,j,bi,bj) = SIskinH(i,j,bi,bj)
c           ENDDO
c         ENDDO

C-        end bi,bj loops
         ENDDO
        ENDDO

#ifdef ALLOW_DIAGNOSTICS
        IF ( useDiagnostics ) THEN
         CALL DIAGNOSTICS_FILL( oceWeight,'CPLoWGHT',0,1,0,1,1,myThid )
         CALL DIAGNOSTICS_FILL( TICES,    'SItIces ',
     I                          0, nITD, 0,1,1, myThid )
         n = nIceLayers*nITD
         CALL DIAGNOSTICS_FILL( SIqIce,   'SIqIce  ',
     I                          0, n ,   0,1,1, myThid )
         n = nSnowLayers*nITD
         CALL DIAGNOSTICS_FILL( SIqSnow,  'SIqSnow ',
     I                          0, n ,   0,1,1, myThid )
         CALL DIAGNOSTICS_FILL( SImeltPd, 'SImeltPd',
     I                          0, nITD, 0,1,1, myThid )
         CALL DIAGNOSTICS_FILL( SIiceAge, 'SIiceAge',
     I                          0, nITD, 0,1,1, myThid )
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */

      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( seqFlag.EQ.2 ) THEN
C--   Calculate Advective Increments

        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)

          DO n=1,nITD
C-        start loop on ice-category "n"
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
              SIadv_Area  (i,j,n,bi,bj) = AREAITD (i,j,n,bi,bj)
     &                              - SIadv_Area  (i,j,n,bi,bj)
              SIadv_Heff  (i,j,n,bi,bj) = HEFFITD (i,j,n,bi,bj)
     &                              - SIadv_Heff  (i,j,n,bi,bj)
              SIadv_Hsnow (i,j,n,bi,bj) = HSNOWITD(i,j,n,bi,bj)
     &                              - SIadv_Hsnow (i,j,n,bi,bj)
              SIadv_meltPd(i,j,n,bi,bj) = SImeltPd(i,j,n,bi,bj)
     &                              - SIadv_meltPd(i,j,n,bi,bj)
              SIadv_iceAge(i,j,n,bi,bj) = SIiceAge(i,j,n,bi,bj)
     &                              - SIadv_iceAge(i,j,n,bi,bj)
              SIadv_tIces (i,j,n,bi,bj) = TICES   (i,j,n,bi,bj)
     &                              - SIadv_tIces (i,j,n,bi,bj)
            ENDDO
           ENDDO
           DO l=1,nIceLayers
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              SIadv_qIce (i,j,l,n,bi,bj) = SIqIce (i,j,l,n,bi,bj)
     &                               - SIadv_qIce (i,j,l,n,bi,bj)
             ENDDO
            ENDDO
           ENDDO
           DO l=1,nSnowLayers
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              SIadv_qSnow(i,j,l,n,bi,bj) = SIqSnow(i,j,l,n,bi,bj)
     &                               - SIadv_qSnow(i,j,l,n,bi,bj)
             ENDDO
            ENDDO
           ENDDO
C-        end loop on ice-category "n"
          ENDDO

c         DO j=1-OLy,sNy+OLy
c           DO i=1-OLx,sNx+OLx
c             SIadv_skinS(i,j,bi,bj) = SIskinS(i,j,bi,bj)
c    &                           - SIadv_skinS(i,j,bi,bj)
c             SIadv_skinH(i,j,bi,bj) = SIskinH(i,j,bi,bj)
c    &                           - SIadv_skinH(i,j,bi,bj)
c           ENDDO
c         ENDDO

C-        end bi,bj loops
         ENDDO
        ENDDO

#ifdef ALLOW_DIAGNOSTICS
        IF ( useDiagnostics ) THEN
         CALL DIAGNOSTICS_FILL( SIadv_Area,   'SI_dArea',
     I                          0, nITD, 0,1,1, myThid )
         CALL DIAGNOSTICS_FILL( SIadv_Heff,   'SI_dHeff',
     I                          0, nITD, 0,1,1, myThid )
         CALL DIAGNOSTICS_FILL( SIadv_Hsnow,  'SI_dHsnw',
     I                          0, nITD, 0,1,1, myThid )
         CALL DIAGNOSTICS_FILL( SIadv_tIces,  'SI_dTIce',
     I                          0, nITD, 0,1,1, myThid )
         n = nIceLayers*nITD
         CALL DIAGNOSTICS_FILL( SIadv_qIce,   'SI_dQIce',
     I                          0, n ,   0,1,1, myThid )
         n = nSnowLayers*nITD
         CALL DIAGNOSTICS_FILL( SIadv_qSnow,  'SI_dQSnw',
     I                          0, n ,   0,1,1, myThid )
         CALL DIAGNOSTICS_FILL( SIadv_meltPd, 'SI_dMPnd',
     I                          0, nITD, 0,1,1, myThid )
         CALL DIAGNOSTICS_FILL( SIadv_iceAge, 'SI_dIcAg',
     I                          0, nITD, 0,1,1, myThid )
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */

      ENDIF

# endif /* ALLOW_SEAICE */
#endif /* HACK_FOR_GMAO_CPL */

      RETURN
      END
