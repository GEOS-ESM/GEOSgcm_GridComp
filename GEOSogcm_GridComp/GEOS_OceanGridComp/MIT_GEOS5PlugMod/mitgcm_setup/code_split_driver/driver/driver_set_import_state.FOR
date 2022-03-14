#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

      MODULE DRIVER_SET_IMPORT_STATE_MOD
C     A GEOS-5 compatible version of set import state.
C     Imports variables from MAPL that are used as
C     boundary/forcing conditions of MITgcm.
C     Mapping from MAPL data structure to MITgcm 
C     five-dimensional structures takes place here.
C
C     o called from "plug" code as
C     CALL DRIVER_SET_IMPORT_STATE( wrap%ptr,   'TAUX',   TAUX )
C
C     Temp until BNDY_STATE code is written
      USE IMPORT_STATE_MOD

#define G5KIND      4
#define REAL_       real(kind=G5KIND)

      INTERFACE DRIVER_SET_IMPORT_STATE
       MODULE PROCEDURE DRIVER_SET_IMPORT_STATE_2D
       MODULE PROCEDURE DRIVER_SET_IMPORT_STATE_3D
      END INTERFACE

      CONTAINS
CBOP
C     !ROUTINE: DRIVER_SET_IMPORT_STATE
C     !INTERFACE:
      SUBROUTINE DRIVER_SET_IMPORT_STATE_2D( mitgcmIState, code, fld)
C     !USES:
      USE MITGCM_STATE_MOD , ONLY :
     &    MITGCM_ISTATE
#include "SIZE.h"

C     == Routine arguments ==
      TYPE(MITGCM_ISTATE), POINTER :: mitgcmIState
      CHARACTER*(*) code
      REAL_, pointer :: fld(:,:)

C     == Local variables ==
!      INTEGER i, j, bi, bj, iSLo, jSLo

      IF ( code .EQ. 'TAUX' ) THEN
!      Copy fld to TAUX
         mitgcmIState%import%TAUX = fld
      ELSEIF ( code .EQ. 'TAUY' ) THEN
!      Copy fld to TAUY
         mitgcmIState%import%TAUY = fld
      ELSEIF ( code .EQ. 'TAUXI' ) THEN
!      Copy fld to TAUXI
         mitgcmIState%import%TAUXI = fld
      ELSEIF ( code .EQ. 'TAUYI' ) THEN
!      Copy fld to TAUYI
         mitgcmIState%import%TAUYI = fld
      ELSEIF ( code .EQ. 'PS' ) THEN
!      Copy fld to PS
         mitgcmIState%import%PS = fld
      ELSEIF ( code .EQ. 'QFLX' ) THEN
!      Copy fld to QFLX
         mitgcmIState%import%QFLX = fld
      ELSEIF ( code .EQ. 'DISCHARGE' ) THEN
!      Copy fld to DISCHARGE
         mitgcmIState%import%DISCHARGE = fld
      ELSEIF ( code .EQ. 'HFLX' ) THEN
!      Copy fld to HFLX
         mitgcmIState%import%HFLX = fld
      ELSEIF ( code .EQ. 'SFLX' ) THEN
!      Copy fld to SFLX
         mitgcmIState%import%SFLX = fld
      ELSEIF ( code .EQ. 'LATS' ) THEN
!      Copy fld to LAT
         mitgcmIState%import%LAT = fld
      ELSEIF ( code .EQ. 'LONS' ) THEN
!      Copy fld to LON
         mitgcmIState%import%LON = fld
      ELSEIF ( code .EQ. 'WGHT' ) THEN
!      Copy fld to WGHT
         mitgcmIState%import%WGHT = fld
      ELSEIF ( code .EQ. 'SI' ) THEN
!      Copy fld to SI
         mitgcmIState%import%SI = fld
      ELSEIF ( code .EQ. 'HI' ) THEN
!      Copy fld to HI
         mitgcmIState%import%HI = fld
      ELSEIF ( code .EQ. 'TS' ) THEN
!     Copy fld to TS
         mitgcmIState%import%TS = fld
      ELSEIF ( code .EQ. 'SWHEAT' ) THEN
!      Copy fld to SWHEAT
         mitgcmIState%import%SWHEAT = fld
      ENDIF

C
      RETURN
      END SUBROUTINE

      SUBROUTINE DRIVER_SET_IMPORT_STATE_3D( mitgcmIState, code, fld)
C     !USES:
      USE MITGCM_STATE_MOD , ONLY :
     &    MITGCM_ISTATE
      IMPLICIT NONE
#include "SIZE.h"

C     == Routine arguments ==
      TYPE(MITGCM_ISTATE), POINTER :: mitgcmIState
      CHARACTER*(*) code
      REAL_, pointer :: fld(:,:,:)

C     == Local variables ==
      INTEGER i, j, k, bi, bj, iSLo, jSLo

#ifdef ALLOW_SEAICE
      IF ( code .EQ. 'FRACICE' ) THEN
!      Copy fld to FRACICE
         mitgcmIState%import%FRAICE = fld
      ELSEIF ( code .EQ. 'TI' ) THEN
!      Copy fld to TI
         mitgcmIState%import%TI = fld
      ELSEIF ( code .EQ. 'VOLICE' ) THEN
!      Copy fld to VOLICE
         mitgcmIState%import%VOLICE = fld
      ELSEIF ( code .EQ. 'VOLSNO' ) THEN
!      Copy fld to VOLSNO
         mitgcmIState%import%VOLSNO = fld
      ELSEIF ( code .EQ. 'ERGICE' ) THEN
!      Copy fld to ERGICE
         mitgcmIState%import%ERGICE = fld
      ELSEIF ( code .EQ. 'ERGSNO' ) THEN
!      Copy fld to ERGSNO
         mitgcmIState%import%ERGSNO = fld
      ELSEIF ( code .EQ. 'MPOND' ) THEN
!      Copy fld to MPOND
         mitgcmIState%import%MPOND = fld
      ELSEIF ( code .EQ. 'TAUAGE' ) THEN
!      Copy fld to TAUAGE
         mitgcmIState%import%TAUAGE = fld
      ENDIF
#endif /* ALLOW_SEAICE */

      RETURN
      END SUBROUTINE
      END MODULE

