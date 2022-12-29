#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

      MODULE DRIVER_GET_EXPORT_STATE_MOD
C     A GEOS-5 compatible version of get export state.
C     Exports variable from MITgcm that are used by other
C     components in a MAPL tree.
C     Mapping from MITgcm five-dimensional to MAPL data
C     structures takes place here.

#define G5KIND      4
#define REAL_       real(kind=G5KIND)

      INTERFACE DRIVER_GET_EXPORT_STATE
       MODULE PROCEDURE DRIVER_GET_EXPORT_STATE_2D
       MODULE PROCEDURE DRIVER_GET_EXPORT_STATE_3D
      END INTERFACE

      CONTAINS
CBOP
C     !ROUTINE: DRIVER_GET_EXPORT_STATE
C     !INTERFACE:
      SUBROUTINE DRIVER_GET_EXPORT_STATE_2D( mitgcmIState, code, fld)
C     !USES:
      USE MITGCM_STATE_MOD , ONLY :
     &    MITGCM_ISTATE
      IMPLICIT NONE
#include "SIZE.h"

C     == Routine arguments ==
      TYPE(MITGCM_ISTATE), POINTER :: mitgcmIState
      CHARACTER*(*) code
      REAL_, pointer :: fld(:,:)

C     == Local variables ==
      INTEGER i, j, bi, bj, iSLo, jSLo

      IF ( code .EQ. 'US' ) THEN
!     Copy US to fld
         fld = mitgcmIState%export%US

      ELSEIF ( code .EQ. 'VS' ) THEN
!     Copy VS to fld
         fld = mitgcmIState%export%VS

      ELSEIF ( code .EQ. 'TS' ) THEN
!     Copy TS to fld
         fld = mitgcmIState%export%TS

      ELSEIF ( code .EQ. 'SS' ) THEN
!     Copy SS to fld
         fld = mitgcmIState%export%SS

      ELSEIF ( code .EQ. 'TAUXBOT' ) THEN
!     Copy TAUXBOT to fld
         fld = 0.0 !@@@@ mitgcmIState%export%TAUXBOT
      ELSEIF ( code .EQ. 'TAUYBOT' ) THEN
!     Copy TAUYBOT to fld
         fld = 0.0 !@@@@ mitgcmIState%export%TAUYBOT

#ifdef ALLOW_SEAICE
      ELSEIF ( code .EQ. 'UI' ) THEN
!     Copy UI to fld
         fld = mitgcmIState%export%UI
      ELSEIF ( code .EQ. 'VI' ) THEN
!     Copy VI to fld
         fld = mitgcmIState%export%VI

      ELSEIF ( code .EQ. 'DELSI' ) THEN
!     Copy DELSI to fld
         fld = mitgcmIState%export%DELSI
      ELSEIF ( code .EQ. 'DELHI' ) THEN
!     Copy DELHI to fld
         fld = mitgcmIState%export%DELHI
      ELSEIF ( code .EQ. 'FRZMLT' ) THEN
!     Copy FRZMLT to fld
         fld = mitgcmIState%export%FRZMLT
#else  /* ALLOW_SEAICE */
      ELSEIF ( code .EQ. 'UI'
     &    .OR. code .EQ. 'VI'
     &    .OR. code .EQ. 'DELSI'
     &    .OR. code .EQ. 'FRZMLT'
     &    .OR. code .EQ. 'DELHI' ) THEN
         fld = 0.0
#endif /* ALLOW_SEAICE */
      ELSE

       print *,'DRIVER_GET_EXPORT_STATE_2D:'//trim(code)
     &      //'. This code is stubbed!'

      ENDIF

      RETURN
      END SUBROUTINE

      SUBROUTINE DRIVER_GET_EXPORT_STATE_3D( mitgcmIState, code, fld)
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

      IF ( code .EQ. 'MASK' ) THEN
       IF ( ASSOCIATED(fld) ) THEN
             fld = mitgcmIState%export%MASK
       ELSE
        print *,'DRIVER_GET_EXPORT_STATE_3D:'//trim(code)
     &      //' : Pointer not associated!'
       ENDIF

#ifdef ALLOW_SEAICE
      ELSEIF ( code .EQ. 'DELFRACICE' ) THEN
!     Copy DELFRAICE to fld
         fld = mitgcmIState%export%DELFRAICE
      ELSEIF ( code .EQ. 'DELVOLICE' ) THEN
!     Copy DELVOLICE to fld
         fld = mitgcmIState%export%DELVOLICE
      ELSEIF ( code .EQ. 'DELVOLSNO' ) THEN
!     Copy DELVOLSNO to fld
         fld = mitgcmIState%export%DELVOLSNO
      ELSEIF ( code .EQ. 'DELERGICE' ) THEN
!     Copy DELERGICE to fld
         fld = mitgcmIState%export%DELERGICE
      ELSEIF ( code .EQ. 'DELERGSNO' ) THEN
!     Copy DELERGSNO to fld
         fld = mitgcmIState%export%DELERGSNO
      ELSEIF ( code .EQ. 'DELMPOND' ) THEN
!     Copy DELMPOND to fld
         fld = mitgcmIState%export%DELMPOND
      ELSEIF ( code .EQ. 'DELTAUAGE' ) THEN
!     Copy DELTAUAGE to fld
         fld = mitgcmIState%export%DELTAUAGE
      ELSEIF ( code .EQ. 'DELTI' ) THEN
!     Copy DELTI to fld
         fld = mitgcmIState%export%DELTI
#else  /* ALLOW_SEAICE */
      ELSEIF ( code .EQ. 'DELFRACICE'
     &    .OR. code .EQ. 'DELVOLICE'
     &    .OR. code .EQ. 'DELVOLSNO'
     &    .OR. code .EQ. 'DELERGICE'
     &    .OR. code .EQ. 'DELERGSNO'
     &    .OR. code .EQ. 'DELMPOND'
     &    .OR. code .EQ. 'DELTAUAGE'
     &    .OR. code .EQ. 'DELTI' ) THEN
         fld = 0.0
#endif /* ALLOW_SEAICE */
      ELSE

       print *,'DRIVER_GET_EXPORT_STATE_3D:'//trim(code)
     &      //'. This code is stubbed!'

      ENDIF

      RETURN
      END SUBROUTINE

      END MODULE
