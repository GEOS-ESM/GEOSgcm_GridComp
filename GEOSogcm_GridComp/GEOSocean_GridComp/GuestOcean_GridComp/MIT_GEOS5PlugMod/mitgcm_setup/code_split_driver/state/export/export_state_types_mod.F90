! $Id: $

#include "CPP_OPTIONS.h"

      MODULE EXPORT_STATE_TYPES_MOD
!     Dynamically allocated type for the export information in MITgcm

      IMPLICIT NONE

!     US - surface currents (cell-centered, along latitude circles)
!     VS - surface currents (cell-centered, along longitude lines)
!     SA - sin of angle between local coordinate lines and lat-lon coordinate lines
!     CA - cos of angle between local coordinate lines and lat-lon coordinate lines
!     TS - sea surface temperature in K
!     SS - sea surface salinty     in g/kg
!     MASK - Ocean (=1) / Land (=0) 3-D mask
      TYPE MITGCM_EXPORT  
       SEQUENCE
       _RL , POINTER :: US(   :,:) => NULL()
       _RL , POINTER :: VS(   :,:) => NULL()
       _RL , POINTER :: SA(   :,:) => NULL()
       _RL , POINTER :: CA(   :,:) => NULL()
       _RL , POINTER :: TS(   :,:) => NULL()
       _RL , POINTER :: SS(   :,:) => NULL()
       _RL , POINTER :: MASK(:,:,:) => NULL()
!ALT: the following variable are added to couple the sea-ice
!      _RL, POINTER  :: FRACICE(:,:,:) => NULL()
       _RL, POINTER  :: UI(       :,:) => NULL()
       _RL, POINTER  :: VI(       :,:) => NULL()
!      _RL, POINTER  :: TAUXI(    :,:) => NULL()
!      _RL, POINTER  :: TAUYI(    :,:) => NULL()
       _RL, POINTER  :: TAUXBOT(  :,:) => NULL()
       _RL, POINTER  :: TAUYBOT(  :,:) => NULL()
       _RL, POINTER  :: DELFRAICE(:,:,:) => NULL()
       _RL, POINTER  :: DELVOLICE(:,:,:) => NULL()
       _RL, POINTER  :: DELVOLSNO(:,:,:) => NULL()
       _RL, POINTER  :: DELERGICE(:,:,:) => NULL()
       _RL, POINTER  :: DELERGSNO(:,:,:) => NULL()
       _RL, POINTER  :: DELMPOND (:,:,:) => NULL()
       _RL, POINTER  :: DELTAUAGE(:,:,:) => NULL()
       _RL, POINTER  :: DELTI    (:,:,:) => NULL()
       _RL, POINTER  :: DELSI    (  :,:) => NULL()
       _RL, POINTER  :: DELHI    (  :,:) => NULL()
      END TYPE

      END MODULE
