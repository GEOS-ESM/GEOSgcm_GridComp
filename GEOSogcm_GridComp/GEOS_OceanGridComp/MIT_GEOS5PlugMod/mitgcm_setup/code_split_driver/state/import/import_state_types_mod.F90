! $Id: $

#include "CPP_OPTIONS.h"

      MODULE IMPORT_STATE_TYPES_MOD
!     Dynamically allocated type for the import information in MITgcm
!     ** No halos **

      IMPLICIT NONE

!     TAUX   - surface stress (cell-centered, for now
!              along EW latitude, +ve E)
!     TAUY   - surface stress (cell-centered, for now
!              along NS longitude, +ve N => a mess at the pole)
!     PS     - surface pressure (cell-centered)
!     SWHEAT - shortwave heating (spread over several levels)
!     QFLX   - freshwater flux (from skin layer)
!     HFLX   - atmos boundary layer heat flux (sensible and latent I think).
!     SFLX   - salt flux (from skin layer)
!     LAT    - cell center latitudes for grid that imports have been interpolated to
!     LON    - cell center longitudes for grid that imports have been interpolated to
      TYPE MITGCM_IMPORT
       SEQUENCE
       _RL , POINTER :: TAUX(     :,:  ) => NULL()
       _RL , POINTER :: TAUY(     :,:  ) => NULL()
       _RL , POINTER :: PS(       :,:  ) => NULL()
       _RL , POINTER :: SWHEAT(   :,:  ) => NULL()
       _RL , POINTER :: QFLX(     :,:  ) => NULL()
       _RL , POINTER :: DISCHARGE(:,:  ) => NULL()
       _RL , POINTER :: HFLX(     :,:  ) => NULL()
       _RL , POINTER :: SFLX(     :,:  ) => NULL()
       _RL , POINTER :: LAT(      :,:  ) => NULL()
       _RL , POINTER :: LON(      :,:  ) => NULL()
       _RL , POINTER :: WGHT(     :,:  ) => NULL()
!ALT: the following variable are added to couple the sea-ice
       _RL , POINTER :: TAUXI(  :,:  ) => NULL()
       _RL , POINTER :: TAUYI(  :,:  ) => NULL()
       _RL, POINTER  :: FRAICE( :,:,:) => NULL()
       _RL, POINTER  :: VOLICE( :,:,:) => NULL()
       _RL, POINTER  :: VOLSNO( :,:,:) => NULL()
       _RL, POINTER  :: ERGICE( :,:,:) => NULL()
       _RL, POINTER  :: ERGSNO( :,:,:) => NULL()
       _RL, POINTER  :: MPOND ( :,:,:) => NULL()
       _RL, POINTER  :: TAUAGE( :,:,:) => NULL()
       _RL, POINTER  :: TI(     :,:,:) => NULL()
       _RL, POINTER  :: SI(     :,:  ) => NULL()
       _RL, POINTER  :: HI(     :,:  ) => NULL()
       _RL, POINTER  :: TS(     :,:  ) => NULL()
      END TYPE

      END MODULE
