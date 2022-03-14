      MODULE MITGCM_STATE_TYPES_MOD
!     -- MITgcm state type --
      USE DYNVARS_H_STATE_MOD , ONLY :    &
          MITGCM_DYNVARS_H

      USE TIMEVARS_STATE_MOD , ONLY :    &
          MITGCM_TIMEVARS

      USE STACKVARS_STATE_MOD , ONLY :    &
          MITGCM_STACKVARS

      USE STACKVARS_STATE_MOD , ONLY :    &
          MITGCM_STACKVARS

      USE FFIELDS_H_STATE_MOD , ONLY :    &
          MITGCM_FFIELDS_H

      USE EXPORT_STATE_MOD , ONLY :       &
          MITGCM_EXPORT    

      USE IMPORT_STATE_MOD , ONLY :       &
          MITGCM_IMPORT    

      IMPLICIT NONE

!     -- Type that holds internal state --
!     -- DYNVARS   :: Copy of DYNVARS.h internal state.
!     -- TIMEVARS  :: Copy of PARAMS.h information on time.
!     -- STACKVARS :: Copy of stack variable state informtaion.
!     -- FFIELDS   :: Copy of FFIELDS.h information.
!     -- EXPORT    :: Export fields exchanged with a coupling layer.
!     -- IMPORT    :: Import fields exchanged with a coupling layer.
      TYPE MITGCM_ISTATE
       SEQUENCE
       TYPE(MITGCM_DYNVARS_H), POINTER :: dynvars_h => NULL()
       TYPE(MITGCM_TIMEVARS),  POINTER :: timevars  => NULL()
       TYPE(MITGCM_STACKVARS), POINTER :: stackvars => NULL()
       TYPE(MITGCM_FFIELDS_H), POINTER :: ffields_h => NULL()
       TYPE(MITGCM_EXPORT),    POINTER :: export    => NULL()
       TYPE(MITGCM_IMPORT),    POINTER :: import    => NULL()
      END TYPE

      TYPE MITGCM_MULTI_ISTATE
       SEQUENCE
       TYPE(MITGCM_ISTATE), POINTER :: ISTATE_LIST(:)
      END TYPE

      TYPE MITGCM_ISTATE_CONTAINER
       SEQUENCE
       TYPE(MITGCM_ISTATE), POINTER :: P => NULL()
      END TYPE

      TYPE MITGCM_ISTATE_WRAP_TYPE
       TYPE(MITGCM_ISTATE), pointer :: Ptr
      END TYPE


      END MODULE
