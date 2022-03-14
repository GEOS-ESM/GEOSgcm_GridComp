! $Id: $

#include "CPP_OPTIONS.h"
      MODULE DYNVARS_H_STATE_TYPES_MOD

!     Dynamically allocated type for the DYNVARS information in MITgcm
      IMPLICIT NONE

      TYPE MITGCM_DYNVARS
       SEQUENCE
       _RL , POINTER :: uVel(   :,:,:,:,:) => NULL()
       _RL , POINTER :: vVel(   :,:,:,:,:) => NULL()
       _RL , POINTER :: wVel(   :,:,:,:,:) => NULL()
       _RL , POINTER :: theta(  :,:,:,:,:) => NULL()
       _RL , POINTER :: salt(   :,:,:,:,:) => NULL()
       _RL , POINTER :: guNM1(  :,:,:,:,:) => NULL()
       _RL , POINTER :: gvNM1(  :,:,:,:,:) => NULL()
       _RL , POINTER :: gtNM1(  :,:,:,:,:) => NULL()
       _RL , POINTER :: gsNM1(  :,:,:,:,:) => NULL()
       _RL , POINTER :: gU(     :,:,:,:,:) => NULL()
       _RL , POINTER :: gv(     :,:,:,:,:) => NULL()
       _RL , POINTER :: uVelD(  :,:,:,:,:) => NULL()
       _RL , POINTER :: vVelD(  :,:,:,:,:) => NULL()
       _RL , POINTER :: uNM1(   :,:,:,:,:) => NULL()
       _RL , POINTER :: vNM1(   :,:,:,:,:) => NULL()
       _RL , POINTER :: etaN(   :,:,:,:)   => NULL()
       _RL , POINTER :: etaNM1( :,:,:,:)   => NULL()
      END TYPE

      TYPE MITGCM_DYNVARS_A_TILE
       SEQUENCE
       _RL , POINTER :: uVel(   :,:,:) => NULL()
       _RL , POINTER :: vVel(   :,:,:) => NULL()
      END TYPE

      TYPE MIGCM_DYNVARS_H_TILESET
       SEQUENCE
       TYPE(MITGCM_DYNVARS_A_TILE), POINTER :: dynvars(:,:)
      END TYPE

      TYPE MITGCM_DYNVARS_R
       SEQUENCE
       _RL , POINTER :: etaH(:,:,:,:) => NULL()
      END TYPE

      TYPE MITGCM_DYNVARS_H
       SEQUENCE
       TYPE(MITGCM_DYNVARS) ,   POINTER :: DYNVARS   => NULL()
       TYPE(MITGCM_DYNVARS_R) , POINTER :: DYNVARS_R => NULL()
      END TYPE

      END MODULE
