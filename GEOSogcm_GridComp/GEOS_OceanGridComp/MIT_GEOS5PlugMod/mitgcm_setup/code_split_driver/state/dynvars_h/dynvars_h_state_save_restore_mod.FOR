! $Id: $

#include "CPP_OPTIONS.h"

      MODULE DYNVARS_H_STATE_SAVE_RESTORE_MOD

      USE DYNVARS_H_STATE_TYPES_MOD 
      IMPLICIT NONE

      INTERFACE SAVE_STATE
       MODULE PROCEDURE DYNVARS_H_SAVE_STATE
      END INTERFACE
      INTERFACE RESTORE_STATE
       MODULE PROCEDURE DYNVARS_H_RESTORE_STATE
      END INTERFACE

      CONTAINS

      SUBROUTINE DYNVARS_H_SAVE_STATE( dynVarsPtr )          
!     -- Save DYNVARS values to internal state --
#include "SIZE.h"
#include "DYNVARS.h"
!     -- Routine arguments --
      TYPE(MITGCM_DYNVARS_H), POINTER :: dynVarsPtr
!     -- Local variables --
      TYPE(MITGCM_DYNVARS), POINTER :: p
      TYPE(MITGCM_DYNVARS_R),  POINTER :: p_dr

      p => dynVarsPtr%dynvars
      p%uVel    = uVel
      p%vVel    = vVel
      p%wVel    = wVel
      p%theta   = theta
      p%salt    = salt 
      p%gUNM1   = gUNM1
      p%gVNM1   = gVNM1
      p%gTNM1   = gTNM1
      p%gSNM1   = gSNM1
      p%gU      = gU   
      p%gV      = gV   
      p%etaN    = etaN 
      p_dr => dynVarsPtr%dynvars_r
      p_dr%etaH = etaH

      RETURN
      END SUBROUTINE

      SUBROUTINE DYNVARS_H_RESTORE_STATE( dynVarsPtr )          
!     -- Restore DYNVARS values from internal state --
#include "SIZE.h"
#include "DYNVARS.h"
!     -- Routine arguments --
      TYPE(MITGCM_DYNVARS_H), POINTER :: dynVarsPtr
!     -- Local variables --
      TYPE(MITGCM_DYNVARS), POINTER :: p
      TYPE(MITGCM_DYNVARS_R),  POINTER :: p_dr

      p => dynVarsPtr%dynvars
      uVel  = p%uVel
      vVel  = p%vVel
      wVel  = p%wVel
      theta = p%theta
      salt  = p%salt 
      gUNM1 = p%gUNM1
      gVNM1 = p%gVNM1
      gTNM1 = p%gTNM1
      gSNM1 = p%gSNM1
      gU    = p%gU
      gV    = p%gV
      etaN  = p%etaN
      p_dr => dynVarsPtr%dynvars_r
      etaH  = p_dr%etaH
     
      RETURN
      END SUBROUTINE

      END MODULE
