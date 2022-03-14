#include "MEM_MANAGE_MACROS.h"

      MODULE FFIELDS_H_STATE_MOD
#include "CPP_OPTIONS.h"

      TYPE FFIELDS_FU
       SEQUENCE
       _RS, POINTER :: FU(:,:,:,:)  => NULL()
      END TYPE

      TYPE FFIELDS_FV
       SEQUENCE
       _RS, POINTER :: FV(:,:,:,:)  => NULL()
      END TYPE

      TYPE FFIELDS_QNET
       SEQUENCE
       _RS, POINTER :: QNET(:,:,:,:)  => NULL()
      END TYPE

      TYPE FFIELDS_QSW
       SEQUENCE
       _RS, POINTER :: QSW(:,:,:,:)  => NULL()
      END TYPE

      TYPE FFIELDS_EMPMR
       SEQUENCE
       _RS, POINTER :: EMPMR(:,:,:,:)  => NULL()
      END TYPE

      TYPE FFIELDS_SALTFLUX
       SEQUENCE
       _RS, POINTER :: SALTFLUX(:,:,:,:)  => NULL()
      END TYPE

      TYPE FFIELDS_SST
       SEQUENCE
       _RS, POINTER :: SST(:,:,:,:)  => NULL()
      END TYPE

      TYPE FFIELDS_SSS
       SEQUENCE
       _RS, POINTER :: SSS(:,:,:,:)  => NULL()
      END TYPE

      TYPE FFIELDS_LAMBDATHETACLIMRELAX
       SEQUENCE
       _RS, POINTER :: LAMBDATHETACLIMRELAX(:,:,:,:)  => NULL()
      END TYPE

      TYPE FFIELDS_LAMBDASALTCLIMRELAX
       SEQUENCE
       _RS, POINTER :: LAMBDASALTCLIMRELAX(:,:,:,:)  => NULL()
      END TYPE

#ifdef ATMOSPHERIC_LOADING
      TYPE FFIELDS_PLOAD
       SEQUENCE
       _RS, POINTER :: PLOAD(:,:,:,:)  => NULL()
      END TYPE

      TYPE FFIELDS_SICELOAD
       SEQUENCE
       _RS, POINTER :: SICELOAD(:,:,:,:)  => NULL()
      END TYPE
#endif

#ifndef ALLOW_EXF
      TYPE TDFIELDS
       SEQUENCE
       _RS , POINTER :: taux0    (:,:,:,:) => NULL()
       _RS , POINTER :: tauy0    (:,:,:,:) => NULL()
       _RS , POINTER :: Qnet0    (:,:,:,:) => NULL()
       _RS , POINTER :: EmPmR0   (:,:,:,:) => NULL()
       _RS , POINTER :: saltFlux0(:,:,:,:) => NULL()
       _RS , POINTER :: SST0     (:,:,:,:) => NULL()
       _RS , POINTER :: SSS0     (:,:,:,:) => NULL()
       _RS , POINTER :: taux1    (:,:,:,:) => NULL()
       _RS , POINTER :: tauy1    (:,:,:,:) => NULL()
       _RS , POINTER :: Qnet1    (:,:,:,:) => NULL()
       _RS , POINTER :: EmPmR1   (:,:,:,:) => NULL()
       _RS , POINTER :: saltFlux1(:,:,:,:) => NULL()
       _RS , POINTER :: SST1     (:,:,:,:) => NULL()
       _RS , POINTER :: SSS1     (:,:,:,:) => NULL()
#ifdef ATMOSPHERIC_LOADING
       _RS , POINTER :: pload0   (:,:,:,:) => NULL()
       _RS , POINTER :: pload1   (:,:,:,:) => NULL()
#endif
#ifdef SHORTWAVE_HEATING
       _RS , POINTER :: Qsw1     (:,:,:,:) => NULL()
       _RS , POINTER :: Qsw0     (:,:,:,:) => NULL()
#endif
      END TYPE
#endif

      TYPE SURFACE_FORCING
       SEQUENCE
       _RL , POINTER :: surfaceForcingU(:,:,:,:) => NULL()
       _RL , POINTER :: surfaceForcingV(:,:,:,:) => NULL()
       _RL , POINTER :: surfaceForcingT(:,:,:,:) => NULL()
       _RL , POINTER :: surfaceForcingS(:,:,:,:) => NULL()
       _RL , POINTER :: surfaceForcingTice(:,:,:,:) => NULL()
      END TYPE

      TYPE MITGCM_FFIELDS_H
       SEQUENCE
       TYPE(FFIELDS_FU)       , POINTER :: FFIELDS_FU => NULL()
       TYPE(FFIELDS_FV)       , POINTER :: FFIELDS_FV => NULL()
       TYPE(FFIELDS_QNET)     , POINTER :: FFIELDS_QNET => NULL()
       TYPE(FFIELDS_QSW)      , POINTER :: FFIELDS_QSW   => NULL()
       TYPE(FFIELDS_EMPMR)    , POINTER :: FFIELDS_EMPMR => NULL()
       TYPE(FFIELDS_SALTFLUX) , POINTER :: FFIELDS_SALTFLUX => NULL()
       TYPE(FFIELDS_SST)      , POINTER :: FFIELDS_SST       => NULL()
       TYPE(FFIELDS_SSS)      , POINTER :: FFIELDS_SSS       => NULL()
       TYPE(FFIELDS_LAMBDATHETACLIMRELAX) , POINTER :: FFIELDS_LAMBDATHETACLIMRELAX => NULL()
       TYPE(FFIELDS_LAMBDASALTCLIMRELAX)  , POINTER :: FFIELDS_LAMBDASALTCLIMRELAX => NULL()
#ifdef ALLOW_ATMOSPHERIC_LOADING
       TYPE(FFIELDS_PLOAD)    , POINTER :: FFIELDS_PLOAD => NULL()
       TYPE(FFIELDS_SICELOAD) , POINTER :: FFIELDS_SICELOAD => NULL()
#endif
#ifndef ALLOW_EXF
       TYPE(TDFIELDS)         , POINTER :: TDFIELDS => NULL()
#endif
       TYPE(SURFACE_FORCING)  , POINTER :: SURFACE_FORCING => NULL()
      END TYPE

      INTERFACE CREATE
       MODULE PROCEDURE MITGCM_FFIELDS_H_CREATE
      END INTERFACE
      INTERFACE DESTROY
       MODULE PROCEDURE MITGCM_FFIELDS_H_DESTROY
      END INTERFACE

      CONTAINS

      SUBROUTINE MITGCM_FFIELDS_H_CREATE( ffieldshPtr,      &
                 sNx, sNy, OLx, OLy, n3ds, nSx, nSy         &
                 )
!--   Create a forcing fields type, allocating space based on the
!--   per instance mesh size.
      TYPE(MITGCM_FFIELDS_H), POINTER :: ffieldshPtr
      INTEGER                :: sNx, sNy, OLx, OLy
      INTEGER, dimension(:)  :: n3ds
      INTEGER                :: nSx, nSy

      RETURN
      END SUBROUTINE

      SUBROUTINE MITGCM_FFIELDS_H_DESTROY( ffieldshPtr )
!--   Free memory for a forcing fields type.
      TYPE(MITGCM_FFIELDS_H), POINTER :: ffieldshPtr

      IF ( ASSOCIATED(ffieldshPtr) ) THEN
       _DEALLOC( ffieldshPtr%ffields_fu%fu )
       _DEALLOC( ffieldshPtr%ffields_fu%fu )
       _DEALLOC( ffieldshPtr%ffields_fu                      )

       _DEALLOC( ffieldshPtr%ffields_fv%fv                   )
       _DEALLOC( ffieldshPtr%ffields_fv                      )

       _DEALLOC( ffieldshPtr%ffields_qnet%qnet               )
       _DEALLOC( ffieldshPtr%ffields_qnet                    )

       _DEALLOC( ffieldshPtr%ffields_qsw%qsw                 )
       _DEALLOC( ffieldshPtr%ffields_qsw                     )

       _DEALLOC( ffieldshPtr%ffields_empmr%empmr             )
       _DEALLOC( ffieldshPtr%ffields_empmr                   )

       _DEALLOC( ffieldshPtr%ffields_saltflux%saltflux       )
       _DEALLOC( ffieldshPtr%ffields_saltflux                )

       _DEALLOC( ffieldshPtr%ffields_sst%sst                 )
       _DEALLOC( ffieldshPtr%ffields_sst                     )

       _DEALLOC( ffieldshPtr%ffields_sss%sss                 )
       _DEALLOC( ffieldshPtr%ffields_sss                     )

       _DEALLOC( ffieldshPtr%ffields_lambdathetaclimrelax%lambdathetaclimrelax )
       _DEALLOC( ffieldshPtr%ffields_lambdathetaclimrelax    )

       _DEALLOC( ffieldshPtr%ffields_lambdasaltclimrelax%lambdasaltclimrelax )
       _DEALLOC( ffieldshPtr%ffields_lambdasaltclimrelax     )

#ifdef ALLOW_ATMOSPHERIC_LOADING
       _DEALLOC( ffieldshPtr%ffields_pload%pload             )
       _DEALLOC( ffieldshPtr%ffields_pload                   )
       _DEALLOC( ffieldshPtr%ffields_sIceLoad%IceLoad        )
       _DEALLOC( ffieldshPtr%ffields_sIceLoad                )
#endif

#ifndef ALLOW_EXF
       _DEALLOC( ffieldshPtr%tdfields%taux0                  )
       _DEALLOC( ffieldshPtr%tdfields%tauy0                  )
       _DEALLOC( ffieldshPtr%tdfields%qnet0                  )
       _DEALLOC( ffieldshPtr%tdfields%empmr0                 )
       _DEALLOC( ffieldshPtr%tdfields%saltflux0              )
       _DEALLOC( ffieldshPtr%tdfields%sst0                   )
       _DEALLOC( ffieldshPtr%tdfields%sss0                   )
       _DEALLOC( ffieldshPtr%tdfields%taux1                  )
       _DEALLOC( ffieldshPtr%tdfields%tauy1                  )
       _DEALLOC( ffieldshPtr%tdfields%qnet1                  )
       _DEALLOC( ffieldshPtr%tdfields%empmr1                 )
       _DEALLOC( ffieldshPtr%tdfields%saltflux1              )
       _DEALLOC( ffieldshPtr%tdfields%sst1                   )
       _DEALLOC( ffieldshPtr%tdfields%sss1                   )
#ifdef ATMOSPHERIC_LOADING
       _DEALLOC( ffieldshPtr%tdfields%pload0                 )
       _DEALLOC( ffieldshPtr%tdfields%pload1                 )
#endif
#ifdef SHORTWAVE_HEATING
       _DEALLOC( ffieldshPtr%tdfields%qsw0                   )
       _DEALLOC( ffieldshPtr%tdfields%qsw1                   )
#endif
       _DEALLOC( ffieldshPtr%tdfields                        )
#endif

       _DEALLOC( ffieldshPtr%surface_forcing%surfaceForcingU )
       _DEALLOC( ffieldshPtr%surface_forcing%surfaceForcingV )
       _DEALLOC( ffieldshPtr%surface_forcing%surfaceForcingT )
       _DEALLOC( ffieldshPtr%surface_forcing%surfaceForcingS )
       _DEALLOC( ffieldshPtr%surface_forcing%surfaceForcingT )
       _DEALLOC( ffieldshPtr%surface_forcing                 )

       _DEALLOC( ffieldshPtr                                 )
      ENDIF

      RETURN
      END SUBROUTINE

      END MODULE
