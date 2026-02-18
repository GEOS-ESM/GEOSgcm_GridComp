MODULE routing_model

  use MAPL_ConstantsMod,     ONLY: rho => MAPL_RHOWTR

  IMPLICIT NONE

  private

  real, parameter :: RRM_M        = 0.45      ! Parameter in hydraulic geometry formula
  real, parameter :: RRM_mm       = 0.35      ! Parameter in hydraulic geometry formula 

  public :: river_routing_hyd, SEARCH_DNST, RRM_M, RRM_mm

CONTAINS

!!  REAL FUNCTION RESCONST (LS, P1, P2)
!!
!!    IMPLICIT NONE
!!
!!    REAL, INTENT (IN)    :: LS, P1, P2
!!
!!    RESCONST  = P1 * ((1./LS)**P2) 
!!
!!  END FUNCTION RESCONST
!!
!!! -------------------------------------------------------------------------------------------------------
!!  
!!  RECURSIVE SUBROUTINE SEARCH_DNST (K, NCAT_G, DNST, Pfaf_all, DNST_OUT)
!!
!!    ! objective:  [best guess by rreichle as of 3 Feb 2026]
!!    ! find Pfafstetter index of downstream (DNST) catchment, where "downstream" is really "downriver" in the
!!    !   jargon of the routing model with "main rivers" and "local streams"  
!!    
!!    implicit none
!!    
!!    integer, intent (in)                     :: NCAT_G, K
!!    integer, intent (in), dimension (NCAT_G) :: Pfaf_all, DNST 
!!    integer, intent (inout)                  :: DNST_OUT
!!    
!!    if (DNST(K) == -1) then 
!!       DNST_OUT = -1
!!    else
!!       
!!       if(Pfaf_all(DNST(K)) >= 1) then
!!          DNST_OUT = Pfaf_all(DNST(K)) 
!!       else
!!          if(DNST(DNST(K)) == -1) then
!!             DNST_OUT = -1
!!          else  
!!             call SEARCH_DNST (DNST(DNST(K)), NCAT_G, DNST, Pfaf_all, DNST_OUT)
!!          endif
!!       endif
!!    endif
!!    
!!    RETURN
!!    
!!  END SUBROUTINE SEARCH_DNST

  
  ! ======================================================================================
  !
  ! HYDRAULIC GEOMETRY ROUTING MODEL
  !
  ! --------------------------------------------------------------------------------------
  ! Routing Model Input Parameters
  ! ------------------------------
  !**** NCAT          = NUMBER OF CATCHMENTS IN THE STUDY DOMAIN
  !**** ROUTE_DT      = TIME STEP FOR ROUTING MODEL                                       [s]
  !**** Qrunf0        = RUNOFF PRODUCED BY LAND SURFACE MODEL IN THE CATCHMENT            [m^3/s]              
  !**** RRM_ALPHA_RIV = ALPHA PARAMETER FOR MAIN RIVER                                        
  !**** RRM_ALPHA_STR = ALPHA PARAMETER FOR LOCAL STREAM                               
                                                                                       
  ! Routing Model Prognostics                                                          
  ! -------------------------                                                          
  !**** Ws0           = AMOUNT OF WATER IN "LOCAL STREAM"                                 [m^3]
  !**** Wr0           = AMOUNT OF WATER IN RIVER                                          [m^3]
  
  ! Routing Model Diagnostics
  ! -------------------------
  !**** QS            = TRANSFER OF MOISTURE FROM STREAM VARIABLE TO RIVER VARIABLE       [m^3/s]
  !**** QOUT          = TRANSFER OF RIVER WATER TO THE DOWNSTREAM (DOWNRIVER) CATCHMENT   [m^3/s]  
  
  SUBROUTINE RIVER_ROUTING_HYD (             &
       NCAT,ROUTE_DT,                        &
       Qrunf0, RRM_ALPHA_RIV, RRM_ALPHA_STR, &
       Ws0,Wr0,                              &
       Qs,Qout)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)                     :: NCAT,ROUTE_DT
    REAL,    INTENT(IN),   DIMENSION (NCAT) :: Qrunf0
    REAL,    INTENT(IN),   DIMENSION (NCAT) :: RRM_ALPHA_RIV, RRM_ALPHA_STR
    REAL,    INTENT(INOUT),DIMENSION (NCAT) :: Ws0,Wr0
    REAL,    INTENT(OUT),  DIMENSION (NCAT) :: Qs,Qout

    real, parameter                         :: small    = 1.e-20 

    real, dimension(NCAT)                   :: Qrunf,Ws,Wr
    real, dimension(NCAT)                   :: Qs0,ks,Ws_last
    
    real                                    :: dt

    ! convert volume units to mass
    Qrunf     = Qrunf0     * rho          ! m3/s -> kg/s  
    Ws        = Ws0        * rho          ! m3   -> kg
    Wr        = Wr0        * rho          ! m3   -> kg

    dt        = ROUTE_DT                  ! integer -> real                                                                      ! If river input is too small, set alp_r to 0

    ! Update state variables: ks, Ws, and Qs 
    where(Qrunf<=small)Qrunf=0.                                            ! Set runoff to zero if it's too small
    Qs0=max(0.,RRM_ALPHA_STR * Ws**(1./(1.-RRM_mm)))                       ! Initial flow from local stream storage (kg/s)
    ks =max(0.,(RRM_ALPHA_STR/(1.-RRM_mm)) * Ws**(RRM_mm/(1.-RRM_mm)))     ! Flow coefficient (s^-1)
    Ws_last=Ws                                                             ! Store the current water storage 
    where(ks>small)  Ws=Ws + (Qrunf-Qs0)/ks*(1.-exp(-ks*dt))               ! Update storage (kg)
    where(ks<=small) Ws=Ws + (Qrunf-Qs0)*dt                                ! Simplified update if ks is small
    Ws=max(0.,Ws)                                                          ! Ensure storage is non-negative
    Qs=max(0.,Qrunf-(Ws-Ws_last)/dt)                                       ! Calculate the local stream flow (kg/s)

    ! Calculate variables related to river routing: Qr0, kr
    Wr=Wr+Qs*dt
    Qout=max(0.,RRM_ALPHA_RIV * Wr**(1./(1.-RRM_mm)))                      ! River flow based on water storage (kg/s)
    Qout=min(Qout,Wr/dt)
    Wr=max(0.,Wr-Qout*dt) 

    ! convert mass units back to volume
    Ws0  = Ws  /rho        ! kg   -> m3
    Wr0  = Wr  /rho        ! kg   -> m3
    Qs   = Qs  /rho        ! kg/s -> m3/s
    Qout = Qout/rho        ! kg/s -> m3/s

    RETURN

  END SUBROUTINE RIVER_ROUTING_HYD

  ! -------------------------------------------------------------------------------------------------------

END MODULE routing_model

! ======================== EOF =========================================================
