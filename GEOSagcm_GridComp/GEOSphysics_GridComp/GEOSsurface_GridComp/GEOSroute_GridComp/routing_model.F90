MODULE routing_model

  use MAPL_ConstantsMod,     ONLY: rho => MAPL_RHOWTR

  IMPLICIT NONE

  private

  real, parameter :: M        = 0.45      ! Parameter in hydraulic geometry formula
  real, parameter :: mm       = 0.35      ! Parameter in hydraulic geometry formula 

  public :: river_routing_hyd, SEARCH_DNST, M, mm

CONTAINS

! -------------------------------------------------------------------------------------------------------

  REAL FUNCTION RESCONST (LS, P1, P2)

    IMPLICIT NONE

    REAL, INTENT (IN)    :: LS, P1, P2

    RESCONST  = P1 * ((1./LS)**P2) 

  END FUNCTION RESCONST

! -------------------------------------------------------------------------------------------------------
  
  RECURSIVE SUBROUTINE SEARCH_DNST (K, NCAT_G, DNST, Pfaf_all, DNST_OUT)
    
    implicit none
    
    integer, intent (in)                     :: NCAT_G, K
    integer, intent (in), dimension (NCAT_G) :: Pfaf_all, DNST 
    integer, intent (inout)                  :: DNST_OUT
    
    if (DNST(K) == -1) then 
       DNST_OUT = -1
    else
       
       if(Pfaf_all(DNST(K)) >= 1) then
          DNST_OUT = Pfaf_all(DNST(K)) 
       else
          if(DNST(DNST(K)) == -1) then
             DNST_OUT = -1
          else  
             call SEARCH_DNST (DNST(DNST(K)), NCAT_G, DNST, Pfaf_all, DNST_OUT)
          endif
       endif
    endif
    
    RETURN
    
  END SUBROUTINE SEARCH_DNST

! -------------------------------------------------------------------------------------------------------

  ! ======================================================================================
  !
  ! HYDRAULIC GEOMETRY ROUTING MODEL
  !
  ! --------------------------------------------------------------------------------------
  ! Routing Model Input Parameters
  ! ------------------------------
  !**** NCAT       = NUMBER OF CATCHMENTS IN THE STUDY DOMAIN
  !**** Qrunf0     = RUNOFF PRODUCED BY LAND SURFACE MODEL IN THE CATCHMENT        [m^3/s]
  !**** llc_ori    = MAIN RIVER LENGTH SCALE                                       [m]
  !**** lstr       = LOCAL STREAMS LENGTH SCALE                                    [m]
  !**** qstr_clmt0 = CLIMATOLOGY RUNOFF                                            [m^3/s]  
  !**** qri_clmt0  = CLIMATOLOGY DISCHAR                                           [m^3/s]
  !**** qin_clmt0  = CLIMATOLOGY INFLOW                                            [m^3/s]               
  !**** K          = K PARAMETER FOR MAIN RIVER
  !**** Kstr0      = K PARAMETER FOR LOCAL STREAM                                  [m^3/s]
  
  ! Routing Model Prognostics
  ! -------------------------
  !**** Ws0        = AMOUNT OF WATER IN "LOCAL STREAM"                             [m^3]
  !**** Wr0        = AMOUNT OF WATER IN RIVER                                      [m^3]
  
  ! Routing Model Diagnostics
  ! -------------------------
  !**** QS         = TRANSFER OF MOISTURE FROM STREAM VARIABLE TO RIVER VARIABLE   [m^3/s]
  !**** QOUT       = TRANSFER OF RIVER WATER TO THE DOWNSTREAM CATCHMENT           [m^3/s]  
  
  SUBROUTINE RIVER_ROUTING_HYD (          &
       NCAT,ROUTE_DT,                     &
       Qrunf0,llc_ori,lstr,               &
       qstr_clmt0, qri_clmt0, qin_clmt0,  &
       K, Kstr0,                          &
       Ws0,Wr0,                           &
       Qs,Qout)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)                     :: NCAT,ROUTE_DT
    REAL,    INTENT(IN),   DIMENSION (NCAT) :: Qrunf0,llc_ori,lstr
    REAL,    INTENT(IN),   DIMENSION (NCAT) :: qstr_clmt0,qri_clmt0,qin_clmt0
    REAL,    INTENT(IN),   DIMENSION (NCAT) :: K, Kstr0
    REAL,    INTENT(INOUT),DIMENSION (NCAT) :: Ws0,Wr0
    REAL,    INTENT(OUT),  DIMENSION (NCAT) :: Qs,Qout



    real, parameter :: small    = 1.e-20 
    real, parameter :: fac_kstr = 0.01      ! Factor for local stream scaling

    real,dimension(NCAT) :: Qrunf,qstr_clmt,qri_clmt,qin_clmt,Ws,Wr,Kstr
    real,dimension(NCAT) :: nume,deno,llc,alp_s,alp_r,Qs0,ks,Ws_last
    
    real    :: dt, rho

    integer :: i,j 

    ! convert volume units to mass
    Qrunf     = Qrunf0     * rho          ! m3/s -> kg/s  
    !llc_ori  = llc_ori0   * 1.e3         ! km   -> m
    !lstr     = lstr0      * 1.e3         ! km   -> m
    qstr_clmt = qstr_clmt0 * rho          ! m3/s -> kg/s
    qri_clmt  = qri_clmt0  * rho          ! m3/s -> kg/s
    qin_clmt  = qin_clmt0  * rho          ! m3/s -> kg/s
    Ws        = Ws0        * rho          ! m3   -> kg
    Wr        = Wr0        * rho          ! m3   -> kg

    Kstr = fac_kstr * Kstr0
    dt   = ROUTE_DT

    ! Adjust llc (length of river channel)
    nume = qri_clmt**(2.-M) - qin_clmt**(2.-M)                                                                  ! Numerator for the llc calculation
    deno = (2.-M) * (qri_clmt - qin_clmt) * (qri_clmt**(1.-M))                                                  ! Denominator for the llc calculation
    where(abs(deno) >  small) llc = llc_ori * (nume / deno)                                                     ! Compute llc where denominator is not too small
    where(abs(deno) <= small) llc = llc_ori * 0.5                                                               ! Set llc to half of original value if denominator is small

    ! Calculate alp_s (stream coefficient) and alp_r (river coefficient)
    where(qstr_clmt >  small) alp_s = (rho**(-M) * qstr_clmt**(M-mm) * Kstr * (0.5*lstr)**(-1.))**(1./(1.-mm))  ! For non-zero streamflow
    where(qstr_clmt <= small) alp_s = 0.                                                                        ! If streamflow is too small, set alp_s to 0
    where(qri_clmt  >  small) alp_r = (rho**(-M) * qri_clmt**(M-mm)  * K * llc**(-1.)          )**(1./(1.-mm))  ! For non-zero river input
    where(qri_clmt  <= small) alp_r = 0.                                                                        ! If river input is too small, set alp_r to 0

    ! Update state variables: ks, Ws, and Qs 
    where(Qrunf<=small)Qrunf=0.                                ! Set runoff to zero if it's too small
    Qs0=max(0.,alp_s * Ws**(1./(1.-mm)))                       ! Initial flow from stream storage (kg/s)
    ks =max(0.,(alp_s/(1.-mm)) * Ws**(mm/(1.D0-mm)))           ! Flow coefficient (s^-1)
    Ws_last=Ws                                                 ! Store the current water storage 
    where(ks>small)  Ws=Ws + (Qrunf-Qs0)/ks*(1.-exp(-ks*dt))   ! Update storage (kg)
    where(ks<=small) Ws=Ws + (Qrunf-Qs0)*dt                    ! Simplified update if ks is small
    Ws=max(0.,Ws)                                              ! Ensure storage is non-negative
    Qs=max(0.,Qrunf-(Ws-Ws_last)/dt)                           ! Calculate the stream flow (kg/s)

    ! Calculate variables related to river routing: Qr0, kr
    Wr=Wr+Qs*dt
    Qout=max(0.,alp_r * Wr**(1./(1.-mm)))                      ! River flow based on water storage (kg/s)
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
