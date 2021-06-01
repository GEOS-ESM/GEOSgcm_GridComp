MODULE routing_model
  
  IMPLICIT NONE

  private
  public :: river_routing, SEARCH_DNST, RRM_TIMESTEP
  integer,    parameter :: RRM_TIMESTEP = 3600 ! timestep used for calibrating P1,P2,P3,P4,P5,P6

  CONTAINS

  ! ------------------------------------------------------------------------

  SUBROUTINE RIVER_ROUTING (                &
       NCAT,RRM_DT,                         &
       RUNCATCH,AREACAT,LENGSC,             &
       WSTREAM,WRIVER,                      &
       SFLOW,DISCHARGE)
    
    IMPLICIT NONE
    INTEGER, INTENT(IN)                     :: NCAT
    REAL,    INTENT(IN)                     :: RRM_DT
    REAL,    INTENT(IN),   DIMENSION (NCAT) :: RUNCATCH,AREACAT,LENGSC
    REAL,    INTENT(INOUT),DIMENSION (NCAT) :: WSTREAM, WRIVER
    REAL,    INTENT(OUT),  DIMENSION (NCAT) :: SFLOW,DISCHARGE

    REAL,   PARAMETER    :: K_SIMPLE = 0.111902, K_RES_MAX = 0.8                       ! m1_r2com_c1
    REAL,   PARAMETER    :: P1 = 0.010611, P2 = 0.188556, P3 = 0.096864,   &
                            P4 = 0.691310, P5 = 0.365747, P6 = 0.009831    ! m5_calib_240

    INTEGER :: N,I,J 
    REAL    :: COEFF, LS, COEFF1, COEFF2,ROFF 

    ! Routing Model Input Parameters
    ! ------------------------------
    !**** NCAT      = NUMBER OF CATCHMENTS IN THE STUDY DOMAIN
    !**** RUNCATCH  = RUNOFF PRODUCED BY LAND SURFACE MODEL IN THE CATCHMENT [m3/s]
    !**** AREACAT   = AREA OF CATCHMENT [km^2]
    !**** LENGSC    = LENGTHSCALE OF CATCHMENT FOR RIVER CALCULATION [km]
    !                 Note: We assume LENGSC for stream to river calculation as AREACAT/LENGSC

    ! Routing Model Prognostics
    ! -------------------------
    !**** WSTREAM   = AMOUNT OF WATER IN "LOCAL STREAM"  [m^3]
    !**** WRIVER    = AMOUNT OF WATER IN RIVER           [m^3]

    ! Routing Model Diagnostics
    ! -------------------------
    !**** SFLOW     = WATER TRANSFER RATE FROM STREAMS TO THE LOCAL RIVER CHANNEL [m^3/s]
    !**** INFLOW    = TRANSFER OF RIVER WATER FROM UPSTREAM CATCHMENTS [m^3/s]
    !                 i.e. sum of DISCHARGEs  from all upstream catchments.
    !                 This is computed outside this subroutine.
    !**** DISCHARGE = TRANSFER OF RIVER WATER TO THE DOWNSTREAM CATCHMENT  [m^3/s]

    SFLOW     = 0.
    DISCHARGE = 0.       
 
    DO N=1,NCAT
              
       ! Updating WSTREAM
       
       WSTREAM(N)    = WSTREAM(N)  + RUNCATCH(N) * RRM_DT
       LS            = AREACAT(N) / (AMAX1(1.,LENGSC (N)))
       ROFF          = RUNCATCH(N) * AREACAT(N)
       IF(ROFF < 2. ) THEN
          COEFF = RESCONST (LS, P1, P2)
       ELSEIF(ROFF > 10.) THEN
          COEFF = RESCONST (LS, P3, P4)
       ELSE
          COEFF1 = RESCONST (LS, P1, P2)    
          COEFF2 = RESCONST (LS, P3, P4)   
          COEFF  = COEFF1 + (ROFF - 2.)*(COEFF2 - COEFF1)/8.
       ENDIF

       IF(COEFF > K_RES_MAX) COEFF = K_SIMPLE
 
       SFLOW(N)     = COEFF * WSTREAM(N)
       WSTREAM(N)   = WSTREAM(N) - SFLOW(N)
       WRIVER(N)    = WRIVER(N)  + SFLOW(N)
       SFLOW(N)     = SFLOW(N) / RRM_DT 

       ! Updating WRIVER
       
       LS            = AMAX1(1.,LENGSC (N)) 
       COEFF         = RESCONST (LS, P5, P6)
       IF(COEFF > K_RES_MAX) COEFF = K_SIMPLE 

       DISCHARGE(N)  = COEFF * WRIVER(N)
       WRIVER(N)     = WRIVER(N)   - DISCHARGE(N)
       DISCHARGE(N)  = DISCHARGE(N) / RRM_DT 
       
    ENDDO
   
    RETURN
    
  END SUBROUTINE RIVER_ROUTING

! -------------------------------------------------------------------------------------------

  REAL FUNCTION RESCONST (LS, P1, P2)

    IMPLICIT NONE

    REAL, INTENT (IN)    :: LS, P1, P2

    RESCONST  = P1 * ((1./LS)**P2) 

  END FUNCTION RESCONST

! -------------------------------------------------------------------------------------------
  
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

END MODULE routing_model
