MODULE routing_model
  
  IMPLICIT NONE

  private
  public :: river_routing, SEARCH_DNST, ROUTE_DT
  integer                   ,    parameter :: ROUTE_DT = 3600

  CONTAINS


  ! ------------------------------------------------------------------------

  SUBROUTINE RIVER_ROUTING (                &
       NCAT,                                &
       RUNCATCH,AREACAT,LENGSC,             &
       WSTREAM,WRIVER,                      &
       QSFLOW,QOUTFLOW)
    
    IMPLICIT NONE
    INTEGER, INTENT(IN)                     :: NCAT
    REAL,    INTENT(IN),   DIMENSION (NCAT) :: RUNCATCH,AREACAT,LENGSC
    REAL,    INTENT(INOUT),DIMENSION (NCAT) :: WSTREAM, WRIVER
    REAL,    INTENT(OUT),  DIMENSION (NCAT) :: QSFLOW,QOUTFLOW

    REAL,   PARAMETER    :: K_SIMPLE = 0.111902, K_RES_MAX = 0.8                       ! m1_r2com_c1
    REAL                 :: P1 = 0.136305, P2 = 0.965619, P3 = 0.364546, P4 = 0.000578 ! m4_r2com_c3  

    INTEGER :: N,I,J 
    REAL    :: COEFF, LS 

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
    !**** QSFLOW    = TRANSFER OF MOISTURE FROM STREAM VARIABLE TO RIVER VARIABLE [m^3/s]
    !**** QINFLOW   = TRANSFER OF RIVER WATER FROM UPSTREAM CATCHMENTS [m^3/s] - i.e. sum of
    !                 QOUTFLOWs from all upstream catchments. This is computed outside this subroutine
    !**** QOUTFLOW  = TRANSFER OF RIVER WATER TO THE DOWNSTREAM CATCHMENT  [m^3/s]

    QSFLOW   = 0.
    QOUTFLOW = 0.       
 
    DO N=1,NCAT
              
       ! Updating WSTREAM
       
       WSTREAM(N)    = WSTREAM(N)  + RUNCATCH(N) * REAL (ROUTE_DT)
       LS            = AREACAT(N) / (AMAX1(1.,LENGSC (N)))
       COEFF         = RESCONST (LS, P1, P2)   
       
       IF(COEFF > K_RES_MAX) COEFF = K_SIMPLE
 
       QSFLOW(N)     = COEFF * WSTREAM(N)
       WSTREAM(N)    = WSTREAM(N) - QSFLOW(N)
       WRIVER(N)     = WRIVER(N)  + QSFLOW(N)
       QSFLOW(N)     = QSFLOW(N) / REAL (ROUTE_DT) 

       ! Updating WRIVER
       
       LS            = AMAX1(1.,LENGSC (N)) 
       COEFF         = RESCONST (LS, P3, P4)
       IF(COEFF > K_RES_MAX) COEFF = K_SIMPLE 

       QOUTFLOW(N)   = COEFF * WRIVER(N)
       WRIVER(N)     = WRIVER(N)   - QOUTFLOW(N)
       QOUTFLOW(N)   = QOUTFLOW(N) / REAL (ROUTE_DT) 
       
    ENDDO
   
    RETURN
    
  END SUBROUTINE RIVER_ROUTING

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

END MODULE routing_model
