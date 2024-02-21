#include "MAPL_Generic.h"

!==========================================================
module SurfParams
  
  ! Justin, 12 Apr 2018  - Created - Replaces LAND_UPD ifdefs , added SurfParams_init,
  !                        called from GEOS_LandGridCompMod, GEOS_LandiceGridCompMod
  ! jkolassa, 11 Jun 2020 - added LSM_CHOICE as input; introduced separate land parameter collections
  !
  use MAPL_ExceptionHandling
  implicit none
  
  private
  
  public :: SurfParams_init
  
  ! ---------------------------------------------------------------------------
  ! Switch for catchment physics parameters
  ! ---------------------------------------------------------------------------
  
  
  ! Variables set by SurfParams_init:
  REAL, PUBLIC, SAVE :: CSOIL_2      		 ! J/K - heat capacity associated w/ tsurf
  REAL, PUBLIC, SAVE :: WEMIN        		 ! kg/m^2  minimum SWE in areal fraction
  REAL, PUBLIC, SAVE :: AICEV, AICEN 		 ! Used in StieglitzSnow
  REAL, PUBLIC, SAVE :: FLWALPHA     		 ! SRFLW multiplier with SRFLW < 0
  REAL, PUBLIC, SAVE :: ASTRFR, STEXP		 ! stress parameters in energy2
  REAL, PUBLIC, SAVE :: RSWILT       		 ! parameters in rsurfp
  
  LOGICAL, PUBLIC, SAVE :: LAND_FIX  		 ! Used for fixes and init changes that
  	   	   	   	     		 ! are still default in Icarus GCM

contains
  
  ! Call to get "constants" that really are variables changeable during land/landice initialization
  
  subroutine SurfParams_init(LAND_PARAMS,LSM_CHOICE, rc)
    
    implicit none
    
    CHARACTER(*), INTENT(IN) :: LAND_PARAMS
    INTEGER, INTENT(IN)      :: LSM_CHOICE
    INTEGER, OPTIONAL, INTENT(OUT) :: rc
    
    LOGICAL, SAVE :: init_called = .FALSE. ! Flag if SurfParams_init has been called
    
    ! ---------------------------------------------------------------------------
    
    if (init_called) then ! already called
       !write (*,*) "SurfParams_init being called again"
       return
    end if
    
    if (LSM_CHOICE==1) then
       
       select case (LAND_PARAMS)
       case ("Icarus")  ! "Old" LDASsa physics, current default for Icarus GCM
          LAND_FIX = .FALSE.
          CSOIL_2  = 200.
          WEMIN    = 26.
          AICEV    = 0.149
          AICEN    = 19.851
          FLWALPHA = 1.      ! i.e., FLWALPHA unchanged
          ASTRFR   = 0.333
          STEXP    = 1.
          RSWILT   = 500.
          
       case ("V24_C05")  ! V24_C05 changes, default for LDAS m4-17-0
          LAND_FIX = .TRUE.
          CSOIL_2  = 70000.  ! Post H5_0
          WEMIN    = 13.
          AICEV    = 0.107
          AICEN    = 19.893
          FLWALPHA = 0.01
          ASTRFR   = 1.
          STEXP    = 2.
          RSWILT   = 2000.
          
       case ("NRv7.2")  ! SMAP NRv7.2 changes, default for after LDAS m4-17-6
          LAND_FIX = .TRUE.
          CSOIL_2  = 70000. ! Post H5_0
          WEMIN    = 13.
          AICEV    = 0.149
          AICEN    = 19.851
          FLWALPHA = 0.04
          ASTRFR   = 0.333  ! reverted
          STEXP    = 1.     ! reverted
          RSWILT   = 500.   ! reverted
          
       case DEFAULT
          _ASSERT(.FALSE.,'LAND_PARAMS not valid or incompatible with LSM_CHOICE')
       end select
       
    else if (LSM_CHOICE==2) then
       
       select case (LAND_PARAMS)
       case ("CN_CLM40")  ! parameters to reproduce Fanwei Zeng's Catchment-CN.4.0 runs (e0004s_transientCO2_05) done with build /gpfsm/dnb31/fzeng/LDASsa_m3-16_0_p2_CatchCatchCN_for_MERRA3
          LAND_FIX        = .TRUE.
          CSOIL_2         = 70000. ! Post H5_0
          WEMIN           = 13.
          AICEV           = 0.149
          AICEN           = 19.851
          FLWALPHA        = 1.
          ASTRFR          = 0.333  ! reverted
          STEXP           = 1.     ! reverted
          RSWILT          = 1500.
          
       case DEFAULT
          _ASSERT(.FALSE.,'LAND_PARAMS not valid or incompatible with LSM_CHOICE')
       end select
       
!    else if (LSM_CHOICE==3) then
!       select case (LAND_PARAMS)      
!          
!       case ("CN_CLM45")  ! parameters to reproduce Eunjee Lee's Catchment-CN4.5 fire carbon emission simulations
!          LAND_FIX 	= .TRUE.
!          CSOIL_2  	= 70000. ! Post H5_0
!          WEMIN    	= 13.
!          AICEV    	= 0.107
!          AICEN    	= 19.893
!          FLWALPHA 	= 0.005
!          ASTRFR   	= 0.333  ! reverted
!          STEXP    	= 1.     ! reverted
!          RSWILT   	= 2000.
!          
!       case DEFAULT
!          _ASSERT(.FALSE.,'LAND_PARAMS not valid or incompatible with LSM_CHOICE')
!       end select

    else
       _ASSERT(.FALSE.,'land model choice not valid')
    end if ! LSM_CHOICE
    
    init_called = .TRUE.
    _RETURN(_SUCCESS)
    
  end subroutine SurfParams_init
  
endmodule SurfParams
