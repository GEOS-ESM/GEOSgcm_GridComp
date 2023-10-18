
module catch_constants
  
  ! reichle+koster, 2007
  ! reichle, 10 Oct 2008 - added "echo_catch_constants()"
  ! reichle, 28 Oct 2010 - moved DZ, SHR, PHI, FSN from subroutines gndtp0() and gndtmp()
  !                      - moved FWETL, FWETC from subroutine interc()
  !                      - renamed N_gndtmp -> N_gt
  ! reichle, 23 Nov 2010 - replaced PHIGT with POROS(N), ALHMGT with ALHM
  !                      - replaced constants with values from MAPL_Constants.F90
  !                        where possible
  ! reichle, 30 Nov 2010 - zero-diff revisions and clean-up for off-line (land-only) MERRA
  !                        replay capability
  !                         - restored PHIGT, ALHMGT 
  !                         - moved MIN_SNOW_MASS->MINSWE, DZ1MAX, and SATCAPFR to 
  !                            catch_constants
  !                         - moved "small" back from catch_constants() into snowrt()
  ! reichle, 14 Aug 2014 - moved constants that are only needed in subroutine catchment() 
  !                         to catchment.F90 where they are now private to the F90 module 
  !                         "catchment_model"
  !                      - removed all constants that come straight from MAPL_Constants
  !                      - renamed all remaining public constants by prefacing with "catch_*"
  !                      - moved "echo_catch_constants()" to catchment.F90 (and renamed)
  ! Sarith,  10 Nov 2015 - moved SHR, EPSILON, SCONST, CSOIL_1,  CSOIL_2, N_sm, and SATCAPFR here
  !                        from land models.
  !                      - added routing model constants  N_Pfafs_LandCatchs and N_Pfaf_Catchs
  ! Justin,  12 Apr 2018 - removed ifdef LAND_UPD, moved CSOIL_2 to SurfParams
  ! reichle, 27 Jan 2022 - cleanup: moved "public" constants from lsm_routines to here; added "CATCH_" prefix

  ! ---------------------------------------------------------------------------
  !
  
  USE MAPL_ConstantsMod, ONLY:          &
       MAPL_ALHF

  
  ! use constants from SURFPARAMS in echo_catch_contants()
  
  USE SURFPARAMS,        ONLY:                   &
       LAND_FIX, CSOIL_2, WEMIN, AICEV, AICEN,   &
       FLWALPHA, ASTRFR, STEXP, RSWILT
    
  implicit none
  private  

  public :: echo_catch_constants
  
  ! ---------------------------------------------------------------------------
  
  INTEGER, PARAMETER, PUBLIC :: CATCH_N_PFAFS       = 291284   ! # of Pfafstetter hydrological catchments (global)
  INTEGER, PARAMETER, PUBLIC :: CATCH_N_PFAFSROUTE  = 290188   ! # of catchments used in runoff routing model
  !                                                                (excl. "submerged" catchments)
  
  ! ---------------------------------------------------------------------------
  
  INTEGER, PARAMETER, PUBLIC :: CATCH_N_SNOW        = 3        ! # layers in snow model
  INTEGER, PARAMETER, PUBLIC :: CATCH_N_GT          = 6        ! # layers in ground temperature model
  INTEGER, PARAMETER, PUBLIC :: CATCH_N_ZONES       = 3        ! # hydrological regimes considered  
  
  ! ---------------------------------------------------------------------------
  !
  ! constants for use with snowrt() and snow_albedo() of module StieglitzSnow
  
  REAL,    PARAMETER, PUBLIC :: CATCH_SNWALB_RHOFS  = 150.     ! kg/m^3  
  REAL,    PARAMETER, PUBLIC :: CATCH_SNWALB_VISMAX = 0.7      !  
  REAL,    PARAMETER, PUBLIC :: CATCH_SNWALB_NIRMAX = 0.5      ! 
  REAL,    PARAMETER, PUBLIC :: CATCH_SNWALB_SLOPE  = -0.0006  ! 
  REAL,    PARAMETER, PUBLIC :: CATCH_MAXSNDEPTH    = 1.e20    ! 
  
  ! ---------------------------------------------------------------------------
  !
  ! layer depth associated with snow-free land surface soil temperatures
  !
  ! ==================
  ! Note by Randy Koster & Rolf Reichle when CSOIL=200. was still used (~2018):
  ! DZTC = .05 is a hardwired setting of the depth of the bottom of
  ! the surface soil layer.  It should be made a parameter that is tied to
  ! the heat capacity CSOIL, which had been set to either CSOIL_1 or
  ! CSOIL_2 based on vegetation type.  For now we leave
  ! it set to 0.05 despite an apparent inconsistency with CSOIL as
  ! currently used.  We do this (again, for now) because the effects of the
  ! inconsistency are drowned out by our arbitrary assumption, in computing
  ! the thermal conductivities, that the unsaturated soil has a degree of
  ! saturation of 50%.  For the flux calculation, setting the depth to .05m
  ! here provides approximately the same fluxes as setting the depth to much
  ! closer to 0 (as the value of CSOIL_2 suggests) and assuming a degree of
  ! saturation of about 25%, which is no less realistic an assumption.  There
  ! are other impacts in wet climates regarding the effect of
  ! the depth of the water table on the thermal conductivity; these impacts
  ! are presumably very small.
  ! ==================
  !
  ! DZTSURF (formerly DZTC) is the layer depth associated w/ surface soil temperatures:
  !   Catchment:    tc1, tc2, tc4  (CSOIL also includes veg canopy)
  !   CatchmentCN:  tg1, tg2, tg4  (tc[X] are separate canopy temperatures)
  
  REAL,    PARAMETER, PUBLIC :: CATCH_DZTSURF       = 0.05     ! m  layer depth for tc[X] or tg[X]

  ! ---------------------------------------------------------------------------
  !
  ! layer depths and other associated with ground heat diffusion model (gndtp0() and gndtmp())
  
  REAL,    PARAMETER, DIMENSION(CATCH_N_GT), PUBLIC :: CATCH_DZGT = &  ! m  layer depths
       (/ 0.0988, 0.1952, 0.3859, 0.7626, 1.5071, 10.0 /)
  
  ! PHIGT and ALHMGT are needed for backward compatibility with
  !  off-line (land-only) MERRA replay:
  !
  ! PHIGT = porosity used in gndtp0() and gndtmp()
  !         if neg,  POROS(n) from soil moisture submodel will be used
  !
  !               |   PHIGT      ALHMGT
  ! ------------------------------------------------
  !  MERRA        |      0.45    3.34e+5
  !  Fortuna-2_3  |  -9999.       ALHM
  
  REAL,    PARAMETER, PUBLIC :: CATCH_PHIGT         = -9999.
  REAL,    PARAMETER, PUBLIC :: CATCH_ALHMGT        = MAPL_ALHF

  REAL,    PARAMETER, PUBLIC :: CATCH_FSN           = 1.e3*CATCH_ALHMGT ! unit change J/kg/K -> J/m/K

  ! miscellaneous Catchment model constants
  
  REAL,    PARAMETER, PUBLIC :: CATCH_SHR           = 2400.    ! J/kg/K  spec heat of rock 
  !                                                              [where "per kg" is something like 
  !                                                               "per kg of water equiv. density"]
  
  REAL,    PARAMETER, PUBLIC :: CATCH_SCONST        = 1.9E6/920.   ! some snow constant
  
  REAL,    PARAMETER, PUBLIC :: CATCH_CSOIL_1       = 70000.   ! J/K - heat capacity associated w/ tsurf

  REAL,    PARAMETER, PUBLIC :: CATCH_C_CANOP       = 200.     ! J/K - heat capacity associated w/ tc (CatchCN)

  REAL,    PARAMETER, PUBLIC :: CATCH_SATCAPFR      = 0.2      ! SATCAP = SATCAPFR * LAI

  
  ! peatCLSM implementation smahanam  3-16-2021
  !
  ! Use of peat-specific hydrology (PEATCLSM) is triggered by a porosity threshold.
  ! Porosity of peat tiles depends on bcs version.
  !
  !     bcs version           | source of peat info       | porosity
  !     -----------------------------------------------------------------
  !     NLv3, NLv4            | HWSD                      | poros=0.80
  !     NLv5                  | PEATMAP                   | poros=0.93
  !
  ! - reichle, 26 Jan 2022
  
  REAL, PARAMETER, PUBLIC :: PEATCLSM_POROS_THRESHOLD  = 0.90    ! [m3/m3]

  ! max zbar for specific yield calc in PEATCLSM
  
  REAL, PARAMETER, PUBLIC :: PEATCLSM_ZBARMAX_4_SYSOIL = 0.45    ! [m] 

contains

  ! ****************************************************************************************
  
  subroutine echo_catch_constants(logunit)
    
    implicit none
    
    integer, intent(in) :: logunit
    
    write (logunit,*)
    write (logunit,*) '-----------------------------------------------------------'
    write (logunit,*)
    write (logunit,*) 'echo_catch_constants():'
    write (logunit,*)
    write (logunit,*) 'CATCH_N_PFAFS              = ', CATCH_N_PFAFS       
    write (logunit,*) 'CATCH_N_PFAFSROUTE         = ', CATCH_N_PFAFSROUTE  
    write (logunit,*) 'CATCH_N_SNOW               = ', CATCH_N_SNOW        
    write (logunit,*) 'CATCH_N_GT                 = ', CATCH_N_GT          
    write (logunit,*) 'CATCH_N_ZONES              = ', CATCH_N_ZONES       
    write (logunit,*) 'CATCH_SNWALB_RHOFS         = ', CATCH_SNWALB_RHOFS  
    write (logunit,*) 'CATCH_SNWALB_VISMAX        = ', CATCH_SNWALB_VISMAX 
    write (logunit,*) 'CATCH_SNWALB_NIRMAX        = ', CATCH_SNWALB_NIRMAX 
    write (logunit,*) 'CATCH_SNWALB_SLOPE         = ', CATCH_SNWALB_SLOPE  
    write (logunit,*) 'CATCH_MAXSNDEPTH           = ', CATCH_MAXSNDEPTH    
    write (logunit,*) 'CATCH_DZTSURF              = ', CATCH_DZTSURF          
    write (logunit,*) 'CATCH_DZGT                 = ', CATCH_DZGT          
    write (logunit,*) 'CATCH_PHIGT                = ', CATCH_PHIGT         
    write (logunit,*) 'CATCH_ALHMGT               = ', CATCH_ALHMGT        
    write (logunit,*) 'CATCH_FSN                  = ', CATCH_FSN           
    write (logunit,*) 'CATCH_SHR                  = ', CATCH_SHR           
    write (logunit,*) 'CATCH_SCONST               = ', CATCH_SCONST        
    write (logunit,*) 'CATCH_CSOIL_1              = ', CATCH_CSOIL_1       
    write (logunit,*) 'CATCH_C_CANOP              = ', CATCH_C_CANOP       
    write (logunit,*) 'CATCH_SATCAPFR             = ', CATCH_SATCAPFR
    write (logunit,*)
    write (logunit,*) 'PEATCLSM_POROS_THRESHOLD   = ', PEATCLSM_POROS_THRESHOLD
    write (logunit,*) 'PEATCLSM_ZBARMAX_4_SYSOIL  = ', PEATCLSM_ZBARMAX_4_SYSOIL
    write (logunit,*)
    write (logunit,*) 'Constants from SURFPARAMS:'
    write (logunit,*)
    write (logunit,*) 'LAND_FIX                   = ', LAND_FIX 
    write (logunit,*) 'CSOIL_2                    = ', CSOIL_2 
    write (logunit,*) 'WEMIN                      = ', WEMIN   
    write (logunit,*) 'AICEV                      = ', AICEV   
    write (logunit,*) 'AICEN                      = ', AICEN   
    write (logunit,*) 'FLWALPHA                   = ', FLWALPHA
    write (logunit,*) 'ASTRFR                     = ', ASTRFR  
    write (logunit,*) 'STEXP                      = ', STEXP   
    write (logunit,*) 'RSWILT                     = ', RSWILT  
    write (logunit,*)    
    write (logunit,*) 'end echo_catch_constants()'
    write (logunit,*)
    write (logunit,*) '-----------------------------------------------------------'
    write (logunit,*)

  end subroutine echo_catch_constants

  ! *******************************************************************************
  
end module catch_constants

! ==================== EOF ========================================================
