
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
  ! Sarith, 10 Nov 2015  - moved SHR, EPSILON, SCONST, CSOIL_1,  CSOIL_2, N_sm, and SATCAPFR here
  !                        from land models.
  !                      - added routing model constants  N_Pfafs_LandCatchs and N_Pfaf_Catchs
  ! Justin, 12 Apr 2018  - removed ifdef LAND_UPD, moved CSOIL_2 to SurfParams

  USE MAPL_ConstantsMod, ONLY:          &
       PIE               => MAPL_PI,     &  ! -                       
       ALHE              => MAPL_ALHL,   &  ! J/kg  @15C              
       ALHM              => MAPL_ALHF,   &  ! J/kg                    
       ALHS              => MAPL_ALHS,   &  ! J/kg                    
       TF                => MAPL_TICE,   &  ! K                       
       RGAS              => MAPL_RGAS,   &  ! J/(kg K)                
       SHW               => MAPL_CAPWTR, &  ! J/kg/K  spec heat of wat
       SHI               => MAPL_CAPICE, &  ! J/kg/K  spec heat of ice
       EPSILON           => MAPL_EPSILON
 
  implicit none
  private  

  ! ---------------------------------------------------------------------------

  INTEGER, PARAMETER, PUBLIC :: N_Pfaf_Catchs      = 291284 ! # of Pfafstetter hydrological catchements in the globe
  INTEGER, PARAMETER, PUBLIC :: N_Pfaf_LandCatchs  = 290188 ! # of Pfafstetter hydrological catchments used within
                                                        ! the runoff routing model (excluding submerged catchments)

  ! ---------------------------------------------------------------------------

  INTEGER, PARAMETER, PUBLIC :: CATCH_N_SNOW   = 3      ! # layers in snow model
  INTEGER, PARAMETER, PUBLIC :: CATCH_N_GT     = 6      ! # layers in ground temperature model
  INTEGER, PARAMETER, PUBLIC :: N_sm           = 3      ! # hydrological regimes considered  

  ! ---------------------------------------------------------------------------
  !
  ! constants for use with snowrt() and snow_albedo() of module StieglitzSnow
  
  REAL,    PARAMETER, PUBLIC :: CATCH_SNWALB_RHOFS  = 150.    ! kg/m^3  
  REAL,    PARAMETER, PUBLIC :: CATCH_SNWALB_VISMAX = 0.7     !  
  REAL,    PARAMETER, PUBLIC :: CATCH_SNWALB_NIRMAX = 0.5     ! 
  REAL,    PARAMETER, PUBLIC :: CATCH_SNWALB_SLOPE  = -0.0006 ! 
  REAL,    PARAMETER, PUBLIC :: CATCH_MAXSNDEPTH    = 1.e20   ! 
  REAL,    PARAMETER, PUBLIC :: CATCH_DZ1MAX        = 0.08    ! m   
  
  ! ---------------------------------------------------------------------------
  ! moved back from CLSM and CLSM-CN 
  ! - Sarith, 10 Nov 2015

  REAL,    PARAMETER, PUBLIC :: SHR      = 2400.  ! J/kg/K  spec heat of rock 
  !                                                     [where "per kg" is something like 
  !                                                      "per kg of water equiv. density"]
  REAL,    PARAMETER, PUBLIC :: SCONST   = 1.9E6/920.
  ! Changed CSOIL_2 back to pre-MERRA value of 70,000 J/K.
  ! This amounts to a change in the very definition of the surface temperature and
  !  also shifts the layers for ground heat content (soil temperature) downward. 
  ! - reichle, 16 Nov 2015
  REAL,    PARAMETER, PUBLIC :: CSOIL_1  = 70000. ! J/K - heat capacity associated w/ tsurf
! #ifdef LAND_UPD
!   REAL,    PARAMETER, PUBLIC :: CSOIL_2  = 70000.   ! J/K - heat capacity associated w/ tsurf ! Post H5_0  
! #else
!   REAL,    PARAMETER, PUBLIC :: CSOIL_2  =   200.   ! J/K - heat capacity associated w/ tsurf 
! #endif
  REAL,    PARAMETER, PUBLIC :: C_CANOP  = 200.   ! J/K - heat capacity associated w/ tc
  REAL,    PARAMETER, PUBLIC :: SATCAPFR = 0.2    ! SATCAP = SATCAPFR * LAI



end module catch_constants
