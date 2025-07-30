#include "MAPL_Generic.h"

MODULE IRRIGATION_MODULE

  USE MAPL
  USE ESMF
  
  IMPLICIT NONE

  ! Irrigation Module
  !
  ! First Version: March 21, 2021 (Sarith Mahanama)
  ! Second Version: June 25, 2024 (Stefano Casirati)
  !
  ! This module computes irrigation rates by 4 different methods: sprinkler, flood, furrow, and drip.
  ! The computed irrigation rates (exports) are imports to Catchment[CN], where the irrigation water 
  ! is added to water balance as follows:
  !   1) sprinkler irrigation rate added to large scale precipitation (for irrigated tiles fractions);
  !   2) drip irrigation volume    added to rootzone excess (for irrigated tiles fractions); 
  !   3) furrow irrigation volume  added to rootzone excess (for irrigated tiles fractions);
  !   4) flood irrigation volume   added to surface excess (only for paddy tiles fractions).

  ! The model uses the rootzone soil moisture state at the local start time of irrigation to compute 
  ! irrigation rates for the day and maintains the same rate throughout the irrigation duration.
  ! 
  ! Sprinkler and Flood/Furrow Irrigation methods were adapted from LIS CLSMF2.5 irrigation module 
  ! (Rodell et al., 2024 (Under Review) 
  ! Drip irrigation method calculation is similar to that of sprinkler, albeit the drip irrigation 
  ! method assumes a 10% water loss. (Source FAO)
  !
  ! (1) EXPORTS - MODEL OUTPUTS TO THE LAND MODEL (IRRIGATION RATES):
  !     1) IRRG_RATE_SPR  [kg m-2 s-1]
  !     2) IRRG_RATE_DRP  [kg m-2 s-1]
  !     3) IRRG_RATE_FRW  [kg m-2 s-1]   
  !     4) IRRG_RATE_PDY  [kg m-2 s-1]
  !     5) IRRG_RATE_TOT  [kg m-2 s-1]    (diagnostic only)
  !
  ! (2) IRRIGATED AND PADDY TILES:
  ! During land BC's generation, the fraction of irrigated crops and paddy is set to zero 
  ! if their sum is below an irrigation threshold (default 1%).
  ! Irrigated fractions can be irrigated with sprinkler, drip, and furrow,
  ! while paddy fractions can only be irrigated using the flood irrigation method.
  ! Vegetation characteristics and vegetation dynamic parameters for irrigated
  ! crops and paddy tiles were taken from the nearest grass or cropland tile.
  !
  ! (3) MODEL INPUTS:
  !     SMWP    : rootzone soil moisture content at wilting point                       [mm] 
  !     SMSAT   : rootzone soil moisture content at saturation                          [mm] 
  !     SMREF   : rootzone soil moisture is at lower tercile of RZ soil moisture range  [mm]
  !     SMCNT   : currrent rootzone soil moisture content                               [mm] 
  !     RZDEF   : rootzone moisture deficit to reach complete soil saturation for paddy [mm]
  !     LOCAL_HOUR to set irrigation switch.
  !
  ! (4) SEASONAL CYLCE OF CROP WATER DEMAND:
  !     The module provides 2 options to determine the seasonal cycle of crop water demand:
  !     4.1) IRRG_TRIGGER: 0 - SUBROUTINE irrigrate_lai_trigger
  !          The LAI-based trigger (Default and the current LIS implementation)
  !          uses precomputed minimum and maximum LAI on irrigateed pixels to determine
  !          beginning and end of crop growing seasons.
  !
  !          This LAI-based trigger is also equipped with an additional control parameter, IRRG_METHOD, 
  !          which is good to choose the method of irrigation that would run on corresponding fraction:
  !              0: (Default) All 4 methods (sprinkler/drip/furrow/paddy) concurrently, according to specified area fracs
  !              1: Only sprinkler    irrigation on *entire* tile (regardless of IRRIGFRAC or PADDYFRAC)
  !              2: Only drip         irrigation on *entire* tile (regardless of IRRIGFRAC or PADDYFRAC)
  !              3: Only furrow/flood irrigation on *entire* tile (regardless of IRRIGFRAC or PADDYFRAC)
  !
  !          IRRG_TRIGGER: 0 SPECIFIC INPUTS:
  !              IRRG_IRRIGFRAC     : fraction of tile covered by sprinkler/drip/furrow-irrigated crops;
  !                                     ranges between 0 and 1 (if IRRG_IRRIGFRAC + IRRG_PADDYFRAC > Irrigation Threshold)
  !              IRRG_PADDYFRAC     : fraction of tile covered by paddy;
  !                                     ranges between 0 and 1 (if IRRG_IRRIGFRAC + IRRG_PADDYFRAC > Irrigation Threshold)
  !
  !              Within IRRIGFRAC, allocation by irrigation method:
  !
  !                  Note: IRRG_IRRIGFRAC_SPR + IRRG_IRRIGFRAC_DRP + IRRG_IRRIGFRAC_FRW = 1.
  !
  !              IRRG_IRRIGFRAC_SPR : fraction of IRRG_IRRIGFRAC equipped for sprinkler irrigation
  !              IRRG_IRRIGFRAC_DRP : fraction of IRRG_IRRIGFRAC equipped for drip irrigation
  !              IRRG_IRRIGFRAC_FRW : fraction of IRRG_IRRIGFRAC equipped for flood/furrow irrigation
  !
  !              LAI          : time varying Leaf Area Index from the model
  !              IRRG_LAIMIN  : Minimum LAI spatially averaged over the irrigated tile fraction 
  !              IRRG_LAIMAX  : Maximum LAI spatially averaged over the irrigated tile fraction
  !
  !     4.2) IRRG_TRIGGER: 1 - SUBROUTINE irrigrate_crop_calendar
  !          Uses 26 crop calendars based on monthly crop growing areas of below crops.
  !               1    Wheat                      14    Oil palm                              
  !               2    Maize                      15    Rape seed / Canola    
  !               3    Rice                       16    Groundnuts / Peanuts  
  !               4    Barley                     17    Pulses                
  !               5    Rye                        18    Citrus                
  !               6    Millet                     19    Date palm             
  !               7    Sorghum                    20    Grapes / Vine         
  !               8    Soybeans                   21    Cotton                
  !               9    Sunflower                  22    Cocoa                 
  !              10    Potatoes                   23    Coffee                
  !              11    Cassava                    24    Others perennial      
  !              12    Sugar cane                 25    Fodder grasses        
  !              13    Sugar beet                 26    Others annual         
  !
  !          IRRG_TRIGGER: 1 SPECIFIC INPUTS:
  !              DOFYR        : day of year
  !              IRRG_TYPE    : Preferred Irrig method (NTILES, 26) -
  !                             0   CONCURRENT (default), 
  !                             1   SPRINKLER ONLY 
  !                             2   DRIP ONLY 
  !                             3   FLOOD/FURROW ONLY, and 
  !                             <0  AVOID this method 
  !              IRRG_CROPIRRIGFRAC: Crop irrigated fraction (NTILES, 26) (per Section 2, fractions have been 
  !                             adjusted such that IRRG_CROPIRRIGFRAC=1. on paddy tiles; the sum of available 
  !                             crop fractions equals 1. on irrigated crop tiles;
  !                             and is zero on non-irrigated tiles.
  !              IRRG_DOY_PLANT   : DOY start planting (NTILES, 2, 26) - up to two seasons
  !              IRRG_DOY_HARVEST : DOY end harvesting (NTILES, 2, 26) - up to two seasons
  !              If IRRG_DOY_PLANT = IRRG_DOY_HARVEST = 998, the crop is not grown on that tile 
  !
  ! (5) MODEL UPDATES (OPTIONAL INTERNALS):
  !     SRATE, DRATE, and FRATE contain irrigation rates applied on individual fractions at any given time.
  !     The second dimensions of 2D arrays is for different crop fractions i.e. the second dimension is 2 for above
  !     IRRG_TRIGGER: 0 to separately store irrigation rates in irrigated crop and paddy fractions.
  !     It would be 26 for IRRG_TRIGGER: 1.
  !     The crop calendar implemetation (IRRG_TRIGGER: 1) computes IRRG_RATE_SPR, IRRG_RATE_DRP, 
  !     IRRG_RATE_FRW, and IRRG_RATE_PDY as weighted averages of irrigation rates from
  !     all active crops in SRATE, DRATE and FRATE arrays.
  
  PRIVATE
  
  INTEGER, PARAMETER, PUBLIC :: IRRG_NCROPS = 26, IRRG_NSEASONS = 2

  type, public :: irrig_params
     
     ! Below parameters can be set via RC file.

     ! IRRGRR: Do we really want to hardwire defaults here *and* in GEOS_SurfaceGridComp.rc ???

     REAL :: irrig_thres      =   0.01 ! threshold of tile fraction to turn the irrigation model on. 
     REAL :: lai_thres        =   0.6  ! threshold of LAI range to turn irrigation on
     REAL :: efcor            =  25.0  ! efficiency correction (% water loss: efcor = 0% denotes 100% efficient water use)
     REAL :: MIDS_LNGTH       =   0.6  ! Mid-season length as a fraction of crop growing season length (to be used with IRRG_TRIGGER: 1)
     
     ! Sprinkler parameters
     ! --------------------
     REAL :: sprinkler_stime  =   6.0  ! sprinkler irrigation start time (local)                                         [hours]
     REAL :: sprinkler_dur    =   4.0  ! sprinkler irrigation duration                                                   [hours]
     REAL :: sprinkler_thres  =   0.7  ! threshold for soil moisture availability ("MA") to trigger sprinkler irrigation [dim-less]
     
     ! Drip parameters 
     ! ---------------
     REAL :: drip_stime       =   8.0  ! drip irrigation start time (local)                                              [hours] 
     REAL :: drip_dur         =   8.0  ! drip irrigation duration                                                        [hours]
     REAL :: drip_thres       =   0.7  ! threshold for soil moisture availability ("MA") to trigger drip      irrigation [dim-less]
     
     ! Flood parameters
     ! ----------------
     REAL :: flood_stime      =   6.0  ! flood irrigation start time (local)                                             [hours]
     REAL :: flood_dur        =   8.0  ! flood irrigation duration                                                       [hours]
     REAL :: flood_thres      =   0.6  ! threshold for soil moisture availability ("MA") to trigger flood     irrigation [dim-less]


     
  end type irrig_params
  
  type, public, extends (irrig_params) :: irrigation_model
     
   contains
     
     ! public
     procedure, public :: init_model
     generic,   public :: run_model     => irrigrate_lai_trigger, irrigrate_crop_calendar
     generic,   public :: update_irates => update_irates_lai,     update_irates_ccalendar
     
     ! private
     procedure, private :: irrigrate_lai_trigger
     procedure, private :: irrigrate_crop_calendar
     procedure, private :: cwd => crop_water_deficit
     procedure, private :: irrig_by_method
     procedure, private :: update_irates_lai
     procedure, private :: update_irates_ccalendar
     
  end type irrigation_model

contains

  ! ----------------------------------------------------------------------------
  
  SUBROUTINE init_model (IP, SURFRC)

    implicit none
    class (irrigation_model), intent(inout) :: IP
    type(irrig_params)                      :: DP
    CHARACTER(*), INTENT(IN)                :: SURFRC
    type(ESMF_Config)                       :: SCF
    integer                                 :: status, RC
    character(len=ESMF_MAXSTR)              :: Iam

    Iam='IRRIGATION_MODULE: init_model'
    
    SCF = ESMF_ConfigCreate(__RC__) 

    CALL ESMF_ConfigLoadFile     (SCF,SURFRC,rc=status) ; VERIFY_(STATUS)

    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_SPR_STIME:' , VALUE=IP%sprinkler_stime, DEFAULT=DP%sprinkler_stime, __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_SPR_DUR:'   , VALUE=IP%sprinkler_dur,   DEFAULT=DP%sprinkler_dur  , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_SPR_THRES:' , VALUE=IP%sprinkler_thres, DEFAULT=DP%sprinkler_thres, __RC__ )

    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_DRP_STIME:' , VALUE=IP%drip_stime,      DEFAULT=DP%drip_stime     , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_DRP_DUR:'   , VALUE=IP%drip_dur,        DEFAULT=DP%drip_dur       , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_DRP_THRES:' , VALUE=IP%drip_thres,      DEFAULT=DP%drip_thres     , __RC__ )

    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_FLD_STIME:' , VALUE=IP%flood_stime,     DEFAULT=DP%flood_stime    , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_FLD_DUR:'   , VALUE=IP%flood_dur,       DEFAULT=DP%flood_dur      , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_FLD_THRES:' , VALUE=IP%flood_thres,     DEFAULT=DP%flood_thres    , __RC__ )

    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_EFCOR:'     , VALUE=IP%efcor,           DEFAULT=DP%efcor          , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_LAI_THRES:' , VALUE=IP%lai_thres,       DEFAULT=DP%lai_thres      , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_MIDS_LNGTH:', VALUE=IP%MIDS_LNGTH,      DEFAULT=DP%MIDS_LNGTH     , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='IRRG_FRAC_THRES:', VALUE=IP%irrig_thres,     DEFAULT=DP%irrig_thres    , __RC__ )
    CALL ESMF_ConfigDestroy      (SCF, __RC__)

  END SUBROUTINE init_model

  ! ----------------------------------------------------------------------------

  SUBROUTINE irrigrate_lai_trigger (this,IRRG_METHOD, local_hour,                                        &
            IRRG_IRRIGFRAC, IRRG_PADDYFRAC, IRRG_IRRIGFRAC_SPR, IRRG_IRRIGFRAC_DRP, IRRG_IRRIGFRAC_FRW,  &           
            SMWP, SMSAT, SMREF, SMCNT, LAI, IRRG_LAIMIN,IRRG_LAIMAX, RZDEF,                              &
            IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY,                                  &
            SRATE, DRATE, FRATE)

    implicit none
    class (irrigation_model), intent(inout) :: this
    integer, intent (in)                    :: IRRG_METHOD
    real, dimension (:), intent (in)        :: local_hour
    real, dimension (:), intent (in)        :: IRRG_IRRIGFRAC, IRRG_PADDYFRAC, IRRG_IRRIGFRAC_SPR, &
         IRRG_IRRIGFRAC_DRP, IRRG_IRRIGFRAC_FRW, SMWP, SMSAT, SMREF, SMCNT, LAI, IRRG_LAIMIN, IRRG_LAIMAX, RZDEF
    real, dimension (:), intent (inout)     :: IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY
    real, dimension (:,:),intent (inout)    :: SRATE, DRATE, FRATE
    INTEGER                                 :: NTILES, N, crop
    REAL                                    :: ma, H1, H2, HC, IT, ROOTFRAC, LAITHRES
    logical                                 :: season_end
    
    NTILES = SIZE (IRRG_IRRIGFRAC)
    TILE_LOOP : DO N = 1, NTILES
       IF(IRRG_LAIMAX (N) > IRRG_LAIMIN (N)) THEN
          LAITHRES = IRRG_LAIMIN (N) + this%lai_thres * (IRRG_LAIMAX (N) - IRRG_LAIMIN (N))          
          ROOTFRAC = MIN((LAI(N) - IRRG_LAIMIN (N)) / (IRRG_LAIMAX(N) - IRRG_LAIMIN(N)) ,1.0)          
       ELSE
          ROOTFRAC = 0.
       ENDIF
       HC = local_hour(n)

       season_end        =  .true.
       
       CHECK_LAITHRES : IF (LAI(N) >= LAITHRES) THEN
          season_end = .false.
          CHECK_IRRIGFRACS: IF ((IRRG_IRRIGFRAC(N) > 0.).OR.(IRRG_PADDYFRAC(N)>0.)) THEN

             !-----------------------------------------------------------------------------
             !     Get the rootzone moisture availability to the plant
             !-----------------------------------------------------------------------------
             if (IRRG_IRRIGFRAC(N) > 0.) then
                if(SMREF(N) > SMWP(N))then
                        ma = (SMCNT(N) - SMWP(N)) /(SMREF(N) - SMWP(N))
                else
                        ma = -1.
                endif
             
                if(ma >= 0) then
                                
                        SELECT CASE (IRRG_METHOD)                
                        CASE (0)  ! CONCURRENTLY SPRINKER + FLOOD + FURROW + DRIP on corresponding fractions

                        call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), &
                                SRATE = SRATE (N,1), &
                                DRATE = DRATE (N,1), &
                                FRATE = FRATE (N,1))

                        SRATE (N,1) =  SRATE (N,1)*IRRG_IRRIGFRAC_SPR(N)
                        DRATE (N,1) =  DRATE (N,1)*IRRG_IRRIGFRAC_DRP (N)
                        FRATE (N,1) =  FRATE (N,1)*IRRG_IRRIGFRAC_FRW (N) 
                   
                        CASE (1)  ! SPRINKLER only

                        call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), &
                                SRATE = SRATE (N,1))
                   
                        DRATE (N,1) = 0.
                        FRATE (N,1) = 0.
                                      
                        CASE (2)  ! DRIP only

                        call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), &
                                DRATE = DRATE (N,1))
                   
                        SRATE (N,1) = 0.
                        FRATE (N,1) = 0.

                        CASE (3)  ! FLOOD only

                        call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), &
                                FRATE = FRATE (N,1))

                        SRATE (N,1) = 0.
                        DRATE (N,1) = 0.
                   
                        CASE DEFAULT
                        PRINT *, 'irrigrate_lai_trigger: IRRG_METHOD can be 0,1,2, or3'
                        CALL EXIT(1)
                        END SELECT
                endif
             endif
                     
             if (IRRG_PADDYFRAC (N) > 0.) then
             
                H1 = this%flood_stime
                H2 = this%flood_stime + this%flood_dur
                if ((HC >= H1).AND.(HC < H2)) then
                        ! use RZDEF at H1 during H1 <= HC < H2 to compute irrigrate for paddy. 
                        if(H1 == HC) FRATE (N,2) = RZDEF(N) *0.1/(H2 - H1)/ 3600.
                        else
                        FRATE (N,2) = 0.
                endif
                SRATE (N,2) = 0.
                DRATE (N,2) = 0.
             endif 
          
          ELSE
             
             SRATE (N,:) = 0.
             DRATE (N,:) = 0.
             FRATE (N,:) = 0.
             
          ENDIF CHECK_IRRIGFRACS
       ENDIF CHECK_LAITHRES

       ! turn off irrigation if LAI is smaller than the LAI trigger marking end of the season
       if(season_end) then
          DO crop = 1, 2                         ! With LAI trigger Crop 1 = Irrigated Fraction, Crop 2 = Paddy Fraction 
             SRATE (N,crop) = 0.
             DRATE (N,crop) = 0.
             FRATE (N,crop) = 0.
          END DO
       endif
       
    END DO TILE_LOOP

    ! Update IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY EXPORTS to be sent to land models.
        
    call this%update_irates( IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY, &
         IRRG_IRRIGFRAC, IRRG_PADDYFRAC, SRATE, DRATE, FRATE )
    
  END SUBROUTINE irrigrate_lai_trigger

  ! ----------------------------------------------------------------------------

  SUBROUTINE irrigrate_crop_calendar(this,dofyr,local_hour, &
       IRRG_IRRIGFRAC_SPR, IRRG_IRRIGFRAC_DRP, IRRG_IRRIGFRAC_FRW,                        &
       IRRG_CROPIRRIGFRAC,IRRG_DOY_PLANT, IRRG_DOY_HARVEST, IRRG_TYPE ,  &
       SMWP,SMSAT,SMREF,SMCNT, RZDEF,                       &  
       IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY, SRATE, DRATE, FRATE)

    implicit none
    class(irrigation_model),intent(inout):: this
    integer, intent (in)                 :: dofyr
    real, dimension (:),   intent (in)   :: local_hour, IRRG_IRRIGFRAC_SPR, IRRG_IRRIGFRAC_DRP, IRRG_IRRIGFRAC_FRW
    real, dimension (:),   intent (in)   :: SMWP, SMSAT, SMREF, SMCNT, RZDEF
    real, dimension(:,:),  intent (in)   :: IRRG_CROPIRRIGFRAC                                             ! IRRG_NCROPS
    real, dimension(:,:),  intent (in)   :: IRRG_TYPE                                                      ! IRRG_NCROPS
    real, dimension(:,:,:),intent (in)   :: IRRG_DOY_PLANT                                                 ! IRRG_NSEASONS, IRRG_NCROPS
    real, dimension(:,:,:),intent (in)   :: IRRG_DOY_HARVEST                                               ! IRRG_NSEASONS, IRRG_NCROPS
    real, dimension (:),intent (inout)   :: IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY
    real, dimension (:,:),intent (inout) :: SRATE, DRATE, FRATE
    INTEGER                              :: NTILES, N, crop, sea, ITYPE, I
    REAL                                 :: ma, H1, H2, HC, IT, ROOTFRAC, void_frac
    logical                              :: season_end (IRRG_NCROPS)
    NTILES = SIZE (local_hour)
         
    TILE_LOOP : DO N = 1, NTILES
       HC = local_hour(n)
       IF_IRR: if(SUM(IRRG_CROPIRRIGFRAC(N,:)) > 0.) then
          ! the tile is irrigated crop or paddy
          season_end =  .true.
          CROP_LOOP: DO crop = 1, IRRG_NCROPS
             CROP_IN_TILE: if(IRRG_CROPIRRIGFRAC(N,crop) > 0.) then
                ! crop is grown in this tile
                TWO_SEASONS: do sea = 1, IRRG_NSEASONS
                   IS_CROP: IF(IRRG_DOY_PLANT(N, sea, crop) /= 998) THEN
                      ! crop is grown in sea
                      IS_SEASON: IF(IS_WITHIN_SEASON(dofyr,NINT(IRRG_DOY_PLANT(N, sea, crop)),NINT(IRRG_DOY_HARVEST(N, sea, crop)))) THEN
                         ! dofyr falls within the crop season
                         season_end(crop) = .false.
                         PADDY_OR_CROP: if (crop == 3) then                           
                            ! PADDY TILE
                            H1 = this%flood_stime
                            H2 = this%flood_stime + this%flood_dur
                            if ((HC >= H1).AND.(HC < H2)) then
                               ! use RZDEF at H1 during H1 <= HC < H2 to compute irrigrate. 
                               if(H1 == HC) FRATE (N,crop) = RZDEF(N) /(H2 - H1)/ 3600.          
                            else
                               FRATE (N,crop) = 0.
                            endif
                            SRATE (N,crop) = 0.
                            DRATE (N,crop) = 0.

                         else
                            
                            ! IRRIGATED CROP: compute sum of irrigrates from 25 crops.
                            
                            ROOTFRAC = CROP_SEASON_STAGE (this%MIDS_LNGTH, dofyr,NINT(IRRG_DOY_PLANT(N, sea, crop)),NINT(IRRG_DOY_HARVEST(N, sea, crop)))
                            if(SMREF(N) > SMWP(N))then
                               ma = (SMCNT(N) - SMWP(N)) /(SMREF(N) - SMWP(N))
                            else
                               ma = -1.
                            endif
                            
                            SOILM: if(ma >= 0) then
                               
                               ITYPE = NINT(IRRG_TYPE(N,crop))

                               CROP_IMETHOD: if (ITYPE == 0) then
                                  
                                  ! concurrently on sprinkler, drip and flood fractions
                                  call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), &
                                       SRATE = SRATE (N,crop), &
                                       DRATE = DRATE (N,crop), &
                                       FRATE = FRATE (N,crop))
                                  
                                  SRATE (N,crop) =  SRATE (N,crop)*IRRG_IRRIGFRAC_SPR(N)
                                  DRATE (N,crop) =  DRATE (N,crop)*IRRG_IRRIGFRAC_DRP (N)
                                  FRATE (N,crop) =  FRATE (N,crop)*IRRG_IRRIGFRAC_FRW (N)

                               elseif (ITYPE > 0) then
                                  
                                  ! only this method
                                  if (ITYPE == 1) call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), SRATE = SRATE (N,crop))
                                  if (ITYPE == 2) call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), DRATE = DRATE (N,crop))
                                  if (ITYPE == 3) call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), FRATE = FRATE (N,crop))

                               elseif (ITYPE < 0) then
                                  
                                  ! crop does not use IRRG_METHOD -(ITYPE)
                                  void_frac = 0.                                  
                                  DO I = 1,3
                                     if(I == ABS(ITYPE))then
                                        ! this itype isn't used by this crop other 2 fractions equally share this fraction
                                        if (I == 1) then
                                           void_frac = IRRG_IRRIGFRAC_SPR(N)/2.
                                           SRATE(N,crop) = 0.
                                        endif
                                        if (I == 2) then
                                           void_frac = IRRG_IRRIGFRAC_DRP (N)/2.
                                           DRATE(N,crop) = 0.
                                        endif
                                        if (I == 3)then
                                           void_frac = IRRG_IRRIGFRAC_FRW (N)/2.
                                           FRATE(N,crop) = 0.
                                        endif
                                     else
                                        if (I == 1) call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), SRATE = SRATE (N,crop))
                                        if (I == 2) call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), DRATE = DRATE (N,crop))
                                        if (I == 3) call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), FRATE = FRATE (N,crop))
                                     endif
                                  END DO
                                  DO I = 1,3
                                     if(I /= ABS(ITYPE))then
                                        if (I == 1) SRATE (N,crop) =  SRATE (N,crop)*(IRRG_IRRIGFRAC_SPR(N) + void_frac)
                                        if (I == 2) DRATE (N,crop) =  DRATE (N,crop)*(IRRG_IRRIGFRAC_DRP (N)     + void_frac)
                                        if (I == 3) FRATE (N,crop) =  FRATE (N,crop)*(IRRG_IRRIGFRAC_FRW (N)    + void_frac)
                                     endif
                                  ENDDO

                               endif CROP_IMETHOD
                            endif SOILM
                         ENDIF PADDY_OR_CROP
                      ENDIF IS_SEASON
                   end IF IS_CROP
                end do TWO_SEASONS
             endif CROP_IN_TILE
          END DO CROP_LOOP

          ! turn off irrigation for crops that ended the season
          DO crop = 1, IRRG_NCROPS
             if(season_end(crop)) then
                SRATE (N,crop) = 0.
                DRATE (N,crop) = 0.
                FRATE (N,crop) = 0.
             endif
          END DO
                             
       endif IF_IRR
    END DO TILE_LOOP
    
    ! Update IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY EXPORTS to be sent to land models
    ! They are weighted averaged over 26 crop fractions.

    call this%update_irates( IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY, &
         IRRG_CROPIRRIGFRAC, SRATE, DRATE, FRATE )
    
  END SUBROUTINE irrigrate_crop_calendar

  ! ----------------------------------------------------------------------------

  SUBROUTINE update_irates_lai (this,IRRG_RATE_SPR,IRRG_RATE_DRP,IRRG_RATE_FRW,IRRG_RATE_PDY, &
       IRRG_IRRIGFRAC,IRRG_PADDYFRAC,SRATE,DRATE,FRATE)
    
    implicit none

    class(irrigation_model),  intent(inout) :: this

    real,    dimension(:),    intent(in)    :: IRRG_IRRIGFRAC, IRRG_PADDYFRAC
    real,    dimension(:,:),  intent(in)    :: SRATE, DRATE, FRATE
    real,    dimension(:),    intent(inout) :: IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY

    integer                                 :: N, NT

    ! INITIALIZE EXPORTS
    IRRG_RATE_SPR = 0.
    IRRG_RATE_DRP = 0.
    IRRG_RATE_FRW = 0.
    IRRG_RATE_PDY = 0.

    NT = size (IRRG_IRRIGFRAC)

    !_ASSERT(size (SRATE,2)==IRRG_NCROPS,'Irrigation model LAI trigger irrig tile types mismatch')
    
    DO N = 1, NT
       IF ((IRRG_IRRIGFRAC(N) + IRRG_PADDYFRAC(N))  > 0.) THEN
          IRRG_RATE_SPR (N) = SRATE (N,1)
          IRRG_RATE_DRP (N) = DRATE (N,1)
          IRRG_RATE_FRW (N) = FRATE (N,1)
          IRRG_RATE_PDY (N) = FRATE (N,2) 
       ENDIF
    END DO
    
  END SUBROUTINE update_irates_lai

  !...............................................................................
  
  SUBROUTINE update_irates_ccalendar(this,IRRG_RATE_SPR,IRRG_RATE_DRP,IRRG_RATE_FRW,IRRG_RATE_PDY, &
       IRRG_CROPIRRIGFRAC,SRATE,DRATE,FRATE)

    implicit none
    
    class(irrigation_model), intent(inout)  :: this

    real,    dimension(:,:), intent (in)    :: IRRG_CROPIRRIGFRAC ! IRRG_NCROPS
    real,    dimension(:,:), intent (in)    :: SRATE, DRATE, FRATE
    real,    dimension(:),   intent (inout) :: IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY

    integer                                 :: N, NT, crop

    ! INITIALIZE EXPORTS
    IRRG_RATE_SPR = 0.
    IRRG_RATE_DRP = 0.
    IRRG_RATE_FRW = 0.
    IRRG_RATE_PDY = 0.

    !_ASSERT(size (SRATE,2)==IRRG_NCROPS,'Irrigation model crop calendar trigger IRRG_NCROPS mismatch')
    NT =  size (IRRG_RATE_SPR)
    DO N = 1, NT
       if(SUM(IRRG_CROPIRRIGFRAC(N,:)) > 0.) then
          DO crop = 1, IRRG_NCROPS
             IRRG_RATE_SPR(N)    = IRRG_RATE_SPR(N) + SRATE (N,crop)*IRRG_CROPIRRIGFRAC(N,crop)/SUM(IRRG_CROPIRRIGFRAC(N,:))
             IRRG_RATE_DRP(N)    = IRRG_RATE_DRP(N) + DRATE (N,crop)*IRRG_CROPIRRIGFRAC(N,crop)/SUM(IRRG_CROPIRRIGFRAC(N,:))
             if (crop==3) then
                ! If crop is rice (crop ==3) then use flood irrigation. Otherwise use furrow irrigation.
                IRRG_RATE_PDY(N) = IRRG_RATE_PDY(N) + FRATE (N,crop)*IRRG_CROPIRRIGFRAC(N,crop)/SUM(IRRG_CROPIRRIGFRAC(N,:))
             else 
                IRRG_RATE_FRW(N) = IRRG_RATE_FRW(N) + FRATE (N,crop)*IRRG_CROPIRRIGFRAC(N,crop)/SUM(IRRG_CROPIRRIGFRAC(N,:))
             endif
          END DO
       endif
    END DO

  END SUBROUTINE update_irates_ccalendar

  ! ----------------------------------------------------------------------------

  SUBROUTINE irrig_by_method (this, HC, ma, ROOTFRAC, SMCNT, SMREF, SRATE, DRATE, FRATE)

    implicit none
    class (irrigation_model), intent(inout) :: this
    REAL, intent (in)                       :: HC, ma, ROOTFRAC,SMCNT, SMREF
    REAL, optional, intent (inout)          :: SRATE, DRATE, FRATE 
    REAL                                    :: H1, H2, IT

    if(present (SRATE)) then
       ! SPRINKLER IRRIGATION
       H1 = this%sprinkler_stime
       H2 = this%sprinkler_stime + this%sprinkler_dur
       IT = this%sprinkler_thres
       if ((HC >= H1).AND.(HC < H2)) then
          ! The model uses rootzone soil moisture state at H1 to compute irrigation
          ! rates for the day and maintains the same rate through out the irrigation
          ! duration (H1 <= HC < H2).
          if((ma <= IT).AND.(H1 == HC)) &
               SRATE = this%cwd (ROOTFRAC,SMCNT,SMREF,this%efcor)/(H2 - H1)/3600.
       else
          SRATE = 0.
       endif
    endif

    if(present (DRATE)) then
       ! DRIP IRRIGATION
       H1 = this%drip_stime
       H2 = this%drip_stime + this%drip_dur
       IT = this%drip_thres 
       if ((HC >= H1).AND.(HC < H2)) then
          ! use SMCNT at H1 during H1 <= HC < H2 to compute irrigrate.
          ! Notice drip uses the same soil moisture threshold of sprinkler but with 10.% efficiency correction.
          if((ma <= IT).AND.(H1 == HC)) &
               DRATE = this%cwd(ROOTFRAC,SMCNT,SMREF,10.)/(H2 - H1)/3600.
       else
          DRATE = 0.
       endif
    endif

    if(present (FRATE)) then
       ! FURROW IRRIGATION
       H1 = this%flood_stime
       H2 = this%flood_stime + this%flood_dur
       IT = this%flood_thres
       if ((HC >= H1).AND.(HC < H2)) then
          ! use SMCNT at H1 during H1 <= HC < H2 to compute irrigrate.
          ! Notice Furrow irrigation uses the same soil moisture threshold of sprinkler but with 
          ! the efficiency correction increased by 15 (e.g., Field application efficiency Sprinkler 75%, Surface Irrigation 60%.
          ! Source FAO)
          if((ma <= IT).AND.(H1 == HC)) &
               FRATE = this%cwd (ROOTFRAC,SMCNT,SMREF,this%efcor+15.)/(H2 - H1)/3600.
       else
          FRATE = 0.
       endif
    endif
    
  END SUBROUTINE irrig_by_method
  
  ! ----------------------------------------------------------------------------

  REAL FUNCTION crop_water_deficit (this, rootfrac, asmc, smcref, efcor)

    implicit none
    class(irrigation_model),intent(inout):: this
    real, intent (in)                    :: rootfrac, asmc, smcref, efcor
    
    crop_water_deficit = 0.
    if(smcref > asmc) crop_water_deficit = rootfrac*(smcref - asmc)*100.0/(100.0-efcor)   
    
  END FUNCTION crop_water_deficit
  
  ! ----------------------------------------------------------------------------

  logical FUNCTION IS_WITHIN_SEASON (DOY,DP, DH)

    implicit none
    integer, intent (in) :: DOY,DP, DH

    IS_WITHIN_SEASON = .false.
    if(DH > DP) then
       if((DOY >= DP).AND.(DOY <=  DH)) IS_WITHIN_SEASON = .true.
    elseif (DH < DP) then
       if((DOY >= DP).AND.(DOY <= 366)) IS_WITHIN_SEASON = .true.
       if((DOY >=  1).AND.(DOY <=  DH)) IS_WITHIN_SEASON = .true.
    endif
      
  end FUNCTION IS_WITHIN_SEASON

  ! ----------------------------------------------------------------------------
  
  real FUNCTION CROP_SEASON_STAGE (MSL, DOY,DP, DH)
    
    ! MSL : mid season length [-] as a fraction of the length of the season
    ! DOY : doy of year
    ! DP  : plant date
    ! DH  : harvest date
    !                                MSL
    !               1.0 ___________________________
    !                      /|                     |\
    !                     / |                     | \
    !                    /  |                     |  \
    !                   /   |                     |   \
    !          --------------------- DOY ---------------------------->
    !                  t0  t1                   t2   t3
    !                  DP                             DH  
    !                  |<----       SEAL            -->|
   
    implicit none
    real, intent    (in) :: MSL
    integer, intent (in) :: DOY,DP, DH
    real                 :: seal, t0, t1, t2, t3, CTime
    
    CROP_SEASON_STAGE = 0.
    
    if(DH > DP) then
       seal  = real (DH - DP + 1)
       CTime = real (DOY)
    else
       ! the crop season is fall-to-spring
       seal = real(366 - DP + 1 + DH)      
       if((DOY >= DP).AND.(DOY <= 366)) CTIME = real (DOY) 
       if((DOY >=  1).AND.(DOY <=  DH)) CTIME = real (DOY) + 365.
    endif

    t0 = real (DP)
    t1 = t0 + seal * (1. - MSL)/2. ! assumes equal development and late periods
    t2 = t1 + seal*MSL
    t3 = t0 + seal

    if (ctime < t1)                     CROP_SEASON_STAGE = (CTIME -t0)/(t1 - t0)
    if ((t1 <= ctime).and.(ctime < t2)) CROP_SEASON_STAGE = 1.
    if (ctime >= t2)                    CROP_SEASON_STAGE = (t3 - ctime)/(t3 - t2)  

    if(CROP_SEASON_STAGE > 1.) CROP_SEASON_STAGE = 1.
    
  end FUNCTION CROP_SEASON_STAGE
  
END MODULE IRRIGATION_MODULE
