#include "MAPL_Generic.h"

MODULE IRRIGATION_MODULE

  USE MAPL
  USE ESMF
  
  IMPLICIT NONE

  ! This module computes irrigation rates by 3 different methods: sprinkler, flood and drip.
  ! Computed irrigation rates return to the land model as rates that water is added to
  ! the hydrological cycle by irrigation. Subsequently, land models add irrigation feedback: 
  !      sprinkler irrigation rate to large scale precipitation;
  !      drip irrigation volume to rootzone excess; and
  !      flood irrigation volume to rootzone excess.
  ! The model uses rootzone soil moisture state at the local start time of irrigation to compute 
  ! irrigation rates for the day and maintains the same rate through out the irrigation duration.
  ! 
  ! Sprinkler and Flood Irrigation methods were adapted from LIS CLSMF2.5 irrigatrion module:
  ! https://github.com/NASA-LIS/LISF/blob/master/lis/surfacemodels/land/clsm.f2.5/irrigation/clsmf25_getirrigationstates.F90 
  ! Drip irrigation method calculation is similar to that of sprinkler, albeit the drip irrigation method assumes a 0% water loss.
  !
  ! March 21, 2021 (Sarith Mahanama) - First Version
  !
  ! (1) EXPORTS - MODEL OUTPUTS TO THE LAND MODEL (IRRIGATION RATES):
  !    1) SPRINKLERRATE [kg m-2 s-1]
  !    2) DRIPRATE [kg m-2 s-1]
  !    3) FLOODRATE [kg m-2 s-1]    
  !
  ! (2) COMPUTATIONAL TILES:
  !    During land BCs generation, land tiles where irrigated crops or paddy is present were further subdivided to a
  !    mosaic of upto 3 fractions: i) non-irrigated land, ii) irrigated crop, or iii) paddy.
  !    The model treats each fraction as a separate computational tile and runs on each individual fraction with own parameters and prognostics.
  !    All fractions inherited model and soil parameters from the main land tile that they belong to. A special treatment of setting BF3 to a high
  !    value (25.) was applied to paddy/crop tiles to account for the uniquely flat nature of farmlands. Vegetation characteristics and vegetation dynamic
  !    parameters for irrigated crop and paddy tiles were taken from the nearest grass or crops land tile. 
  !    During tiling and BCs data preparation, computed fractional coverages for land tiles were also adjusted
  !    to reflect each computational tile under the land grid component represents entirely one of the 3 irrigation surface types: a non-irrigated land,
  !    OR a irrigated-crops OR a paddy tile.
  !
  ! (3) MODEL INPUTS:
  !     SMWP    : rootzone soil moisture content at wilting point [mm] 
  !     SMSAT   : rootzone soil moisture content at saturation [mm] 
  !     SMREF   : rootzone soil moisture is at lower tercile of RZ soil moisture range [mm]
  !     SMCNT   : currrent root zone soil moisture content [mm] 
  !     RZDEF   : rootzone moisture deficit to reach complete soil saturation for paddy [mm]
  !     LOCAL_HOUR to set irrigation switch.
  !
  ! (4) SEASONAL CYLCE OF CROP WATER DEMAND:
  ! The module provides 2 options to determine the seasonal cycle of crop water demand:
  !    4.1) IRRIG_TRIGGER: 0 - SUBROUTINE irrigrate_lai_trigger
  !          The LAI-based trigger (Default and the current LIS implementation)
  !          uses precomputed minimum and maximum LAI on irrigateed pixels to determine
  !          beginning and end of crop growing seasons.
  !
  !          This LAI-based trigger is also equipped with an additional control parameter, IRRIG_METHOD, which is good to choose the method of irrigation
  !          that woould run on corresponding fractions       
  !          i)  0: (Default) All 3 methods (sprinkler/flood/drip) concurrently.
  !          ii) 1: Sprinkler irrigation on entire tile.
  !          iv) 2: Drip irrigation on entire tile.
  !          iii)3: Flood irrigation on entire tile.
  !
  !     IRRIG_TRIGGER: 0 SPECIFIC INPUTS:
  !          IRRIGFRAC    : fraction of tile covered by irrigated crops (per Section 2, values will be 0. or 1.)
  !          PADDYFRAC    : fraction of tile covered by paddy (per Section 2, values will be 0. or 1.)
  !          SPRINKLERFR  : fraction of tile equipped for sprinkler irrigation
  !          DRIPFR       : fraction of tile equipped for drip irrigation
  !          FLOODFR      : fraction of tile equipped for flood irrigation
  !          LAI          : time varying Leaf Area Index from the model
  !          LAIMIN       : Minimum LAI spatially averaged over the irrigated tile fraction 
  !          LAIMAX       : Maximum LAI spatially averaged over the irrigated tile fraction
  !
  !    4.2) IRRIG_TRIGGER: 1 - SUBROUTINE irrigrate_crop_calendar
  !          uses 26 crop calendars based on monthly crop growing areas of below crops.
  !          1       Wheat                      14      Oil palm                              
  !          2       Maize                      15      Rape seed / Canola    
  !          3       Rice                       16      Groundnuts / Peanuts  
  !          4       Barley                     17      Pulses                
  !          5       Rye                        18      Citrus                
  !          6       Millet                     19      Date palm             
  !          7       Sorghum                    20      Grapes / Vine         
  !          8       Soybeans                   21      Cotton                
  !          9       Sunflower                  22      Cocoa                 
  !          10      Potatoes                   23      Coffee                
  !          11      Cassava                    24      Others perennial      
  !          12      Sugar cane                 25      Fodder grasses        
  !          13      Sugar beet                 26      Others annual         
  !
  !     IRRIG_TRIGGER: 1 SPECIFIC INPUTS:
  !          DOFYR        : day of year
  !          IRRIGTYPE    : Preferred Irrig method (NTILES, 26) -
  !                         (0)CONCURRENT (default), (1)SPRINKLER ONLY (2)DRIP ONLY (3)FLOOD ONLY, and (-negative) AVOID this method 
  !          CROPIRRIGFRAC: Crop irrigated fraction (NTILES, 26) (per Section 2, fractions have been adjusted such that
  !                         CROPIRRIGFRAC is 1. on paddy tiles; the sum of available crop fractions is 1. on irrigated crop tiles;
  !                         and is zero on non-irrigated tiles.
  !          IRRIGPLANT   : DOY start planting (NTILES, 2, 26) - up to two seasons
  !          IRRIGHARVEST : DOY end harvesting (NTILES, 2, 26) - up to two seasons
  !          If IRRIGPLANT/IRRIGHARVEST = 998, the crop is not grown on that tile 
  !
  ! (5) MODEL UPDATES (OPTIONAL INTERNALS):
  !     SRATE, DRATE, and FRATE contain irrigation rates applied on individual fractions at any given time.
  !     The second dimensions of 2D arrays is for different crop fractions i.e. the second dimension is 2 for above
  !     IRRIG_TRIGGER: 0 to separately store irrigation rates in irrigated crop and paddy fractions.
  !     It would be 26 for IRRIG_TRIGGER: 1.
  !     Note also that runnning the irrigation model on subtiling mode (specific irrigated crop and paddy tiles with their own land prognostics)
  !     is preferred (Section 2) to running the irrigation model on fractions of typical land tiles. For the subtiling mode with own soil moisture
  !     prognostics, IRRIGFRAC, PADDYFRAC and CROPIRRIGFRAC fractions have been adjusted to represent irrigation type of the subtile in question. i.e.  
  !     IRRIFRAC has been set to 1. on irrigated-crop subtiles; PADDYFRAC and CROPIRRIGFRAC(N,3) are set to 1. on paddy subtiles; and  SUM of
  !     CROPIRRIGFRAC(N,:) excluding the 3rd element is set to 1. on irrigated crop tiles.
  !
  !     The crop calendar implemetation (IRRIG_TRIGGER: 1) computes SPRINKLERRATE, DRIPRATE, and FLOODRATE as weighted averages of irrigation rates from
  !     all active crops in SRATE, DRATE and FRATE arrays.
  
  PRIVATE
  
  INTEGER, PARAMETER, PUBLIC :: NUM_CROPS = 26, NUM_SEASONS = 2

  type, public :: irrig_params
     
     ! Below parameters can be set via RC file.

     REAL :: irrig_thres      =  0.5 ! threshold of tile fraction to turn the irrigation model on. 
     REAL :: lai_thres        =  0.6 ! threshold of LAI range to turn irrigation on
     REAL :: efcor            = 30.0 ! Efficiency Correction (% water loss: efcor = 0% denotes 100% efficient water use)
     REAL :: MIDS_LENGTH      =  0.6 ! Mid-season length as a fraction of crop growing season length (to be used with IRRIG_TRIGGER: 1)
     
     ! Sprinkler parameters
     ! --------------------
     REAL :: sprinkler_stime  =  6.0 ! sprinkler irrigatrion start time [hours]
     REAL :: sprinkler_dur    =  4.0 ! sprinkler irrigation duration [hours]
     REAL :: sprinkler_thres  =  0.7 ! soil moisture threshhold to trigger sprinkler irrigation
     
     ! Drip parameters 
     ! ---------------
     REAL :: drip_stime       =  8.0 ! drip irrigatrion start time [hours] 
     REAL :: drip_dur         =  8.0 ! drip irrigation duration [hours]
     
     ! Flood parameters
     ! ----------------
     REAL :: flood_stime      =  6.0 ! flood irrigatrion start time [hours]
     REAL :: flood_dur        =  1.0 ! flood irrigation duration [hours]
     REAL :: flood_thres      =  0.6 ! soil moisture threshhold to trigger flood irrigation

     
  end type irrig_params
  
  type, public, extends (irrig_params) :: irrigation_model
     
   contains
     
     ! public
     procedure, public :: init_model
     generic,   public :: run_model => irrigrate_lai_trigger, irrigrate_crop_calendar
     generic,   public :: update_irates => update_irates_lai, update_irates_ccalendar
     
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
    CALL ESMF_ConfigGetAttribute (SCF, label='SPRINKLER_STIME:', VALUE=IP%sprinkler_stime, DEFAULT=DP%sprinkler_stime, __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='SPRINKLER_DUR:'  , VALUE=IP%sprinkler_dur,   DEFAULT=DP%sprinkler_dur  , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='SPRINKLER_THRES:', VALUE=IP%sprinkler_thres, DEFAULT=DP%sprinkler_thres, __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='DRIP_STIME:'     , VALUE=IP%drip_stime,      DEFAULT=DP%drip_stime     , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='DRIP_DUR:'       , VALUE=IP%drip_dur,        DEFAULT=DP%drip_dur       , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='FLOOD_STIME:'    , VALUE=IP%flood_stime,     DEFAULT=DP%flood_stime    , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='FLOOD_DUR:'      , VALUE=IP%flood_dur,       DEFAULT=DP%flood_dur      , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='FLOOD_THRES:'    , VALUE=IP%flood_thres,     DEFAULT=DP%flood_thres    , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='IRR_EFCOR:'      , VALUE=IP%efcor,           DEFAULT=DP%efcor          , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='LAI_THRES:'      , VALUE=IP%lai_thres,       DEFAULT=DP%lai_thres      , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='MIDS_LENGTH:'    , VALUE=IP%MIDS_LENGTH,     DEFAULT=DP%lai_thres      , __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='IRRIG_THRES:'    , VALUE=IP%irrig_thres,     DEFAULT=DP%irrig_thres    , __RC__ )
    CALL ESMF_ConfigDestroy      (SCF, __RC__)

  END SUBROUTINE init_model

  ! ----------------------------------------------------------------------------

  SUBROUTINE irrigrate_lai_trigger (this,IRRIG_METHOD, local_hour,         &
            IRRIGFRAC, PADDYFRAC, SPRINKLERFR, DRIPFR, FLOODFR,            &           
            SMWP, SMSAT, SMREF, SMCNT, LAI, LAIMIN,LAIMAX, RZDEF,          &
            SPRINKLERRATE, DRIPRATE, FLOODRATE, SRATE, DRATE, FRATE)

    implicit none
    class (irrigation_model), intent(inout) :: this
    integer, intent (in)                    :: IRRIG_METHOD
    real, dimension (:), intent (in)        :: local_hour
    real, dimension (:), intent (in)        :: IRRIGFRAC, PADDYFRAC, SPRINKLERFR, &
         DRIPFR, FLOODFR, SMWP, SMSAT, SMREF, SMCNT, LAI, LAIMIN, LAIMAX, RZDEF
    real, dimension (:), intent (inout)     :: SPRINKLERRATE, DRIPRATE, FLOODRATE
    real, dimension (:,:),intent (inout)    :: SRATE, DRATE, FRATE
    INTEGER                                 :: NTILES, N, crop
    REAL                                    :: ma, H1, H2, HC, IT, ROOTFRAC, LAITHRES
    logical                                 :: season_end
    
    NTILES = SIZE (IRRIGFRAC)
    TILE_LOOP : DO N = 1, NTILES
       IF(LAIMAX (N) > LAIMIN (N)) THEN
          LAITHRES = LAIMIN (N) + this%lai_thres * (LAIMAX (N) - LAIMIN (N))          
          ROOTFRAC = MIN((LAI(N) - LAIMIN (N)) / (LAIMAX(N) - LAIMIN(N)) ,1.0)          
       ELSE
          ROOTFRAC = 0.
       ENDIF
       HC = local_hour(n)

       season_end        =  .true.
       
       CHECK_LAITHRES : IF (LAI(N) >= LAITHRES) THEN
          season_end = .false.
          CHECK_IRRIGFRACS: IF (IRRIGFRAC(N) > 0.) THEN

             !-----------------------------------------------------------------------------
             !     Get the root zone moisture availability to the plant
             !-----------------------------------------------------------------------------

             if(SMREF(N) > SMWP(N))then
                ma = (SMCNT(N) - SMWP(N)) /(SMREF(N) - SMWP(N))
             else
                ma = -1.
             endif
             
             if(ma >= 0) then
                                
                SELECT CASE (IRRIG_METHOD)                
                CASE (0)  ! CONCURRENTLY SPRINKER + FLOOD + DRIP on corresponding fractions

                   call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), &
                        SRATE = SRATE (N,1), &
                        DRATE = DRATE (N,1), &
                        FRATE = FRATE (N,1))

                   SRATE (N,1) =  SRATE (N,1)*SPRINKLERFR(N)
                   DRATE (N,1) =  DRATE (N,1)*DRIPFR (N)
                   FRATE (N,1) =  FRATE (N,1)*FLOODFR (N) 
                   
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
                   PRINT *, 'irrigrate_lai_trigger: IRRIG_METHOD can be 0,1,2, or3'
                   CALL EXIT(1)
                END SELECT
             endif
                          
          IF (PADDYFRAC (N) > 0.) THEN
             
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
             
          ELSE
             
             SRATE (N,:) = 0.
             DRATE (N,:) = 0.
             FRATE (N,:) = 0.
             
          ENDIF CHECK_IRRIGFRACS
       ENDIF CHECK_LAITHRES

       ! turn off irrigation if LAI is smaller than the LAI trigger marking end of the season
       if(season_end) then
          DO crop = 1, 2
             SRATE (N,crop) = 0.
             DRATE (N,crop) = 0.
             FRATE (N,crop) = 0.
          END DO
       endif
       
    END DO TILE_LOOP

    ! Update SPRINKLERRATE, DRIPRATE,  FLOODRATE EXPORTS to be sent to land models.
    ! FLOODRATE is weighted averaged over irrigated crops + paddy fractions.
        
    call this%update_irates (SPRINKLERRATE,DRIPRATE,FLOODRATE, &
         IRRIGFRAC,PADDYFRAC,SRATE,DRATE,FRATE)
    
  END SUBROUTINE irrigrate_lai_trigger

  ! ----------------------------------------------------------------------------

  SUBROUTINE irrigrate_crop_calendar(this,dofyr,local_hour, &
       SPRINKLERFR, DRIPFR, FLOODFR,                        &
       CROPIRRIGFRAC,IRRIGPLANT, IRRIGHARVEST, IRRIGTYPE ,  &
       SMWP,SMSAT,SMREF,SMCNT, RZDEF,                       &  
       SPRINKLERRATE, DRIPRATE, FLOODRATE, SRATE, DRATE, FRATE)

    implicit none
    class(irrigation_model),intent(inout):: this
    integer, intent (in)                 :: dofyr
    real, dimension (:),   intent (in)   :: local_hour, SPRINKLERFR, DRIPFR, FLOODFR
    real, dimension (:),   intent (in)   :: SMWP, SMSAT, SMREF, SMCNT, RZDEF
    real, dimension(:,:),  intent (in)   :: CROPIRRIGFRAC ! NUM_CROPS
    real, dimension(:,:),  intent (in)   :: IRRIGTYPE     ! NUM_CROPS
    real, dimension(:,:,:),intent (in)   :: IRRIGPLANT    ! NUM_SEASONS, NUM_CROPS
    real, dimension(:,:,:),intent (in)   :: IRRIGHARVEST  ! NUM_SEASONS, NUM_CROPS
    real, dimension (:),intent (inout)   :: SPRINKLERRATE, DRIPRATE, FLOODRATE
    real, dimension (:,:),intent (inout) :: SRATE, DRATE, FRATE
    INTEGER                              :: NTILES, N, crop, sea, ITYPE, I
    REAL                                 :: ma, H1, H2, HC, IT, ROOTFRAC, void_frac
    logical                              :: season_end (NUM_CROPS)
    NTILES = SIZE (local_hour)
         
    TILE_LOOP : DO N = 1, NTILES
       HC = local_hour(n)
       IF_IRR: if(SUM(CROPIRRIGFRAC(N,:)) > 0.) then
          ! the tile is irrigated crop or paddy
          season_end =  .true.
          CROP_LOOP: DO crop = 1, NUM_CROPS
             CROP_IN_TILE: if(CROPIRRIGFRAC(N,crop) > 0.) then
                ! crop is grown in this tile
                TWO_SEASONS: do sea = 1, NUM_SEASONS
                   IS_CROP: IF(IRRIGPLANT(N, sea, crop) /= 998) THEN
                      ! crop is grown in sea
                      IS_SEASON: IF(IS_WITHIN_SEASON(dofyr,NINT(IRRIGPLANT(N, sea, crop)),NINT(IRRIGHARVEST(N, sea, crop)))) THEN
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
                            
                            ROOTFRAC = CROP_SEASON_STAGE (this%MIDS_LENGTH, dofyr,NINT(IRRIGPLANT(N, sea, crop)),NINT(IRRIGHARVEST(N, sea, crop)))
                            if(SMREF(N) > SMWP(N))then
                               ma = (SMCNT(N) - SMWP(N)) /(SMREF(N) - SMWP(N))
                            else
                               ma = -1.
                            endif
                            
                            SOILM: if(ma >= 0) then
                               
                               ITYPE = NINT(IRRIGTYPE(N,crop))

                               CROP_IMETHOD: if (ITYPE == 0) then
                                  
                                  ! concurrently on sprinkler, drip and flood fractions
                                  call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), &
                                       SRATE = SRATE (N,crop), &
                                       DRATE = DRATE (N,crop), &
                                       FRATE = FRATE (N,crop))
                                  
                                  SRATE (N,crop) =  SRATE (N,crop)*SPRINKLERFR(N)
                                  DRATE (N,crop) =  DRATE (N,crop)*DRIPFR (N)
                                  FRATE (N,crop) =  FRATE (N,crop)*FLOODFR (N)

                               elseif (ITYPE > 0) then
                                  
                                  ! only this method
                                  if (ITYPE == 1) call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), SRATE = SRATE (N,crop))
                                  if (ITYPE == 2) call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), DRATE = DRATE (N,crop))
                                  if (ITYPE == 3) call this%irrig_by_method (HC, ma, ROOTFRAC, SMCNT(N), SMREF(N), FRATE = FRATE (N,crop))

                               elseif (ITYPE < 0) then
                                  
                                  ! crop does not use IRRIG_METHOD -(ITYPE)
                                  void_frac = 0.                                  
                                  DO I = 1,3
                                     if(I == ABS(ITYPE))then
                                        ! this itype isn't used by this crop other 2 fractions equally share this fraction
                                        if (I == 1) then
                                           void_frac = SPRINKLERFR(N)/2.
                                           SRATE(N,crop) = 0.
                                        endif
                                        if (I == 2) then
                                           void_frac = DRIPFR (N)/2.
                                           DRATE(N,crop) = 0.
                                        endif
                                        if (I == 3)then
                                           void_frac = FLOODFR (N)/2.
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
                                        if (I == 1) SRATE (N,crop) =  SRATE (N,crop)*(SPRINKLERFR(N) + void_frac)
                                        if (I == 2) DRATE (N,crop) =  DRATE (N,crop)*(DRIPFR (N)     + void_frac)
                                        if (I == 3) FRATE (N,crop) =  FRATE (N,crop)*(FLOODFR (N)    + void_frac)
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
          DO crop = 1, NUM_CROPS
             if(season_end(crop)) then
                SRATE (N,crop) = 0.
                DRATE (N,crop) = 0.
                FRATE (N,crop) = 0.
             endif
          END DO
                             
       endif IF_IRR
    END DO TILE_LOOP
    
    ! Update SPRINKLERRATE, DRIPRATE,  FLOODRATE EXPORTS to be sent to land models
    ! They are weighted averaged over 26 crop fractions.

    call this%update_irates (SPRINKLERRATE,DRIPRATE,FLOODRATE, &
       CROPIRRIGFRAC,SRATE,DRATE,FRATE)
    
  END SUBROUTINE irrigrate_crop_calendar

  ! ----------------------------------------------------------------------------

  SUBROUTINE update_irates_lai (this,SPRINKLERRATE,DRIPRATE,FLOODRATE, &
       IRRIGFRAC,PADDYFRAC,SRATE,DRATE,FRATE)
    
    implicit none

    class(irrigation_model),intent(inout):: this
    real, dimension (:), intent (in)     :: IRRIGFRAC, PADDYFRAC
    real, dimension (:,:), intent (in)   :: SRATE, DRATE, FRATE
    real, dimension (:),intent (inout)   :: SPRINKLERRATE, DRIPRATE, FLOODRATE
    integer                              :: N, NT

    ! INITIALIZE EXPORTS
    SPRINKLERRATE = 0.
    DRIPRATE      = 0.
    FLOODRATE     = 0.

    NT = size (IRRIGFRAC)

    !_ASSERT(size (SRATE,2)==NUM_CROPS,'Irrigation model LAI trigger irrig tile types mismatch')
    
    DO N = 1, NT
       IF ((IRRIGFRAC(N) + PADDYFRAC(N))  > 0.) THEN
          SPRINKLERRATE (N) = SRATE (N,1)
          DRIPRATE (N)      = DRATE (N,1)
          FLOODRATE (N)     = (IRRIGFRAC(N)* FRATE (N,1) + PADDYFRAC(N)*FRATE (N,2)) &
               /(IRRIGFRAC(N) + PADDYFRAC(N)) 
       ENDIF
    END DO
    
  END SUBROUTINE update_irates_lai

  !...............................................................................
  
  SUBROUTINE update_irates_ccalendar(this,SPRINKLERRATE,DRIPRATE,FLOODRATE, &
       CROPIRRIGFRAC,SRATE,DRATE,FRATE)

    implicit none
    class(irrigation_model),intent(inout):: this
    real, dimension(:,:),  intent (in)   :: CROPIRRIGFRAC ! NUM_CROPS
    real, dimension (:,:), intent (in)   :: SRATE, DRATE, FRATE
    real, dimension (:),intent (inout)   :: SPRINKLERRATE, DRIPRATE, FLOODRATE
    integer                              :: N, NT, crop

    ! INITIALIZE EXPORTS
    SPRINKLERRATE = 0.
    DRIPRATE      = 0.
    FLOODRATE     = 0.

    !_ASSERT(size (SRATE,2)==NUM_CROPS,'Irrigation model crop calandar trigger NUM_CROPS mismatch')

    NT =  size (SPRINKLERRATE)
    DO N = 1, NT
       if(SUM(CROPIRRIGFRAC(N,:)) > 0.) then
          DO crop = 1, NUM_CROPS
             SPRINKLERRATE(N) = SPRINKLERRATE(N) + SRATE (N,crop)*CROPIRRIGFRAC(N,crop)/SUM(CROPIRRIGFRAC(N,:))
             DRIPRATE(N)      = DRIPRATE(N)      + DRATE (N,crop)*CROPIRRIGFRAC(N,crop)/SUM(CROPIRRIGFRAC(N,:))
             FLOODRATE(N)     = FLOODRATE(N)     + FRATE (N,crop)*CROPIRRIGFRAC(N,crop)/SUM(CROPIRRIGFRAC(N,:))
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
       IT = this%sprinkler_thres 
       if ((HC >= H1).AND.(HC < H2)) then
          ! use SMCNT at H1 during H1 <= HC < H2 to compute irrigrate.
          ! Notice drip uses the same soil moisture threshold of sprinkler but with 0.% efficiency correction.
          if((ma <= IT).AND.(H1 == HC)) &
               DRATE = this%cwd(ROOTFRAC,SMCNT,SMREF,0.)/(H2 - H1)/3600.
       else
          DRATE = 0.
       endif
    endif

    if(present (FRATE)) then
       ! FLOOD IRRIGATION
       H1 = this%flood_stime
       H2 = this%flood_stime + this%flood_dur
       IT = this%flood_thres
       if ((HC >= H1).AND.(HC < H2)) then
          ! use SMCNT at H1 during H1 <= HC < H2 to compute irrigrate.
          if((ma <= IT).AND.(H1 == HC)) &
               FRATE = this%cwd (ROOTFRAC,SMCNT,SMREF,this%efcor)/(H2 - H1)/3600.
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
