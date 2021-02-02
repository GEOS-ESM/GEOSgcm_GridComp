#include "MAPL_Generic.h"

MODULE IRRIGATION_MODULE

  USE MAPL
  USE ESMF
  
  IMPLICIT NONE

  ! This module computes irrigation rates by 3 different methods: sprinkler, flood and drip.
  ! Computed irrigation rates return to the land model as rates that water is added to
  ! the hydrological cycle by irrigation. Thus, land models add irrigation feedback: 
  !      sprinkler irrigation rate to large scale precipitation;
  !      drip irrigation volume to rootzone excess; and
  !      flood irrigation volume to rootzone excess.
  ! Computed irrigation rates are treated as internals to ensure stop-start zerodiff of the modeling system.
  ! The model uses rootzone soil moisture state at the local start time of irrigation to compute 
  ! irrigation rates for the day and maintains the same rate through out the irrigation duration.
  ! 
  ! Sprinkler and Flood Irrigation methods were adapted from LIS CLSMF2.5 irrigatrion module:
  ! https://github.com/NASA-LIS/LISF/blob/master/lis/surfacemodels/land/clsm.f2.5/irrigation/clsmf25_getirrigationstates.F90 
  ! Drip irrigation method calculation is similar to that of sprinkler, albeit the drip irrigation method assumes a 0% water loss.
  !
  ! December 24, 2020 (Sarith) - First Version
  !
  ! (1) MODEL OUTPUTS:
  !    1) SPRINKLERRATE [kg m-2 s-1]
  !    2) DRIPRATE [kg m-2 s-1]
  !    3) FLOODRATE [kg m-2 s-1]    
  !
  ! (2) COMPUTATIONAL TILES:
  !    During land BCs generation, land tiles where irrigated crops or paddy is present were further subdivided to a
  !    mosaic of upto 3 fractions: i) non-irrigated land, ii) irrigated crop, or iii) paddy.
  !    The model treats each fraction as a separate computational tile and runs on each individual fraction with own parameters and prognostics.
  !    All fractions inherited model and soil parameters from the main land tile that they belong to. A special treatment of setting BF3 to a high
  !    value (25.) was applied to paddy tiles to account for the uniquely flat nature of paddy fields. Vegetation characteristics and vegetation dynamic
  !    parameters for irrigated crop and paddy tiles were taken from the nearest grass or crops land tile. 
  !    During tiling and BCs data preparation, computed fractional coverages for land tiles were also adjusted
  !    to reflect each computational tile under the land grid component represents entirely one of the 3 irrigation surface types: a non-irrigated land,
  !    OR a irrigated-crops OR a paddy tile.
  !
  ! (3) MODEL INPUTS:
  !     SMWP    : rootzone soil moisture content at wilting point [mm] 
  !     SMSAT   : rootzone soil moisture content at saturation [mm] 
  !     SMREF   : rootzone soil moisture at lower tercile [mm]
  !     SMCNT   : currrent root zone soil moisture content [mm] 
  !     RZDEF   : rootzone moisture deficit to reach complete soil saturation for paddy [mm]
  !     LOCAL_HOUR to set irrigation switch.
  !
  ! (3) SEASONAL CYLCE OF CROP WATER DEMAND:
  ! The module provides 2 options to determine the seasonal cycle of crop water demand:
  !    3.1) IRRIG_TRIGGER: 0 - SUBROUTINE irrigrate_lai_trigger
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
  !    3.2) IRRIG_TRIGGER: 1 - SUBROUTINE irrigrate_crop_calendar
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
  !          IRRIGTYPE    : Preferred Irrig method (NTILES, 26) - (1)SPRINKLER_(2)DRIP_(3)FLOOD
  !          CROPIRRIGFRAC: Crop irrigated fraction (NTILES, 26) (per Section 2, fractions have been adjusted such that
  !                         CROPIRRIGFRAC is 1. on paddy tiles; the sum of available crop fractions is 1. on irrigated crop tiles;
  !                         and is zero on non-irrigated tiles.
  !          IRRIGPLANT   : DOY start planting (NTILES, 2, 26) - up to two seasons
  !          IRRIGHARVEST : DOY end harvesting (NTILES, 2, 26) - up to two seasons
  !          If IRRIGPLANT/IRRIGHARVEST = 998, the crop is not grown on that tile 
  
  PRIVATE
  
  INTEGER, PARAMETER, PUBLIC :: NUM_CROPS = 26, NUM_SEASONS = 2

  type, public :: irrig_params
     
     ! Below parameters can be set via RC file.
     
     REAL :: lai_thres        =  0.6 ! threshold of LAI range to turn irrigation on
     REAL :: efcor            = 24.0 ! Efficiency Correction (% water loss: efcor = 0% denotes 100% efficient water use)
     
     ! Sprinkler parameters
     ! --------------------
     REAL :: sprinkler_stime  =  6.0 ! sprinkler irrigatrion start time [hours]
     REAL :: sprinkler_dur    =  4.0 ! sprinkler irrigation duration [hours]
     REAL :: sprinkler_thres  =  0.5 ! soil moisture threshhold to trigger sprinkler irrigation
     
     ! Drip parameters 
     ! ---------------
     REAL :: drip_stime       = 10.0 ! drip irrigatrion start time [hours] 
     REAL :: drip_dur         =  8.0 ! drip irrigation duration [hours]
     
     ! Flood parameters
     ! ----------------
     REAL :: flood_stime      =  6.0 ! flood irrigatrion start time [hours]
     REAL :: flood_dur        =  1.0 ! flood irrigation duration [hours]
     REAL :: flood_thres      = 0.25 ! soil moisture threshhold to trigger flood irrigation

     
  end type irrig_params
  
  type, public, extends (irrig_params) :: irrigation_model
     
   contains
     
     ! public
     procedure, public :: init_model
     generic,   public :: run_model => irrigrate_lai_trigger, irrigrate_crop_calendar

     ! private
     procedure, private :: irrigrate_lai_trigger
     procedure, private :: irrigrate_crop_calendar
     procedure, private :: cwd => crop_water_deficit
     
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
    CALL ESMF_ConfigDestroy      (SCF, __RC__)

  END SUBROUTINE init_model

  ! ----------------------------------------------------------------------------

  SUBROUTINE irrigrate_lai_trigger (this,IRRIG_METHOD, local_hour,         &
            IRRIGFRAC, PADDYFRAC, SPRINKLERFR, DRIPFR, FLOODFR,            &           
            SMWP, SMSAT, SMREF, SMCNT, LAI, LAIMIN,LAIMAX, RZDEF,          &
            SPRINKLERRATE, DRIPRATE, FLOODRATE)

    implicit none
    class (irrigation_model), intent(inout) :: this
    integer, intent (in)                    :: IRRIG_METHOD
    real, dimension (:), intent (in)        :: local_hour
    real, dimension (:), intent (in)        :: IRRIGFRAC, PADDYFRAC, SPRINKLERFR, &
         DRIPFR, FLOODFR, SMWP, SMSAT, SMREF, SMCNT, LAI, LAIMIN, LAIMAX, RZDEF
    real, dimension (:), intent (inout)     :: SPRINKLERRATE, DRIPRATE, FLOODRATE
    INTEGER                                 :: NTILES, N
    REAL                                    :: ma, H1, H2, HC, IT, LF, LAITHRES

    NTILES = SIZE (IRRIGFRAC)
    TILE_LOOP : DO N = 1, NTILES
       IF(LAIMAX (N) > LAIMIN (N)) THEN
          LAITHRES = LAIMIN (N) + this%lai_thres * (LAIMAX (N) - LAIMIN (N))          
          LF       = (LAI(N) - LAIMIN (N)) / (LAIMAX(N) - LAIMIN(N))
       ELSE
          LF = 0.
       ENDIF
       HC = local_hour(n)
       CHECK_LAITHRES : IF (LAI(N) >= LAITHRES) THEN          
          CHECK_IRRIGFRACS: IF (IRRIGFRAC(N) > 0.) THEN

             !-----------------------------------------------------------------------------
             !     Get the root zone moisture availability to the plant
             !-----------------------------------------------------------------------------

             if(SMREF(N) > SMWP(N))then
                ma = (LF*SMCNT(N) - SMWP(N)) /(SMREF(N) - SMWP(N))
             else
                ma = -1.
             endif
             
             if(ma >= 0) then
                                
                SELECT CASE (IRRIG_METHOD)                
                CASE (0)  ! CONCURRENTLY SPRINKER + FLOOD + DRIP on corresponding fractions
                   
                   H1 = this%sprinkler_stime
                   H2 = this%sprinkler_stime + this%sprinkler_dur
                   IT = this%sprinkler_thres                   
                   if ((HC >= H1).AND.(HC < H2)) then
                      if((ma <= IT).AND.(H1 == HC)) &
                           SPRINKLERRATE (N) = this%cwd (LF*SMCNT(N),SMREF(N),this%efcor) &
                                               * SPRINKLERFR(N)/(H2 - H1)/3600.
                   else
                      SPRINKLERRATE (N) = 0
                   endif
                   
                   H1 = this%flood_stime
                   H2 = this%flood_stime + this%flood_dur
                   IT = this%flood_thres
                   if ((HC >= H1).AND.(HC < H2)) then
                      if((ma <= IT).AND.(H1 == HC)) &
                           FLOODRATE (N) = this%cwd (LF*SMCNT(N),SMREF(N),this%efcor) &
                                           * FLOODFR (N)/(H2 - H1)/3600.
                   else
                      FLOODRATE (N) = 0.
                   endif

                   H1 = this%drip_stime
                   H2 = this%drip_stime + this%drip_dur
                   IT = this%sprinkler_thres 
                   if ((HC >= H1).AND.(HC < H2)) then
                      if((ma <= IT).AND.(H1 == HC)) DRIPRATE (N) = this%cwd(SMCNT(N),SMREF(N),0.) * DRIPFR (N) /(H2 - H1)/3600.
                   else
                      DRIPRATE (N) = 0.
                   endif
                   
                CASE (1)  ! SPRINKLER only
                   
                   H1 = this%sprinkler_stime
                   H2 = this%sprinkler_stime + this%sprinkler_dur
                   IT = this%sprinkler_thres    
                   if ((HC >= H1).AND.(HC < H2)) then
                      if((ma <= IT).AND.(H1 == HC)) &
                           SPRINKLERRATE (N) = this%cwd(LF*SMCNT(N),SMREF(N),this%efcor) &
                                               /(H2 - H1)/3600.
                   else
                      SPRINKLERRATE (N) = 0
                   endif
                   
                   DRIPRATE (N) = 0.
                   FLOODRATE (N)= 0.
                                      
                CASE (2)  ! DRIP only
                   
                   H1 = this%drip_stime
                   H2 = this%drip_stime + this%drip_dur
                   IT = this%sprinkler_thres 
                   if ((HC >= H1).AND.(HC < H2)) then
                      if((ma <= IT).AND.(H1 == HC)) DRIPRATE (N) = this%cwd(LF*SMCNT(N),SMREF(N),0.) &
                                               /(H2 - H1)/3600.
                   else
                      DRIPRATE (N) = 0.
                   endif
                   SPRINKLERRATE (N) = 0.
                   FLOODRATE (N)     = 0.

                CASE (3)  ! FLOOD only
                   
                   H1 = this%flood_stime
                   H2 = this%flood_stime + this%flood_dur
                   IT = this%flood_thres
                   if ((HC >= H1).AND.(HC < H2)) then
                      if((ma <= IT).AND.(H1 == HC)) &
                           FLOODRATE (N) = this%cwd(LF*SMCNT(N),SMREF(N),this%efcor) &
                                           /(H2 - H1)/3600.
                   else
                      FLOODRATE (N) = 0.
                   endif
                   SPRINKLERRATE (N) = 0.
                   DRIPRATE (N)      = 0.
                   
                CASE DEFAULT
                   PRINT *, 'irrigrate_lai_trigger: IRRIG_METHOD can be 0,1,2, or3'
                   CALL EXIT(1)
                END SELECT
             endif
             
          ELSEIF (PADDYFRAC (N) > 0.) THEN
             
             H1 = this%flood_stime
             H2 = this%flood_stime + this%flood_dur
             if ((HC >= H1).AND.(HC < H2)) then
                if(H1 == HC) FLOODRATE (N) = RZDEF(N) /(H2 - H1)/ 3600.
             else
                 FLOODRATE (N) = 0.
             endif
             SPRINKLERRATE (N) = 0.
             DRIPRATE (N)      = 0.
             
          ELSE
             
             SPRINKLERRATE (N) = 0.
             DRIPRATE (N)      = 0.
             FLOODRATE (N)     = 0.
             
          ENDIF CHECK_IRRIGFRACS
       ENDIF CHECK_LAITHRES
    END DO TILE_LOOP

  END SUBROUTINE irrigrate_lai_trigger

  ! ----------------------------------------------------------------------------

  SUBROUTINE irrigrate_crop_calendar(this,dofyr,local_hour, &
       CROPIRRIGFRAC,IRRIGPLANT, IRRIGHARVEST, IRRIGTYPE ,  &
       SMWP,SMSAT,SMREF,SMCNT, RZDEF,                       &  
       SPRINKLERRATE, DRIPRATE, FLOODRATE)

    implicit none
    class(irrigation_model),intent(inout):: this
    integer, intent (in)                 :: dofyr
    real, dimension (:),   intent (in)   :: local_hour
    real, dimension (:),   intent (in)   :: SMWP, SMSAT, SMREF, SMCNT, RZDEF
    real, dimension(:,:),  intent (in)   :: CROPIRRIGFRAC ! NUM_CROPS
    real, dimension(:,:),  intent (in)   :: IRRIGTYPE     ! NUM_CROPS
    real, dimension(:,:,:),intent (in)   :: IRRIGPLANT    ! NUM_SEASONS, NUM_CROPS
    real, dimension(:,:,:),intent (in)   :: IRRIGHARVEST  ! NUM_SEASONS, NUM_CROPS
    real, dimension (:),intent (inout)   :: SPRINKLERRATE, DRIPRATE, FLOODRATE
    INTEGER                              :: NTILES, N, crop, sea
    REAL                                 :: ma, H1, H2, HC, IT, ICFRAC
    REAL, ALLOCATABLE, DIMENSION (:)     :: FRATE, SRATE, DRATE

    NTILES = SIZE (local_hour)
    ALLOCATE (FRATE (NUM_CROPS))
    ALLOCATE (SRATE (NUM_CROPS))
    ALLOCATE (DRATE (NUM_CROPS))
    
    TILE_LOOP : DO N = 1, NTILES
       HC = local_hour(n)
       IRR_OR_NOT: if(SUM(CROPIRRIGFRAC(N,:)) > 0.) then
          
          if(SMREF(N) > SMWP(N))then
             ma = (SMCNT(N) - SMWP(N)) /(SMREF(N) - SMWP(N))
          else
             ma = -1.
          endif
          ICFRAC = SUM (CROPIRRIGFRAC(N,:)) - CROPIRRIGFRAC(N,3)
          IF (ICFRAC > 0.) THEN
             srate = 0.
             drate = 0.
             frate = 0.
             DO crop = 1, NUM_CROPS
                if(crop /= 3) then
                   if ((IRRIGTYPE (N,crop) == 1.).AND.(SUM (CROPIRRIGFRAC(N,:),mask=IRRIGTYPE (N,:) == 1.) > 0)) &
                        SRATE (crop) = SPRINKLERRATE(N)*CROPIRRIGFRAC(N,crop)/SUM (CROPIRRIGFRAC(N,:),mask=IRRIGTYPE (N,:) == 1.)
                   if ((IRRIGTYPE (N,crop) == 2.).AND.(SUM (CROPIRRIGFRAC(N,:),mask=IRRIGTYPE (N,:) == 2.) > 0)) &
                        DRATE (crop) = DRIPRATE(N)     *CROPIRRIGFRAC(N,crop)/SUM (CROPIRRIGFRAC(N,:),mask=IRRIGTYPE (N,:) == 2.)
                   if ((IRRIGTYPE (N,crop) == 3.).AND.(SUM (CROPIRRIGFRAC(N,:),mask=IRRIGTYPE (N,:) == 3.) > 0)) &
                        FRATE (crop) = FLOODRATE(N)    *CROPIRRIGFRAC(N,crop)/SUM (CROPIRRIGFRAC(N,:),mask=IRRIGTYPE (N,:) == 3.)
                endif
             END DO
           END IF
          
          CROP_LOOP: DO crop = 1, NUM_CROPS
             CROP_IN_TILE: if(CROPIRRIGFRAC(N,crop) > 0.) then
                do sea = 1, NUM_SEASONS
                      IS_CROP: IF(IRRIGPLANT(N, sea, crop) /= 998) THEN
                         IS_SEASON: IF(IS_WITHIN_SEASON(dofyr,NINT(IRRIGPLANT(N, sea, crop)),NINT(IRRIGHARVEST(N, sea, crop)))) THEN
                            PADDY_OR_CROP: if (crop == 3) then
                               
                               ! a paddy tile
                               H1 = this%flood_stime
                               H2 = this%flood_stime + this%flood_dur
                               if ((HC >= H1).AND.(HC < H2)) then
                                  if(H1 == HC) FLOODRATE (N) = RZDEF(N) /(H2 - H1)/ 3600.
                               else
                                  FLOODRATE (N) = 0.
                               endif
                               SPRINKLERRATE (N) = 0.
                               DRIPRATE (N)      = 0.

                            else
                               
                               ! irrigated crop - compute sum of irrigrates from 25 crops   
                               SELECT CASE (NINT(IRRIGTYPE(N,crop)))
                               CASE (1)
                                  ! (1)SPRINKLER
                                  H1 = this%sprinkler_stime
                                  H2 = this%sprinkler_stime + this%sprinkler_dur
                                  IT = this%sprinkler_thres
                                  if ((HC >= H1).AND.(HC < H2)) then
                                     if((ma >= 0.).AND.(ma < IT).AND.(H1 == HC)) &
                                          SRATE (crop) = this%cwd(SMCNT(N),SMREF(N),this%efcor) &
                                          /(H2 - H1)/3600.
                                  else
                                     SRATE (crop) = 0
                                  endif
                                   
                               CASE (2)
                                  ! (2)DRIP
                                  H1 = this%drip_stime
                                  H2 = this%drip_stime + this%drip_dur
                                  IT = this%sprinkler_thres 
                                  if ((HC >= H1).AND.(HC < H2)) then
                                     if((ma >= 0.).AND.(ma < IT).AND.(H1 == HC)) DRATE (crop) = this%cwd(SMCNT(N),SMREF(N),0.) &
                                          /(H2 - H1)/3600.
                                  else
                                     DRATE (crop) = 0.
                                  endif
                                  
                               CASE (3)
                                  ! (3)FLOOD
                                  H1 = this%flood_stime
                                  H2 = this%flood_stime + this%flood_dur
                                  IT = this%flood_thres
                                  if ((HC >= H1).AND.(HC < H2)) then
                                     if((ma >= 0.).AND.(ma < IT).AND.(H1 == HC)) &
                                          FRATE (crop) = this%cwd(SMCNT(N),SMREF(N),this%efcor) &
                                           /(H2 - H1)/3600.
                                  else
                                     FRATE (crop) = 0.
                                  endif
                                  
                               CASE DEFAULT
                                  PRINT *, 'irrigrate_crop_calendar: IRRIG_METHOD can be 1,2, or3'
                                  CALL EXIT(1)
                               END SELECT
                            ENDIF PADDY_OR_CROP
                         ENDIF IS_SEASON
                      end IF IS_CROP
                   end do
             endif CROP_IN_TILE
          END DO CROP_LOOP

          IF (ICFRAC > 0.) THEN             
             SPRINKLERRATE(N) = 0.
             DRIPRATE(N)      = 0.
             FLOODRATE(N)     = 0.
             DO crop = 1, NUM_CROPS
                if(crop /= 3) then
                   if (IRRIGTYPE (N,crop) == 1.) SPRINKLERRATE(N) = SPRINKLERRATE(N) + SRATE (crop)*CROPIRRIGFRAC(N,crop)!/SUM (CROPIRRIGFRAC(N,:),mask=IRRIGTYPE (N,:) == 1.)
                   if (IRRIGTYPE (N,crop) == 2.) DRIPRATE(N)      = DRIPRATE(N)      + DRATE (crop)*CROPIRRIGFRAC(N,crop)!/SUM (CROPIRRIGFRAC(N,:),mask=IRRIGTYPE (N,:) == 2.)
                   if (IRRIGTYPE (N,crop) == 3.) FLOODRATE(N)     = FLOODRATE(N)     + FRATE (crop)*CROPIRRIGFRAC(N,crop)!/SUM (CROPIRRIGFRAC(N,:),mask=IRRIGTYPE (N,:) == 3.)
                endif
             END DO
          ENDIF
          
       else
          
          SPRINKLERRATE (N) = 0.
          DRIPRATE (N)      = 0.
          FLOODRATE (N)     = 0.
          
       endif IRR_OR_NOT
    END DO TILE_LOOP

    deallocate (FRATE, SRATE, DRATE)
    
  END SUBROUTINE irrigrate_crop_calendar

  ! ----------------------------------------------------------------------------

  logical FUNCTION IS_WITHIN_SEASON (DOY,DP, DH)

    implicit none
    integer, intent (in) :: DOY,DP, DH

    IS_WITHIN_SEASON = .false.
    if(DH > DP) then
       if((DOY >= DP).AND.(DOY <=  DH)) IS_WITHIN_SEASON = .true.
    else
       if((DOY >= DP).AND.(DOY <= 366)) IS_WITHIN_SEASON = .true.
       if((DOY >=  1).AND.(DOY <=  DH)) IS_WITHIN_SEASON = .true.
    endif
      
  end FUNCTION IS_WITHIN_SEASON
  
  ! ----------------------------------------------------------------------------

  REAL FUNCTION crop_water_deficit (this, asmc, smcref, efcor)

    implicit none
    class(irrigation_model),intent(inout):: this
    real, intent (in)                    :: asmc, smcref, efcor

    crop_water_deficit = (smcref - asmc)*100.0/(100.0-efcor)   

  END FUNCTION crop_water_deficit
  
  ! ----------------------------------------------------------------------------

END MODULE IRRIGATION_MODULE
