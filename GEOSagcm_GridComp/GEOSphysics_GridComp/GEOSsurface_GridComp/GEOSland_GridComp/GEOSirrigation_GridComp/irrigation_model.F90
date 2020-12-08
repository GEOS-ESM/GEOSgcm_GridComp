MODULE IRRIGATION_MODULE

  USE ESMF
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC IRRIGATION, NUM_CROPS, NUM_SEASONS

  INTEGER, PARAMETER :: NUM_CROPS = 26, NUM_SEASONS = 2

CONTAINS
  
  SUBROUTINE IRRIGATION (SURFRC)

    IMPLICIT NONE

    CHARACTER(:), INTENT(IN) :: SURFRC
    
    ! Sprinkler parameters
    ! --------------------
    REAL :: sprinkler_stime  = 6.0 ! sprinkler irrigatrion start time [hours]
    REAL :: sprinkler_dur    = 4.0 ! sprinkler irrigation duration [hours]
    REAL :: sprinkler_thres  = 0.5 ! soil moisture threshhold to trigger sprinkler irrigation
      
    ! Drip parameters 
    ! ---------------
    REAL :: drip_stime       =  6.0 ! drip irrigatrion start time [hours]
    REAL :: drip_dur         = 12.0 ! drip irrigation duration [hours]
      
    ! Flood parameters
    ! ----------------
    REAL :: flood_stime      =  6.0  ! flood irrigatrion start time [hours]
    REAL :: flood_dur        =  1.0  ! flood irrigation duration [hours]
    REAL :: flood_thres      =  0.25 ! soil moisture threshhold to trigger flood irrigation
    REAL :: efcor            = 76.0  ! Efficiency Correction (%)    
    TYPE(ESMF_Config)       :: SCF 

    SCF = ESMF_ConfigCreate(__RC__) 
    CALL ESMF_ConfigLoadFile     (SCF,SURFRC,rc=status) ; VERIFY_(STATUS)
    CALL ESMF_ConfigGetAttribute (SCF, label='SPRINKLER_STIME:', VALUE=sprinkler_stime, DEFAULT= 6.00, __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='SPRINKLER_DUR:'  , VALUE=sprinkler_dur,   DEFAULT= 4.00, __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='SPRINKLER_THRES:', VALUE=sprinkler_thres, DEFAULT= 0.50, __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='DRIP_STIME:'     , VALUE=drip_stime,      DEFAULT= 6.00, __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='DRIP_DUR:'       , VALUE=drip_dur,        DEFAULT=12.00, __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='FLOOD_STIME:'    , VALUE=flood_stime,     DEFAULT= 6.00, __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='FLOOD_DUR:'      , VALUE=flood_dur,       DEFAULT= 1.00, __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='FLOOD_THRES:'    , VALUE=flood_thres,     DEFAULT= 0.25, __RC__ )
    CALL ESMF_ConfigGetAttribute (SCF, label='IRR_EFCOR:'      , VALUE=efcor,           DEFAULT=76.00, __RC__ )
    CALL ESMF_ConfigDestroy      (SCF, __RC__)
    
  CONTAINS
    
    ! ********************************************************************
  
    SUBROUTINE irrigation_rate (IRRIG_METHOD,                                    &
         NTILES, AGCM_HH, AGCM_MI, AGCM_S, lons, IRRIGFRAC, PADDYFRAC,  &
         CLMPT,CLMST, CLMPF, CLMSF, LAIMAX, LAIMIN, LAI,                &
         POROS, WPWET, VGWMAX, RZMC, IRRIGRATE)
      
      ! !DESCRIPTION:
      !
      ! NOTE: This is an experimental feature under development.
      !
      ! Calculate water requirement and apply the amount to precipitation.
      !
      ! Irrigate when available root zone soil moisture falls below tunable
      !        irrigation threshold parameter.
      ! Below GRIPC irrigated data provide fractions of croplands and paddy croplands.
      ! The irrigation model is applied on a tile if:
      ! (1) the irrigated fraction of the tile is greater than 0. AND
      ! (2) primary or secondary type in the tile is CLM4 type 16 (cropland) AND
      ! (3) LAI exceeds the LAI theshhold  (60% of LAI range)
      !
      ! GRIPC croplands and paddy croplands fractions determine whether to apply
      !    either sprinkler or flood OR both irrigation methods. Each method has
      !    its own local start times, durations and irrigation threshold parameters.
      !
      ! We assume plants need available soil moisture stay above 1/3 of soil moisture range
      !    [ wilting - saturation]
      ! Irrigation amount is scaled to grid total crop fraction when intensity
      ! is less than the fraction.  Irrigation is expanded to non-crop, non-forest,
      ! non-baresoil/urban tiles if intensity exceeds grid total crop fraction.
      ! In latter case, scaled irrigation is applied to grassland first,
      ! then further applied over the rest of tiles equally if the intensity
      ! exceeds grassland fraction as well.
      !
      ! Optionally efficiency correction is applied to account for field loss.
      !
      ! REVISION HISTORY:
      !
      ! Aug 2018: Sarith Mahanama ; Version 1 adapted from LIS subroutine clsmf25_getirrigationstates.F90
      
      IMPLICIT NONE
      
      ! INPUTS
      ! ------
      INTEGER, INTENT (in)                       :: IRRIG_METHOD, NTILES, AGCM_HH, AGCM_MI, AGCM_S
      REAL   , INTENT (in), DIMENSION (ntiles)   :: lons, IRRIGFRAC, PADDYFRAC, LAIMAX,  &
           LAIMIN, LAI, CLMPT,CLMST, CLMPF, CLMSF, POROS, WPWET, VGWMAX, RZMC
      ! IRRIG_METHOD : 0 sprinkler and flood combined; 1 sprinkler irrigation ; 2 flood irrigation
      ! AGCM_HH / AGCM_MI / AGCM_S/ lons : Current hour, minute, second (UTC) and longitude
      
      ! Irrigation hotspots : Using the Global Rain-Fed, Irrigated, and Paddy Croplands (GRIPC) Dataset  (Salmon et al., 2015)
      !       Salmon JM, Friedl MA, Frolking S, Wisser D and Douglas EM: Global rain-fed, irrigated,
      !             and paddy croplands: A new high resolution map derived from remote sensing, crop
      !             inventories and climate data, Int. J. Appl. Earth Obs. Geoinf, 38, 321–334,
      !             doi:10.1016/j.jag.2015.01.014, 2015.
      
      ! IRRIGFRAC                        : Fraction of irrigated croplands [-] = total number of 500m irrigated croplands pixels in the tile /
      !                                                                          total number of 500m pixels in the tile
      ! PADDYFRAC                        : Fraction of paddy croplands [-]     = total number of 500m paddy croplands pixels in the tile /
      !                                                                          total number of 500m pixels in the tile
      ! LAIMAX / LAIMIN / LAI            : Maximum, minimum and current Leaf Area Index
      ! CLMPT / CLMST                    : CLM4 primary and secondary types (Note type 16 is cropland)
      ! CLMPF / CLMSF                    : CLM4 fractions of primary and secondary types
      ! POROS / WPWET / VGWMAX / RZMC    : porosity [m3/m3], wilting point wetness [-], maximum and current root zone soil moisture content [m3/m3]
      
      ! ONLY output
      ! -----------
      REAL   , INTENT (out), DIMENSION (ntiles)  :: IRRIGRATE
      

      
      ! local vars
      ! ----------
      REAL    :: smcwlt, smcref, smcmax, asmc, laithresh, laifac, RZDEP, vfrac, ma,  &
           otimee, irrig_thresh, IrrigScale, s_irate, f_irate, local_long, local_hour
      INTEGER :: n, t, vtyp
      
      IRRIGRATE (:) = 0.
      
      TILE_LOOP : DO N = 1, NTILES
         
         local_long = 180. * lons(n) / PIE                                             ! local logitude [degrees]
         local_hour = AGCM_HH + AGCM_MI / 60. + AGCM_S / 3600. + 12.* local_long /180. ! local time [hours]
         IF (local_hour >= 24.) local_hour =  local_hour - 24.
         IF (local_hour <   0.) local_hour =  local_hour + 24.
         
         laithresh  = laimin (n) + 0.60 * (laimax (n) - laimin (n))
         IF(laimax (n) /= laimin (n)) THEN
            laifac = (lai(n) - laimin (n)) / (laimax(n) - laimin(n))
         ELSE
            laifac = 0.
         ENDIF
         
         RZDEP      = laifac * VGWMAX (n) / poros (n)                                  ! root zone depth [mm]
         smcwlt     = RZDEP * wpwet (n) * poros (n)                                    ! RZ soil moisture content at wilting point [mm]
         smcref     = RZDEP * (wpwet (n) + 0.333 * (1. - wpwet (n))) * poros(n)        ! RZ reference soil moisture content [mm]
         smcmax     = RZDEP * poros (n)                                                ! RZ soil moisture at saturatopm [mm]
         asmc       = RZDEP * rzmc (n)                                                 ! actual RZ soil moisture content [mm]
         
         CHECK_IRRIG_INTENSITY : IF ((IRRIGFRAC(N) + PADDYFRAC(N)) > 0.) THEN
            
            s_irate = 0.
            f_irate = 0.
            
            TWO_CLMTYPS : DO t = 1, 2
               
               IF (t == 1) THEN
                  ! Primary CLM fraction
                  vtyp  = NINT (CLMPT (n))
                  vfrac = CLMPF (n)
               ENDIF
               
               IF (t == 2) THEN
                  ! Secondary CLM fraction
                  vtyp  = NINT (CLMST (n))
                  vfrac = CLMSF (n)
               ENDIF
               
               CHECK_CROP_LAITHRESH : IF ((vtyp == 16).AND.(vfrac > 0.).AND.(lai(n) >= laithresh).AND.(laifac > 0.)) THEN
                  
                  !-----------------------------------------------------------------------------
                  !     Compute irrigation scale parameter :
                  !     Scale the irrigation intensity to the crop % when intensity < crop%.
                  !     Expand irrigation for non-crop, non-forest when intensity > crop %
                  !     in preference order of grassland first then rest.
                  !-----------------------------------------------------------------------------
                  
                  IF ((IRRIGFRAC(N) + PADDYFRAC(N)) <  vfrac) THEN
                     IrrigScale = vfrac / (IRRIGFRAC(N) + PADDYFRAC(N))
                  ELSE
                     IrrigScale = 1.
                  ENDIF
                  
                  !-----------------------------------------------------------------------------
                  !     Get the root zone moisture availability to the plant
                  !-----------------------------------------------------------------------------
                  
                  IF(smcref.GE.smcwlt) THEN
                     ma = (asmc - smcwlt) /(smcref - smcwlt)
                  ELSE
                     ma = -1
                  ENDIF
                  
                  SELECT CASE (IRRIG_METHOD)
                     
                     !--------------------------------------------------------------------------------------------------------------------------
                     !     IRRIGRATE : irrigation rate required to fill up water deficit before END OF IRRIGATION PERIOD  (otimee - local_hour)
                     !--------------------------------------------------------------------------------------------------------------------------
                     
                  CASE (0)
                     ! SPRINKLER AND FLOOD IRRIGATION COMBINED
                     ! ---------------------------------------
                     C_SPRINKLER : IF((IRRIGFRAC (N) > 0.).AND.(local_hour >= sprinkler_stime).AND. (local_hour < sprinkler_stime + sprinkler_dur)) THEN
                        otimee = sprinkler_stime + sprinkler_dur ;  irrig_thresh = sprinkler_thres
                        IF ((ma  <= irrig_thresh).AND.(ma.GE.0)) THEN
                           s_irate = crop_water_deficit (IRRIGFRAC (N) * irrigScale, asmc, smcref, efcor) / (otimee - local_hour) /3600.0
                        ENDIF
                     ENDIF C_SPRINKLER
                     
                     C_FLOOD : IF((PADDYFRAC (N) > 0.).AND.(local_hour >= flood_stime).AND. (local_hour <= flood_stime + flood_dur)) THEN
                        otimee = flood_stime + flood_dur ; irrig_thresh = flood_thres
                        IF ((ma  <= irrig_thresh).AND.(ma.GE.0)) THEN
                           f_irate = crop_water_deficit (PADDYFRAC (N) * irrigScale, asmc, smcref, efcor) / (otimee - local_hour) /3600.0
                        ENDIF
                     ENDIF C_FLOOD
                     
                     IRRIGRATE (N) = (s_irate * IRRIGFRAC (N) + f_irate * PADDYFRAC (N)) / (IRRIGFRAC(N) + PADDYFRAC(N)) ! weighted averaged sprinkler + flood
                     
                  CASE (1)
                     ! SPRINKLER IRRIGATION ONLY
                     ! -------------------------
                     SPRINKLER : IF(((IRRIGFRAC (N) + PADDYFRAC (N)) > 0.).AND.(local_hour >= sprinkler_stime).AND. (local_hour < sprinkler_stime + sprinkler_dur)) THEN
                        otimee = sprinkler_stime + sprinkler_dur ;  irrig_thresh = sprinkler_thres
                        IF ((ma  <= irrig_thresh).AND.(ma.GE.0)) THEN
                           IRRIGRATE (N) = crop_water_deficit ((IRRIGFRAC (N) + PADDYFRAC (N)) * irrigScale, asmc, smcref, efcor) / &
                                (otimee - local_hour) /3600.0
                        ENDIF
                     ENDIF SPRINKLER
                     
                  CASE (2)
                     ! FLOOD IRRIGATION ONLY
                     ! ---------------------
                     FLOOD : IF(((IRRIGFRAC (N) + PADDYFRAC (N)) > 0.).AND.(local_hour >= flood_stime).AND. (local_hour <= flood_stime + flood_dur)) THEN
                        otimee = flood_stime + flood_dur ; irrig_thresh = flood_thres
                        IF ((ma  <= irrig_thresh).AND.(ma.GE.0)) THEN
                           IRRIGRATE (N) = crop_water_deficit ((IRRIGFRAC (N) +PADDYFRAC (N)) * irrigScale, asmc, smcref, efcor) / &
                                (otimee - local_hour) /3600.0
                        ENDIF
                     ENDIF FLOOD
                     
                  CASE DEFAULT
                     PRINT *, 'IN IRRIGATION_RATE : IRRIGATION_METHOD can  be 0, 1, or 2'
                     CALL EXIT(1)
                  END SELECT
               END IF CHECK_CROP_LAITHRESH
            END DO TWO_CLMTYPS
         END IF CHECK_IRRIG_INTENSITY
      END DO TILE_LOOP
      
    END SUBROUTINE irrigation_rate
    
    ! ********************************************************************
    
    REAL FUNCTION sprinkler_irrig_water (lroot,  smc, rdpth, smcref, irrigScale)
      
      IMPLICIT NONE
      
      INTEGER, INTENT (in)                 :: lroot
      REAL, DIMENSION (lroot), INTENT (in) :: smc, rdpth
      REAL, INTENT (in)                    :: smcref, irrigScale
      REAL                                 :: water
      INTEGER                              :: k
      
      water = 0.
      DO k=1,lroot
         water = water + (smcref- smc(k))*rdpth(k)*1000.0
      ENDDO
      
      !-----------------------------------------------------------------------------
      !     Scale the irrigation intensity to the crop % when intensity < crop%.
      !     Expand irrigation for non-crop, non-forest when intensity > crop %
      !     in preference order of grassland first then rest.
      !     *scale is pre-computed for each tile in getirrpmapetc module in a way
      !     that is transparent to every tile irrigated or non-irrigated
      !-----------------------------------------------------------------------------
      
      water = water * irrigScale
      
      !-----------------------------------------------------------------------------
      !     Apply efficiency correction
      !-----------------------------------------------------------------------------
      
      sprinkler_irrig_water = water * (100.0/(100.0-efcor))
      
    END FUNCTION sprinkler_irrig_water
    
    ! ********************************************************************
    
    REAL FUNCTION drip_irrig_water
      
      IMPLICIT NONE
      
      REAL :: twater, twater2
      
      ! Need to get crop coefficient so that we can caculate unstressed Transp
      !       RC=RSMIN/(XLAI*RCS*RCT*RCQ)
      !       PCIRR=(RR+DELTA)/(RR*(1.+RC*CH)+DELTA)
      ! CALL TRANSP (with PCIRR)
      
      ! Then add enough water to get from actual Transp to unstressed Transp
      
      twater = 0.0
      !-----------------------------------------------------------------------------
      !     Apply efficiency correction
      !-----------------------------------------------------------------------------
      twater2 = twater
      drip_irrig_water = twater*(100.0/(100.0-efcor))
      
    END FUNCTION drip_irrig_water
    
    ! ********************************************************************
    
    REAL FUNCTION flood_irrig_water (mxsoildpth, smc, soildep, smcmax, irrigScale)
      
      IMPLICIT NONE
      
      INTEGER, INTENT (in)                      :: mxsoildpth
      REAL, DIMENSION (mxsoildpth), INTENT (in) :: smc, soildep
      REAL, INTENT (in)                         :: smcmax, irrigScale
      REAL                                      :: water
      INTEGER                                   :: l
      
      water = 0.
      
      DO l = 1, mxsoildpth
         IF( l == 1 ) THEN
            water = (smcmax - smc(l))*soildep(l)*1000.0
         ELSE
            ! BZ modification 4/2/2015 to saturate entire column and apply ippix 
            water = water + (smcmax - smc(l))*soildep(l)*1000.0
            !   twater = twater + (smcmax - noah33_struc(n)%noah(t)%smc(2))*sldpth(2)*1000.0
            !   twater = twater + (smcmax - noah33_struc(n)%noah(t)%smc(3))*sldpth(3)*1000.0
            !   twater = twater + (smcmax - noah33_struc(n)%noah(t)%smc(4))*sldpth(4)*1000.0
         ENDIF
      END DO
      
      !-----------------------------------------------------------------------------
      !     Scale the irrigation intensity to the crop % when intensity < crop%.
      !     Expand irrigation for non-crop, non-forest when intensity > crop %
      !     in preference order of grassland first then rest.
      !     *scale is pre-computed for each tile in getirrpmapetc module in a way
      !     that is transparent to every tile irrigated or non-irrigated
      !-----------------------------------------------------------------------------
      
      water = water * irrigScale
      
      !-----------------------------------------------------------------------------
      !     Apply efficiency correction
      !-----------------------------------------------------------------------------
      
      flood_irrig_water = water * (100.0/(100.0-efcor))
      
    END FUNCTION flood_irrig_water
  END SUBROUTINE IRRIGATION
END MODULE IRRIGATION_MODULE
