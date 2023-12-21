
module StieglitzSnow
  
  ! This is a merge of Teppei's snow tracer code with the Heracles (H52) version
  
  ! reichle, 12 Aug 2014 - moved "*_calc_asnow()" to here from catchment.F90
  !                      - renamed "get_tf_nd()" to "StieglitzSnow_calc_tpsnow()"
  !                      - removed "relayer()" (obsolete)
  !                      - renamed "N_sm" ("soil moisture") to "N_zones" (could be LandIce...)
  !                      - cleaned up MAPL constants
  !                      - additional minor clean-up
  ! reichle, 19 Nov 2015 - changed WEMIN back to 13 kg/m2 
  !                          (snow cover fraction better agrees w/ MODIS over land)
  ! Sarith, 4/21/2016    - WEMIN made public, sibalb uses it
  ! Justin, 4/16/2018    - moved WEMIN, AICEV, AICEN to SurfParams,
  !  	    		   removed LAND_UPD ifdef's
  
  USE MAPL_ConstantsMod, ONLY:  &
       PIE  => MAPL_PI,         &  ! -          
       ALHE => MAPL_ALHL,       &  ! J/kg  @15C 
       ALHM => MAPL_ALHF,       &  ! J/kg       
       TF   => MAPL_TICE,       &  ! K          
       RHOW => MAPL_RHOWTR         ! kg/m^3     
  
  USE MAPL_BaseMod,      ONLY: MAPL_LANDICE
 
  USE SurfParams,        ONLY: WEMIN, AICEV, AICEN
  
  implicit none

  public :: StieglitzSnow_snowrt           ! used by LandIce, Catchment[CN]
  public :: StieglitzSnow_trid             ! used by LandIce
  public :: StieglitzSnow_snow_albedo      ! used by LandIce, Catchment[CN], LDAS
  public :: StieglitzSnow_calc_asnow       ! used by          Catchment[CN], LDAS
  public :: StieglitzSnow_calc_tpsnow      ! used by          Catchment[CN], LDAS
  public :: StieglitzSnow_echo_constants   ! used by                         LDAS
  public :: StieglitzSnow_relayer          ! used by                         LDAS, land-atm DAS

  public :: StieglitzSnow_RHOMA            ! used by                         LDAS, land-atm DAS
  public :: StieglitzSnow_MINSWE           ! used by LandIce
  public :: StieglitzSnow_CPW              ! used by LandIce

  interface StieglitzSnow_calc_asnow

     module procedure StieglitzSnow_calc_asnow_1
     module procedure StieglitzSnow_calc_asnow_2
     module procedure StieglitzSnow_calc_asnow_3
     
  end interface StieglitzSnow_calc_asnow

  interface StieglitzSnow_calc_tpsnow
     
     module procedure StieglitzSnow_calc_tpsnow_scalar       ! replicates original get_tf0d()
     module procedure StieglitzSnow_calc_tpsnow_vector       ! replicates original get_tf_nd()

  end interface StieglitzSnow_calc_tpsnow

  ! constants specific to StieglitzSnow
  
  real,          parameter :: StieglitzSnow_RHOMA  = 500.    ! kg/m^3  maximum snow density
  real,          parameter :: StieglitzSnow_MINSWE = 0.013   ! kg/m^2  min SWE to avoid immediate melt
  real,          parameter :: StieglitzSnow_CPW    = 2065.22 ! J/kg/K  specific heat of ice at 0 deg C (??) [=MAPL_CAPICE??]

  real, private, parameter :: SNWALB_VISMIN        = 0.5    
  real, private, parameter :: SNWALB_NIRMIN        = 0.3

  !================================ Added by Teppei Yasunari ==================================
  !--------------------------------------------------------------------------------------------
  ! Teppei J. Yasunari produced this file, 23 May 2014
  ! Teppei, 23 May 2014 - Moved the constants from StieglitzSnow.F90 to here and revised comments.
  ! Teppei, 19 August 2014 - no longer considered goswim_constants.F90 
  !                          and put all the GOSWIM-related constants here
  ! Teppei, 27 August 2014 - if condition was revised for ALB_WITH_IMPURITY
  !                        - ABVIS and ABNIR are not needed in the subroutine because spcified below
  !
  !--------------------------------------------------------------------------------------------
  ! snow albedo related constants
  !--------------------------------------------------------------------------------------------
  !   See details in: 
  !       Yasunari, T. J., K.-M. Lau, S. P. P. Mahanama, P. R. Colarco, A. M. da Silva, 
  !       T. Aoki, K. Aoki, N. Murao, S. Yamagata, and Y. Kodama, 2014: 
  !       GOddard SnoW Impurity Module (GOSWIM) for the NASA GEOS-5 Earth System Model: 
  !       Preliminary comparisons with observations in Sapporo, Japan. SOLA, 10, 50-56, 
  !       doi:10.2151/sola.2014-011.
  !   The URL of the paper: https://www.jstage.jst.go.jp/article/sola/10/0/10_2014-011/_article
  !--------------------------------------------------------------------------------------------
  
  !  **************** IMPORTANT  IMPORTANT  IMPORTANT  IMPORTANT ***********************
  !  Below number of constituents in each group must be consistent with the AERO_DP (expChem) bundle
  !  in GEOSchem_GridComp/GOCART_GridComp/GOCART_GridCompMod.F90
  
  integer, parameter, public :: NUM_DUDP = 5, NUM_DUSV = 5, NUM_DUWT = 5, NUM_DUSD = 5
  integer, parameter, public :: NUM_BCDP = 2, NUM_BCSV = 2, NUM_BCWT = 2, NUM_BCSD = 2
  integer, parameter, public :: NUM_OCDP = 2, NUM_OCSV = 2, NUM_OCWT = 2, NUM_OCSD = 2
  integer, parameter, public :: NUM_SUDP = 1, NUM_SUSV = 1, NUM_SUWT = 1, NUM_SUSD = 1
  integer, parameter, public :: NUM_SSDP = 5, NUM_SSSV = 5, NUM_SSWT = 5, NUM_SSSD = 5 
  
  integer, public, parameter :: N_constit = 9         ! Number of constituents in snow
  
  ! (for riv, rin,aicev, aicen, and denice, instead use Teppei-defined 
  !  values below)
  !--------------------------------------------------------------------------------------------
  !    Spectrally integrated values for VIS and NIR using 
  !    the updated ice refractive indices by Warren and Brandt, (JGR, 2008).
  !    were used. 
  !    VIS: 300-690 nm; NIR: 690-3847 nm
  real, private, parameter :: RIV=0.018, RIN=0.017
  
  !--------------------------------------------------------------------------------------------
  ! So as to explain Abs. Co. for ice of 10 [m-1] by Kondo et al. (1988)
  ! the spectrally integrated Abs. Co. in VIS was used to get the one in NIR
  ! as a tuning parameter (see Yasunari et al., JGR, 2011).
  
  real, private, parameter :: DENICE=917.
  
  !--------------------------------------------------------------------------------------------
  
  !     Mass Absorption Coefficient or Mass Absorption Cross-section (MAC) [m2 g-1]
  !     Then the representative MACs for VIS & NIR was estimated from
  !     the GOCART/GEOS-5 optical properties.
  !     VIS: 300-690 nm; NIR: 690-3850 nm
  !     (Spectrally integrated with the standard surface radiation:
  !      ASTM G173-03 Tables; http://rredc.nrel.gov/solar/spectra/am1.5/)
  !     Updated on May 10, 2012
  
  ! ---------------------------------------------------------------------------
  !
  ! constants for snow constituents (dust, carbon, etc.)
  
  ! MAC, visible (VIS)
  real, private, parameter, dimension(N_constit) :: ABVIS = (/                              & 
       0.148,      &   ! Dust1 
       0.106,      &   ! Dust2
       0.076,      &   ! Dust3
       0.051,      &   ! Dust4
       0.032,      &   ! Dust5
       7.747,      &   ! Black carbon hydrophobic
       11.227,     &   ! Black carbon hydrophilic
       0.103,      &   ! Organic carbon hydrophobic
       0.114   /)      ! Organic carbon hydrophic
  
  ! MAC, near-infrared (NIR)
  real, private, parameter, dimension(N_constit) :: ABNIR = (/                              &
       0.095,      &   ! Dust1   	
       0.080,      &   ! Dust2
       0.062,      &   ! Dust3
       0.043,      &   ! Dust4
       0.029,      &   ! Dust5
       4.621,      &   ! Black carbon hydrophobic
       6.528,      &   ! Black carbon hydrophilic
       0.092,      &   ! Organic carbon hydrophobic
       0.127   /)      ! Organic carbon hydrophic
  
  !--------------------------------------------------------------------------------------------
  
  !     Scavenging coefficients for flushing effect in snow for constituents:
  !     Based on GOCART/GEOS-5 particle size;
  !     Tuning parameters so as to satisfy NCAR/CLM based scavenging efficiencies;
  !     See more in Yasunari et al. (SOLA, 2014)
  
  real, private, parameter, dimension(N_constit) :: SCAV = (/                               &
       0.065442,  &   ! Dust 1
       0.077829,  &   ! Dust 2
       0.306841,  &   ! Dust 3
       0.      ,  &   ! Dust 4
       0.      ,  &   ! Dust 5
       0.074361,  &   ! Black carbon hydrophobic
       0.502814,  &   ! Black carbon hydrophilic
       0.075855,  &   ! Organic carbon hydrophobic
       0.535225 /)    ! Organic carbon hydrophic
  
  !  Representative particle size in diameter 
  !  based on effective radius GOCART/GEOS-5 (dust 1-5 bins, BC, and OC) [um]

  real, private, parameter, dimension(N_constit) :: PSIZE = (/                              &
       1.272,     &   ! Dust 1
       2.649,     &   ! Dust 2
       4.602,     &   ! Dust 3
       8.334,     &   ! Dust 4
       15.341,    &   ! Dust 5
       0.078,     &   ! Black carbon hydrophobic
       0.148,     &   ! Black carbon hydrophilic
       0.175,     &   ! Organic carbon hydrophobic
       0.441 /)       ! Organic carbon hydrophic

  !============================================================================================
  
contains
  
  subroutine StieglitzSnow_snowrt(N_zones, N_snow, tileType,                     &  ! in 
       maxsndepth, rhofs, targetthick,                                           &  ! in 
       t1, area, tkgnd, precip, snowf, ts, dts, eturb, dedtc, hsturb, dhsdtc,    &  ! in 
       hlwtc, dhlwtc, raddn, zc1, totdepos,                                      &  ! in 
       wesn, htsnn, sndz, rconstit,                                              &  ! inout
       hlwout, fices, tpsn, rmelt,                                               &  ! out
       areasc, areasc0, pre, fhgnd, evap, shflux, lhflux, hcorr, ghfluxsno,      &  ! out
       sndzsc, wesnprec, sndzprec, sndz1perc,                                    &  ! out
       wesnperc, wesndens, wesnrepar, mltwtr,                                    &  ! out  
       excs, drho0, wesnbot, tksno, dtss        )                                   ! out
    
    !*********************************************************************
    ! AUTHORS:  M. Stieglitz, M. Suarez, R. Koster & S. Dery.
    ! VERSION:  2003b - This version last updated:  05/30/03.
    !*********
    ! INPUTS:
    !*********
    !  N_zones     : number of zones in the horizontal dimension (eg, 3 for Catchment, 1 for LandIce)
    !  N_snow      : number of snow layers  
    !  N_constit   : Number of constituent tracers in snow
    !  tileType    : surface type of the tile
    !  t1          : Temperature of catchment zones  [C]
    !  ts          : Air temperature [K]
    !  area        : Fraction of snow-free area in each catchment zone [0-1]
    !  precip      : Precipitation (Rain+snowfall) [kg/m^2/s == mm/s]
    !  snowf       : Snowfall per unit area of catchment [kg/m^2/s == mm/s]
    !  dts         : Time step  [s]
    !  eturb       : Evaporation per unit area of snow [kg/m^2/s == mm/s]
    !  dedtc       : d(eturb)/d(ts) [kg/m^2/s/K]
    !  hsturb      : Sensible heat flux per unit area of snow  [W/m^2]
    !  dhsdtc      : d(hsturb)/d(ts)  [W/m^2/K]
    !  hlwtc       : Emitted IR per unit area of snow  [W/m^2]
    !  dhlwtc      : d(hlwtc)/d(ts)  [W/m^2/K]
    !  raddn       : Net solar + incident terrestrial per unit area of snow [W/m^2]
    !  tkgnd       : Thermal diffusivity of soil in catchment zones [W/m/K]
    !  zc1         : Half-thickness (mid-point) of top soil layer [m]
    !***  Bin Zhao added *************************************************
    !  maxsndepth  :  Maximum snow depth beyond which snow gets thrown away
    !  rhofs       :  fresh snow density 
    !  targetthick :  the target thickness distribution relayer redistribute mass 
    !                and energy to; currently its value is surface type dependent           
    !                for catchment, the 1st array element the target thickness
    !                               the rest define a sigma distribution;
    !                for landice,   it is an array with specified thicknesses       
    !*********
    ! UPDATES:
    !*********
    !  wesn        : Layer water contents per unit area of catchment [kg/m^2]
    !  htsnn       : Layer heat contents relative to liquid water at 0 C [J/m^2]
    !  sndz        : Layer depths [m]
    !  rconstit    :  Mass of constituents in snow layer [kg] (i.e., [kg m-2])
    !  rmelt       : Flushed mass amount of constituents from the bottom snow layer [kg m-2 s-1 (kg/m^2/s)]
    !*********
    ! OUTPUTS: 
    !*********
    !  tpsn        : Layer temperatures [C]
    !  fices       : Layer frozen fraction [0-1]
    !  areasc      : Areal snow coverage at beginning of step [0-1]
    !  areasc0     : Areal snow coverage at end of step [0-1]
    !  pre         : Liquid water outflow from snow base [kg/m^2/s]
    !  fhgnd       : Heat flux at snow base at catchment zones  [W/m^2]
    !  hlwout      : Final emitted IR flux per unit area of snow [W/m^2]
    !  lhflux      : Final latent heat flux per unit area of snow [W/m^2]
    !  shflux      : Final sensible heat flux per unit area of snow   [W/m^2]
    !  evap        : Final evaporation per unit area of snow   [kg/m^2/s]
    !***  Bin Zhao added *************************************************
    !  sndzsc      :  top layer thickness change due to sublimation/condensation
    !  wesnprec    :  top layer water content change due to precip (different from precip itself)
    !  sndzprec    :  top layer thickness change due to precip 
    !  sndz1perc   :  top layer thickness change due to percolation
    !  wesnperc    :  layer water content change due to percolation
    !  wesndens    :  layer water content change due to densification
    !  wesnrepar   :  layer water content change due to relayer
    !  mltwtr      :  total melt water production rate
    !  excs        :  frozen part of water content from densification excess
    !  drho0       :  layer density change due to densification
    !  wesnbot     :  excessive water content due to thickness exceeding maximum depth
    !  tksno       :  layer conductivity
    !  dtss        :  top layer temperature change
    !
    !******************************************************************************
    ! NOTA:  By convention, wesn is representative for a catchment area
    !        equal to 1 whereas sndz is relative to the area covered by snow only.
    !******************************************************************************
    
    implicit none
    
    !      real, parameter :: lhv    = 2.4548E6 !  2.5008e6  !  @ 0 C [J/kg]
    !      real, parameter :: lhs    = 2.8368E6 !  2.8434e6  !  @ 0 C [J/kg]
    !      real, parameter :: lhf    = (lhs-lhv)             !  @ 0 C [J/kg]
    
    !rr      real, parameter :: cpw_liquid = 4185. ! [J/kg/K]
    
    !      real, parameter :: tfrz   = 273.16     !  @ 0 C [K]
    !      real, parameter :: rhofs  = 150.       !  [kg/m^3]
    !      real, parameter :: rhoma  = 500.       !  [kg/m^3]
    !      real, parameter :: rhow   = 1000.      !  [kg/m^3]
    !      real, parameter :: wemin  = 13.        !  [kg/m^2]
    
    real, parameter :: snfr   = 0.01       !  holding capacity
    real, parameter :: small  = 1.e-6      !  small number 

    integer,                               intent(in)    :: N_zones, N_snow, tileType

    real,    dimension(N_zones),           intent(in)    :: t1, area, tkgnd
    real,    dimension(N_constit),         intent(in)    :: totdepos
    real,                                  intent(in)    :: ts, precip, snowf, dts, dedtc, raddn, hlwtc
    real,                                  intent(in )   :: dhsdtc, dhlwtc, eturb, hsturb, zc1

    real,    dimension(N_snow),            intent(inout) :: wesn, htsnn, sndz
    real,    dimension(N_snow, N_constit), intent(inout) :: rconstit

    real,    dimension(N_snow),            intent(out)   :: tpsn, fices
    real,    dimension(N_zones),           intent(out)   :: fhgnd
    real,                                  intent(out)   :: hlwout, lhflux, shflux, areasc0, evap, areasc, pre
    real,                                  intent(out)   :: hcorr
    real,    dimension(N_constit),         intent(out)   :: rmelt
    real,                                  intent(out)   :: ghfluxsno
    
    real,                                  intent(out)   :: wesnprec
    real,                                  intent(out)   :: sndzsc, sndzprec
    real,                                  intent(out)   :: sndz1perc
    real,    dimension(N_snow),            intent(out)   :: wesnperc
    real,    dimension(N_snow),            intent(out)   :: wesndens
    real,    dimension(N_snow),            intent(out)   :: wesnrepar
    real,                                  intent(out)   :: mltwtr
    real,    dimension(N_snow),            intent(out)   :: excs
    real,    dimension(N_snow),            intent(out)   :: drho0
    real,                                  intent(out)   :: wesnbot
    real,    dimension(N_snow),            intent(out)   :: tksno
    real,                                  intent(out)   :: dtss

    real,                                  intent(in)    :: maxsndepth
    real,                                  intent(in)    :: rhofs
    real,    dimension(N_snow),            intent(in)    :: targetthick

    ! ----------------------------------------
    !
    ! Locals
    
    real :: tsx, mass,snowd,rainf,denom,alhv,lhturb,dlhdtc,                 &
         enew,eold,tdum,fnew,tnew,icedens,densfac,hnew,scale,t1ave,         &
         flxnet,fdum,dw,waterin,waterout,snowin,snowout, mtwt,              &
         waterbal,precision,flow,term,dz,w(0:N_snow),HTSPRIME,              &
         wlossfrac
    real :: excsdz, excswe, sndzsum, mtwt0
    
    real, dimension(size(wesn)  ) :: cmpc,dens
    real, dimension(size(wesn)  ) :: tksn
    real, dimension(size(wesn)  ) :: dtc,q,cl,cd,cr
    real, dimension(size(wesn)+1) :: fhsn,df
    real, dimension(size(wesn)  ) :: htest,ttest,ftest
    

    ! by Teppei --- GOSWIM related variables
    !============================================================================================
    
    real, dimension(size(wesn)  ) :: denblk          ! bulk snow density
    real, dimension(size(wesn)  ) :: po              ! snow porosity
    
    !============================================================================================

    logical, dimension(size(wesn)  ) :: ice1, tzero
    real,    dimension(size(wesn)  ) :: dens0
    real,    dimension(N_constit)    :: flow_r,rconc
    
    integer :: i,izone,k
    logical :: logdum
    
    snowd = sum(wesn)
    snowin = snowd
    ghfluxsno = 0.
    
    !rr   correction for "cold" snow
    tsx   = min(ts-tf,0.)*StieglitzSnow_CPW
    
    !rr   correction for heat content of rain
    !rr       tsx_rain = max(ts-tf,0.)*cpw_liquid
    
    df     = 0.
    dtc    = 0.
    tpsn   = 0.
    fices  = 0.
    areasc = 0.
    areasc0= 0.
    pre    = 0.
    fhgnd  = 0.
    hlwout = 0.
    shflux = 0.
    lhflux = 0.
    evap   = 0.
    excs   = 0.
    hcorr  = 0.
    dens   = rhofs
    rainf  = precip - snowf   ! [kg/m^2/s]
    
    sndzsc    = 0.
    wesnprec  = 0.
    sndzprec  = 0.
    sndz1perc = 0.
    wesnperc  = 0.
    wesndens  = 0.
    wesnrepar = 0.
    wesnbot   = 0.
    dtss      = 0. 
    excswe    = 0.
    
    rmelt  = 0.0
    mltwtr = 0.0
    drho0  = 0.0
    tksno  = 0.0
    
    if(snowd <= StieglitzSnow_MINSWE) then ! initial snow mass is negligible
       
       ! Melt off initial (very small) snowpack; new snow pack is based
       !   on new snowfall only (if any)
       
       call StieglitzSnow_calc_asnow( snowd, areasc )
       areasc0 = 0.
       pre = snowd/dts + areasc*rainf       ! pre = melted snowpack plus rainfall
       wesn  = 0.
       hcorr = hcorr + raddn*areasc + sum(htsnn)/dts
       htsnn = 0.
       sndz  = 0.
       mltwtr = snowd/dts                   ! mltwtr = melted snowpack
       do k=1,N_constit
          rmelt(k)=sum(rconstit(:,k))/dts
       enddo
       rconstit(:,:) = 0.
       
       if(snowf > 0.) then  ! only initialize with non-liquid part of precip
                            ! liquid precip (rainf) is part of outflow from snow base (see "pre" above)
          
          wesn    = snowf*dts/float(N_snow)  
          htsnn   = (tsx-alhm)*wesn
          call StieglitzSnow_calc_asnow( snowf*dts, areasc0 )

          !*** should have fractional snow cover taken into account
          sndz = wesn/(max(areasc0,small)*rhofs)

          hcorr  = hcorr - tsx*snowf           ! randy
          
          ! Add constituent to top snow layer, in area covered by snow.
          do k=1,N_constit
             rconstit(1,k)=rconstit(1,k)+areasc0*totdepos(k)*dts
          enddo

          ! call relayer without heat content adjustment

          call StieglitzSnow_relayer( N_snow, N_constit, tileType, targetthick, &
               htsnn, wesn, sndz, rconstit )
          
          call StieglitzSnow_calc_tpsnow(N_snow, htsnn, wesn, tpsn, fices)
          
       endif   ! (snowf > 0.)
       
       return  ! if there was no snow at start of time step
       
    endif      ! (snowd <= StieglitzSnow_MINSWE)
    
    ! ---------------------------------------------------------------
    !
    ! derive new snow pack from existing snow pack and new snowfall: 
    
    !**** Determine the fractional snow coverage
    
    call StieglitzSnow_calc_asnow( snowd, areasc )
    
    !**** Set the mean density & diffusivity of the layers
    
    do i=1,N_snow
       if(sndz(i) > 0) dens(i) = max(wesn(i)/(areasc*sndz(i)),rhofs)
    enddo
    tksn  = 3.2217e-06*dens**2
    tksno = tksn
    dens0 = dens
    
    !**** Determine temperature & frozen fraction of snow layers
    
    call StieglitzSnow_calc_tpsnow(N_snow, htsnn, wesn, tpsn, fices)
    
    mtwt  = sum(wesn*(1.-fices)) 
    
    !**** Calculate the ground-snow energy flux at 3 zones
    
    denom = 1./(sndz(N_snow)*0.5-zc1)
    fhgnd = -sqrt(tkgnd*tksn(N_snow))*area*denom*(tpsn(N_snow)-t1)
    fhsn(N_snow+1) = sum(fhgnd)
    do i=1,N_zones
       df(N_snow+1)=df(N_snow+1)-sqrt(tkgnd(i)*tksn(N_snow))*area(i)*denom
    enddo
    
    !**** Ensure against excessive heat flux between ground and snow:
    !**** if heat flux found to cause the lowest snow layer temperature
    !**** to "overshoot" (e.g. to become higher than the ground temperature
    !**** when it had been lower), reduce the heat flux.  If the lowest 
    !**** snow layer starts off at zero and the new temperature is greater
    !**** than zero, reduce the heat flux to melt only half of the lowest
    !**** layer snow.

    t1ave=sum(t1*area)/sum(area)
    htest=htsnn
    htest(N_snow)=htest(N_snow)+fhsn(N_snow+1)*dts*areasc
    
    call StieglitzSnow_calc_tpsnow(N_snow, htest, wesn, ttest, ftest)
    
    scale=1.
    if((t1ave-tpsn(N_snow))*(t1ave-ttest(N_snow)) .lt. 0.) then
       scale=0.5*(tpsn(N_snow)-t1ave)/(tpsn(N_snow)-ttest(N_snow))
    endif
    if(tpsn(N_snow) .eq. 0. .and. ttest(N_snow) .gt. 0. .and.                &
         abs(fhsn(N_snow+1)) .gt. 1.e-10) then
       scale=(-0.5*htsnn(N_snow)/(dts*areasc))/fhsn(N_snow+1)
    endif
    
    fhsn(N_snow+1)=fhsn(N_snow+1)*scale
    df(N_snow+1)=df(N_snow+1)*scale
    fhgnd=fhgnd*scale
    
    !**** Calculate heat fluxes between snow layers.
    
    do i=2,N_snow
       df(i) =  -sqrt(tksn(i-1)*tksn(i))/((sndz(i-1)+sndz(i))*0.5)
       fhsn(i)= df(i)*(tpsn(i-1)-tpsn(i))
    enddo
    
    ghfluxsno = fhsn(2)
    
    !**** Effective heat of vaporization includes bringing snow to 0 C
    
    alhv   = alhe + alhm                            ! randy
    
    !**** Initial estimate of latent heat flux change with Tc
    
    lhturb = alhv*eturb
    dlhdtc = alhv*dedtc
    
    !**** Initial estimate of net surface flux & its change with Tc

    fhsn(1) = lhturb + hsturb + hlwtc - raddn
    df(1)   = -(dlhdtc + dhsdtc + dhlwtc)
    
    !**** Prepare array elements for solution & coefficient matrices.
    !**** Terms are as follows:  left (cl), central (cd) & right (cr)
    !**** diagonal terms in coefficient matrix & solution (q) terms.
    
    do i=1,N_snow
       
       call StieglitzSnow_calc_tpsnow(htsnn(i),wesn(i),tdum,fdum, ice1(i),tzero(i), .true.)
       
       if(ice1(i)) then
          cl(i) = df(i)
          cd(i) = StieglitzSnow_CPW*wesn(i)/dts - df(i) - df(i+1)
          cr(i) = df(i+1)
          q(i)  = fhsn(i+1)-fhsn(i)
       else
          cl(i) = 0.
          cd(i) = 1.
          cr(i) = 0.
          q(i)  = 0.
       endif
       
    enddo
    
    cl(1)      = 0.
    cr(N_snow) = 0.
    
    do i=1,N_snow-1
       if(.not.ice1(i)) cl(i+1) = 0.
    enddo
    
    do i=2,N_snow
       if(.not.ice1(i)) cr(i-1) = 0.
    enddo
    
    
    !**** Solve the tri-diagonal matrix for implicit change in Tc.
    
    call STIEGLITZSNOW_TRID(dtc,cl,cd,cr,q,N_snow)
    
    !**** Check temperature changes for passages across critical points,i.e.
    !**** If implicit change has taken layer past melting/freezing, correct.
    
    do i=1,N_snow
       if(tpsn(i)+dtc(i) > 0. .or. htsnn(i)+wesn(i)*StieglitzSnow_CPW*dtc(i) > 0.) then
          dtc(i)=-tpsn(i)
       endif
       if(.not.ice1(i)) dtc(i)=0.
    enddo
    
    !**** Further adjustments; compute new values of h associated with
    !**** all adjustments.
    
    eold=sum(htsnn)
    
    do i=1,N_snow
       
       !**** Quick check for "impossible" condition:
       
       if(.not.tzero(i) .and. .not.ice1(i)) then
          write(*,*) 'bad snow condition: fice,tpsn =',fices(i),tpsn(i)
          stop
       endif
       
       !****  Condition 1: layer starts fully frozen (temp < 0.)
       
       if(.not.tzero(i)) then
          tnew=tpsn(i)+dtc(i)
          fnew=1.
          
       endif
       
       !****  Condition 2: layer starts with temp = 0, fices < 1.
       !      Corrections for flxnet calculation: Koster, March 18, 2003.
       
       if(.not.ice1(i)) then
          tnew=0.
          if(i==1) flxnet= fhsn(i+1)+df(i+1)*(dtc(i)-dtc(i+1))              &
               -fhsn(i)-df(i)*dtc(i)
          if(i > 1 .and. i < N_snow) flxnet=                                &
               fhsn(i+1)+df(i+1)*(dtc(i)-dtc(i+1))                          &
               -fhsn(i)-df(i)*(dtc(i-1)-dtc(i))
          if(i==N_snow) flxnet=fhsn(i+1)+df(i+1)*dtc(i)                     &
               -fhsn(i)-df(i)*(dtc(i-1)-dtc(i))
          HTSPRIME=HTSNN(I)+AREASC*FLXNET*DTS
          call StieglitzSnow_calc_tpsnow( HTSPRIME, wesn(i), tdum, fnew, logdum, logdum, .true. )
          fnew=amax1(0.,  amin1(1.,  fnew))
          
       endif
       
       !****  Condition 3: layer starts with temp = 0, fices = 1.
       !      Corrections for flxnet calculation: Koster, March 18, 2003.
       
       if(ice1(i) .and. tzero(i)) then
          if(dtc(i) < 0.) then
             tnew=tpsn(i)+dtc(i)
             fnew=1.
          endif
          if(dtc(i) >= 0.) then
             tnew=0.
             if(i==1) flxnet=fhsn(i+1)+df(i+1)*(dtc(i)-dtc(i+1))            &
                  -fhsn(i)-df(i)*dtc(i)
             if(i > 1 .and. i < N_snow) flxnet=                             &
                  fhsn(i+1)+df(i+1)*(dtc(i)-dtc(i+1))                       &
                  -fhsn(i)-df(i)*(dtc(i-1)-dtc(i))
             if(i==N_snow) flxnet=fhsn(i+1)+df(i+1)*dtc(i)                  &
                  -fhsn(i)-df(i)*(dtc(i-1)-dtc(i))
             
             HTSPRIME=HTSNN(I)+AREASC*FLXNET*DTS
             call StieglitzSnow_calc_tpsnow( HTSPRIME, wesn(i), tdum, fnew, logdum, logdum, .true. )
             fnew=amax1(0.,  amin1(1.,  fnew))
          endif
       endif
       
       !**** Now update heat fluxes & compute sublimation or deposition.
       !**** Note: constituents (dust, etc.) do not evaporate.
       
       if(i == 1) then
          dtss   = dtc(1)
          lhflux = lhturb + dlhdtc*dtc(1)
          shflux = hsturb + dhsdtc*dtc(1)
          hlwout = hlwtc  + dhlwtc*dtc(1)
          evap = lhflux/alhv
          dw = -evap*dts*areasc
          if(-dw > wesn(1) ) then
             dw = -wesn(1)
             evap = -dw/(dts*areasc)
             hcorr=hcorr+(lhflux-evap*alhv)*areasc
             lhflux=evap*alhv
          endif
          wesn(1)  = wesn(1) + dw
          denom = 1./dens(1)
          if(dw > 0.) denom = 1./StieglitzSnow_RHOMA
          sndz(1) = sndz(1) + dw*denom
          sndzsc = dw*denom
       endif
       
       if(i == N_snow) then
          do izone=1,N_zones
             fhgnd(izone)=fhgnd(izone)+area(izone)*df(N_snow+1)*dtc(N_snow)
          enddo
       endif
       
       !**** Now update thermodynamic quantities.
       
       htsnn(i)=(StieglitzSnow_CPW*tnew-fnew*alhm)*wesn(i)
       tpsn(i) = tnew    
       fices(i)= fnew

    enddo     ! (i=1,N_snow)
    
    ! -----------------------------------------------------------------------------
    !
    !**** Store excess heat in hcorr.
    
    enew=sum(htsnn)
    hcorr=hcorr-((enew-eold)/dts+areasc*(lhflux+shflux+hlwout-raddn)        &
         -areasc*(fhsn(N_snow+1)+df(N_snow+1)*dtc(N_snow))                  &
         )
    
    call StieglitzSnow_calc_tpsnow(N_snow, htsnn, wesn, tpsn, fices)
    
    mltwtr = max(0., sum(wesn*(1.-fices)) - mtwt)
    mltwtr = mltwtr / dts
    
    mtwt0 = sum(wesn*fices)
    
    !rr!**** Add rainwater and snow at ts., bal. budget with shflux.
    !rr   (tried and failed 19 Jun 2003, reichle)
    !rr
    !rr       wesn (1) = wesn (1) + (rainf*areasc+snowf)*dts
    !rr       htsnn(1) = htsnn(1) + (tsx -alhm)*(snowf*dts) + tsx_rain*rainf*dts
    !rr       sndz (1) = sndz (1) + (snowf/rhofs)*dts
    !rr       !  shflux   = shflux   + tsx*snowf                    ! randy
    !rr       hcorr   = hcorr   - (tsx-alhm)*snowf - tsx_rain*rainf ! randy
    
    !**** Add rainwater at 0 C, snow at ts., bal. budget with shflux.
    
    wesn (1) = wesn (1) + (rainf*areasc+snowf)*dts
    htsnn(1) = htsnn(1) + (tsx -alhm)*(snowf*dts)
    sndz (1) = sndz (1) + (snowf/rhofs)*dts
    wesnprec = (rainf*areasc+snowf)*dts
    sndzprec = (snowf/rhofs)*dts
    hcorr    = hcorr   -  tsx*snowf          ! randy
    
    snowd=sum(wesn)
    
    call StieglitzSnow_calc_tpsnow(N_snow, htsnn, wesn, tpsn, fices)
    
    !**** Constituent deposition: Add to top snow layer, in area covered by snow.
    do k=1,N_constit
       rconstit(1,k)=rconstit(1,k)+areasc*totdepos(k)*dts
    enddo

    !**** Move meltwater through the pack.
    !**** Updated by Koster, August 27, 2002.
    
    pre = 0.
    rmelt(:) = 0.
    flow = 0.
    flow_r(:) = 0.
    
    wesnperc = wesn
    
    do i=1,N_snow
       
       if(flow > 0.) then
          wesn (i) =  wesn(i) + flow               ! add "flow" [kg/m2] from layer i-1 to wesn(i)
          do k=1,N_constit
             rconstit(i,k)=rconstit(i,k)+flow_r(k)
          enddo
          call StieglitzSnow_calc_tpsnow(N_snow, htsnn, wesn, tpsn, fices)  
       endif
       
       pre  = max((1.-fices(i))*wesn(i), 0.)
       flow = 0.
       flow_r(:) = 0.
       rconc(:) = 0.
       
       if(snowd > wemin) then
          
          icedens = wesn(i)*fices(i)/(sndz(i)+1.e-20)
          densfac = amax1(0., amin1(1., icedens/rhofs))
          term    = densfac*snfr*(sndz(i)*rhow-wesn(i)*fices(i))
          
          if(pre > term) then
             pre = min(pre - term, wesn(i))        ! when asnow=1, retain some liquid water in snow pack
             do k=1,N_constit
                rconc(k)=rconstit(i,k)/wesn(i)
             enddo
             wesn(i) = wesn(i) - pre
             flow = pre
          endif
       else
          do k=1,N_constit
             rconc(k)=rconstit(i,k)/wesn(i)         
          enddo
          wesn(i) = wesn(i) - pre                  ! when asnow<1, remove all liquid water from snow pack
          flow = pre
       endif
       
       !------------------------------- by Teppei ---------------------------------------

       if(areasc.ge.0.999) then
          
          do k=1,N_constit
             ! mass loss by excess water
             ! To calculate bulk snow density in snow layers
             
             denblk(i)=( wesn(i)/(areasc+1.e-20) )/(sndz(i)+1.e-20)
             
             ! porosity of snow layers
             po(i)=-7.2E-07*(denblk(i)**2.0)-(0.00063*denblk(i))+0.967073
             if(denblk(i) > 800.) po(i)=0.
             
             ! constituents flushing
             flow_r(k)=(po(i)*(1.0-(min(PSIZE(k),5.0)/5.0))*SCAV(k))             &
                  *(rconc(k)+1.0e-20)*flow
             if ( (flow < 1.0e-20).or.(denblk(i) > 800.) ) flow_r(k)=0. 
             
             rconstit(i,k)=rconstit(i,k)-flow_r(k)
             rconstit(i,k)=amax1(0.,rconstit(i,k)) ! guard against truncation error
             
          enddo
          
       endif
       
       if(areasc.lt.0.999) then
          
          do k=1,N_constit
             flow_r(k)=rconc(k)*flow
             rconstit(i,k)=rconstit(i,k)-flow_r(k)
             rconstit(i,k)=amax1(0.,rconstit(i,k)) ! guard against truncation error
          enddo
          
       endif

       !---------------------------------------------------------------------------------
       
       !**** Adjust top layer snow depth to get proper density values
       !**** But limit this change for large throughflow (STEPH 06/19/03)
       
       if (i==1) then
          dz=min(flow/dens(i),0.5*sndz(i))
          sndz(i)=sndz(i)-dz
          sndz1perc = -dz
       endif

    enddo  ! (i=1,N_snow)
    
    ! ----------------------------------------------------------------------------------------

    wesnperc = wesn - wesnperc
    
    pre = flow/dts                                 ! convert outflow to flux units [kg/m2/s]
    do k=1,N_constit
       rmelt(k)=rmelt(k)+flow_r(k)/dts
    enddo
    snowd=sum(wesn)
    
    !**** Update snow density by compaction (Pitman et al. 1991)
    
    mass = 0.
    w    = 0.
    
    wesndens = wesn
    
    if(snowd > wemin) then ! Compaction only after full coverage.
       
       do i=1,N_snow
          dens(i) = rhofs
          if(sndz(i)>0.) dens(i) = max(wesn(i)/(sndz(i)),rhofs)
       enddo
       
       drho0 = dens 
       
       cmpc    = exp(14.643 - (4000./min(tpsn+tf,tf))-.02*dens)
       
       do i=1,N_snow

          w(i)    = wesn(i)
          mass    = mass + 0.5*(w(i)+w(i-1))
          dens(i) = dens(i)*(1. + (dts*0.5e-7*9.81)*mass*cmpc(i))
          
          !**** Clip densities below maximum value, adjust quantities accordingly
          !**** while conserving heat & mass (STEPH 06/21/03).
          
          if(dens(i) > StieglitzSnow_RHOMA) then
             
             if (tileType==MAPL_LANDICE) then               ! restrict SWE adjustment to LANDICE tiles
                
                ! excs = SWE in excess of max density given fixed snow depth

                excs(i) = (dens(i)-StieglitzSnow_RHOMA)*sndz(i)           ! solid + liquid
                wlossfrac=excs(i)/wesn(i)
                wesn(i) = wesn(i) - excs(i)                               ! remove EXCS from SWE
                do k=1,N_constit
                   rmelt(k)=rmelt(k)+rconstit(i,k)*wlossfrac/dts
                   rconstit(i,k)=rconstit(i,k)*(1.-wlossfrac)
                   rconstit(i,k)=amax1(0.,rconstit(i,k))                  ! guard against truncation error
                enddo
                hnew = (StieglitzSnow_CPW*tpsn(i)-fices(i)*alhm)*wesn(i)  ! adjust heat content accordingly
                hcorr= hcorr+(htsnn(i)-hnew)/dts                          ! add excess heat content into residual accounting term
                htsnn(i)= hnew

             end if
             
             dens(i) = StieglitzSnow_RHOMA
          endif
       enddo
       drho0 = dens - drho0
    endif
    
    wesndens = wesn - wesndens
    
    if (tileType==MAPL_LANDICE) then                        ! finish SWE adjustment for LANDICE tiles
       
       pre  = pre + sum(excs*max(1.-fices,0.0))/dts
       excs = excs * fices / dts
       
    end if

    snowd=sum(wesn)
    call StieglitzSnow_calc_asnow( snowd, areasc0 )
    areasc0 = max(small, areasc0 )
    sndz = (wesn/areasc0)/dens
    
    sndzsum = sum(sndz)
    if(sndzsum > maxsndepth) then
       excsdz  = sndzsum - maxsndepth
       excswe  = dens(N_snow) * excsdz 
       wlossfrac=excswe/wesn(N_snow)
       wesn(N_snow) = wesn(N_snow) - excswe
       do k=1,N_constit
          rmelt(k)=rmelt(k)+rconstit(N_snow,k)*wlossfrac/dts
          rconstit(N_snow,k)=rconstit(N_snow,k)*(1.-wlossfrac)
          rconstit(N_snow,k)=amax1(0.,rconstit(N_snow,k)) ! guard against truncation error
       enddo
       hnew = (StieglitzSnow_CPW*tpsn(N_snow)-fices(N_snow)*alhm)*wesn(N_snow)
       htsnn(N_snow)= hnew
       sndz(N_snow) = sndz(N_snow) - excsdz
       wesnbot = excswe
    endif
    
    !**** Restore layers to sigma values.
    
    wesnrepar = wesn
        
    ! call relayer with adjustment of heat content and hcorr accounting

    call StieglitzSnow_relayer( N_snow, N_constit, tileType, targetthick, &
         htsnn, wesn, sndz, rconstit, tpsn, fices, dts, hcorr  )
    
    wesnrepar = wesn - wesnrepar
    
    !**** Reset fractional area coverage.
    
    call StieglitzSnow_calc_asnow( sum(wesn), areasc0 )
    
    !**** Final check for water balance.
    
    waterin   = (rainf*areasc+snowf)*dts + max(dw,0.)
    waterout  = pre*dts - min(dw,0.)
    snowout   = sum(wesn) + sum(excs) + excswe
    waterbal  = snowin + waterin - waterout - snowout
    precision = snowout*small
    
#if 0
    if((waterbal > precision).and.(waterbal > small) .or. pre < -1.e-13  ) then
       write(*,*) 'Warning: Imbalance in snow water budget!', waterbal
       write(*,*) 'waterin    = ', waterin
       write(*,*) 'snowin     = ', snowin
       write(*,*) 'waterout   = ', waterout
       write(*,*) 'snowout    = ', snowout
       write(*,*) 'dw         = ', dw
       write(*,*) 'excswe     = ', excswe
       write(*,*) 'sum(excs)  = ', sum(excs)
       write(*,*) 'snowf*dts  = ', snowf*dts
       write(*,*) 'sum(wesn)  = ', sum(wesn)
       write(*,*) (wesn(i), i=1,N_snow)
       write(*,*) 'sum(sndz)  = ', sum(sndz)
       write(*,*) (sndz(i), i=1,N_snow)
       write(*,*) 'dens0      = '
       write(*,*) (dens0(i), i=1,N_snow)
       !write(*,*) 'sum(wesn0) = ', sum(wesn0)
       !write(*,*) (wesn0(i), i=1,N_snow)
       write(*,*) 'sum(wesn1) = ', sum(wesn1)
       write(*,*) (wesn1(i), i=1,N_snow)
       write(*,*) 'sum(wesn2) = ', sum(wesn2)
       write(*,*) (wesn2(i), i=1,N_snow)
       !write(*,*) 'sum(sndz0) = ', sum(sndz0)
       !write(*,*) (sndz0(i), i=1,N_snow)
       write(*,*) 'sum(sndz1) = ', sum(sndz1)
       write(*,*) (sndz1(i), i=1,N_snow)
       write(*,*) 'sum(sndz2) = ', sum(sndz2)
       write(*,*) (sndz2(i), i=1,N_snow)
       !stop
    endif
#endif
    
    return  !  end snow
    
  end subroutine StieglitzSnow_snowrt
  
  ! **********************************************************************
  
  subroutine FindTargetThickDist_Landice(N_snow, sndz, dzmax, topthick, thickdist)
    
    ! get snow layer target thicknesses to be used with relayer for *landice*

    integer, intent(in)                       :: N_snow
    real,    intent(in)                       :: sndz(N_snow)
    real,    intent(in)                       :: dzmax(N_snow)
    real,    intent(out)                      :: topthick
    real,    intent(out), dimension(N_snow-1) :: thickdist
    
    real,                 dimension(N_snow)   :: sndzt
    real                                      :: totald, dzdiff, restthick
    integer                                   :: i
    integer,              dimension(N_snow)   :: mark
    logical                                   :: lth_satisfy
    
    totald = sum(sndz)
    sndzt = totald/float(N_snow)
    
    if (sndzt(1) < dzmax(1)) then
       
       topthick = dzmax(1)
       do i=2,N_snow
          thickdist(i-1) = 1.0/real(N_snow-1,kind=4)
       enddo
       
    else
       
       mark = 0
       do
          lth_satisfy = .true.
          do i=1,N_snow
             if(mark(i) == 0 .and. sndzt(i) > dzmax(i)) then
                sndzt(i) = dzmax(i)
                mark(i)  = 1
                lth_satisfy = .false.
             endif
          enddo
          if(lth_satisfy) exit
          dzdiff = 0.0
          do i=1,N_snow
             if(mark(i) == 1) then
                dzdiff = dzdiff + sndzt(i)
             endif
          enddo
          restthick = (totald-dzdiff)/float(N_snow-sum(mark))
          do i=1,N_snow
             if(mark(i) == 0) then
                sndzt(i) = restthick
             endif
          enddo
       enddo
          
       topthick = sndzt(1)
       totald = totald - topthick
       do i=2,N_snow
          thickdist(i-1) = sndzt(i)/totald
       enddo
          
    endif
       
    return
       
  end subroutine FindTargetThickDist_Landice
  
  ! **********************************************************************
  
  subroutine StieglitzSnow_relayer(N_snow, N_constit, tileType, targetthick, &
       htsnn, wesn, sndz, rconstit, tpsn, fices, dts, hcorr )
    
    ! relayer for land and landice tiles

    ! revised to included processing of target thickness parameters and 
    !   optional snow heat content adjustment
    !   
    ! optional arguments        action
    ! -----------------------------------------------------------
    ! none                      original relayer() (redistribution only)
    ! tpsn, fices               + adjust heat content (originally done externally)
    ! tpsn, fices, dts, hcorr   + account for heat content adjustment in correction term

    implicit none
    
    integer, intent(in)                                 :: N_snow, N_constit, tileType

    real,    intent(in),    dimension(N_snow)           :: targetthick

    real,    intent(inout), dimension(N_snow)           :: htsnn, wesn, sndz
    real,    intent(inout), dimension(N_snow,N_constit) :: rconstit
    
    real,    intent(out),   dimension(N_snow), optional :: tpsn, fices

    real,    intent(in),                       optional :: dts
    real,    intent(inout),                    optional :: hcorr

    ! ----------------------------
    !
    ! local variables:
    
    character(len=*), parameter              :: Iam = 'StieglitzSnow_relayer'

    real                                     :: thick_toplayer
    real,    dimension(N_snow-1)             :: thickdist

    real,    dimension(N_snow,  2+N_Constit) :: h, s
    
    integer                                  :: i, k, ilow, ihigh
    
    real                                     :: dz, hnew
    real                                     :: totalthick, tdum, fdum
    real,    dimension(N_snow)               :: tol_old, bol_old, tol_new, bol_new
    real,    dimension(N_snow)               :: thickness

    logical                                  :: adjust_htsnn, update_hcorr, kflag

    logical, dimension(N_snow)               :: ice10, tzero0
    
    !**** thickness(1) : final thickness of topmost snow layer (m)
    !**** h            : array holding specific heat, water, and constituent contents
    !**** s            : array holding the total and final heat, water, and constit. contents
    !**** ilow         : first layer used in a particular relayering calculation
    !**** ihigh        : final layer used in a particular relayering calculation 
    !**** totalthick   : total thickness of layers 2 through N_snow
    !**** thickness    : array holding final thicknesses (m) of the snow layers
    !**** tol_old(i)   : depth (from surface) of the top of layer i, before          &
    !****                  relayering
    !**** bol_old(i)   : depth (from surface) of the bottom of layer i, before       &
    !****                  relayering
    !**** tol_old(i)   : depth (from surface) of the top of layer i, after           &
    !****                  relayering
    !**** bol_old(i)   : depth (from surface) of the bottom of layer i, after        &
    !****                  relayering
   
    ! ---------------------------------------
    !
    ! process optional arguments (required to maintain 0-diff; reichle 13 Oct 2023)
    
    if     ( present(tpsn) .and. present(fices) ) then
       
       adjust_htsnn = .true.
       
    elseif ( present(tpsn) .or.  present(fices) ) then
       
       write(*,*) Iam, '(): bad optional arguments (tpsn, fices)'
       stop
       
    else
       
       adjust_htsnn = .false.
       
    end if
    
    if (adjust_htsnn) then
       
       if     ( present(dts) .and. present(hcorr) ) then
          
          update_hcorr = .true.
          
       elseif ( present(dts) .or.  present(hcorr) ) then
          
          write(*,*) Iam, '(): bad optional arguments (dts, hcorr)'
          stop
          
       else
          
          update_hcorr = .false.
          
       end if
       
       ! determine frozen fraction and temperature before relayering
       
       do i=1,N_snow
          call StieglitzSnow_calc_tpsnow(htsnn(i),wesn(i),tdum,fdum,ice10(i),tzero0(i), .true. )
       enddo

    end if
    
    ! process "targetthick" snow depth parameters:
    !
    ! targetthick:     contents depend on tileType (Catch[CN] or Landice)
    !
    ! thick_toplayer:  *target* thickness of top snow layer (m);
    !                    NOTE: thickness(1) [see below] is final thickness of top layer (m)
    !
    ! thickdist:       assigned (final) distribution of thickness in layers 2:N_snow,
    !                    in terms of fraction
    
    select case (tileType)
    case (MAPL_LANDICE)
       call FindTargetThickDist_Landice(N_snow, sndz, targetthick, thick_toplayer, thickdist)
    case default
       thick_toplayer  = targetthick(1)
       thickdist       = targetthick(2:N_snow)    
    end select

    ! ----------------------------------------------------------------------------------------
    !
    ! start of original relayer()

    totalthick   = sum(sndz)     ! total snow depth

    ! make sure thickness of top layer does not exceed total thickness

    thickness(1) = amin1(totalthick*0.9, thick_toplayer)

    totalthick   = totalthick-thickness(1)

    do i=1,N_snow-1
       thickness(i+1)=thickdist(i)*totalthick
    enddo
    
    !**** Initialize some variables.
    
    h  = 0.
    s  = 0.
    dz = 0.
    
    !**** Compute specific heat & water contents of old layers.
    
    do i=1,N_snow
       if (sndz(i) > 0.) then
          h(i,1) = htsnn(i)/sndz(i)
          h(i,2) =  wesn(i)/sndz(i)
          do k=1,N_Constit
             h(i,2+k)=rconstit(i,k)/sndz(i)
          enddo
       endif
    enddo
    
    !**** Determine old and new boundary depths (cumulative from top)
    !**** (tol refers to "top of layer", bol refers to "bottom of layer"
    
    tol_old(1)=0.
    bol_old(1)=sndz(1)
    tol_new(1)=0.
    bol_new(1)=thickness(1)
    
    do i=2,N_snow
       tol_old(i)=bol_old(i-1)
       bol_old(i)=bol_old(i-1)+sndz(i)
       tol_new(i)=bol_new(i-1)
       bol_new(i)=bol_new(i-1)+thickness(i)
    enddo
    
    !**** Redistribute quantities
    
    !**** Step 1: Do top layer
    ihigh=1
    do k=1,N_snow
       if(bol_old(k) .lt. bol_new(1)) ihigh=k+1
    enddo
    
    do k=1,ihigh
       if(k .lt. ihigh) dz=sndz(k)
       if(k .eq. ihigh) dz=bol_new(1)-tol_old(k)
       s(1,:)=s(1,:)+h(k,:)*dz
    enddo
    
    !**** Step 2: Do remaining layers
    do i=2,N_snow
       
       ilow=ihigh
       do k=ilow,N_snow
          if(bol_old(k) .lt. bol_new(i)) ihigh=k+1
       enddo
       
       if(ihigh .eq. N_snow+1)  ihigh=N_snow ! Account for potential truncation problem 
       
       do k=ilow,ihigh
          if(k .eq. ilow .and. k .lt. ihigh) dz=bol_old(k)-tol_new(i)
          if(k .eq. ilow .and. k .eq. ihigh) dz=bol_new(i)-tol_new(i)
          if(k .gt. ilow .and. k .lt. ihigh) dz=bol_old(k)-tol_old(k)
          if(k .gt. ilow .and. k .eq. ihigh) dz=bol_new(i)-tol_old(k)
          s(i,:)=s(i,:)+h(k,:)*dz
       enddo
       
    enddo
    
    htsnn = s(:,1)
    wesn  = s(:,2)
    do k=1,N_Constit
       rconstit(:,k)=s(:,2+k)
    enddo
    sndz=thickness
    
    ! end of original relayer()
    !
    ! ----------------------------------------------------------------------------------------
    
    if (adjust_htsnn) then
       
       call StieglitzSnow_calc_tpsnow(N_snow, htsnn, wesn, tpsn, fices)
       
       !**** Check that (ice10,tzero) conditions are conserved through
       !**** relayering process (or at least that (fices,tpsn) conditions don't 
       !**** go through the (1,0) point); excess goes to hcorr.
       
       ! for each layer, check snow conditions (partially/fully frozen, temp at/below zero) 
       !   before and after relayer; in select cases, adjust snow heat content and temp
       !
       ! NOTE: logicals before relayer were computed with    "buffer" (use_threshold_fac=.true. )
       !       reals    after  relayer were computed without "buffer" (use_threshold_fac=.false.)
       
       do i=1,N_snow                                                          
          
          kflag = .false.                                                   ! default: do nothing
          
          ! set klfag to .true. under certain conditions:
          
          if(     ice10(i) .and.      tzero0(i) .and.                   &   ! if     before relayer: fully     frozen and at    0 deg
               (fices(i) .ne. 1. .or.  tpsn(i) .ne. 0.) ) kflag=.true.      !    and after  relayer: partially frozen or  below 0 deg (or above 0 deg?)
          
          if(.not.ice10(i) .and.      tzero0(i) .and.                   &   ! if     before relayer: partially frozen and at    0 deg
               (fices(i) .eq. 1. .and. tpsn(i) .lt. 0.) ) kflag=.true.      !    and after  relayer: fully     frozen and below 0 deg
          
          if(     ice10(i) .and. .not.tzero0(i) .and.                   &   ! if     before relayer: fully frozen     and below 0 deg
               (fices(i) .ne. 1. .and. tpsn(i) .eq. 0.) ) kflag=.true.      !    and after  relayer: partially frozen and at    0 deg
          
          if (kflag) then                                                    
             
             ! make fully frozen and at 0 deg

             hnew    = -alhm*wesn(i)
             
             ! add "excess" heat content to hcorr
             
             if (update_hcorr) hcorr = hcorr+(htsnn(i)-hnew)/dts          

             htsnn(i)= hnew
             tpsn(i) = 0.
             fices(i)= 1.
          endif
       
       enddo
    
    end if  ! (adjust_htsnn)
    
    return
    
  end subroutine StieglitzSnow_relayer
  
  ! **********************************************************************
  
  subroutine StieglitzSnow_calc_tpsnow_scalar( h, w, t, f, ice1, tzero,  &
       use_threshold_fac )
    
    ! diagnose snow temperature and frozen fraction from snow mass and snow heat content
    !
    ! scalar version of StieglitzSnow_calc_tpsnow() with two differences:
    !   1) contains hardcoded multiplier 1.+eps in "if (hbw < -ALHM)" condition
    !                                and 1.-eps in "if (hbw > -ALHM)" condition
    !   2) additional logical outputs ice1 and tzero:
    !         ice1  = .true.  -->  frozen fraction (fice) equal to 1.
    !         tzero = .true.  -->  snow temperature at 0 deg C
    
    ! reichle,  6 Oct 2023:
    !   modified to have single subroutine that can maintain the above-mentioned differences
    !   (and thus 0-diff test results) and can provide a single interface with the science 
    !   calculations being in just one place

    implicit none
    
    !RR      real, parameter :: cpw    = 2065.22   !  @ 0 C [J/kg/K] -- ALREADY DEFINED ABOVE
    !      real, parameter :: lhv    = 2.4548E6 !  2.5008e6   !  @ 0 C [J/kg]
    !      real, parameter :: lhs    = 2.8368E6 !  2.8434e6 !  @ 0 C [J/kg]
    !rr   real, parameter :: lhv    = 2.5008e6  !  @ 0 C [J/kg]
    !rr   real, parameter :: lhs    = 2.8434e6  !  @ 0 C [J/kg]
    !      real, parameter :: lhf    = (lhs-lhv) !  @ 0 C [J/kg]
    
    real,    intent(in )   :: w, h          ! snow mass (SWE), snow heat content 
    real,    intent(out)   :: t, f          ! snow temperature, frozen ("ice") fraction
    
    logical, intent(out)   :: ice1, tzero   ! frozen fraction==1?, snow temp at 0 deg C?
    
    logical, intent(in)    :: use_threshold_fac
    
    ! ------------------------------------------------------------
    
    real,    parameter     :: tfac=1./StieglitzSnow_CPW
    real,    parameter     :: ffac=1./alhm
    
    real                   :: hbw
    
    real                   :: threshold1, threshold2
    
    ! ------------------------------------------------------------------------------
    
    if (use_threshold_fac) then
       
       ! replicates original get_tf0d()
       
       threshold1 = -1.00001*alhm     
       threshold2 = -0.99999*alhm

    else

       ! replicates original get_tf_nd() / StieglitzSnow_calc_tpsnow[_vector]()
       
       threshold1 = -alhm                      
       threshold2 = -alhm                      
              
    end if
    
    ! -------------------------------------------------------------------    

    hbw=0.
    
    if (w > 0.) hbw = h/w
    
    if     (hbw < threshold1) then    ! fully     frozen, temp below 0 deg

       t     = (hbw+alhm)*tfac
       f     = 1.
       ice1  = .true.
       tzero = .false.

    elseif (hbw > threshold2) then    ! partially frozen, temp at    0 deg

       t     = 0.
       f     = -hbw*ffac
       ice1  = .false.
       tzero = .true.

    else                              ! fully     frozen, temp at    0 deg
       t     = 0.
       f     = 1.
       ice1  = .true.
       tzero = .true.

    endif
    
    if (f <  0.) then                 ! (i.e., h>0 and f=-hbw/alhm via "partially frozen, temp at 0 deg")

       t = hbw*tfac                   ! t>0.  ??????
       f = 0.

    endif
    
    if (w == 0.) then                 ! no snow

       t = 0.
       f = 0.

    endif
    
    return
    
  end subroutine StieglitzSnow_calc_tpsnow_scalar
  
  ! **********************************************************************
  
  subroutine StieglitzSnow_calc_tpsnow_vector( N, h, w, t, f )
    
    ! renamed for clarity:   get_tf_nd() --> StieglitzSnow_calc_tpsnow()
    ! reichle, 12 Aug 2014
    
    !     n-dimensional version of get_tf
    !     
    !     avoid slow "where" statements
    !     
    !     can be called for any number of layers or catchments, for example
    
    !     1.) call StieglitzSnow_calc_tpsnow( ncatm, htsnn1(1:ncatm), wesn1(1:ncatm),
    !                         tpsn(1:ncatm),f(1:ncatm) )
    !     
    !     2.) call StieglitzSnow_calc_tpsnow(N_snow, h, w, t, f)
    
    !     reichle, 22 Aug 2002
    !     reichle, 29 Apr 2003 (updated parameter values)
    
    ! modified to call StieglitzSnow_calc_tpsnow_scalar() while maintaining 0-diff
    !   [avoiding same science equations in two different places]

    integer,               intent(in)  :: N
    
    real,    dimension(N), intent(in)  :: h, w
    real,    dimension(N), intent(out) :: t, f
    
    ! -----------------------------------
    
    integer            :: ii      
    
    logical            :: ice1, tzero   
    
    logical, parameter :: use_threshold_fac = .false.

    ! ----------------------------------
    
    do ii=1,N
       
       call StieglitzSnow_calc_tpsnow_scalar( h(ii), w(ii), t(ii), f(ii), ice1, tzero,  &
            use_threshold_fac )
       
    end do

  end subroutine StieglitzSnow_calc_tpsnow_vector
  
  ! ********************************************************************
  
  subroutine StieglitzSnow_calc_asnow_1( N_snow, NTILES, wesnn, asnow )
    
    ! Calculate diagnostic snow area from prognostic SWE
    !
    ! *_1(): input SWE for multiple snow layers at multiple tiles
    !
    ! reichle, Nov 3, 2004
    ! reichle,  2 Apr 2012 - revised for use without catch_types structures
    ! reichle, 12 Aug 2014 - moved to here from catchment.F90 
    !                          [formerly known as catch_calc_asnow()]
    !
    ! ----------------------------------------------------------------
    
    implicit none
    
    integer,                        intent(in)  :: N_snow, NTILES    
    real, dimension(N_snow,NTILES), intent(in)  :: wesnn    
    real, dimension(       NTILES), intent(out) :: asnow
    
    ! -----------------------------------------------------------
    
    asnow = max( min( sum(wesnn,1)/wemin, 1. ), 0. )
    
  end subroutine StieglitzSnow_calc_asnow_1
  
  ! *************************
  
  subroutine StieglitzSnow_calc_asnow_2( N_snow, wesnn, asnow )
    
    ! Calculate diagnostic snow area from prognostic SWE 
    !
    ! *_2(): input SWE for multiple snow layers at single tile
    
    implicit none
    
    integer,                 intent(in)  :: N_snow
    real, dimension(N_snow), intent(in)  :: wesnn
    real,                    intent(out) :: asnow
    
    ! -----------------------------------------------------------
    
    asnow = max( min( sum(wesnn)/wemin, 1. ), 0. )
    
  end subroutine StieglitzSnow_calc_asnow_2

  ! *************************
  
  subroutine StieglitzSnow_calc_asnow_3( totswe, asnow )
    
    ! Calculate diagnostic snow area from prognostic SWE
    !
    ! *_3(): input total SWE at single tile
    
    implicit none
    
    real,    intent(in)  :: totswe
    real,    intent(out) :: asnow
    
    ! -----------------------------------------------------------
    
    asnow = max( min( totswe/wemin, 1. ), 0. )
    
  end subroutine StieglitzSnow_calc_asnow_3

  ! ********************************************************************
  
  SUBROUTINE StieglitzSnow_trid(X,DD,D,RD,B,N)
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: N
      REAL*4, INTENT(IN), DIMENSION(N) :: DD, RD
      REAL*4, INTENT(INOUT), DIMENSION(N) :: D, B
      REAL*4, INTENT(OUT),DIMENSION(N) :: X

      integer I,J
      real*4  RSF
      RSF=0.
      do I=2,N
         J=N+1-I
         if(D(J+1).ne.0.) RSF=RD(J)/D(J+1)
         D(J)=D(J)-DD(J+1)*RSF
         B(J)=B(J)- B(J+1)*RSF
      enddo
      if(D(1).ne.0.) X(1)=B(1)/D(1)
      do J=2,N
         if(D(J).ne.0.) X(J)=(B(J)-DD(J)*X(J-1))/D(J)
      enddo
      RETURN

  END SUBROUTINE StieglitzSnow_trid

  !=======================================================================
  !      Version 5.0.2  by Teppei J. Yasuanari on 02/14/2011
  
  SUBROUTINE StieglitzSnow_snow_albedo(                                  &
       NCH, N_snow, N_constit_type, ITYP, VLAI, ZTH,                     &
       RHOFRESH,                                                         &   
       SNWALB_VISMAX, SNWALB_NIRMAX, SLOPE,                              & 
       WESN, HTSNN, SNDZ,                                                &
       AVISDR, ANIRDR, AVISDF, ANIRDF,                                   &
       ASNVDR, ASNNDR, ASNVDF, ASNNDF,                                   &
       RCONSTIT, UM, RTS, PARDIR, PARDIF                                 &
       )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NCH, N_snow, N_constit_type
    
    ! --------------------------------------
    ! INPUTS:
    
    ! NCH:      number of tiles considered
    ! N_snow:   number of snow layers
    ! N_constit: number of constituents
    ! ITYP:     vegetation type
    ! VLAI:     the leaf area index.
    ! ZTH:      The cosine of the solar zenith angle.
    ! RHOFRESH: density of fresh snow
    ! SNWALB_VISMAX:   max of visible snow albedo
    ! SNWALB_NIRMAX:   max of NIR snow albedo
    ! SLOPE:           slope the albedo decreases higher density
    !                  SLOPE > 0:  it gets recomputed
    !                  SLOPE < 0:  it is directly used in the computation; 
    !                              for example,to recover the old formulation,
    !                              set SLOPE=-0.0006, SNWALB_VISMAX=0.7, SNWALB_NIRMAX=0.5 
    ! WESN:     snow water equivalent in each layer
    ! HTSNN:    heat content of each layer
    ! SNDZ:     depth of each layer
    ! UM:       wind speed
    ! RTS:      surface temperature
    ! PARDIR:   photosynthetically active radiation, direct
    ! PARDIF:   photosynthetically active radiation, diffuse
    ! RCONSTIT: array of constituent masses
    ! ABVIS:    const array (not used when there is no constituents)    
    ! ABNIR:    const array (not used when there is no constituents) 
    ! AVISDR:   visible, direct albedo (snow-free).
    ! ANIRDR:   near infra-red, direct albedo (snow-free).
    ! AVISDF:   visible, diffuse albedo (snow-free).
    ! ANIRDF:   near infra-red, diffuse albedo (snow-free).
    
    INTEGER, INTENT(IN), DIMENSION(NCH)           :: ITYP
    REAL,    INTENT(IN)                           :: RHOFRESH 
    REAL,    INTENT(IN)                           :: SNWALB_VISMAX, SNWALB_NIRMAX, SLOPE 
    REAL,    INTENT(IN), DIMENSION(NCH)           :: AVISDR, ANIRDR, AVISDF 
    REAL,    INTENT(IN), DIMENSION(NCH)           :: ANIRDF, VLAI, ZTH
    REAL,    INTENT(IN), DIMENSION(N_Snow,NCH)    :: WESN, HTSNN, SNDZ
    
    REAL,    INTENT(IN), DIMENSION(NCH),                  OPTIONAL :: UM, RTS, PARDIR, PARDIF
    
    REAL,    INTENT(IN), DIMENSION(NCH,N_snow,N_Constit), OPTIONAL :: RCONSTIT
    
    ! --------------------------------------
    ! OUTPUTS:
    
    ! ASNVDR: snow albedo, visible direct
    ! ASNNDR: snow albedo, near-infrared direct
    ! ASNVDF: snow albedo, visible diffuse
    ! ASNNDF: snow albedo, near-infrared diffuse
    
    REAL, INTENT(OUT), DIMENSION(NCH) :: ASNVDR, ASNNDR, ASNVDF, ASNNDF
    
    ! --------------------------------------
    ! Other variables as needed.  Includes:
    ! SSA: snow specific surface area
    ! EFFG: effective ice thickness (m)
    
    INTEGER            :: I
    INTEGER, PARAMETER :: NTYPS_SIB=9
    
    REAL ::                                             &
         FAC, FVEG, TOTDEP, SWE, DENS_EXC, AREASC,      &
         DENSITY, ASNVDR_VEG, ASNNDR_VEG, ASNVDF_VEG,   &
         ASNNDF_VEG, GK_B
    
    REAL, DIMENSION(NTYPS_SIB) :: SNWMID
    
    DATA SNWMID /50.,30.,45.,20.,30.,20.,2.,2.,2./
    
    ! *********************************************************************

    if(SLOPE < 0.0) then
       GK_B = SLOPE
    else
       GK_B = (0.85808-0.6)/(RHOFRESH-StieglitzSnow_RHOMA)
    endif
    
    DO I=1,NCH
       
       SWE=SUM(WESN(:,I))
       
       TOTDEP=SNDZ(1,I)
       call StieglitzSnow_calc_asnow( SWE, AREASC )
       !DENSITY=(SWE/(AREASC+1.e-20)) / (TOTDEP+1.e-20)
       !*** only use top layer density to dentermine albedo 
       DENSITY=(WESN(1,I)/(AREASC+1.e-20)) / (TOTDEP+1.e-20)
       DENS_EXC=MAX(0., DENSITY-RHOFRESH)
       
       !*********************************************************************
       
       !     Using snow tracer albedo scheme only when N_constit_type > 0
       
       if(N_constit_type > 0) then
          call ALB_WITH_IMPURITY (N_snow, ZTH(I), SNWALB_VISMAX, SNWALB_NIRMAX,        & 
               WESN(:,I),HTSNN(:,I),SNDZ(:,I), UM(I), RTS(I), PARDIR(I), PARDIF(I),    &
               ASNVDR(I), ASNNDR(I), ASNVDF(I), ASNNDF(I),RCONSTIT(I,:,:))
       else
          
          !       Use these when you use the original snow albedo model (comment out for alb)
          !ASNVDR(I) = MAX(SNWALB_VISMIN, SNWALB_VISMAX - DENS_EXC*.0006)
          !ASNNDR(I) = MAX(SNWALB_NIRMIN, SNWALB_NIRMAX - DENS_EXC*.0006)
          !ASNVDF(I) = MAX(SNWALB_VISMIN, SNWALB_VISMAX - DENS_EXC*.0006)
          !ASNNDF(I) = MAX(SNWALB_NIRMIN, SNWALB_NIRMAX - DENS_EXC*.0006)
          ASNVDR(I) = MAX(SNWALB_VISMIN, SNWALB_VISMAX + GK_B*DENS_EXC)
          ASNNDR(I) = MAX(SNWALB_NIRMIN, SNWALB_NIRMAX + GK_B*DENS_EXC)
          ASNVDF(I) = MAX(SNWALB_VISMIN, SNWALB_VISMAX + GK_B*DENS_EXC)
          ASNNDF(I) = MAX(SNWALB_NIRMIN, SNWALB_NIRMAX + GK_B*DENS_EXC)
          
       endif
       
       ! ACCOUNT FOR VEGETATION MASKING, FOR EACH COMPONENT
       
       ! A) FIRST DO MASKING IN VEGETATED FRACTION:
       FAC = SWE / (SWE + SNWMID(ITYP(I)))
       ASNVDR_VEG=AVISDR(I) + (ASNVDR(I)-AVISDR(I))*FAC
       ASNNDR_VEG=ANIRDR(I) + (ASNNDR(I)-ANIRDR(I))*FAC
       ASNVDF_VEG=AVISDF(I) + (ASNVDF(I)-AVISDF(I))*FAC
       ASNNDF_VEG=ANIRDF(I) + (ASNNDF(I)-ANIRDF(I))*FAC
       
       ! B) NOW ACCOUNT FOR SUBGRID VEGETATION FRACTION
       FVEG=AMIN1( 1., VLAI(I)/2. )
       ASNVDR(I)=ASNVDR(I)*(1.-FVEG)+ASNVDR_VEG*FVEG
       ASNNDR(I)=ASNNDR(I)*(1.-FVEG)+ASNNDR_VEG*FVEG
       ASNVDF(I)=ASNVDF(I)*(1.-FVEG)+ASNVDF_VEG*FVEG
       ASNNDF(I)=ASNNDF(I)*(1.-FVEG)+ASNNDF_VEG*FVEG
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE StieglitzSnow_Snow_Albedo
  
  ! ****************************************************************************
  
  SUBROUTINE ALB_WITH_IMPURITY (N_snow, ZTH,               &
       SNWALB_VISMAX, SNWALB_NIRMAX,                       &   
       WESN, HTSNN, SNDZ, UM, RTS, PARDIR, PARDIF,         &
       ASNVDR, ASNNDR, ASNVDF, ASNNDF,                     &
       RCONSTIT                                            &
       )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N_snow
    
    
    ! --------------------------------------
    ! INPUTS:
    
    ! ZTH:      The cosine of the solar zenith angle.
    ! WESN:     snow water equivalent in each layer
    ! HTSNN:    heat content of each layer
    ! SNDZ:     depth of each layer
    ! UM:       wind speed
    ! RTS:      surface temperature
    ! PARDIR:   photosynthetically active radiation, direct
    ! PARDIF:   photosynthetically active radiation, diffuse
    ! RCONSTIT: array of constituent masses
    ! AVISDR:   visible, direct albedo (snow-free).
    ! ANIRDR:   near infra-red, direct albedo (snow-free).
    ! AVISDF:   visible, diffuse albedo (snow-free).
    ! ANIRDF:   near infra-red, diffuse albedo (snow-free).
    
    REAL, INTENT(IN)                                        :: SNWALB_VISMAX, SNWALB_NIRMAX
    REAL, INTENT(IN)                                        :: ZTH, UM, RTS, PARDIR, PARDIF
    REAL, INTENT(IN), DIMENSION(N_snow,N_Constit), OPTIONAL :: RCONSTIT
    REAL, INTENT(IN), DIMENSION(N_Snow)                     :: WESN, HTSNN, SNDZ
    
    ! --------------------------------------
    ! OUTPUTS:
    
    ! ASNVDR: snow albedo, visible direct
    ! ASNNDR: snow albedo, near-infrared direct
    ! ASNVDF: snow albedo, visible diffuse
    ! ASNNDF: snow albedo, near-infrared diffuse
    
    REAL, INTENT(OUT) :: ASNVDR, ASNNDR, ASNVDF, ASNNDF
    
    ! --------------------------------------
    ! Other variables as needed.  Includes:
    ! SSA: snow specific surface area
    ! RHO_FS: fresh snow density
    ! EFFG: effective ice thickness (m)
    
    INTEGER :: M,K
    INTEGER, PARAMETER :: NTYPS_SIB=9
    
    REAL :: rho_fs, DEGSZA, SZASIN, COS50,                 &
         WSS, TS, FAC, TOTDEP, SWE, DENS_EXC, AREASC,      &
         DENSITY, SUM1, SUM2
    
    REAL :: SZTH
    
    REAL, DIMENSION(N_snow)           :: ABSCOV, ABSCON, EFFG, SSA, DENEL
    
    REAL, DIMENSION(N_snow,N_Constit) :: CONCENT 
    
    REAL, DIMENSION(N_snow)           :: CWESN, CHTSNN, CSNDZ
    REAL, DIMENSION(N_snow)           :: TPSN, FICES
    REAL, DIMENSION(N_snow)           :: CTPSN, CFICES
    
    ! *********************************************************************

    SZTH=ZTH
    DEGSZA=ACOS(SZTH)*180./PIE
    SZASIN=SQRT(1.-(SZTH**2.0))
    COS50=COS(2.*PIE*50./360.)
    
    !    When it is cloud-covered, SZTA is set to 50 degree.  THE VALUE
    !    USED HERE (0.1) CAN BE TUNED!
    IF(pardir/(pardir+pardif+1.e-20) < 0.1) SZTH=COS50
    
    ! CTPSN: LAYER TEMPERATURE [degree C]
    DO M=1,N_Snow
       CWESN(M)=WESN(M)
       CHTSNN(M)=HTSNN(M)
    END DO
    
    CALL StieglitzSnow_calc_tpsnow(N_snow,CHTSNN,CWESN,CTPSN,CFICES)
    
    DO M=1,N_Snow
       TPSN(M)=CTPSN(M)
       FICES(M)=CFICES(M)
    END DO
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
    ! SNOW ALBEDOES
    
    SWE=SUM(WESN(:))
    TOTDEP=SUM(SNDZ(:))
    call StieglitzSnow_calc_asnow( SWE, AREASC )
    DENSITY=(SWE/(AREASC+1.e-20)) / (TOTDEP+1.e-20)
    
    WSS=UM
    TS=RTS
    
    CALL WFSDEN(WSS,TS,RHO_FS)
    DENS_EXC=MAX(0., DENSITY-RHO_FS)
    
    !*********************************************************************
        
    IF(SWE > 0.01) THEN
       
       !========== SNOW CASE ==========
           
       DO M=1,N_SNOW
          
          ! Dry snow density in each snow layer [kg m-3] (N_Snow layers)
          
          DENEL(M)=(FICES(M)*WESN(M)/(AREASC+1.e-20))/(SNDZ(M) + 1.e-20)
          DENEL(M)=MIN(StieglitzSnow_RHOMA,DENEL(M))
          
          ! Mass concentrations of dust,BC, and OC [kg kg-1] (3 layers)
          
          IF(PRESENT(RCONSTIT)) THEN
             DO K=1,N_Constit
                CONCENT(M,K)=RCONSTIT(M,K)/(WESN(M) + 1.e-20)
             ENDDO
          ELSE
             DO K=1,N_Constit
                CONCENT(M,K)=0.
             ENDDO
          ENDIF
          
          !     Snow specific surface area (SSA) [m2 kg-1] Equation based on Yamazaki et al. [1991,1993]
          !     The equation is also based on the observed data by Narita [1971], Teionkagaku      
          
          SSA(M)=10.0**(-15.32*((DENEL(M)/1000.)**3.0)              &
               +16.65*((DENEL(M)/1000.)**2.0)                       &
               -7.3*(DENEL(M)/1000.)+2.23)
          
          !      Another Equation for SSA based on the equation (1) in Domine et al., JGR, vol.112, 2007
          !       SSA(M,I)=( -308.2D0*LOG(DENEL(M,I)/1000.D0)-206.D0 )/10.D0
          !     When each snow layer has excess water (snow gets wet)
          !       SSA is decreased to 60% values of the original value.
          
          IF(TPSN(M).GE.-0.001) SSA(M)=0.6*SSA(M)
          
          !     effective ice thickness comparable to effective snow grain radius
          !     if EFFG is multiplied by 3 
          
          EFFG(M)=2.D0/(DENICE*SSA(M) + 1.e-20)
          
          !   Calculation of Absorption coefficients [m-1] Equations for mass con. [kg kg-1]
          
          SUM2=0.
          do K=1,N_Constit
             SUM2=SUM2+CONCENT(M,K)
          ENDDO
          
          ! VIS range
          SUM1=0.
          do K=1,N_Constit
             SUM1=SUM1+CONCENT(M,K)*1.E03*ABVIS(K)
          ENDDO
          ABSCOV(M)=DENICE*SUM1+(1.-SUM2)*AICEV
              
          ! NIR range
          SUM1=0.
          do K=1,N_Constit
             SUM1=SUM1+CONCENT(M,K)*1.E03*ABNIR(K)
          ENDDO
          ABSCON(M)=DENICE*SUM1+(1.-SUM2)*AICEN 
          
       ENDDO
       
       ! VIS & NIR CASES BY NEW SNOW ALBEDO MODEL
       
       ASNVDR=ALB(N_Snow,RIV,ABSCOV(:),EFFG(:),DENEL(:),SNDZ(:),DENICE)
       ASNNDR=ALB(N_Snow,RIN,ABSCON(:),EFFG(:),DENEL(:),SNDZ(:),DENICE)
       
       ! Diffuse components are calculated with the assumed SZTH of 50 degrees
       ASNVDF = ASNVDR 
       ASNNDF = ASNNDR
       
       !+++++ Taking the effect of solar zenith angle (SZA) into account ++++++            
       ! The second terms for VIS and NIR in equations (6) & (7)
       ! [Marks and Dozier, Water Resources Research, Vol. 28, 3043-3054]
       
       ! For VIS [Note: unit for EFFG is um]
       ASNVDR = ASNVDR                                                    &
            -( (SQRT(1.5*(EFFG(1)*1.0E06))*1.375E-03)*(1.-COS50) )        & 
            +( (SQRT(1.5*(EFFG(1)*1.0E06))*1.375E-03)*(1.-SZTH) )   
       
       ! For NIR [Note: unit for EFFG is um]
       ASNNDR = ASNNDR                                                    &
            -( (( SQRT(1.5*(EFFG(1)*1.0E06))*2.0E-03)+0.1)*(1.-COS50))    & 
            +( (( SQRT(1.5*(EFFG(1)*1.0E06))*2.0E-03)+0.1)*(1.-SZTH) )
       
    ELSE
       
       !========== NO SNOW CASE ==========            
       !           Use these when you use the original snow albedo model 
       
       ASNVDR = MAX(SNWALB_VISMIN, SNWALB_VISMAX - DENS_EXC*.0006)
       ASNNDR = MAX(SNWALB_NIRMIN, SNWALB_NIRMAX - DENS_EXC*.0006)
       ASNVDF = MAX(SNWALB_VISMIN, SNWALB_VISMAX - DENS_EXC*.0006)
       ASNNDF = MAX(SNWALB_NIRMIN, SNWALB_NIRMAX - DENS_EXC*.0006)
       
    ENDIF
    
    RETURN
    
  END SUBROUTINE ALB_WITH_IMPURITY
      
  ! **********************************************************************

  !**** ------------------------------------------------------------------
  !**** //////////////// Added by Teppei J. Yasunari ///////////////START1
  !**** ------------------------------------------------------------------
  
  !=======================================================================
  !=                         To determine RHOFS                          =
  !=======================================================================
  !      Version 5.1.0  by Teppei J. Yasuanari on 09/15/2014
  
  SUBROUTINE WFSDEN(UM,RTS,RHO_FS)
    
    
    REAL, INTENT(IN)  :: UM,RTS
    REAL, INTENT(OUT) :: RHO_FS
    REAL              :: ARHOFS
    INTEGER           :: IFG
    
    !--------------------------------------------------------------!
    !          YOU CAN CHOOSE ONE OF 10 TYPES OF "RHOFS"           !
    !--------------------------------------------------------------!
    
    !     ESAT : Saturated water vapor pressure [hPa]
    !     E : Water vapor pressure [hPa]
    !     TW : Wet bulb temperature [degree C]
    !
    !      V1=UM
    !      T2=RTS-273.15
    !      SR=TSNOW*3600.
    !
    !     Saturated water vapor pressure

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !     Fresh snow density by Melloh et al. (2002)
    !     Hydrological Processes, 16, 3571-3584
    !
    ! RHOFS TYPE 1
    !
    !      IF(T2.GE.0.) THEN
    !       TW=T2-((1013.25/(0.667*RPRES))*(ESAT-E))
    !      ELSE IF (T2.LT.0.) THEN
    !       TW=T2-((1013.25/(0.589*RPRES))*(ESAT-E))
    !      END IF
    !
    !    TWW=TW
    !      ARHOFS=(0.05+(0.0017*(((TW+273.15)-258.16)**1.5)))*1000.
    !      ARHOFS=MAX(ARHOFS,90.)
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !     Multiple regression equation for new snow density
    !     based on Kajikawa et al. (2004), Sepyo, 66, p.561-565.
    !
    ! RHOFS TYPE 2
    !     Spatial dendrite (R**2=0.948; 10% confidence limit)
    !      ARHOFS=13.3+(53.9*SR)+(6.54*V1)
    !
    ! RHOFS TYPE 3
    !     Stellar crystal (R**2=0.817; 5% confidence limit)
    !      ARHOFS=23.4+(37.5*SR)+(7.32*V1)+(0.579*T2)
    !
    ! RHOFS TYPE 4
    !     Rimed stellar crystal (R**2=0.442; 5% confidence limit)
    !      ARHOFS=41.2+(8.26*SR)+(5.16*V1)+(0.422*T2)
    !
    ! RHOFS TYPE 5
    !     Rimed spatial dendrites (R**2=0.369; 5% confidence limit)
    !      ARHOFS=67.5+(23.4*SR)-(1.29*V1)+(3.65*T2)
    !
    ! RHOFS TYPE 6
    !     Rimed spatial dendrites 2 (R**2=0.321; 5% confidence limit)
    !      ARHOFS=43.1*EXP(0.106*V1)
    !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 2001 J. Hydrometeorology (Boone and Etchevers)
    ! Fresh snow density based on wind speed and air temeperature
    !
    ! RHOFS TYPE 7
    !        ARHOFS=109.+(6.*(RTS-273.16))+(26.*(V1**0.5))
    !
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Equation 18  by Yamazaki et al. [1998]
    ! Seppyo, 60, 131-141.
    
    ! RHOFS TYPE 8
    ARHOFS=67.+(13.*UM)
    IFG=2
    IF(RTS-TF >= 0.) ARHOFS=200.
    
    !+++ USE Constant value of fresh snow density [kg m-3]+++
    ! RHOFS TYPE 9
    !      ARHOFS=110.                 ! Used in Yasunari et al. (2010)
    
    ! RHOFS TYPE 10
    !      ARHOFS=150.                 ! Used in original snowpack model
    
    
    ! RHOFS does not exceed the maximum of 200. [kg m-3].
    RHO_FS=MIN(ARHOFS,200.)
    
    !      write(*,*) 'FRESH SNOW DENSITY',RHOFS
    
    ! (comment out for den)
    rho_fs=150.
    
    RETURN

  END SUBROUTINE WFSDEN
  
  ! **************************************************************************************
  
  
  !=======================================================================
  !=                   New snow albedo functions                         =
  !=======================================================================
  !      Version 5.1.0  by Teppei J. Yasuanari on 09/15/2014
  
  FUNCTION ALB(L,RI0,ABSCO0,EFFG0,DENEL0,SNDZ0,DENICE)
    
    IMPLICIT NONE
    
    INTEGER :: J,K,L2
    INTEGER, INTENT(IN) :: L
    REAL :: ALB
    DOUBLE PRECISION :: ALBD,DZ,CDZ
    REAL, INTENT(IN) :: RI0,DENICE
    REAL, INTENT(IN), DIMENSION(L) :: EFFG0,ABSCO0,SNDZ0,DENEL0
    DOUBLE PRECISION :: RI,DDENICE
    DOUBLE PRECISION, DIMENSION(L) :: EFFG,ABSCO,SNDZ,DENEL
    DOUBLE PRECISION, DIMENSION(L) :: REF,TRANS,A,B,AMU,TJALPHA,TJBETA,DI
    
    ! From single to double precision
    
    RI=DBLE(RI0)
    DDENICE=DBLE(DENICE)
    EFFG(:)=DBLE(EFFG0(:))
    ABSCO(:)=DBLE(ABSCO0(:))
    SNDZ(:)=DBLE(SNDZ0(:)) 
    DENEL(:)=DBLE(DENEL0(:))
    
    ! Only the top 3m of snow contributes to albedo 
    
    CDZ=0.d0
    
    DO J=1,L
       CDZ=CDZ+SNDZ(J)
       IF(CDZ > 3.d0) EXIT
       L2=J
    END DO
    
    DO J=1,L2
       
       !      Reflectance in one ice layer
       
       REF(J)= RI +                                                         &
            ( ((1.d0-RI)**2.d0)*RI*EXP(-2.d0*ABSCO(J)*EFFG(J))/             &
            MAX(1.d-20,(1.d0-( (RI*EXP(-ABSCO(J)*EFFG(J)))**2.d0 )) ) )
       
       !      Transparency in one ice layer
       
       TRANS(J)=                                                            &
            ( ((1.d0-RI)**2.d0)*EXP(-ABSCO(J)*EFFG(J))/                     &
            MAX(1.d-20,(1.d0-( (RI*EXP(-ABSCO(J)*EFFG(J)))**2.d0 )) ) )
       
       
       A(J)=( (1.d0-TRANS(J))/(EFFG(J)+1.d-20) )*( DENEL(J)/DDENICE )
       B(J)=( REF(J)/(EFFG(J)+1.d-20) )*( DENEL(J)/DDENICE )
       
       AMU(J)=SQRT( MAX(0.d0,(A(J)**2.d0) - (B(J)**2.d0)) )
       TJALPHA(J)=( A(J)-AMU(J) )/(B(J)+1.d-20)
       TJBETA(J)=( A(J)+AMU(J) )/(B(J)+1.d-20)
       
    END DO
    
    DI(L2)=0.d0
    
    DO K=L2-1,1,-1
       
       !----- Revision for multi-layers (Teppei, Sep. 15, 2014) -----
       DZ=SUM(SNDZ(1:K))
       
       DI(K)=(TJALPHA(K)-TJBETA(K+1))*DI(K+1)*EXP(2.d0*DZ*AMU(K+1))+TJALPHA(K)-TJALPHA(K+1)
       DI(K)=DI(K)/( (TJBETA(K+1)-TJBETA(K))*DI(K+1)*EXP(2.d0*DZ*AMU(K+1))+TJALPHA(K+1)-TJBETA(K) )
       DI(K)=DI(K)*EXP(-2.d0*DZ*AMU(K))
       
    END DO
    
    !     Surface Snow Albedo
    
    ALBD= RI +                                                               &
         ( ((1.d0-RI)**2.d0)*(TJBETA(1)*DI(1)+TJALPHA(1))/                   &
         MAX(1.d-20,(1.d0+DI(1)-(RI*(TJBETA(1)*DI(1)+TJALPHA(1))) ) ) )
    
    ALB=REAL(ALBD)
    
    RETURN
  END FUNCTION ALB
  
  ! **********************************************************************
  
  subroutine StieglitzSnow_echo_constants(logunit)
    
    ! reichle, 12 Aug 2014
    
    implicit none
    
    integer, intent(in) :: logunit
    
    write (logunit,*)
    write (logunit,*) '-----------------------------------------------------------'
    write (logunit,*)
    write (logunit,*) 'StieglitzSnow_echo_constants():'
    write (logunit,*)
    write (logunit,*) 'PIE                  = ', PIE      
    write (logunit,*) 'ALHE                 = ', ALHE     
    write (logunit,*) 'ALHM                 = ', ALHM     
    write (logunit,*) 'TF                   = ', TF    
    write (logunit,*) 'RHOW                 = ', RHOW     
    write (logunit,*)
    write (logunit,*) 'StieglitzSnow_MINSWE = ', StieglitzSnow_MINSWE
    write (logunit,*) 'WEMIN                = ', WEMIN
    write (logunit,*) 'StieglitzSnow_CPW    = ', StieglitzSnow_CPW
    write (logunit,*) 'StieglitzSnow_RHOMA  = ', StieglitzSnow_RHOMA    
    write (logunit,*) 
    write (logunit,*) 'SNWALB_VISMIN        = ', SNWALB_VISMIN
    write (logunit,*) 'SNWALB_NIRMIN        = ', SNWALB_NIRMIN
    write (logunit,*) 
    write (logunit,*) 'end StieglitzSnow_echo_constants()'
    write (logunit,*)
    write (logunit,*) '-----------------------------------------------------------'
    write (logunit,*)
    
  end subroutine StieglitzSnow_echo_constants
  
end module StieglitzSnow

! ============================ EOF =========================================================
