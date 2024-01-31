      MODULE bl_seaspray_mod
 
!-------------------------------------------------------------------
!     From: 
!          Jian-Wen Bao [jian-wen.bao@noaa.gov] & Chris Fairall [chris.fairall@noaa.gov] (NOAA)
!     Reference:
!          Bao et al., Mon Wea Rev, 2011. 10.1175/MWR-D-11-00007.1
!
!     Code was provided by Jian-Wen Bao on 02/17/2015
!-------------------------------------------------------------------
!  Input variables
!
!  UU    ----------- wind speed at the lowest atm. level [m/s]
!  ZU    ----------- height of the lowest atm. model level [m]
!  Ts    ----------- sea surface temperature [C]
!  Ta    ----------- air temperature at the lowest model level [C]
!  MR    ----------- Mixing ratio at the lowest model level
!  Patm  ----------- Pressure at the lowest model level [mb]
!  hwave ----------- significant wave height [m/s]
!  cwave ----------- Phase speed of breaking wave [m]
!  p     ----------- Wave energy dissipation [(m/s)^3]
!  usr   ----------- Friction velocity [m/s]
!  hss   ----------- Sensible heat flux [W/m^2]
!  hll   ----------- Latent heat flux [W/m^2]
!
!  Output variable
!
!  usr_new --------- New friction velocity [m/s]
!  fxh_seaspray ---- New sensible heat flux [W/m^2]
!  fxe_seaspray ---- New latent heat flux [W/m^2]
!
!  Implemented by Laura Bianco
!--------------------------------------------------------
     

!--------------------------------------------------------
! Example:
!
!     real UU,ZU,Ts,Ta,MR,Patm,hwave,cwave,p,usr,hss_GFDL,hll_GFDL
!     real usr_new,fxh_seaspray,fxe_seaspray
!     real s,hss,hll,massf,hs_tot,hl_tot
!     real S_bar1
!     real z_r,omega,alpha,vfm
!     real sourcestrength_tune, feedback_tune
!     
!     ! Introducing the variables
!     
!     sourcestrength_tune = 0.4 ! Scales droplet mass flux, 0.5 = midway between Anreas/Fairall and de Leeuw
!     feedback_tune = 0.2       ! Doesn't change droplet contribution to enthalphy flux, but changes balance
!                               ! between droplet Hs and Hl
!
!     UU = 50.                  ! Wind speed
!     ZU = 25.                  ! Reference height for bulk data
!     Ts = 29.                  ! SST
!     Ta = 27.                  ! Air T
!     Patm = 980.               ! Atmospheric pressure
!     s = 0.90                  ! Saturatoin ratio (Relative humidity % /100)
!     hss = -999.0              ! Sensible heat flux
!     hll = -999.0              ! Latent heat flux

!     ! Start of spray-mediated fluxes
!     call execute_spray_param (sourcestrength_tune,feedback_tune,UU,ZU,
!    &             Ts,Ta,s,Patm,hss,hll,hwave,cwave,p,usr,massf,hs_tot,
!    &             hl_tot,usr_new,S_bar1,z_r,omega,alpha,vfm)
!
!     print *, 'hss (W/m^2) = ', hss, 'hss_spray (W/m^2) = ', hs_tot
!     print *, 'hll (W/m^2) = ', hll, 'hll_spray (W/m^2) = ', hl_tot
!     print *, 'usr (m/s)   = ', usr, 'usr_spray (m/s)   = ', usr_new
!
!     end
!--------------------------------------------------------

      implicit none

      private
      
#ifdef GEOS
      public :: online_spray
#else
      public :: execute_spray_param
#endif

      contains

      ! End of spray-mediated fluxes.



#ifdef GEOS
!********************************************************
! ESRL/PSD spray parameterization

      subroutine online_spray(sourcestrength_tune,feedback_tune, &
                         u,zu,ts,ta,s,Patm,hss,hll,hwave,cwave,p,usr,   &
                         massf,hs_tot,hl_tot,usr_new,S_bar1,z_r,        &
                         omega,alpha,vfm)

        implicit none

        real, intent(in) :: sourcestrength_tune,feedback_tune,u,zu,ts,ta, &
                s,Patm,hwave,cwave,p
        real, intent(inout) :: hss,hll   
        real, intent(out)::usr,massf,hs_tot,hl_tot,usr_new,S_bar1, &
                z_r,omega,alpha
        ! Output from spray_param subroutine.
        real, intent(out) :: vfm

        ! Set constants and calculations
        real :: Q, Le, gamma1, beta, wetdep
        real :: rhoa                       ! Air density
        !real, parameter :: tdk = 273.15    ! Celsius to Kelvin
        real, parameter :: tdk = 273.6     ! Celsius to Kelvin
        !real, parameter :: Rdgas = 287.05  ! Gas constant, dry air
        real, parameter :: Rdgas = 287.1   ! Gas constant, dry air
        real, parameter :: grav = 9.83     ! Gravitational acceleration
        real, parameter :: kon = 0.40      ! von Karman constant
        real, parameter :: rhow = 1e3      ! Water density
        real, parameter :: cpa = 1004.67   ! Specific heat air
        real, parameter :: cpw = 4.2e3     ! Specific heat water

        ! Wave properties
        real :: h

        real, parameter :: nit = 3         ! Number of iterations of spray param


        ! Gunn and Kinzer variables to compute iterated ustar.
        real, parameter :: rho_0 = 1.20    ! Reference density
        real, parameter :: Beta1 = 4.7
        real :: z
        real :: sigma
        real :: hsource
        real :: k2,S_r
        real :: cdzn,z0,U10

        integer :: i                       ! Iteration of droplet buoyancy

        ! Calculated constants
        Q = qsat(ta,Patm)*s/1e3                       ! Saturation ratio (RH%/100)
        Le = (2.501-0.00237*ts)*1e6                   ! Latent heat vaporation
        rhoa = Patm*100/(Rdgas*(ta+tdk)*(1+0.61*Q))   ! Air density
        sigma = (rhow-rhoa)/rhoa
        gamma1 = 241*17.501/(ta+241)**2
        beta = 1/(1+1e-3*qsat(ta,Patm)*Le/cpa*gamma1)
        wetdep = (1-s)*(1-beta)/gamma1


        call find_ust(u,kon,zu,rhoa,rhow,hwave,cwave,p,usr,h,z0)

        z_r = z0
        z = zu

        ! Iterate spray_param until usr is stabilized (depends on massflux).

        call spray_param(sourcestrength_tune,feedback_tune,u,zu,hss, &
                         hll,ts,ta,s,Patm,usr,hwave,cwave,p,h,       &
                         grav,kon,Q,tdk,Rdgas,Le,cpa,rhoa,cpw,rhow,  &
                         gamma1,beta,massf,U10,hs_tot,hl_tot,vfm)

        do i = 1,nit

           ! Droplet buoyancy

           ! Fall speed calculation

           k2 = 2.2*1e3*(rho_0/rhoa)**0.5         ! [cm^(1/2) s^(-1)]
           k2 = k2*1e-2                           ! [(m/s^2)^0.5] 

           omega = vfm/(kon*usr)                  ! From Lykossov (2001)

           !hsource = 0.05                         ! Set a height scale for the source (m)
           hsource = 0.01*h                       ! Set a height scale for the source (m)

           sigma = (rhow-rhoa)/rhoa
           S_r = massf/sigma/vfm/rhoa*(h/hsource)**omega ! Compute surface source from parameterization at wave height
           alpha = (Beta1*grav*kon**2*z_r*sigma*S_r)/(usr**2)

           if (omega /= 1) then

              S_bar1 = (omega**-1)*log(1+((alpha*omega**2)/ &
                       (1-omega))*((z/z_r)**(1-omega)-1))

           else

             S_bar1 = log(1+alpha*log(z/z_r))        

           endif

           usr_new = u*kon/(log(z/z0)+S_bar1)     ! u*
           cdzn = (usr_new/u)**2                  ! Drag coefficient

!!         hss = -999.0                           ! Sensible heat flux
!!         hll = -999.0                           ! Latent heat flux
          
           call spray_param(sourcestrength_tune,feedback_tune,u,zu,hss, &
                            hll,ts,ta,s,Patm,usr_new,hwave,cwave,p,h,   &
                            grav,kon,Q,tdk,Rdgas,Le,cpa,rhoa,cpw,rhow,  &
                            gamma1,beta,massf,U10,hs_tot,hl_tot,vfm)
        
        end do

        return

      end subroutine online_spray
#endif


!********************************************************
! ESRL/PSD spray parameterization

      subroutine execute_spray_param(sourcestrength_tune,feedback_tune, &
                         u,zu,ts,ta,s,Patm,hss,hll,hwave,cwave,p,usr,   &
                         massf,hs_tot,hl_tot,usr_new,S_bar1,z_r,        &
                         omega,alpha,vfm)

        implicit none

        real, intent(in) :: sourcestrength_tune,feedback_tune,u,zu,ts,ta, &
                s,Patm
        real, intent(inout) :: hss,hll   
        real, intent(out)::usr,massf,hs_tot,hl_tot,usr_new,S_bar1, &
                z_r,omega,alpha
        ! Output from spray_param subroutine.
        real, intent(out) :: vfm

        ! Set constants and calculations
        real :: Q, Le, gamma1, beta, wetdep
        real :: rhoa                       ! Air density
        !real, parameter :: tdk = 273.15    ! Celsius to Kelvin
        real, parameter :: tdk = 273.6     ! Celsius to Kelvin
        !real, parameter :: Rdgas = 287.05  ! Gas constant, dry air
        real, parameter :: Rdgas = 287.1   ! Gas constant, dry air
        real, parameter :: grav = 9.83     ! Gravitational acceleration
        real, parameter :: kon = 0.40      ! von Karman constant
        real, parameter :: rhow = 1e3      ! Water density
        real, parameter :: cpa = 1004.67   ! Specific heat air
        real, parameter :: cpw = 4.2e3     ! Specific heat water

        ! Wave properties
        real :: hwave,cwave,p,h

        real, parameter :: nit = 3         ! Number of iterations of spray param


        ! Gunn and Kinzer variables to compute iterated ustar.
        real, parameter :: rho_0 = 1.20    ! Reference density
        real, parameter :: Beta1 = 4.7
        real :: z
        real :: sigma
        real :: hsource
        real :: k2,S_r
        real :: cdzn,z0,U10

        integer :: i                       ! Iteration of droplet buoyancy

        ! Calculated constants
        Q = qsat(ta,Patm)*s/1e3                       ! Saturation ratio (RH%/100)
        Le = (2.501-0.00237*ts)*1e6                   ! Latent heat vaporation
        rhoa = Patm*100/(Rdgas*(ta+tdk)*(1+0.61*Q))   ! Air density
        sigma = (rhow-rhoa)/rhoa
        gamma1 = 241*17.501/(ta+241)**2
        beta = 1/(1+1e-3*qsat(ta,Patm)*Le/cpa*gamma1)
        wetdep = (1-s)*(1-beta)/gamma1

        call find_ust(u,kon,zu,rhoa,rhow,hwave,cwave,p,usr,h,z0)

        z_r = z0
        z = zu

        ! Iterate spray_param until usr is stabilized (depends on massflux).

        call spray_param(sourcestrength_tune,feedback_tune,u,zu,hss, &
                         hll,ts,ta,s,Patm,usr,hwave,cwave,p,h,       &
                         grav,kon,Q,tdk,Rdgas,Le,cpa,rhoa,cpw,rhow,  &
                         gamma1,beta,massf,U10,hs_tot,hl_tot,vfm)

        do i = 1,nit

           ! Droplet buoyancy

           ! Fall speed calculation

           k2 = 2.2*1e3*(rho_0/rhoa)**0.5         ! [cm^(1/2) s^(-1)]
           k2 = k2*1e-2                           ! [(m/s^2)^0.5] 

           omega = vfm/(kon*usr)                  ! From Lykossov (2001)

           !hsource = 0.05                         ! Set a height scale for the source (m)
           hsource = 0.01*h                       ! Set a height scale for the source (m)

           sigma = (rhow-rhoa)/rhoa
           S_r = massf/sigma/vfm/rhoa*(h/hsource)**omega ! Compute surface source from parameterization at wave height
           alpha = (Beta1*grav*kon**2*z_r*sigma*S_r)/(usr**2)

           if (omega /= 1) then

              S_bar1 = (omega**-1)*log(1+((alpha*omega**2)/ &
                       (1-omega))*((z/z_r)**(1-omega)-1))

           else

             S_bar1 = log(1+alpha*log(z/z_r))        

           endif

           usr_new = u*kon/(log(z/z0)+S_bar1)     ! u*
           cdzn = (usr_new/u)**2                  ! Drag coefficient

           hss = -999.0                           ! Sensible heat flux
           hll = -999.0                           ! Latent heat flux
          
           call spray_param(sourcestrength_tune,feedback_tune,u,zu,hss, &
                            hll,ts,ta,s,Patm,usr_new,hwave,cwave,p,h,   &
                            grav,kon,Q,tdk,Rdgas,Le,cpa,rhoa,cpw,rhow,  &
                            gamma1,beta,massf,U10,hs_tot,hl_tot,vfm)
        
        end do

        return

        end subroutine execute_spray_param

!*****************************************************************

        subroutine find_ust(u,kon,zu,rhoa,rhow,hwave,cwave,p, &
                            usr,h,z0)

          implicit none
          real :: z0
          real :: cd10n
          real :: cdhf     ! sqrt(drag coefficient)
          real :: cdzn     ! Drag coefficient
          real :: usr      ! Ustar (m/s)
          real :: hwave    ! Wave significant height (m)
          real :: cwave    ! Wave phase speed (m/s)
          real :: h        ! Wave top height above mean sea level (m)
          real :: U10      ! 10m wind speed
          real :: p        ! Energy input by breaking of waves (m/s)^3 right behind the breaker
          real :: u,kon,zu,rhoa,rhow,wab,cw10,aa,ab,ac
          integer :: i
          real, parameter :: grav = 9.83 
          real, parameter :: pi = 3.14159265359
          integer :: ustx
          real :: cEwave,twave,a,b,charn
          real :: ustw

#ifdef GEOS
       ustx = -4
#else
       ustx = -1
#endif
       ! -1 uses parametrized wave state (cwave, hwave, and p all parameterized)
       !    It also finds initial ustar to feed into spray parameterization
       ! -2 uses cwave from input and hwave and p parameterized
       ! -3 uses cwave and hwave from input and p parameterized
       ! -4 uses cwave, hwave, and p from input
       ! -5 uses cwave and p from input and hwave parameterized
       
       select case (ustx)

       case (-1)
       
         U10 = u*log(10/1e-3)/log(zu/1e-3)
         do i = 1,3  
            cwave = 7.+10.*U10/75.               ! Parameterized cwave
            wab = min(2.,cwave/U10)
            cw10 = U10
            aa = 0.241+1.9e-3*cw10-7e-6*cw10**2-0.123/exp(cw10/10)
            ab = 0.0797-10.5e-3*cw10+9.4e-5*cw10**2-0.181/exp(cw10/10)
            ac = -0.159+9.3e-3*cw10-9.1e-5*cw10**2+0.229/exp(cw10/10)
            cd10n = (0.94+0.06*(1.-exp(-U10/10.))) &
                   *(aa+ab*wab+ac*wab*wab)*1e-2
            z0 = 10./exp(kon/sqrt(cd10n))
            cdzn = (kon/log(zu/z0))**2
            cdhf = sqrt(cdzn)
            usr = u*cdhf
            U10 = usr/sqrt(cd10n)
         enddo
         hwave = 5.+10.*U10/80.                  ! Parameterized hwave
         h = hwave/2.                            ! Wave top height above mean sea level
         cEwave = 1.*(-0.4+0.25*U10)
         p = rhoa*cEwave*usr**2/rhow             ! Parameterized p
         if (p .lt. 0.0) then
            p = rhoa/rhow*3.5*usr**3.5           ! If p<0 we use the old parameterization
         endif

       case (-2)

         z0 = zu/exp(kon*u/usr)
         cd10n = (kon/log(10./z0))**2
         U10 = usr/sqrt(cd10n)

         hwave = 5.+10.*U10/80.                  ! parameterized hwave
         h = hwave/2
         cEwave = 1.*(-0.4+0.25*U10)
         p = rhoa*cEwave*usr**2/rhow             ! parameterized p
         if (p .lt. 0.0) then
            p = rhoa/rhow*3.5*usr**3.5
         endif

       case (-3)

         z0 = zu/exp(kon*u/usr)
         cd10n = (kon/log(10./z0))**2
         U10 = usr/sqrt(cd10n)

         h = hwave/2
         cEwave = 1.*(-0.4+0.25*U10)
         p = rhoa*cEwave*usr**2/rhow             ! parameterized p
         if (p .lt. 0.0) then
            p = rhoa/rhow*3.5*usr**3.5
         endif

       case (-4)

         z0 = zu/exp(kon*u/usr)
         cd10n = (kon/log(10./z0))**2
         U10 = usr/sqrt(cd10n)

         h = hwave/2

       case (-5)
     
         z0 = zu/exp(kon*u/usr)
         cd10n = (kon/log(10./z0))**2
         U10 = usr/sqrt(cd10n)
         hwave = 5.+10.*U10/80.                  ! parameterized hwave

         h = hwave/2
      
      end select
      
      return

      end subroutine find_ust

!************************************************************

        subroutine spray_param(sourcestrength,feedtune,Uz,zu,hss,   &
                         hll,ts,ta,s,Patm,ust,hwave,cwave,p,h,      &
                         grav,kon,Q,tdk,Rdgas,Le,cpa,rhoa,cpw,rhow, &
                         gamma1,beta,massf,U10,hs_tot,hl_tot,vfm)

        !Input
          !sourcestrength  1.0 = original Andreas and Fairall; reccommended ~ 0.3
          !Uz              Wind speed (m/s) at height zu
          !zu              Height (m) of bulk met data
          !hss             Sensible heat carried by turbulence from bulk algorithm (set to -999 if not available)
          !hll             Latent heat flux carried by water vapor from bulk algorithm (set to -999 if not available)
          !ts              Sea surface temperature (C) at zu
          !ta              Air temperature (C) at zu
          !s               Water vapor Saturation ratio (= RH (%)/100) at zu
          !Patm            Pressure (mb)
          !ust             Friction velocity (m/s) from wave model
          !hwave           Significant wave height (m)
          !cwave           Phase speed of breaking waves (m/s)
          !p               Wave energy dissipation (m/s)^3
          !h               Wave top height (m)
          !grav            Gravitational acceleration constant (m/s^2)
          !kon             von Karman constant
          !Q               Specific humidity (RH%/100)
          !tdk             Celsius to Kelvin
          !Rdgas           Gas constant, dry air
          !Le              Latent heat vaporization
          !cpa             Air specific heat
          !rhoa            Air density
          !cpw             Water specific heat
          !rhow            Water density
          !visa            Air kinematic viscosity

        !Output
          !hsd             Sensible heat carried by droplet mass flux
          !hld             Latent heat flux carried by droplet evaporation, including feedback effects
          !hss             Sensible heat carried by turbulence if not passed by input
          !hll             Latent heat flux carried by water vapor, if not passed by input
          !feed            Feedback coefficient (0 to 1) computed in the model, a small value means no feedback
          !massf           Mass flux of sea spray in kg/m^2/s
          !U10             Wind speed at 10-m height, used for other calculations
          !hs_tot          Total sensible heat flux realized above dropletlayer=hss+hsd_hs_epsilon+rho*cp*U*Ch*dtf
          !hl_tot          Total latent heat flux realized above dropletlayer=hll+hld-rho*cp*U*Ch*dtf
          !hs_eps          Heat energy dissipated in the atmospheric surface layer
          !hsde            Total droplet heat transfer before feedback
          !rmx             1-third power of third moment of droplet size spectral density (micron)
          !vmx             Mass weighted droplet fall velocity (m/s)

        !Comments
          !Evaporation in the droplet layer cools air by dtf and increases saturation by dqf
          !dtf is determined integrating the flux profile through the evaporation layer
          !Total enthalpy transfer to the atmosphere (neglecting dissipation) is hss+hll+hsdr
          !Total sensible tansfer is hss+hsdr-hldr+rho*Le*U*Ch*dtq
          !Total latent transfer is hll+hldr-rho*cp*U*Ch*dtf

        !References
          !Fairall et al, 1994, GLobal Atmos Ocean Syst., 2, 121-142
          !Bao et al., 2000, Mon Wea Rev., 128, 2190-2210

        implicit none

        !Input
        real :: sourcestrength, feedtune
        real :: Uz,zu,hss,hll,ts,ta,s,Patm,ust,hwave,cwave,p,h,grav,kon, &
                Q,tdk,Rdgas,Le,cpa,rhoa,cpw,rhow,gamma1,beta
        !Output
        real :: hsd,hld,tau,feed,massf,U10, &
                hs_tot,hl_tot,dtf,hs_eps,hsde

        !For U10 and bulk fluxes
        real :: ustr,z0,charn,cd10n,cdzn2,dti,dt
        real :: zot,rmy,visa,Rr
        
        !From drop_source
        real :: rm,vfm
        
        !Source area and volume
        real :: sfcc,sfcc2,sfcb,Sa,Sv
        
        !Fluxes
        real :: wetdep,qll0

        !For feedback
        !Flux gradient model
        real :: zl,zm,Gh4,Gq4,Zo4,qqs0_4,bow,qllp0,sp4,Hso4,Hlo4,hsde4, &
                Dtf0,Dqf0,hld4,Hst4,Hlt4,Dtf4,Dqf4,Qs4,Q4,bow4,hsd4,    &
                hsdl4,qllp4,feed4,hs_tot4,hl_tot4,dqf
        !Bao enthalpy model
        real :: baotime,dtime,hsde5,tw,hsde50,hld5,fe0,fe,numx05,     &
                denx05,tap0,dtf5,sp5,bow5,hsd5,hsdl5,numx5,denx5,tap, &
                qllp5,dqf5,feed5,hs_tot5,hl_tot5

        real, parameter :: pi = 3.14159265359
        real    :: ctznh,hh
        integer :: i,j
        integer :: feedt

!       print *, sourcestrength,feedtune,Uz,zu,hss
!       print *, hll,ts,ta,s,Patm,ust,hwave,cwave,p,h
!       print *, grav,kon,Q,tdk,Rdgas,Le,cpa,rhoa,cpw,rhow
!       print *, gamma1,beta,massf,U10,hs_tot,hl_tot,vfm

        visa = 1.326e-5*(1+6.542e-3*ta+8.301e-6*ta**2-4.84e-9*ta**3)

        !Estimate u10 and bulk fluxes
        if (ust > 0) then
          ustr = ust
          z0 = zu/exp(kon*Uz/ustr)
          z0 = max(z0,1e-4)
        else
          charn = 0.018
          ustr = 0.04*Uz
          do i = 1,10
            z0 = charn*ustr**2/grav
            cd10n = min(3.5e-3,(kon/log(10/z0))**2) ! cd is not allowed to exceed 3.5e-3
            z0 = 10/exp(kon/sqrt(cd10n))
            cdzn2 = (kon/log(zu/z0))**2
            ustr = Uz*sqrt(cdzn2)
          enddo
        endif
        tau = rhoa*ustr**2
        U10 = ustr*log(10/z0)/kon
        dt = ts-0.3-ta
        Rr = ustr*z0/visa
        zot = min(1.1e-5,5.5e-5*Rr**-0.6)                     ! coare 3.0 scalar roughtnes
        ctznh = 0.4/log(zu/zot)

        if (hll .lt. -900) then
           hll = rhoa*Le*ctznh*ustr*1e-3* &
                 (0.98*qsat(ts-0.3,Patm)-qsat(ta,Patm)*s)     ! Direct turbulent moisture flux, no feedback
        endif

        if (hss .lt. -900) then
           hss = rhoa*cpa*ctznh*dt*ustr                       ! Direct turbulent sensible, no feedback
        endif
!       print *, 'hss =', hss, 'hll = ', hll
        h = hwave/2 
        hh = max(0.03*h,10*z0)                                ! Gust level height
        vfm = (Uz+cwave/2-ustr/0.4*log(zu/hh))*0.07*1.15-0.3  ! Estimate mean fall velocity
#ifdef GEOS
        vfm = max(vfm,1e-6)
#else
        vfm = max(vfm,0.0)
#endif
        rm = (50*vfm**0.7+14)                                 ! Estimate mass mean radius (micron)
        rm = (150*vfm**0.8+5)                                 ! Estimate mass mean radius (micron)
        rmy = 55*vfm**0.7+20.
        Sv = 1.e-5*(1+(h/3)**0.1)*(rmy/50)**2.5               ! Normalized source volume m/s
        Sa = 1.2*(p/6e-4)**0.15*(rmy/73)**-1                  ! Normalized source area 1/s
        massf = sourcestrength*rhow*p*Sv

        ! Droplet and turbulent fluxes
        wetdep = (1-s)*(1-beta)/gamma1                        ! Wet bulb depression
        qll0 = sourcestrength*Sa*p*h*rhoa*Le*1e-3*qsat(ta,Patm)*beta & ! Estimate of qll from parameterization
               *(1-s)*(1-0.27*(1+1/(1-s))**(1.0/3.0))
        hsde = 0.92*cpw*massf*(dt+wetdep)                     ! Droplet enthalpy flux; 0.92 is loss of heat not
                                                              ! transfered from very large droplets
        hs_eps = 0.5*rhoa*ustr**3/kon*(log(h/10)+kon*U10/ustr)! Dissipation heating in droplet layer

        feedt = int(feedtune)

        select case (feedt)

          case (-999)
            zl = 0.1*h
            zm = h
            Gh4 = rhoa*cpa*0.4*ustr
            Gq4 = rhoa*Le*0.4*ustr
            Zo4 = 10/exp(kon**2/1.1e-3/log(10/z0))
            qqs0_4 = 0.92*sourcestrength*p*Sv*rhow*cpw*(dt+wetdep)
            bow = dt/wetdep
            qllp0 = qll0 + bow/(1+bow)*hsde
            sp4 = s
            Hso4 = hss
            Hlo4 = hll
            hsde4 = hsde
            Dtf0 = 0
            Dqf0 = 0
            do i = 1,20
              hld4 = sourcestrength*Sa*p*h*rhoa*Le*1e-3*qsat(ta,Patm)* &
                     beta*(1-sp4)*(1-0.27*(1+1/(1-sp4))**(1.0/3.0))
              Hst4 = Hso4+qqs0_4-hld4
              Hlt4 = Hlo4+hld4
              Dtf4 = -(Hst4-Hso4)*(1-zl/(zm-zl)*log(zm/zl))/Gh4
              Dqf4 = (Hlt4-Hlo4)*(1-zl/(zm-zl)*log(zm/zl))/Gq4
              Dtf0 = Dtf0-0.2*(Dtf0-Dtf4)
              Dqf0 = Dqf0-0.2*(Dqf0-Dqf4)
              Qs4 = qsat(ta-Dtf0,Patm)/1e3  ! Specific humidity (kg/kg)
              Q4 = Q+Dqf4
              sp4 = Q4/Qs4
              bow4 = (ts-(ta-Dtf4))/(wetdep-Dtf4)
              hsd4 = hsde4*bow4/(1+bow4)
              hsdl4 = hsde4-hsd4
            end do
            qllp4 = hld4+hsdl4
            feed4 = (qllp4-rhoa*Le*U10*1.1e-3*Dqf4)/qllp0
            hs_tot4 = hss+hsd4-hld4+hs_eps+rhoa*cpa*U10*1e-3*Dtf4
            hld4 = hld4+hsdl4
            hl_tot4 = hll+hld4-rhoa*Le*U10*1e-3*Dqf4

            hld = hld4
            hsd = hsd4
            hs_tot = hs_tot4
            hl_tot = hl_tot4
            feed = feed4
            dtf = Dtf4
            dqf = Dqf4

          case default
            baotime = 40*feedtune
            dtime = 0.5*h/vfm*baotime
            hsde5 = hsde
            bow = dt/wetdep
            tw = ta-wetdep
            hsde50 = bow/(1+bow)*hsde
            hld5 = sourcestrength*Sa*p*h*rhoa*Le*1e-3*qsat(ta,Patm)* &
                   beta*(1-s)*(1-0.27*(1+1/(1-s))**(1.0/3.0))
            qllp0 = qll0+bow/(1+bow)*hsde
            fe0 = hld5/Le/massf
            fe = fe0
            numx05 = cpw*massf/h*(dt+(1-fe)*(ta-tw))-Le*fe*massf/h
            denx05 = cpa*rhoa+cpa*fe*massf/h*dtime
            tap0 = -dtime*numx05/denx05
            tap = tap0
            dtf5 = 0
            do i = 1,20                           ! Feedback parameterization
              dtf5 = dtf5+0.1*(tap-dtf5)          ! Change in air temperature caused by droplet evaporation
              dtf5 = max(dtf5,-0.1) 
              sp5 = s+dtf5*gamma1/(1-beta)        ! Saturation ratio, after feedback effects
              !print *, s, dtf5, gamma1, beta
              sp5 = min(sp5,0.98)                 ! Doesn't allow to exceed seawater saturation
              bow5 = (ts-(ta-dtf5))/(wetdep-dtf5) ! Ratio of droplet sensible to evap heats to get to wet bulb temp
              hsd5 = hsde5*bow5/(1.0+bow5)
              hsdl5 = hsde5-hsd5
              hld5 = sourcestrength*Sa*p*h*rhoa*Le*1e-3* &
                     qsat(ta,Patm)*beta*(1-sp5)*         &
                     (1.0-0.27*(1.0+1.0/(1.0-sp5))**(1.0/3.0))
              fe = hld5/Le/massf
              numx5 = cpw*massf/h*(dt+dtf5+(1.0-fe)* &
                      (ta-dtf5-tw))-Le*fe*massf/h
              denx5 = cpa*rhoa+cpa*fe*massf/h*dtime
              tap = -dtime*numx5/denx5
            enddo
            qllp5 = hld5+hsdl5
            dqf5 = cpa/Le*dtf5
            feed5 = (qllp5-rhoa*Le*U10*1.1e-3*dqf5)/qllp0
            hs_tot5 = hss+hsd5-hld5+hs_eps+rhoa*cpa*U10*1e-3*dtf5 ! With this formulation the term -hld5 might
                                                                  ! make hs_tot5 too negative.
                                                                  ! The next line is the new formulation with -hld5 removed.
            ! hs_tot5 = hss+hsd5+hs_eps+rhoa*cpa*U10*1e-3*dtf5
            hld5 = hld5 + hsdl5           !lb note: Chirs's code is hld4 = hld5 + hsdl5 but I think it is wrong
            hl_tot5 = hll+hld5-rhoa*Le*U10*1e-3*dqf5

            hld = hld5
            hsd = hsd5
            hs_tot = hs_tot5
            hl_tot = hl_tot5
            feed = feed5
            dtf = dtf5
            dqf = dqf5
        end select

        return

        end subroutine spray_param

!**********************************************************

        real function qsat(ta,Patm)

          implicit none
          real :: ta, Patm
          real :: es
          es = 6.112*exp(17.502*ta/(ta+241.0))*(1.0007+3.46e-6*Patm)
          qsat = es*622/(Patm-0.378*es)

          return

        end function qsat

!************************************************************

end module bl_seaspray_mod
