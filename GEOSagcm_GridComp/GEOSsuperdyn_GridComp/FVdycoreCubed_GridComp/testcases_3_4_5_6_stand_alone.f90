MODULE testcases_3_4_5_6

  !=======================================================================
  !
  !  Functions for setting up initial conditions for the dynamical core test cases 3-6
  !  The input parameters depend on the test case (see comments for each test case below).
  !
  !  Given longitude (radians), latitude (radians), eta (pressure) and rotation_angle (degrees)
  !  the functions will return temperature, surface geopotential, zonal and meridional wind
  !  components, respectively.
  !
  !  Author: Christiane Jablonowski (University of Michigan, cjablono@umich.edu)
  !
  !  May/5/2008
  !
  !=======================================================================
  
implicit none
!=======================================================================
!  physical constants
!=======================================================================
  integer, parameter :: r8 = SELECTED_REAL_KIND(12) ! 8 byte real

  real(r8), parameter ::           &
       Rd         = 287.04_r8,     &             ! gas constant J/(K kg)
       cp         = 1004.64_r8,    &             ! specific heat at constant pressure J/(K kg)
       kappa      = Rd/cp,         &             ! kappa = 2/7
       g          = 9.80616_r8,    &             ! gravitational acceleration (m/s^2)
       a          = 6371229._r8,   &             ! Earth's radius in m
       pi         = 3.14159265358979323846_r8,&  ! pi
       omega      = 2._r8*pi/86164._r8, &        ! Earth's angular velocity 1/s
       pih        = pi*0.5_r8,     &             ! pi/2
       deg2rad    = pi/180._r8

CONTAINS

!==========================================================================================
! pure 3D advection, time-dependent
!==========================================================================================
  SUBROUTINE advection (tracer_variant, lon, lat, height, rotation_angle,  &
                        u_wind, v_wind, temperature, surface_geopotential, &
                        surface_pressure, q5, q6)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      character*2, intent(in) :: tracer_variant                          ! identifies test variant 'yy', here tracers
                                                                         ! e.g. 0 : no tracer, set to zero
                                                                         !      5 : tracer q5 only
                                                                         !     56 : both tracers q5 and q6
      real(r8), intent(in)  :: lon,                                     & ! longitude in radians
                               lat,                                     & ! latitude in radians
                               height,                                  & ! height of the level in m
                               rotation_angle                             ! alpha in degrees
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind,                                  & ! meridional wind in m/s
                               temperature,                             & ! temperature in K
                               surface_geopotential,                    & ! surface geopotential in m^2/s^2
                               surface_pressure,                        & ! surface pressure in Pa
                               q5,                                      & ! tracer q5
                               q6                                         ! tracer q6
!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
      real(r8), parameter :: p0      = 100000._r8,                      & ! reference pressure 
                             u0      = (2._r8*pi*a)/(12._r8*86400._r8), & ! circumference / 12 days
                             tau     = 3._r8 * 86400._r8,               & ! period: 3 days expressed in s
                             omega_0 = pi*40000._r8 / tau,              & ! 0.4848 Pa/s
                             T0      = 300._r8,                         & ! constant temperature
                             H       = Rd * T0 / g,                     & ! scale height
                             RR      = 1/3._r8,                         & ! horizontal half width divided by 'a'
                             ZZ      = 1000._r8,                        & ! vertical half width
                             z0      = 4500._r8,                        & ! center point in z
                             lambda0 = 1.5_r8*pi,                       & ! center point in longitudes
                             phi0    = 0._r8,                           & ! center point in latitudes
                             slot    = 1._r8/8._r8                        ! half width of the slot in radians
!----------------------------------------------------------------------- 
!     local variables
!-----------------------------------------------------------------------                             
      real(r8) :: alpha
      real(r8) :: sin_tmp, cos_tmp
      real(r8) :: d1, d2, s, r
      
      alpha = rotation_angle*deg2rad
!-----------------------------------------------------------------------
!    initialize the wind components
!-----------------------------------------------------------------------
      u_wind =  u0 * (cos(lat)*cos(alpha) + sin(lat)*cos(lon)*sin(alpha))
      v_wind = -u0 *  sin(lon) * sin(alpha)
!-----------------------------------------------------------------------
!     initialization of the vertical velocity: 
!     must be implemented in the dynamical core
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     initialize T (temperature)
!-----------------------------------------------------------------------
      temperature = T0
!-----------------------------------------------------------------------
!     initialize surface geopotential
!-----------------------------------------------------------------------
      surface_geopotential = 0._r8
!-----------------------------------------------------------------------
!     initialize PS (surface pressure)
!-----------------------------------------------------------------------
      surface_pressure = p0
!-----------------------------------------------------------------------
!     Tracer variables
!-----------------------------------------------------------------------
      q5 = 0._r8   ! default
      q6 = 0._r8   ! default
!-----------------------------------------------------------------------
!     tracer q5
!-----------------------------------------------------------------------
      if (tracer_variant(1:1) == '5' .or. tracer_variant(2:2) == '5') then
        sin_tmp = sin(lat) * sin(phi0)
        cos_tmp = cos(lat) * cos(phi0)
        r  = ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))       ! great circle distance without 'a'
        d1 = min( 1._r8, (r/RR)**2 + ((height-z0)/ZZ)**2 )
        q5 = 0.5_r8 * (1._r8 + cos(pi*d1))
      endif
!-----------------------------------------------------------------------
!     tracer q6
!-----------------------------------------------------------------------
      if (tracer_variant(1:1) == '6' .or. tracer_variant(2:2) == '6') then
        sin_tmp = sin(lat) * sin(phi0)
        cos_tmp = cos(lat) * cos(phi0)
        r  = ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))       ! great circle distance without 'a'
        d2 = (r/RR)**2 + ((height-z0)/ZZ)**2
        if (d2 <= 1._r8) then
          q6 = 1._r8
        else
          q6 = 0._r8
        endif
        if ((height > z0) .and. ((phi0-slot) < lat .and. lat < (phi0+slot)) ) q6 = 0._r8   ! slotted ellipse               
      endif
  end subroutine advection
 
!==========================================================================================
! Rossby_Haurwitz wave, wavenumber 4
!==========================================================================================
  SUBROUTINE Rossby_Haurwitz (lon, lat, pressure,                                &
                              u_wind, v_wind, temperature, surface_geopotential, &
                              surface_pressure)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(in)  :: lon,                                     & ! longitude in radians
                               lat,                                     & ! latitude in radians
                               pressure                                   ! pressure at full model level in Pa
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind,                                  & ! meridional wind in m/s
                               temperature,                             & ! temperature in K
                               surface_geopotential,                    & ! surface geopotential in m^2/s^2
                               surface_pressure                           ! surface pressure in Pa

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
      real(r8),parameter :: u0      = 50._r8,                          &   ! reference wind
                            T0      = 288._r8,                         &   ! reference temperature
                            n       = 4._r8,                           &   ! wavenumber
                            MM      = u0/(n*a),                        &   ! parameter M and p_ref=95500.
                            KK      = u0/(n*a),                        &   ! parameter K
                            gamma   = 0.0065_r8,                       &   ! lapse rate in K/m
                            p_ref   = 95500._r8                            ! reference pressure
                            
!----------------------------------------------------------------------- 
!     local
!----------------------------------------------------------------------- 
      real(r8) :: tmp1, tmp2, tmp3
      real(r8) :: sin_lat, cos_lat, sin_slat, cos_slat
      real(r8) :: exponent_1, exponent_2
      real(r8) :: AA, BB, CC
      real(r8) :: phis_perturb
      
!-----------------------------------------------------------------------
!     initialize the wind components
!-----------------------------------------------------------------------
      cos_lat = cos(lat)
      sin_lat = sin(lat)
      tmp1 = a * MM * cos_lat
      tmp2 = a * KK * cos_lat**(n-1._r8)*(n*sin_lat**2 - cos_lat**2)
      tmp3 = -a * KK * n * cos_lat**(n-1._r8) * sin_lat
      u_wind = tmp1 + tmp2 * cos(n*lon)
      v_wind = tmp3 * sin(n*lon)
!-----------------------------------------------------------------------
!     initialize surface geopotential
!-----------------------------------------------------------------------
      surface_geopotential = 0._r8     
!-----------------------------------------------------------------------
!     initialize surface pressure and temperature
!-----------------------------------------------------------------------
      tmp1       = gamma/(g*T0)
      tmp2       = a*a
      exponent_1 = g/(gamma*Rd)
      exponent_2 = (gamma*Rd)/g
      
      cos_lat = cos(lat)
      AA = tmp2 * (0.5_r8 * MM*(2._r8*omega+MM) * cos_lat**2 + 0.25_r8 * KK**2 * cos_lat**(2._r8*n) * &
                  ( (n+1._r8)*cos_lat**2 + (2._r8*n*n - n - 2._r8)) - 0.5_r8*n*n*KK**2 * cos_lat**(2._r8*(n-1)))
      BB = tmp2 * (2._r8*(omega+MM)*KK/((n+1._r8)*(n+2._r8)) * cos_lat**n * &
                   ( (n*n + 2._r8*n +2._r8) - (n+1._r8)**2 * cos_lat**2 ))
      CC = tmp2 * (0.25_r8 * KK**2 * cos_lat**(2._r8*n) * ( (n+1._r8)*cos_lat**2 - (n+2._r8)))
      phis_perturb = AA + BB * cos(n*lon) + CC * cos(2._r8*n*lon)
      surface_pressure = p_ref * (1._r8 + tmp1*phis_perturb)**exponent_1   ! surface pressure
      temperature      = T0 * (pressure/p_ref)**exponent_2                 ! temperature
      
  end subroutine Rossby_Haurwitz

!==========================================================================================
! Mountain induced Rossby wave
!==========================================================================================
  SUBROUTINE mountain_Rossby (choice, lon, lat, pressure,                        &
                              u_wind, v_wind, temperature, surface_geopotential, &
                              surface_pressure)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      integer, intent(in)   :: choice                                     ! identifies test variant 'x'
                                                                          ! e.g. 0 : gaussian hill mountian
                                                                          !      1 : cylindrical step-function hill
      real(r8), intent(in)  :: lon,                                     & ! longitude in radians
                               lat,                                     & ! latitude in radians
                               pressure                                   ! pressure at full model level in Pa
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind,                                  & ! meridional wind in m/s
                               temperature,                             & ! temperature in K
                               surface_geopotential,                    & ! surface geopotential in m^2/s^2
                               surface_pressure                           ! surface pressure in Pa
!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
      real(r8),parameter :: u0      = 20._r8,                          &   ! 20 m/s
                            T0      = 288._r8,                         &   ! temperature
                            N2      = g*g/(cp*T0),                     &   ! squared Brunt Vaisala frequency N^2
                            h0      = 2000._r8,                        &   ! amplitude of the mountain, 2km
                            d       = 1500.e3_r8,                      &   ! half width 1500 km
                            lambda0 = 0.5_r8*pi,                       &   ! center point in longitudes
                            phi0    = pi/6._r8,                        &   ! center point in latitudes
                            p_sp    = 93000._r8                            ! pressure at the South Pole in Pa
!-----------------------------------------------------------------------
!   local variables
!-----------------------------------------------------------------------
      real(r8) :: sin_tmp, cos_tmp
      real(r8) :: tmp1, tmp2, tmp3
      real(r8) :: r

!-----------------------------------------------------------------------
!    initialize the wind components
!-----------------------------------------------------------------------
      u_wind = u0 * cos(lat)
      v_wind = 0._r8
!-----------------------------------------------------------------------
!     initialize T (temperature)
!-----------------------------------------------------------------------
      temperature = T0
!-----------------------------------------------------------------------
!     initialize surface geopotential and surface pressure
!-----------------------------------------------------------------------
      tmp1 = (a * N2 * u0)/(2._r8 * g*g * kappa) * (u0/a + 2._r8 * omega)
      tmp2 = N2 / (g*g * kappa)
   
      sin_tmp = sin(lat) * sin(phi0)
      cos_tmp = cos(lat) * cos(phi0)
      tmp3 = tmp1*((sin(lat))**2 - 1._r8)
      r = a * ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))   ! great circle distance with 'a'
      surface_geopotential = g*h0 * exp(-(r/d)**2)        ! Gaussian profile of the mountain

      if (choice==1) then
        surface_geopotential = 0.0
        sin_tmp = sin(lat) * sin(phi0)
        cos_tmp = cos(lat) * cos(phi0)
        tmp3 = tmp1*((sin(lat))**2 - 1._r8)
        r = a * ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))   ! great circle distance with 'a'
        if (r<=d) surface_geopotential = g*h0               ! Cylindrical step-function profile of the mountain
      endif 

      surface_pressure     = p_sp * exp( -tmp3 - tmp2*surface_geopotential)

  end subroutine mountain_Rossby
  
!==========================================================================================
! gravity waves
!==========================================================================================
  SUBROUTINE gravity_wave (choice, lon, lat, height,                          &
                           u_wind, v_wind, temperature, surface_geopotential, &
                           surface_pressure, pressure)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      integer, intent(in)   :: choice                                     ! identifies test variant 'x'
                                                                          ! e.g. 0 : no Coriolis, N=0.01 1/s, u0=0  m/s
                                                                          !      1 : no Coriolis, N=0.01 1/s, u0=40 m/s
      real(r8), intent(in)  :: lon,                                     & ! longitude in radians
                               lat,                                     & ! latitude in radians
                               height                                     ! height of the layer in meters
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(r8), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind,                                  & ! meridional wind in m/s
                               temperature,                             & ! temperature in K
                               surface_geopotential,                    & ! surface geopotential in m^2/s^2
                               surface_pressure
      real(r8), optional, intent(out) :: pressure
!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
      real(r8),parameter :: p0      = 100000._r8,                      & ! reference pressure 
                            T0      = 300._r8,                         & ! reference temperature
                            RR      = a/3._r8,                         & ! half width   
                            Lz      = 20.e3_r8,                        & ! vertical wave length, 20 km 
                            delta_theta = 10._r8                         ! potential temperature perturbation amplitude                   

!----------------------------------------------------------------------- 
!     local variables
!----------------------------------------------------------------------- 
      real(r8) :: sin_tmp, cos_tmp
      real(r8) :: tmp1, tmp2, tmp3
      real(r8) :: theta                                                    ! potential temperature
      real(r8) :: theta_mean
      real(r8) :: pres                                                     ! pressure
      real(r8) :: r                                                        ! great circle distance

!-----------------------------------------------------------------------
!     more test case parameters
!----------------------------------------------------------------------- 
      real(r8) :: N2,                                                  &   ! squared Brunt-Vaisala frequency
                  S,                                                   &   ! parameter
                  u0,                                                  &   ! background wind speed
                  lambda0,                                             &   ! center point in longitudes 
                  phi0,                                                &   ! center point in latitudes
                  gw_omega,                                            &   ! rotation
                  p_eq
                  
!-----------------------------------------------------------------------
!    initialize parameters
!-----------------------------------------------------------------------
      lambda0 = pi                                                         ! center point in longitudes
      p_eq    = p0                                                         ! surface pressure at the equator

      select case (choice)
      case (0) 
        N2 = 1.e-4_r8                                                      ! squared Brunt Vaisala frequency N^2
        u0 = 0._r8                                                         ! background wind speed
        phi0 = 0._r8                                                       ! center point in latitudes (0 deg)
        gw_omega = 0._r8                                                   ! no rotation
      case (1) 
        N2 = (g*g)/(cp*T0)                                                 ! squared Brunt Vaisala frequency N^2
        u0 = 0._r8                                                         ! background wind speed
        phi0 = 0._r8                                                       ! center point in latitudes (0 deg)
        gw_omega = 0._r8                                                   ! no rotation
      case (2) 
        N2 = (g*g)/(cp*T0)                                                 ! squared Brunt Vaisala frequency N^2
        u0 = 40._r8                                                        ! background wind speed
        phi0 = 0._r8                                                       ! center point in latitudes (0 deg)
        gw_omega = 0._r8                                                   ! no rotation
      case (3) 
        N2 = (g*g)/(cp*T0)                                                 ! squared  Brunt Vaisala frequency N^2
        u0 = 0._r8                                                         ! background wind speed
        phi0 = pi/4._r8                                                    ! center point in latitudes (45 deg)
        gw_omega = omega                                                   ! Earth's rotation
      end select
      S = g*g/(cp*N2)                                                      ! parameter

!-----------------------------------------------------------------------
!     initialize the wind components
!-----------------------------------------------------------------------
      u_wind = u0 * cos(lat)
      v_wind = 0._r8
!-----------------------------------------------------------------------
!     initialize surface geopotential
!-----------------------------------------------------------------------
      surface_geopotential = 0._r8
!-----------------------------------------------------------------------
!     initialize surface pressure
!-----------------------------------------------------------------------
      tmp1 = (a * N2 * u0)/(2._r8 * g*g * kappa) * (u0/a + 2._r8 * gw_omega)
      surface_pressure = p_eq * exp( - tmp1*((sin(lat))**2) )
!-----------------------------------------------------------------------
!     initialize temperature
!-----------------------------------------------------------------------                 
    ! height  = -(g/(N2)) * log( (T0/S) * ( (pressure/p0)**kappa - 1._r8 ) + 1._r8 )
      pres    = p0 * ( (1._r8 - S/T0) + S/T0 * exp(- (N2*height)/g) )**(cp/Rd)
      if ( present(pressure) ) pressure = pres
      sin_tmp = sin(lat) * sin(phi0)
      cos_tmp = cos(lat) * cos(phi0)
      theta_mean = T0 /( T0/S * ((pres/p0)**kappa - 1._r8) + 1._r8 )
      r = a * ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))     ! great circle distance with radius
      if (r < RR) then
        tmp1 = 0.5_r8 * (1._r8 + cos(pi*r/RR))
      else
        tmp1 = 0._r8
      endif
      theta = theta_mean + delta_theta * tmp1 * sin(2._r8*pi*height/Lz)
      temperature = theta * (pres/p0)**kappa

  end subroutine gravity_wave


END MODULE testcases_3_4_5_6
