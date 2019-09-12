MODULE jw
  IMPLICIT NONE
  !
  !  Functions for setting up initial conditions for the Jablonowski-Williamson test case.
  !
  !  Given longitude (radians), latitude (radians), eta (pressure) and rotation_angle (degrees)
  !  the functions will return temperature, surface geopotential, zonal and meridional wind
  !  components, respectively.
  !
  !  lperturb=.FALSE. result in initial conditions for the steady-state test case.
  !  lperturb=.TRUE.  result in initial conditions for the baroclinic wave test case.
  !
  !     T   : FUNCTION temperature         (lon,lat,eta,rotation_angle)
  !     PHIS: FUNCTION surface_geopotential(lon,lat,rotation_angle)
  !     U   : FUNCTION u_wind              (lon,lat,eta,lperturb,rotation_angle)
  !     V   : FUNCTION v_wind              (lon,lat,eta,lperturb,rotation_angle)
  !
  !  The non-rotated (rotation_angle=0) version of the test cases is described in: 
  !
  !                 Jablonowski, C., and D. L. Williamson, 2006: A baroclinic instability 
  !                 test case for atmospheric model dynamical cores. 
  !                 Quart. J. Roy. Meteor. Soc., 132, 2943-2975.
  !
  !                 Jablonowski, C., and D. L. Williamson, 2006: A Baroclinic Wave Test Case 
  !                 for Dynamical Cores of General Circulation Models: Model Intercomparisons, 
  !                 NCAR Technical Note, NCAR/TN-469+STR, 89 pp. 
  !  
  !  The rotated version simply rotates the initial conditions so that the spherical coordinate
  !  poles do not conicide with the earth's rotation axis. Thereby the Coriolis parameter is
  !  a function of latitude and longitude:
  !
  !      f = 2*Omega*(-cos(lon)*cos(lat)*sin(rotation_angle)+sin(lat)*cos(rotation_angle))
  !
  !  where Omega = 7.292 x 10E-5/s and rotation_angle is the angle between the flow direction
  !  and equator.
  !
  !  Author: Peter Hjort Lauritzen (NCAR, pel@ucar.edu)
  !
  INTEGER,PARAMETER :: r8 = SELECTED_REAL_KIND(12) ! 8 byte real

  REAL(r8), PARAMETER ::       &
       eta_strato = 0.2d0     ,&! tropopause level
       u0         = 35.d0     ,&! 35 m/s
       T0         = 288.d0    ,&! global mean T at surface
       p0         = 100000.d0 ,&! global mean surface pressure
       eta0       = 0.252d0   ,&! center of jets (hybrid)
       !
       radius                 = 10.d0,             & ! radius of the perturbation
       perturbation_amplitude =  1.d0,             & ! amplitude of u perturbation 1 m/s
       perturbation_longitude = 20.d0,             & ! longitudinal position, 20E
       perturbation_latitude  = 40.d0,             & ! latitudinal position, 40N
       perturbation_latitude_tracer = 55.d0,       &
       !
       !
       !
       Rd         = 287.d0    ,&! gas constant J/(K kg)
       g          = 9.80616d0 ,&! gravitational acceleration (m/s^2)
       a          = 6371229.d0,&! Earth's radius in m
       omega      = 7.29212d-5,&! angular velocity 1/s
       gamma      = 0.005d0   ,&! lapse rate
       pi         = 3.14159265358979323846_R8,&  ! pi
       deg2rad    = pi/180.d0

CONTAINS
  !
  !********************************************************************
  !
  ! Temperature (equation (6) in Jablonowski and Williamson, 2006)
  !
  !********************************************************************
  !
  REAL(r8) FUNCTION temperature(lon,lat,eta,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta, lon, lat, rotation_angle
    REAL(r8)             :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF

    temperature  = t_mean(eta) + t_deviation(rot_lon,rot_lat,eta)
  END FUNCTION temperature
  !
  ! Horizontally averaged temperature (equation (4) and (5) in Jablonowski and Williamson (2006))
  !
  REAL(r8) FUNCTION t_mean(eta)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta
    REAL(r8)             :: exponent, delta_T
    
    exponent = Rd*gamma/g
    delta_T  = 480000.d0  ! in K, for T mean calculation
    IF (eta.gt.(eta_strato)) THEN
       t_mean = T0*eta**exponent                                ! mean temperature at each level (troposphere)
    ELSE
       t_mean = T0*eta**exponent + delta_T*(eta_strato-eta)**5  ! mean temperature at each level (stratosphere)
    ENDIF
  END FUNCTION t_mean
  !
  ! Temperature deviation from the horizontal mean 
  ! (equation (6) minus horizontally averaged temperature)
  !
  REAL(r8) FUNCTION t_deviation(lon,lat,eta)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta, lon, lat
    REAL(r8)             :: factor, phi_vertical, rot_lon, rot_lat, a_omega

    factor       = eta*pi*u0/Rd             
    phi_vertical = (eta - eta0) * 0.5d0*pi
    a_omega      = a*omega

    rot_lon = lon
    rot_lat = lat

    t_deviation = factor * 1.5d0 * SIN(phi_vertical) * (cos(phi_vertical))**0.5d0 *                        &
                  ((-2.d0*(SIN(rot_lat))**6 * ((COS(rot_lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*              &
                  u0 * (COS(phi_vertical))**1.5d0  +                                                       &
                  (8.d0/5.d0*(COS(rot_lat))**3 * ((SIN(rot_lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega*0.5d0 )
  END FUNCTION t_deviation
  !
  !**************************************************************************
  !
  ! Surface geopotential (equaiton (7) in Jablonowski and Williamson, 2006)
  !
  !**************************************************************************
  !  
  REAL(r8) FUNCTION surface_geopotential(lon,lat,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: lon, lat, rotation_angle
    REAL(r8)             :: cos_tmp, eta_sfc, rot_lon, rot_lat, a_omega

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF

    eta_sfc    = 1.d0
    cos_tmp    = u0 * (cos((eta_sfc-eta0)*pi*0.5d0))**1.5d0
    a_omega    = a*omega

    surface_geopotential = ((-2.d0*(SIN(rot_lat))**6 * ((COS(rot_lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*COS_tmp   &
                 + (8.d0/5.d0*(COS(rot_lat))**3 * ((SIN(rot_lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega)*COS_tmp
  END FUNCTION surface_geopotential
  !
  !********************************************************************
  !
  ! wind components (equation 2 in Jablonowski and Williamson, 2006)
  !
  !********************************************************************
  !  
  REAL(r8) FUNCTION u_wind(lon,lat,eta,lperturb,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: lon,lat,eta,rotation_angle
    LOGICAL, INTENT(IN)  :: lperturb
    REAL(r8) :: cos_lat, u_lat, phi_vertical, rot_lon, rot_lat, sin_tmp, cos_tmp, r, u_perturb, v_lat
    REAL(r8) :: perturb_lon, perturb_lat, v_tmp

    perturb_lon = perturbation_longitude*deg2rad
    perturb_lat = perturbation_latitude*deg2rad

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF

    phi_vertical = (eta - eta0) *0.5d0*pi
    u_lat = (COS(phi_vertical))**1.5d0 * 4.d0 * u0 * (sin(rot_lat))**2 * (cos(rot_lat))**2
    u_wind = u_lat

    IF (lperturb) THEN

       sin_tmp = SIN(perturb_lat)*SIN(rot_lat)
       cos_tmp = COS(perturb_lat)*COS(rot_lat)
                  
       r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-perturb_lon) )    ! great circle distance
       u_perturb = perturbation_amplitude*EXP(- (r*radius)**2 )
       u_lat     = u_perturb + u_lat
    ENDIF
    IF (ABS(rotation_angle)<1.0E-8) THEN
       u_wind = u_lat
    ELSE
       v_lat = 0.0d0
       !
       ! rotate wind components
       !
       CALL turnwi(u_lat,v_lat, u_wind,v_tmp,lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,-1)
       IF (ABS(u_wind)<1.0E-10) u_wind=0.0d0
    ENDIF
  END FUNCTION u_wind

  REAL(r8) FUNCTION v_wind(lon,lat,eta,lperturb,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: lon,lat,eta,rotation_angle
    LOGICAL, INTENT(IN)  :: lperturb
    REAL(r8) :: cos_lat, u_lat, phi_vertical, rot_lon, rot_lat, sin_tmp, cos_tmp, r, u_perturb, v_lat
    REAL(r8) :: perturb_lon, perturb_lat, u_tmp

    perturb_lon = perturbation_longitude*deg2rad
    perturb_lat = perturbation_latitude*deg2rad

    IF (ABS(rotation_angle)<1.0E-8) THEN
       v_wind = 0.0d0
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)


       phi_vertical = (eta - eta0) *0.5d0*pi
       u_lat = (COS(phi_vertical))**1.5d0 * 4.d0 * u0 * (sin(rot_lat))**2 * (cos(rot_lat))**2
 
       IF (lperturb) THEN
          
          sin_tmp = SIN(perturb_lat)*SIN(rot_lat)
          cos_tmp = COS(perturb_lat)*COS(rot_lat)
          
          r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-perturb_lon) )    ! great circle distance
          u_perturb = perturbation_amplitude*EXP(- (r*radius)**2 )
          u_lat     = u_perturb + u_lat
       ENDIF

       v_lat = 0.0d0
       !
       ! pole point velocities are not well-defined
       !
       IF (ABS(pi*0.5d0-lat)<1.0E-8.OR.ABS(pi*0.5d0+lat)<1.0E-8) THEN
          v_wind = 0.0d0
       ELSE
          !
          ! rotate wind components
          !
          CALL turnwi(u_lat,v_lat, u_tmp,v_wind,lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,-1)
       ENDIF
    ENDIF
  END FUNCTION v_wind
  !
  !******************************************************************************
  !
  ! Subroutines for rotation
  !
  !******************************************************************************
  !

  SUBROUTINE regrot(pxreg,pyreg,pxrot,pyrot,pxcen,pycen,kcall)
    IMPLICIT NONE
!
!----------------------------------------------------------------------
!
!*    conversion between regular and rotated spherical coordinates.
!*
!*    pxreg     longitudes of the regular coordinates
!*    pyreg     latitudes of the regular coordinates
!*    pxrot     longitudes of the rotated coordinates
!*    pyrot     latitudes of the rotated coordinates
!*              all coordinates given in degrees n (negative for s)
!*              and degrees e (negative values for w)
!*    pxcen     regular longitude of the south pole of the rotated grid
!*    pycen     regular latitude of the south pole of the rotated grid
!*
!*    kcall=-1: find regular as functions of rotated coordinates.
!*    kcall= 1: find rotated as functions of regular coordinates.
!
!-----------------------------------------------------------------------
!
      integer kxdim,kydim,kx,ky,kcall
      real(r8) :: pxreg,pyreg,&
                  pxrot,pyrot,&
                  pxcen,pycen
!
!-----------------------------------------------------------------------
!
      real(r8) zsycen,zcycen,zxmxc,zsxmxc,zcxmxc,zsyreg,zcyreg, &
               zsyrot,zcyrot,zcxrot,zsxrot,zpi,zpih
      integer jy,jx

      zpih = pi*0.5d0
!
      !----------------------------------------------------------------------
!
      zsycen = SIN((pycen+zpih))
      zcycen = COS((pycen+zpih))
!
      IF (kcall.eq.1) then
!
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         zsyrot = zcycen*zsyreg - zsycen*zcyreg*zcxmxc
         zsyrot = max(zsyrot,-1.d0)
         zsyrot = min(zsyrot,+1.d0)
         !
         pyrot = ASIN(zsyrot)
         !
         zcyrot = COS(pyrot)
         zcxrot = (zcycen*zcyreg*zcxmxc +zsycen*zsyreg)/zcyrot
         zcxrot = max(zcxrot,-1.d0)
         zcxrot = min(zcxrot,+1.d0)
         zsxrot = zcyreg*zsxmxc/zcyrot
         !
         pxrot = ACOS(zcxrot)
         !
         IF (zsxrot<0.0) pxrot = -pxrot
               !
      ELSEIF (kcall.eq.-1) then
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         zsyreg = zcycen*zsyrot + zsycen*zcyrot*zcxrot
         zsyreg = max(zsyreg,-1.d0)
         zsyreg = min(zsyreg,+1.d0)
         !
         pyreg = ASIN(zsyreg)
         !
         zcyreg = COS(pyreg)
         zcxmxc = (zcycen*zcyrot*zcxrot -&
              zsycen*zsyrot)/zcyreg
         zcxmxc = max(zcxmxc,-1.d0)
         zcxmxc = min(zcxmxc,+1.d0)
         zsxmxc = zcyrot*zsxrot/zcyreg
         zxmxc  = ACOS(zcxmxc)
         IF (zsxmxc<0.0) zxmxc = -zxmxc
         !
         pxreg = zxmxc + pxcen
         !
      ELSE
         WRITE(6,'(1x,''invalid kcall in regrot'')')
         STOP
      ENDIF
    END SUBROUTINE regrot

    SUBROUTINE turnwi(puarg,pvarg,pures,pvres,         &
                      pxreg,pyreg,pxrot,pyrot,   &
                      pxcen,pycen,kcall)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!*    turn horizontal velocity components between regular and
!*    rotated spherical coordinates.
!
!*    puarg : input u components
!*    pvarg : input v components
!*    pures : output u components
!*    pvres : output v components
!*    pa    : transformation coefficients
!*    pb    :    -"-
!*    pc    :    -"-
!*    pd    :    -"-
!*    pxreg : regular longitudes
!*    pyreg : regular latitudes
!*    pxrot : rotated longitudes
!*    pyrot : rotated latitudes
!*    kxdim              : dimension in the x (longitude) direction
!*    kydim              : dimension in the y (latitude) direction
!*    kx                 : number of gridpoints in the x direction
!*    ky                 : number of gridpoints in the y direction
!*    pxcen              : regular longitude of the south pole of the
!*                         transformed grid
!*    pycen              : regular latitude of the south pole of the
!*                         transformed grid
!*
!*    kcall < 0          : find wind components in regular coordinates
!*                         from wind components in rotated coordinates
!*    kcall > 0          : find wind components in rotated coordinates
!*                         from wind components in regular coordinates
!*    note that all coordinates are given in degrees n and degrees e.
!*       (negative values for s and w)
!
!-----------------------------------------------------------------------

      integer kxdim,kydim,kx,ky,kcall
      real(r8) puarg,pvarg,    &
               pures,pvres,    &
               pa,   pb,       &
               pc,   pd,       &
               pxreg,pyreg,    &
               pxrot,pyrot
      real(r8) pxcen,pycen
!
!-----------------------------------------------------------------------
!
      integer jy,jx
      real(r8) zpih,zsyc,zcyc,zsxreg,zcxreg,zsyreg,zcyreg,zxmxc,&
               zsxmxc,zcxmxc,zsxrot,zcxrot,zsyrot,zcyrot
!
!-----------------------------------------------------------------------
!
      IF (kcall.eq.1) then
         zpih = pi*0.5d0
         zsyc = SIN(pycen+zpih)
         zcyc = COS(pycen+zpih)
         !
         zsxreg = SIN(pxreg)
         zcxreg = COS(pxreg)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         !
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         !
         pa = zcyc*zsxmxc*zsxrot + zcxmxc*zcxrot
         pb = zcyc*zcxmxc*zsyreg*zsxrot - zsyc*zcyreg*zsxrot - &
              zsxmxc*zsyreg*zcxrot
         pc = zsyc*zsxmxc/zcyrot
         pd = (zsyc*zcxmxc*zsyreg + zcyc*zcyreg)/zcyrot
         !
         pures = pa*puarg + pb*pvarg
         pvres = pc*puarg + pd*pvarg
      ELSEIF (kcall.eq.-1) then
         zpih = pi*0.5d0
         zsyc = SIN(pycen+zpih)
         zcyc = COS(pycen+zpih)
         !
         zsxreg = SIN(pxreg)
         zcxreg = COS(pxreg)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         !
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         !
         pa = zcxmxc*zcxrot + zcyc*zsxmxc*zsxrot
         pb = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot -&
              zcxmxc*zsxrot*zsyrot
         pc =-zsyc*zsxrot/zcyreg
         pd = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg
         !
         pures = pa*puarg + pb*pvarg
         pvres = pc*puarg + pd*pvarg
      ELSE
         write(6,'(1x,''invalid kcall in turnwi'')')
         STOP
      ENDIF
    END SUBROUTINE turnwi  

!********************************************************************
!
! Tracers
!
!********************************************************************

!-----------------------------------------------------------------------
! Tracer q1 and q2
!-----------------------------------------------------------------------
  REAL(r8) FUNCTION tracer_q1_q2(lon,lat,eta,rotation_angle, eta_c)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta, lon, lat, rotation_angle, eta_c
    REAL(r8) :: rot_lon, rot_lat, sin_tmp, cos_tmp, r
    REAL(r8) :: rot_perturb_lon, rot_perturb_lat, tmp

    rot_perturb_lon = perturbation_longitude*deg2rad
    rot_perturb_lat = perturbation_latitude_tracer *deg2rad

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF
    sin_tmp = SIN(rot_perturb_lat)*SIN(rot_lat)
    cos_tmp = COS(rot_perturb_lat)*COS(rot_lat)
    r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-rot_perturb_lon) )    ! great circle distance

    tmp = EXP(- ((r*radius)**2 + ((eta-eta_c)/0.1_r8)**2))
    IF (ABS(tmp)<1.0E-8) tmp = 0.0
    tracer_q1_q2 = tmp
  END FUNCTION tracer_q1_q2

!-----------------------------------------------------------------------
! Tracer q3
!-----------------------------------------------------------------------
  REAL(r8) FUNCTION tracer_q3(lon,lat,eta,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta, lon, lat, rotation_angle
    REAL(r8) :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF
    tracer_q3 = 0.5_r8 * ( tanh( 3._r8*abs(rot_lat)-pi ) + 1._r8)

  END FUNCTION tracer_q3

!-----------------------------------------------------------------------
! Tracer q, absolute value of the relative vorticity of the unperturbed initial state
!           multiplied by 10^5
!-----------------------------------------------------------------------
  REAL(r8) FUNCTION tracer_q(lon,lat,eta,rotation_angle)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: eta, lon, lat, rotation_angle
    REAL(r8) :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF
    tracer_q = abs(-4._r8 * u0/a * (cos((eta-eta0)*pi*0.5_r8))**1.5_r8 * sin(rot_lat) * &
               cos(rot_lat) * (2._r8-5._r8*(sin(rot_lat))**2)) * 1.e5_r8
    if (tracer_q < 1.e-9_r8) tracer_q = 0._r8  !  otherwise error in netcdf file

  END FUNCTION tracer_q

END MODULE jw
