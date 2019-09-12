MODULE sw
  IMPLICIT NONE
  INTEGER,PARAMETER :: r8 = SELECTED_REAL_KIND(12) ! 8 byte real

  REAL(r8), PARAMETER ::       &
       alpha      = 0.        ,&! angle of axis rotation about the poles
       !
       !
       g          = 9.80616d0 ,&! gravitational acceleration (m/s^2)
       a          = 6371229.d0,&! Earth's radius in m
       omega      = 7.29212d-5,&! angular velocity 1/s
       pi         = 3.14159265358979323846_R8,&  ! pi
       deg2rad    = pi/180.d0

CONTAINS
  !
  !**************************************************************************
  !
  ! Surface geopotential
  !
  !**************************************************************************
  !  
  REAL(r8) FUNCTION surface_geopotential(lon,lat,test_case)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: lon, lat
    INTEGER , INTENT(IN) :: test_case
    REAL(r8) :: r, r0, p1(2), p2(2) 

      select case (test_case)

      case(1)

         surface_geopotential = 0.0

      case(2)

         surface_geopotential = 0.0

      case(5)

         r0 = pi/9.
         p1(1) = pi/2.
         p1(2) = pi/6.
         p2(1) = lon
         p2(2) = lat
         r = SQRT(MIN(r0*r0, (p2(1)-p1(1))*(p2(1)-p1(1)) + (p2(2)-p1(2))*(p2(2)-p1(2))))
         surface_geopotential = 2000.0*g*(1.0-(r/r0))

      end select

  END FUNCTION surface_geopotential
  !
  !**************************************************************************
  !
  ! heights
  !
  !**************************************************************************
  !  
  REAL(r8) FUNCTION height(lon,lat,test_case)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: lon, lat
    INTEGER , INTENT(IN) :: test_case
    REAL(r8) :: Ubar
    REAL(r8) :: r, r0, p1(2), p2(2)

      select case (test_case)

      case(1)

         r0 = a/3.
         p1(1) = pi/2.
         p1(2) = 0.0
         p2(1) = lon
         p2(2) = lat
         r = a*acos( sin(p1(2))*sin(p2(2)) + cos(p1(2))*cos(p2(2))*cos(p2(1)-p1(1)) )
         if (r < r0) then
            height = (1000.0/2.0)*(1.0+cos(pi*r/r0))
         else
            height = 0.0
         endif

      case(2)

         Ubar = (2.0*pi*a)/(12.0*86400.0)
         height = 2.94e4 - (a*omega*Ubar + (Ubar*Ubar)/2.) * &
                   ( -1.*cos(lon)*cos(lat)*sin(alpha) + &
                                  sin(lat)*cos(alpha) ) ** 2

      case(5)

         Ubar = (2.0*pi*a)/(12.0*86400.0)
         height = 5400.*g - (a*omega*Ubar + (Ubar*Ubar)/2.) * &
                   ( -1.*cos(lon)*cos(lat)*sin(alpha) + &
                                  sin(lat)*cos(alpha) ) ** 2

      end select

  END FUNCTION height
  !
  !********************************************************************
  !
  ! wind components
  !
  !********************************************************************
  !  
  REAL(r8) FUNCTION u_wind(lon,lat,test_case)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: lon,lat
    INTEGER, INTENT(IN)  :: test_case
    REAL(r8) :: Ubar

         Ubar = (2.0*pi*a)/(12.0*86400.0)
         u_wind = Ubar*( cos(lat)*cos(alpha) + cos(lon)*sin(lat)*sin(alpha))

  END FUNCTION u_wind

  REAL(r8) FUNCTION v_wind(lon,lat,test_case)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: lon,lat
    INTEGER, INTENT(IN)  :: test_case
    REAL(r8) :: Ubar

         Ubar = (2.0*pi*a)/(12.0*86400.0)
         v_wind = -Ubar*sin(lon)*sin(alpha)

  END FUNCTION v_wind

END MODULE sw
