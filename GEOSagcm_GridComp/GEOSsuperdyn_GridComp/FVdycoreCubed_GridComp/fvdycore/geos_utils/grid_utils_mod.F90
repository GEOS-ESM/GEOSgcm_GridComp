!-*- F90 -*-
module GRID_UTILS_mod
  !--------------------------------------------------------------------!
  ! author:  Michael Herzog                                            !
  ! email:   Michael.Herzog@noaa.gov                                   !
  ! date:    Feb 2007                                                  !
  ! version: 0.1                                                       !
  !                                                                    !
  ! routines for grid calculations on the sphere                       !
  !--------------------------------------------------------------------!

  use fv_arrays_mod, only: REAL8

  implicit none

  private
  public :: latlon2xyz,   xyz2latlon,                                   &
            dist2side,    spherical_angle,                              &
            great_circle, unit_vect_latlon
  public :: get_dx, get_dy, get_dxa, get_dya,                           &
            get_center_vect, get_west_vect, get_south_vect,             &
            get_cosa_center, vect_cross, normalize_vect

  real(REAL8), parameter :: pi = 3.141592653589793,                            &
                     big_number  = 1.e+30,                              &
                     tiny_number = 1.e-30

contains
  !====================================================================!
  subroutine get_dx(xyz_corner, isd, ied, jsd, jed, is, ie, js, je, dx, rdx)
    !------------------------------------------------------------------!
    ! calculate normalized dx                                          !
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed, is, ie, js, je
    real(REAL8), dimension(3,isd:ied+1,jsd:jed+1), intent(in)  :: xyz_corner
    real(REAL8), dimension  (isd:ied  ,jsd:jed+1), intent(out) :: dx
    real(REAL8), dimension  (isd:ied  ,jsd:jed+1), optional, intent(out) :: rdx
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: i, j

    do j=js,je+1
       do i=is,ie
          dx(i,j)=great_circle(xyz_corner(:,i,j), xyz_corner(:,i+1,j))
       enddo
    enddo
    if (present(rdx)) then
       do j=js,je+1
          do i=is,ie
             rdx(i,j)=1./dx(i,j)
          enddo
       enddo
    endif

  end subroutine get_dx
  !====================================================================!
  subroutine get_dy(xyz_corner, isd, ied, jsd, jed, is, ie, js, je, dy, rdy)
    !------------------------------------------------------------------!
    ! calculate normalized dy                                          !
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed, is, ie, js, je
    real(REAL8), dimension(3,isd:ied+1,jsd:jed+1), intent(in)  :: xyz_corner
    real(REAL8), dimension  (isd:ied+1,jsd:jed  ), intent(out) :: dy
    real(REAL8), dimension  (isd:ied+1,jsd:jed  ), optional, intent(out) :: rdy
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: i, j

    do j=js,je
       do i=is,ie+1
          dy(i,j)=great_circle(xyz_corner(:,i,j), xyz_corner(:,i,j+1))
       enddo
    enddo
    if (present(rdy)) then
       do j=js,je
          do i=is,ie+1
             rdy(i,j)=1./dy(i,j)
          enddo
       enddo
    endif

  end subroutine get_dy
  !====================================================================!
  subroutine get_dxa(xyz_corner, isd, ied, jsd, jed, is, ie, js, je, dxa, rdxa)
    !------------------------------------------------------------------!
    ! calculate normalized dx                                          !
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed, is, ie, js, je
    real(REAL8), dimension(3,isd:ied+1,jsd:jed+1), intent(in)  :: xyz_corner
    real(REAL8), dimension  (isd:ied  ,jsd:jed  ), intent(out) :: dxa
    real(REAL8), dimension  (isd:ied  ,jsd:jed  ), optional, intent(out) :: rdxa
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: i, j
    real(REAL8), dimension(3,isd:ied+1,jsd:jed) :: mid_pt

    do j=js,je
       do i=is,ie+1
          mid_pt(:,i,j)=xyz_corner(:,i,j)+xyz_corner(:,i,j+1)
       enddo
    enddo
    do j=js,je
       do i=is,ie
          dxa(i,j)=great_circle(mid_pt(:,i,j), mid_pt(:,i+1,j))
       enddo
    enddo
    if (present(rdxa)) then
       do j=js,je
          do i=is,ie
             rdxa(i,j)=1./dxa(i,j)
          enddo
       enddo
    endif

  end subroutine get_dxa
  !====================================================================!
  subroutine get_dya(xyz_corner, isd, ied, jsd, jed, is, ie, js, je, dya, rdya)
    !------------------------------------------------------------------!
    ! calculate normalized dx                                          !
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed, is, ie, js, je
    real(REAL8), dimension(3,isd:ied+1,jsd:jed+1), intent(in)  :: xyz_corner
    real(REAL8), dimension  (isd:ied  ,jsd:jed  ), intent(out) :: dya
    real(REAL8), dimension  (isd:ied  ,jsd:jed  ), optional, intent(out) :: rdya
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: i, j
    real(REAL8), dimension(3,isd:ied,jsd:jed+1) :: mid_pt

    do j=js,je+1
       do i=is,ie
          mid_pt(:,i,j)=xyz_corner(:,i,j)+xyz_corner(:,i+1,j)
       enddo
    enddo
    do j=js,je
       do i=is,ie
          dya(i,j)=great_circle(mid_pt(:,i,j), mid_pt(:,i,j+1))
       enddo
    enddo
    if (present(rdya)) then
       do j=js,je
          do i=is,ie
             rdya(i,j)=1./dya(i,j)
          enddo
       enddo
    endif

  end subroutine get_dya
  !====================================================================!
  subroutine get_center_vect(xyz_corner, isd, ied, jsd, jed, is, ie, js, je, ec1, ec2)
    !------------------------------------------------------------------!
    ! calculate unity coordinate vectors for cell center               !
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed, is, ie, js, je
    real(REAL8), dimension(3,isd:ied+1,jsd:jed+1), intent(in)  :: xyz_corner
    real(REAL8), dimension(3,isd:ied  ,jsd:jed  ), intent(out) :: ec1, ec2
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: i, j, k

    do j=js,je
       do i=is,ie
          do k=1,3
            ec1(k,i,j)=xyz_corner(k,i+1,j)+xyz_corner(k,i+1,j+1)          &
                      -xyz_corner(k,i  ,j)-xyz_corner(k,i  ,j+1)
            ec2(k,i,j)=xyz_corner(k,i,j+1)+xyz_corner(k,i+1,j+1)          &
                      -xyz_corner(k,i,j  )-xyz_corner(k,i+1,j  ) 
          enddo
          call normalize_vect(ec1(:,i,j))
          call normalize_vect(ec2(:,i,j))
       enddo
    enddo
    
  end subroutine get_center_vect
  !====================================================================!
  subroutine get_cosa_center(ec1, ec2, isd, ied, jsd, jed, is, ie, js, je, cosa_s, sina_s)
    !------------------------------------------------------------------!
    ! cosine of angle between coordinate axes at cell center           !
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed, is, ie, js, je
    real(REAL8), dimension(3,isd:ied,jsd:jed), intent(in)  :: ec1, ec2
    real(REAL8), dimension  (isd:ied,jsd:jed), intent(out) :: cosa_s, sina_s
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8) :: sin2
    integer :: i, j

    do j=js,je
       do i=is,ie
          cosa_s(i,j)=inner_prod(ec1(:,i,j), ec2(:,i,j))
          sin2=1.-cosa_s(i,j)*cosa_s(i,j)
          if (sin2>0.) then
             sina_s(i,j)=sqrt(sin2)
          else
             cosa_s(i,j)=1.
             sina_s(i,j)=tiny_number
          endif
       enddo
    enddo

  end subroutine get_cosa_center
  !====================================================================!
  subroutine get_west_vect(xyz_corner, isd, ied, jsd, jed, is, ie, js, je, &
                           west_edge, east_edge, iwest, ieast, ew1, ew2)
    !------------------------------------------------------------------!
    ! calculate unity coordinate vectors for west cell boundary        !
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed, is, ie, js, je
    integer, intent(in) :: iwest, ieast
    real(REAL8), dimension(3,isd:ied+1,jsd:jed+1), intent(in)  :: xyz_corner
    real(REAL8), dimension(3,isd:ied+1,jsd:jed  ), intent(out) :: ew1, ew2
    logical, intent(in) :: west_edge, east_edge
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8), dimension(3) :: mid_pt
    integer :: i, j, iss

    !------------------------------------------------------------------!
    ! interior cells                                                   !
    !------------------------------------------------------------------!
    iss=max(is,isd+1)
    do j=js,je
       do i=iss,ie+1
          ew1(:,i,j) = - xyz_corner(:,i-1,j  ) + xyz_corner(:,i+1,j  )  &
                       - xyz_corner(:,i-1,j+1) + xyz_corner(:,i+1,j+1)
          mid_pt(:) = xyz_corner(:,i,j) + xyz_corner(:,i,j+1)
          call project_sphere_v(ew1(:,i,j), mid_pt(:))
          call normalize_vect(ew1(:,i,j))

          ew2(:,i,j) = - xyz_corner(:,i,j) + xyz_corner(:,i,j+1)
          call normalize_vect(ew2(:,i,j))
       enddo
    enddo
    !------------------------------------------------------------------!
    ! fix edges: one-sided difference to avoid kink                    !
    !------------------------------------------------------------------!
    if (west_edge) then
       i=iwest
       do j=js,je
          ew1(:,i,j) = - xyz_corner(:,i,j  ) + xyz_corner(:,i+1,j  )  &
                       - xyz_corner(:,i,j+1) + xyz_corner(:,i+1,j+1)
          mid_pt(:) = xyz_corner(:,i,j) + xyz_corner(:,i,j+1)
          call project_sphere_v(ew1(:,i,j), mid_pt(:))
          call normalize_vect(ew1(:,i,j))
       enddo
    endif
    if (east_edge) then
       i=ieast
       do j=js,je
          ew1(:,i,j) = - xyz_corner(:,i-1,j  ) + xyz_corner(:,i,j  )  &
                       - xyz_corner(:,i-1,j+1) + xyz_corner(:,i,j+1)
          mid_pt(:) = xyz_corner(:,i,j) + xyz_corner(:,i,j+1)
          call project_sphere_v(ew1(:,i,j), mid_pt(:))
          call normalize_vect(ew1(:,i,j))
       enddo
    endif

  end subroutine get_west_vect
  !====================================================================!
  subroutine get_south_vect(xyz_corner, isd, ied, jsd, jed, is, ie, js, je, &
                           south_edge, north_edge, jsouth, jnorth, es1, es2)
    !------------------------------------------------------------------!
    ! calculate unity coordinate vectors for south cell boundary       !
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed, is, ie, js, je
    integer, intent(in) :: jsouth, jnorth
    real(REAL8), dimension(3,isd:ied+1,jsd:jed+1), intent(in)  :: xyz_corner
    real(REAL8), dimension(3,isd:ied  ,jsd:jed+1), intent(out) :: es1, es2
    logical, intent(in) :: south_edge, north_edge
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8), dimension(3) :: mid_pt
    integer :: i, j, jss

    !------------------------------------------------------------------!
    ! interior cells                                                   !
    !------------------------------------------------------------------!
    jss=max(js,jsd+1)
    do j=jss,je+1
       do i=is,ie
          es2(:,i,j) = - xyz_corner(:,i  ,j-1) + xyz_corner(:,i  ,j+1)  &
                       - xyz_corner(:,i+1,j-1) + xyz_corner(:,i+1,j+1)
          mid_pt(:) = xyz_corner(:,i,j) + xyz_corner(:,i+1,j)
          call project_sphere_v(es2(:,i,j), mid_pt(:))
          call normalize_vect(es2(:,i,j))

          es1(:,i,j) = - xyz_corner(:,i,j) + xyz_corner(:,i+1,j)
          call normalize_vect(es1(:,i,j))
       enddo
    enddo
    !------------------------------------------------------------------!
    ! fix edges: one-sided difference to avoid kink                    !
    !------------------------------------------------------------------!
    if (south_edge) then
       j=jsouth
       do i=is,ie
          es2(:,i,j) = - xyz_corner(:,i  ,j) + xyz_corner(:,i  ,j+1)    &
                       - xyz_corner(:,i+1,j) + xyz_corner(:,i+1,j+1)
          mid_pt(:) = xyz_corner(:,i,j) + xyz_corner(:,i+1,j)
          call project_sphere_v(es2(:,i,j), mid_pt(:))
          call normalize_vect(es2(:,i,j))
       enddo
    endif
    if (north_edge) then
       j=jnorth
       do i=is,ie
          es2(:,i,j) = - xyz_corner(:,i  ,j-1) + xyz_corner(:,i  ,j)    &
                       - xyz_corner(:,i+1,j-1) + xyz_corner(:,i+1,j)
          mid_pt(:) = xyz_corner(:,i,j) + xyz_corner(:,i+1,j)
          call project_sphere_v(es2(:,i,j), mid_pt(:))
          call normalize_vect(es2(:,i,j))
       enddo
    endif

  end subroutine get_south_vect
  !====================================================================!
  subroutine unit_vect_latlon(sph, elon, elat)
    !------------------------------------------------------------------!
    ! calculate unit vector for latlon in cartesian coordinates        !
    !------------------------------------------------------------------!
    real(REAL8), intent(in)  :: sph(2)
    real(REAL8), intent(out) :: elon(3), elat(3)

    real(REAL8) :: sin_lon, cos_lon, sin_lat, cos_lat
    
    sin_lon = sin(sph(1))
    cos_lon = cos(sph(1))
    sin_lat = sin(sph(2))
    cos_lat = cos(sph(2))
    
    elon(1) = -sin_lon
    elon(2) =  cos_lon
    elon(3) =  0.
    
    elat(1) = -sin_lat*cos_lon
    elat(2) = -sin_lat*sin_lon
    elat(3) =  cos_lat
    
  end subroutine unit_vect_latlon
  !====================================================================!
  subroutine project_sphere_v(v1, v2)
    !------------------------------------------------------------------!
    ! project v1 such that it is perpendicular to v2                   !
    !------------------------------------------------------------------!
    real(REAL8), dimension(3), intent(inout) :: v1
    real(REAL8), dimension(3), intent(in)    :: v2
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8) :: prod
    
    prod=inner_prod(v1, v2)
    v1(:)=v1(:)-prod*v2(:)

  end subroutine project_sphere_v
  !====================================================================!
  subroutine normalize_vect(vector)
    !------------------------------------------------------------------!
    ! normalize vector                                                 !
    !------------------------------------------------------------------!
    real(REAL8), intent(inout) :: vector(3)
    real(REAL8) :: vec2

    vec2=vector(1)*vector(1)+vector(2)*vector(2)+vector(3)*vector(3)
    if (vec2>0.) then
       vec2=1./sqrt(vec2)
       vector(:)=vec2*vector(:)
    endif

  end subroutine normalize_vect
  !====================================================================!
  subroutine latlon2xyz(sph_coor, xyz_coor)
    !------------------------------------------------------------------!
    ! calculate cartesian coordinates from spherical coordinates       !
    !                                                                  !
    ! input:                                                           !
    ! sph_coor [rad]   latlon coordinate                               !
    !                                                                  !
    ! output:                                                          !
    ! xyz_coor [1]     normalized cartesian vector                     !
    !------------------------------------------------------------------!
    real(REAL8), dimension(2), intent(in)    :: sph_coor
    real(REAL8), dimension(3), intent(inout) :: xyz_coor

    xyz_coor(1) = cos(sph_coor(2)) * cos(sph_coor(1))
    xyz_coor(2) = cos(sph_coor(2)) * sin(sph_coor(1))
    xyz_coor(3) = sin(sph_coor(2))
    
  end subroutine latlon2xyz
  !====================================================================!
  subroutine xyz2latlon(xyz_coor, sph_coor)
    !------------------------------------------------------------------!
    ! calculate spherical coordinates from cartesian coordinates       !
    !                                                                  !
    ! input:                                                           !
    ! xyz_coor [1]     normalized cartesian vector                     !
    !                                                                  !
    ! output:                                                          !
    ! sph_coor [rad]   latlon coordinate                               !
    !------------------------------------------------------------------!

    real(REAL8), dimension(3), intent(in)    :: xyz_coor
    real(REAL8), dimension(2), intent(inout) :: sph_coor
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8) :: coslat, radius, pi2
    real(REAL8), parameter :: epsilon=1.e-7
    integer :: i,j

    pi2=8.*atan(1.)

    radius=sqrt(xyz_coor(1)*xyz_coor(1) &
               +xyz_coor(2)*xyz_coor(2) &
               +xyz_coor(3)*xyz_coor(3))

    sph_coor(1)=atan2(xyz_coor(2),xyz_coor(1))
    if (sph_coor(1)<0.) sph_coor(1)=sph_coor(1)+pi2
    sph_coor(2)=asin(xyz_coor(3)/radius)

  end subroutine xyz2latlon
  !====================================================================!
  function dist2side(v1, v2, point)
    !------------------------------------------------------------------!
    ! calculate shortest normalized distance on sphere                 !
    ! from point to straight line defined by v1 and v2                 !
    !------------------------------------------------------------------!
    real(REAL8) :: dist2side
    real(REAL8), dimension(3), intent(in) :: v1, v2, point
    
    real(REAL8) :: angle, side

    angle = spherical_angle(v1, v2, point)
    side  = great_circle(v1, point)
    dist2side = asin(sin(side)*sin(angle))

  end function dist2side
  !====================================================================!
  function spherical_angle(v1, v2, v3)
    !------------------------------------------------------------------!
    ! calculate spherical angle of a triangle formed by v1, v2 and v3  !
    ! at v1                                                            !
    !------------------------------------------------------------------!
    real(REAL8) :: spherical_angle
    real(REAL8), dimension(3), intent(in) :: v1, v2, v3
    real(REAL8) :: px, py, pz, qx, qy, qz, abs_p, abs_q

    ! vector product between v1 and v2
    px = v1(2)*v2(3) - v1(3)*v2(2)
    py = v1(3)*v2(1) - v1(1)*v2(3)
    pz = v1(1)*v2(2) - v1(2)*v2(1)
    ! vector product between v1 and v3
    qx = v1(2)*v3(3) - v1(3)*v3(2) 
    qy = v1(3)*v3(1) - v1(1)*v3(3) 
    qz = v1(1)*v3(2) - v1(2)*v3(1) 
    
    ! angle between p and q
    abs_p=px**2+py**2+pz**2
    abs_q=qx**2+qy**2+qz**2
    if (abs_p*abs_q==0.) then
       spherical_angle=0.
    else
       spherical_angle = (px*qx+py*qy+pz*qz)/sqrt(abs_p*abs_q)
       spherical_angle = sign(min(1.,abs(spherical_angle)),spherical_angle)
       spherical_angle = acos(spherical_angle)
    endif

  end function spherical_angle
  !====================================================================!
  function great_circle(v1, v2)
    !------------------------------------------------------------------!
    ! calculate normalized great circle distance between v1 and v2     ! 
    !------------------------------------------------------------------!
    real(REAL8) :: great_circle
    real(REAL8), dimension(3), intent(in) :: v1, v2

    great_circle=(v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3))                  &
           /sqrt((v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))                  &
                *(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3)))
    great_circle = sign(min(1.,abs(great_circle)),great_circle)
    great_circle=acos(great_circle)

  end function great_circle
  !====================================================================!
  function inner_prod(v1, v2)
    !------------------------------------------------------------------!
    ! calculate inner product between v1 and v2                        !
    !------------------------------------------------------------------!
    real(REAL8) :: inner_prod
    real(REAL8), dimension(3), intent(in) :: v1, v2

    inner_prod=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)

  end function inner_prod
  !====================================================================!

 subroutine vect_cross(e, p1, p2)
 real(REAL8), intent(in) :: p1(3), p2(3)
 real(REAL8), intent(out):: e(3)
!
! Perform cross products of 3D vectors: e = P1 X P2
! 
      e(1) = p1(2)*p2(3) - p1(3)*p2(2)
      e(2) = p1(3)*p2(1) - p1(1)*p2(3)
      e(3) = p1(1)*p2(2) - p1(2)*p2(1)

 end subroutine vect_cross

end module GRID_UTILS_mod
