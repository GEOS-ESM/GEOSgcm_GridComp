!-*- F90 -*-
module FLOW_PROJ_mod
  !--------------------------------------------------------------------!
  ! author:  Michael Herzog                                            !
  ! email:   Michael.Herzog@noaa.gov                                   !
  ! date:    Feb 2007                                                  !
  ! version: 0.1                                                       !
  !                                                                    !
  ! routines for projection of flow vectors and components             !
  !--------------------------------------------------------------------!
  use fv_arrays_mod, only: REAL8

  implicit none

  private
  public :: d2a, d2a_vect, a2d_vect

contains
  !====================================================================!
  subroutine d2a(u, v, dx, dy, rdxa, rdya, cosa_s, &
       isd, ied, jsd, jed, ksd, ked, is, ie, js, je, ks, ke, ua, va)
    !------------------------------------------------------------------!
    ! D -> A: co-variant d-grid u,d to a-grid                          !
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed, ksd, ked, is, ie, js, je, ks, ke
    real(REAL8), dimension(isd:ied  ,jsd:jed+1,ksd:ked), intent(in) :: u
    real(REAL8), dimension(isd:ied+1,jsd:jed  ,ksd:ked), intent(in) :: v
    real(REAL8), dimension(isd:ied  ,jsd:jed+1), intent(in) :: dx
    real(REAL8), dimension(isd:ied+1,jsd:jed  ), intent(in) :: dy
    real(REAL8), dimension(isd:ied  ,jsd:jed  ), intent(in) :: rdxa, rdya, cosa_s
    real(REAL8), dimension(is :ie ,js :je ,ks :ke ), intent(out) :: ua, va
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8), dimension(is:ie  ,js:je+1) :: vt
    real(REAL8), dimension(is:ie+1,js:je  ) :: ut
    real(REAL8) :: utmp, vtmp, rsin2
    integer :: i,j,k
    !------------------------------------------------------------------!
    ! D -> A: co-variant u,d to contra-variant ua, va                  !
    ! vorticity preserving                                             !
    !------------------------------------------------------------------!
    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             vt(i,j) = u(i,j,k)*dx(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             ut(i,j) = v(i,j,k)*dy(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie
             utmp = 0.5*(vt(i,j) + vt(i,j+1))*rdxa(i,j)
             vtmp = 0.5*(ut(i,j) + ut(i+1,j))*rdya(i,j)
             rsin2 = 1. / (1.-cosa_s(i,j)*cosa_s(i,j))
             ua(i,j,k) = (utmp-vtmp*cosa_s(i,j)) * rsin2
             va(i,j,k) = (vtmp-utmp*cosa_s(i,j)) * rsin2
          enddo
       enddo
    enddo

  end subroutine d2a
  !====================================================================!
  !====================================================================!
  subroutine d2a_vect(u, v, dx, dy, rdxa, rdya, cosa_s, ec1, ec2,       &
       isd, ied, jsd, jed, ksd, ked, is, ie, js, je, ks, ke, va_xyz)
    !------------------------------------------------------------------!
    ! D -> A: co-variant d-grid u,d to flow vector va_xyz on a-grid    !
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed, ksd, ked, is, ie, js, je, ks, ke
    real(REAL8), dimension(isd:ied  ,jsd:jed+1,ksd:ked), intent(in) :: u
    real(REAL8), dimension(isd:ied+1,jsd:jed  ,ksd:ked), intent(in) :: v
    real(REAL8), dimension(isd:ied  ,jsd:jed+1), intent(in) :: dx
    real(REAL8), dimension(isd:ied+1,jsd:jed  ), intent(in) :: dy
    real(REAL8), dimension(isd:ied  ,jsd:jed  ), intent(in) :: rdxa, rdya, cosa_s
    real(REAL8), dimension(3,isd:ied,jsd:jed)  , intent(in) :: ec1, ec2

    real(REAL8), dimension(3,isd:ied,jsd:jed,ksd:ked), intent(out) :: va_xyz
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8), dimension(is:ie  ,js:je+1) :: vt
    real(REAL8), dimension(is:ie+1,js:je  ) :: ut
    real(REAL8) :: utmp, vtmp, ua, va, rsin2
    integer :: i,j,k
    !------------------------------------------------------------------!
    ! D -> A: co-variant u,d to contra-variant ua, va                  !
    ! vorticity preserving                                             !
    !------------------------------------------------------------------!
    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             vt(i,j) = u(i,j,k)*dx(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             ut(i,j) = v(i,j,k)*dy(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie
             utmp = 0.5*(vt(i,j) + vt(i,j+1))*rdxa(i,j)
             vtmp = 0.5*(ut(i,j) + ut(i+1,j))*rdya(i,j)
             rsin2 = 1. / (1.-cosa_s(i,j)*cosa_s(i,j))
             ua = (utmp-vtmp*cosa_s(i,j)) * rsin2
             va = (vtmp-utmp*cosa_s(i,j)) * rsin2
             !---------------------------------------------------------!
             ! va_xyz = ua * e1 + va * e2                              !
             !---------------------------------------------------------!
             va_xyz(1,i,j,k) = ua * ec1(1,i,j) + va * ec2(1,i,j)
             va_xyz(2,i,j,k) = ua * ec1(2,i,j) + va * ec2(2,i,j)
             va_xyz(3,i,j,k) = ua * ec1(3,i,j) + va * ec2(3,i,j)
          enddo
       enddo
    enddo

  end subroutine d2a_vect
  !====================================================================!
  subroutine a2d_vect(va_xyz, ew1, ew2, es1, es2,                       &
       isd, ied, jsd, jed, ksd, ked, is, ie, js, je, ks, ke,            &
       edge_interp, edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n, &
       west_edge, east_edge, south_edge, north_edge,                    &
       sw_corner, se_corner, nw_corner, ne_corner,                      &
       iwest, ieast, jsouth, jnorth, u, v)
    !------------------------------------------------------------------!
    ! A -> D: flow vector va_xyz on a-grid to co-variant u,d on d-grid !
    !------------------------------------------------------------------!
    use FILL_CORNER_mod, only: fill_xyz_corner, XDir, YDir
    
    integer, intent(in) :: isd, ied, jsd, jed, ksd, ked, is, ie, js, je, ks, ke
    integer, intent(in) :: iwest, ieast, jsouth, jnorth

    real(REAL8), dimension(3,isd:ied,jsd:jed,ksd:ked), intent(inout) :: va_xyz
    real(REAL8), dimension(3,isd:ied+1,jsd:jed), intent(in) :: ew1, ew2
    real(REAL8), dimension(3,isd:ied,jsd:jed+1), intent(in) :: es1, es2
    real(REAL8), dimension(isd:ied), intent(in) :: edge_vect_s, edge_vect_n
    real(REAL8), dimension(jsd:jed), intent(in) :: edge_vect_w, edge_vect_e

    real(REAL8), dimension(isd:ied  ,jsd:jed+1,ksd:ked), intent(out) :: u
    real(REAL8), dimension(isd:ied+1,jsd:jed  ,ksd:ked), intent(out) :: v

    logical, intent(in) ::  edge_interp,                                  &
                            west_edge, east_edge, south_edge, north_edge, &
                            sw_corner, se_corner, nw_corner, ne_corner
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(REAL8), dimension(:,:), allocatable :: v1d, vt1d
    real(REAL8) :: vx, vy, vz
    integer :: i, j, k, nghost, im2, jm2
    !------------------------------------------------------------------!
    ! calculate d-grid u                                               !
    !------------------------------------------------------------------!
    nghost=max(iwest-is, ie-ieast, jsouth-js, je-jnorth)
    call fill_xyz_corner(va_xyz, YDir, nghost,                          &
                         isd, ied, jsd, jed, ksd, ked, ks, ke,          &
                         sw_corner, se_corner, nw_corner, ne_corner,    &
                         iwest, ieast, jsouth, jnorth)
    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             vx = va_xyz(1,i,j-1,k) + va_xyz(1,i,j,k)
             vy = va_xyz(2,i,j-1,k) + va_xyz(2,i,j,k)
             vz = va_xyz(3,i,j-1,k) + va_xyz(3,i,j,k)
             
             u(i,j,k) = 0.5*(vx*es1(1,i,j) + vy*es1(2,i,j) + vz*es1(3,i,j))
          enddo
       enddo
    enddo
    if (edge_interp) then
       !---------------------------------------------------------------!
       ! fix south edge                                                !
       !---------------------------------------------------------------!
       allocate(v1d(3,isd:ied), vt1d(3,isd:ied))
       im2=(ied-isd)/2
       if (south_edge) then
          j=jsouth
          do k=ks,ke
             do i=is-1,ie+1
                v1d(1,i) = va_xyz(1,i,j-1,k) + va_xyz(1,i,j,k)
                v1d(2,i) = va_xyz(2,i,j-1,k) + va_xyz(2,i,j,k)
                v1d(3,i) = va_xyz(3,i,j-1,k) + va_xyz(3,i,j,k)
             enddo
             do i=is,ie
                if ( i<iwest .or. (i>im2 .and. i<=ieast) ) then
                   vt1d(1,i) = edge_vect_s(i)*v1d(1,i-1)+(1.-edge_vect_s(i))*v1d(1,i)
                   vt1d(2,i) = edge_vect_s(i)*v1d(2,i-1)+(1.-edge_vect_s(i))*v1d(2,i)
                   vt1d(3,i) = edge_vect_s(i)*v1d(3,i-1)+(1.-edge_vect_s(i))*v1d(3,i)
                else
                   vt1d(1,i) = edge_vect_s(i)*v1d(1,i+1)+(1.-edge_vect_s(i))*v1d(1,i)
                   vt1d(2,i) = edge_vect_s(i)*v1d(2,i+1)+(1.-edge_vect_s(i))*v1d(2,i)
                   vt1d(3,i) = edge_vect_s(i)*v1d(3,i+1)+(1.-edge_vect_s(i))*v1d(3,i)
                endif
             enddo
             ! Projection:
             do i=is,ie
                u(i,j,k) = 0.5*(vt1d(1,i)*es1(1,i,j)                       &
                               +vt1d(2,i)*es1(2,i,j)                       &
                               +vt1d(3,i)*es1(3,i,j))
             enddo
          enddo
       endif
       !---------------------------------------------------------------!
       ! fix the north edge                                            !
       !---------------------------------------------------------------!
       if (north_edge) then
          j=jnorth+1
          do k=ks,ke
             do i=is-1,ie+1
                v1d(1,i) = va_xyz(1,i,j-1,k) + va_xyz(1,i,j,k)
                v1d(2,i) = va_xyz(2,i,j-1,k) + va_xyz(2,i,j,k)
                v1d(3,i) = va_xyz(3,i,j-1,k) + va_xyz(3,i,j,k)
             enddo
             do i=is,ie
                if ( i<iwest .or. (i>im2 .and. i<=ieast) ) then
                   vt1d(1,i) = edge_vect_n(i)*v1d(1,i-1)+(1.-edge_vect_n(i))*v1d(1,i)
                   vt1d(2,i) = edge_vect_n(i)*v1d(2,i-1)+(1.-edge_vect_n(i))*v1d(2,i)
                   vt1d(3,i) = edge_vect_n(i)*v1d(3,i-1)+(1.-edge_vect_n(i))*v1d(3,i)
                else
                   vt1d(1,i) = edge_vect_n(i)*v1d(1,i+1)+(1.-edge_vect_n(i))*v1d(1,i)
                   vt1d(2,i) = edge_vect_n(i)*v1d(2,i+1)+(1.-edge_vect_n(i))*v1d(2,i)
                   vt1d(3,i) = edge_vect_n(i)*v1d(3,i+1)+(1.-edge_vect_n(i))*v1d(3,i)
                endif
             enddo
             ! Projection:
             do i=is,ie
                u(i,j,k) = 0.5*(vt1d(1,i)*es1(1,i,j)                       &
                               +vt1d(2,i)*es1(2,i,j)                       &
                               +vt1d(3,i)*es1(3,i,j))
             enddo
          enddo
       endif
       deallocate(v1d, vt1d)
    endif
    !------------------------------------------------------------------!
    ! calculate d-grid v                                               !
    !------------------------------------------------------------------!
    call fill_xyz_corner(va_xyz, XDir, nghost,                          &
                         isd, ied, jsd, jed, ksd, ked, ks, ke,          &
                         sw_corner, se_corner, nw_corner, ne_corner,    &
                         iwest, ieast, jsouth, jnorth)
    do k=ks,ke
       do j=js,je
          do i=is,ie+1
             vx = va_xyz(1,i-1,j,k) + va_xyz(1,i,j,k)
             vy = va_xyz(2,i-1,j,k) + va_xyz(2,i,j,k)
             vz = va_xyz(3,i-1,j,k) + va_xyz(3,i,j,k)
             
             v(i,j,k) = 0.5*(vx*ew2(1,i,j) + vy*ew2(2,i,j) + vz*ew2(3,i,j))
          enddo
       enddo
    enddo
    if (edge_interp) then
       !---------------------------------------------------------------!
       ! fix west edge                                                 !
       !---------------------------------------------------------------!
       allocate(v1d(3,jsd:jed), vt1d(3,jsd:jed))
       jm2=(jed-jsd)/2
       if (west_edge) then
          i=iwest
          do k=ks,ke
             do j=js-1,je+1
                v1d(1,j) = va_xyz(1,i-1,j,k) + va_xyz(1,i,j,k)
                v1d(2,j) = va_xyz(2,i-1,j,k) + va_xyz(2,i,j,k)
                v1d(3,j) = va_xyz(3,i-1,j,k) + va_xyz(3,i,j,k)
             enddo
             do j=js,je
                if ( j<jsouth .or. (j>jm2 .and. j<=jnorth) ) then
                   vt1d(1,j) = edge_vect_w(j)*v1d(1,j-1)+(1.-edge_vect_w(j))*v1d(1,j)
                   vt1d(2,j) = edge_vect_w(j)*v1d(2,j-1)+(1.-edge_vect_w(j))*v1d(2,j)
                   vt1d(3,j) = edge_vect_w(j)*v1d(3,j-1)+(1.-edge_vect_w(j))*v1d(3,j)
                else
                   vt1d(1,j) = edge_vect_w(j)*v1d(1,j+1)+(1.-edge_vect_w(j))*v1d(1,j)
                   vt1d(2,j) = edge_vect_w(j)*v1d(2,j+1)+(1.-edge_vect_w(j))*v1d(2,j)
                   vt1d(3,j) = edge_vect_w(j)*v1d(3,j+1)+(1.-edge_vect_w(j))*v1d(3,j)
                endif
             enddo
             ! Projection:
             do j=js,je
                v(i,j,k) = 0.5*(vt1d(1,j)*ew2(1,i,j)                       &
                               +vt1d(2,j)*ew2(2,i,j)                       &
                               +vt1d(3,j)*ew2(3,i,j))
             enddo
          enddo
       endif
       !---------------------------------------------------------------!
       ! fix east edge                                                 !
       !---------------------------------------------------------------!
       if (east_edge) then
          i=iwest+1
          do k=ks,ke
             do j=js-1,je+1
                v1d(1,j) = va_xyz(1,i-1,j,k) + va_xyz(1,i,j,k)
                v1d(2,j) = va_xyz(2,i-1,j,k) + va_xyz(2,i,j,k)
                v1d(3,j) = va_xyz(3,i-1,j,k) + va_xyz(3,i,j,k)
             enddo
             do j=js,je
                if ( j<jsouth .or. (j>jm2 .and. j<=jnorth) ) then
                   vt1d(1,j) = edge_vect_e(j)*v1d(1,j-1)+(1.-edge_vect_e(j))*v1d(1,j)
                   vt1d(2,j) = edge_vect_e(j)*v1d(2,j-1)+(1.-edge_vect_e(j))*v1d(2,j)
                   vt1d(3,j) = edge_vect_e(j)*v1d(3,j-1)+(1.-edge_vect_e(j))*v1d(3,j)
                else
                   vt1d(1,j) = edge_vect_e(j)*v1d(1,j+1)+(1.-edge_vect_e(j))*v1d(1,j)
                   vt1d(2,j) = edge_vect_e(j)*v1d(2,j+1)+(1.-edge_vect_e(j))*v1d(2,j)
                   vt1d(3,j) = edge_vect_e(j)*v1d(3,j+1)+(1.-edge_vect_e(j))*v1d(3,j)
                endif
             enddo
             ! Projection:
             do j=js,je
                v(i,j,k) = 0.5*(vt1d(1,j)*ew2(1,i,j)                       &
                               +vt1d(2,j)*ew2(2,i,j)                       &
                               +vt1d(3,j)*ew2(3,i,j))
             enddo
          enddo
       endif
       deallocate(v1d, vt1d)
    endif

  end subroutine a2d_vect
  !====================================================================!
end module FLOW_PROJ_mod
