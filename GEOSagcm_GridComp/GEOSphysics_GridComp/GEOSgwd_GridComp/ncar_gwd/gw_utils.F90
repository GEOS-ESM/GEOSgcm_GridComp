module gw_utils

!
! This module contains utility code for the gravity wave modules.
!

implicit none
private
save

! Real kind for gravity wave parameterization.
!!!integer, public, parameter :: GW_PRC = selected_real_kind(12)
integer,public,parameter :: GW_PRC = kind(1.d0) ! build precision

! Public interface
public :: get_unit_vector
public :: dot_2d
public :: midpoint_interp

contains

! Take two components of a vector, and find the unit vector components and
! total magnitude.
subroutine get_unit_vector(u, v, u_n, v_n, mag)
  real(GW_PRC), intent(in) :: u(:)
  real(GW_PRC), intent(in) :: v(:)
  real(GW_PRC), intent(out) :: u_n(:)
  real(GW_PRC), intent(out) :: v_n(:)
  real(GW_PRC), intent(out) :: mag(:)

  integer :: i

  mag = sqrt(u*u + v*v)

  ! Has to be a loop/if instead of a where, because floating point
  ! exceptions can trigger even on a masked divide-by-zero operation
  ! (especially on Intel).
  do i = 1, size(mag)
     if (mag(i) > 0._GW_PRC) then
        u_n(i) = u(i)/mag(i)
        v_n(i) = v(i)/mag(i)
     else
        u_n(i) = 0._GW_PRC
        v_n(i) = 0._GW_PRC
     end if
  end do

end subroutine get_unit_vector

! Vectorized version of a 2D dot product (since the intrinsic dot_product
! is more suitable for arrays of contiguous vectors).
function dot_2d(u1, v1, u2, v2)
  real(GW_PRC), intent(in) :: u1(:), v1(:)
  real(GW_PRC), intent(in) :: u2(:), v2(:)

  real(GW_PRC) :: dot_2d(size(u1))

  dot_2d = u1*u2 + v1*v2

end function dot_2d

! Pure function that interpolates the values of the input array along
! dimension 2. This is obviously not a very generic routine, unlike, say,
! CAM's lininterp. But it's used often enough that it seems worth providing
! here.
pure function midpoint_interp(arr) result(interp)
  real(GW_PRC), intent(in) :: arr(:,:)
  real(GW_PRC) :: interp(size(arr,1),size(arr,2)-1)

  integer :: i

  do i = 1, size(interp,2)
     interp(:,i) = 0.5_GW_PRC * (arr(:,i)+arr(:,i+1))
  end do

end function midpoint_interp

end module gw_utils
