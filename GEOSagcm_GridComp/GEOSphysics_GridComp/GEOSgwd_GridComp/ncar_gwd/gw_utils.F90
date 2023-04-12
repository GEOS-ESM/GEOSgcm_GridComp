module gw_utils

!
! This module contains utility code for the gravity wave modules.
!

implicit none
private

! Real kind for gravity wave parameterization.
integer,public,parameter :: GW_R4 = SELECTED_REAL_KIND(6,37)
integer,public,parameter :: GW_R8 = SELECTED_REAL_KIND(15,307)
integer,public,parameter :: GW_PRC = GW_R4

! Public interface

interface get_unit_vector
   module procedure get_unit_vector_r4
   module procedure get_unit_vector_r8
end interface get_unit_vector
public :: get_unit_vector

interface dot_2d
   module procedure dot_2d_r4
   module procedure dot_2d_r8
end interface dot_2d
public :: dot_2d

interface midpoint_interp
   module procedure midpoint_interp_r4
   module procedure midpoint_interp_r8
end interface midpoint_interp
public :: midpoint_interp

contains

! Take two components of a vector, and find the unit vector components and
! total magnitude.
subroutine get_unit_vector_r4(u, v, u_n, v_n, mag)
  real(GW_R4), intent(in) :: u(:)
  real(GW_R4), intent(in) :: v(:)
  real(GW_R4), intent(out) :: u_n(:)
  real(GW_R4), intent(out) :: v_n(:)
  real(GW_R4), intent(out) :: mag(:)

  integer :: i

  mag = sqrt(u*u + v*v)

  ! Has to be a loop/if instead of a where, because floating point
  ! exceptions can trigger even on a masked divide-by-zero operation
  ! (especially on Intel).
  do i = 1, size(mag)
     if (mag(i) > 0._GW_R4) then
        u_n(i) = u(i)/mag(i)
        v_n(i) = v(i)/mag(i)
     else
        u_n(i) = 0._GW_R4
        v_n(i) = 0._GW_R4
     end if
  end do

end subroutine get_unit_vector_r4

subroutine get_unit_vector_r8(u, v, u_n, v_n, mag)
  real(GW_R8), intent(in) :: u(:)
  real(GW_R8), intent(in) :: v(:)
  real(GW_R8), intent(out) :: u_n(:)
  real(GW_R8), intent(out) :: v_n(:)
  real(GW_R8), intent(out) :: mag(:)

  integer :: i

  mag = sqrt(u*u + v*v)

  ! Has to be a loop/if instead of a where, because floating point
  ! exceptions can trigger even on a masked divide-by-zero operation
  ! (especially on Intel).
  do i = 1, size(mag)
     if (mag(i) > 0._GW_R8) then
        u_n(i) = u(i)/mag(i)
        v_n(i) = v(i)/mag(i)
     else
        u_n(i) = 0._GW_R8
        v_n(i) = 0._GW_R8
     end if
  end do

end subroutine get_unit_vector_r8

! Vectorized version of a 2D dot product (since the intrinsic dot_product
! is more suitable for arrays of contiguous vectors).
function dot_2d_r4(u1, v1, u2, v2) result(dot_2d)
  real(GW_R4), intent(in) :: u1(:), v1(:)
  real(GW_R4), intent(in) :: u2(:), v2(:)

  real(GW_R4) :: dot_2d(size(u1))

  dot_2d = u1*u2 + v1*v2

end function dot_2d_r4

function dot_2d_r8(u1, v1, u2, v2) result(dot_2d)
  real(GW_R8), intent(in) :: u1(:), v1(:)
  real(GW_R8), intent(in) :: u2(:), v2(:)

  real(GW_R8) :: dot_2d(size(u1))

  dot_2d = u1*u2 + v1*v2

end function dot_2d_r8

! Pure function that interpolates the values of the input array along
! dimension 2. This is obviously not a very generic routine, unlike, say,
! CAM's lininterp. But it's used often enough that it seems worth providing
! here.
pure function midpoint_interp_r4(arr) result(interp)
  real(GW_R4), intent(in) :: arr(:,:)
  real(GW_R4) :: interp(size(arr,1),size(arr,2)-1)

  integer :: i

  do i = 1, size(interp,2)
     interp(:,i) = 0.5_GW_R4 * (arr(:,i)+arr(:,i+1))
  end do

end function midpoint_interp_r4

pure function midpoint_interp_r8(arr) result(interp)
  real(GW_R8), intent(in) :: arr(:,:)
  real(GW_R8) :: interp(size(arr,1),size(arr,2)-1)

  integer :: i

  do i = 1, size(interp,2)
     interp(:,i) = 0.5_GW_R8 * (arr(:,i)+arr(:,i+1))
  end do

end function midpoint_interp_r8

end module gw_utils
