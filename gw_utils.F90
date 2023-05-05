module gw_utils

!
! This module contains utility code for the gravity wave modules.
!

implicit none
private
save

! Real kind for gravity wave parameterization.
integer, public, parameter :: r8 = selected_real_kind(12)

! Public interface
public :: get_unit_vector
public :: dot_2d
public :: midpoint_interp
public :: midpoint_interp_s

contains

! Take two components of a vector, and find the unit vector components and
! total magnitude.
subroutine get_unit_vector(u, v, u_n, v_n, mag)
!$acc routine gang
  real(r8), intent(in) :: u(:)
  real(r8), intent(in) :: v(:)
  real(r8), intent(out) :: u_n(:)
  real(r8), intent(out) :: v_n(:)
  real(r8), intent(out) :: mag(:)

  integer :: i

!   mag = sqrt(u*u + v*v)
!$acc loop gang vector
  do i = 1, size(u)
      mag(i) = sqrt(u(i)*u(i) + v(i)*v(i))
  enddo

  ! Has to be a loop/if instead of a where, because floating point
  ! exceptions can trigger even on a masked divide-by-zero operation
  ! (especially on Intel).
!$acc loop gang vector
  do i = 1, size(mag)
     if (mag(i) > 0._r8) then
        u_n(i) = u(i)/mag(i)
        v_n(i) = v(i)/mag(i)
     else
        u_n(i) = 0._r8
        v_n(i) = 0._r8
     end if
  end do

end subroutine get_unit_vector

! Vectorized version of a 2D dot product (since the intrinsic dot_product
! is more suitable for arrays of contiguous vectors).
function dot_2d(u1, v1, u2, v2)
!$acc routine seq
  real(r8), intent(in) :: u1(:), v1(:)
  real(r8), intent(in) :: u2(:), v2(:)

  real(r8) :: dot_2d(size(u1))

  dot_2d = u1*u2 + v1*v2

end function dot_2d

! Pure function that interpolates the values of the input array along
! dimension 2. This is obviously not a very generic routine, unlike, say,
! CAM's lininterp. But it's used often enough that it seems worth providing
! here.
pure function midpoint_interp(arr) result(interp)
!$acc routine seq
  real(r8), intent(in) :: arr(:,:)
  real(r8) :: interp(size(arr,1),size(arr,2)-1)

  integer :: i

!!$acc loop vector
  do i = 1, size(interp,2)
     interp(:,i) = 0.5_r8 * (arr(:,i)+arr(:,i+1))
  end do
!!$acc end loop

end function midpoint_interp

subroutine midpoint_interp_s(arr, interp)
!$acc routine gang
   real(r8), intent(in)  :: arr(:,:)
   real(r8), intent(out) :: interp(:,:)

   integer :: i , j
!$acc loop gang vector collapse(2)
   do i = 1, size(interp,2)
      do j = 1, size(interp,1)
        interp(j,i) = 0.5_r8 * (arr(j,i)+arr(j,i+1))
      enddo
   enddo

end subroutine
end module gw_utils

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”

! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”

! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at

! http://www.apache.org/licenses/LICENSE-2.0

! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.