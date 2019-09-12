!-*- F90 -*-
module FILL_CORNER_mod

  use fv_arrays_mod, only: REAL8

  implicit none

  private
  public :: XDir, YDir
  public :: fill_xyz_corner

  integer, parameter :: XDir=1, YDir=2

contains
  !====================================================================!
  subroutine fill_xyz_corner(va_xyz, dir, nghost,                        &
                             isd, ied, jsd, jed, ksd, ked, ks, ke,       &
                             sw_corner, se_corner, nw_corner, ne_corner, &
                             iwest, ieast, jsouth, jnorth)
    !------------------------------------------------------------------!
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed, ksd, ked, ks, ke
    integer, intent(in) :: dir, nghost, iwest, ieast, jsouth, jnorth
    logical, intent(in) :: sw_corner, se_corner, nw_corner, ne_corner
    real(REAL8), dimension(3,isd:ied,jsd:jed,ksd:ked), intent(inout) :: va_xyz
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: i, k

    do k=ks,ke
       do i=1,3
          call fill_2d_corner(va_xyz(i,:,:,k), dir, nghost, isd, ied, jsd, jed, &
                              sw_corner, se_corner, nw_corner, ne_corner,       &
                              iwest, ieast, jsouth, jnorth)
       enddo
    enddo

  end subroutine fill_xyz_corner
  !====================================================================!
  subroutine fill_2d_corner(var, dir, nghost, isd, ied, jsd, jed,       &
                            sw_corner, se_corner, nw_corner, ne_corner, &
                            iwest, ieast, jsouth, jnorth)
    !------------------------------------------------------------------!
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed
    integer, intent(in) :: dir, nghost, iwest, ieast, jsouth, jnorth
    logical, intent(in) :: sw_corner, se_corner, nw_corner, ne_corner
    real(REAL8), dimension(isd:ied,jsd:jed), intent(inout) :: var
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: i, j

    select case(dir)
    case(XDir)
       if ( sw_corner ) then
          do j=1,nghost
             do i=1,nghost
                var(iwest-i, jsouth-j) = var(iwest-j, i+jsouth-1)
             enddo
          enddo
       endif
       if ( se_corner ) then
          do j=1,nghost
             do i=1,nghost
                var(ieast+i, jsouth-j) = var(ieast+j, i+jsouth-1)
             enddo
          enddo
       endif
       if ( nw_corner ) then
          do j=1,nghost
             do i=1,nghost
                var(iwest-i, jnorth+j) = var(iwest-j, jnorth+1-i)
             enddo
          enddo
       endif
       if ( ne_corner ) then
          do j=1,nghost
             do i=1,nghost
                var(ieast+i, jnorth+j) = var(ieast+j, jnorth+1-i)
             enddo
          enddo
       endif
       
    case(YDir)
       
       if ( sw_corner ) then
          do j=1,nghost
             do i=1,nghost
                var(iwest-j, jsouth-i) = var(iwest-1+i, jsouth-j)
             enddo
          enddo
       endif
       if ( se_corner ) then
          do j=1,nghost
             do i=1,nghost
                var(ieast+j, jsouth-i) = var(ieast+1-i, jsouth-j)
             enddo
          enddo
       endif
       if ( nw_corner ) then
          do j=1,nghost
             do i=1,nghost
                var(iwest-j, jnorth+i) = var(iwest-1+i, jnorth+j)
             enddo
          enddo
       endif
       if ( ne_corner ) then
          do j=1,nghost
             do i=1,nghost
                var(ieast+j, jnorth+i) = var(ieast+1-i, jnorth+j)
             enddo
          enddo
       endif
       
    end select
  end subroutine fill_2d_corner
  !====================================================================!
end module FILL_CORNER_mod
