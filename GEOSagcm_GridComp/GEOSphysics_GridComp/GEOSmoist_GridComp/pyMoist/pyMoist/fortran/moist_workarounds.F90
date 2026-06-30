module moist_dsl_workarounds

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! This module provides workarounds for memory passing when
! beyond the scope of the NDSL python bridge.
! Fortran API:
!  - CNV_Tracers_To_SOA: move the ProcessLibrary CNV_Tracers to a SOA for DSL consumption
!  - CNV_Tracers_To_AOS: copy back the tracers value part to the AOS structure the original code expects
! C API:
!  - C__get_CNV_Tracers_SOA__size: reach for CNV_Tracers 4th dimensions size from C
!  - C__get_CNV_Tracers_SOA__Q   : reach for Q as a 4D array
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

use GEOSmoist_Process_Library
use iso_c_binding
use MAPL

implicit none
private

real, allocatable, dimension (:,:,:,:), target :: CNV_Tracers_SOA__Q
real, allocatable, dimension (:), target       :: CNV_Tracers_SOA__fscav
real, allocatable, dimension (:, :), target    :: CNV_Tracers_SOA__Vect_Hcts
logical, allocatable, dimension (:), target    :: CNV_Tracers_SOA__use_gcc_washout

public :: CNV_Tracers_To_SOA, CNV_Tracers_To_AOS

contains

subroutine CNV_Tracers_To_SOA()
    integer :: IM, JM, LM, n, ncnst

    ncnst = size(CNV_Tracers)

    IM = size(CNV_Tracers(1)%Q,1)
    JM = size(CNV_Tracers(1)%Q,2)
    LM = size(CNV_Tracers(1)%Q,3)

    if ( .not. allocated(CNV_Tracers_SOA__Q) ) then
        allocate(CNV_Tracers_SOA__Q(IM, JM, LM, ncnst))
        allocate(CNV_Tracers_SOA__fscav(ncnst))
        allocate(CNV_Tracers_SOA__Vect_Hcts(ncnst, 4))
        allocate(CNV_Tracers_SOA__use_gcc_washout(ncnst))
        do n = 1, ncnst
            CNV_Tracers_SOA__fscav(n) = CNV_Tracers(n)%fscav
            CNV_Tracers_SOA__Vect_Hcts(n, :) = CNV_Tracers(n)%Vect_Hcts(:)
            CNV_Tracers_SOA__use_gcc_washout(n) = CNV_Tracers(n)%use_gcc_washout
        enddo
        call WRITE_PARALLEL("CNV_Tracers_SOA ready")
    endif

    do n = 1, ncnst
        CNV_Tracers_SOA__Q(:, :, :, n) = CNV_Tracers(n)%Q(:,:,:)
    enddo
end subroutine

subroutine CNV_Tracers_To_AOS()
    integer :: IM, JM, LM, n, ncnst

    ncnst = size(CNV_Tracers)

    IM = size(CNV_Tracers(1)%Q,1)
    JM = size(CNV_Tracers(1)%Q,2)
    LM = size(CNV_Tracers(1)%Q,3)

    do n = 1, ncnst
        CNV_Tracers(n)%Q(:,:,:)  = CNV_Tracers_SOA__Q(:, :, :, n)
    enddo
end subroutine

function C__get_CNV_Tracers_SOA__size() result(size_) bind(c, name="get_CNV_Tracers_SOA__size")
    integer(c_int) :: size_
    size_ = size(CNV_Tracers)
end function

function C__get_CNV_Tracers_SOA__Q() result(Q_as_C_pointer) bind(c, name="get_CNV_Tracers_SOA__Q")
    type(c_ptr) :: Q_as_C_pointer
    Q_as_C_pointer=c_loc(CNV_Tracers_SOA__Q)
end function

function C__get_CNV_Tracers_SOA__fscav() result(fscav_as_C_pointer) bind(c, name="get_CNV_Tracers_SOA__fscav")
    type(c_ptr) :: fscav_as_C_pointer
    fscav_as_C_pointer=c_loc(CNV_Tracers_SOA__fscav)
end function

function C__get_CNV_Tracers_SOA__Vect_Hcts() result(Vect_Hcts_as_C_pointer) bind(c, name="get_CNV_Tracers_SOA__Vect_Hcts")
    type(c_ptr) :: Vect_Hcts_as_C_pointer
    Vect_Hcts_as_C_pointer=c_loc(CNV_Tracers_SOA__Vect_Hcts)
end function

function C__get_CNV_Tracers_SOA__use_gcc_washout() result(use_gcc_washout_as_C_pointer) bind(c, name="get_CNV_Tracers_SOA__use_gcc_washout")
    type(c_ptr) :: use_gcc_washout_as_C_pointer
    use_gcc_washout_as_C_pointer=c_loc(CNV_Tracers_SOA__use_gcc_washout)
end function

end module moist_dsl_workarounds
