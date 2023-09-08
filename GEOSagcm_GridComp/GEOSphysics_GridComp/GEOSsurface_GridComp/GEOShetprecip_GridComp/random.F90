module random
  use MAPL, only: MAPL_PI
  implicit none
  private
!  public init_random_seed
!  public random_stduniform
  public random_normal

contains

#if 0
!ALT: the subroutine below is a cut-and-paste from this link:
!https://stackoverflow.com/questions/37304793/random-number-generator-produces-same-sequence-even-though-its-seeded

 subroutine init_random_seed()
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if

       pid = getpid()
       t = ieor(t, int(pid, kind(t)))

       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
  end subroutine init_random_seed

! the next 2 subroutines are cut-and-paste from 
! https://masuday.github.io/fortran_tutorial/random.html

subroutine random_stduniform(u)
   implicit none
   real,intent(out) :: u
   real :: r
   call random_number(r)
   u = 1 - r
end subroutine random_stduniform

subroutine random_stdnormal(x)
   implicit none
   real,intent(out) :: x
   real,parameter :: pi=3.14159265
   real :: u1,u2
   call random_stduniform(u1)
   call random_stduniform(u2)
   x = sqrt(-2*log(u1))*cos(2*pi*u2)
end subroutine random_stdnormal

#else

!AT adopted from https://rosettacode.org/wiki/Random_numbers#Fortran
subroutine random_normal(array)
  implicit none
  real :: array(:)

  integer :: n

  integer :: i
  real, parameter :: pi=MAPL_PI
  real :: temp

  n = size(array)
  call random_number(array) ! Uniform distribution
  array = 1.0 - array ! make it 0 < x <= 1
 
! Now convert to normal distribution
  DO i = 1, n-1, 2
    temp = SQRT(-2.0*LOG(array(i))) * COS(2*pi*array(i+1))
    array(i+1) = SQRT(-2.0*LOG(array(i))) * SIN(2*pi*array(i+1))
    array(i) = temp
  END DO
 
end subroutine random_normal
#endif

end module random
