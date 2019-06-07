
#include "Raster.h"

module leap_year
  
  implicit none
  
contains
  
  integer function days_in_month(year, month)
    
    ! return the number of days in a given month
    
    implicit none
    
    integer :: year, month
    
    integer, dimension(12), parameter :: days_in_month_leap = &
         (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    
    integer, dimension(12), parameter :: days_in_month_nonleap = &
         (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    
    if (is_leap_year(year)) then
       days_in_month = days_in_month_leap(month) 
    else
       days_in_month = days_in_month_nonleap(month) 
    end if
    
  end function days_in_month
    
  ! ------------------------------------------------------------------
  
  integer function days_in_year(year)
    
    ! return the number of days in a given year
    
    implicit none
    
    integer :: year
    
    if (is_leap_year(year)) then
       days_in_year = 366
    else
       days_in_year = 365
    end if
    
  end function days_in_year
  
  ! ------------------------------------------------------------------
  
  logical function is_leap_year(year)
    
    implicit none

    integer :: year
    
    if (mod(year,4) /= 0) then 
       is_leap_year = .false.
    else if (mod(year,400) == 0) then
       is_leap_year = .true.
    else if (mod(year,100) == 0) then 
       is_leap_year = .false.
    else
       is_leap_year = .true.
    end if
    
  end function is_leap_year

  ! ------------------------------------------------------------------
  
  integer function pentad_of_year(day_of_year, year)
    
    implicit none
    
    integer :: day_of_year, year
    
    ! determine pentad
    
    if ((is_leap_year(year)) .and. day_of_year>=59) then
       
       pentad_of_year = (day_of_year-2)/5+1
       
    else
       
       pentad_of_year = (day_of_year-1)/5+1
       
    end if
    
  end function pentad_of_year
  
end module leap_year


! ========= EOF =========================================================
  

