module clm_time_manager

   implicit none
   private

! Public methods

! gkw: this is just to get code to compile

   public ::&
      get_step_size,            &! return step size in seconds
      get_rad_step_size,        &! return radiation step size in seconds
      get_nstep,                &! return timestep number
      get_days_per_year          ! return the days per year for current year

contains

integer function get_step_size( dt )

  ! Return the step size in seconds.

 integer, optional, intent(in) :: dt  ! set to this time step

 integer, save :: dt_default = -999

 if ( present(dt) ) then
   dt_default = dt
 end if

 if(dt_default < 0) stop 'CN: dt_default = 0'
 get_step_size = dt_default 
  
end function get_step_size

integer function get_nstep()

  ! Return the timestep number.

 get_nstep     =    0 ! for FireMod 
  
end function get_nstep

integer function get_rad_step_size()

  ! Return the step size in seconds.

 get_rad_step_size = -999999999  ! gkw: to make sure this is not used
  
end function get_rad_step_size

integer function get_days_per_year( year )

 integer, optional, intent(in) :: year  ! current year

 integer, save :: curr_year = 1999
 logical :: is_leap_year

 if ( present(year) ) then
   curr_year = year
 end if

 if (mod(curr_year,4) /= 0) then 
    is_leap_year = .false.
  else if (mod(curr_year,400) == 0) then
    is_leap_year = .true.
  else if (mod(curr_year,100) == 0) then 
    is_leap_year = .false.
  else
    is_leap_year = .true.
 end if

!!!is_leap_year = .false. ! gkw: 71l test 20110920

 if(is_leap_year) then
   get_days_per_year = 366
  else
   get_days_per_year = 365
 endif
  
end function get_days_per_year

end module clm_time_manager
