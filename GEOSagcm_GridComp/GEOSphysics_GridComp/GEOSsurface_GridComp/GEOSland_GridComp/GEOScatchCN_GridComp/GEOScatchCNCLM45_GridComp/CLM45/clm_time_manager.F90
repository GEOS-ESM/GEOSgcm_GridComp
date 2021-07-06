module clm_time_manager

   use update_model_para4cn, only: curr_year,curr_month,curr_day,curr_dofyr,curr_hour,curr_min,curr_sec
   use clm_varctl  , only: iulog

   implicit none
   private

! Public methods

! gkw: this is just to get code to compile
   
   public ::&
      get_step_size,            &! return step size in seconds
      get_rad_step_size,        &! return radiation step size in seconds
      get_nstep,                &! return CN timestep number

      get_curr_date,            &! return date components at end of current timestep
!      get_start_date,           &! return components of the start date
!      get_driver_start_ymd,     &! return year/month/day (as integer in YYYYMMDD format) of driver start date
!      get_ref_date,             &! return components of the reference date
!      get_curr_time,            &! return components of elapsed time since reference date at end of current timestep
      get_curr_calday,          &! return calendar day at end of current timestep
      get_calday,               &! return calendar day from input date
!      get_calendar,             &! return calendar
      
      get_days_per_year,        &! return the days per year for current year
      
      is_end_curr_day,          &! return true on last timestep in current day
      is_restart                 ! return true if this is a restart run

contains

!=========================================================================================

integer function get_step_size( dt )

  ! Return the step size in seconds.

 integer, optional, intent(in) :: dt  ! set to this time step

 integer, save :: dt_default = -999

 if ( present(dt) ) then
   dt_default = dt
 end if

 if(dt_default < 0) stop 'CN: dt_default < 0'
 get_step_size = dt_default 
  
end function get_step_size

!=========================================================================================

integer function get_nstep(istep)

  ! Return the timestep number.
  
  integer*8, optional, intent(in) :: istep
  
  integer, save :: istep_default = -999
  
  if ( present(istep) ) then
   istep_default = istep
  end if

 if(istep_default < 0) stop 'CN: istep_default < 0'
 get_nstep = istep_default  ! for FireMod
  
end function get_nstep

!=========================================================================================

integer function get_rad_step_size()

  ! Return the step size in seconds.

 get_rad_step_size = -999999999  ! gkw: to make sure this is not used
  
end function get_rad_step_size

!=========================================================================================

subroutine get_curr_date(yr, mon, day, tod)

  ! Return date components valid at end of current timestep with an optional
  ! offset (positive or negative) in seconds.  
  
  implicit none
  integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)
       
  yr  = curr_year
  mon = curr_month
  day = curr_day
  tod = 3600*curr_hour + 60*curr_min + curr_sec
  
end subroutine get_curr_date

!=========================================================================================

function get_curr_calday()

  ! Return calendar day at end of current timestep with optional offset.
  ! Calendar day 1.0 = 0Z on Jan 1.
  
 real :: get_curr_calday
  
 get_curr_calday = curr_dofyr

end function get_curr_calday

!=========================================================================================

function get_calday(ymd, tod)

! Return calendar day corresponding to specified time instant.
! Calendar day 1.0 = 0Z on Jan 1.

! fzeng: 
! combined info from 
! (1) subroutine get_dofyr_pentad in Catchment date_time_util.F90: the method 
! (2) subroutine ESMF_TimeGetDayOfYearInteger in CLM4.5 ESMF_TimeMod.F90: 
!     output day of the year ranges from 1 to 365
! (3) function get_calday and function TimeSetymd in CLM4.5 clm_time_manager.F90
! (4) function days_in_month in GEOSsurface_GridComp/Shared/Raster/src/leap_year.F90

! Arguments
   integer, intent(in) :: &
      ymd,   &! date in yearmmdd format
      tod     ! time of day (seconds past 0Z)

! Return value
   real :: get_calday

! Local variables   
   integer :: yr, mon, day          ! Year, month, day as integers
   integer :: i
   integer, dimension(12), parameter :: days_in_month_nonleap = &
         (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
   
   yr  = ymd / 10000
   mon = (ymd - yr*10000) / 100
   day =  ymd - yr*10000 - mon*100
   
   get_calday = day
   do i=1,mon-1 
      get_calday = get_calday + days_in_month_nonleap(i)
   end do
   
   if ( (get_calday > 366.0) .and. (get_calday <= 367.0) )then
       get_calday = get_calday - 1.0
   end if

   if ( (get_calday < 1.0) .or. (get_calday > 366.0) )then
      write(iulog,*) 'clm::get_calday = ', get_calday
      stop 'clm::get_calday: error calday out of range'
   end if

end function get_calday

!=========================================================================================

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

!=========================================================================================

function is_end_curr_day( )

  ! Return true if current timestep is last timestep in current day.

  ! Return value
 logical :: is_end_curr_day

  ! Local variables
 integer ::&
    yr,    &! year
    mon,   &! month
    day,   &! day of month
    tod     ! time of day (seconds past 0Z)

 call get_curr_date(yr, mon, day, tod)
 is_end_curr_day = (tod == 0)

end function is_end_curr_day

!=========================================================================================

logical function is_restart( )

  ! Determine if it's a restart run

 is_restart = .false.

end function is_restart

end module clm_time_manager
