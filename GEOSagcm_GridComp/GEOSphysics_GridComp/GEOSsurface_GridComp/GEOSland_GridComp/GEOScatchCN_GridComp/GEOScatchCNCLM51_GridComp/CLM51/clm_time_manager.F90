module clm_time_manager

#include "MAPL_Generic.h"
#include "shr_assert.h"

   use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
   use update_model_para4cn, only: curr_year,curr_month,curr_day,curr_dofyr,curr_hour,curr_min,curr_sec, &
                                   prev_year,prev_month,prev_day,prev_dofyr,prev_hour,prev_min,prev_sec
   use clm_varctl  , only: iulog
   use MAPL_ExceptionHandling
   use ESMF

   implicit none
   private

! Public methods
   
   public ::&
      get_step_size,            &! return step size in seconds
      get_step_size_real,       &! return step size in seconds, real-valued
      get_rad_step_size,        &! return radiation step size in seconds
      get_nstep,                &! return CN timestep number

      get_curr_date,            &! return date components at end of current timestep
      get_prev_date,            &! return date components at beginning of current timestep
!      get_curr_ESMF_Time,       &! get current time in terms of the ESMF_Time
!      get_start_date,           &! return components of the start date
!      get_driver_start_ymd,     &! return year/month/day (as integer in YYYYMMDD format) of driver start date
!      get_ref_date,             &! return components of the reference date
!      get_curr_time,            &! return components of elapsed time since reference date at end of current timestep
      get_curr_calday,          &! return calendar day at end of current timestep
      get_calday,               &! return calendar day from input date
!      get_calendar,             &! return calendar
      
      get_days_per_year,        &! return the days per year for current year
      get_local_timestep_time,  &! return the local time for the input longitude to the nearest time-step
      get_local_time,           &! return the local time for the input longitude
      get_curr_yearfrac,        &! return the fractional position in the current year, as of the end of the current timestep
      get_prev_yearfrac,        &! return the fractional position in the current year, as of the beginning of the current timestep


      is_end_curr_day,          &! return true on last timestep in current day
      is_beg_curr_year,         &! return true on first timestep in current year
      is_restart,               &! return true if this is a restart run
      is_first_step,            &! dummy function here, because it is loaded, but not used
      is_near_local_noon,       &! return true if near local noon
      update_rad_dtime          ! track radiation interval via nstep

   integer,  parameter :: uninit_int = -999999999
   integer, save ::&
        dtime          = -999999999,  &! timestep in seconds
        dtime_rad      = -999999999,  &! radiation interval in seconds
        nstep_rad_prev = -999999999    ! radiation interval in seconds

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

real(r8) function get_step_size_real()

    ! Return the step size in seconds, as a real value

    get_step_size_real = real(get_step_size(), r8)

  end function get_step_size_real

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

 get_nstep = get_nstep - 1  
end function get_nstep

!=========================================================================================
  subroutine update_rad_dtime(doalb)
    !---------------------------------------------------------------------------------
    ! called only on doalb timesteps to save off radiation nsteps
    ! 
    ! Local Arguments
    logical,intent(in) ::  doalb
    integer :: dtime,nstep

    if (doalb) then

       dtime=get_step_size()
       nstep = get_nstep()

       if (nstep_rad_prev == uninit_int ) then 
          dtime_rad = dtime
          nstep_rad_prev = nstep
       else 
          dtime_rad = (nstep - nstep_rad_prev) * dtime
          nstep_rad_prev = nstep
       endif
    end if
  end subroutine update_rad_dtime

  !=========================================================================================

integer function get_rad_step_size()

    character(len=*), parameter :: sub = 'clm::get_rad_step_size'

!    if ( .not. check_timemgr_initialized(sub) ) return

    if (nstep_rad_prev == uninit_int ) then 
       get_rad_step_size=get_step_size()
    else 
       get_rad_step_size=dtime_rad
    end if

end function get_rad_step_size

!=========================================================================================

  subroutine get_curr_date(yr, mon, day, tod, offset)

  ! Return date components valid at end of current timestep with an optional
  ! offset (positive or negative) in seconds.  

  implicit none
  integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

  integer, optional, intent(in) :: offset  ! Offset from current time in seconds (not used)

  yr  = curr_year
  mon = curr_month
  day = curr_day
  tod = 3600*curr_hour + 60*curr_min + curr_sec

  end subroutine get_curr_date
!=========================================================================================

  subroutine get_prev_date(yr, mon, day, tod)

    ! Return date components valid at beginning of current timestep.

  implicit none
  integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)
  !---------------------------------------------

  yr  = prev_year
  mon = prev_month
  day = prev_day
  tod = 3600*prev_hour + 60*prev_min + prev_sec

  end subroutine get_prev_date
!=========================================================================================

function get_curr_calday(offset)

  ! Return calendar day at end of current timestep with optional offset.
  ! Calendar day 1.0 = 0Z on Jan 1.
 integer, optional, intent(in) :: offset  ! Offset from current time in seconds.(not used)
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

logical function is_first_step(first)

  ! Return value
  logical, optional, intent(in) :: first  ! set to this time step

  logical, save :: is_first_default = .true.

 if ( present(first) ) then
   is_first_default = first
 end if

 is_first_step = is_first_default

end function is_first_step

!=========================================================================================

logical function is_restart( )

  ! Determine if it's a restart run

 is_restart = .false.

end function is_restart

!=========================================================================================

  logical function is_near_local_noon( londeg, deltasec )

    !---------------------------------------------------------------------------------
    ! Is this longitude near it's local noon?
    !
    ! uses
    use clm_varcon, only: degpsec, isecspday
    ! Arguments
    real(r8), intent(in) :: londeg   ! Longitude in degrees
    integer , intent(in) :: deltasec ! Number of seconds before or after local noon

    ! Local variables
    integer :: local_secs                         ! Local time in seconds
    integer, parameter :: noonsec = isecspday / 2 ! seconds at local noon
    !---------------------------------------------------------------------------------
    SHR_ASSERT( deltasec < noonsec, "deltasec must be less than 12 hours" )
    local_secs = get_local_timestep_time( londeg )

    if ( local_secs >= (noonsec - deltasec) .and. local_secs <= (noonsec + deltasec)) then
       is_near_local_noon = .true.
    else
       is_near_local_noon = .false.
    end if

    !---------------------------------------------------------------------------------
  end function is_near_local_noon

  !=========================================================================================

  integer function get_local_timestep_time( londeg, offset, rc )

    !---------------------------------------------------------------------------------
    ! Get the local time for this longitude that is evenly divisible by the time-step
    !
    ! uses
    use clm_varcon, only: degpsec, isecspday
    ! Arguments
    real(r8)         , intent(in) :: londeg  ! Longitude in degrees
    integer, optional, intent(in) :: offset  ! Offset from current time in seconds (either sign)
    integer, optional, intent(out) :: rc

    ! Local variables
    integer  :: yr, mon, day    ! year, month, day, unused
    integer  :: secs            ! seconds into the day
    real(r8) :: lon             ! positive longitude
    integer  :: offset_sec      ! offset seconds (either 0 for current time or -dtime for previous time)
    !---------------------------------------------------------------------------------
    if ( present(offset) ) then
       offset_sec = offset
       _ASSERT(.FALSE.,"offset function not enabled")
    else
       offset_sec = 0
    end if
    SHR_ASSERT( londeg >= -180.0_r8, "londeg must be greater than -180" )
    SHR_ASSERT( londeg <= 360.0_r8,  "londeg must be less than 360" )
    call  get_curr_date(yr, mon, day, secs )
    lon = londeg
    if ( lon < 0.0_r8 ) lon = lon + 360.0_r8
    get_local_timestep_time  = secs + nint((lon/degpsec)/real(dtime,r8))*dtime
    get_local_timestep_time  = mod(get_local_timestep_time,isecspday)
  end function get_local_timestep_time

  !=========================================================================================

!  function get_curr_ESMF_Time( )
!
!    ! Return the current time as ESMF_Time
!
!    type(ESMF_Time) :: get_curr_ESMF_Time
!    character(len=*), parameter :: sub = 'clm::get_curr_ESMF_Time'
!    integer :: rc, status
!
!   ! if ( .not. check_timemgr_initialized(sub) ) return
!
!    call ESMF_ClockGet( tm_clock, currTime=get_curr_ESMF_Time, rc=STATUS )
!    VERIFY_(STATUS)
!
!  end function get_curr_ESMF_Time

  integer function get_local_time( londeg, starttime, offset )

    !---------------------------------------------------------------------------------
    ! Get the local time for this longitude
    !
    ! uses
    use clm_varcon, only: degpsec, isecspday
    ! Arguments
    real(r8)         , intent(in) :: londeg       ! Longitude in degrees
    integer, optional, intent(in) :: starttime    ! Start time (sec)
    integer, optional, intent(in) :: offset       ! Offset from current time in seconds (either sign)

    ! Local variables
    integer  :: yr, mon, day    ! year, month, day, unused
    integer  :: secs            ! seconds into the day
    integer  :: start           ! start seconds
    integer  :: offset_sec      ! offset seconds (either 0 for current time or -dtime for previous time)
    real(r8) :: lon             ! positive longitude
    !---------------------------------------------------------------------------------
    if ( present(starttime) ) then
       start = starttime
    else
       start = 0
    end if
    if ( present(offset) ) then
       offset_sec = offset
    else
       offset_sec = 0
    end if
    SHR_ASSERT( start >= 0,            "starttime must be greater than or equal to zero" )
    SHR_ASSERT( start <= isecspday,    "starttime must be less than or equal to number of seconds in a day" )
    SHR_ASSERT( londeg >= -180.0_r8,   "londeg must be greater than -180" )
    SHR_ASSERT( londeg <= 360.0_r8,    "londeg must be less than 360" )
    SHR_ASSERT( (offset_sec == 0) .or. (offset_sec == -dtime), "offset must be zero or negative time-step" )
    call  get_curr_date(yr, mon, day, secs, offset=offset_sec )
    lon = londeg
    if ( lon < 0.0_r8 ) lon = lon + 360.0_r8
    get_local_time  = modulo(secs + nint(londeg/degpsec), isecspday)
    get_local_time  = modulo(get_local_time - start,isecspday)
  end function get_local_time

  !-----------------------------------------------------------------------
  logical function is_beg_curr_year()
    !
    ! !DESCRIPTION:
    ! Return true if current timestep is first timestep in current year.
    !
    ! !LOCAL VARIABLES:
    integer ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)

    character(len=*), parameter :: subname = 'is_beg_curr_year'
    !-----------------------------------------------------------------------

   ! if ( .not. check_timemgr_initialized(subname) ) return

    call get_curr_date(yr, mon, day, tod)
    is_beg_curr_year = (mon == 1 .and. day == 1 .and. tod == dtime)

  end function is_beg_curr_year

  !=========================================================================================

  function get_curr_yearfrac( offset )

    !---------------------------------------------------------------------------------
    ! Get the fractional position in the current year, as of the end of the current
    ! timestep. This is 0 at midnight on Jan 1, and 1 at the end of Dec 31.

    !
    ! Arguments
    real(r8) :: get_curr_yearfrac  ! function result

    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
    ! Positive for future times, negative 
    ! for previous times.

    character(len=*), parameter :: sub = 'clm::get_curr_yearfrac'
    real(r8) :: cday               ! current calendar day (1.0 = 0Z on Jan 1)
    real(r8) :: days_per_year      ! days per year

  !  if ( .not. check_timemgr_initialized(sub) ) return

    cday          = get_curr_calday(offset=offset)
    days_per_year = get_days_per_year()

    get_curr_yearfrac = (cday - 1._r8)/days_per_year

  end function get_curr_yearfrac

  !=========================================================================================

  function get_prev_yearfrac()

    !---------------------------------------------------------------------------------
    ! Get the fractional position in the current year, as of the beginning of the current
    ! timestep. This is 0 at midnight on Jan 1, and 1 at the end of Dec 31.

    !
    ! Arguments
    real(r8) :: get_prev_yearfrac  ! function result

    character(len=*), parameter :: sub = 'clm::get_curr_yearfrac'

  !  if ( .not. check_timemgr_initialized(sub) ) return

    get_prev_yearfrac = get_curr_yearfrac(offset = -dtime)

  end function get_prev_yearfrac

  !=========================================================================================

!  function get_curr_calday(offset)
!
!    ! Return calendar day at end of current timestep with optional offset.
!    ! Calendar day 1.0 = 0Z on Jan 1.
!
!    ! Arguments
!    integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
!    ! Positive for future times, negative 
!    ! for previous times.
!    ! Return value
!    real(r8) :: get_curr_calday
!
!    ! Local variables
!    character(len=*), parameter :: sub = 'clm::get_curr_calday'
!    integer :: rc
!    type(ESMF_Time) :: date
!    type(ESMF_TimeInterval) :: off, diurnal
!    integer :: year, month, day, tod
!    !-----------------------------------------------------------------------------------------
!
!!    if ( .not. check_timemgr_initialized(sub) ) return
!
!    call ESMF_ClockGet( tm_clock, currTime=date, rc=rc )
!    call chkrc(rc, sub//': error return from ESMF_ClockGet')
!
!    if (present(offset)) then
!       if (offset > 0) then
!          call ESMF_TimeIntervalSet( off, s=offset, rc=rc )
!          call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
!          date = date + off
!       else if (offset < 0) then
!          call ESMF_TimeIntervalSet( off, s=-offset, rc=rc )
!          call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
!          date = date - off
!       end if
!    end if
!
!    if ( tm_perp_calendar ) then
!       call ESMF_TimeGet(date, yy=year, mm=month, dd=day, s=tod, rc=rc)
!       call chkrc(rc, sub//': error return from ESMF_TimeGet')
!       call ESMF_TimeIntervalSet( diurnal, s=tod, rc=rc )
!       call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
!       date = tm_perp_date + diurnal
!    end if
!
!    call ESMF_TimeGet( date, dayOfYear_r8=get_curr_calday, rc=rc )
!    call chkrc(rc, sub//': error return from ESMF_TimeGet')
!    !----------------------------------------------------------------------------------------!
!    !!!!!!!!!!!!!! WARNING HACK TO ENABLE Gregorian CALENDAR WITH SHR_ORB !!!!!!!!!!!!!!!!!!!!
!    !!!! The following hack fakes day 366 by reusing day 365. This is just because the  !!!!!!
!    !!!! current shr_orb_decl calculation can't handle days > 366.                      !!!!!!
!    !!!!       Dani Bundy-Coleman and Erik Kluzek Aug/2008                              !!!!!!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    if ( (get_curr_calday > 366.0) .and. (get_curr_calday <= 367.0) .and. &
!         (trim(calendar) == GREGORIAN_C) )then
!       get_curr_calday = get_curr_calday - 1.0_r8
!    end if
!    !!!!!!!!!!!!!! END HACK TO ENABLE Gregorian CALENDAR WITH SHR_ORB !!!!!!!!!!!!!!!!!!!!!!!!
!    !----------------------------------------------------------------------------------------!
!    if ( (get_curr_calday < 1.0) .or. (get_curr_calday > 366.0) )then
!       write(iulog,*) sub, ' = ', get_curr_calday
!       if ( present(offset) ) write(iulog,*) 'offset = ', offset
!       call shr_sys_abort( sub//': error get_curr_calday out of bounds' )
!    end if
!
!  end function get_curr_calday
end module clm_time_manager
