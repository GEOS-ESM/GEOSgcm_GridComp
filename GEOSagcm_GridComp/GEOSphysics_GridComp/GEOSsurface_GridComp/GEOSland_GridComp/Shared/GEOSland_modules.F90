#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

MODULE GEOSland_modules

  use ESMF
  use MAPL
  
  implicit none

  private

  type, public :: MODISReader
   contains
     procedure, public :: modis_date
     procedure, public :: read_modis_data
  end type MODISReader

  type, public :: TimeUtil
   contains
     generic,   public  :: TimeFrac => TimeFrac_esmf, TimeFrac_yyyymmdd
     procedure, public  :: MidTime
     procedure, private :: TimeFrac_esmf
     procedure, private :: TimeFrac_yyyymmdd
  end type TimeUtil

contains
  
  ! ---------------------------------------------------------------------------
    
  integer function modis_date (this,DOY, interval) result (MOD_DOY)
    
    implicit none
    class (MODISReader), intent(inout)         :: this
    integer, intent(in) :: DOY, interval
    integer, parameter  :: N_MODIS_DAYS8 = 46, N_MODIS_DAYS5 = 74
    integer, dimension (N_MODIS_DAYS8), target ::     &
         MODIS_DOYS8 = (/                             &
         1  ,  9, 17, 25, 33, 41, 49, 57, 65,         &
         73 , 81, 89, 97,105,113,121,129,137,         &
         145,153,161,169,177,185,193,201,209,         &
         217,225,233,241,249,257,265,273,281,         &
         289,297,305,313,321,329,337,345,353,361/)
    
    integer, dimension (N_MODIS_DAYS5), target ::     &
         MODIS_DOYS5 = (/                             &
         1  ,  6, 11, 16, 21, 26, 31, 36, 41, 46,     &
         51 , 56, 61, 66, 71, 76, 81, 86, 91, 96,     &
         101,106,111,116,121,126,131,136,141,146,     &
         151,156,161,166,171,176,181,186,191,196,     &
         201,206,211,216,221,226,231,236,241,246,     &
         251,256,261,266,271,276,281,286,291,296,     &
         301,306,311,316,321,326,331,336,341,346,     &
         351,356,361,366/)
    integer, dimension(:), pointer :: MODIS_DOYS
    integer :: i,N_MODIS_DATES
    
    select case (interval)
    case (8)
       MODIS_DOYS => MODIS_DOYS8
       N_MODIS_DATES = N_MODIS_DAYS8
    case (5)
       MODIS_DOYS => MODIS_DOYS5
       N_MODIS_DATES = N_MODIS_DAYS5
    end select
    
    if (DOY < MODIS_DOYS(N_MODIS_DATES)) then
       do i = 1, N_MODIS_DATES
          if (MODIS_DOYS(i) > DOY) exit
       end do
       MOD_DOY = MODIS_DOYS(i-1)
    else
       MOD_DOY = MODIS_DOYS(N_MODIS_DATES)
    endif
    
  end function modis_date
  
  ! ---------------------------------------------------------------------------
  
  subroutine read_modis_data (this,MAPL,CUR_YY, MOD_DOY, b4_modis_date, &
       GRIDNAME, MODIS_PATH, label, MODIS_DATA, MODIS_NIR) 

    implicit none
    class (MODISReader), intent(inout)       :: this
    type(MAPL_MetaComp), intent(in),pointer  :: MAPL
    character (*), intent (in)               :: GRIDNAME, MODIS_PATH, label
    integer, intent (in)                     :: CUR_YY, MOD_DOY
    logical, intent (in)                     :: b4_modis_date
    real, dimension(:), intent(out)          :: MODIS_DATA
    real, dimension(:), intent(out),optional :: MODIS_NIR
    type(ESMF_Grid)                          :: TILEGRID
    type(MAPL_LocStream)                     :: LOCSTREAM
    integer, pointer                         :: mask(:)
    integer                                  :: status, unit
    character*300                            :: filename
    CHARACTER(len=7)                         :: YYYYDoY
    logical                                  :: file_exists
    
    if(b4_modis_date) then
       WRITE (YYYYDoY,'(a4,i3.3)') 'YYYY',MOD_DOY
    else
       WRITE (YYYYDoY,'(i4.4,i3.3)') CUR_YY,MOD_DOY
    endif
    
    filename = trim(MODIS_PATH)//'/'//trim(GRIDNAME)//'/'//trim(label)//'_data.'//YYYYDoY
    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
       ASSERT_(.FALSE.)
    endif
    call MAPL_Get(MAPL, LocStream=LOCSTREAM, RC=STATUS)            ; VERIFY_(STATUS)
    call MAPL_LocStreamGet(LOCSTREAM, TILEGRID=TILEGRID, RC=STATUS); VERIFY_(STATUS)
    call MAPL_TileMaskGet(tilegrid,  mask, rc=status)              ; VERIFY_(STATUS)
    unit = GETFILE(trim(filename), form="unformatted", RC=STATUS)  ; VERIFY_(STATUS)
    call MAPL_VarRead(unit,tilegrid,MODIS_DATA,mask=mask,RC=STATUS); VERIFY_(STATUS)
    if(trim(label) == 'alb') &
         call MAPL_VarRead(unit,tilegrid,MODIS_NIR,mask=mask,RC=STATUS); VERIFY_(STATUS)
    call FREE_FILE(unit, RC=STATUS)                                ; VERIFY_(STATUS)
    
  end subroutine read_modis_data

  ! ---------------------------------------------------------------------------

  real function TimeFrac_ESMF (this, PTIME, CTIME, NTIME) result (tf)
    
    implicit none
    class (TimeUtil), intent(inout) :: this
    type(ESMF_Time),INTENT(IN) :: PTIME, CTIME, NTIME 
    real(ESMF_KIND_R8)         :: tdiff_21, tdiff_20, tfrac      
    type(ESMF_TimeInterval)    :: T21, T20
    integer                    :: status
    
    T21 = NTIME - PTIME
    T20 = NTIME - CTIME
    
    call ESMF_TimeIntervalGet(T21, s_r8=tdiff_21,rc=status) ; VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(T20,s_r8=tdiff_20,rc=status)  ; VERIFY_(STATUS)
    tfrac = tdiff_20/tdiff_21
    tf    = tfrac
    
  end function TimeFrac_ESMF

  ! ---------------------------------------------------------------------------
  
  real function TimeFrac_YYYYMMDD (this, dayp, dayc, dayn) result (tf)
    
    implicit none
    class (TimeUtil), intent(inout) :: this
    integer, intent (in)       :: dayp,dayc, dayn
    type(ESMF_Time)            :: PTIME, CTIME, NTIME 
    integer                    :: status, YY, MO, DA
    
    YY = dayp / 10000
    MO = (dayp - YY*10000) / 100
    DA = dayp - (YY*10000 + MO*100)
    call ESMF_TimeSet (ptime, yy=yy, mm=mo, dd=da, rc=status)     ; VERIFY_(STATUS)

    YY = dayc / 10000
    MO = (dayc - YY*10000) / 100
    DA = dayc - (YY*10000 + MO*100)
    call ESMF_TimeSet (ctime, yy=yy, mm=mo, dd=da, rc=status)     ; VERIFY_(STATUS)

    YY = dayn / 10000
    MO = (dayn - YY*10000) / 100
    DA = dayn - (YY*10000 + MO*100)
    call ESMF_TimeSet (ntime, yy=yy, mm=mo, dd=da, rc=status)     ; VERIFY_(STATUS)

    tf    = this%TimeFrac (PTIME, CTIME, NTIME)
    
  end function TimeFrac_YYYYMMDD

  ! ---------------------------------------------------------------------------

  function MidTime (this, DOY1, DOY2, CYEAR, start, last)  result (MT)

    implicit none
    class (TimeUtil), intent(inout) :: this
    integer, intent (in)    :: DOY1, DOY2, CYEAR
    type(ESMF_Time)         :: TIME1, MT, CYB
    type(ESMF_TimeInterval) :: DAY_SHIFT
    integer                 :: status, yyyy, ref_year
    logical, intent(in), optional :: start, last

    ref_year = 2002
    if (mod(cyear,4) /= 0) then 
       ref_year = 2002
    else if (mod(cyear,400) == 0) then
       ref_year = 2004
    else if (mod(cyear,100) == 0) then 
       ref_year = 2002
    else
       ref_year = 2004
    end if
    
    if (present(start)) then
       yyyy = ref_year -1
    else
       yyyy = ref_year
    endif
    
    call ESMF_TimeSet (CYB, yy=yyyy, mm=1, dd=1, rc=status)                  ; VERIFY_(STATUS)
    
    if(present(last)) then
       yyyy = ref_year + 1
       call ESMF_TimeSet (CYB, yy=yyyy, mm=1, dd=1, rc=status)               ; VERIFY_(STATUS)
    endif
    
    call ESMF_TimeIntervalSet (DAY_SHIFT, h=24*(DOY1-1), rc=status )         ; VERIFY_(STATUS)
    TIME1 = CYB + DAY_SHIFT
    
    if (DOY2 > DOY1) then
       call ESMF_TimeIntervalSet (DAY_SHIFT, h=24*(DOY2-DOY1)/2, rc=status ) ; VERIFY_(STATUS)
       MT = TIME1 + DAY_SHIFT
    else
       ! climatologies don't account for leap years)
       call ESMF_TimeIntervalSet (DAY_SHIFT, h=24*(366-DOY1)/2, rc=status )  ; VERIFY_(STATUS)
       MT = TIME1 + DAY_SHIFT
       if (DOY2 > 1) then
          call ESMF_TimeIntervalSet (DAY_SHIFT, h=24*(DOY2-1)/2, rc=status ) ; VERIFY_(STATUS)
          MT = MT + DAY_SHIFT
       endif
    endif
        
  end function MidTime
      
END MODULE GEOSland_modules


 
