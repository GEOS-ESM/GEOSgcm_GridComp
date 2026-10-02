! GEOS-owned stubs for symbols that LISF sources reference but that are
! never actually called with our LIS_misc.h feature-switch settings.
! Kept here instead of editing the vendored @LISF sources.

! LISF's noahmp36_wrf_routines.F90 only provides wrf_message/wrf_error_fatal
! when COUPLED is undefined. LISF_cmake/CMakeLists.txt defines COUPLED (to
! exclude lis/offline's conflicting `program` main), but NoahMP.3.6/4.0.1
! still call these routines unconditionally, so they must be provided here.
subroutine wrf_error_fatal (string)
   character (len=*) :: string
   print *,string
   stop
end subroutine wrf_error_fatal

subroutine wrf_message(message)
    implicit none
    character(len=*), intent(in) :: message
    write(0,*) trim(message)
end subroutine wrf_message

subroutine wrf_error_fatal3(file, line, message)
    implicit none
    character(len=*), intent(in) :: file
    integer,          intent(in) :: line
    character(len=*), intent(in) :: message
    write(0,*) trim(file), 'line: ', line, ': ', trim(message)
    stop (1)
end subroutine wrf_error_fatal3

! LIS_misc.h undefs USE_GRIBAPI, so retrieve_GLDAS1data/retrieve_NLDAS2data
! (in HYMAP3_router's runoff readers) never execute their call
! grib_get(igrib,"values",var(index,:),ios) at runtime -- but that call
! isn't itself guarded by USE_GRIBAPI, so the symbol is still referenced at
! link time. Stub it out to match that one call site; never invoked.
subroutine grib_get(igrib, key, values, status)
    implicit none
    integer,       intent(in)  :: igrib
    character(*),  intent(in)  :: key
    real,          intent(out) :: values(*)
    integer,       intent(out) :: status
    status = -1
end subroutine grib_get
