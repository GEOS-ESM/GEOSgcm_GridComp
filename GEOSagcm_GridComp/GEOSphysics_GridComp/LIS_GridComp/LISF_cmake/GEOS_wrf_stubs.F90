! GEOS-owned replacement for the wrf_message/wrf_error_fatal stubs that
! LISF's noahmp36_wrf_routines.F90 only provides when COUPLED is undefined.
! LISF_cmake/CMakeLists.txt defines COUPLED (to exclude lis/offline's
! conflicting `program` main), but NoahMP.3.6/4.0.1 still call these
! routines unconditionally, so they must be provided from here instead.
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
