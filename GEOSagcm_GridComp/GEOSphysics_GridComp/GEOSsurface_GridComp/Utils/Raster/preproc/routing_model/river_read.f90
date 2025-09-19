module river_read
!module for reading river routing-related netcdf data

  implicit none
  include 'netcdf.inc'
  
  public :: read_ncfile_int1d
  public :: read_ncfile_real1d
  public :: read_ncfile_double1d
  
  public :: read_ncfile_int2d
  public :: read_ncfile_int3d  
  public :: read_ncfile_real2d
  public :: read_ncfile_real3d
  public :: read_ncfile_double2d
  public :: read_ncfile_double3d
  
  contains
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_int1d(filename,varname,var,n)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in)          :: n
  integer, intent(inout)       :: var(n)
  
  character(len=4)             :: subname="read"
  integer                      :: ncid, varid  

  call check_ret(nf_open(filename,0,ncid),subname)
  call check_ret(nf_inq_varid(ncid,varname,varid),subname)
  call check_ret(nf_get_var_int(ncid,varid,var),subname)
  call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_int1d    
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_real1d(filename,varname,var,n)
  character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in)          :: n
	real, intent(inout)          :: var(n)
	
  character(len=4)             :: subname="read"
	integer                      :: ncid, varid	

  call check_ret(nf_open(filename,0,ncid),subname)
  call check_ret(nf_inq_varid(ncid,varname,varid),subname)
	call check_ret(nf_get_var_real(ncid,varid,var),subname)
	call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_real1d  
!------------------------------------------------------------------------------------------  
  subroutine read_ncfile_double1d(filename,varname,var,n)
  character(len=*), intent(in)   :: filename
	character(len=*), intent(in)   :: varname
	integer, intent(in)            :: n
	real*8, intent(inout)          :: var(n)
	
  character(len=4)               :: subname="read"
	integer                        :: ncid, varid	

    call check_ret(nf_open(filename,0,ncid),subname)
    call check_ret(nf_inq_varid(ncid,varname,varid),subname)
	call check_ret(nf_get_var_double(ncid,varid,var),subname)
	call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_double1d
!------------------------------------------------------------------------------------------  
  subroutine read_ncfile_int2d(filename,varname,var,nlon,nlat)
  character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in)          :: nlon, nlat
	integer, intent(inout)       :: var(nlon,nlat)
	
  character(len=4)             :: subname="read"
	integer                      :: ncid, varid	

    call check_ret(nf_open(filename,0,ncid),subname)
    call check_ret(nf_inq_varid(ncid,varname,varid),subname)
	call check_ret(nf_get_var_int(ncid,varid,var),subname)
	call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_int2d
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_int3d(filename,varname,var,nlon,nlat,nlev)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in)          :: nlon, nlat, nlev
  integer, intent(inout)       :: var(nlon,nlat,nlev)
  
  character(len=4)             :: subname="read"
  integer                      :: ncid, varid  

  call check_ret(nf_open(filename,0,ncid),subname)
  call check_ret(nf_inq_varid(ncid,varname,varid),subname)
  call check_ret(nf_get_var_int(ncid,varid,var),subname)
  call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_int3d   
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_real2d(filename,varname,var,nlon,nlat)
  character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in)          :: nlon, nlat
	real, intent(inout)          :: var(nlon,nlat)
	
  character(len=4)             :: subname="read"
	integer                      :: ncid, varid	

  call check_ret(nf_open(filename,0,ncid),subname)
  call check_ret(nf_inq_varid(ncid,varname,varid),subname)
	call check_ret(nf_get_var_real(ncid,varid,var),subname)
	call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_real2d
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_real3d(filename,varname,var,nlon,nlat,nlev)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in)          :: nlon, nlat, nlev
  real, intent(inout)          :: var(nlon,nlat,nlev)
  
  character(len=4)             :: subname="read"
  integer                      :: ncid, varid  

  call check_ret(nf_open(filename,0,ncid),subname)
  call check_ret(nf_inq_varid(ncid,varname,varid),subname)
  call check_ret(nf_get_var_real(ncid,varid,var),subname)
  call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_real3d  
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_double2d(filename,varname,var,nlon,nlat)
  character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in)          :: nlon, nlat
	real*8, intent(inout)        :: var(nlon,nlat)
	
  character(len=4)             :: subname="read"
	integer                      :: ncid, varid	

  call check_ret(nf_open(filename,0,ncid),subname)
  call check_ret(nf_inq_varid(ncid,varname,varid),subname)
	call check_ret(nf_get_var_double(ncid,varid,var),subname)
	call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_double2d
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_double3d(filename,varname,var,nlon,nlat,nlev)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in)          :: nlon, nlat, nlev
  real*8, intent(inout)        :: var(nlon,nlat,nlev)

  character(len=4)             :: subname="read"
  integer                      :: ncid, varid

  call check_ret(nf_open(filename,0,ncid),subname)
  call check_ret(nf_inq_varid(ncid,varname,varid),subname)
  call check_ret(nf_get_var_double(ncid,varid,var),subname)
  call check_ret(nf_close(ncid), subname)

  end subroutine read_ncfile_double3d
!------------------------------------------------------------------------------------------
  subroutine check_ret(ret, calling)
  integer, intent(in) :: ret
  character(len=*)    :: calling

  if (ret /= NF_NOERR) then 
    write(6,*)'netcdf error from ',trim(calling)
    call endrun(nf_strerror(ret))
  end if
  end subroutine check_ret
!-----------------------------------------------------------------------
  subroutine endrun(msg,subname)
  character(len=*), intent(in), optional :: msg    
  character(len=*), intent(in), optional :: subname   

  if (present (subname)) then 
    write(6,*) 'ERROR in subroutine :', trim(subname)
  end if

  if (present (msg)) then
      write(6,*)'ENDRUN:', msg
  else
      write(6,*) 'ENDRUN: called without a message string'
  end if

  stop 
  end subroutine endrun  
!-----------------------------------------------------------------------

end module river_read
	
