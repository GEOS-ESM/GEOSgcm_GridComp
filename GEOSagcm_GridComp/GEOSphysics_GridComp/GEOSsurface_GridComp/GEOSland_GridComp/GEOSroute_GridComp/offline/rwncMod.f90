module rwncfile

  use ncdio
  implicit none
  
  public :: read_ncfile_int1d
  public :: read_ncfile_real1d
  public :: read_ncfile_double1d
  
  public :: read_ncfile_int2d
  public :: read_ncfile_int3d  
  public :: read_ncfile_real2d
  public :: read_ncfile_real3d
  public :: read_ncfile_double2d
  public :: read_ncfile_double3d
  
  public :: write_ncfile_int2d
  public :: write_ncfile_real2d
  public :: write_ncfile_double2d
  
  public :: create_ncfile_byte2d
  public :: create_ncfile_short2d 
  public :: create_ncfile_short3d  
  public :: create_ncfile_int3d  
  public :: create_ncfile_int2d

  public :: create_ncfile_long2d
  public :: create_ncfile_real2d
  public :: create_ncfile_real3d
  public :: create_ncfile_double2d
  
  contains
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_int1d(filename,varname,var,n)
    character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in) :: n
  integer, intent(inout) :: var(n)
  
    character(len=4) :: subname="read"
  integer :: ncid, varid  

    call check_ret(nf_open(filename,0,ncid),subname)
    call check_ret(nf_inq_varid(ncid,varname,varid),subname)
  call check_ret(nf_get_var_int(ncid,varid,var),subname)
  call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_int1d    
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_real1d(filename,varname,var,n)
    character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in) :: n
	real, intent(inout) :: var(n)
	
    character(len=4) :: subname="read"
	integer :: ncid, varid	

    call check_ret(nf_open(filename,0,ncid),subname)
    call check_ret(nf_inq_varid(ncid,varname,varid),subname)
	call check_ret(nf_get_var_real(ncid,varid,var),subname)
	call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_real1d  
!------------------------------------------------------------------------------------------  
  subroutine read_ncfile_double1d(filename,varname,var,n)
    character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in) :: n
	real*8, intent(inout) :: var(n)
	
    character(len=4) :: subname="read"
	integer :: ncid, varid	

    call check_ret(nf_open(filename,0,ncid),subname)
    call check_ret(nf_inq_varid(ncid,varname,varid),subname)
	call check_ret(nf_get_var_double(ncid,varid,var),subname)
	call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_double1d
!------------------------------------------------------------------------------------------  
  subroutine read_ncfile_int2d(filename,varname,var,nlon,nlat)
    character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in) :: nlon, nlat
	integer, intent(inout) :: var(nlon,nlat)
	
    character(len=4) :: subname="read"
	integer :: ncid, varid	

    call check_ret(nf_open(filename,0,ncid),subname)
    call check_ret(nf_inq_varid(ncid,varname,varid),subname)
	call check_ret(nf_get_var_int(ncid,varid,var),subname)
	call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_int2d
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_int3d(filename,varname,var,nlon,nlat,nlev)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nlon, nlat, nlev
  integer, intent(inout) :: var(nlon,nlat,nlev)
  
  character(len=4) :: subname="read"
  integer :: ncid, varid  

  call check_ret(nf_open(filename,0,ncid),subname)
  call check_ret(nf_inq_varid(ncid,varname,varid),subname)
  call check_ret(nf_get_var_int(ncid,varid,var),subname)
  call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_int3d   
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_real2d(filename,varname,var,nlon,nlat)
    character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in) :: nlon, nlat
	real, intent(inout) :: var(nlon,nlat)
	
    character(len=4) :: subname="read"
	integer :: ncid, varid	

    call check_ret(nf_open(filename,0,ncid),subname)
    call check_ret(nf_inq_varid(ncid,varname,varid),subname)
	call check_ret(nf_get_var_real(ncid,varid,var),subname)
	call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_real2d
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_real3d(filename,varname,var,nlon,nlat,nlev)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nlon, nlat, nlev
  real, intent(inout) :: var(nlon,nlat,nlev)
  
  character(len=4) :: subname="read"
  integer :: ncid, varid  

  call check_ret(nf_open(filename,0,ncid),subname)
  call check_ret(nf_inq_varid(ncid,varname,varid),subname)
  call check_ret(nf_get_var_real(ncid,varid,var),subname)
  call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_real3d  
!------------------------------------------------------------------------------------------
  subroutine read_ncfile_double2d(filename,varname,var,nlon,nlat)
    character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in) :: nlon, nlat
	real*8, intent(inout) :: var(nlon,nlat)
	
    character(len=4) :: subname="read"
	integer :: ncid, varid	

    call check_ret(nf_open(filename,0,ncid),subname)
    call check_ret(nf_inq_varid(ncid,varname,varid),subname)
	call check_ret(nf_get_var_double(ncid,varid,var),subname)
	call check_ret(nf_close(ncid), subname)  
  
  end subroutine read_ncfile_double2d


  subroutine read_ncfile_double3d(filename,varname,var,nlon,nlat,nlev)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nlon, nlat, nlev
  real*8, intent(inout) :: var(nlon,nlat,nlev)

  character(len=4) :: subname="read"
  integer :: ncid, varid

  call check_ret(nf_open(filename,0,ncid),subname)
  call check_ret(nf_inq_varid(ncid,varname,varid),subname)
  call check_ret(nf_get_var_double(ncid,varid,var),subname)
  call check_ret(nf_close(ncid), subname)

  end subroutine read_ncfile_double3d
!------------------------------------------------------------------------------------------
  subroutine write_ncfile_int2d(filename,varname,var,nlon,nlat)
    character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in) :: nlon, nlat
	integer, intent(inout) :: var(nlon,nlat)
	
    character(len=4) :: subname="write"
	integer :: ncid, varid, omode		
	
    call check_ret(nf_open(filename, nf_write, ncid), subname)
    call check_ret(nf_set_fill(ncid, nf_nofill, omode), subname)
    call ncd_ioglobal(varname=varname, data=var, ncid=ncid, flag='write')
    call check_ret(nf_sync(ncid),  subname)	  
    call check_ret(nf_close(ncid), subname) 
  end subroutine write_ncfile_int2d
!------------------------------------------------------------------------------------------
  subroutine write_ncfile_real2d(filename,varname,var,nlon,nlat)
    character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in) :: nlon, nlat
	real, intent(inout) :: var(nlon,nlat)
	
    character(len=4) :: subname="write"
	integer :: ncid, varid, omode		
	
    call check_ret(nf_open(filename, nf_write, ncid), subname)
    call check_ret(nf_set_fill(ncid, nf_nofill, omode), subname)
    call ncd_ioglobal(varname=varname, data=var, ncid=ncid, flag='write')
    call check_ret(nf_sync(ncid),  subname)	  
    call check_ret(nf_close(ncid), subname) 
  end subroutine write_ncfile_real2d
!------------------------------------------------------------------------------------------
  subroutine write_ncfile_double2d(filename,varname,var,nlon,nlat)
    character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in) :: nlon, nlat
	real*8, intent(inout) :: var(nlon,nlat)
	
    character(len=4) :: subname="write"
	integer :: ncid, varid, omode		
	
    call check_ret(nf_open(filename, nf_write, ncid), subname)
    call check_ret(nf_set_fill(ncid, nf_nofill, omode), subname)
    call ncd_ioglobal(varname=varname, data=var, ncid=ncid, flag='write')
    call check_ret(nf_sync(ncid),  subname)	  
    call check_ret(nf_close(ncid), subname) 
  end subroutine write_ncfile_double2d
!------------------------------------------------------------------------------------------
  subroutine create_ncfile_int2d(filename,varname,var,lon,lat,nlon,nlat)
    character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in) :: nlon, nlat
	integer, intent(inout) :: var(nlon,nlat)
	real*8, intent(in) :: lon(nlon),lat(nlat)
	
    character(len=4) :: subname="create"
	integer :: ncid, varid, dimid	
	real*8 :: lon1(nlon), lat1(nlat)
	  
	  lon1=lon
	  lat1=lat
      call check_ret(nf_create(trim(filename), nf_clobber, ncid), subname)
	  call check_ret(nf_def_dim(ncid,'lon',nlon, dimid), subname)
	  call check_ret(nf_def_dim(ncid,'lat',nlat, dimid), subname)
      call ncd_defvar(ncid=ncid, varname='lon', xtype=nf_double, dim1name='lon', &
                  long_name='longtitude', units='degrees_east')
      call ncd_defvar(ncid=ncid, varname='lat', xtype=nf_double, dim2name='lat', &
                  long_name='latitude', units='degrees_north')
      call ncd_defvar(ncid=ncid, varname=varname, xtype=nf_int, dim1name='lon', &
                  dim2name='lat', long_name=varname, units='unitless', fill_value=-9999.)						  
      call check_ret(nf_enddef(ncid), subname)
      call ncd_ioglobal(varname='lon', data=lon1, flag='write',ncid=ncid)
      call ncd_ioglobal(varname='lat', data=lat1, flag='write',ncid=ncid)
	  call ncd_ioglobal(varname=varname, data=var, flag='write',ncid=ncid)
	  call check_ret(nf_close(ncid), subname)	
  end subroutine create_ncfile_int2d

  subroutine create_ncfile_long2d(filename,varname,var,lon,lat,nlon,nlat)
    character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: varname
        integer, intent(in) :: nlon, nlat
        integer*8, intent(inout) :: var(nlon,nlat)
        real*8, intent(in) :: lon(nlon),lat(nlat)
        
    character(len=4) :: subname="create"
        integer :: ncid, varid, dimid   
        real*8 :: lon1(nlon), lat1(nlat)
        
          lon1=lon
          lat1=lat
      call check_ret(nf_create(trim(filename), NF_NETCDF4, ncid), subname)
          call check_ret(nf_def_dim(ncid,'lon',nlon, dimid), subname)
          call check_ret(nf_def_dim(ncid,'lat',nlat, dimid), subname)
      call ncd_defvar(ncid=ncid, varname='lon', xtype=nf_double, dim1name='lon',&
                  long_name='longtitude', units='degrees_east')
      call ncd_defvar(ncid=ncid, varname='lat', xtype=nf_double, dim2name='lat',&
                  long_name='latitude', units='degrees_north')
      call ncd_defvar(ncid=ncid, varname=varname, xtype=nf_int64, dim1name='lon',&
                  dim2name='lat', long_name=varname, units='unitless',fill_value=-9999.)                        
      call check_ret(nf_enddef(ncid), subname)
      call ncd_ioglobal(varname='lon', data=lon1, flag='write',ncid=ncid)
      call ncd_ioglobal(varname='lat', data=lat1, flag='write',ncid=ncid)
          call ncd_ioglobal(varname=varname, data=var, flag='write',ncid=ncid)
          call check_ret(nf_close(ncid), subname)       
  end subroutine create_ncfile_long2d

!------------------------------------------------------------------------------------------
  subroutine create_ncfile_byte2d(filename,varname,var,lon,lat,nlon,nlat)
    character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nlon, nlat
  byte, intent(inout) :: var(nlon,nlat)
  real*8, intent(in) :: lon(nlon),lat(nlat)
  
    character(len=4) :: subname="create"
  integer :: ncid, varid, dimid 
  real*8 :: lon1(nlon), lat1(nlat)
    
    lon1=lon
    lat1=lat
      call check_ret(nf_create(trim(filename), nf_clobber, ncid), subname)
    call check_ret(nf_def_dim(ncid,'lon',nlon, dimid), subname)
    call check_ret(nf_def_dim(ncid,'lat',nlat, dimid), subname)
      call ncd_defvar(ncid=ncid, varname='lon', xtype=nf_double, dim1name='lon', &
                  long_name='longtitude', units='degrees_east')
      call ncd_defvar(ncid=ncid, varname='lat', xtype=nf_double, dim2name='lat', &
                  long_name='latitude', units='degrees_north')
      call ncd_defvar(ncid=ncid, varname=varname, xtype=nf_byte, dim1name='lon', &
                  dim2name='lat', long_name=varname, units='unitless',fill_value=-128. )             
      call check_ret(nf_enddef(ncid), subname)
      call ncd_ioglobal(varname='lon', data=lon1, flag='write',ncid=ncid)
      call ncd_ioglobal(varname='lat', data=lat1, flag='write',ncid=ncid)
    call ncd_ioglobal(varname=varname, data=var, flag='write',ncid=ncid)
    call check_ret(nf_close(ncid), subname) 
  end subroutine create_ncfile_byte2d 

!------------------------------------------------------------------------------------------
  subroutine create_ncfile_short2d(filename,varname,var,lon,lat,nlon,nlat)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nlon, nlat
  integer*2, intent(inout) :: var(nlon,nlat)
  real*8, intent(in) :: lon(nlon),lat(nlat)
  
  character(len=4) :: subname="create"
  integer :: ncid, varid, dimid 
  real*8 :: lon1(nlon), lat1(nlat)
    
    lon1=lon
    lat1=lat
    call check_ret(nf_create(trim(filename), nf_clobber, ncid), subname)
    call check_ret(nf_def_dim(ncid,'lon',nlon, dimid), subname)
    call check_ret(nf_def_dim(ncid,'lat',nlat, dimid), subname)
      call ncd_defvar(ncid=ncid, varname='lon', xtype=nf_double, dim1name='lon', &
                  long_name='longtitude', units='degrees_east')
      call ncd_defvar(ncid=ncid, varname='lat', xtype=nf_double, dim2name='lat', &
                  long_name='latitude', units='degrees_north')
      call ncd_defvar(ncid=ncid, varname=varname, xtype=nf_short, dim1name='lon', &
                  dim2name='lat', long_name=varname, units='unitless',fill_value=-9999. )             
      call check_ret(nf_enddef(ncid), subname)
      call ncd_ioglobal(varname='lon', data=lon1, flag='write',ncid=ncid)
      call ncd_ioglobal(varname='lat', data=lat1, flag='write',ncid=ncid)
    call ncd_ioglobal(varname=varname, data=var, flag='write',ncid=ncid)
    call check_ret(nf_close(ncid), subname) 
  end subroutine create_ncfile_short2d  


!------------------------------------------------------------------------------------------
  subroutine create_ncfile_real2d(filename,varname,var,lon,lat,nlon,nlat)
    character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in) :: nlon, nlat
	real, intent(inout) :: var(nlon,nlat)
	real*8, intent(in) :: lon(nlon),lat(nlat)
	
    character(len=4) :: subname="create"
	integer :: ncid, varid, dimid	
	real*8 :: lon1(nlon), lat1(nlat)
	  
	  lon1=lon
	  lat1=lat	
      call check_ret(nf_create(trim(filename), nf_clobber, ncid), subname)
	  call check_ret(nf_def_dim(ncid,'lon',nlon, dimid), subname)
	  call check_ret(nf_def_dim(ncid,'lat',nlat, dimid), subname)
      call ncd_defvar(ncid=ncid, varname='lon', xtype=nf_double, dim1name='lon', &
                  long_name='longtitude', units='degrees_east')
      call ncd_defvar(ncid=ncid, varname='lat', xtype=nf_double, dim2name='lat', &
                  long_name='latitude', units='degrees_north')
      call ncd_defvar(ncid=ncid, varname=varname, xtype=nf_float, dim1name='lon', &
                  dim2name='lat', long_name=varname, units='unitless', fill_value=-9999.)						  
      call check_ret(nf_enddef(ncid), subname)
      call ncd_ioglobal(varname='lon', data=lon1, flag='write',ncid=ncid)
      call ncd_ioglobal(varname='lat', data=lat1, flag='write',ncid=ncid)
	  call ncd_ioglobal(varname=varname, data=var, flag='write',ncid=ncid)
	  call check_ret(nf_close(ncid), subname)	
  end subroutine create_ncfile_real2d

!------------------------------------------------------------------------------------------
  subroutine create_ncfile_short3d(filename,varname,var,lon,lat,lev,nlon,nlat,nlev)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nlon, nlat, nlev
  integer*2, intent(inout) :: var(nlon,nlat,nlev)
  real*8, intent(in) :: lon(nlon),lat(nlat),lev(nlev)
  
    character(len=4) :: subname="create"
    integer :: ncid, varid, dimid 
    real*8 :: lon1(nlon), lat1(nlat), lev1(nlev)
    
      lon1=lon
      lat1=lat  
      lev1=lev
      call check_ret(nf_create(trim(filename), nf_clobber, ncid), subname)
      call check_ret(nf_def_dim(ncid,'lon',nlon, dimid), subname)
      call check_ret(nf_def_dim(ncid,'lat',nlat, dimid), subname)
      call check_ret(nf_def_dim(ncid,'lev',nlev, dimid), subname)

      call ncd_defvar(ncid=ncid, varname='lon', xtype=nf_double, dim1name='lon', &
                  long_name='longtitude', units='degrees_east')
      call ncd_defvar(ncid=ncid, varname='lat', xtype=nf_double, dim2name='lat', &
                  long_name='latitude', units='degrees_north')
      call ncd_defvar(ncid=ncid, varname='lev', xtype=nf_double, dim2name='lev', &
                  long_name='level', units='unitless')

      call ncd_defvar(ncid=ncid, varname=varname, xtype=nf_short, dim1name='lon', &
                  dim2name='lat', dim3name='lev', long_name=varname, units='unitless', fill_value=-9999.)             
      call check_ret(nf_enddef(ncid), subname)
      call ncd_ioglobal(varname='lon', data=lon1, flag='write',ncid=ncid)
      call ncd_ioglobal(varname='lat', data=lat1, flag='write',ncid=ncid)
      call ncd_ioglobal(varname='lev', data=lev1, flag='write',ncid=ncid)      
      call ncd_ioglobal(varname=varname, data=var, flag='write',ncid=ncid)
      call check_ret(nf_close(ncid), subname) 
  end subroutine create_ncfile_short3d
!------------------------------------------------------------------------------------------
  subroutine create_ncfile_int3d(filename,varname,var,lon,lat,lev,nlon,nlat,nlev)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nlon, nlat, nlev
  integer, intent(inout) :: var(nlon,nlat,nlev)
  real*8, intent(in) :: lon(nlon),lat(nlat),lev(nlev)
  
    character(len=4) :: subname="create"
    integer :: ncid, varid, dimid 
    real*8 :: lon1(nlon), lat1(nlat), lev1(nlev)
    
      lon1=lon
      lat1=lat  
      lev1=lev
      call check_ret(nf_create(trim(filename), nf_clobber, ncid), subname)
      call check_ret(nf_def_dim(ncid,'lon',nlon, dimid), subname)
      call check_ret(nf_def_dim(ncid,'lat',nlat, dimid), subname)
      call check_ret(nf_def_dim(ncid,'lev',nlev, dimid), subname)

      call ncd_defvar(ncid=ncid, varname='lon', xtype=nf_double, dim1name='lon', &
                  long_name='longtitude', units='degrees_east')
      call ncd_defvar(ncid=ncid, varname='lat', xtype=nf_double, dim2name='lat', &
                  long_name='latitude', units='degrees_north')
      call ncd_defvar(ncid=ncid, varname='lev', xtype=nf_double, dim2name='lev', &
                  long_name='level', units='unitless')

      call ncd_defvar(ncid=ncid, varname=varname, xtype=nf_int, dim1name='lon', &
                  dim2name='lat', dim3name='lev', long_name=varname, units='unitless', fill_value=-9999.)             
      call check_ret(nf_enddef(ncid), subname)
      call ncd_ioglobal(varname='lon', data=lon1, flag='write',ncid=ncid)
      call ncd_ioglobal(varname='lat', data=lat1, flag='write',ncid=ncid)
      call ncd_ioglobal(varname='lev', data=lev1, flag='write',ncid=ncid)      
      call ncd_ioglobal(varname=varname, data=var, flag='write',ncid=ncid)
      call check_ret(nf_close(ncid), subname) 
  end subroutine create_ncfile_int3d
!------------------------------------------------------------------------------------------
  subroutine create_ncfile_real3d(filename,varname,var,lon,lat,lev,nlon,nlat,nlev)
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nlon, nlat, nlev
  real, intent(inout) :: var(nlon,nlat,nlev)
  real*8, intent(in) :: lon(nlon),lat(nlat),lev(nlev)
  
    character(len=4) :: subname="create"
    integer :: ncid, varid, dimid 
    real*8 :: lon1(nlon), lat1(nlat), lev1(nlev)
    
      lon1=lon
      lat1=lat  
      lev1=lev
      call check_ret(nf_create(trim(filename), nf_clobber, ncid), subname)
      call check_ret(nf_def_dim(ncid,'lon',nlon, dimid), subname)
      call check_ret(nf_def_dim(ncid,'lat',nlat, dimid), subname)
      call check_ret(nf_def_dim(ncid,'lev',nlev, dimid), subname)

      call ncd_defvar(ncid=ncid, varname='lon', xtype=nf_double, dim1name='lon', &
                  long_name='longtitude', units='degrees_east')
      call ncd_defvar(ncid=ncid, varname='lat', xtype=nf_double, dim2name='lat', &
                  long_name='latitude', units='degrees_north')
      call ncd_defvar(ncid=ncid, varname='lev', xtype=nf_double, dim2name='lev', &
                  long_name='level', units='unitless')

      call ncd_defvar(ncid=ncid, varname=varname, xtype=nf_float, dim1name='lon', &
                  dim2name='lat', dim3name='lev', long_name=varname, units='unitless', fill_value=-9999.)             
      call check_ret(nf_enddef(ncid), subname)
      call ncd_ioglobal(varname='lon', data=lon1, flag='write',ncid=ncid)
      call ncd_ioglobal(varname='lat', data=lat1, flag='write',ncid=ncid)
      call ncd_ioglobal(varname='lev', data=lev1, flag='write',ncid=ncid)      
      call ncd_ioglobal(varname=varname, data=var, flag='write',ncid=ncid)
      call check_ret(nf_close(ncid), subname) 
  end subroutine create_ncfile_real3d

!------------------------------------------------------------------------------------------
  subroutine create_ncfile_double2d(filename,varname,var,lon,lat,nlon,nlat)
    character(len=*), intent(in) :: filename
	character(len=*), intent(in) :: varname
	integer, intent(in) :: nlon, nlat
	real*8, intent(inout) :: var(nlon,nlat)
	real*8, intent(in) :: lon(nlon),lat(nlat)
	
    character(len=4) :: subname="create"
	integer :: ncid, varid, dimid	
	real*8 :: lon1(nlon), lat1(nlat)
	  
	  lon1=lon
	  lat1=lat	
      call check_ret(nf_create(trim(filename), nf_clobber, ncid), subname)
	  call check_ret(nf_def_dim(ncid,'lon',nlon, dimid), subname)
	  call check_ret(nf_def_dim(ncid,'lat',nlat, dimid), subname)
      call ncd_defvar(ncid=ncid, varname='lon', xtype=nf_double, dim1name='lon', &
                  long_name='longtitude', units='degrees_east')
      call ncd_defvar(ncid=ncid, varname='lat', xtype=nf_double, dim2name='lat', &
                  long_name='latitude', units='degrees_north')
      call ncd_defvar(ncid=ncid, varname=varname, xtype=nf_double, dim1name='lon', &
                  dim2name='lat', long_name=varname, units='unitless', fill_value=-9999.)						  
      call check_ret(nf_enddef(ncid), subname)
      call ncd_ioglobal(varname='lon', data=lon1, flag='write',ncid=ncid)
      call ncd_ioglobal(varname='lat', data=lat1, flag='write',ncid=ncid)
	  call ncd_ioglobal(varname=varname, data=var, flag='write',ncid=ncid)
	  call check_ret(nf_close(ncid), subname)	
  end subroutine create_ncfile_double2d
!------------------------------------------------------------------------------------------
end module rwncfile
	
