
module ncdio
 use netcdf
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdioMod
!
! !DESCRIPTION:
! Generic interfaces to write fields to netcdf files
!
! !USES:
!
! !PUBLIC TYPES:
  implicit none
  include 'netcdf.inc' !netcdf库文件
  save
  public :: check_ret   ! checks return status of netcdf calls
  public :: check_var   ! determine if variable is on netcdf file
  public :: check_dim   ! validity check on dimension
  public :: ncd_defvar
!
! !REVISION HISTORY:
!
!EOP
!
! !PRIVATE METHODS:
!
  interface ncd_iolocal
     module procedure ncd_iolocal_int_1d
     module procedure ncd_iolocal_real_1d
     module procedure ncd_iolocal_double_1d	 
     module procedure ncd_iolocal_int_2d
     module procedure ncd_iolocal_real_2d
     module procedure ncd_iolocal_double_2d		 
  end interface

  interface ncd_ioglobal
     module procedure ncd_ioglobal_int_var
     module procedure ncd_ioglobal_real_var
     module procedure ncd_ioglobal_double_var	 
     module procedure ncd_ioglobal_int_1d
     module procedure ncd_ioglobal_real_1d
     module procedure ncd_ioglobal_double_1d	 
     module procedure ncd_ioglobal_byte_2d 
     module procedure ncd_ioglobal_short_2d  
     module procedure ncd_ioglobal_int_2d
     module procedure ncd_ioglobal_long_2d
     module procedure ncd_ioglobal_real_2d
     module procedure ncd_ioglobal_double_2d	 
     module procedure ncd_ioglobal_int_3d
     module procedure ncd_ioglobal_short_3d      
     module procedure ncd_ioglobal_real_3d
     module procedure ncd_ioglobal_double_3d	 
  end interface

  private :: endrun
  logical, public, parameter :: nc_masterproc = .true. ! proc 0 logical for printing msgs

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dim
!
! !INTERFACE:
  subroutine check_dim(ncid, dimname, value)
!
! !DESCRIPTION:
! Validity check on dimension
! 判断nc文件中指定维数dimname的长度与指定值value相等
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer, intent(in) :: value
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: dimid, dimlen    ! temporaries
!-----------------------------------------------------------------------

    call check_ret(nf_inq_dimid (ncid, trim(dimname), dimid), 'check_dim') !查询维数的代码
    call check_ret(nf_inq_dimlen (ncid, dimid, dimlen), 'check_dim') !查询维数的大小
    if (dimlen /= value) then
       write (6,*) 'CHECK_DIM error: mismatch of input dimension ',dimlen, &
            ' with expected value ',value,' for variable ',trim(dimname)
       call endrun()
    end if

  end subroutine check_dim

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_var
!
! !INTERFACE:
  subroutine check_var(ncid, varname, varid, readvar)
! 判断NC文件中是否含有名为varname的变量，如有则返回readvar=true且返回变量号varid，否则报错。
! !DESCRIPTION:
! Check if variable is on netcdf file
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)          :: ncid
    character(len=*), intent(in) :: varname
    integer, intent(out)         :: varid
    logical, intent(out)         :: readvar
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ret     ! return value
!-----------------------------------------------------------------------

    readvar = .true.
    if (nc_masterproc) then
       ret = nf_inq_varid (ncid, varname, varid)
       if (ret/=NF_NOERR) then
          write(6,*)'CHECK_VAR: variable ',trim(varname),' is not on initial dataset'
          readvar = .false.
       end if
    end if
  end subroutine check_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_ret
!
! !INTERFACE:
  subroutine check_ret(ret, calling)
! 返回NC文件操作是否正确
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ret
    character(len=*) :: calling
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

    if (ret /= NF_NOERR) then !如果打开nc文件出错，则提示出错信息
       write(6,*)'netcdf error from ',trim(calling)
       call endrun(nf_strerror(ret))
    end if

  end subroutine check_ret

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_defvar
!
! !INTERFACE:
  subroutine ncd_defvar(ncid, varname, xtype, &
       dim1name, dim2name, dim3name, dim4name, dim5name, &
       long_name, units, cell_method, missing_value, fill_value, &
       imissing_value, ifill_value)
! 定义NC变量，
! ncid--NC文件号
! varname--变量名称
! xtype--变量类型
! dim1name--第一维的名称
! dim2name--第二维的名称
! dim3name--第三维的名称
! dim4name--第四维的名称
! dim5name--第五维的名称
! long_name--属性-变量的完整名称
! units--属性-变量的单位
! cell_method--属性-值的来源说明
! missing_value--属性-实型缺测值
! fill_value--属性-实型的缺省值
! imissing_value--属性-整型的缺测值
! ifill_value--属性-整型的缺省值
! !DESCRIPTION:
!  Define a netcdf variable
!
! !ARGUMENTS:
    implicit none
    integer         , intent(in)  :: ncid                    ! input unit
    character(len=*), intent(in)  :: varname                 ! variable name
    integer         , intent(in)  :: xtype                   ! external type
    character(len=*), intent(in), optional :: dim1name       ! dimension name
    character(len=*), intent(in), optional :: dim2name       ! dimension name
    character(len=*), intent(in), optional :: dim3name       ! dimension name
    character(len=*), intent(in), optional :: dim4name       ! dimension name
    character(len=*), intent(in), optional :: dim5name       ! dimension name
    character(len=*), intent(in), optional :: long_name      ! attribute
    character(len=*), intent(in), optional :: units          ! attribute
    character(len=*), intent(in), optional :: cell_method    ! attribute
    real     , intent(in), optional :: missing_value  ! attribute for real
    real     , intent(in), optional :: fill_value     ! attribute for real
    integer         , intent(in), optional :: imissing_value ! attribute for int
    integer         , intent(in), optional :: ifill_value    ! attribute for int
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n              ! indices
    integer :: ndims          ! dimension counter
    integer :: dimid(5)       ! dimension ids
    integer :: varid          ! variable id
    integer :: itmp           ! temporary
    character(len=256) :: str ! temporary
    character(len=32) :: subname='NCD_DEFVAR_REAL' ! subroutine name
!-----------------------------------------------------------------------

    if (.not. nc_masterproc) return

    ! Determine dimension ids for variable

    dimid(:) = 0
	ndims=0
    if (present(dim1name)) then
	   ndims=ndims+1
       call check_ret(nf_inq_dimid(ncid, dim1name, dimid(ndims)), subname)
    end if
    if (present(dim2name)) then
	   ndims=ndims+1
       call check_ret(nf_inq_dimid(ncid, dim2name, dimid(ndims)), subname)
    end if
    if (present(dim3name)) then
	   ndims=ndims+1
       call check_ret(nf_inq_dimid(ncid, dim3name, dimid(ndims)), subname)
    end if
    if (present(dim4name)) then
	   ndims=ndims+1
       call check_ret(nf_inq_dimid(ncid, dim4name, dimid(ndims)), subname)
    end if
    if (present(dim5name)) then
	   ndims=ndims+1
       call check_ret(nf_inq_dimid(ncid, dim5name, dimid(ndims)), subname)
    end if


    ! Define variable

    if (present(dim1name) .or. present(dim2name) .or. present(dim3name) .or. &
		present(dim4name) .or. present(dim5name)) then
	   call check_ret(nf_def_var(ncid, trim(varname), xtype, ndims, dimid(1:ndims), varid), subname)
    else
       call check_ret(nf_def_var(ncid, varname, xtype, 0, 0, varid), subname)
    end if
    if (present(long_name)) then
       call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
    end if
    if (present(units)) then
       call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
    end if
    if (present(cell_method)) then
       str = 'time: ' // trim(cell_method)
       call check_ret(nf_put_att_text(ncid, varid, 'cell_method', len_trim(str), trim(str)), subname)
    end if
    if (present(fill_value)) then
       call check_ret(nf_put_att_real(ncid, varid, '_FillValue', xtype, 1, fill_value), subname)
    end if
    if (present(missing_value)) then
       call check_ret(nf_put_att_real(ncid, varid, 'missing_value', xtype, 1, missing_value), subname)
    end if
    if (present(ifill_value)) then
       call check_ret(nf_put_att_int(ncid, varid, '_FillValue', xtype, 1, ifill_value), subname)
    end if
    if (present(imissing_value)) then
       call check_ret(nf_put_att_int(ncid, varid, 'missing_value', xtype, 1, imissing_value), subname)
    end if

  end subroutine ncd_defvar

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_int_1d
!
! !INTERFACE:

  subroutine ncd_iolocal_int_1d(varname, data, flag, ncid, &
		lb_lon, lb_lat, lb_lvl, lb_t, ub_lon, ub_lat, ub_lvl, ub_t, &
		long_name, units, readvar)
! 读/写局部一维实型变量:将一笔资料阵列写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! lb_lon--经向的起始号
! lb_lat--纬向的起始号
! lb_lvl--层次的起始号
! lb_t--时间的起始号
! ub_lon--经向的起始号
! ub_lat--纬向的起始号
! ub_lvl--层次的起始号
! ub_t--时间的起始号
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! I/O for 1d int field
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname            ! variable name
    integer         , intent(inout) :: data(:)            ! local decomposition data
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    integer         , optional, intent(in) :: lb_lon    ! start for longitude
    integer         , optional, intent(in) :: lb_lat    ! start for latitute sizes
    integer         , optional, intent(in) :: lb_lvl    ! start for level size
    integer         , optional, intent(in) :: lb_t      ! start for time size
    integer         , optional, intent(in) :: ub_lon    ! start for longitude
    integer         , optional, intent(in) :: ub_lat    ! start for latitute sizes
    integer         , optional, intent(in) :: ub_lvl    ! start for level size
    integer         , optional, intent(in) :: ub_t      ! start for time size
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! variable id
    integer :: ndim                     ! dimension counter
    integer :: start(4)                 ! starting indices for netcdf field
    integer :: count(4)                 ! count values for netcdf field
	character(len=32) :: inq_name		! inquid variable name
	character(len=8) :: inq_xtype		! inquid variable type
	integer :: inq_ndims				! inquid variable dimention
	integer :: inq_dimids(4)			! inquid variable dimention id
	character(len=255) :: inq_natts		! inquid variable attachment
    character(len=32) :: subname='NCD_IOLOCAL_INT_1D' ! subroutine name
    logical :: varpresent               ! if true, variable is on tape
!-----------------------------------------------------------------------

    ! Write field as 1d field 
	if (flag == 'write') then
		if (nc_masterproc) then
			call check_ret(nf_inq_varid(ncid, varname, varid), subname)
			! Write 1d field
			ndim=0
			count=1
			if (present(lb_lon) .and. present(ub_lon)) then
				ndim=ndim+1
				start(ndim)=lb_lon
				count(ndim)=ub_lon-lb_lon+1
			else if(present(lb_lon) .neqv. present(ub_lon))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_lat)) then
				ndim=ndim+1
				start(ndim)=lb_lat
				count(ndim)=ub_lat-lb_lat+1
			else if(present(lb_lat) .neqv. present(ub_lat))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_lvl)) then
				ndim=ndim+1
				start(ndim)=lb_lvl
				count(ndim)=ub_lvl-lb_lvl+1
			else if(present(lb_lvl) .neqv. present(ub_lvl))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_t)) then
				ndim=ndim+1
				start(ndim)=lb_t
				count(ndim)=ub_t-lb_t+1
			else if(present(lb_t) .neqv. present(lb_t))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif

			if ((.not. present(lb_lon)) .and. (.not. present(lb_lat)) .and. &
				(.not. present(lb_lvl)) .and. (.not. present(lb_t))) then
				call endrun('must specify one dimention!',subname)
			endif

			call check_ret(nf_put_vara_int(ncid, varid, start(1:ndim), count(1:ndim), data), subname)
			if (present(long_name)) then
				call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
			end if
			if (present(units)) then
				call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
			end if
		end if   ! end of if-nc_masterproc block
    ! Read field as 1d field 
	else if (flag == 'read') then
		if (nc_masterproc) then
			call check_var(ncid, varname, varid, varpresent)
			if (varpresent) then
				ndim=0
				count=1
				if (present(lb_lon) .and. present(ub_lon)) then
					ndim=ndim+1
					start(ndim)=lb_lon
					count(ndim)=ub_lon-lb_lon+1
				else if(present(lb_lon) .neqv. present(ub_lon))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_lat)) then
					ndim=ndim+1
					start(ndim)=lb_lat
					count(ndim)=ub_lat-lb_lat+1
				else if(present(lb_lat) .neqv. present(ub_lat))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_lvl)) then
					ndim=ndim+1
					start(ndim)=lb_lvl
					count(ndim)=ub_lvl-lb_lvl+1
				else if(present(lb_lvl) .neqv. present(ub_lvl))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_t)) then
					ndim=ndim+1
					start(ndim)=lb_t
					count(ndim)=ub_t-lb_t+1
				else if(present(lb_t) .neqv. present(lb_t))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif

				if ((.not. present(lb_lon)) .and. (.not. present(lb_lat)) .and. &
					(.not. present(lb_lvl)) .and. (.not. present(lb_t))) then
					call endrun('must specify one dimention!',subname)
				endif

				!read data
				call check_ret(nf_get_vara_int(ncid, varid, start(1:ndim), count(1:ndim), data), subname)
			else
				call endrun('the varibal does not difined!',subname)
			end if
		end if
		if (present(readvar)) readvar = varpresent
	end if

  end subroutine ncd_iolocal_int_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_real_1d
!
! !INTERFACE:
  subroutine ncd_iolocal_real_1d(varname, data, flag, ncid, &
		lb_lon, lb_lat, lb_lvl, lb_t, ub_lon, ub_lat, ub_lvl, ub_t, &
		long_name, units, readvar)
! 读/写局部一维实型变量:将一笔资料阵列写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! lb_lon--经向的起始号
! lb_lat--纬向的起始号
! lb_lvl--层次的起始号
! lb_t--时间的起始号
! ub_lon--经向的起始号
! ub_lat--纬向的起始号
! ub_lvl--层次的起始号
! ub_t--时间的起始号
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! I/O for 1d int field
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname            ! variable name
    real, intent(inout) :: data(:)            ! local decomposition data
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    integer         , optional, intent(in) :: lb_lon    ! start for longitude
    integer         , optional, intent(in) :: lb_lat    ! start for latitute sizes
    integer         , optional, intent(in) :: lb_lvl    ! start for level size
    integer         , optional, intent(in) :: lb_t      ! start for time size
    integer         , optional, intent(in) :: ub_lon    ! start for longitude
    integer         , optional, intent(in) :: ub_lat    ! start for latitute sizes
    integer         , optional, intent(in) :: ub_lvl    ! start for level size
    integer         , optional, intent(in) :: ub_t      ! start for time size
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! variable id
    integer :: ndim                     ! dimension counter
    integer :: start(4)                 ! starting indices for netcdf field
    integer :: count(4)                 ! count values for netcdf field
	character(len=32) :: inq_name		! inquid variable name
	character(len=8) :: inq_xtype		! inquid variable type
	integer :: inq_ndims				! inquid variable dimention
	integer :: inq_dimids(4)			! inquid variable dimention id
	character(len=255) :: inq_natts		! inquid variable attachment
    character(len=32) :: subname='NCD_IOLOCAL_REAL_1D' ! subroutine name
    logical :: varpresent               ! if true, variable is on tape
!-----------------------------------------------------------------------

    ! Write field as 1d field 
	if (flag == 'write') then
		if (nc_masterproc) then
			call check_ret(nf_inq_varid(ncid, varname, varid), subname)
			! Write 1d field
			ndim=0
			count=1
			if (present(lb_lon) .and. present(ub_lon)) then
				ndim=ndim+1
				start(ndim)=lb_lon
				count(ndim)=ub_lon-lb_lon+1
			else if(present(lb_lon) .neqv. present(ub_lon))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_lat)) then
				ndim=ndim+1
				start(ndim)=lb_lat
				count(ndim)=ub_lat-lb_lat+1
			else if(present(lb_lat) .neqv. present(ub_lat))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_lvl)) then
				ndim=ndim+1
				start(ndim)=lb_lvl
				count(ndim)=ub_lvl-lb_lvl+1
			else if(present(lb_lvl) .neqv. present(ub_lvl))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_t)) then
				ndim=ndim+1
				start(ndim)=lb_t
				count(ndim)=ub_t-lb_t+1
			else if(present(lb_t) .neqv. present(lb_t))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif

			if ((.not. present(lb_lon)) .and. (.not. present(lb_lat)) .and. &
				(.not. present(lb_lvl)) .and. (.not. present(lb_t))) then
				call endrun('must specify one dimention!',subname)
			endif

			call check_ret(nf_put_vara_real(ncid, varid, start(1:ndim), count(1:ndim), data), subname)
			if (present(long_name)) then
				call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
			end if
			if (present(units)) then
				call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
			end if
		end if   ! end of if-nc_masterproc block
    ! Read field as 1d field 
	else if (flag == 'read') then
		if (nc_masterproc) then
			call check_var(ncid, varname, varid, varpresent)
			if (varpresent) then
				ndim=0
				count=1
				if (present(lb_lon) .and. present(ub_lon)) then
					ndim=ndim+1
					start(ndim)=lb_lon
					count(ndim)=ub_lon-lb_lon+1
				else if(present(lb_lon) .neqv. present(ub_lon))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_lat)) then
					ndim=ndim+1
					start(ndim)=lb_lat
					count(ndim)=ub_lat-lb_lat+1
				else if(present(lb_lat) .neqv. present(ub_lat))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_lvl)) then
					ndim=ndim+1
					start(ndim)=lb_lvl
					count(ndim)=ub_lvl-lb_lvl+1
				else if(present(lb_lvl) .neqv. present(ub_lvl))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_t)) then
					ndim=ndim+1
					start(ndim)=lb_t
					count(ndim)=ub_t-lb_t+1
				else if(present(lb_t) .neqv. present(lb_t))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif

				if ((.not. present(lb_lon)) .and. (.not. present(lb_lat)) .and. &
					(.not. present(lb_lvl)) .and. (.not. present(lb_t))) then
					call endrun('must specify one dimention!',subname)
				endif

				!read data
				call check_ret(nf_get_vara_real(ncid, varid, start(1:ndim), count(1:ndim), data), subname)
			else
				call endrun('the varibal does not difined!',subname)
			end if
		end if
		if (present(readvar)) readvar = varpresent
	end if

  end subroutine ncd_iolocal_real_1d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_real_1d
!
! !INTERFACE:
  subroutine ncd_iolocal_double_1d(varname, data, flag, ncid, &
		lb_lon, lb_lat, lb_lvl, lb_t, ub_lon, ub_lat, ub_lvl, ub_t, &
		long_name, units, readvar)
! 读/写局部一维实型变量:将一笔资料阵列写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! lb_lon--经向的起始号
! lb_lat--纬向的起始号
! lb_lvl--层次的起始号
! lb_t--时间的起始号
! ub_lon--经向的起始号
! ub_lat--纬向的起始号
! ub_lvl--层次的起始号
! ub_t--时间的起始号
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! I/O for 1d int field
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname            ! variable name
    real*8, intent(inout) :: data(:)            ! local decomposition data
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    integer         , optional, intent(in) :: lb_lon    ! start for longitude
    integer         , optional, intent(in) :: lb_lat    ! start for latitute sizes
    integer         , optional, intent(in) :: lb_lvl    ! start for level size
    integer         , optional, intent(in) :: lb_t      ! start for time size
    integer         , optional, intent(in) :: ub_lon    ! start for longitude
    integer         , optional, intent(in) :: ub_lat    ! start for latitute sizes
    integer         , optional, intent(in) :: ub_lvl    ! start for level size
    integer         , optional, intent(in) :: ub_t      ! start for time size
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! variable id
    integer :: ndim                     ! dimension counter
    integer :: start(4)                 ! starting indices for netcdf field
    integer :: count(4)                 ! count values for netcdf field
	character(len=32) :: inq_name		! inquid variable name
	character(len=8) :: inq_xtype		! inquid variable type
	integer :: inq_ndims				! inquid variable dimention
	integer :: inq_dimids(4)			! inquid variable dimention id
	character(len=255) :: inq_natts		! inquid variable attachment
    character(len=32) :: subname='NCD_IOLOCAL_REAL_1D' ! subroutine name
    logical :: varpresent               ! if true, variable is on tape
!-----------------------------------------------------------------------

    ! Write field as 1d field 
	if (flag == 'write') then
		if (nc_masterproc) then
			call check_ret(nf_inq_varid(ncid, varname, varid), subname)
			! Write 1d field
			ndim=0
			count=1
			if (present(lb_lon) .and. present(ub_lon)) then
				ndim=ndim+1
				start(ndim)=lb_lon
				count(ndim)=ub_lon-lb_lon+1
			else if(present(lb_lon) .neqv. present(ub_lon))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_lat)) then
				ndim=ndim+1
				start(ndim)=lb_lat
				count(ndim)=ub_lat-lb_lat+1
			else if(present(lb_lat) .neqv. present(ub_lat))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_lvl)) then
				ndim=ndim+1
				start(ndim)=lb_lvl
				count(ndim)=ub_lvl-lb_lvl+1
			else if(present(lb_lvl) .neqv. present(ub_lvl))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_t)) then
				ndim=ndim+1
				start(ndim)=lb_t
				count(ndim)=ub_t-lb_t+1
			else if(present(lb_t) .neqv. present(lb_t))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif

			if ((.not. present(lb_lon)) .and. (.not. present(lb_lat)) .and. &
				(.not. present(lb_lvl)) .and. (.not. present(lb_t))) then
				call endrun('must specify one dimention!',subname)
			endif

			call check_ret(nf_put_vara_double(ncid, varid, start(1:ndim), count(1:ndim), data), subname)
			if (present(long_name)) then
				call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
			end if
			if (present(units)) then
				call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
			end if
		end if   ! end of if-nc_masterproc block
    ! Read field as 1d field 
	else if (flag == 'read') then
		if (nc_masterproc) then
			call check_var(ncid, varname, varid, varpresent)
			if (varpresent) then
				ndim=0
				count=1
				if (present(lb_lon) .and. present(ub_lon)) then
					ndim=ndim+1
					start(ndim)=lb_lon
					count(ndim)=ub_lon-lb_lon+1
				else if(present(lb_lon) .neqv. present(ub_lon))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_lat)) then
					ndim=ndim+1
					start(ndim)=lb_lat
					count(ndim)=ub_lat-lb_lat+1
				else if(present(lb_lat) .neqv. present(ub_lat))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_lvl)) then
					ndim=ndim+1
					start(ndim)=lb_lvl
					count(ndim)=ub_lvl-lb_lvl+1
				else if(present(lb_lvl) .neqv. present(ub_lvl))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_t)) then
					ndim=ndim+1
					start(ndim)=lb_t
					count(ndim)=ub_t-lb_t+1
				else if(present(lb_t) .neqv. present(lb_t))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif

				if ((.not. present(lb_lon)) .and. (.not. present(lb_lat)) .and. &
					(.not. present(lb_lvl)) .and. (.not. present(lb_t))) then
					call endrun('must specify one dimention!',subname)
				endif

				!read data
				call check_ret(nf_get_vara_double(ncid, varid, start(1:ndim), count(1:ndim), data), subname)
			else
				call endrun('the varibal does not difined!',subname)
			end if
		end if
		if (present(readvar)) readvar = varpresent
	end if

  end subroutine ncd_iolocal_double_1d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_int_2d
!
! !INTERFACE:
  subroutine ncd_iolocal_int_2d(varname, data, flag, ncid, &
		lb_lon, lb_lat, lb_lvl, lb_t, ub_lon, ub_lat, ub_lvl, ub_t, &
		long_name, units, readvar)
! 读/写局部二维整型变量:将一笔资料阵列写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! lb_lon--经向的起始号
! lb_lat--纬向的起始号
! lb_lvl--层次的起始号
! lb_t--时间的起始号
! ub_lon--经向的起始号
! ub_lat--纬向的起始号
! ub_lvl--层次的起始号
! ub_t--时间的起始号
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! I/O for 2d real field
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname            ! variable name
    integer         , intent(inout) :: data(:,:)            ! local decomposition data
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    integer         , optional, intent(in) :: lb_lon    ! start for longitude
    integer         , optional, intent(in) :: lb_lat    ! start for latitute sizes
    integer         , optional, intent(in) :: lb_lvl    ! start for level size
    integer         , optional, intent(in) :: lb_t      ! start for time size
    integer         , optional, intent(in) :: ub_lon    ! start for longitude
    integer         , optional, intent(in) :: ub_lat    ! start for latitute sizes
    integer         , optional, intent(in) :: ub_lvl    ! start for level size
    integer         , optional, intent(in) :: ub_t      ! start for time size
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! variable id
    integer :: ndim                     ! dimension counter
    integer :: start(4)                 ! starting indices for netcdf field
    integer :: count(4)                 ! count values for netcdf field
	character(len=32) :: inq_name		! inquid variable name
	character(len=8) :: inq_xtype		! inquid variable type
	integer :: inq_ndims				! inquid variable dimention
	integer :: inq_dimids(4)			! inquid variable dimention id
	character(len=255) :: inq_natts		! inquid variable attachment
    character(len=32) :: subname='NCD_IOLOCAL_INT_2D' ! subroutine name
    logical :: varpresent               ! if true, variable is on tape
!-----------------------------------------------------------------------

    ! Write field as 2d field 
	if (flag == 'write') then
		if (nc_masterproc) then
			call check_ret(nf_inq_varid(ncid, varname, varid), subname)
			! Write 2d field
			ndim=0
			count=1
			if (present(lb_lon) .and. present(ub_lon)) then
				ndim=ndim+1
				start(ndim)=lb_lon
				count(ndim)=ub_lon-lb_lon+1
			else if(present(lb_lon) .neqv. present(ub_lon))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_lat)) then
				ndim=ndim+1
				start(ndim)=lb_lat
				count(ndim)=ub_lat-lb_lat+1
			else if(present(lb_lat) .neqv. present(ub_lat))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_lvl)) then
				ndim=ndim+1
				start(ndim)=lb_lvl
				count(ndim)=ub_lvl-lb_lvl+1
			else if(present(lb_lvl) .neqv. present(ub_lvl))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_t)) then
				ndim=ndim+1
				start(ndim)=lb_t
				count(ndim)=ub_t-lb_t+1
			else if(present(lb_t) .neqv. present(lb_t))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif

			if ((.not. present(lb_lon)) .and. (.not. present(lb_lat)) .and. &
				(.not. present(lb_lvl)) .and. (.not. present(lb_t))) then
				call endrun('must specify one dimention!',subname)
			endif

			call check_ret(nf_put_vara_int(ncid, varid, start(1:ndim), count(1:ndim), data), subname)
			if (present(long_name)) then
				call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
			end if
			if (present(units)) then
				call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
			end if
		end if   ! end of if-nc_masterproc block
    ! Read field as 1d field 
	else if (flag == 'read') then
		if (nc_masterproc) then
			call check_var(ncid, varname, varid, varpresent)
			if (varpresent) then
				ndim=0
				count=1
				if (present(lb_lon) .and. present(ub_lon)) then
					ndim=ndim+1
					start(ndim)=lb_lon
					count(ndim)=ub_lon-lb_lon+1
				else if(present(lb_lon) .neqv. present(ub_lon))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_lat)) then
					ndim=ndim+1
					start(ndim)=lb_lat
					count(ndim)=ub_lat-lb_lat+1
				else if(present(lb_lat) .neqv. present(ub_lat))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_lvl)) then
					ndim=ndim+1
					start(ndim)=lb_lvl
					count(ndim)=ub_lvl-lb_lvl+1
				else if(present(lb_lvl) .neqv. present(ub_lvl))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_t)) then
					ndim=ndim+1
					start(ndim)=lb_t
					count(ndim)=ub_t-lb_t+1
				else if(present(lb_t) .neqv. present(lb_t))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif

				if ((.not. present(lb_lon)) .and. (.not. present(lb_lat)) .and. &
					(.not. present(lb_lvl)) .and. (.not. present(lb_t))) then
					call endrun('must specify one dimention!',subname)
				endif

				call check_ret(nf_get_vara_int(ncid, varid, start(1:ndim), count(1:ndim), data), subname)
			else
				call endrun('the varibal does not difined!',subname)
			end if
		end if
		if (present(readvar)) readvar = varpresent
	end if

  end subroutine ncd_iolocal_int_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_real_2d
!
! !INTERFACE:
  subroutine ncd_iolocal_real_2d(varname, data, flag, ncid, &
		lb_lon, lb_lat, lb_lvl, lb_t, ub_lon, ub_lat, ub_lvl, ub_t, &
		long_name, units, readvar)
! 读/写局部二维实型变量:将一笔资料阵列写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! lb_lon--经向的起始号
! lb_lat--纬向的起始号
! lb_lvl--层次的起始号
! lb_t--时间的起始号
! ub_lon--经向的起始号
! ub_lat--纬向的起始号
! ub_lvl--层次的起始号
! ub_t--时间的起始号
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! I/O for 2d real field
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname            ! variable name
    real, intent(inout) :: data(:,:)            ! local decomposition data
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    integer         , optional, intent(in) :: lb_lon    ! start for longitude
    integer         , optional, intent(in) :: lb_lat    ! start for latitute sizes
    integer         , optional, intent(in) :: lb_lvl    ! start for level size
    integer         , optional, intent(in) :: lb_t      ! start for time size
    integer         , optional, intent(in) :: ub_lon    ! start for longitude
    integer         , optional, intent(in) :: ub_lat    ! start for latitute sizes
    integer         , optional, intent(in) :: ub_lvl    ! start for level size
    integer         , optional, intent(in) :: ub_t      ! start for time size
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! variable id
    integer :: ndim                     ! dimension counter
    integer :: start(4)                 ! starting indices for netcdf field
    integer :: count(4)                 ! count values for netcdf field
	character(len=32) :: inq_name		! inquid variable name
	character(len=8) :: inq_xtype		! inquid variable type
	integer :: inq_ndims				! inquid variable dimention
	integer :: inq_dimids(4)			! inquid variable dimention id
	character(len=255) :: inq_natts		! inquid variable attachment
    character(len=32) :: subname='NCD_IOLOCAL_REAL_2D' ! subroutine name
    logical :: varpresent               ! if true, variable is on tape
!-----------------------------------------------------------------------

    ! Write field as 2d field 
	if (flag == 'write') then
		if (nc_masterproc) then
			call check_ret(nf_inq_varid(ncid, varname, varid), subname)
			! Write 2d field
			ndim=0
			count=1
			if (present(lb_lon) .and. present(ub_lon)) then
				ndim=ndim+1
				start(ndim)=lb_lon
				count(ndim)=ub_lon-lb_lon+1
			else if(present(lb_lon) .neqv. present(ub_lon))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_lat)) then
				ndim=ndim+1
				start(ndim)=lb_lat
				count(ndim)=ub_lat-lb_lat+1
			else if(present(lb_lat) .neqv. present(ub_lat))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_lvl)) then
				ndim=ndim+1
				start(ndim)=lb_lvl
				count(ndim)=ub_lvl-lb_lvl+1
			else if(present(lb_lvl) .neqv. present(ub_lvl))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_t)) then
				ndim=ndim+1
				start(ndim)=lb_t
				count(ndim)=ub_t-lb_t+1
			else if(present(lb_t) .neqv. present(lb_t))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif

			if ((.not. present(lb_lon)) .and. (.not. present(lb_lat)) .and. &
				(.not. present(lb_lvl)) .and. (.not. present(lb_t))) then
				call endrun('must specify one dimention!',subname)
			endif

			call check_ret(nf_put_vara_real(ncid, varid, start(1:ndim), count(1:ndim), data), subname)
			if (present(long_name)) then
				call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
			end if
			if (present(units)) then
				call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
			end if
		end if   ! end of if-nc_masterproc block
    ! Read field as 1d field 
	else if (flag == 'read') then
		if (nc_masterproc) then
			call check_var(ncid, varname, varid, varpresent)
			if (varpresent) then
				ndim=0
				count=1
				if (present(lb_lon) .and. present(ub_lon)) then
					ndim=ndim+1
					start(ndim)=lb_lon
					count(ndim)=ub_lon-lb_lon+1
				else if(present(lb_lon) .neqv. present(ub_lon))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_lat)) then
					ndim=ndim+1
					start(ndim)=lb_lat
					count(ndim)=ub_lat-lb_lat+1
				else if(present(lb_lat) .neqv. present(ub_lat))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_lvl)) then
					ndim=ndim+1
					start(ndim)=lb_lvl
					count(ndim)=ub_lvl-lb_lvl+1
				else if(present(lb_lvl) .neqv. present(ub_lvl))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_t)) then
					ndim=ndim+1
					start(ndim)=lb_t
					count(ndim)=ub_t-lb_t+1
				else if(present(lb_t) .neqv. present(lb_t))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif

				if ((.not. present(lb_lon)) .and. (.not. present(lb_lat)) .and. &
					(.not. present(lb_lvl)) .and. (.not. present(lb_t))) then
					call endrun('must specify one dimention!',subname)
				endif

				call check_ret(nf_get_vara_real(ncid, varid, start(1:ndim), count(1:ndim), data), subname)
			else
				call endrun('the varibal does not difined!',subname)
			end if
		end if
		if (present(readvar)) readvar = varpresent
	end if

  end subroutine ncd_iolocal_real_2d


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_iolocal_real_2d
!
! !INTERFACE:
  subroutine ncd_iolocal_double_2d(varname, data, flag, ncid, &
		lb_lon, lb_lat, lb_lvl, lb_t, ub_lon, ub_lat, ub_lvl, ub_t, &
		long_name, units, readvar)
! 读/写局部二维实型变量:将一笔资料阵列写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! lb_lon--经向的起始号
! lb_lat--纬向的起始号
! lb_lvl--层次的起始号
! lb_t--时间的起始号
! ub_lon--经向的起始号
! ub_lat--纬向的起始号
! ub_lvl--层次的起始号
! ub_t--时间的起始号
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! I/O for 2d real field
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname            ! variable name
    real*8, intent(inout) :: data(:,:)            ! local decomposition data
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    integer         , optional, intent(in) :: lb_lon    ! start for longitude
    integer         , optional, intent(in) :: lb_lat    ! start for latitute sizes
    integer         , optional, intent(in) :: lb_lvl    ! start for level size
    integer         , optional, intent(in) :: lb_t      ! start for time size
    integer         , optional, intent(in) :: ub_lon    ! start for longitude
    integer         , optional, intent(in) :: ub_lat    ! start for latitute sizes
    integer         , optional, intent(in) :: ub_lvl    ! start for level size
    integer         , optional, intent(in) :: ub_t      ! start for time size
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! variable id
    integer :: ndim                     ! dimension counter
    integer :: start(4)                 ! starting indices for netcdf field
    integer :: count(4)                 ! count values for netcdf field
	character(len=32) :: inq_name		! inquid variable name
	character(len=8) :: inq_xtype		! inquid variable type
	integer :: inq_ndims				! inquid variable dimention
	integer :: inq_dimids(4)			! inquid variable dimention id
	character(len=255) :: inq_natts		! inquid variable attachment
    character(len=32) :: subname='NCD_IOLOCAL_REAL_2D' ! subroutine name
    logical :: varpresent               ! if true, variable is on tape
!-----------------------------------------------------------------------

    ! Write field as 2d field 
	if (flag == 'write') then
		if (nc_masterproc) then
			call check_ret(nf_inq_varid(ncid, varname, varid), subname)
			! Write 2d field
			ndim=0
			count=1
			if (present(lb_lon) .and. present(ub_lon)) then
				ndim=ndim+1
				start(ndim)=lb_lon
				count(ndim)=ub_lon-lb_lon+1
			else if(present(lb_lon) .neqv. present(ub_lon))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_lat)) then
				ndim=ndim+1
				start(ndim)=lb_lat
				count(ndim)=ub_lat-lb_lat+1
			else if(present(lb_lat) .neqv. present(ub_lat))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_lvl)) then
				ndim=ndim+1
				start(ndim)=lb_lvl
				count(ndim)=ub_lvl-lb_lvl+1
			else if(present(lb_lvl) .neqv. present(ub_lvl))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif
			if (present(lb_t)) then
				ndim=ndim+1
				start(ndim)=lb_t
				count(ndim)=ub_t-lb_t+1
			else if(present(lb_t) .neqv. present(lb_t))then
				call endrun('must specify the low and up boundary at the same time!',subname)
			endif

			if ((.not. present(lb_lon)) .and. (.not. present(lb_lat)) .and. &
				(.not. present(lb_lvl)) .and. (.not. present(lb_t))) then
				call endrun('must specify one dimention!',subname)
			endif

			call check_ret(nf_put_vara_double(ncid, varid, start(1:ndim), count(1:ndim), data), subname)
			if (present(long_name)) then
				call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
			end if
			if (present(units)) then
				call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
			end if
		end if   ! end of if-nc_masterproc block
    ! Read field as 1d field 
	else if (flag == 'read') then
		if (nc_masterproc) then
			call check_var(ncid, varname, varid, varpresent)
			if (varpresent) then
				ndim=0
				count=1
				if (present(lb_lon) .and. present(ub_lon)) then
					ndim=ndim+1
					start(ndim)=lb_lon
					count(ndim)=ub_lon-lb_lon+1
				else if(present(lb_lon) .neqv. present(ub_lon))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_lat)) then
					ndim=ndim+1
					start(ndim)=lb_lat
					count(ndim)=ub_lat-lb_lat+1
				else if(present(lb_lat) .neqv. present(ub_lat))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_lvl)) then
					ndim=ndim+1
					start(ndim)=lb_lvl
					count(ndim)=ub_lvl-lb_lvl+1
				else if(present(lb_lvl) .neqv. present(ub_lvl))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif
				if (present(lb_t)) then
					ndim=ndim+1
					start(ndim)=lb_t
					count(ndim)=ub_t-lb_t+1
				else if(present(lb_t) .neqv. present(lb_t))then
					call endrun('must specify the low and up boundary at the same time!',subname)
				endif

				if ((.not. present(lb_lon)) .and. (.not. present(lb_lat)) .and. &
					(.not. present(lb_lvl)) .and. (.not. present(lb_t))) then
					call endrun('must specify one dimention!',subname)
				endif

				call check_ret(nf_get_vara_double(ncid, varid, start(1:ndim), count(1:ndim), data), subname)
			else
				call endrun('the varibal does not difined!',subname)
			end if
		end if
		if (present(readvar)) readvar = varpresent
	end if

  end subroutine ncd_iolocal_double_2d


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_var
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_var(varname, data, flag, ncid, long_name, units, nt, readvar)
! 读/写全局零维整型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! I/O of integer variable
!

! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname            ! variable name
	integer         , intent(inout)  :: data               ! local decomposition data
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
	integer         , intent(in)  :: ncid               ! input unit
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    integer         , optional, intent(in) :: nt        ! time sample index
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            ! error status
    integer :: dimid(1)                       ! dimension id
    integer :: start(1), count(1)             ! output bounds
    integer :: varid                          ! variable id
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_INT_VAR' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = nt; count(1) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
		  if (present(long_name)) then
		     call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
		  end if
		  if (present(units)) then
			 call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
		  end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
		      call check_ret(nf_get_var_int(ncid, varid, data), subname)
		  else
			  call endrun('the varibal does not difined!',subname)
		  endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_var
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_var(varname, data, flag, ncid, long_name, units, nt, readvar)
! 读/写全局零维实型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! I/O of real variable
!

! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname            ! variable name
    real        , intent(inout)  :: data               ! local decomposition data
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            ! error status
    integer :: dimid(1)                       ! dimension id
    integer :: start(1), count(1)             ! output bounds
    integer :: varid                          ! variable id
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_VAR' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = nt; count(1) = 1
             call check_ret(nf_put_vara_real(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_real(ncid, varid, data), subname)
          end if
		  if (present(long_name)) then
		     call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
		  end if
		  if (present(units)) then
			 call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
		  end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
		      call check_ret(nf_get_var_real(ncid, varid, data), subname)
		  else
			  call endrun('the varibal does not difined!',subname)
		  endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_var
!
! !INTERFACE:
  subroutine ncd_ioglobal_double_var(varname, data, flag, ncid, long_name, units, nt, readvar)
! 读/写全局零维实型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! I/O of real variable
!

! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: varname            ! variable name
    real*8        , intent(inout)  :: data               ! local decomposition data
    character(len=*), intent(in)  :: flag               ! 'read' or 'write'
    integer         , intent(in)  :: ncid               ! input unit
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ier                            ! error status
    integer :: dimid(1)                       ! dimension id
    integer :: start(1), count(1)             ! output bounds
    integer :: varid                          ! variable id
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_VAR' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = nt; count(1) = 1
             call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_double(ncid, varid, data), subname)
          end if
		  if (present(long_name)) then
		     call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
		  end if
		  if (present(units)) then
			 call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
		  end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
		      call check_ret(nf_get_var_double(ncid, varid, data), subname)
		  else
			  call endrun('the varibal does not difined!',subname)
		  endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_double_var

!----------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_1d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_1d(varname, data, flag, ncid, long_name, units, nt, readvar)
! 读/写全局一维整型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! Master I/O for 1d integer data
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data(:)          ! local decomposition data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: dimid(2), ndims                ! dimension ids
    integer :: start(2), count(2)             ! output bounds
    integer :: ier                            ! error code
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_INT_1D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data)
             start(2) = nt; count(2) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
		  if (present(long_name)) then
		     call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
		  end if
		  if (present(units)) then
			 call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
		  end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
		      call check_ret(nf_get_var_int(ncid, varid, data), subname)
		  else
			  call endrun('the varibal does not difined!',subname)
		  endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_1d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_1d(varname, data, flag, ncid, long_name, units, nt, readvar)
! 读/写全局一维实型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! Master I/O for 1d real data
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    real        , intent(inout) :: data(:)          ! local decomposition input data
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: ier                            ! error code
    integer :: dimid(2), ndims                ! dimension ids
    integer :: start(2), count(2)             ! output bounds
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_1D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data)
             start(2) = nt; count(2) = 1
             call check_ret(nf_put_vara_real(ncid, varid, start, count, data), subname)
          else
!             call check_ret(nf_put_var_real(ncid, varid, data), subname)
call check_ret(nf_put_var_real(ncid, varid, data), subname)
          end if
		  if (present(long_name)) then
		     call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
		  end if
		  if (present(units)) then
			 call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
		  end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
		      call check_ret(nf_get_var_real(ncid, varid, data), subname)
		  else
			  call endrun('the varibal does not difined!',subname)
		  endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_1d
!
! !INTERFACE:
  subroutine ncd_ioglobal_double_1d(varname, data, flag, ncid, long_name, units, nt, readvar)
! 读/写全局一维实型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! Master I/O for 1d real data
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    real*8        , intent(inout) :: data(:)          ! local decomposition input data
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: ier                            ! error code
    integer :: dimid(2), ndims                ! dimension ids
    integer :: start(2), count(2)             ! output bounds
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_1D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data)
             start(2) = nt; count(2) = 1
             call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
          else
!             call check_ret(nf_put_var_double(ncid, varid, data), subname)
call check_ret(nf_put_var_double(ncid, varid, data), subname)
          end if
		  if (present(long_name)) then
		     call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
		  end if
		  if (present(units)) then
			 call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
		  end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
		      call check_ret(nf_get_var_double(ncid, varid, data), subname)
		  else
			  call endrun('the varibal does not difined!',subname)
		  endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_double_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_2d(varname, data, flag, ncid, long_name, units, nt, readvar)
! 读/写全局二维整型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! netcdf I/O of global 2d integer array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data(:,:)        ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: dimid(3), ndims                ! dimension ids
    integer :: start(3), count(3)             ! output bounds
    integer :: ier                            ! error code
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_2D_INT_IO' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = nt; count(3) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
		  if (present(long_name)) then
		     call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
		  end if
		  if (present(units)) then
			 call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
		  end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
		      call check_ret(nf_get_var_int(ncid, varid, data), subname)
		  else
			  call endrun('the varibal does not difined!',subname)
		  endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_2d

!-----------------------------------------------------------------------

!BOP
!
! !IROUTINE: ncd_ioglobal_int_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_long_2d(varname, data, flag, ncid, long_name, units, nt, readvar)
! 读/写全局二维整型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! netcdf I/O of global 2d integer array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer*8         , intent(inout) :: data(:,:)        ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: dimid(3), ndims                ! dimension ids
    integer :: start(3), count(3)             ! output bounds
    integer :: ier                            ! error code
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_2D_INT_IO' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = nt; count(3) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
      if (present(long_name)) then
         call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
      end if
      if (present(units)) then
       call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
      end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
          call check_ret(nf_get_var_int(ncid, varid, data), subname)
      else
        call endrun('the varibal does not difined!',subname)
      endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_long_2d

!-----------------------------------------------------------------------

!BOP
!
! !IROUTINE: ncd_ioglobal_byte_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_byte_2d(varname, data, flag, ncid, long_name, units, nt, readvar)
! 读/写全局二维整型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! netcdf I/O of global 2d integer array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    byte, intent(inout) :: data(:,:)        ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: dimid(3), ndims                ! dimension ids
    integer :: start(3), count(3)             ! output bounds
    integer :: ier                            ! error code
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_2D_INT1_IO' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = nt; count(3) = 1
             call check_ret(nf_put_vara_int1(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int1(ncid, varid, data), subname)
          end if
      if (present(long_name)) then
         call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
      end if
      if (present(units)) then
       call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
      end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
          call check_ret(nf_get_var_int1(ncid, varid, data), subname)
      else
        call endrun('the varibal does not difined!',subname)
      endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_byte_2d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_short_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_short_2d(varname, data, flag, ncid, long_name, units, nt, readvar)
! 读/写全局二维整型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! netcdf I/O of global 2d integer array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer*2, intent(inout) :: data(:,:)        ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: dimid(3), ndims                ! dimension ids
    integer :: start(3), count(3)             ! output bounds
    integer :: ier                            ! error code
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_2D_INT2_IO' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = nt; count(3) = 1
             call check_ret(nf_put_vara_int2(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int2(ncid, varid, data), subname)
          end if
      if (present(long_name)) then
         call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
      end if
      if (present(units)) then
       call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
      end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
          call check_ret(nf_get_var_int2(ncid, varid, data), subname)
      else
        call endrun('the varibal does not difined!',subname)
      endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_short_2d  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_2d(varname, data, flag, &
                                  ncid, long_name, units, nt, readvar)
! 读/写全局二维实型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! long_name--属性-变量全称
! units--属性-变量单位
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! netcdf I/O of global 2d real array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real        , intent(inout) :: data(:,:)        ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: ier                            ! error code
    integer :: dimid(3), ndims                ! dimension ids
    integer :: start(3), count(3)             ! output bounds
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_2D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = nt; count(3) = 1
!             call check_ret(nf_put_vara_real(ncid, varid, start, count, data), subname)
call check_ret(nf_put_vara_real(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_real(ncid, varid, data), subname)
          end if
		  if (present(long_name)) then
		     call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
		  end if
		  if (present(units)) then
			 call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
		  end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
		      call check_ret(nf_get_var_real(ncid, varid, data), subname)
		  else
			  call endrun('the varibal does not difined!',subname)
		  endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_2d
!
! !INTERFACE:
  subroutine ncd_ioglobal_double_2d(varname, data, flag, &
                                  ncid, long_name, units, nt, readvar)
! 读/写全局二维实型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! long_name--属性-变量全称
! units--属性-变量单位
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! netcdf I/O of global 2d real array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real*8        , intent(inout) :: data(:,:)        ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                          ! netCDF variable id
    integer :: ier                            ! error code
    integer :: dimid(3), ndims                ! dimension ids
    integer :: start(3), count(3)             ! output bounds
    logical :: varpresent                     ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_2D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = nt; count(3) = 1
!             call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_double(ncid, varid, data), subname)
          end if
		  if (present(long_name)) then
		     call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
		  end if
		  if (present(units)) then
			 call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
		  end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
          if (varpresent) then
		      call check_ret(nf_get_var_double(ncid, varid, data), subname)
		  else
			  call endrun('the varibal does not difined!',subname)
		  endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_double_2d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_short_3d
!
! !INTERFACE:
  subroutine ncd_ioglobal_short_3d(varname, data, flag, &
                                 ncid, long_name, units, nt, readvar)
! 读/写全局三维整型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! long_name--属性-变量全称
! units--属性-变量单位
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! netcdf I/O of global 3d integer array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer*2       , intent(inout) :: data(:,:,:)      ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! netCDF variable id
    integer :: dimid(4), ndims          ! dimension ids
    integer :: start(4), count(4)       ! output bounds
    integer :: ier                      ! error code
    logical :: varpresent               ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_3D_INT2_IO' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = 1;  count(3) = size(data, dim=3)
             start(4) = nt; count(4) = 1
             call check_ret(nf_put_vara_int2(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int2(ncid, varid, data), subname)
          end if
      if (present(long_name)) then
         call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
      end if
      if (present(units)) then
       call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
      end if
      end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
      if (varpresent) then
        if (present(nt)) then
         start(1) = 1;  count(1) = size(data, dim=1)
         start(2) = 1;  count(2) = size(data, dim=2)
         start(3) = 1;  count(3) = size(data, dim=3)
         start(4) = nt; count(4) = 1
         call check_ret(nf_get_vara_int2(ncid, varid, start, count, data), subname)
        else
         call check_ret(nf_get_var_int2(ncid, varid, data), subname)
        end if
      else
        call endrun('the varibal does not difined!',subname)
      endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_short_3d
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_int_3d
!
! !INTERFACE:
  subroutine ncd_ioglobal_int_3d(varname, data, flag, &
                                 ncid, long_name, units, nt, readvar)
! 读/写全局三维整型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! long_name--属性-变量全称
! units--属性-变量单位
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! netcdf I/O of global 3d integer array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    integer         , intent(inout) :: data(:,:,:)      ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! netCDF variable id
    integer :: dimid(4), ndims          ! dimension ids
    integer :: start(4), count(4)       ! output bounds
    integer :: ier                      ! error code
    logical :: varpresent               ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_3D_INT_IO' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = 1;  count(3) = size(data, dim=3)
             start(4) = nt; count(4) = 1
             call check_ret(nf_put_vara_int(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_int(ncid, varid, data), subname)
          end if
 		  if (present(long_name)) then
		     call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
		  end if
		  if (present(units)) then
			 call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
		  end if
      end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
		  if (varpresent) then
			  if (present(nt)) then
				 start(1) = 1;  count(1) = size(data, dim=1)
				 start(2) = 1;  count(2) = size(data, dim=2)
				 start(3) = 1;  count(3) = size(data, dim=3)
				 start(4) = nt; count(4) = 1
				 call check_ret(nf_get_vara_int(ncid, varid, start, count, data), subname)
			  else
				 call check_ret(nf_get_var_int(ncid, varid, data), subname)
			  end if
		  else
			  call endrun('the varibal does not difined!',subname)
		  endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_int_3d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_3d
!
! !INTERFACE:
  subroutine ncd_ioglobal_real_3d(varname, data, flag, &
                                  ncid, long_name, units, nt, readvar)
! 读/写全局三维实型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! long_name--属性-变量全称
! units--属性-变量单位
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! netcdf I/O of global 3d real array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real, intent(inout) :: data(:,:,:)      ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! netCDF variable id
    integer :: ier                      ! error code
    integer :: dimid(4), ndims          ! dimension ids
    integer :: start(4), count(4)       ! output bounds
    logical :: varpresent               ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_3D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = 1;  count(3) = size(data, dim=3)
             start(4) = nt; count(4) = 1
             call check_ret(nf_put_vara_real(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_real(ncid, varid, data), subname)
          end if
		  if (present(long_name)) then
		     call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
		  end if
		  if (present(units)) then
			 call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
		  end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
		  if (varpresent) then
			  if (present(nt)) then
				 start(1) = 1;  count(1) = size(data, dim=1)
				 start(2) = 1;  count(2) = size(data, dim=2)
				 start(3) = 1;  count(3) = size(data, dim=3)
				 start(4) = nt; count(4) = 1
				 call check_ret(nf_get_vara_real(ncid, varid, start, count, data), subname)
			  else
				 call check_ret(nf_get_var_real(ncid, varid, data), subname)
			  end if
		  else
			  call endrun('the varibal does not difined!',subname)
		  endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_real_3d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_ioglobal_real_3d
!
! !INTERFACE:
  subroutine ncd_ioglobal_double_3d(varname, data, flag, &
                                  ncid, long_name, units, nt, readvar)
! 读/写全局三维实型变量:将所有的资料阵列均写入档案中
! varname--变量名
! data--变量存储数组
! long_name--属性-变量全称
! units--属性-变量单位
! flag--读/写的标记
! ncid--NC文件对应的文件号
! nt--时间步
! readvar--所读取的变量是否存在于该NC文件中
! !DESCRIPTION:
! netcdf I/O of global 3d real array
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    integer         , intent(in)    :: ncid             ! input unit
    character(len=*), intent(in)    :: varname          ! variable name
    real*8, intent(inout) :: data(:,:,:)      ! local decomposition input data
    character(len=*), optional, intent(in) :: long_name ! variable long name
    character(len=*), optional, intent(in) :: units     ! variable units
    logical         , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: varid                    ! netCDF variable id
    integer :: ier                      ! error code
    integer :: dimid(4), ndims          ! dimension ids
    integer :: start(4), count(4)       ! output bounds
    logical :: varpresent               ! if true, variable is on tape
    character(len=32) :: subname='NCD_IOGLOBAL_REAL_3D' ! subroutine name
!-----------------------------------------------------------------------

    if (flag == 'write') then

       if (nc_masterproc) then
          call check_ret(nf_inq_varid(ncid, varname, varid), subname)
          if (present(nt)) then
             start(1) = 1;  count(1) = size(data, dim=1)
             start(2) = 1;  count(2) = size(data, dim=2)
             start(3) = 1;  count(3) = size(data, dim=3)
             start(4) = nt; count(4) = 1
             call check_ret(nf_put_vara_double(ncid, varid, start, count, data), subname)
          else
             call check_ret(nf_put_var_double(ncid, varid, data), subname)
          end if
		  if (present(long_name)) then
		     call check_ret(nf_put_att_text(ncid, varid, 'long_name', len_trim(long_name), trim(long_name)), subname)
		  end if
		  if (present(units)) then
			 call check_ret(nf_put_att_text(ncid, varid, 'units', len_trim(units), trim(units)), subname)
		  end if
       end if

    else if (flag == 'read') then

       if (nc_masterproc) then
          call check_var(ncid, varname, varid, varpresent)
		  if (varpresent) then
			  if (present(nt)) then
				 start(1) = 1;  count(1) = size(data, dim=1)
				 start(2) = 1;  count(2) = size(data, dim=2)
				 start(3) = 1;  count(3) = size(data, dim=3)
				 start(4) = nt; count(4) = 1
				 call check_ret(nf_get_vara_double(ncid, varid, start, count, data), subname)
			  else
				 call check_ret(nf_get_var_double(ncid, varid, data), subname)
			  end if
		  else
			  call endrun('the varibal does not difined!',subname)
		  endif
       end if
       if (present(readvar)) readvar = varpresent

    end if

  end subroutine ncd_ioglobal_double_3d

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: endrun
!
! !INTERFACE:
subroutine endrun(msg,subname)
!
! !DESCRIPTION:
! Abort the model for abnormal termination
   implicit none
! !ARGUMENTS:
   character(len=*), intent(in), optional :: msg    ! string to be printed
   character(len=*), intent(in), optional :: subname    ! subname

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

end module ncdio






