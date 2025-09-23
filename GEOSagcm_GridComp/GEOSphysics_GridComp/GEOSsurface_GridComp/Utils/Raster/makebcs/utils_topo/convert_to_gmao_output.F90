program create_example
use netcdf
implicit none


character(len=512) :: fin,str,fonc4,fobin,fsnc4,fsbin
integer :: im_world,jm_world
integer :: varid, lonid, latid
integer :: i, j, nc
integer :: zid,gwdid,trbid
integer :: ncid,rc
integer :: dimids(2)
integer :: status
character(len=:), allocatable :: Istr, Jstr
real :: g
real, allocatable :: z1d(:),turb1d(:),gwd1d(:)

real, allocatable :: a(:,:)
real  :: asum

logical :: isLL
integer :: nargs

nargs = command_argument_count()

isLL = .false.
do i=1,nargs
     call get_command_argument(i,str)
     select case(trim(str))
     case ('-i','--input')
         call get_command_argument(i+1,fin)
     case ('--im')
          call get_command_argument(i+1,str)
          read(str,'(I10)')im_world
          jm_world=im_world*6
     case ('--jm')
          call get_command_argument(i+1,str)
          read(str,'(I10)')jm_world
          isLL=.true.
     end select
enddo

g = 9.80616

allocate(z1d(im_world*jm_world),turb1d(im_world*jm_world),gwd1d(im_world*jm_world))
allocate(a(im_world,jm_world))

call check(nf90_open(fin,NF90_NOWRITE,ncid),"open")
call check(nf90_inq_varid(ncid,'PHIS',varid),"find phis")
call check(nf90_get_var(ncid,varid,z1d),"read phis")
call check(nf90_inq_varid(ncid,'SGH',varid),"find sgh")
call check(nf90_get_var(ncid,varid,gwd1d),"read sgh")
call check(nf90_inq_varid(ncid,'SGH30',varid),"find sgh30")
call check(nf90_get_var(ncid,varid,turb1d),"read sgh30")
call check(nf90_close(ncid),"close")

Istr  = i_to_string(im_world)
Jstr  = i_to_string(jm_world)
fsnc4 = Istr//'x'//Jstr//'.nc4'
fsbin = Istr//'x'//Jstr//'.data'
!if (isLL) then
!   if (im_world < 288) then
!      write(fsnc4,"(i3.3,'x',i2.2,'.nc4')")im_world,jm_world
!      write(fsbin,"(i3.3,'x',i2.2,'.data')")im_world,jm_world
!   else if (im_world < 1152) then
!      write(fsnc4,"(i3.3,'x',i3.3,'.nc4')")im_world,jm_world
!      write(fsbin,"(i3.3,'x',i3.3,'.data')")im_world,jm_world
!   else
!      write(fsnc4,"(i4.4,'x',i3.3,'.nc4')")im_world,jm_world
!      write(fsbin,"(i4.4,'x',i3.3,'.data')")im_world,jm_world
!   end if
!else
!   if (im_world < 168) then
!      write(fsnc4,"(i2.2,'x',i3.3,'.nc4')")im_world,jm_world
!      write(fsbin,"(i2.2,'x',i3.3,'.data')")im_world,jm_world
!   else if (im_world >= 168 .and. im_world < 1000) then
!      write(fsnc4,"(i3.3,'x',i4.4,'.nc4')")im_world,jm_world
!      write(fsbin,"(i3.3,'x',i4.4,'.data')")im_world,jm_world
!   else if (im_world >= 1000 .and. im_world < 1666) then
!      write(fsnc4,"(i4.4,'x',i4.4,'.nc4')")im_world,jm_world
!      write(fsbin,"(i4.4,'x',i4.4,'.data')")im_world,jm_world
!   else if (im_world >= 1666) then
!      write(fsnc4,"(i4.4,'x',i5.5,'.nc4')")im_world,jm_world
!      write(fsbin,"(i4.4,'x',i5.5,'.data')")im_world,jm_world
!   endif
!end if
   
! --- DYN (mean elevation in meters) -----------------------------------------
fonc4 = "gmted_DYN_ave_"//trim(fsnc4)
fobin = "gmted_DYN_ave_"//trim(fsbin)

call check(nf90_create(fonc4, IOR(NF90_CLOBBER, NF90_NETCDF4), ncid), "create dyn")
open(file=fobin, unit=21, form='unformatted')

call check(nf90_def_dim(ncid, "lon", im_world, lonid), "defdim lon dyn")
call check(nf90_def_dim(ncid, "lat", jm_world, latid), "defdim lat dyn")
call check(nf90_def_var(ncid, "z", NF90_FLOAT, (/lonid, latid/), varid), "defvar z")

call check(nf90_put_att(ncid, varid, "long_name",     "mean elevation"), "att z long_name")
call check(nf90_put_att(ncid, varid, "units",         "m"),               "att z units")
! optional nicety:
 call check(nf90_put_att(ncid, varid, "standard_name","surface_altitude"), "att z stdname")

call check(nf90_enddef(ncid), "enddef dyn")

nc = 0
do j = 1, jm_world
  do i = 1, im_world
    nc = nc + 1
    ! NCAR writes PHIS (m^2 s^-2); convert to meters. Guard sentinels if present.
    if (z1d(nc) < 0.9e36) then
      a(i,j) = z1d(nc) / g
    else
        a(i,j) = 0.0   ! or 1.0e36 if you decide to add _FillValue  
    end if
  end do
end do

if (isLL) then
  asum = 0.0
  do i = 1, im_world
    asum = asum + a(i,1)
  end do
  a(:,1) = asum / float(im_world)

  asum = 0.0
  do i = 1, im_world
    asum = asum + a(i,jm_world)
  end do
  a(:,jm_world) = asum / float(im_world)
end if

call check(nf90_put_var(ncid, varid, a), "putvar z")
call check(nf90_close(ncid), "close dyn")
write(21) a
close(21)

! --- GWD (write variance m^2) ----------------------------------------------
fonc4 = "gmted_GWD_var_"//trim(fsnc4)
fobin = "gmted_GWD_var_"//trim(fsbin)

call check(nf90_create(fonc4, IOR(NF90_CLOBBER, NF90_NETCDF4), ncid), "create gwd")
open(file=fobin, unit=21, form='unformatted')

call check(nf90_def_dim(ncid, "lon", im_world, lonid), "defdim lon gwd")
call check(nf90_def_dim(ncid, "lat", jm_world, latid), "defdim lat gwd")
call check(nf90_def_var(ncid, "gwd", NF90_FLOAT, (/lonid, latid/), varid), "defvar gwd")

call check(nf90_put_att(ncid, varid, "long_name", "variance of subgrid orography for GWD"), "att gwd long_name")
call check(nf90_put_att(ncid, varid, "units",     "m2"),                                   "att gwd units")

call check(nf90_enddef(ncid), "enddef gwd")

nc = 0
do j = 1, jm_world
  do i = 1, im_world
    nc = nc + 1
    ! Square std dev (m) from PE file to variance (m^2) here.
    ! Guard against sentinels/negatives so we don't square 1e36.
    if (gwd1d(nc) >= 0.0 .and. gwd1d(nc) < 0.9e36) then
      a(i,j) = gwd1d(nc) * gwd1d(nc)
    else
      a(i,j) = 1.0e36
    end if
  end do
end do

if (isLL) then
  asum = 0.0
  do i = 1, im_world
    asum = asum + a(i,1)
  end do
  a(:,1) = asum / float(im_world)

  asum = 0.0
  do i = 1, im_world
    asum = asum + a(i,jm_world)
  end do
  a(:,jm_world) = asum / float(im_world)
end if

call check(nf90_put_var(ncid, varid, a), "putvar gwd")
call check(nf90_close(ncid), "close gwd")
write(21) a
close(21)

! --- TURB (write variance m^2) ----------------------------------------------
fonc4 = "gmted_TRB_var_"//trim(fsnc4)
fobin = "gmted_TRB_var_"//trim(fsbin)

call check(nf90_create(fonc4, IOR(NF90_CLOBBER, NF90_NETCDF4), ncid), "create trb")
open(file=fobin, unit=21, form='unformatted')

call check(nf90_def_dim(ncid, "lon", im_world, lonid), "defdim lon trb")
call check(nf90_def_dim(ncid, "lat", jm_world, latid), "defdim lat trb")
call check(nf90_def_var(ncid, "trb", NF90_FLOAT, (/lonid, latid/), varid), "defvar trb")

call check(nf90_put_att(ncid, varid, "long_name", "variance of subgrid orography for TRB"), "att trb long_name")
call check(nf90_put_att(ncid, varid, "units",     "m2"),                                    "att trb units")

call check(nf90_enddef(ncid), "enddef trb")

nc = 0
do j = 1, jm_world
  do i = 1, im_world
    nc = nc + 1
    ! Square std dev (m) from PE file to variance (m^2) here.
    ! Guard against sentinels/negatives so we don't square 1e36.
    if (turb1d(nc) >= 0.0 .and. turb1d(nc) < 0.9e36) then
      a(i,j) = turb1d(nc) * turb1d(nc)
    else
      a(i,j) = 1.0e36
    end if
  end do
end do

if (isLL) then
  asum = 0.0
  do i = 1, im_world
    asum = asum + a(i,1)
  end do
  a(:,1) = asum / float(im_world)

  asum = 0.0
  do i = 1, im_world
    asum = asum + a(i,jm_world)
  end do
  a(:,jm_world) = asum / float(im_world)
end if

call check(nf90_put_var(ncid, varid, a), "putvar trb")
call check(nf90_close(ncid), "close trb")
write(21) a
close(21)
! ---------------------------------------------------------------------------

contains

   subroutine check(status,loc)

      integer, intent ( in) :: status
      character(len=*), intent ( in) :: loc

      if(status /= NF90_noerr) then
         write (*,*) "Error at ", loc
         write (*,*) nf90_strerror(status)
         stop "Stopped"
      end if

   end subroutine check

   function i_to_string(count) result(str)
      character(len=:), allocatable :: str
      integer, intent(in) :: count
      character(len=9)    :: buffer
      write(buffer,'(i0)') count
      str = trim(buffer)
   end function i_to_string

end program create_example

