subroutine latlon2cube(npx, npy, nlon, nlat, data_ll, data_cs)

 use ESMF
 use MAPL_Mod, only : MAPL_UNDEF
 use MAPL_IOMod, only : GETFILEUNIT, FREE_FILE
 use MAPL_ConstantsMod, only : pi=> MAPL_PI_R8
 use fv_grid_utils_mod, only : gnomonic_grids, cell_center2
 use fv_grid_tools_mod, only : mirror_grid
 use fv_arrays_mod,     only : REAL4, REAL8, R_GRID


 implicit none

 integer, intent(in) :: npx, npy, nlon, nlat
 real, dimension(npx , npy ), intent(out) :: data_cs
 real, dimension(nlon, nlat), intent(in ) :: data_ll

 integer :: ntiles=6
 integer :: npts
 integer :: ndims=2
 real(R_GRID), allocatable :: grid_global(:,:,:,:)
 real(R_GRID), allocatable :: agrid(:,:,:)

 real(R_GRID), allocatable :: xlon(:), ylat(:)

 real, dimension(nlon, nlat)           :: dll_flipped

 integer                               :: l2c_unit
 character(len=ESMF_MAXSTR)            :: l2c_fname
 logical                               :: l2c_file_exists

 real(ESMF_KIND_R4), allocatable :: l2c(:,:,:)
 integer           , allocatable :: id1(:,:), id2(:,:), jdc(:,:)
 logical, save :: do_init=.true.

 real(R_GRID) :: dlon, dlat

 logical :: found
 integer :: grid_type = 0

 integer :: i,j,n,l,itile, i1,i2, j1,j2

  npts = npx+1

 if (do_init) then

   write(l2c_fname,'(i5.5,"x",i5.5,"_l2c_",i5.5,"x",i5.5,".bin")') nlon, nlat, npx,npy

   inquire(FILE=TRIM(l2c_fname), EXIST=l2c_file_exists)

   if (.not. l2c_file_exists) then

   print*, 'Computing weights for ', TRIM(l2c_fname)
!--------------------------------------------------------------------!
! initialize cubed sphere grid                                       !
!--------------------------------------------------------------------!
  allocate( grid_global(npts,npts,ndims,ntiles) )
  call gnomonic_grids(grid_type, npts-1, grid_global(:,:,1,1), grid_global(:,:,2,1))
  call mirror_grid(grid_global, 0, npts, npts, 2, 6)
  do n=1,ntiles
     do j=1,npts
        do i=1,npts
!---------------------------------
! Shift the corner away from Japan
!---------------------------------
! This will result in the corner close to east coast of China
           grid_global(i,j,1,n) = grid_global(i,j,1,n) - pi/18.
           if ( grid_global(i,j,1,n) < 0. )              &
                grid_global(i,j,1,n) = grid_global(i,j,1,n) + 2.*pi
           if (ABS(grid_global(i,j,1,1)) < 1.e-10) grid_global(i,j,1,1) = 0.0
           if (ABS(grid_global(i,j,2,1)) < 1.e-10) grid_global(i,j,2,1) = 0.0
        enddo
     enddo
  enddo
!---------------------------------
! Clean Up Corners
!---------------------------------
  grid_global(  1,1:npts,:,2)=grid_global(npts,1:npts,:,1)
  grid_global(  1,1:npts,:,3)=grid_global(npts:1:-1,npts,:,1)
  grid_global(1:npts,npts,:,5)=grid_global(1,npts:1:-1,:,1)
  grid_global(1:npts,npts,:,6)=grid_global(1:npts,1,:,1)
  grid_global(1:npts,  1,:,3)=grid_global(1:npts,npts,:,2)
  grid_global(1:npts,  1,:,4)=grid_global(npts,npts:1:-1,:,2)
  grid_global(npts,1:npts,:,6)=grid_global(npts:1:-1,1,:,2)
  grid_global(  1,1:npts,:,4)=grid_global(npts,1:npts,:,3)
  grid_global(  1,1:npts,:,5)=grid_global(npts:1:-1,npts,:,3)
  grid_global(npts,1:npts,:,3)=grid_global(1,1:npts,:,4)
  grid_global(1:npts,  1,:,5)=grid_global(1:npts,npts,:,4)
  grid_global(1:npts,  1,:,6)=grid_global(npts,npts:1:-1,:,4)
  grid_global(  1,1:npts,:,6)=grid_global(npts,1:npts,:,5)

  !------------------------------------------------------------------!
  ! define the agrid at cell centers  0:360  -90:90                  !
  !------------------------------------------------------------------!
  allocate( agrid(npx,npy,ndims) )
  agrid(:,:,:) = -1.e25
  do n=1,ntiles
  do j=1,npts-1
     j1 = (npts-1)*(n-1) + j
     do i=1,npts-1
         call cell_center2(grid_global(i,j,  1:2,n), grid_global(i+1,j,  1:2,n),   &
                           grid_global(i,j+1,1:2,n), grid_global(i+1,j+1,1:2,n),   &
                           agrid(i,j1,1:2) )
         agrid(i,j1,1) = agrid(i,j1,1)
     enddo
  enddo
  enddo
  deallocate ( grid_global )

  !------------------------------------------------------------------!
  ! initialize latlon grid                                           !
  !------------------------------------------------------------------!
  allocate(xlon(nlon), ylat(nlat))
  !------------------------------------------------------------------!
  ! LatLon locations should be   0:360  -90:90                       !
  !------------------------------------------------------------------!
  dlon=(pi+pi)/real(nlon)
  dlat=pi/real(nlat-1)
  do i=1,nlon
     xlon(i)=(real(i)-1)*dlon
  enddo
  ylat(1)   =-0.5*pi
  ylat(nlat)= 0.5*pi
  do j=2,nlat-1
     ylat(j)=ylat(1)+(real(j)-1.)*dlat
  enddo

  !--------------------------------------------------------------------!
  ! allocate storage for l2c weights                                   !
  !--------------------------------------------------------------------!
  allocate( l2c(npx, npy, 4) )
  allocate( id1(npx, npy)    )
  allocate( id2(npx, npy)    )
  allocate( jdc(npx, npy)    )

  !--------------------------------------------------------------------!
  ! calculate weights for cubed sphere to latlon grid  interpolation   !
  !--------------------------------------------------------------------! 
  call remap_coef( npx, npy, nlon, nlat, agrid, xlon, ylat, id1, id2, jdc, l2c )

  deallocate ( agrid )
  deallocate ( xlon )
  deallocate ( ylat )
 !do_init = .false.

  l2c_unit = GETFILEUNIT(TRIM(l2c_fname))
  OPEN (UNIT=l2c_unit,FILE=TRIM(l2c_fname), FORM='unformatted')
  WRITE(UNIT=l2c_unit) l2c
  WRITE(UNIT=l2c_unit) id1
  WRITE(UNIT=l2c_unit) id2
  WRITE(UNIT=l2c_unit) jdc
  CLOSE(UNIT=l2c_unit)
  call FREE_FILE(l2c_unit)
  print*, 'Wrote weights for ', TRIM(l2c_fname)

  else

  !--------------------------------------------------------------------!
  ! allocate storage for l2c weights & read them in                    !
  !--------------------------------------------------------------------!
  allocate( l2c(npx, npy, 4) )
  allocate( id1(npx, npy)    )
  allocate( id2(npx, npy)    )
  allocate( jdc(npx, npy)    )
  l2c_unit = GETFILEUNIT(TRIM(l2c_fname))
 !print*, 'Reading weights from ', TRIM(l2c_fname)
  OPEN (UNIT=l2c_unit, FILE=TRIM(l2c_fname), FORM='unformatted', STATUS='OLD')
  READ (UNIT=l2c_unit) l2c
  READ (UNIT=l2c_unit) id1
  READ (UNIT=l2c_unit) id2
  READ (UNIT=l2c_unit) jdc
  CLOSE(UNIT=l2c_unit)
  call FREE_FILE(l2c_unit)
 !do_init = .false.

  endif

 endif ! do_init

  ! Latlon data needs to flip East/West hemispheres
  do j=1,nlat
     do i=1,nlon/2
        dll_flipped(i,j) = data_ll((nlon/2)+i,j)
        dll_flipped((nlon/2)+i,j) = data_ll(i,j)
     enddo
  enddo

  !--------------------------------------------------------------------!
  ! perform interpolation                                              !
  !--------------------------------------------------------------------!
  do j=1,npy
     do i=1,npx
        i1 = id1(i,j)
        i2 = id2(i,j)
        j1 = jdc(i,j)
        data_cs(i,j) = l2c(i,j,1)*dll_flipped(i1,j1  ) + l2c(i,j,2)*dll_flipped(i2,j1  ) +  &
                       l2c(i,j,3)*dll_flipped(i2,j1+1) + l2c(i,j,4)*dll_flipped(i1,j1+1)
     enddo
  enddo

  deallocate( l2c )
  deallocate( id1 )
  deallocate( id2 )
  deallocate( jdc )


contains


 subroutine remap_coef( npx, npy, im, jm, agrid, lon, lat, id1, id2, jdc, l2c )

  integer,   intent(in):: npx, npy
  integer,   intent(in):: im, jm
  real(R_GRID), intent(in):: agrid(npx,npy,2)
  real(R_GRID), intent(in):: lon(im), lat(jm)
  real(ESMF_KIND_R4), intent(out):: l2c(npx,npy,4)
  integer,   intent(out), dimension(npx,npy):: id1, id2, jdc
! local:
  real(R_GRID) :: rdlon(im)
  real(R_GRID) :: rdlat(jm)
  real(R_GRID) :: a1, b1
  integer i,j, i1, i2, jc, i0, j0

  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

  i1 = -999
  i2 = -999
  jc = -999

! * Interpolate to cubed sphere cell center
  do 5000 j=1,npy
     do i=1,npx

       if ( agrid(i,j,1)>lon(im) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
       elseif ( agrid(i,j,1)<lon(1) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
       else
            do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
               i1 = i0;  i2 = i0+1
               a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif
111    continue

       if ( agrid(i,j,2)<lat(1) ) then
            jc = 1
            b1 = 0.
       elseif ( agrid(i,j,2)>lat(jm) ) then
            jc = jm-1
            b1 = 1.
       else
          do j0=1,jm-1
          if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
               jc = j0
               b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
       endif
222    continue

       l2c(i,j,1) = (1.-a1) * (1.-b1)
       l2c(i,j,2) =     a1  * (1.-b1)
       l2c(i,j,3) =     a1  *     b1
       l2c(i,j,4) = (1.-a1) *     b1
       id1(i,j) = i1
       id2(i,j) = i2
       jdc(i,j) = jc
     enddo    !i-loop
5000 continue !j-loop

 end subroutine remap_coef


end subroutine latlon2cube

