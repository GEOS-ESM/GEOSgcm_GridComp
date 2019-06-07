  program dyn_interp_rst

  !--------------------------------------------------------------------!
  ! purpose: driver for interpolation between two cubed sphere grids   !
  !          with different spatial resolution for GEOS restarts       !
  !--------------------------------------------------------------------!
  use ESMF

! Cube to Cube Utilities
  use CUB2CUB_mod,       only : get_c2c_weight, do_c2c_interpolation
  use MAPL_ConstantsMod, only : pi=> MAPL_PI
  use fv_grid_utils_mod, only : gnomonic_grids
  use fv_grid_tools_mod, only : mirror_grid
  use GHOST_CUBSPH_mod,  only : B_grid, A_grid, ghost_cubsph_update

! Lat-lon to Cube Utilities
  use fms_mod,        only: fms_init, fms_end
  use fv_arrays_mod,  only: fv_atmos_type
  use fv_control_mod, only: npx,npy,npz, ntiles
  use fv_control_mod, only: fv_init, fv_end
  use fv_mp_mod,      only: gid, masterproc, tile, mp_gather
  use fv_grid_tools_mod,  only: grid, agrid, area, dx, dy, dxc, dyc
  use fv_surf_map_mod,    only: map_to_cubed_simple
  use external_ic_mod, only: cubed_a2d

  implicit none

  integer :: npx_in
  integer :: npy_in
  integer :: npx_out
  integer :: npy_out

  integer :: npts
  integer :: ndims=2
  real(ESMF_KIND_R8), allocatable :: xs(:,:), ys(:,:)
  real(ESMF_KIND_R8), allocatable :: grid_in(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: grid_out(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: corner_in(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: corner_out(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: weight_c2c(:,:,:,:)
  integer,            allocatable :: index_c2c(:,:,:,:)

  real(ESMF_KIND_R8) :: dlon, dlat
  real(ESMF_KIND_R4), allocatable :: r4latlon(:,:)
  real(ESMF_KIND_R8), allocatable :: r8latlon(:,:)
  real(ESMF_KIND_R8), allocatable :: r8tmp(:,:)
  real(ESMF_KIND_R8), allocatable :: r8_global(:,:,:)
  real(ESMF_KIND_R8), allocatable :: varo(:,:)
  type(fv_atmos_type) :: Atm(1)

  real(ESMF_KIND_R8), allocatable :: ua(:,:,:)
  real(ESMF_KIND_R8), allocatable :: va(:,:,:)
  real(ESMF_KIND_R8), allocatable :: ud(:,:,:)
  real(ESMF_KIND_R8), allocatable :: vd(:,:,:)

  real(ESMF_KIND_R8), allocatable :: ak(:), bk(:)
  integer :: IUNIT=15
  integer :: OUNIT=16

  integer :: header(6)
  integer :: grid_type = 0
  integer :: i,j,l,k,j1,j2,status

  character(30) :: myName

  character(30) :: str_arg

  external :: getarg, iargc
  integer iargc

  if (IARGC() /= 5) then
     print*, 'ABORT: need 5 arguments input_res_x,input_res_y and output_res_x,output_res_y and levels (No vertical interp supported yet)'
     stop
  endif
 
  CALL GETARG(1, str_arg)
  read (str_arg,'(I10)') npx_in
  CALL GETARG(2, str_arg)
  read (str_arg,'(I10)') npy_in
  CALL GETARG(3, str_arg)
  read (str_arg,'(I10)') npx_out
  CALL GETARG(4, str_arg)
  read (str_arg,'(I10)') npy_out
  CALL GETARG(5, str_arg)
  read (str_arg,'(I10)') npz
  print*, npx_in, npx_out, npz


 if (npx_in == npy_in) then  ! Cube to Cube
  !--------------------------------------------------------------------!
  ! initialize Input cubed sphere grid                                 !
  !--------------------------------------------------------------------!
  ntiles=6
  npts = npx_in+1
  allocate( xs(npts,npts) )
  allocate( ys(npts,npts) )
  allocate( grid_in(npts,npts,ndims,ntiles) )
  call gnomonic_grids(grid_type, npts-1, xs, ys)
  do j=1,npts
     do i=1,npts
        grid_in(i,j,1,1) = xs(i,j)
        grid_in(i,j,2,1) = ys(i,j)
     enddo
  enddo
  deallocate ( xs )
  deallocate ( ys )
  ! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi]
  call mirror_grid(grid_in, npts, npts, 2, 6)
  allocate( corner_in(ndims,0:npts+1,0:npts+1,ntiles) )
  corner_in(1,1:npts,1:npts,:) = grid_in(:,:,1,:)
  corner_in(2,1:npts,1:npts,:) = grid_in(:,:,2,:)
  !------------------------------------------------------------------!
  ! do halo update                                                   !
  !------------------------------------------------------------------!
  do l=1,ntiles
     corner_in(1:2,0     ,0     ,l)=0.
     corner_in(1:2,npts+1,0     ,l)=0.
     corner_in(1:2,0     ,npts+1,l)=0.
     corner_in(1:2,npts+1,npts+1,l)=0.
     call ghost_cubsph_update(corner_in(1,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1, 1, &
                              1, ntiles, 1, 1, l, B_grid)
     call ghost_cubsph_update(corner_in(2,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1, 1, &
                              1, ntiles, 1, 1, l, B_grid)
  enddo
  deallocate (grid_in)
 
  !--------------------------------------------------------------------!
  ! initialize Output cubed sphere grid                                !
  !--------------------------------------------------------------------!
  npts = npx_out+1
  allocate( xs(npts,npts) )
  allocate( ys(npts,npts) )
  allocate( grid_out(npts,npts,ndims,ntiles) )
  call gnomonic_grids(grid_type, npts-1, xs, ys)
  do j=1,npts
     do i=1,npts
        grid_out(i,j,1,1) = xs(i,j)
        grid_out(i,j,2,1) = ys(i,j)
     enddo
  enddo
  deallocate ( xs )
  deallocate ( ys )
  ! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi]
  call mirror_grid(grid_out, npts, npts, 2, 6)
  allocate( corner_out(ndims,0:npts+1,0:npts+1,ntiles) )
  corner_out(1,1:npts,1:npts,:) = grid_out(:,:,1,:)
  corner_out(2,1:npts,1:npts,:) = grid_out(:,:,2,:)
  !------------------------------------------------------------------!
  ! do halo update                                                   !
  !------------------------------------------------------------------!
  do l=1,ntiles
     corner_out(1:2,0     ,0     ,l)=0.
     corner_out(1:2,npts+1,0     ,l)=0.
     corner_out(1:2,0     ,npts+1,l)=0.
     corner_out(1:2,npts+1,npts+1,l)=0.
     call ghost_cubsph_update(corner_out(1,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1, 1, &
                              1, ntiles, 1, 1, l, B_grid)
     call ghost_cubsph_update(corner_out(2,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1, 1, &
                              1, ntiles, 1, 1, l, B_grid)
  enddo
  deallocate (grid_out)

  !--------------------------------------------------------------------!
  ! calculate weights and indices from bilinear interpolation          !
  ! from grid_in to grid_out                                           !
  !--------------------------------------------------------------------!
  allocate(index_c2c (3, npx_out, npy_out, ntiles))
  allocate(weight_c2c(4, npx_out, npy_out, ntiles))
  call get_c2c_weight(ntiles, npx_in+1, npy_in+1, corner_in,                &
                      npx_out+1, npy_out+1, corner_out,                     &
                      index_c2c,  weight_c2c)

  open(IUNIT,file='fvcore_internal_restart_in' ,access='sequential',form='unformatted',status='old')
  open(OUNIT,file='fvcore_internal_restart_out',access='sequential',form='unformatted')

! Headers
  print*, ' '
  print*, '-----------------------------'
  print*, ' '
  read (IUNIT, IOSTAT=status) header
  write(OUNIT) header
  print*, header
  print*, ' '
  print*, '-----------------------------'
  print*, ' '
  read (IUNIT, IOSTAT=status) header(1:5)
  print*, header(1:5)  
  header(1) = npx_out
  header(2) = npy_out*6
  write(OUNIT) header(1:5)
  print*, header(1:5)

! AK and BK
  allocate( ak(npz+1) )
  allocate( bk(npz+1) )
  read (IUNIT, IOSTAT=status) ak
  read (IUNIT, IOSTAT=status) bk
  write(OUNIT) ak
  write(OUNIT) bk
  deallocate( ak )
  deallocate( bk )

! U and V
#define FLOW
#if defined(FLOW)
  print*, 'U and V'
  call read_interp_write_flow(IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz, ntiles, &
                              index_c2c, weight_c2c, corner_in, corner_out)
#else
! U
  print*, 'U'
  call read_interp_write(IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz, ntiles, &
                              index_c2c, weight_c2c)
! V
  print*, 'V'
  call read_interp_write(IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz, ntiles, &
                              index_c2c, weight_c2c)
#endif
! PT
  print*, 'PT'
  call read_interp_write(IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz, ntiles, &
                              index_c2c, weight_c2c)
! PE
  print*, 'PE'
  call read_interp_write(IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz+1, ntiles, &
                              index_c2c, weight_c2c)
! PKZ
  print*, 'PKZ'
  call read_interp_write(IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz, ntiles, &
                              index_c2c, weight_c2c)

  close(IUNIT)
  close(OUNIT)

  deallocate(corner_in, corner_out, index_c2c, weight_c2c)
 
 else ! LAT-LON to Cube

!#define GRADS_READABLE

#ifndef GRADS_READABLE
! Start up FMS/MPP
  call fms_init()
  ntiles = 6
  npx = npx_out+1
  npy = npy_out+1
  call fv_init(Atm, 1800.d0)
#endif

! Init latlon grid
  allocate( grid_in(npx_in,npy_in,ndims,1) )
  dlon=(pi+pi)/real(npx_in)
  dlat=pi/real(npy_in)
  do j=1,npy_in
     do i=1,npx_in
        grid_in(i,j,1,1) = real(i)*dlon
        grid_in(i,j,2,1) = -0.5*pi + (real(j)-1.)*dlat
     enddo
  enddo

! Open Files
  print*, ' Open Files '
  open(IUNIT,file='fvcore_internal_restart_in' ,access='sequential',form='unformatted',status='old')
  open(OUNIT,file='fvcore_internal_restart_out',access='sequential',form='unformatted')

#ifdef GRADS_READABLE
  read (IUNIT, IOSTAT=status) header
  print*, header
  read (IUNIT, IOSTAT=status) header(1:5)
  print*, header(1:5)
! AK and BK
  allocate( ak(npz+1) )
  allocate( bk(npz+1) )
  read (IUNIT, IOSTAT=status) ak
  read (IUNIT, IOSTAT=status) bk
     do k=1,npz+1
        print*, ak(k), bk(k)
     enddo
  deallocate( ak )
  deallocate( bk )
  allocate ( r8latlon(npx_in,npy_in) )
  allocate ( r4latlon(npx_in,npy_in) )
  do i=1,5
  if (i==4) then
     npz = 73
  else
     npz=72
  endif
  do k=1,npz
    read (IUNIT, IOSTAT=status) r8latlon
    if (i==1) r8latlon(:,1) = r8latlon(:,2)*0.25
! Regrid from -180:180 to 0:360
    r4latlon(1           :npx_in/2,:) = r8latlon(npx_in/2 + 1 :npx_in  , :)
    r4latlon(npx_in/2 + 1:npx_in  ,:) = r8latlon(1            :npx_in/2, :)
    write(OUNIT) r4latlon
  enddo
  enddo
  deallocate( r8latlon )
  deallocate( r4latlon )

#else

! Headers
  if (gid==masterproc) print*, ' '
  if (gid==masterproc) print*, '-----------------------------'
  if (gid==masterproc) print*, ' '
  read (IUNIT, IOSTAT=status) header
  if (gid==masterproc) print*, header
  if (gid==masterproc) write(OUNIT) header
  if (gid==masterproc) print*, header
  if (gid==masterproc) print*, ' '
  if (gid==masterproc) print*, '-----------------------------'
  if (gid==masterproc) print*, ' '
  read (IUNIT, IOSTAT=status) header(1:5)
  if (gid==masterproc) print*, header(1:5)
  header(1) = npx_out
  header(2) = npy_out*6
  if (gid==masterproc) write(OUNIT) header(1:5)
  if (gid==masterproc) print*, header(1:5)

! AK and BK
  allocate( ak(npz+1) )
  allocate( bk(npz+1) )
  read (IUNIT, IOSTAT=status) ak
  read (IUNIT, IOSTAT=status) bk
  if (gid==masterproc) print*, ' '
  if (gid==masterproc) print*, '-----------------------------'
  if (gid==masterproc) print*, ' '
  if (gid==masterproc) then
     do k=1,npz+1
        print*, ak(k), bk(k)
     enddo
  endif
  if (gid==masterproc) print*, ' '
  if (gid==masterproc) print*, '-----------------------------'
  if (gid==masterproc) print*, ' '
  if (gid==masterproc) write(OUNIT) ak
  if (gid==masterproc) write(OUNIT) bk
  deallocate( ak )
  deallocate( bk )

! Read and regrid
  allocate ( r8latlon(npx_in,npy_in) )
  allocate ( r4latlon(npx_in,npy_in) )
  allocate ( r8tmp(Atm(1)%isd:Atm(1)%ied,Atm(1)%jsd:Atm(1)%jed) )
  allocate ( r8_global(npx_out,npy_out,ntiles) ) 
  allocate ( varo(npx_out,npy_out*ntiles) )

! State Vars

! Need to do U/V and move to cubed-sphere orientation
  allocate ( ua(Atm(1)%isd:Atm(1)%ied  ,Atm(1)%jsd:Atm(1)%jed  ,npz) )
  allocate ( va(Atm(1)%isd:Atm(1)%ied  ,Atm(1)%jsd:Atm(1)%jed  ,npz) )
  allocate ( ud(Atm(1)%isd:Atm(1)%ied  ,Atm(1)%jsd:Atm(1)%jed+1,npz) )
  allocate ( vd(Atm(1)%isd:Atm(1)%ied+1,Atm(1)%jsd:Atm(1)%jed  ,npz) )
  do i=1,2
  do k=1,npz
    if (gid==masterproc) print*, 'Working on ', i, ' of 5, Level=',k
    read (IUNIT, IOSTAT=status) r8latlon
! If U-Wind treat the Southern Pole as a reduced magnitude copy of next latitude
    if (i==1) r8latlon(:,1) = r8latlon(:,2)*0.25
! Regrid from -180:180 to 0:360
    r4latlon(1           :npx_in/2,:) = r8latlon(npx_in/2 + 1 :npx_in  , :)
    r4latlon(npx_in/2 + 1:npx_in  ,:) = r8latlon(1            :npx_in/2, :)
    call map_to_cubed_simple(npx_in, npy_in, grid_in(1,:,2,1), grid_in(:,1,1,1), r4latlon, grid, agrid, r8tmp, npx_out, npy_out)
    if (i==1) ua(:,:,k) = r8tmp
    if (i==2) va(:,:,k) = r8tmp
  enddo
  enddo
! Re-orient a2d3d lat-lon to cubed orientation and write
  call cubed_a2d(npx, npy, npz, ua, va, ud, vd )
  do i=1,2
  do k=1,npz
    if (i==1) then
    r8_global(Atm(1)%isc:Atm(1)%iec,Atm(1)%jsc:Atm(1)%jec,tile) = &
           ud(Atm(1)%isc:Atm(1)%iec,Atm(1)%jsc:Atm(1)%jec,k)
    else
    r8_global(Atm(1)%isc:Atm(1)%iec,Atm(1)%jsc:Atm(1)%jec,tile) = &
           vd(Atm(1)%isc:Atm(1)%iec,Atm(1)%jsc:Atm(1)%jec,k)
    endif
    call mp_gather(r8_global, Atm(1)%isc, Atm(1)%iec, Atm(1)%jsc, Atm(1)%jec, npx_out, npy_out, ntiles)
    write(myName, "(i1.1,'x',i2.2,'_Var')") i,k
    call Write_Profile(r8_global, npx_out, npy_out, ntiles, myName)
    do l=1,6
       j1 = (npy_out)*(l-1) + 1
       j2 = (npy_out)*(l-1) + npy_out
       varo(:,j1:j2)=r8_global(:,:,l)
    enddo
    if (gid==masterproc) write(OUNIT) varo
  enddo
  enddo
  deallocate( ua )
  deallocate( va )
  deallocate( ud )
  deallocate( vd )

 do i=3,5
  if (i==4) then
     npz = 73
  else
     npz=72
  endif
  do k=1,npz
    if (gid==masterproc) print*, 'Working on ', i, ' of 5, Level=',k
    read (IUNIT, IOSTAT=status) r8latlon
! If U-Wind treat the Southern Pole as a reduced magnitude copy of next latitude
    if (i==1) r8latlon(:,1) = r8latlon(:,2)*0.25
! Regrid from -180:180 to 0:360
    r4latlon(1           :npx_in/2,:) = r8latlon(npx_in/2 + 1 :npx_in  , :)
    r4latlon(npx_in/2 + 1:npx_in  ,:) = r8latlon(1            :npx_in/2, :)
    call map_to_cubed_simple(npx_in, npy_in, grid_in(1,:,2,1), grid_in(:,1,1,1), r4latlon, grid, agrid, r8tmp, npx_out, npy_out)
    r8_global(Atm(1)%isc:Atm(1)%iec,Atm(1)%jsc:Atm(1)%jec,tile) = &
        r8tmp(Atm(1)%isc:Atm(1)%iec,Atm(1)%jsc:Atm(1)%jec)
    call mp_gather(r8_global, Atm(1)%isc, Atm(1)%iec, Atm(1)%jsc, Atm(1)%jec, npx_out, npy_out, ntiles)

    write(myName, "(i1.1,'x',i2.2,'_Var')") i,k
    call Write_Profile(r8_global, npx_out, npy_out, ntiles, myName)

    do l=1,6
       j1 = (npy_out)*(l-1) + 1
       j2 = (npy_out)*(l-1) + npy_out
       varo(:,j1:j2)=r8_global(:,:,l)
    enddo
    if (gid==masterproc) write(OUNIT) varo
  enddo
 enddo

  deallocate ( r8latlon  )
  deallocate ( r4latlon  )
  deallocate ( r8tmp     )
  deallocate ( r8_global )
  deallocate ( varo      )


#endif

  close(IUNIT)
  close(OUNIT)

 call fms_end()

 endif


 contains
 
 subroutine read_interp_write(IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz, ntiles, &
                              index_c2c, weight_c2c)
  integer, intent(IN) :: IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz, ntiles
  integer           , intent(IN) ::  index_c2c(3, npx_out, npy_out, ntiles)
  real(ESMF_KIND_R8), intent(IN) :: weight_c2c(4, npx_out, npy_out, ntiles)

  real(ESMF_KIND_R8), allocatable :: vari(:,:)
  real(ESMF_KIND_R8), allocatable :: varo(:,:)
  real(ESMF_KIND_R8), allocatable :: var_in(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: var_out(:,:,:,:)

  integer :: i,j,l,k,j1,j2,status

  allocate( vari(npx_in,npy_in*ntiles) )
  allocate( varo(npx_out,npy_out*ntiles) )
  allocate( var_in(0:npx_in+1,0:npy_in+1,npz,ntiles) )
  allocate( var_out(npx_out,npy_out,npz,ntiles) )
! 
  do k=1,npz
    read (IUNIT, IOSTAT=status) vari
    do l=1,ntiles
       j1 = (npy_in)*(l-1) + 1
       j2 = (npy_in)*(l-1) + npx_in
       var_in(1:npx_in,1:npy_in,k,l)=vari(:,j1:j2)
    enddo

!if (k==1) then
!   print*, vari(:,1)
!endif

  enddo
  do l=1,ntiles
    call ghost_cubsph_update(var_in, 0, npx_in+1, 0, npy_in+1, npz, 1, ntiles,  &
                             1, npz, l, A_grid)
  enddo
  call do_c2c_interpolation(var_in, 0, npx_in+1, 0, npy_in+1, npz, ntiles, &
                            index_c2c, weight_c2c, npx_out, npy_out, var_out)
  do k=1,npz
    do l=1,ntiles
       j1 = (npy_out)*(l-1) + 1
       j2 = (npy_out)*(l-1) + npy_out
       varo(:,j1:j2)=var_out(:,:,k,l)
    enddo
    write(OUNIT) varo

!if (k==1) then
!   print*, '--------------'
!   print*, varo(:,1)
!endif
!stop

  enddo

  deallocate ( vari, varo, var_in, var_out )

 end subroutine read_interp_write

 subroutine read_interp_write_flow(IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz, ntiles, &
                              index_c2c, weight_c2c, corner_in, corner_out)

  use GRID_UTILS_mod, only: latlon2xyz
  use GRID_UTILS_mod,   only: get_dx, get_dxa, get_dy, get_dya,     &
                              get_center_vect, get_west_vect,       &
                              get_south_vect, get_cosa_center
  use FLOW_PROJ_mod,    only: d2a_vect, a2d_vect

  integer, intent(IN) :: IUNIT, OUNIT, npx_in, npy_in, npx_out, npy_out, npz, ntiles
  integer           , intent(IN) ::  index_c2c(3, npx_out, npy_out, ntiles)
  real(ESMF_KIND_R8), intent(IN) :: weight_c2c(4, npx_out, npy_out, ntiles)
  real(ESMF_KIND_R8), intent(IN) :: corner_in(2,0:npx_in+2,0:npy_in+2,ntiles)
  real(ESMF_KIND_R8), intent(IN) :: corner_out(2,0:npx_out+2,0:npy_out+2,ntiles)

  logical :: west_edge  = .true., east_edge  = .true.,                  &
             south_edge = .true., north_edge = .true.
  logical :: sw_corner = .true., se_corner = .true.,                    &
             nw_corner = .true., ne_corner = .true.
  logical :: edge_interp = .false.

  real(ESMF_KIND_R8), dimension(:,:,:,:), allocatable :: xyz_corner_in, xyz_corner_out

  real(ESMF_KIND_R8), allocatable :: vari(:,:)
  real(ESMF_KIND_R8), allocatable :: varo(:,:)
  real(ESMF_KIND_R8), allocatable :: u_in(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: u_out(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: v_in(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: v_out(:,:,:,:)

  integer :: i,j,l,k,n,j1,j2,itile,status
  integer :: nx_in, ny_in, nx_out, ny_out, nz
  real(ESMF_KIND_R8), dimension(:,:,:,:,:), allocatable :: vxyz_in, vxyz_out
  real(ESMF_KIND_R8), dimension(:,:,:), allocatable :: u, v ,ec1, ec2, ew1, ew2, es1, es2
  real(ESMF_KIND_R8), dimension(:,:), allocatable :: dx, dy, dxa, dya, rdxa, rdya, cosa_s, sina_s
  real(ESMF_KIND_R8), dimension(:), allocatable :: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n


            nx_in = npx_in+1
            ny_in = npy_in+1
            nx_out = npx_out+1
            ny_out = npy_out+1

    !------------------------------------------------------------------!
    ! calculate xyz cell corners and cell centers                      !
    !------------------------------------------------------------------!
    allocate(xyz_corner_in (3, 0:nx_in+1,  0:ny_in+1,  ntiles),        &
             xyz_corner_out(3, 0:nx_out+1, 0:ny_out+1, ntiles))
    do l=1,ntiles
       do j=0,ny_in+1
          do i=0,nx_in+1
            call latlon2xyz(corner_in(:,i,j,l), xyz_corner_in(:,i,j,l))
          enddo
       enddo
    enddo
    do l=1,ntiles
       do j=0,ny_out+1
          do i=0,nx_out+1
            call latlon2xyz(corner_out(:,i,j,l), xyz_corner_out(:,i,j,l))
          enddo
       enddo
    enddo


  allocate( vari(npx_in,npy_in*ntiles) )
  allocate( varo(npx_out,npy_out*ntiles) )
  allocate( u_in(npx_in,npy_in,npz,ntiles) )
  allocate( v_in(npx_in,npy_in,npz,ntiles) )
  allocate( u_out(npx_out,npy_out,npz,ntiles) )
  allocate( v_out(npx_out,npy_out,npz,ntiles) )


! Read U
  do k=1,npz
    read (IUNIT, IOSTAT=status) vari
    do l=1,ntiles
       j1 = (npy_in)*(l-1) + 1
       j2 = (npy_in)*(l-1) + npx_in
       u_in(1:npx_in,1:npy_in,k,l)=vari(:,j1:j2)
    enddo
  enddo
! Read U
  do k=1,npz
    read (IUNIT, IOSTAT=status) vari
    do l=1,ntiles
       j1 = (npy_in)*(l-1) + 1
       j2 = (npy_in)*(l-1) + npx_in
       v_in(1:npx_in,1:npy_in,k,l)=vari(:,j1:j2)
    enddo
  enddo


   ! Flow interpolation for U and V components
    do k=1,npz
            nz=1

            !----------------------------------------------------------!
            ! horizontal flow                                          !
            !----------------------------------------------------------!
            allocate(u(0:nx_in,0:ny_in+1,nz), v(0:nx_in+1,0:ny_in,nz), &
                     vxyz_in(3,0:nx_in,0:ny_in,nz,ntiles))
            allocate(dx(0:nx_in,0:ny_in+1), dxa(0:nx_in,0:ny_in), rdxa(0:nx_in,0:ny_in))
            allocate(dy(0:nx_in+1,0:ny_in), dya(0:nx_in,0:ny_in), rdya(0:nx_in,0:ny_in))
            allocate(ec1(3,0:nx_in,0:ny_in), ec2(3,0:nx_in,0:ny_in))
            allocate(cosa_s(0:nx_in,0:ny_in), sina_s(0:nx_in,0:ny_in))
            !----------------------------------------------------------!
            ! loop over tiles                                          !
            !----------------------------------------------------------!
            do itile=1,ntiles
               !-------------------------------------------------------!
               ! read horizontal flow                                  !
               !-------------------------------------------------------!
               u(1:nx_in-1,1:ny_in-1,1:)=u_in(1:nx_in-1,1:ny_in-1,k:k,itile)
               v(1:nx_in-1,1:ny_in-1,1:)=v_in(1:nx_in-1,1:ny_in-1,k:k,itile)
               !-------------------------------------------------------!
               ! fill shared edges on D-Grid                           !
               !-------------------------------------------------------!
               if (itile==1) u(1:nx_in-1,ny_in,1:) = -REVERSE( v_in(1        ,1:ny_in-1,k:k,3) ) 
               if (itile==2) u(1:nx_in-1,ny_in,1:) =         ( u_in(1:nx_in-1,1        ,k:k,3) )
               if (itile==3) u(1:nx_in-1,ny_in,1:) = -REVERSE( v_in(1        ,1:ny_in-1,k:k,5) )
               if (itile==4) u(1:nx_in-1,ny_in,1:) =         ( u_in(1:nx_in-1,1        ,k:k,5) )
               if (itile==5) u(1:nx_in-1,ny_in,1:) = -REVERSE( v_in(1        ,1:ny_in-1,k:k,1) )
               if (itile==6) u(1:nx_in-1,ny_in,1:) =         ( u_in(1:nx_in-1,1        ,k:k,1) )

               if (itile==1) v(nx_in,1:ny_in-1,1:) =         ( v_in(1        ,1:ny_in-1,k:k,2) )
               if (itile==2) v(nx_in,1:ny_in-1,1:) = -REVERSE( u_in(1:nx_in-1,1        ,k:k,4) )
               if (itile==3) v(nx_in,1:ny_in-1,1:) =         ( v_in(1        ,1:ny_in-1,k:k,4) )
               if (itile==4) v(nx_in,1:ny_in-1,1:) = -REVERSE( u_in(1:nx_in-1,1        ,k:k,6) )
               if (itile==5) v(nx_in,1:ny_in-1,1:) =         ( v_in(1        ,1:ny_in-1,k:k,6) )
               if (itile==6) v(nx_in,1:ny_in-1,1:) = -REVERSE( u_in(1:nx_in-1,1        ,k:k,2) )

               !-------------------------------------------------------!
               ! geometrical properties of input grid                  !
               !-------------------------------------------------------!
               call get_dx (xyz_corner_in(:,:,:,itile), 0, nx_in, 0, ny_in,             &
                                                        0, nx_in, 0, ny_in, dx)
               call get_dxa(xyz_corner_in(:,:,:,itile), 0, nx_in, 0, ny_in,             &
                                                        0, nx_in, 0, ny_in, dxa, rdxa=rdxa)
               call get_dy (xyz_corner_in(:,:,:,itile), 0, nx_in, 0, ny_in,             &
                                                        0, nx_in, 0, ny_in, dy)
               call get_dya(xyz_corner_in(:,:,:,itile), 0, nx_in, 0, ny_in,             &
                                                        0, nx_in, 0, ny_in, dya, rdya=rdya)
               call get_center_vect(xyz_corner_in(:,:,:,itile), 0, nx_in, 0, ny_in,     &
                                                                0, nx_in, 0, ny_in, ec1, ec2)
               call get_cosa_center(ec1, ec2, 0, nx_in, 0, ny_in,     &
                                              0, nx_in, 0, ny_in, cosa_s, sina_s)
               !-------------------------------------------------------!
               ! calculate flow vector for a-grid                      !
               !-------------------------------------------------------!
               call d2a_vect(u, v, dx, dy, rdxa, rdya, cosa_s, ec1, ec2, &
                             0, nx_in  , 0, ny_in  , 1, nz,            &
                             1, nx_in-1, 1, ny_in-1, 1, nz,            &
                             vxyz_in(:,:,:,:,itile))
            enddo 
            deallocate(u, v, dx, dy, dxa, dya, rdxa, rdya, ec1, ec2, cosa_s, sina_s)
            allocate(vxyz_out(3,0:nx_out,0:ny_out,nz,ntiles))
            !----------------------------------------------------------!
            ! ghost cell update of vxyz_in                               !
            !----------------------------------------------------------!
            vxyz_in(:,0,    0,    :,:)=0.
            vxyz_in(:,nx_in,0    ,:,:)=0.
            vxyz_in(:,nx_in,ny_in,:,:)=0.
            vxyz_in(:,0,    ny_in,:,:)=0.
            do n=1,3
               do itile=1,ntiles
                  call ghost_cubsph_update(vxyz_in(n,:,:,:,:), 0, nx_in, 0, ny_in, nz, 1, ntiles, &
                                           1, nz, itile, A_grid)
               enddo
            enddo
            !----------------------------------------------------------!
            ! do interpolation of flow vector                          !
            !----------------------------------------------------------!
            do n=1,3
               call do_c2c_interpolation(vxyz_in(n,:,:,:,:), 0, nx_in, 0, ny_in, nz, ntiles, &
                                         index_c2c, weight_c2c, nx_out-1, ny_out-1,          &
                                         vxyz_out(n,1:nx_out-1,1:ny_out-1,:,:))
            enddo
            deallocate(vxyz_in)
            !----------------------------------------------------------!
            ! loop over tiles                                          !
            !----------------------------------------------------------!
            allocate(u(0:nx_out,0:ny_out+1,nz), v(0:nx_out+1,0:ny_out,nz))
            allocate(ew1(3,0:nx_out,0:ny_out+1), ew2(3,0:nx_out,0:ny_out+1))
            allocate(es1(3,0:nx_out+1,0:ny_out), es2(3,0:nx_out+1,0:ny_out))
            allocate(edge_vect_w(0:ny_out), edge_vect_e(0:ny_out),    &
                     edge_vect_s(0:nx_out), edge_vect_n(0:nx_out))
            do itile=1,ntiles
               !-------------------------------------------------------!
               ! ghost cell update of vxyz_out                           !
               !-------------------------------------------------------!
               do n=1,3
                  call ghost_cubsph_update(vxyz_out(n,:,:,:,:), 0, nx_out, 0, ny_out, nz, 1, ntiles, &
                                           1, nz, itile, A_grid)
               enddo
               !-------------------------------------------------------!
               ! geometrical properties of output grid                 !
               !-------------------------------------------------------!
               call get_west_vect(xyz_corner_out(:,:,:,itile),                       &
                    0, nx_out, 0, ny_out, 1, nx_out-1, 1, ny_out-1,              &
                    west_edge, east_edge, 1, nx_out-1, ew1, ew2)
               call get_south_vect(xyz_corner_out(:,:,:,itile),                      &
                    0, nx_out, 0, ny_out, 1, nx_out-1, 1, ny_out-1,              &
                    west_edge, east_edge, 1, ny_out-1, es1, es2)
               !-------------------------------------------------------!
               ! calculate co-variant flow components on d_grid        !
               !-------------------------------------------------------!
               call a2d_vect(vxyz_out(:,:,:,:,itile), ew1, ew2, es1, es2,            &
                    0, nx_out,   0, ny_out,   1, nz,                               &
                    1, nx_out-1, 1, ny_out-1, 1, nz,                               &
                    edge_interp, edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n, &
                    west_edge, east_edge, south_edge, north_edge,                    &
                    sw_corner, se_corner, nw_corner, ne_corner,                      &
                    1, nx_out-1, 1, ny_out-1, u, v)
               !-------------------------------------------------------!
               ! write flow components                                 !
               !-------------------------------------------------------!
               u_out(1:nx_out-1, 1:ny_out-1, k:k,itile)=u(1:nx_out-1, 1:ny_out-1, 1:nz)
               v_out(1:nx_out-1, 1:ny_out-1, k:k,itile)=v(1:nx_out-1, 1:ny_out-1, 1:nz)
            enddo
            deallocate(vxyz_out, u, v, ew1, ew2, es1, es2,                 &
                       edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n)


  enddo ! npz


  do k=1,npz
    do l=1,ntiles
       j1 = (npy_out)*(l-1) + 1
       j2 = (npy_out)*(l-1) + npy_out
       varo(:,j1:j2)=u_out(:,:,k,l)
    enddo
    write(OUNIT) varo
  enddo
  do k=1,npz
    do l=1,ntiles
       j1 = (npy_out)*(l-1) + 1
       j2 = (npy_out)*(l-1) + npy_out
       varo(:,j1:j2)=v_out(:,:,k,l)
    enddo
    write(OUNIT) varo
  enddo

! print*, u_in(:,1,1,1)
! print*, '---'
! print*, u_out(:,1,1,1)

  deallocate ( vari, varo, u_in, v_in, u_out, v_out )
  deallocate ( xyz_corner_in, xyz_corner_out )

 end subroutine read_interp_write_flow


  Function REVERSE(A) Result(B)
    real(ESMF_KIND_R8), Intent(In) :: A(:,:)
    real(ESMF_KIND_R8) :: B(Size(A,1),Size(A,2))

    Integer :: i, n

    n = Size(A, 1)

    Do i = 1, n
       B(i,:) = A(1+n-i,:)
    End Do

  End Function REVERSE

  Subroutine Write_Profile(arr_global, npx, npy, npz, name)
    integer,            intent(IN) :: npx, npy, npz
    real(ESMF_KIND_R8), intent(IN) :: arr_global(npx,npy,npz)
    character(len=*),   intent(IN) :: name

    integer :: k
    real(ESMF_KIND_R8) :: rng(3,npz)
    real(ESMF_KIND_R8) :: GSUM

    IF (gid==masterproc) Then
       rng(1,:) = MINVAL(MINVAL(arr_global,DIM=1),DIM=1)
       rng(2,:) = MAXVAL(MAXVAL(arr_global,DIM=1),DIM=1)
       rng(3,:) = SUM(SUM(arr_global,DIM=1),DIM=1)/(npx*npy)
       GSUM     = SUM(SUM(SUM(arr_global,DIM=1),DIM=1),DIM=1)

       print*,'***********'
       print*,'stats for ',trim(name)

       Do k = 1, npz
          Write(*,'(a,i4.0,3(f21.9,1x))')'k:',k,rng(:,k)
       End Do
  !    Write(*,"('GlobalSum: ',f21.9)") GSUM
       print*,'***********'
       print*,' '
    End IF

  End Subroutine Write_Profile

  end

