      PROGRAM CreateTopography

      use ESMF
      use constants_mod,  only: pi, grav

! Shared Utilities
      use fms_mod,        only: fms_init, fms_end, file_exist
      use mpp_mod,        only: mpp_error, FATAL, NOTE
      use fv_arrays_mod,  only: fv_atmos_type, FVPRC, REAL4, REAL8
      use fv_control_mod, only: npx,npy,npz, ntiles
      use fv_control_mod, only: fv_init, fv_end
      use fv_mp_mod,      only: gid, masterproc, tile, mp_gather, mp_barrier
      use fv_grid_tools_mod,  only: grid, agrid, area, dx, dy, dxc, dyc

! Surf Map Utilities
      use fv_surf_map_mod,    only: surfdrv, map_to_cubed_simple

      IMPLICIT NONE

      integer:: im
      integer:: jm
      real(REAL8), allocatable :: phis_m(:,:)
      real(REAL8), allocatable :: r8tmp(:,:)

      real(REAL8), allocatable :: phis_global(:,:,:)
      real(REAL4), allocatable :: rtopo(:,:)

      real(REAL8), allocatable :: gwd_global(:,:,:)
      real(REAL4), allocatable :: rgwd(:,:)

      real(REAL8), allocatable :: trb_global(:,:,:)
      real(REAL4), allocatable :: rtrb(:,:)

      character*125 :: fname1
      character*125 :: fname2
      character*125 :: fname3

      real(REAL4) :: rmax,rmin
      real(REAL8) :: dlon, dlat
      integer :: nlon, nlat
      integer :: c2c_interp_npts
 
      real(REAL8), allocatable :: lat1(:)
      real(REAL8), allocatable :: lon1(:)
      real(REAL8), allocatable :: r8latlon1(:,:)
      real(REAL8), allocatable :: r8latlon2(:,:)
      real(REAL8) :: dt

      integer :: npts

      integer :: i,j,itile, j1,j2

      character(30) :: str_arg
#ifndef __GFORTRAN__
      external :: getarg, iargc
      integer iargc
#endif

      type(fv_atmos_type) :: Atm(1)

! Start up FMS/MPP
      call fms_init()

      if (IARGC() /= 5) then
        call mpp_error(FATAL, 'ABORT: need 4 argument npx,npy, nlon,nlat, c2c_interp_npts')
      endif
      CALL GETARG(1, str_arg)
      read (str_arg,'(I10)') im
      CALL GETARG(2, str_arg)
      read (str_arg,'(I10)') jm
      jm = jm*6
      CALL GETARG(3, str_arg)
      read (str_arg,'(I10)') nlon
      CALL GETARG(4, str_arg)
      read (str_arg,'(I10)') nlat
      CALL GETARG(5, str_arg)
      read (str_arg,'(I10)') c2c_interp_npts

! allocate lat/lon arrays
      allocate ( lat1(nlat+1) )
      allocate ( lon1(nlon+1) )
      allocate ( r8latlon1(nlon,nlat) )
      allocate ( r8latlon2(nlon,nlat) )

! Get Topography
      npx = IM+1
      npy = (JM/6) + 1
      npz = 32
      ntiles = 6
      dt = 1800.0
      call fv_init(Atm, dt)

    if (npx-1 <= 1440) then
      allocate ( phis_m(Atm(1)%isd:Atm(1)%ied,Atm(1)%jsd:Atm(1)%jed) )
      allocate ( phis_global(npx-1,npy-1,6) )
      allocate ( rtopo(npx-1,(npy-1)*6) )
      call surfdrv(npx, npy, grid, agrid, area, dx, dy, dxc, dyc,  &
                   phis_m, gid==masterproc)
      phis_global(Atm(1)%isc:Atm(1)%iec,Atm(1)%jsc:Atm(1)%jec,tile) = &
           phis_m(Atm(1)%isc:Atm(1)%iec,Atm(1)%jsc:Atm(1)%jec)/grav
      call mp_gather(phis_global, Atm(1)%isc, Atm(1)%iec, Atm(1)%jsc, Atm(1)%jec, npx-1, npy-1, ntiles)
      if (gid==masterproc) then
        do itile=1,ntiles
           j1 = (npy-1)*(itile-1) + 1
           j2 = (npy-1)*(itile-1) + npy-1
           rtopo(:,j1:j2) = phis_global(1:npx-1,1:npy-1,itile)
        enddo
        write(fname1, "('topo_DYN_ave_',i4.4,'x',i5.5,'.data')") IM,JM
        open(44,file=fname1,form='unformatted',status='unknown')
        write(44) rtopo
        close(44)
      endif
      deallocate( rtopo )
      deallocate( phis_global )
      deallocate( phis_m )
    endif

! Get GWD/TRB Variance 
  ! if (npx-1 <= 180) then
      allocate( r8tmp(Atm(1)%isd:Atm(1)%ied,Atm(1)%jsd:Atm(1)%jed) )
      allocate ( gwd_global(im,im,6) )
      allocate ( rgwd(im,jm) )
! Init latlon grid
      dlon=(pi+pi)/real(nlon)
      dlat=pi/real(nlat)
      do i=1,nlon+1
         lon1(i)=real(i-1)*dlon
      enddo
      lat1(1)   =-0.5*pi
      lat1(nlat+1)= 0.5*pi
      do j=2,nlat
         lat1(j)=lat1(1)+(real(j)-1.)*dlat
      enddo
! Read GWD
      write(fname1, "('topo_GWD_var_',i3.3,'x',i2.2,'_DC.data')") nlon,nlat
      if (.not. file_exist(fname1)) then
         write(fname1, "('topo_GWD_var_',i3.3,'x',i3.3,'_DC.data')") nlon,nlat
         if (.not. file_exist(fname1)) then
            write(fname1, "('topo_GWD_var_',i4.4,'x',i3.3,'_DC.data')") nlon,nlat
            if (.not. file_exist(fname1)) then
               write(fname1, "('topo_GWD_var_',i4.4,'x',i4.4,'_DC.data')") nlon,nlat
               if (.not. file_exist(fname1)) call mpp_error(FATAL,fname1)
            endif
         endif
      endif
      open(40,file=fname1,form='unformatted',status='old')
      read(40) r8latlon1
      close(40)
      if (gid==masterproc) then
          rmax =  vmax(r8latlon1,rmin,nlon,nlat,1)
          write(6,*) 'gwd_max=', rmax*1.E3
          write(6,*) 'gwd_min=', rmin*1.E3
      endif
      call mp_barrier()
! Regrid from -180:180 to 0:360
      r8latlon2(1         :nlon/2,:) = r8latlon1(nlon/2 + 1 :nlon  , :)
      r8latlon2(nlon/2 + 1:nlon  ,:) = r8latlon1(1          :nlon/2, :)
! Interp GWD
      call map_to_cubed_simple(nlon, nlat, lat1, lon1, r8latlon2, grid, agrid, r8tmp, npx, npy)
       gwd_global(Atm(1)%isc:Atm(1)%iec,Atm(1)%jsc:Atm(1)%jec,tile) = &
           r8tmp(Atm(1)%isc:Atm(1)%iec,Atm(1)%jsc:Atm(1)%jec)**2.0
      call mp_gather(gwd_global, Atm(1)%isc, Atm(1)%iec, Atm(1)%jsc, Atm(1)%jec, npx-1, npy-1, ntiles)
      if (gid==masterproc) then
        do itile=1,ntiles
           j1 = (npy-1)*(itile-1) + 1
           j2 = (npy-1)*(itile-1) + npy-1
            rgwd(:,j1:j2) =  gwd_global(1:npx-1,1:npy-1,itile)
        enddo
        write(fname2, "('topo_GWD_var_',i4.4,'x',i5.5,'.data')") IM,JM
        open(45,file=fname2,form='unformatted',status='unknown')
        write(45) rgwd
        close(45)
      endif
      deallocate( rgwd )
      deallocate( gwd_global )

      allocate ( trb_global(im,im,6) )
      allocate ( rtrb(im,jm) )
! Read TRB
      write(fname1, "('topo_TRB_var_',i3.3,'x',i2.2,'_DC.data')") nlon,nlat
      if (.not. file_exist(fname1)) then
         write(fname1, "('topo_TRB_var_',i3.3,'x',i3.3,'_DC.data')") nlon,nlat
         if (.not. file_exist(fname1)) then
            write(fname1, "('topo_TRB_var_',i4.4,'x',i3.3,'_DC.data')") nlon,nlat
            if (.not. file_exist(fname1)) then
               write(fname1, "('topo_TRB_var_',i4.4,'x',i4.4,'_DC.data')") nlon,nlat
               if (.not. file_exist(fname1)) call mpp_error(FATAL,fname1)
            endif
         endif
      endif
      open(40,file=fname1,form='unformatted',status='old')
      read(40) r8latlon1
      close(40)
      if (gid==masterproc) then
          rmax =  vmax(r8latlon1,rmin,nlon,nlat,1)
          write(6,*) 'trb_max=', rmax*1.E3
          write(6,*) 'trb_min=', rmin*1.E3
      endif
      call mp_barrier()
! Regrid from -180:180 to 0:360
      r8latlon2(1         :nlon/2,:) = r8latlon1(nlon/2 + 1 :nlon  , :)
      r8latlon2(nlon/2 + 1:nlon  ,:) = r8latlon1(1          :nlon/2, :)
! Interp TRB
      call map_to_cubed_simple(nlon, nlat, lat1, lon1, r8latlon2, grid, agrid, r8tmp, npx, npy)
       trb_global(Atm(1)%isc:Atm(1)%iec,Atm(1)%jsc:Atm(1)%jec,tile) = &
           r8tmp(Atm(1)%isc:Atm(1)%iec,Atm(1)%jsc:Atm(1)%jec)
      call mp_gather(trb_global, Atm(1)%isc, Atm(1)%iec, Atm(1)%jsc, Atm(1)%jec, npx-1, npy-1, ntiles)
      if (gid==masterproc) then
        do itile=1,ntiles
           j1 = (npy-1)*(itile-1) + 1
           j2 = (npy-1)*(itile-1) + npy-1
            rtrb(:,j1:j2) =  trb_global(1:npx-1,1:npy-1,itile)
        enddo
        write(fname3, "('topo_TRB_var_',i4.4,'x',i5.5,'.data')") IM,JM
        open(46,file=fname3,form='unformatted',status='unknown')
        write(46) rtrb
        close(46)
      endif
      deallocate( rtrb )
      deallocate( trb_global )
      deallocate( r8tmp )
  ! else
  !   call cube2cube_gwd_trb(c2c_interp_npts,npx-1)
  ! endif

    deallocate( lat1 )
    deallocate( lon1 )
    deallocate( r8latlon1 )
    deallocate( r8latlon2 )

      call fv_end(Atm)
      call fms_end()

    contains

#ifdef DO_NOT_SKIP
 subroutine cube2cube_gwd_trb(npx_in,npx_out)
! Cube to Cube Utilities
  use CUB2CUB_mod,       only : get_c2c_weight, do_c2c_interpolation
  use fv_grid_utils_mod, only : gnomonic_grids
  use fv_grid_tools_mod, only : mirror_grid
  use GHOST_CUBSPH_mod,  only : B_grid, A_grid, ghost_cubsph_update
      IMPLICIT NONE
  integer, intent(in) :: npx_in
  integer, intent(in) :: npx_out

  integer :: npts
  integer :: ndims=2
  real(ESMF_KIND_R8), allocatable :: xs(:,:), ys(:,:)
  real(ESMF_KIND_R8), allocatable :: grid_in(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: grid_out(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: corner_in(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: corner_out(:,:,:,:)
  real(ESMF_KIND_R8), allocatable :: weight_c2c(:,:,:,:)
  integer,            allocatable :: index_c2c(:,:,:,:)

  real(ESMF_KIND_R4), allocatable :: vari(:,:)
  real(ESMF_KIND_R4), allocatable :: varo(:,:)
  real(ESMF_KIND_R8), allocatable :: var_in(:,:,:)
  real(ESMF_KIND_R8), allocatable :: var_out(:,:,:)

  integer :: grid_type = 0
  integer :: i,j,l,j1,j2,status

  integer :: IUNIT=15
  integer :: OUNIT=16

  character*125 :: fname1
  character*125 :: fname2

  if (gid==masterproc) then

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
  do l=1,ntiles
     do j=1,npts
        do i=1,npts
!---------------------------------
! Shift the corner away from Japan
!---------------------------------
! This will result in the corner close to east coast of China
           grid_in(i,j,1,l) = grid_in(i,j,1,l) - pi/18.
           if ( grid_in(i,j,1,l) < 0. )              &
                grid_in(i,j,1,l) = grid_in(i,j,1,l) + 2.*pi
           if (ABS(grid_in(i,j,1,l)) < 1.e-10) grid_in(i,j,1,l) = 0.0
           if (ABS(grid_in(i,j,2,l)) < 1.e-10) grid_in(i,j,2,l) = 0.0
        enddo
     enddo
  enddo
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
                               ntiles,  1, l, B_grid)
     call ghost_cubsph_update(corner_in(2,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1, 1, &
                               ntiles,  1, l, B_grid)
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
  do l=1,ntiles
     do j=1,npts
        do i=1,npts
!---------------------------------
! Shift the corner away from Japan
!---------------------------------
! This will result in the corner close to east coast of China
           grid_out(i,j,1,l) = grid_out(i,j,1,l) - pi/18.
           if ( grid_out(i,j,1,l) < 0. )              &
                grid_out(i,j,1,l) = grid_out(i,j,1,l) + 2.*pi
           if (ABS(grid_out(i,j,1,l)) < 1.e-10) grid_out(i,j,1,l) = 0.0
           if (ABS(grid_out(i,j,2,l)) < 1.e-10) grid_out(i,j,2,l) = 0.0
        enddo
     enddo
  enddo
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
     call ghost_cubsph_update(corner_out(1,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1,  &
                              1, ntiles,  1, l, B_grid)
     call ghost_cubsph_update(corner_out(2,0:npts+1,0:npts+1,:), 0, npts+1, 0, npts+1,  &
                              1, ntiles,  1, l, B_grid)
  enddo
  deallocate (grid_out)

  !--------------------------------------------------------------------!
  ! calculate weights and indices from bilinear interpolation          !
  ! from grid_in to grid_out                                           !
  !--------------------------------------------------------------------!
  allocate(index_c2c (3, npx_out, npx_out, ntiles))
  allocate(weight_c2c(4, npx_out, npx_out, ntiles))
  call get_c2c_weight(ntiles, npx_in+1, npx_in+1, corner_in,                &
                      npx_out+1, npx_out+1, 1,npx_out, 1,npx_out, 1,ntiles, &
                      corner_out, index_c2c,  weight_c2c)

  allocate( vari(npx_in,ntiles*(npx_in)) )
  allocate( varo(npx_out,ntiles*(npx_out)) )
  allocate( var_in(0:npx_in+1,0:npx_in+1,ntiles) )

! Do DYN_ave interp
  allocate( var_out(npx_out,npx_out,ntiles) )
  write(fname1, "('topo_DYN_ave_',i4.4,'x',i5.5,'.data')") (npx_in),ntiles*(npx_in)
  if (.not. file_exist(fname1)) call mpp_error(FATAL,fname1)
  open(IUNIT,file=fname1,form='unformatted',status='old')
  read(IUNIT) vari
  close(IUNIT)
  do l=1,ntiles
     j1 = (npx_in)*(l-1) + 1
     j2 = (npx_in)*(l-1) + npx_in
     var_in(1:npx_in,1:npx_in,l)=vari(:,j1:j2)
  enddo
  do l=1,ntiles
    call ghost_cubsph_update(var_in, 0, npx_in+1, 0, npx_in+1,  1, ntiles,  &
                              1, l, A_grid)
  enddo
  call do_c2c_interpolation(var_in, 0, npx_in+1, 0, npx_in+1, 1, ntiles, &
                            index_c2c, weight_c2c, npx_out, npx_out, var_out)
  write(fname2, "('topo_DYN_ave_',i4.4,'x',i5.5,'_interp.data')") (npx_out),ntiles*(npx_out)
  open(OUNIT,file=fname2,form='unformatted',status='unknown')
  do l=1,ntiles
     j1 = (npx_out)*(l-1) + 1
     j2 = (npx_out)*(l-1) + npx_out
     varo(:,j1:j2)=var_out(:,:,l)
  enddo
  write(OUNIT) varo
  close(OUNIT)

! Do GWD interp
  if (npx_out-1 <= 360) then
  write(fname1, "('topo_GWD_var_',i4.4,'x',i5.5,'.data')") (npx_in),ntiles*(npx_in)
  if (.not. file_exist(fname1)) call mpp_error(FATAL,fname1)
  open(IUNIT,file=fname1,form='unformatted',status='old')
  read(IUNIT) vari
  close(IUNIT)
  do l=1,ntiles
     j1 = (npx_in)*(l-1) + 1
     j2 = (npx_in)*(l-1) + npx_in
     var_in(1:npx_in,1:npx_in,l)=vari(:,j1:j2)
  enddo
  do l=1,ntiles
    call ghost_cubsph_update(var_in, 0, npx_in+1, 0, npx_in+1,  1, ntiles,  &
                              1, l, A_grid)
  enddo
  call do_c2c_interpolation(var_in, 0, npx_in+1, 0, npx_in+1, 1, ntiles, &
                            index_c2c, weight_c2c, npx_out, npx_out, var_out)
  else
    var_out(:,:,:) = 0.0
  endif
  write(fname2, "('topo_GWD_var_',i4.4,'x',i5.5,'.data')") (npx_out),ntiles*(npx_out)
  open(OUNIT,file=fname2,form='unformatted',status='unknown')
  do l=1,ntiles
     j1 = (npx_out)*(l-1) + 1
     j2 = (npx_out)*(l-1) + npx_out
     varo(:,j1:j2)=var_out(:,:,l)
  enddo
  write(OUNIT) varo
  close(OUNIT)

! Do TRB interp
  if (npx_out-1 <= 360) then
  write(fname1, "('topo_TRB_var_',i4.4,'x',i5.5,'.data')") (npx_in),ntiles*(npx_in)
  if (.not. file_exist(fname1)) call mpp_error(FATAL,fname1)
  open(IUNIT,file=fname1,form='unformatted',status='old')
  read(IUNIT) vari
  close(IUNIT)
  do l=1,ntiles
     j1 = (npx_in)*(l-1) + 1
     j2 = (npx_in)*(l-1) + npx_in
     var_in(1:npx_in,1:npx_in,l)=vari(:,j1:j2)
  enddo
  do l=1,ntiles
    call ghost_cubsph_update(var_in, 0, npx_in+1, 0, npx_in+1,  1, ntiles,  &
                              1, l, A_grid)
  enddo
  call do_c2c_interpolation(var_in, 0, npx_in+1, 0, npx_in+1, 1, ntiles, &
                            index_c2c, weight_c2c, npx_out, npx_out, var_out)
  else
    var_out(:,:,:) = 0.0
  endif
  write(fname2, "('topo_TRB_var_',i4.4,'x',i5.5,'.data')") (npx_out),ntiles*(npx_out) 
  open(OUNIT,file=fname2,form='unformatted',status='unknown')
  do l=1,ntiles
     j1 = (npx_out)*(l-1) + 1 
     j2 = (npx_out)*(l-1) + npx_out
     varo(:,j1:j2)=var_out(:,:,l)
  enddo
  write(OUNIT) varo
  close(OUNIT)

  deallocate( var_out )
  deallocate( var_in )
  deallocate( varo )
  deallocate( vari )

  deallocate( corner_out )
  deallocate( corner_in )
  deallocate( weight_c2c )
  deallocate( index_c2c )

  endif  ! masterproc

 end subroutine cube2cube_gwd_trb
#endif

 real function vmax(a,pmin,m,n,z)
      IMPLICIT NONE
      integer m,n,z, i,j,k
      real(REAL4) pmin, pmax
      real(REAL8) a(m,n,z)

      pmax = a(1,1,1)
      pmin = a(1,1,1)

      do k=1,z
      do j=1,n
      do i=1,m
         pmax = max(pmax,a(i,j,k))
         pmin = min(pmin,a(i,j,k))
      enddo
      enddo
      enddo

      vmax = pmax
 end function vmax

      END

