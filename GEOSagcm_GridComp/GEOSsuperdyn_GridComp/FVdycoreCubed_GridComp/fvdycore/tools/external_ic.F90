module external_ic_mod

#ifdef MAPL_MODE
#define DEALLOCGLOB_(A) if(associated(A)) then;A=0;if(MAPL_ShmInitialized) then; call MAPL_DeAllocNodeArray(A,rc=status);else; deallocate(A);endif;NULLIFY(A);endif
#endif

#ifndef DYCORE_SOLO
   use amip_interp_mod,    only: i_sst, j_sst, sst_ncep
#endif
   use fv_arrays_mod,      only: REAL4, REAL8, FVPRC, R_GRID
   use fms_mod,            only: file_exist, read_data, field_exist
   use fms_io_mod,         only: get_tile_string, field_size
   use mpp_mod,            only: mpp_error, FATAL, NOTE, mpp_broadcast,mpp_npes
   use mpp_parameter_mod,  only: AGRID_PARAM=>AGRID
   use mpp_domains_mod,    only: mpp_get_tile_id, domain2d, mpp_update_domains, mpp_get_boundary, DGRID_NE
   use tracer_manager_mod, only: get_tracer_names, get_number_tracers, get_tracer_index
   use field_manager_mod,  only: MODEL_ATMOS

   use constants_mod,     only: pi=>pi_8, omega, grav, kappa, rdgas, rvgas, cp_air
#ifdef MAPL_MODE
   use MAPL_IOMod
   use MAPL_ShmemMod
   use pFIO
#endif
   use, intrinsic :: iso_fortran_env, only: REAL64, REAL32

   use fv_arrays_mod,     only: fv_atmos_type, fv_grid_type, fv_grid_bounds_type
   use fv_diagnostics_mod,only: prt_maxmin
   use fv_mp_mod,         only: is_master, ng, mp_barrier, mp_gather, mp_bcst, &
         is,js,ie,je, isd,jsd,ied,jed, fill_corners, YDir
   use fv_grid_utils_mod, only: ptop_min
   use fv_grid_utils_mod, only: normalize_vect
   use fv_mapz_mod,       only: mappm
   use fv_surf_map_mod,   only: surfdrv
   use fv_timing_mod,     only: timing_on, timing_off
   use init_hydro_mod,    only: p_var
#ifndef MAPL_MODE
   use fv_io_mod,         only: fv_io_read_tracers
   use sim_nc_mod,        only: open_ncfile, close_ncfile, get_ncdim1, get_var1_double, get_var2_real,   &
         get_var3_r4
#endif
   use fv_timing_mod,      only: timing_on, timing_off
   use mpp_mod,           only: mpp_pe
   use memutils_mod, only: print_memuse_stats
   use fv_nwp_nudge_mod,  only: T_is_Tv

   implicit none

#include "mpif.h"

   private

   real(FVPRC), parameter:: zvir = rvgas/rdgas - 1.
   real(FVPRC) :: deg2rad

   public get_external_ic,get_cubed_sphere_terrain

   integer, SAVE :: tile, npes_x, npes_y

   integer :: status
   integer :: IUNIT=15
   integer :: OUNIT=16

contains

   subroutine get_external_ic( Atm, fv_domain, use_geos_latlon_restart, use_geos_cubed_restart, ntracers )

      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      logical, optional,   intent(in)    :: use_geos_latlon_restart
      logical, optional,   intent(in)    :: use_geos_cubed_restart
      integer, optional,   intent(in)    :: ntracers(3)
      real(FVPRC):: alpha = 0.
      integer i,j,k,nq
      logical :: cubed_sphere,fv_diag_ic

! * Initialize coriolis param:

      do j=jsd,jed+1
         do i=isd,ied+1
            Atm(1)%gridstruct%fc(i,j) = 2.*omega*( -1.*cos(Atm(1)%gridstruct%grid(i,j,1))*cos(Atm(1)%gridstruct%grid(i,j,2))*sin(alpha) + &
                  sin(Atm(1)%gridstruct%grid(i,j,2))*cos(alpha) )
         enddo
      enddo

      do j=jsd,jed
         do i=isd,ied
            Atm(1)%gridstruct%f0(i,j) = 2.*omega*( -1.*cos(Atm(1)%gridstruct%agrid(i,j,1))*cos(Atm(1)%gridstruct%agrid(i,j,2))*sin(alpha) + &
                  sin(Atm(1)%gridstruct%agrid(i,j,2))*cos(alpha) )
         enddo
      enddo

      cubed_sphere=.true.
      fv_diag_ic=.false.
      call mpp_update_domains( Atm(1)%gridstruct%f0, Atm(1)%domain )
      if ( cubed_sphere ) call fill_corners(Atm(1)%gridstruct%f0, Atm(1)%npx, Atm(1)%npy, YDir)

! Read in cubed_sphere terrain
      if ( Atm(1)%flagstruct%mountain ) then
         call get_cubed_sphere_terrain(Atm, fv_domain)
      else
         Atm(1)%phis = 0.
      endif

! Read in the specified external dataset and do all the needed transformation
      if ( present(use_geos_latlon_restart)) then
         if (allocated(Atm(1)%q)) deallocate( Atm(1)%q )
         allocate  ( Atm(1)%q(isd:ied,jsd:jed,Atm(1)%npz,Atm(1)%ncnst) )
         call get_geos_latlon_ic( Atm, fv_domain, Atm(1)%ncnst, ntracers )
      elseif ( present(use_geos_cubed_restart)) then
         if (allocated(Atm(1)%q)) deallocate( Atm(1)%q )
         allocate  ( Atm(1)%q(isd:ied,jsd:jed,Atm(1)%npz,Atm(1)%ncnst) )
         call get_geos_cubed_ic( Atm, fv_domain, Atm(1)%ncnst, ntracers )
      else
         if ( Atm(1)%flagstruct%ncep_ic ) then
            nq = 1
            call timing_on('NCEP_IC')
            call get_ncep_ic( Atm, fv_domain, nq )
            call timing_off('NCEP_IC')
#ifndef MAPL_MODE          
#ifndef NO_FV_TRACERS
            call fv_io_read_tracers( fv_domain, Atm )
            if(is_master()) write(6,*) 'All tracers except sphum replaced by FV IC'
#endif
#endif
         elseif ( fv_diag_ic ) then
! Interpolate/remap diagnostic output from a FV model diagnostic output file on uniform lat-lon A grid:
            nq = size(Atm(1)%q,4)
! Needed variables: lon, lat, pfull(dim), zsurf, ps, ucomp, vcomp, w, temp, and all q
! delz not implemnetd yet; set make_nh = .true.
            call get_diag_ic( Atm, fv_domain, nq )
         else
! The following is to read in lagacy lat-lon FV core restart file
!  is Atm%q defined in all cases?

            nq = size(Atm(1)%q,4)
            call get_fv_ic( Atm, fv_domain, nq )
         endif
      endif

      call prt_maxmin('T', Atm(1)%pt, is, ie, js, je, ng, Atm(1)%npz, 1.0_FVPRC)

      call p_var(Atm(1)%npz,  is, ie, js, je, Atm(1)%ak(1),  ptop_min,         &
            Atm(1)%delp, Atm(1)%delz, Atm(1)%pt, Atm(1)%ps,               &
            Atm(1)%pe,   Atm(1)%peln, Atm(1)%pk, Atm(1)%pkz,              &
            kappa, Atm(1)%q, ng, Atm(1)%ncnst, dble(Atm(1)%gridstruct%area),Atm(1)%flagstruct%dry_mass,           &
            Atm(1)%flagstruct%adjust_dry_mass, Atm(1)%flagstruct%mountain, Atm(1)%flagstruct%moist_phys,   &
            Atm(1)%flagstruct%hydrostatic, Atm(1)%flagstruct%nwat, Atm(1)%domain,Atm(1)%flagstruct%make_nh)

   end subroutine get_external_ic



   subroutine get_cubed_sphere_terrain( Atm, fv_domain )
      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      integer              :: ntileMe
      integer, allocatable :: tile_id(:)
      character(len=64)    :: fname
      integer              ::  n
      integer              ::  jbeg, jend
      real(FVPRC) ftop

#ifndef MAPL_MODE

      ntileMe = size(Atm(:))  ! This will have to be modified for mult tiles per PE
! always one at this point

      allocate( tile_id(ntileMe) )
      tile_id = mpp_get_tile_id( fv_domain )

      do n=1,ntileMe

         call get_tile_string(fname, 'INPUT/fv_core.res.tile', tile_id(n), '.nc' )

         if( file_exist(fname) ) then
            call read_data(fname, 'phis', Atm(n)%phis(is:ie,js:je),      &
                  domain=fv_domain, tile_count=n)
         else
            call surfdrv(  Atm(n)%npx, Atm(n)%npy, grid, agrid,   &
                  area, dx, dy, dxc, dyc, Atm(n)%phis, is_master())
            call mpp_error(NOTE,'terrain datasets generated using USGS data')
         endif

      end do

      call mpp_update_domains( Atm(1)%phis, domain )
      ftop = g_sum(Atm(1)%phis(is:ie,js:je), is, ie, js, je, ng, area, 1)

      call prt_maxmin('ZS', Atm(1)%phis,  is, ie, js, je, ng, 1, 1./grav)
      if(is_master()) write(6,*) 'mean terrain height (m)=', ftop/grav

      deallocate( tile_id )

#endif

   end subroutine get_cubed_sphere_terrain


   subroutine get_diag_ic( Atm, fv_domain, nq )
      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      integer, intent(in):: nq
! local:
      character(len=128) :: fname, tracer_name
      real(kind=4), allocatable:: wk1(:), wk2(:,:), wk3(:,:,:)
      real(FVPRC), allocatable:: tp(:,:,:), qp(:,:,:,:)
      real(FVPRC), allocatable:: ua(:,:,:), va(:,:,:), wa(:,:,:)
      real(FVPRC), allocatable:: lat(:), lon(:), ak0(:), bk0(:)
      real(FVPRC):: s2c(is:ie,js:je,4)
      integer, dimension(is:ie,js:je):: id1, id2, jdc
      real(FVPRC) psc(is:ie,js:je)
      real(FVPRC) gzc(is:ie,js:je)
      integer:: i, j, k, im, jm, km, npz, npt
      integer:: i1, i2, j1, ncid
      integer:: jbeg, jend
      integer tsize(4), tr_ind
      logical:: found

      integer  sphum, liq_wat, rainwat, ice_wat, snowwat, graupel
#ifndef MAPL_MODE

      deg2rad = pi/180.

      npz = Atm(1)%npz

! Zero out all initial tracer fields:
      Atm(1)%q = 0.

      fname = Atm(1)%res_latlon_dynamics

      if( file_exist(fname) ) then
         call open_ncfile( fname, ncid )        ! open the file
         call get_ncdim1( ncid, 'lon',   tsize(1) )
         call get_ncdim1( ncid, 'lat',   tsize(2) )
         call get_ncdim1( ncid, 'pfull', tsize(3) )

         im = tsize(1); jm = tsize(2); km = tsize(3)

         if(is_master())  write(*,*) fname, ' FV_diag IC dimensions:', tsize

         allocate (  lon(im) )
         allocate (  lat(jm) )

         call get_var1_double (ncid, 'lon', im, lon )
         call get_var1_double (ncid, 'lat', jm, lat )

! Convert to radian
         do i=1,im
            lon(i) = lon(i) * deg2rad  ! lon(1) = 0.
         enddo
         do j=1,jm
            lat(j) = lat(j) * deg2rad
         enddo

         allocate ( ak0(km+1) )
         allocate ( bk0(km+1) )
! npz
         if ( npz /= km ) then
            call mpp_error(FATAL,'==>Error in get_diag_ic: vertical dim must be the same')
         else
            ak0(:) = Atm(1)%ak(:)
            bk0(:) = Atm(1)%bk(:)
         endif
      else
         call mpp_error(FATAL,'==> Error in get_diag_ic: Expected file '//trim(fname)//' does not exist')
      endif

! Initialize lat-lon to Cubed bi-linear interpolation coeff:
      call remap_coef( im, jm, lon, lat, id1, id2, jdc, s2c, Atm(1)%gridstruct%agrid, Atm(1)%bd )

! Find bounding latitudes:
      jbeg = jm-1;         jend = 2
      do j=js,je
         do i=is,ie
            j1 = jdc(i,j)
            jbeg = min(jbeg, j1) 
            jend = max(jend, j1+1)
         enddo
      enddo

! remap surface pressure and height:
      allocate ( wk2(im,jbeg:jend) )
      call get_var3_r4( ncid, 'ps', 1,im, jbeg,jend, 1,1, wk2 )
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            psc(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                  s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
         enddo
      enddo

      call get_var3_r4( ncid, 'zsurf', 1,im, jbeg,jend, 1,1, wk2 )
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            gzc(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                  s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
         enddo
      enddo
      deallocate ( wk2 )

! Read in temperature:
      allocate ( wk3(1:im,jbeg:jend, 1:km) )
      call get_var3_r4( ncid, 'temp', 1,im, jbeg,jend, 1,km, wk3 )
      allocate (  tp(is:ie,js:je,km) )
      do k=1,km
         do j=js,je
            do i=is,ie
               i1 = id1(i,j)
               i2 = id2(i,j)
               j1 = jdc(i,j)
               tp(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                     s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
            enddo
         enddo
      enddo

! Read in all tracers:
      allocate ( qp(is:ie,js:je,km,nq) )
      qp = 0.
      do tr_ind=1, nq
         call get_tracer_names(MODEL_ATMOS, tr_ind, tracer_name)
         if (field_exist(fname,tracer_name)) then
            call get_var3_r4( ncid, tracer_name, 1,im, jbeg,jend, 1,km, wk3 )
            do k=1,km
               do j=js,je
                  do i=is,ie
                     i1 = id1(i,j)
                     i2 = id2(i,j)
                     j1 = jdc(i,j)
                     qp(i,j,k,tr_ind) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                           s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
                  enddo
               enddo
            enddo
         endif
      enddo
      call remap_scalar(im, jm, km, npz, nq, nq, ak0, bk0, psc, gzc, tp, qp, Atm(1))
      deallocate ( tp )
      deallocate ( qp )

! Winds:
      call get_var3_r4( ncid, 'ucomp', 1,im, jbeg,jend, 1,km, wk3 )
      allocate ( ua(is:ie,js:je,km) )
      do k=1,km
         do j=js,je
            do i=is,ie
               i1 = id1(i,j)
               i2 = id2(i,j)
               j1 = jdc(i,j)
               ua(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                     s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
            enddo
         enddo
      enddo

      call get_var3_r4( ncid, 'vcomp', 1,im, jbeg,jend, 1,km, wk3 )
      allocate ( va(is:ie,js:je,km) )
      do k=1,km
         do j=js,je
            do i=is,ie
               i1 = id1(i,j)
               i2 = id2(i,j)
               j1 = jdc(i,j)
               va(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                     s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
            enddo
         enddo
      enddo
      call remap_winds(im, jm, km, npz, ak0, bk0, psc, ua, va, Atm)

      deallocate ( ua )
      deallocate ( va )

      if ( .not. Atm(1)%hydrostatic ) then
         if (field_exist(fname,'w')) then
            allocate ( wa(is:ie,js:je,km) )
            call get_var3_r4( ncid, 'w', 1,im, jbeg,jend, 1,km, wk3 )
            do k=1,km
               do j=js,je
                  do i=is,ie
                     i1 = id1(i,j)
                     i2 = id2(i,j)
                     j1 = jdc(i,j)
                     wa(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                           s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
                  enddo
               enddo
            enddo
            call remap_wz(im, jm, km, npz, ng, ak0, bk0, psc, wa, Atm(1)%w, Atm)
            deallocate ( wa )
         else    ! use "w = - fac * omega" ?
            Atm(1)%w(:,:,:) = 0.
         endif
! delz:
         if (field_exist(fname,'delz')) then
            allocate ( wa(is:ie,js:je,km) )
            call get_var3_r4( ncid, 'delz', 1,im, jbeg,jend, 1,km, wk3 )
            do k=1,km
               do j=js,je
                  do i=is,ie
                     i1 = id1(i,j)
                     i2 = id2(i,j)
                     j1 = jdc(i,j)
                     wa(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                           s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
                  enddo
               enddo
            enddo
            call remap_wz(im, jm, km, npz, 0,  ak0, bk0, psc, wa, Atm(1)%delz, Atm)
            deallocate ( wa )
         else    ! Force make = T
            Atm(1)%make_nh = .true.
         endif

      endif   ! hydrostatic test

      call close_ncfile ( ncid )
      deallocate ( wk3 )
      deallocate ( ak0 )
      deallocate ( bk0 )
      deallocate ( lat )
      deallocate ( lon )
#endif

   end subroutine get_diag_ic

   subroutine get_geos_cubed_ic( Atm, fv_domain, nq, ntracers )
      use GHOST_CUBSPH_mod,  only : A_grid, ghost_cubsph_update
      use CUB2CUB_mod,    only: get_c2c_weight,                 &
            interpolate_c2c
      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      integer, intent(in):: nq, ntracers(3)

      character(len=128) :: fname, fname1
      real(FVPRC), allocatable:: pkz0(:,:)
      real(FVPRC), allocatable:: ps0(:,:), gz0(:,:), t0(:,:,:), q0(:,:,:)
      real(FVPRC), allocatable:: u0(:,:,:), v0(:,:,:)
      real(FVPRC), allocatable:: ak0(:), bk0(:)
      integer :: i, j, k, l, iq, j1, j2, im, jm, km, npts, npx, npy, npz
      integer :: iq_moist0 , iq_moist1
      integer :: iq_gocart0, iq_gocart1
      integer :: iq_pchem0 , iq_pchem1
      integer :: ntiles=6
      logical found
      integer :: header(6)
      character (len=8) :: imc, jmc

      real(FVPRC) :: qtmp

      integer:: i1, i2
      integer:: ic, jc, jc0
      real(FVPRC) psc(is:ie,js:je)
      real(FVPRC) gzc(is:ie,js:je)
      real(FVPRC), allocatable:: tp(:,:,:), qp(:,:,:,:)
      real(FVPRC), allocatable:: ua(:,:,:), va(:,:,:)

      real(REAL64), allocatable :: akbk_r8(:)

      real(REAL64), dimension(:,:,:,:), allocatable :: corner_in, corner_out, weight_c2c
      integer, dimension(:,:,:,:), allocatable :: index_c2c

      integer :: is_i,ie_i, js_i,je_i
      integer :: isd_i,ied_i, jsd_i,jed_i
      integer :: ng_i = 5
      type(domain2d) :: domain_i

      real(FVPRC), allocatable :: ebuffer(:,:)
      real(FVPRC), allocatable :: nbuffer(:,:)
      real(FVPRC), allocatable :: wbuffer(:,:)
      real(FVPRC), allocatable :: sbuffer(:,:)

      integer (kind=MPI_OFFSET_KIND) :: slice_2d
      integer (kind=MPI_OFFSET_KIND) :: offset
      character(len=64) :: strTxt
      integer :: nmoist, ngocart, npchem, nqmap

      integer            :: filetype
      logical            :: isNC4
      type(Netcdf4_Fileformatter) :: formatter
      type(FileMetadata), allocatable :: cfg(:)
      type(StringIntegerMap) :: dims
      integer            :: nDims, nVars, ivar, dimSizes(3)
      character(len=128) :: vname
      real(FVPRC),   allocatable  :: gslice_r4(:,:)
      real*8, allocatable  :: gslice_r8(:,:)
      integer            :: tileoff,lvar_cnt
      character(len=128), allocatable :: vnames(:)
      character(len=:), pointer :: var_name
      type(StringVariableMap), pointer :: vars
      type(StringVariableMapIterator) :: iter


!bma added
      character(len=128) :: moist_order(9) = (/"Q   ","QLLS","QLCN","CLLS","CLCN","QILS","QICN","NCPL","NCPI"/)

      iq_moist0 =           1
      iq_moist1 =           ntracers(1)
      iq_gocart0=iq_moist1 +1
      iq_gocart1=iq_moist1 +ntracers(2)
      iq_pchem0 =iq_gocart1+1
      iq_pchem1 =iq_gocart1+ntracers(3)
      nmoist  = ntracers(1)
      ngocart = ntracers(2)
      npchem  = ntracers(3)

      npx = Atm(1)%npx
      npy = Atm(1)%npy
      npz = Atm(1)%npz

! Zero out all initial tracer fields:
      Atm(1)%q = 0.

! Read input FV core restart file
      fname = "fvcore_internal_restart_in"

      if( file_exist(fname) ) then

         call MAPL_NCIOGetFileType(fname,filetype)
         if (filetype >=0 ) then
            isNC4 = .true.
         else
            isNC4 = .false.
         end if

         if (isNC4) then

            allocate(cfg(1))
            call formatter%open(fname,pFIO_READ,rc=status)
            cfg(1) = formatter%read(rc=status)
            im =cfg(1)%get_dimension('lon',rc=status)
            jm =cfg(1)%get_dimension('lat',rc=status)
            km =cfg(1)%get_dimension('lev',rc=status)
             
            allocate(gslice_r8(im,jm))

         else

            open(IUNIT,file=fname ,access='sequential',form='unformatted',status='old')
            read (IUNIT, IOSTAT=status) header
            if (is_master()) print*, header
            read (IUNIT, IOSTAT=status) header(1:5)
            if (is_master()) print*, header(1:5)

            im=header(1)
            jm=header(2)
            km=header(3)

         end if

         if(is_master()) write(*,*) 'Using GEOS restart:', fname

         if ( file_exist(fname) ) then
            if(is_master())  write(*,*) 'External IC dimensions:', im   , jm       , km
            if(is_master())  write(*,*) 'Interpolating to      :', npx-1, (npy-1)*6, npz
         else
            call mpp_error(FATAL,'==> Error from get_external_ic:        &
                  & field not found')
         endif

!--------------------------------------------------------------------!
! setup input cubed-sphere domain                                    !
!--------------------------------------------------------------------!
         call mpp_domain_decomp(domain_i,im+1,(jm/ntiles)+1,ntiles,ng_i,ntiles, &
               is_i,ie_i,js_i,je_i,isd_i,ied_i,jsd_i,jed_i,tile)
!--------------------------------------------------------------------!
! initialize cubed sphere grid: in                                   !
!--------------------------------------------------------------------!
         allocate(corner_in(2,is_i:ie_i+1,js_i:je_i+1,tile:tile))
         call init_cubsph_grid(im+1, is_i,ie_i, js_i,je_i, ntiles, corner_in)
         call print_memuse_stats('get_geos_cubed_ic: init corner_in')
!--------------------------------------------------------------------!
! initialize cubed sphere grid: out                                  !
!--------------------------------------------------------------------!
         allocate(corner_out(2,is:ie+1,js:je+1,tile:tile))
         do l=tile,tile
            do j=js,je+1
               do i=is,ie+1
                  corner_out(1,i,j,l) = Atm(1)%gridstruct%grid(i,j,1)
                  corner_out(2,i,j,l) = Atm(1)%gridstruct%grid(i,j,2)
               enddo
            enddo
         enddo
         call print_memuse_stats('get_geos_cubed_ic: init corner_out')
!--------------------------------------------------------------------!
! calculate weights and indices from bilinear interpolation          !
! from grid_in to grid_out                                           !
!--------------------------------------------------------------------!
         allocate(index_c2c (3, is:ie, js:je, tile:tile))
         allocate(weight_c2c(4, is:ie, js:je, tile:tile))
         call get_c2c_weight(ntiles, im+1, (jm/6)+1, &
               is_i, ie_i, js_i, je_i, isd_i, ied_i, jsd_i, jed_i, &
               corner_in(:,is_i:ie_i+1,js_i:je_i+1,tile), &
               npx, npy, is,ie, js,je, tile,tile, &
               corner_out(:,is:ie+1,js:je+1,tile:tile), &
               index_c2c,  weight_c2c, domain_i)
         npts = im
         call print_memuse_stats('get_geos_cubed_ic: get_c2c_weight')

         allocate ( ak0(km+1) )
         allocate ( bk0(km+1) )
         allocate ( akbk_r8(km+1) )
         if (isNC4) then
            call MAPL_VarRead(formatter,"AK",akbk_r8)
         else
            read (IUNIT, IOSTAT=status) akbk_r8
         end if
         ak0 = akbk_r8
         if (isNC4) then
            call MAPL_VarRead(formatter,"BK",akbk_r8)
         else
            read (IUNIT, IOSTAT=status) akbk_r8
         end if
         bk0 = akbk_r8
         deallocate ( akbk_r8 )
         call print_memuse_stats('get_geos_cubed_ic: read ak/bk')
         close (IUNIT)

! Read U
         allocate (  u0(isd_i:ied_i,jsd_i:jed_i+1,km) )
         u0(:,:,:) = 0.0 
         if (isNC4) then
            tileoff = (tile-1)*(jm/ntiles)
            do k=1,km
               call MAPL_VarRead(formatter,"U",gslice_r8,lev=k)
               u0(is_i:ie_i,js_i:je_i,k) = gslice_r8(is_i:ie+i,tileoff+js_i:tileoff+je_i)
            enddo
         else
!offset = sequential access: 4 + INT(6) + 8 + INT(5) + 8 + DBL(NPZ+1) + 8 + DBL(NPZ+1) + 8
            offset =                    4 + 24     + 8 + 20     + 8 + (km+1)*8   + 8 + (km+1)*8   + 8
            if (is_master()) print*, offset
            call parallel_read_file_r8(fname, npts, is_i,ie_i, js_i,je_i, km, offset, u0(is_i:ie_i,js_i:je_i,:))
         end if
         call print_memuse_stats('get_geos_cubed_ic: read U')
! Read V
         allocate (  v0(isd_i:ied_i+1,jsd_i:jed_i,km) )
         v0(:,:,:) = 0.0
         if (isNC4) then
            tileoff = (tile-1)*(jm/ntiles)
            do k=1,km
               call MAPL_VarRead(formatter,"V",gslice_r8,lev=k)
               v0(is_i:ie_i,js_i:je_i,k) = gslice_r8(is_i:ie+i,tileoff+js_i:tileoff+je_i)
            enddo
         else
            if (is_master()) print*, offset
            call parallel_read_file_r8(fname, npts, is_i,ie_i, js_i,je_i, km, offset, v0(is_i:ie_i,js_i:je_i,:))
         end if
         call print_memuse_stats('get_geos_cubed_ic: read V')
         allocate ( sbuffer(is_i:ie_i,km) )
         allocate ( wbuffer(js_i:je_i,km) )
         allocate ( nbuffer(is_i:ie_i,km) )
         allocate ( ebuffer(js_i:je_i,km) )
         call mpp_get_boundary(u0, v0, domain_i, &
               wbuffery=wbuffer, ebuffery=ebuffer, &
               sbufferx=sbuffer, nbufferx=nbuffer, &
               gridtype=DGRID_NE )
         do k=1,km
            do i=is_i,ie_i    
               u0(i,je_i+1,k) = nbuffer(i,k)
            enddo
            do j=js_i,je_i
               v0(ie_i+1,j,k) = ebuffer(j,k)
            enddo
         enddo
         deallocate ( sbuffer )
         deallocate ( wbuffer )
         deallocate ( nbuffer )
         deallocate ( ebuffer )
         call mpp_update_domains( u0, v0, domain_i, gridtype=DGRID_NE, complete=.true. )
         call prt_maxmin(' U_geos', u0, is_i, ie_i, js_i, je_i, ng_i, km, 1.0_FVPRC)
         call prt_maxmin(' V_geos', v0, is_i, ie_i, js_i, je_i, ng_i, km, 1.0_FVPRC)
         allocate ( ua(is:ie,js:je,km) )
         allocate ( va(is:ie,js:je,km) )
         call interp_c2c_vect(npts, npts, npx-1, npy-1, km, ntiles, domain_i, &
               is,ie, js,je, isd_i,ied_i, jsd_i,jed_i, is_i,ie_i, js_i,je_i, &
               u0, v0, ua, va, index_c2c, weight_c2c, corner_in(:,is_i:ie_i+1,js_i:je_i+1,tile), corner_out, Atm(1)%gridstruct)
         deallocate ( v0 )
         deallocate ( u0 )
         deallocate ( corner_in )
         deallocate ( corner_out )
! Read T
         allocate (  t0(isd_i:ied_i,jsd_i:jed_i,km) )
         t0(:,:,:) = 0.0
         if (isNC4) then
            tileoff = (tile-1)*(jm/ntiles)
            do k=1,km
               call MAPL_VarRead(formatter,"PT",gslice_r8,lev=k)
               t0(is_i:ie_i,js_i:je_i,k) = gslice_r8(is_i:ie+i,tileoff+js_i:tileoff+je_i)
            enddo
         else
            if (is_master()) print*, offset
            call parallel_read_file_r8(fname, npts, is_i,ie_i, js_i,je_i, km, offset, t0(is_i:ie_i,js_i:je_i,:))
         end if
         call print_memuse_stats('get_geos_cubed_ic: read T')
! Read PE at Surface only
         allocate ( ps0(isd_i:ied_i,jsd_i:jed_i) )
         ps0(:,:) = 0.0
         if (isNC4) then
            tileoff = (tile-1)*(jm/ntiles)
            call MAPL_VarRead(formatter,"PE",gslice_r8,lev=km+1)
            ps0(is_i:ie_i,js_i:je_i) = gslice_r8(is_i:ie+i,tileoff+js_i:tileoff+je_i)
         else
            slice_2d = npts*npts*ntiles
            do k=1,km
               offset = offset + slice_2d*8 + 8  ! skip first KM levels of Edge Pressure to find surface pressure
            enddo
            if (is_master()) print*, offset, slice_2d
            call parallel_read_file_r8(fname, npts, is_i,ie_i, js_i,je_i, 1, offset, ps0(is_i:ie_i,js_i:je_i))
         end if
         call mpp_update_domains(ps0, domain_i)
! Read PKZ
         allocate ( pkz0(isd_i:ied_i,jsd_i:jed_i) )
         pkz0(:,:) = 0.0
         if (isNC4) then
            tileoff = (tile-1)*(jm/ntiles)
            do k=1,km
               call MAPL_VarRead(formatter,"PKZ",gslice_r8,lev=k)
               pkz0(is_i:ie_i,js_i:je_i) = gslice_r8(is_i:ie+i,tileoff+js_i:tileoff+je_i)
               t0(is_i:ie_i,js_i:je_i,k) = t0(is_i:ie_i,js_i:je_i,k)*pkz0(is_i:ie_i,js_i:je_i)
            enddo
         else
            if (is_master()) print*, offset
            do k=1,km
               call parallel_read_file_r8(fname, npts, is_i,ie_i, js_i,je_i, 1, offset, pkz0(is_i:ie_i,js_i:je_i))
! t0 needs to be just temperature with no virtual effect
               t0(is_i:ie_i,js_i:je_i,k) = t0(is_i:ie_i,js_i:je_i,k)*pkz0(is_i:ie_i,js_i:je_i)
            enddo
         end if
         call print_memuse_stats('get_geos_cubed_ic: converted T')
         deallocate ( pkz0 )

         if (isNC4) then
            call formatter%close()
            deallocate(cfg)
            deallocate(gslice_r8)
         end if

         allocate ( gz0(isd_i:ied_i,jsd_i:jed_i) )
         gz0(:,:) = 0.0

         write(imc, "(i8)") im
         write(jmc, "(i8)") jm
         imc = adjustl(imc)
         jmc = adjustl(jmc)

         write(fname1, "('topo_DYN_ave_',a,'x',a,'.data')") trim(imc), trim(jmc)
         if (.not. file_exist(fname1)) then
            call mpp_error(FATAL,'get_geos_cubed_ic: cannot find topo_DYN_ave file')
         endif
         call print_memuse_stats('get_geos_cubed_ic: '//TRIM(fname1)//' being read')
         offset = 4
         call parallel_read_file_r4(fname1, npts, is_i,ie_i, js_i,je_i, 1, offset, gz0(is_i:ie_i,js_i:je_i))
         call mpp_update_domains(gz0, domain_i)
         gz0 = gz0*grav

! Read cubed-sphere phis from file since IMPORT is not ready yet
         offset = 4
         call parallel_read_file_r4('topo_dynave.data', Atm(1)%npx-1, is,ie, js,je, 1, offset, Atm(1)%phis(is:ie,js:je))
         call mpp_update_domains(Atm(1)%phis, Atm(1)%domain)
         Atm(1)%phis = Atm(1)%phis*grav
         call print_memuse_stats('get_geos_cubed_ic: phis')

! Horiz Interp for surface pressure 
         call prt_maxmin('PS_geos', ps0, is_i, ie_i, js_i, je_i, ng_i, 1, 1.0_FVPRC)
         do j=js,je
            do i=is,ie
               ic=index_c2c(1,i,j,tile)
               jc=index_c2c(2,i,j,tile)
               psc(i,j)=weight_c2c(1,i,j,tile)*ps0(ic  ,jc  )  &      
                     +weight_c2c(2,i,j,tile)*ps0(ic  ,jc+1)  &
                     +weight_c2c(3,i,j,tile)*ps0(ic+1,jc+1)  &
                     +weight_c2c(4,i,j,tile)*ps0(ic+1,jc  )
            enddo
         enddo
         deallocate ( ps0 )
! Horiz Interp for surface height
         call prt_maxmin('GZ_geos', gz0, is_i, ie_i, js_i, je_i, ng_i, 1, 1.0/grav)
         do j=js,je
            do i=is,ie
               ic=index_c2c(1,i,j,tile)
               jc=index_c2c(2,i,j,tile)
               gzc(i,j)=weight_c2c(1,i,j,tile)*gz0(ic  ,jc  )  &       
                     +weight_c2c(2,i,j,tile)*gz0(ic  ,jc+1)  &
                     +weight_c2c(3,i,j,tile)*gz0(ic+1,jc+1)  &
                     +weight_c2c(4,i,j,tile)*gz0(ic+1,jc  )
            enddo
         enddo
         deallocate ( gz0 )

! Horiz Interp for Q
         allocate ( q0(isd_i:ied_i,jsd_i:jed_i,km) )
         allocate ( qp(is:ie,js:je,km,nq) )
         q0(:,:,:) = 0.0
         qp(:,:,:,:) = 0.0

! Horiz Interp for moist tracers
! is there a moist restart file to interpolate?
! Read in tracers: only sphum at this point
         if( file_exist("moist_internal_restart_in") .and. ntracers(1) > 0 ) then
            if (is_master()) print*, 'Trying to interpolate moist_internal_restart_in'

            call MAPL_NCIOGetFileType("moist_internal_restart_in",filetype)

            if (filetype /= 0) then
               offset=4
            else
               lvar_cnt = 0
               allocate(gslice_r4(im,jm))
               allocate(cfg(1))
               call formatter%open("moist_internal_restart_in",pFIO_READ,rc=status)
               cfg(1) = formatter%read(rc=status)
               call MAPL_IOCountNonDimVars(cfg(1),nvars,rc=status)
               if (nVars /= iq_moist1-iq_moist0+1) call mpp_error(FATAL,'Wrong number of variables in moist file') 
               tileoff = (tile-1)*(jm/ntiles)
            end if

            do ivar=iq_moist0,iq_moist1
               if (filetype ==0) lvar_cnt=lvar_cnt+1
               do k=1,km
                  if (filetype /= 0) then
                     call parallel_read_file_r4('moist_internal_restart_in', npts, is_i,ie_i, js_i,je_i, 1, offset, q0(is_i:ie_i,js_i:je_i,k))
                     call mpp_update_domains(q0(:,:,k), domain_i)
                  else
                     vname = trim(moist_order(lvar_cnt))
                     call MAPL_VarRead(formatter,vname,gslice_r4,lev=k)
                     q0(is_i:ie_i,js_i:je_i,k)=gslice_r4(is_i:ie+i,tileoff+js_i:tileoff+je_i)
                  end if
                  call mpp_update_domains(q0(:,:,k), domain_i)
                  do j=js,je
                     do i=is,ie
                        ic=index_c2c(1,i,j,tile)
                        jc=index_c2c(2,i,j,tile)
                        qp(i,j,k,ivar)=weight_c2c(1,i,j,tile)*q0(ic  ,jc  ,k)  &
                              +weight_c2c(2,i,j,tile)*q0(ic  ,jc+1,k)  &
                              +weight_c2c(3,i,j,tile)*q0(ic+1,jc+1,k)  &
                              +weight_c2c(4,i,j,tile)*q0(ic+1,jc  ,k)
                     enddo
                  enddo

               enddo
          !WMP if (ivar == 1) t0(is_i:ie_i,js_i:je_i,:) = (t0(is_i:ie_i,js_i:je_i,:)/(1.0 + zvir*q0(is_i:ie_i,js_i:je_i,:)))
               call prt_maxmin( 'Q_geos_moist', q0, is_i, ie_i, js_i, je_i, ng_i, km, 1._FVPRC)
            enddo

            if (filetype == 0) then
               call formatter%close()
               deallocate(cfg)
               deallocate(gslice_r4)
            end if

         end if

! Horiz Interp for GOCART tracers
! is there a gocart restart file to interpolate?
! Read in tracers: only sphum at this point
         if( file_exist("gocart_internal_restart_in") .and. ntracers(2) > 0 ) then
            if (is_master()) print*, 'Trying to interpolate gocart_internal_restart_in'

            call MAPL_NCIOGetFileType("gocart_internal_restart_in",filetype)

            if (filetype /= 0) then
               offset=4
            else
               lvar_cnt = 0 
               allocate(gslice_r4(im,jm))
               allocate(cfg(1))
               call formatter%open("gocart_internal_restart_in",pFIO_READ,rc=status)
               cfg(1) = formatter%read(rc=status)
               call MAPL_IOCountNonDimVars(cfg(1),nvars,rc=status)
               if (nVars /= iq_gocart1-iq_gocart0+1) call mpp_error(FATAL,'Wrong number of variables in gocart file') 
               tileoff = (tile-1)*(jm/ntiles)

               allocate(vnames(nVars))
               vars => cfg(1)%get_variables()
               iter = vars%begin()
        
               lvar_cnt=0 
               do while(iter /= vars%end())
                  var_name => iter%key()
                  if (.not.cfg(1)%is_coordinate_variable(var_name)) then
                     lvar_cnt=lvar_cnt+1
                     vnames(lvar_cnt)=var_name
                  end if
                  call iter%next()
               enddo
            end if

            lvar_cnt = 0
            do ivar=iq_gocart0,iq_gocart1
               if (filetype ==0) lvar_cnt=lvar_cnt+1
               do k=1,km
                  if (filetype /= 0) then
                     call parallel_read_file_r4('gocart_internal_restart_in', npts, is_i,ie_i, js_i,je_i, 1, offset, q0(is_i:ie_i,js_i:je_i,k))
                     call mpp_update_domains(q0(:,:,k), domain_i)
                  else
                     call MAPL_VarRead(formatter,vnames(lvar_cnt),gslice_r4,lev=k)
                     q0(is_i:ie_i,js_i:je_i,k)=gslice_r4(is_i:ie+i,tileoff+js_i:tileoff+je_i)
                  end if
                  call mpp_update_domains(q0(:,:,k), domain_i)
                  do j=js,je
                     do i=is,ie
                        ic=index_c2c(1,i,j,tile)
                        jc=index_c2c(2,i,j,tile)
                        qp(i,j,k,ivar)=weight_c2c(1,i,j,tile)*q0(ic  ,jc  ,k)  &
                              +weight_c2c(2,i,j,tile)*q0(ic  ,jc+1,k)  &
                              +weight_c2c(3,i,j,tile)*q0(ic+1,jc+1,k)  &
                              +weight_c2c(4,i,j,tile)*q0(ic+1,jc  ,k)
                     enddo
                  enddo

               enddo
               call prt_maxmin( 'Q_geos_gocart', q0, is_i, ie_i, js_i, je_i, ng_i, km, 1._FVPRC)
            enddo

            if (filetype == 0) then
               call formatter%close()
               deallocate(cfg)
               deallocate(vnames)
               deallocate(gslice_r4)
            end if

         end if
! Horiz Interp for pchem tracers
! is there a gocart restart file to interpolate?
! Read in tracers: only sphum at this point
         if( file_exist("pchem_internal_restart_in") .and. ntracers(3) > 0 ) then
            if (is_master()) print*, 'Trying to interpolate pchem_internal_restart_in'

            call MAPL_NCIOGetFileType("pchem_internal_restart_in",filetype)

            if (filetype /= 0) then
               offset=4
            else
               lvar_cnt = 0 
               allocate(gslice_r4(im,jm))
               allocate(cfg(1))
               call formatter%open("pchem_internal_restart_in",pFIO_READ,rc=status)
               cfg(1) = formatter%read(rc=status)
               call MAPL_IOCountNonDimVars(cfg(1),nvars,rc=status)
               if (nVars /= iq_pchem1-iq_pchem0+1) call mpp_error(FATAL,'Wrong number of variables in pchem file') 
               tileoff = (tile-1)*(jm/ntiles)

               allocate(vnames(nVars))
               vars => cfg(1)%get_variables()
               iter = vars%begin()
        
               lvar_cnt=0 
               do while(iter /= vars%end())
                  var_name => iter%key()
                  if (.not.cfg(1)%is_coordinate_variable(var_name)) then
                     lvar_cnt=lvar_cnt+1
                     vnames(lvar_cnt)=var_name
                  end if
                  call iter%next()
               enddo
            end if

            lvar_cnt=0
            do ivar=iq_pchem0,iq_pchem1
               if (filetype == 0) lvar_cnt=lvar_cnt+1
               do k=1,km
                  if (filetype /= 0) then
                     call parallel_read_file_r4('pchem_internal_restart_in', npts, is_i,ie_i, js_i,je_i, 1, offset, q0(is_i:ie_i,js_i:je_i,k))
                     call mpp_update_domains(q0(:,:,k), domain_i)
                  else
                     call MAPL_VarRead(formatter,vnames(lvar_cnt),gslice_r4,lev=k)
                     q0(is_i:ie_i,js_i:je_i,k)=gslice_r4(is_i:ie+i,tileoff+js_i:tileoff+je_i)
                  end if
                  call mpp_update_domains(q0(:,:,k), domain_i)
                  do j=js,je
                     do i=is,ie
                        ic=index_c2c(1,i,j,tile)
                        jc=index_c2c(2,i,j,tile)
                        qp(i,j,k,ivar)=weight_c2c(1,i,j,tile)*q0(ic  ,jc  ,k)  &
                              +weight_c2c(2,i,j,tile)*q0(ic  ,jc+1,k)  &
                              +weight_c2c(3,i,j,tile)*q0(ic+1,jc+1,k)  &
                              +weight_c2c(4,i,j,tile)*q0(ic+1,jc  ,k)
                     enddo
                  enddo

               enddo
               call prt_maxmin( 'Q_geos_pchem', q0, is_i, ie_i, js_i, je_i, ng_i, km, 1._FVPRC)
            enddo

            if (filetype == 0) then
               call formatter%close()
               deallocate(cfg)
               deallocate(vnames)
               deallocate(gslice_r4)
            end if

         end if
                   
! Horiz Interp for T
         deallocate ( q0 )
         call mpp_update_domains(t0, domain_i)
         call prt_maxmin( 'T_geos', t0, is_i, ie_i, js_i, je_i, ng_i, km, 1.0_FVPRC)
         allocate (  tp(is:ie,js:je,km) )
         do k=1,km
            do j=js,je
               do i=is,ie
                  ic=index_c2c(1,i,j,tile)
                  jc=index_c2c(2,i,j,tile)
                  tp(i,j,k)=weight_c2c(1,i,j,tile)*t0(ic  ,jc  ,k)  &
                        +weight_c2c(2,i,j,tile)*t0(ic  ,jc+1,k)  &
                        +weight_c2c(3,i,j,tile)*t0(ic+1,jc+1,k)  &
                        +weight_c2c(4,i,j,tile)*t0(ic+1,jc  ,k)
               enddo
            enddo
         enddo
         deallocate ( t0 )
         deallocate( index_c2c )
         deallocate( weight_c2c )

! Horz/Vert remap for scalars
         nqmap =  nmoist + ngocart + npchem 

         call remap_scalar(im, jm, km, npz, nqmap, nqmap, ak0, bk0, psc, gzc, tp, qp, Atm(1))

         deallocate ( tp )
         deallocate ( qp )
         call print_memuse_stats('get_geos_cubed_ic: remap_scalar')
! Horz/Vert remap for U/V
         call remap_winds(im, jm, km, npz, ak0, bk0, psc, ua, va, Atm)
         deallocate ( ua )
         deallocate ( va )
         call print_memuse_stats('get_geos_cubed_ic: remap_winds')

      else
         call mpp_error(FATAL,'==> Error from get_external_ic:        &
               & Expected file '//trim(fname)//' does not exist')
      endif

      if (allocated(bk0)) deallocate ( bk0 )
      if (allocated(ak0)) deallocate ( ak0 )

! Finished, let's check the results !

      call prt_maxmin('GZ_model', Atm(1)%phis, is, ie, js, je, ng, 1, 1.0/grav)
      call prt_maxmin('PS_model', Atm(1)%ps, is, ie, js, je, ng, 1, 0.01_FVPRC)
      call prt_maxmin('DP_model', Atm(1)%delp, is, ie, js, je, ng, npz, 1.0_FVPRC)
      call prt_maxmin(' U_model', Atm(1)%u, is, ie, js, je, ng, npz, 1.0_FVPRC)
      call prt_maxmin(' V_model', Atm(1)%v, is, ie, js, je, ng, npz, 1.0_FVPRC)
      call prt_maxmin('PT_model', Atm(1)%pt, is, ie, js, je, ng, npz, 1.0_FVPRC)
! Range check the MOIST tracers
      do iq=iq_moist0,iq_moist1
         do k=1,npz
            do j=js,je
               do i=is,ie
                  Atm(1)%q(i,j,k,iq) = MIN(Atm(1)%q(i,j,k,iq),1.d0)
                  Atm(1)%q(i,j,k,iq) = MAX(Atm(1)%q(i,j,k,iq),0.d0)
               enddo
            enddo
         enddo
      enddo
      do iq=1,nmoist+npchem+ngocart
         call prt_maxmin('QP_model', Atm(1)%q(is:ie,js:je,1:npz,iq), is, ie, js, je, 0, npz, 1._FVPRC)
      enddo

   end subroutine get_geos_cubed_ic

   subroutine get_geos_latlon_ic( Atm, fv_domain, nq, ntracers)
      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      integer, intent(in):: nq,ntracers(3)

      character(len=128) :: fname, fname1
      real(FVPRC), allocatable:: pkz0(:,:)
      real(FVPRC), allocatable:: ps0(:,:), gz0(:,:), t0(:,:,:), q0(:,:,:)
      real(FVPRC), allocatable:: u0(:,:,:), v0(:,:,:), ua0(:,:), va0(:,:)
      real(FVPRC), allocatable:: lat(:), lon(:), ak0(:), bk0(:)
      integer :: i, j, k, l, iq, j1, j2, im, jm, km, npz
      logical found
      real(FVPRC) dak, dbk
      integer :: header(6)
      character (len=8) :: imc, jmc

      integer:: i1, i2, nmoist, ngocart, npchem, nqmap, offset
      real(FVPRC):: s2c(is:ie,js:je,4)
      integer, dimension(is:ie,js:je):: id1, id2, jdc
      real(FVPRC) psc(is:ie,js:je)
      real(FVPRC) gzc(is:ie,js:je)          
      real(FVPRC), allocatable:: tp(:,:,:), qp(:,:,:,:)
      real(FVPRC), allocatable:: ua(:,:,:), va(:,:,:)

      real(REAL4), allocatable :: phis_r4(:,:)
      real(REAL64), allocatable :: r8latlon(:,:)
      real(REAL4), allocatable :: r4latlon(:,:)
      real(REAL64), allocatable :: akbk_r8(:)

      character(len=128) :: tag
      character(len=3)   :: cid

      integer            :: filetype
      logical            :: isNC4
      integer            :: nDims, nVars, ivar, dimSizes(3)
      character(len=128) :: vname
      integer :: iq_moist0 , iq_moist1
      integer :: iq_gocart0, iq_gocart1
      integer :: iq_pchem0 , iq_pchem1
      integer :: lvar_cnt
      type(Netcdf4_Fileformatter) :: formatter
      type(FileMetadata), allocatable :: cfg(:)
      character(len=128), allocatable :: vnames(:)
      character(len=:), pointer :: var_name
      type(StringVariableMap), pointer :: vars
      type(StringVariableMapIterator) :: iter


!bma added
      character(len=128) :: moist_order(9) = (/"Q   ","QLLS","QLCN","CLLS","CLCN","QILS","QICN","NCPL","NCPI"/)

      iq_moist0 =           1
      iq_moist1 =           ntracers(1)
      iq_gocart0=iq_moist1 +1
      iq_gocart1=iq_moist1 +ntracers(2)
      iq_pchem0 =iq_gocart1+1
      iq_pchem1 =iq_gocart1+ntracers(3)
      nmoist  = ntracers(1)
      ngocart = ntracers(2)
      npchem  = ntracers(3)

      npz = Atm(1)%npz

! Zero out all initial tracer fields:
      Atm(1)%q = 0.

! Read in lat-lon FV core restart file
      fname = "fvcore_internal_restart_in"

      if( file_exist(fname) ) then


         call MAPL_NCIOGetFileType(fname,filetype)
         if (filetype >=0 ) then
            isNC4 = .true.
         else
            isNC4 = .false.
         end if

         if (isNC4) then

            allocate(cfg(1))
            call formatter%open(fname,pFIO_READ,rc=status)
            cfg(1) = formatter%read(rc=status)
            im =cfg(1)%get_dimension('lon',rc=status)
            jm =cfg(1)%get_dimension('lat',rc=status)
            km =cfg(1)%get_dimension('lev',rc=status)

         else

            open(IUNIT,file=fname ,access='sequential',form='unformatted',status='old')
            read (IUNIT, IOSTAT=status) header
            if (is_master()) print*, header
            read (IUNIT, IOSTAT=status) header(1:5)
            if (is_master()) print*, header(1:5)

            im=header(1)
            jm=header(2)
            km=header(3)

         end if

         if(is_master()) write(*,*) 'Using GEOS restart:', fname
         if(is_master())  write(*,*) 'External IC dimensions:', im, jm, km

         allocate (  lon(im) )
         do i=1,im
            lon(i) = (0.5 + real(i-1)) * 2.*pi/real(im)
         enddo
         allocate (  lat(jm) )
         do j=1,jm
            lat(j) = -0.5*pi + real(j-1)*pi/real(jm-1)   ! SP to NP 
         enddo

         call remap_coef( im, jm, lon, lat, id1, id2, jdc, s2c , Atm(1)%gridstruct%agrid, Atm(1)%bd)

         allocate ( ak0(km+1) )
         allocate ( bk0(km+1) )
         allocate ( akbk_r8(km+1) )
         if (isNC4) then
            call MAPL_VarRead(formatter,"AK",akbk_r8)
         else
            read (IUNIT, IOSTAT=status) akbk_r8
         end if
         ak0 = akbk_r8
         if (isNC4) then
            call MAPL_VarRead(formatter,"BK",akbk_r8)
         else
            read (IUNIT, IOSTAT=status) akbk_r8
         end if
         bk0 = akbk_r8
         deallocate ( akbk_r8 )

         call print_memuse_stats('get_geos_latlon_ic: read ak/bk')

         allocate ( r8latlon(im,jm) )
! Read U
         allocate (  u0(im,jm,km) )
         do k=1,km
            if (isNC4) then
               call MAPL_VarRead(formatter,"U",r8latlon,lev=k)
            else
               read (IUNIT, IOSTAT=status) r8latlon
            end if
! Regrid from -180:180 to 0:360
            u0(1       :im/2,:,k) = r8latlon(im/2 + 1 :im  , :)
            u0(im/2 + 1:im  ,:,k) = r8latlon(1        :im/2, :)
         enddo
         call print_memuse_stats('get_geos_latlon_ic: read u')
! Read V
         allocate (  v0(im,jm,km) )
         do k=1,km
            if (isNC4) then
               call MAPL_VarRead(formatter,"V",r8latlon,lev=k)
            else
               read (IUNIT, IOSTAT=status) r8latlon
            end if
! Regrid from -180:180 to 0:360
            v0(1       :im/2,:,k) = r8latlon(im/2 + 1 :im  , :)
            v0(im/2 + 1:im  ,:,k) = r8latlon(1        :im/2, :)
         enddo
         call print_memuse_stats('get_geos_latlon_ic: read v')
         if(is_master()) call pmaxmin( 'U_geos',   u0(:,2:jm,:), im*(jm-1), km, 1.0_FVPRC)
         if(is_master()) call pmaxmin( 'V_geos',   v0,           im*jm    , km, 1.0_FVPRC)
         allocate ( ua(is:ie,js:je,km) )
         allocate ( va(is:ie,js:je,km) )
         allocate ( ua0(im,jm) )
         allocate ( va0(im,jm) )
         do k=1,km
! Move latlon D winds to cell centers (A-grid)
            call d2a3d(u0(:,:,k), v0(:,:,k),  ua0(:,:),  va0(:,:), im, jm, 1, lon)
! Horiz Interp for U
            do j=js,je
               do i=is,ie
                  i1 = id1(i,j)
                  i2 = id2(i,j)
                  j1 = jdc(i,j)
                  ua(i,j,k) = s2c(i,j,1)*ua0(i1,j1  ) + s2c(i,j,2)*ua0(i2,j1  ) +  &
                        s2c(i,j,3)*ua0(i2,j1+1) + s2c(i,j,4)*ua0(i1,j1+1)
               enddo
            enddo
! Horiz Interp for V
            do j=js,je
               do i=is,ie
                  i1 = id1(i,j)
                  i2 = id2(i,j)
                  j1 = jdc(i,j)
                  va(i,j,k) = s2c(i,j,1)*va0(i1,j1  ) + s2c(i,j,2)*va0(i2,j1  ) +  &
                        s2c(i,j,3)*va0(i2,j1+1) + s2c(i,j,4)*va0(i1,j1+1)
               enddo
            enddo
         enddo
         call print_memuse_stats('get_geos_latlon_ic: d2a3d')
         deallocate ( v0 )
         deallocate ( u0 )
         deallocate ( ua0 )
         deallocate ( va0 )
! Read T
         allocate (  t0(im,jm,km) )
         do k=1,km
            if (isNC4) then
               call MAPL_VarRead(formatter,"PT",r8latlon,lev=k)
            else
               read (IUNIT, IOSTAT=status) r8latlon
            end if
! Regrid from -180:180 to 0:360
            t0(1       :im/2,:,k) = r8latlon(im/2 + 1 :im  , :)
            t0(im/2 + 1:im  ,:,k) = r8latlon(1        :im/2, :)
         enddo
         call print_memuse_stats('get_geos_latlon_ic: read t')
! Read PE
         do k=1,km+1 
            if (isNC4) then
               call MAPL_VarRead(formatter,"PE",r8latlon,lev=k)
            else
               read (IUNIT, IOSTAT=status) r8latlon
            end if
         enddo
! Regrid from -180:180 to 0:360
         allocate ( ps0(im,jm) )
         ps0(1       :im/2,:) = r8latlon(im/2 + 1 :im  , :)
         ps0(im/2 + 1:im  ,:) = r8latlon(1        :im/2, :)
         allocate ( pkz0(im,jm) )
         do k=1,km
            if (isNC4) then
               call MAPL_VarRead(formatter,"PKZ",r8latlon,lev=k)
            else
               read (IUNIT, IOSTAT=status) r8latlon
            end if
! Regrid from -180:180 to 0:360
            pkz0(1       :im/2,:) = r8latlon(im/2 + 1 :im  , :)
            pkz0(im/2 + 1:im  ,:) = r8latlon(1        :im/2, :)
! t0 needs to be just temperature with no virtual effect
            t0(:,:,k) = t0(:,:,k)*pkz0
         enddo
         deallocate ( r8latlon )
         call print_memuse_stats('get_geos_latlon_ic: converted T')
         deallocate ( pkz0 )
         if (isNC4) then
            call formatter%close()
            deallocate(cfg)
         else
            close (IUNIT)
         end if

         write(imc, "(i8)") im
         write(jmc, "(i8)") jm
         imc = adjustl(imc)
         jmc = adjustl(jmc)

         write(fname1, "('topo_DYN_ave_',a,'x',a,'_DC.data')") trim(imc), trim(jmc)
         if (.not. file_exist(fname1)) then
            CALL mpp_error(FATAL,'get_geos_latlon_ic: cannot find topo_DYN_ave file') 
         endif
         call print_memuse_stats('get_geos_latlon_ic: '//TRIM(fname1)//' being read')
         allocate ( r4latlon(im,jm) )
         open(IUNIT,file=fname1,form='unformatted',status='old')
         read(IUNIT) r4latlon
         close(IUNIT)
! Regrid from -180:180 to 0:360
         allocate ( gz0(im,jm) )
         gz0(1       :im/2,:) = r4latlon(im/2 + 1 :im  , :)
         gz0(im/2 + 1:im  ,:) = r4latlon(1        :im/2, :)
         gz0 = gz0*grav
         deallocate ( r4latlon )

! Read cubed-sphere phis from file since IMPORT is not ready yet
         write(imc, "(i8)")    Atm(1)%npx-1
         write(jmc, "(i8)") 6*(Atm(1)%npy-1)
         imc = adjustl(imc)
         jmc = adjustl(jmc)

         write(fname1, "('topo_DYN_ave_',a,'x',a,'.data')") trim(imc), trim(jmc)
         if (.not. file_exist(fname1)) then
            call mpp_error(FATAL,'get_geos_latlon_ic: cannot find topo_DYN_ave file')
         endif
         allocate( phis_r4(Atm(1)%npx-1,6*(Atm(1)%npy-1)) )
         open(IUNIT,file=fname1,form='unformatted',status='old')
         read(IUNIT) phis_r4
         close(IUNIT)
         Atm(1)%phis(is:ie,js:je) = phis_r4(is:ie,js+(tile-1)*(Atm(1)%npy-1):je+(tile-1)*(Atm(1)%npy-1))*grav
         call mpp_update_domains(Atm(1)%phis, Atm(1)%domain)
         deallocate( phis_r4 )
         call print_memuse_stats('get_geos_latlon_ic: phis')

! Horiz Interp for surface pressure 
         if(is_master()) call pmaxmin( 'PS_geos', ps0, im,    jm, 0.01_FVPRC)
         do j=js,je
            do i=is,ie
               i1 = id1(i,j)
               i2 = id2(i,j)
               j1 = jdc(i,j)
               psc(i,j) = s2c(i,j,1)*ps0(i1,j1  ) + s2c(i,j,2)*ps0(i2,j1  ) +  &
                     s2c(i,j,3)*ps0(i2,j1+1) + s2c(i,j,4)*ps0(i1,j1+1)
            enddo
         enddo
         deallocate ( ps0 )
! Horiz Interp for surface height
         if(is_master()) call pmaxmin( 'ZS_geos', gz0, im,    jm, 1.0/grav)
         do j=js,je
            do i=is,ie
               i1 = id1(i,j)
               i2 = id2(i,j)
               j1 = jdc(i,j)
               gzc(i,j) = s2c(i,j,1)*gz0(i1,j1  ) + s2c(i,j,2)*gz0(i2,j1  ) +  &
                     s2c(i,j,3)*gz0(i2,j1+1) + s2c(i,j,4)*gz0(i1,j1+1)
            enddo
         enddo
         deallocate ( gz0 )

! Horiz Interp for MOIST
! ----------------------
         allocate (  q0(im,jm,km) )
         allocate ( qp(is:ie,js:je,km,nq) )
         qp = 0.0

! Horiz Interp for moist tracers
! is there a moist restart file to interpolate?
! Read in tracers: only sphum at this point
         if( file_exist("moist_internal_restart_in") .and. ntracers(1) > 0 ) then
            if (is_master()) print*, 'Trying to interpolate moist_internal_restart_in'
            allocate ( r4latlon(im,jm) )

            call MAPL_NCIOGetFileType("moist_internal_restart_in",filetype)

            if (filetype /= 0) then
               open(IUNIT,file="moist_internal_restart_in" ,access='sequential',form='unformatted',status='old')
            else
               lvar_cnt = 0
               call formatter%open("moist_internal_restart_in",pFIO_READ,rc=status)
               cfg = formatter%read(rc=status)
               call MAPL_IOCountNonDimVars(cfg(1),nvars,rc=status)
               if (nVars /= iq_moist1-iq_moist0+1) call mpp_error(FATAL,'Wrong number of variables in moist file')
            end if

            do ivar=iq_moist0,iq_moist1
               if (filetype ==0) lvar_cnt=lvar_cnt+1
               do k=1,km
                  if (filetype /= 0) then
                     read (IUNIT, IOSTAT=status) r4latlon
                  else
                     vname = trim(moist_order(lvar_cnt))
                     call MAPL_VarRead(formatter,vname,r4latlon,lev=k)
                  end if
                  q0(1       :im/2,:,k) = r4latlon(im/2 + 1 :im  , :) ! Regrid from -180:180 to 0:360
                  q0(im/2 + 1:im  ,:,k) = r4latlon(1        :im/2, :) ! Regrid from -180:180 to 0:360
                  do j=js,je
                     do i=is,ie
                        i1 = id1(i,j)
                        i2 = id2(i,j)
                        j1 = jdc(i,j)
                        qp(i,j,k,ivar) = s2c(i,j,1)*q0(i1,j1  ,k) + s2c(i,j,2)*q0(i2,j1  ,k) +  &
                              s2c(i,j,3)*q0(i2,j1+1,k) + s2c(i,j,4)*q0(i1,j1+1,k)
                     enddo
                  enddo
               enddo
               if (ivar == 1) t0 = (t0/(1.0 + zvir*q0(:,:,:)))
               if (is_master()) call pmaxmin( 'MOIST_Q_',  q0(:,:,:), im*jm, km, 1.0_FVPRC)
            enddo
            if (filetype == 0) then
               call formatter%close()
               deallocate(cfg)
            else
               close(IUNIT)
            end if
            deallocate(r4latlon)

         end if

! Horiz Interp for GOCART tracers
! is there a gocart restart file to interpolate?
! Read in tracers: only sphum at this point
         if( file_exist("gocart_internal_restart_in") .and. ntracers(2) > 0 ) then
            if (is_master()) print*, 'Trying to interpolate gocart_internal_restart_in'
            allocate ( r4latlon(im,jm) )
            call MAPL_NCIOGetFileType("gocart_internal_restart_in",filetype)

            if (filetype /= 0) then
               open(IUNIT,file="gocart_internal_restart_in" ,access='sequential',form='unformatted',status='old')
            else
               lvar_cnt = 0
               allocate(cfg(1))
               call formatter%open("gocart_internal_restart_in",pFIO_READ,rc=status)
               cfg(1) = formatter%read(rc=status)
               call MAPL_IOCountNonDimVars(cfg(1),nvars,rc=status)
               if (nVars /= iq_gocart1-iq_gocart0+1) call mpp_error(FATAL,'Wrong number of variables in gocart file')

               allocate(vnames(nVars))
               vars => cfg(1)%get_variables()
               iter = vars%begin()

               do while(iter /= vars%end())
                  var_name => iter%key()
                  if (.not.cfg(1)%is_coordinate_variable(var_name)) then
                     lvar_cnt=lvar_cnt+1
                     vnames(lvar_cnt)=var_name
                  end if
               enddo

            end if

            lvar_cnt=0
            do ivar=iq_gocart0,iq_gocart1
               if (filetype ==0) lvar_cnt=lvar_cnt+1
               do k=1,km
                  if (filetype /= 0) then
                     read (IUNIT, IOSTAT=status) r4latlon
                  else
                     call MAPL_VarRead(formatter,vnames(lvar_cnt),r4latlon,lev=k)
                  end if
                  q0(1       :im/2,:,k) = r4latlon(im/2 + 1 :im  , :) ! Regrid from -180:180 to 0:360
                  q0(im/2 + 1:im  ,:,k) = r4latlon(1        :im/2, :) ! Regrid from -180:180 to 0:360
                  do j=js,je
                     do i=is,ie
                        i1 = id1(i,j)
                        i2 = id2(i,j)
                        j1 = jdc(i,j)
                        qp(i,j,k,ivar) = s2c(i,j,1)*q0(i1,j1  ,k) + s2c(i,j,2)*q0(i2,j1  ,k) +  &
                              s2c(i,j,3)*q0(i2,j1+1,k) + s2c(i,j,4)*q0(i1,j1+1,k)
                     enddo
                  enddo

               enddo
               if (is_master()) call pmaxmin( 'GOCART_Q_',  q0(:,:,:), im*jm, km, 1.0_FVPRC)
            enddo

            if (filetype == 0) then
               call formatter%close()
               deallocate(cfg)
               deallocate(vnames)
            else
               close(IUNIT)
            end if
            deallocate(r4latlon)
         end if

! Horiz Interp for pchem tracers
! is there a pchem restart file to interpolate?
! Read in tracers: only sphum at this point
         if( file_exist("pchem_internal_restart_in") .and. ntracers(3) > 0 ) then
            if (is_master()) print*, 'Trying to interpolate pchem_internal_restart_in'
            allocate ( r4latlon(im,jm) )
            call MAPL_NCIOGetFileType("pchem_internal_restart_in",filetype)

            if (filetype /= 0) then
               open(IUNIT,file="pchem_internal_restart_in" ,access='sequential',form='unformatted',status='old')
            else
               lvar_cnt = 0
               allocate(cfg(1))
               call formatter%open("pchem_internal_restart_in",pFIO_READ,rc=status)
               cfg(1) = formatter%read(rc=status)
               call MAPL_IOCountNonDimVars(cfg(1),nvars,rc=status)
               if (nVars /= iq_pchem1-iq_pchem0+1) call mpp_error(FATAL,'Wrong number of variables in pchem file')

               allocate(vnames(nVars))
               vars => cfg(1)%get_variables()
               iter = vars%begin()

               do while(iter /= vars%end())
                  var_name => iter%key()
                  if (.not.cfg(1)%is_coordinate_variable(var_name)) then
                     lvar_cnt=lvar_cnt+1
                     vnames(lvar_cnt)=var_name
                  end if
               enddo

            end if

            do ivar=iq_pchem0,iq_pchem1
               if (filetype ==0) lvar_cnt=lvar_cnt+1
               do k=1,km
                  if (filetype /= 0) then
                     read (IUNIT, IOSTAT=status) r4latlon
                  else
                     call MAPL_VarRead(formatter,vnames(lvar_cnt),r4latlon,lev=k)
                  end if
                  q0(1       :im/2,:,k) = r4latlon(im/2 + 1 :im  , :) ! Regrid from -180:180 to 0:360
                  q0(im/2 + 1:im  ,:,k) = r4latlon(1        :im/2, :) ! Regrid from -180:180 to 0:360
                  do j=js,je
                     do i=is,ie
                        i1 = id1(i,j)
                        i2 = id2(i,j)
                        j1 = jdc(i,j)
                        qp(i,j,k,ivar) = s2c(i,j,1)*q0(i1,j1  ,k) + s2c(i,j,2)*q0(i2,j1  ,k) +  &
                              s2c(i,j,3)*q0(i2,j1+1,k) + s2c(i,j,4)*q0(i1,j1+1,k)
                     enddo
                  enddo

               enddo
               if (is_master()) call pmaxmin( 'PCHEM_Q_',  q0(:,:,:), im*jm, km, 1.0_FVPRC)
            enddo

            if (filetype == 0) then
               call formatter%close()
               deallocate(cfg)
               deallocate(vnames)
            else
               close(IUNIT)
            end if
            deallocate(r4latlon)
         end if

         call print_memuse_stats('get_geos_latlon_ic: remap_tracers')
         deallocate ( q0 )

! Horiz Interp for T
         if(is_master()) call pmaxmin( 'T_geos',   t0, im*jm, km, 1.0_FVPRC) 
         allocate (  tp(is:ie,js:je,km) )
         do k=1,km
            do j=js,je
               do i=is,ie
                  i1 = id1(i,j)
                  i2 = id2(i,j)
                  j1 = jdc(i,j)
                  tp(i,j,k) = s2c(i,j,1)*t0(i1,j1  ,k) + s2c(i,j,2)*t0(i2,j1  ,k) +  &
                        s2c(i,j,3)*t0(i2,j1+1,k) + s2c(i,j,4)*t0(i1,j1+1,k)
               enddo
            enddo
         enddo
         deallocate ( t0 )
         call print_memuse_stats('get_geos_latlon_ic: remap_t')

! Horz/Vert remap for MOIST, GOCART, and PCHEM scalars (Assuming Total Number is divisible by KM)
! -----------------------------------------------------------------------------------------------
         nqmap = nmoist + ngocart + npchem 

         call remap_scalar(im, jm, km, npz, nqmap, nqmap, ak0, bk0, psc, gzc, tp, qp, Atm(1))

         deallocate ( tp )
         deallocate ( qp ) 
         call print_memuse_stats('get_geos_latlon_ic: remap_scalar')

! Horz/Vert remap for U/V
         call remap_winds(im, jm, km, npz, ak0, bk0, psc, ua, va, Atm)
         deallocate ( ua )
         deallocate ( va )
         call print_memuse_stats('get_geos_latlon_ic: remap_winds')

      else
         call mpp_error(FATAL,'==> Error from get_external_ic:        &
               & Expected file '//trim(fname)//' does not exist')
      endif

      if (allocated(bk0)) deallocate ( bk0 )
      if (allocated(ak0)) deallocate ( ak0 )
      if (allocated(lat)) deallocate ( lat )
      if (allocated(lon)) deallocate ( lon )

! Finished, let's check the results !

      call prt_maxmin('GZ_model', Atm(1)%phis, is, ie, js, je, ng,   1, 1.0/grav)
      call prt_maxmin('PS_model', Atm(1)%ps  , is, ie, js, je, ng,   1, 0.01_FVPRC)
      call prt_maxmin('DP_model', Atm(1)%delp, is, ie, js, je, ng, npz, 1.0_FVPRC)
      call prt_maxmin(' U_model', Atm(1)%u   , is, ie, js, je, ng, npz, 1.0_FVPRC)
      call prt_maxmin(' V_model', Atm(1)%v   , is, ie, js, je, ng, npz, 1.0_FVPRC)
      call prt_maxmin('PT_model', Atm(1)%pt  , is, ie, js, je, ng, npz, 1.0_FVPRC)

! Range check the MOIST tracers
      do iq=iq_moist0,iq_moist1
         do k=1,npz
            do j=js,je
               do i=is,ie
                  Atm(1)%q(i,j,k,iq) = MIN(Atm(1)%q(i,j,k,iq),1.d0)
                  Atm(1)%q(i,j,k,iq) = MAX(Atm(1)%q(i,j,k,iq),0.d0)
               enddo
            enddo
         enddo
      enddo

      if (is_master()) print*
      do iq=iq_moist0,iq_moist1
         call prt_maxmin('QP_MOIST_Q', Atm(1)%q(is:ie,js:je,1:npz,iq), is, ie, js, je, 0, npz, 1._FVPRC)
      enddo
      if (is_master()) print*
      do iq=iq_gocart0,iq_gocart1
         call prt_maxmin('QP_GOCART_Q', Atm(1)%q(is:ie,js:je,1:npz,iq), is, ie, js, je, 0, npz, 1._FVPRC)
      enddo
      if (is_master()) print*
      do iq=iq_pchem0,iq_pchem1
         call prt_maxmin('QP_PCHEM_Q', Atm(1)%q(is:ie,js:je,1:npz,iq), is, ie, js, je, 0, npz, 1._FVPRC)
      enddo

   end subroutine get_geos_latlon_ic

   subroutine get_ncep_ic( Atm, fv_domain, nq )
      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      integer, intent(in):: nq
! local:
      character(len=128) :: fname
      real(kind=4), allocatable:: wk1(:), wk2(:,:), wk3(:,:,:)
      real(FVPRC), allocatable:: tp(:,:,:), qp(:,:,:)
      real(FVPRC), allocatable:: ua(:,:,:), va(:,:,:)
      real(FVPRC), allocatable:: lat(:), lon(:), ak0(:), bk0(:)
      real(FVPRC):: s2c(is:ie,js:je,4)
      integer, dimension(is:ie,js:je):: id1, id2, jdc
      real(FVPRC) psc(is:ie,js:je)
      real(FVPRC) gzc(is:ie,js:je)
      real(FVPRC) tmean
      integer:: i, j, k, im, jm, km, npz, npt
      integer:: i1, i2, j1, ncid
      integer:: jbeg, jend
      integer tsize(4)
      logical:: read_ts = .true.
      logical:: land_ts = .false.
      logical:: found
#ifndef MAPL_MODE

      deg2rad = pi/180.

      npz = Atm(1)%npz

! Zero out all initial tracer fields:
      Atm(1)%q = 0.

      fname = Atm(1)%res_latlon_dynamics

      if( file_exist(fname) ) then
         call open_ncfile( fname, ncid )        ! open the file
         call get_ncdim1( ncid, 'lon', tsize(1) )
         call get_ncdim1( ncid, 'lat', tsize(2) )
         call get_ncdim1( ncid, 'lev', tsize(3) )

         im = tsize(1); jm = tsize(2); km = tsize(3)

         if(is_master())  write(*,*) fname, ' NCEP IC dimensions:', tsize

         allocate (  lon(im) )
         allocate (  lat(jm) )

         call get_var1_double (ncid, 'lon', im, lon )
         call get_var1_double (ncid, 'lat', jm, lat )

! Convert to radian
         do i=1,im
            lon(i) = lon(i) * deg2rad  ! lon(1) = 0.
         enddo
         do j=1,jm
            lat(j) = lat(j) * deg2rad
         enddo

         allocate ( ak0(km+1) )
         allocate ( bk0(km+1) )
         call get_var1_double (ncid, 'hyai', km+1, ak0, found )
         if ( .not. found )  ak0(:) = 0.

         call get_var1_double (ncid, 'hybi', km+1, bk0 )

! Note: definition of NCEP hybrid is p(k) = a(k)*1.E5 + b(k)*ps
         ak0(:) = ak0(:) * 1.E5

! Limiter to prevent NAN at top during remapping
         if ( bk0(1) < 1.E-9 ) ak0(1) = max(1.e-9, ak0(1))
      else
         call mpp_error(FATAL,'==> Error in get_external_ic: Expected file '//trim(fname)//' does not exist')
      endif

! Initialize lat-lon to Cubed bi-linear interpolation coeff:
      call remap_coef( im, jm, lon, lat, id1, id2, jdc, s2c, Atm(1)%gridstruct%agrid, Atm(1)%bd )

! Find bounding latitudes:
      jbeg = jm-1;         jend = 2
      do j=js,je
         do i=is,ie
            j1 = jdc(i,j)
            jbeg = min(jbeg, j1) 
            jend = max(jend, j1+1)
         enddo
      enddo

! remap surface pressure and height:

      allocate ( wk2(im,jbeg:jend) )
      call get_var3_r4( ncid, 'PS', 1,im, jbeg,jend, 1,1, wk2 )

      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            psc(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                  s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
         enddo
      enddo

      call get_var3_r4( ncid, 'PHIS', 1,im, jbeg,jend, 1,1, wk2 )
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            gzc(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                  s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
         enddo
      enddo

      deallocate ( wk2 )
      allocate ( wk2(im,jm) )

      if ( read_ts ) then       ! read skin temperature; could be used for SST

         call get_var2_real( ncid, 'TS', im, jm, wk2 )

         if ( .not. land_ts ) then
            allocate ( wk1(im) )

            do j=1,jm
! Read NCEP ORO (1; land; 0: ocean; 2: sea_ice)
               call get_var3_r4( ncid, 'ORO', 1,im, j,j, 1,1, wk1 )
               tmean = 0.
               npt = 0
               do i=1,im
                  if( abs(wk1(i)-1.) > 0.99 ) then   ! ocean or sea ice
                     tmean = tmean + wk2(i,j)
                     npt = npt + 1
                  endif
               enddo
!------------------------------------------------------
! Replace TS over interior land with zonal mean SST/Ice
!------------------------------------------------------
               if ( npt /= 0 ) then
                  tmean= tmean / real(npt)
                  do i=1,im
                     if( abs(wk1(i)-1.) <= 0.99 ) then  ! Land points
                        if ( i==1 ) then
                           i1 = im;     i2 = 2
                        elseif ( i==im ) then
                           i1 = im-1;   i2 = 1
                        else
                           i1 = i-1;    i2 = i+1
                        endif
                        if ( abs(wk1(i2)-1.)>0.99 ) then     ! east side has priority
                           wk2(i,j) = wk2(i2,j)
                        elseif ( abs(wk1(i1)-1.)>0.99 ) then ! west side
                           wk2(i,j) = wk2(i1,j)
                        else
                           wk2(i,j) = tmean
                        endif
                     endif
                  enddo
               endif
            enddo   ! j-loop
            deallocate ( wk1 )
         endif   !(.not.land_ts)

         do j=js,je
            do i=is,ie
               i1 = id1(i,j)
               i2 = id2(i,j)
               j1 = jdc(i,j)
               Atm(1)%ts(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                     s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
            enddo
         enddo
         call prt_maxmin('SST_model', Atm(1)%ts, is, ie, js, je, 0, 1, 1._FVPRC)

! Perform interp to FMS SST format/grid
#ifndef DYCORE_SOLO
         call ncep2fms(im, jm, lon, lat, wk2)
         if( is_master() ) then
            write(*,*) 'External_ic_mod: i_sst=', i_sst, ' j_sst=', j_sst
            call pmaxmin( 'SST_ncep_fms',  sst_ncep, i_sst, j_sst, 1._FVPRC)
         endif
#endif
      endif  !(read_ts)

      deallocate ( wk2 )

! Read in temperature:
      allocate ( wk3(1:im,jbeg:jend, 1:km) )
      call get_var3_r4( ncid, 'T', 1,im, jbeg,jend, 1,km, wk3 )

      allocate (  tp(is:ie,js:je,km) )
      do k=1,km
         do j=js,je
            do i=is,ie
               i1 = id1(i,j)
               i2 = id2(i,j)
               j1 = jdc(i,j)
               tp(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                     s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
            enddo
         enddo
      enddo

! Read in tracers: only sphum at this point
      call get_var3_r4( ncid, 'Q', 1,im, jbeg,jend, 1,km, wk3 )

      allocate ( qp(is:ie,js:je,km) )
      do k=1,km
         do j=js,je
            do i=is,ie
               i1 = id1(i,j)
               i2 = id2(i,j)
               j1 = jdc(i,j)
               qp(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                     s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
            enddo
         enddo
      enddo

      call remap_scalar(im, jm, km, npz, nq, nq, ak0, bk0, psc, gzc, tp, qp, Atm)
      deallocate ( tp )
      deallocate ( qp )

! Winds:
      call get_var3_r4( ncid, 'U', 1,im, jbeg,jend, 1,km, wk3 )

      allocate ( ua(is:ie,js:je,km) )
      do k=1,km
         do j=js,je
            do i=is,ie
               i1 = id1(i,j)
               i2 = id2(i,j)
               j1 = jdc(i,j)
               ua(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                     s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
            enddo
         enddo
      enddo

      call get_var3_r4( ncid, 'V', 1,im, jbeg,jend, 1,km, wk3 )
      call close_ncfile ( ncid )

      allocate ( va(is:ie,js:je,km) )
      do k=1,km
         do j=js,je
            do i=is,ie
               i1 = id1(i,j)
               i2 = id2(i,j)
               j1 = jdc(i,j)
               va(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                     s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
            enddo
         enddo
      enddo
      deallocate ( wk3 )

      call remap_winds(im, jm, km, npz, ak0, bk0, psc, ua, va, Atm)

      deallocate ( ua )
      deallocate ( va )

      deallocate ( ak0 )
      deallocate ( bk0 )
      deallocate ( lat )
      deallocate ( lon )
#endif

   end subroutine get_ncep_ic



   subroutine get_fv_ic( Atm, fv_domain, nq )
      type(fv_atmos_type), intent(inout) :: Atm(:)
      type(domain2d),      intent(inout) :: fv_domain
      integer, intent(in):: nq

      character(len=128) :: fname, tracer_name
      real(FVPRC), allocatable:: ps0(:,:), gz0(:,:), u0(:,:,:), v0(:,:,:), t0(:,:,:), dp0(:,:,:), q0(:,:,:,:)
      real(FVPRC), allocatable:: ua(:,:,:), va(:,:,:)
      real(FVPRC), allocatable:: lat(:), lon(:), ak0(:), bk0(:)
      integer :: i, j, k, im, jm, km, npz, tr_ind
      integer tsize(4)
!     integer sphum, liq_wat, ice_wat, cld_amt       ! GFDL AM2 physics
      logical found

#ifndef MAPL_MODE
      npz = Atm(1)%npz

! Zero out all initial tracer fields:
      Atm(1)%q = 0.

! Read in lat-lon FV core restart file
      fname = Atm(1)%res_latlon_dynamics

      if( file_exist(fname) ) then
         call field_size(fname, 'T', tsize, field_found=found)
         if(is_master()) write(*,*) 'Using lat-lon FV restart:', fname 

         if ( found ) then
            im = tsize(1); jm = tsize(2); km = tsize(3)
            if(is_master())  write(*,*) 'External IC dimensions:', tsize
         else
            call mpp_error(FATAL,'==> Error in get_external_ic: field not found')
         endif

! Define the lat-lon coordinate:
         allocate (  lon(im) )
         allocate (  lat(jm) )

         do i=1,im
            lon(i) = (0.5 + real(i-1)) * 2.*pi/real(im)
         enddo

         do j=1,jm
            lat(j) = -0.5*pi + real(j-1)*pi/real(jm-1)   ! SP to NP 
         enddo

         allocate ( ak0(1:km+1) )
         allocate ( bk0(1:km+1) )
         allocate ( ps0(1:im,1:jm) )
         allocate ( gz0(1:im,1:jm) )
         allocate (  u0(1:im,1:jm,1:km) )
         allocate (  v0(1:im,1:jm,1:km) )
         allocate (  t0(1:im,1:jm,1:km) )
         allocate ( dp0(1:im,1:jm,1:km) )

         call read_data (fname, 'ak', ak0)
         call read_data (fname, 'bk', bk0)
         call read_data (fname, 'Surface_geopotential', gz0)
         call read_data (fname, 'U',     u0)
         call read_data (fname, 'V',     v0)
         call read_data (fname, 'T',     t0)
         call read_data (fname, 'DELP', dp0)

! Share the load
         if(is_master()) call pmaxmin( 'ZS_data', gz0, im,    jm, 1./grav)
         if(is_master()) call pmaxmin( 'U_data',   u0, im*jm, km, 1._FVPRC)
         if(is_master()) call pmaxmin( 'V_data',   v0, im*jm, km, 1._FVPRC)
         if(is_master()) call pmaxmin( 'T_data',   t0, im*jm, km, 1._FVPRC)
         if(is_master()) call pmaxmin( 'DEL-P',   dp0, im*jm, km, 0.01_FVPRC)


      else
         call mpp_error(FATAL,'==> Error in get_external_ic: Expected file '//trim(fname)//' does not exist')
      endif

! Read in tracers: only AM2 "physics tracers" at this point
      fname = Atm(1)%res_latlon_tracers

      if( file_exist(fname) ) then
         if(is_master()) write(*,*) 'Using lat-lon tracer restart:', fname 

         allocate ( q0(im,jm,km,Atm(1)%ncnst) )
         q0 = 0.

         do tr_ind = 1, nq
            call get_tracer_names(MODEL_ATMOS, tr_ind, tracer_name)
            if (field_exist(fname,tracer_name)) then
               call read_data(fname, tracer_name, q0(1:im,1:jm,1:km,tr_ind))
               call mpp_error(NOTE,'==>  Have read tracer '//trim(tracer_name)//' from '//trim(fname))
               cycle
            endif
         enddo
      else
         call mpp_error(FATAL,'==> Error in get_external_ic: Expected file '//trim(fname)//' does not exist')
      endif

! D to A transform on lat-lon grid:
      allocate (  ua(im,jm,km) )
      allocate (  va(im,jm,km) )

      call d2a3d(u0, v0,  ua,  va, im, jm, km, lon)

      deallocate ( u0 ) 
      deallocate ( v0 ) 

      if(is_master()) call pmaxmin( 'UA', ua, im*jm, km, 1._FVPRC)
      if(is_master()) call pmaxmin( 'VA', va, im*jm, km, 1._FVPRC)

      do j=1,jm
         do i=1,im
            ps0(i,j) = ak0(1)
         enddo
      enddo

      do k=1,km
         do j=1,jm
            do i=1,im
               ps0(i,j) = ps0(i,j) + dp0(i,j,k)
            enddo
         enddo
      enddo

      if (is_master()) call pmaxmin( 'PS_data (mb)', ps0, im, jm, 0.01_FVPRC)

! Horizontal interpolation to the cubed sphere grid center
! remap vertically with terrain adjustment

      call remap_xyz( im, 1, jm, jm, km, npz, nq, Atm(1)%ncnst, lon, lat, ak0, bk0,   &
            ps0,  gz0, ua, va, t0, q0, Atm )

      deallocate ( ak0 ) 
      deallocate ( bk0 ) 
      deallocate ( ps0 ) 
      deallocate ( gz0 ) 
      deallocate ( t0 ) 
      deallocate ( q0 ) 
      deallocate ( dp0 ) 
      deallocate ( ua ) 
      deallocate ( va ) 
      deallocate ( lat ) 
      deallocate ( lon ) 
#endif

   end subroutine get_fv_ic


#ifndef DYCORE_SOLO
   subroutine ncep2fms(im, jm, lon, lat, wk)

      integer, intent(in):: im, jm
      real(FVPRC),    intent(in):: lon(im), lat(jm)
      real(kind=4),    intent(in):: wk(im,jm)
! local:
      real(FVPRC) :: rdlon(im)
      real(FVPRC) :: rdlat(jm)
      real(FVPRC):: a1, b1
      real(FVPRC):: delx, dely
      real(FVPRC):: xc, yc    ! "data" location
      real(FVPRC):: c1, c2, c3, c4
      integer i,j, i1, i2, jc, i0, j0, it, jt

      do i=1,im-1
         rdlon(i) = 1. / (lon(i+1) - lon(i))
      enddo
      rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

      do j=1,jm-1
         rdlat(j) = 1. / (lat(j+1) - lat(j))
      enddo

! * Interpolate to "FMS" 1x1 SST data grid
! lon: 0.5, 1.5, ..., 359.5
! lat: -89.5, -88.5, ... , 88.5, 89.5

      delx = 360./real(i_sst)
      dely = 180./real(j_sst)

      jt = 1
      do 5000 j=1,j_sst

         yc = (-90. + dely * (0.5+real(j-1)))  * deg2rad
         if ( yc<lat(1) ) then
            jc = 1
            b1 = 0.
         elseif ( yc>lat(jm) ) then
            jc = jm-1
            b1 = 1.
         else
            do j0=jt,jm-1
               if ( yc>=lat(j0) .and. yc<=lat(j0+1) ) then
                  jc = j0
                  jt = j0
                  b1 = (yc-lat(jc)) * rdlat(jc)
                  go to 222
               endif
            enddo
         endif
222      continue
         it = 1

         do i=1,i_sst
            xc = delx * (0.5+real(i-1)) * deg2rad
            if ( xc>lon(im) ) then
               i1 = im;     i2 = 1
               a1 = (xc-lon(im)) * rdlon(im)
            elseif ( xc<lon(1) ) then
               i1 = im;     i2 = 1
               a1 = (xc+2.*pi-lon(im)) * rdlon(im)
            else
               do i0=it,im-1
                  if ( xc>=lon(i0) .and. xc<=lon(i0+1) ) then
                     i1 = i0;  i2 = i0+1
                     it = i0
                     a1 = (xc-lon(i1)) * rdlon(i0)
                     go to 111
                  endif
               enddo
            endif
111         continue

            if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
               write(*,*) 'gid=', mpp_pe(), i,j,a1, b1
            endif

            c1 = (1.-a1) * (1.-b1)
            c2 =     a1  * (1.-b1)
            c3 =     a1  *     b1
            c4 = (1.-a1) *     b1
! Interpolated surface pressure
            sst_ncep(i,j) = c1*wk(i1,jc  ) + c2*wk(i2,jc  ) +    &
                  c3*wk(i2,jc+1) + c4*wk(i1,jc+1)
         enddo   !i-loop
5000     continue   ! j-loop

      end subroutine ncep2fms
#endif


 subroutine remap_coef( im, jm, lon, lat, id1, id2, jdc, s2c, agrid, bd )

   type(fv_grid_bounds_type), intent(IN) :: bd
  integer, intent(in):: im, jm
  real(FVPRC),    intent(in):: lon(im), lat(jm)
  real(FVPRC),    intent(out):: s2c(bd%is:bd%ie,bd%js:bd%je,4)
  integer, intent(out), dimension(bd%is:bd%ie,bd%js:bd%je):: id1, id2, jdc
  real(FVPRC),    intent(in):: agrid(bd%isd:bd%ied,bd%jsd:bd%jed,2)
! local:
  real(FVPRC) :: rdlon(im)
  real(FVPRC) :: rdlat(jm)
  real(FVPRC):: a1, b1
  integer i,j, i1, i2, jc, i0, j0

  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = bd%is
  ie  = bd%ie
  js  = bd%js
  je  = bd%je
  isd = bd%isd
  ied = bd%ied
  jsd = bd%jsd
  jed = bd%jed
  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

! * Interpolate to cubed sphere cell center
  do 5000 j=js,je

     do i=is,ie

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

       if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
            write(*,*) 'gid=', mpp_pe(), i,j,a1, b1
       endif

       s2c(i,j,1) = (1.-a1) * (1.-b1)
       s2c(i,j,2) =     a1  * (1.-b1)
       s2c(i,j,3) =     a1  *     b1
       s2c(i,j,4) = (1.-a1) *     b1
       id1(i,j) = i1
       id2(i,j) = i2
       jdc(i,j) = jc
     enddo   !i-loop
5000 continue   ! j-loop

 end subroutine remap_coef

 subroutine remap_scalar(im, jm, km, npz, nq, ncnst, ak0, bk0, psc, gzc, ta, qa, Atm)
  type(fv_atmos_type), intent(inout) :: Atm
  integer, intent(in):: im, jm, km, npz, nq, ncnst
  real(FVPRC),    intent(in):: ak0(km+1), bk0(km+1)
  real(FVPRC),    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je):: psc, gzc
  real(FVPRC),    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km):: ta
  real(FVPRC),    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km,ncnst):: qa
! local:
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,km):: tp
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,km+1):: pe0, pn0
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,npz):: qn1
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pe1, pn1
  real(FVPRC) pt0(km), gz(km+1), pk0(km+1)
  real(FVPRC) qp(Atm%bd%is:Atm%bd%ie,km,ncnst)
  real(FVPRC) pst, p1, p2, alpha, rdg
  integer i,j,k, iq
  integer  sphum
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = Atm%bd%is
  ie  = Atm%bd%ie
  js  = Atm%bd%js
  je  = Atm%bd%je
  isd = Atm%bd%isd
  ied = Atm%bd%ied
  jsd = Atm%bd%jsd
  jed = Atm%bd%jed


#ifdef MAPL_MODE
            sphum   = 1
#else
            sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
#endif
  if ( sphum/=1 ) then
       call mpp_error(FATAL,'SPHUM must be 1st tracer')
  endif

  do 5000 j=js,je
     do i=is,ie

       do iq=1,ncnst
          do k=1,km
             qp(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo

    if ( T_is_Tv ) then
! The "T" field in NCEP analysis is actually virtual temperature (Larry H. post processing)
! BEFORE 20051201
       do k=1,km
          tp(i,k) = ta(i,j,k)
       enddo
    else
       do k=1,km
          tp(i,k) = ta(i,j,k)*(1.+zvir*qp(i,k,sphum))
       enddo
    endif
! Tracers:

       do k=1,km+1
          pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
          pn0(i,k) = log(pe0(i,k))
            pk0(k) = pe0(i,k)**kappa
       enddo

#ifdef USE_DATA_ZS
       Atm%  ps(i,j) = psc(i,j)
       Atm%phis(i,j) = gzc(i,j)
#else

! * Adjust interpolated ps to model terrain
       gz(km+1) = gzc(i,j)
       do k=km,1,-1
           gz(k) = gz(k+1) + rdgas*tp(i,k)*(pn0(i,k+1)-pn0(i,k))
       enddo
! Only lowest layer potential temp is needed
          pt0(km) = tp(i,km)/(pk0(km+1)-pk0(km))*(kappa*(pn0(i,km+1)-pn0(i,km)))
       if( Atm%phis(i,j)>gzc(i,j) ) then
           do k=km,1,-1
              if( Atm%phis(i,j) <  gz(k)  .and.    &
                  Atm%phis(i,j) >= gz(k+1) ) then
                  pst = pk0(k) + (pk0(k+1)-pk0(k))*(gz(k)-Atm%phis(i,j))/(gz(k)-gz(k+1))
                  go to 123
              endif
           enddo
       else
! Extrapolation into the ground
           pst = pk0(km+1) + (gzc(i,j)-Atm%phis(i,j))/(cp_air*pt0(km))
       endif

123    Atm%ps(i,j) = pst**(1./kappa)
#endif
     enddo   !i-loop


     do i=is,ie
        pe1(i,1) = Atm%ak(1)
        pn1(i,1) = log(pe1(i,1))
     enddo
     do k=2,npz+1
       do i=is,ie
          pe1(i,k) = Atm%ak(k) + Atm%bk(k)*Atm%ps(i,j)
          pn1(i,k) = log(pe1(i,k))
       enddo
     enddo

! * Compute delp
     do k=1,npz
        do i=is,ie
           Atm%delp(i,j,k) = pe1(i,k+1) - pe1(i,k)
        enddo
     enddo

!---------------
! map tracers
!----------------
      do iq=1,ncnst
         call mappm(km, pe0, qp(is,1,iq), npz, pe1,  qn1, is,ie, 0, 11, Atm%ptop)
         if ( iq==sphum .and. Atm%flagstruct%ncep_ic ) then
             p1 = 200.E2
             p2 =  75.E2
! Blend model sphum with NCEP data
         do k=1,npz
            do i=is,ie
               pst = 0.5*(pe1(i,k)+pe1(i,k+1))
               if ( pst > p1 ) then
                    Atm%q(i,j,k,iq) = qn1(i,k)
               elseif( pst > p2 ) then            ! p2 < pst < p1
                    alpha = (pst-p2)/(p1-p2)
                    Atm%q(i,j,k,1) = qn1(i,k)*alpha + Atm%q(i,j,k,1)*(1.-alpha)
               endif
            enddo
         enddo
         else
         do k=1,npz
            do i=is,ie
               Atm%q(i,j,k,iq) = qn1(i,k)
            enddo
         enddo
         endif
      enddo

!-------------------------------------------------------------
! map virtual temperature using geopotential conserving scheme.
!-------------------------------------------------------------
      call mappm(km, pn0, tp, npz, pn1, qn1, is,ie, 1, 9, Atm%ptop)
      do k=1,npz
         do i=is,ie
            Atm%pt(i,j,k) = qn1(i,k)/(1.+zvir*Atm%q(i,j,k,sphum))
         enddo
      enddo

      if ( .not. Atm%flagstruct%hydrostatic .and. Atm%flagstruct%ncep_ic ) then
! Replace delz with NCEP hydrostatic state
         rdg = -rdgas / grav
         do k=1,npz
            do i=is,ie
               atm%delz(i,j,k) = rdg*qn1(i,k)*(pn1(i,k+1)-pn1(i,k))
            enddo
         enddo
      endif

5000 continue

  call prt_maxmin('PS_model', Atm%ps, is, ie, js, je, ng, 1, 0.01_FVPRC)

  if (is_master()) write(*,*) 'done remap_scalar'

 end subroutine remap_scalar

            subroutine remap_winds(im, jm, km, npz, ak0, bk0, psc, ua, va, Atm)
               type(fv_atmos_type), intent(inout) :: Atm(:)
               integer, intent(in):: im, jm, km, npz
               real(FVPRC),    intent(in):: ak0(km+1), bk0(km+1)
               real(FVPRC),    intent(in):: psc(is:ie,js:je)
               real(FVPRC),    intent(in), dimension(is:ie,js:je,km):: ua, va
! local:
               real(FVPRC), dimension(isd:ied,jsd:jed,npz):: ut, vt   ! winds
               real(FVPRC), dimension(is:ie, km+1):: pe0
               real(FVPRC), dimension(is:ie,npz+1):: pe1
               real(FVPRC), dimension(is:ie,npz):: qn1
               integer i,j,k

               ut = 0.0
               vt = 0.0
               do 5000 j=js,je

                  do k=1,km+1
                     do i=is,ie
                        pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
                     enddo
                  enddo

                  do k=1,npz+1
                     do i=is,ie
                        pe1(i,k) = Atm(1)%ak(k) + Atm(1)%bk(k)*Atm(1)%ps(i,j)
                     enddo
                  enddo

!------
! map u
!------
                  call mappm(km, pe0, ua(is:ie,j,1:km), npz, pe1, qn1, is,ie, -1, 4, Atm(1)%ptop)
                  do k=1,npz
                     do i=is,ie
                        ut(i,j,k) = qn1(i,k)
                     enddo
                  enddo
!------
! map v
!------
                  call mappm(km, pe0, va(is:ie,j,1:km), npz, pe1, qn1, is,ie, -1, 4, Atm(1)%ptop)
                  do k=1,npz
                     do i=is,ie
                        vt(i,j,k) = qn1(i,k)
                     enddo
                  enddo

5000              continue

                  call prt_maxmin('UT', ut, is, ie, js, je, ng, npz, 1.0_FVPRC)
                  call prt_maxmin('VT', vt, is, ie, js, je, ng, npz, 1.0_FVPRC)

!----------------------------------------------
! winds: lat-lon ON A to Cubed-D transformation:
!----------------------------------------------
                  call cubed_a2d(Atm(1)%npx, Atm(1)%npy, npz, ut, vt, Atm(1)%u, Atm(1)%v, Atm(1)%gridstruct, &
                  Atm(1)%domain, Atm(1)%bd )

                  if (is_master()) write(*,*) 'done remap_winds'

               end subroutine remap_winds


               subroutine remap_wz(im, jm, km, npz, mg, ak0, bk0, psc, wa, wz, Atm)
                  type(fv_atmos_type), intent(inout) :: Atm(:)
                  integer, intent(in):: im, jm, km, npz
                  integer, intent(in):: mg     ! mg = 0 for delz; mg=3 for w
                  real(FVPRC),    intent(in):: ak0(km+1), bk0(km+1)
                  real(FVPRC),    intent(in):: psc(is:ie,js:je)
                  real(FVPRC),    intent(in), dimension(is:ie,js:je,km):: wa
                  real(FVPRC),   intent(out):: wz(is-mg:ie+mg,js-mg:je+mg,npz)
! local:
                  real(FVPRC), dimension(is:ie, km+1):: pe0
                  real(FVPRC), dimension(is:ie,npz+1):: pe1
                  real(FVPRC), dimension(is:ie,npz):: qn1
                  integer i,j,k

                  do 5000 j=js,je

                     do k=1,km+1
                        do i=is,ie
                           pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
                        enddo
                     enddo

                     do k=1,npz+1
                        do i=is,ie
                           pe1(i,k) = Atm(1)%ak(k) + Atm(1)%bk(k)*Atm(1)%ps(i,j)
                        enddo
                     enddo

!------
! map w
!------
                     call mappm(km, pe0, wa(is:ie,j,1:km), npz, pe1, qn1, is,ie, -1, 4, Atm(1)%ptop)
                     do k=1,npz
                        do i=is,ie
                           wz(i,j,k) = qn1(i,k)
                        enddo
                     enddo

5000                 continue

! call prt_maxmin('WZ', wz, is, ie, js, je, mg, npz, 1._FVPRC, is_master())
! if (is_master()) write(*,*) 'done remap_wz'

                  end subroutine remap_wz

  subroutine remap_xyz( im, jbeg, jend, jm, km, npz, nq, ncnst, lon, lat, ak0, bk0, ps0, gz0,   &
                        ua, va, ta, qa, Atm )

  type(fv_atmos_type), intent(inout), target :: Atm
  integer, intent(in):: im, jm, km, npz, nq, ncnst
  integer, intent(in):: jbeg, jend
  real(FVPRC),    intent(in):: lon(im), lat(jm), ak0(km+1), bk0(km+1)
  real(FVPRC),    intent(in):: gz0(im,jbeg:jend), ps0(im,jbeg:jend)
  real(FVPRC),    intent(in), dimension(im,jbeg:jend,km):: ua, va, ta
  real(FVPRC),    intent(in), dimension(im,jbeg:jend,km,ncnst):: qa

  real(FVPRC), pointer, dimension(:,:,:) :: agrid

! local:
  real(FVPRC), dimension(Atm%bd%isd:Atm%bd%ied,Atm%bd%jsd:Atm%bd%jed,npz):: ut, vt   ! winds 
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,km):: up, vp, tp
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,km+1):: pe0, pn0
  real(FVPRC) pt0(km), gz(km+1), pk0(km+1)
  real(FVPRC) qp(Atm%bd%is:Atm%bd%ie,km,ncnst)
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,npz):: qn1
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pe1, pn1
  real(FVPRC) :: rdlon(im)
  real(FVPRC) :: rdlat(jm)
  real(FVPRC):: a1, b1, c1, c2, c3, c4
  real(FVPRC):: gzc, psc, pst
  integer i,j,k, i1, i2, jc, i0, j0, iq
! integer  sphum, liq_wat, ice_wat, cld_amt
  integer  sphum
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = Atm%bd%is
  ie  = Atm%bd%ie
  js  = Atm%bd%js
  je  = Atm%bd%je
  isd = Atm%bd%isd
  ied = Atm%bd%ied
  jsd = Atm%bd%jsd
  jed = Atm%bd%jed

  !!NOTE: Only Atm is used in this routine.
  agrid => Atm%gridstruct%agrid

  sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
! liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
! ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
! cld_amt = get_tracer_index(MODEL_ATMOS, 'cld_amt')

   if ( sphum/=1 ) then
        call mpp_error(FATAL,'SPHUM must be 1st tracer')
   endif
  pk0(1) = ak0(1)**kappa

  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

! * Interpolate to cubed sphere cell center
  do 5000 j=js,je

     do i=is,ie
        pe0(i,1) = ak0(1)
        pn0(i,1) = log(ak0(1))
     enddo

     do i=is,ie

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

#ifndef DEBUG_REMAP
       if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
            write(*,*) i,j,a1, b1
       endif
#endif
       c1 = (1.-a1) * (1.-b1)
       c2 =     a1  * (1.-b1)
       c3 =     a1  *     b1
       c4 = (1.-a1) *     b1

! Interpolated surface pressure
       psc = c1*ps0(i1,jc  ) + c2*ps0(i2,jc  ) +    &
             c3*ps0(i2,jc+1) + c4*ps0(i1,jc+1)

! Interpolated surface geopotential
       gzc = c1*gz0(i1,jc  ) + c2*gz0(i2,jc  ) +    &
             c3*gz0(i2,jc+1) + c4*gz0(i1,jc+1)

! 3D fields:
       do iq=1,ncnst
!          if ( iq==sphum .or. iq==liq_wat .or. iq==ice_wat .or. iq==cld_amt ) then
          do k=1,km
             qp(i,k,iq) = c1*qa(i1,jc,  k,iq) + c2*qa(i2,jc,  k,iq) +  &
                          c3*qa(i2,jc+1,k,iq) + c4*qa(i1,jc+1,k,iq)
          enddo
!          endif
       enddo
       do k=1,km
          up(i,k) = c1*ua(i1,jc,  k) + c2*ua(i2,jc,  k) +  &
                    c3*ua(i2,jc+1,k) + c4*ua(i1,jc+1,k)
          vp(i,k) = c1*va(i1,jc,  k) + c2*va(i2,jc,  k) +  &
                    c3*va(i2,jc+1,k) + c4*va(i1,jc+1,k)
          tp(i,k) = c1*ta(i1,jc,  k) + c2*ta(i2,jc,  k) +  &
                    c3*ta(i2,jc+1,k) + c4*ta(i1,jc+1,k)
! Virtual effect:
          tp(i,k) = tp(i,k)*(1.+zvir*qp(i,k,sphum))
       enddo
! Tracers:

       do k=2,km+1
          pe0(i,k) = ak0(k) + bk0(k)*psc
          pn0(i,k) = log(pe0(i,k))
          pk0(k) = pe0(i,k)**kappa
       enddo

#ifdef USE_DATA_ZS
       Atm%  ps(i,j) = psc
       Atm%phis(i,j) = gzc
#else

! * Adjust interpolated ps to model terrain
       gz(km+1) = gzc
       do k=km,1,-1
           gz(k) = gz(k+1) + rdgas*tp(i,k)*(pn0(i,k+1)-pn0(i,k))
       enddo
! Only lowest layer potential temp is needed
          pt0(km) = tp(i,km)/(pk0(km+1)-pk0(km))*(kappa*(pn0(i,km+1)-pn0(i,km)))
       if( Atm%phis(i,j)>gzc ) then
           do k=km,1,-1
              if( Atm%phis(i,j) <  gz(k)  .and.    &
                  Atm%phis(i,j) >= gz(k+1) ) then
                  pst = pk0(k) + (pk0(k+1)-pk0(k))*(gz(k)-Atm%phis(i,j))/(gz(k)-gz(k+1))
                  go to 123
              endif
           enddo
       else
! Extrapolation into the ground
           pst = pk0(km+1) + (gzc-Atm%phis(i,j))/(cp_air*pt0(km))
       endif

123    Atm%ps(i,j) = pst**(1./kappa)
#endif
     enddo   !i-loop


! * Compute delp from ps
     do i=is,ie
        pe1(i,1) = Atm%ak(1)
        pn1(i,1) = log(pe1(i,1))
     enddo
     do k=2,npz+1
       do i=is,ie
          pe1(i,k) = Atm%ak(k) + Atm%bk(k)*Atm%ps(i,j)
          pn1(i,k) = log(pe1(i,k))
       enddo
     enddo

     do k=1,npz
        do i=is,ie
           Atm%delp(i,j,k) = pe1(i,k+1) - pe1(i,k)
        enddo
     enddo

! Use kord=9 for winds; kord=11 for tracers
!------
! map u
!------
      call mappm(km, pe0, up, npz, pe1, qn1, is,ie, -1, 9, Atm%ptop)
      do k=1,npz
         do i=is,ie
            ut(i,j,k) = qn1(i,k)
         enddo
      enddo
!------
! map v
!------
      call mappm(km, pe0, vp, npz, pe1, qn1, is,ie, -1, 9, Atm%ptop)
      do k=1,npz
         do i=is,ie
            vt(i,j,k) = qn1(i,k)
         enddo
      enddo

!---------------
! map tracers
!----------------
      do iq=1,ncnst
! Note: AM2 physics tracers only
!         if ( iq==sphum .or. iq==liq_wat .or. iq==ice_wat .or. iq==cld_amt ) then
         call mappm(km, pe0, qp(is,1,iq), npz, pe1,  qn1, is,ie, 0, 11, Atm%ptop)
         do k=1,npz
            do i=is,ie
               Atm%q(i,j,k,iq) = qn1(i,k)
            enddo
         enddo
!         endif
      enddo

!-------------------------------------------------------------
! map virtual temperature using geopotential conserving scheme.
!-------------------------------------------------------------
      call mappm(km, pn0, tp, npz, pn1, qn1, is,ie, 1, 9, Atm%ptop)
      do k=1,npz
         do i=is,ie
            Atm%pt(i,j,k) = qn1(i,k)/(1.+zvir*Atm%q(i,j,k,sphum))
         enddo
      enddo

5000 continue

  call prt_maxmin('PS_model', Atm%ps, is, ie, js, je, ng, 1, 0.01_FVPRC)
  call prt_maxmin('UT', ut, is, ie, js, je, ng, npz, 1._FVPRC)
  call prt_maxmin('VT', vt, is, ie, js, je, ng, npz, 1._FVPRC)

!----------------------------------------------
! winds: lat-lon ON A to Cubed-D transformation:
!----------------------------------------------
  call cubed_a2d(Atm%npx, Atm%npy, npz, ut, vt, Atm%u, Atm%v, Atm%gridstruct, Atm%domain, Atm%bd )

  if (is_master()) write(*,*) 'done remap_xyz'

 end subroutine remap_xyz



                     subroutine init_cubsph_grid(npts, is,ie, js,je, ntiles, sph_corner)  
!------------------------------------------------------------------!
! read/generate cubed sphere grid                                  !
! calculate cell center from cell corners                          !
!                                                                  !
! input:                                                           !
! npts, is,ie, js,je, ntiles       number of grid points and tiles !
!                                                                  !
! output:                                                          !
! sph_corner             cell corners in spherical coor            !
!------------------------------------------------------------------!
                        use GHOST_CUBSPH_mod, only: B_grid, ghost_cubsph_update             
                        use fv_grid_utils_mod, only : gnomonic_grids
                        use fv_grid_tools_mod, only : mirror_grid

                        integer, intent(in) :: npts, is,ie, js,je, ntiles
                        real*8, dimension(2,is:ie+1,js:je+1), intent(out) :: sph_corner
!------------------------------------------------------------------!
! local variables                                                  !
!------------------------------------------------------------------!
                        integer :: i, j, l, n
                        real*8, pointer :: xs(:,:), ys(:,:)
                        real*8, pointer :: grid_in(:,:,:,:)
                        integer :: grid_type = 0
!------------------------------------------------------------------!
! create sph_corner                                                !
!------------------------------------------------------------------!
#ifdef SMEM_MAPL_MODE
! allocate global arrays (preferable in shared memory)
                        if(MAPL_ShmInitialized) then
                           if (is_master()) write(*,*) 'Using MAPL_Shmem in external_ic: init_cubsph_grid' 
                           call MAPL_AllocNodeArray(grid_in,Shp=(/npts,npts,2,ntiles/),rc=STATUS)
                        else
                           if (is_master()) write(*,*) 'WARNING... in external_ic: Global grid allocate'
                           allocate( grid_in(npts,npts,2,ntiles) )
                        endif
                        if (is_master()) then
                          allocate( xs(npts,npts) )
                          allocate( ys(npts,npts) )
                          call gnomonic_grids(grid_type, npts-1, xs, ys)
                          do j=1,npts
                             do i=1,npts
                                grid_in(i,j,1,1) = xs(i,j)
                                grid_in(i,j,2,1) = ys(i,j)
                             enddo
                          enddo
                          deallocate ( xs )
                          deallocate ( ys )
                        endif
                        if(MAPL_ShmInitialized) then
                          call MAPL_SyncSharedMemory(rc=STATUS)
                          call MAPL_BroadcastToNodes( grid_in, N=size(grid_in), ROOT=masterproc, RC=status)
                          call MAPL_SyncSharedMemory(rc=STATUS)
                        else
                          call mpp_broadcast(grid_in, size(grid_in), masterproc)
                        endif
#else
                        allocate( xs(npts,npts) )
                        allocate( ys(npts,npts) )
                        allocate( grid_in(npts,npts,2,ntiles) )
                        call gnomonic_grids(grid_type, npts-1, xs, ys)
                        do j=1,npts
                           do i=1,npts
                              grid_in(i,j,1,1) = xs(i,j)
                              grid_in(i,j,2,1) = ys(i,j)
                           enddo
                        enddo
                        deallocate ( xs )
                        deallocate ( ys )
#endif

! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi]
                        call mirror_grid(grid_in, 0, npts, npts, 2, ntiles)
                        sph_corner(1,is:ie+1,js:je+1) = grid_in(is:ie+1,js:je+1,1,tile)
                        sph_corner(2,is:ie+1,js:je+1) = grid_in(is:ie+1,js:je+1,2,tile)
                        do j=js,je+1
                           do i=is,ie+1
!---------------------------------
! Shift the corner away from Japan
!---------------------------------
! This will result in the corner close to east coast of China
                              sph_corner(1,i,j) = sph_corner(1,i,j) - pi/18.
                              if ( sph_corner(1,i,j) < 0. )              &
                                    sph_corner(1,i,j) = sph_corner(1,i,j) + 2.*pi
                              if (ABS(sph_corner(1,i,j)) < 1.e-10) sph_corner(1,i,j) = 0.0
                              if (ABS(sph_corner(2,i,j)) < 1.e-10) sph_corner(2,i,j) = 0.0
                           enddo
                        enddo
#ifdef SMEM_MAPL_MODE
                        call MAPL_SyncSharedMemory(rc=STATUS)
                        DEALLOCGLOB_(grid_in)
                        call MAPL_SyncSharedMemory(rc=STATUS)
#else
                        deallocate ( grid_in )
#endif
!------------------------------------------------------------------!
! do halo update                                                   !
!------------------------------------------------------------------!

                     end subroutine init_cubsph_grid

                     subroutine interp_c2c_vect(npx_in, npy_in, npx_out, npy_out, npz, ntiles, domain_i, &
                           is,ie, js,je, isd_i,ied_i, jsd_i,jed_i, is_i,ie_i, js_i,je_i, &
                           u_in, v_in, u_out, v_out, index_c2c, weight_c2c, corner_in, corner_out, gridstruct)
                        use GRID_UTILS_mod, only: latlon2xyz
                        use GRID_UTILS_mod,   only: get_dx, get_dxa, get_dy, get_dya,     &
                              get_center_vect, get_west_vect,       &
                              get_south_vect, get_cosa_center
                        use FLOW_PROJ_mod,    only: d2a, d2a_vect, a2d_vect
                        use GHOST_CUBSPH_mod,  only : A_grid, ghost_cubsph_update
                        use CUB2CUB_mod,    only: do_c2c_interpolation
                        integer, intent(IN) :: npx_in, npy_in, npx_out, npy_out, npz, ntiles
                        type(domain2d), intent(INOUT) :: domain_i
                        integer, intent(IN) :: is,ie, js,je, isd_i,ied_i, jsd_i,jed_i, is_i,ie_i, js_i,je_i
                        integer, intent(IN) ::  index_c2c(3, is:ie, js:je )
                        real(REAL64),  intent(IN) :: weight_c2c(4, is:ie, js:je )
                        real(REAL64),  intent(IN) :: corner_in(2,is_i:ie_i+1,js_i:je_i+1)
                        real(REAL64),  intent(IN) :: corner_out(2,is:ie+1,js:je+1)
                        real(FVPRC),  intent(IN) ::  u_in(isd_i:ied_i,jsd_i:jed_i+1,npz)
                        real(FVPRC),  intent(IN) ::  v_in(isd_i:ied_i+1,jsd_i:jed_i,npz)
                        real(FVPRC),  intent(OUT):: u_out(is:ie,js:je,npz)
                        real(FVPRC),  intent(OUT):: v_out(is:ie,js:je,npz)
                        type(fv_grid_type), intent(IN), target :: gridstruct

                        real(FVPRC) :: tmp(isd_i:ied_i,jsd_i:jed_i)

                        integer :: i,j,l,k,n
                        real(REAL64), dimension(:,:,:), allocatable :: vxyz_in, vxyz_out
                        real(REAL64), dimension(:,:,:), allocatable :: ec1, ec2, ew1, ew2, es1, es2
                        real(REAL64), dimension(:,:)  , allocatable :: dx, dy, dxa, dya, rdxa, rdya, cosa_s, sina_s
                        real(FVPRC) :: u1, v1, vx, vy, vz
                        real(REAL64) :: pc1(3), pc2(3)
                        integer :: ic, jc, lc

                        real(REAL64), dimension(:,:,:), allocatable :: xyz_corner_in, xyz_corner_out

                        character(len=64) :: strTxt

!------------------------------------------------------------------!
! calculate xyz cell corners and cell centers                      !
!------------------------------------------------------------------!
                        tmp = 0.0
                        allocate(xyz_corner_in (3, isd_i:ied_i+1, jsd_i:jed_i+1), &
                              xyz_corner_out(3, is   :ie   +1, js   :je   +1))
                        do j=js_i,je_i+1
                           do i=is_i,ie_i+1
                              call latlon2xyz(corner_in(:,i,j), xyz_corner_in(:,i,j))
                           enddo
                        enddo
                        do j=js,je+1
                           do i=is,ie+1
                              call latlon2xyz(corner_out(:,i,j), xyz_corner_out(:,i,j))
                           enddo
                        enddo
                        call print_memuse_stats('interp_c2c_vect: GRIDS')

!----------------------------------------------------------!
! allocate horizontal flow variables                       !
!----------------------------------------------------------!
                        allocate(vxyz_in(3,isd_i:ied_i,jsd_i:jed_i))
                        allocate(dx(isd_i:ied_i,jsd_i:jed_i+1), dxa(isd_i:ied_i,jsd_i:jed_i), rdxa(isd_i:ied_i,jsd_i:jed_i))
                        allocate(dy(isd_i:ied_i+1,jsd_i:jed_i), dya(isd_i:ied_i,jsd_i:jed_i), rdya(isd_i:ied_i,jsd_i:jed_i))
                        allocate(ec1(3,isd_i:ied_i,jsd_i:jed_i), ec2(3,isd_i:ied_i,jsd_i:jed_i))
                        allocate(cosa_s(isd_i:ied_i,jsd_i:jed_i), sina_s(isd_i:ied_i,jsd_i:jed_i))
                        allocate(vxyz_out(3,is:ie,js:je))

!-------------------------------------------------------!
! geometrical properties of input grid                  !
!-------------------------------------------------------!
                        call get_dx (xyz_corner_in(:,isd_i:ied_i+1,jsd_i:jed_i+1), &
                              isd_i,ied_i  ,jsd_i,jed_i, &
                              is_i ,ie_i   ,js_i ,je_i , dx)
                        call get_dxa(xyz_corner_in(:,isd_i:ied_i+1,jsd_i:jed_i+1), &
                              isd_i,ied_i  ,jsd_i,jed_i, &
                              is_i ,ie_i   ,js_i ,je_i , dxa, rdxa=rdxa)
                        call get_dy (xyz_corner_in(:,isd_i:ied_i+1,jsd_i:jed_i+1), &
                              isd_i,ied_i  ,jsd_i,jed_i, &
                              is_i ,ie_i   ,js_i ,je_i , dy)
                        call get_dya(xyz_corner_in(:,isd_i:ied_i+1,jsd_i:jed_i+1), &
                              isd_i,ied_i  ,jsd_i,jed_i, &
                              is_i ,ie_i   ,js_i ,je_i , dya, rdya=rdya)
                        call get_center_vect(xyz_corner_in(:,isd_i:ied_i+1,jsd_i:jed_i+1), &
                              isd_i,ied_i  ,jsd_i,jed_i, &
                              is_i ,ie_i   ,js_i ,je_i , ec1, ec2)
                        call get_cosa_center(ec1, ec2, isd_i,ied_i  ,jsd_i,jed_i, &
                              is_i ,ie_i   ,js_i ,je_i , cosa_s, sina_s)

! Flow interpolation for U and V components
                        do k=1,npz
!-------------------------------------------------------!
! calculate flow vector for a-grid                      !
!-------------------------------------------------------!
                           call d2a_vect(DBLE(u_in(:,:,k)), DBLE(v_in(:,:,k)), DBLE(dx), DBLE(dy), DBLE(rdxa), DBLE(rdya), DBLE(cosa_s), DBLE(ec1), DBLE(ec2), &
                                 isd_i, ied_i, jsd_i, jed_i, 1, 1,           &
                                 is_i , ie_i , js_i , je_i , 1, 1,           &
                                 vxyz_in(:,isd_i:ied_i,jsd_i:jed_i))
                           write(strTxt,'(A,i3.3)') 'interp_c2c_vect: INPUT D2A level:', k
                           if (k==npz) call print_memuse_stats(strTxt)
!----------------------------------------------------------!
! ghost cell update of vxyz_in                               !
!----------------------------------------------------------!
                           do n=1,3
                              tmp(is_i:ie_i,js_i:je_i) = vxyz_in(n,is_i:ie_i,js_i:je_i)
                              call mpp_update_domains(tmp, domain_i)
                              vxyz_in(n,:,:) = tmp
                           enddo
                           do n=1,3
                              do j=js,je
                                 do i=is,ie
!----------------------------------------------------------!
! do interpolation of flow vector on A-Grid                !
!----------------------------------------------------------!
                                    ic=index_c2c(1,i,j)
                                    jc=index_c2c(2,i,j)
                                    vxyz_out(n,i,j)=weight_c2c(1,i,j)*vxyz_in(n,ic  ,jc  )  &
                                          +weight_c2c(2,i,j)*vxyz_in(n,ic  ,jc+1)  &
                                          +weight_c2c(3,i,j)*vxyz_in(n,ic+1,jc+1)  &
                                          +weight_c2c(4,i,j)*vxyz_in(n,ic+1,jc  )
                                 enddo
                              enddo
                           enddo
                           do j=js,je
                              do i=is,ie
                                 vx = vxyz_out(1,i,j)
                                 vy = vxyz_out(2,i,j)
                                 vz = vxyz_out(3,i,j)
!----------------------------------------------------------!
! convert flow vector to wind vectors on A-Grid            !
!----------------------------------------------------------!
                                 pc1(:)=xyz_corner_out(:,i+1,j)+xyz_corner_out(:,i+1,j+1)          &
                                       -xyz_corner_out(:,i  ,j)-xyz_corner_out(:,i  ,j+1)
                                 pc2(:)=xyz_corner_out(:,i,j+1)+xyz_corner_out(:,i+1,j+1)          &
                                       -xyz_corner_out(:,i,j  )-xyz_corner_out(:,i+1,j  )
                                 call normalize_vect(pc1(:))
                                 call normalize_vect(pc2(:))
                                 u1 = vx*pc1(1) + vy*pc1(2) + vz*pc1(3)
                                 v1 = vx*pc2(1) + vy*pc2(2) + vz*pc2(3)
!----------------------------------------------------------!
! rotate wind vectors from cubed to latlon orientation     !
!----------------------------------------------------------!
                                 u_out(i,j,k) = 2.0*(gridstruct%a11(i,j)*u1 + gridstruct%a12(i,j)*v1)
                                 v_out(i,j,k) = 2.0*(gridstruct%a21(i,j)*u1 + gridstruct%a22(i,j)*v1)
                              enddo
                           enddo
                           write(strTxt,'(A,i3.3)') 'interp_c2c_vect: OUTPUT A-grid C2L level:', k
                           if (k==npz) call print_memuse_stats(strTxt)
                        enddo ! npz

                        deallocate(dx, dy, dxa, dya, rdxa, rdya, ec1, ec2, cosa_s, sina_s)
                        deallocate(vxyz_in)
                        deallocate(vxyz_out)
                        deallocate ( xyz_corner_in, xyz_corner_out )

                     end subroutine interp_c2c_vect

                     Function REVERSE(A) Result(B)
                        real(FVPRC), Intent(In) :: A(:,:)
                        real(FVPRC) :: B(Size(A,1),Size(A,2))

                        Integer :: i, n

                        n = Size(A, 1)

                        Do i = 1, n
                           B(i,:) = A(1+n-i,:)
                        End Do

                     End Function REVERSE

 subroutine cubed_a2d( npx, npy, npz, ua, va, u, v, gridstruct, fv_domain, bd )

! Purpose; Transform wind on A grid to D grid

  use mpp_domains_mod,    only: mpp_update_domains

  type(fv_grid_bounds_type), intent(IN) :: bd
  integer, intent(in):: npx, npy, npz
  real(FVPRC), intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz):: ua, va
  real(FVPRC), intent(out):: u(bd%isd:bd%ied,  bd%jsd:bd%jed+1,npz)
  real(FVPRC), intent(out):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz)
  type(fv_grid_type), intent(IN), target :: gridstruct
  type(domain2d), intent(INOUT) :: fv_domain
! local:
  real(FVPRC) v3(3,bd%is-1:bd%ie+1,bd%js-1:bd%je+1)
  real(FVPRC) ue(3,bd%is-1:bd%ie+1,bd%js:bd%je+1)    ! 3D winds at edges
  real(FVPRC) ve(3,bd%is:bd%ie+1,bd%js-1:bd%je+1)    ! 3D winds at edges
  real(FVPRC), dimension(bd%is:bd%ie):: ut1, ut2, ut3
  real(FVPRC), dimension(bd%js:bd%je):: vt1, vt2, vt3
  integer i, j, k, im2, jm2

  real(REAL64), pointer, dimension(:,:,:)   :: vlon, vlat
  real(REAL64), pointer, dimension(:)       :: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n
  real(REAL64), pointer, dimension(:,:,:,:) :: ew, es

  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = bd%is
  ie  = bd%ie
  js  = bd%js
  je  = bd%je
  isd = bd%isd
  ied = bd%ied
  jsd = bd%jsd
  jed = bd%jed
  vlon => gridstruct%vlon
  vlat => gridstruct%vlat

  edge_vect_w => gridstruct%edge_vect_w
  edge_vect_e => gridstruct%edge_vect_e
  edge_vect_s => gridstruct%edge_vect_s
  edge_vect_n => gridstruct%edge_vect_n

  ew => gridstruct%ew
  es => gridstruct%es

  call mpp_update_domains(ua, fv_domain, complete=.false.)
  call mpp_update_domains(va, fv_domain, complete=.true.)

    im2 = (npx-1)/2
    jm2 = (npy-1)/2
    do k=1, npz
! Compute 3D wind on A grid
       do j=js-1,je+1
          do i=is-1,ie+1
             !v3(1,i,j) = ua(i,j,k)*vlon(1,i,j) + va(i,j,k)*vlat(1,i,j)
             !v3(2,i,j) = ua(i,j,k)*vlon(2,i,j) + va(i,j,k)*vlat(2,i,j)
             !v3(3,i,j) = ua(i,j,k)*vlon(3,i,j) + va(i,j,k)*vlat(3,i,j)
             v3(1,i,j) = ua(i,j,k)*vlon(i,j,1) + va(i,j,k)*vlat(i,j,1)
             v3(2,i,j) = ua(i,j,k)*vlon(i,j,2) + va(i,j,k)*vlat(i,j,2)
             v3(3,i,j) = ua(i,j,k)*vlon(i,j,3) + va(i,j,k)*vlat(i,j,3)
          enddo
       enddo

! A --> D
! Interpolate to cell edges
       do j=js,je+1
          do i=is-1,ie+1
             ue(1,i,j) = 0.5*(v3(1,i,j-1) + v3(1,i,j))
             ue(2,i,j) = 0.5*(v3(2,i,j-1) + v3(2,i,j))
             ue(3,i,j) = 0.5*(v3(3,i,j-1) + v3(3,i,j))
          enddo
       enddo

       do j=js-1,je+1
          do i=is,ie+1
             ve(1,i,j) = 0.5*(v3(1,i-1,j) + v3(1,i,j))
             ve(2,i,j) = 0.5*(v3(2,i-1,j) + v3(2,i,j))
             ve(3,i,j) = 0.5*(v3(3,i-1,j) + v3(3,i,j))
          enddo
       enddo

! --- E_W edges (for v-wind):
     if (.not. gridstruct%nested) then
     if ( is==1) then
       i = 1
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_w(j)*ve(1,i,j-1)+(1.-edge_vect_w(j))*ve(1,i,j)
             vt2(j) = edge_vect_w(j)*ve(2,i,j-1)+(1.-edge_vect_w(j))*ve(2,i,j)
             vt3(j) = edge_vect_w(j)*ve(3,i,j-1)+(1.-edge_vect_w(j))*ve(3,i,j)
        else
             vt1(j) = edge_vect_w(j)*ve(1,i,j+1)+(1.-edge_vect_w(j))*ve(1,i,j)
             vt2(j) = edge_vect_w(j)*ve(2,i,j+1)+(1.-edge_vect_w(j))*ve(2,i,j)
             vt3(j) = edge_vect_w(j)*ve(3,i,j+1)+(1.-edge_vect_w(j))*ve(3,i,j)
        endif
       enddo
       do j=js,je
          ve(1,i,j) = vt1(j)
          ve(2,i,j) = vt2(j)
          ve(3,i,j) = vt3(j)
       enddo
     endif

     if ( (ie+1)==npx ) then
       i = npx
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = edge_vect_e(j)*ve(1,i,j-1)+(1.-edge_vect_e(j))*ve(1,i,j)
             vt2(j) = edge_vect_e(j)*ve(2,i,j-1)+(1.-edge_vect_e(j))*ve(2,i,j)
             vt3(j) = edge_vect_e(j)*ve(3,i,j-1)+(1.-edge_vect_e(j))*ve(3,i,j)
        else
             vt1(j) = edge_vect_e(j)*ve(1,i,j+1)+(1.-edge_vect_e(j))*ve(1,i,j)
             vt2(j) = edge_vect_e(j)*ve(2,i,j+1)+(1.-edge_vect_e(j))*ve(2,i,j)
             vt3(j) = edge_vect_e(j)*ve(3,i,j+1)+(1.-edge_vect_e(j))*ve(3,i,j)
        endif
       enddo
       do j=js,je
          ve(1,i,j) = vt1(j)
          ve(2,i,j) = vt2(j)
          ve(3,i,j) = vt3(j)
       enddo
     endif

! N-S edges (for u-wind):
     if ( js==1 ) then
       j = 1
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_s(i)*ue(1,i-1,j)+(1.-edge_vect_s(i))*ue(1,i,j)
             ut2(i) = edge_vect_s(i)*ue(2,i-1,j)+(1.-edge_vect_s(i))*ue(2,i,j)
             ut3(i) = edge_vect_s(i)*ue(3,i-1,j)+(1.-edge_vect_s(i))*ue(3,i,j)
        else
             ut1(i) = edge_vect_s(i)*ue(1,i+1,j)+(1.-edge_vect_s(i))*ue(1,i,j)
             ut2(i) = edge_vect_s(i)*ue(2,i+1,j)+(1.-edge_vect_s(i))*ue(2,i,j)
             ut3(i) = edge_vect_s(i)*ue(3,i+1,j)+(1.-edge_vect_s(i))*ue(3,i,j)
        endif
       enddo
       do i=is,ie
          ue(1,i,j) = ut1(i)
          ue(2,i,j) = ut2(i)
          ue(3,i,j) = ut3(i)
       enddo
     endif

     if ( (je+1)==npy ) then
       j = npy
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = edge_vect_n(i)*ue(1,i-1,j)+(1.-edge_vect_n(i))*ue(1,i,j)
             ut2(i) = edge_vect_n(i)*ue(2,i-1,j)+(1.-edge_vect_n(i))*ue(2,i,j)
             ut3(i) = edge_vect_n(i)*ue(3,i-1,j)+(1.-edge_vect_n(i))*ue(3,i,j)
        else
             ut1(i) = edge_vect_n(i)*ue(1,i+1,j)+(1.-edge_vect_n(i))*ue(1,i,j)
             ut2(i) = edge_vect_n(i)*ue(2,i+1,j)+(1.-edge_vect_n(i))*ue(2,i,j)
             ut3(i) = edge_vect_n(i)*ue(3,i+1,j)+(1.-edge_vect_n(i))*ue(3,i,j)
        endif
       enddo
       do i=is,ie
          ue(1,i,j) = ut1(i)
          ue(2,i,j) = ut2(i)
          ue(3,i,j) = ut3(i)
       enddo
     endif

     endif ! .not. nested
     do j=js,je+1
        do i=is,ie
           u(i,j,k) =  ue(1,i,j)*es(1,i,j,1) +  &
                       ue(2,i,j)*es(2,i,j,1) +  &
                       ue(3,i,j)*es(3,i,j,1)
        enddo
     enddo
     do j=js,je
        do i=is,ie+1
           v(i,j,k) = ve(1,i,j)*ew(1,i,j,2) +  &
                      ve(2,i,j)*ew(2,i,j,2) +  &
                      ve(3,i,j)*ew(3,i,j,2)
        enddo
     enddo

   enddo         ! k-loop

 end subroutine cubed_a2d


                     subroutine d2a3d(u, v,  ua,   va,  im,  jm, km, lon)
                        integer, intent(in):: im, jm, km           ! Dimensions
                        real(FVPRC), intent(in ) :: lon(im)
                        real(FVPRC), intent(in ), dimension(im,jm,km):: u, v
                        real(FVPRC), intent(out), dimension(im,jm,km):: ua, va
! local
                        real(FVPRC) :: coslon(im),sinlon(im)    ! Sine and cosine in longitude
                        integer i, j, k
                        integer imh
                        real(FVPRC) un, vn, us, vs

                        integer :: ks, ke

                        imh = im/2

                        do i=1,im
                           sinlon(i) = sin(lon(i))
                           coslon(i) = cos(lon(i))
                        enddo

                        do k=1,km
                           do j=2,jm-1
                              do i=1,im
                                 ua(i,j,k) = 0.5*(u(i,j,k) + u(i,j+1,k))
                              enddo
                           enddo

                           do j=2,jm-1
                              do i=1,im-1
                                 va(i,j,k) = 0.5*(v(i,j,k) + v(i+1,j,k))
                              enddo
                              va(im,j,k) = 0.5*(v(im,j,k) + v(1,j,k))
                           enddo

! Projection at SP
                           us = 0.
                           vs = 0.
                           do i=1,imh
                              us = us + (ua(i+imh,2,k)-ua(i,2,k))*sinlon(i)      &
                                    + (va(i,2,k)-va(i+imh,2,k))*coslon(i)
                              vs = vs + (ua(i+imh,2,k)-ua(i,2,k))*coslon(i)      &
                                    + (va(i+imh,2,k)-va(i,2,k))*sinlon(i)
                           enddo
                           us = us/im
                           vs = vs/im
                           do i=1,imh
                              ua(i,1,k)   = -us*sinlon(i) - vs*coslon(i)
                              va(i,1,k)   =  us*coslon(i) - vs*sinlon(i)
                              ua(i+imh,1,k)   = -ua(i,1,k)
                              va(i+imh,1,k)   = -va(i,1,k)
                           enddo

! Projection at NP
                           un = 0.
                           vn = 0.
                           do i=1,imh
                              un = un + (ua(i+imh,jm-1,k)-ua(i,jm-1,k))*sinlon(i)    &
                                    + (va(i+imh,jm-1,k)-va(i,jm-1,k))*coslon(i)
                              vn = vn + (ua(i,jm-1,k)-ua(i+imh,jm-1,k))*coslon(i)    &
                                    + (va(i+imh,jm-1,k)-va(i,jm-1,k))*sinlon(i)
                           enddo

                           un = un/im
                           vn = vn/im
                           do i=1,imh
                              ua(i,jm,k) = -un*sinlon(i) + vn*coslon(i)
                              va(i,jm,k) = -un*coslon(i) - vn*sinlon(i)
                              ua(i+imh,jm,k) = -ua(i,jm,k)
                              va(i+imh,jm,k) = -va(i,jm,k)
                           enddo
                        enddo

                     end subroutine d2a3d



                     subroutine pmaxmin( qname, a, im, jm, fac )

                        integer, intent(in):: im, jm
                        character(len=*) :: qname
                        integer i, j
                        real(FVPRC) a(im,jm)

                        real(FVPRC) qmin(jm), qmax(jm)
                        real(FVPRC) pmax, pmin
                        real(FVPRC) fac                     ! multiplication factor

                        do j=1,jm
                           pmax = a(1,j)
                           pmin = a(1,j)
                           do i=2,im
                              pmax = max(pmax, a(i,j))
                              pmin = min(pmin, a(i,j))
                           enddo
                           qmax(j) = pmax
                           qmin(j) = pmin
                        enddo
!
! Now find max/min of amax/amin
!
                        pmax = qmax(1)
                        pmin = qmin(1)
                        do j=2,jm
                           pmax = max(pmax, qmax(j))
                           pmin = min(pmin, qmin(j))
                        enddo

                        write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

                     end subroutine pmaxmin

                     subroutine pmaxmin4d( qname, a, im, jm, km, lm, fac )

                        character*(*)  qname 
                        integer, intent(in):: im, jm, km, lm
                        integer i, j, k, l
                        real(FVPRC) a(im,jm,km,lm)

                        real(FVPRC) qmin(jm), qmax(jm)
                        real(FVPRC) pmax, pmin
                        real(FVPRC) fac                     ! multiplication factor

                        pmax = a(1,1,1,1)
                        pmin = a(1,1,1,1)
                        do l=1,lm
                           do k=1,km
                              do j=1,jm
                                 do i=1,im
                                    pmax = max(pmax, a(i,j,k,l))
                                    pmin = min(pmin, a(i,j,k,l))
                                 enddo
                              enddo
                           enddo
                        enddo

                        write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

                     end subroutine pmaxmin4d


                     subroutine mpp_domain_decomp(domain,npx,npy,nregions,ng,grid_type, &
                           is,ie,js,je,isd,ied,jsd,jed,tile)
                        use mpp_domains_mod, only: mpp_domains_init, MPP_DOMAIN_TIME, &
                              mpp_define_mosaic, mpp_get_compute_domain, mpp_get_data_domain, &
                              mpp_domains_set_stack_size, mpp_define_layout
                        use mpp_mod,         only : mpp_pe
                        type(domain2D), intent(OUT) :: domain
                        integer, intent(IN)  :: npx,npy,nregions,ng,grid_type
                        integer, intent(OUT) :: is,ie,js,je,isd,ied,jsd,jed,tile

                        integer :: layout(2)
                        integer, allocatable :: pe_start(:), pe_end(:)

                        integer :: num_contact, ntiles, npes_per_tile
                        integer, allocatable, dimension(:)       :: npes_tile, tile1, tile2
                        integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
                        integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
                        integer, allocatable, dimension(:,:)     :: layout2D, global_indices

                        character*80 :: evalue
                        integer :: ios,nx,ny,n,num_alloc
                        character(len=32) :: type = "unknown"
                        logical :: is_symmetry
                        logical :: debug=.false.
                        integer :: npes

                        npes = mpp_npes()

                        nx = npx-1
                        ny = npy-1

                        call print_memuse_stats('external_ic:mpp_domain_decomp: top')

                        call mpp_domains_init(MPP_DOMAIN_TIME)

                        call print_memuse_stats('external_ic:mpp_domain_decomp: mpp_domains_init')

! call mpp_domains_set_stack_size(10000)
! call mpp_domains_set_stack_size(900000)
! call mpp_domains_set_stack_size(1500000)
                        call mpp_domains_set_stack_size(3000000)

                        select case(nregions)
                        case ( 1 )  ! Lat-Lon "cyclic"

                           select case (grid_type)
                           case (3)   ! Lat-Lon "cyclic"
                              type="Lat-Lon: cyclic"
                              ntiles = 4
                              num_contact = 8
                              if( mod(npes,ntiles) .NE. 0 ) then
                                 call mpp_error(NOTE,'TEST_MPP_DOMAINS: for Cyclic mosaic, npes should be multiple of ntiles. ' // &
                                       'No test is done for Cyclic mosaic. ' )
                                 return
                              end if
                              npes_per_tile = npes/ntiles
                              call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
                              layout = (/1,npes_per_tile/) ! force decomp only in Lat-Direction
                           case (4)   ! Cartesian, double periodic
                              type="Cartesian: double periodic"
                              ntiles = 1
                              num_contact = 2
                              npes_per_tile = npes/ntiles
                              call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
                           case (5)   ! latlon patch
                              type="Lat-Lon: patch"
                              ntiles = 1
                              num_contact = 0
                              npes_per_tile = npes/ntiles
                              call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
                           case (6)   ! latlon strip
                              type="Lat-Lon: strip"
                              ntiles = 1
                              num_contact = 1
                              npes_per_tile = npes/ntiles
                              call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
                           case (7)   ! Cartesian, channel
                              type="Cartesian: channel"
                              ntiles = 1
                              num_contact = 1
                              npes_per_tile = npes/ntiles
                              call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
                           end select

                        case ( 6 )  ! Cubed-Sphere
                           type="Cubic: cubed-sphere"
                           ntiles = 6
                           num_contact = 12
!--- cubic grid always have six tiles, so npes should be multiple of 6
                           if( mod(npes,ntiles) .NE. 0 .OR. npx-1 .NE. npy-1) then
                              call mpp_error(NOTE,'mpp_domain_decomp: for Cubic_grid mosaic, npes should be multiple of ntiles(6) ' // &
                                    'and npx-1 should equal npy-1, mpp_domain_decomp is NOT done for Cubic-grid mosaic. ' )
                              return
                           end if
                           npes_per_tile = npes/ntiles
                           call  mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
                           npes_x = layout(1)
                           npes_y = layout(2)
                           if ( (npx/npes_x < ng) .or. (npy/npes_y < ng) ) then
                              write(*,310) npes_x, npes_y, npx/npes_x, npy/npes_y
310                           format('Invalid layout, NPES_X:',i4.4,'NPES_Y:',i4.4,'ncells_X:',i4.4,'ncells_Y:',i4.4)
                              call mpp_error(FATAL, 'mpp_domain_decomp: bad decomp')
                           endif
                        case default
                           call mpp_error(FATAL, 'mpp_domain_decomp: no such test: '//type)
                        end select

                        call print_memuse_stats('external_ic:mpp_domain_decomp: mpp_define_layout')

                        allocate(layout2D(2,ntiles), global_indices(4,ntiles), npes_tile(ntiles) )
                        allocate(pe_start(ntiles),pe_end(ntiles))
                        npes_tile = npes_per_tile
                        do n = 1, ntiles
                           global_indices(:,n) = (/1,npx-1,1,npy-1/)
                           layout2D(:,n)         = layout
                           pe_start(n) = (n-1)*layout(1)*layout(2)
                           pe_end(n)   = pe_start(n) + layout(1)*layout(2) -1
                        end do
                        num_alloc=max(1,num_contact)
                        allocate(tile1(num_alloc), tile2(num_alloc) )
                        allocate(istart1(num_alloc), iend1(num_alloc), jstart1(num_alloc), jend1(num_alloc) )
                        allocate(istart2(num_alloc), iend2(num_alloc), jstart2(num_alloc), jend2(num_alloc) )

                        is_symmetry = .true.

                        call print_memuse_stats('external_ic:mpp_domain_decomp: allocates 1')

                        select case(nregions)
                        case ( 1 )

                           select case (grid_type)
                           case (3)   ! Lat-Lon "cyclic"
!--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
                              tile1(1) = 1; tile2(1) = 2
                              istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
                              istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
!--- Contact line 2, between tile 1 (SOUTH) and tile 3 (NORTH)  --- cyclic
                              tile1(2) = 1; tile2(2) = 3
                              istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
                              istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
!--- Contact line 3, between tile 1 (WEST) and tile 2 (EAST) --- cyclic
                              tile1(3) = 1; tile2(3) = 2
                              istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
                              istart2(3) = nx; iend2(3) = nx; jstart2(3) = 1;  jend2(3) = ny
!--- Contact line 4, between tile 1 (NORTH) and tile 3 (SOUTH)
                              tile1(4) = 1; tile2(4) = 3
                              istart1(4) = 1;  iend1(4) = nx; jstart1(4) = ny;  jend1(4) = ny
                              istart2(4) = 1;  iend2(4) = nx; jstart2(4) = 1;   jend2(4) = 1
!--- Contact line 5, between tile 2 (SOUTH) and tile 4 (NORTH) --- cyclic
                              tile1(5) = 2; tile2(5) = 4
                              istart1(5) = 1;  iend1(5) = nx; jstart1(5) = 1;  jend1(5) = 1
                              istart2(5) = 1;  iend2(5) = nx; jstart2(5) = ny; jend2(5) = ny
!--- Contact line 6, between tile 2 (NORTH) and tile 4 (SOUTH)
                              tile1(6) = 2; tile2(6) = 4
                              istart1(6) = 1;  iend1(6) = nx; jstart1(6) = ny;  jend1(6) = ny
                              istart2(6) = 1;  iend2(6) = nx; jstart2(6) = 1;   jend2(6) = 1
!--- Contact line 7, between tile 3 (EAST) and tile 4 (WEST)
                              tile1(7) = 3; tile2(7) = 4
                              istart1(7) = nx; iend1(7) = nx; jstart1(7) = 1;  jend1(7) = ny
                              istart2(7) = 1;  iend2(7) = 1;  jstart2(7) = 1;  jend2(7) = ny
!--- Contact line 8, between tile 3 (WEST) and tile 4 (EAST) --- cyclic
                              tile1(8) = 3; tile2(8) = 4
                              istart1(8) = 1;  iend1(8) = 1;  jstart1(8) = 1;  jend1(8) = ny
                              istart2(8) = nx; iend2(8) = nx; jstart2(8) = 1;  jend2(8) = ny
                              is_symmetry = .false.
                           case (4)   ! Cartesian, double periodic
!--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
                              tile1(1) = 1; tile2(1) = 1
                              istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
                              istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
!--- Contact line 2, between tile 1 (SOUTH) and tile 1 (NORTH)  --- cyclic
                              tile1(2) = 1; tile2(2) = 1
                              istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
                              istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
                           case (5)   ! latlon patch

                           case (6)   !latlon strip
!--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
                              tile1(1) = 1; tile2(1) = 1
                              istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
                              istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
                           case (7)   ! Cartesian, channel
!--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
                              tile1(1) = 1; tile2(1) = 1
                              istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
                              istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
                           end select

                        case ( 6 )  ! Cubed-Sphere
!--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
                           tile1(1) = 1; tile2(1) = 2
                           istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
                           istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
!--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
                           tile1(2) = 1; tile2(2) = 3
                           istart1(2) = 1;  iend1(2) = nx; jstart1(2) = ny; jend1(2) = ny
                           istart2(2) = 1;  iend2(2) = 1;  jstart2(2) = ny; jend2(2) = 1
!--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
                           tile1(3) = 1; tile2(3) = 5
                           istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
                           istart2(3) = nx; iend2(3) = 1;  jstart2(3) = ny; jend2(3) = ny
!--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
                           tile1(4) = 1; tile2(4) = 6
                           istart1(4) = 1;  iend1(4) = nx; jstart1(4) = 1;  jend1(4) = 1
                           istart2(4) = 1;  iend2(4) = nx; jstart2(4) = ny; jend2(4) = ny
!--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
                           tile1(5) = 2; tile2(5) = 3
                           istart1(5) = 1;  iend1(5) = nx; jstart1(5) = ny; jend1(5) = ny
                           istart2(5) = 1;  iend2(5) = nx; jstart2(5) = 1;  jend2(5) = 1
!--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
                           tile1(6) = 2; tile2(6) = 4
                           istart1(6) = nx; iend1(6) = nx; jstart1(6) = 1;  jend1(6) = ny
                           istart2(6) = nx; iend2(6) = 1;  jstart2(6) = 1;  jend2(6) = 1
!--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
                           tile1(7) = 2; tile2(7) = 6
                           istart1(7) = 1;  iend1(7) = nx; jstart1(7) = 1;  jend1(7) = 1
                           istart2(7) = nx; iend2(7) = nx; jstart2(7) = ny; jend2(7) = 1
!--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
                           tile1(8) = 3; tile2(8) = 4
                           istart1(8) = nx; iend1(8) = nx; jstart1(8) = 1;  jend1(8) = ny
                           istart2(8) = 1;  iend2(8) = 1;  jstart2(8) = 1;  jend2(8) = ny
!--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
                           tile1(9) = 3; tile2(9) = 5
                           istart1(9) = 1;  iend1(9) = nx; jstart1(9) = ny; jend1(9) = ny
                           istart2(9) = 1;  iend2(9) = 1;  jstart2(9) = ny; jend2(9) = 1
!--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
                           tile1(10) = 4; tile2(10) = 5
                           istart1(10) = 1;  iend1(10) = nx; jstart1(10) = ny; jend1(10) = ny
                           istart2(10) = 1;  iend2(10) = nx; jstart2(10) = 1;  jend2(10) = 1
!--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
                           tile1(11) = 4; tile2(11) = 6
                           istart1(11) = nx; iend1(11) = nx; jstart1(11) = 1;  jend1(11) = ny
                           istart2(11) = nx; iend2(11) = 1;  jstart2(11) = 1;  jend2(11) = 1
!--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
                           tile1(12) = 5; tile2(12) = 6
                           istart1(12) = nx; iend1(12) = nx; jstart1(12) = 1;  jend1(12) = ny
                           istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;  jend2(12) = ny
                        end select

                        call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              pe_start=pe_start, pe_end=pe_end, symmetry=is_symmetry,              &
                              shalo = ng, nhalo = ng, whalo = ng, ehalo = ng, name = type)
                        call print_memuse_stats('external_ic:mpp_domain_decomp: mpp_define_mosaic')

                        deallocate(pe_start,pe_end)

!--- find the tile number
                        tile = mpp_pe()/npes_per_tile+1
                        call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
                        call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

                        call print_memuse_stats('external_ic:mpp_domain_decomp: mpp_get domains')

                     end subroutine mpp_domain_decomp
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

                     subroutine parallel_read_file_r8(fname, npts, is,ie, js,je, km, offset, var)
                        character(len=*), intent(IN) :: fname
                        integer,            intent(IN) :: npts, is,ie, js,je, km
                        integer (kind=MPI_OFFSET_KIND), intent(INOUT) :: offset
                        real(FVPRC),               intent(INOUT) :: var(is:ie, js:je, km)

                        integer :: ntiles=6
                        real(REAL64) :: var_r8(is:ie, js:je)
                        integer :: k

                        integer :: MUNIT=17
                        integer :: lsize, gsizes(2), distribs(2), dargs(2), psizes(2)
                        integer :: filetype
                        integer :: mcol, mrow, irow, jcol, mpiio_rank
                        integer :: rank, total_pes
                        integer :: mpistatus(MPI_STATUS_SIZE)
                        integer (kind=MPI_OFFSET_KIND) :: slice_2d

                        real(FVPRC) :: xmod, ymod
                        character(128) :: strErr

                        xmod = mod(npts,npes_x)
                        write(strErr, "(i4.4,' not evenly divisible by ',i4.4)") npts, npes_x
                        if (xmod /= 0) call mpp_error(FATAL, strErr)
                        ymod = mod(npts*6,npes_y)
                        write(strErr, "(i4.4,' not evenly divisible by ',i4.4)") npts*6, npes_y
                        if (ymod /= 0) call mpp_error(FATAL, strErr)

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, MUNIT, STATUS)
                        gsizes(1) = npts
                        gsizes(2) = npts * 6
                        distribs(1) = MPI_DISTRIBUTE_BLOCK
                        distribs(2) = MPI_DISTRIBUTE_BLOCK
                        dargs(1) = MPI_DISTRIBUTE_DFLT_DARG
                        dargs(2) = MPI_DISTRIBUTE_DFLT_DARG
                        psizes(1) = npes_x
                        psizes(2) = npes_y * 6
                        call MPI_COMM_SIZE(MPI_COMM_WORLD, total_pes, STATUS)
                        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, STATUS)
                        mcol = npes_x
                        mrow = npes_y*ntiles
                        irow = rank/mcol       !! logical row number
                        jcol = mod(rank, mcol) !! logical column number
                        mpiio_rank = jcol*mrow + irow
                        call MPI_TYPE_CREATE_DARRAY(total_pes, mpiio_rank, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, STATUS)
                        call MPI_TYPE_COMMIT(filetype, STATUS)
                        lsize = (ie-is+1)*(je-js+1)
                        slice_2d = npts*npts*ntiles
                        do k=1,km
                           call MPI_FILE_SET_VIEW(MUNIT, offset, MPI_DOUBLE_PRECISION, filetype, "native", MPI_INFO_NULL, STATUS)
                           call MPI_FILE_READ_ALL(MUNIT, var_r8, lsize, MPI_DOUBLE_PRECISION, mpistatus, STATUS)
                           var(:,:,k) = var_r8
                           offset = offset + slice_2d*8 + 8
                        enddo
                        call MPI_FILE_CLOSE(MUNIT, STATUS) 

                     end subroutine parallel_read_file_r8

                     subroutine parallel_read_file_r4(fname, npts, is,ie, js,je, km, offset, var)
                        character(len=*), intent(IN) :: fname
                        integer,            intent(IN) :: npts, is,ie, js,je, km
                        integer (kind=MPI_OFFSET_KIND), intent(INOUT) :: offset
                        real(FVPRC),               intent(INOUT) :: var(is:ie, js:je, km)

                        integer :: ntiles=6
                        real(REAL4) :: var_r4(is:ie, js:je)
                        integer :: k

                        integer :: MUNIT=17
                        integer :: lsize, gsizes(2), distribs(2), dargs(2), psizes(2)
                        integer :: filetype
                        integer :: mcol, mrow, irow, jcol, mpiio_rank
                        integer :: rank, total_pes
                        integer :: mpistatus(MPI_STATUS_SIZE)
                        integer (kind=MPI_OFFSET_KIND) :: slice_2d

                        real(FVPRC) :: xmod, ymod
                        character(128) :: strErr

                        xmod = mod(npts,npes_x)
                        write(strErr, "(i4.4,' not evenly divisible by ',i4.4)") npts, npes_x
                        if (xmod /= 0) call mpp_error(FATAL, strErr)
                        ymod = mod(npts*6,npes_y)
                        write(strErr, "(i4.4,' not evenly divisible by ',i4.4)") npts*6, npes_y
                        if (ymod /= 0) call mpp_error(FATAL, strErr)

                        call MPI_FILE_OPEN(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, MUNIT, STATUS)
                        gsizes(1) = npts
                        gsizes(2) = npts * 6
                        distribs(1) = MPI_DISTRIBUTE_BLOCK
                        distribs(2) = MPI_DISTRIBUTE_BLOCK
                        dargs(1) = MPI_DISTRIBUTE_DFLT_DARG
                        dargs(2) = MPI_DISTRIBUTE_DFLT_DARG
                        psizes(1) = npes_x
                        psizes(2) = npes_y * 6
                        call MPI_COMM_SIZE(MPI_COMM_WORLD, total_pes, STATUS)
                        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, STATUS)
                        mcol = npes_x
                        mrow = npes_y*ntiles
                        irow = rank/mcol       !! logical row number
                        jcol = mod(rank, mcol) !! logical column number
                        mpiio_rank = jcol*mrow + irow
                        call MPI_TYPE_CREATE_DARRAY(total_pes, mpiio_rank, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_FORTRAN, MPI_REAL, filetype, STATUS)
                        call MPI_TYPE_COMMIT(filetype, STATUS)
                        lsize = (ie-is+1)*(je-js+1)
                        slice_2d = npts*npts*ntiles
                        do k=1,km
                           call MPI_FILE_SET_VIEW(MUNIT, offset, MPI_REAL, filetype, "native", MPI_INFO_NULL, STATUS)
                           call MPI_FILE_READ_ALL(MUNIT, var_r4, lsize, MPI_REAL, mpistatus, STATUS)
                           var(:,:,k) = var_r4
                           offset = offset + slice_2d*4 + 8
                        enddo
                        call MPI_FILE_CLOSE(MUNIT, STATUS)

                     end subroutine parallel_read_file_r4

                  end module external_ic_mod

