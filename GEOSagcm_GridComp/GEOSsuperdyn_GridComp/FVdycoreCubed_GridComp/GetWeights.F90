   
!  $Id$
   
!!!#define REAL8 8
   
   subroutine GetWeights_init (in_ntiles,in_ncnst,in_npx,in_npy,in_npz,&
         in_nx,in_ny,in_hydro,in_mknh,comm)
      use fms_mod,           only: fms_init, set_domain
      use fv_control_mod,    only: fv_init1, fv_init2
      use fv_arrays_mod,     only: REAL4, REAL8, FVPRC
      use FV_StateMod,       only : FV_Atm  
      implicit none
      integer,intent(in) :: in_ntiles,in_ncnst
      integer,intent(in) :: in_npx,in_npy,in_npz
      integer,intent(in) :: in_nx,in_ny
      logical,intent(in) :: in_hydro,in_mknh
      integer            :: comm
!#ifdef SINGLE_FV
      real(FVPRC),    parameter :: dt_who_cares = 1800.
!#else
!      real*8,  parameter :: dt_who_cares = 1800.d0
!#endif
      integer :: p_split
      logical, allocatable :: grids_on_my_pe(:)
   
      p_split = 1

      call fms_init(comm)

      call fv_init1(FV_atm, dt_who_cares, grids_on_my_pe, p_split)

      FV_Atm(1)%flagstruct%ntiles = in_ntiles
      FV_Atm(1)%flagstruct%hydrostatic = in_hydro
      FV_Atm(1)%flagstruct%Make_NH     = in_mknh
      FV_Atm(1)%flagstruct%ncnst  = in_ncnst ! nq
      FV_Atm(1)%flagstruct%npx    = in_npx
      FV_Atm(1)%flagstruct%npy    = in_npy
      FV_Atm(1)%flagstruct%npz    = in_npz

      FV_Atm(1)%flagstruct%npx=FV_Atm(1)%flagstruct%npx+1
      FV_Atm(1)%flagstruct%npy=FV_Atm(1)%flagstruct%npy+1

      call fv_init2(FV_atm, dt_who_cares, grids_on_my_pe, p_split)

      call set_domain(FV_Atm(1)%domain)

   end subroutine GetWeights_init

   subroutine GetWeights(npx, npy, nlat, nlon, index, weight, id1, id2, jdc, l2c, &
         ee1, ee2, ff1, ff2, gg1, gg2, e1, e2, f1, f2, g1, g2, sublons, sublats, AmNodeRoot, WriteNetcdf) 
#include "MAPL_Generic.h"

      use MAPL_BaseMod
      use MAPL_GenericMod
      use MAPL_ShmemMod
      use MAPL_ErrorHandlingMod
      use fv_grid_utils_mod, only : gnomonic_grids, cell_center2, mid_pt_sphere
      use fv_grid_tools_mod, only : mirror_grid
      use fv_grid_tools_mod, only : get_unit_vector
      use fv_grid_utils_mod, only : inner_prod
      use fv_arrays_mod,     only : REAL4, REAL8, FVPRC
      use CUB2LATLON_mod,    only : init_latlon_grid, get_c2l_weight
      use GHOST_CUBSPH_mod,  only : B_grid, A_grid, ghost_cubsph_update
      use FV_StateMod,       only : FV_Atm  
      use fv_mp_mod,         only : is,js,ie,je, is_master

      include "netcdf.inc"

      integer,                     intent(in   ) :: npx,  npy
      integer,                     intent(in   ) :: nlon, nlat
      integer,                     intent(  out) :: index(3,nlon,nlat)
      real(REAL8),                 intent(  out) :: weight(4,nlon,nlat)
      integer,                     intent(  out) :: id1(npx,npy)
      integer,                     intent(  out) :: id2(npx,npy)
      integer,                     intent(  out) :: jdc(npx,npy)
      real(REAL8),                 intent(  out) :: l2c(4,npx,npy)
      real(REAL8),                 intent(  out) :: ee1(npx,npy,3)
      real(REAL8),                 intent(  out) :: ee2(npx,npy,3) 
      real(REAL8),                 intent(  out) :: ff1(npx,npy,3) 
      real(REAL8),                 intent(  out) :: ff2(npx,npy,3) 
      real(REAL8),                 intent(  out) :: gg1(npx,npy,3) 
      real(REAL8),                 intent(  out) :: gg2(npx,npy,3)
      real(REAL8), pointer                       ::  e1(:,:,:)
      real(REAL8), pointer                       ::  e2(:,:,:)
      real(REAL8), pointer                       ::  f1(:,:,:)
      real(REAL8), pointer                       ::  f2(:,:,:)
      real(REAL8), pointer                       ::  g1(:,:,:)
      real(REAL8), pointer                       ::  g2(:,:,:)
      real(REAL8), optional                      :: sublons(:)
      real(REAL8), optional                      :: sublats(:)
      logical,     optional                      :: AmNodeRoot
      logical,     optional                      :: WriteNetcdf
! Locals
!-------

! NETCDF stuff for weights file.
!-----------------------------------
      integer                               :: c2l_unit
      character(len=128)                    :: c2l_fname
      logical                               :: c2l_file_exists
      integer :: LONDIM,LATDIM,NPXDIM,NPYDIM,IN3DIM,IN4DIM
      integer :: INDEX_ID,WEIGTH_ID,L2C_ID
      integer :: ID1_ID,ID2_ID,JDC_ID
      integer :: EE1_ID,EE2_ID,FF1_ID,FF2_ID,GG1_ID,GG2_ID
      integer :: STATUS

      integer :: npts, n, l, j, j1

      integer, parameter :: ntiles=6
      integer, parameter :: ndims=2

      character(len=128), parameter :: Iam="GetWeights"
      integer :: rc

      real(REAL8), parameter :: PI=3.14159265358979323846

! Real*8 are needed to make fv calls.
!-----------------------------------

      real(REAL8), allocatable :: grid_global(:,:,:,:)
      real(REAL8), allocatable :: sph_corner (:,:,:,:)
      real(REAL8), allocatable :: xlon(:), ylat(:)
      real(REAL8), allocatable :: agrid(:,:,:)
      real(REAL8), allocatable :: slon(:), slat(:), clon(:), clat(:)



      if (AmNodeRoot) then
         write(c2l_fname,'("PE",i0,"x",i0,"-CF_c2l_PC",i0,"x",i0,"-DC.nc4")') npx,npy,nlon,nlat
         inquire(FILE=TRIM(c2l_fname), EXIST=c2l_file_exists)
         if (.not. c2l_file_exists) then

            if (is_master()) print *, 'Computing weights for ', TRIM(c2l_fname)

            npts = npx + 1

            allocate( grid_global(npts,npts,ndims,ntiles))

            call gnomonic_grids(A_grid, npx, grid_global(:,:,1,1), grid_global(:,:,2,1))

! mirror_grid assumes that the tile=1 is centered 
!   on equator and greenwich meridian Lon[-pi,pi]
!------------------------------------------------

            call mirror_grid(grid_global, 0, npts, npts, ndims, ntiles)

! Shift the corner away from Japan.
!  This will result in the corner 
!  close to the east coast of China.
!-----------------------------------

            grid_global(:,:,1,:) = grid_global(:,:,1,:) - PI/18.

            where(grid_global(:,:,1,:) < 0.) &
                  grid_global(:,:,1,:) = grid_global(:,:,1,:) + 2.* PI

! Keep Equator and Greenwich exact
!---------------------------------

            where(abs(grid_global(:,:,:,1)) < 1.e-10) grid_global(:,:,:,1) = 0.0


! Clean Up Corners
!---------------------------------

            grid_global(1   ,   :,:,2)=grid_global(npts     ,:        ,:,1)
            grid_global(1   ,   :,:,3)=grid_global(npts:1:-1,npts     ,:,1)
            grid_global(:   ,npts,:,5)=grid_global(1        ,npts:1:-1,:,1)
            grid_global(:   ,npts,:,6)=grid_global(:        ,1        ,:,1)
            grid_global(:   ,   1,:,3)=grid_global(:        ,npts     ,:,2)
            grid_global(:   ,   1,:,4)=grid_global(npts     ,npts:1:-1,:,2)
            grid_global(npts,   :,:,6)=grid_global(npts:1:-1,1        ,:,2)
            grid_global(1   ,   :,:,4)=grid_global(npts     ,:        ,:,3)
            grid_global(1   ,   :,:,5)=grid_global(npts:1:-1,npts     ,:,3)
            grid_global(npts,   :,:,3)=grid_global(1        ,:        ,:,4)
            grid_global(:   ,   1,:,5)=grid_global(:        ,npts     ,:,4)
            grid_global(:   ,   1,:,6)=grid_global(npts     ,npts:1:-1,:,4)
            grid_global(1   ,   :,:,6)=grid_global(npts     ,:        ,:,5)

! This is for C2L
!----------------

            allocate( sph_corner(ndims,0:npts+1,0:npts+1,ntiles))

            sph_corner(1,1:npts,1:npts,:) = grid_global(:,:,1,:)
            sph_corner(2,1:npts,1:npts,:) = grid_global(:,:,2,:)

! This is for L2C
!----------------

            allocate( agrid(npx,npy,ndims))
            agrid = 1.e+25
            do n=1,ntiles
               do j=1,npx
                  j1 = npx*(n-1) + j
                  do i=1,npx
                     call cell_center2(grid_global(i,j,  1:2,n), grid_global(i+1,j,  1:2,n),   &
                           grid_global(i,j+1,1:2,n), grid_global(i+1,j+1,1:2,n),   &
                           agrid(i,j1,1:2) )
                  enddo
               enddo
            enddo

            do n=1,ntiles
               do j=1,npx
                  j1 = npx*(n-1) + j
                  do i=1,npx
                     call CreateCube2LatLonRotation( &
                           grid_global(i:i+1,j:j+1,:,n), agrid(i,j1,1:2), &
                           ee1(i,j1,:),ee2(i,j1,:),ff1(i,j1,:),ff2(i,j1,:),gg1(i,j1,:),gg2(i,j1,:))
                  enddo
               enddo
            enddo

            deallocate ( grid_global )

! do halo update                                                   !
!------------------------------------------------------------------!

            do n=1,ntiles
               sph_corner(1:2,0     ,0     ,n)=0.
               sph_corner(1:2,npts+1,0     ,n)=0.
               sph_corner(1:2,0     ,npts+1,n)=0.
               sph_corner(1:2,npts+1,npts+1,n)=0.
               do L=1,2
                  call ghost_cubsph_update(var=sph_corner(L,0:npts+1,0:npts+1,:), &
                        n1x=0, nx=npts+1, n1y=0, ny=npts+1,           &
                        ng=1, ntiles=ntiles, ng_update=1, tile=n, grid_type=B_grid         )
               end do
            enddo


! initialize latlon grid
!-----------------------
            allocate(xlon(nlon), ylat(nlat))

            if (present(sublons) .and. present(sublats)) then
               do n=1,nlon
                  xlon(n) = sublons(n) !+ 180.0_8
               enddo
               do n=1,nlat
                  ylat(n) = sublats(n)
               enddo
            else
               call init_latlon_grid(xlon, ylat, nlon, nlat)
            endif

! calculate weights for bilinear interpolation
! from cubed sphere to latlon grid            
!---------------------------------------------
            if (present(sublons) .and. present(sublats)) then
               call get_c2l_weight(sph_corner, npts, npts, ntiles, &
                     xlon, ylat, nlon, nlat,  &
                     index, weight, .true.)
            else
               call get_c2l_weight(sph_corner, npts, npts, ntiles, &
                     xlon, ylat, nlon, nlat,  &
                     index, weight, .false.)
            endif

            deallocate ( sph_corner )

! calculate weights for bilinear interpolation
! from cubed sphere to latlon grid            
!---------------------------------------------

            call remap_coef( agrid, xlon, ylat, id1, id2, jdc, l2c )

            deallocate ( xlon, ylat )
            deallocate ( agrid      )

! write out NETCDF weights file            
!---------------------------------------------
            if (present(WriteNetcdf)) then
               if (WriteNetcdf) then

                  if (is_master()) print *, 'Writing weights to ', TRIM(c2l_fname)
                  STATUS = NF_CREATE (trim(c2l_fname), IOR(NF_CLOBBER,NF_NETCDF4), c2l_unit)

                  STATUS = NF_DEF_DIM(c2l_unit, 'lat', nlat, LATDIM)
                  STATUS = NF_DEF_DIM(c2l_unit, 'lon', nlon, LONDIM)
                  STATUS = NF_DEF_DIM(c2l_unit, 'npx', npx , NPXDIM)
                  STATUS = NF_DEF_DIM(c2l_unit, 'npy', npy , NPYDIM)
                  STATUS = NF_DEF_DIM(c2l_unit, 'in3',    3, IN3DIM)
                  STATUS = NF_DEF_DIM(c2l_unit, 'in4',    4, IN4DIM)

                  STATUS = NF_DEF_VAR (c2l_unit,  'index', NF_INT    , 3, (/IN3DIM,LONDIM,LATDIM/),  INDEX_ID)
                  STATUS = NF_DEF_VAR (c2l_unit, 'weight', NF_DOUBLE , 3, (/IN4DIM,LONDIM,LATDIM/), WEIGTH_ID)
                  STATUS = NF_DEF_VAR (c2l_unit,    'l2c', NF_DOUBLE , 3, (/IN4DIM,NPXDIM,NPYDIM/),    L2C_ID)

                  STATUS = NF_DEF_VAR (c2l_unit,    'id1', NF_INT    , 2, (/NPXDIM,NPYDIM/)       ,    ID1_ID)
                  STATUS = NF_DEF_VAR (c2l_unit,    'id2', NF_INT    , 2, (/NPXDIM,NPYDIM/)       ,    ID2_ID)
                  STATUS = NF_DEF_VAR (c2l_unit,    'jdc', NF_INT    , 2, (/NPXDIM,NPYDIM/)       ,    JDC_ID)

                  STATUS = NF_DEF_VAR (c2l_unit,    'ee1', NF_DOUBLE , 3, (/NPXDIM,NPYDIM,IN3DIM/),    EE1_ID)
                  STATUS = NF_DEF_VAR (c2l_unit,    'ee2', NF_DOUBLE , 3, (/NPXDIM,NPYDIM,IN3DIM/),    EE2_ID)
                  STATUS = NF_DEF_VAR (c2l_unit,    'ff1', NF_DOUBLE , 3, (/NPXDIM,NPYDIM,IN3DIM/),    FF1_ID)
                  STATUS = NF_DEF_VAR (c2l_unit,    'ff2', NF_DOUBLE , 3, (/NPXDIM,NPYDIM,IN3DIM/),    FF2_ID)
                  STATUS = NF_DEF_VAR (c2l_unit,    'gg1', NF_DOUBLE , 3, (/NPXDIM,NPYDIM,IN3DIM/),    GG1_ID)
                  STATUS = NF_DEF_VAR (c2l_unit,    'gg2', NF_DOUBLE , 3, (/NPXDIM,NPYDIM,IN3DIM/),    GG2_ID)

                  STATUS = NF_ENDDEF(c2l_unit)

                  STATUS = NF_PUT_VARA_INT     (c2l_unit,  INDEX_ID, (/1,1,1/), (/3,nlon,nlat/),  index)
                  STATUS = NF_PUT_VARA_DOUBLE  (c2l_unit, WEIGTH_ID, (/1,1,1/), (/4,nlon,nlat/), weight)
                  STATUS = NF_PUT_VARA_DOUBLE  (c2l_unit,    L2C_ID, (/1,1,1/), (/4,npx ,npy /),    l2c)

                  STATUS = NF_PUT_VARA_INT     (c2l_unit,    ID1_ID, (/1,1/)  , (/npx,npy/)    ,    id1)
                  STATUS = NF_PUT_VARA_INT     (c2l_unit,    ID2_ID, (/1,1/)  , (/npx,npy/)    ,    id2)
                  STATUS = NF_PUT_VARA_INT     (c2l_unit,    JDC_ID, (/1,1/)  , (/npx,npy/)    ,    jdc)

                  STATUS = NF_PUT_VARA_DOUBLE  (c2l_unit,    EE1_ID, (/1,1,1/), (/npx,npy,3/)  ,    ee1)
                  STATUS = NF_PUT_VARA_DOUBLE  (c2l_unit,    EE2_ID, (/1,1,1/), (/npx,npy,3/)  ,    ee2)
                  STATUS = NF_PUT_VARA_DOUBLE  (c2l_unit,    FF1_ID, (/1,1,1/), (/npx,npy,3/)  ,    ff1)
                  STATUS = NF_PUT_VARA_DOUBLE  (c2l_unit,    FF2_ID, (/1,1,1/), (/npx,npy,3/)  ,    ff2)
                  STATUS = NF_PUT_VARA_DOUBLE  (c2l_unit,    GG1_ID, (/1,1,1/), (/npx,npy,3/)  ,    gg1)
                  STATUS = NF_PUT_VARA_DOUBLE  (c2l_unit,    GG2_ID, (/1,1,1/), (/npx,npy,3/)  ,    gg2)

                  STATUS = NF_CLOSE (c2l_unit)
               endif
            endif

         else  ! NOT WriteNetcdf, so read in the weights

            if (is_master()) print *, 'Reading weights for ', TRIM(c2l_fname)

! read NETCDF weights file            
!---------------------------------------------
            STATUS = NF_OPEN (trim(c2l_fname), NF_NOWRITE, c2l_unit)

            STATUS = NF_INQ_VARID (c2l_unit,  'index',  INDEX_ID)
            STATUS = NF_INQ_VARID (c2l_unit, 'weight', WEIGTH_ID)
            STATUS = NF_INQ_VARID (c2l_unit,    'l2c',    L2C_ID)

            STATUS = NF_INQ_VARID (c2l_unit,    'id1',    ID1_ID)
            STATUS = NF_INQ_VARID (c2l_unit,    'id2',    ID2_ID)
            STATUS = NF_INQ_VARID (c2l_unit,    'jdc',    JDC_ID)

            STATUS = NF_INQ_VARID (c2l_unit,    'ee1',    EE1_ID)
            STATUS = NF_INQ_VARID (c2l_unit,    'ee2',    EE2_ID)
            STATUS = NF_INQ_VARID (c2l_unit,    'ff1',    FF1_ID)
            STATUS = NF_INQ_VARID (c2l_unit,    'ff2',    FF2_ID)
            STATUS = NF_INQ_VARID (c2l_unit,    'gg1',    GG1_ID)
            STATUS = NF_INQ_VARID (c2l_unit,    'gg2',    GG2_ID)

            STATUS = NF_GET_VARA_INT     (c2l_unit,  INDEX_ID, (/1,1,1/), (/3,nlon,nlat/),  index)
            STATUS = NF_GET_VARA_DOUBLE  (c2l_unit, WEIGTH_ID, (/1,1,1/), (/4,nlon,nlat/), weight)
            STATUS = NF_GET_VARA_DOUBLE  (c2l_unit,    L2C_ID, (/1,1,1/), (/4,npx ,npy /),    l2c)

            STATUS = NF_GET_VARA_INT     (c2l_unit,    ID1_ID, (/1,1/)  , (/npx,npy/)    ,    id1)
            STATUS = NF_GET_VARA_INT     (c2l_unit,    ID2_ID, (/1,1/)  , (/npx,npy/)    ,    id2)
            STATUS = NF_GET_VARA_INT     (c2l_unit,    JDC_ID, (/1,1/)  , (/npx,npy/)    ,    jdc)

            STATUS = NF_GET_VARA_DOUBLE  (c2l_unit,    EE1_ID, (/1,1,1/), (/npx,npy,3/)  ,    ee1)
            STATUS = NF_GET_VARA_DOUBLE  (c2l_unit,    EE2_ID, (/1,1,1/), (/npx,npy,3/)  ,    ee2)
            STATUS = NF_GET_VARA_DOUBLE  (c2l_unit,    FF1_ID, (/1,1,1/), (/npx,npy,3/)  ,    ff1)
            STATUS = NF_GET_VARA_DOUBLE  (c2l_unit,    FF2_ID, (/1,1,1/), (/npx,npy,3/)  ,    ff2)
            STATUS = NF_GET_VARA_DOUBLE  (c2l_unit,    GG1_ID, (/1,1,1/), (/npx,npy,3/)  ,    gg1)
            STATUS = NF_GET_VARA_DOUBLE  (c2l_unit,    GG2_ID, (/1,1,1/), (/npx,npy,3/)  ,    gg2)

            STATUS = NF_CLOSE (c2l_unit)

         endif ! WriteNetcdf

      endif ! AmNodeRoot

! Everyone Needs their copy of the Vector Rotation arrays
!--------------------------------------------------------
      if(.not. allocated(FV_Atm)) return

      call MAPL_SyncSharedMemory(rc=STATUS)
      _VERIFY(STATUS)
      allocate(e1(is:ie,js:je,3))
      allocate(e2(is:ie,js:je,3))
      allocate(f1(is:ie,js:je,3))
      allocate(f2(is:ie,js:je,3))
      allocate(g1(is:ie,js:je,3))
      allocate(g2(is:ie,js:je,3))
      n = FV_Atm(1)%tile
      do j=1,npx
         j1 = npx*(n-1) + j
         do i=1,npx
            if(i>=is.and.j>=js.and.i<=ie.and.j<=je) then
               e1(i,j,:) = ee1(i,j1,:)
               e2(i,j,:) = ee2(i,j1,:)
               f1(i,j,:) = ff1(i,j1,:)
               f2(i,j,:) = ff2(i,j1,:)
               g1(i,j,:) = gg1(i,j1,:)
               g2(i,j,:) = gg2(i,j1,:)
            end if
         enddo
      enddo

      return

   contains

      subroutine remap_coef( agrid, lon, lat, id1, id2, jdc, l2c )

         real(REAL8), intent(in)  :: agrid(:,:,:)
         real(REAL8), intent(in)  :: lon(:), lat(:)
         real(REAL8), intent(out) :: l2c(:,:,:)
         integer,  intent(out) :: id1(:,:), id2(:,:), jdc(:,:)

! local:

         real(REAL8) :: a1, b1
         integer  :: i,j, im, jm, i1, i2, jc, i0, j0
         real(REAL8) :: slt(size(lat))

! Interpolate to cubed sphere cell center

         im = size(lon)
         jm = size(lat)

         slt = sin(lat)

         do j=1,size(agrid,2)
            do i=1,size(agrid,1)

               do i1= 1, im
                  i2 = mod(i1,im) + 1
                  d1 = modulo(agrid(i,j,1) - lon(i1), 2*pi)
                  d2 = modulo(lon(i2) - lon(i1), 2*pi)

                  if (d2 >= d1) then
                     a1 = d1 / d2
                     exit
                  end if
               enddo

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
                        b1 = (sin(agrid(i,j,2))-slt(jc)) / (slt(jc+1) - slt(jc))
                        exit
                     endif
                  enddo
               endif

               l2c(1,i,j) = (1.-a1) * (1.-b1)
               l2c(2,i,j) =     a1  * (1.-b1)
               l2c(3,i,j) =     a1  *     b1
               l2c(4,i,j) = (1.-a1) *     b1

               id1(i,j)   = i1
               id2(i,j)   = i2
               jdc(i,j)   = jc

            enddo
         enddo

      end subroutine remap_coef


      subroutine CreateCube2LatLonRotation(grid, center, ee1, ee2, ff1, ff2, gg1, gg2)

         real(REAL8), intent(IN)    :: grid(0:1,0:1,2), center(2)
         real(REAL8), intent(OUT), dimension(3) :: ee1, ee2, ff1, ff2, gg1, gg2

         real(REAL8), dimension(2) :: p1, p2, p3, p4
         real :: H, F


         call mid_pt_sphere(grid(0,0,:),grid(0,1,:), p1)
         call mid_pt_sphere(grid(0,0,:),grid(1,0,:), p2)
         call mid_pt_sphere(grid(1,0,:),grid(1,1,:), p3)
         call mid_pt_sphere(grid(0,1,:),grid(1,1,:), p4)

         call get_unit_vector(p3, center, p1, ee1)
         call get_unit_vector(p4, center, p2, ee2)

         H   = dot_product(ee1,ee2)
         F   = 1.0/(H**2-1.0)

         ff1 = F*(ee2*H-ee1)
         ff2 = F*(ee1*H-ee2)

         gg1(1) = -SIN(center(1) - PI)
         gg1(2) =  COS(center(1) - PI)
         gg1(3) = 0.0
         gg2(1) = -SIN(center(2))*gg1(2)
         gg2(2) =  SIN(center(2))*gg1(1)
         gg2(3) =  COS(center(2))

         return
      end subroutine CreateCube2LatLonRotation

   end subroutine GetWeights

!  version of a2d3d routine from FV_StateMod

   subroutine A2D2C(U,V,npz,getC)

! Move A-Grid winds/tendencies oriented on lat/lon to the D-grid, or Optionally C-grid cubed-sphere orientation
! will return d-grid winds unless user asks for c-grid winds

      use mpp_domains_mod,   only: mpp_update_domains, mpp_get_boundary, DGRID_NE
      use mpp_parameter_mod, only: AGRID
      use FV_StateMod,       only: fv_atm
      use fv_arrays_mod,     only: REAL4, REAL8, FVPRC
      use fv_mp_mod,         only: is,js,ie,je,isd,jsd,ied,jed,ng
      use sw_core_mod,       only: d2a2c_vect
      implicit none

      real,    intent(INOUT)           :: U(:,:,:)
      real,    intent(INOUT)           :: V(:,:,:)
      integer, intent(   IN)           :: npz
      logical, intent(   IN)           :: getC

!local variables
      integer :: i,j,k, im2,jm2

      real(REAL8) :: ud(is:ie,js:je+1,npz)
      real(REAL8) :: vd(is:ie+1,js:je,npz)

      real(FVPRC) :: ut(isd:ied, jsd:jed)
      real(FVPRC) :: vt(isd:ied, jsd:jed)

      real(REAL8) :: v3(is-1:ie+1,js-1:je+1,3)
      real(REAL8) :: ue(is-1:ie+1,js  :je+1,3)    ! 3D winds at edges
      real(REAL8) :: ve(is  :ie+1,js-1:je+1,3)    ! 3D winds at edges
      real(REAL8), dimension(is:ie):: ut1, ut2, ut3
      real(REAL8), dimension(js:je):: vt1, vt2, vt3

      real(FVPRC) :: uctemp(isd:ied+1,jsd:jed  ,npz)
      real(FVPRC) :: vctemp(isd:ied  ,jsd:jed+1,npz)
      real(FVPRC) ::  utemp(isd:ied  ,jsd:jed+1,npz)
      real(FVPRC) ::  vtemp(isd:ied+1,jsd:jed  ,npz)
      real(FVPRC) :: uatemp(isd:ied,jsd:jed,npz)
      real(FVPRC) :: vatemp(isd:ied,jsd:jed,npz)
!#else
!      real*8 :: uctemp(isd:ied+1,jsd:jed  ,npz)
!!      real*8 :: vctemp(isd:ied  ,jsd:jed+1,npz)
!#endif

      real(FVPRC) :: wbuffer(js:je,npz)
      real(FVPRC) :: sbuffer(is:ie,npz)
      real(FVPRC) :: ebuffer(js:je,npz)
      real(FVPRC) :: nbuffer(is:ie,npz)
      integer     :: npx, npy
      integer :: STATUS

      npx = FV_Atm(1)%npx
      npy = FV_Atm(1)%npy

      uatemp = 0.0d0
      vatemp = 0.0d0

      uatemp(is:ie,js:je,:) = U
      vatemp(is:ie,js:je,:) = V

      im2 = (npx-1)/2
      jm2 = (npy-1)/2

! Cubed-Sphere
      call mpp_update_domains(uatemp, FV_Atm(1)%domain, complete=.false.)
      call mpp_update_domains(vatemp, FV_Atm(1)%domain, complete=.true.)
      do k=1, npz
! Compute 3D wind tendency on A grid
         do j=js-1,je+1
            do i=is-1,ie+1
               v3(i,j,1) = uatemp(i,j,k)*fv_atm(1)%gridstruct%vlon(i,j,1) + vatemp(i,j,k)*fv_atm(1)%gridstruct%vlat(i,j,1)
               v3(i,j,2) = uatemp(i,j,k)*fv_atm(1)%gridstruct%vlon(i,j,2) + vatemp(i,j,k)*fv_atm(1)%gridstruct%vlat(i,j,2)
               v3(i,j,3) = uatemp(i,j,k)*fv_atm(1)%gridstruct%vlon(i,j,3) + vatemp(i,j,k)*fv_atm(1)%gridstruct%vlat(i,j,3)
            enddo
         enddo
! A --> D
! Interpolate to cell edges
         do j=js,je+1
            do i=is-1,ie+1
               ue(i,j,1) = v3(i,j-1,1) + v3(i,j,1)
               ue(i,j,2) = v3(i,j-1,2) + v3(i,j,2)
               ue(i,j,3) = v3(i,j-1,3) + v3(i,j,3)
            enddo
         enddo

         do j=js-1,je+1
            do i=is,ie+1
               ve(i,j,1) = v3(i-1,j,1) + v3(i,j,1)
               ve(i,j,2) = v3(i-1,j,2) + v3(i,j,2)
               ve(i,j,3) = v3(i-1,j,3) + v3(i,j,3)
            enddo
         enddo
! --- E_W edges (for v-wind):
         if ( is==1 ) then
            i = 1
            do j=js,je
               if ( j>jm2 ) then
                  vt1(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,1)
                  vt2(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,2)
                  vt3(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,3)
               else
                  vt1(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,1)
                  vt2(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,2)
                  vt3(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,3)
               endif
            enddo
            do j=js,je
               ve(i,j,1) = vt1(j)
               ve(i,j,2) = vt2(j)
               ve(i,j,3) = vt3(j)
            enddo
         endif
         if ( (ie+1)==npx ) then
            i = npx
            do j=js,je
               if ( j>jm2 ) then
                  vt1(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j-1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,1)
                  vt2(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j-1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,2)
                  vt3(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j-1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,3)
               else
                  vt1(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j+1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,1)
                  vt2(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j+1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,2)
                  vt3(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j+1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,3)
               endif
            enddo
            do j=js,je
               ve(i,j,1) = vt1(j)
               ve(i,j,2) = vt2(j)
               ve(i,j,3) = vt3(j)
            enddo
         endif
! N-S edges (for u-wind):
         if ( js==1 ) then
            j = 1
            do i=is,ie
               if ( i>im2 ) then
                  ut1(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i-1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,1)
                  ut2(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i-1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,2)
                  ut3(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i-1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,3)
               else
                  ut1(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i+1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,1)
                  ut2(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i+1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,2)
                  ut3(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i+1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,3)
               endif
            enddo
            do i=is,ie
               ue(i,j,1) = ut1(i)
               ue(i,j,2) = ut2(i)
               ue(i,j,3) = ut3(i)
            enddo
         endif

         if ( (je+1)==npy ) then
            j = npy
            do i=is,ie
               if ( i>im2 ) then
                  ut1(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,1)
                  ut2(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,2)
                  ut3(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,3)
               else
                  ut1(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,1)
                  ut2(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,2)
                  ut3(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,3)
               endif
            enddo
            do i=is,ie
               ue(i,j,1) = ut1(i)
               ue(i,j,2) = ut2(i)
               ue(i,j,3) = ut3(i)
            enddo
         endif
! Update:
         do j=js,je+1
            do i=is,ie
               ud(i,j,k) = 0.5*( ue(i,j,1)*fv_atm(1)%gridstruct%es(1,i,j,1) +  &
                     ue(i,j,2)*fv_atm(1)%gridstruct%es(2,i,j,1) +  &
                     ue(i,j,3)*fv_atm(1)%gridstruct%es(3,i,j,1) )
            enddo
         enddo
         do j=js,je
            do i=is,ie+1
               vd(i,j,k) = 0.5*( ve(i,j,1)*fv_atm(1)%gridstruct%ew(1,i,j,2) +  &
                     ve(i,j,2)*fv_atm(1)%gridstruct%ew(2,i,j,2) +  &
                     ve(i,j,3)*fv_atm(1)%gridstruct%ew(3,i,j,2) )
            enddo
         enddo

      enddo         ! k-loop

      if (getC) then

         ! Now we have D-Grid winds, need to make call to d2a2c_vect

         utemp = 0.0d0
         vtemp = 0.0d0
         utemp(is:ie,js:je,:) = ud(is:ie,js:je,:)
         vtemp(is:ie,js:je,:) = vd(is:ie,js:je,:)

         ! update shared edges
         call mpp_get_boundary(utemp, vtemp, FV_Atm(1)%domain, &
                               wbuffery=wbuffer, ebuffery=ebuffer, &
                               sbufferx=sbuffer, nbufferx=nbuffer, &
                               gridtype=DGRID_NE, complete=.true. )
         do k=1,npz
            do i=is,ie
               utemp(i,je+1,k) = nbuffer(i,k)
            enddo
            do j=js,je
               vtemp(ie+1,j,k) = ebuffer(j,k)
            enddo
         enddo

         call mpp_update_domains(utemp, vtemp, FV_Atm(1)%domain, gridtype=DGRID_NE, complete=.true.)
         do k=1,npz
            call d2a2c_vect(utemp(:,:,k),  vtemp(:,:,k), &
                      uatemp(:,:,k), vatemp(:,:,k), &
                      uctemp(:,:,k), vctemp(:,:,k), ut, vt, .true., &
                      fv_atm(1)%gridstruct,fv_atm(1)%bd,npx,npy,.false.,0)
         enddo

         U(:,:,:) = uctemp(is:ie,js:je,:)
         V(:,:,:) = vctemp(is:ie,js:je,:)

      else

         ! return d-grid winds
         u = ud(is:ie,js:je,:)
         v = vd(is:ie,js:je,:)

      end if

   end subroutine A2D2C
