!  $Id$

#define SUCCESS 0
#define VERIFY_(A) if((A)/=0) then; if(present(rc)) rc=A; PRINT *, Iam, __LINE__; return; endif
#define ASSERT_(A) if(.not.(A)) then; if(present(rc)) rc=1; PRINT *, Iam, __LINE__; return; endif
#define RETURN_(A) if(present(rc)) rc=A; return

#define DEALLOCGLOB_(A) call deallocGlob(A,status);VERIFY_(status)

#define DEALLOCLOCL_(A) if(associated(A)) then; deallocate(A, stat=STATUS); VERIFY_(STATUS); NULLIFY(A); endif

#ifdef TAU_PROFILE
#undef ASSERT_
#define ASSERT_(A)

#undef VERIFY_
#define VERIFY_(A)

#undef RETURN_
#define RETURN_(A)
#endif

Module CubeLatLonTransformMod
  
  use ESMF
  use MAPL_BaseMod
  use MAPL_LocStreamMod
  use MAPL_IOMod
  use MAPL_CommsMod
  use MAPL_ShmemMod
  use, intrinsic :: iso_fortran_env, only: REAL64

  implicit none
  private

  public T_CubeLatLonTransform
  public T_CubeCubeTransform
  public CubeLatLonIsCreated
  public CubeLatLonSubset
  public CubeLatLonCreate
  public CubeLatLonDestroy
  public CubeToLatLon
  public CubeCubeIsCreated
  public CubeCubeCreate
  public CubeCubeDestroy
  public CubeToCube
  public LatLonToCube
  public CartesianToSpherical
  public SphericalToCartesian
  public RunTileTransform
  public MAPL_RegridLSCreate
  public TileTransformTranspose
  public get_conservative_weights !JK patch for conservative interp
  public RestaggerWindsCube

  include 'mpif.h'

  character*30 :: Iam="CubeToLatLon in line "

 type T_CubeLatLonTransform
     private
     real(REAL64),pointer    :: weight(:,:,:) => null()
     real(REAL64),pointer    :: l2c(:,:,:)   => null()
     integer, pointer    :: index (:,:,:) => null()
     integer, pointer    :: id1(:,:) => null()
     integer, pointer    :: id2(:,:) => null()
     integer, pointer    :: jdc(:,:) => null()
     logical             :: Created=.false.
     character(len=120)  :: name
     integer             :: npx, npy, nlon, nlat
! global
     real(REAL64), pointer   :: ee1(:,:,:)  => null()
     real(REAL64), pointer   :: ee2(:,:,:)  => null()
     real(REAL64), pointer   :: ff1(:,:,:)  => null()
     real(REAL64), pointer   :: ff2(:,:,:)  => null()
     real(REAL64), pointer   :: gg1(:,:,:)  => null()
     real(REAL64), pointer   :: gg2(:,:,:)  => null()
! local
     real(REAL64), pointer   :: e1(:,:,:)  => null()
     real(REAL64), pointer   :: e2(:,:,:)  => null()
     real(REAL64), pointer   :: f1(:,:,:)  => null()
     real(REAL64), pointer   :: f2(:,:,:)  => null()
     real(REAL64), pointer   :: g1(:,:,:)  => null()
     real(REAL64), pointer   :: g2(:,:,:)  => null()
!
     real(REAL64), pointer   :: elon(:,:,:) => null()
     real(REAL64), pointer   :: elat(:,:,:) => null()
     real(REAL64), pointer   :: elon_local(:,:,:) => null()
     real(REAL64), pointer   :: elat_local(:,:,:) => null()
!
     logical             :: lsCreated = .false.
     type(MAPL_LocStream) :: locStIn
     type(MAPL_LocStream) :: locStOut
     type(MAPL_LocStreamXform) :: XformInOut
     type(MAPL_LocStreamXform) :: XformOutIn
     logical             :: subset = .false.
  end type T_CubeLatLonTransform

 type T_CubeCubeTransform
     private
     real(REAL64),pointer    :: weight(:,:,:,:) => NULL()
     integer, pointer    :: index (:,:,:,:) => NULL()
     real(REAL64), pointer   :: ee1(:,:,:) => NULL()
     real(REAL64), pointer   :: ee2(:,:,:) => NULL()
     real(REAL64), pointer   :: ff1(:,:,:) => NULL()
     real(REAL64), pointer   :: ff2(:,:,:) => NULL()
     logical             :: Created=.false.
     character(len=120)  :: name
     integer             :: npx, npy, npxout, npyout
  end type T_CubeCubeTransform

  interface CubeToLatLon
     module procedure CubeToLatLonr8
     module procedure CubeToLatLonr4
  end interface

  interface LatLonToCube
     module procedure LatLonToCuber8
     module procedure LatLonToCuber4
  end interface

  interface SphericalToCartesian
     module procedure SphericalToCartesianR4
     module procedure SphericalToCartesianREAL64
     module procedure SphericalToCartesianR4C2C
  end interface

  interface CartesianToSpherical
     module procedure CartesianToSphericalR4
     module procedure CartesianToSphericalREAL64
     module procedure CartesianToSphericalR4C2C
  end interface

  integer, parameter :: ntiles=6
  integer, parameter :: ndims=2
  integer, parameter :: r8=REAL64
  integer, parameter :: maxstring=120

  real(REAL64), parameter :: PI=3.14159265358979323846d0

! This EXTERNAL subroutine is in the fv directory
!  and has real*8 interfaces

  interface
     subroutine GetWeights(npx, npy, nlat, nlon, &
          index, weight, id1, id2, jdc, l2c,     &
          ee1, ee2, ff1, ff2, gg1, gg2,          &
           e1,  e2,  f1,  f2,  g1,  g2,          &
          sublons, sublats, AmNodeRoot, WriteNetcdf) 
       use, intrinsic :: iso_fortran_env, only: REAL64
       integer,  intent(in   ) :: npx,  npy
       integer,  intent(in   ) :: nlon, nlat
       integer,  intent(  out) :: index(3,nlon,nlat)
       real(REAL64), intent(  out) :: weight(4,nlon,nlat)
       integer,  intent(  out) :: id1(npx,npy)
       integer,  intent(  out) :: id2(npx,npy)
       integer,  intent(  out) :: jdc(npx,npy)
       real(REAL64), intent(  out) :: l2c(4,npx,npy)
       real(REAL64), intent(  out) :: ee1(npx,npy,3)
       real(REAL64), intent(  out) :: ee2(npx,npy,3)
       real(REAL64), intent(  out) :: ff1(npx,npy,3)
       real(REAL64), intent(  out) :: ff2(npx,npy,3)
       real(REAL64), intent(  out) :: gg1(npx,npy,3)
       real(REAL64), intent(  out) :: gg2(npx,npy,3)
       real(REAL64), pointer       ::  e1(:,:,:)
       real(REAL64), pointer       ::  e2(:,:,:)
       real(REAL64), pointer       ::  f1(:,:,:)
       real(REAL64), pointer       ::  f2(:,:,:)
       real(REAL64), pointer       ::  g1(:,:,:)
       real(REAL64), pointer       ::  g2(:,:,:)
       real(REAL64), intent(in)    :: sublons(:)
       real(REAL64), intent(in)    :: sublats(:)
       logical , optional      :: AmNodeRoot
       logical , optional      :: WriteNetcdf
     end subroutine GetWeights
  end interface

  interface
     subroutine GetWeightsC2C(npx, npy, npxout, npyout, index, weight, &
         ee1, ee2, ff1, ff2)
       use, intrinsic :: iso_fortran_env, only: REAL64
       integer,  intent(in   ) :: npx,  npy
       integer,  intent(in   ) :: npxout, npyout
       integer,  intent(  out) :: index(:,:,:,:)
       real(REAL64), intent(  out) :: weight(:,:,:,:)
       real(REAL64), intent(  out) :: ee1(:,:,:)
       real(REAL64), intent(  out) :: ee2(:,:,:)
       real(REAL64), intent(  out) :: ff1(:,:,:)
       real(REAL64), intent(  out) :: ff2(:,:,:)
     end subroutine GetWeightsC2C
  end interface

!JK patch for conservative interp----------------------------
integer , save   :: NT_tiles=-1
integer , allocatable, dimension(:) :: Tile_LL_ii, Tile_LL_jj
integer , allocatable, dimension(:) :: Tile_CS_nx, Tile_CS_ny
real(REAL64), allocatable, dimension(:) :: Tile_LL_wfrac, Tile_CS_wfrac
real(REAL64), allocatable, dimension(:) :: Tile_area, Tile_x, Tile_y
logical, save                 :: DO_CONSERVATIVE=.false.
!JK patch for conservative interp----------------------------


interface deallocGlob
   module procedure deallocGlob_i4_2
   module procedure deallocGlob_i4_3
   module procedure deallocGlob_i4_4
   module procedure deallocGlob_r8_3
   module procedure deallocGlob_r8_4
end interface deallocGlob

interface deallocLoc
   module procedure deallocLocl_r8_3
end interface deallocLoc

contains

  subroutine CubeLatLonDestroy( Tr, rc)
    type(T_CubeLatLonTransform), intent(inout) :: Tr
    integer, optional,           intent(  out) :: rc

    integer :: status

    call MAPL_SyncSharedMemory(rc=STATUS)
    VERIFY_(STATUS)
    DEALLOCGLOB_(Tr%index)
    DEALLOCGLOB_(Tr%weight)
    DEALLOCGLOB_(Tr%l2c)
    DEALLOCGLOB_(Tr%id1)
    DEALLOCGLOB_(Tr%id2)
    DEALLOCGLOB_(Tr%jdc)
    DEALLOCGLOB_(Tr%ee1)
    DEALLOCGLOB_(Tr%ee2)
    DEALLOCGLOB_(Tr%ff1)
    DEALLOCGLOB_(Tr%ff2)
    DEALLOCGLOB_(Tr%gg1)
    DEALLOCGLOB_(Tr%gg2)

    DEALLOCLOCL_(Tr%elon)
    DEALLOCLOCL_(Tr%elat)
    DEALLOCLOCL_(Tr%e1)
    DEALLOCLOCL_(Tr%e2)
    DEALLOCLOCL_(Tr%f1)
    DEALLOCLOCL_(Tr%f2)
    DEALLOCLOCL_(Tr%g1)
    DEALLOCLOCL_(Tr%g2)

    Tr%Created = .false.

    if (Tr%lsCreated) then
!ALT: if we created LocStream tile transforms we should destroy them
!     unfortunately MAPL does not destroy LocationStreams yet 
    end if

    RETURN_(SUCCESS)
  end subroutine CubeLatLonDestroy

  logical function CubeLatLonIsCreated(Tr)

    type(T_CubeLatLonTransform), intent(in ) :: Tr

    CubeLatLonIsCreated = Tr%Created

  end function CubeLatLonIsCreated

  subroutine CubeLatLonSubset(Tr,doSubset)
    logical,                     intent(in ) :: doSubset
    type(T_CubeLatLonTransform)              :: Tr
    Tr%subset=dosubset
  end subroutine CubeLatLonSubset

  function CubeLatLonCreate( npx, npy, nlon, nlat, lons, lats, local_ij, rc ) result(Tr)

    integer,                     intent(in ) :: npx,  npy
    integer,                     intent(in ) :: nlon, nlat
    real(REAL64),                    intent(in ) :: lons(:), lats(:)
    integer, optional, intent(in) :: local_ij(2,2)
    type(T_CubeLatLonTransform)              :: Tr
    integer, optional,           intent(out) :: rc


! npx      : inner dimension of global cube arrays (number of cells along cube edge) 
! npy      : outer dimension of global cube arrays ( 6*npx )
! nlon     : inner dimension of global LL arrays (Number of longitude points)
! nlat     : outer dimension of global LL arrays (Number of latitude points)
! lons     : the local nlon longitudes of LL grid, in radians
! lats     : the local nlat latitudes of LL grid, in radians
! local_ij : bounding box in (i,j) for local portion of lat-lon grid
! Tr       : The structure that holds the output transform
! rc       : return code

! Creates all necessary data to transform fields between Cube and LatLon
! grids, in bothe directions.  Data is stored in the output transform Tr. The
! transforms, Tr,is transposable by the MAPL transforming routines.

! Locals
!-------

    integer :: npts, status
    integer :: i, j
    real(REAL64), allocatable :: slon(:), slat(:)
    real(REAL64), allocatable :: clon(:), clat(:)

    integer :: i0, i1, j0, j1

! global vector rotations to be copied into local Tr versions
!    real(REAL64), pointer   :: ee1(:,:,:) => null()
!    real(REAL64), pointer   :: ee2(:,:,:) => null()
!    real(REAL64), pointer   :: ff1(:,:,:) => null()
!    real(REAL64), pointer   :: ff2(:,:,:) => null()
!    real(REAL64), pointer   :: gg1(:,:,:) => null()
!    real(REAL64), pointer   :: gg2(:,:,:) => null()

! Begin
!------

    ASSERT_(.not.Tr%Created)

    npts = npx + 1

    write(Tr%name,'(i5.5,"x",i5.5,"_c2l_",i5.5,"x",i5.5)') npx,npy,nlon,nlat
!    write(*,'(i5.5,"x",i5.5,"_c2l_",i5.5,"x",i5.5)') npx,npy,nlon,nlat

    Tr%npx  = npx
    Tr%npy  = npy
    Tr%nlon = nlon
    Tr%nlat = nlat

  ! allocate storage for weights and indeces for C2L
  !-------------------------------------------------

    DEALLOCGLOB_(Tr%index)
    DEALLOCGLOB_(Tr%weight)
    DEALLOCGLOB_(Tr%l2c)
    DEALLOCGLOB_(Tr%id1)
    DEALLOCGLOB_(Tr%id2)
    DEALLOCGLOB_(Tr%jdc)

    DEALLOCLOCL_(Tr%elon)
    DEALLOCLOCL_(Tr%elat)
    DEALLOCLOCL_(Tr%e1)
    DEALLOCLOCL_(Tr%e2)
    DEALLOCLOCL_(Tr%f1)
    DEALLOCLOCL_(Tr%f2)
    DEALLOCLOCL_(Tr%g1)
    DEALLOCLOCL_(Tr%g2)

    call MAPL_AllocNodeArray(Tr%index,(/3,nlon,nlat/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%index(3,nlon,nlat),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%weight,(/4,nlon,nlat/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%weight(4,nlon,nlat),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%l2c,(/4,npx,npy/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%l2c(4,npx,npy),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%id1,(/npx,npy/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%id1(npx,npy),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%id2,(/npx,npy/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%id2(npx,npy),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%jdc,(/npx,npy/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%jdc(npx,npy),stat=status)
    VERIFY_(STATUS)

    allocate(Tr%elon(size(lons),size(lats),3),stat=STATUS)
    VERIFY_(STATUS)

    if (present(local_ij)) then
      i0 = local_ij(1,1)
      i1 = local_ij(2,1)
      j0 = local_ij(1,2)
      j1 = local_ij(2,2)
    else
      i0 = 1
      i1 = size(lons)
      j0 = 1
      j1 = size(lats)
    end if

    tr%elon_local => tr%elon(i0:i1,j0:j1,:)

    allocate(Tr%elat(size(lons),size(lats),3),stat=STATUS)
    VERIFY_(STATUS)
    tr%elat_local => tr%elat(i0:i1,j0:j1,:)

    call MAPL_AllocNodeArray(Tr%ee1,(/npx,npy,3/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%ee1(npx,npy,3),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%ee2,(/npx,npy,3/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%ee2(npx,npy,3),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%ff1,(/npx,npy,3/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%ff1(npx,npy,3),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%ff2,(/npx,npy,3/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%ff2(npx,npy,3),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%gg1,(/npx,npy,3/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%gg1(npx,npy,3),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%gg2,(/npx,npy,3/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%gg2(npx,npy,3),stat=status)
    VERIFY_(STATUS)

! Argument AmNodeRoot passed to GetWeights identifies if we're using SHMEM
! and then only the NodeRoot gets weights, otherwise everyone does
!-------------------------------------------------------------------------
    call GetWeights(npx, npy, nlat, nlon, Tr%index, Tr%weight, &
         Tr%id1, Tr%id2, Tr%jdc, Tr%l2c,  &
         Tr%ee1, Tr%ee2, Tr%ff1, Tr%ff2,  Tr%gg1, Tr%gg2, &
         Tr%e1, Tr%e2, Tr%f1, Tr%f2, Tr%g1, Tr%g2, lons, lats, &
         AmNodeRoot = (MAPL_AmNodeRoot .or. (.not. MAPL_ShmInitialized))   &
#ifdef WRITE_WEIGHTS_TO_FILE
         , WriteNetcdf = MAPL_am_I_root() &
#endif
         )

! Deallocate large global vector rotation transforms

    call MAPL_SyncSharedMemory(rc=STATUS)
    VERIFY_(STATUS)

!    DEALLOCGLOB_(Tr%ee1)
!    DEALLOCGLOB_(Tr%ee2)
!    DEALLOCGLOB_(Tr%ff1)
!    DEALLOCGLOB_(Tr%ff2)
!    DEALLOCGLOB_(Tr%gg1)
!    DEALLOCGLOB_(Tr%gg2)

! Cartesian to latlon spherical on latlon grid

       allocate(slat(size(lats)),clat(size(lats)))
       allocate(slon(size(lons)),clon(size(lons)))

       do j=1,size(lats)
          SLAT(j) = SIN(lats(j))
          CLAT(j) = COS(lats(j))
       end do

       do I=1,size(lons)
          SLON(I) = sin(lons(i) - PI)
          CLON(I) = cos(lons(i) - PI)
       end DO

       do j=1,size(lats)
          do I=1,size(lons)
             Tr%elon(I,J,1) = -SLON(I)
             Tr%elon(I,J,2) =  CLON(I)
             Tr%elon(I,J,3) =  0.0
             Tr%elat(I,J,1) = -SLAT(J)*CLON(I)
             Tr%elat(I,J,2) = -SLAT(J)*SLON(I)
             Tr%elat(I,J,3) =  CLAT(J)
          end do
       end do

       deallocate(slon,clon,slat,clat)

    Tr%Created=.true.

    RETURN_(SUCCESS)
  end function CubeLatLonCreate

  subroutine CubeCubeDestroy( Tr, rc)
    type(T_CubeCubeTransform), intent(inout) :: Tr
    integer, optional,         intent(  out) :: rc

    integer :: status

    DEALLOCGLOB_(Tr%weight)
    DEALLOCGLOB_(Tr%index)
    DEALLOCGLOB_(Tr%ee1)
    DEALLOCGLOB_(Tr%ee2)
    DEALLOCGLOB_(Tr%ff1)
    DEALLOCGLOB_(Tr%ff2)

    Tr%Created = .false.

    RETURN_(SUCCESS)
  end subroutine CubeCubeDestroy

  logical function CubeCubeIsCreated(Tr)

    type(T_CubeCubeTransform), intent(in ) :: Tr

    CubeCubeIsCreated = Tr%Created

  end function CubeCubeIsCreated

  function CubeCubeCreate( npx, npy, npxout, npyout, rc ) result(Tr)

    integer,                     intent(in ) :: npx,  npy
    integer,                     intent(in ) :: npxout,  npyout
    type(T_CubeCubeTransform)                :: Tr
    integer, optional,           intent(out) :: rc

! Locals
!-------

    integer :: npts, status
    integer, parameter :: ntiles=6

! Real*8 are needed to make fv calls.
!-----------------------------------

! Begin
!------

    ASSERT_(.not.Tr%Created)

!ALT    npts = npx + 1
    npts = npxout ! + 1

    write(Tr%name,'(i5.5,"x",i5.5,"_c2c_",i5.5,"x",i5.5)') npx,npy,npxout,npyout

    Tr%npx  = npx
    Tr%npy  = npy
    Tr%npxout  = npxout
    Tr%npyout  = npyout

  ! allocate storage for weights and indeces for C2C
  !-------------------------------------------------

    DEALLOCGLOB_(Tr%weight)
    DEALLOCGLOB_(Tr%index)
    DEALLOCGLOB_(Tr%ee1)
    DEALLOCGLOB_(Tr%ee2)
    DEALLOCGLOB_(Tr%ff1)
    DEALLOCGLOB_(Tr%ff2)
    
    ! ALT: index and weight are allocated at the output grid resolution
    call MAPL_AllocNodeArray(Tr%weight,(/4,npxout,npyout/6,6/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%weight(4,npxout,npyout/6,6),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%index,(/3,npxout,npyout/6,6/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%index(3,npxout,npyout/6,6),stat=status)
    VERIFY_(STATUS)

    ! ALT: ff1 and ff2 are allocated at the input grid resolution
    call MAPL_AllocNodeArray(Tr%ff1,(/npx,npy,3/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%ff1(npx,npy,3),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%ff2,(/npx,npy,3/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%ff2(npx,npy,3),stat=status)
    VERIFY_(STATUS)

    ! ALT: ee1 and ee2 are allocated at the output grid resolution
    call MAPL_AllocNodeArray(Tr%ee1,(/npxout,npyout,3/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%ee1(npxout,npyout,3),stat=status)
    VERIFY_(STATUS)

    call MAPL_AllocNodeArray(Tr%ee2,(/npxout,npyout,3/),rc=STATUS)
    if(STATUS==MAPL_NoShm) allocate(Tr%ee2(npxout,npyout,3),stat=status)
    VERIFY_(STATUS)

    if (MAPL_AmNodeRoot .or. (.not. MAPL_ShmInitialized)) then
       call GetWeightsC2C(npx, npy, npxout, npyout, Tr%index, Tr%weight, &
            Tr%ee1, Tr%ee2, Tr%ff1, Tr%ff2)
    end if
    call MAPL_SyncSharedMemory(rc=STATUS)
    VERIFY_(STATUS)

       
    Tr%Created=.true.

    RETURN_(SUCCESS)
  end function CubeCubeCreate

  subroutine CubeToCube(Tr, data_cs_in, data_cs_out, rc)

    type(T_CubeCubeTransform),    intent(in )   :: Tr
    real,                         intent(inout) :: data_cs_in(:,:)
    real,                         intent(inout) :: data_cs_out(:,:)
    integer, optional,            intent(out)   :: rc

! Locals
!-------

    integer               :: nx,j1,j2,status,itile
    real(REAL64), allocatable :: var_cs_in(:,:,:), var_cs_out(:,:,:)

    ASSERT_(Tr%Created)

    nx   = Tr%npx

    !--------------------------------------------------------------------!
    ! perform interpolation                                              !
    !--------------------------------------------------------------------!

    allocate ( var_cs_in(0:nx+1,0:nx+1,ntiles),stat=status)
    VERIFY_(STATUS)

    var_cs_in=0.0

    do itile=1,ntiles
       j1 = nx*(itile-1) + 1
       j2 = nx*(itile-1) + nx
       var_cs_in(1:nx,1:nx,itile)=data_cs_in(:,j1:j2)
    end do

    nx   = Tr%npxout
    allocate ( var_cs_out(0:nx+1,0:nx+1,ntiles),stat=status)
    VERIFY_(STATUS)
    var_cs_out=0.0

    call C2CInterp(var_cs_in, var_cs_out, Tr%index, Tr%weight)

    do itile=1,ntiles
       j1 = nx*(itile-1) + 1
       j2 = nx*(itile-1) + nx
       data_cs_out(:,j1:j2) = var_cs_out(1:nx,1:nx,itile)
    end do
       
    deallocate ( var_cs_in, var_cs_out, stat=status)
    VERIFY_(STATUS)

    RETURN_(SUCCESS)
  end subroutine CubeToCube

  subroutine C2CInterp(var_in, var_out, index_c2c, weight_c2c)

    !------------------------------------------------------------------!
    ! do bilinear interpolation from cubed sphere to latlon grid       !
    ! using precalculated weights from get_weight                  !
    !------------------------------------------------------------------!

    real(REAL64), dimension(0:,0:,:), intent(inout) :: var_in
    real(REAL64), dimension(0:,0:,:), intent(inout) :: var_out
    real(REAL64), dimension(:,:,:,:), intent(in   ) :: weight_c2c
    integer,      dimension(:,:,:,:), intent(in   ) :: index_c2c

    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!

    integer           :: i, j, l, jx, ic, jc, lc, nx, ny

    nx   = size(var_out,1)-2
    ny   = size(var_out,2)-2


    call GhostCube(var_in)

    FACES: do l=1,ntiles
       JLOOP: do j=1,ny
          ILOOP: do i=1,nx
             ic=index_c2c(1,i,j,l)
             jc=index_c2c(2,i,j,l)
             lc=index_c2c(3,i,j,l)

             var_out(i,j,l)=weight_c2c(1,i,j,l)*var_in(ic  ,jc  , lc)  &
                           +weight_c2c(2,i,j,l)*var_in(ic  ,jc+1, lc)  &
                           +weight_c2c(3,i,j,l)*var_in(ic+1,jc+1, lc)  &
                           +weight_c2c(4,i,j,l)*var_in(ic+1,jc  , lc)
          enddo ILOOP
       enddo JLOOP
    enddo FACES

    return
  end subroutine C2CInterp


  subroutine get_conservative_weights(name)

    real(REAL64), allocatable :: FF_LL(:), FF_CS(:), CS_data(:,:), LL_data(:,:)
    character(len=120)  :: name
    integer             :: n, n_temp, nlon_, nlat_, nx_g, ny_g, UNIT
    integer             :: STATUS
    type (ESMF_VM)      :: vm

    DO_CONSERVATIVE = .true.
    call ESMF_VMGetCurrent(vm, rc=status)

!---read the weight fractions and connectivities---------------------------
    UNIT = GETFILE(name, DO_OPEN=0, ALL_PES=.true., RC=STATUS)
    open (UNIT=UNIT, FILE=NAME)
    if ( MAPL_am_I_root() ) read(UNIT,*) NT_tiles, n_temp, n_temp
    if ( MAPL_am_I_root() ) read(UNIT,*) n_temp
    if ( MAPL_am_I_root() ) read(UNIT,*) name
    if ( MAPL_am_I_root() ) read(UNIT,*) nlon_
    if ( MAPL_am_I_root() ) read(UNIT,*) nlat_
    if ( MAPL_am_I_root() ) read(UNIT,*) name
    if ( MAPL_am_I_root() ) read(UNIT,*) nx_g
    if ( MAPL_am_I_root() ) read(UNIT,*) ny_g
    call MAPL_CommsBcast(vm, DATA=NT_tiles, N=1, ROOT=0, RC=status)
    call MAPL_CommsBcast(vm, DATA=nlon_, N=1, ROOT=0, RC=status)
    call MAPL_CommsBcast(vm, DATA=nlat_, N=1, ROOT=0, RC=status)
    call MAPL_CommsBcast(vm, DATA=nx_g , N=1, ROOT=0, RC=status)
    call MAPL_CommsBcast(vm, DATA=ny_g , N=1, ROOT=0, RC=status)
    allocate (Tile_LL_ii(NT_tiles), Tile_LL_jj(NT_tiles))
    allocate (Tile_CS_nx(NT_tiles), Tile_CS_ny(NT_tiles))
    allocate (Tile_LL_wfrac(NT_tiles), Tile_CS_wfrac(NT_tiles))
    allocate (Tile_area(NT_tiles), Tile_x(NT_tiles), Tile_y(NT_tiles))
    do n=1,NT_tiles
       if ( MAPL_am_I_root() ) then
       read(UNIT,*) n_temp, Tile_area(n),  Tile_x(n),     Tile_y(n),       &
                            Tile_LL_ii(n), Tile_LL_jj(n), Tile_LL_wfrac(n),&
                    n_temp, Tile_CS_nx(n), Tile_CS_ny(n), Tile_CS_wfrac(n),&
                    n_temp 
       endif
    enddo
    call MAPL_CommsBcast(vm, DATA=Tile_LL_ii, N=NT_tiles, ROOT=0, RC=status)
    call MAPL_CommsBcast(vm, DATA=Tile_LL_jj, N=NT_tiles, ROOT=0, RC=status)
    call MAPL_CommsBcast(vm, DATA=Tile_CS_nx, N=NT_tiles, ROOT=0, RC=status)
    call MAPL_CommsBcast(vm, DATA=Tile_CS_ny, N=NT_tiles, ROOT=0, RC=status)
    call MAPL_CommsBcast(vm, DATA=Tile_LL_wfrac, N=NT_tiles, ROOT=0, RC=status)
    call MAPL_CommsBcast(vm, DATA=Tile_CS_wfrac, N=NT_tiles, ROOT=0, RC=status)
    call FREE_FILE(UNIT)

!---normalizing the weight fractions--------------------------------------
    allocate(LL_data(nlon_, nlat_), CS_data(nx_g, ny_g))
    allocate(FF_CS(nt_tiles), FF_LL(nt_tiles))
    CS_data(:,:)=0.; LL_data(:,:)=0.
    do n=1,ntiles
      CS_data(Tile_CS_nx(n),Tile_CS_ny(n))=&
      CS_data(Tile_CS_nx(n),Tile_CS_ny(n))+Tile_CS_wfrac(n)
      LL_data(Tile_LL_ii(n),Tile_LL_jj(n))=&
      LL_data(Tile_LL_ii(n),Tile_LL_jj(n))+Tile_LL_wfrac(n)
    enddo
    FF_CS=0.; FF_LL=0.
    do n=1,ntiles
      FF_CS(n)=CS_data(Tile_CS_nx(n),Tile_CS_ny(n))
      FF_LL(n)=LL_data(Tile_LL_ii(n),Tile_LL_jj(n))
    enddo
    do n=1,ntiles
      Tile_CS_wfrac(n)=Tile_CS_wfrac(n)/FF_CS(n)
      Tile_LL_wfrac(n)=Tile_LL_wfrac(n)/FF_LL(n)
    enddo
    deallocate(FF_CS, FF_LL)
    deallocate(LL_data, CS_data)

  end subroutine get_conservative_weights


  subroutine CubeToLatLonr8( Tr, data_cs, data_ll, transpose, misval, rc)

    type(T_CubeLatLonTransform),  intent(in )   :: Tr
    real(REAL64),                     intent(inout) :: data_cs(:,:)
    real(REAL64),                     intent(inout) :: data_ll(:,:)
    logical, optional,            intent(in )   :: transpose
    real,    optional,            intent(in )   :: misval
    integer, optional,            intent(out)   :: rc

! Locals
!-------

    integer               :: nx,j1,j2,status,itile
    real(REAL64), allocatable :: var_cs(:,:,:)
    
    ASSERT_(Tr%Created)

    nx   = Tr%npx

    !--------------------------------------------------------------------!
    ! perform interpolation                                              !
    !--------------------------------------------------------------------!

    allocate ( var_cs(0:nx+1,0:nx+1,ntiles),stat=status)
    VERIFY_(STATUS)

    ASSERT_(.not. (transpose .and. present(misval)))

    var_cs=0.0

    if(.not.transpose) then
       data_ll=0.0
       do itile=1,ntiles
          j1 = nx*(itile-1) + 1
          j2 = nx*(itile-1) + nx
          var_cs(1:nx,1:nx,itile)=data_cs(:,j1:j2)
       enddo
    end if

    call C2LInterp(var_cs, data_ll, Tr%index, Tr%weight,&
          misval, Tr%subset, transpose)

    if(transpose) then
       do itile=1,ntiles
          j1 = nx*(itile-1) + 1
          j2 = nx*(itile-1) + nx
          data_cs(:,j1:j2) = var_cs(1:nx,1:nx,itile)
       enddo
    end if

    deallocate ( var_cs ,stat=status)
    VERIFY_(STATUS)

    RETURN_(SUCCESS)
  end subroutine CubeToLatLonr8

  subroutine CubeToLatLonr4( Tr, data_cs, data_ll, transpose, misval, rc)

    type(T_CubeLatLonTransform),  intent(in )   :: Tr
    real,                         intent(inout) :: data_cs(:,:)
    real,                         intent(inout) :: data_ll(:,:)
    logical, optional,            intent(in )   :: transpose
    real,    optional,            intent(in )   :: misval
    integer, optional,            intent(out)   :: rc

! Locals
!-------

    integer               :: nx,j1,j2,status,itile
    real(REAL64), allocatable :: var_cs(:,:,:), data_ll8(:,:)

!JK patch for conservative interp---------------
    real(REAL64), allocatable :: data_cs8(:,:)


    ASSERT_(Tr%Created)

    nx   = Tr%npx

    !--------------------------------------------------------------------!
    ! perform interpolation                                              !
    !--------------------------------------------------------------------!

    allocate ( var_cs(0:nx+1,0:nx+1,ntiles),stat=status )
    VERIFY_(STATUS)
    allocate ( data_ll8(size(data_ll,1),size(data_ll,2)),stat=status)
    VERIFY_(STATUS)

    if (DO_CONSERVATIVE) then
       allocate ( data_cs8(size(data_cs,1),size(data_cs,2)),stat=status )
       VERIFY_(STATUS)
    endif

    ASSERT_(.not. (transpose .and. present(misval)))

    var_cs=0.0

    if(.not.transpose) then
       do itile=1,ntiles
          j1 = nx*(itile-1) + 1
          j2 = nx*(itile-1) + nx
          var_cs(1:nx,1:nx,itile)=data_cs(:,j1:j2)
       enddo
       if (DO_CONSERVATIVE) data_cs8=data_cs  !JK for conservative interp---
    else
       data_ll8=data_ll
    end if

    if (DO_CONSERVATIVE) then  !JK for conservative interp---

    if(.not.transpose) then
      call CToL_interp  &
         (data_ll8, data_cs8, NT_Tiles, Tr%nlat, Tr%nlon, Tr%npx, Tr%npy)
    else    
      call CToL_interp_b&
         (data_ll8, data_cs8, NT_Tiles, Tr%nlat, Tr%nlon, Tr%npx, Tr%npy)
    endif

    else                       !JK for conservative interp---

    call C2LInterp(var_cs, data_ll8, Tr%index, Tr%weight,&
                   misval,  Tr%subset, transpose)

    endif                      !JK for conservative interp---

    if(transpose) then
       do itile=1,ntiles
          j1 = nx*(itile-1) + 1
          j2 = nx*(itile-1) + nx
          data_cs(:,j1:j2) = var_cs(1:nx,1:nx,itile)
       enddo
       if (DO_CONSERVATIVE) data_cs=data_cs8  !JK for conservative interp---
    else
       data_ll=data_ll8
    end if

    if (DO_CONSERVATIVE) deallocate ( data_cs8 ) !JK for conservative interp---

    deallocate ( var_cs, data_ll8 )

    RETURN_(SUCCESS)
  end subroutine CubeToLatLonr4

  subroutine LatLonToCuber8( Tr, data_ll, data_cs, transpose, misval, rc)

    type(T_CubeLatLonTransform),  intent(in )   :: Tr
    real(REAL64),                     intent(inout) :: data_ll(:,:)
    real(REAL64),                     intent(inout) :: data_cs(:,:)
    logical, optional,            intent(in )   :: transpose
    real,    optional,            intent(in )   :: misval
    integer, optional,            intent(out)   :: rc

! Locals
!-------

    integer               :: nx,j1,j2,status,itile
    real(REAL64), allocatable :: var_cs(:,:,:)

    ASSERT_(Tr%Created)

    nx   = Tr%npx

    !--------------------------------------------------------------------!
    ! perform interpolation                                              !
    !--------------------------------------------------------------------!

    allocate ( var_cs(0:nx+1,0:nx+1,ntiles),stat=status)
    VERIFY_(STATUS)

    var_cs=0.

    if(transpose) then
       data_ll=0.
       do itile=1,ntiles
          j1 = nx*(itile-1) + 1
          j2 = nx*(itile-1) + nx
          var_cs(1:nx,1:nx,itile) = data_cs(:,j1:j2)
       enddo
    end if

    call L2CInterp(data_ll, var_cs, Tr%id1,  Tr%id2,  Tr%jdc, &
                   Tr%l2c, misval, transpose)

    if(.not.transpose) then
       do itile=1,ntiles
          j1 = nx*(itile-1) + 1
          j2 = nx*(itile-1) + nx
          data_cs(:,j1:j2) = var_cs(1:nx,1:nx,itile)
       enddo
    end if

    deallocate ( var_cs ,STAT=STATUS)
    VERIFY_(STATUS)

    RETURN_(SUCCESS)
  end subroutine LatLonToCuber8

  subroutine LatLonToCuber4( Tr, data_ll, data_cs, transpose, misval, rc)

    type(T_CubeLatLonTransform),  intent(in )   :: Tr
    real,                         intent(inout) :: data_ll(:,:)
    real,                         intent(inout) :: data_cs(:,:)
    logical, optional,            intent(in )   :: transpose
    real,    optional,            intent(in )   :: misval
    integer, optional,            intent(out)   :: rc

! Locals
!-------

    integer               :: nx,j1,j2,status,itile
    real(REAL64), allocatable :: data_cs8(:,:,:), data_ll8(:,:)
!JK for conservative interp--------------
    real(REAL64), allocatable :: cs8_data(:,:)
!JK for conservative interp--------------


    ASSERT_(Tr%Created)

    nx   = Tr%npx

    allocate ( data_ll8(size(data_ll,1),size(data_ll,2)),stat=status)
    VERIFY_(STATUS)
    allocate ( data_cs8(0:nx+1,0:nx+1,ntiles),stat=status)
    VERIFY_(STATUS)

    if (DO_CONSERVATIVE) then   !JK for conservative interp--------------
      allocate ( cs8_data(size(data_cs,1),size(data_cs,2)),stat=status )
      VERIFY_(STATUS)
    endif

    !--------------------------------------------------------------------!
    ! perform interpolation                                              !
    !--------------------------------------------------------------------!

    data_cs8=0.

    if(.not.transpose) then
       data_ll8 = data_ll
    end if

    if(transpose) then
       data_ll8=0.
       do itile=1,ntiles
          j1 = nx*(itile-1) + 1
          j2 = nx*(itile-1) + nx
          data_cs8(1:nx,1:nx,itile) = data_cs(:,j1:j2)
       enddo

       if (DO_CONSERVATIVE) cs8_data=data_cs !JK for conservative interp---

    end if

    if (DO_CONSERVATIVE) then !JK for conservative interp---

    if(.not.transpose) then
      call LToC_interp &
         (cs8_data, data_ll8, NT_Tiles, Tr%nlat, Tr%nlon, Tr%npx, Tr%npy)
    else
      call LToC_interp_b&
         (cs8_data, data_ll8, NT_Tiles, Tr%nlat, Tr%nlon, Tr%npx, Tr%npy)
    endif

    else                     !JK for conservative interp---

    call L2CInterp(data_ll8, data_cs8, Tr%id1,  Tr%id2,  Tr%jdc, &
                   Tr%l2c, misval, transpose)

    endif                    !JK for conservative interp---

    if(.not.transpose) then
       do itile=1,ntiles
          j1 = nx*(itile-1) + 1
          j2 = nx*(itile-1) + nx
          data_cs(:,j1:j2) = data_cs8(1:nx,1:nx,itile)
       enddo

       if (DO_CONSERVATIVE) data_cs=cs8_data !JK for conservative interp---

    end if

    if(transpose) then
       data_ll = data_ll8
    end if

    if (DO_CONSERVATIVE) deallocate ( cs8_data )

    deallocate ( data_cs8 )
    deallocate ( data_ll8 )

    RETURN_(SUCCESS)
  end subroutine LatLonToCuber4

  subroutine C2LInterp(cubsph, latlon, index, weight, misval, subset, transpose)

    !------------------------------------------------------------------!
    ! do bilinear interpolation from cubed sphere to latlon grid       !
    ! using precalculated weights from get_weight                  !
    !------------------------------------------------------------------!

    real(REAL64), dimension(0:,0:,:), intent(inout) :: cubsph
    real(REAL64), dimension(:,:),     intent(inout) :: latlon
    real(REAL64), dimension(:,:,:),   intent(in)    :: weight
    integer,  dimension(:,:,:),   intent(in)    :: index
    real, optional, intent(in)    :: misval
    logical,                      intent(in)    :: subset
    logical,                      intent(in)    :: transpose

    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!

    integer           :: i, j, ic, jc, nx, ny, nlon, nlat, tile, ii
    real(REAL64)          :: ww

    nx   = size(cubsph,1)-2
    ny   = size(cubsph,2)-2


    if(.not.transpose) then
       
       call GhostCube(cubsph)

    else

       cubsph = 0.0

    endif ! not transpose

    nlon = size(latlon,1)
    nlat = size(latlon,2)

    if( (mod(nlon,2)/=0) .and. (.not.subset) ) then
       print *, "NLON not even in cubetolatlon. Stopping."
       stop
    endif

    FACES: do tile = 1,ntiles
       ILOOP: do i=1,nlon
          ! if subsetting do not flip as this can make no sense
          if (.not.subset) then 
             ii = mod(i - 1 + nlon/2,nlon) + 1
          else
             ii = i
          end if
          JLOOP: do j=1,nlat

             HAVE_POINT: if (tile==index(3,i,j)) then

                ic=index(1,i,j)
                jc=index(2,i,j)

                ADJOINT: if(.not.transpose) then
                   UNDEF: if( .not. present(misval)) then
                      latlon(ii,j) = weight(1,i,j)*cubsph(ic  ,jc  ,tile)  &
                                   + weight(2,i,j)*cubsph(ic  ,jc+1,tile)  &
                                   + weight(3,i,j)*cubsph(ic+1,jc+1,tile)  &
                                   + weight(4,i,j)*cubsph(ic+1,jc  ,tile)
                   else
                      ww          = 0.0
                      latlon(ii,j) = 0.0

                      if(cubsph(ic  ,jc  ,tile)/=misval) then
                         latlon(ii,j) = latlon(ii,j) + weight(1,i,j)*cubsph(ic  ,jc  ,tile)
                         ww           = ww           + weight(1,i,j)
                      end if

                      if(cubsph(ic  ,jc+1,tile)/=misval) then
                         latlon(ii,j) = latlon(ii,j) + weight(2,i,j)*cubsph(ic  ,jc+1,tile)
                         ww           = ww           + weight(2,i,j)
                      end if

                      if(cubsph(ic+1,jc+1,tile)/=misval) then
                         latlon(ii,j) = latlon(ii,j) + weight(3,i,j)*cubsph(ic+1,jc+1,tile)
                         ww           = ww           + weight(3,i,j)
                      end if

                      if(cubsph(ic+1,jc  ,tile)/=misval) then
                         latlon(ii,j) = latlon(ii,j) + weight(4,i,j)*cubsph(ic+1,jc  ,tile)
                         ww           = ww           + weight(4,i,j)
                      end if

                      if(ww==0.0) then
                         latlon(ii,j) = misval
                      else
                         latlon(ii,j) = latlon(ii,j) / ww
                      end if
                   end if UNDEF

                else  ! Transpose
                   
                   cubsph(ic  ,jc  ,tile)=cubsph(ic  ,jc  ,tile)+weight(1,i,j)*latlon(ii,j)
                   cubsph(ic  ,jc+1,tile)=cubsph(ic  ,jc+1,tile)+weight(2,i,j)*latlon(ii,j)
                   cubsph(ic+1,jc+1,tile)=cubsph(ic+1,jc+1,tile)+weight(3,i,j)*latlon(ii,j)
                   cubsph(ic+1,jc  ,tile)=cubsph(ic+1,jc  ,tile)+weight(4,i,j)*latlon(ii,j)

                end if ADJOINT
             endif HAVE_POINT

          enddo JLOOP
       enddo ILOOP
    end do FACES


    if(transpose) then
       call GhostCubeT(cubsph)
    endif

    return
  end subroutine C2LInterp

  subroutine L2CInterp(latlon, cubsph, id1, id2, jdc, weight, misval, transpose)

    !------------------------------------------------------------------!
    ! do bilinear interpolation from cubed sphere to latlon grid       !
    ! using precalculated weights from get_weight                  !
    !------------------------------------------------------------------!

    real(REAL64), dimension(0:,0:,:), intent(inout) :: cubsph
    real(REAL64), dimension(:,:),   intent(inout) :: latlon
    real(REAL64), dimension(:,:,:), intent(in)    :: weight
    integer,  dimension(:,:),   intent(in)    :: id1, id2, jdc
    real, optional, intent(in)    :: misval
    logical,                    intent(in)    :: transpose

    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!

    integer       :: i, j, i1, i2, j1, nx, ny, nz, jx, k, nlon
    real(REAL64)      :: ww

    nx   = size(cubsph,1)-2
    ny   = size(cubsph,2)-2
    nz   = size(cubsph,3)
    nlon = size(latlon,1)

    if(transpose) latlon = 0.0    

    FACES: do k=1,nz
       FACE_Y: do jx=1,ny
          FACE_X: do i=1,nx

             j  = (k-1)*ny + jx
             j1 = jdc(i,j)

             i1 = id1(i,j)
             i2 = id2(i,j)

             ADJOINT: if(.not.transpose) then    

                UNDEF: if(.not. present(misval)) then
                   cubsph(i,jx,k) = weight(1,i,j)*latlon(i1,j1  )   &
                                  + weight(2,i,j)*latlon(i2,j1  )   &
                                  + weight(3,i,j)*latlon(i2,j1+1)   &
                                  + weight(4,i,j)*latlon(i1,j1+1)
                else
                   ww          = 0.0
                   cubsph(i,jx,k) = 0.0

                   if(latlon(i1,j1  )/=misval) then
                      cubsph(i,jx,k) = cubsph(i,jx,k) + weight(1,i,j)*latlon(i1,j1  )
                      ww             = ww             + weight(1,i,j)
                   end if
                   if(latlon(i2,j1  )/=misval) then
                      cubsph(i,jx,k) = cubsph(i,jx,k) + weight(2,i,j)*latlon(i2,j1  )
                      ww             = ww             + weight(2,i,j)
                   end if
                   if(latlon(i2,j1+1)/=misval) then
                      cubsph(i,jx,k) = cubsph(i,jx,k) + weight(3,i,j)*latlon(i2,j1+1)
                      ww             = ww             + weight(3,i,j)
                   end if
                   if(latlon(i1,j1+1)/=misval) then
                      cubsph(i,jx,k) = cubsph(i,jx,k) + weight(4,i,j)*latlon(i1,j1+1)
                      ww             = ww             + weight(4,i,j)
                   end if

                   if(ww==0.0) then
                      cubsph(i,jx,k) = misval
                   else
                      cubsph(i,jx,k) = cubsph(i,jx,k) / ww
                   end if

                end if UNDEF

             else

                latlon(i1,j1  ) = latlon(i1,j1  ) + weight(1,i,j)*cubsph(i,jx,k)
                latlon(i2,j1  ) = latlon(i2,j1  ) + weight(2,i,j)*cubsph(i,jx,k)
                latlon(i2,j1+1) = latlon(i2,j1+1) + weight(3,i,j)*cubsph(i,jx,k)
                latlon(i1,j1+1) = latlon(i1,j1+1) + weight(4,i,j)*cubsph(i,jx,k)

             end if ADJOINT

          enddo FACE_X
       enddo FACE_Y
    end do FACES

    return
  end subroutine L2CInterp

  
  subroutine GhostCube(x)
    real(REAL64), intent(INOUT) :: x(0:,0:,:)

    integer :: nx, ny
    
    nx = size(x,1)-2
    ny = nx

    x(1:nx,0,   1) = x(1:nx,ny,   6)
    x(1:nx,ny+1,1) = x(1,ny:1:-1, 3)
    x(0,1:ny,   1) = x(nx:1:-1,ny,5)
    x(nx+1,1:ny,1) = x(1,1:ny,    2)

    x(1:nx,0,   2) = x(nx,ny:1:-1,6)
    x(1:nx,ny+1,2) = x(1:nx,1,    3)
    x(0,1:ny,   2) = x(nx,1:ny,   1)
    x(nx+1,1:ny,2) = x(nx:1:-1,1, 4)

    x(1:nx,0,   3) = x(1:nx,ny,   2)
    x(1:nx,ny+1,3) = x(1,ny:1:-1, 5)
    x(0,1:ny,   3) = x(nx:1:-1,ny,1)
    x(nx+1,1:ny,3) = x(1,1:ny,    4)

    x(1:nx,0,   4) = x(nx,ny:1:-1,2)
    x(1:nx,ny+1,4) = x(1:nx,1,    5)
    x(0,1:ny,   4) = x(nx,1:ny,   3)
    x(nx+1,1:ny,4) = x(nx:1:-1,1, 6)

    x(1:nx,0,   5) = x(1:nx,ny,   4)
    x(1:nx,ny+1,5) = x(1,ny:1:-1, 1)
    x(0,1:ny,   5) = x(nx:1:-1,ny,3)
    x(nx+1,1:ny,5) = x(1,1:ny,    6)

    x(1:nx,0,   6) = x(nx,ny:1:-1,4)
    x(1:nx,ny+1,6) = x(1:nx,1,    1)
    x(0,1:ny,   6) = x(nx,1:ny,   5)
    x(nx+1,1:ny,6) = x(nx:1:-1,1, 2)

!  Set ghost corners to zero (to prevent "uninitialized" error). 
!  These values should not be used in any calculation

    x(0,   ny+1,:) = 0.0
    x(0,      0,:) = 0.0
    x(nx+1,   0,:) = 0.0
    x(nx+1,ny+1,:) = 0.0

  end subroutine GhostCube

  subroutine GhostCubeT(x)
    real(REAL64), intent(INOUT) :: x(0:,0:,:)

    integer :: nx, ny
    
    nx = size(x,1)-2
    ny = nx

    x(1:nx,ny,   6) = x(1:nx,ny,   6) + x(1:nx,0,   1)
    x(1,ny:1:-1, 3) = x(1,ny:1:-1, 3) + x(1:nx,ny+1,1)
    x(nx:1:-1,ny,5) = x(nx:1:-1,ny,5) + x(0,1:ny,   1)
    x(1,1:ny,    2) = x(1,1:ny,    2) + x(nx+1,1:ny,1)
                                                                            
    x(nx,ny:1:-1,6) = x(nx,ny:1:-1,6) + x(1:nx,0,   2)
    x(1:nx,1,    3) = x(1:nx,1,    3) + x(1:nx,ny+1,2)
    x(nx,1:ny,   1) = x(nx,1:ny,   1) + x(0,1:ny,   2)
    x(nx:1:-1,1, 4) = x(nx:1:-1,1, 4) + x(nx+1,1:ny,2)
                                                                            
    x(1:nx,ny,   2) = x(1:nx,ny,   2) + x(1:nx,0,   3)
    x(1,ny:1:-1, 5) = x(1,ny:1:-1, 5) + x(1:nx,ny+1,3)
    x(nx:1:-1,ny,1) = x(nx:1:-1,ny,1) + x(0,1:ny,   3)
    x(1,1:ny,    4) = x(1,1:ny,    4) + x(nx+1,1:ny,3)
                                                                            
    x(nx,ny:1:-1,2) = x(nx,ny:1:-1,2) + x(1:nx,0,   4)
    x(1:nx,1,    5) = x(1:nx,1,    5) + x(1:nx,ny+1,4)
    x(nx,1:ny,   3) = x(nx,1:ny,   3) + x(0,1:ny,   4)
    x(nx:1:-1,1, 6) = x(nx:1:-1,1, 6) + x(nx+1,1:ny,4)
                                                                            
    x(1:nx,ny,   4) = x(1:nx,ny,   4) + x(1:nx,0,   5)
    x(1,ny:1:-1, 1) = x(1,ny:1:-1, 1) + x(1:nx,ny+1,5)
    x(nx:1:-1,ny,3) = x(nx:1:-1,ny,3) + x(0,1:ny,   5)
    x(1,1:ny,    6) = x(1,1:ny,    6) + x(nx+1,1:ny,5)
                                                                            
    x(nx,ny:1:-1,4) = x(nx,ny:1:-1,4) + x(1:nx,0,   6)
    x(1:nx,1,    1) = x(1:nx,1,    1) + x(1:nx,ny+1,6)
    x(nx,1:ny,   5) = x(nx,1:ny,   5) + x(0,1:ny,   6)
    x(nx:1:-1,1, 2) = x(nx:1:-1,1, 2) + x(nx+1,1:ny,6)

! Zero Halo

    x(   :,   0,:) = 0.0
    x(   :,ny+1,:) = 0.0
    x(   0,   :,:) = 0.0
    x(nx+1,   :,:) = 0.0

  end subroutine GhostCubeT

! Transforms a 2D vector on the surface of the sphere to
! a 3D Cartesian vector or 3D to 2D when (Inverse == True).
! It can also apply the transpose of either operation (Transpose==True).
! It can deal with 2D vectors defined either along Lat-Lon coordinates
! or along Cube-sphere grid lines (Cube==True) as defined by
! the FVCube dynamical core. The transforming coefficients for 
! all of these operations are kept in Trans and are precomputed when
! it is initialized for a particular pair of grids when one grid is cube
! and the other is lat-lon.


  subroutine SphericalToCartesianR4(Tr, U, V, Uxyz, Transpose, SphIsLL, Rotate, RC)
    type(T_CubeLatLonTransform), intent(IN ) :: Tr
    real,                        intent(IN ) :: U(:,:,:), V(:,:,:)
    real,                        intent(OUT) :: Uxyz(:,:,:)
    logical,                     intent(IN ) :: Transpose
    logical,                     intent(IN ) :: SphIsLL
    logical,                     intent(IN ) :: Rotate
    integer, optional,           intent(OUT) :: RC

    integer           :: K, LM
    real(REAL64), pointer :: e1(:,:,:), e2(:,:,:) 

    if(.not.Rotate) then
       ASSERT_(.not.Transpose .and. .not.SphIsLL)
    end if

    if(SphIsLL) then
       if (size(U,1) == size(Tr%elon,1) .and. &
            size(U,2) == size(Tr%elon,2) .and. &
            size(V,1) == size(Tr%elat,1) .and. &
            size(V,2) == size(Tr%elat,2)) then
          e1=>Tr%elon
          e2=>Tr%elat
       else
          e1=>Tr%elon_local
          e2=>Tr%elat_local
       end if
    else
       if(.not.Rotate) then
          if (size(U,1) == size(Tr%gg1,1) .and. &
              size(U,2) == size(Tr%gg1,2) .and. &
              size(V,1) == size(Tr%gg2,1) .and. &
              size(V,2) == size(Tr%gg2,2)) then
             e1=>Tr%gg1
             e2=>Tr%gg2
          else
             e1=>Tr%g1
             e2=>Tr%g2
          end if
       elseif(.not.Transpose) then
          if (size(U,1) == size(Tr%ff1,1) .and. &
              size(U,2) == size(Tr%ff1,2) .and. &
              size(V,1) == size(Tr%ff2,1) .and. &
              size(V,2) == size(Tr%ff2,2)) then
             e1=>Tr%ff1
             e2=>Tr%ff2
          else
             e1=>Tr%f1
             e2=>Tr%f2
          endif
       else
          e1=>Tr%e1
          e2=>Tr%e2
       end if
    end if

    LM = size(U,3)

    do k=1,LM
       where(U(:,:,k) == MAPL_UNDEF .or. V(:,:,k) == MAPL_UNDEF)
          Uxyz(:,:,k     ) = MAPL_UNDEF
          Uxyz(:,:,k+  LM) = MAPL_UNDEF
          Uxyz(:,:,k+2*LM) = MAPL_UNDEF
       elsewhere
          Uxyz(:,:,k     ) = U(:,:,k)*e1(:,:,1) + V(:,:,k)*e2(:,:,1)
          Uxyz(:,:,k+  LM) = U(:,:,k)*e1(:,:,2) + V(:,:,k)*e2(:,:,2)
          Uxyz(:,:,k+2*LM) = U(:,:,k)*e1(:,:,3) + V(:,:,k)*e2(:,:,3)
       end where
    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine SphericalToCartesianR4

  subroutine SphericalToCartesianREAL64(Tr, U, V, Uxyz, Transpose, SphIsLL, Rotate, RC)
    type(T_CubeLatLonTransform), intent(IN ) :: Tr
    real(REAL64),                    intent(IN ) :: U(:,:,:), V(:,:,:)
    real(REAL64),                    intent(OUT) :: Uxyz(:,:,:)
    logical,                     intent(IN ) :: Transpose
    logical,                     intent(IN ) :: SphIsLL
    logical,                     intent(IN ) :: Rotate
    integer, optional,           intent(OUT) :: RC

    integer           :: K, LM
    real(REAL64), pointer :: e1(:,:,:), e2(:,:,:) 

    if(.not.Rotate) then
       ASSERT_(.not.Transpose .and. .not.SphIsLL)
    end if

    if(SphIsLL) then
       if (size(U,1) == size(Tr%elon,1) .and. &
             size(U,2) == size(Tr%elon,2) .and. &
             size(V,1) == size(Tr%elat,1) .and. &
             size(V,2) == size(Tr%elat,2)) then
          e1=>Tr%elon
          e2=>Tr%elat
       else
          e1=>Tr%elon_local
          e2=>Tr%elat_local
       end if
     else
       if(.not.Rotate) then
          if (size(U,1) == size(Tr%gg1,1) .and. &
              size(U,2) == size(Tr%gg1,2) .and. &
              size(V,1) == size(Tr%gg2,1) .and. &
              size(V,2) == size(Tr%gg2,2)) then
             e1=>Tr%gg1
             e2=>Tr%gg2
          else
             e1=>Tr%g1
             e2=>Tr%g2
          end if
       elseif(.not.Transpose) then
          if (size(U,1) == size(Tr%ff1,1) .and. &
              size(U,2) == size(Tr%ff1,2) .and. &
              size(V,1) == size(Tr%ff2,1) .and. &
              size(V,2) == size(Tr%ff2,2)) then
             e1=>Tr%ff1
             e2=>Tr%ff2
          else
             e1=>Tr%f1
             e2=>Tr%f2
          endif
       else
          e1=>Tr%e1
          e2=>Tr%e2
       end if
    end if

    LM = size(U,3)

    do k=1,LM
       where(U(:,:,k) == MAPL_UNDEF .or. V(:,:,k) == MAPL_UNDEF)
          Uxyz(:,:,k     ) = MAPL_UNDEF
          Uxyz(:,:,k+  LM) = MAPL_UNDEF
          Uxyz(:,:,k+2*LM) = MAPL_UNDEF
       elsewhere
          Uxyz(:,:,k     ) = U(:,:,k)*e1(:,:,1) + V(:,:,k)*e2(:,:,1)
          Uxyz(:,:,k+  LM) = U(:,:,k)*e1(:,:,2) + V(:,:,k)*e2(:,:,2)
          Uxyz(:,:,k+2*LM) = U(:,:,k)*e1(:,:,3) + V(:,:,k)*e2(:,:,3)
       end where
    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine SphericalToCartesianREAL64

  subroutine CartesianToSphericalR4(Tr, Uxyz, U, V, Transpose, SphIsLL, Rotate, RC)
    type(T_CubeLatLonTransform), intent(IN ) :: Tr
    real,                        intent(OUT) :: U(:,:,:), V(:,:,:)
    real,                        intent(IN ) :: Uxyz(:,:,:)
    logical,                     intent(IN ) :: Transpose
    logical,                     intent(IN ) :: SphIsLL
    logical,                     intent(IN ) :: Rotate
    integer, optional,           intent(OUT) :: RC

    integer           :: K, LM
    real(REAL64), pointer :: e1(:,:,:), e2(:,:,:) 

    if(.not.Rotate) then
       ASSERT_(.not.Transpose .and. .not.SphIsLL)
    end if

    if(SphIsLL) then
       if (size(Uxyz,1) == size(Tr%elon,1) &
            .and. size(Uxyz,2) == size(Tr%elon,2) &
            .and. size(Uxyz,1) == size(Tr%elat,1) &
            .and. size(Uxyz,2) == size(Tr%elat,2)) then
          e1=>Tr%elon
          e2=>Tr%elat
       else
          e1=>Tr%elon_local
          e2=>Tr%elat_local
       end if
    else
       if(.not.Rotate) then
          if (size(Uxyz,1) == size(Tr%gg1,1) &
               .and. size(Uxyz,2) == size(Tr%gg1,2) &
               .and. size(Uxyz,1) == size(Tr%gg2,1) &
               .and. size(Uxyz,2) == size(Tr%gg2,2)) then
             e1=>Tr%gg1
             e2=>Tr%gg2
          else
             e1=>Tr%g1
             e2=>Tr%g2
          endif
       elseif(Transpose) then
          e1=>Tr%f1
          e2=>Tr%f2
       else
          if (size(Uxyz,1) == size(Tr%ee1,1) &
               .and. size(Uxyz,2) == size(Tr%ee1,2) &
               .and. size(Uxyz,1) == size(Tr%ee2,1) &
               .and. size(Uxyz,2) == size(Tr%ee2,2)) then
             e1=>Tr%ee1
             e2=>Tr%ee2
          else
             e1=>Tr%e1
             e2=>Tr%e2
          end if
       end if
    end if

    LM = size(U,3)
    
    do k=1,LM
       where (Uxyz(:,:,k) == MAPL_UNDEF .or. Uxyz(:,:,k+LM) == MAPL_UNDEF &
            .or. Uxyz(:,:,k+2*LM) == MAPL_UNDEF)
          U(:,:,k) = MAPL_UNDEF
          V(:,:,k) = MAPL_UNDEF
       elsewhere
 
          U(:,:,k) = Uxyz(:,:,k     )*e1(:,:,1) + &
                     Uxyz(:,:,k+  LM)*e1(:,:,2) + &
                     Uxyz(:,:,k+2*LM)*e1(:,:,3)

          V(:,:,k) = Uxyz(:,:,k     )*e2(:,:,1) + &
                     Uxyz(:,:,k+  LM)*e2(:,:,2) + &
                     Uxyz(:,:,k+2*LM)*e2(:,:,3)
       end where
    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine CartesianToSphericalR4

  subroutine CartesianToSphericalREAL64(Tr, Uxyz, U, V, Transpose, SphIsLL, Rotate, RC)
    type(T_CubeLatLonTransform), intent(IN ) :: Tr
    real(REAL64),                    intent(OUT) :: U(:,:,:), V(:,:,:)
    real(REAL64),                    intent(IN ) :: Uxyz(:,:,:)
    logical,                     intent(IN ) :: Transpose
    logical,                     intent(IN ) :: SphIsLL
    logical, optional,           intent(IN ) :: Rotate
    integer, optional,           intent(OUT) :: RC

    integer           :: K, LM
    real(REAL64), pointer :: e1(:,:,:), e2(:,:,:) 

    if(.not.Rotate) then
       ASSERT_(.not.Transpose .and. .not.SphIsLL)
    end if

    if(SphIsLL) then
       if (size(Uxyz,1) == size(Tr%elon,1) &
            .and. size(Uxyz,2) == size(Tr%elon,2) &
            .and. size(Uxyz,1) == size(Tr%elat,1) &
            .and. size(Uxyz,2) == size(Tr%elat,2)) then
          e1=>Tr%elon
          e2=>Tr%elat
       else
          e1=>Tr%elon_local
          e2=>Tr%elat_local
       end if
    else
       if(.not.Rotate) then
          if (size(Uxyz,1) == size(Tr%gg1,1) &
               .and. size(Uxyz,2) == size(Tr%gg1,2) &
               .and. size(Uxyz,1) == size(Tr%gg2,1) &
               .and. size(Uxyz,2) == size(Tr%gg2,2)) then
             e1=>Tr%gg1
             e2=>Tr%gg2
          else
             e1=>Tr%g1
             e2=>Tr%g2
          endif
       elseif(Transpose) then
          e1=>Tr%f1
          e2=>Tr%f2
       else
          if (size(Uxyz,1) == size(Tr%ee1,1) &
               .and. size(Uxyz,2) == size(Tr%ee1,2) &
               .and. size(Uxyz,1) == size(Tr%ee2,1) &
               .and. size(Uxyz,2) == size(Tr%ee2,2)) then
             e1=>Tr%ee1
             e2=>Tr%ee2
          else
             e1=>Tr%e1
             e2=>Tr%e2
          end if
       end if
    end if

    LM = size(U,3)
    
    do k=1,LM
       where (Uxyz(:,:,k) == MAPL_UNDEF .or. Uxyz(:,:,k+LM) == MAPL_UNDEF &
            .or. Uxyz(:,:,k+2*LM) == MAPL_UNDEF)
          U(:,:,k) = MAPL_UNDEF
          V(:,:,k) = MAPL_UNDEF
       elsewhere
 
          U(:,:,k) = Uxyz(:,:,k     )*e1(:,:,1) + &
                     Uxyz(:,:,k+  LM)*e1(:,:,2) + &
                     Uxyz(:,:,k+2*LM)*e1(:,:,3)

          V(:,:,k) = Uxyz(:,:,k     )*e2(:,:,1) + &
                     Uxyz(:,:,k+  LM)*e2(:,:,2) + &
                     Uxyz(:,:,k+2*LM)*e2(:,:,3)
       end where
    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine CartesianToSphericalREAL64

! create LS and grindIn
! create LS and grindOut
! create Xform
  subroutine MAPL_RegridLSCreate(TR, GridIn, GridOut, TILINGFILE, RC)
    type(T_CubeLatLonTransform), intent(INOUT) :: Tr
    type(ESMF_Grid),             intent(INout) :: GridIn
    type(ESMF_Grid),             intent(INout) :: GridOut
    character(len=*),            intent(IN   ) :: TILINGFILE
    integer, optional,           intent(  OUT) :: RC


    type(MAPL_LocStream) :: LocStIn
    type(MAPL_LocStream) :: LocStOut
    type(ESMF_DistGrid)  :: distgrid
    type(ESMF_DELayout)  :: layout
    integer :: status


    call ESMF_GridGet(GridIn, DistGrid=distgrid, rc=status)
    VERIFY_(STATUS)
    call ESMF_DistGridGet(distGRID, deLayout=layout, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamCreate(locStIn, LAYOUT=layout, FILENAME=TILINGFILE, &
                              NAME='LocStIn',                           &
                              grid=GridIn, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamCreate(locStOut, LAYOUT=layout, FILENAME=TILINGFILE, &
                              NAME='LocStOut',                           &
                              grid=GridOut, RC=STATUS)
    VERIFY_(STATUS)

   call MAPL_LocStreamCreateXform ( XFORM=TR%XFORMinout, &
                                    LocStreamOut=locStOut, &
                                    LocStreamIn=locstIn, &
                                    NAME='XFORM_InOut', &
                                    RC=STATUS )
   VERIFY_(STATUS)

   call MAPL_LocStreamCreateXform ( XFORM=TR%XFORMoutin, &
                                    LocStreamOut=locStIn, &
                                    LocStreamIn=locstOut, &
                                    NAME='XFORM_OutIn', &
                                    RC=STATUS )
   VERIFY_(STATUS)

   TR%lsCreated = .true.
   TR%locStIn = locStIn
   TR%locStOut = locStOut

   RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_RegridLSCreate


  subroutine RunTileTransform(TR, PTR2d_IN, PTR2d_OUT, TRANSPOSE, RC)
    type(T_CubeLatLonTransform), target, intent(IN ) :: Tr
    real                                     :: ptr2d_in(:,:)
    real                                     :: ptr2d_out(:,:)
    logical                                  :: transpose
    integer, optional        , intent(  OUT) :: RC
    

! Errlog vars
    integer :: status

! Local vars
    type(MAPL_LocStreamXform), pointer :: XFORM
    type(MAPL_LocStream),      pointer :: LS_IN, LS_OUT
    integer                            :: NTILES_IN, NTILES_OUT
    real, allocatable, dimension(:)    :: tile_in, tile_out

    LS_IN => TR%LocStIn
    LS_OUT => TR%LocStOut
    XFORM => TR%XFORMinout

! query ntiles
    call MAPL_LocStreamGet(LS_In, NT_LOCAL = ntiles_in, rc=status)
    VERIFY_(STATUS)
    call MAPL_LocStreamGet(LS_Out, NT_LOCAL = ntiles_out, rc=status)
    VERIFY_(STATUS)

    allocate(tile_in (ntiles_in ), stat=status)
    VERIFY_(STATUS)
    allocate(tile_out(ntiles_out), stat=status)
    VERIFY_(STATUS)

    if (.not.transpose) then
! forward run
! G2T
       call MAPL_LocStreamTransform(LS_IN, TILE_IN, PTR2d_IN, RC=STATUS)
       VERIFY_(STATUS)
! T2T
       call MAPL_LocStreamTransform(tile_out, XFORM, tile_in, RC=STATUS ) 
       VERIFY_(STATUS)
! T2G
       call MAPL_LocStreamTransform(LS_OUT, PTR2d_OUT, TILE_OUT, RC=STATUS)
       VERIFY_(STATUS)
    else 
! adjoint run
! T2G adj
       call MAPL_LocStreamTransform(LS_IN, PTR2d_IN, TILE_IN, &
            TRANSPOSE=transpose, RC=STATUS)
       VERIFY_(STATUS)
! T2T adj
       call MAPL_LocStreamTransform(tile_out, XFORM, tile_in, RC=STATUS ) 
       VERIFY_(STATUS)
! G2T adj
       call MAPL_LocStreamTransform(LS_OUT, TILE_OUT, PTR2d_OUT, &
            TRANSPOSE=transpose, RC=STATUS)
       VERIFY_(STATUS)
    endif

    deallocate(tile_out)
    deallocate(tile_in )

    RETURN_(ESMF_SUCCESS)
  end subroutine RunTileTransform

  subroutine TileTransformTranspose(TR, RC)
    type(T_CubeLatLonTransform), intent(INOUT) :: TR
    integer, optional          , intent(  OUT) :: RC

    type(MAPL_LocStream) :: ls
    type(MAPL_LocStreamXform) :: Xform

    TR%lsCreated = .false.

    ls = TR%locStIn
    TR%locStIn = TR%locStOut
    TR%locStOut = ls

    xform = TR%xformInOut
    TR%xformInOut = TR%xformOutIn
    TR%xformOutIn = xform

    RETURN_(ESMF_SUCCESS)
  end subroutine TileTransformTranspose

subroutine LToC_interp(CS_data, LL_data, NT_Tiles, nlat, nlon, nxg, nyg)
integer n, NT_Tiles, nlat, nlon, nxg, nyg
real (REAL64) :: LL_data(nlon,nlat), CS_data(nxg,nyg)
real (REAL64), allocatable :: T_data(:)
allocate(T_data(NT_tiles))
!LL2Tile-----------------------------------------------------
do n=1,NT_tiles
  T_data(n)=LL_data(Tile_LL_ii(n),Tile_LL_jj(n))
enddo
!Tile2Cube---------------------------------------------------
CS_data(:,:)=0.
do n=1,NT_tiles
  CS_data(Tile_CS_nx(n),Tile_CS_ny(n))=CS_data(Tile_CS_nx(n),Tile_CS_ny(n))+&
                                       T_data(n)*Tile_CS_wfrac(n)
enddo
deallocate(T_data)
end subroutine LToC_interp

subroutine CToL_interp(LL_data, CS_data, NT_Tiles, nlat, nlon, nxg, nyg)
integer n, NT_Tiles, nlat, nlon, nxg, nyg
real (REAL64) :: LL_data(nlon,nlat), CS_data(nxg,nyg)
real (REAL64), allocatable :: T_data(:)

allocate(T_data(NT_tiles))
!Cube2Tile---------------------------------------------------
T_data=0.
do n=1,NT_tiles
  T_data(n)=CS_data(Tile_CS_nx(n),Tile_CS_ny(n))
enddo
!Tile2LL-----------------------------------------------------
LL_data(:,:)=0.
do n=1,NT_tiles
  LL_data(Tile_LL_ii(n),Tile_LL_jj(n))=LL_data(Tile_LL_ii(n),Tile_LL_jj(n))+&
                                       T_data(n)*Tile_LL_wfrac(n)
enddo
deallocate(T_data)
end subroutine CToL_interp
SUBROUTINE LToC_INTERP_B(cs_datab, ll_datab, NT_tiles, nlat, nlon, nxg, nyg)
  IMPLICIT NONE
  INTEGER :: NT_tiles, nlat, nlon, nxg, nyg
  REAL (REAL64) :: ll_datab(nlon, nlat), cs_datab(nxg, nyg)
  REAL (REAL64) :: t_datab(NT_tiles)
  INTEGER :: n
  t_datab = 0.0
  DO n=1,NT_tiles
    t_datab(n)=t_datab(n)+tile_cs_wfrac(n)*cs_datab(tile_cs_nx(n),tile_cs_ny(n))
  END DO
  ll_datab = 0.0
  DO n=1,NT_tiles
    ll_datab(tile_ll_ii(n),tile_ll_jj(n))=ll_datab(tile_ll_ii(n),tile_ll_jj(n))+&
                                          t_datab(n)
    t_datab(n) = 0.0
  END DO
  cs_datab = 0.0
END SUBROUTINE LToC_INTERP_B
SUBROUTINE CToL_INTERP_B(ll_datab, cs_datab, NT_tiles, nlat, nlon, nxg, nyg)
  IMPLICIT NONE
  INTEGER :: NT_tiles, nlat, nlon, nxg, nyg
  REAL (REAL64) :: ll_datab(nlon, nlat), cs_datab(nxg, nyg)
  REAL (REAL64) :: t_datab(NT_tiles)
  INTEGER :: n
  t_datab = 0.0
  DO n=1,NT_tiles
    t_datab(n)=t_datab(n)+tile_ll_wfrac(n)*ll_datab(tile_ll_ii(n),tile_ll_jj(n))
  END DO
  cs_datab = 0.0
  DO n=1,NT_tiles
    cs_datab(tile_cs_nx(n),tile_cs_ny(n))=cs_datab(tile_cs_nx(n),tile_cs_ny(n))+&
                                          t_datab(n)
    t_datab(n) = 0.0
  END DO
  ll_datab = 0.0
END SUBROUTINE CToL_INTERP_B

subroutine ReStaggerWindsCube(U, V, D2A)
  real,                      intent(INOUT) :: U(:,:,:)
  real,                      intent(INOUT) :: V(:,:,:)
  logical,                   intent(IN   ) :: D2A
  
  integer :: im, jm, lm
  integer :: tile, l, j1, j2, status
  integer, parameter :: ntiles=6
  real, allocatable, dimension(:,:,:) :: UAP, VAP

  IM = size(U,1)
  JM = size(V,2)/6
  LM = size(U,3)
  
  allocate (UAP(0:IM+1,0:JM+1,ntiles), stat=status)
  allocate (VAP(0:IM+1,0:JM+1,ntiles), stat=status)


  do L = 1, LM
     do tile = 1,ntiles
        j1 = (tile-1)*JM + 1
        j2 = tile*JM 
        UAP(1:IM,1:JM,tile) = U(1:IM,j1:j2,L)
        VAP(1:IM,1:JM,tile) = V(1:IM,j1:j2,L)
     enddo
     call GhostCubeVector(UAP, VAP, D2A)

     if(D2A) then
        do tile = 1,ntiles
           j1 = (tile-1)*IM + 1
           j2 = tile*IM 
           U(1:IM,j1:j2,L) = 0.5*(UAP(1:IM,1:JM,tile) + UAP(1:IM,2:JM+1,tile))
           V(1:IM,j1:j2,L) = 0.5*(VAP(1:IM,1:JM,tile) + VAP(2:IM+1,1:JM,tile))
        end do
     else
        do tile = 1,ntiles
           j1 = (tile-1)*JM + 1
           j2 = tile*JM 

           U(1:IM,j1:j2,L) = 0.5*(UAP(1:IM,1:JM,tile) + UAP(1:IM,0:JM-1,tile))
           V(1:IM,j1:j2,L) = 0.5*(VAP(1:IM,1:JM,tile) + VAP(0:IM-1,1:JM,tile))
        end do
     end if
  end do

  deallocate(VAP)
  deallocate(UAP)

  return
end subroutine ReStaggerWindsCube

  subroutine SphericalToCartesianR4C2C(Tr, U, V, Uxyz)
    type(T_CubeCubeTransform), intent(IN ) :: Tr
    real,                      intent(IN ) :: U(:,:,:), V(:,:,:)
    real,                      intent(OUT) :: Uxyz(:,:,:)
    
    integer           :: K, LM
    real(REAL64), pointer :: e1(:,:,:), e2(:,:,:) 

    e1=>Tr%ff1
    e2=>Tr%ff2

    LM = size(U,3)
    
    do k=1,LM
       Uxyz(:,:,k     ) = U(:,:,k)*e1(:,:,1) + V(:,:,k)*e2(:,:,1)
       Uxyz(:,:,k+  LM) = U(:,:,k)*e1(:,:,2) + V(:,:,k)*e2(:,:,2)
       Uxyz(:,:,k+2*LM) = U(:,:,k)*e1(:,:,3) + V(:,:,k)*e2(:,:,3)
    end do

    return
  end subroutine SphericalToCartesianR4C2C

  subroutine CartesianToSphericalR4C2C(Tr, Uxyz, U, V)
    type(T_CubeCubeTransform),   intent(IN ) :: Tr
    real,                        intent(OUT) :: U(:,:,:), V(:,:,:)
    real,                        intent(IN ) :: Uxyz(:,:,:)

    integer           :: K, LM
    real(REAL64), pointer :: e1(:,:,:), e2(:,:,:) 

    e1=>Tr%ee1
    e2=>Tr%ee2

    LM = size(U,3)
    
    do k=1,LM
       U(:,:,k) = Uxyz(:,:,k     )*e1(:,:,1) + &
                  Uxyz(:,:,k+  LM)*e1(:,:,2) + &
                  Uxyz(:,:,k+2*LM)*e1(:,:,3)

       V(:,:,k) = Uxyz(:,:,k     )*e2(:,:,1) + &
                  Uxyz(:,:,k+  LM)*e2(:,:,2) + &
                  Uxyz(:,:,k+2*LM)*e2(:,:,3)
    end do

    return
  end subroutine CartesianToSphericalR4C2C

  subroutine GhostCubeVector(u,v,d2a)
    real,    intent(INOUT) :: u(0:,0:,:)
    real,    intent(INOUT) :: v(0:,0:,:)
    logical, intent(IN   ) :: d2a ! .true. when input is on D-grid

    integer :: nx, ny

! Note that the mapping between the "flipped" edges differs for
! points on the A and D grid. 
! For "center" vars, 1 <=> n; for "edge" vars, 1 <=> n+1

    nx = size(u,1)-2
    ny = nx

    ! FACE 1
    u(1:nx,0,   1) =  u(1:nx,ny,   6)
    u(1:nx,ny+1,1) = -v(1,ny:1:-1, 3)
    u(nx+1,1:ny,1) =  u(1,1:ny,    2)

    v(1:nx,0,   1) =  v(1:nx,ny,   6)
    v(0,1:ny,   1) = -u(nx:1:-1,ny,5)
    v(nx+1,1:ny,1) =  v(1,1:ny,    2)

    ! FACE 2
    u(1:nx,0,   2) = -v(nx,ny:1:-1,6)
    u(1:nx,ny+1,2) =  u(1:nx,1,    3)
    u(0,1:ny,   2) =  u(nx,1:ny,   1)

    v(1:nx,ny+1,2) =  v(1:nx,1,    3)
    v(0,1:ny,   2) =  v(nx,1:ny,   1)
    v(nx+1,1:ny,2) = -u(nx:1:-1,1, 4)

    ! FACE 3
    u(1:nx,0,   3) =  u(1:nx,ny,   2)
    u(1:nx,ny+1,3) = -v(1,ny:1:-1, 5)
    u(nx+1,1:ny,3) =  u(1,1:ny,    4)

    v(1:nx,0,   3) =  v(1:nx,ny,   2)
    v(0,1:ny,   3) = -u(nx:1:-1,ny,1)
    v(nx+1,1:ny,3) =  v(1,1:ny,    4)

    ! FACE 4
    u(1:nx,0,   4) = -v(nx,ny:1:-1,2)
    u(1:nx,ny+1,4) =  u(1:nx,1,    5)
    u(0,1:ny,   4) =  u(nx,1:ny,   3)

    v(1:nx,ny+1,4) =  v(1:nx,1,    5)
    v(0,1:ny,   4) =  v(nx,1:ny,   3)
    v(nx+1,1:ny,4) = -u(nx:1:-1,1, 6)

    ! FACE 5
    u(1:nx,0,   5) =  u(1:nx,ny,   4)
    u(1:nx,ny+1,5) = -v(1,ny:1:-1, 1)
    u(nx+1,1:ny,5) =  u(1,1:ny,    6)

    v(1:nx,0,   5) =  v(1:nx,ny,   4)
    v(0,1:ny,   5) = -u(nx:1:-1,ny,3)
    v(nx+1,1:ny,5) =  v(1,1:ny,    6)

    ! FACE 6
    u(1:nx,0,   6) = -v(nx,ny:1:-1,4)
    u(1:nx,ny+1,6) =  u(1:nx,1,    1)
    u(0,1:ny,   6) =  u(nx,1:ny,   5)

    v(1:nx,ny+1,6) =  v(1:nx,1,    1)
    v(0,1:ny,   6) =  v(nx,1:ny,   5)
    v(nx+1,1:ny,6) = -u(nx:1:-1,1, 2)

! Deal with "flipped" edges
    if (D2A) then
       u(0,2:ny,   1) =  v(nx:2:-1,ny,5)
       v(2:nx,ny+1,1) =  u(1,ny:2:-1, 3)
       u(nx+1,2:ny,2) =  v(nx:2:-1,1, 4)
       v(2:nx,0,   2) =  u(nx,ny:2:-1,6)
       u(0,2:ny,   3) =  v(nx:2:-1,ny,1)
       v(2:nx,ny+1,3) =  u(1,ny:2:-1, 5)
       u(nx+1,2:ny,4) =  v(nx:2:-1,1, 6)
       v(2:nx,0,   4) =  u(nx,ny:2:-1,2)
       u(0,2:ny,   5) =  v(nx:2:-1,ny,3)
       v(2:nx,ny+1,5) =  u(1,ny:2:-1, 1)
       u(nx+1,2:ny,6) =  v(nx:2:-1,1, 2)
       v(2:nx,0,   6) =  u(nx,ny:2:-1,4)
       ! fill these edges (on first cell/row) around corner
       u(0,   1,   1) =  v(1,ny,6)
       v(1,   ny+1,1) = -v(1,ny,5)
       u(nx+1,   1,2) = -u(nx,1,6)
       v(1,   0,   2) =  u(nx,1,1)
       u(0,   1,   3) =  v(1,ny,2)
       v(1,   ny+1,3) = -v(1,ny,1)
       u(nx+1,1,   4) = -u(nx,1,2)
       v(1,   0,   4) =  u(nx,1,3)
       u(0,   1,   5) =  v(1,ny,4)
       v(1,   ny+1,5) = -v(1,ny,3)
       u(nx+1,1,   6) = -u(nx,1,4)
       v(1,   0,   6) =  u(nx,1,5)
    else
       u(0,1:ny,   1) =  v(nx:1:-1,ny,5)
       v(1:nx,ny+1,1) =  u(1,ny:1:-1, 3)
       u(nx+1,1:ny,2) =  v(nx:1:-1,1, 4)
       v(1:nx,0,   2) =  u(nx,ny:1:-1,6)
       u(0,1:ny,   3) =  v(nx:1:-1,ny,1)
       v(1:nx,ny+1,3) =  u(1,ny:1:-1, 5)
       u(nx+1,1:ny,4) =  v(nx:1:-1,1, 6)
       v(1:nx,0,   4) =  u(nx,ny:1:-1,2)
       u(0,1:ny,   5) =  v(nx:1:-1,ny,3)
       v(1:nx,ny+1,5) =  u(1,ny:1:-1, 1)
       u(nx+1,1:ny,6) =  v(nx:1:-1,1, 2)
       v(1:nx,0,   6) =  u(nx,ny:1:-1,4)
    end if

!  Set ghost corners to zero (to prevent "uninitialized" error). 
!  These values should not be used in any calculation

    u(0,   ny+1,:) = 0.0
    u(0,      0,:) = 0.0
    u(nx+1,   0,:) = 0.0
    u(nx+1,ny+1,:) = 0.0

    v(0,   ny+1,:) = 0.0
    v(0,      0,:) = 0.0
    v(nx+1,   0,:) = 0.0
    v(nx+1,ny+1,:) = 0.0

  end subroutine GhostCubeVector

  subroutine deallocGlob_i4_4(a, status)
     integer, pointer :: a(:,:,:,:)
     integer, intent(out) :: status

     if (associated(a)) then
        a=0
        call MAPL_DeAllocNodeArray(a, rc=status)
        if (status==MAPL_NoShm) deallocate(a, stat=status)
        if (status /= 0) return
        nullify(a)
     else
        status = 0
     end if

  end subroutine deallocGlob_i4_4

  subroutine deallocGlob_i4_3(a, status)
     integer, pointer :: a(:,:,:)
     integer, intent(out) :: status

     if (associated(a)) then
        a=0
        call MAPL_DeAllocNodeArray(a, rc=status)
        if (status==MAPL_NoShm) deallocate(a, stat=status)
        if (status /= 0) return
        nullify(a)
     else
        status = 0
     end if
     
  end subroutine deallocGlob_i4_3

  subroutine deallocGlob_i4_2(a, status)
     integer, pointer :: a(:,:)
     integer, intent(out) :: status
     
     if (associated(a)) then
        a=0
        call MAPL_DeAllocNodeArray(a, rc=status)
        if (status==MAPL_NoShm) deallocate(a, stat=status)
        if (status /= 0) return
        nullify(a)
     else
        status = 0
     end if
     
  end subroutine deallocGlob_i4_2

  subroutine deallocGlob_r8_3(a, status)
     real(REAL64), pointer :: a(:,:,:)
     integer, intent(out) :: status
     
     if (associated(a)) then
        a=0
        call MAPL_DeAllocNodeArray(a, rc=status)
        if (status==MAPL_NoShm) deallocate(a, stat=status)
        if (status /= 0) return
        nullify(a)
     else
        status = 0
     end if
                
  end subroutine deallocGlob_r8_3

  subroutine deallocGlob_r8_4(a, status)
     real(REAL64), pointer :: a(:,:,:,:)
     integer, intent(out) :: status

     if (associated(a)) then
        a=0
        call MAPL_DeAllocNodeArray(a, rc=status)
        if (status==MAPL_NoShm) deallocate(a, stat=status)
        if (status /= 0) return
        nullify(a)
     else
        status = 0
     end if

  end subroutine deallocGlob_r8_4

  subroutine deallocLocl_r8_3(a, status)
     real(REAL64), pointer :: a(:,:,:)
     integer, intent(out) :: status

     if (associated(a)) then
        deallocate(a)
        if (status /= 0) return
        nullify(a)
     else
        status = 0
     end if
  end subroutine deallocLocl_r8_3

end Module CubeLatLonTransformMod
