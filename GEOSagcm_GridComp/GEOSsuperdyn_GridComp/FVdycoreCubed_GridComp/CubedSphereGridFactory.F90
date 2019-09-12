!-----------------------------------------------------
! Note that this implementation only supports
! "square" faces on the cube.   I.e. the number of
! cells along each axis (of each face) are the same.
! IM_WORLD is used for this quantity, and there is no
! equivalent for the "other" axis.
!-----------------------------------------------------

#define _SUCCESS      0
#define _FAILURE     1
#define _VERIFY(A)   if(  A/=0) then; if(present(rc)) rc=A; PRINT *, Iam, __LINE__; return; endif
#define _ASSERT(A)   if(.not.(A)) then; if(present(rc)) rc=_FAILURE; PRINT *, Iam, __LINE__; return; endif
#define _RETURN(A)   if(present(rc)) rc=A; return
#include "unused_dummy.H"


module CubedSphereGridFactoryMod
   use fv_grid_utils_mod, only: gnomonic_grids, cell_center2
   use fv_grid_tools_mod, only: mirror_grid
   use MAPL_AbstractGridFactoryMod
   use MAPL_MinMaxMod
   use MAPL_KeywordEnforcerMod
   use ESMF
   use MAPL_CommsMod
   use MAPL_IOMod, only : GETFILE, FREE_FILE 
   use fms_mod, only: fms_init
   use mpp_domains_mod,   only: domain2d, mpp_update_domains
   use mpp_mod,           only: mpp_error, FATAL, NOTE
   !use CubedSphereA2D2CMod
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
   private

   public :: CubedSphereGridFactory

   integer, parameter :: ndims = 2
   integer, parameter :: UNDEFINED_INTEGER = 1-huge(1)
   real, parameter :: UNDEFINED_REAL = huge(1.)
   character(len=*), parameter :: UNDEFINED_CHAR = '**'

   integer, parameter :: FV_GRID_TYPE_DEFAULT = 0
   character(len=*), parameter :: GRID_NAME_DEFAULT = 'UNKNOWN'

   integer, parameter :: NUM_CUBE_FACES = 6

   interface
      subroutine A2D2C(U,V,npz,getC)
      implicit none
      real,    intent(INOUT)           :: U(:,:,:)
      real,    intent(INOUT)           :: V(:,:,:)
      integer, intent(   IN)           :: npz
      logical, optional, intent(   IN)           :: getC
      end subroutine A2D2C
   end interface

   type, extends(AbstractGridFactory) :: CubedSphereGridFactory
      private

      
      character(len=:), allocatable :: grid_name
      integer :: grid_type = UNDEFINED_INTEGER

      ! Grid dimensions - Note that we only support "square" grids
      integer :: im_world = UNDEFINED_INTEGER
      integer :: lm = UNDEFINED_INTEGER
      integer :: ntiles = NUM_CUBE_FACES

      ! Domain decomposition: - note that we only support "square" dec
      integer :: nx = UNDEFINED_INTEGER
      integer :: ny = UNDEFINED_INTEGER
      integer, allocatable :: ims(:)
      integer, allocatable :: jms(:)
      ! rectangle decomposition
      integer, allocatable :: jms_2d(:,:)

      ! For halo
      type (domain2d) :: domain

      logical :: halo_initialized = .false.

   contains

      procedure :: make_new_grid
      procedure :: create_basic_grid

      procedure :: initialize_from_config_with_prefix
      procedure :: initialize_from_esmf_distGrid

      procedure :: halo_init
      procedure :: halo

      procedure :: check_and_fill_consistency
      procedure :: equals
      procedure :: generate_grid_name
      procedure :: to_string

      procedure :: a2d2c => a2d2c_factory
   end type CubedSphereGridFactory
   
   character(len=*), parameter :: MOD_NAME = 'CubedSphereGridFactory::'
   
   interface CubedSphereGridFactory
      module procedure CubedSphereGridFactory_from_parameters
   end interface CubedSphereGridFactory

   interface set_with_default
      module procedure set_with_default_integer
      module procedure set_with_default_real
      module procedure set_with_default_character
      module procedure set_with_default_bounds
   end interface set_with_default


contains


   function CubedSphereGridFactory_from_parameters(unusable, grid_name, grid_type, &
        & im_world, lm, nx, ny, ims, jms, &
        & rc) result(factory)
      type (CubedSphereGridFactory) :: factory
      class (KeywordEnforcer), optional, intent(in) :: unusable
      character(len=*), optional, intent(in) :: grid_name
      integer, optional, intent(in) :: grid_type

      ! grid details:
      integer, optional, intent(in) :: im_world
      integer, optional, intent(in) :: lm

      ! decomposition:
      integer, optional, intent(in) :: nx
      integer, optional, intent(in) :: ny
      integer, optional, intent(in) :: ims(:)
      integer, optional, intent(in) :: jms(:)

      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME // 'CubedSphereGridFactory_from_parameters'

      if (present(unusable)) print*,shape(unusable)

      call set_with_default(factory%grid_name, grid_name, GRID_NAME_DEFAULT)
      call set_with_default(factory%grid_type, grid_type, FV_GRID_TYPE_DEFAULT)

      call set_with_default(factory%nx, nx, UNDEFINED_INTEGER)
      call set_with_default(factory%ny, ny, UNDEFINED_INTEGER)

      call set_with_default(factory%im_world, im_world, UNDEFINED_INTEGER)
      call set_with_default(factory%lm, lm, UNDEFINED_INTEGER)

      ! default is unallocated
      if (present(ims)) factory%ims = ims
      if (present(jms)) factory%jms = jms

      call factory%check_and_fill_consistency(rc=status)

      _VERIFY(status)

      _RETURN(_SUCCESS)

   end function CubedSphereGridFactory_from_parameters


   function make_new_grid(this, unusable, rc) result(grid)
      type (ESMF_Grid) :: grid
      class (CubedSphereGridFactory), intent(in) :: this
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME // 'make_grid'

      grid = this%create_basic_grid(rc=status)
      _VERIFY(status)

      _RETURN(_SUCCESS)

   end function make_new_grid


   
   function create_basic_grid(this, unusable, rc) result(grid)
      type (ESMF_Grid) :: grid
      class (CubedSphereGridFactory), intent(in) :: this
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      integer :: i,nTile
      integer, allocatable :: ims(:,:), jms(:,:)
      real(kind=ESMF_KIND_R8), pointer :: centers(:,:)
      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME // 'create_basic_grid'

      if (this%grid_type <=3) then
         nTile=6
      else
         nTile=1
      end if

      allocate(ims(this%nx,nTile))
      do i=1,nTile
         ims(:,i)=this%ims
      enddo

      if(allocated(this%jms_2d)) then
         _ASSERT(size(this%jms_2d,2) == 6) 
         allocate(jms, source = this%jms_2d)
      else
         allocate(jms(this%ny,nTile))
         do i=1,nTile
            jms(:,i)=this%jms
         end do
      endif

      if (this%grid_type <= 3) then
         grid = ESMF_GridCreateCubedSPhere(this%im_world,countsPerDEDim1PTile=ims, &
                   countsPerDEDim2PTile=jms ,name=this%grid_name, &
                   staggerLocList=[ESMF_STAGGERLOC_CENTER,ESMF_STAGGERLOC_CORNER], coordSys=ESMF_COORDSYS_SPH_RAD, rc=status)
         _VERIFY(status)
         call ESMF_AttributeSet(grid, 'GridType', 'Cubed-Sphere', rc=status)
      else
         grid = ESMF_GridCreateNoPeriDim( &
              & name = this%grid_name, &
              & countsPerDEDim1=this%ims, &
              & countsPerDEDim2=this%jms, &
              & indexFlag=ESMF_INDEX_DELOCAL, &
              & gridEdgeLWidth=[0,0], &
              & gridEdgeUWidth=[1,1], &
              & coordDep1=[1,2], &
              & coordDep2=[1,2], &
              & coordSys=ESMF_COORDSYS_SPH_RAD, &
              & rc=status)
         _VERIFY(status)
         call ESMF_AttributeSet(grid, 'GridType', 'Doubly-Periodic', rc=status)
         _VERIFY(status)
         call ESMF_GridAddCoord(grid,rc=status)
         _VERIFY(status)
         call ESMF_GridGetCoord(grid, coordDim=1, localDE=0, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           farrayPtr=centers, rc=status)
         _VERIFY(status)
         centers=0.0
         call ESMF_GridGetCoord(grid, coordDim=2, localDE=0, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           farrayPtr=centers, rc=status)
         _VERIFY(status)
         centers=0.0

      end if

      deallocate(ims,jms)

      if (this%lm /= UNDEFINED_INTEGER) then
         call ESMF_AttributeSet(grid, name='GRID_LM', value=this%lm, rc=status)
         _VERIFY(status)
      end if

      call ESMF_AttributeSet(grid, name='NEW_CUBE', value=1,rc=status)
      _VERIFY(status)

      _RETURN(_SUCCESS)
   end function create_basic_grid
   
   subroutine initialize_from_config_with_prefix(this, config, prefix, unusable, rc)
      use esmf
      class (CubedSphereGridFactory), intent(inout) :: this
      type (ESMF_Config), intent(inout) :: config
      character(len=*), intent(in) :: prefix  ! effectively optional due to overload without this argument
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME//'make_geos_grid_from_config'
      character(len=ESMF_MAXSTR) :: tmp
      type (ESMF_VM) :: vm
      integer :: vmcomm, ndes

      if (present(unusable)) print*,shape(unusable)

      call ESMF_VMGetCurrent(vm, rc=status)
      _VERIFY(status)

      call ESMF_ConfigGetAttribute(config, tmp, label=prefix//'GRIDNAME:', default=GRID_NAME_DEFAULT)
      this%grid_name = trim(tmp)

      call ESMF_ConfigGetAttribute(config, this%grid_type, label=prefix//'CS_GRID_TYPE:', default=FV_GRID_TYPE_DEFAULT)

      call ESMF_ConfigGetAttribute(config, this%nx, label=prefix//'NX:', default=UNDEFINED_INTEGER)
      call ESMF_ConfigGetAttribute(config, this%ny, label=prefix//'NY:', default=UNDEFINED_INTEGER)

      call ESMF_ConfigGetAttribute(config, this%im_world, label=prefix//'IM_WORLD:', default=UNDEFINED_INTEGER)

      call get_multi_integer(this%ims, 'IMS:', rc=status)
      _VERIFY(status)

      call ESMF_ConfigGetAttribute(config, tmp, label=prefix//'JMS_FILE:', rc=status)
      if (status == _SUCCESS) then
         call get_jms_from_file(this%jms_2d, trim(tmp),this%ny, rc=status)
         _VERIFY(status)
      else
         call get_multi_integer(this%jms, 'JMS:', rc=status)
         _VERIFY(status)
      endif

      call ESMF_ConfigGetAttribute(config, this%lm, label=prefix//'LM:', default=UNDEFINED_INTEGER)

      call this%check_and_fill_consistency(rc=status)
      _VERIFY(status)

      ! halo initialization
        
      call ESMF_VmGet(VM, mpicommunicator=vmcomm, petCount=ndes, rc=status)
      _VERIFY(status)


      _RETURN(_SUCCESS)

   contains
      
      subroutine get_multi_integer(values, label, rc)
         integer, allocatable, intent(out) :: values(:)
         character(len=*) :: label
         integer, optional, intent(out) :: rc

         integer :: i
         integer :: n
         integer :: tmp
         integer :: status
         
         call ESMF_ConfigFindLabel(config, label=prefix//label,rc=status)
         if (status /= _SUCCESS) then
            _RETURN(_SUCCESS)
         end if

         ! First pass:  count values
         n = 0
         do
            call ESMF_ConfigGetAttribute(config, tmp, default=UNDEFINED_INTEGER, rc=status)
            if (status /= _SUCCESS) then
               exit
            else
               n = n + 1
            end if
         end do

         ! Second pass: allocate and fill
         allocate(values(n), stat=status) ! no point in checking status
         _VERIFY(status)
         call ESMF_ConfigFindLabel(config, label=prefix//label,rc=status)
         do i = 1, n
            call ESMF_ConfigGetAttribute(config, values(i), rc=status)
            _VERIFY(status)
         end do

         _RETURN(_SUCCESS)

      end subroutine get_multi_integer

      subroutine get_jms_from_file(values, file_name, n, rc)
         integer, allocatable, intent(out) :: values(:,:)
         character(len=*), intent(in) :: file_name
         integer, intent(in) :: n
         integer, optional, intent(out) :: rc

         logical :: FileExists
         integer :: i,k,face,total, unit, max_procs
         integer :: status, N_proc,NF
         integer, allocatable :: values_tmp(:), values_(:,:)

    
         N_proc = n*6 ! it has been devided by 6. get back the original NY
         allocate(values_tmp(N_proc), stat=status) ! no point in checking status
         _VERIFY(status)

         inquire(FILE = trim(file_name), EXIST=FileExists)
         if ( .not. FileExists) then
            print*, file_name //  " does not exist"
             _RETURN(_FAILURE)

         elseif (MAPL_AM_I_Root(VM)) then

            UNIT = GETFILE ( trim(file_name), form="formatted", rc=status )
            _VERIFY(STATUS)
            read(UNIT,*) total, max_procs
            if (total /= N_proc) then
                print*, "n /= total"
                _RETURN(_FAILURE)
            endif
            do i = 1,total
                read(UNIT,*) values_tmp(i)
            enddo
            call FREE_FILE(UNIT)
         endif

         call MAPL_CommsBcast(VM, max_procs,  n=1, ROOT=MAPL_Root, rc=status)
         call MAPL_CommsBcast(VM, values_tmp, n=N_proc, ROOT=MAPL_Root, rc=status)
         _VERIFY(STATUS)

         ! distributed to 6 faces
         allocate(values_(max_procs,6))
         values_ = 0
         k = 1
         do NF = 1, 6
            face = 0
            do i = 1, max_procs
               values_(i,NF) = values_tmp(k)
               face = face + values_tmp(k)
               k = k+1
               if (face == this%im_world) exit
            enddo             
          enddo
          values = values_

         _RETURN(_SUCCESS)

      end subroutine get_jms_from_file

      subroutine get_bounds(bounds, label, rc)
         type(RealMinMax), intent(out) :: bounds
         character(len=*) :: label
         integer, optional, intent(out) :: rc

         integer :: i
         integer :: n
         integer :: status
         
         call ESMF_ConfigFindLabel(config, label=prefix//label,rc=status)
         if (status /= _SUCCESS) then
            _RETURN(_SUCCESS)
         end if

         ! Must be 2 values: min and max
         call ESMF_ConfigGetAttribute(config, bounds%min, rc=status)
         _VERIFY(status)
         call ESMF_ConfigGetAttribute(config, bounds%max, rc=status)
         _VERIFY(status)

         _RETURN(_SUCCESS)

      end subroutine get_bounds

      
   end subroutine initialize_from_config_with_prefix
   
   subroutine halo_init(this, halo_width)
      use CubeHaloMod, only: mpp_domain_decomp
      class (CubedSphereGridFactory), intent(inout) :: this
      integer, optional, intent(in) :: halo_width
      integer        :: npx, npy
      integer        :: ng, nregions, grid_type

      nregions = 6
      npx = this%im_world + 1
      npy = npx
      ng = 1        ! halowidth
      if(present(halo_width)) ng = halo_width
      grid_type = 0 ! cubed-sphere

      call mpp_domain_decomp(this%domain, npx, npy, nregions, ng, grid_type, &
           & this%nx*this%ny, this%nx, this%ny)
      
   end subroutine halo_init

   function to_string(this) result(string)
      character(len=:), allocatable :: string
      class (CubedSphereGridFactory), intent(in) :: this

      _UNUSED_DUMMY(this)
      string = 'CubedSphereGridFactory'

   end function to_string



   subroutine check_and_fill_consistency(this, unusable, rc)
      use MAPL_BaseMod, only: MAPL_DecomposeDim
      class (CubedSphereGridFactory), intent(inout) :: this
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME // 'check_and_fill_consistency'


      if (.not. allocated(this%grid_name)) then
         this%grid_name = GRID_NAME_DEFAULT
      end if

      if (this%grid_type == UNDEFINED_INTEGER) then
         this%grid_type = FV_GRID_TYPE_DEFAULT ! fv default
      end if

      ! Check decomposition/bounds
      ! WY notes: not necessary for this assert
      !_ASSERT(allocated(this%ims) .eqv. allocated(this%jms))
      call verify(this%nx, this%im_world, this%ims, rc=status)
      if (allocated(this%jms_2d)) then
        _ASSERT(size(this%jms_2d,2)==6) 
        _ASSERT(sum(this%jms_2d) == 6*this%im_world)
      else
         call verify(this%ny, this%im_world, this%jms, rc=status)
      endif
      
      _RETURN(_SUCCESS)
         
   contains

      subroutine verify(n, m_world, ms, rc)
         integer, intent(inout) :: n
         integer, intent(inout) :: m_world
         integer, allocatable, intent(inout) :: ms(:)
         integer, optional, intent(out) :: rc

         integer :: status

         if (allocated(ms)) then
            _ASSERT(size(ms) > 0)

            if (n == UNDEFINED_INTEGER) then
               n = size(ms)
            else
               _ASSERT(n == size(ms))
            end if

            if (m_world == UNDEFINED_INTEGER) then
               m_world = sum(ms)
            else
               _ASSERT(m_world == sum(ms))
            end if

         else

            _ASSERT(n /= UNDEFINED_INTEGER)
            _ASSERT(m_world /= UNDEFINED_INTEGER)
            allocate(ms(n), stat=status)
            _VERIFY(status)

            call MAPL_DecomposeDim (m_world, ms, n, symmetric=.true.)

         end if

         _RETURN(_SUCCESS)

      end subroutine verify

   end subroutine check_and_fill_consistency


   elemental subroutine set_with_default_integer(to, from, default)
      integer, intent(out) :: to
      integer, optional, intent(in) :: from
      integer, intent(in) :: default
      
      if (present(from)) then
         to = from
      else
         to = default
      end if
      
   end subroutine set_with_default_integer
   
   
   elemental subroutine set_with_default_real(to, from, default)
      real, intent(out) :: to
      real, optional, intent(in) :: from
      real, intent(in) :: default
      
      if (present(from)) then
         to = from
      else
         to = default
      end if
      
   end subroutine set_with_default_real
   
   subroutine set_with_default_character(to, from, default)
      character(len=:), allocatable, intent(out) :: to
      character(len=*), optional, intent(in) :: from
      character(len=*), intent(in) :: default
      
      if (present(from)) then
         to = from
      else
         to = default
      end if
      
   end subroutine set_with_default_character


   elemental subroutine set_with_default_bounds(to, from, default)
      type (RealMinMax), intent(out) :: to
      type (RealMinMax), optional, intent(in) :: from
      type (RealMinMax), intent(in) :: default
      
      if (present(from)) then
         to = from
      else
         to = default
      end if
      
   end subroutine set_with_default_bounds
   

   logical function equals(a, b)
      class (CubedSphereGridFactory), intent(in) :: a
      class (AbstractGridFactory), intent(in) :: b

      select type (b)
      class default
         equals = .false.
         return
      class is (CubedSphereGridFactory)
         equals = .true.

         equals = (a%im_world == b%im_world)
         if (.not. equals) return
         
         equals = (a%lm == b%lm)
         if (.not. equals) return
         
         ! same decomposition
         equals = all(a%ims == b%ims) 
         if (.not. equals) return

         if ( allocated(a%jms) .and. allocated(b%jms)) then
            equals = all(a%jms == b%jms)
            if (.not. equals) return
         else
            equals = all(a%jms_2d == b%jms_2d)
            if (.not. equals) return
         endif
         
      end select
         
   end function equals

   subroutine initialize_from_esmf_distGrid(this, dist_grid, lon_array, lat_array, unusable, rc)
      class (CubedSphereGridFactory), intent(inout)  :: this
      type (ESMF_DistGrid), intent(in) :: dist_grid
      type (ESMF_LocalArray), intent(in) :: lon_array
      type (ESMF_LocalArray), intent(in) :: lat_array
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME // 'CubedSphereGridFactory_initialize_from_esmf_distGrid'

      ! not implemented
      _ASSERT(.false.)

   end subroutine initialize_from_esmf_distGrid

   subroutine halo(this, array, unusable, halo_width, rc)
      use, intrinsic :: iso_fortran_env, only: REAL32
      class (CubedSphereGridFactory), intent(inout) :: this
      real(kind=REAL32), intent(inout) :: array(:,:)
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(in) :: halo_width
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME // 'halo'

      if (.not. this%halo_initialized) then
         call this%halo_init(halo_width = halo_width)
         this%halo_initialized = .true.
      end if

      call mpp_update_domains(array, this%domain)

      _RETURN(_SUCCESS)

   end subroutine halo

   function generate_grid_name(this) result(name)
      class (CubedSphereGridFactory), intent(in) :: this
      character(len=:), allocatable :: name

      character(len=4) :: im_string

      write(im_string,'(i4.4)') this%im_world

      name = 'CF' // im_string //'x6C'

   end function generate_grid_name

   subroutine a2d2c_factory(this, u, v, lm, getC, unusable, rc)
      class (CubedSphereGridFactory), intent(in) :: this
      real, intent(inout) :: u(:,:,:)
      real, intent(inout) :: v(:,:,:)
      integer, intent(in) :: lm
      logical, intent(in) :: getC
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      character(len=*), parameter :: Iam= MOD_NAME // 'A2D2C'
      integer :: status

      ! From GetWeights
      call a2d2c(u,v,lm,getC)

      _RETURN(_SUCCESS)

   end subroutine A2D2C_FACTORY

end module CubedSphereGridFactoryMod
