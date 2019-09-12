#define _SUCCESS      0
#define _FAILURE     1
#define _VERIFY(A)   if(  A/=0) then; if(present(rc)) rc=A; PRINT *, Iam, __LINE__; return; endif
#define _ASSERT(A)   if(.not.A) then; if(present(rc)) rc=_FAILURE; PRINT *, Iam, __LINE__; return; endif
#define _RETURN(A)   if(present(rc)) rc=A; return

module CubeToLatLonRegridderMod
   use MAPL_AbstractRegridderMod
   use CubeLatLonTransformMod
   use MAPL_GridSpecMod
   use MAPL_RegridderSpecMod
   use, intrinsic :: iso_fortran_env, only: REAL32
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none

   type, extends (AbstractRegridder) :: CubeToLatLonRegridder
      private
      type (T_CubeLatLonTransform) :: transform
   contains
      procedure :: initialize_subclass
      procedure :: regrid_scalar_2d_real32
      procedure :: regrid_scalar_3d_real32
      procedure :: regrid_vector_3d_real32
      procedure :: transpose_regrid_scalar_2d_real32
      procedure :: transpose_regrid_scalar_3d_real32
      procedure :: transpose_regrid_vector_3d_real32
   end type CubeToLatLonRegridder


   interface CubeToLatLonRegridder
      module procedure newCubeToLatLonRegridder
   end interface
   
   character(len=*), parameter :: MOD_NAME = 'MAPL_CubeToLatLonRegridder::'
   integer, parameter :: NUM_DIMS = 2 

contains


   function newCubeToLatLonRegridder(regridder_spec, rc) result(regridder)
      type (CubeToLatLonRegridder) :: regridder
      type (RegridderSpec), intent(inout) :: regridder_spec
      integer, optional, intent(out) :: rc

      call regridder%initialize(regridder_spec)

   end function newCubeToLatLonRegridder

   subroutine initialize_subclass(this, unusable, rc)
      use MAPL_KeywordEnforcerMod
      use MAPL_BaseMod
      use MAPL_AbstractGridFactoryMod
      use MAPL_LatLonGridFactoryMod
      use MAPL_GridManagerMod
      class (CubeToLatLonRegridder), intent(inout) :: this
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc
         
      real(kind=REAL64), allocatable :: xll(:), yll(:)
      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME//'initialize_subclass'
      integer, parameter :: NUM_DIMS = 2
      type (RegridderSpec) :: spec

      integer :: N_in(NUM_DIMS), N_out(NUM_DIMS)
      integer :: dims(5)
      integer :: i0, i1, j0, j1
      integer :: local_ij(2,2)

      class (AbstractGridFactory), pointer :: factory
      
      spec = this%get_spec()

      associate ( grid_in => spec%grid_in, grid_out => spec%grid_out )

          call MAPL_GridGet(grid_in, globalCellCountPerDim=dims, rc=status)
          _VERIFY(status)
          N_in = dims(1:2)
          call MAPL_GridGet(grid_out, globalCellCountPerDim=dims, rc=status)
          _VERIFY(status)
          N_out = dims(1:2)

          ! Use factory to get global lat/lon 1d arrays
          
          factory => grid_manager%get_factory(grid_out)

          select type (latlon_factory => factory)
          type is (LatLonGridFactory)
             xll = latlon_factory%get_longitudes()
             yll = latlon_factory%get_latitudes()
          class default
             ! factory must be lat lon for regrid to lat-lon
             _ASSERT(.false.)
          end select

          call MAPL_grid_interior(grid_out, i0, i1, j0, j1)
          local_ij = reshape([i0,i1,j0,j1],[2,2])
          this%transform = CubeLatLonCreate(N_in(1), N_in(2), N_out(1), N_out(2), &
               & xll(:), yll(:), local_ij, rc=status)
          _VERIFY(status)
          call CubeLatLonSubset(this%transform,.true.)

      end associate

      _RETURN(_SUCCESS)
   end subroutine initialize_subclass

   subroutine regrid_scalar_2d_real32(this, q_in, q_out, rc)
      class (CubeToLatLonRegridder), intent(in) :: this
      real (kind=REAL32), intent(in) :: q_in(:,:)
      real (kind=REAL32), intent(out) :: q_out(:,:)
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME//'regrid_scalar_2d_real32'
      real (kind=REAL32) :: undef
      real (kind=REAL32), allocatable :: cs_tmp(:,:)


      cs_tmp = q_in ! need to copy because of mixed intent in CubeToLatLon().
      if (this%has_undef_value()) then
         undef = this%get_undef_value()
         call CubeToLatLon(this%transform, cs_tmp, q_out, transpose=.false., misval=undef, rc=status)
         _VERIFY(status)
      else
         call CubeToLatLon(this%transform, cs_tmp, q_out, transpose=.false., rc=status)
         _VERIFY(status)
      end if
      deallocate(cs_tmp,stat=status)
      _VERIFY(status)

      _RETURN(_SUCCESS)
      
   end subroutine regrid_scalar_2d_real32


   subroutine regrid_scalar_3d_real32(this, q_in, q_out, rc)
      use MAPL_CommsMod
      use MAPL_BaseMod

      class (CubeToLatLonRegridder), intent(in) :: this
      real (kind=REAL32), intent(in) :: q_in(:,:,:)
      real (kind=REAL32), intent(out) :: q_out(:,:,:)
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME//'regrid_scalar_2d_real32'
      integer :: k

      type (RegridderSpec) :: spec

      logical :: redistribute

      _ASSERT(size(q_in,3) == size(q_out,3))

      spec = this%get_spec()

      block
        integer :: N_in(NUM_DIMS)
        integer :: dims(5)
        
        call MAPL_GridGet(spec%grid_in, globalCellCountPerDim=dims, rc=status)
        _VERIFY(status)
        N_in = dims(1:2)

        if ((N_in(1) /= size(q_in,1)) .or. (N_in(2) /= size(q_in,2))) then
           redistribute = .true.
        else
           redistribute = .false.
        end if

      end block

      if (redistribute) then
         block
           real (kind=REAL32), pointer :: q_in_global(:,:,:)
           real (kind=REAL32), allocatable :: q_out_global(:,:,:)
           integer :: N_out(NUM_DIMS)
           integer :: dims(5)
        
           q_in_global=> null()

           call MAPL_CollectiveGather3D(spec%grid_in, q_in, q_in_global, rc=status)
           _VERIFY(status)

           call MAPL_GridGet(spec%grid_out, globalCellCountPerDim=dims, rc=status)
           _VERIFY(status)
           N_out = dims(1:2)

           allocate(q_out_global(n_out(1), n_out(2), size(q_in_global,3)))

           if (size(q_in_global) > 1) then
              do k = 1, size(q_in_global,3)
                 call this%regrid(q_in_global(:,:,k), q_out_global(:,:,k), rc=status)
                 _VERIFY(status)
              end do
           end if

           deallocate(q_in_global)

           call MAPL_CollectiveScatter3D(spec%grid_out, q_out_global, q_out, rc=status)
           _VERIFY(status)

           deallocate(q_out_global)
         end block
      else
         if (size(q_in) > 1) then
            do k = 1, size(q_in,3)
               call this%regrid(q_in(:,:,k), q_out(:,:,k), rc=status)
               _VERIFY(status)
            end do
         end if
      end if

      _RETURN(_SUCCESS)
      
   end subroutine regrid_scalar_3d_real32


   subroutine regrid_vector_3d_real32(this, u_in, v_in, u_out, v_out, rotate, rc)
      use MAPL_CommsMod
      use MAPL_BaseMod
      use, intrinsic :: iso_fortran_env, only: REAL32
      class (CubeToLatLonRegridder), intent(in) :: this
      real(kind=REAL32), intent(in) :: u_in(:,:,:)
      real(kind=REAL32), intent(in) :: v_in(:,:,:)
      real(kind=REAL32), intent(out) :: u_out(:,:,:)
      real(kind=REAL32), intent(out) :: v_out(:,:,:)
      logical, optional, intent(in ) :: rotate
      integer, optional, intent(out) :: rc
      character(len=*), parameter :: Iam = MOD_NAME//'regrid_vector_3d_real32'


      logical :: rotateBefore, rotateAfter
      integer :: status

      logical, parameter :: inputIsLL = .false.
      type (RegridderSpec) :: spec

      real (kind=REAL32), allocatable :: uvw_in(:,:,:)
      real (kind=REAL32), allocatable :: uvw_out(:,:,:)

      _ASSERT(size(u_in,3) == size(u_out,3))
      _ASSERT(size(v_in,3) == size(v_out,3))
      _ASSERT(size(u_in,3) == size(v_in,3))

      RotateBefore = InputIsLL
      RotateAfter  = .not.InputIsLL
      if (present(rotate)) then
         if (rotate) then
            RotateBefore = .true.
            RotateAfter  = .true.
         end if
      end if

      allocate(uvw_in(size(u_in,1),size(u_in,2),3*size(u_in,3)))
        
      call SphericalToCartesian(this%transform, u_in , v_in , UVW_in , &
           & transpose=.false., SphIsLL=InputIsLL,  Rotate=RotateBefore, rc=status)
      _VERIFY(status)
      
      allocate(uvw_out(size(u_out,1),size(u_out,2),3*size(u_out,3)))
      call this%regrid(uvw_in, uvw_out, rc=status)
      _VERIFY(status)
      deallocate(uvw_in)
      
      call CartesianToSpherical(this%transform, uvw_out, u_out, v_out, &
           transpose=.false., SphIsLL=.not.InputIsLL, Rotate=RotateAfter, rc=status)
      _VERIFY(status)
      deallocate(uvw_out)
      _RETURN(_SUCCESS)

   end subroutine regrid_vector_3d_real32


   subroutine transpose_regrid_scalar_2d_real32(this, q_in, q_out, rc)
      class (CubeToLatLonRegridder), intent(in) :: this
      real (kind=REAL32), intent(in) :: q_in(:,:)
      real (kind=REAL32), intent(out) :: q_out(:,:)
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME//'regrid_scalar_2d_real32'
      real (kind=REAL32) :: undef

      real (kind=REAL32), allocatable :: ll_tmp(:,:)

      ll_tmp = q_in ! need to copy because of mixed intent in CubeToLatLon().

      if (this%has_undef_value()) then
         undef = this%get_undef_value()
         call CubeToLatLon(this%transform, q_out, ll_tmp, transpose=.true., misval=undef, rc=status)
         _VERIFY(status)
      else
         call CubeToLatLon(this%transform, q_out, ll_tmp, transpose=.true., rc=status)
         _VERIFY(status)
      end if

      deallocate(ll_tmp,stat=status)
      _VERIFY(status)

      _RETURN(_SUCCESS)

   end subroutine transpose_regrid_scalar_2d_real32

   subroutine transpose_regrid_scalar_3d_real32(this, q_in, q_out, rc)
      use MAPL_CommsMod
      use MAPL_BaseMod

      class (CubeToLatLonRegridder), intent(in) :: this
      real (kind=REAL32), intent(in) :: q_in(:,:,:)
      real (kind=REAL32), intent(out) :: q_out(:,:,:)
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME//'transpose_regrid_scalar_2d_real32'
      integer :: k

      type (RegridderSpec) :: spec

      logical :: redistribute

      _ASSERT(size(q_in,3) == size(q_out,3))

      spec = this%get_spec()

      block
        integer :: N_in(NUM_DIMS)
        integer :: dims(5)
        
        call MAPL_GridGet(spec%grid_out, globalCellCountPerDim=dims, rc=status)
        _VERIFY(status)
        N_in = dims(1:2)

        if ((N_in(1) /= size(q_in,1)) .or. (N_in(2) /= size(q_in,2))) then
           redistribute = .true.
        else
           redistribute = .false.
        end if

      end block

      if (redistribute) then
         block
           real (kind=REAL32), pointer :: q_in_global(:,:,:)
           real (kind=REAL32), allocatable :: q_out_global(:,:,:)
           integer :: N_out(NUM_DIMS)
           integer :: dims(5)
        
           q_in_global=> null()

           call MAPL_CollectiveGather3D(spec%grid_out, q_in, q_in_global, rc=status)
           _VERIFY(status)

           call MAPL_GridGet(spec%grid_in, globalCellCountPerDim=dims, rc=status)
           _VERIFY(status)
           N_out = dims(1:2)

           allocate(q_out_global(n_out(1), n_out(2), size(q_in_global,3)))

           if (size(q_in_global) > 1) then
              do k = 1, size(q_in_global,3)
                 call this%transpose_regrid(q_in_global(:,:,k), q_out_global(:,:,k), rc=status)
                 _VERIFY(status)
              end do
           end if

           deallocate(q_in_global)

           call MAPL_CollectiveScatter3D(spec%grid_in, q_out_global, q_out, rc=status)
           _VERIFY(status)

           deallocate(q_out_global)
         end block
      else
         if (size(q_in) > 1) then
            do k = 1, size(q_in,3)
               call this%transpose_regrid(q_in(:,:,k), q_out(:,:,k), rc=status)
               _VERIFY(status)
            end do
         end if
      end if

      _RETURN(_SUCCESS)
      
   end subroutine transpose_regrid_scalar_3d_real32
   
   subroutine transpose_regrid_vector_3d_real32(this, u_in, v_in, u_out, v_out, rotate, rc)
      use MAPL_CommsMod
      use MAPL_BaseMod
      use, intrinsic :: iso_fortran_env, only: REAL32
      class (CubeToLatLonRegridder), intent(in) :: this
      real(kind=REAL32), intent(in) :: u_in(:,:,:)
      real(kind=REAL32), intent(in) :: v_in(:,:,:)
      real(kind=REAL32), intent(out) :: u_out(:,:,:)
      real(kind=REAL32), intent(out) :: v_out(:,:,:)
      logical, optional, intent(in ) :: rotate
      integer, optional, intent(out) :: rc
      character(len=*), parameter :: Iam = MOD_NAME//'transpose_regrid_vector_3d_real32'


      logical :: rotateBefore, rotateAfter
      integer :: status

      logical, parameter :: inputIsLL = .true.
      type (RegridderSpec) :: spec

      real (kind=REAL32), allocatable :: uvw_in(:,:,:)
      real (kind=REAL32), allocatable :: uvw_out(:,:,:)

      _ASSERT(size(u_in,3) == size(u_out,3))
      _ASSERT(size(v_in,3) == size(v_out,3))
      _ASSERT(size(u_in,3) == size(v_in,3))

      RotateBefore = InputIsLL
      RotateAfter  = .not.InputIsLL
      if (present(rotate)) then
         if (rotate) then
            RotateBefore = .true.
            RotateAfter  = .true.
         end if
      end if

      spec = this%get_spec()

      allocate(uvw_in(size(u_in,1),size(u_in,2),3*size(u_in,3)))
        
      call SphericalToCartesian(this%transform, u_in , v_in , UVW_in , &
           & transpose=.true., SphIsLL=InputIsLL,  Rotate=RotateBefore, rc=status)
      _VERIFY(status)
      
      allocate(uvw_out(size(u_out,1),size(u_out,2),3*size(u_out,3)))
      call this%transpose_regrid(uvw_in, uvw_out, rc=status)
      _VERIFY(status)
      deallocate(uvw_in)
      
      call CartesianToSpherical(this%transform, uvw_out, u_out, v_out, &
           transpose=.true., SphIsLL=.not.InputIsLL, Rotate=RotateAfter, rc=status)
      _VERIFY(status)
      deallocate(uvw_out)
      _RETURN(_SUCCESS)

   end subroutine transpose_regrid_vector_3d_real32

end module CubeToLatLonRegridderMod
