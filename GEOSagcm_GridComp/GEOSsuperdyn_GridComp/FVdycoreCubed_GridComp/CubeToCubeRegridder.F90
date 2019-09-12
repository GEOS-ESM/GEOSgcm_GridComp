#define _SUCCESS      0
#define _FAILURE     1
#define _VERIFY(A)   if(  A/=0) then; if(present(rc)) rc=A; PRINT *, Iam, __LINE__; return; endif
#define _ASSERT(A)   if(.not.A) then; if(present(rc)) rc=_FAILURE; PRINT *, Iam, __LINE__; return; endif
#define _RETURN(A)   if(present(rc)) rc=A; return

module CubeToCubeRegridderMod
   use MAPL_AbstractRegridderMod
   use CubeLatLonTransformMod
   use MAPL_GridSpecMod
   use MAPL_RegridderSpecMod
   use, intrinsic :: iso_fortran_env, only: REAL32
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none

   type, extends (AbstractRegridder) :: CubeToCubeRegridder
      private
      type (T_CubeCubeTransform), allocatable :: transform
   contains
      procedure :: initialize_subclass
      procedure :: regrid_scalar_2d_real32
      procedure :: regrid_scalar_3d_real32
      procedure :: regrid_vector_3d_real32
      procedure :: transpose_regrid_scalar_2d_real32
   end type CubeToCubeRegridder


   interface CubeToCubeRegridder
      module procedure newCubeToCubeRegridder
   end interface
   
   character(len=*), parameter :: MOD_NAME = 'MAPL_CubeToCubeRegridder::'
   integer, parameter :: NUM_DIMS = 2 
   
contains


   function newCubeToCubeRegridder(regridder_spec, rc) result(regridder)
      use ESMF
      type (CubeToCubeRegridder) :: regridder
      type (RegridderSpec), intent(in) :: regridder_spec
      integer, optional, intent(out) :: rc

      call regridder%initialize(regridder_spec)

   end function newCubeToCubeRegridder


   subroutine initialize_subclass(this, unusable, rc)
      use MAPL_KeywordEnforcerMod
      use MAPL_BaseMod
      use MAPL_CommsMod
      use ESMF
      class (CubeToCubeRegridder), intent(inout) :: this
      class (KeywordEnforcer), optional, intent(in) :: unusable
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME//'initialize'
      type (RegridderSpec) :: spec
      
      spec = this%get_spec()

      associate ( grid_in => spec%grid_in, grid_out => spec%grid_out )

        block
          integer :: dims(5)
          integer :: N_in(NUM_DIMS), N_out(NUM_DIMS)

          call MAPL_GridGet(grid_in, globalCellCountPerDim=dims, rc=status)
          _VERIFY(status)
          N_in = dims(1:2)
          call MAPL_GridGet(grid_out, globalCellCountPerDim=dims, rc=status)
          _VERIFY(status)
          N_out = dims(1:2)

          this%transform = CubeCubeCreate(N_in(1), N_in(2), N_out(1), N_out(2), rc=status)
          _VERIFY(status)

        end block

      end associate

      _RETURN(_SUCCESS)

   end subroutine initialize_subclass


   subroutine regrid_scalar_2d_real32(this, q_in, q_out, rc)
      class (CubeToCubeRegridder), intent(in) :: this
      real (kind=REAL32), intent(in) :: q_in(:,:)
      real (kind=REAL32), intent(out) :: q_out(:,:)
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME//'regrid_scalar_2d_real32'
      real (kind=REAL32), allocatable :: cs_tmp(:,:)

      cs_tmp = q_in

!!$      if (this%has_undef_value()) then
!!$         print*,'CubeToCube does not support misval'
!!$         _ASSERT(.false.)
!!$      else
         call CubeToCube(this%transform, cs_tmp, q_out, rc=status)
         _VERIFY(status)
!!$      end if

      deallocate(cs_tmp,stat=status)
      _VERIFY(status)

      _RETURN(_SUCCESS)
      
   end subroutine regrid_scalar_2d_real32


   subroutine regrid_scalar_3d_real32(this, q_in, q_out, rc)
      use MAPL_CommsMod
      use MAPL_BaseMod

      class (CubeToCubeRegridder), intent(in) :: this
      real (kind=REAL32), intent(in) :: q_in(:,:,:)
      real (kind=REAL32), intent(out) :: q_out(:,:,:)
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME//'regrid_scalar_2d_real32'
      integer :: k

      type (RegridderSpec) :: spec
      logical :: redistribute

      _ASSERT(size(q_in,3) == size(q_out,3))


      block
        integer :: N_in(NUM_DIMS)
        integer :: dims(5)
        
        spec = this%get_spec()
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

           call MAPL_GridGet(spec%grid_out, globalCellCountPerDim=dims, rc=status)
           _VERIFY(status)
           N_out = dims(1:2)

           allocate(q_out_global(n_out(1), n_out(2), size(q_in_global,3)))

           call MAPL_CollectiveGather3D(spec%grid_in, q_in, q_in_global, rc=status)
           _VERIFY(status)

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
      use, intrinsic :: iso_fortran_env, only: REAL32
      class (CubeToCubeRegridder), intent(in) :: this
      real(kind=REAL32), intent(in) :: u_in(:,:,:)
      real(kind=REAL32), intent(in) :: v_in(:,:,:)
      real(kind=REAL32), intent(out) :: u_out(:,:,:)
      real(kind=REAL32), intent(out) :: v_out(:,:,:)
      logical, optional, intent(in) :: rotate
      integer, optional, intent(out) :: rc
      character(len=*), parameter :: Iam = MOD_NAME//'regrid_vector_3d_real32'

      real(kind=REAL32) :: UVW_in (size(u_in,1),size(u_in,2),size(u_in,3)*3)
      real(kind=REAL32) :: UVW_out (size(u_out,1),size(u_out,2),size(u_out,3)*3)
      integer :: status

      real(kind=REAL32), allocatable :: u_in_restaggered(:,:,:)
      real(kind=REAL32), allocatable :: v_in_restaggered(:,:,:)
      associate (Trans => this%transform)

        u_in_restaggered = u_in
        v_in_restaggered = v_in
        call RestaggerWindsCube(u_in_restaggered, v_in_restaggered, D2A=.true.)

        call SphericalToCartesian(Trans, u_in_restaggered, v_in_restaggered, UVW_in)

        call this%regrid(UVW_in, UVW_out, rc=status)
        _VERIFY(status)

        call CartesianToSpherical(Trans, UVW_out, u_out, v_out)

        call RestaggerWindsCube(u_out, v_out, D2A=.false.)

      end associate

      _RETURN(_SUCCESS)
      
   end subroutine regrid_vector_3d_real32


   subroutine transpose_regrid_scalar_2d_real32(this, q_in, q_out, rc)
      class (CubeToCubeRegridder), intent(in) :: this
      real (kind=REAL32), intent(in) :: q_in(:,:)
      real (kind=REAL32), intent(out) :: q_out(:,:)
      integer, optional, intent(out) :: rc

      integer :: status
      character(len=*), parameter :: Iam = MOD_NAME//'regrid_scalar_2d_real32'

      print*,'CubeToCube does not support transpose operations.'
      _ASSERT(.false.)

   end subroutine transpose_regrid_scalar_2d_real32
   

end module CubeToCubeRegridderMod
