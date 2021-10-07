module MAPL_MaplComponent
   implicit none

   type, abstract :: MaplComponent
      private
   contains
      procedure, deferred :: fill_private_state
      procedure, deferred :: create_entry_points
      procedure, deferred :: create_state_specs
      procedure :: creat_grid_spec
      procedure :: create_child_specs
      procedure :: create_connection_specs
   end type MaplComponent

   type :: EntryPoint
      procedure(), nopass :: method
      type(ESMF_METHOD_FLAG) :: method_type
      character(:), allocatable :: phase
   contains
      procedure :: execute
   end type EntryPoint

   type, extends(MaplComponent) :: FallibleComponent
      logical :: is_valid
      class(MaplComponent), allocatable :: reference
   contains
      procedure :: fill_private_state
   end type FallibleComponent

contains

   function create_grid_spec(this) result(spec)
      type(GridSpec) :: spec
      class(MaplComponent), intent(in) :: this
      spec = GridSpec(MAPL_INHERIT_FROM_PARENT)
   end function create_grid_spec

   function create_child_specs(this) result(specs)
      type(ChildSpecification), allocatable :: specs(:)
      class(MaplComponent), intent(in) :: this
      specs = [ChildSpecification :: ]
   end function create_child_specs

   function create_connection_specs(this) result(specs)
      type(ConnectionSpecification) :: specs(:)
      specs = [ConnectionSpecification :: ]
   end function create_connection_specs

   function new_EntryPoint(proc, method_type, phase)
      procedure(
   end function new_EntryPoint

end module MAPL_MaplComponent
