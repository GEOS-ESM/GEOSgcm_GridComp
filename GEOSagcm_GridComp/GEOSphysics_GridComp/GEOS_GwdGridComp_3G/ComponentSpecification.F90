module mapl_ComponentSpecification
   implicit none
   private

   public :: ComponentSpecificiton
   type :: ComponentSpecificiton
      private
      class(MaplComponent), allocatable :: component
      type(StateSpecification), allocatable :: state_specs(:)
      type(EntryPoint), allocatable :: entry_points(:)
      type(ChildSpecification), allocatable :: child_specs(:)
      type(ConnectionSpecification), allocatable :: connection_specs(:)
      type(GridSpecification), allocatable :: grid_spec
   contains
      procedure :: get_component
      procedure :: get_state_specs
      procedure :: get_entry_points
      procedure :: get_child_specs
      procedure :: get_connection_Specs
      procedure :: get_grid_spec
   end type ComponentSpecificiton

contains

   function new_ComponentSpecifaction(component) result(spec)
      type(ComponentSpecification) :: spec
      class(MaplComponent), intent(in) :: component

      spec = ComponentSpecification( &
           component, &
           component%get_state_specs(), &
           component%get_entry_points(), &
           component%get_child_specs(), &
           component%get_connection_specs(), &
           component%get_grid_spec())

   end function new_ComponentSpecifaction

end module mapl_ComponentSpecification
