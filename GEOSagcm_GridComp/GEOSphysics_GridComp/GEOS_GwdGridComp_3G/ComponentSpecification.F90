module mapl_ComponentSpecification
   implicit none
   private

   public :: ComponentSpecifaciton

contains

   function new_ComponentSpecifaction(component, config) result(spec)
      type(ComponentSpecification) :: spec
      class(MaplComponent), intent(in) :: component
      type(Configuration), intent(in) :: config

      if (component%is_valid()) then
         spec%component = component
         spec%state_specification = component%get_state_spec()
         spec%entry_points = component%get_entry_points()
         spec%child_specs = component%get_child_specs(config)
         spec%connections = component%get_connections(config)
      end if

   end function new_ComponentSpecifaction

end module mapl_ComponentSpecification
