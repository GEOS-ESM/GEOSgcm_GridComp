module mapl_ComponentSpecification
   implicit none
   private

   public :: ComponentSpecifaciton

contains

   function new_ComponentSpecifaction(component)
      class(MaplComponent), intent(in) :: component

      if (component%is_valid()) then
         spec%component = component
         spec%state_specification = component%get_state_spec()
         spec%entry_points = component%get_entry_points()
         spec%child_specs = component%get_child_specs()
         spec%connecitons = component%get_connections()
      end if

   end function new_ComponentSpecifaction
end module mapl_ComponentSpecification
