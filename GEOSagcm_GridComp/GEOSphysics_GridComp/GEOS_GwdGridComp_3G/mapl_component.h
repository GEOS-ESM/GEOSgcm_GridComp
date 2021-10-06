! The following multiline macro creates a BIND(C) procedure that
! invokes the specification factory method of a gridded component.
! This supports DSO access.

#define DEFINE_COMPONENT(module_name, binding_name)	\
  function define_component(p_config, rc) result(p_spec)  bind(name=binding_name) ; \
      use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer ;\
      use module_name, only: specification;\
      use MAPL, only: ComponentSpecification ;\
      use yaFyaml, only: Configuration ;\
      type(c_ptr) :: p_spec ;\
      type(c_ptr), intent(in) :: p_config ;\
      integer, intent(in) :: rc ;\
      type(Configuration), pointer :: config ;\
      type(ComponentSpecification), pointer :: spec ;\
      integer :: status ;\
      call c_f_pointer(p_config, config) ;\
      call c_f_pointer(p_spec, spec) ;\
      spec = configure_component(config, rc=status) ;	\
   end function define_component
