module GEOS_ChemGridComp
   implicit none
   private


   public :: get_specification

   type, extends(MaplComponent) :: ChemGridComp
      private
      logical :: enable_PCHEM
      logical :: enable_ACHEM
      ...
   contains
      procedure :: get_state_specifications
      !      procedure :: get_children
      procedure :: get_connections
   end type ChemGridComp

contains

   function new_ChemGridComp(config) result(comp)
      type(GwdComponent) :: gwd
      type(Configuration), intent(in) :: config
      integer, optional, intent(out) :: rc

      integer :: status


      call config%get(comp%enable_PCHEM, "ENABLE_PCHEM", default=.false., _RC)
      call config%get(comp%enable_ACHEM, "ENABLE_ACHEM", default=.false., _RC)

      ...
      
      call comp%set_valid()

   end function new_ChemGridComp

   function get_specification(config) result(spec)
      type(ComponentSpecification) :: spec
      type(Configuration), intent(in) :: config ! yaml

      spec = ComponentSpecification(ChemGridComponent(config))

   end function get_specification

   ! Uses ACG to produce spec files to include.
   ! Entire function could become a macro.
   function make_state_specification(this) result(state_specifications)
      class(GhemGridComponent), intent(in) :: this
      type(StateSpecifications) :: state_specifications

!#include "import_specs.h"
#include "export_specs.h"
#include "internal_specs.h"

   end function make_state_specification


   function get_entry_points(this) result(entry_points)
      type(EntryPoint), allocatable :: entry_points(:)
      class(ChemGridComp), intent(in) :: this

      entry_points = [ &
           & EntryPoint(init, ESMF_METHOD_INITIALIZE), &
           & EntryPoint(run1, ESMF_METHOD_RUN, phase='1'), &
           & EntryPoint(run2, ESMF_METHOD_RUN, phase='2')
      ]

   end function get_entry_points

   function get_connections(this, config)


      ! handle the simple cases automatically
      connections = this%MaplComponent%get_connections(config)
      ! append the hardwired special cases
      connections = [connections, get_hardwired_connections()]

   contains

      function get_hardwired_connections()
         integer :: n
         ! hardwire the logic for any others
         Search_SO4v: do n = this%chemReg%i_SU, this%chemReg%j_SU
            if( trim(this%chemReg%vname(n)) == "SO4v" )  then
               call MAPL_AddConnectivity ( GC, &
                    SHORT_NAME  = (/ "GOCART::SO4v" /), &
                    DST_ID=GMICHEM, SRC_ID=GOCART, __RC__)
            end if
         end do Search_SO4v
      end function get_hardwired_connections

   end function get_connections
   


   
