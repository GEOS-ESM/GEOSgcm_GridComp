module pyMLINC_interface_mod

   use iso_c_binding, only: c_int, c_float, c_double, c_bool, c_ptr

   implicit none

   private
   public :: pyMLINC_interface_setservice_f, pyMLINC_interface_run_f
   public :: a_pod_struct_type

   !-----------------------------------------------------------------------
   ! See `interface.h` for explanation of the POD-strict struct
   !-----------------------------------------------------------------------
   type, bind(c) :: a_pod_struct_type
      integer(kind=c_int) :: npx
      integer(kind=c_int) :: npy
      integer(kind=c_int) :: npz
      ! Magic number
      integer(kind=c_int) :: make_flags_C_interop = 123456789
   end type


   interface

      subroutine pyMLINC_interface_setservice_f() bind(c, name='pyMLINC_interface_setservice_c')
      end subroutine pyMLINC_interface_setservice_f

      subroutine pyMLINC_interface_run_f(options, in_buffer, out_buffer) bind(c, name='pyMLINC_interface_run_c')

         import c_float, a_pod_struct_type

         implicit none
         ! This is an interface to a C function, the intent ARE NOT enforced
         ! by the compiler. Consider them developer hints
         type(a_pod_struct_type), intent(in) :: options
         real(kind=c_float), dimension(*), intent(in) :: in_buffer
         real(kind=c_float), dimension(*), intent(out) :: out_buffer

      end subroutine pyMLINC_interface_run_f
   
   end interface

end module pyMLINC_interface_mod
