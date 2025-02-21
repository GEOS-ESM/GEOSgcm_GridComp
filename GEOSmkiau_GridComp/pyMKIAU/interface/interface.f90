module pyMKIAU_interface_mod

   use iso_c_binding, only: c_int, c_float, c_double, c_bool, c_ptr

   implicit none

   private
   public :: pyMKIAU_interface_f_setservice, pyMKIAU_interface_f_run
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

      subroutine pyMKIAU_interface_f_setservice() bind(c, name='pyMKIAU_interface_c_setservice')
      end subroutine pyMKIAU_interface_f_setservice

      subroutine pyMKIAU_interface_f_run(options, in_buffer, out_buffer) bind(c, name='pyMKIAU_interface_c_run')

         import c_float, a_pod_struct_type

         implicit none
         ! This is an interface to a C function, the intent ARE NOT enforced
         ! by the compiler. Consider them developer hints
         type(a_pod_struct_type), intent(in) :: options
         real(kind=c_float), dimension(*), intent(in) :: in_buffer
         real(kind=c_float), dimension(*), intent(out) :: out_buffer

      end subroutine pyMKIAU_interface_f_run
   
   end interface

end module pyMKIAU_interface_mod
