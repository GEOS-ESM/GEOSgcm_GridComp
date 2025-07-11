module pyMLINC_interface_mod

   use iso_c_binding, only: c_int, c_float, c_double, c_bool, c_ptr

   implicit none

   private

   public :: pyMLINC_interface_init_f, pyMLINC_interface_run_f

   interface

      subroutine pyMLINC_interface_init_f(magic_number) bind(c, name="pyMLINC_interface_init_c")
         import c_int
         implicit none
         integer(kind=c_int), value, intent(in) :: magic_number
      end subroutine pyMLINC_interface_init_f

      subroutine pyMLINC_interface_run_f( &
           ! Input
           xdim, ydim, zdim, &
           u, v, t, &
           qv, ql, qi, qr, qs, qg, &
           ps, &
           ! Output
           dtdt, &
           ! LAST ARGUMENT - input
           magic_number) bind(c, name="pyMLINC_interface_run_c")
         import c_int, c_float
         implicit none
         ! This is an interface to a C function, the intent is NOT enforced
         ! by the compiler. Consider them developer hints
         integer(kind=c_int), value, intent(in) :: xdim, ydim, zdim
         real(kind=c_float), dimension(*), intent(in) :: u
         real(kind=c_float), dimension(*), intent(in) :: v
         real(kind=c_float), dimension(*), intent(in) :: t
         real(kind=c_float), dimension(*), intent(in) :: qv
         real(kind=c_float), dimension(*), intent(in) :: ql
         real(kind=c_float), dimension(*), intent(in) :: qi
         real(kind=c_float), dimension(*), intent(in) :: qr
         real(kind=c_float), dimension(*), intent(in) :: qs
         real(kind=c_float), dimension(*), intent(in) :: qg
         real(kind=c_float), dimension(*), intent(in) :: ps
         real(kind=c_float), dimension(*), intent(out) :: dtdt
         integer(kind=c_int), value, intent(in) :: magic_number
      end subroutine pyMLINC_interface_run_f
   
   end interface

end module pyMLINC_interface_mod
