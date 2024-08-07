module pymoist_interface_mod

   use iso_c_binding, only: c_int, c_float, c_double, c_bool
   
   implicit none

   private
   public :: pymoist_interface_f_init
   public :: pymoist_interface_f_run_AerActivation
   public :: pymoist_interface_f_finalize
   public :: make_moist_flags_C_interop
   public :: moist_flags_interface_type

   !-----------------------------------------------------------------------
   ! Shadow C interoperable config struct for FV. See `fv_arrays.f90` for
   ! the original structure, docs and default values
   !-----------------------------------------------------------------------
   type, bind(c) :: moist_flags_interface_type
      ! Grid information
      integer(kind=c_int) :: npx
      integer(kind=c_int) :: npy
      integer(kind=c_int) :: npz
      integer(kind=c_int) :: layout_x
      integer(kind=c_int) :: layout_y
      integer(kind=c_int) :: n_tiles
      ! Aer Activation
      integer(kind=c_int) :: n_modes
      ! Magic number
      integer(kind=c_int) :: make_flags_C_interop = 123456789
   end type


   interface

      subroutine pymoist_interface_f_init( moist_flags) bind(c, name='pymoist_interface_c_init')

         import c_int, c_float, c_double, moist_flags_interface_type

         implicit none
         type(moist_flags_interface_type), intent(in) :: moist_flags

      end subroutine pymoist_interface_f_init

      subroutine pymoist_interface_f_run_AerActivation( &
         aero_dgn, aero_num, aero_hygroscopicity, aero_sigma, &
         frland, nn_ocean, nn_land, &
         t, plo, &
         qicn, qils, qlcn, qlls, &
         vvel, tke, &
         nacti, nwfa, nact ) bind(c, name='pymoist_interface_c_run_AerActivation')

         import c_int, c_float, c_double

         implicit none

         ! Input
         real(kind=c_float), dimension(*), intent(in) :: aero_dgn
         real(kind=c_float), dimension(*), intent(in) :: aero_num
         real(kind=c_float), dimension(*), intent(in) :: aero_hygroscopicity
         real(kind=c_float), dimension(*), intent(in) :: aero_sigma

         real(kind=c_float), dimension(*), intent(in) :: frland
         real(kind=c_float), value, intent(in) :: nn_ocean
         real(kind=c_float), value, intent(in) :: nn_land

         real(kind=c_float), dimension(*), intent(in) :: t
         real(kind=c_float), dimension(*), intent(in) :: plo

         real(kind=c_float), dimension(*), intent(in) :: qicn
         real(kind=c_float), dimension(*), intent(in) :: qils
         real(kind=c_float), dimension(*), intent(in) :: qlcn
         real(kind=c_float), dimension(*), intent(in) :: qlls

         real(kind=c_float), dimension(*), intent(in) :: vvel
         real(kind=c_float), dimension(*), intent(in) :: tke

         ! Output
         real(kind=c_float), dimension(*), intent(in) :: nacti
         real(kind=c_float), dimension(*), intent(in) :: nwfa
         real(kind=c_float), dimension(*), intent(in) :: nact

      end subroutine pymoist_interface_f_run_AerActivation

      subroutine pymoist_interface_f_finalize() bind(c, name='pymoist_interface_c_finalize')
      end subroutine pymoist_interface_f_finalize

   end interface

contains

   subroutine make_moist_flags_C_interop(npx, npy, npz , nx, ny, n_tiles, n_modes, moist_flags)
      integer, intent(in) :: nx, ny, npx, npy, npz, n_tiles
      integer, intent(in) :: n_modes ! Aer Activation
      type(moist_flags_interface_type), intent(out) :: moist_flags

      moist_flags%npx = npx
      moist_flags%npy = npy
      moist_flags%npz = npz
      moist_flags%layout_x = nx
      moist_flags%layout_y = ny
      moist_flags%n_tiles = n_tiles
      moist_flags%n_modes = n_modes
   end subroutine make_moist_flags_C_interop

end module pymoist_interface_mod
