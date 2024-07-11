! a module that holds the information about pso_params_type
module pso_params_type_mod

  implicit none
  save
  private
  public :: pso_params_type

  ! create the type pso_params_type
  type :: pso_params_type
    ! param_vals: first dimension is each parameter, second dimension
    ! is each particle
    ! in current case, that is the 5 different g1's then the five
    ! different ksat parameters
    real, allocatable, dimension(:,:) :: param_vals
    ! particle_num: which particle number are we?
    integer                           :: particle_num
    real, allocatable, dimension(:)   :: map_vals
    real, allocatable, dimension(:)   :: lai_vals
    real, allocatable, dimension(:)   :: sand_vals
    real, allocatable, dimension(:)   :: k0_vals
    real, allocatable, dimension(:)   :: canopy_vals
    integer, allocatable, dimension(:)   :: local_tile_nums
    integer, allocatable, dimension(:)   :: all_tile_nums
    integer                           :: total_ens
  end type pso_params_type
  
end module pso_params_type_mod

module pso_params
    use pso_params_type_mod, only: pso_params_type

    implicit none
    save
    private
    public :: pso_vals

    type(pso_params_type) :: pso_vals
end module pso_params
