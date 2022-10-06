module pso_params_types_landshared

  implicit none
  private
  public :: pso_params_type_landshared

  type :: pso_params_type_landshared
    real, allocatable, dimension(:,:) :: param_vals
    integer                           :: ens_num
    real, allocatable, dimension(:)   :: map_vals
    integer                           :: total_ens
  end type pso_params_type_landshared
  
end module pso_params_types_landshared

module pso_params_mod_landshared
    use pso_params_types_landshared, only: pso_params_type_landshared

    implicit none
    private
    public :: pso_params

    type(pso_params_type_landshared) :: pso_params
end module pso_params_mod_landshared
