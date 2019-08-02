module CubeGridPrototype
  implicit none
  
  private
  public :: register_grid_and_regridders
contains

  subroutine register_grid_and_regridders()
    use MAPL_RegridderManagerMod, only: regridder_manager
    use MAPL_RegridderSpecMod, only: REGRID_METHOD_BILINEAR
    use LatLonToCubeRegridderMod
    use CubeToLatLonRegridderMod
    use CubeToCubeRegridderMod

    type (CubeToLatLonRegridder) :: cube_to_latlon_prototype
    type (LatLonToCubeRegridder) :: latlon_to_cube_prototype
    type (CubeToCubeRegridder) :: cube_to_cube_prototype

    associate (method => REGRID_METHOD_BILINEAR, mgr => regridder_manager)
      call mgr%add_prototype('Cubed-Sphere', 'LatLon', method, cube_to_latlon_prototype)
      call mgr%add_prototype('LatLon', 'Cubed-Sphere', method, latlon_to_cube_prototype)
      call mgr%add_prototype('Cubed-Sphere', 'Cubed-Sphere', method, cube_to_cube_prototype)
    end associate

  end subroutine register_grid_and_regridders
  
end module CubeGridPrototype
