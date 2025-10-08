# CMakeLists.txt for Utils/Raster/preproc/topography
if(COMMAND esma_add_subdirectories)
  esma_add_subdirectories(@ncar_topo utils_topo)
else()
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/@ncar_topo/CMakeLists.txt")
    add_subdirectory(@ncar_topo)
  endif()
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/utils_topo/CMakeLists.txt")
    add_subdirectory(utils_topo)
  endif()
endif()

