esma_set_this()

set (srcs
  GEOS_HsGridComp.F90
)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL
  TYPE SHARED
)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

mapl_acg (
  ${this} HS_StateSpecs.rc
  IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS
)

add_subdirectory(regression)
