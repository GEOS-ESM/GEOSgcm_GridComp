esma_set_this()

message("-- Build MAPLpyish")

set (srcs ./MAPLpyish/c_bridge_to_MAPL.F90)

esma_add_library (${this}
  TYPE SHARED
  SRCS ${srcs}
  DEPENDENCIES MAPL esmf)
