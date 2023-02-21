esma_set_this(OVERRIDE preprocess_soil_bcs)

esma_add_library (${this}
  SRCS m_loss_during_routines.f90 )


ecbuild_add_executable (
  TARGET loss_during_day.x
  SOURCES loss_during_day.f90
  LIBS ${this} OpenMP::OpenMP_Fortran)

ecbuild_add_executable (
  TARGET loss_surf_5cm_gensoil.x
  SOURCES loss_surf_5cm_gensoil.f90
  LIBS ${this} OpenMP::OpenMP_Fortran)
