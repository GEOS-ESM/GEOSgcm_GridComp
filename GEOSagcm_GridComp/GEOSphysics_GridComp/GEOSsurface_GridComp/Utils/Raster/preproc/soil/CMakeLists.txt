# build without installation

add_library(m_loss_during_routines m_loss_during_routines.f90)
add_executable(loss_during_day.x loss_during_day.f90)
target_link_libraries(loss_during_day.x m_loss_during_routines)
add_executable(loss_surf_5cm_gensoil.x loss_surf_5cm_gensoil.f90)
target_link_libraries(loss_surf_5cm_gensoil.x m_loss_during_routines)
