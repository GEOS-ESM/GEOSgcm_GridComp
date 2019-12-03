set (this RRTMG_mods)

set (srcs
  parrrtm.F90 rrlw_cld.F90 rrlw_con.F90 rrlw_kg01.F90
  rrlw_kg02.F90 rrlw_kg03.F90 rrlw_kg04.F90 rrlw_kg05.F90 rrlw_kg06.F90
  rrlw_kg07.F90 rrlw_kg08.F90 rrlw_kg09.F90 rrlw_kg10.F90 rrlw_kg11.F90
  rrlw_kg12.F90 rrlw_kg13.F90 rrlw_kg14.F90 rrlw_kg15.F90 rrlw_kg16.F90
  rrlw_ncpar.F90 rrlw_ref.F90 rrlw_tbl.F90 rrlw_vsn.F90 rrlw_wvn.F90
  dbmemory.F90
  )

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES GEOS_Shared MAPL)

file (COPY _gpudef.inc DESTINATION ${esma_include}/RRTMG)
file (COPY _memdef.inc DESTINATION ${esma_include}/RRTMG)

