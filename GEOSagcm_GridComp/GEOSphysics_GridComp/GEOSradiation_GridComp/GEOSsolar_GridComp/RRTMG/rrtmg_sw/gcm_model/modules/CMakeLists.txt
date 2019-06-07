# Non-default name
esma_set_this(OVERRIDE RRTMG_SW_mods)

set (srcs
  parrrsw.F90 rrsw_aer.F90 rrsw_cld.F90 rrsw_con.F90
  rrsw_kg16.F90 rrsw_kg17.F90 rrsw_kg18.F90 rrsw_kg19.F90
  rrsw_kg20.F90 rrsw_kg21.F90 rrsw_kg22.F90 rrsw_kg23.F90
  rrsw_kg24.F90 rrsw_kg25.F90 rrsw_kg26.F90 rrsw_kg27.F90
  rrsw_kg28.F90 rrsw_kg29.F90 rrsw_ref.F90  rrsw_tbl.F90
  rrsw_vsn.F90 rrsw_wvn.F90
  )

esma_add_library (${this} SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL_Base)

