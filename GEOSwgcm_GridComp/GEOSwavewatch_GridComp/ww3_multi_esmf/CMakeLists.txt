esma_set_this ()

list( APPEND WW3ESMF_SRCS
    constants.F90
    w3adatmd.F90
    w3arrymd.F90
    w3cspcmd.F90
    w3dispmd.F90
    w3fldsmd.F90
    w3gdatmd.F90
    w3flx1md.F90
    w3gsrumd.F90
    w3idatmd.F90
    w3initmd.F90
    w3iobcmd.F90
    w3iogomd.F90
    w3iogrmd.F90
    w3iopomd.F90
    w3iorsmd.F90
    w3iosfmd.F90
    w3iotrmd.F90
    w3nmlmultimd.F90
    w3odatmd.F90
    w3parall.F90
    w3partmd.F90
    w3pro3md.F90
    w3profsmd.F90
    w3sbt1md.F90
    w3sdb1md.F90
    w3servmd.F90
    w3snl1md.F90
    w3sln1md.F90
    w3fld1md.F90
    w3fld2md.F90
    w3src4md.F90
    w3srcemd.F90
    w3timemd.F90
    w3triamd.F90
    w3updtmd.F90
    w3uqckmd.F90
    w3wavemd.F90
    w3wdasmd.F90
    w3wdatmd.F90
    wmmaplmd.F90
    wmfinlmd.F90
    wmgridmd.F90
    wminiomd.F90
    wminitmd.F90
    wmiopomd.F90
    wmmdatmd.F90
    wmscrpmd.F90
    wmunitmd.F90
    wmupdtmd.F90
    wmwavemd.F90
    SCRIP/scrip_constants.f
    SCRIP/scrip_errormod.f90
    SCRIP/scrip_grids.f
    SCRIP/scrip_interface.F90
    SCRIP/scrip_iounitsmod.f90
    SCRIP/scrip_kindsmod.f90
    SCRIP/scrip_netcdfmod.f90
    SCRIP/scrip_remap_conservative.f
    SCRIP/scrip_remap_read.f
    SCRIP/scrip_remap_vars.f
    SCRIP/scrip_remap_write.f
    SCRIP/scrip_timers.f
)

esma_mepo_style(ww3 WW3_rel_path REL_PATH ..)
set (WW3_path ${CMAKE_CURRENT_SOURCE_DIR}/${WW3_rel_path})

set (SRCS)
foreach(src_file ${WW3ESMF_SRCS})
    list(APPEND SRCS ${WW3_path}/model/src/${src_file})
#    message(STATUS "WW3::file  ${src_file}")
#    message(STATUS "WW3::file  ${WW3_path}/model/src/${src_file}")
#    message(STATUS "")
endforeach()

esma_add_library (${this}
   SRCS ${SRCS}
   DEPENDENCIES MAPL esmf NetCDF::NetCDF_Fortran
  )


set (switch_strings   "NOGRB DIST MPI SCRIP PR3 UQ FLX0 LN1 FLD2 ST4 STAB0 NL1 BT1 DB1 MLIM TR0 BS0 XX0 WNX0 WNT0 CRX0 CRT0 O0 O1 O2 O3 O4 O5 O6 O7 IC0 IS0 REF0" )

set_property (SOURCE ${WW3_path}/model/src/w3initmd.F90
  APPEND
  PROPERTY COMPILE_DEFINITIONS
  "__WW3_SWITCHES__=\'${switch_strings}\'"
  )



list (APPEND defs 
    W3_NOGRB
    W3_DIST
    W3_MPI
    W3_SCRIP
    W3_PR3
    W3_UQ
    W3_FLX0
    W3_LN1
    W3_FLD2
    W3_ST4
    W3_STAB0
    W3_NL1
    W3_BT1
    W3_DB1
    W3_MLIM
    W3_TR0
    W3_BS0
    W3_XX0
    W3_WNX0
    W3_WNT0
    W3_CRX0
    W3_CRT0
    W3_O0
    W3_O1
    W3_O2
    W3_O3
    W3_O4
    W3_O5
    W3_O6
    W3_O7
    W3_IC0
    W3_IS0
    W3_REF0
)

target_compile_definitions(${this} PRIVATE ${defs})


