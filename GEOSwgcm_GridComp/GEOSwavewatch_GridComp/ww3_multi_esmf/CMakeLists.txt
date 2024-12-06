esma_set_this ()


list( APPEND WW3ESMF_FTN
    constants.ftn
    w3adatmd.ftn
    w3arrymd.ftn
    w3cspcmd.ftn
    w3dispmd.ftn
    w3fldsmd.ftn
    w3flx1md.ftn
    w3gdatmd.ftn
    w3gsrumd.ftn
    w3idatmd.ftn
    w3initmd.ftn
    w3iobcmd.ftn
    w3iogomd.ftn
    w3iogrmd.ftn
    w3iopomd.ftn
    w3iorsmd.ftn
    w3iosfmd.ftn
    w3iotrmd.ftn
    w3nmlmultimd.ftn
    w3odatmd.ftn
    w3parall.ftn
    w3partmd.ftn
    w3pro3md.ftn
    w3profsmd.ftn
    w3sbt1md.ftn
    w3sdb1md.ftn
    w3servmd.ftn
    w3snl1md.ftn
    w3sln1md.ftn
    w3fld1md.ftn
    w3fld2md.ftn
    w3src4md.ftn
    w3srcemd.ftn
    w3timemd.ftn
    w3triamd.ftn
    w3updtmd.ftn
    w3uqckmd.ftn
    w3wavemd.ftn
    w3wdasmd.ftn
    w3wdatmd.ftn
    wmmaplmd.ftn
    wmfinlmd.ftn
    wmgridmd.ftn
    wminiomd.ftn
    wminitmd.ftn
    wmiopomd.ftn
    wmmdatmd.ftn
    wmscrpmd.ftn
    wmunitmd.ftn
    wmupdtmd.ftn
    wmwavemd.ftn
    SCRIP/scrip_constants.f
    SCRIP/scrip_errormod.f90
    SCRIP/scrip_grids.f
    SCRIP/scrip_interface.ftn
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

get_filename_component(aux_dir ${WW3_path}/model/aux ABSOLUTE)
get_filename_component(ftn_dir ${WW3_path}/model/ftn ABSOLUTE)

message(DEBUG "WW3 aux_dir ${aux_dir}")
message(DEBUG "WW3 ftn_dir ${ftn_dir}")

add_executable(w3adc "${aux_dir}/w3adc.f")

set (WW3ESMF_F90)
foreach(src_file ${WW3ESMF_FTN})
    STRING(REGEX REPLACE ".ftn" ".F90" gen_src_file ${src_file})
    STRING(REGEX REPLACE "/" "_" gen_log_file ${gen_src_file})
    # Testing shows that we only want BYPRODUCTS here if our CMAKE_GENERATOR is Unix Makefiles
    # If we use Ninja, we don't want BYPRODUCTS
    if (CMAKE_GENERATOR STREQUAL "Unix Makefiles")
      set(BYPRODUCTS ${gen_src_file})
    else()
      set(BYPRODUCTS "")
    endif()
    add_custom_command(
        OUTPUT  ${gen_src_file}
        BYPRODUCTS  ${BYPRODUCTS}
        DEPENDS w3adc ${ftn_dir}/${src_file}
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_w3adc.sh ${ftn_dir} ${src_file} > ${gen_log_file}.w3adc.log 2>&1
        COMMENT "Running w3adc ${src_file}")
    list(APPEND WW3ESMF_F90 ${gen_src_file})
endforeach()


esma_add_library (${this}
   SRCS ${WW3ESMF_F90}
   DEPENDENCIES MAPL ESMF::ESMF NetCDF::NetCDF_Fortran
  )

target_include_directories (${this} PRIVATE
  $<BUILD_INTERFACE:${ftn_dir}>
  )
