esma_set_this ()

# model/ftn/wmesmfmd.ftn~
# model/ftn/w3iogrmd.ftn~
# model/ftn/wminitmd.ftn~
# model/ftn/w3timemd.ftn~
# model/ftn/w3iorsmd.ftn~
# 

# endiannes
# --------------------------------
#          GEOS    WGCM   WW3esmf
# --------------------------------
# CVS/J3   big     big    big
# git      little  little little 
# --------------------------------

set (ESMFMKFILE "${BASEDIR}/lib/esmf.mk")


set (alldirs
    #   @ww3-release/model/esmf
    ww3_multi_esmf
  )


if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_WaveWatchGridComp.F90)

  esma_add_library (${this}
    SRCS GEOS_WaveWatchGridComp.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL
    INCLUDES ${INC_ESMF}
    )

else ()

  esma_add_subdirectories (${alldirs})

endif()


file (GLOB_RECURSE rc_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
