esma_set_this()

set (alldirs
  GEOSmkiau_GridComp
  )


#if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_GcmGridComp.F90)

 #ecbuild_declare_project()
  #esma_add_library(${this}
    #SRCS GEOS_GcmGridComp.F90
    ##SRCS gcm_stub.F90
    #SUBCOMPONENTS ${alldirs}
    #DEPENDENCIES MAPL modelE
    #INCLUDES ${INC_ESMF} $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include/modelE>)

  #ecbuild_install_project( NAME GEOSgcm_GridComp)

#else ()

  #esma_add_subdirectories (${alldirs})

#endif()

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include/modelE)
include(ExternalProject)
if(DEFINED ENV{ModelE_Template})
   set(ModelE_Template $ENV{ModelE_Template})
else()
   set(ModelE_Template geos_setup)
endif()
if(DEFINED ENV{ModelE_Run})
   set(ModelE_Run $ENV{ModelE_Run})
else()
   set(ModelE_Run geos_run)
endif()
if(DEFINED ENV{MODELERC})
   set(MODELERC $ENV{MODELERC})
else()
   set(MODELERC /home/bmauer/.modelErc)
endif()

message(INFO " MODELE TEMPLATE: ${ModelE_Template}")
message(INFO " MODELE RUN: ${ModelE_Run}")
message(INFO " MODELERC:  ${MODELERC}")
ExternalProject_Add(modelE
   PREFIX "modelE"
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/@Model5Eagcm_GridComp/decks
   DOWNLOAD_COMMAND ""
   CONFIGURE_COMMAND ""
   UPDATE_COMMAND ""
   BUILD_ALWAYS 1
   #BUILD_COMMAND make clean OVERWRITE=YES && make rundeck RUN=${ModelE_Run} RUNSRC=${ModelE_Template} OVERWRITE=YES && $(MAKE) setup RUN=${ModelE_Run} MAPL=YES OVERWRITE=YES GEOS_BINARY_DIR=${CMAKE_BINARY_DIR}  VERBOSE_OUTPUT=YES MODELERC=${MODELERC} MPI=YES
   BUILD_COMMAND make rundeck RUN=${ModelE_Run} RUNSRC=${ModelE_Template} OVERWRITE=YES 
   #BUILD_COMMAND make clean OVERWRITE=YES
         COMMAND make rundeck RUN=${ModelE_Run} RUNSRC=${ModelE_Template} OVERWRITE=YES 
         COMMAND $(MAKE) setup RUN=${ModelE_Run} MAPL=YES OVERWRITE=YES GEOS_BINARY_DIR=${CMAKE_BINARY_DIR}  VERBOSE_OUTPUT=YES MODELERC=${MODELERC} MPI=YES COMPILE_WITH_TRAPS=YES
         COMMAND find ../model/ -type f -name "*.a" -exec /bin/cp {} ${CMAKE_BINARY_DIR}/lib $<SEMICOLON>
         COMMAND find ../model/ -type f -name "*.mod" -exec /bin/cp {} ${CMAKE_BINARY_DIR}/include/modelE $<SEMICOLON>
   BUILD_IN_SOURCE 1
   INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
   INSTALL_COMMAND cmake -E make_directory <INSTALL_DIR>/lib
           COMMAND cmake -E make_directory <INSTALL_DIR>/bin
           COMMAND find <SOURCE_DIR>/../model/ -type f -name "*.a" -exec /bin/cp {} <INSTALL_DIR>/lib  $<SEMICOLON>
           COMMAND find <SOURCE_DIR>/../model/ -type f -name "*.mod" -exec /bin/cp {} <INSTALL_DIR>/include/modelE  $<SEMICOLON>
           COMMAND find <SOURCE_DIR>/../decks/ -type f -name "geos_run.exe" -exec /bin/cp {} <INSTALL_DIR>/bin  $<SEMICOLON>
   DEPENDS MAPL
)
set(modele_lib_dir ${CMAKE_BINARY_DIR}/lib CACHE PATH "ModelE library dir" FORCE)
set(modele_inc_dir ${CMAKE_BINARY_DIR}/include/modelE CACHE STRING "ModelE include dir" FORCE)
set(modele_model_lib ${modele_lib_dir}/libmodel.a CACHE PATH "ModelE library" FORCE)
#set(modele_model_lib "${modele_lib_dir}/libmodel.a ${modele_lib_dir}/libprofiler.a ${modele_lib_dir}/libdd2d.a ${modele_lib_dir}/libgiss_LSM.a ${modele_lib_dir}/libEnt.a ${modele_lib_dir}/libsolvers.a ${modele_lib_dir}/libMPI_Support.a ${modele_lib_dir}/libshared.a"  CACHE PATH "ModelE library" FORCE)
set(modele_dd2d_lib ${modele_lib_dir}/libdd2d.a CACHE PATH "ModelE library" FORCE)
set(modele_solvers_lib ${modele_lib_dir}/libsolvers.a CACHE PATH "ModelE library" FORCE)
set(modele_prof_lib ${modele_lib_dir}/libprofiler.a CACHE PATH "ModelE library" FORCE)
set(modele_ent_lib ${modele_lib_dir}/libEnt.a CACHE PATH "ModelE library" FORCE)
set(modele_gisslsm_lib ${modele_lib_dir}/libgiss_LSM.a CACHE PATH "ModelE library" FORCE)
set(modele_mpi_lib ${modele_lib_dir}/libMPI_Support.a CACHE PATH "ModelE library" FORCE)
set(modele_shared_lib ${modele_lib_dir}/libshared.a CACHE PATH "ModelE library" FORCE)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GEOS_GcmGridComp.F90)

 ecbuild_declare_project()
  esma_add_library(${this}
    SRCS GEOS_GcmGridComp.F90
    #SRCS gcm_stub.F90
    SUBCOMPONENTS ${alldirs}
    DEPENDENCIES MAPL
    INCLUDES ${INC_ESMF} $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include/modelE>)

  add_dependencies(${this} modelE)

  ecbuild_install_project( NAME GEOSgcm_GridComp)

else ()

  esma_add_subdirectories (${alldirs})

endif()
