esma_set_this ()

set (SRCS
  GEOS_MITDynaGridComp.F90 
  )

if (NOT DEFINED MIT_CONFIG_ID)
   message("You need to specify MIT_CONFIG_ID")
   return()
endif() 

SET (MITGCM_ROOTDIR ${CMAKE_CURRENT_SOURCE_DIR}/../../GEOS_OceanGridComp/MIT_GEOS5PlugMod/@mit)
SET (MITGCM_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup)
SET (MITGCM_CONFIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/configs/${MIT_CONFIG_ID})

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
     SET (MITGCM_BUILD_OPTIONS build_options.gnu)
     set (FC gfortran)
     set (CC gcc)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
     SET (MITGCM_BUILD_OPTIONS build_options.intel)
     # Path to Intel compiler to overwrite defaults for mpif90
     set (FC ifort)
     set (CC icc)
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
     SET (MITGCM_EXTRA_BUILD_OPTIONS "-ieee")
elseif()
     SET (MITGCM_EXTRA_BUILD_OPTIONS "")
endif()


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mitgcm_setup
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${MITGCM_CONFIG_DIR}/code
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup)

file (GLOB ICEOBJS CONFIGURE_DEPENDS ${MITGCM_BUILD_DIR}/build/*.o
                ${MITGCM_BUILD_DIR}/code_split_driver/*.o)

SET_SOURCE_FILES_PROPERTIES(
  ${ICEOBJS}
  PROPERTIES
  EXTERNAL_OBJECT true
  GENERATED true
)

esma_add_library (${this}
  SRCS ${SRCS} ${ICEOBJS}
  DEPENDENCIES MAPL CICE4 $<BUILD_INTERFACE:mitgcmseaice_driver> $<BUILD_INTERFACE:mitgcmseaice>
  INCLUDES $<BUILD_INTERFACE:${MITGCM_BUILD_DIR}/code_split_driver>
  TYPE SHARED
)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

target_link_directories(${this} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/code_split_driver>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/build>
)

add_custom_target(
	make_ice_local_build_dir
	COMMAND test -e local_build || mkdir local_build
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/
)

add_custom_target(
	make_ice_mit_src
	COMMENT "building MITOCEAN ..."
	# this runs genmake2 and "make depend"
	COMMAND  ${MITGCM_ROOTDIR}/tools/genmake2
                 -rootdir ${MITGCM_ROOTDIR}
                 -mods ../code -mpi
		 -of ../${MITGCM_BUILD_OPTIONS}
		 ${MITGCM_EXTRA_BUILD_OPTIONS}
	COMMAND make depend
	DEPENDS make_ice_local_build_dir
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/local_build
)

add_custom_target(
	fix_ice_mit_inc
	COMMAND ./mk_local ../local_build
	DEPENDS make_ice_mit_src
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/inc
)

add_custom_target(
	mitgcmseaice
	COMMAND ./mk_local ../local_build
	COMMAND make install
	DEPENDS fix_ice_mit_inc
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/build
)

add_custom_target(
	mitgcmseaice_driver
	COMMAND make install FC=${FC} CC=${CC}
	DEPENDS mitgcmseaice
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/code_split_driver
)



