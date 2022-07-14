esma_set_this ()

set (SRCS
  MIT_GEOS5PlugMod.F90
  )

if (NOT DEFINED MIT_CONFIG_ID)
   message("You need to specify MIT_CONFIG_ID")
   return()
endif() 

SET (MITGCM_ROOTDIR ${CMAKE_CURRENT_SOURCE_DIR}/@mit)
SET (MITGCM_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup)
SET (MITGCM_CONFIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/configs/${MIT_CONFIG_ID})

#SET (MITLIBS ${MITGCM_BUILD_DIR}/code_split_driver/libmitgcmuv_driver.a)
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

#file (GLOB OBJS CONFIGURE_DEPENDS ${MITGCM_BUILD_DIR}/build/*.o
#                ${MITGCM_BUILD_DIR}/code_split_driver/*.o)


#SET_SOURCE_FILES_PROPERTIES(
#  ${MITLIBS}
#  PROPERTIES
#  EXTERNAL_OBJECT true
#  GENERATED true
#)

#message(STATUS "OBJS set to ${OBJS}")

esma_add_library (${this}
  SRCS ${SRCS} #${OBJS}
#  DEPENDENCIES MAPL mitgcm_setup/code_split_driver/driverdone
  DEPENDENCIES MAPL $<BUILD_INTERFACE:mitgcmuv_driver> $<BUILD_INTERFACE:mitgcmuv>
  INCLUDES $<BUILD_INTERFACE:${MITGCM_BUILD_DIR}/code_split_driver>
#  TYPE SHARED
)

#add_dependencies(${this} 
#	${MITLIBS})

#target_link_libraries(${this} PUBLIC
#  libmitgcmuv_driver.so
#)

#target_link_libraries(${this} PUBLIC
#  libmitgcmuv.a
#)

target_link_directories(${this} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/code_split_driver>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/build>
)

add_custom_command(
	OUTPUT mitgcm_setup/local_build
	COMMAND mkdir local_build
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/
)

add_custom_command(
	OUTPUT mitgcm_setup/local_build/Makefile
	COMMENT "building MITOCEAN ..."
	# this runs genmake2 and "make depend"
	COMMAND  ${MITGCM_ROOTDIR}/tools/genmake2
                 -rootdir ${MITGCM_ROOTDIR}
                 -mods ../code -mpi
		 -of ../${MITGCM_BUILD_OPTIONS}
		 ${MITGCM_EXTRA_BUILD_OPTIONS}
	COMMAND make depend
	DEPENDS mitgcm_setup/local_build
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/local_build
)

add_custom_command(
	OUTPUT mitgcm_setup/inc/incdone
	COMMAND ./mk_local ../local_build
	COMMAND touch incdone
	DEPENDS mitgcm_setup/local_build/Makefile
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/inc
)

#add_custom_command(
#	OUTPUT mitgcm_setup/build/builddone
add_custom_target(mitgcmuv
	COMMAND ./mk_local ../local_build
	COMMAND make install
#	COMMAND touch builddone
	DEPENDS mitgcm_setup/inc/incdone
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/build
)

#add_custom_command(
#	OUTPUT ${MITLIBS}
add_custom_target(mitgcmuv_driver
	COMMAND make install FC=${FC} CC=${CC}
#	DEPENDS mitgcm_setup/build/builddone
	DEPENDS mitgcmuv
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/code_split_driver
)

