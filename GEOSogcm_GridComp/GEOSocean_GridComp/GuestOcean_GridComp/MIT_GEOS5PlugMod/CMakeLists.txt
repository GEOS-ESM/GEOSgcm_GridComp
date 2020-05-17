esma_set_this ()

set (SRCS
  MIT_GEOS5PlugMod.F90
  )

SET (MIT_EXP_ID  c90_llc90_05)
SET (MITGCM_ROOTDIR ${CMAKE_CURRENT_SOURCE_DIR}/@mit)
SET (MITGCM_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup)
SET (MITGCM_EXP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/experiments/${MIT_EXP_ID})
#SET (MIT_LIB_FILE libmitgcmuv.a libmitgcmuv_driver.a)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mitgcm_setup
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${MITGCM_EXP_DIR}/code
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup)

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL lib_mit lib_mit_interfaces
  INCLUDES (${MITGCM_BUILD_DIR}/code_split_driver)
  )

add_custom_target(
	make_local_build_dir
	COMMAND test -e local_build || mkdir local_build
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/
)

add_custom_target(
	make_mit_src
	COMMENT "building MITOCEAN ..."
	# this runs genmake2 and "make depend"
	COMMAND  ${MITGCM_ROOTDIR}/tools/genmake2
                 -rootdir ${MITGCM_ROOTDIR}
                 -mods ../code -mpi
		 -of ../build_options
		 #-of ${MITGCM_ROOTDIR}/tools/build_options/linux_amd64_ifort+mpi_ice_nas
		 -ieee
	COMMAND make depend
	DEPENDS make_local_build_dir
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/local_build
)

add_custom_target(
	fix_mit_inc
	COMMAND ./mk_local ../local_build
	DEPENDS make_mit_src
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/inc
)

add_custom_target(
	lib_mit
	COMMAND ./mk_local ../local_build
	COMMAND make install
	DEPENDS fix_mit_inc
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/build
)

add_custom_target(
	lib_mit_interfaces
	COMMAND make install
	DEPENDS lib_mit
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mitgcm_setup/code_split_driver
)
