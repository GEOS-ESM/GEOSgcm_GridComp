# Builds LISF's `LIS` library natively under CMake (replacing the previous
# approach of shelling out to LISF's own GNU-make build via `make nuopc`).
#
# This is a first-pass port and has NOT been build-verified in this
# environment (no Fortran compiler available here) -- treat it as a
# starting point to iterate on once it's actually compiled.
#
# BASEDIR must point to a Baselibs install (e.g. -DBASEDIR=$BASEDIR/Linux,
# matching how Baselibs is found elsewhere in this build) -- it is used
# here to locate HDF4/HDF-EOS2, which (unlike NetCDF/HDF5/ESMF) do not
# already have CMake targets set up by @cmake/external_libraries/FindBaselibs.cmake.
if (NOT BASEDIR)
  if (DEFINED ENV{BASEDIR})
    set (BASEDIR $ENV{BASEDIR})
  else ()
    message (FATAL_ERROR "BASEDIR is not set. Pass -DBASEDIR=<path> or set the BASEDIR environment variable.")
  endif ()
endif ()

find_package (Python3 COMPONENTS Interpreter REQUIRED)

# @LISF lives one directory up from here (GEOSlis_GridComp/@LISF), not under
# LISF_cmake/ itself.
esma_mepo_style (LISF mepo_LISF REL_PATH ..)

set (LISF_LIS  ${CMAKE_CURRENT_SOURCE_DIR}/${mepo_LISF}/lis)
set (LISF_MAKE ${LISF_LIS}/make)

#-----------------------------------------------------------------------
# user.cfg: plugin selection. VIC is disabled (restricted/unsupported).
# CABLE is disabled too -- cable_canopy.f90 fails to compile here with
# array shape-mismatch errors (rbw/poolcoef1*) unrelated to our Plug,
# which only drives Noah/NoahMP LSMs anyway. Noah.3.9 is disabled -- its
# noah39_main.F90 calls SFCDIF_OFF with more actual than dummy arguments
# (a vendored bug); we use Noah.3.3/NoahMP.3.6/NoahMP.4.0.1, not Noah.3.9.
# DA OBS pildas is disabled -- read_pildassmobs.F90 fails to compile with
# ifort (name conflicts/reuse in specification expressions, a vendored
# bug), unrelated to our Plug, which doesn't use this DA obs plugin.
# Everything else in LISF's default.cfg is left at its default (mostly On).
#-----------------------------------------------------------------------
file (WRITE ${LISF_MAKE}/user.cfg
"
VIC.4.1.1: Off
VIC.4.1.2: Off
CABLE: Off
Noah.3.9: Off
DA OBS pildas: Off
")

#-----------------------------------------------------------------------
# LIS_misc.h: build-time feature switches read by plugins.py and by LISF's
# Fortran sources. COUPLED must be #define'd here (unlike the original
# make-based build's setting) so plugins.py excludes lis/offline -- that
# directory's lisdrv.F90 is a `program` main entry point that would
# collide with GEOS.x's own main when linked into GEOS.
#-----------------------------------------------------------------------
file (WRITE ${LISF_MAKE}/LIS_misc.h
"
#define USE_NETCDF4
#undef USE_GRIBAPI
#undef AFWA_GRIB_CONFIGS
#define USE_HDF4
#define USE_HDFEOS2
#define USE_HDF5
#undef USE_MINPACK
#undef RTMS
#undef PETSc
#undef INC_WATER_PTS
#define COUPLED
#undef USAF_LIS75_SMDA
")

file (WRITE ${LISF_MAKE}/LIS_NetCDF_inc.h
"
#define NETCDF_shuffle  1
#define NETCDF_deflate  1
#define NETCDF_deflate_level  1
")

#-----------------------------------------------------------------------
# Generate Filepath + LIS_plugins.h via LISF's own plugin selector
# (plugins.py reads default.cfg/user.cfg and picks ~180 of LISF's ~200
# plugin directories by default) -- this is LISF-owned logic, so we defer
# to it here rather than re-deriving the same directory list by hand.
# Re-run on every configure so editing user.cfg above always takes effect.
#-----------------------------------------------------------------------
execute_process (
  COMMAND ${Python3_EXECUTABLE} plugins.py
  WORKING_DIRECTORY ${LISF_MAKE}
  RESULT_VARIABLE _plugins_result
)
if (NOT _plugins_result EQUAL 0)
  message (FATAL_ERROR "LISF's plugins.py failed (exit code ${_plugins_result})")
endif ()

#-----------------------------------------------------------------------
# Parse Filepath ("dirs := . ../core ../plugins ...") into real source
# lists, mirroring the GNU Makefile's own non-recursive per-directory
# file(GLOB) (lis/make/Makefile's FIND_FILES/FIND_HEADERS).
#-----------------------------------------------------------------------
file (STRINGS ${LISF_MAKE}/Filepath _filepath_lines)
set (_lis_dirs)
foreach (_line ${_filepath_lines})
  if (_line MATCHES "^dirs[ \t]*:=(.*)$")
    separate_arguments (_tokens UNIX_COMMAND "${CMAKE_MATCH_1}")
    foreach (_tok ${_tokens})
      get_filename_component (_abs_dir "${LISF_MAKE}/${_tok}" ABSOLUTE)
      list (APPEND _lis_dirs ${_abs_dir})
    endforeach ()
  endif ()
endforeach ()
list (REMOVE_DUPLICATES _lis_dirs)

set (LIS_SOURCES)
set (LIS_HEADER_DIRS ${LISF_MAKE})
# Some same-named files appear in more than one Filepath directory (e.g.
# get_cdf_params.F90 under both metforcing/mogreps_g and metforcing/galwem_ge).
# The GNU Makefile's VPATH silently shadows later duplicates, only building
# the first match; replicate that here to avoid "multiple definition" link
# errors from compiling and linking both copies.
set (_lis_source_names)
foreach (_dir ${_lis_dirs})
  file (GLOB _srcs ${_dir}/*.F90 ${_dir}/*.f90 ${_dir}/*.F ${_dir}/*.f ${_dir}/*.c ${_dir}/*.cc)
  foreach (_src ${_srcs})
    get_filename_component (_src_name ${_src} NAME)
    if (NOT _src_name IN_LIST _lis_source_names)
      list (APPEND _lis_source_names ${_src_name})
      list (APPEND LIS_SOURCES ${_src})
    endif ()
  endforeach ()
  file (GLOB _hdrs ${_dir}/*.h)
  if (_hdrs)
    list (APPEND LIS_HEADER_DIRS ${_dir})
  endif ()
endforeach ()

if (NOT LIS_SOURCES)
  message (FATAL_ERROR "No LISF sources found under any Filepath directory -- did plugins.py run correctly?")
endif ()

#-----------------------------------------------------------------------
# HDF4 / HDF-EOS2 -- Baselibs ships these, but (unlike NetCDF/HDF5/ESMF)
# @cmake/external_libraries/FindBaselibs.cmake does not create targets for
# them since MAPL itself does not need them.
#-----------------------------------------------------------------------
find_library (HDF4_MFHDF_LIBRARY NAMES mfhdf PATHS ${BASEDIR}/lib REQUIRED)
find_library (HDF4_DF_LIBRARY NAMES df PATHS ${BASEDIR}/lib REQUIRED)
find_library (HDFEOS_LIBRARY NAMES hdfeos PATHS ${BASEDIR}/lib REQUIRED)
find_library (GCTP_LIBRARY NAMES Gctp PATHS ${BASEDIR}/lib REQUIRED)
find_library (JPEG_LIBRARY NAMES jpeg PATHS ${BASEDIR}/lib REQUIRED)

add_library (hdf4::hdf4 STATIC IMPORTED)
set_target_properties (hdf4::hdf4 PROPERTIES
  IMPORTED_LOCATION ${HDF4_MFHDF_LIBRARY}
  INTERFACE_INCLUDE_DIRECTORIES ${BASEDIR}/include/hdf
  INTERFACE_LINK_LIBRARIES "${HDF4_DF_LIBRARY};${JPEG_LIBRARY}"
)

add_library (hdfeos::hdfeos STATIC IMPORTED)
set_target_properties (hdfeos::hdfeos PROPERTIES
  IMPORTED_LOCATION ${HDFEOS_LIBRARY}
  INTERFACE_INCLUDE_DIRECTORIES ${BASEDIR}/include/hdfeos
  INTERFACE_LINK_LIBRARIES ${GCTP_LIBRARY}
)

#-----------------------------------------------------------------------
# libLIS
# NOTE: on a from-scratch build, do the first build of this target with
# -j1 (or plain `make`/`cmake --build .` without --parallel). With ~thousands
# of sources in one target, a fully-parallel first pass can start compiling
# a file before the .F90 that `USE`s it has finished/been scanned, causing
# spurious "Error in opening the compiled module file" failures. Once built
# once, .mod files + depend.make are complete and later -jN rebuilds are fine.
#-----------------------------------------------------------------------
esma_set_this (OVERRIDE LIS)

# hdf4::hdf4/hdfeos::hdfeos are only needed for the default-enabled DA OBS
# plugins that call HDF4/HDF-EOS directly (#include "hdf.f90" + gdopen/
# gdattach/gdrdfld/gddetach/gdclose): AMSRE_SWE, ANSA_SNWD,
# GCOMW_AMSR2L3SND, GLASS_Albedo, GLASS_LAI, MODISsca, NASA_AMSREsm,
# PMW_snow, SSMI_SNWD. Drop these two libs if all 9 are ever disabled.
esma_add_library (${this}
  SOURCES ${LIS_SOURCES}
  PUBLIC_LIBS esmf NetCDF::NetCDF_Fortran hdf5::hdf5 hdf5::hdf5_hl hdfeos::hdfeos hdf4::hdf4
  TYPE SHARED
)

target_include_directories (${this} PRIVATE ${LIS_HEADER_DIRS})

# CMake would otherwise pick CXX (icpx) as the linker driver here since the
# sources include some .c/.cc files; but ESMF's esmf.mk link flags/libs
# (e.g. -threads -cxxlib, -lifport -lifcoremt) are classic ifort/icpc-only
# and icpx rejects them. Force the Fortran (ifort) driver, which understands them.
set_target_properties (${this} PROPERTIES LINKER_LANGUAGE Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  target_compile_options (${this} PRIVATE -fp-model precise -traceback)
  # -names/-convert/-assume/-nomixed-str-len-arg are ifort/ifx-only flags;
  # icx (the C compiler in this toolchain) rejects them, so restrict them
  # to Fortran sources.
  target_compile_options (${this} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:-nomixed-str-len-arg -names lowercase -convert big_endian -assume byterecl>
  )
endif ()

target_compile_definitions (${this} PRIVATE
  SPMD HIDE_SHR_MSG NO_SHR_VMATH IFC LINUX USE_INCLUDE_MPI LIS_JULES
)
