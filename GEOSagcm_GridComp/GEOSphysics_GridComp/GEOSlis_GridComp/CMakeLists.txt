esma_set_this ()

# GEOS_LisGridComp.F90 wraps LISF's core Fortran API directly
# (LIS_domain_init, LIS_surfaceModel_run, ...); no NUOPC/ESMF coupling
# layer is involved.

# LisFieldCatalog.F90 holds a self-contained, static copy of
# LIS_NUOPC_Gluecode's field stdName/stateName/units metadata and role
# (adImport/adExport/directConn) classification, so LisFields.F90
# doesn't need to build/link LISF's NUOPC cap (or LIS_FORC_AttributesMod)
# just to get this static data.
# LisFields.F90 holds the shared Import/Export/Internal field
# catalog derived from it.
set (SRCS LisFieldCatalog.F90 LisFields.F90 GEOS_LisGridComp.F90)

# LISF_cmake/ builds the `LIS` library target natively under CMake from
# @LISF's own source tree (see LISF_cmake/CMakeLists.txt for details).
add_subdirectory (LISF_cmake)

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL LIS
  TYPE SHARED
)

# Force the Fortran compiler driver as the linker -- otherwise CMake picks
# CXX (icpx) here, same as for LIS itself, and ESMF's link flags (e.g.
# -cxxlib, -threads from esmf.mk) are meant for the Fortran/MPI driver.
set_target_properties (${this} PROPERTIES LINKER_LANGUAGE Fortran)

# Baselibs' HDF4/HDF-EOS2 static archives (libmfhdf.a, libdf.a, libhdfeos.a,
# libhdf5_fortran.a, ...), pulled in transitively via LIS, weren't built
# with -fPIC, so linking them into this SHARED library fails with
# "relocation R_X86_64_32 ... can not be used when making a shared object".
# -z notext tells the linker to allow the resulting text relocations
# instead of requiring PIC everywhere.
target_link_options (${this} PRIVATE -Wl,-z,notext)
