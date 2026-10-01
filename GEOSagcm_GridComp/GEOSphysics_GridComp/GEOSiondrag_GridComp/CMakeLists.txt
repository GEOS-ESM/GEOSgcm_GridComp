esma_set_this ()

set (srcs
  ion_drag_module.F90
  IRI_input_module.F90
  GEOS_IonDragGridComp.F90
  IRI/irisub.for
  IRI/irifun.for
  IRI/iritec.for
  IRI/iridreg.for
  IRI/iriflip.for
  IRI/irirtam.for
  IRI/cira.for
  IRI/igrf.for
  IRI/rocdrift.for
  )

set (resource_files
  IonDrag_GridComp.rc
  )

install( FILES ${resource_files}
  DESTINATION etc
  )

# IRI opens its coefficient files by bare filename at runtime.
# GEOSgcm_App copies ${GEOSDIR}/etc/* into the experiment RC directory and
# later copies ${EXPDIR}/RC/* into the run directory without recursion.
# Install IRI runtime data directly into etc so each file reaches the run
# directory as a regular file rather than remaining under an IRI subdirectory.
set (iri_data_files
  IRI/ig_rz.dat
  IRI/apf107.dat
  IRI/ccir11.asc IRI/ccir12.asc IRI/ccir13.asc IRI/ccir14.asc IRI/ccir15.asc
  IRI/ccir16.asc IRI/ccir17.asc IRI/ccir18.asc IRI/ccir19.asc IRI/ccir20.asc
  IRI/ccir21.asc IRI/ccir22.asc
  IRI/ursi11.asc IRI/ursi12.asc IRI/ursi13.asc IRI/ursi14.asc IRI/ursi15.asc
  IRI/ursi16.asc IRI/ursi17.asc IRI/ursi18.asc IRI/ursi19.asc IRI/ursi20.asc
  IRI/ursi21.asc IRI/ursi22.asc
  IRI/dgrf1945.dat IRI/dgrf1950.dat IRI/dgrf1955.dat IRI/dgrf1960.dat
  IRI/dgrf1965.dat IRI/dgrf1970.dat IRI/dgrf1975.dat IRI/dgrf1980.dat
  IRI/dgrf1985.dat IRI/dgrf1990.dat IRI/dgrf1995.dat IRI/dgrf2000.dat
  IRI/dgrf2005.dat IRI/dgrf2010.dat IRI/dgrf2015.dat IRI/dgrf2020.dat
  IRI/igrf2025.dat IRI/igrf2025s.dat
  IRI/mcsat11.dat IRI/mcsat12.dat IRI/mcsat13.dat IRI/mcsat14.dat
  IRI/mcsat15.dat IRI/mcsat16.dat IRI/mcsat17.dat IRI/mcsat18.dat
  IRI/mcsat19.dat IRI/mcsat20.dat IRI/mcsat21.dat IRI/mcsat22.dat
  )

install( FILES ${iri_data_files}
  DESTINATION etc
  )

# Install the GEOS-MLT ML ion-velocity runtime files.
#
# Expected source-tree layout:
#
#   @GEOSmlionvel_GridComp/
#     geos_mlionvel_driver.py
#     mlionvel_uivi_colfilm_2020-2022.pth
#     NormStats_IonVel_UIVI_2020-2022_train.npz
#
# The Python driver is installed into lib/Python so MAPL_PythonBridge can
# import it by module name. Model and normalization files are installed
# under share/GEOSmlionvel, matching geos_mlionvel_driver.py.
set (MLIONVEL_RUNTIME_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/@GEOSmlionvel_GridComp
  )

install(
  FILES
    ${MLIONVEL_RUNTIME_DIR}/geos_mlionvel_driver.py
  DESTINATION
    lib/Python
  )

install(
  FILES
    ${MLIONVEL_RUNTIME_DIR}/mlionvel_uivi_colfilm_2020-2022.pth
  DESTINATION
    share/GEOSmlionvel/models
  )

install(
  FILES
    ${MLIONVEL_RUNTIME_DIR}/NormStats_IonVel_UIVI_2020-2022_train.npz
  DESTINATION
    share/GEOSmlionvel/norms
  )

esma_add_library (
  ${this}
  SRCS ${srcs}
  DEPENDENCIES GEOS_Shared MAPL ESMF::ESMF NetCDF::NetCDF_Fortran fvdycore TYPE SHARED
  )

mapl_acg (
  ${this} IonDrag_StateSpecs.rc
  IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS
  )
