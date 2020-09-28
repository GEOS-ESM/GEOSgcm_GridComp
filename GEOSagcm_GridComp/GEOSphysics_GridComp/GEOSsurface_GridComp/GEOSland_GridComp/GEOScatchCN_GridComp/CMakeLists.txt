# Note that directory and library name do not correspond!
# GEOS CatchCN model is essentially a two-part model: 1) CLM's Carbon-Nirtrogen module and 2) GMAO's unified land surface model.
# Subdiectories CLM40, CLM45 etc contain CN modules from different CLM versions while the common unified module is located in this directory. 
# As of 9/28/2020, CLM40 is the default CN module of the GEOS catchCN model. Note, the GEOS model uses the default CN module to compile/build the GEOScatchCN_GridComp library.
# SWITCH to change the CN module: Change below line 8 and line 18 of GEOSland_GridComp/CMakeLists.txt to the desired version of the CN module.
     
esma_set_this (OVERRIDE GEOS_CatchCNRoot)
esma_add_subdirectories (CLM45)

set (srcs
  catchmentCN.F90
  catchcn_iau.F90)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL GEOS_Shared GEOS_LandShared)
