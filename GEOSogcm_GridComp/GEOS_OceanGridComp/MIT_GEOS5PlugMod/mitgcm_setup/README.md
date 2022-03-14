Contains code to interface between standard MITgcm and GEOS-5/MAPL world.

* `inc/' This directory contains a copy of some MITgcm header files, with
    comment "C" markers changed to "!" instead - to allow for fixed form and
    free form compile use (we think). This needs to be cleaned up. The files
    should be generated from a script somehow.

* `code/` This is a directory of standard MITgcm code/ customization for a
    specific experiment (global_ocean_c32 in this case).

* `code_split_driver/` This is core F90 code that connects MITgcm data
    structutures and MAPL data structures. The code is in a tree of files but
    Makefile dependencies are currently for a flat set of files in the
    `code_split_driver/` directory.

* code needs to compatible with whatever configuration/checkpoint of MITgcm has
    been selected. The MITgcm checkpoint source are linked in with soft links
    under gmao_mitgcm_couplng/mitgcm_setup/build by the Makefile
