# Standalone of `CUP_gf` version 2020 subroutine in Moist Physics

## Brief Description

This standalone replicates the computations performed by the `CUP_gf` subroutine in [`ConvPar_GF2020.F90`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/2.1.3%2Bmoist-zero-divide-fix%2Baqua-planet-def/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/ConvPar_GF2020.F90#L1859).  The codebase was developed based on the `GEOSgcm_GridComp` branch tagged [`v2.1.2`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/tree/v2.1.2).  MAPL data calls are substituted with file IO from a c180 GEOS aquaplanet run.  

The standalone can be built for execution on CPUs.

## Building and Running

1. Edit the [`Makefile`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/orphan/openacc/moist/cup_gf_gf2020/Makefile)
    - Set `FC` in the `Makefile` to a Fortran compiler.  During development, the standalone has been compiled with `ifort`, `gfortran`, and `nvfortran`.
    - Choose the appropriate optimization `OPT` based on the Fortran compiler.  There are comments in the `Makefile` to guide the choice for `OPT`.

2. Run `make` to build.  This will create a binary called `TEST_MOIST`.

3. Run `./TEST_MOIST /discover/nobackup/projects/geosongpu/physics_standalone_data/moist/cup_gf_GF2020/ <Dataset Number>`
    - `/discover/nobackup/projects/geosongpu/physics_standalone_data/moist/cup_gf_GF2020/` contains the input and comparison dataset for the standalone.
    - `<Dataset Number>` can be set as an integer from `0` to `5`.

## Other Notes
- The standalone verifies with the comparison dataset.
- The stack size may have to be set to unlimited (ex: `ulimit -s unlimited`) to avoid seg faults.
