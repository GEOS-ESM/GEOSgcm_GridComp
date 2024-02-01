# Standalone of `fillq2zero` subroutine in Moist Physics

## Brief Description

This standalone replicates the computations performed by the `fillq2zero` subroutine in [`aer_actv_single_moment.F90`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/v2.1.2/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/aer_actv_single_moment.F90#L39).  The codebase was developed based on the `GEOSgcm_GridComp` branch tagged [`v2.1.2`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/tree/v2.1.2).  MAPL data calls are substituted with file IO from a c180 GEOS aquaplanet run.  

The standalone can be built for execution on both CPUs and GPUs.

## Building and Running

1. Edit the [`Makefile`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/orphan/openacc/moist/radcoup_loop/Makefile)
    - Set `FC` in the `Makefile` to a Fortran compiler.  During development, the standalone has been compiled with `ifort`, `gfortran`, and `nvfortran`.
    - Choose the appropriate optimization `OPT` based on the Fortran compiler.  There are comments in the `Makefile` to guide the choice for `OPT`.

2. Run `make` to build.  This will create a binary called `TEST_MOIST`.

3. Run `./TEST_MOIST /discover/nobackup/projects/geosongpu/physics_standalone_data/moist/radcoup_loop/ <Dataset Number>`
    - `/discover/nobackup/projects/geosongpu/physics_standalone_data/moist/radcoup_loop/` contains the input and comparison dataset for the standalone.
    - `<Dataset Number>` can be set as an integer from `0` to `5`.

## Other Notes
- GPU execution is performed via OpenACC directives.
- The `nvfortran` compiler has been tested up to version 23.11 and can build and execute the standalone using OpenACC directives.
- The `gfortran` compiler has been tested up to version 11.4.0 and can build and execute the standalone using OpenACC directives.
- The standalone verifies with the comparison dataset.