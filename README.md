# Standalone of `compute_uwshcu_inv` routine in Moist Physics

## Brief Description

This standalone replicates the computations performed by the `compute_uwshcu_inv` subroutine in [`uwshcu.F90`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/2.1.3%2Bmoist-zero-divide-fix%2Baqua-planet-def/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/uwshcu.F90#L74).  The codebase was developed based on the `GEOSgcm_GridComp` branch [`2.1.3+moist-zero-divide-fix+aqua-planet-def`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/tree/2.1.3%2Bmoist-zero-divide-fix%2Baqua-planet-def).  MAPL data calls are substituted with file IO from a c180 GEOS aquaplanet run.

The standalone can be built for execution on CPUs.

## Building and Running

1. Edit the [`Makefile`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/orphan/openacc/moist/uw_shal_conv/Makefile)
    - Set `FC` in the `Makefile` to a Fortran compiler.  During development, the standalone has been compiled with `ifort`, `gfortran`, and `nvfortran`.
    - Choose the appropriate optimization `OPT` based on the Fortran compiler.  There are comments in the `Makefile` to guide the choice for `OPT`.

2. Run `make` to build.  This will create a binary called `TEST_MOIST`.

3. Run `./TEST_MOIST /discover/nobackup/projects/geosongpu/physics_standalone_data/moist/uwshcu/ <Dataset Number>`
    - `/discover/nobackup/projects/geosongpu/physics_standalone_data/moist/uwshcu/` contains the input and comparison data for the standalone.
    - `<Dataset Number>` can be set as an integer from `0` to `5`.

## Other Notes
- The standalone contains OpenACC directives that can compile for either CPU and GPU using `nvfortran`.  The standalone does not execute successfully on GPUs; however, it does exhibit performance improvements when executing with OpenACC on CPUs.
