# Standalone of `gfdl_cloud_microphys_driver` routine in Moist Physics

## Brief Description

This standalone replicates the computations performed by the `gfdl_cloud_microphys_driver` subroutine in [`gfdl_cloud_microphys.F90`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/2.1.3%2Bmoist-zero-divide-fix%2Baqua-planet-def/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/gfdl_cloud_microphys.F90#L345).  The codebase was developed based on the `GEOSgcm_GridComp` branch [`2.1.3+moist-zero-divide-fix+aqua-planet-def`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/tree/2.1.3%2Bmoist-zero-divide-fix%2Baqua-planet-def).  MAPL data calls are substituted with file IO from a c180 GEOS aquaplanet run.  

Currently, the standalone can be built for execution on both CPUs and GPUs.

## Building and Running

1. Edit the [`Makefile`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/orphan/openacc/moist/gfdl_cloud_microphysics/Makefile)
    - Set `FC` in the `Makefile` to a Fortran compiler.  During development, the standalone has been compiled with `ifort`, `gfortran`, and `nvfortran`.
    - Choose the appropriate optimization `OPT` based on the Fortran compiler.  There are comments in the `Makefile` to guide the choice for `OPT`.

2. Run `make` to build.  This will create a binary called `TEST_MOIST`.

3. Run `./TEST_MOIST /discover/nobackup/projects/geosongpu/physics_standalone_data/moist/gfdl_cloud_microphysics/ <Dataset Number>`
    - `/discover/nobackup/projects/geosongpu/physics_standalone_data/moist/gfdl_cloud_microphysics/` contains the input and comparison data for the standalone.
    - `<Dataset Number>` can be set as an integer from `0` to `5`.

## Other Notes
- GPU execution is performed via OpenACC directives.
- The `nvfortran` compiler has been tested on the standalone up to version 23.11 and can build and execute the standalone via OpenACC directives.
- The `gfortran` compiler has been tested on the standalone up to version 11.4.0 and can build the standalone with OpenACC directives; however, a `gfortran` compiled binary cannot successfully execute on device.
- The code "verifies" with the comparison dataset, but this may have to be examined more in detail.
- When compiling with `ifort`, the stack size may have to be set to unlimited (ex: With a bash shell, run `ulimit -s unlimited`) so that the standalone does not produce a segmentation fault.
