# OpenACC Port of `GFDL_1M_Run` routine in Moist Physics

## Brief Description

This standalone replicates the computations performed by the `GFDL_1M_Run` subroutine in [`GEOS_GFDL_1M_InterfaceMod.F90`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/2.1.3%2Bmoist-zero-divide-fix%2Baqua-planet-def/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90).  The codebase was developed based on the `GEOSgcm_GridComp` branch [`2.1.3+moist-zero-divide-fix+aqua-planet-def`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/tree/2.1.3%2Bmoist-zero-divide-fix%2Baqua-planet-def).  MAPL data calls are substituted with file IO from a c180 GEOS aquaplanet run.  Currently, the standalone can be built for execution on CPUs, but its output has not verified with the `discover` dataset.

## Building and Running

1. Edit the [`Makefile`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/orphan/openacc/moist/GFDL_1M_Run/Makefile)
    - Set `FC` in the `Makefile` to a Fortran compiler.  During development, the standalone has been compiled with `ifort`, `gfortran`, and `nvfortran`.
    - Choose the appropriate optimization `OPT` based on the Fortran compiler.  There are comments in the `Makefile` to guide the choice for `OPT`.

2. Run `make` to build.  This will create a binary called `TEST_MOIST`.

3. Run `./TEST_MOIST /discover/nobackup/projects/geosongpu/physics_standalone_data/moist/gfdl_1m_run/ <Dataset Number>`
    - `/discover/nobackup/projects/geosongpu/physics_standalone_data/moist/gfdl_1m_run/` contains the input and comparison data for the standalone.
    - `<Dataset Number>` can be set from `0` to `5`.

## Other Notes
- Currently, this standalone does not verify when executing on the CPU.
- There is no OpenACC / OpenMP pragmas in the codebase.
