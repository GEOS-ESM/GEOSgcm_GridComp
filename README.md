# Standalone of Cloud Macrophysics routines in Moist Physics

## Brief Description

This standalone replicates the computations performed within the Cloud Macrophysics region in [`GEOS_GFDL_1M_InterfaceMod.F90`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/v2.1.2/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90#L530).  The codebase was developed based on the `GEOSgcm_GridComp` branch tagged [`v2.1.2`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/tree/v2.1.2).  MAPL data calls are substituted with file IO from a c180 GEOS aquaplanet run.

Currently, the standalone can be built for execution on CPUs and GPUs and has been verified with a `discover` dataset.

## Building and Running

1. Edit the [`Makefile`](https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/orphan/openacc/moist/evap_subl_pdf_loop/Makefile)
    - Set `FC` in the `Makefile` to a Fortran compiler.  During development, the standalone has been compiled with `ifort`, `gfortran`, and `nvfortran`.
    - Choose the appropriate optimization `OPT` based on the Fortran compiler.  There are comments in the `Makefile` to guide the choice for `OPT`.

2. Run `make` to build.  This will create a binary called `TEST_MOIST`.

3. Run `./TEST_MOIST /discover/nobackup/projects/geosongpu/physics_standalone_data/moist/evap_subl_pdf_loop/ <Dataset Number>`
    - `/discover/nobackup/projects/geosongpu/physics_standalone_data/moist/evap_subl_pdf_loop/` contains the input and comparison data for the standalone.
    - `<Dataset Number>` can be set as an integer from `0` to `5`.

## Other Notes
- The standalone contains OpenACC directives that can compile and execute successfully for either CPU and GPU using `nvfortran`.
- `gfortran` can build the standalone with OpenACC offloading to the GPU; however, executing the binary results in the following error: `libgomp: cuStreamSynchronize error: an illegal memory access was encountered`
- When compiling with `ifort`, the stack size may have to be set to unlimited (ex: With a bash shell, run `ulimit -s unlimited`) so that the standalone does not produce a segmentation fault.
