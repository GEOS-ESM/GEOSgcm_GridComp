# pyMoist <> GEOS Fortran bindings

The `pyMoist` package is internally built to work as a standalone science code.

This module ties the Python/DSL version of `pyMoist` to the inner working of the larger GEOS model (Fortran side).

We use the `pyGEOSBridge` and insert init and run hooks in `GEOS_UW_InterfaceMod.F90`
