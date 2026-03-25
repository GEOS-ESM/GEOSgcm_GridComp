# pyMoist

`pyMoist` is the [NDSL](https://github.com/NOAA-GFDL/NDSL) version of the NASA GMAO's GEOS Moist physics and it's required interface to GEOS.

The status of the ported schemes goes as follows:

- GFDL One Moment Microphysics (GFDL 1M): [ported](https://github.com/FlorianDeconinck/GEOSgcm_GridComp/tree/develop/UW_and_GFDL1M/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/pyMoist/pyMoist/microphysics/GFDL_1M), [numerical regression done](https://github.com/FlorianDeconinck/GEOSgcm_GridComp/tree/develop/UW_and_GFDL1M/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/pyMoist/tests/translate_tests/microphysics/GFDL_1M), scientific validation on single column in progress.

- University of Washington Shallow Convection (UW): [ported](https://github.com/FlorianDeconinck/GEOSgcm_GridComp/tree/develop/UW_and_GFDL1M/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/pyMoist/pyMoist/convection/UW), [numerical regression done](https://github.com/FlorianDeconinck/GEOSgcm_GridComp/tree/develop/UW_and_GFDL1M/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/pyMoist/tests/translate_tests/convection/UW), scientific validationon single column in progress.

- Greil-Freitas Deep Convection (GF2020): porting done, integration ongoing.

We have no plans of porting the BACM_1M, MGB2_3M, RAS or THOM_1M schemes.

## Testing

The `tests/savepoint` folder contains the numerical regression testing that served to port the code from the original Fortran. It tests that the difference from original Fortran is minimal.

Note: bit to bit equivalence is not possible considering the targeted hardware and the change of compiled language (NDSL uses C).

## Interface to Fortran

The `fortran/param_interfaces` contains the hook to the MAPL-powered Fortan-Python bridge

## Contributing

Couple things to consider when contributing:

- Your branch name must contain the word `dsl`, lower case. somewhere. Good rule of thumb is to call them `dsl/xxxxxx`.
- Run `pre-commit run --all-files` before committing for code guidelines coherency.
- The convention for PR names is `[DSL] xxxx` to distinguish from other/general GEOS PRs.

