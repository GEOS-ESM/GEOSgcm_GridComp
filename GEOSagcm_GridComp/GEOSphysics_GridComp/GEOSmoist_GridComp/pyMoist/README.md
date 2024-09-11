# pyMoist

`pyMoist` is the [NDSL](https://github.com/NOAA-GFDL/NDSL) version of the NASA GMAO's GEOS Moist physics and it's required interface to GEOS.

## Testing

The `tests/savepoint` folder contains the numerical regression testing that served to port the code from the original Fortran. It tests that the difference from original Fortran is minimal.

Note: bit to bit equivalence is not possible considering the targeted hardware and the change of compiled language (NDSL uses C).

## Interface

The `interface` sub-directory contain the three-way bridge to and from GEOS.

## Develop

- Run `pre-commit run --all-files` before comitting for code guidelines coherency.

This is a new line
