# pyMoist

`pyMoist` is the [NDSL](https://github.com/NOAA-GFDL/NDSL) version of the NASA GMAO's GEOS Moist physics and it's required interface to GEOS.

|  *Model code*   |                  |   Numerical    | Integration | Scientific |      Orchestration      |                 CPU                 | GPU |
| :-------------: | :--------------: | :------------: | :---------: | :--------: | :---------------------: | :---------------------------------: | --- |
|  GFDL_1M  |                  |       🟢       |     🟢      |            |            ❌ (doesn't validate)            |                                     |     |    |
|                 |      Driver      |       🟢       |     N/A     |    N/A     |           🟢            |   📈 1.2 (c24)   |     |
|    UW     |                  |       🟢       |     🟢      |     🟢     |           🟢            | 📉 *0.67 (c24)* |     |     |
|    GF     |                  |       *Ongoing*       |             |            |                         |                                     |     |     |
| Chemistry PCHEM |                  |       *Ongoing*       |             |            |                         |                                     |     |

Legend:

- Numerical: numerical regression vs Fortran on a C24 resolution. See testing.
- Integration: integrating the NDSL code (python) in the current GEOS (fortran), either via direct C-binding and/or via an experimental MAPLpy layer
- Scientific: running 7 days at C180 on Discover/PRISM with expert review
- CPU: speed up (>1) or slow down (<1) versus the Fortran, with the resolution noted
- GPU: idem CPU but on the device, with a "generation to generation" comparison

## Testing

The `tests/savepoint` folder contains the numerical regression testing that served to port the code from the original Fortran. It tests that the difference from original Fortran is minimal.

Note: bit to bit equivalence is not possible considering the targeted hardware and the change of compiled language (NDSL uses C).

## Interface

The `interface` sub-directory contain the three-way bridge to and from GEOS.

## Contributing

Couple things to consider when contributing:

- Your branch name must contain the word `dsl`, lower case. somewhere. Good rule of thumb is to call them `dsl/xxxxxx`.
- Run `pre-commit run --all-files` before committing for code guidelines coherency.
- The convention for PR names is `[NDSL] xxxx` to distinguish from other/general GEOS PRs.
