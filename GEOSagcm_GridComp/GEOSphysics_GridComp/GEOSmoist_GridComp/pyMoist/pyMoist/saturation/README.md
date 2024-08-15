# Saturation formulations

## In Fortran

### GEOS_QSatX

The code lists 3 ways, with each a specific set of formulation.

- `GEOS_QsatLqu` & `GEOS_QsatIce`:
    - Latest code, carrying for every subsequent code the "exact" formulations for the different saturations schemes
    - Can run with in exact mode or in table mode
- `Qsat` & `DQsat` are the "traditional" called, which are flagged as deprecated in docs but still in use
    - Can run in exact mode (and ping back to GEOS_QSatLqu/Ice) or in table mode 

The table (ESINIT) is a constant computation that leverages the GEOS_QsatLqu/GEOS_QSatIce to freeze results in increment
of 0.1 Pa (per documenttion)

A lot of the complexity of the code is due to micro-optimization to re-use inlined code instead of using function calls.

### MAPL_EQSatX

TBD

## Our implementation

### Estimated Saturation Table

:warning: This inplements the GEOS_QSatX only

The class `SaturationVaporPressureTable` in `table.py` computes on-demand the table of estimated saturation based on 0.1 Pa
increments. The `get_table` function in conjection with the `SaturationFormulation` in `fomulation.py` returns the correct table
to be sampled across a few fields.
