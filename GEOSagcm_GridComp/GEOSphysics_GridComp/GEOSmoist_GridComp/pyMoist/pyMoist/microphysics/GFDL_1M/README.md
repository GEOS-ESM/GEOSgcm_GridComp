# GFDL One Moment Microphysics - GEOS Flavored

This is DSL port of the GFDL_1M scheme called from `GEOS_GFDL_1M_InterfaceMod.F90`.

Last version ported: v11.5.2

Numerical regression tests can be found in `./tests/translate_tests/GFDL_1M` with a script to run them in `./tests/scripts/run_tests`.

## Porting status

Last update: 2026/02/12

### Scientific validation

SCM runs TBD.

### Numerical validation

All code has been ported. Last check on numerical tests:

✅: Validates. 🟢: Validates with threshold.  🔴: Doesn't validate. 🟣: Crashes.

| Test               |                |      `debug`      |              `numpy`              | `dace:cpu_KJI` |
| ------------------ | -------------- | :---------------: | :-------------------------------: | :------------: |
| Setup              |                |        🔴         |         🟣 \| Needs .at()         |       ✅        |
| PhaseChange        |                | 🟣 \| crash > esx |         🟣 \| Needs .at()         |       🟢       |
|                    | RHCalculations |        🔴         |         🟣 \| Needs .at()         |       ✅        |
|                    | HydrostaticPDF | 🟣 \| crash > esx | 🟣 \| Data indices must be scalar |       🔴       |
|                    | MeltFreeze     |         ✅         |                 ✅                 |       ✅        |
|                    | Evaporate      |        🔴         |                🔴                 |       ✅        |
|                    | Sublimate      |        🔴         |                🔴                 |       ✅        |
| Driver             |                |        🔴         | 🟣 \| Data indices must be scalar |       🔴       |
|                    | DriverSetup    |         ✅         |                 ✅                 |       ✅        |
|                    | FallSpeed      |         ✅         |                 ✅                 |       ✅        |
|                    | TerminalFall   |         ✅         |                 ✅                 |       ✅        |
|                    | WarmRain       |        🔴         | 🟣 \| Data indices must be scalar |       🔴       |
|                    | IceCloud       |         ✅         | 🟣 \| Data indices must be scalar |       ✅        |
|                    | DriverFinish   |         ✅         |                 ✅                 |       ✅        |
|                    | ---            |                   |                                   |                |
|                    | driver_tables  |   🟣 \| No data   |           🟣 \| No data           | 🟣 \| No data  |
| RedistributeClouds |                |         ✅         |                                   |       ✅        |
| RadiationCoupling  |                |        🔴         | 🟣 \| Data indices must be scalar |       ✅        |
| Finalize           |                |        🔴         | 🟣 \| Data indices must be scalar |       🔴       |
