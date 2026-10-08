# Building LISF's `LIS` library under CMake

`LISF_cmake/CMakeLists.txt` builds LISF's `LIS` library natively under CMake,
replacing the previous approach of shelling out to LISF's own GNU-make build
via `make nuopc`.

## Checkout

```sh
git clone -b v11.10.2+LISFv7.8.0-public git@github.com:GEOS-ESM/GEOSgcm.git
cd GEOSgcm
mepo clone
```

`mepo clone` reads `components.yaml` and checks out all nested sub-repos,
including:

- `GEOSgcm_GridComp` on branch `feature/pchakrab/integrate-lis-into-v11` —
  all of our GEOS-side integration code lives in
  `@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/LIS_GridComp`.
- LISF itself, unmodified, from <https://github.com/NASA-LIS/LISF.git>, checked
  out under `LIS_GridComp/@LISF`.

## Build

```sh
source @env/g5_modules.sh # works with ifort stack, not ifx
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install -DBUILD_LIS=On
cmake --build build -j8
cmake --install build
```

`LIS_GridComp` is opt-in: without `-DBUILD_LIS=On` it is not added to the
build at all, and GEOS builds exactly as it did before.

The two libraries can also be built on their own:

```sh
make LIS            # libLIS.so  (~2000 sources, one target)
make LIS_GridComp   # libLIS_GridComp.so, the GEOS-side wrapper
```

Parallel (`-jN`) builds work from scratch. Verified with CMake 4.4.3 /
`ifort` 2021.13.0 / Baselibs 8.33.0 (`ifort-stack`), Unix Makefiles generator.

The stack matters: HDF4/HDF-EOS2 are enabled unconditionally (see below), and
only the `ifort` Baselibs ships HDF4's Fortran interface. Building against the
`ifx` stack currently fails in the HDF4-reading DA OBS plugins.

## How the source list is derived

We defer to LISF's own plugin selector rather than re-deriving its directory
list by hand:

1. `user.cfg`, `LIS_misc.h` and `LIS_NetCDF_inc.h` are written into
   `lisf_generated/` in the build tree, alongside a copy of LISF's
   `default.cfg`.
2. LISF's `plugins.py` runs there on every configure, reading
   `default.cfg` + `user.cfg` and emitting `Filepath` and `LIS_plugins.h`.
   (Re-running every time means edits to `user.cfg` always take effect.)
3. `Filepath`'s `dirs := . ../core ../plugins ...` line is parsed and each
   directory non-recursively globbed, mirroring `lis/make/Makefile`'s
   `FIND_FILES`/`FIND_HEADERS`. The entries are relative to `lis/make`,
   not to where `Filepath` was written.

Nothing is written into the `@LISF` clone, so it stays clean in
`mepo status`. LISF's own `configure` would instead drop these files
directly into `lis/make/`.

Some file names appear in more than one `Filepath` directory (e.g.
`get_cdf_params.F90` under both `metforcing/mogreps_g` and
`metforcing/galwem_ge`). The GNU Makefile's `VPATH` silently shadows later
duplicates and only builds the first match; we replicate that, otherwise both
copies get compiled and linked and the link fails with "multiple definition".

## Disabled plugins (`user.cfg`)

Everything else in LISF's `default.cfg` is left at its default (mostly On).

| Plugin | Why disabled |
| --- | --- |
| `VIC.4.1.1`, `VIC.4.1.2` | Restricted / unsupported. |
| `CABLE` | `cable_canopy.f90` fails to compile with array shape-mismatch errors (`rbw`/`poolcoef1*`). |
| `Noah.3.9` | `noah39_main.F90` calls `SFCDIF_OFF` with more actual than dummy arguments (vendored bug). We use Noah.3.3 / NoahMP.3.6 / NoahMP.4.0.1. |
| `RUC.3.7` | `LIS_lsm_pluginMod.F90` declares `external ruc37_reset` but no `RUC37_reset.F90` exists anywhere in LISF's `ruc.3.7` plugin (vendored gap), leaving `ruc37_reset_` undefined at link time. |

None of these are used by the Plug, which only drives Noah/NoahMP LSMs.

## Feature switches (`LIS_misc.h`)

One setting differs from the stock make-based build:

**`#define COUPLED`** — makes `plugins.py` exclude `lis/offline`, whose
`lisdrv.F90` is a `program` main entry point that would collide with GEOS.x's
own main at link time.

As a side effect, LISF's `noahmp36_wrf_routines.F90` only provides
`wrf_message`/`wrf_error_fatal` when `COUPLED` is *un*defined, yet NoahMP.3.6
and 4.0.1 still call them unconditionally. `GEOS_LIS_link_stubs.F90` supplies
those (and `grib_get`, referenced from an unguarded call site despite
`USE_GRIBAPI` being undef'd) so we don't have to patch vendored sources.

## HDF4 / HDF-EOS2

`USE_HDF4` and `USE_HDFEOS2` are both `#define`d, so LISF's HDF4-reading DA
OBS plugins are compiled for real rather than to no-ops.

Baselibs ships HDF4 and HDF-EOS2, but — unlike NetCDF/HDF5/ESMF —
`@cmake/external_libraries/FindBaselibs.cmake` creates no targets for them,
since MAPL itself does not need them. `CMakeLists.txt` therefore builds its
own `hdf4::hdf4` and `hdfeos::hdfeos` imported targets from `find_library`
results (`mfhdf`, `df`, `hdfeos`, `Gctp`, `jpeg`).

These have to stay `IMPORTED` rather than becoming an `INTERFACE` library:
`LIS` is installed via the `GEOSgcm-targets` export set, and an `INTERFACE`
library in its link interface fails generation with "requires target ... that
is not in any export set". Namespaced `IMPORTED` targets are exempt.

**These plugins need HDF4's Fortran interface** (`#include "hdf.f90"`), which
not every Baselibs build provides. Against a stack without it:

```
read_PMW_snow.F90(256): #error: can't find include file: hdf.f90
```

HDF4/HDF-EOS use in LISF is confined to 9 DA OBS readers. Three of them —
`AMSRE_SWE`, `MODISsca`, `NASA_AMSREsm` — have no entry in LISF's
`default.cfg`, so `plugins.py` never puts them on `Filepath` and they are not
built here (LISF's own make build does not build them either). The remaining
six are: `ANSA_SNWD`, `GCOMW_AMSR2L3SND`, `GLASS_Albedo`, `GLASS_LAI`,
`PMW_snow`, `SSMI_SNWD`.

If support for an HDF4-less stack is ever needed, gate the `find_library`
block and the two `USE_HDF*` defines on `EXISTS ${BASEDIR}/include/hdf/hdf.f90`
— the macros already guard every HDF4 call site, so the readers degrade to
no-ops cleanly.

## CMake Fortran dependency-scanner bug

`cmFortranParser` loses sync on a line-continuation `&` immediately followed
by a preprocessor directive — the NoahMP dummy-argument-list pattern:

```fortran
SUBROUTINE NOAHMP_SFLX (..., FLDFRC   &
#ifdef WRF_HYDRO
                       ,SFCHEADRT     &
#endif
#ifdef PARFLOW
                       ,QINSUR,ETRANI &
#endif
                       )
```

The scanner then treats the remainder of the file as a continuation, so every
`MODULE` defined after that point is silently dropped from the target's
`provides` list in `CMakeFiles/LIS.dir/fortran.internal`. With no provider
recorded, no ordering edge is emitted into `depend.make`, and parallel make
compiles the consumers first:

```
module_sf_noahmp_groundwater_401.F90(23): error #7002: Error in opening the
  compiled module file.  Check INCLUDE paths.   [NOAHMP_TABLES_401]
```

followed by a cascade of `#6404`/`#6580` errors for every symbol from the
missing module.

`lis_scanner_fixup.py` runs at configure time and restores the lost edges.
It scans every source for the trigger pattern, collects the `MODULE`s each
suspect file defines, matches them against every `USE` in the source list, and
emits one `<consumer>|<provider>,...` record per affected consumer;
`CMakeLists.txt` turns each into an `OBJECT_DEPENDS` on the providers' object
files.

The scan deliberately over-approximates — a file is suspect if it contains the
pattern at all, even though the scanner sometimes recovers — so some redundant
edges are emitted. They are harmless, and the alternative (a hand-maintained
table) silently rots: the original two-entry table missed `AC72_MODULE` and
eight other providers, which only surfaced as an intermittent `-jN` failure.

Because the scan re-runs on every configure, no manual step is needed after
bumping `@LISF`. To inspect what it found:

```sh
# providers it is compensating for
cut -d'|' -f2 <build>/.../LISF_cmake/lis_scanner_fixup.txt | tr ',' '\n' | sort -u
```

In the current plugin set this is 10 providers across 611 consumer records.

## Other build-system notes

**Linker driver.** Both `LIS` and `LIS_GridComp` set
`LINKER_LANGUAGE Fortran`. CMake would otherwise pick CXX (`icpx`) for `LIS`
because the source list includes some `.c`/`.cc` files, and ESMF's `esmf.mk`
link flags (`-threads`, `-cxxlib`, `-lifport`, `-lifcoremt`) are classic
ifort/icpc-only — `icpx` rejects them.

**Fortran-only compile options.** `-nomixed-str-len-arg`, `-names lowercase`,
`-convert big_endian` and `-assume byterecl` are gated behind
`$<COMPILE_LANGUAGE:Fortran>`; `icx` rejects them on the `.c` sources.
