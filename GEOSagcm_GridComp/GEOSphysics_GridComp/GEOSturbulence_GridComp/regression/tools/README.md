# Tools to rebuild a Turbulence regression case from a MAPL2 capture

ctest does not use anything in this directory. These scripts are for the person who creates a case or regenerates its data (a new configuration, a new grid, a change of the TURBULENCE import or tracer list).

## Why they are needed

A case compares the MAPL3 TURBULENCE with a MAPL2 run, but a MAPL2 capture cannot be used as it is:

1. **Format.** The capture is in the old cubed-sphere layout (`lon`, `lat` with six faces stacked). MAPL3 reads the new layout (`Xdim`, `Ydim`, `nf`).
2. **Which component holds what.** MAPL2's TURBULENCE import state is one file. In the standalone case it is split over components: `DataForTurbulence` holds the 44 imports, `DataTracers` holds the 47 diffused tracers (MAPL2's `TR` bundle, the `TR_*` fields, renamed), `DataPhase2` holds the 9 imports Surface changes between the phases (`run2`). Each is restarted from its own file, so the capture has to be cut into those files.
3. **Names and metadata.** The data components are YAML files that declare every field with its standard name and units; they have to match the capture and `Turbulence_StateSpecs.rc`.
4. **Which exports to request.** What is requested changes what TURBULENCE computes (`ALLOC=` is emulated, see `../README.md` and `../../mapl3-porting-notes.md`, step 6), so the History list has to reproduce the MAPL2 run.

Doing this by hand is slow and easy to get wrong (47 tracers, 44 imports, renaming `.`), so it is scripted.

## Scripts

| Script | Does |
|---|---|
| `convert_capture.sh` | Converts the capture files to the new layout with `GEOS_Restart_Utilities/convertOldCStoNewCS.py` (no example file, so no `lons`/`lats`). |
| `find_phase2_fields.py` | Lists the imports that differ between "after Run1" and "before Run2": the `phase2_fields` of `run2/root.yaml` and the fields of `data-phase2.yaml`. |
| `gen_case_yaml.py` | From `Turbulence_StateSpecs.rc` and the converted capture, writes the fragments of `data-for-turbulence.yaml`, `data-tracers.yaml` and `root.yaml`, `history_fields.txt` (the exports to activate) and `tracer_names.txt` (capture name to MAPL3 name). |
| `gen_history.py` | Writes a case's `history.yaml` from `history_fields.txt`, less the exports you list, and (`run1`) the 9 that only Run2 writes. |
| `subset_nc.py` | Writes a restart file from a converted file: the grid variables and a chosen set of variables, optionally renamed. |
| `run_case_by_hand.sh` | Runs a case outside ctest, in a scratch directory, to compare with the capture (`cmpchk.py`) and to produce `checkpoints/last`. |

Python scripts need `numpy` and `netCDF4`.

## 1. The MAPL2 capture

Run GEOSgcm on `develop` (a MAPL2 build; the Turbulence sources must be the ones the port was made from, see `../../mapl3-porting-notes.md`, step 0) at the case's grid and configuration, with these keys in `AGCM.rc`:

```
COMPONENT_TO_RECORD: TURBULENCE
TEST_FRAMEWORK: .true.
TARGET_TIME: '2000-04-14 22:00:00'
MAPL_GridCapture: .true.
```

The run starts at `2000-04-14T21:00:00` with 1200 s steps; the capture is of the fourth step. It writes `TURBULENCE_{import,internal,export}_{before,after}_runPhase{1,2}`, 12 files with no `.nc` suffix, in the old layout. Keep `DT:` equal to the run step (refresh every step). Request nothing special in `HISTORY.rc`: the export capture holds the exports MAPL2 had allocated in a normal run, and that is what the case must reproduce.

Record every Turbulence key `AGCM.rc` sets: they go into the case's `turbulence.yaml`. Keys that are commented out are at their code default (some defaults depend on `LM`: see the case README).

## 2. Rebuilding a case

Let `CAP` be the directory with the 12 capture files, `NEW` the converted capture (for example `~/input/turbulence/<case>/C12L181/new-style`), `CASE` the case directory (for example `regression/turbulence-default-sa`), and `TOOLS` this directory. Copy an existing case to start from: `cp -r turbulence-shocmf-sa turbulence-<name>-sa`.

```
# 1. convert
GEOS_RESTART_UTILITIES=~/GEOS_Restart_Utilities  $TOOLS/convert_capture.sh $CAP $NEW

# 2. what changes between the phases (run2): the answer is the phase2_fields list
$TOOLS/find_phase2_fields.py $NEW

# 3. data components, root connections, exports
mkdir gen
$TOOLS/gen_case_yaml.py ../Turbulence_StateSpecs.rc $NEW gen
```

Compare `gen/data-for-turbulence.yaml` and `gen/data-tracers.yaml` with the case's files of the same names: they must be the same apart from the header comments and the `mapl:`, `misc:`, `import: TRDIFF:` and `items:` lines around them, which are in the case's files but not in the fragments (`gen/data-tracers.yaml` ends with the `items:` list, for `TRDIFF`). If the import or tracer list of the capture changed, edit the case's files (and `root.yaml`'s connections from `gen/connections_imports.yaml`; `data-phase2.yaml` and `phase2_fields` from step 2).

Next, `turbulence.yaml` (the keys of step 1), the comments at the top of the case's files, and `cmpchk.py` (`baseline_dir` is `$NEW`). Then the restart files, into `<regression data>/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSturbulence_GridComp/<case>/<run>/checkpoints/2000-04-14T22:00:00/` for `run1` and `run2`:

```
D=<regression data>/.../GEOSturbulence_GridComp/<case>/run1/checkpoints/2000-04-14T22:00:00
mkdir -p $D
# the imports, without the tracers
$TOOLS/subset_nc.py $NEW/TURBULENCE_import_before_runPhase1.nc $D/DataForTurbulence_internal.nc --all-but-prefix TR_
# the tracers: TR_<name> becomes <name>, '.' becomes '_' (gen/tracer_names.txt: "TR_<name> <name>")
$TOOLS/subset_nc.py $NEW/TURBULENCE_import_before_runPhase1.nc $D/DataTracers_internal.nc @gen/tracer_names.txt
# the internal state (all of it)
$TOOLS/subset_nc.py $NEW/TURBULENCE_internal_before_runPhase1.nc $D/TURBULENCE_internal.nc --all-but-prefix ZZZZ
```

Do the same for `run2` (the same three files), and for `run2` also `DataPhase2_internal.nc` from the import state *before Run2*, with the fields of step 2:

```
$TOOLS/subset_nc.py $NEW/TURBULENCE_import_before_runPhase2.nc $D2/DataPhase2_internal.nc EVAP QDEL QHAT S SDEL SH T TH USTAR
```

(`--all-but-prefix ZZZZ` copies every variable: no variable starts with `ZZZZ`.)

## 3. Which exports to request (`history.yaml`)

`gen/history_fields.txt` has every captured export that is a MAPL3 export. Request all of them first, then take out what must not be requested:

```
$TOOLS/gen_history.py $CASE/run1 run1 gen/history_fields.txt [EXPORT ...]
$TOOLS/gen_history.py $CASE/run2 run2 gen/history_fields.txt [EXPORT ...]
```

An export that this configuration never writes comes out uninitialized in MAPL3 (MAPL2 zeroed it), which shows as `nan` in `cmpchk.py`. Find them by running (next section), then list them as `EXPORT ...` and regenerate. Also take out any export the code reads unguarded when unrequested, for example `ISOTROPY` with SHOC off and `TKE` requested (a floating invalid trap, `../../mapl3-porting-notes.md`, step 9b). A requested export that fails to compare, with no `nan`, is a real difference to investigate, not to exclude.

## 4. Running and checking

```
export GEOS_BUILD=<fixture>/build-debug
export DATA_ROOT=<regression data>/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSturbulence_GridComp
export WORK_DIR=<scratch>
$TOOLS/run_case_by_hand.sh <case> run1
cd $WORK_DIR/<case>-run1 && python3 cmpchk.py     # "Variables that differ: none"
```

Same for `run2`. When both match, keep the output as the ctest baseline and check that it reproduces:

```
mkdir -p $DATA_ROOT/<case>/run1/checkpoints/last
cp $WORK_DIR/<case>-run1/checkpoints/last/*.nc $DATA_ROOT/<case>/run1/checkpoints/last/
$TOOLS/run_case_by_hand.sh <case> run1           # a second run
nccmp -dmgfsN -x lons,lats,corner_lons,corner_lats <baseline file> <second run's file>    # each .nc: identical
```

Add the case to `../cases.txt` and re-run CMake. ctest compares exactly with `nccmp` against `checkpoints/last`, and only for ifx builds (`../run_case.cmake`). Then update `../README.md` and `../../mapl3-porting-notes.md`, and copy the data to the shared regression-data tree and S3 (until then the cases pass without running: `../README.md`).

## Notes

- Build with the same environment as the fixture (the porting notes name the ifx modules); the case runs with `OMP_NUM_THREADS=1` on 6 tasks.
- `gen_case_yaml.py` asserts that every unconditional import is in the capture and that every captured export is a TURBULENCE export; an error there means the capture and `Turbulence_StateSpecs.rc` disagree (different source versions).
- `convertOldCStoNewCS.py` is in `GEOS_Restart_Utilities`, not in this repository.
