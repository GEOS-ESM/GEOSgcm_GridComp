# Topography Pipeline — Quickstart & Reference

This repository builds **GWD-ready topography** on cubed‑sphere grids (uniform and Schmidt‑stretched). The entry point for most users is **`make_topo.py`**, which interactively generates a Slurm job script and runs the full pipeline end‑to‑end via compiled utilities in `bin/`.

---

## Quickstart

```bash
# 1) Ensure binaries and modules are available
ls bin/{bin_to_cube.x,generate_scrip_cube_topo.x,cube_to_target.x,convert_to_gmao_output_topo.x}  # should exist

source g5_modules

# 2) Run the driver
./make_topo.py

# 3) Answer prompts (bin dir, output dir, input GMTED path, and resolutions)
#    The script writes a Slurm job like topo_{res}.j in your chosen out_dir.

# 4) Submit
cd <out_dir>
sbatch topo_c{res}.j
```

Outputs (per‑resolution) in `out_dir/output_<IM>/` as NetCDF files and binary alongside a GMAO restart.

---

## What the pipeline does

1. **Prepare a high‑res intermediate cube (`c3000`)**

   * Converts the global GMTED lat‑lon input
   * Produces: `c3000.gmted_fixedanarticasuperior.nc`

2. **Generate a SCRIP descriptor** (`generate_scrip_cube_topo.x`)

   * Uniform or Schmidt‑stretched (SG001/SG002) per `GenScrip.yaml` built by the job
   * Writes `PE<IM>x<JM>-CF.nc4` + `c<IM>_coords.nc4`

3. **Select `rrfac_max` (stretched grids only)** using `cdo` + `awk`

   * For **uniform grids**, `rrfac ≡ 1` everywhere, so `rrfac_max=1` and the flag is effectively a no‑op.

4. **Remap to target grid** (`cube_to_target.x`)

   * Uses the intermediate cube + `--smoothing_scale` + optional `--rrfac_max` and Laplacian iterations
   * Produces PE* NetCDF in `output_<IM>/`

5. **Convert to restart + GMAO outputs**

   * `scrip_to_restart_topo.py` (adds `-g sg001/sg002` as needed)
   * `convert_to_gmao_output_topo.x -i PE* --im <IM>`

---

## Inputs & dependencies

* **Topography source (lat‑lon):** `gmted_fixed_anartica_superior_caspian.nc4`
* **Binaries (in `bin/`):**

  * `bin_to_cube.x`
  * `generate_scrip_cube_topo.x`
  * `cube_to_target.x`
  * `convert_to_gmao_output_topo.x`
  * `scrip_to_restart_topo.py`
* **Environment**

  * `bin/g5_modules` (loads platform compilers/MPI)
  * Modules: `nco`, `cdo` , `python`
  * `mpirun` available in `$PATH`
* **Scheduler:** Slurm (`sbatch`)

> The driver writes a **csh** job file

---

## Running `make_topo.py`

`make_topo.py` prompts for:

* **bin_dir**: Absolute path to the `bin/` directory containing the compiled tools.
* **out_dir**: Where to put the generated Slurm script and outputs. Default: `/discover/nobackup/<user>/BCS_TOPO/`.
* **path_latlon**: Directory containing `gmted_fixed_anartica_superior_caspian.nc4`.
* **resolutions**: Pick from uniform (`C12 … C5760`) and stretched families (`SG001`, `SG002`).

  * If you choose a stretched family, you will be asked for its concrete IM sizes (e.g., `C270`, `C540`).

The script writes **`topo_<tag>.j`** in `out_dir` and sets sane defaults for time, nodes, and job name.

---

## Smoothing & GWD knobs

* **`smoothing_scale`** (per grid) lives in the Python map inside `make_topo.py`.

  * *Bigger* `smoothing_scale` ⇒ *stronger* smoothing ⇒ *smaller* GWD (less roughness).
* **`ALPHA`** (stretched only) scales the Laplacian smooth step. Injected into `GenScrip.yaml` when applicable.
* **`RRFAC`** (regional refinement factor) **exists only on stretched grids**. Uniform grids have `rrfac=1` everywhere; you can ignore any `rrfac_max` logic for uniform runs.
* **Extra Laplacian cycles**: For `IM=5760`, the job adds `-l 13` to `cube_to_target.x`.

### Defaults table

| Grid          | smoothing_scale | alpha (stretched only) |
| ------------- | --------------: | ---------------------: |
| C12           |           512.0 |                      – |
| C24           |           305.0 |                      – |
| C48           |           166.0 |                      – |
| C90           |            96.2 |                      – |
| C180          |           51.21 |                      – |
| C360          |           28.95 |                      – |
| C720          |            19.5 |                      – |
| C1120         |            8.26 |                      – |
| C1440         |            12.0 |                      – |
| C2880         |           3.285 |                      – |
| C5760         |             3.0 |                      – |
| C270 (SG001)  |           100.0 |                    2.9 |
| C540 (SG001)  |            53.3 |                   2.73 |
| C1080 (SG001) |            17.0 |                    7.0 |
| C1536 (SG002) |            26.8 |                   12.1 |
| C2160 (SG001) |            2.98 |                  14.25 |

> To add/retune a grid, edit the `smoothmap`/`alpha` dicts in `make_topo.py` and re‑run.

---

## Runtime guidance

> Wall estimates assume Discover nodes; adjust for your system.

* `c180 … c2880, c1536, c1120, c2160`: **1 node**, ~**1h**
* `c5760, c540, c270, c48`: **2 nodes**, ~**4h**
* `c90`: **1 node**, ~**3h**
* `c24`: **2 nodes**, ~**8h**
* `c12`: **2 nodes**, ~**19h**, use `qos=long`

> The driver currently sets 1 node in the Slurm header by default; edit the generated job for heavy cases.

---

## Outputs & where to find them

For each `IM`:

```
<out_dir>/
  topo_<tag>.j                 # the job script you submit
  output_<IM>/
    PE<IM>x<JM>* .nc           # remapped topography (latest picked by timestamp)
    gwd_internal_rst.*         # restart written by scrip_to_restart_topo.py
    *_gmao_*.nc                # GMAO‑formatted outputs from convert_to_gmao_output_topo.x
```

The job also leaves the descriptor in `<out_dir>/PE<IM>x<JM>-CF.nc4`.

---

## Stretching presets

* **SG001** ⇒ `TARGET_LON=-98.35`, `TARGET_LAT=39.5`, `STRETCH_FACTOR=2.5` and `IM ∈ {270,540,1080,2160}`
* **SG002** ⇒ `TARGET_LON=-98.35`, `TARGET_LAT=39.5`, `STRETCH_FACTOR=3.0` and `IM = 1536`

These are injected into `GenScrip.yaml` only for matching resolutions.
Users may add new stretching options with different target centers and stretch factors. Note, however, that the stretch factor must be less than 3.0 (3.0 is the maximum supported).

---

## How to add a new resolution

1. Edit `smoothmap` (and `alpha` if stretched) in `make_topo.py`.
2. If it’s a *stretched* grid, extend the SG001/SG002 lists in the job template section.
3. Re‑run `make_topo.py` and select the new resolution.

---

## Design notes

* The pipeline **builds `c3000` once** and **reuses** it across resolutions for consistency and speed.
* For **stretched** grids, `rrfac_max` is computed via `cdo infon` + `awk` on the SCRIP file to inform `cube_to_target.x`. 
* For **uniform** grids, `rrfac=1`, so `rrfac_max` is 1 and can be omitted.
* For `IM=5760`, the driver adds `-l 13` Laplacian cycles for additional smoothing.

---

## Example for one session

```
$ ./make_topo.py
Enter the root path of the bin: /home/me/topo/bin
Enter the path of the output directory: /discover/nobackup/me/BCS_TOPO
Enter the path contains gmted_fixed_anartica_superior_caspian.nc4 (confirm selection):
  /discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/topo/v1/
Select resolutions: [x] C360 [x] SG001
Select resolution of SG001 grid: [x] C540 [ ] C1080 [ ] C2160 [ ] C270
# => writes topo_c360_c540.j under out_dir
$ cd /discover/nobackup/me/BCS_TOPO
$ sbatch topo_c360_c540.j
# or submit one job at the time to have clean run log for each job: topo_c360.j. Save it is res dependent dir to have whole run isolated: /discover/nobackup/me/BCS_TOPO/c360/
```

---

## FAQ

**Q: Can I reuse an existing `c3000`?**
A: Yes. The job checks and reuses `c3000.gmted_fixedanarticasuperior.nc` if present.

**Q: Where do I change the stretch center/factor?**
A: Inside the job template section (look for `TARGET_LON/LAT` and `STRETCH_FACTOR`).

**Q: Do uniform grids use RRFAC?**
A: No. Uniform grids have `rrfac=1` everywhere. The `rrfac_max` step is only meaningful for stretched grids.

**Q: Can I tune GWD magnitude or general smoothing?**
A: Yes. smoothing_scale (all grids) and alpha (stretched only) control smoothing strength and GWD amplitude. Be cautious: aggressive changes can destabilize the Laplacian smoother. In such cases you may need to recode the smoother or carefully retune these knobs.

**Q: Why does this code differ from the NCAR_Topo repository?**
A: There are several reasons:
This repository started from an outdated fork of NCAR_Topo.
We introduced major improvements in job submission and processing, especially through dynamic segments.
We cannot push changes back upstream (no test platform there), so their repo may follow different design choices.

**Q: Can I widen or narrow the stretched‑grid refinement “dome”?**
A: Yes. The dome’s footprint is set by the half‑power radius (half_power_radius_deg) inside the stretched‑grid rrfac logic. By default it’s ~40°/sqrt(max(1, stretch_factor)), so for SF≈2.5 the half‑power radius is ~25° (covers most of CONUS). Increase the constant (e.g., 45–50) to broaden the dome; decrease to tighten it. This only affects stretched grids.

## List of Files in utils_topo

Below are the main files in this directory, with one-line summaries.  
Each file also has its own inline header with more detail.

- **CMakeLists.txt** — builds the utilities.
- **convert_bin_to_netcdf.F90** — small helper to convert raw binary topography into NetCDF (intermediate or testing).
- **convert_to_gmao_output.F90** — final step: produces GMAO deliverables (`gmted_DYN_ave_*.nc4`, `gmted_GWD_var_*.nc4`, `gmted_TRB_var_*.nc4`).
- **geompack.F90** — bundled geometry library (Burkardt routines: convex hull, triangle quality, etc.), needed for SCRIP generation.
- **generate_scrip_cube.F90** — builds SCRIP descriptors (uniform or Schmidt-stretched).
- **make_topo.py** — interactive driver, generates the Slurm job script to run the full pipeline.
- **scrip_to_cube_topo.py** — converts flat `ncol` → cube layout, using example file geometry.
- **scrip_to_restart_topo.py** — converts PE outputs into GWD restart format.

*Last updated: 2025-10-03*
