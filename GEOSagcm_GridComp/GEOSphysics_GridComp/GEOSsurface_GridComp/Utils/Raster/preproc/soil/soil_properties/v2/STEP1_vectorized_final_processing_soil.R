library(dplyr)
library(RSQLite)
library(raster)
library(ncdf4)

# ----------------------------------------------------------------------
# HWSD2 Soil Processing (R)
#
# Inputs:
#   - HWSD2.bil (+ HWSD2.hdr, HWSD2.prj): global HWSD2 raster containing HWSD2_SMU_ID per pixel
#   - HWSD2.sqlite: SQLite version of the HWSD2 attribute database (converted from HWSD2.mdb)
#
# Example source files present in the HWSD2 package:
#   HWSD2.bil, HWSD2.hdr, HWSD2.prj, HWSD2.stx, HWSD2.mdb
#   HWSD2_RASTER.zip, HWSD2_DB.zip
#
# Run environment (Discover):
#   $ module load R/4.3.1
#   $ Rscript STEP1_vectorized_final_processing_soil.R
#
# Outputs:
#   - NetCDF tiles written as 10° x 10° files per layer (D1..D7), named like:
#       Lon_-180_-170_Lat_70_80_D1.nc
#   - Each tile contains 2-D variables on (lat, lon):
#       SOC, Sand, Silt, Clay, Coarse, Bulk_Density, Ref_Bulk_Density, PH,
#       Share, Porosity, Organic_Matter, Sequence_Used
#
# ----------------------------------------------------------------------
# DATA MODEL / ORIENTATION NOTES
#
# 1) SQLite database (HWSD2.sqlite)
#    - Contains soil attributes keyed by HWSD2_SMU_ID (mapping unit ID),
#      depth layer (D1..D7), and component SEQUENCE/SHARE.
#    - The database has NO spatial grid information (no lat/lon).
#
# 2) HWSD2 raster (HWSD2.bil + .hdr/.prj)
#    - Provides the spatial grid (lat/lon, resolution, and row/column layout).
#    - Each raster pixel stores a HWSD2_SMU_ID, which links to the SQLite attributes.
#
# 3) Raster -> matrices in R
#    - raster::as.matrix() returns [row, col] = [lat, lon] in raster row order
#      (typically row 1 = northernmost).
#    - We flip matrix rows once so that latitude is ascending (south -> north),
#
# 4) NetCDF output convention
#    - Variables are written with dimensions (lat, lon), e.g., Sand(lat, lon).
#    - If we ever transpose a matrix, it is ONLY to match the declared NetCDF
#      dimension order. No spatial meaning is "undone"—it is just array layout.
#
# ----------------------------------------------------------------------
# COMPONENT SELECTION POLICY (HWSD2_LAYERS table)
#
# For each HWSD2_SMU_ID and layer:
#   1) Drop any component with invalid values (any selected field < 0).
#   2) Select the remaining component with the largest SHARE (dominant component).
#
# Derived fields:
#   - OrganicMatter = min(SOC * 1.72, 100)  [Van Bemmelen factor]
#   - ParticleDensity: peat/volcanic/mineral heuristic
#   - Porosity = clamp(1 - BULK_DENSITY / ParticleDensity, 0..1)
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
raster_path <- "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/HWSD2.bil"
sqlite_path <- "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/HWSD2.sqlite"
# All output for all layers (D1..D7) is in this dir , for now I kept only D2 
output_directory <- "/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/output_STEP1/"

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
  cat("Created output directory:", output_directory, "\n")
}

# ----------------------------------------------------------------------
# Layers to process
# NOTE: depth_m is metadata here; this script writes 2-D (lat,lon) fields.
# ----------------------------------------------------------------------
depth_ranges <- list(
  "D1" = list(depth = "0-20cm",    depth_m = 0.2),
  "D2" = list(depth = "20-40cm",   depth_m = 0.2),
  "D3" = list(depth = "40-60cm",   depth_m = 0.2),
  "D4" = list(depth = "60-80cm",   depth_m = 0.2),
  "D5" = list(depth = "80-100cm",  depth_m = 0.2),
  "D6" = list(depth = "100-150cm", depth_m = 0.5),
  "D7" = list(depth = "150-200cm", depth_m = 0.5)
)

# ----------------------------------------------------------------------
# Load full HWSD2 raster once
# Raster stores HWSD2_SMU_ID per pixel.
# ----------------------------------------------------------------------
cat("Loading full HWSD2 raster:", raster_path, "\n")
raster_data <- raster(raster_path)
cat("Raster loaded.\n")

# Derive pixel resolution from raster rather than hardcoding
xdim <- res(raster_data)[1]
ydim <- res(raster_data)[2]

# ----------------------------------------------------------------------
# Helper: fetch + select one representative component per SMU for a layer
#
# Selection policy (IMPORTANT):
#   For each HWSD2_SMU_ID, HWSD2 may have multiple components (SEQUENCE).
#   We:
#     1) Drop any component with invalid soil properties (negative values).
#     2) Select the remaining component with the largest SHARE.
#
# Rationale:
#   - negative values are treated as invalid/missing
#   - SHARE selects the dominant component for each mapping unit
#
# Output:
#   One row per HWSD2_SMU_ID (plus NA rows for SMUs with no valid component)
# ----------------------------------------------------------------------
fetch_soil_data_layer <- function(db_conn, layer, smu_ids) {

  smu_ids <- unique(na.omit(smu_ids))
  if (length(smu_ids) == 0) return(NULL)

  # Force integer to avoid scientific notation in SQL IN() list
  smu_ids <- as.integer(smu_ids)
  smu_ids <- smu_ids[!is.na(smu_ids)]
  if (length(smu_ids) == 0) return(NULL)

  query <- paste0(
    "SELECT HWSD2_SMU_ID, SEQUENCE, SHARE, COARSE, SAND, SILT, CLAY, ",
    "BULK AS BULK_DENSITY, REF_BULK AS REF_BULK_DENSITY, ORG_CARBON AS SOC, PH_WATER ",
    "FROM HWSD2_LAYERS ",
    "WHERE LAYER = '", layer, "' ",
    "AND HWSD2_SMU_ID IN (", paste(smu_ids, collapse = ","), ")"
  )

  soil_data <- dbGetQuery(db_conn, query)
  if (nrow(soil_data) == 0) return(NULL)

  # QC + select dominant SHARE component
  soil_selected <- soil_data %>%
    mutate(HWSD2_SMU_ID = as.integer(HWSD2_SMU_ID)) %>%
    filter(
      COARSE >= 0, SAND >= 0, SILT >= 0, CLAY >= 0,
      BULK_DENSITY >= 0, REF_BULK_DENSITY >= 0,
      SOC >= 0, PH_WATER >= 0
    ) %>%
    arrange(HWSD2_SMU_ID, desc(SHARE)) %>%
    group_by(HWSD2_SMU_ID) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      # Organic matter from SOC using Van Bemmelen factor (1.72), clamp to 100%
      OrganicMatter = pmin(SOC * 1.72, 100),

      # Particle density heuristic (g/cm^3)
      ParticleDensity = case_when(
        OrganicMatter > 20 ~ 100 / ((OrganicMatter / 1.3) + ((100 - OrganicMatter) / 2.65)), # peat-like
        BULK_DENSITY < 1.0 & CLAY > 30 & SAND < 30 ~ 2.3,                                    # volcanic-like
        TRUE ~ 2.65                                                                           # mineral
      ),

      # Porosity derived from bulk density and particle density, clamp [0,1]
      Porosity = pmax(0, pmin(1, 1 - (BULK_DENSITY / ParticleDensity)))
    )

  # Ensure all requested SMUs appear (missing -> NA row)
  missing_smus <- setdiff(smu_ids, unique(soil_selected$HWSD2_SMU_ID))
  if (length(missing_smus) > 0) {
    soil_selected <- bind_rows(
      soil_selected,
      data.frame(
        HWSD2_SMU_ID = missing_smus,
        SEQUENCE = NA_integer_,
        SHARE = NA_real_,
        COARSE = NA_real_,
        SAND = NA_real_,
        SILT = NA_real_,
        CLAY = NA_real_,
        BULK_DENSITY = NA_real_,
        REF_BULK_DENSITY = NA_real_,
        SOC = NA_real_,
        PH_WATER = NA_real_,
        OrganicMatter = NA_real_,
        ParticleDensity = NA_real_,
        Porosity = NA_real_
      )
    )
  }

  # One row per SMU_ID (guaranteed)
  soil_selected %>% distinct(HWSD2_SMU_ID, .keep_all = TRUE)
}

# ----------------------------------------------------------------------
# Helper: vectorized mapping from raster SMU_ID values -> soil properties
#
# Instead of looping over SMU_IDs, we:
#   - take the raster cell SMU_ID vector v (length ncell)
#   - build idx = match(v, table$HWSD2_SMU_ID)
#   - index soil property vectors by idx (fast)
#
# IMPORTANT ABOUT ORIENTATION:
#   raster row 1 is typically north.
#   We output lat ascending (south->north), so after creating matrices
#   we flip rows once before writing NetCDF.
# ----------------------------------------------------------------------
vector_map_to_matrices <- function(raster_chunk, soil_tbl) {

  v <- getValues(raster_chunk)              # HWSD2_SMU_ID per cell (length ncell)
  v_int <- as.integer(v)                    # ensure integer IDs

  idx <- match(v_int, soil_tbl$HWSD2_SMU_ID) # NA where SMU missing/NA

  # Create vectors for every cell
  soc_vec    <- soil_tbl$SOC[idx]
  sand_vec   <- soil_tbl$SAND[idx]
  silt_vec   <- soil_tbl$SILT[idx]
  clay_vec   <- soil_tbl$CLAY[idx]
  coarse_vec <- soil_tbl$COARSE[idx]
  bd_vec     <- soil_tbl$BULK_DENSITY[idx]
  rbd_vec    <- soil_tbl$REF_BULK_DENSITY[idx]
  ph_vec     <- soil_tbl$PH_WATER[idx]
  share_vec  <- soil_tbl$SHARE[idx]
  por_vec    <- soil_tbl$Porosity[idx]
  om_vec     <- soil_tbl$OrganicMatter[idx]
  seq_vec    <- soil_tbl$SEQUENCE[idx]

  # reshape: assign cell values into a raster, then as.matrix()
  # This guarantees we follow raster's internal cell ordering.
  tmp <- raster(raster_chunk)

  to_mat <- function(vec) {
    rr <- tmp
    values(rr) <- vec
    as.matrix(rr)  # [row, col] => [lat, lon] in raster's row order (north->south)
  }

  mats <- list(
    SOC = to_mat(soc_vec),
    Sand = to_mat(sand_vec),
    Silt = to_mat(silt_vec),
    Clay = to_mat(clay_vec),
    Coarse = to_mat(coarse_vec),
    Bulk_Density = to_mat(bd_vec),
    Ref_Bulk_Density = to_mat(rbd_vec),
    PH = to_mat(ph_vec),
    Share = to_mat(share_vec),
    Porosity = to_mat(por_vec),
    Organic_Matter = to_mat(om_vec),
    Sequence_Used = to_mat(seq_vec)
  )

  mats
}

# ----------------------------------------------------------------------
# Process one 10° tile for one layer
# Outputs a NetCDF with dims (lat, lon) and lat ascending.
# ----------------------------------------------------------------------
process_chunk <- function(db_conn, layer, lon_min, lon_max, lat_min, lat_max) {

  cat("Tile Lon:", lon_min, lon_max, "Lat:", lat_min, lat_max, "Layer:", layer, "\n")

  chunk_extent <- extent(lon_min, lon_max, lat_min, lat_max)
  raster_chunk <- crop(raster_data, chunk_extent)

  if (is.null(raster_chunk) || ncell(raster_chunk) == 0) {
    cat("  -> Empty tile, skipping.\n")
    return(invisible(NULL))
  }

  ncols_chunk <- ncol(raster_chunk)
  nrows_chunk <- nrow(raster_chunk)

  # lon centers (ascending west->east)
  lon <- seq(from = xmin(raster_chunk) + xdim / 2, by = xdim, length.out = ncols_chunk)

  # lat centers: build ascending south->north
  # (raster row order is typically north->south; we will flip data rows later)
  lat <- seq(from = ymin(raster_chunk) + ydim / 2, by = ydim, length.out = nrows_chunk)

  # Create NetCDF dimensions in the same order as our matrices: (lat, lon)
  lat_dim <- ncdim_def("lat", "degrees_north", lat)
  lon_dim <- ncdim_def("lon", "degrees_east", lon)

  # Fill values 
  fill_f <- 1e15
  fill_i <- -9999L

  var_defs <- list(
    SOC = ncvar_def("SOC", "%", list(lon_dim, lat_dim), missval = fill_f, longname = "Soil Organic Carbon", compression = 4),
    Sand = ncvar_def("Sand", "%", list(lon_dim, lat_dim), missval = fill_f, longname = "Sand", compression = 4),
    Silt = ncvar_def("Silt", "%", list(lon_dim, lat_dim), missval = fill_f, longname = "Silt", compression = 4),
    Clay = ncvar_def("Clay", "%", list(lon_dim, lat_dim), missval = fill_f, longname = "Clay", compression = 4),
    Coarse = ncvar_def("Coarse", "%", list(lon_dim, lat_dim), missval = fill_f, longname = "Coarse Fragments", compression = 4),
    Bulk_Density = ncvar_def("Bulk_Density", "g cm-3", list(lon_dim, lat_dim), missval = fill_f, longname = "Bulk Density", compression = 4),
    Ref_Bulk_Density = ncvar_def("Ref_Bulk_Density", "g cm-3", list(lon_dim, lat_dim), missval = fill_f, longname = "Reference Bulk Density", compression = 4),
    PH = ncvar_def("PH", "-", list(lon_dim, lat_dim), missval = fill_f, longname = "pH in Water (-log(H+))", compression = 4),
    Share = ncvar_def("Share", "%", list(lon_dim, lat_dim), missval = fill_f, longname = "Selected Component Share", compression = 4),
    Porosity = ncvar_def("Porosity", "m3/m3", list(lon_dim, lat_dim), missval = fill_f, longname = "Porosity", compression = 4),
    Organic_Matter = ncvar_def("Organic_Matter", "%", list(lon_dim, lat_dim), missval = fill_f, longname = "Organic Matter Content", compression = 4),
    Sequence_Used = ncvar_def("Sequence_Used", "-", list(lon_dim, lat_dim), missval = fill_i, longname = "Sequence Used for Selected Component", compression = 4, prec = "integer")
  )

  output_file <- file.path(output_directory, paste0(
    "Lon_", lon_min, "_", lon_max, "_Lat_", lat_min, "_", lat_max, "_", layer, ".nc"
  ))
  if (file.exists(output_file)) file.remove(output_file)

  nc <- nc_create(output_file, var_defs)

  # Unique SMUs in this tile
  smu_ids <- unique(getValues(raster_chunk))
  smu_ids <- smu_ids[!is.na(smu_ids)]

  if (length(smu_ids) == 0) {
    cat("  -> No SMU IDs found; writing empty tile.\n")
    nc_close(nc)
    return(invisible(NULL))
  }

  # Fetch soil table for these SMUs and this layer
  soil_tbl <- fetch_soil_data_layer(db_conn, layer, smu_ids)
  if (is.null(soil_tbl) || nrow(soil_tbl) == 0) {
    cat("  -> No soil records returned; writing empty tile.\n")
    nc_close(nc)
    return(invisible(NULL))
  }

  # Vectorized mapping from raster values to matrices
  mats <- vector_map_to_matrices(raster_chunk, soil_tbl)

  # Flip rows so matrix rows correspond to lat ascending (south->north)
  for (nm in names(mats)) {
    mats[[nm]] <- mats[[nm]][nrow(mats[[nm]]):1, , drop = FALSE]
  }

  put_f <- function(name, mat) {
    mat[is.na(mat)] <- fill_f
    ncvar_put(nc, var_defs[[name]], t(mat))
  }
  put_i <- function(name, mat) {
    mat[is.na(mat)] <- fill_i
    ncvar_put(nc, var_defs[[name]], t(mat))
  }
  put_f("SOC", mats$SOC)
  put_f("Sand", mats$Sand)
  put_f("Silt", mats$Silt)
  put_f("Clay", mats$Clay)
  put_f("Coarse", mats$Coarse)
  put_f("Bulk_Density", mats$Bulk_Density)
  put_f("Ref_Bulk_Density", mats$Ref_Bulk_Density)
  put_f("PH", mats$PH)
  put_f("Share", mats$Share)
  put_f("Porosity", mats$Porosity)
  put_f("Organic_Matter", mats$Organic_Matter)
  put_i("Sequence_Used", mats$Sequence_Used)

  nc_close(nc)
  cat("  -> Wrote:", output_file, "\n")
  invisible(NULL)
}

# ----------------------------------------------------------------------
# Process all 10-degree tiles for one layer
# ----------------------------------------------------------------------
process_layer_chunks <- function(layer) {
  cat("\n=== Processing layer:", layer, "===\n")
  t0 <- Sys.time()

  db_conn <- dbConnect(SQLite(), sqlite_path)
  on.exit(dbDisconnect(db_conn), add = TRUE)

  lon_intervals <- seq(-180, 170, by = 10)
  lat_intervals <- seq(-90, 80, by = 10)

  for (lon_min in lon_intervals) {
    lon_max <- lon_min + 10
    for (lat_min in lat_intervals) {
      lat_max <- lat_min + 10
      tryCatch(
        process_chunk(db_conn, layer, lon_min, lon_max, lat_min, lat_max),
        error = function(e) {
          cat("ERROR tile Lon:", lon_min, lon_max, "Lat:", lat_min, lat_max, "Layer:", layer, "\n")
          cat("  ", e$message, "\n")
        }
      )
    }
  }

  cat("Layer", layer, "done. Time:", Sys.time() - t0, "\n")
}

# ----------------------------------------------------------------------
# Run all layers requested
# ----------------------------------------------------------------------
for (layer in names(depth_ranges)) {
  process_layer_chunks(layer)
}

