#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(terra); library(landscapemetrics) })

# ===================== USER INPUTS =====================
WINDOW <- "2016_2018"  # change as needed
IN_RASTERS <- list(
  `2004_2006` = "outputs/rasters/class_2004-2006.tif",
  `2010_2012` = "outputs/rasters/class_2010-2012.tif",
  `2016_2018` = "outputs/rasters/class_2016-2018.tif",
  `2022_2024` = "outputs/rasters/class_2022-2024.tif"
)
IN_PATH      <- IN_RASTERS[[WINDOW]]

FOREST_CODE  <- 1     # land-cover class value for forest
EDGE_WIDTH_M <- 60    # two 30 m pixels (use '>' to drop two rings)
MMU_HA       <- 5     # keep patches ≥ 5 ha

# AOI vector (recommended)
USE_AOI_VECTOR <- TRUE
AOI_PATH <- "vectors/aoi.gpkg"

DIR_RAS <- "outputs/rasters"
DIR_TMP <- "outputs/tmp"
DIR_TAB <- "outputs/tables"
dir.create(DIR_RAS, FALSE, TRUE)
dir.create(DIR_TMP, FALSE, TRUE)
dir.create(DIR_TAB, FALSE, TRUE)
terraOptions(progress = 3, memfrac = 0.5, todisk = TRUE, tempdir = DIR_TMP)

# ===================== READ INPUT ======================
lc <- terra::rast(IN_PATH)
if (is.na(terra::crs(lc, proj = TRUE)) || grepl("degree", terra::crs(lc, proj = TRUE), ignore.case = TRUE)) {
  stop("Input land-cover must be projected in meters (not degrees).")
}

# ===================== AOI (vector) =====================
aoi_v <- NULL
if (USE_AOI_VECTOR && file.exists(AOI_PATH)) {
  aoi_v <- try(terra::vect(AOI_PATH), silent = TRUE)
  if (!inherits(aoi_v, "SpatVector")) stop("Could not read AOI vector: ", AOI_PATH)
  
  # Ensure AOI has a CRS
  if (is.na(terra::crs(aoi_v, proj = TRUE)) || terra::crs(aoi_v, proj = TRUE) == "") {
    stop("AOI vector has no CRS. Define it before running (e.g., crs(aoi) <- 'EPSG:3763').")
  }
  
  # Compare CRS robustly (prefer terra::same.crs when available)
  same_crs <- FALSE
  if ("same.crs" %in% getNamespaceExports("terra")) {
    same_crs <- terra::same.crs(aoi_v, lc)
  } else {
    same_crs <- identical(terra::crs(aoi_v, proj = TRUE), terra::crs(lc, proj = TRUE))
  }
  if (!same_crs) aoi_v <- terra::project(aoi_v, terra::crs(lc))
  
  # Crop LC to AOI bbox (aligned to cell edges), then mask outside polygon to NA
  lc <- terra::crop(lc, aoi_v, snap = "out")
  lc <- terra::mask(lc, aoi_v)
}

# ===================== BUILD AOI MASK + TIGHT EXTENT  =====================

one <- lc; terra::values(one) <- 1
if (!is.null(aoi_v)) {
  aoi_mask <- terra::mask(one, aoi_v)   # NA outside AOI polygon
} else {
  aoi_mask <- terra::mask(one, lc)      # NA where lc is NA (no explicit AOI)
}
aoi_mask <- terra::mask(aoi_mask, lc)   # enforce lc NA

# Trim NA-only borders (no 'values='!)
aoi_ext <- terra::ext(terra::trim(aoi_mask))


# =====================  FOREST MASK  =====================
forest_mask <- terra::ifel(lc == FOREST_CODE, 1, 0)

# =====================  CORE (60 m) =====================
# Treat any non-forest (0) OR outside AOI (NA) as matrix targets
matrix_targets <- terra::ifel((forest_mask != 1) | is.na(lc), 1, NA)
d_m   <- terra::distance(matrix_targets)                    # Euclidean distance (meters)
core60 <- terra::ifel((forest_mask == 1) & (d_m > EDGE_WIDTH_M), 1, 0)
core60 <- terra::mask(core60, aoi_mask, updatevalue = NA)

# =====================  MMU ≥ 5 ha  =====================
cell_area_m2 <- prod(terra::res(lc))
min_cells <- ceiling((MMU_HA * 10000) / cell_area_m2)

core_1na <- terra::classify(core60, rcl = cbind(1,1), others = NA)  # 1=core, NA=else
lab_path <- file.path(DIR_TMP, paste0("labels_", WINDOW, ".tif"))
if (file.exists(lab_path)) unlink(lab_path, force = TRUE)
lab <- terra::patches(core_1na, directions = 8, zeroAsNA = TRUE,
                      filename = lab_path, overwrite = TRUE)

out_final <- file.path(DIR_RAS, paste0("forest_core_", WINDOW, ".tif"))

write_core <- function(r_bin, out_path) {
  # Preserve NA outside AOI and crop to tight AOI extent for alignment
  r_out <- terra::mask(r_bin, aoi_mask, updatevalue = NA)
  r_out <- terra::crop(r_out, aoi_ext, snap = "out")
  if (file.exists(out_path)) unlink(out_path, force = TRUE)
  terra::writeRaster(
    r_out, out_path, overwrite = TRUE,
    datatype = "INT1U", NAflag = 255,
    wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"))
  )
  invisible(r_out)
}

if (is.null(lab)) {
  final0 <- terra::classify(core60, rcl = cbind(1,0), others = 0)
  r_core_written <- write_core(final0, out_final)
  cat(sprintf("No core pixels. Saved: %s\n", out_final))
} else {
  ft <- suppressWarnings(terra::freq(lab)) |> as.data.frame()
  ft <- ft[!is.na(ft$value), , drop = FALSE]
  keep_ids <- if (nrow(ft)) ft$value[ft$count >= min_cells] else integer(0)
  
  kept_1na <- if (length(keep_ids)) terra::classify(lab, rcl = cbind(keep_ids, 1), others = NA) else lab * NA
  final_bin <- terra::ifel(!is.na(kept_1na), 1, 0)
  
  r_core_written <- write_core(final_bin, out_final)
  
  n_core <- suppressWarnings(as.numeric(terra::global(r_core_written == 1, "sum", na.rm = TRUE)[1,1]))
  if (is.na(n_core)) n_core <- 0
  cat(sprintf("Saved final core: %s\nKept core pixels (≥ %.1f ha; ≥ %d cells): %s\n",
              out_final, MMU_HA, min_cells, n_core))
}

# Optional: clean temp label file
try(unlink(lab_path, force = TRUE), silent = TRUE)

# ===================== WRITE FRAG METRICS CSV =====================
OUT_CSV <- file.path(DIR_TAB, paste0("frag_metrics_", WINDOW, ".csv"))

r_core <- terra::rast(out_final)
r_core <- terra::classify(r_core, rcl = matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))

cell_area_m2 <- prod(terra::res(r_core))
n_core_cells <- suppressWarnings(as.numeric(terra::global(r_core == 1, "sum", na.rm = TRUE)))
if (is.na(n_core_cells)) n_core_cells <- 0
core_area_ha <- (n_core_cells * cell_area_m2) / 10000

if (n_core_cells == 0) {
  out <- data.frame(
    core_area_ha = 0,
    edge_density_m_per_ha = 0,
    mean_patch_area_ha = NA_real_,
    mean_nn_dist_m = NA_real_
  )
  write.csv(out, OUT_CSV, row.names = FALSE)
  cat(sprintf("✓ Saved metrics to %s (no core present)\n", OUT_CSV))
} else {
  r_metrics <- terra::classify(r_core, rcl = matrix(c(0, NA, 1, 1), ncol = 2, byrow = TRUE))
  
  pull_class1 <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NA_real_)
    idx <- which(df$class == 1)
    if (length(idx) == 0) return(NA_real_)
    as.numeric(df$value[idx[1]])
  }
  
  ed_tbl  <- try(lsm_c_ed(r_core, directions = 8, count_boundary = TRUE), silent = TRUE)
  edge_density_m_per_ha <- if (!inherits(ed_tbl, "try-error")) pull_class1(ed_tbl) else NA_real_
  
  if (is.na(edge_density_m_per_ha) || edge_density_m_per_ha == 0) {
    te_tbl <- try(lsm_c_te(r_core, directions = 8, count_boundary = TRUE), silent = TRUE)
    if (!inherits(te_tbl, "try-error")) {
      te_m <- pull_class1(te_tbl)
      land_cells <- suppressWarnings(as.numeric(terra::global(!is.na(r_core), "sum", na.rm = TRUE)))
      land_area_ha <- (land_cells * cell_area_m2) / 10000
      if (!is.na(te_m) && land_area_ha > 0) edge_density_m_per_ha <- te_m / land_area_ha
    }
  }
  
  area_mn_tbl <- try(lsm_c_area_mn(r_metrics), silent = TRUE)
  mean_patch_area_ha <- if (!inherits(area_mn_tbl, "try-error")) pull_class1(area_mn_tbl) else NA_real_
  
  enn_mn_tbl <- try(lsm_c_enn_mn(r_metrics, directions = 8), silent = TRUE)
  mean_nn_dist_m <- if (!inherits(enn_mn_tbl, "try-error")) pull_class1(enn_mn_tbl) else NA_real_
  
  out <- data.frame(
    core_area_ha = round(core_area_ha, 2),
    edge_density_m_per_ha = round(edge_density_m_per_ha, 6),
    mean_patch_area_ha = round(mean_patch_area_ha, 6),
    mean_nn_dist_m = round(mean_nn_dist_m, 3)
  )
  write.csv(out, OUT_CSV, row.names = FALSE)
  cat(sprintf("✓ Saved metrics to %s\n", OUT_CSV))
}

# ===================== CLEAN TMP FOLDER =====================
tmp_items <- list.files(DIR_TMP, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
if (length(tmp_items)) {
  try(unlink(tmp_items, recursive = TRUE, force = TRUE), silent = TRUE)
  cat(sprintf("✓ Cleared temporary folder: %s\n", DIR_TMP))
}
