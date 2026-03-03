# ============================================================
# Land-cover classification (Random Forest)
# 
# Saves:
#   - outputs/tables/model_params_[WINDOW]_[STACK_TAG].csv
#   - outputs/tables/predictor_screening_[WINDOW]_[STACK_TAG].csv
#   - outputs/tables/variable_importance_[WINDOW]_[STACK_TAG].csv
#   - outputs/rasters/class_[WINDOW].tif
# ============================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(readr); library(yaml)
  library(glue);  library(purrr); library(tibble); library(stringr)
  library(ranger); library(tidyr)
})
set.seed(4401)

# ----- GDAL / terra temp + cache (adjust paths/sizes as needed)
dir.create("temp_tiles", showWarnings = FALSE, recursive = TRUE)
terraOptions(
  memfrac = 0.5,          # conservative
  progress = 1,
  todisk   = TRUE,
  tempdir  = "temp_tiles" # fast disk with free space
)
Sys.setenv(
  GDAL_CACHEMAX = "512",  # MB for GDAL block cache
  R_MAX_VSIZE   = "20Gb"  # Windows: allow larger virtual memory
)

# ----------------------- USER KNOBS -----------------------
WINDOW         <- "2022-2024" # 2004-2006; 2010-2012; 2016-2018; 2022-2024
STACK          <- "allbands"    # "allbands" or "fullextended"
COR_THRESH     <- 0.97
TREES          <- 1000
PRED_CORES     <- 1             
BATCH_EXTRACT  <- 50000L        # rows per extraction batch 
MAX_PER_CLASS  <- 150000L       # cap training rows per class 

# ----------------------- PATHS ----------------------------
paths        <- yaml::read_yaml("config/paths.yml")
REF_GRID     <- paths$reference_grid
AOI_PATH     <- paths$aoi

# ----------------------- DERIVED KEYS ---------------------
win_key <- switch(gsub("_","-",WINDOW),
                  "2004-2006"="04-06","2010-2012"="10-12","2016-2018"="16-18","2022-2024"="22-24",
                  stop("Bad WINDOW value"))
STACK_TAG <- switch(tolower(STACK),
                    "allbands"="all_bands","fullextended"="full_extended",
                    stop("STACK must be 'allbands' or 'fullextended'"))
COMPOSITE <- switch(tolower(STACK),
                    "allbands"     = glue("rasters/composites/AllBands_{win_key}.tif"),
                    "fullextended" = glue("rasters/composites/Composite_{win_key}_FullExtended.tif"))

TRAIN_GPKG_CANDIDATES <- c(
  glue("outputs/vectors/train_{WINDOW}.gpkg"),
  glue("outputs/vectors/train_{gsub('_','-',WINDOW)}.gpkg"),
  glue("outputs/vectors/train_{gsub('-','_',WINDOW)}.gpkg")
)

# ----------------------- OUTPUTS --------------------------
OUT_RASTERS  <- "outputs/rasters"
OUT_VECTORS  <- "outputs/vectors"
OUT_TABLES   <- "outputs/tables"
OUT_LOGS     <- "logs"
dir.create(OUT_RASTERS, TRUE, TRUE)
dir.create(OUT_VECTORS, TRUE, TRUE)
dir.create(OUT_TABLES,  TRUE, TRUE)
dir.create(OUT_LOGS,    TRUE, TRUE)

# ----------------------- CLASSES (Base-7) -----------------
classes <- tibble::tribble(
  ~code, ~name,
  1L, "Forest",
  2L, "Agrocanopy",
  3L, "Paddy",
  4L, "GrassScrub",
  5L, "Built",
  6L, "Water",
  7L, "Rock"
)

# ----------------------- LOAD AOI + REF + COMPOSITE -------
aoi     <- sf::st_read(AOI_PATH, quiet = TRUE) |> sf::st_make_valid()
ref     <- terra::rast(REF_GRID)
crs_ref <- terra::crs(ref)

x <- terra::rast(COMPOSITE)
x <- terra::mask(terra::crop(x, ref), ref)  # align to ref grid (streamed)
names(x) <- make.unique(make.names(names(x), allow_ = TRUE))
all_band_names <- names(x)
gc()

# ----------------------- PICK TRAIN LAYER -----------------
vec_path <- TRAIN_GPKG_CANDIDATES[file.exists(TRAIN_GPKG_CANDIDATES)][1]
if (is.na(vec_path)) {
  stop("Could not find a training file. Tried:\n- ",
       paste(TRAIN_GPKG_CANDIDATES, collapse = "\n- "))
}
li <- sf::st_layers(vec_path)
need_cols <- c("base7_code","split")
spatial_layers <- li$name[li$geomtype != "None" & li$features > 0]
if (!length(spatial_layers)) stop("No spatial layers found in ", vec_path)

layer_to_read <- NULL
for (ly in spatial_layers) {
  tmp <- try(suppressWarnings(sf::st_read(vec_path, layer = ly, quiet = TRUE, n_max = 1)),
             silent = TRUE)
  if (inherits(tmp, "try-error")) next
  if (all(need_cols %in% names(tmp))) { layer_to_read <- ly; break }
}
if (is.null(layer_to_read)) {
  layer_to_read <- spatial_layers[1]
  message("No layer explicitly matched required fields; using '", layer_to_read, "'.")
}

pts <- sf::st_read(vec_path, layer = layer_to_read, quiet = TRUE)
if (!all(need_cols %in% names(pts))) {
  stop("Layer '", layer_to_read, "' lacks required fields: ",
       paste(setdiff(need_cols, names(pts)), collapse = ", "))
}
pts <- sf::st_transform(pts, crs_ref)
pts <- suppressWarnings(sf::st_intersection(pts, aoi))

# ----------------------- EXTRACT PREDICTORS (BATCHED) -----
batch_extract <- function(rst, pts_sf, batch = 50000L) {
  n <- nrow(pts_sf)
  idx <- split(seq_len(n), ceiling(seq_len(n)/batch))
  out <- vector("list", length(idx))
  for (i in seq_along(idx)) {
    ii <- idx[[i]]
    v  <- terra::extract(rst, terra::vect(pts_sf[ii, , drop = FALSE]))
    out[[i]] <- dplyr::bind_cols(
      sf::st_drop_geometry(pts_sf[ii, , drop = FALSE]),
      tibble::as_tibble(v)[, -1, drop = FALSE]
    )
    rm(v); gc()
    cat(sprintf(" - extracted batch %d/%d (%d rows)\n", i, length(idx), length(ii)))
  }
  dplyr::bind_rows(out)
}

dat <- batch_extract(x, pts, batch = BATCH_EXTRACT) %>%
  dplyr::filter(!is.na(base7_code)) %>%
  dplyr::mutate(base7_code = factor(base7_code, levels = classes$code))

train_df <- dat %>% dplyr::filter(split == "train")
cat("n_train=", nrow(train_df), "\n", sep = "")
if (nrow(train_df) == 0) stop("No training rows with split=='train'.")
y_col <- "base7_code"

# ----------------------- MULTICOLLINEARITY SCREENING -------
cand <- base::intersect(names(train_df), all_band_names)
cand <- cand[vapply(train_df[, cand, drop = FALSE], is.numeric, TRUE)]
if (!length(cand)) stop("No usable numeric band columns were found in training data.")

prune_by_corr <- function(df, cols, thr = 0.95) {
  if (length(cols) < 2) {
    return(list(keep = cols,
                dropped = tibble::tibble(variable=character(),
                                         dropped_because_of=character(),
                                         abs_r=numeric())))
  }
  X <- as.data.frame(df[, cols, drop = FALSE])
  is_all_na <- vapply(X, function(v) all(is.na(v)), logical(1))
  is_const  <- vapply(X, function(v) stats::sd(v, na.rm = TRUE) == 0, logical(1))
  drop0     <- names(which(is_all_na | is_const))
  keep0     <- base::setdiff(cols, drop0)
  if (length(keep0) < 2) {
    return(list(keep=keep0,
                dropped=tibble::tibble(variable=drop0, dropped_because_of="all_NA_or_const", abs_r=NA_real_)))
  }
  Xk <- as.data.frame(df[, keep0, drop = FALSE])
  Xk <- Xk[, vapply(Xk, function(v) !all(is.na(v)), TRUE), drop = FALSE]
  if (ncol(Xk) < 2) {
    return(list(keep=colnames(Xk),
                dropped=tibble::tibble(variable=drop0, dropped_because_of="all_NA_or_const", abs_r=NA_real_)))
  }
  m <- try(stats::cor(Xk, use = "pairwise.complete.obs"), silent = TRUE)
  if (inherits(m, "try-error") || is.null(m) || ncol(m) < 2) {
    return(list(keep=colnames(Xk),
                dropped=tibble::tibble(variable=drop0, dropped_because_of="all_NA_or_const", abs_r=NA_real_)))
  }
  m[is.na(m)] <- 0
  keep    <- colnames(m)
  dropped <- tibble::tibble(variable=character(), dropped_because_of=character(), abs_r=numeric())
  repeat {
    A <- abs(m); A[lower.tri(A, diag = TRUE)] <- 0
    mx <- suppressWarnings(max(A))
    if (!is.finite(mx) || mx < thr) break
    ind   <- which(A == mx, arr.ind = TRUE)[1, ]
    v1    <- colnames(A)[ind[1]]
    v2    <- colnames(A)[ind[2]]
    mean1 <- mean(abs(m[v1, base::setdiff(colnames(m), v1)]))
    mean2 <- mean(abs(m[v2, base::setdiff(colnames(m), v2)]))
    drop  <- if (mean1 >= mean2) v1 else v2
    keep    <- base::setdiff(keep, drop)
    dropped <- dplyr::bind_rows(
      dropped,
      tibble::tibble(variable=drop,
                     dropped_because_of=if (drop == v1) v2 else v1,
                     abs_r=as.numeric(mx))
    )
    if (length(keep) < 2) break
    m <- m[keep, keep, drop = FALSE]
  }
  if (length(drop0)) {
    dropped <- dplyr::bind_rows(
      tibble::tibble(variable=drop0, dropped_because_of="all_NA_or_const", abs_r=NA_real_), dropped
    )
  }
  list(keep = keep, dropped = dropped)
}

mc     <- prune_by_corr(train_df, cand, thr = COR_THRESH)
x_cols <- mc$keep

pred_screen <- dplyr::bind_rows(
  tibble::tibble(variable = x_cols, status="kept",
                 dropped_because_of = NA_character_, abs_r = NA_real_),
  dplyr::mutate(mc$dropped, status = "dropped_by_corr")
) %>%
  dplyr::mutate(correlation_threshold = COR_THRESH) %>%
  dplyr::arrange(dplyr::desc(status), variable)

readr::write_csv(pred_screen,
                 file.path(OUT_TABLES, glue::glue("predictor_screening_{WINDOW}_{STACK_TAG}.csv"))
)

if (length(x_cols) == 0) {
  stop("After screening, no predictors remain. Raise COR_THRESH or verify band extraction/types.")
}
cat("Predictors from raster:", length(cand),
    " | kept after NA/const & corr screening:", length(x_cols), "\n")
gc()

# ----------------------- OPTIONAL CAP PER CLASS ------------
if (is.finite(MAX_PER_CLASS)) {
  set.seed(4401)
  train_df <- train_df %>%
    dplyr::select(base7_code, dplyr::all_of(x_cols)) %>%
    dplyr::mutate(.rand = runif(n())) %>%                 # randomize rows first
    dplyr::arrange(base7_code, .rand) %>%
    dplyr::group_by(base7_code) %>%
    dplyr::mutate(.rk = dplyr::row_number()) %>%          # rank within class
    dplyr::ungroup() %>%
    dplyr::filter(.rk <= MAX_PER_CLASS) %>%               # keep top-K per class
    dplyr::select(-.rand, -.rk)
  
  cat("n_train (capped)=", nrow(train_df), "\n", sep = "")
} else {
  train_df <- train_df %>%
    dplyr::select(base7_code, dplyr::all_of(x_cols))
}

# ----------------------- FIT FINAL MODEL (ranger lean) ----
mtry_val  <- max(1L, floor(sqrt(length(x_cols))))  # classic RF default
min_n_val <- 5L
rf_final <- ranger::ranger(
  dependent.variable.name = y_col,
  data           = train_df[, c(y_col, x_cols), drop = FALSE],
  num.trees      = TREES,
  mtry           = mtry_val,
  min.node.size  = min_n_val,
  classification = TRUE,
  probability    = FALSE,
  importance     = "impurity",
  write.forest   = TRUE,
  num.threads    = 1,
  save.memory    = TRUE
)

# Save model params + used bands
readr::write_csv(tibble::tibble(
  window     = gsub("_","-", WINDOW),
  algorithm  = "ranger_RF",
  trees      = TREES,
  mtry       = mtry_val,
  min_n      = min_n_val,
  stack      = STACK,
  stack_tag  = STACK_TAG,
  kept_bands = paste(x_cols, collapse = ","),
  all_bands  = paste(all_band_names, collapse = ",")
), file.path(OUT_TABLES, glue::glue("model_params_{WINDOW}_{STACK_TAG}.csv")))
gc()

# Variable importance
vi <- tibble::tibble(
  variable   = names(rf_final$variable.importance),
  importance = as.numeric(rf_final$variable.importance)
) %>% dplyr::arrange(dplyr::desc(importance))
readr::write_csv(vi,
                 file.path(OUT_TABLES, glue::glue("variable_importance_{WINDOW}_{STACK_TAG}.csv"))
)

# ----------------------- RASTER NA IMPUTE → DISK ----------
band_medians <- vapply(
  x_cols, function(b) stats::median(train_df[[b]], na.rm = TRUE), FUN.VALUE = numeric(1)
)
band_medians[!is.finite(band_medians)] <- 0

impute_dir <- file.path("logs", paste0("imputed_", WINDOW, "_", STACK_TAG))
dir.create(impute_dir, showWarnings = FALSE, recursive = TRUE)

imputed_paths <- character(length(x_cols))
for (i in seq_along(x_cols)) {
  b   <- x_cols[i]
  lyr <- x[[b]]            # streamed read
  med <- band_medians[[b]]
  lyr2 <- terra::ifel(is.na(lyr), med, lyr)
  outb <- file.path(impute_dir, paste0(b, ".tif"))
  terra::writeRaster(
    lyr2, outb, overwrite = TRUE,
    gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","PREDICTOR=2")
  )
  imputed_paths[i] <- outb
  rm(lyr, lyr2); gc()
}
x_kept <- terra::rast(imputed_paths)
names(x_kept) <- x_cols
rm(x); gc()

# ----------------------- RASTER PREDICTION → TILED --------
TILE_PX    <- 512   # smaller tiles = safer on 8 GB
OVERLAP_PX <- 0     # no overlap

tmp_dir <- file.path(OUT_RASTERS, paste0("tiles_", WINDOW, "_", STACK_TAG))
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# Build tile extents in MAP coordinates (uses x/y resolution)
make_tile_extents <- function(r, tile_px = 512, overlap_px = 0) {
  ex  <- terra::ext(r)
  dx  <- terra::xres(r); dy <- terra::yres(r)
  tw  <- tile_px * dx
  th  <- tile_px * dy
  ovx <- overlap_px * dx
  ovy <- overlap_px * dy
  
  xmin <- ex[1]; xmax <- ex[2]; ymin <- ex[3]; ymax <- ex[4]
  
  # Steps (subtract overlap; with overlap=0 this equals tw/th)
  x_stp <- max(1e-9, tw - ovx)
  y_stp <- max(1e-9, th - ovy)
  
  exts <- list(); k <- 0L
  y_top <- ymax
  while (y_top > ymin) {
    y_bot <- max(ymin, y_top - th)
    x_left <- xmin
    while (x_left < xmax) {
      x_right <- min(xmax, x_left + tw)
      k <- k + 1L
      exts[[k]] <- terra::ext(x_left, x_right, y_bot, y_top)
      x_left <- x_left + x_stp
    }
    y_top <- y_top - y_stp
  }
  exts
}

tile_exts <- make_tile_extents(x_kept, TILE_PX, OVERLAP_PX)

# Predict each tile directly to its own GeoTIFF (streamed)
tile_files <- character(length(tile_exts))
for (i in seq_along(tile_exts)) {
  ex <- tile_exts[[i]]
  
  # IMPORTANT: snap='in' to align to the pixel grid → avoids extent mismatch
  x_tile <- terra::crop(x_kept, ex, snap = "in")
  if (is.null(x_tile) || any(dim(x_tile)[1:2] == 0)) next
  
  fout <- file.path(tmp_dir, sprintf("pred_tile_%04d.tif", i))
  terra::predict(
    x_tile, rf_final, type = "response",
    filename = fout, overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER")),
    cores = 1, na.rm = FALSE, progress = 1
  )
  tile_files[i] <- fout
  rm(x_tile); gc()
}
tile_files <- tile_files[file.exists(tile_files) & nzchar(tile_files)]
stopifnot(length(tile_files) > 0)

# Mosaic tiles directly to disk (collection is a bit more forgiving)
tmp_pred <- file.path(OUT_RASTERS, glue::glue("class_{WINDOW}_{STACK_TAG}__tmp.tif"))
terra::mosaic(
  terra::sprc(tile_files),
  filename = tmp_pred, overwrite = TRUE,
  wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"))
)
gc()

# ----------------------- ALIGN TO REF_GRID + CLIP TO AOI → ONE OUTPUT ----------------
# Input: tmp_pred  (the mosaic you just wrote)
# Output: outputs/rasters/class_[WINDOW].tif  (aligned to REF_GRID + clipped to AOI)

pred_tmp <- terra::rast(tmp_pred)
aoi_v    <- terra::vect(aoi)   # AOI as SpatVector

# Align to REF grid (CRS + resolution + origin)
aligned_path <- file.path(OUT_RASTERS, glue::glue("class_{WINDOW}__aligned_tmp.tif"))
if (!terra::same.crs(pred_tmp, ref)) {
  pred_aligned <- terra::project(
    pred_tmp, ref, method = "near",
    filename = aligned_path, overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"))
  )
} else {
  pred_aligned <- terra::resample(
    pred_tmp, ref, method = "near",
    filename = aligned_path, overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"))
  )
}
rm(pred_tmp); gc()

# Crop to REF extent and clip to AOI
cropped_path <- file.path(OUT_RASTERS, glue::glue("class_{WINDOW}__crop_tmp.tif"))
pred_cropped <- terra::crop(
  pred_aligned, ref,
  filename = cropped_path, overwrite = TRUE,
  wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"))
)
rm(pred_aligned); gc()

# Final single file (no STACK_TAG in the name)
out_ras <- file.path(OUT_RASTERS, glue::glue("class_{WINDOW}.tif"))
terra::mask(
  pred_cropped, aoi_v,
  filename = out_ras, overwrite = TRUE,
  wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"))
)
rm(pred_cropped); gc()

# Cleanup intermediates (optional)
try(unlink(tmp_pred),      silent = TRUE)
try(unlink(aligned_path),  silent = TRUE)
try(unlink(cropped_path),  silent = TRUE)

cat(glue::glue("✅ wrote final aligned & AOI-clipped raster: {out_ras}\n"))


# Cleanup intermediates (optional)
try(unlink(tmp_pred), silent = TRUE)
try(unlink(tile_files,  force = TRUE), silent = TRUE)
try(unlink(tmp_dir,     recursive = TRUE, force = TRUE), silent = TRUE)
try(unlink(terra::tmpFiles(), force = TRUE), silent = TRUE)

cat(glue::glue("✅ wrote raster (tiled): {out_ras}\n"))


# ----------------------- LOG SESSION ----------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_LOGS, "sessionInfo.txt"))
cat("📝 logged R session to logs/sessionInfo.txt\n")
