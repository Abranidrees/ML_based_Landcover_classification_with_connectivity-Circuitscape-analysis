#!/usr/bin/env Rscript
# Frag summary figure — ultra-stable + low-memory
# ------------------------------------------------------------------

suppressPackageStartupMessages({ library(terra) })

# ========= USER SETTINGS =========
OUT_PNG    <- "outputs/figures/frag_summary.png"
CORE_PATHS <- list(
  `2004_2006` = "outputs/rasters/forest_core_2004_2006.tif",
  `2010_2012` = "outputs/rasters/forest_core_2010_2012.tif",
  `2016_2018` = "outputs/rasters/forest_core_2016_2018.tif",
  `2022_2024` = "outputs/rasters/forest_core_2022_2024.tif"
)

# Area bins (ha) and legend (edit to taste)
BIN_BREAKS_HA <- c(5, 10, 50, 200, Inf)
BIN_LABELS    <- c("5–10 ha", "10–50 ha", "50–200 ha", "≥200 ha")
BIN_COLORS    <- c("#c6dbef", "#6baed6", "#3182bd", "#08519c")  # sequential blues

# Downsample factor for FIGURE ONLY 
# Try 2 (mild), 4 (good), 8 (very safe). Larger = coarser visualization.
AGG_FAC <- 4

# Temp dir + terra options
DIR_TMP <- "outputs/tmp"
dir.create(DIR_TMP, FALSE, TRUE)
terraOptions(tempdir = DIR_TMP, todisk = TRUE, memfrac = 0.30, progress = 1)

# ========= HELPER: safe temp file path =========
.tmp_file <- function(prefix) {
  file.path(DIR_TMP, paste0(prefix, "_", sprintf("%09d", as.integer(runif(1, 1, 1e9))), ".tif"))
}

# ========= HELPER: compute binned raster for ONE panel  =========
compute_bins_to_disk <- function(core_path) {
  if (!file.exists(core_path)) stop("Missing core raster: ", core_path)
  
  # read & (optionally) downsample for plotting only
  r0 <- rast(core_path)
  # Ensure core=1, background=0; then convert 0 -> NA after aggregate
  # (aggregate 'max' keeps any core presence within the block)
  r_aggr_path <- .tmp_file("core_aggr")
  r_aggr <- aggregate(r0, fact = AGG_FAC, fun = "max",
                      filename = r_aggr_path, overwrite = TRUE,
                      wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER")))
  
  # keep core only (0 -> NA), write to disk
  core1_path <- .tmp_file("core1")
  core1 <- classify(r_aggr,
                    rcl = matrix(c(0, NA, 1, 1), ncol = 2, byrow = TRUE),
                    filename = core1_path, overwrite = TRUE,
                    wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"),
                                datatype = "INT1U", NAflag = 255))
  
  # label patches (8-neighbour) to disk — the memory heavy step
  lab_path <- .tmp_file("labels")
  lab <- patches(core1, directions = 8, zeroAsNA = TRUE,
                 filename = lab_path, overwrite = TRUE,
                 wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"),
                             datatype = "INT32"))  # labels need INT32
  
  # freq -> area (ha)
  cell_area_m2 <- prod(res(lab))
  ft <- as.data.frame(freq(lab))
  ft <- ft[!is.na(ft$value), , drop = FALSE]
  
  # no core present
  if (!nrow(ft)) {
    bins_path <- .tmp_file("bins_empty")
    bins <- lab * NA
    writeRaster(bins, bins_path, overwrite = TRUE,
                wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"),
                            datatype = "INT1U", NAflag = 255))
    bins <- rast(bins_path)
    bins <- as.factor(bins)
    levels(bins) <- data.frame(value = seq_along(BIN_LABELS), label = BIN_LABELS)
    # cleanup big intermediates
    try(unlink(c(r_aggr_path, core1_path, lab_path), force = TRUE), silent = TRUE)
    gc()
    return(bins)
  }
  
  ft$area_ha <- ft$count * cell_area_m2 / 10000
  
  # assign bin index per label
  bin_idx <- findInterval(ft$area_ha, BIN_BREAKS_HA,
                          left.open = FALSE, rightmost.closed = TRUE)
  keep <- !is.na(bin_idx)
  rcl  <- cbind(ft$value[keep], bin_idx[keep])
  
  # labels -> bin codes, write to disk
  bins_path <- .tmp_file("bins")
  bins <- classify(lab, rcl = rcl, others = NA,
                   filename = bins_path, overwrite = TRUE,
                   wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"),
                               datatype = "INT1U", NAflag = 255))
  bins <- as.factor(bins)
  levels(bins) <- data.frame(value = seq_along(BIN_LABELS), label = BIN_LABELS)
  
  # cleanup big intermediates 
  try(unlink(c(r_aggr_path, core1_path, lab_path), force = TRUE), silent = TRUE)
  gc()
  bins
}

# =========  vectorize =========
plot_panel_vector <- function(core_path, main_title) {
  r <- rast(core_path)
  if (AGG_FAC > 1) {
    r_aggr_path <- .tmp_file("core_aggr_vec")
    r <- aggregate(r, fact = AGG_FAC, fun = "max",
                   filename = r_aggr_path, overwrite = TRUE,
                   wopt = list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER")))
  }
  r <- classify(r, rcl = matrix(c(0, NA, 1, 1), ncol = 2, byrow = TRUE))  # 0->NA
  # dissolve contiguous core pixels into polygons (touches=TRUE merges corner-touching)
  poly <- as.polygons(r, dissolve = TRUE, trunc = FALSE, touches = TRUE)
  if (is.null(poly) || nrow(poly) == 0) {
    plot.new(); title(main_title, line = 0.2, cex.main = 1)
    return(invisible())
  }
  poly$area_ha <- expanse(poly, unit = "ha")
  # bin by area
  # Build breaks that exactly mirror BIN_BREAKS_HA (last is Inf already)
  brks <- BIN_BREAKS_HA
  # labels via cut; right=FALSE means intervals [a,b)
  poly$bin <- cut(poly$area_ha, breaks = brks, labels = BIN_LABELS, right = FALSE, include.lowest = TRUE)
  
  # Plot polygons (no legend here; shared legend later)
  plot(poly["bin"], axes = FALSE, legend = FALSE, border = NA,
       col = BIN_COLORS, mar = c(2,2,2,2))
  title(main_title, line = 0.2, cex.main = 1)
}

# ========= DRAW FIGURE  =========
dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)
png(OUT_PNG, width = 2000, height = 1600, res = 300)
op <- par(no.readonly = TRUE)
layout(matrix(c(1,2,3,4,5,5), nrow = 3, byrow = TRUE), heights = c(1,1,0.28))

for (nm in names(CORE_PATHS)) {
  par(mar = c(2,2,2,2))
  core_path <- CORE_PATHS[[nm]]
  if (!file.exists(core_path)) {
    plot.new(); title(paste0(nm, " (missing)"), line = 0.2, cex.main = 1)
    next
  }
  
  # Try raster-based bins; on OOM fall back to vector plotting
  ok <- TRUE
  bins <- try(compute_bins_to_disk(core_path), silent = TRUE)
  if (inherits(bins, "try-error")) ok <- FALSE
  
  if (ok) {
    # Plot the binned raster
    plot(bins, col = BIN_COLORS, legend = FALSE, axes = FALSE, box = FALSE, asp = NA)
    title(nm, line = 0.2, cex.main = 1)
    # cleanup the bins on disk immediately
    try({
      srcs <- sources(bins)
      rm(bins); gc()
      if (length(srcs)) for (s in srcs) if (!is.na(s) && file.exists(s)) unlink(s, force = TRUE)
    }, silent = TRUE)
  } else {
    # Fallback vector rendering
    message("Raster labeling failed for ", nm, " – falling back to vector plotting.")
    plot_panel_vector(core_path, nm)
  }
}

# Shared legend
par(mar = c(1,1,1,1))
plot.new()
legend("center",
       legend = BIN_LABELS,
       fill   = BIN_COLORS,
       border = NA, bty = "n", cex = 0.95,
       title  = "Core patch area (ha)")

par(op); dev.off()
cat(sprintf("✓ Wrote %s\n", OUT_PNG))

# Final tidy of temp folder contents (keeps folder)
tmp_items <- list.files(DIR_TMP, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
if (length(tmp_items)) try(unlink(tmp_items, recursive = TRUE, force = TRUE), silent = TRUE)
