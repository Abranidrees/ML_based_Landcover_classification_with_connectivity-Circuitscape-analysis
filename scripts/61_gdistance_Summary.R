#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(glue)
})

# ================================================================
# USER SETTINGS
# ================================================================
species <- "PFL"      # choose: "PFL" or "GEN"
proj_root <- getwd()

cat("\n=====================================================\n")
cat(glue("     LCP SUMMARY — SPECIES = {species} (LATE ONLY)\n"))
cat("=====================================================\n\n")

# ================================================================
# LOAD ACC_COST (LATE ONLY)
# ================================================================
acc_path <- file.path(
  proj_root, "outputs/rasters",
  glue("acc_cost_{species}_Late.tif")
)

cat(glue("📌 Loading raster:\n   {acc_path}\n"))
acc_ras <- rast(acc_path)
acc_crs <- st_crs(crs(acc_ras))
cat("✔ Raster CRS loaded.\n\n")

# ================================================================
# LOAD LCP LINES
# ================================================================
lcp_path <- file.path(
  proj_root, "outputs/vectors",
  glue("LCP_{species}_Late.gpkg")
)

cat(glue("📌 Loading LCP vectors:\n   {lcp_path}\n"))
lcp_sf <- st_read(lcp_path, quiet = TRUE)
lcp_crs <- st_crs(lcp_sf)
cat("✔ LCP CRS loaded.\n\n")

# ================================================================
# FIX CRS IF MISMATCH
# ================================================================
cat("📌 Checking CRS match...\n")

if (acc_crs$wkt != lcp_crs$wkt) {
  cat("⚠️ CRS mismatch — reprojecting LCP to raster CRS...\n")
  lcp_sf <- st_transform(lcp_sf, acc_crs)
  cat("   ✔ Reprojection complete.\n\n")
} else {
  cat("✔ CRS matches.\n\n")
}

# ================================================================
# GEOMETRY CHECK
# ================================================================
cat("📌 Checking geometry validity...\n")

geom <- st_geometry(lcp_sf)
is_null  <- sapply(geom, is.null)
is_empty <- st_is_empty(lcp_sf)
geom_type <- as.character(st_geometry_type(lcp_sf))

valid_type <- geom_type %in% c("LINESTRING", "MULTILINESTRING")
valid_idx <- which(!is_null & !is_empty & valid_type)

cat(glue("Total LCPs: {nrow(lcp_sf)}\n"))
cat(glue("Valid geometries: {length(valid_idx)}\n\n"))

lcp_sf <- lcp_sf[valid_idx, ]

# ================================================================
# COST FUNCTION (FULLY SAFE)
# ================================================================
calc_lcp_cost <- function(geom) {
  
  if (is.null(geom)) return(NA_real_)
  
  # safely extract geometry
  if (!inherits(geom, "sfg")) {
    geom <- st_geometry(geom)[[1]]
    if (is.null(geom)) return(NA_real_)
  }
  
  # merge MULTILINESTRING
  if (st_geometry_type(geom) == "MULTILINESTRING") {
    merged <- try(st_line_merge(geom), silent = TRUE)
    if (!inherits(merged, "try-error") && !is.null(merged)) {
      geom <- merged
    }
  }
  
  # densify safely
  seg_length <- res(acc_ras)[1] / 3
  
  dense <- try(st_segmentize(geom, dfMaxLength = seg_length), silent = TRUE)
  if (inherits(dense, "try-error") || is.null(dense)) return(NA_real_)
  
  xy <- try(st_coordinates(dense), silent = TRUE)
  if (inherits(xy, "try-error")) return(NA_real_)
  if (!all(c("X","Y") %in% colnames(xy))) return(NA_real_)
  
  # SAFE EXTRACT FOR SINGLE-BAND RASTER
  ext <- terra::extract(acc_ras, xy[,1:2])
  
  if (is.null(ext) || nrow(ext) == 0)
    return(NA_real_)
  
  # cost values are always in the LAST column
  vals <- as.numeric(ext[, ncol(ext)])
  
  vals <- vals[is.finite(vals)]
  
  if (length(vals) == 0)
    return(NA_real_)
  
  sum(vals)
  
}

# ================================================================
# COMPUTE SUMMARY
# ================================================================
cat("📌 Computing LCP length + cost...\n\n")

summary_list <- lapply(seq_len(nrow(lcp_sf)), function(i) {
  
  ln <- lcp_sf[i, ]
  
  geom <- st_geometry(ln)[[1]]
  
  if (is.null(geom))
    return(data.frame(
      from_node = ln$from_id,
      to_node   = ln$to_id,
      lcp_length_m = NA,
      lcp_cost = NA
    ))
  
  # length (meters)
  len_m <- as.numeric(st_length(geom))
  
  # cost
  cst <- calc_lcp_cost(geom)
  
  data.frame(
    from_node    = ln$from_id,
    to_node      = ln$to_id,
    lcp_length_m = len_m,
    lcp_cost     = cst
  )
})

summary_df <- bind_rows(summary_list)

cat(glue("\n✔ Done. Total rows: {nrow(summary_df)}"))
cat(glue("\n   NA costs: {sum(is.na(summary_df$lcp_cost))}\n\n"))

# ================================================================
# SAVE CSV
# ================================================================
out_csv <- file.path(
  proj_root, "outputs/tables",
  glue("lcp_summary_{species}.csv")
)

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(summary_df, out_csv, row.names = FALSE)

cat("🎉 LCP SUMMARY SAVED:\n")
cat(glue("    {out_csv}\n\n"))
cat("Columns: from_node, to_node, lcp_length_m, lcp_cost\n")
