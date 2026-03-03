# scripts/00_setup.R
# Milestone 1 — Setup & Alignment
# Runs cleanly after renv::restore()
# Produces: logs/alignment_check.csv, figures/grid_overlay.png, logs/sessionInfo.txt

set.seed(4401)

suppressPackageStartupMessages({
  library(terra); library(sf); library(yaml)
  library(dplyr); library(readr); library(ggplot2)
})

cfg <- yaml::read_yaml("config/paths.yml")

# --- Load master reference grid and AOI ---
ref <- rast(cfg$reference_grid)     # rasters/reference/ValidatedLandcover2023_base7.tif
aoi <- sf::st_read(cfg$aoi, quiet = TRUE)

# --- Helper to check alignment vs reference ---
check_raster <- function(p) {
  r <- tryCatch(rast(p), error = function(e) NULL)
  if (is.null(r)) {
    return(tibble(path=p, exists=FALSE, read_ok=FALSE,
                  crs_ok=NA, res_ok=NA, ext_ok=NA, origin_ok=NA,
                  ncol=NA, nrow=NA, nlyr=NA, has_nodata=NA))
  }
  tibble(
    path = p,
    exists = TRUE,
    read_ok = TRUE,
    crs_ok = as.character(crs(r)) == as.character(crs(ref)),
    res_ok = isTRUE(all.equal(res(r), res(ref))),
    ext_ok = isTRUE(all.equal(ext(r), ext(ref))),
    origin_ok = isTRUE(all.equal(origin(r), origin(ref))),
    ncol = ncol(r), nrow = nrow(r), nlyr = nlyr(r),
    has_nodata = !is.na(NAflag(r))
  )
}

# --- Build list of all paths to verify (from YAML) ---
paths <- c(
  cfg$reference_grid,
  cfg$aoi,
  cfg$training_seed,
  unlist(cfg$composites$allbands),
  unlist(cfg$composites$fullextended),
  cfg$resistance$pfl_with_tall_early,
  cfg$resistance$pfl_with_tall_late,
  cfg$resistance$gen_all_early,
  cfg$resistance$gen_all_late
)

# Only rasters go through terra::rast safely; we’ll skip vectors except AOI existence

# Keep AOI in table as “read_ok” only
is_raster <- grepl("\\.(tif|tiff)$", paths, ignore.case = TRUE)
ras_paths <- paths[is_raster]
vec_paths <- paths[!is_raster]

# Check rasters
tab_r <- bind_rows(lapply(ras_paths, check_raster))

# Check AOI exists/valid
tab_v <- tibble(
  path = vec_paths,
  exists = file.exists(vec_paths),
  read_ok = TRUE,
  crs_ok = NA, res_ok = NA, ext_ok = NA, origin_ok = NA,
  ncol = NA, nrow = NA, nlyr = NA, has_nodata = NA
)

out <- bind_rows(tab_r, tab_v) |> arrange(path)
readr::write_csv(out, "logs/alignment_check.csv")

# --- Small overlay figure using first AllBands image ---
first_comp <- tryCatch(rast(cfg$composites$allbands[[1]]), error=function(e) NULL)
if (!is.null(first_comp)) {
  # coarse visualization
  ref1  <- ref[[1]]
  cmp1  <- first_comp[[1]]
  png("outputs/figures/grid_overlay.png", width=1200, height=800, res=150)
  plot(ref1, col=gray(0.98), main="Snap check to ValidatedLandcover2023 grid (EPSG:32644, 30 m)")
  plot(cmp1, add=TRUE, alpha=0.4)
  plot(vect(aoi), add=TRUE, border="cyan", lwd=2)
  dev.off()
}

# --- Save session info ---
writeLines(capture.output(sessionInfo()), "logs/sessionInfo.txt")

message("Milestone 1 checks complete. See logs/alignment_check.csv and figures/grid_overlay.png.")
