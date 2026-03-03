#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(yaml)
  library(glue)
})

# ============================================================
# SPEED / GENERAL SETTINGS
# ============================================================
agg_factor <- 15L      # 30 m * 15 = 450 m
max_nodes  <- 100L
period       <- "Late"              # RUN TWICE: "Early" and then "Late"
species_list <- c("PFL", "GEN")
min_patch_ha <- 5

# Avoid inflated core during aggregation

core_prop_thr <- 0.5

# Reproducibility  
delete_tmp <- FALSE

# Use 4 Julia threads (adjust to your CPU)
Sys.setenv(JULIA_NUM_THREADS = "4")

# ============================================================
# LOAD CONFIG + AOI + LOG SETUP
# ============================================================
proj_root <- getwd()
paths     <- yaml::read_yaml(file.path(proj_root, "config", "paths.yml"))

check_file <- function(rel){
  fp <- if (grepl("^(?:[A-Za-z]:|/)", rel)) rel else file.path(proj_root, rel)
  if (!file.exists(fp)) stop(glue("❌ File not found: {fp}"))
  fp
}

cat("\n=======================================\n")
cat(glue(" CIRCUITSCAPE SETUP — PERIOD = {period}\n"))
cat(glue(" agg_factor     = {agg_factor}\n"))
cat(glue(" max_nodes      = {max_nodes}\n"))
cat(glue(" min_patch_ha   = {min_patch_ha}\n"))
cat(glue(" core_prop_thr  = {core_prop_thr}\n"))
cat("=======================================\n\n")

# AOI
aoi_path <- check_file(paths$aoi)
cat(glue("📌 Using AOI from: {aoi_path}\n"))
aoi <- vect(aoi_path)

# Run-log setup
log_dir  <- file.path(proj_root, "outputs", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp    <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_path <- file.path(log_dir, glue("circuitscape_run_log_{stamp}.txt"))

run_summary <- data.frame(
  species      = character(),
  period       = character(),
  n_nodes      = integer(),
  n_pairs      = numeric(),
  start_time   = character(),
  end_time     = character(),
  duration_min = numeric(),
  ini_path     = character(),
  stringsAsFactors = FALSE
)

# ============================================================
# MAP inputs
# ============================================================
core_map <- list(
  Early = "forest_core_2004_2006.tif",
  Late  = "forest_core_2022_2024.tif"
)


# Period-specific resistance keys under paths$resistance
res_map <- list(
  PFL = list(Early = "pfl_Early", Late = "pfl_Late"),
  GEN = list(Early = "gen_Early", Late = "gen_Late")
)

# Validate resistance keys exist in paths.yml
for (sp in names(res_map)) {
  for (per in names(res_map[[sp]])) {
    key <- res_map[[sp]][[per]]
    if (is.null(paths$resistance[[key]])) {
      stop(glue("❌ paths.yml missing resistance entry for key: {key}"))
    }
  }
}

# ============================================================
# TEMP ROOT FOR ALL CIRCUITSCAPE FILES (INCLUDING GPKG)
# ============================================================
tmp_root <- file.path(proj_root, "outputs", "tmp_circuitscape")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

# Node outputs 
nodes_vec_out <- file.path(
  tmp_root,
  glue("forest_nodes_{period}_agg{agg_factor}_thr{core_prop_thr}_max{max_nodes}.gpkg")
)
nodes_ras_out <- file.path(
  tmp_root,
  glue("forest_nodes_{period}_agg{agg_factor}_thr{core_prop_thr}_max{max_nodes}.tif")
)

n_nodes <- NA_integer_
nodes_sf <- NULL
nodes_ras <- NULL

# ============================================================
# BUILD / LOAD NODE RASTER FOR THIS PERIOD
# ============================================================
if (file.exists(nodes_ras_out)) {
  cat(glue("📌 Existing node raster found (temp): {nodes_ras_out}\n"))
  nodes_ras <- rast(nodes_ras_out)
  
  vals    <- values(nodes_ras, mat = FALSE)
  n_nodes <- length(unique(na.omit(vals)))
  cat(glue("  ✔ Nodes (patches) in raster: {n_nodes}\n\n"))
  
  if (file.exists(nodes_vec_out)) {
    nodes_sf <- sf::st_read(nodes_vec_out, quiet = TRUE)
    cat(glue("  ✔ Nodes vector found (temp): {nodes_vec_out} ({nrow(nodes_sf)} points)\n\n"))
  } else {
    cat("  ℹ No vector nodes file found (only raster will be used by Circuitscape).\n\n")
  }
  
} else {
  
  # ----------------------------------------
  # Compute nodes from forest core and save (TEMP)
  # ----------------------------------------
  core_path <- check_file(file.path("outputs/rasters", core_map[[period]]))
  core_ras  <- rast(core_path)
  cat(glue("📌 Loaded forest-core raster for {period}: {core_path}\n"))
  
  
  # Aggregate using MEAN (proportion) + fixed threshold
  cat(glue("📌 Aggregating forest core by MEAN (proportion) factor = {agg_factor}\n"))
  core_ras_agg <- aggregate(core_ras, fact = agg_factor, fun = "mean", na.rm = TRUE)
  
  cat(glue("📌 Applying fixed proportion threshold ≥ {core_prop_thr}\n"))
  core_bin <- ifel(core_ras_agg >= core_prop_thr, 1, NA)
  
  cat(glue("📌 Deriving patches ≥ {min_patch_ha} ha on aggregated grid...\n"))
  
  # Patches with identical neighbor rule across periods
  patch_r    <- terra::patches(core_bin, directions = 8)
  patch_poly <- terra::as.polygons(patch_r, dissolve = TRUE)
  names(patch_poly)[1] <- "patch_id"
  
  patch_poly$area_ha <- terra::expanse(patch_poly, unit = "ha")
  patch_sel <- patch_poly[
    !is.na(patch_poly$patch_id) & patch_poly$area_ha >= min_patch_ha
  ]
  
  n_patches_raw <- nrow(patch_sel)
  cat(glue("  ✔ Patches ≥ {min_patch_ha} ha (before capping): {n_patches_raw}\n"))
  
  if (n_patches_raw < 2) {
    stop("Not enough patches ≥ 5 ha for connectivity analysis.")
  }
  
  # Cap nodes by area
  if (n_patches_raw > max_nodes) {
    cat(glue("  ⚠ Too many nodes ({n_patches_raw}) → keeping largest {max_nodes} by area\n"))
    ord <- order(patch_sel$area_ha, decreasing = TRUE)
    patch_sel <- patch_sel[ord[seq_len(max_nodes)], ]
  }
  
  n_nodes <- nrow(patch_sel)
  cat(glue("  ✔ Nodes used for Circuitscape: {n_nodes}\n"))
  
  
  # Safer point placement
  patch_sf <- sf::st_as_sf(patch_sel)
  nodes_sf <- sf::st_point_on_surface(patch_sf)
  
  cat(glue("  ✔ Extracted {nrow(nodes_sf)} node points\n"))
  sf::st_write(nodes_sf, nodes_vec_out, delete_dsn = TRUE, quiet = TRUE)
  cat(glue("  ✔ Saved nodes vector (temp) to: {nodes_vec_out}\n"))
  
  # Build node raster
  keep_ids <- patch_sel$patch_id
  nodes_ras <- patch_r
  v_nodes <- values(nodes_ras, mat = FALSE)
  v_nodes[!(v_nodes %in% keep_ids)] <- NA
  nodes_ras[] <- v_nodes
  
  writeRaster(nodes_ras, nodes_ras_out, overwrite = TRUE)
  cat(glue("  ✔ Saved node raster (temp) to: {nodes_ras_out}\n\n"))
}

# INFO: pairwise number of solves
if (!is.na(n_nodes)) {
  n_pairs <- n_nodes * (n_nodes - 1) / 2
  cat(glue("ℹ Pairwise mode: about {format(n_pairs, big.mark=',')} solves for {n_nodes} nodes.\n\n"))
}

# ============================================================
# CREATE CIRCUITSCAPE .ini FILES (USING TEMP OUTPUT FOLDER)
# ============================================================
ini_dir <- file.path(proj_root, "config")
dir.create(ini_dir, showWarnings = FALSE)

create_ini <- function(species, period, res_path_cs, nodes_ras_path) {
  
  ini_path <- file.path(
    ini_dir,
    glue("circuitscape_{species}_{period}_agg{agg_factor}_max{max_nodes}.ini")
  )
  
  out_base <- file.path(
    tmp_root,
    glue("cs_{species}_{period}_agg{agg_factor}_max{max_nodes}")
  )
  
  dir.create(dirname(out_base), recursive = TRUE, showWarnings = FALSE)
  
  lines <- c(
    "[Options]",
    "data_type = raster",
    "scenario = pairwise",
    "write_cur_maps = True",
    "write_cum_cur_map_only = True",
    "write_volt_maps = False",
    "",
    "[Pairwise]",
    glue("habitat_file = {res_path_cs}"),
    glue("point_file   = {nodes_ras_path}"),
    glue("output_file  = {out_base}"),
    "",
    "[Raster]",
    "habitat_map_is_resistances = True",
    "connect_four_neighbors      = False",
    "connect_eight_neighbors     = True",
    "connect_using_avg_resistances = True"
  )
  
  writeLines(lines, ini_path)
  return(ini_path)
}

# ============================================================
# PREPARE RESISTANCE RASTERS (AGGREGATED, TEMP)
# ============================================================
ini_paths <- list()

for (sp in species_list) {
  
  # ---- Period-specific resistance selection (ENFORCED) ----
  res_key <- res_map[[sp]][[period]]
  if (is.null(res_key)) {
    stop(glue("❌ No resistance mapping for species={sp}, period={period}"))
  }
  
  res_path_orig <- check_file(paths$resistance[[res_key]])
  
  cat(glue("\n📌 Resistance selected: species={sp}, period={period} → key={res_key}\n"))
  cat(glue("   Path: {res_path_orig}\n"))
  
  # Aggregated version for Circuitscape
  res_base    <- sub("\\.tif$", "", basename(res_path_orig))
  res_cs_path <- file.path(
    tmp_root,
    glue("{res_base}_cs_agg{agg_factor}.tif")
  )
  
  if (!file.exists(res_cs_path)) {
    cat(glue("📌 Creating aggregated resistance for {sp} ({period}) (factor {agg_factor})\n"))
    r30 <- rast(res_path_orig)
    
    if (agg_factor > 1L) {
      r_cs <- aggregate(r30, fact = agg_factor, fun = "mean", na.rm = TRUE)
    } else {
      r_cs <- r30
    }
    
    writeRaster(r_cs, res_cs_path, overwrite = TRUE)
    cat(glue("  ✔ Saved (temp): {res_cs_path}\n"))
  } else {
    cat(glue("📌 Using existing aggregated resistance (temp): {res_cs_path}\n"))
  }
  
  # --------------------------------------------------------
  # Use aggregated resistance as template for nodes
  # --------------------------------------------------------
  template <- rast(res_cs_path)
  
  if (is.null(nodes_ras)) {
    nodes_ras <- rast(nodes_ras_out)
  }
  
  nodes_ras_aligned <- crop(nodes_ras, template)
  nodes_ras_aligned <- resample(nodes_ras_aligned, template, method = "near")
  nodes_ras_aligned <- as.int(nodes_ras_aligned)
  
  if (!compareGeom(nodes_ras_aligned, template, stopOnError = FALSE)) {
    stop(glue("❌ Hnodes grid does not match resistance grid for {sp} {period}"))
  }
  
  writeRaster(nodes_ras_aligned, nodes_ras_out, overwrite = TRUE)
  nodes_ras <- nodes_ras_aligned
  
  cat(glue("  ✔ Nodes aligned to resistance template for {sp} {period}\n"))
  
  # Create INI
  ini_paths[[sp]] <- create_ini(
    species        = sp,
    period         = period,
    res_path_cs    = res_cs_path,
    nodes_ras_path = nodes_ras_out
  )
}

cat("\n=======================================\n")
cat(" CIRCUITSCAPE SETUP COMPLETE\n")
cat(" INI FILES CREATED:\n")
print(ini_paths)
cat("=======================================\n\n")

# Also create copies with required names (overwritten on each run)
if (!is.null(ini_paths$PFL)) {
  file.copy(ini_paths$PFL,
            file.path(ini_dir, "circuitscape_pfl.ini"),
            overwrite = TRUE)
}
if (!is.null(ini_paths$GEN)) {
  file.copy(ini_paths$GEN,
            file.path(ini_dir, "circuitscape_gen.ini"),
            overwrite = TRUE)
}

# ============================================================
# RUN CIRCUITSCAPE FROM R (JULIA FROM paths.yml)
# + CONVERT CUMULATIVE MAP TO FINAL CLIPPED GeoTIFF
# ============================================================
get_julia_exe <- function() {
  if (!is.null(paths$julia)) {
    return(check_file(paths$julia))
  }
  julia <- Sys.which("julia")
  if (nzchar(julia)) return(julia)
  
  stop(
    "Could not find Julia.\n",
    "Set 'julia' in config/paths.yml, e.g.\n",
    '  julia: "tools/Julia-1.12.1/bin/julia.exe"'
  )
}

run_circuitscape <- function(ini_path, species) {
  cat(glue("\n🚀 Running Circuitscape for: {ini_path}\n"))
  
  julia_exe <- get_julia_exe()
  ini_norm  <- chartr("\\", "/", ini_path)
  expr <- sprintf('using Circuitscape; compute("%s")', ini_norm)
  
  start_time <- Sys.time()
  
  system2(
    command = julia_exe,
    args    = c("-t", "4", "-e", shQuote(expr)),
    stdout  = "",
    stderr  = ""
  )
  
  end_time <- Sys.time()
  dur_min  <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  run_summary <<- rbind(
    run_summary,
    data.frame(
      species      = species,
      period       = period,
      n_nodes      = n_nodes,
      n_pairs      = if (!is.na(n_nodes)) n_nodes * (n_nodes - 1) / 2 else NA_real_,
      start_time   = format(start_time, "%Y-%m-%d %H:%M:%S"),
      end_time     = format(end_time, "%Y-%m-%d %H:%M:%S"),
      duration_min = dur_min,
      ini_path     = ini_path,
      stringsAsFactors = FALSE
    )
  )
  
  # After Circuitscape: convert cumulative ASC to final TIFF
  base_name <- file.path(
    tmp_root,
    glue("cs_{species}_{period}_agg{agg_factor}_max{max_nodes}")
  )
  cum_asc <- paste0(base_name, "_cum_curmap.asc")
  
  if (!file.exists(cum_asc)) {
    warning(glue("❗ Cumulative current map not found: {cum_asc}"))
    return(invisible(NULL))
  }
  
  species_lower <- tolower(species)
  final_tif <- file.path(
    proj_root, "outputs", "rasters",
    glue("cs_current_{species_lower}_{period}.tif")
  )
  
  cat(glue("📌 Converting {basename(cum_asc)} → {basename(final_tif)} (clipped to AOI)\n"))
  cs_r <- rast(cum_asc)
  
  cs_r <- crop(cs_r, aoi)
  cs_r <- mask(cs_r, aoi)
  
  writeRaster(cs_r, final_tif, overwrite = TRUE)
  cat(glue("  ✔ Final clipped Circuitscape raster written: {final_tif}\n"))
  
  invisible(NULL)
}

for (sp in names(ini_paths)) {
  run_circuitscape(ini_paths[[sp]], species = sp)
}

cat("\n✅ Circuitscape runs completed for period = ", period, "\n")
cat("   In outputs/rasters you now have:\n")
cat(glue("   - cs_current_pfl_{period}.tif (clipped to AOI)\n"))
cat(glue("   - cs_current_gen_{period}.tif (clipped to AOI)\n\n"))

# ============================================================
# WRITE SHORT RUN LOG
# ============================================================

get_n_nodes_for_period <- function(per) {
  nodes_vec <- file.path(
    tmp_root,
    glue("forest_nodes_{per}_agg{agg_factor}_thr{core_prop_thr}_max{max_nodes}.gpkg")
  )
  nodes_ras <- file.path(
    tmp_root,
    glue("forest_nodes_{per}_agg{agg_factor}_thr{core_prop_thr}_max{max_nodes}.tif")
  )
  
  if (file.exists(nodes_vec)) {
    n_nodes <- nrow(sf::st_read(nodes_vec, quiet = TRUE))
    return(list(n_nodes = n_nodes, source = nodes_vec))
  }
  
  if (file.exists(nodes_ras)) {
    r <- rast(nodes_ras)
    n_nodes <- length(unique(na.omit(values(r, mat = FALSE))))
    return(list(n_nodes = n_nodes, source = nodes_ras))
  }
  
  return(list(n_nodes = NA_integer_, source = NA_character_))
}

collect_existing_tmp_summary <- function(periods = c("Early", "Late"),
                                         species_list = c("PFL", "GEN")) {
  
  rows <- list()
  
  for (per in periods) {
    nn <- get_n_nodes_for_period(per)
    n_nodes <- nn$n_nodes
    n_pairs <- if (!is.na(n_nodes) && n_nodes >= 2) n_nodes * (n_nodes - 1) / 2 else NA_real_
    
    for (sp in species_list) {
      ini_path <- file.path(
        ini_dir,
        glue("circuitscape_{sp}_{per}_agg{agg_factor}_max{max_nodes}.ini")
      )
      
      base_name <- file.path(tmp_root, glue("cs_{sp}_{per}_agg{agg_factor}_max{max_nodes}"))
      cum_asc   <- paste0(base_name, "_cum_curmap.asc")
      
      rows[[length(rows) + 1]] <- data.frame(
        species         = sp,
        period          = per,
        n_nodes         = n_nodes,
        n_pairs         = n_pairs,
        nodes_source    = nn$source,
        ini_path        = if (file.exists(ini_path)) ini_path else NA_character_,
        ini_exists      = file.exists(ini_path),
        cum_asc         = cum_asc,
        cum_asc_exists  = file.exists(cum_asc),
        cum_asc_mtime   = if (file.exists(cum_asc)) format(file.info(cum_asc)$mtime, "%Y-%m-%d %H:%M:%S") else NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  
  do.call(rbind, rows)
}

run_summary_all <- collect_existing_tmp_summary(
  periods = c("Early", "Late"),
  species_list = species_list
)

# Write TXT log only
con <- file(log_path, "w")
writeLines("Circuitscape existing-output log (read from outputs/tmp_circuitscape)", con)
writeLines(glue("Created: {Sys.time()}"), con)
writeLines(glue("agg_factor: {agg_factor}"), con)
writeLines(glue("max_nodes: {max_nodes}"), con)
writeLines(glue("min_patch_ha: {min_patch_ha}"), con)
writeLines(glue("core_prop_thr: {core_prop_thr}"), con)
writeLines("", con)

writeLines("Summary table:", con)
writeLines(paste(capture.output(print(run_summary_all, row.names = FALSE)), collapse = "\n"), con)

# Quick warnings in the log (helpful for QC)
missing_ini <- run_summary_all[!run_summary_all$ini_exists, ]
missing_asc <- run_summary_all[!run_summary_all$cum_asc_exists, ]

writeLines("\nChecks:", con)
if (nrow(missing_ini) > 0) {
  writeLines("⚠ Missing INI files for:", con)
  writeLines(paste(capture.output(print(missing_ini[, c("species","period","ini_path")], row.names = FALSE)), collapse = "\n"), con)
} else {
  writeLines("✔ All expected INI files exist.", con)
}

if (nrow(missing_asc) > 0) {
  writeLines("\n⚠ Missing cumulative ASC files for:", con)
  writeLines(paste(capture.output(print(missing_asc[, c("species","period","cum_asc")], row.names = FALSE)), collapse = "\n"), con)
} else {
  writeLines("\n✔ All expected cumulative ASC files exist.", con)
}

close(con)

cat(glue("\n📝 Combined TXT log written: {log_path}\n"))


# ============================================================
# Normalize cumulative current by node pairs
# ============================================================

normalize_cum_current <- function(
    periods = c("Early", "Late"),
    species_list = c("PFL", "GEN")
) {
  suppressPackageStartupMessages({
    library(terra)
    library(sf)
    library(glue)
    library(yaml)
  })
  
  proj_root <- getwd()
  paths     <- yaml::read_yaml(file.path(proj_root, "config", "paths.yml"))
  
  check_file <- function(rel){
    fp <- if (grepl("^(?:[A-Za-z]:|/)", rel)) rel else file.path(proj_root, rel)
    if (!file.exists(fp)) stop(glue("❌ File not found: {fp}"))
    fp
  }
  
  # AOI
  aoi_path <- check_file(paths$aoi)
  aoi <- vect(aoi_path)
  
  # Same tmp_root as your main script
  tmp_root <- file.path(proj_root, "outputs", "tmp_circuitscape")
  if (!dir.exists(tmp_root)) stop(glue("❌ tmp_root not found: {tmp_root}"))
  
  out_dir <- file.path(proj_root, "outputs", "rasters")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (per in periods) {
    
    # Use the period node vector (fast + stable) to get n_nodes
    nodes_vec_out <- file.path(
      tmp_root,
      glue("forest_nodes_{per}_agg{agg_factor}_thr{core_prop_thr}_max{max_nodes}.gpkg")
    )
    
    if (file.exists(nodes_vec_out)) {
      nodes_sf <- sf::st_read(nodes_vec_out, quiet = TRUE)
      n_nodes  <- nrow(nodes_sf)
    } else {
      # Fallback: try nodes raster (slower)
      nodes_ras_out <- file.path(
        tmp_root,
        glue("forest_nodes_{per}_agg{agg_factor}_thr{core_prop_thr}_max{max_nodes}.tif")
      )
      if (!file.exists(nodes_ras_out)) {
        warning(glue("⚠ Nodes file not found for {per}. Skipping period. Expected: {nodes_vec_out}"))
        next
      }
      nr <- rast(nodes_ras_out)
      n_nodes <- length(unique(na.omit(values(nr, mat = FALSE))))
    }
    
    if (is.na(n_nodes) || n_nodes < 2) {
      warning(glue("⚠ Not enough nodes for {per} (n_nodes={n_nodes}). Skipping."))
      next
    }
    
    n_pairs <- n_nodes * (n_nodes - 1) / 2
    cat(glue("\n🧮 Normalizing period={per}: n_nodes={n_nodes}, n_pairs={format(n_pairs, big.mark=',')}\n"))
    
    for (sp in species_list) {
      base_name <- file.path(tmp_root, glue("cs_{sp}_{per}_agg{agg_factor}_max{max_nodes}"))
      cum_asc   <- paste0(base_name, "_cum_curmap.asc")
      
      if (!file.exists(cum_asc)) {
        warning(glue("❗ Missing cumulative ASC: {cum_asc}"))
        next
      }
      
      sp_lower <- tolower(sp)
      final_tif <- file.path(out_dir, glue("cs_current_{sp_lower}_{per}_mean.tif"))
      
      cat(glue("📌 {basename(cum_asc)} → {basename(final_tif)} (AOI clip + / n_pairs)\n"))
      
      r <- rast(cum_asc)
      
      # Clip to AOI 
      r <- crop(r, aoi)
      r <- mask(r, aoi)
      
      # Normalize
      r_mean <- r / n_pairs
      
      writeRaster(r_mean, final_tif, overwrite = TRUE)
      cat(glue("  ✔ Wrote: {final_tif}\n"))
    }
  }
  
  cat("\n✅ Done. Normalized mean-current rasters are in outputs/rasters/*_mean.tif\n")
}

# ---- run it (after you have ASC outputs) ----
normalize_cum_current(
  periods = c("Early", "Late"),
  species_list = c("PFL", "GEN")
)

# ============================================================
# Comparison of normalized Circuitscape current (mean per node-pair)
#
# Design: 2×2 factorial
#   Factor 1: Period  (Early vs Late)
#   Factor 2: Species (PFL vs GEN)
# Response: log10(mean_current + 1e-6)
#
# Model: log_current ~ period * species + block_id
#   - period         : Early vs Late
#   - species        : PFL vs GEN
#   - period:species : interaction (difference-in-differences)
#
# Sampling: forest-core only (aggregated to the Circuitscape grid)
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(glue)
})

# ---- Inputs----
mean_paths <- list(
  pfl_Early = file.path(proj_root, "outputs", "rasters", "cs_current_pfl_Early_mean.tif"),
  pfl_Late  = file.path(proj_root, "outputs", "rasters", "cs_current_pfl_Late_mean.tif"),
  gen_Early = file.path(proj_root, "outputs", "rasters", "cs_current_gen_Early_mean.tif"),
  gen_Late  = file.path(proj_root, "outputs", "rasters", "cs_current_gen_Late_mean.tif")
)

for (nm in names(mean_paths)) {
  if (!file.exists(mean_paths[[nm]])) stop(glue("❌ Missing raster: {mean_paths[[nm]]}"))
}

r_pfl_E <- rast(mean_paths$pfl_Early); names(r_pfl_E) <- "PFL_Early"
r_pfl_L <- rast(mean_paths$pfl_Late);  names(r_pfl_L) <- "PFL_Late"
r_gen_E <- rast(mean_paths$gen_Early); names(r_gen_E) <- "GEN_Early"
r_gen_L <- rast(mean_paths$gen_Late);  names(r_gen_L) <- "GEN_Late"

# Use one layer as the template grid for alignment
template <- r_pfl_E

# Ensure all layers share identical geometry (extent/resolution/origin/CRS)
align_to_template <- function(r, template) {
  if (!compareGeom(r, template, stopOnError = FALSE)) {
    r <- resample(r, template, method = "bilinear")  # continuous surface
  }
  r
}

r_pfl_L <- align_to_template(r_pfl_L, template)
r_gen_E <- align_to_template(r_gen_E, template)
r_gen_L <- align_to_template(r_gen_L, template)

r4 <- c(r_pfl_E, r_pfl_L, r_gen_E, r_gen_L)

# ---- Build forest-core-only sampling mask on the Circuitscape grid ----
core_E_path <- file.path(proj_root, "outputs", "rasters", "forest_core_2004_2006.tif")
core_L_path <- file.path(proj_root, "outputs", "rasters", "forest_core_2022_2024.tif")
if (!file.exists(core_E_path) || !file.exists(core_L_path)) {
  stop("❌ Forest-core rasters not found; cannot build forest-only sampling mask.")
}

core_E_30 <- rast(core_E_path)
core_L_30 <- rast(core_L_path)

# Aggregate forest core to the Circuitscape grid using mean proportion, then threshold.

core_E_prop <- aggregate(core_E_30, fact = agg_factor, fun = "mean", na.rm = TRUE)
core_L_prop <- aggregate(core_L_30, fact = agg_factor, fun = "mean", na.rm = TRUE)

core_E_prop <- align_to_template(core_E_prop, template)
core_L_prop <- align_to_template(core_L_prop, template)

# Forest-core definition:

forest_core_mask <- (core_E_prop >= core_prop_thr) | (core_L_prop >= core_prop_thr)

# Require valid values in all four current layers so all contrasts are comparable
mask_all4 <- !is.na(r_pfl_E) & !is.na(r_pfl_L) & !is.na(r_gen_E) & !is.na(r_gen_L)

# Restrict to AOI
forest_core_mask <- crop(forest_core_mask, aoi) |> mask(aoi)
mask_all4        <- crop(mask_all4, aoi) |> mask(aoi)

sample_mask <- forest_core_mask & mask_all4
cat(glue("🌲 Sampling only forest-core cells (thr={core_prop_thr}) with valid values in all four rasters.\n"))

# ---- Spatial blocks ----
block_size_m <- 5000  
block_id_r <- rast(template)
res(block_id_r) <- block_size_m
ext(block_id_r) <- ext(template)
crs(block_id_r) <- crs(template)
block_id_r <- init(block_id_r, "cell")  # unique block IDs

# ---- Sample points ----
set.seed(4401)
n_samples <- 20000L
sample_method <- "regular"  # "regular" (stable) or "random"

pts <- spatSample(sample_mask, size = n_samples, method = sample_method, as.points = TRUE, na.rm = TRUE)
cat(glue("📌 Sampled {nrow(pts)} points using method='{sample_method}'.\n"))

# ---- Extract values for the four rasters + block IDs ----
vals <- terra::extract(r4, pts, xy = TRUE, ID = TRUE)
blk  <- terra::extract(block_id_r, pts, ID = TRUE)
vals$block_id <- as.factor(blk[, 2])

# ---- Wide → long: one row per (point × scenario) ----
long <- reshape(
  vals,
  varying   = c("PFL_Early", "PFL_Late", "GEN_Early", "GEN_Late"),
  v.names   = "current",
  timevar   = "layer",
  times     = c("PFL_Early", "PFL_Late", "GEN_Early", "GEN_Late"),
  direction = "long"
)

# Split “PFL_Early” into species and period factors
sp_per <- do.call(rbind, strsplit(as.character(long$layer), "_"))
long$species <- factor(sp_per[, 1], levels = c("PFL", "GEN"))
long$period  <- factor(sp_per[, 2], levels = c("Early", "Late"))

# ---- Transform response (stabilize skew and handle zeros) ----
eps <- 1e-6
long$log_current <- log10(long$current + eps)

# ---- Fit factorial model: tests main effects + interaction ----
m <- lm(log_current ~ period * species + block_id, data = long)

# Extract the three key tests from the ANOVA table
anova_tbl <- anova(m)
tests_tbl <- anova_tbl[c("period", "species", "period:species"), c("Df", "F value", "Pr(>F)")]
colnames(tests_tbl) <- c("Df", "F_value", "p_value")

# ---- Interpretable contrasts (log10 scale) ----
co <- coef(m)

# Reference levels: period=Early, species=PFL
delta_PFL_EarlyToLate <- unname(co["periodLate"])                         # PFL: Late - Early
did_GENvsPFL          <- unname(co["periodLate:speciesGEN"])              # interaction (diff-in-diff)
delta_GEN_EarlyToLate <- unname(co["periodLate"] + co["periodLate:speciesGEN"])  # GEN: Late - Early

diff_Early_GENminusPFL <- unname(co["speciesGEN"])                        # Early: GEN - PFL
diff_Late_GENminusPFL  <- unname(co["speciesGEN"] + co["periodLate:speciesGEN"]) # Late: GEN - PFL

contrasts_tbl <- data.frame(
  contrast = c(
    "PFL: Late - Early",
    "GEN: Late - Early",
    "Early: GEN - PFL",
    "Late:  GEN - PFL",
    "Interaction (difference-in-differences)"
  ),
  estimate_log10 = c(
    delta_PFL_EarlyToLate,
    delta_GEN_EarlyToLate,
    diff_Early_GENminusPFL,
    diff_Late_GENminusPFL,
    did_GENvsPFL
  ),
  stringsAsFactors = FALSE
)

# ---- Write results to TXT ----
out_dir_tbl <- file.path(proj_root, "outputs", "tables")
dir.create(out_dir_tbl, recursive = TRUE, showWarnings = FALSE)

out_txt <- file.path(out_dir_tbl, "cs_current_comparison_forestcore_period_species.txt")

con <- file(out_txt, "w")
writeLines("Normalized Circuitscape current comparison (forest-core only)", con)
writeLines(glue("Created: {Sys.time()}"), con)
writeLines(glue("Sampling: method={sample_method}, n_samples={n_samples}, forest_core_thr={core_prop_thr}"), con)
writeLines("Mask: forest core (Early OR Late) AND valid in all four rasters", con)
writeLines(glue("Transform: log10(current + {eps})"), con)
writeLines("", con)

writeLines("Tests from ANOVA (model: log_current ~ period * species + block_id):", con)
writeLines(paste(capture.output(print(tests_tbl)), collapse = "\n"), con)

writeLines("\nContrasts (log10 scale):", con)
writeLines(paste(capture.output(print(contrasts_tbl, row.names = FALSE)), collapse = "\n"), con)

writeLines("\nModel summary:", con)
writeLines(paste(capture.output(print(summary(m))), collapse = "\n"), con)

close(con)

cat(glue("\n✅ Comparison finished. Results written to: {out_txt}\n"))


# =======================================================================
# DIFF RASTERS (Late − Early) from NORMALIZED mean-current maps
# =======================================================================

suppressPackageStartupMessages({
  library(terra)
  library(glue)
})

proj_root <- getwd()
ras_dir   <- file.path(proj_root, "outputs", "rasters")

species_list <- c("PFL", "GEN")

load_r <- function(path) {
  if (!file.exists(path)) stop(glue("❌ Missing raster: {path}"))
  rast(path)
}

align_to <- function(r, template) {
  if (!compareGeom(r, template, stopOnError = FALSE)) {
    r <- resample(r, template, method = "bilinear")  # continuous surface
  }
  r
}

make_diff <- function(sp) {
  sp_l <- tolower(sp)
  
  # Inputs: normalized mean rasters (per-node-pair)
  early_path <- file.path(ras_dir, glue("cs_current_{sp_l}_Early_mean.tif"))
  late_path  <- file.path(ras_dir, glue("cs_current_{sp_l}_Late_mean.tif"))
  
  cat(glue("\n📌 {sp}: Late − Early (using normalized mean-current rasters)\n"))
  cat(glue("   Early: {early_path}\n"))
  cat(glue("   Late : {late_path}\n"))
  
  r_early <- load_r(early_path)
  r_late  <- load_r(late_path)
  
  # Use Early grid as the template for subtraction
  r_late <- align_to(r_late, r_early)
  
  # Difference
  r_diff <- r_late - r_early
  names(r_diff) <- glue("{sp}_LateMinusEarly")
  
  # Output filename exactly as requested
  out_path <- file.path(ras_dir, glue("cs_diff_{sp_l}_LateMinusEarly.tif"))
  writeRaster(r_diff, out_path, overwrite = TRUE)
  
  cat(glue("  ✔ Wrote: {out_path}\n"))
  invisible(out_path)
}

out_files <- lapply(species_list, make_diff)

cat("\n✅ Difference rasters created:\n")
cat(paste0(" - ", unlist(out_files), collapse = "\n"), "\n")


# ==============================================
# CONNECTIVITY CHANGE TABLES (forest-core only)
# - Area (ha) of top-5% current by period
# - outputs/tables/connectivity_change_[pfl|gen].csv
# ==============================================

suppressPackageStartupMessages({
  library(terra)
  library(glue)
})

proj_root <- getwd()
ras_dir   <- file.path(proj_root, "outputs", "rasters")
tab_dir   <- file.path(proj_root, "outputs", "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

species_list <- c("PFL", "GEN")
top_p <- 0.95  # 95th percentile => top 5%


stopifnot(exists("agg_factor"), exists("core_prop_thr"), exists("aoi"))

load_r <- function(path) {
  if (!file.exists(path)) stop(glue("❌ Missing raster: {path}"))
  rast(path)
}

align_to <- function(r, template) {
  if (!compareGeom(r, template, stopOnError = FALSE)) {
    r <- resample(r, template, method = "bilinear")  # continuous surface
  }
  r
}

build_forest_core_mask <- function(template) {
  core_E_path <- file.path(ras_dir, "forest_core_2004_2006.tif")
  core_L_path <- file.path(ras_dir, "forest_core_2022_2024.tif")
  if (!file.exists(core_E_path) || !file.exists(core_L_path)) {
    stop("❌ forest_core rasters not found; cannot compute forest-core-only mask.")
  }
  
  core_E_30 <- rast(core_E_path)
  core_L_30 <- rast(core_L_path)
  
  core_E_prop <- aggregate(core_E_30, fact = agg_factor, fun = "mean", na.rm = TRUE)
  core_L_prop <- aggregate(core_L_30, fact = agg_factor, fun = "mean", na.rm = TRUE)
  
  core_E_prop <- align_to(core_E_prop, template)
  core_L_prop <- align_to(core_L_prop, template)
  
  forest_mask <- (core_E_prop >= core_prop_thr) | (core_L_prop >= core_prop_thr)
  
  forest_mask <- crop(forest_mask, aoi)
  forest_mask <- mask(forest_mask, aoi)
  
  forest_mask
}

top_area_ha <- function(r, mask_r, prob = 0.95) {
  # apply mask
  r_m <- crop(r, mask_r)
  r_m <- mask(r_m, mask_r)
  
  # pull values (non-NA)
  v <- values(r_m, mat = FALSE)
  v <- v[!is.na(v)]
  if (length(v) == 0) {
    return(list(threshold = NA_real_, n_cells = 0, area_ha = 0))
  }
  
  # threshold for "top (1-prob)" based on all valid cells
  thr <- as.numeric(stats::quantile(v, probs = prob, na.rm = TRUE, names = FALSE))
  
  # count cells >= threshold
  top_mask <- r_m >= thr
  n_cells  <- as.numeric(terra::global(top_mask, fun = "sum", na.rm = TRUE)[1, 1])
  
  # pixel area in ha (assumes meters)
  pix_area_ha <- abs(res(r_m)[1] * res(r_m)[2]) / 10000
  
  list(threshold = thr, n_cells = n_cells, area_ha = n_cells * pix_area_ha)
}


compute_connectivity_change <- function(sp) {
  sp_l <- tolower(sp)
  
  early_path <- file.path(ras_dir, glue("cs_current_{sp_l}_Early_mean.tif"))
  late_path  <- file.path(ras_dir, glue("cs_current_{sp_l}_Late_mean.tif"))
  
  cat(glue("\n📌 Connectivity change table for {sp} (top 5% area, forest-core only)\n"))
  cat(glue("   Early: {early_path}\n"))
  cat(glue("   Late : {late_path}\n"))
  
  r_early <- load_r(early_path)
  r_late  <- load_r(late_path)
  
  # Keep a single grid for consistent masking + quantiles
  template <- r_early
  r_late   <- align_to(r_late, template)
  
  forest_mask <- build_forest_core_mask(template)
  
  cat("  🔎 Computing Early top-5% area...\n")
  tE <- top_area_ha(r_early, forest_mask, prob = top_p)
  
  cat("  🔎 Computing Late top-5% area...\n")
  tL <- top_area_ha(r_late, forest_mask, prob = top_p)
  
  df <- data.frame(
    species           = sp,
    period            = c("Early", "Late"),
    threshold_current = c(tE$threshold, tL$threshold),
    top5_cells        = c(tE$n_cells,   tL$n_cells),
    top5_area_ha      = c(tE$area_ha,   tL$area_ha),
    delta_area_ha     = c(NA_real_,     tL$area_ha - tE$area_ha),
    delta_area_pct    = c(
      NA_real_,
      if (tE$area_ha > 0) 100 * (tL$area_ha - tE$area_ha) / tE$area_ha else NA_real_
    ),
    notes = "",
    stringsAsFactors = FALSE
  )
  
  out_csv <- file.path(tab_dir, glue("connectivity_change_{sp_l}.csv"))
  write.csv(df, out_csv, row.names = FALSE)
  
  cat(glue("  ✔ Wrote: {out_csv}\n"))
  invisible(out_csv)
}

out_files <- lapply(species_list, compute_connectivity_change)

cat("\n✅ Connectivity tables created:\n")
cat(paste0(" - ", unlist(out_files), collapse = "\n"), "\n")
