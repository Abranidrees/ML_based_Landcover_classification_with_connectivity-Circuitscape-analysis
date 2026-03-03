#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(raster)
  library(glue)
  library(yaml)
  library(gdistance)
})

# ============================================================================
# LOAD CONFIG + MASTER GRID
# ============================================================================
proj_root <- getwd()
paths     <- yaml::read_yaml(file.path(proj_root, "config", "paths.yml"))

check_file <- function(rel){
  fp <- if(grepl("^(?:[A-Za-z]:|/)", rel)) rel else file.path(proj_root, rel)
  if(!file.exists(fp)) stop(glue("❌ File not found: {fp}"))
  fp
}

ref_grid <- rast(check_file("outputs/rasters/forest_core_2004_2006.tif"))
aoi      <- vect(check_file(paths$aoi))

# ============================================================================
# USER SELECTION
# ============================================================================
species <- "PFL"     # "PFL" or "GEN"
period  <- "Late"    # "Early" or "Late"

stopifnot(species %in% c("PFL","GEN"))
stopifnot(period  %in% c("Early","Late"))

# ============================================================================
# LOAD RESISTANCE + CORE  
# ============================================================================
core_map <- list(
  Early = "forest_core_2004_2006.tif",
  Late  = "forest_core_2022_2024.tif"
)


#   pfl_Early, pfl_Late, gen_Early, gen_Late
res_key <- switch(species,
                  PFL = paste0("pfl_", period),
                  GEN = paste0("gen_", period)
)

if (is.null(paths$resistance[[res_key]])) {
  stop(glue("❌ paths.yml is missing resistance key: {res_key} (expected one of pfl_Early/pfl_Late/gen_Early/gen_Late)"))
}

res_path <- check_file(paths$resistance[[res_key]])
res_raw  <- rast(res_path)
core_raw <- rast(check_file(file.path("outputs/rasters", core_map[[period]])))

cat(glue("📌 Resistance selected: {res_key}\n"))
cat("📌 Aligning rasters...\n")

res  <- resample(res_raw, ref_grid, "near")
core <- resample(core_raw, ref_grid, "near")

aoi  <- project(aoi, crs(res))
res  <- mask(crop(res, aoi), aoi)
core <- mask(crop(core, aoi), aoi)

# ============================================================================
#  FIX NA / ZERO / NEGATIVE IN RESISTANCE
# ============================================================================
cat("📌 Fixing NA values in resistance...\n")

vals <- values(res)
mx   <- max(vals, na.rm = TRUE)
vals[is.na(vals)] <- mx * 5
vals[vals <= 0]   <- mx * 10
values(res) <- vals

px_size <- res(ref_grid)[1]
cat(glue("  ✔ Pixel size: {px_size} m\n\n"))

# ============================================================================
# NODES FROM FOREST CORE PATCHES ≥ 5 ha
# ============================================================================
min_patch_ha <- 5
cat(glue("📌 Deriving nodes from forest-core patches ≥ {min_patch_ha} ha...\n"))

core_bin <- core
v_core <- values(core_bin, mat = FALSE)
v_core[is.na(v_core)] <- 0
v_core[v_core != 1]   <- NA
core_bin[] <- v_core

patch_r    <- terra::patches(core_bin, directions = 8)
patch_poly <- terra::as.polygons(patch_r, dissolve = TRUE)
names(patch_poly)[1] <- "patch_id"

patch_poly$area_ha <- terra::expanse(patch_poly, unit="ha")
patch_sel <- patch_poly[
  !is.na(patch_poly$patch_id) & patch_poly$area_ha >= min_patch_ha
]

n_patches <- nrow(patch_sel)
cat(glue("  ✔ Patches ≥ {min_patch_ha} ha: {n_patches}\n\n"))
if(n_patches < 2) stop("Not enough patches for connectivity.")

patch_sf <- st_as_sf(patch_sel)
nodes_sf <- st_centroid(patch_sf)
nodes_xy <- st_coordinates(nodes_sf)[,1:2, drop=FALSE]
colnames(nodes_xy) <- c("x","y")

cat(glue("  ✔ Extracted {nrow(nodes_xy)} patch-centroid nodes\n\n"))

# code later references `nodes` but only created `nodes_xy`
nodes <- nodes_xy

# ============================================================================
# ACC COST TILE-BASED
# ============================================================================
tile_px <- 1500
overlap <- 500

ncols <- ncol(res)
nrows <- nrow(res)

step_x <- tile_px - overlap
step_y <- tile_px - overlap

tile_cols <- seq(1, ncols, by = step_x)
tile_rows <- seq(1, nrows, by = step_y)

cat(glue("📌 Tiles: {length(tile_cols)} × {length(tile_rows)}\n"))
cat(glue("    Tile size = {tile_px}px, overlap = {overlap}px\n\n"))

out_tiles <- file.path(proj_root, "outputs/temp_tiles")
dir.create(out_tiles, recursive = TRUE, showWarnings = FALSE)

process_tile <- function(x0, y0, tid){
  
  x1 <- min(x0 + tile_px - 1, ncols)
  y1 <- min(y0 + tile_px - 1, nrows)
  
  cells <- cellFromRowColCombine(res, y0:y1, x0:x1)
  ext_tile <- ext(res, cells = cells)
  
  cat(glue("➡ Tile {tid}...\n"))
  
  r_tile <- crop(res, ext_tile)
  
  if(all(is.na(values(r_tile)))){
    cat("   Empty tile → NA\n")
    o <- r_tile; values(o) <- NA
    outfile <- file.path(out_tiles, glue("tile_{tid}.tif"))
    writeRaster(o, outfile, overwrite=TRUE)
    return(outfile)
  }
  
  r_tile_r <- raster::raster(r_tile)
  raster::crs(r_tile_r) <- crs(r_tile)
  
  # COST = resistance × pixel size
  res_local <- raster::values(r_tile_r)
  res_local[res_local <= 0 | is.na(res_local)] <- max(res_local, na.rm=TRUE) * 10
  
  cost <- res_local * px_size
  speed <- 1 / cost
  
  r_speed <- raster::setValues(r_tile_r, speed)
  
  # Transition matrix
  tr <- gdistance::transition(r_speed, transitionFunction = mean, directions=8)
  tr <- gdistance::geoCorrection(tr, type="c", multpl=TRUE)
  
  # Local nodes 
  buf <- overlap * px_size
  nodes_local <- nodes[
    nodes[,1] >= xmin(ext_tile)-buf &
      nodes[,1] <= xmax(ext_tile)+buf &
      nodes[,2] >= ymin(ext_tile)-buf &
      nodes[,2] <= ymax(ext_tile)+buf,
    , drop = FALSE]
  
  if(nrow(nodes_local) == 0){
    cat("   No nodes → NA\n")
    o <- r_tile; values(o) <- NA
    outfile <- file.path(out_tiles, glue("tile_{tid}.tif"))
    writeRaster(o, outfile, overwrite=TRUE)
    return(outfile)
  }
  
  acc_r <- gdistance::accCost(tr, nodes_local)
  raster::crs(acc_r) <- crs(r_tile)
  
  acc_t <- rast(acc_r)
  crs(acc_t) <- crs(r_tile)
  
  outfile <- file.path(out_tiles, glue("tile_{tid}.tif"))
  writeRaster(acc_t, outfile, overwrite=TRUE)
  
  return(outfile)
}

all_tiles <- list()
tid <- 1

for(cx in tile_cols){
  for(cy in tile_rows){
    all_tiles[[tid]] <- process_tile(cx, cy, tid)
    tid <- tid + 1
  }
}

tile_paths <- unlist(all_tiles)

# Check valid tiles
is_valid_tile <- function(path) {
  r <- rast(path)
  rv <- values(r)
  if (all(is.na(rv))) return(FALSE)
  if (!all(is.finite(rv), na.rm = TRUE)) return(FALSE)
  TRUE
}

valid_tiles <- tile_paths[sapply(tile_paths, is_valid_tile)]

cat("\n✔ VALID TILES FOUND:", length(valid_tiles), "of", length(tile_paths), "\n")
cat("🔄 Creating harmonized mosaic...\n")

rs <- lapply(valid_tiles, rast)
for (i in seq_along(rs)) crs(rs[[i]]) <- crs(ref_grid)

raw_mosaic <- do.call(mosaic, c(rs, fun = min))
names(raw_mosaic) <- "acc_cost"

cat("   ✔ Initial mosaic created\n")

# HARMONIZATION
hp <- raw_mosaic - focal(raw_mosaic, w = matrix(1,3,3), fun=mean, na.policy="omit")
hp_smooth <- focal(hp, w = matrix(1,5,5), fun=mean, na.policy="omit")
final_harmonized <- raw_mosaic - hp + hp_smooth
aoi_proj <- project(aoi, crs(final_harmonized))
final_clipped <- mask(crop(final_harmonized, aoi_proj), aoi_proj)

# SAVE OUTPUT
outfile <- file.path(
  proj_root,
  "outputs/rasters",
  glue("acc_cost_{species}_{period}.tif")
)
writeRaster(final_clipped, outfile, overwrite = TRUE)
cat(glue("\n🎯 FINAL SEAMLESS OUTPUT SAVED:\n    {outfile}\n\n"))

# ============================================================================
# CLEAN UP TEMPORARY TILE FOLDER
# ============================================================================
gd_tiles_dir <- file.path(proj_root, "outputs/temp_tiles")
if (dir.exists(gd_tiles_dir)) {
  unlink(gd_tiles_dir, recursive = TRUE, force = TRUE)
  cat("   ✔ Deleted: gd_tiles folder\n")
} else {
  cat("   ✔ No gd_tiles folder found, nothing to delete\n")
}
cat("🧩 Cleanup complete.\n")

# ============================================================================
#  LEAST-COST PATH VECTORS BETWEEN CORE PATCHES 
# ============================================================================
k_nn            <- 1
max_dist_km     <- 15
local_buffer_km <- 5

coords <- nodes_xy
n_nodes <- nrow(coords)
cat(glue("📌 Using {n_nodes} nodes for LCP.\n"))

pair_i <- integer(0)
pair_j <- integer(0)
max_dist_m <- max_dist_km * 1000

cat("📌 Building node pairs...\n")

for(i in seq_len(n_nodes)){
  dx <- coords[,1] - coords[i,1]
  dy <- coords[,2] - coords[i,2]
  d  <- sqrt(dx^2 + dy^2)
  d[i] <- Inf
  
  cand <- which(d <= max_dist_m)
  if(!length(cand)) next
  
  k_use <- min(k_nn, length(cand))
  neighs <- cand[order(d[cand])[1:k_use]]
  
  pair_i <- c(pair_i, rep(i, k_use))
  pair_j <- c(pair_j, neighs)
}

pair_df <- unique(data.frame(
  i = pmin(pair_i, pair_j),
  j = pmax(pair_i, pair_j)
))
pair_i <- pair_df$i
pair_j <- pair_df$j

n_pairs <- length(pair_i)
cat(glue("  ✔ {n_pairs} unique pairs.\n\n"))
if(n_pairs==0) stop("No node pairs generated.")

build_lcp_local <- function(i, j){
  p1 <- coords[i,]
  p2 <- coords[j,]
  
  dist_ij <- sqrt(sum((p1-p2)^2))
  buf_m <- max(local_buffer_km*1000, dist_ij*0.25)
  
  xmin <- min(p1[1],p2[1]) - buf_m
  xmax <- max(p1[1],p2[1]) + buf_m
  ymin <- min(p1[2],p2[2]) - buf_m
  ymax <- max(p1[2],p2[2]) + buf_m
  
  e <- ext(res)
  xmin <- max(xmin, e$xmin)
  xmax <- min(xmax, e$xmax)
  ymin <- max(ymin, e$ymin)
  ymax <- min(ymax, e$ymax)
  
  if(xmin>=xmax || ymin>=ymax) return(NULL)
  
  bbox <- ext(xmin, xmax, ymin, ymax)
  res_sub <- crop(res, bbox)
  
  if(is.null(res_sub)) return(NULL)
  if(all(is.na(values(res_sub, mat=FALSE)))) return(NULL)
  
  res_sub_r <- raster(res_sub)
  raster::crs(res_sub_r) <- crs(res_sub)
  
  cost_vec <- values(res_sub_r) * px_size
  cond_vec <- 1/cost_vec
  cond_vec[!is.finite(cond_vec)] <- 0
  
  cond_r <- setValues(res_sub_r, cond_vec)
  tr <- transition(cond_r, mean, directions=8)
  tr <- geoCorrection(tr, type="c", multpl=TRUE)
  
  # compute VALID pairwise cost distance and store it
  cd <- try(gdistance::costDistance(tr, matrix(p1,1), matrix(p2,1)), silent = TRUE)
  cd_val <- if (inherits(cd, "try-error") || length(cd) == 0) NA_real_ else as.numeric(cd[1,1])
  
  sp_line <- try(
    shortestPath(tr, matrix(p1,1), matrix(p2,1), output="SpatialLines"),
    silent=TRUE
  )
  if(inherits(sp_line,"try-error") || is.null(sp_line)) return(NULL)
  
  sf_line <- st_as_sf(sp_line)
  sf_line$from_id  <- nodes_sf$patch_id[i]
  sf_line$to_id    <- nodes_sf$patch_id[j]
  sf_line$costDist <- cd_val
  sf_line$length_m <- as.numeric(st_length(sf_line))
  
  sf_line
}

cat("📌 Tracing LCPs...\n")
lcp_list <- list()

for(k in seq_len(n_pairs)){
  out <- build_lcp_local(pair_i[k], pair_j[k])
  if(!is.null(out)) lcp_list[[length(lcp_list)+1]] <- out
  
  if(k %% 20 == 0) cat(glue("   ... {k}/{n_pairs}\n"))
}

lcp_list <- lcp_list[lengths(lcp_list)>0]
if(!length(lcp_list)) stop("No valid LCPs produced.")

lcp_sf <- do.call(rbind, lcp_list)
st_crs(lcp_sf) <- st_crs(nodes_sf)

outfile_lcp <- file.path(
  proj_root, "outputs/vectors",
  glue("LCP_{species}_{period}.gpkg")
)
dir.create(dirname(outfile_lcp), recursive=TRUE, showWarnings = FALSE)
st_write(lcp_sf, outfile_lcp, delete_dsn=TRUE, quiet=TRUE)

cat("\n🎯 LCP NETWORK SAVED:\n")
cat("    ", outfile_lcp, "\n\n")
