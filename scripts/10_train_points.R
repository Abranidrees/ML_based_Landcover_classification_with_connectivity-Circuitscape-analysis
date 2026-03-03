# ===============================
# Spatial K-fold split by hex blocks (≈70/30 per chosen fold)
# + De-cluster training points: cap per HEX×CLASS (+ optional min spacing)
# ===============================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(readr)
  library(stringr); library(ggplot2)
})

# ---------- Window & Inputs
WINDOW      <- "2022-2024" # 2004-2006; 2010-2012; 2016-2018; 2022-2024
aoi_path    <- "vectors/aoi.gpkg"
points_in   <- "vectors/Training_points_2022-24.gpkg"

# ---------- Outputs
seed_out    <- file.path("outputs", "vectors", paste0("train_seed_", WINDOW, ".gpkg"))
train_out   <- file.path("outputs", "vectors", paste0("train_", WINDOW, ".gpkg"))
counts_out  <- file.path("outputs", "tables", paste0("train_counts_", WINDOW, ".csv"))
notes_out   <- file.path("outputs", "tables", paste0("train_notes_",  WINDOW, ".csv"))
folds_png   <- file.path("outputs", "figures", paste0("training_folds_", WINDOW, ".png"))

dir.create(dirname(seed_out),  TRUE, TRUE)
dir.create(dirname(train_out), TRUE, TRUE)
dir.create(dirname(counts_out), TRUE, TRUE)
dir.create(dirname(folds_png), TRUE, TRUE)

# ---------- Parameters
BLOCK_KM              <- 3.0
KFOLDS                <- 3
CHOSEN_FOLD           <- 3
VALID_PROP            <- 0.30
GAP_M                 <- 0
N_PER_HEX_PER_CLASS   <- 6
MIN_DIST_M            <- 20
SEED                  <- 4401
set.seed(SEED)

# ---------- Helpers
auto_project_utm <- function(g){
  ctr  <- st_centroid(st_union(st_geometry(g))) |> st_coordinates()
  lon  <- ctr[1]; lat <- ctr[2]
  zone <- floor((lon + 180)/6) + 1
  epsg <- if (lat >= 0) 32600 + zone else 32700 + zone
  st_transform(g, epsg)
}

thin_min_dist <- function(sf_pts, min_dist_m){
  if (min_dist_m <= 0 || nrow(sf_pts) <= 1) return(sf_pts)
  keep_idx <- integer(0)
  cand_idx <- sample(seq_len(nrow(sf_pts)))
  for (i in cand_idx){
    if (!length(keep_idx)) { keep_idx <- i; next }
    d <- as.numeric(sf::st_distance(sf_pts[i,], sf_pts[keep_idx,]))
    if (all(d >= min_dist_m)) keep_idx <- c(keep_idx, i)
  }
  sf_pts[keep_idx,]
}

# ---------- Read AOI & points
aoi <- st_read(aoi_path, quiet = TRUE) |> st_make_valid()
pts <- st_read(points_in, quiet = TRUE) |> st_make_valid()
if (!"base7_code" %in% names(pts)) stop("Points need 'base7_code'")
if (!"notes" %in% names(pts)) stop("Points need 'notes'")

if (is.na(st_is_longlat(aoi)) || st_is_longlat(aoi)) aoi <- auto_project_utm(aoi)
pts <- st_transform(pts, st_crs(aoi))

# ---------- Seeds
pts_seed <- pts |>
  mutate(notes_chr = as.character(notes)) |>
  filter(grepl("\\bseed\\b", notes_chr, ignore.case = TRUE)) |>
  select(-notes_chr)
if (!"window" %in% names(pts_seed)) pts_seed$window <- WINDOW
st_write(pts_seed, seed_out, delete_layer = TRUE, quiet = TRUE)
cat("✅ Seed preview → ", seed_out, " (n = ", nrow(pts_seed), ")\n", sep = "")

# ---------- Hex grid + diversify + folds
block_m <- BLOCK_KM * 1000
hex_raw <- st_make_grid(aoi, cellsize = block_m, square = FALSE)
blocks  <- st_sf(block_id = seq_along(hex_raw), geometry = hex_raw) |>
  st_intersection(aoi) |>
  mutate(area_km2 = as.numeric(st_area(geometry))/1e6) |>
  filter(area_km2 > 0)
rm(hex_raw)

pts_b <- st_join(pts, blocks["block_id"], left = TRUE, largest = TRUE) |>
  filter(!is.na(block_id))

cap_group <- function(.x, key){
  x <- thin_min_dist(.x, MIN_DIST_M)
  n_keep <- min(nrow(x), N_PER_HEX_PER_CLASS)
  if (n_keep == 0L) return(x[0,, drop = FALSE])
  dplyr::slice_sample(x, n = n_keep, replace = FALSE)
}

set.seed(SEED)
pts_cap <- pts_b %>%
  dplyr::group_by(block_id, base7_code) %>%
  dplyr::group_modify(cap_group) %>%
  dplyr::ungroup()

blk_cls  <- pts_cap |> st_drop_geometry() |> count(block_id, base7_code, name = "n")
blk_load <- blk_cls |> group_by(block_id) |> summarise(total = sum(n), .groups = "drop") |> arrange(desc(total))

fold_vec   <- integer(nrow(blk_load)); fold_sizes <- rep(0L, KFOLDS)
for (i in seq_len(nrow(blk_load))) { j <- which.min(fold_sizes); fold_vec[i] <- j; fold_sizes[j] <- fold_sizes[j] + blk_load$total[i] }
blk_fold <- blk_load |> mutate(fold = fold_vec) |> select(block_id, fold)
blocks   <- blocks   |> left_join(blk_fold, by = "block_id")
pts_folded <- pts_cap |> left_join(blk_fold, by = "block_id")

# ---------- Split + GAP
pts_split <- pts_folded |> mutate(split = if_else(fold == CHOSEN_FOLD, "valid", "train"))
if (GAP_M > 0){
  valid_poly <- blocks |> filter(fold == CHOSEN_FOLD) |> st_union()
  buf <- st_buffer(valid_poly, GAP_M)
  inside_buf <- lengths(st_intersects(pts_split, buf)) > 0
  removed <- sum(inside_buf & pts_split$split == "train")
  pts_split <- pts_split |> filter(!(split == "train" & inside_buf))
  message(sprintf("Applied %d m gap; removed %d train points near valid edges.", GAP_M, removed))
}

# ---------- Final schema (KEEP geometry), write tables & file
keep_fields <- base::intersect(c("id","base7_code","window","notes"), names(pts_split))

# Make sure `tmp` is sf, then do a geometry-safe select:
tmp <- sf::st_as_sf(pts_split) %>%
  dplyr::mutate(window = if (!"window" %in% names(.)) WINDOW else window)

final <- {
  atts <- dplyr::select(tmp, dplyr::any_of(keep_fields), fold, split)
  # explicitly carry geometry column through select:
  sf::st_as_sf(dplyr::mutate(atts, geometry = sf::st_geometry(tmp)), sf_column_name = "geometry")
}

st_write(final, train_out, delete_layer = TRUE, quiet = TRUE)
cat("✅ Train/valid (by blocks, fold=", CHOSEN_FOLD, ") → ",
    train_out, " (n = ", nrow(final), ")\n", sep = "")

# Counts & notes
counts_wide <- final |> st_drop_geometry() |>
  count(base7_code, split, name = "N") |>
  complete(base7_code, split = c("train","valid"), fill = list(N = 0)) |>
  group_by(base7_code) |>
  summarise(total = sum(N), train = N[split=="train"], valid = N[split=="valid"], .groups="drop") |>
  arrange(base7_code)
write_csv(counts_wide, counts_out)

notes_counts <- final |> st_drop_geometry() |>
  mutate(notes = as.character(notes)) |>
  count(base7_code, notes, name = "N") |>
  arrange(base7_code, desc(N), notes)
write_csv(notes_counts, notes_out)

# ----------  Map
p <- ggplot() +
  geom_sf(data = blocks, aes(fill = factor(fold)), color = NA, alpha = 0.8) +
  geom_sf(data = sf::st_geometry(aoi), fill = NA, color = "black", linewidth = 0.3) +
  geom_sf(data = final |> dplyr::mutate(is_valid = split == "valid"),
          aes(shape = is_valid), size = 0.6, alpha = 0.7, color = "black") +
  scale_fill_brewer(palette = "Set2", name = "Fold") +
  scale_shape_manual(values = c(16,4), name = "Chosen split",
                     labels = c("Train", paste0("Valid (fold ", CHOSEN_FOLD, ")"))) +
  guides(shape = guide_legend(override.aes = list(color = "black", size = 2))) +
  labs(
    title    = paste0("Spatial K-folds by Hex Blocks — Window ", WINDOW),
    subtitle = paste0(
      "Hex ~", BLOCK_KM, " km · K=", KFOLDS,
      " · Cap=", N_PER_HEX_PER_CLASS, "/hex/class",
      if (MIN_DIST_M > 0) paste0(" · MinDist=", MIN_DIST_M, " m") else "",
      " · Gap=", GAP_M, " m · Seed=", SEED
    ),
    caption = "Points diversified by (hex × class) cap; blocks assigned entirely to a fold."
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

ggsave(folds_png, p, width = 8.5, height = 7, dpi = 300)
cat("🗺️  Folds map → ", folds_png, "\n", sep = "")

# ---------- QA: proportions
achieved <- final |> st_drop_geometry() |>
  count(base7_code, split, name = "N") |>
  group_by(base7_code) |>
  mutate(total = sum(N),
         prop_valid = ifelse(total > 0, N[split == "valid"]/total, NA_real_),
         prop_train = ifelse(total > 0, N[split == "train"]/total, NA_real_)) |>
  ungroup() |>
  arrange(base7_code, split)
message("---- Achieved per-class proportions (target valid ≈ ", round(100*VALID_PROP), "%) ----")
print(achieved |> select(base7_code, split, N, total, prop_valid) |> arrange(base7_code, split))
