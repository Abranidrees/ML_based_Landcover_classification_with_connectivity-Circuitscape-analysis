# =========================
# Reclassify 2023 raster → Base-7 and extract to training points (clipped + Int16)
# Requires: terra, sf, dplyr, readr (optional)
# =========================

library(terra)
library(sf)
library(dplyr)
library(readr)

# ---------- Inputs (edit paths as needed)
aoi_path      <- file.path("vectors", "aoi.gpkg")                     # your AOI polygon
in_raster     <- file.path("rasters", "reference", "ValidatedLandcover2023.tif")
seed_points   <- file.path("vectors", "TrainingPoints.gpkg")          # your training seeds (points)

# ---------- Outputs
out_raster    <- file.path("rasters", "reference", "ValidatedLandcover2023_base7.tif")
out_points    <- file.path("vectors", "train_seed_2023_base7.gpkg")
out_counts    <- file.path("tables", "train_counts_seed2023_base7.csv")

# ---------- Read data
r   <- rast(in_raster)                               # keeps original CRS/resolution/grid
aoi <- st_read(aoi_path, quiet = TRUE) |> st_make_valid()
aoi <- st_transform(aoi, crs(r))                     # align AOI to raster CRS

# ---------- Define reclassification (original_code → Base-7)
# Mapping:
# 0 TallF → 1, 7 OpenF → 1, 9 ShortF → 1,
# 4 HomeG → 2, 5 Plant → 2,
# 11 Paddy → 3,
# 1 Scrub → 4, 3 Grass → 4,
# 2 Built → 5,
# 10 Water → 6, 8 Wetlands → 6,
# 6 Rock → 7,
# 255 NoData → mask (NA)

lut <- data.frame(
  from = c(0, 7, 9, 4, 5, 11, 1, 3, 2, 10, 8, 6, 255),
  to   = c(1, 1, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7, NA)
)

# ---------- Reclassify
# others=NA ensures any unexpected codes become NA.
r_base7 <- classify(r, rcl = lut, others = NA)

# ---------- Clip to AOI (crop extent + mask edges), keep pixel size & alignment
r_base7 <- crop(r_base7, vect(aoi)) |> mask(vect(aoi))

# ---------- Write reclassified raster as Int16 (INT2S)
# Values are small (1..7) but you requested Int16.
dir.create(dirname(out_raster), showWarnings = FALSE, recursive = TRUE)
writeRaster(
  r_base7,
  filename = out_raster,
  overwrite = TRUE,
  gdal = c("COMPRESS=LZW", "TILED=YES"),
  datatype = "INT2S"   # Int16 (signed)
)

cat("✅ Wrote clipped + reclassified raster (Int16):", out_raster, "\n")

# ---------- Extract Base-7 to training seed points
pts <- st_read(seed_points, quiet = TRUE) |> st_make_valid()
pts <- st_transform(pts, crs(r))  # ensure same CRS as raster (uses source raster CRS)

# Extract returns a data.frame with an ID column and the raster value column
ex <- terra::extract(r_base7, vect(pts))
# Column 2 is the extracted value (column name equals layer name)
base7_vals <- ex[[2]]

# Append as new attribute
pts$base7_code <- as.integer(base7_vals)

# ---------- (Optional) quick counts per class for QC
counts <- pts |>
  st_drop_geometry() |>
  count(base7_code, name = "N") |>
  arrange(base7_code)

# Save counts table
dir.create(dirname(out_counts), showWarnings = FALSE, recursive = TRUE)
write_csv(counts, out_counts)

# ---------- Write updated points (new file to keep provenance clean)
st_write(pts, out_points, delete_layer = TRUE, quiet = TRUE)
cat("✅ Wrote updated points with base7_code:", out_points, "\n")
cat("📊 Counts saved to:", out_counts, "\n")
