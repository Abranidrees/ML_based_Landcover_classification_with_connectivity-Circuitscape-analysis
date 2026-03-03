#!/usr/bin/env Rscript
# =========================================================
# 40_change.R — Change analysis (tables + figures)
# Spec: Base-7; Rock (7) masked in summaries/plots
# Inputs (aligned CRS/res; extents may differ):
#   outputs/rasters/class_2004-2006.tif
#   outputs/rasters/class_2010-2012.tif
#   outputs/rasters/class_2016-2018.tif
#   outputs/rasters/class_2022-2024.tif
# Config:
#   config/colors_base7.csv  (columns: Code,Class,Hex)
# Outputs:
#   outputs/tables/transition_*.csv
#   outputs/tables/change_summary.csv
#   outputs/tables/change_qc.csv
#   outputs/figures/change_gains_losses.png
#   outputs/figures/change_net.png
#   outputs/figures/change_heatmap_*.png
# =========================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(purrr)
  library(grid)  # for unit()
})

set.seed(4401)

# ----------------------------
# 1) Paths & config
# ----------------------------
IN_RASTERS <- list(
  `2004_2006` = "outputs/rasters/class_2004-2006.tif",
  `2010_2012` = "outputs/rasters/class_2010-2012.tif",
  `2016_2018` = "outputs/rasters/class_2016-2018.tif",
  `2022_2024` = "outputs/rasters/class_2022-2024.tif"
)

DIR_TBL <- "outputs/tables"
DIR_FIG <- "outputs/figures"
dir.create(DIR_TBL, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

# OPTIONAL: AOI mask raster (preferred by client)
# - If it exists, it will be applied before paired-mask intersection.
# - Expected: 1 inside AOI, 0 outside OR NA outside.
AOI_MASK <- "outputs/rasters/aoi_mask.tif"  # optional, safe if missing

# Read class labels & palette from config
COLORS_CSV <- "config/colors_base7.csv"   # columns: Code,Class,Hex
stopifnot(file.exists(COLORS_CSV))
cls_raw <- read_csv(COLORS_CSV, show_col_types = FALSE)
names(cls_raw) <- tolower(names(cls_raw))
req_cols <- c("code","class","hex")
if (!all(req_cols %in% names(cls_raw))) {
  stop("colors_base7.csv must have columns: Code, Class, Hex (case-insensitive).")
}
cls_tbl <- cls_raw |>
  transmute(code = as.integer(code),
            name = as.character(class),
            hex  = as.character(hex)) |>
  arrange(code)

code_to_name <- setNames(cls_tbl$name, cls_tbl$code)
code_to_hex  <- setNames(cls_tbl$hex,  cls_tbl$code)

# ----------------------------
# 2) Helpers
# ----------------------------

# Mask Rock (7) from summaries/plots
mask_rock <- function(r) terra::ifel(r == 7, NA, r)

# Hectares per pixel (projected or lon/lat)
pixel_ha <- function(r) {
  crs_txt <- tryCatch(terra::crs(r, proj = TRUE), error = function(e) NA_character_)
  rs <- terra::res(r)
  message(sprintf("Raster CRS: %s", ifelse(is.na(crs_txt), "<unknown>", crs_txt)))
  message(sprintf("Pixel size (CRS units): %.6f × %.6f", rs[1], rs[2]))
  if (terra::is.lonlat(r)) {
    message("Geographic CRS detected (lon/lat). Using geodesic mean cell area.")
    cs <- terra::cellSize(r, unit = "m")
    mean_m2 <- as.numeric(terra::global(cs, "mean", na.rm = TRUE)[1,1])
    return(mean_m2 / 10000.0)
  } else {
    return((rs[1] * rs[2]) / 10000.0)
  }
}

# Short-scale labels (version-safe)
lab_short <- function(accuracy = 1) {
  if (utils::packageVersion("scales") >= "1.2.0") {
    scales::label_number(accuracy = accuracy, scale_cut = scales::cut_short_scale())
  } else {
    function(x) scales::number(x, accuracy = accuracy, big.mark = ",")
  }
}

# Align 'to' to 'from' grid if needed (CRS/res must already match)
harmonize_to_template <- function(r_from, r_to) {
  if (crs(r_from) != crs(r_to)) {
    stop("CRS mismatch between rasters in a change pair.")
  }
  if (!all(res(r_from) == res(r_to))) {
    stop(sprintf("Resolution mismatch: from=%.6f×%.6f vs to=%.6f×%.6f",
                 res(r_from)[1], res(r_from)[2], res(r_to)[1], res(r_to)[2]))
  }
  if (!all(ext(r_from) == ext(r_to)) || ncol(r_from) != ncol(r_to) || nrow(r_from) != nrow(r_to)) {
    message("Harmonizing grids: resampling 'to' to match 'from' (nearest).")
    r_to <- terra::resample(r_to, r_from, method = "near")
  }
  r_to
}

# Apply AOI mask (preferred by client); safe no-op if missing
apply_aoi_if_present <- function(r) {
  if (!is.null(AOI_MASK) && file.exists(AOI_MASK)) {
    aoi <- rast(AOI_MASK)
    if (crs(aoi) != crs(r)) stop("AOI mask CRS does not match raster CRS.")
    if (!all(res(aoi) == res(r)) || !all(ext(aoi) == ext(r)) ||
        ncol(aoi) != ncol(r) || nrow(aoi) != nrow(r)) {
      message("AOI mask grid differs; harmonizing AOI to raster grid (nearest).")
      aoi <- terra::resample(aoi, r, method = "near")
      aoi <- terra::crop(aoi, r)
    }
    
    # If AOI looks like 0/1, mask out 0. Otherwise, use NA footprint.
    rng <- tryCatch(terra::global(aoi, "range", na.rm = TRUE), error = function(e) NULL)
    if (!is.null(rng)) {
      mn <- as.numeric(rng[1,1]); mx <- as.numeric(rng[1,2])
      if (!is.na(mn) && !is.na(mx) && mn >= 0 && mx <= 1) {
        r <- terra::mask(r, aoi, maskvalues = 0)
      } else {
        r <- terra::mask(r, aoi)
      }
    } else {
      r <- terra::mask(r, aoi)
    }
  }
  r
}

# Crosstab → transition (ha), with paired valid-pixel intersection (client requirement)
compute_transition <- function(r_from, r_to, pair_label) {
  r_to <- harmonize_to_template(r_from, r_to)
  
  # Apply AOI first (preferred), then rock masking
  r_from <- apply_aoi_if_present(r_from)
  r_to   <- apply_aoi_if_present(r_to)
  
  r1 <- mask_rock(r_from)
  r2 <- mask_rock(r_to)
  
  # Restrict analysis to pixels valid in BOTH rasters (publication-safe)
  common <- (!is.na(r1)) & (!is.na(r2))
  common01 <- terra::ifel(common, 1, 0)
  
  r1 <- terra::mask(r1, common01, maskvalues = 0)
  r2 <- terra::mask(r2, common01, maskvalues = 0)
  
  # QC: paired analysis area
  px_ha <- pixel_ha(r1)
  paired_cells <- as.numeric(terra::global(common01, "sum", na.rm = TRUE)[1,1])
  paired_area_ha <- paired_cells * px_ha
  message(sprintf("Paired analysis area (%s): %.2f ha", pair_label, paired_area_ha))
  
  ct <- terra::crosstab(c(r1, r2), useNA = FALSE, long = TRUE)
  stopifnot(ncol(ct) >= 3)
  names(ct)[1:3] <- c("from", "to", "npx")
  ct$from <- as.integer(ct$from); ct$to <- as.integer(ct$to); ct$npx <- as.numeric(ct$npx)
  
  ct$area_ha <- ct$npx * px_ha
  ct$pair    <- pair_label
  
  ct <- ct[order(ct$from, ct$to), c("pair","from","to","npx","area_ha")]
  
  # Store QC on the object for later change_qc.csv
  attr(ct, "paired_cells") <- paired_cells
  attr(ct, "paired_area_ha") <- paired_area_ha
  attr(ct, "px_ha") <- px_ha
  
  ct
}

# ✅ Client replacement: gross gain/loss + persistence + net (publication-safe)
summarise_change <- function(ct_list) {
  all_ct <- dplyr::bind_rows(ct_list)
  
  # totals at each time (still useful)
  area_from <- all_ct |>
    dplyr::group_by(pair, code = from) |>
    dplyr::summarise(area_from_ha = sum(area_ha), .groups="drop")
  
  area_to <- all_ct |>
    dplyr::group_by(pair, code = to) |>
    dplyr::summarise(area_to_ha = sum(area_ha), .groups="drop")
  
  # persistence (from == to)
  pers <- all_ct |>
    dplyr::filter(from == to) |>
    dplyr::group_by(pair, code = from) |>
    dplyr::summarise(persist_ha = sum(area_ha), .groups="drop")
  
  # gross gain / gross loss (exclude persistence)
  gain <- all_ct |>
    dplyr::filter(from != to) |>
    dplyr::group_by(pair, code = to) |>
    dplyr::summarise(gain_ha = sum(area_ha), .groups="drop")
  
  loss <- all_ct |>
    dplyr::filter(from != to) |>
    dplyr::group_by(pair, code = from) |>
    dplyr::summarise(loss_ha = sum(area_ha), .groups="drop")
  
  out <- purrr::reduce(
    list(area_from, area_to, gain, loss, pers),
    dplyr::full_join,
    by = c("pair","code")
  ) |>
    dplyr::mutate(
      dplyr::across(c(area_from_ha, area_to_ha, gain_ha, loss_ha, persist_ha),
                    ~tidyr::replace_na(.x, 0)),
      net_ha = area_to_ha - area_from_ha
    ) |>
    dplyr::arrange(pair, code)
  
  out
}

pair_nice <- function(p) {
  p |>
    str_replace_all("_", " ") |>
    str_replace("2004 2006", "2004–06") |>
    str_replace("2010 2012", "2010–12") |>
    str_replace("2016 2018", "2016–18") |>
    str_replace("2022 2024", "2022–24")
}

# ----------------------------
# 3) Load class rasters; log CRS/px size once
# ----------------------------
lc <- lapply(IN_RASTERS, rast)
invisible(pixel_ha(lc[[1]]))

# ----------------------------
# 4) Change pairs
# ----------------------------
pairs <- list(
  `2004_2006__2010_2012` = c("2004_2006", "2010_2012"),
  `2010_2012__2016_2018` = c("2010_2012", "2016_2018"),
  `2016_2018__2022_2024` = c("2016_2018", "2022_2024")
)

# ----------------------------
# 5) Transition tables per pair → CSV
# ----------------------------
transitions <- list()
pair_qc <- list()

for (lab in names(pairs)) {
  from_nm <- pairs[[lab]][1]; to_nm <- pairs[[lab]][2]
  ct <- compute_transition(lc[[from_nm]], lc[[to_nm]], pair_label = lab)
  
  # store QC for this pair
  pair_qc[[lab]] <- tibble(
    pair = lab,
    paired_cells = as.numeric(attr(ct, "paired_cells")),
    paired_area_ha = as.numeric(attr(ct, "paired_area_ha"))
  )
  
  ct_out <- ct |>
    mutate(from_name = code_to_name[as.character(from)],
           to_name   = code_to_name[as.character(to)]) |>
    select(pair, from, from_name, to, to_name, npx, area_ha)
  
  write_csv(ct_out, file.path(DIR_TBL, paste0("transition_", lab, ".csv")))
  transitions[[lab]] <- ct
}

pair_qc_tbl <- bind_rows(pair_qc)

# ----------------------------
# 6) Change summary (gross gain/loss + persistence + net) → CSV + QC
# ----------------------------
chg <- summarise_change(transitions) |>
  mutate(class_name = code_to_name[as.character(code)],
         pair_label = pair_nice(pair)) |>
  select(pair, pair_label, code, class_name,
         area_from_ha, area_to_ha, persist_ha, gain_ha, loss_ha, net_ha)

write_csv(chg, file.path(DIR_TBL, "change_summary.csv"))

qc <- chg |>
  group_by(pair, pair_label) |>
  summarise(
    sum_net_ha = sum(net_ha),
    total_gain_ha = sum(gain_ha),
    total_loss_ha = sum(loss_ha),
    total_persist_ha = sum(persist_ha),
    total_from_ha = sum(area_from_ha),
    total_to_ha   = sum(area_to_ha),
    .groups = "drop"
  ) |>
  left_join(pair_qc_tbl, by = "pair")

write_csv(qc, file.path(DIR_TBL, "change_qc.csv"))

# ----------------------------
# 7) Figures
# ----------------------------

# Palette
pal_vec <- code_to_hex[as.character(sort(unique(chg$code)))]

# 7a) Gains & Losses (NOW gross gain/loss; publication-safe)
bars_long <- chg |>
  select(pair_label, code, class_name, gain_ha, loss_ha) |>
  pivot_longer(c(gain_ha, loss_ha), names_to = "type", values_to = "area_ha") |>
  mutate(type = ifelse(type == "gain_ha", "Gain", "Loss"),
         code = as.integer(code))

p_bars <- ggplot(bars_long,
                 aes(x = reorder(class_name, code),
                     y = area_ha,
                     fill = factor(code, levels = sort(unique(code))))) +
  geom_col(position = "stack") +
  facet_grid(type ~ pair_label, scales = "free_y") +
  scale_y_continuous(labels = lab_short(accuracy = 1)) +
  scale_fill_manual(values = pal_vec,
                    breaks = sort(unique(bars_long$code)),
                    labels = code_to_name[as.character(sort(unique(bars_long$code)))],
                    name = "Class") +
  labs(x = "Class", y = "Area (ha)",
       title = "Gross Gains and Gross Losses by Class (ha)",
       subtitle = "Persistence excluded from gains/losses; Rock (7) masked; paired-valid footprint enforced") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(DIR_FIG, "change_gains_losses.png"),
       p_bars, width = 12, height = 7.2, dpi = 350)

# 7b) Net change (headline)
net_fig <- chg |>
  select(pair_label, code, class_name, net_ha) |>
  mutate(code = as.integer(code))

p_net <- ggplot(net_fig,
                aes(x = reorder(class_name, code),
                    y = net_ha,
                    fill = factor(code, levels = sort(unique(code))))) +
  geom_hline(yintercept = 0) +
  geom_col() +
  scale_y_continuous(labels = lab_short(accuracy = 1)) +
  scale_fill_manual(values = pal_vec, guide = "none") +
  facet_wrap(~ pair_label, ncol = 1, scales = "free_y") +
  labs(x = "Class", y = "Net change (ha)",
       title = "Net change by class (ha)",
       subtitle = "Positive = net gain; Negative = net loss; Rock (7) masked; paired-valid footprint enforced") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(DIR_FIG, "change_net.png"),
       p_net, width = 8.5, height = 8.8, dpi = 350)

# 7c) Transition heatmaps (fixed class order + fixed legend)
pair_ids <- c("2004_2006__2010_2012",
              "2010_2012__2016_2018",
              "2016_2018__2022_2024")

# Desired display order (must match names in colors_base7.csv)
desired_order <- c("Forest","Agrocanopy","Paddy","GrassScrub","Built","Water")

# Map codes we will keep (drop Rock)
keep_tbl <- cls_tbl |>
  filter(name %in% desired_order) |>
  arrange(factor(name, levels = desired_order)) |>
  select(code, name)

# One global max so all heatmaps share the same colorbar
global_max <- max(vapply(pair_ids, function(lab) {
  fp <- file.path(DIR_TBL, paste0("transition_", lab, ".csv"))
  if (!file.exists(fp)) return(0)
  suppressMessages(read_csv(fp, show_col_types = FALSE)$area_ha) |> max(na.rm = TRUE)
}, numeric(1)), na.rm = TRUE)

for (lab in pair_ids) {
  ct_path <- file.path(DIR_TBL, paste0("transition_", lab, ".csv"))
  if (!file.exists(ct_path)) {
    warning("Missing transition CSV: ", ct_path, " — skipping.")
    next
  }
  
  ct <- suppressMessages(read_csv(ct_path, show_col_types = FALSE)) |>
    # keep only our 6 classes
    semi_join(keep_tbl, by = c("from" = "code")) |>
    semi_join(keep_tbl, by = c("to"   = "code")) |>
    mutate(area_ha = as.numeric(area_ha))
  
  # Build a full 6x6 grid (even if some combos are missing)
  full_grid <- expand.grid(from = keep_tbl$code,
                           to   = keep_tbl$code,
                           KEEP.OUT.ATTRS = FALSE)
  
  hm <- full_grid |>
    left_join(ct, by = c("from","to")) |>
    mutate(area_ha   = replace_na(area_ha, 0),
           from_name = factor(code_to_name[as.character(from)],
                              levels = desired_order),
           to_name   = factor(code_to_name[as.character(to)],
                              levels = desired_order),
           pair_label = dplyr::case_when(
             lab == "2004_2006__2010_2012" ~ "2004–06 → 2010–12",
             lab == "2010_2012__2016_2018" ~ "2010–12 → 2016–18",
             TRUE                          ~ "2016–18 → 2022–24"
           )) |>
    filter(!is.na(from_name), !is.na(to_name))
  
  p_hm <- ggplot(hm, aes(x = to_name, y = from_name, fill = area_ha)) +
    geom_tile() +
    coord_equal() +
    scale_fill_gradientn(
      colours = c("#f7fbff","#6baed6","#08306b"),
      trans   = scales::pseudo_log_trans(sigma = 1),
      limits  = c(0, global_max),
      oob     = scales::squish,
      labels  = lab_short(accuracy = 1),
      name    = "ha",
      guide   = guide_colorbar(
        direction      = "vertical",
        title.position = "top",
        label.position = "right",
        barheight      = unit(10, "cm"),
        barwidth       = unit(0.5, "cm"),
        ticks = TRUE
      )
    ) +
    labs(x = "To class", y = "From class",
         title = "Transition matrix (ha, pseudo-log scale)",
         subtitle = paste0("Pair: ", unique(hm$pair_label), " — Rock (7) masked; paired-valid footprint enforced")) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position   = "right",
      legend.box.margin = margin(4, 6, 4, 4),
      legend.title      = element_text(size = 10, margin = margin(b = 4)),
      legend.text       = element_text(size = 9),
      axis.text.x       = element_text(angle = 35, hjust = 1),
      panel.grid.minor  = element_blank(),
      plot.margin       = margin(6, 10, 6, 6)
    )
  
  ggsave(file.path(DIR_FIG, paste0("change_heatmap_", lab, ".png")),
         p_hm, width = 10.6, height = 7.2, dpi = 350)
}

message("Tables → ", DIR_TBL)
message("Figures → ", DIR_FIG)
