# ============================================================
# Accuracy assessment (standalone) — with buffer, area-adjusted OA,
# F1/BA/Macro-F1, and compact colorblind-safe figures
# Keeps original CSV outputs unchanged and adds new ones.
# ============================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(readr); library(glue); library(tibble)
  library(yardstick); library(stringr); library(ggplot2); library(tidyr)
  library(viridisLite)
})

# ----------------------- USER KNOBS (edit these) -----------------------
WINDOW <- "2016-2018" # 2004-2006; 2010-2012; 2016-2018; 2022-2024
STACK  <- "allbands"   # "allbands" or "fullextended"
buf_m  <- 500          # drop valid points within this distance of any train point
# targets + bootstrap
TARGET_OA <- 0.85      # used for ✓/✗ call-out
# Optionally set per-class targets as named vectors, else leave NULL
TARGET_UA <- NULL      # e.g., c(Forest=0.80, Water=0.90)
TARGET_PA <- NULL

BOOT_B    <- 1000      # bootstrap draws
BOOT_SEED <- 42

# top confusions to report
TOP_K_CONFUSIONS <- 5   # how many off-diagonal confusions to list
# ----------------------- DERIVED KEYS ----------------------------------
STACK_TAG <- switch(tolower(STACK),
                    "allbands"     = "all_bands",
                    "fullextended" = "full_extended",
                    stop("STACK must be 'allbands' or 'fullextended'"))

# ----------------------- PATHS -----------------------------------------
PRED_RASTER <- glue("outputs/rasters/class_{WINDOW}.tif")
TRAIN_GPKG_CANDIDATES <- c(
  glue("outputs/vectors/train_{WINDOW}.gpkg"),
  glue("outputs/vectors/train_{gsub('_','-',WINDOW)}.gpkg"),
  glue("outputs/vectors/train_{gsub('-','_',WINDOW)}.gpkg")
)
OUT_TABLES  <- "outputs/tables"
OUT_FIGS    <- "outputs/figures"
OUT_LOGS    <- "outputs/logs"
COLORS_CSV  <- "config/colors_base7.csv"
dir.create(OUT_TABLES, TRUE, TRUE)
dir.create(OUT_FIGS, TRUE, TRUE)
dir.create(OUT_LOGS, TRUE, TRUE)

# ----------------------- CLASSES ---------------------------------------
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

# Load colorblind-safe class colors (optional but recommended)
# Expect columns like: code,name,hex  (hex = "#RRGGBB")
class_colors <- tryCatch({
  cc <- readr::read_csv(COLORS_CSV, show_col_types = FALSE)
  # harmonize column names
  if (!"code" %in% names(cc)) {
    if ("base7_code" %in% names(cc)) cc <- dplyr::rename(cc, code = base7_code)
  }
  if (!"hex" %in% names(cc)) {
    if ("color" %in% names(cc)) cc <- dplyr::rename(cc, hex = color)
  }
  dplyr::select(cc, code, hex) %>%
    mutate(code = as.integer(code)) %>%
    right_join(classes, by = "code")
}, error = function(e) {
  warning("Could not read ", COLORS_CSV, " — will fallback to viridis/dark greys.")
  classes %>% mutate(hex = NA_character_)
})

# ----------------------- LOAD PREDICTION -------------------------------
pred_num <- terra::rast(PRED_RASTER)
crs_ref  <- terra::crs(pred_num)

# ----------------------- PICK VALIDATION LAYER -------------------------
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

# ----------------------- ENFORCE TRAIN–VALID BUFFER --------------------
train_pts <- pts[pts$split == "train", ]
valid_pts <- pts[pts$split == "valid", ]
if (nrow(valid_pts) == 0) stop("No validation rows (split == 'valid') before buffering.")

# Drop valid points within buffer of any train point
if (nrow(train_pts) > 0 && nrow(valid_pts) > 0) {
  nb <- sf::st_is_within_distance(valid_pts, train_pts, dist = buf_m)
  keep_idx <- lengths(nb) == 0
  n_dropped <- sum(!keep_idx)
  valid_pts <- valid_pts[keep_idx, ]
  message("Buffered validation: dropped ", n_dropped, " valid points within ", buf_m, " m of train pts. ",
          "Remaining valid: ", nrow(valid_pts))
}
stopifnot(nrow(valid_pts) > 0)

# ----------------------- EXTRACT PREDICTIONS AT VALIDATION POINTS ------
pred_at_val <- terra::extract(pred_num, terra::vect(valid_pts))[, 2]

acc_df <- sf::st_drop_geometry(valid_pts) %>%
  dplyr::transmute(
    base7_code = factor(base7_code, levels = classes$code),
    pred       = factor(pred_at_val,  levels = classes$code)
  ) %>%
  dplyr::filter(!is.na(pred))


# ----------------------- HELPERS ---------------------------------------
safe_div <- function(n, d) ifelse(d > 0, n / d, NA_real_)

compute_confusion_and_metrics <- function(df) {
  tab <- as.matrix(table(
    Pred = factor(df$pred,       levels = classes$code),
    Ref  = factor(df$base7_code, levels = classes$code)
  ))
  diagv   <- diag(tab)
  row_tot <- rowSums(tab)        # predicted totals  -> UA denominator
  col_tot <- colSums(tab)        # reference totals  -> PA denominator
  UA <- safe_div(diagv, row_tot) # precision
  PA <- safe_div(diagv, col_tot) # recall
  OA <- sum(diagv) / sum(tab)
  kappa_val <- tryCatch(yardstick::kap(df, base7_code, pred)$.estimate, error = function(e) NA_real_)
  list(tab=tab, UA=UA, PA=PA, OA=OA, kappa=kappa_val, row_tot=row_tot, col_tot=col_tot)
}

# ----------------------- CONFUSION + METRICS ----------
tab <- as.matrix(table(
  Pred = factor(acc_df$pred,       levels = classes$code),
  Ref  = factor(acc_df$base7_code, levels = classes$code)
))

safe_div <- function(n, d) ifelse(d > 0, n / d, NA_real_)
diagv   <- diag(tab)
row_tot <- rowSums(tab)        # predicted totals  -> UA denominator
col_tot <- colSums(tab)        # reference totals  -> PA denominator
UA <- safe_div(diagv, row_tot) # User's accuracy (precision)
PA <- safe_div(diagv, col_tot) # Producer's accuracy (recall)

# Tidy by-class UA/PA + counts
by_class <- dplyr::bind_rows(
  tibble::tibble(code = classes$code, metric = "PA", value = PA),
  tibble::tibble(code = classes$code, metric = "UA", value = UA)
)
counts <- tibble::tibble(code = classes$code, N_valid = as.numeric(col_tot))
acc_tidy <- dplyr::full_join(by_class, counts, by = "code") %>%
  dplyr::left_join(classes, by = "code") %>%
  dplyr::arrange(code, metric)
readr::write_csv(acc_tidy, file.path(OUT_TABLES, glue::glue("accuracy_{WINDOW}.csv")))

# Raw confusion (long)
cm_tbl <- as.data.frame(as.table(tab))    # Pred, Ref, Freq
names(cm_tbl) <- c("Pred", "Ref", "n")
readr::write_csv(cm_tbl, file.path(OUT_TABLES, glue::glue("accuracy_{WINDOW}_confusion_raw.csv")))

# Paper-style matrix with OA and Kappa
mat_df   <- as.data.frame.matrix(tab)
ref_cols <- paste0("Ref_", classes$code)
stopifnot(ncol(mat_df) == length(ref_cols))
colnames(mat_df) <- ref_cols

OA        <- sum(diagv) / sum(tab)
kappa_val <- yardstick::kap(acc_df, base7_code, pred)$.estimate

acc_matrix <- tibble::tibble(
  Pred_code = classes$code,
  Pred_name = classes$name
) %>%
  dplyr::bind_cols(mat_df) %>%
  dplyr::mutate(
    Row_total          = row_tot,
    Users_accuracy     = round(UA, 4),
    Producers_accuracy = round(PA, 4),
    Overall_Accuracy   = NA_real_,
    Kappa              = NA_real_
  )

col_tot_row <- tibble::tibble(
  Pred_code = NA_integer_,
  Pred_name = "Column_totals"
) %>%
  dplyr::bind_cols(
    tibble::as_tibble_row(stats::setNames(as.list(as.numeric(col_tot)), ref_cols))
  ) %>%
  dplyr::mutate(
    Row_total          = sum(row_tot),
    Users_accuracy     = NA_real_,
    Producers_accuracy = NA_real_,
    Overall_Accuracy   = NA_real_,
    Kappa              = NA_real_
  )

overall_row <- tibble::tibble(
  Pred_code = NA_integer_,
  Pred_name = "Overall",
  !!!stats::setNames(as.list(rep(NA_real_, length(ref_cols))), ref_cols),
  Row_total          = sum(tab),
  Users_accuracy     = NA_real_,
  Producers_accuracy = NA_real_,
  Overall_Accuracy   = round(OA, 4),
  Kappa              = round(kappa_val, 4)
)

acc_matrix_full <- dplyr::bind_rows(acc_matrix, col_tot_row, overall_row)
#readr::write_csv(acc_matrix_full, file.path(OUT_TABLES, glue::glue("accuracy_matrix_{WINDOW}.csv")))

# ----------------------- EXTRA METRICS: F1, BA, Macro-F1  --------------
prec <- UA
rec  <- PA
F1   <- 2 * (prec * rec) / (prec + rec)
BA   <- mean(rec, na.rm = TRUE)
MacroF1 <- mean(F1, na.rm = TRUE)

acc_extended <- tibble::tibble(
  code  = classes$code,
  name  = classes$name,
  N_pred = as.numeric(row_tot),
  N_ref  = as.numeric(col_tot),
  UA    = round(prec, 4),
  PA    = round(rec, 4),
  F1    = round(F1, 4)
)

acc_extended_overall <- tibble::tibble(
  code = NA_integer_, name = "Overall",
  N_pred = sum(row_tot), N_ref = sum(col_tot),
  UA = NA_real_, PA = NA_real_, F1 = NA_real_
) %>%
  mutate(
    OA = round(OA, 4),
    Kappa = round(kappa_val, 4),
    BA = round(BA, 4),
    MacroF1 = round(MacroF1, 4)
  )

# write extra metrics ( originals untouched)
readr::write_csv(
  bind_rows(
    acc_extended %>% mutate(OA = NA_real_, Kappa = NA_real_, BA = NA_real_, MacroF1 = NA_real_),
    acc_extended_overall
  ),
  file.path(OUT_TABLES, glue::glue("accuracy_matrix_{WINDOW}.csv"))
)

# ----------------------- AREA-ADJUSTED ACCURACY (Olofsson) -----------------------
# Robust pixel counts per MAP class
pixel_counts <- vapply(classes$code, function(k) {
  v <- try(terra::global(pred_num == k, "sum", na.rm = TRUE)[1,1], silent = TRUE)
  if (inherits(v, "try-error") || is.na(v)) 0 else as.numeric(v)
}, numeric(1))

total_pixels <- sum(pixel_counts)

area_tbl <- tibble::tibble(
  code         = classes$code,
  name         = classes$name,
  pixel_count  = as.integer(pixel_counts),
  pixel_prop_raw = if (total_pixels > 0) pixel_counts / total_pixels else rep(0, length(pixel_counts))
)

# weights w_i aligned to MAP classes (rows of tab) — and normalize (idempotent if already sums to 1)
w_i <- setNames(area_tbl$pixel_prop_raw, as.character(area_tbl$code))
ws <- sum(w_i)
if (ws > 0) w_i <- w_i / ws

# Row-normalize confusion: p_ij = P(Ref=j | Map=i)
p_ij <- sweep(tab, 1, row_tot, "/")

# Per-class UA and area-adjusted OA
UA_row <- diag(p_ij)
OA_adj <- sum(w_i * UA_row, na.rm = TRUE)

# Detailed CSV (your preferred format)
area_debug <- tibble::tibble(
  class_code       = area_tbl$code,
  class_name       = area_tbl$name,
  pixel_count      = as.integer(area_tbl$pixel_count),
  pixel_prop_raw   = round(area_tbl$pixel_prop_raw, 8),
  pixel_prop_final = round(as.numeric(w_i), 8),
  row_total_points = as.integer(row_tot),
  col_total_points = as.integer(col_tot),
  UA_row           = round(UA_row, 8),
  contrib_to_OA    = round(as.numeric(w_i) * UA_row, 8)
)

area_debug_footer <- tibble::tibble(
  class_code       = NA_integer_,
  class_name       = "TOTAL/CHECKS",
  pixel_count      = sum(area_debug$pixel_count, na.rm = TRUE),
  pixel_prop_raw   = round(sum(area_debug$pixel_prop_raw,   na.rm = TRUE), 8),
  pixel_prop_final = round(sum(area_debug$pixel_prop_final, na.rm = TRUE), 8),
  row_total_points = sum(area_debug$row_total_points, na.rm = TRUE),
  col_total_points = sum(area_debug$col_total_points, na.rm = TRUE),
  UA_row           = NA_real_,
  contrib_to_OA    = round(sum(area_debug$contrib_to_OA,    na.rm = TRUE), 8)  # equals OA_adj
)

readr::write_csv(
  dplyr::bind_rows(area_debug, area_debug_footer),
  file.path(OUT_TABLES, glue::glue("accuracy_area_adjusted_{WINDOW}.csv"))
)


# ----------------------- CONFUSION "NOTES" (TOP ERRORS) ---------
# pct_of_ref = n / column total for that Ref
notes_tbl <- as.data.frame(as.table(tab)) %>%
  tibble::as_tibble(.name_repair = "minimal") %>%
  dplyr::rename(Pred = 1, Ref = 2, n = 3) %>%
  mutate(
    ref_code  = as.integer(as.character(Ref)),
    pred_code = as.integer(as.character(Pred)),
    ref_name  = classes$name[match(ref_code, classes$code)],
    pred_name = classes$name[match(pred_code, classes$code)],
    ref_col_tot = col_tot[match(ref_code, classes$code)],
    pct_of_ref = safe_div(n, ref_col_tot)
  ) %>%
  filter(!is.na(pct_of_ref), pred_code != ref_code) %>%
  arrange(desc(pct_of_ref)) %>%
  slice_head(n = TOP_K_CONFUSIONS) %>%
  mutate(
    window     = WINDOW,
    quick_note = "common mix; check boundary/texture; add training tiles"
  ) %>%
  select(window, ref_code, ref_name, pred_code, pred_name, n, pct_of_ref, quick_note)

readr::write_csv(notes_tbl, file.path(OUT_TABLES, glue("accuracy_notes_{WINDOW}.csv")))

# ----------------------- TARGETS CALL-OUT -----------------------
target_ok <- !is.na(TARGET_OA) && is.finite(TARGET_OA) && (OA >= TARGET_OA)
target_mark <- ifelse(target_ok, "\u2713", "\u2717") # ✓ / ✗

# Optional: per-class target checks (reported to log only)
check_targets <- function(vals, targets_named) {
  if (is.null(targets_named) || length(targets_named) == 0) return(NA_character_)
  nm <- names(targets_named)
  cls_idx <- match(nm, classes$name)
  got <- vals[cls_idx]
  ok <- mapply(function(g,t) is.finite(g) && g >= t, got, targets_named)
  paste0(nm, ":", ifelse(ok, "✓","✗"), collapse = " ")
}

UA_named <- setNames(UA, classes$name)
PA_named <- setNames(PA, classes$name)
ua_targets_summary <- check_targets(UA_named, TARGET_UA)
pa_targets_summary <- check_targets(PA_named, TARGET_PA)

# -----------------------  BOOTSTRAP CIs --------------------------
set.seed(BOOT_SEED)
N <- nrow(acc_df)
boot_metrics <- function(idx) {
  dfb <- acc_df[idx, , drop = FALSE]
  cm_b <- compute_confusion_and_metrics(dfb)
  # Return OA + per-class UA/PA
  tibble::tibble(
    metric = c("OA", paste0("UA_", classes$name), paste0("PA_", classes$name)),
    value  = c(cm_b$OA, cm_b$UA, cm_b$PA)
  )
}

boot_list <- replicate(BOOT_B, sample.int(N, N, replace = TRUE), simplify = FALSE)
boot_df <- purrr::map_dfr(boot_list, boot_metrics, .id = "boot_id")

# summarize CIs
summ_ci <- boot_df %>%
  group_by(metric) %>%
  summarise(
    est  = mean(value, na.rm = TRUE),
    lo95 = quantile(value, 0.025, na.rm = TRUE, names = FALSE),
    hi95 = quantile(value, 0.975, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  mutate(B = BOOT_B)

# reshape to requested tidy format: metric (OA|UA|PA), code, name, est, lo95, hi95, B
to_row <- function(prefix, values_named) {
  tibble::tibble(
    metric = prefix,
    code   = classes$code,
    name   = classes$name,
    est    = as.numeric(values_named)
  )
}
# point estimates
ci_point <- bind_rows(
  tibble::tibble(metric="OA", code=NA_integer_, name="Overall", est=OA),
  to_row("UA", UA_named),
  to_row("PA", PA_named)
)

# merge with CI bands
extract_ci <- function(prefix, nm) {
  summ_ci %>% filter(metric == nm) %>% select(est_ci=est, lo95, hi95, B)
}
# Build UA and PA CI rows by matching names
ua_ci_rows <- purrr::map_dfr(classes$name, function(nm){
  row <- extract_ci("UA", paste0("UA_", nm))
  tibble::tibble(metric="UA", code=classes$code[classes$name==nm], name=nm,
                 est_ci=row$est_ci, lo95=row$lo95, hi95=row$hi95, B=row$B)
})
pa_ci_rows <- purrr::map_dfr(classes$name, function(nm){
  row <- extract_ci("PA", paste0("PA_", nm))
  tibble::tibble(metric="PA", code=classes$code[classes$name==nm], name=nm,
                 est_ci=row$est_ci, lo95=row$lo95, hi95=row$hi95, B=row$B)
})
oa_ci_row <- summ_ci %>% filter(metric=="OA") %>%
  transmute(metric="OA", code=NA_integer_, name="Overall",
            est_ci=est, lo95, hi95, B)

# join point est with CI est (keep point as 'est')
acc_ci <- bind_rows(
  ci_point %>% left_join(bind_rows(oa_ci_row, ua_ci_rows, pa_ci_rows),
                         by=c("metric","code","name"))
)

readr::write_csv(acc_ci, file.path(OUT_TABLES, glue("accuracy_ci_{WINDOW}.csv")))

# ----------------------- CONFUSION % BY REF (CSV) ---------------
cm_long <- as.data.frame(as.table(tab)) |>
  tibble::as_tibble(.name_repair = "minimal") |>
  dplyr::rename(Pred = 1, Ref = 2, n = 3)

col_sum_df <- cm_long |>
  dplyr::group_by(Ref) |>
  dplyr::summarise(coln = sum(n), .groups = "drop")

cm_pct_by_ref <- cm_long |>
  dplyr::left_join(col_sum_df, by = "Ref") |>
  dplyr::mutate(pct_by_ref = ifelse(coln > 0, n/coln, NA_real_)) |>
  dplyr::select(Pred, Ref, n, pct_by_ref)

readr::write_csv(cm_pct_by_ref, file.path(OUT_TABLES, glue("accuracy_confusion_pct_by_ref_{WINDOW}.csv")))

# ----------------------- FIGURES (compact, colorblind-safe) ------------
colors_raw <- readr::read_csv(COLORS_CSV, show_col_types = FALSE) %>%
  dplyr::rename_with(tolower) %>%
  dplyr::mutate(code = as.integer(code))

# Palette resolution + sanity check
pal_resolved <- classes %>%
  left_join(select(colors_raw, code, hex), by="code") %>%
  mutate(
    source = ifelse(!is.na(hex) & nzchar(hex), "config", "fallback")
  )
# Fallback greys only where missing
if (any(is.na(pal_resolved$hex) | pal_resolved$hex=="")) {
  missing_idx <- which(is.na(pal_resolved$hex) | pal_resolved$hex=="")
  pal_resolved$hex[missing_idx] <- scales::grey_pal()(length(missing_idx))
  warning("Missing hex in palette for: ",
          paste(pal_resolved$name[missing_idx], collapse=", "),
          " — using greys as fallback.")
}
# Write resolved palette for traceability
readr::write_csv(pal_resolved %>% select(code, name, hex, source),
                 file.path(OUT_TABLES, glue("colors_resolved_{WINDOW}.csv")))

# Join only hex; keep 'name' coming from acc_tidy
acc_plot_df <- acc_tidy %>%
  dplyr::left_join(dplyr::select(pal_resolved, code, hex), by = "code") %>%
  dplyr::mutate(
    name  = dplyr::coalesce(name, classes$name[match(code, classes$code)]),
    label = factor(name, levels = classes$name),
    metric = factor(metric, levels = c("PA","UA"))
  )

pal_vec <- setNames(pal_resolved$hex, pal_resolved$name)

subtitle_txt <- glue::glue(
  "OA={round(OA,3)} | κ={round(kappa_val,3)} | BA={round(BA,3)} | Macro-F1={round(MacroF1,3)} | OA(adj)={round(OA_adj,3)} | target_OA:{ifelse(target_ok,'✓','✗')}"
)

# BAR CHART: PA solid, UA transparent
p_bar <- ggplot(acc_plot_df, aes(x = label, y = value)) +
  geom_col(data = dplyr::filter(acc_plot_df, metric == "PA"),
           aes(fill = label),
           position = position_dodge(width = 0.8), width = 0.7, alpha = 0.95) +
  geom_col(data = dplyr::filter(acc_plot_df, metric == "UA"),
           aes(fill = label),
           position = position_dodge(width = 0.8), width = 0.7, alpha = 0.55) +
  scale_fill_manual(values = pal_vec[levels(acc_plot_df$label)], guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(
    x = NULL, y = "Accuracy",
    title = glue::glue("Per-class PA (solid) & UA (transparent) — {WINDOW}"),
    subtitle = subtitle_txt
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  filename = file.path(OUT_FIGS, glue::glue("accuracy_{WINDOW}.png")),
  plot = p_bar, width = 8, height = 4.5, dpi = 400
)

# HEATMAP (column-normalized)
cm_long_plot <- cm_pct_by_ref %>%
  mutate(
    Pred_name = classes$name[match(as.integer(as.character(Pred)), classes$code)],
    Ref_name  = classes$name[match(as.integer(as.character(Ref)),  classes$code)],
    Pred_lab  = factor(Pred_name, levels = classes$name),
    Ref_lab   = factor(Ref_name,  levels = classes$name)
  )

p_heat <- ggplot(cm_long_plot, aes(x = Ref_lab, y = Pred_lab, fill = pct_by_ref)) +
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(pct_by_ref), "", scales::percent(pct_by_ref, accuracy = 0.1))),
            size = 3) +
  scale_fill_viridis_c(option = "C", limits = c(0,1), na.value = "white") +
  labs(x = "Reference", y = "Predicted",
       title = glue::glue("Confusion Matrix (column-normalized) — {WINDOW}"),
       subtitle = "Cell labels: % of reference class mapped as row class") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  filename = file.path(OUT_FIGS, glue::glue("confusion_{WINDOW}.png")),
  plot = p_heat, width = 7.5, height = 6.5, dpi = 400
)

# ----------------------- REPRODUCIBILITY LOG --------------------
log_path <- file.path(OUT_LOGS, glue("accuracy_run_{WINDOW}.txt"))

log_lines <- c(
  glue("window: {WINDOW}"),
  glue("stack: {STACK_TAG}"),
  glue("buffer_m: {buf_m}"),
  glue("valid_points_kept: {nrow(acc_df)}"),
  glue("OA: {round(OA,4)}"),
  glue("kappa: {round(kappa_val,4)}"),
  glue("BA: {round(BA,4)}"),
  glue("MacroF1: {round(MacroF1,4)}"),
  glue("OA_adjusted: {round(OA_adj,4)}"),
  glue("target_OA_threshold: {TARGET_OA}"),
  glue("target_OA_pass: {ifelse(target_ok,'yes','no')}"),
  glue("UA_targets: {ifelse(is.na(ua_targets_summary), '(none)', ua_targets_summary)}"),
  glue("PA_targets: {ifelse(is.na(pa_targets_summary), '(none)', pa_targets_summary)}"),
  glue("bootstrap_B: {BOOT_B}"),
  glue("bootstrap_seed: {BOOT_SEED}"),
  "---- sessionInfo ----",
  capture.output(utils::sessionInfo())
)

# Direct write; no extra files
writeLines(log_lines, con = log_path, useBytes = TRUE)



# ----------------------- CONSOLE SUMMARY -------------------------------
cat(glue::glue(
  "📊 wrote accuracy tables for {WINDOW} ({STACK_TAG})\n",
  "➕ wrote extended metrics + area-adjusted + CIs + notes + pct-by-ref CSVs\n",
  "🎨 palette resolved → colors_resolved_{WINDOW}.csv\n",
  "🖼️ saved figures: {file.path(OUT_FIGS, glue('accuracy_{WINDOW}.png'))}, ",
  "{file.path(OUT_FIGS, glue('confusion_{WINDOW}.png'))}\n",
  "🎯 target OA {TARGET_OA}: {ifelse(target_ok,'✓','✗')}\n",
  "🗒️ log: {log_path}\n"
))