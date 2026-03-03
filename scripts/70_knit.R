#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rmarkdown)
  library(yaml)
  library(glue)
})

# -----------------------------
# Project root 
# -----------------------------
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA)
proj_root <- if (!is.na(this_file)) normalizePath(file.path(dirname(this_file), "..")) else getwd()
setwd(proj_root)


# -----------------------------
# Paths
# -----------------------------
rmd_in  <- file.path("reports", "methods_results.Rmd")
out_dir <- "reports"
log_dir <- "logs"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rmd_in)) {
  stop(glue("❌ Missing report input: {rmd_in}"))
}

# -----------------------------
# Render HTML 
# -----------------------------
message("🧶 Rendering methods_results.Rmd -> HTML ...")
html_out <- rmarkdown::render(
  input       = rmd_in,
  output_dir  = out_dir,
  output_file = "methods_results.html",
  quiet       = FALSE,
  envir       = new.env(parent = globalenv())
)

message(glue("✅ HTML written: {html_out}"))


# -----------------------------
# Save session info for reproducibility
# -----------------------------
sess_fp <- file.path(log_dir, "sessionInfo.txt")
writeLines(c(
  paste0("Date: ", Sys.time()),
  paste0("Working directory: ", getwd()),
  "",
  capture.output(sessionInfo())
), sess_fp)

message(glue("🧾 sessionInfo saved: {sess_fp}"))
message("🎉 Done.")
