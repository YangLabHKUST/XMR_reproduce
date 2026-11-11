#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else getwd()
default_results_dir <- normalizePath(file.path(script_dir, "..", "results"))
default_output_dir <- normalizePath(file.path(script_dir, "..", "figures"))

results_dir <- if (length(args) >= 1L) args[[1L]] else default_results_dir
output_dir <- if (length(args) >= 2L) args[[2L]] else default_output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

xmr_filename <- "SHBG_T2D_XMR_EUR_selected_FEMA_P_target_clump_noEffectHetQC.tsv"
xmr_path <- file.path(results_dir, "crossPopIntercept0", xmr_filename)
camera_results_subdir <- "CAMeRa_results_clump10000kb_V2"
camera_source_file <- file.path(
  camera_results_subdir, "SHBG_T2D_CAMeRa_all_results.tsv"
)
camera_iv_source_file <- file.path(
  camera_results_subdir, "threshold_5e-08", "FEMA_IV",
  "CAMeRa_harmonised_SHBG_T2D.tsv"
)

mrbee_filename <- "SHBG_T2D_MRBEE_TL_all_results.tsv"
mrbee_path <- file.path(results_dir, "MRBEE_TL_FEMA_source_results_clump10000kb", mrbee_filename)
mrbee_source_file <- if (dirname(mrbee_path) == results_dir) {
  basename(mrbee_path)
} else {
  file.path(basename(dirname(mrbee_path)), basename(mrbee_path))
}

paths <- c(
  XMR = xmr_path,
  CAMERA = file.path(results_dir, camera_source_file),
  CAMERA_IVs = file.path(results_dir, camera_iv_source_file),
  MRBEE = mrbee_path,
  TEMR = file.path(results_dir, "TEMR_results_clump10000kb", "SHBG_T2D_TEMR_all_results.tsv"),
  MRAPSS = file.path(results_dir, "MRAPSS_results_clump10000kb", "SHBG_T2D_MRAPSS_all_results.tsv")
)

missing_files <- paths[!file.exists(paths)]
if (length(missing_files) > 0L) {
  stop("Missing result file(s):\n", paste(missing_files, collapse = "\n"))
}

read_result <- function(path, required_columns) {
  x <- fread(path, na.strings = c("", "NA", "NaN", "Inf", "-Inf"))
  missing_columns <- setdiff(required_columns, names(x))
  if (length(missing_columns) > 0L) {
    stop(
      basename(path),
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
  x
}

canonical_threshold <- function(x) {
  value <- suppressWarnings(as.numeric(x))
  fifelse(
    is.finite(value) & abs(log10(value) - log10(5e-5)) < 1e-6,
    "5e-05",
    fifelse(
      is.finite(value) & abs(log10(value) - log10(5e-8)) < 1e-6,
      "5e-08",
      NA_character_
    )
  )
}

standard_frame <- function(
    method_group,
    population,
    population_label,
    threshold,
    beta,
    se,
    p,
    iv_total,
    iv_effective,
    status,
    error,
    source_file,
    selection_note,
    se_origin = "Reported",
    row_order = NULL) {
  n <- length(population)
  if (is.null(row_order)) row_order <- seq_len(n)
  data.table(
    method_group = rep(method_group, length.out = n),
    population = toupper(trimws(as.character(population))),
    population_label = as.character(population_label),
    row_order = as.numeric(row_order),
    threshold = canonical_threshold(threshold),
    beta = suppressWarnings(as.numeric(beta)),
    se = suppressWarnings(as.numeric(se)),
    p = suppressWarnings(as.numeric(p)),
    iv_total = suppressWarnings(as.numeric(iv_total)),
    iv_effective = suppressWarnings(as.numeric(iv_effective)),
    status = tolower(trimws(as.character(status))),
    error = as.character(error),
    source_file = rep(source_file, length.out = n),
    selection_note = rep(selection_note, length.out = n),
    se_origin = rep(se_origin, length.out = n)
  )
}

placeholder_rows <- function(method_group, populations, threshold, source_file, note) {
  standard_frame(
    method_group = method_group,
    population = populations,
    population_label = populations,
    threshold = rep(threshold, length(populations)),
    beta = rep(NA_real_, length(populations)),
    se = rep(NA_real_, length(populations)),
    p = rep(NA_real_, length(populations)),
    iv_total = rep(NA_real_, length(populations)),
    iv_effective = rep(NA_real_, length(populations)),
    status = rep("not_applicable", length(populations)),
    error = rep(note, length(populations)),
    source_file = source_file,
    selection_note = note,
    se_origin = "Not available"
  )
}

population_levels <- c("EUR", "EAS", "AFR", "CSA")
xmr_thresholds <- "5e-05"
threshold_levels <- c("5e-05", "5e-08")

# XMR: selected EUR exposure hits, two-population signed FEMA P ranking,
# target-population LD clumping at 5e-05. Use nvalid for effective IVs.
xmr_raw <- read_result(
  paths[["XMR"]],
  c(
    "selection_mode", "clump_priority_mode", "population", "threshold",
    "beta", "se", "p", "nIV", "nvalid", "status", "error",
    "effect_heterogeneity_qc", "clump_kb", "clump_r2"
  )
)
xmr_source_file <- if (dirname(paths[["XMR"]]) == results_dir) {
  basename(paths[["XMR"]])
} else {
  file.path(
    basename(dirname(paths[["XMR"]])),
    basename(paths[["XMR"]])
  )
}
xmr_raw[, threshold_key := canonical_threshold(threshold)]
xmr_raw <- xmr_raw[
  threshold_key %chin% xmr_thresholds &
    selection_mode == "EUR_SELECTED_FEMA_P_CLUMP" &
    clump_priority_mode == "FEMA_P"
]
expected_xmr <- CJ(
  population = c("EAS", "AFR", "CSA"),
  threshold_key = xmr_thresholds
)
if (nrow(xmr_raw) != nrow(expected_xmr) ||
    !fsetequal(xmr_raw[, .(population, threshold_key)], expected_xmr)) {
  stop("Expected one XMR FEMA-P row per target at 5e-05")
}
if (any(xmr_raw$effect_heterogeneity_qc != FALSE, na.rm = TRUE) ||
    any(is.na(xmr_raw$effect_heterogeneity_qc)) ||
    any(xmr_raw$clump_kb != 10000, na.rm = TRUE) ||
    any(is.na(xmr_raw$clump_kb)) ||
    any(abs(xmr_raw$clump_r2 - 0.001) > 1e-12, na.rm = TRUE) ||
    any(is.na(xmr_raw$clump_r2))) {
  stop("XMR rows do not match no-effect-heterogeneity-QC 10 Mb/r2=0.001 results")
}
xmr <- standard_frame(
  method_group = "XMR",
  population = xmr_raw$population,
  population_label = xmr_raw$population,
  threshold = xmr_raw$threshold,
  beta = xmr_raw$beta,
  se = xmr_raw$se,
  p = xmr_raw$p,
  iv_total = xmr_raw$nIV,
  iv_effective = xmr_raw$nvalid,
  status = xmr_raw$status,
  error = xmr_raw$error,
  source_file = xmr_source_file,
  selection_note = paste0(
    "EUR exposure P below ", xmr_raw$threshold_key,
    "; cross-population intercept fixed at zero; no effect-heterogeneity QC; ",
    "two-population signed FEMA P ranking; target-population LD clump ",
    "(10 Mb, r2 = 0.001)"
  )
)
xmr <- rbindlist(list(
  placeholder_rows(
    "XMR", "EUR", NA_character_, xmr_source_file,
    "EUR is the auxiliary population; no EUR XMR target estimate"
  ),
  xmr
), use.names = TRUE)

# CAMERA: discovery-stage FEMA regional IVs at P < 5e-8. Count the IVs
# actually used per population, not the size of the shared candidate set.
camera_raw <- read_result(
  paths[["CAMERA"]],
  c(
    "pops", "threshold", "iv_mode", "post_fema_clump", "beta", "se",
    "p", "nIV", "n_harmonised_associations", "status", "error"
  )
)
camera_raw[, threshold_key := canonical_threshold(threshold)]
camera_raw <- camera_raw[
  threshold_key == "5e-08" &
    iv_mode == "FEMA_IV" &
    !post_fema_clump &
    pops %chin% population_levels
]
if (nrow(camera_raw) != 4L || !setequal(camera_raw$pops, population_levels)) {
  stop("Expected four CAMERA discovery/FEMA-IV rows at threshold 5e-08")
}
camera_harmonised <- read_result(
  paths[["CAMERA_IVs"]],
  c("SNP", "pops", "beta.x", "se.x", "beta.y", "se.y")
)
if (anyNA(camera_harmonised$SNP) || any(camera_harmonised$SNP == "") ||
    anyNA(camera_harmonised$pops) ||
    !setequal(camera_harmonised$pops, population_levels) ||
    anyDuplicated(camera_harmonised[, .(pops, SNP)]) ||
    camera_harmonised[, any(
      !is.finite(beta.x) | !is.finite(se.x) | se.x <= 0 |
        !is.finite(beta.y) | !is.finite(se.y) | se.y <= 0
    )]) {
  stop("CAMERA harmonised data must contain unique, valid IVs per population")
}
camera_iv_counts <- camera_harmonised[, .(nIV_used = uniqueN(SNP)), by = pops]
camera_raw[camera_iv_counts, on = "pops", nIV_used := i.nIV_used]
if (anyNA(camera_raw$nIV_used) ||
    any(camera_raw$nIV_used > camera_raw$nIV) ||
    anyNA(camera_raw$n_harmonised_associations) ||
    any(camera_raw$n_harmonised_associations != nrow(camera_harmonised))) {
  stop("CAMERA population-specific IV counts do not match the result metadata")
}
camera <- standard_frame(
  method_group = "CAMERA",
  population = camera_raw$pops,
  population_label = camera_raw$pops,
  threshold = camera_raw$threshold,
  beta = camera_raw$beta,
  se = camera_raw$se,
  p = camera_raw$p,
  iv_total = camera_raw$nIV_used,
  iv_effective = NA_real_,
  status = camera_raw$status,
  error = camera_raw$error,
  source_file = camera_source_file,
  selection_note = paste0(
    "Population-specific discovery hits; ancestry-matched discovery clump; ",
    "regional FEMA IV selection; no post-FEMA clump; ",
    "IV count is the number of unique SNPs used in each population's harmonised data"
  )
)
camera[, iv_count_source_file := camera_iv_source_file]

# MRBEE: combine one EUR source estimate with the EAS, AFR, and CSA
# target-specific estimates. The input repeats the same EUR result for each
# source-target pair, so verify that the repeated records are identical first.
mrbee_raw <- read_result(
  paths[["MRBEE"]],
  c(
    "method", "clump_mode", "population", "target_population", "threshold",
    "beta", "se", "p", "nIV", "iv_selection", "clump_panel", "status", "error"
  )
)
mrbee_raw[, threshold_key := canonical_threshold(threshold)]
mrbee_raw <- mrbee_raw[
  threshold_key == "5e-08" & clump_mode == "EUR_clump"
]

source_raw <- mrbee_raw[method == "MRBEE_source_specific"]
specific_raw <- mrbee_raw[method == "MRBEE_target_specific"]
tl_raw <- mrbee_raw[method == "MRBEE_TL_UV"]
if (nrow(source_raw) != 3L ||
    !setequal(source_raw$target_population, c("EAS", "AFR", "CSA"))) {
  stop("Expected three repeated MRBEE source rows at threshold 5e-08")
}
if (nrow(specific_raw) != 3L ||
    !setequal(specific_raw$population, c("EAS", "AFR", "CSA"))) {
  stop("Expected MRBEE-specific rows for EAS, AFR, and CSA")
}
if (nrow(tl_raw) != 3L || !setequal(tl_raw$population, c("EAS", "AFR", "CSA"))) {
  stop("Expected MRBEE-TL rows for EAS, AFR, and CSA")
}

source_unique <- unique(source_raw[, .(
  population, threshold, beta, se, p, nIV, iv_selection, clump_panel, status, error
)])
if (nrow(source_unique) != 1L || source_unique$population[[1L]] != "EUR") {
  stop("MRBEE source rows are not one identical EUR result")
}

mrbee_source <- standard_frame(
  method_group = "MRBEE",
  population = source_unique$population,
  population_label = source_unique$population,
  threshold = source_unique$threshold,
  beta = source_unique$beta,
  se = source_unique$se,
  p = source_unique$p,
  iv_total = source_unique$nIV,
  iv_effective = NA_real_,
  status = source_unique$status,
  error = source_unique$error,
  source_file = mrbee_source_file,
  selection_note = paste0(
    source_unique$iv_selection, "; IV threshold 5e-08; ",
    source_unique$clump_panel, "; ",
    "source-specific estimator"
  )
)

mrbee_specific <- standard_frame(
  method_group = "MRBEE",
  population = specific_raw$population,
  population_label = specific_raw$population,
  threshold = specific_raw$threshold,
  beta = specific_raw$beta,
  se = specific_raw$se,
  p = specific_raw$p,
  iv_total = specific_raw$nIV,
  iv_effective = NA_real_,
  status = specific_raw$status,
  error = specific_raw$error,
  source_file = mrbee_source_file,
  selection_note = paste0(
    "Four-population joint-chi-square FEMA P below 5e-08; EUR-LD shared clump; ",
    "target-specific estimator"
  )
)
mrbee <- rbindlist(list(mrbee_source, mrbee_specific), use.names = TRUE)

mrbee_tl <- standard_frame(
  method_group = "MRBEE-TL",
  population = tl_raw$population,
  population_label = tl_raw$population,
  threshold = tl_raw$threshold,
  beta = tl_raw$beta,
  se = tl_raw$se,
  p = tl_raw$p,
  iv_total = tl_raw$nIV,
  iv_effective = NA_real_,
  status = tl_raw$status,
  error = tl_raw$error,
  source_file = mrbee_source_file,
  selection_note = paste0(
    "Four-population joint-chi-square FEMA P below 5e-08; EUR-LD shared clump; ",
    "transfer-learning estimator"
  )
)
mrbee_tl <- rbindlist(list(
  placeholder_rows(
    "MRBEE-TL", "EUR", "5e-08", mrbee_source_file,
    "MRBEE transfer-learning estimator is not defined for the EUR auxiliary population"
  ),
  mrbee_tl
), use.names = TRUE)

# TEMR: use target-population results at P < 5e-8 and reconstruct SE from the
# reported two-sided P value, including machine-precision handling for P = 0.
temr_raw <- read_result(
  paths[["TEMR"]],
  c(
    "target_population", "threshold", "beta", "p", "nIV", "status", "error"
  )
)
temr_raw[, threshold_key := canonical_threshold(threshold)]
temr_raw <- temr_raw[threshold_key == "5e-08"]
if (nrow(temr_raw) != 3L ||
    !setequal(temr_raw$target_population, c("EAS", "AFR", "CSA"))) {
  stop("Expected TEMR rows for EAS, AFR, and CSA at threshold 5e-08")
}
temr_p_safe <- pmin(
  pmax(as.numeric(temr_raw$p), .Machine$double.xmin),
  1
)
temr_z <- qnorm(temr_p_safe / 2, lower.tail = FALSE)
temr_se <- abs(as.numeric(temr_raw$beta)) / temr_z
temr <- standard_frame(
  method_group = "TEMR",
  population = temr_raw$target_population,
  population_label = temr_raw$target_population,
  threshold = temr_raw$threshold,
  beta = temr_raw$beta,
  se = temr_se,
  p = temr_raw$p,
  iv_total = temr_raw$nIV,
  iv_effective = NA_real_,
  status = temr_raw$status,
  error = temr_raw$error,
  source_file = file.path("TEMR_results_clump10000kb", basename(paths[["TEMR"]])),
  selection_note = paste0(
    "Target or auxiliary exposure P below 5e-08; unsigned combined exposure P ",
    "ranking; target-population LD clump"
  ),
  se_origin = "Derived from two-sided P value"
)
temr <- rbindlist(list(
  placeholder_rows(
    "TEMR", "EUR", "5e-08", basename(paths[["TEMR"]]),
    "The new TEMR result file reports target populations only"
  ),
  temr
), use.names = TRUE)

# MR-APSS: target-population IV selection and clumping at P < 5e-5.
mrapss_raw <- read_result(
  paths[["MRAPSS"]],
  c(
    "population", "threshold", "beta", "se", "p", "nIV", "nvalid",
    "status", "error"
  )
)
mrapss_raw[, threshold_key := canonical_threshold(threshold)]
mrapss_raw <- mrapss_raw[threshold_key == "5e-05"]
if (nrow(mrapss_raw) != 4L || !setequal(mrapss_raw$population, population_levels)) {
  stop("Expected four MR-APSS population rows at threshold 5e-05")
}
mrapss <- standard_frame(
  method_group = "MR-APSS",
  population = mrapss_raw$population,
  population_label = mrapss_raw$population,
  threshold = mrapss_raw$threshold,
  beta = mrapss_raw$beta,
  se = mrapss_raw$se,
  p = mrapss_raw$p,
  iv_total = mrapss_raw$nIV,
  iv_effective = mrapss_raw$nvalid,
  status = mrapss_raw$status,
  error = mrapss_raw$error,
  source_file = file.path("MRAPSS_results_clump10000kb", basename(paths[["MRAPSS"]])),
  selection_note = paste0(
    "Target-population exposure P below 5e-05; target-population LD clump"
  )
)

plot_data <- rbindlist(
  list(
    xmr,
    camera,
    mrbee,
    mrbee_tl,
    temr,
    mrapss
  ),
  use.names = TRUE,
  fill = TRUE
)

method_levels <- c(
  "XMR",
  "CAMERA",
  "MRBEE",
  "MRBEE-TL",
  "TEMR",
  "MR-APSS"
)
plot_data[, method_order := match(method_group, method_levels)]
plot_data[, population_order := match(population, population_levels)]
plot_data[, threshold_order := match(threshold, threshold_levels)]
plot_data[is.na(threshold) & status == "not_applicable", threshold_order := 0L]
if (anyNA(plot_data$method_order) || anyNA(plot_data$population_order) ||
    anyNA(plot_data$threshold_order)) {
  stop("Unmapped method, population, or threshold order in selected plotting data")
}
setorder(plot_data, method_order, population_order, threshold_order)

plot_data[, valid := status == "ok" & is.finite(beta) & is.finite(se) & se > 0]
plot_data[, `:=`(
  ci_low = fifelse(valid, beta - 1.96 * se, NA_real_),
  ci_high = fifelse(valid, beta + 1.96 * se, NA_real_)
)]

format_p_value <- function(value, valid) {
  if (!isTRUE(valid) || !is.finite(value)) return("--")
  if (value <= 0 || value < 1e-300) return("<1e-300")
  if (value < 0.001) return(formatC(value, format = "e", digits = 2))
  sprintf("%.3f", value)
}

plot_data[, estimate_label := vapply(seq_len(.N), function(i) {
  if (!isTRUE(valid[[i]])) return("--")
  se_text <- if (se[[i]] < 0.001) {
    formatC(se[[i]], format = "e", digits = 2)
  } else {
    sprintf("%.3f", se[[i]])
  }
  if (se_origin[[i]] == "Derived from two-sided P value") {
    se_text <- paste0("~", se_text)
  }
  sprintf("%.3f (%s)", beta[[i]], se_text)
}, character(1L))]
plot_data[, threshold_label := fifelse(is.na(threshold), "--", threshold)]
plot_data[, p_label := vapply(
  seq_len(.N),
  function(i) format_p_value(p[[i]], valid[[i]]),
  character(1L)
)]
plot_data[, iv_label := vapply(seq_len(.N), function(i) {
  total <- iv_total[[i]]
  effective <- iv_effective[[i]]
  if (!is.finite(total)) return("--")
  total_label <- format(round(total), big.mark = ",", scientific = FALSE)
  if (is.finite(effective)) {
    return(paste0(
      total_label,
      " (",
      format(round(effective), big.mark = ",", scientific = FALSE),
      ")"
    ))
  }
  total_label
}, character(1L))]

# Method-first vertical layout: bold method header followed by available
# population rows. Missing estimators remain visible as -- placeholders.
plot_data[, y_raw := NA_real_]
method_headers <- data.table(
  method_group = method_levels,
  y_raw = NA_real_
)
cursor <- 0
for (m in seq_along(method_levels)) {
  method_value <- method_levels[[m]]
  method_headers[method_group == method_value, y_raw := -cursor]
  cursor <- cursor + 0.65
  idx <- which(plot_data$method_group == method_value)
  for (i in idx) {
    plot_data$y_raw[[i]] <- -cursor
    cursor <- cursor + 0.86
  }
  if (m < length(method_levels)) cursor <- cursor + 0.55
}

shift_y <- -min(c(plot_data$y_raw, method_headers$y_raw)) + 1
plot_data[, y := y_raw + shift_y]
method_headers[, y := y_raw + shift_y]

method_fills <- c(
  "XMR" = "#FCEFED",
  "CAMERA" = "#EDF5FA",
  "MRBEE" = "#EDF8F5",
  "MRBEE-TL" = "#FDF7EB",
  "TEMR" = "#F7F0FA",
  "MR-APSS" = "#F2F1F1"
)
# Match simulation-1017/results-V2/fromhxh/figures/plot_beta0_qq_method_thresholds.R.
# Use one MRBEE colour across source and target populations in this figure.
method_colours <- c(
  "XMR" = "#D73027",
  "CAMERA" = "#0072B2",
  "MRBEE" = "#009E73",
  "MRBEE-TL" = "#E69F00",
  "TEMR" = "#8E44AD",
  "MR-APSS" = "#332F2C"
)

method_groups <- split(
  plot_data,
  factor(plot_data$method_group, levels = method_levels)
)
bands <- rbindlist(lapply(method_groups, function(group) {
  data.table(
    method_group = group$method_group[[1L]],
    ymin = min(group$y) - 0.40,
    ymax = max(group$y) + 0.40,
    fill = unname(method_fills[group$method_group[[1L]]])
  )
}))

row_lines <- plot_data[, .(y = y - 0.43)]
separator_data <- rbindlist(lapply(method_groups[-length(method_groups)], function(group) {
  data.table(y = min(group$y) - 0.72)
}))

header_y <- max(method_headers$y) + 1.05
header_rule_y <- max(method_headers$y) + 0.43
forest_limits <- c(-0.55, 0.25)
forest_breaks <- seq(-0.5, 0.2, by = 0.1)
forest_width_scale <- 1.225
# Widen only the forest column while preserving the effect scale and other columns.
forest_x <- function(x) {
  forest_limits[[1L]] + (x - forest_limits[[1L]]) * forest_width_scale
}

x_layout <- list(
  xmin = -1.10,
  method = -1.04,
  population = -1.00,
  estimate = 0.515,
  p = 0.89,
  threshold = 1.13,
  iv = 1.44,
  xmax = 1.69
)

axis_y <- min(plot_data$y) - 0.86
axis_tick_y <- axis_y - 0.14
axis_text_y <- axis_y - 0.31
axis_title_y <- axis_y - 1.35
body_y_limits <- c(axis_title_y - 0.18, header_y + 0.40)

valid_data <- plot_data[valid == TRUE]
valid_data[, colour_key := method_group]
invalid_data <- plot_data[valid == FALSE]
grid_data <- data.table(x = forest_x(forest_breaks[forest_breaks != 0]))
tick_data <- data.table(x = forest_x(forest_breaks), label = sprintf("%.1f", forest_breaks))
text_scale <- 1.35

body <- ggplot() +
  geom_rect(
    data = bands,
    aes(ymin = ymin, ymax = ymax, fill = fill),
    xmin = x_layout$xmin,
    xmax = x_layout$xmax,
    colour = NA
  ) +
  scale_fill_identity() +
  geom_segment(
    data = row_lines,
    aes(x = x_layout$xmin, xend = x_layout$xmax, y = y, yend = y),
    colour = "white",
    linewidth = 0.30
  ) +
  geom_segment(
    data = separator_data,
    aes(x = x_layout$xmin, xend = x_layout$xmax, y = y, yend = y),
    colour = "#3D4349",
    linewidth = 0.58
  ) +
  annotate(
    "segment",
    x = x_layout$xmin,
    xend = x_layout$xmax,
    y = header_rule_y,
    yend = header_rule_y,
    colour = "#3D4349",
    linewidth = 0.65
  ) +
  geom_segment(
    data = grid_data,
    aes(
      x = x,
      xend = x,
      y = min(plot_data$y) - 0.40,
      yend = max(plot_data$y) + 0.40
    ),
    colour = "#D9DEE5",
    linewidth = 0.34
  ) +
  annotate(
    "segment",
    x = forest_x(0),
    xend = forest_x(0),
    y = min(plot_data$y) - 0.40,
    yend = max(plot_data$y) + 0.40,
    colour = "#3F4A54",
    linewidth = 0.72,
    linetype = "dashed"
  ) +
  geom_errorbar(
    data = valid_data,
    aes(
      y = y,
      xmin = forest_x(ci_low),
      xmax = forest_x(ci_high),
      colour = colour_key
    ),
    orientation = "y",
    width = 0.18,
    linewidth = 0.80
  ) +
  geom_point(
    data = valid_data,
    aes(x = forest_x(beta), y = y, colour = colour_key),
    size = 3.3,
    stroke = 0.25
  ) +
  geom_text(
    data = invalid_data,
    aes(x = forest_x(-0.12), y = y, label = "--"),
    size = 4.4 * text_scale,
    colour = "#6B7280"
  ) +
  annotate(
    "rect",
    xmin = forest_x(forest_limits[[1L]]),
    xmax = forest_x(forest_limits[[2L]]),
    ymin = min(plot_data$y) - 0.40,
    ymax = max(plot_data$y) + 0.40,
    fill = NA,
    colour = "#8A929B",
    linewidth = 0.52
  ) +
  geom_text(
    data = method_headers,
    aes(x = x_layout$method, y = y, label = method_group),
    hjust = 0,
    fontface = "bold",
    size = 4.85 * text_scale,
    colour = "#17212B"
  ) +
  geom_text(
    data = plot_data,
    aes(x = x_layout$population, y = y, label = population_label),
    hjust = 0,
    size = 4.25 * text_scale,
    colour = "#1F2933"
  ) +
  geom_text(
    data = plot_data,
    aes(x = x_layout$estimate, y = y, label = estimate_label),
    hjust = 0,
    size = 4.05 * text_scale,
    colour = "#111827"
  ) +
  geom_text(
    data = plot_data,
    aes(x = x_layout$p, y = y, label = p_label),
    hjust = 0,
    size = 4.05 * text_scale,
    colour = "#111827"
  ) +
  geom_text(
    data = plot_data,
    aes(x = x_layout$threshold, y = y, label = threshold_label),
    hjust = 0,
    size = 4.05 * text_scale,
    colour = "#111827"
  ) +
  geom_text(
    data = plot_data,
    aes(x = x_layout$iv, y = y, label = iv_label),
    hjust = 0,
    size = 4.05 * text_scale,
    colour = "#111827"
  ) +
  annotate(
    "text",
    x = x_layout$method,
    y = header_y,
    label = "Method / Population",
    hjust = 0,
    fontface = "bold",
    size = 4.9 * text_scale
  ) +
  annotate(
    "text",
    x = forest_x(mean(forest_limits)),
    y = header_y,
    label = "bold(hat(beta))~bold('(95% CI)')",
    parse = TRUE,
    hjust = 0.5,
    fontface = "bold",
    size = 4.9 * text_scale
  ) +
  annotate(
    "text",
    x = x_layout$estimate,
    y = header_y,
    label = "bold(hat(beta))~bold('(SE)')",
    parse = TRUE,
    hjust = 0,
    fontface = "bold",
    size = 4.9 * text_scale
  ) +
  annotate(
    "text",
    x = x_layout$p,
    y = header_y,
    label = "P-value",
    hjust = 0,
    fontface = "bold",
    size = 4.9 * text_scale
  ) +
  annotate(
    "text",
    x = x_layout$threshold,
    y = header_y,
    label = "IV threshold",
    hjust = 0,
    fontface = "bold",
    size = 4.9 * text_scale
  ) +
  annotate(
    "text",
    x = x_layout$iv,
    y = header_y,
    label = "#IV",
    hjust = 0,
    fontface = "bold",
    size = 4.9 * text_scale
  ) +
  annotate(
    "segment",
    x = forest_x(forest_limits[[1L]]),
    xend = forest_x(forest_limits[[2L]]),
    y = axis_y,
    yend = axis_y,
    colour = "#5F6872",
    linewidth = 0.50
  ) +
  geom_segment(
    data = tick_data,
    aes(x = x, xend = x, y = axis_y, yend = axis_tick_y),
    colour = "#5F6872",
    linewidth = 0.43
  ) +
  geom_text(
    data = tick_data,
    aes(x = x, y = axis_text_y, label = label),
    size = 4.05 * text_scale,
    colour = "#27313A"
  ) +
  annotate(
    "text",
    x = forest_x(mean(forest_limits)),
    y = axis_title_y,
    label = "Effect size",
    size = 4.7 * text_scale,
    colour = "#111827"
  ) +
  scale_colour_manual(values = method_colours, guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_y_continuous(expand = expansion(mult = 0.02)) +
  coord_cartesian(
    xlim = c(x_layout$xmin, x_layout$xmax),
    ylim = body_y_limits,
    clip = "off"
  ) +
  theme_void() +
  theme(plot.margin = margin(1, 3, 1, 3))

title_plot <- ggdraw() +
  draw_label(
    "SHBG-T2D",
    x = 0.5,
    hjust = 0.5,
    fontface = "bold",
    size = 17 * text_scale
  )

final_plot <- plot_grid(
  title_plot,
  body,
  ncol = 1,
  rel_heights = c(0.032, 1)
)

file_stem <- "SHBG_T2D_selected_methods_forest_by_method"
data_path <- file.path(output_dir, paste0(file_stem, "_data.tsv"))
pdf_path <- file.path(output_dir, paste0(file_stem, ".pdf"))
png_path <- file.path(output_dir, paste0(file_stem, ".png"))

export_data <- copy(plot_data)
setcolorder(export_data, c(
  "method_group", "population", "population_label", "threshold",
  "beta", "se", "ci_low", "ci_high", "p", "iv_total", "iv_effective",
  "estimate_label", "p_label", "threshold_label", "iv_label", "status", "error", "valid",
  "se_origin", "selection_note", "source_file", "method_order",
  "population_order", "threshold_order", "row_order", "y_raw", "y"
))
fwrite(export_data, data_path, sep = "\t", na = "NA")

input_manifest <- data.table(
  input = names(paths), path = normalizePath(unname(paths)),
  md5 = unname(tools::md5sum(paths)),
  rows = vapply(paths, function(p) nrow(fread(p)), integer(1L))
)
fwrite(input_manifest, file.path(output_dir, paste0(file_stem, "_inputs.tsv")), sep = "\t")
checks <- data.table(
  check = c("XMR version", "XMR target rows", "Forest rows", "Valid estimates",
            "Placeholder rows", "CAMERA actual IV counts", "TEMR uncertainty"),
  result = c("crossPopIntercept0", nrow(xmr_raw), nrow(plot_data),
             sum(plot_data$valid), sum(plot_data$status == "not_applicable"),
             paste(camera_iv_counts$pops, camera_iv_counts$nIV_used, collapse = "; "),
             "Approximate SE/CI reconstructed from P; P=0 floored at machine precision, not an exact SE")
)
stopifnot(nrow(xmr_raw) == 3L, nrow(plot_data) == 24L,
          sum(plot_data$valid) == 21L,
          all(plot_data[method_group == "XMR" & valid == TRUE, threshold] == "5e-05"),
          all(plot_data$p[plot_data$valid] >= 0 & plot_data$p[plot_data$valid] <= 1))
fwrite(checks, file.path(output_dir, paste0(file_stem, "_checks.tsv")), sep = "\t")
capture.output(sessionInfo(), file = file.path(output_dir, paste0(file_stem, "_sessionInfo.txt")))

figure_width <- 16.2
figure_height <- 17.6
ggsave(
  pdf_path,
  final_plot,
  width = figure_width,
  height = figure_height,
  units = "in",
  device = grDevices::pdf,
  useDingbats = FALSE,
  bg = "white"
)
ggsave(
  png_path,
  final_plot,
  width = figure_width,
  height = figure_height,
  units = "in",
  dpi = 320,
  bg = "white"
)

message("Saved data: ", data_path)
message("Saved PDF: ", pdf_path)
message("Saved PNG: ", png_path)
