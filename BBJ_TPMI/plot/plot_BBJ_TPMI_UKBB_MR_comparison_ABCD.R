#!/usr/bin/env Rscript

# BBJ/TPMI: crossPopIntercept0 XMR without effect-heterogeneity QC; UKBB: MR-APSS.
# Read the three result files directly; do not use cached plot data.
# A: significant associations; B: overlap; C: two BBJ examples; D: method discovery counts.
# The figure summarizes cross-cohort discoveries and method comparisons.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else getwd()
default_results_dir <- normalizePath(file.path(script_dir, "..", "results"))
default_output_dir <- normalizePath(file.path(script_dir, "..", "figures"))

results_dir <- if (length(args) >= 1L) args[[1L]] else default_results_dir
output_dir <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  default_output_dir
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  BBJ = file.path(results_dir, "masterClump5e-05", "BBJ", "crossPopIntercept0", "XMR_EAS_BBJ_EURselected_two_thresholds_noEffectHetQC.tsv"),
  TPMI = file.path(results_dir, "masterClump5e-05", "TPMI", "crossPopIntercept0", "XMR_EAS_TPMI_EURselected_two_thresholds_noEffectHetQC.tsv"),
  EUR = file.path(results_dir, "EUR", "MRAPSS_EUR_exposure_PanUKB_15outcomes_5e-05.tsv")
)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop("Missing required result file(s):\n", paste(missing_files, collapse = "\n"))
}

long_data <- rbindlist(lapply(names(input_files), function(cohort_name) {
  path <- input_files[[cohort_name]]
  dat <- fread(path, na.strings = c("", "NA", "NaN", "Inf", "-Inf"))
  if (cohort_name == "EUR") {
    if (!all(c("beta_se", "beta_pvalue") %in% names(dat))) {
      stop("EUR results require beta_se and beta_pvalue: ", path)
    }
    setnames(dat, c("beta_se", "beta_pvalue"), c("se", "p"))
  }
  if (cohort_name != "EUR") dat <- dat[abs(as.numeric(threshold) - 5e-5) < 1e-12]
  required <- c(
    "analysis_id", "method", "exposure", "exposure_name", "outcome",
    "outcome_name", "threshold", "beta", "se", "p", "n_IV", "pi0", "status"
  )
  if (cohort_name != "EUR") required <- c(required, "nvalid")
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop(basename(path), " is missing: ", paste(missing, collapse = ", "))
  }
  dat[, `:=`(
    cohort = cohort_name,
    source_file = normalizePath(path),
    source_row = seq_len(.N),
    exposure = as.character(exposure),
    outcome = as.character(outcome),
    exposure_name = trimws(as.character(exposure_name)),
    outcome_name = trimws(as.character(outcome_name)),
    status = tolower(trimws(as.character(status))),
    threshold = as.numeric(threshold),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p),
    n_IV = as.numeric(n_IV),
    pi0 = as.numeric(pi0)
  )]
  if (cohort_name == "EUR") {
    # MR-APSS reports the expected valid-IV count as pi0 * nrow(MRdat).
    dat[, `:=`(nvalid = n_IV * pi0, nvalid_source = "n_IV * pi0")]
  } else {
    dat[, `:=`(nvalid = as.numeric(nvalid), nvalid_source = "reported nvalid")]
  }
  if (any(!is.finite(dat$threshold)) ||
      any(abs(log10(dat$threshold) - log10(5e-5)) > 1e-6)) {
    stop(basename(path), " contains an IV threshold other than 5e-05")
  }
  if (anyNA(dat$exposure_name) || anyNA(dat$outcome_name) ||
      any(!nzchar(dat$exposure_name)) || any(!nzchar(dat$outcome_name)) ||
      dat[, uniqueN(.SD), .SDcols = c("exposure_name", "outcome_name")] != nrow(dat)) {
    stop(basename(path), " has missing or duplicated trait pairs")
  }
  dat[, .(
    cohort, method, analysis_id, exposure, exposure_name, outcome, outcome_name,
    threshold, beta, se, p, n_IV, pi0, nvalid, nvalid_source,
    status, source_file, source_row
  )]
}))

expected_cohorts <- names(input_files)
bbj_pairs <- long_data[cohort == "BBJ", .(exposure_name, outcome_name)]
for (cohort_name in c("TPMI", "EUR")) {
  if (!fsetequal(
    bbj_pairs, long_data[cohort == cohort_name, .(exposure_name, outcome_name)]
  )) {
    stop("BBJ, TPMI, and EUR must contain the same trait pairs")
  }
}
exposure_levels <- unique(long_data[cohort == "BBJ", exposure_name])
outcome_levels <- unique(long_data[cohort == "BBJ", outcome_name])
long_data[, `:=`(
  exposure_label = trimws(gsub("\\s*\\([^()]*\\)", "", exposure_name, perl = TRUE)),
  outcome_label = trimws(gsub("\\s*\\([^()]*\\)", "", outcome_name, perl = TRUE)),
  valid_result = status == "ok" & is.finite(beta) & is.finite(se) & se > 0 &
    is.finite(p) & p >= 0 & p <= 1,
  exposure_index = length(exposure_levels) - match(exposure_name, exposure_levels) + 1L,
  outcome_index = match(outcome_name, outcome_levels),
  effect_heterogeneity_qc = fifelse(cohort == "EUR", NA, FALSE)
)]
if (anyNA(long_data$valid_result) || any(!long_data$valid_result)) {
  stop("A result has failed or has invalid beta/SE/P; inspect input files before plotting")
}


long_data[, n_tests := sum(valid_result), by = cohort]
long_data[, bonferroni_threshold := 0.05 / n_tests]
long_data[, significant_bonferroni :=
  valid_result & is.finite(p) & p < bonferroni_threshold
]

test_counts <- unique(long_data[, .(cohort, n_tests, bonferroni_threshold)])
if (test_counts[, uniqueN(n_tests)] != 1L) {
  stop("Cohorts have different numbers of valid tests")
}
n_tests <- test_counts$n_tests[[1L]]
bonferroni_alpha <- 0.05 / n_tests

# Panel A: Bonferroni-significant associations in all three cohorts.
exposure_map <- unique(long_data[, .(exposure_name, exposure_label, exposure_index)])
outcome_map <- unique(long_data[, .(outcome_name, outcome_label, outcome_index)])
if (exposure_map[, .N, by = exposure_index][, any(N != 1L)] ||
    outcome_map[, .N, by = outcome_index][, any(N != 1L)]) {
  stop("Trait names do not map uniquely to matrix coordinates")
}
setorder(exposure_map, exposure_index)
setorder(outcome_map, outcome_index)

exposure_abbreviations <- c(
  "Body mass index" = "BMI", "Hemoglobin" = "Hb",
  "Diastolic blood pressure" = "DBP", "Alanine aminotransferase" = "ALT",
  "Glucose" = "Glu", "HbA1c" = "HbA1c", "LDL-cholesterol" = "LDL-C"
)
# Match outcome abbreviations used by the BBJ/TPMI analysis scripts.
outcome_abbreviations <- c(
  "Breast cancer" = "BC", "Prostate cancer" = "PrC",
  "Iron deficiency anemia" = "IDA", "Uterine fibroid" = "UF",
  "Epilepsy" = "Ep", "Cataract" = "Cat", "Angina pectoris" = "Ang",
  "Myocardial infarction" = "MI", "Peripheral arterial disease" = "PAD",
  "Unstable angina pectoris" = "UAP", "Asthma" = "As",
  "Chronic obstructive pulmonary disease" = "COPD",
  "Cholelithiasis" = "CL", "Type 2 diabetes" = "T2D",
  "Rheumatoid arthritis" = "RA"
)
exposure_map[, abbreviation := unname(exposure_abbreviations[exposure_label])]
outcome_map[, abbreviation := unname(outcome_abbreviations[outcome_label])]
stopifnot(!anyNA(exposure_map$abbreviation), !anyNA(outcome_map$abbreviation),
          uniqueN(exposure_map$abbreviation) == 7L,
          uniqueN(outcome_map$abbreviation) == 15L)

cohort_x_offset <- c("BBJ" = -0.19, "TPMI" = 0.19, "EUR" = 0)
cohort_y_offset <- c("BBJ" = 0.18, "TPMI" = 0.18, "EUR" = -0.18)
# Cohort colours are shared across A/C and distinct from simulation method colours.
cohort_colours <- c(
  "BBJ" = "#2A9DAD",
  "TPMI" = "#9A6B3F",
  "EUR" = "#8F9B45"
)
venn_colours <- setNames(
  unname(cohort_colours[c("BBJ", "TPMI", "EUR")]),
  c("BBJ", "TPMI", "PanUKB")
)

panel_a_data <- long_data[
  cohort %chin% expected_cohorts & significant_bonferroni == TRUE
]
if (!nrow(panel_a_data)) {
  stop("No Bonferroni-significant associations")
}
panel_a_data[, `:=`(
  direction = factor(
    fifelse(beta >= 0, "Positive", "Negative"),
    levels = c("Positive", "Negative")
  ),
  abs_beta = abs(beta),
  x_plot = outcome_index + unname(cohort_x_offset[cohort]),
  y_plot = exposure_index + unname(cohort_y_offset[cohort])
)]

grid_data <- CJ(
  outcome_index = sort(unique(outcome_map$outcome_index)),
  exposure_index = sort(unique(exposure_map$exposure_index))
)

max_abs_beta <- max(panel_a_data$abs_beta, na.rm = TRUE)
effect_size_breaks <- pretty(c(0, max_abs_beta), n = 4)
effect_size_breaks <- effect_size_breaks[
  effect_size_breaks > 0 & effect_size_breaks <= max_abs_beta * 1.001
]
if (!length(effect_size_breaks)) effect_size_breaks <- max_abs_beta

panel_a <- ggplot() +
  geom_tile(
    data = grid_data,
    aes(x = outcome_index, y = exposure_index),
    width = 0.96,
    height = 0.96,
    fill = "#FBFCFD",
    colour = "#D9E0E6",
    linewidth = 0.34
  ) +
  geom_point(
    data = panel_a_data,
    aes(
      x = x_plot,
      y = y_plot,
      fill = cohort,
      shape = direction,
      size = abs_beta
    ),
    colour = "#26323B",
    stroke = 0.34,
    alpha = 0.95
  ) +
  scale_fill_manual(
    name = "Cohort",
    values = cohort_colours,
    breaks = c("BBJ", "TPMI", "EUR"),
    labels = c("BBJ", "TPMI", "UKBB")
  ) +
  scale_shape_manual(
    name = "Effect direction",
    values = c("Positive" = 24, "Negative" = 25),
    labels = c("Positive", "Negative")
  ) +
  scale_size_area(
    name = "|Effect size|",
    max_size = 9.5,
    breaks = effect_size_breaks,
    limits = c(0, max_abs_beta * 1.015),
    labels = label_number(accuracy = 0.01)
  ) +
  scale_x_continuous(
    name = "Outcome",
    breaks = outcome_map$outcome_index,
    labels = outcome_map$abbreviation,
    limits = range(outcome_map$outcome_index) + c(-0.55, 0.55),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    name = "Exposure",
    breaks = exposure_map$exposure_index,
    labels = exposure_map$abbreviation,
    limits = range(exposure_map$exposure_index) + c(-0.55, 0.55),
    expand = expansion(mult = 0)
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(
    title = NULL
  ) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(shape = 24, size = 4.5, alpha = 1)
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(fill = "#6B7280", size = 4.5, alpha = 1)
    ),
    size = guide_legend(
      order = 3,
      override.aes = list(shape = 24, fill = "#6B7280", alpha = 1)
    )
  ) +
  theme_minimal(base_size = 14.5, base_family = "sans") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 25,
      colour = "#202A33",
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(
      size = 14.4,
      colour = "#4B5965",
      margin = margin(b = 6)
    ),
    axis.title = element_text(face = "bold", size = 25, colour = "#26323B"),
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 27)),
    axis.text.x = element_text(
      angle = 52,
      hjust = 1,
      vjust = 1,
      size = 20,
      colour = "#26323B"
    ),
    axis.text.y = element_text(size = 20, colour = "#26323B"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "#6B7785", linewidth = 0.5),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.justification = "center",
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 20),
    legend.key.height = grid::unit(0.6, "lines"),
    legend.spacing.x = grid::unit(0.25, "lines"),
    legend.margin = margin(t = 3),
    plot.margin = margin(6, 8, 5, 5)
  )

# Panel B: area-proportional Bonferroni Venn diagram.
significance_wide <- dcast(
  long_data,
  exposure_name + outcome_name ~ cohort,
  value.var = "significant_bonferroni"
)
significance_wide[, PanUKB := EUR]

region_counts <- data.table(
  region = c(
    "BBJ only",
    "TPMI only",
    "PanUKB only",
    "BBJ + TPMI only",
    "BBJ + PanUKB only",
    "TPMI + PanUKB only",
    "BBJ + TPMI + PanUKB",
    "None"
  ),
  n_pairs = c(
    significance_wide[, sum(BBJ & !TPMI & !PanUKB)],
    significance_wide[, sum(!BBJ & TPMI & !PanUKB)],
    significance_wide[, sum(!BBJ & !TPMI & PanUKB)],
    significance_wide[, sum(BBJ & TPMI & !PanUKB)],
    significance_wide[, sum(BBJ & !TPMI & PanUKB)],
    significance_wide[, sum(!BBJ & TPMI & PanUKB)],
    significance_wide[, sum(BBJ & TPMI & PanUKB)],
    significance_wide[, sum(!BBJ & !TPMI & !PanUKB)]
  )
)

circle_overlap_area <- function(distance, radius1, radius2) {
  if (distance >= radius1 + radius2) return(0)
  if (distance <= abs(radius1 - radius2)) {
    return(pi * min(radius1, radius2)^2)
  }
  angle1 <- acos(
    (distance^2 + radius1^2 - radius2^2) / (2 * distance * radius1)
  )
  angle2 <- acos(
    (distance^2 + radius2^2 - radius1^2) / (2 * distance * radius2)
  )
  lens_term <- sqrt(
    (-distance + radius1 + radius2) *
      (distance + radius1 - radius2) *
      (distance - radius1 + radius2) *
      (distance + radius1 + radius2)
  )
  radius1^2 * angle1 + radius2^2 * angle2 - 0.5 * lens_term
}

distance_for_overlap <- function(radius1, radius2, target_area) {
  maximum_area <- pi * min(radius1, radius2)^2
  if (target_area <= 0) return(radius1 + radius2)
  if (target_area >= maximum_area) return(abs(radius1 - radius2))
  uniroot(
    function(distance) {
      circle_overlap_area(distance, radius1, radius2) - target_area
    },
    interval = c(
      abs(radius1 - radius2) + .Machine$double.eps^0.5,
      radius1 + radius2 - .Machine$double.eps^0.5
    )
  )$root
}

# Set column widths for the association matrix and cohort-overlap panels.
original_column_widths <- 24 * c(1.82, 1) / 2.82
column_widths <- original_column_widths * c(1, 1.20)
canvas_width <- sum(column_widths)
venn_margin_inches <- 16 / 72.27
venn_width_ratio <- (column_widths[2] - venn_margin_inches) /
  (original_column_widths[2] - venn_margin_inches)
# Scale ellipse height to the panel viewport.
ellipse_y_scale <- 0.80 * 1.20 / venn_width_ratio
ellipse_polygon <- function(
    group, center_x, center_y, radius, draw_order, n = 720L) {
  angle <- seq(0, 2 * pi, length.out = n)
  data.table(
    group = group,
    draw_order = draw_order,
    angle = angle,
    x = center_x + radius * cos(angle),
    y = center_y + ellipse_y_scale * radius * sin(angle)
  )
}

bbj_total <- significance_wide[, sum(BBJ)]
tpmi_total <- significance_wide[, sum(TPMI)]
panukb_total <- significance_wide[, sum(PanUKB)]
triple_total <- significance_wide[, sum(BBJ & TPMI & PanUKB)]
bbj_panukb_total <- significance_wide[, sum(BBJ & PanUKB)]
bbj_tpmi_total <- significance_wide[, sum(BBJ & TPMI)]
tpmi_panukb_total <- significance_wide[, sum(TPMI & PanUKB)]

if (tpmi_panukb_total != tpmi_total || bbj_tpmi_total != triple_total) {
  stop("Observed Venn regions do not support the area-proportional layout")
}

radius_bbj <- sqrt(bbj_total / pi)
radius_tpmi <- sqrt(tpmi_total / pi)
radius_panukb <- sqrt(panukb_total / pi)
center_bbj <- 0
center_panukb <- distance_for_overlap(
  radius_bbj,
  radius_panukb,
  bbj_panukb_total
)
center_tpmi <- distance_for_overlap(radius_bbj, radius_tpmi, bbj_tpmi_total)

venn_polygons <- rbindlist(list(
  ellipse_polygon("PanUKB", center_panukb, 0, radius_panukb, 1L),
  ellipse_polygon("BBJ", center_bbj, 0, radius_bbj, 2L),
  ellipse_polygon("TPMI", center_tpmi, 0, radius_tpmi, 3L)
))
setorder(venn_polygons, draw_order, angle)
venn_polygons[, group := factor(group, levels = c("PanUKB", "BBJ", "TPMI"))]

count_positions <- data.table(
  region = c(
    "BBJ only",
    "PanUKB only",
    "BBJ + PanUKB only",
    "TPMI + PanUKB only",
    "BBJ + TPMI + PanUKB"
  ),
  x = c(-2.15, 3.95, 0.55, 3.10, 1.45),
  y = c(0, 0, 2.20, 0, -0.25) * ellipse_y_scale
)
count_positions <- merge(count_positions, region_counts, by = "region", sort = FALSE)

cohort_totals <- data.table(
  cohort = c("BBJ", "TPMI", "UKBB"),
  n_significant = c(bbj_total, tpmi_total, panukb_total),
  x = c(-1.75, center_tpmi, 3.75),
  y = c(3.12, 0.65, 3.12) * ellipse_y_scale
)
none_count <- region_counts[region == "None", n_pairs]

# Compress B horizontally about its centre without changing its height or text size.
venn_x_center <- mean(c(-2.90, 4.95))
for (coordinates in list(venn_polygons, count_positions, cohort_totals)) {
  coordinates[, x := venn_x_center + 0.80 * (x - venn_x_center)]
}

panel_b_plot <- ggplot() +
  geom_polygon(
    data = venn_polygons,
    aes(x = x, y = y, group = group, fill = group, colour = group),
    alpha = 0.32,
    linewidth = 0.95
  ) +
  geom_text(
    data = count_positions,
    aes(x = x, y = y, label = n_pairs),
    size = 25 / ggplot2::.pt,
    fontface = "bold",
    colour = "#202A33"
  ) +
  geom_text(
    data = cohort_totals,
    aes(
      x = x,
      y = y,
      label = paste0(cohort, "\n", n_significant, " significant")
    ),
    size = 23 / ggplot2::.pt,
    lineheight = 1.02,
    fontface = "bold",
    colour = "#202A33"
  ) +
  scale_fill_manual(values = venn_colours) +
  scale_colour_manual(values = venn_colours) +
  coord_equal(
    xlim = c(-2.90, 4.95),
    ylim = c(-3.20, 3.65) * ellipse_y_scale,
    expand = FALSE,
    clip = "off"
  ) +
  labs(
    title = NULL
  ) +
  theme_void(base_size = 15.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 25,
      colour = "#202A33",
      hjust = 0,
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(
      size = 14.4,
      colour = "#4B5965",
      hjust = 0.5,
      margin = margin(b = 4)
    ),
    plot.caption = element_text(
      size = 13.2,
      colour = "#5F6B73",
      hjust = 0.5,
      margin = margin(t = 3)
    ),
    plot.margin = margin(0, 8, 0, 8)
  )

# Reuse the standalone plots in isolated environments, without rewriting their files.
component_script_dir <- script_dir
count_env <- new.env(parent = globalenv())
count_env$render_only <- TRUE
sys.source(file.path(component_script_dir, "plot_BBJ_TPMI_method_significant_pairs.R"),
           envir = count_env)
panel_d_data <- copy(count_env$counts)
panel_d_plot <- count_env$figure + theme(
  text = element_text(size = 19),
  strip.text = element_text(size = 25, face = "bold", family = "Helvetica",
                            colour = "#202A33", vjust = 0.5,
                            margin = margin(t = 3, b = 14)),
  axis.title.y = element_text(size = 21, margin = margin(r = 12)),
  axis.text = element_text(size = 18, colour = "#252525"),
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = 7)),
  plot.margin = margin(16, 18, 18, 14)
)
stopifnot(inherits(panel_d_plot$layers[[2]]$geom, "GeomText"))
panel_d_plot$layers[[2]]$aes_params$size <- 5 + 5 / ggplot2::.pt

case_env <- new.env(parent = globalenv())
case_env$render_only <- TRUE
sys.source(file.path(component_script_dir, "plot_BBJ_DBP_MI_BMI_PAD_shared_columns.R"),
           envir = case_env)
panel_c_data <- copy(case_env$d)
panel_c_plot <- case_env$panel_plot +
  labs(title = "BBJ") +
  theme(
    plot.title = element_text(size = 25, face = "bold", family = "Helvetica",
                              colour = "#202A33", hjust = 0.5, vjust = 0.5,
                              margin = margin(t = 3, b = 14)),
    plot.title.position = "plot",
    plot.margin = margin(16, 7, 0, 7)
  )

stopifnot(
  nrow(panel_c_data) == 10L, nrow(panel_d_data) == 10L,
  !any(panel_c_data$method %in% c("TEMR", "MRBEE ancestry-specific")),
  identical(unname(count_env$method_colours[names(case_env$method_colours)]),
            unname(case_env$method_colours)),
  all(panel_d_data[method == "XMR", n_significant] ==
      c(bbj_total, tpmi_total))
)

panel_a_header <- ggdraw() +
  draw_label(
    "A",
    x = 0.0025,
    y = 1 - 0.275 / 0.95,
    hjust = 0,
    vjust = 0.5,
    fontface = "bold",
    size = 25,
    colour = "#202A33"
  ) +
  draw_label("XMR", x = 0.5, y = 1 - 0.75 / 0.95,
             fontface = "bold", size = 25, colour = "#202A33")

panel_b_header <- ggdraw() +
  draw_label(
    "B",
    x = 0.0025,
    y = 1 - 0.275 / 0.95,
    hjust = 0,
    vjust = 0.5,
    fontface = "bold",
    size = 25,
    colour = "#202A33"
  ) +
  draw_label("XMR", x = 0.5, y = 1 - 0.75 / 0.95,
             fontface = "bold", size = 25, colour = "#202A33")

panel_c_header <- ggdraw() +
  draw_label("C", x = 0.0025, y = 0.5, hjust = 0, vjust = 0.5,
             fontface = "bold", size = 25, colour = "#202A33")
panel_d_header <- ggdraw() +
  draw_label("D", x = 0.0025, y = 0.5, hjust = 0, vjust = 0.5,
             fontface = "bold", size = 25, colour = "#202A33")

top_widths <- column_widths
header_height <- 0.95
bottom_height <- 6.50 * 1.30

# Size the top row to its content, keeping A's square cells and B's dimensions.
grDevices::pdf(NULL, width = canvas_width, height = 16, family = "Helvetica")
panel_a_grob <- ggplotGrob(panel_a)
a_panel_row <- panel_a_grob$layout$t[panel_a_grob$layout$name == "panel"]
a_panel_column <- panel_a_grob$layout$l[panel_a_grob$layout$name == "panel"]
a_outer_width <- grid::convertWidth(
  sum(panel_a_grob$widths[-a_panel_column]), "in", valueOnly = TRUE
)
panel_height <- (canvas_width * top_widths[1] / sum(top_widths) - a_outer_width) *
  (diff(range(exposure_map$exposure_index)) + 1.10) /
  (diff(range(outcome_map$outcome_index)) + 1.10)
panel_a_grob$heights[a_panel_row] <- grid::unit(panel_height, "in")
panel_a_grob$respect <- FALSE
a_content_height <- grid::convertHeight(sum(panel_a_grob$heights), "in", valueOnly = TRUE)
grDevices::dev.off()
b_content_height <- (column_widths[2] - venn_margin_inches) *
  diff(c(-3.20, 3.65) * ellipse_y_scale) / diff(c(-2.90, 4.95))
body_height <- max(a_content_height, b_content_height + 0.30)
top_height <- header_height + body_height
canvas_height <- top_height + bottom_height
panel_a_block <- plot_grid(panel_a_header, panel_a_grob, ncol = 1,
                          rel_heights = c(header_height, body_height))
panel_b_body <- ggdraw() +
  draw_plot(panel_b_plot, x = 0, y = 0, width = 1, height = 1)
panel_b_block <- plot_grid(panel_b_header, panel_b_body, ncol = 1,
                          rel_heights = c(header_height, body_height))

panel_c <- plot_grid(
  panel_c_header, panel_c_plot,
  ncol = 1, rel_heights = c(0.08, 0.92)
)
panel_d <- plot_grid(
  panel_d_header, panel_d_plot, ncol = 1, rel_heights = c(0.08, 0.92)
)
bottom_row <- plot_grid(
  panel_c, panel_d, nrow = 1, rel_widths = top_widths
)

top_row <- plot_grid(
  panel_a_block,
  panel_b_block,
  nrow = 1,
  rel_widths = top_widths
)

figure <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(top_height, bottom_height)
)

if (!exists("render_only", inherits = FALSE) || !render_only) {
file_stem <- "BBJ_TPMI_UKBB_MR_comparison_ABCD"
pdf_path <- file.path(output_dir, paste0(file_stem, ".pdf"))
png_path <- file.path(output_dir, paste0(file_stem, ".png"))
panel_a_data_path <- file.path(
  output_dir,
  paste0(file_stem, "_panelA_data.tsv")
)
panel_b_data_path <- file.path(output_dir, paste0(file_stem, "_panelB_data.tsv"))
panel_c_data_path <- file.path(output_dir, paste0(file_stem, "_panelC_data.tsv"))
panel_d_data_path <- file.path(output_dir, paste0(file_stem, "_panelD_data.tsv"))

ggsave(
  pdf_path,
  figure,
  width = canvas_width,
  height = canvas_height,
  units = "in",
  device = grDevices::pdf,
  title = "MR discoveries and method comparisons across BBJ, TPMI, and UKBB",
  useDingbats = FALSE,
  encoding = "WinAnsi.enc",
  bg = "white",
  limitsize = FALSE
)
ggsave(
  png_path,
  figure,
  width = canvas_width,
  height = canvas_height,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

fwrite(
  panel_a_data[, .(
    cohort,
    exposure_name,
    outcome_name,
    exposure_label,
    outcome_label,
    beta,
    se,
    p,
    n_IV,
    pi0,
    nvalid,
    nvalid_source,
    n_tests,
    bonferroni_threshold,
    significant_bonferroni,
    direction,
    abs_beta,
    exposure_index,
    outcome_index,
    x_plot,
    y_plot
  )],
  panel_a_data_path,
  sep = "\t",
  na = "NA"
)
fwrite(region_counts, panel_b_data_path, sep = "\t")
fwrite(rbindlist(list(
  exposure_map[, .(role = "Exposure", abbreviation, full_name = exposure_name)],
  outcome_map[, .(role = "Outcome", abbreviation, full_name = outcome_name)]
)), file.path(output_dir, paste0(file_stem, "_panelA_abbreviations.tsv")), sep = "\t")
fwrite(panel_c_data, panel_c_data_path, sep = "\t", na = "NA")
fwrite(panel_d_data, panel_d_data_path, sep = "\t", na = "NA")
fwrite(significance_wide, file.path(output_dir, paste0(file_stem, "_panelB_all_data.tsv")),
       sep = "\t", na = "NA")
fwrite(count_env$all_data, file.path(output_dir, paste0(file_stem, "_panelD_all_data.tsv")),
       sep = "\t", na = "NA")
summary_table <- long_data[, .(
  method = unique(method),
  n_pairs = .N,
  n_valid = sum(valid_result),
  bonferroni_threshold = unique(bonferroni_threshold),
  n_significant = sum(significant_bonferroni),
  n_positive_significant = sum(significant_bonferroni & beta > 0),
  n_negative_significant = sum(significant_bonferroni & beta < 0)
), by = cohort]
print(summary_table)
print(region_counts)
stopifnot(n_tests == 105L, nrow(long_data) == 315L)
fwrite(long_data, file.path(output_dir, paste0(file_stem, "_all_data.tsv")), sep = "\t", na = "NA")
fwrite(region_counts, file.path(output_dir, paste0(file_stem, "_panelB_counts.tsv")), sep = "\t")
fwrite(panel_d_data, file.path(output_dir, paste0(file_stem, "_panelD_counts.tsv")), sep = "\t")
fwrite(summary_table, file.path(output_dir, paste0(file_stem, "_checks.tsv")), sep = "\t")
component_inputs <- unique(c(
  unname(input_files), count_env$input_paths, case_env$input_paths
))
fwrite(data.table(
  path = normalizePath(component_inputs),
  md5 = unname(tools::md5sum(component_inputs)),
  rows = vapply(component_inputs, function(p) nrow(fread(p)), integer(1L))
), file.path(output_dir, paste0(file_stem, "_inputs.tsv")), sep = "\t")
capture.output(sessionInfo(), file = file.path(output_dir, paste0(file_stem, "_sessionInfo.txt")))

print(panel_c_data[, .(pair, method, threshold, n_IV, beta, se, p)])
print(panel_d_data[, .(cohort, method, threshold, n_significant, n_included)])
message("Saved PDF: ", normalizePath(pdf_path))
message("Saved PNG: ", normalizePath(png_path))
message("Saved panel A data: ", normalizePath(panel_a_data_path))
message("Saved panel B overlap data: ", normalizePath(panel_b_data_path))
message("Saved panel C shared-column case data: ", normalizePath(panel_c_data_path))
message("Saved panel D count data: ", normalizePath(panel_d_data_path))
message(sprintf("Canvas: %.2f x %.2f inches; A matrix: %.2f inches high; second row: %.2f inches",
                canvas_width, canvas_height, panel_height, bottom_height))
message("Cross-cohort MR comparison; BBJ/UKBB labels aligned; TPMI label centred; panel layout unchanged")
}
