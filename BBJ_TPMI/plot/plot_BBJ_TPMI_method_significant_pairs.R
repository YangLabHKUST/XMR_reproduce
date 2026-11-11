#!/usr/bin/env Rscript

# One combined count-only figure; correction always uses 105 planned pairs.
# Each method uses its own available fits; MR-APSS failures do not exclude other methods.
# XMR / MR-APSS: 5e-05; CAMERA / MRBEE / MRBEE-TL: 5e-08.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else getwd()
project_dir <- normalizePath(file.path(script_dir, ".."))
results_dir <- file.path(project_dir, "results")
output_dir <- file.path(project_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cohorts <- c("BBJ", "TPMI")
n_planned <- 105L
significance_cutoff <- 0.05 / n_planned

input_paths <- file.path(results_dir, c(
  "masterClump5e-05/BBJ/crossPopIntercept0/XMR_EAS_BBJ_EURselected_two_thresholds_noEffectHetQC.tsv",
  "masterClump5e-05/TPMI/crossPopIntercept0/XMR_EAS_TPMI_EURselected_two_thresholds_noEffectHetQC.tsv",
  "sharedIV_masterClump5e-05/MRAPSS/BBJ/MRAPSS_BBJ_results.tsv",
  "sharedIV_masterClump5e-05/MRAPSS/TPMI/MRAPSS_TPMI_results.tsv",
  "masterClump5e-05/MRBEE_TL_EUR_to_BBJ_TPMI.tsv",
  "masterClump5e-05/CAMERA_EUR_to_BBJ_TPMI.tsv"
))
stopifnot(all(file.exists(input_paths)))
method_names <- c(
  XMR = "XMR", MRAPSS = "MR-APSS",
  CAMERA = "CAMERA",
  MRBEE_target_specific = "MRBEE", MRBEE_TL_UV = "MRBEE-TL"
)
method_levels <- unname(method_names)
method_colours <- c(
  "XMR" = "#D73027", "MR-APSS" = "#332F2C",
  "CAMERA" = "#0072B2",
  "MRBEE" = "#009E73", "MRBEE-TL" = "#E69F00"
)
method_labels <- c(
  "XMR" = "XMR", "MR-APSS" = "MR-APSS",
  "CAMERA" = "CAMERA",
  "MRBEE" = "MRBEE", "MRBEE-TL" = "MRBEE-TL"
)
selected_thresholds <- c("XMR" = 5e-5, "MR-APSS" = 5e-5,
                         "CAMERA" = 5e-8,
                         "MRBEE" = 5e-8, "MRBEE-TL" = 5e-8)

# Preserve the reported P values and missing status; do not replace failed fits by P=1.
all_data <- rbindlist(lapply(input_paths, function(path) {
  d <- fread(path, colClasses = c(exposure = "character", outcome = "character"))
  required <- c("method", "target", "exposure", "outcome", "threshold",
                "beta", "se", "p", "status", "error")
  stopifnot(all(required %in% names(d)))
  d[, source_row := seq_len(.N)]
  if (all(d$method == "XMR")) {
    stopifnot(all(d$background_cross_population_intercept == 0),
              all(d$master_clump_threshold == 5e-5))
  }
  d <- d[method != "MRBEE_source_specific"]
  stopifnot(all(d$method %in% names(method_names)))
  d[, method_label := unname(method_names[method])]
  d <- d[threshold == unname(selected_thresholds[method_label])]
  d[, .(
    cohort = fifelse(target == "TW", "TPMI", target),
    method = method_label,
    exposure = as.character(as.numeric(exposure)),
    outcome = as.character(as.numeric(outcome)),
    threshold, beta, se, p, status = as.character(status), error = as.character(error),
    source_file = normalizePath(path), source_row
  )]
}))
stopifnot(all(all_data$cohort %in% cohorts),
          all(all_data$status %in% c("ok", "failed")))
keys <- c("cohort", "method", "exposure", "outcome")
stopifnot(!anyDuplicated(all_data[, ..keys]))
reference_pairs <- all_data[cohort == "BBJ" & method == "XMR", .(exposure, outcome)]
stopifnot(nrow(reference_pairs) == n_planned,
          uniqueN(reference_pairs$exposure) == 7L,
          uniqueN(reference_pairs$outcome) == 15L,
          fsetequal(reference_pairs, CJ(exposure = unique(reference_pairs$exposure),
                                      outcome = unique(reference_pairs$outcome))))
for (pop in cohorts) for (m in method_levels) {
  stopifnot(fsetequal(reference_pairs,
                     all_data[cohort == pop & method == m, .(exposure, outcome)]))
}
all_data[, valid := status == "ok" & is.finite(beta) & is.finite(se) & se > 0 &
           is.finite(p) & p >= 0 & p <= 1]
stopifnot(!anyNA(all_data$valid))
all_data[, significant := fifelse(valid, p < significance_cutoff, NA)]
all_data[, bonferroni_cutoff := significance_cutoff]

# Keep each method's own valid fits. MR-APSS availability is retained for audit only.
mrapss_availability <- all_data[method == "MR-APSS", .(
  cohort, exposure, outcome, mrapss_available = valid
)]
all_data[mrapss_availability, on = .(cohort, exposure, outcome),
         mrapss_available := i.mrapss_available]
stopifnot(!anyNA(all_data$mrapss_available))
all_data[, included := valid]
all_data[, exclusion_reason := fifelse(
  included, NA_character_, "Method-specific fit unavailable"
)]
stopifnot(all(all_data[method != "MR-APSS", included]),
          all_data[cohort == "BBJ" & method == "MR-APSS", sum(included)] == 105L,
          all_data[cohort == "TPMI" & method == "MR-APSS", sum(included)] == 96L)
counts <- all_data[, .(
  n_planned = .N, n_included = sum(included), n_excluded = sum(!included),
  n_valid = sum(included & valid), n_unavailable = sum(!valid),
  n_significant = sum(significant[included], na.rm = TRUE),
  bonferroni_cutoff = significance_cutoff
), by = .(cohort, method, threshold)]
counts[, method := factor(method, levels = method_levels)]
counts[, cohort := factor(cohort, levels = cohorts)]
setorder(counts, cohort, method)
stopifnot(nrow(counts) == 10L, all(counts$n_planned == 105L))

# BBJ is the left panel and TPMI the right; both use the same y scale.
y_break_max <- max(5, 5 * ceiling(max(counts$n_significant) / 5))
figure <- ggplot(counts, aes(method, n_significant, fill = method)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = n_significant), vjust = -0.45,
            size = 5, colour = "#252525", family = "Helvetica") +
  facet_grid(cols = vars(cohort)) +
  scale_fill_manual(values = method_colours, guide = "none") +
  scale_x_discrete(labels = method_labels, expand = expansion(add = 0.6)) +
  scale_y_continuous(limits = c(0, y_break_max + 2.5),
                     breaks = seq(0, y_break_max, 5), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Number of significant pairs") +
  theme_classic(base_size = 14, base_family = "Helvetica") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 19, face = "bold",
                              margin = margin(t = 3, b = 14)),
    panel.spacing.x = grid::unit(0.9, "cm"),
    panel.grid.major.y = element_line(colour = "#E7E7E7", linewidth = 0.3),
    axis.text = element_text(size = 13, colour = "#252525"),
    axis.text.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(size = 14, margin = margin(r = 12)),
    axis.line = element_line(colour = "#666666", linewidth = 0.35),
    axis.ticks = element_line(colour = "#666666", linewidth = 0.35),
    plot.margin = margin(12, 18, 12, 14)
  )
if (!exists("render_only", inherits = FALSE) || !render_only) {
filename <- file.path(output_dir, "BBJ_TPMI_method_significant_pairs_Bonferroni")
ggsave(paste0(filename, ".pdf"), figure, width = 12.5, height = 4.9,
       device = grDevices::pdf, useDingbats = FALSE, bg = "white")
ggsave(paste0(filename, ".png"), figure, width = 12.5, height = 4.9,
       dpi = 220, bg = "white")
stopifnot(all(file.exists(paste0(filename, c(".pdf", ".png")))),
          all(file.info(paste0(filename, c(".pdf", ".png")))$size > 0))

stem <- "BBJ_TPMI_method_significant_pairs"
export_table <- function(d, suffix) {
  fwrite(d, file.path(output_dir, paste0(stem, "_", suffix, ".tsv")), sep = "\t", na = "NA")
}
export_table(counts, "counts")
export_table(all_data, "all_data")
export_table(all_data[included == TRUE & significant == TRUE], "significant_pairs")
export_table(all_data[valid == FALSE], "unavailable_fits")
export_table(all_data[included == FALSE], "excluded_right_group_pairs")
export_table(data.table(
  source_file = normalizePath(input_paths),
  md5 = unname(tools::md5sum(input_paths)), bytes = file.info(input_paths)$size
), "input_manifest")
capture.output(sessionInfo(), file = file.path(output_dir, paste0(stem, "_sessionInfo.txt")))
print(counts[, .(cohort, method, threshold, n_planned, n_included,
                 n_significant, n_unavailable, n_excluded)])
message("Saved combined BBJ and TPMI count-only figure: ", filename, ".pdf")
}
