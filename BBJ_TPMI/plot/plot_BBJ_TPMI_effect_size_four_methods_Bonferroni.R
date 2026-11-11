#!/usr/bin/env Rscript

# Compare BBJ and TPMI estimates for the same 105 exposure--outcome pairs.
# XMR: 5e-05, crossPopIntercept0; CAMERA / MRBEE target / MRBEE-TL: 5e-08.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else getwd()
project_dir <- normalizePath(file.path(script_dir, ".."))
results_dir <- file.path(project_dir, "results", "masterClump5e-05")
output_dir <- file.path(project_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(output_dir, "BBJ_TPMI_effect_size_four_methods_Bonferroni")

input_paths <- file.path(results_dir, c(
  "BBJ/crossPopIntercept0/XMR_EAS_BBJ_EURselected_two_thresholds_noEffectHetQC.tsv",
  "TPMI/crossPopIntercept0/XMR_EAS_TPMI_EURselected_two_thresholds_noEffectHetQC.tsv",
  "CAMERA_EUR_to_BBJ_TPMI.tsv",
  "MRBEE_TL_EUR_to_BBJ_TPMI.tsv"
))
stopifnot(all(file.exists(input_paths)))
input_hashes <- tools::md5sum(input_paths)
method_names <- c(XMR = "XMR", CAMERA = "CAMERA",
                  MRBEE_target_specific = "MRBEE", MRBEE_TL_UV = "MRBEE-TL")
method_levels <- unname(method_names)
iv_thresholds <- c(XMR = 5e-5, CAMERA = 5e-8, MRBEE = 5e-8, "MRBEE-TL" = 5e-8)
n_tests <- 105L
bonferroni_cutoff <- 0.05 / n_tests

# Exact significance-category colours from the supplied four_methods_FDR.pdf.
# Classifications still use Bonferroni correction, not the reference PDF's FDR.
category_levels <- c("Both significant", "BBJ only", "TPMI only", "Neither")
category_colours <- c("Both significant" = "#FF7F50", "BBJ only" = "#009ACD",
                      "TPMI only" = "#66CD00", "Neither" = "#000000")

all_data <- rbindlist(lapply(input_paths, function(path) {
  d <- fread(path, colClasses = c(exposure = "character", outcome = "character"))
  d[, source_row := seq_len(.N)]
  stopifnot(all(c("method", "target", "exposure", "outcome", "threshold",
                  "beta", "se", "p", "status") %in% names(d)))
  if (all(d$method == "XMR")) {
    stopifnot(all(d$background_cross_population_intercept == 0),
              all(d$master_clump_threshold == 5e-5))
    d[, `:=`(CI_low = beta - 1.96 * se, CI_high = beta + 1.96 * se)]
  }
  # Source-specific MRBEE rows are EUR estimates and must not enter this comparison.
  d <- d[method %chin% names(method_names)]
  d[, method_label := unname(method_names[method])]
  d <- d[threshold == unname(iv_thresholds[method_label])]
  d[, .(method = method_label, cohort = fifelse(target == "TW", "TPMI", target),
         exposure = as.character(as.numeric(exposure)),
         outcome = as.character(as.numeric(outcome)),
         threshold, beta, se, p, CI_low, CI_high, status, error,
         source_file = normalizePath(path), source_row)]
}))
keys <- c("method", "cohort", "exposure", "outcome")
stopifnot(!anyDuplicated(all_data[, ..keys]),
          all(all_data$cohort %chin% c("BBJ", "TPMI")))
reference_pairs <- all_data[method == "XMR" & cohort == "BBJ", .(exposure, outcome)]
stopifnot(nrow(reference_pairs) == n_tests,
          uniqueN(reference_pairs$exposure) == 7L,
          uniqueN(reference_pairs$outcome) == 15L,
          fsetequal(reference_pairs, CJ(exposure = unique(reference_pairs$exposure),
                                      outcome = unique(reference_pairs$outcome))))
for (m in method_levels) for (pop in c("BBJ", "TPMI")) {
  stopifnot(fsetequal(reference_pairs,
    all_data[method == m & cohort == pop, .(exposure, outcome)]))
}
all_data[, valid := !is.na(status) & status == "ok" & is.finite(beta) &
           is.finite(se) & se > 0 & is.finite(p) & p >= 0 & p <= 1 &
           is.finite(CI_low) & is.finite(CI_high) & CI_low <= beta & beta <= CI_high]
# All supplied fits at the requested thresholds are valid. Never turn missing fits into P=1.
stopifnot(nrow(all_data) == 8L * n_tests, all(all_data$valid))
all_data[, `:=`(significant = p < bonferroni_cutoff,
                bonferroni_threshold = bonferroni_cutoff)]

match_keys <- c("method", "exposure", "outcome", "threshold")
values <- c("beta", "se", "p", "CI_low", "CI_high", "significant", "source_file", "source_row")
bbj <- all_data[cohort == "BBJ", c(match_keys, values), with = FALSE]
tpmi <- all_data[cohort == "TPMI", c(match_keys, values), with = FALSE]
setnames(bbj, values, paste0(values, "_BBJ"))
setnames(tpmi, values, paste0(values, "_TPMI"))
pairs <- merge(bbj, tpmi, by = match_keys, all = TRUE)
stopifnot(nrow(pairs) == 4L * n_tests, !anyNA(pairs))
pairs[, category := factor(fcase(
  significant_BBJ & significant_TPMI, "Both significant",
  significant_BBJ, "BBJ only", significant_TPMI, "TPMI only",
  default = "Neither"
), levels = category_levels)]
pairs[, bonferroni_threshold := bonferroni_cutoff]
# Trait codes and names follow the BBJ/TPMI XMR input files.
exposure_abbreviations <- c("21001" = "BMI", "30020" = "Hb", "30620" = "ALT",
                            "30740" = "Glu", "30750" = "HbA1c",
                            "30780" = "LDL-C", "4079" = "DBP")
outcome_abbreviations <- c("174.1" = "BC", "185" = "PrC", "218.1" = "UF",
                           "250.2" = "T2D", "280.1" = "IDA", "345.1" = "Ep",
                           "366" = "Cat", "411.1" = "UAP", "411.2" = "MI",
                           "411.3" = "Ang", "443.9" = "PAD", "495" = "As",
                           "496.21" = "COPD", "574.1" = "CL", "714.1" = "RA")
pairs[, `:=`(exposure_abbreviation = unname(exposure_abbreviations[exposure]),
              outcome_abbreviation = unname(outcome_abbreviations[outcome]))]
stopifnot(!anyNA(pairs$exposure_abbreviation), !anyNA(pairs$outcome_abbreviation))
pairs[, pair_label := paste(exposure_abbreviation, outcome_abbreviation, sep = "-")]
pairs[, plot_label := fifelse(category == "Both significant", pair_label, "")]
summary_table <- pairs[, .(
  n_pairs = .N, pearson_r = cor(beta_BBJ, beta_TPMI),
  both = sum(significant_BBJ & significant_TPMI),
  bbj_only = sum(significant_BBJ & !significant_TPMI),
  tpmi_only = sum(!significant_BBJ & significant_TPMI),
  neither = sum(!significant_BBJ & !significant_TPMI),
  significant_BBJ = sum(significant_BBJ), significant_TPMI = sum(significant_TPMI),
  both_significant_opposite_direction = sum(significant_BBJ & significant_TPMI &
                                            beta_BBJ * beta_TPMI < 0)
), by = .(method, threshold)]
summary_table <- summary_table[match(method_levels, method)]
stopifnot(all(summary_table$n_pairs == n_tests),
          all(summary_table$both + summary_table$bbj_only +
              summary_table$tpmi_only + summary_table$neither == n_tests))

# Descriptive OLS fit across all 105 pairs, separately for each method.
# The ribbon is the pointwise 95% confidence interval for the fitted mean,
# not an individual-pair prediction interval or a measurement-error regression.
linear_models <- setNames(lapply(method_levels, function(m) {
  lm(beta_TPMI ~ beta_BBJ, data = pairs[method == m])
}), method_levels)
fit_data <- rbindlist(lapply(method_levels, function(m) {
  d <- pairs[method == m]
  grid <- data.frame(beta_BBJ = seq(min(d$beta_BBJ), max(d$beta_BBJ), length.out = 200))
  prediction <- predict(linear_models[[m]], newdata = grid,
                        interval = "confidence", level = 0.95)
  data.table(method = m, beta_BBJ = grid$beta_BBJ,
             fitted_TPMI = prediction[, "fit"],
             lower = prediction[, "lwr"], upper = prediction[, "upr"])
}))
fit_summary <- rbindlist(lapply(method_levels, function(m) {
  fit <- linear_models[[m]]
  interval <- confint(fit, level = 0.95)
  data.table(method = m, n_pairs = nobs(fit),
             intercept = unname(coef(fit)[1L]), slope = unname(coef(fit)[2L]),
             slope_CI_low = interval[2L, 1L], slope_CI_high = interval[2L, 2L],
             residual_df = df.residual(fit), residual_sd = sigma(fit))
}))
stopifnot(all(fit_summary$n_pairs == n_tests),
          all(fit_data$lower <= fit_data$fitted_TPMI),
          all(fit_data$upper >= fit_data$fitted_TPMI))

# Shared limits include every confidence interval, including the wide MRBEE-TL intervals.
x_limits <- c(-0.40, 1.65)
y_limits <- c(-0.25, 0.73)
stopifnot(all(pairs$CI_low_BBJ > x_limits[1L]), all(pairs$CI_high_BBJ < x_limits[2L]),
          all(pairs$CI_low_TPMI > y_limits[1L]), all(pairs$CI_high_TPMI < y_limits[2L]))

make_panel <- function(m, j) {
  d <- copy(pairs[method == m])
  s <- summary_table[method == m]
  d[, draw_order := fifelse(category == "Neither", 1L,
                            fifelse(category == "Both significant", 3L, 2L))]
  setorder(d, draw_order)
  ggplot(d, aes(beta_BBJ, beta_TPMI, colour = category)) +
    geom_hline(yintercept = 0, colour = "#B5B5B5", linetype = "dotted", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "#B5B5B5", linetype = "dotted", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, colour = "#858585", linetype = "dashed", linewidth = 0.5) +
    geom_ribbon(data = fit_data[method == m],
                 aes(x = beta_BBJ, ymin = lower, ymax = upper),
                 inherit.aes = FALSE, fill = "#FFC0CB", alpha = 0.4,
                 colour = NA, show.legend = FALSE) +
    geom_segment(aes(x = CI_low_BBJ, xend = CI_high_BBJ, yend = beta_TPMI),
                  linewidth = 0.32, alpha = 0.4, show.legend = FALSE) +
    geom_segment(aes(y = CI_low_TPMI, yend = CI_high_TPMI, xend = beta_BBJ),
                  linewidth = 0.32, alpha = 0.4, show.legend = FALSE) +
    geom_point(data = d[category == "Neither"], size = 1.9, alpha = 0.75) +
    geom_point(data = d[category != "Neither"], size = 2.8, alpha = 0.9) +
    geom_line(data = fit_data[method == m], aes(x = beta_BBJ, y = fitted_TPMI),
               inherit.aes = FALSE, colour = "#FF0000", linewidth = 0.9,
               show.legend = FALSE) +
    geom_text_repel(aes(label = plot_label), colour = "#333333", family = "Helvetica",
                     size = 14 / .pt, box.padding = 0.4, point.padding = 0.25,
                     point.size = 2.8, min.segment.length = 0,
                     segment.colour = "#777777", segment.size = 0.35,
                     force = 2, force_pull = 0.3, nudge_y = 0.025,
                     seed = 20261027, max.overlaps = Inf, max.iter = 60000,
                     max.time = Inf, show.legend = FALSE) +
    annotate("text", x = -0.34, y = 0.68, hjust = 0, vjust = 1,
              label = sprintf("Pearson r = %.3f", s$pearson_r),
              size = 16 / .pt, colour = "#222222", family = "Helvetica") +
    scale_colour_manual(values = category_colours, breaks = category_levels,
                         drop = FALSE, name = "Significance") +
    scale_x_continuous(breaks = seq(0, 1.5, 0.5), expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = seq(-0.2, 0.6, 0.2), expand = expansion(mult = 0)) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    labs(title = m,
         x = "Effect size in BBJ", y = if (j == 1L) "Effect size in TPMI" else NULL,
         tag = LETTERS[j]) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3.5, alpha = 1))) +
    theme_classic(base_size = 16, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 23, face = "bold", hjust = 0.5,
                                margin = margin(b = 14)),
      plot.tag = element_text(size = 23, face = "bold"),
      axis.title = element_text(size = 18),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text = element_text(size = 14, colour = "#333333"),
      axis.line = element_line(linewidth = 0.4, colour = "#555555"),
      axis.ticks = element_line(linewidth = 0.4, colour = "#555555"),
      legend.position = "bottom", legend.title = element_text(size = 21),
      legend.text = element_text(size = 21), legend.key.width = grid::unit(0.7, "cm"),
      legend.box.spacing = grid::unit(0, "pt"),
      plot.margin = margin(10, 14, 8, 10)
    )
}

panels <- lapply(seq_along(method_levels), function(j) make_panel(method_levels[j], j))
legend <- get_legend(panels[[1L]])
panel_row <- plot_grid(plotlist = lapply(panels, function(p) p + theme(legend.position = "none")),
                      nrow = 1, align = "hv", axis = "tblr")
figure <- plot_grid(panel_row, legend, ncol = 1, rel_heights = c(1, 0.10))
if (!exists("render_only", inherits = FALSE) || !render_only) {
ggsave(paste0(stem, ".pdf"), figure, width = 24, height = 6.8,
        device = grDevices::pdf, useDingbats = FALSE, bg = "white",
        title = "BBJ and TPMI MR effect estimates with Bonferroni significance")
ggsave(paste0(stem, ".png"), figure, width = 24, height = 6.8, dpi = 220, bg = "white")
fwrite(all_data, paste0(stem, "_all_data.tsv"), sep = "\t", na = "NA")
fwrite(pairs, paste0(stem, "_plot_data.tsv"), sep = "\t", na = "NA")
fwrite(pairs[category == "Both significant"], paste0(stem, "_labelled_pairs.tsv"),
        sep = "\t", na = "NA")
fwrite(summary_table, paste0(stem, "_summary.tsv"), sep = "\t")
fwrite(fit_data, paste0(stem, "_linear_fit.tsv"), sep = "\t")
fwrite(fit_summary, paste0(stem, "_linear_fit_summary.tsv"), sep = "\t")
fwrite(data.table(source_file = normalizePath(input_paths), md5 = unname(input_hashes)),
        paste0(stem, "_inputs.tsv"), sep = "\t")
capture.output(sessionInfo(), file = paste0(stem, "_sessionInfo.txt"))
stopifnot(identical(input_hashes, tools::md5sum(input_paths)))
pdf_file <- paste0(stem, ".pdf")
stopifnot(identical(input_hashes, tools::md5sum(input_paths)))
print(summary_table)
message("Saved: ", pdf_file)
}
