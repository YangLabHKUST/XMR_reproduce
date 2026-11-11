#!/usr/bin/env Rscript

# Compare BBJ XMR with EUR MR-APSS at IV P < 5e-05 for 105 planned pairs.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else getwd()
project_dir <- normalizePath(file.path(script_dir, ".."))
output_dir <- file.path(project_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(output_dir, "BBJ_XMR_EUR_MRAPSS_effect_size_Bonferroni")
input_paths <- file.path(project_dir, "results", c(
  "masterClump5e-05/BBJ/crossPopIntercept0/XMR_EAS_BBJ_EURselected_two_thresholds_noEffectHetQC.tsv",
  "EUR/MRAPSS_EUR_exposure_PanUKB_15outcomes_5e-05.tsv"
))
stopifnot(all(file.exists(input_paths)))
input_hashes <- tools::md5sum(input_paths)
n_tests <- 105L
bonferroni_cutoff <- 0.05 / n_tests

bbj <- fread(input_paths[1L], colClasses = c(exposure = "character", outcome = "character"))
bbj[, source_row := seq_len(.N)]
bbj <- bbj[threshold == 5e-5]
stopifnot(all(bbj$method == "XMR"), all(bbj$target == "BBJ"),
          all(bbj$background_cross_population_intercept == 0))
bbj <- bbj[, .(
  exposure, outcome, exposure_name_BBJ = exposure_name, outcome_name_BBJ = outcome_name,
  beta_BBJ = beta, se_BBJ = se, p_BBJ = p, status_BBJ = status,
  threshold_BBJ = threshold, source_file_BBJ = normalizePath(input_paths[1L]),
  source_row_BBJ = source_row
)]
eur <- fread(input_paths[2L], colClasses = c(exposure = "character", outcome = "character",
                                            panukb_code = "character"))
eur[, source_row := seq_len(.N)]
eur <- eur[threshold == 5e-5]
stopifnot(all(eur$method == "MRAPSS"), all(eur$outcome_dataset == "PanUKB_EUR"))
# EUR stores abbreviations in outcome; match its panukb_code to the BBJ outcome code.
eur <- eur[, .(
  exposure, outcome = panukb_code, exposure_name_EUR = exposure_name,
  outcome_name_EUR = outcome_name, outcome_abbreviation = outcome,
  beta_EUR = beta, se_EUR = beta_se, p_EUR = beta_pvalue, status_EUR = status,
  threshold_EUR = threshold, source_file_EUR = normalizePath(input_paths[2L]),
  source_row_EUR = source_row
)]
for (d in list(bbj, eur)) {
  stopifnot(nrow(d) == n_tests, !anyDuplicated(d, by = c("exposure", "outcome")))
}
stopifnot(fsetequal(bbj[, .(exposure, outcome)], eur[, .(exposure, outcome)]))
pairs <- merge(bbj, eur, by = c("exposure", "outcome"), all = TRUE)
stopifnot(nrow(pairs) == n_tests, !anyNA(pairs),
          uniqueN(pairs$exposure) == 7L, uniqueN(pairs$outcome) == 15L,
          all(pairs$exposure_name_BBJ == pairs$exposure_name_EUR),
          all(pairs$outcome_name_BBJ == pairs$outcome_name_EUR),
          all(pairs$status_BBJ == "ok"), all(pairs$status_EUR == "ok"),
          all(is.finite(pairs$beta_BBJ)), all(is.finite(pairs$beta_EUR)),
          all(is.finite(pairs$se_BBJ) & pairs$se_BBJ > 0),
          all(is.finite(pairs$se_EUR) & pairs$se_EUR > 0),
          all(is.finite(pairs$p_BBJ) & pairs$p_BBJ >= 0 & pairs$p_BBJ <= 1),
          all(is.finite(pairs$p_EUR) & pairs$p_EUR >= 0 & pairs$p_EUR <= 1))
pairs[, `:=`(
  CI_low_BBJ = beta_BBJ - 1.96 * se_BBJ, CI_high_BBJ = beta_BBJ + 1.96 * se_BBJ,
  CI_low_EUR = beta_EUR - 1.96 * se_EUR, CI_high_EUR = beta_EUR + 1.96 * se_EUR,
  significant_BBJ = p_BBJ < bonferroni_cutoff, significant_EUR = p_EUR < bonferroni_cutoff,
  bonferroni_threshold = bonferroni_cutoff
)]
category_levels <- c("Both significant", "BBJ only", "EUR only", "Neither")
category_colours <- c("Both significant" = "#FF7F50", "BBJ only" = "#009ACD",
                      "EUR only" = "#66CD00", "Neither" = "#000000")
pairs[, category := factor(fcase(
  significant_BBJ & significant_EUR, "Both significant",
  significant_BBJ, "BBJ only", significant_EUR, "EUR only", default = "Neither"
), levels = category_levels)]
exposure_abbreviations <- c("21001" = "BMI", "30020" = "Hb", "30620" = "ALT",
                            "30740" = "Glu", "30750" = "HbA1c",
                            "30780" = "LDL-C", "4079" = "DBP")
pairs[, exposure_abbreviation := unname(exposure_abbreviations[exposure])]
stopifnot(!anyNA(pairs$exposure_abbreviation))
pairs[, pair_label := paste(exposure_abbreviation, outcome_abbreviation, sep = "-")]
pairs[, plot_label := fifelse(category == "Both significant", pair_label, "")]
summary_table <- pairs[, .(
  n_pairs = .N, pearson_r = cor(beta_BBJ, beta_EUR),
  significant_BBJ = sum(significant_BBJ), significant_EUR = sum(significant_EUR),
  both = sum(significant_BBJ & significant_EUR),
  bbj_only = sum(significant_BBJ & !significant_EUR),
  eur_only = sum(!significant_BBJ & significant_EUR),
  neither = sum(!significant_BBJ & !significant_EUR),
  both_significant_opposite_direction = sum(significant_BBJ & significant_EUR &
                                            beta_BBJ * beta_EUR < 0)
)]

# Descriptive regression uses every matched pair, not only significant pairs.
fit <- lm(beta_EUR ~ beta_BBJ, data = pairs)
fit_grid <- data.frame(beta_BBJ = seq(min(pairs$beta_BBJ), max(pairs$beta_BBJ), length.out = 200))
fit_prediction <- predict(fit, newdata = fit_grid, interval = "confidence", level = 0.95)
fit_data <- data.table(beta_BBJ = fit_grid$beta_BBJ,
                       fitted_EUR = fit_prediction[, "fit"],
                       lower = fit_prediction[, "lwr"], upper = fit_prediction[, "upr"])
fit_interval <- confint(fit, level = 0.95)
fit_summary <- data.table(n_pairs = nobs(fit), intercept = unname(coef(fit)[1L]),
                          slope = unname(coef(fit)[2L]), slope_CI_low = fit_interval[2L, 1L],
                          slope_CI_high = fit_interval[2L, 2L],
                          residual_df = df.residual(fit), residual_sd = sigma(fit))
x_limits <- c(-0.30, 1.50)
y_limits <- c(-0.20, 0.77)
stopifnot(all(pairs$CI_low_BBJ > x_limits[1L]), all(pairs$CI_high_BBJ < x_limits[2L]),
          all(pairs$CI_low_EUR > y_limits[1L]), all(pairs$CI_high_EUR < y_limits[2L]),
          all(fit_data$lower > y_limits[1L]), all(fit_data$upper < y_limits[2L]))
pairs[, draw_order := fifelse(category == "Neither", 1L,
                              fifelse(category == "Both significant", 3L, 2L))]
setorder(pairs, draw_order)
figure <- ggplot(pairs, aes(beta_BBJ, beta_EUR, colour = category)) +
  geom_hline(yintercept = 0, colour = "#B5B5B5", linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "#B5B5B5", linetype = "dotted", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, colour = "#858585", linetype = "dashed", linewidth = 0.5) +
  geom_ribbon(data = fit_data, aes(x = beta_BBJ, ymin = lower, ymax = upper),
                inherit.aes = FALSE, fill = "#FFC0CB", alpha = 0.4, colour = NA) +
  geom_segment(aes(x = CI_low_BBJ, xend = CI_high_BBJ, yend = beta_EUR),
                 linewidth = 0.35, alpha = 0.4, show.legend = FALSE) +
  geom_segment(aes(y = CI_low_EUR, yend = CI_high_EUR, xend = beta_BBJ),
                 linewidth = 0.35, alpha = 0.4, show.legend = FALSE) +
  geom_point(data = pairs[category == "Neither"], size = 2, alpha = 0.75) +
  geom_point(data = pairs[category != "Neither"], size = 3, alpha = 0.9) +
  geom_line(data = fit_data, aes(x = beta_BBJ, y = fitted_EUR), inherit.aes = FALSE,
              colour = "#FF0000", linewidth = 0.9, show.legend = FALSE) +
  geom_text_repel(aes(label = plot_label), colour = "#333333", family = "Helvetica",
                   size = 14 / .pt, box.padding = 0.45, point.padding = 0.3,
                   point.size = 3, min.segment.length = 0, segment.colour = "#777777",
                   segment.size = 0.35, force = 3, force_pull = 0.15,
                   seed = 20261027, max.overlaps = Inf, max.iter = 100000,
                   max.time = Inf, show.legend = FALSE) +
  annotate("text", x = -0.25, y = 0.72, hjust = 0, vjust = 1,
              label = sprintf("Pearson r = %.3f", summary_table$pearson_r),
              size = 20 / .pt, colour = "#222222", family = "Helvetica") +
  scale_colour_manual(values = category_colours, breaks = category_levels,
                       drop = FALSE, name = "Significance") +
  scale_x_continuous(breaks = seq(0, 1.5, 0.25), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = seq(-0.2, 0.6, 0.2), expand = expansion(mult = 0)) +
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  labs(title = "BBJ XMR versus EUR MR-APSS",
        x = "Effect size in BBJ (XMR)", y = "Effect size in EUR (MR-APSS)") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3.5, alpha = 1))) +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5,
                                  margin = margin(b = 16)),
         axis.title = element_text(size = 21), axis.title.x = element_text(margin = margin(t = 12)),
         axis.title.y = element_text(margin = margin(r = 12)),
         axis.text = element_text(size = 15, colour = "#333333"),
         axis.line = element_line(linewidth = 0.4, colour = "#555555"),
         axis.ticks = element_line(linewidth = 0.4, colour = "#555555"),
         legend.position = "bottom", legend.title = element_text(size = 21),
         legend.text = element_text(size = 21), legend.key.width = grid::unit(0.7, "cm"),
         legend.box.margin = margin(t = 12), plot.margin = margin(16, 18, 14, 12))
if (!exists("render_only", inherits = FALSE) || !render_only) {
ggsave(paste0(stem, ".pdf"), figure, width = 12, height = 9.8,
        device = grDevices::pdf, useDingbats = FALSE, bg = "white",
        title = "BBJ XMR and EUR MR-APSS effect estimates with Bonferroni significance")
ggsave(paste0(stem, ".png"), figure, width = 12, height = 9.8, dpi = 220, bg = "white")
fwrite(pairs, paste0(stem, "_plot_data.tsv"), sep = "\t")
fwrite(pairs[category == "Both significant"], paste0(stem, "_labelled_pairs.tsv"), sep = "\t")
fwrite(summary_table, paste0(stem, "_summary.tsv"), sep = "\t")
fwrite(fit_data, paste0(stem, "_linear_fit.tsv"), sep = "\t")
fwrite(fit_summary, paste0(stem, "_linear_fit_summary.tsv"), sep = "\t")
fwrite(data.table(source_file = normalizePath(input_paths), md5 = unname(input_hashes)),
        paste0(stem, "_inputs.tsv"), sep = "\t")
capture.output(sessionInfo(), file = paste0(stem, "_sessionInfo.txt"))
pdf_file <- paste0(stem, ".pdf")
stopifnot(identical(input_hashes, tools::md5sum(input_paths)))
print(summary_table)
message("Saved: ", pdf_file)
}
