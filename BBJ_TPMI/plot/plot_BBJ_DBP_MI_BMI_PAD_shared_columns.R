#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

# Two BBJ examples with shared method and IV-threshold columns.
script_file <- normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(FALSE), value = TRUE
)))
analysis_dir <- dirname(dirname(script_file))
results_dir <- file.path(analysis_dir, "results")
output_dir <- file.path(analysis_dir, "figures")
stem <- "BBJ_DBP_MI_BMI_PAD_shared_columns"
cutoff <- 0.05 / 105
method_order <- c("XMR", "CAMERA", "MRBEE", "MRBEE-TL", "MR-APSS")
method_colours <- c(
  "XMR" = "#D73027", "CAMERA" = "#0072B2", "MRBEE" = "#009E73",
  "MRBEE-TL" = "#E69F00", "MR-APSS" = "#332F2C"
)
input_paths <- file.path(results_dir, c(
  "masterClump5e-05/BBJ/crossPopIntercept0/XMR_EAS_BBJ_EURselected_two_thresholds_noEffectHetQC.tsv",
  "masterClump5e-05/CAMERA_EUR_to_BBJ_TPMI.tsv",
  "masterClump5e-05/MRBEE_TL_EUR_to_BBJ_TPMI.tsv",
  "sharedIV_masterClump5e-05/MRAPSS/BBJ/MRAPSS_BBJ_results.tsv"
))
stopifnot(all(file.exists(input_paths)))
d <- rbindlist(lapply(input_paths, function(path) {
  z <- fread(path)
  z[, source_row := seq_len(.N)]
  z <- z[target == "BBJ" & ((exposure == 4079 & outcome == 411.2) |
                            (exposure == 21001 & outcome == 443.9))]
  if (all(z$method == "XMR")) {
    stopifnot(all(z$background_cross_population_intercept == 0),
              all(z$master_clump_threshold == 5e-05))
  }
  z <- z[method != "MRBEE_source_specific"]
  z <- z[(method %chin% c("XMR", "MRAPSS") & threshold == 5e-05) |
           (!method %chin% c("XMR", "MRAPSS") & threshold == 5e-08)]
  stopifnot(nrow(z) > 0L, all(z$status == "ok"))
  labels <- c(XMR = "XMR", CAMERA = "CAMERA", MRAPSS = "MR-APSS",
              MRBEE_target_specific = "MRBEE", MRBEE_TL_UV = "MRBEE-TL")
  stopifnot(all(z$method %in% names(labels)))
  iv_col <- if ("nIV" %in% names(z)) "nIV" else "n_IV"
  reported_ci <- all(c("CI_low", "CI_high") %in% names(z))
  data.table(
    cohort = "BBJ", pair = fifelse(z$exposure == 4079, "DBP_MI", "BMI_PAD"),
    exposure = z$exposure, outcome = z$outcome,
    method = unname(labels[z$method]), source_method = z$method,
    threshold = z$threshold, n_IV = z[[iv_col]], beta = z$beta, se = z$se, p = z$p,
    ci_lower = if (reported_ci) z$CI_low else z$beta - qnorm(0.975) * z$se,
    ci_upper = if (reported_ci) z$CI_high else z$beta + qnorm(0.975) * z$se,
    ci_source = if (reported_ci) "Reported CI" else "Normal approximation: beta +/- 1.96 SE",
    status = z$status, source_file = normalizePath(path), source_row = z$source_row,
    source_analysis_id = z$analysis_id
  )
}))
stopifnot(
  nrow(d) == 10L, !anyDuplicated(d[, .(pair, method)]),
  all(is.finite(d$beta)), all(is.finite(d$se) & d$se > 0),
  all(is.finite(d$p) & d$p >= 0 & d$p <= 1),
  all(is.finite(d$n_IV) & d$n_IV > 0 & d$n_IV == floor(d$n_IV)),
  all(d$ci_lower < d$beta & d$beta < d$ci_upper),
  all(d[, .(n_thresholds = uniqueN(threshold)), by = method]$n_thresholds == 1L)
)
d[, `:=`(
  y = 6L - match(method, method_order), significant = p < cutoff,
  bonferroni_cutoff = cutoff, threshold_label = sprintf("%.0e", threshold),
  p_label = ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
                   formatC(p, format = "f", digits = 3))
)]
d[, `:=`(
  point_fill = ifelse(significant, unname(method_colours[method]), "white"),
  estimate_label = sprintf("paste('%.3f (', %s('%s'), ')')", beta,
                           ifelse(significant, "bold", "plain"), p_label)
)]
setorder(d, pair, -y)

# Different exposure-outcome scales retain separate, fully labelled axes.
layout <- data.table(
  pair = c("DBP_MI", "BMI_PAD"),
  title = c("bold('DBP' %->% 'Myocardial infarction (MI)')",
            "bold('BMI' %->% 'Peripheral arterial disease (PAD)')"),
  title_x = c(6.55, 14.30), iv_x = c(3.90, 10.75),
  forest_left = c(4.40, 11.30), forest_right = c(6.85, 14.60),
  estimate_x = c(8.60, 16.60),
  effect_min = c(-0.20, -0.08), effect_max = c(0.65, 0.13)
)
d[layout, on = "pair", `:=`(
  iv_x = i.iv_x, forest_left = i.forest_left, forest_right = i.forest_right,
  estimate_x = i.estimate_x,
  effect_min = i.effect_min, effect_max = i.effect_max
)]
stopifnot(all(d$ci_lower > d$effect_min), all(d$ci_upper < d$effect_max))
map_effect <- function(value, lower, upper, left, right) {
  left + (value - lower) / (upper - lower) * (right - left)
}
d[, `:=`(
  effect_x = map_effect(beta, effect_min, effect_max, forest_left, forest_right),
  lo_x = map_effect(ci_lower, effect_min, effect_max, forest_left, forest_right),
  hi_x = map_effect(ci_upper, effect_min, effect_max, forest_left, forest_right)
)]
ticks <- rbindlist(lapply(seq_len(nrow(layout)), function(i) {
  l <- layout[i]
  values <- if (l$pair == "DBP_MI") seq(-0.2, 0.6, 0.2) else seq(-0.05, 0.10, 0.05)
  data.table(x = map_effect(values, l$effect_min, l$effect_max,
                           l$forest_left, l$forest_right),
             label = sprintf(if (l$pair == "DBP_MI") "%.1f" else "%.2f", values))
}))
layout[, zero_x := map_effect(0, effect_min, effect_max, forest_left, forest_right)]
shared <- unique(d[, .(method, y, threshold_label)])
ink <- "#202A33"
font_pt <- 19.5
text_size <- font_pt / ggplot2::.pt
bands <- data.table(y = c(2, 4))

panel_plot <- ggplot() +
  geom_rect(data = bands, aes(ymin = y - 0.48, ymax = y + 0.48),
            xmin = 0, xmax = 18.10, fill = "#F4F6F8", colour = NA) +
  geom_text(data = layout, aes(x = title_x, label = title), y = 7.35,
            parse = TRUE, size = 21 / ggplot2::.pt, colour = ink) +
  annotate("segment", x = 0, xend = 18.10, y = 5.65, yend = 5.65,
           colour = "#B8BEC4", linewidth = 0.4) +
  annotate("segment", x = 10.35, xend = 10.35, y = 0.55, yend = 7.70,
           colour = "#D9DEE2", linewidth = 0.35) +
  geom_segment(data = layout, aes(x = zero_x, xend = zero_x), y = 0.55, yend = 5.43,
               colour = "#92999F", linetype = "dashed", linewidth = 0.45) +
  geom_segment(data = d, aes(x = lo_x, xend = hi_x, y = y, yend = y, colour = method),
               linewidth = 0.85) +
  geom_segment(data = d, aes(x = lo_x, xend = lo_x, y = y - 0.1,
                            yend = y + 0.1, colour = method), linewidth = 0.6) +
  geom_segment(data = d, aes(x = hi_x, xend = hi_x, y = y - 0.1,
                            yend = y + 0.1, colour = method), linewidth = 0.6) +
  geom_point(data = d, aes(x = effect_x, y = y, colour = method, fill = point_fill),
             shape = 21, size = 3.4, stroke = 0.95) +
  geom_text(data = shared, aes(y = y, label = method), x = 0.10,
            hjust = 0, size = text_size, colour = ink) +
  geom_text(data = shared, aes(y = y, label = threshold_label), x = 2.60,
            size = text_size, colour = ink) +
  geom_text(data = d, aes(x = iv_x, y = y, label = n_IV), size = text_size, colour = ink) +
  geom_text(data = d, aes(x = estimate_x, y = y, label = estimate_label),
            parse = TRUE, size = text_size, colour = ink) +
  annotate("text", x = 0.10, y = 6.18, label = "Method", hjust = 0,
           fontface = "bold", size = text_size, colour = ink) +
  annotate("text", x = 2.60, y = 6.18, label = "IV\nthreshold", lineheight = 0.95,
           fontface = "bold", size = text_size, colour = ink) +
  geom_text(data = layout, aes(x = iv_x), y = 6.18, label = "IVs",
            fontface = "bold", size = text_size, colour = ink) +
  geom_text(data = layout, aes(x = (forest_left + forest_right) / 2),
            y = 6.18, label = "Effect estimate\n(95% CI)", lineheight = 0.95,
            fontface = "bold", size = text_size, colour = ink) +
  geom_text(data = layout, aes(x = estimate_x), y = 6.18,
            label = "hat(beta)~'('*italic(P)*')'", parse = TRUE,
            size = text_size, colour = ink) +
  geom_segment(data = layout, aes(x = forest_left, xend = forest_right),
               y = 0.30, yend = 0.30, colour = "#92999F", linewidth = 0.4) +
  geom_segment(data = ticks, aes(x = x, xend = x), y = 0.30, yend = 0.17,
               colour = "#92999F", linewidth = 0.4) +
  geom_text(data = ticks, aes(x = x, label = label), y = -0.10,
            size = 16.5 / ggplot2::.pt, colour = ink) +
  scale_colour_manual(values = method_colours, guide = "none") +
  scale_fill_identity() +
  scale_x_continuous(limits = c(-0.05, 18.20), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(-0.65, 7.90), expand = expansion(mult = 0)) +
  theme_void(base_family = "Helvetica") +
  theme(plot.margin = margin(0, 7, 0, 7))

if (!exists("render_only", inherits = FALSE) || !render_only) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  figure <- ggdraw() +
    draw_label("BBJ", x = 0.5, y = 0.96, fontface = "bold", size = 25, colour = ink) +
    draw_plot(panel_plot, x = 0.01, y = 0.01, width = 0.98, height = 0.90)
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  ggsave(pdf_path, figure, width = 15.5, height = 6.5,
         device = grDevices::pdf, family = "Helvetica", useDingbats = FALSE, bg = "white")
  ggsave(file.path(output_dir, paste0(stem, ".png")), figure,
         width = 15.5, height = 6.5, dpi = 220, bg = "white")
  fwrite(d, file.path(output_dir, paste0(stem, "_data.tsv")), sep = "\t")
  fwrite(data.table(source_file = input_paths, md5 = unname(tools::md5sum(input_paths))),
         file.path(output_dir, paste0(stem, "_inputs.tsv")), sep = "\t")
  print(d[, .(pair, method, threshold, n_IV, beta, se, p)])
  message("Saved: ", pdf_path)
}
