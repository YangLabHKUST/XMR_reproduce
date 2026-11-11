#!/usr/bin/env Rscript

# BBJ target effects with FinnGen versus PanUKB auxiliary data.
# Define a common trait-pair set using XMR availability across all thresholds.
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
output_dir <- file.path(project_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stem <- "Auxiliary_sensitivity_four_methods_two_thresholds_effect_overlap"
thresholds <- c("5e-05", "5e-08")
all_xmr_thresholds <- c("5e-05", "5e-06", "5e-07", "5e-08")
method_names <- c(XMR = "XMR", CAMeRa = "CAMERA",
                  MRBEE_target_specific = "MRBEE", MRBEE_TL_UV = "MRBEE-TL")
method_levels <- unname(method_names)
# Use the same significance-category palette for every method.
category_colours <- c("Both" = "#332288", "FinnGen only" = "#6B1636",
                      "PanUKB only" = "#E6D84A", "Neither" = "#C3C7CD")
category_alpha <- 0.75
legend_size <- 30
n_tests <- 168L
significance_cutoff <- 0.05 / n_tests
category_levels <- c("Both", "FinnGen only", "PanUKB only", "Neither")

# Read the experiment outputs for the four methods.
input_paths <- file.path(project_dir, "results", c(
  "crossPopIntercept0/XMR_Fin_PanUKB_to_BBJ_68_traits_FEMA_P_target_clump_noEffectHetQC.tsv",
  "masterClump5e-05/CAMeRa_Fin_PanUKB_to_BBJ_68_traits.tsv",
  "masterClump5e-05/MRBEE_TL_Fin_PanUKB_to_BBJ_68_traits.tsv"
))
stopifnot(all(file.exists(input_paths)))
input_hashes <- tools::md5sum(input_paths)
raw <- rbindlist(lapply(input_paths, function(path) {
  d <- fread(path, na.strings = c("NA", "", "NaN"))
  stopifnot(all(c("method", "auxiliary", "exposure", "outcome", "threshold",
                  "beta", "se", "p", "status") %in% names(d)))
  d[, `:=`(source_file = normalizePath(path), source_row = seq_len(.N))]
  if (!"CI_low" %in% names(d)) {
    d[, `:=`(CI_low = beta - 1.96 * se, CI_high = beta + 1.96 * se)]
  }
  d
}), fill = TRUE)
raw[, threshold_key := sprintf("%.0e", threshold)]
raw[, valid := !is.na(status) & status == "ok" & is.finite(beta) &
      is.finite(se) & se > 0 & is.finite(p) & p >= 0 & p <= 1 &
      is.finite(CI_low) & is.finite(CI_high) & CI_low <= beta & beta <= CI_high]
stopifnot(!anyNA(raw$valid), !anyDuplicated(raw,
  by = c("method", "auxiliary", "exposure", "outcome", "threshold_key")))

xmr_raw <- raw[method == "XMR"]
stopifnot(setequal(xmr_raw$threshold_key, all_xmr_thresholds),
          setequal(xmr_raw$auxiliary, c("Fin", "PanUKB")),
          setequal(xmr_raw$exposure, c("BMI", "BW", "Hei")))
outcome_checks <- xmr_raw[, .(
  n_rows = .N, n_exposures = uniqueN(exposure), n_thresholds = uniqueN(threshold_key),
  n_auxiliary = uniqueN(auxiliary), all_valid = all(valid)
), by = outcome]
outcome_checks[, included := n_rows == 24L & n_exposures == 3L &
                 n_thresholds == 4L & n_auxiliary == 2L & all_valid]
reference_pairs <- unique(xmr_raw[outcome %chin% outcome_checks[included == TRUE, outcome],
                                  .(exposure, outcome)])
stopifnot(sum(outcome_checks$included) == 56L, nrow(reference_pairs) == n_tests)

selected <- merge(raw[method %chin% names(method_names) & threshold_key %chin% thresholds],
                  reference_pairs, by = c("exposure", "outcome"))
selected[, method_label := unname(method_names[method])]
availability <- selected[, .(n_reported = .N, n_valid = sum(valid)),
                         by = .(method_label, auxiliary, threshold_key)]
stopifnot(nrow(availability) == 16L, all(availability$n_reported == n_tests),
          all(availability$n_valid == n_tests), all(selected$target == "BBJ"))
keys <- c("method_label", "threshold_key", "exposure", "outcome")
values <- c("beta", "se", "p", "CI_low", "CI_high", "source_file", "source_row")
fin <- selected[auxiliary == "Fin", c(keys, values), with = FALSE]
pan <- selected[auxiliary == "PanUKB", c(keys, values), with = FALSE]
setnames(fin, values, paste0(values, "_fin"))
setnames(pan, values, paste0(values, "_panukb"))
pairs <- merge(fin, pan, by = keys, all = TRUE)
setnames(pairs, "method_label", "method")
stopifnot(nrow(pairs) == 8L * n_tests, !anyNA(pairs))
pairs[, `:=`(significant_fin = p_fin < significance_cutoff,
             significant_panukb = p_panukb < significance_cutoff,
             pair_label = paste(exposure, outcome, sep = "-"),
             bonferroni_threshold = significance_cutoff)]
pairs[, category := factor(fcase(
  significant_fin & significant_panukb, "Both",
  significant_fin, "FinnGen only", significant_panukb, "PanUKB only",
  default = "Neither"
), levels = category_levels)]
stats <- pairs[, .(
  n_pairs = .N, pearson_r = cor(beta_fin, beta_panukb),
  both = sum(significant_fin & significant_panukb),
  fin_only = sum(significant_fin & !significant_panukb),
  panukb_only = sum(!significant_fin & significant_panukb),
  neither = sum(!significant_fin & !significant_panukb)
), by = .(method, threshold_key)]
stats[, `:=`(fin_total = both + fin_only, panukb_total = both + panukb_only,
             union = both + fin_only + panukb_only)]
stats <- stats[order(match(threshold_key, thresholds), match(method, method_levels))]
stopifnot(all(stats$n_pairs == n_tests), all(stats$union + stats$neither == n_tests))

bar_data <- rbindlist(list(
  stats[, .(method, threshold_key, category = "Both", n = both, total = union)],
  stats[, .(method, threshold_key, category = "FinnGen only", n = fin_only, total = union)],
  stats[, .(method, threshold_key, category = "PanUKB only", n = panukb_only, total = union)]
))
bar_data[, category := factor(category, levels = category_levels)]
bar_data <- bar_data[order(match(threshold_key, thresholds), match(method, method_levels), category)]
bar_data[, bottom := cumsum(n) - n, by = .(method, threshold_key)]
stopifnot(all(bar_data[, .(valid = sum(n) == unique(total)),
                       by = .(method, threshold_key)]$valid))
effect_limits <- c(-0.22, 0.52)
effect_span <- diff(effect_limits)
effect_breaks <- c(-0.2, 0, 0.2, 0.4)
all_ci <- unlist(pairs[, .(CI_low_fin, CI_high_fin, CI_low_panukb, CI_high_panukb)])
stopifnot(all(all_ci > effect_limits[1L] & all_ci < effect_limits[2L]))
bar_max <- max(5, ceiling((max(bar_data$total) + 2) / 5) * 5)

make_panel <- function(m, th, tag) {
  show_axis_titles <- m == method_levels[1L]
  d <- copy(pairs[method == m & threshold_key == th])
  a <- stats[method == m & threshold_key == th]
  palette <- category_colours
  d[, draw_order := fifelse(category == "Neither", 1L, fifelse(category == "Both", 3L, 2L))]
  setorder(d, draw_order)
  # Mark only shared discoveries with labels, avoiding crowded labels for TL-only hits.
  both_points <- d[category == "Both"]
  bar_min <- effect_limits[2L] + 0.015 * effect_span
  bar_end <- effect_limits[2L] + 0.25 * effect_span
  segments <- copy(bar_data[method == m & threshold_key == th])
  segments[, `:=`(bar_x = bar_min + 0.115 * effect_span,
                  fill = unname(palette[as.character(category)]))]
  segments[, `:=`(
    xmin = bar_x - 0.045 * effect_span, xmax = bar_x + 0.045 * effect_span,
    ymin = effect_limits[1L] + bottom / bar_max * effect_span,
    ymax = effect_limits[1L] + (bottom + n) / bar_max * effect_span,
    ymid = effect_limits[1L] + (bottom + n / 2) / bar_max * effect_span
  )]
  # Put small-count labels beside thin segments so enlarged labels do not overlap.
  thin_segments <- copy(segments[n > 0 & n <= 2])
  thin_segments[, `:=`(
    label_x = bar_x + fifelse(category == "FinnGen only", -0.10, 0.10) * effect_span,
    edge_x = bar_x + fifelse(category == "FinnGen only", -0.045, 0.045) * effect_span
  )]
  b <- unique(segments[, .(bar_x, total)])

  p <- ggplot(d, aes(beta_fin, beta_panukb, colour = category)) +
    geom_rect(data = segments[n > 0], aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                         ymax = ymax, fill = I(fill)),
               inherit.aes = FALSE, colour = "white", linewidth = 0.5,
               alpha = category_alpha) +
    geom_label(data = segments[n > 2], aes(x = bar_x, y = ymid, label = n),
                inherit.aes = FALSE, fill = "white", colour = "#222222", linewidth = 0,
                label.padding = grid::unit(0.04, "lines"), size = 25 / .pt,
                fontface = "bold", show.legend = FALSE) +
    geom_segment(data = thin_segments,
                 aes(x = edge_x, xend = label_x, y = ymid, yend = ymid),
                 inherit.aes = FALSE, colour = "#666666", linewidth = 0.25) +
    geom_label(data = thin_segments, aes(x = label_x, y = ymid, label = n),
               inherit.aes = FALSE, fill = "white", colour = "#222222", linewidth = 0,
               label.padding = grid::unit(0.04, "lines"), size = 25 / .pt,
               fontface = "bold", show.legend = FALSE) +
    geom_text(data = b, aes(x = bar_x, y = effect_limits[1L] + total / bar_max * effect_span,
                           label = total), inherit.aes = FALSE, vjust = -0.8,
                colour = "#222222", size = 25 / .pt, fontface = "bold") +
    annotate("segment", x = effect_limits[1L], xend = effect_limits[2L],
             y = effect_limits[1L], yend = effect_limits[1L], linewidth = 0.45) +
    annotate("segment", x = effect_limits[1L], xend = effect_limits[2L],
             y = 0, yend = 0, colour = "#CCCCCC", linewidth = 0.3, linetype = "dotted") +
    annotate("segment", x = 0, xend = 0, y = effect_limits[1L], yend = effect_limits[2L],
             colour = "#CCCCCC", linewidth = 0.3, linetype = "dotted") +
    annotate("segment", x = effect_limits[1L], xend = effect_limits[2L],
             y = effect_limits[1L], yend = effect_limits[2L],
             colour = "#555555", linewidth = 0.45, linetype = "dashed") +
    geom_segment(aes(x = CI_low_fin, xend = CI_high_fin, yend = beta_panukb),
                  alpha = 0.20, linewidth = 0.23, show.legend = FALSE) +
    geom_segment(aes(y = CI_low_panukb, yend = CI_high_panukb, xend = beta_fin),
                  alpha = 0.20, linewidth = 0.23, show.legend = FALSE) +
    geom_point(data = d[category == "Neither"], size = 1.9,
               alpha = category_alpha, stroke = 0) +
    geom_point(data = d[category %in% c("FinnGen only", "PanUKB only")], size = 3.9,
               alpha = category_alpha, stroke = 0) +
    geom_point(data = both_points, size = 4.6, alpha = category_alpha, stroke = 0) +
    geom_text_repel(data = both_points, aes(label = pair_label), colour = "#303030",
                   size = 20 / .pt, box.padding = 0.28, point.padding = 0.15,
                   point.size = 4.6,
                   min.segment.length = 0, segment.colour = "#888888", segment.size = 0.2,
                   force = 2, force_pull = 0.2, max.overlaps = Inf,
                   max.time = Inf, max.iter = 60000, seed = 20261027,
                   xlim = effect_limits, ylim = c(effect_limits[1L], effect_limits[2L] - 0.09),
                   show.legend = FALSE) +
    annotate("text", x = effect_limits[1L] + 0.03 * effect_span,
             y = effect_limits[2L] - 0.03 * effect_span,
             label = sprintf("Pearson r = %.3f", a$pearson_r),
             hjust = 0, vjust = 1, size = 25 / .pt, colour = "#222222") +
    scale_colour_manual(values = palette, guide = "none") +
    scale_x_continuous(breaks = effect_breaks, expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = effect_breaks, expand = expansion(mult = 0)) +
    coord_fixed(ratio = 1, xlim = c(effect_limits[1L], bar_end), ylim = effect_limits,
                 expand = FALSE, clip = "off") +
    # Empty two-line labels preserve identical panel geometry in columns 2--4.
    labs(x = if (show_axis_titles) "Effect estimate\n(FinnGen as auxiliary)" else "\n",
         y = if (show_axis_titles) "Effect estimate\n(PanUKB as auxiliary)" else "\n") +
    theme_classic(base_size = 22, base_family = "Helvetica") +
    theme(axis.title = element_text(size = 25, lineheight = 0.95),
          axis.title.x = element_text(hjust = 0.5 / 1.25, margin = margin(t = 8)),
          axis.title.y.left = element_text(margin = margin(r = 8)),
          axis.text = element_text(size = 20, colour = "#333333"),
          axis.line.x = element_blank(), axis.line.y = element_line(linewidth = 0.45),
          axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks = element_line(linewidth = 0.4),
          plot.margin = margin(6, 12, 24, 10), legend.position = "none")
  if (tag %in% c("A", "E")) {
    p <- p + annotate("text", x = b$bar_x,
                       y = effect_limits[1L] + (b$total / bar_max + 0.16) * effect_span,
                       label = "Significant\npairs", vjust = 0, lineheight = 0.95,
                       size = 20 / .pt, colour = "black")
  }
  header <- ggdraw()
  if (tag %in% c("A", "E")) {
    header <- header + draw_label(paste0("IV threshold: ", th), x = 0.50, y = 0.5,
                                   fontface = "bold", size = 31, colour = "black")
  }
  plot_grid(header, p, ncol = 1, rel_heights = c(0.12, 1))
}

# A--D identify method columns; both thresholds share the same column label.
column_headers <- plot_grid(plotlist = lapply(seq_along(method_levels), function(j) {
  ggdraw() +
    draw_label(LETTERS[j], x = 0.02, y = 0.65, hjust = 0,
               fontface = "bold", size = 29, colour = "black") +
    draw_label(method_levels[j], x = 0.5, y = 0.5, fontface = "bold", size = 30,
               colour = "black")
}), nrow = 1)
rows <- lapply(seq_along(thresholds), function(i) {
  ps <- lapply(seq_along(method_levels), function(j) {
    make_panel(method_levels[j], thresholds[i], LETTERS[(i - 1L) * 4L + j])
  })
  plot_grid(plotlist = ps, nrow = 1, align = "hv", axis = "tblr")
})

# All methods use the same significance-category colours; compact, centred legend.
legend <- ggdraw() +
  draw_label("Significance", x = 0.195, y = 0.54, hjust = 0,
              fontface = "bold", size = legend_size) +
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = category_colours["Both"],
                                          alpha = category_alpha, col = NA)),
            x = 0.286, y = 0.36, width = 0.018, height = 0.36) +
  draw_label("Both", x = 0.310, y = 0.54, hjust = 0, size = legend_size) +
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = category_colours["FinnGen only"],
                                          alpha = category_alpha, col = NA)),
            x = 0.359, y = 0.36, width = 0.018, height = 0.36) +
  draw_label("FinnGen as auxiliary only", x = 0.383, y = 0.54, hjust = 0, size = legend_size) +
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = category_colours["PanUKB only"],
                                          alpha = category_alpha, col = NA)),
            x = 0.552, y = 0.36, width = 0.018, height = 0.36) +
  draw_label("PanUKB as auxiliary only", x = 0.576, y = 0.54, hjust = 0, size = legend_size) +
  draw_grob(grid::circleGrob(gp = grid::gpar(fill = "#C3C7CD", col = NA)),
            x = 0.744, y = 0.36, width = 0.014, height = 0.36) +
  draw_label("Neither", x = 0.764, y = 0.54, hjust = 0, size = legend_size)
figure <- plot_grid(column_headers, rows[[1L]], rows[[2L]], legend,
                    ncol = 1, rel_heights = c(0.10, 1, 1, 0.12))
figure <- ggdraw(figure)
for (column_boundary in c(0.25, 0.5, 0.75)) {
  figure <- figure + draw_line(x = rep(column_boundary, 2), y = c(0.065, 0.952),
                                colour = "#A8ADB3", linetype = "dashed", linewidth = 0.6)
}

pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
png_file <- file.path(output_dir, paste0(stem, ".png"))
ggsave(pdf_file, figure, width = 32, height = 17, device = grDevices::pdf,
        useDingbats = FALSE, bg = "white", title = "Auxiliary sensitivity across four MR methods")
ggsave(png_file, figure, width = 32, height = 17, dpi = 220, bg = "white")
fwrite(pairs, file.path(output_dir, paste0(stem, "_data.tsv")), sep = "\t")
fwrite(stats, file.path(output_dir, paste0(stem, "_summary.tsv")), sep = "\t")
fwrite(bar_data, file.path(output_dir, paste0(stem, "_bar_counts.tsv")), sep = "\t")
fwrite(outcome_checks, file.path(output_dir, paste0(stem, "_outcome_checks.tsv")), sep = "\t")
fwrite(data.table(source_file = normalizePath(input_paths), md5 = unname(input_hashes)),
        file.path(output_dir, paste0(stem, "_inputs.tsv")), sep = "\t")
capture.output(sessionInfo(), file = file.path(output_dir, paste0(stem, "_sessionInfo.txt")))
stopifnot(identical(tools::md5sum(input_paths), input_hashes))
print(stats)
message("Saved: ", pdf_file)
