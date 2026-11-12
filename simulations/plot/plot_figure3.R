#!/usr/bin/env Rscript
# Draw Figure 3 from the tables written by prepare_figure3.R.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(cowplot)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript plot_figure3.R TABLE_DIR FIGURE_DIR")
table_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
qq <- fread(file.path(table_dir, "qq.tsv"))
power <- fread(file.path(table_dir, "power_bars.tsv"))
curves <- fread(file.path(table_dir, "power_curves.tsv"))
validity <- fread(file.path(table_dir, "fit_counts.tsv"))
methods <- c("XMR_FEMA", "CAMeRa_FEMA", "MRBEE_SOURCE", "MRBEE_TARGET",
             "MRBEE_TL", "TEMR", "MRAPSS_TARGET")
method_labels <- c(XMR_FEMA = "XMR", CAMeRa_FEMA = "CAMERA",
    MRBEE_SOURCE = "MRBEE (EUR)", MRBEE_TARGET = "MRBEE (EAS)",
    MRBEE_TL = "MRBEE-TL", TEMR = "TEMR", MRAPSS_TARGET = "MRAPSS")
method_colours <- c(XMR_FEMA = "#D73027", CAMeRa_FEMA = "#0072B2",
    MRBEE_SOURCE = "#707070", MRBEE_TARGET = "#009E73", MRBEE_TL = "#E69F00",
    TEMR = "#8E44AD", MRAPSS_TARGET = "#332F2C")
method_shapes <- c(XMR_FEMA = 16, CAMeRa_FEMA = 17, MRBEE_SOURCE = 1,
    MRBEE_TARGET = 15, MRBEE_TL = 4, TEMR = 18, MRAPSS_TARGET = 8)
series_levels <- c("XMR_5e05", "MRAPSS_TARGET", "CAMeRa_FEMA", "MRBEE_TARGET")
power_levels <- list(Baseline = series_levels[1:2], Reduced = series_levels)
series_labels <- c(XMR_5e05 = "XMR (5e-05)",
    MRAPSS_TARGET = "MRAPSS", CAMeRa_FEMA = "CAMERA (5e-08)",
    MRBEE_TARGET = "MRBEE (EAS; 5e-08)")
series_colours <- c(XMR_5e05 = "#D73027",
    MRAPSS_TARGET = "#332F2C", CAMeRa_FEMA = "#0072B2", MRBEE_TARGET = "#009E73")

base_theme <- theme_classic(base_size = 12, base_family = "sans") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 11, colour = "#333333"),
        axis.line = element_line(linewidth = .4, colour = "#606060"),
        axis.ticks = element_line(linewidth = .35, colour = "#606060"),
        plot.subtitle = element_text(size = 10, colour = "#666666", hjust = .5,
                                     margin = margin(b = 9)),
        legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 17),
        legend.key.width = grid::unit(.75, "cm"), plot.margin = margin(10, 15, 10, 10))
n <- 50L; ranks <- seq_len(n)
envelope <- data.table(x = -log10((ranks - .5) / n),
    lo = -log10(qbeta(.975, ranks, n - ranks + 1)),
    hi = pmin(-log10(qbeta(.025, ranks, n - ranks + 1)), 3))
setorder(envelope, x)

make_qq <- function(g, nt, th, legend = FALSE) {
  d <- qq[profile == g & n_target == nt & threshold == th]
  shown_methods <- if (th == 5e-8) setdiff(methods, "MRAPSS_TARGET") else methods
  counts <- validity[profile == g & n_target == nt & threshold == th & beta_true == 0 &
                       method %in% shown_methods]
  stopifnot(nrow(counts) == length(shown_methods), all(counts$n_valid == 50L))
  p <- ggplot(d, aes(expected_log10, observed_display, colour = method, shape = method)) +
    geom_ribbon(data = envelope, aes(x = x, ymin = lo, ymax = hi),
                inherit.aes = FALSE, fill = "#E4E7EA") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#454545", linewidth = .5) +
    geom_point(size = 2.1, stroke = .65) +
    scale_colour_manual(values = method_colours, limits = shown_methods, labels = method_labels[shown_methods], drop = FALSE) +
    scale_shape_manual(values = method_shapes, limits = shown_methods, labels = method_labels[shown_methods], drop = FALSE) +
    scale_x_continuous(breaks = seq(0, 2.5, .5), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 3, .5), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 3), expand = FALSE, clip = "off") +
    labs(
         x = expression(Expected~~-log[10](italic(P))),
         y = expression(Observed~~-log[10](italic(P)))) + base_theme
  if (legend) p <- p + theme(legend.position = "bottom", legend.text = element_text(size = 18)) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3)),
           shape = guide_legend(nrow = 1))
  p
}
make_power <- function(g, nt, legend = FALSE) {
  d <- copy(power[profile == g & n_target == nt])
  lev <- power_levels[[g]]
  d[, series := factor(series, levels = lev)]
  pd <- position_dodge(width = .8)
  p <- ggplot(d, aes(factor(beta_true, levels=c(0,.1,.2,.3,.4)), power, fill=series, group=series)) +
    geom_col(position=pd, width=.72) +
    geom_hline(yintercept=.05,linetype="dotted",colour="#888888",linewidth=.4) +
    geom_errorbar(aes(ymin=lower,ymax=upper),position=pd,width=.15,linewidth=.4,colour="#333333",show.legend=FALSE) +
    scale_fill_manual(values=series_colours,limits=lev,labels=series_labels[lev],drop=FALSE,name=NULL) +
    scale_x_discrete(labels=sprintf("%.1f",c(0,.1,.2,.3,.4)),drop=FALSE) +
    scale_y_continuous(breaks=seq(0,1,.2),limits=c(0,1),expand=expansion(mult=c(.025,.025))) +
    labs(x="True causal effect",y="Power") + base_theme
  if(legend) p <- p + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2,byrow=TRUE))
  p
}

# Fourth column: the same target-only MR-APSS reference appears in all three panels.
curve_levels <- c("XMR (N aux = 20K)", "XMR (N aux = 40K)", "XMR (N aux = 60K)", "MRAPSS")
curve_labels <- c("XMR (N aux = 20K, N target = 6K)", "XMR (N aux = 40K, N target = 6K)",
                  "XMR (N aux = 60K, N target = 6K)", "MRAPSS")
curves[, series := factor(fifelse(method == "XMR_FEMA", sprintf("XMR (N aux = %sK)", n_aux / 1000),
                                 "MRAPSS"), levels = curve_levels)]
make_curve <- function(correlation, legend = FALSE) {
  d <- copy(curves[rho == correlation | method == "MRAPSS_TARGET"])
  setorder(d, series, beta_true)
  ggplot(d, aes(beta_true, power, colour = series, shape = series, linetype = series, group = series)) +
    geom_hline(yintercept = .05, linetype = "dotted", colour = "#888888", linewidth = .4) +
    geom_line(linewidth = .7, alpha = 1) +
    geom_point(size = 3.5, stroke = .8, alpha = 1) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .009, linewidth = .45,
                  alpha = 1, show.legend = FALSE) +
    scale_colour_manual(values = setNames(c(rep("#D73027", 3), "#332F2C"), curve_levels),
                        limits = curve_levels, labels = curve_labels, drop = FALSE) +
    scale_shape_manual(values = setNames(c(17, 16, 15, 8), curve_levels),
                       limits = curve_levels, labels = curve_labels, drop = FALSE) +
    scale_linetype_manual(values = setNames(c("dashed", "solid", "dotdash", "dotted"), curve_levels),
                          limits = curve_levels, labels = curve_labels, drop = FALSE) +
    scale_x_continuous(breaks = seq(0, .4, .1), labels = sprintf("%.1f", seq(0, .4, .1)),
                       limits = c(0, .4), expand = expansion(mult = c(.06, .06))) +
    scale_y_continuous(breaks = seq(0, 1, .2), limits = c(0, 1),
                       expand = expansion(mult = c(.025, .025))) +
    labs(x = "True causal effect", y = "Power", colour = NULL, shape = NULL, linetype = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3.5, alpha = 1)),
           shape = guide_legend(ncol = 1, byrow = TRUE),
           linetype = guide_legend(ncol = 1, byrow = TRUE)) +
    theme_classic(base_size = 11, base_family = "sans") +
    theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10, colour = "#333333"),
          axis.title.x = element_text(margin = margin(t = 8)),
          axis.title.y = element_text(margin = margin(r = 8)),
          axis.line = element_line(linewidth = .4, colour = "#606060"),
          axis.ticks = element_line(linewidth = .35, colour = "#606060"),
          legend.position = if (legend) "inside" else "none",
          legend.position.inside = c(.01, .98), legend.justification = c(0, 1),
          legend.background = element_rect(fill = scales::alpha("white", .85), colour = NA),
          legend.text = element_text(size = 16), legend.key.width = grid::unit(.9, "cm"),
          legend.spacing.x = grid::unit(.15, "cm"), plot.margin = margin(14, 22, 5, 9))
}

# Assemble only the twelve panels used in Figure 3, in reading order.
settings <- data.table(profile = c("Baseline", "Baseline", "Reduced"),
                       n_target = c(6000L, 8000L, 6000L), background = c(.4, .4, .2),
                       curve_rho = c(.2, .5, .8))
plots <- unlist(lapply(seq_len(nrow(settings)), function(i) {
  row <- settings[i]
  list(make_qq(row$profile, row$n_target, 5e-5), make_qq(row$profile, row$n_target, 5e-8),
       make_power(row$profile, row$n_target), make_curve(row$curve_rho, legend = i == 1L))
}), recursive = FALSE)
grobs <- lapply(plots, ggplotGrob)
widths <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
grobs <- lapply(grobs, function(g) {g$widths <- widths; g$respect <- FALSE; g})
make_row <- function(i) {
  ix <- (i - 1L) * 4L + 1:4
  plot_grid(plotlist = grobs[ix], nrow = 1, labels = LETTERS[ix],
            label_size = 20, label_x = .005, label_y = .995, hjust = 0, vjust = 1)
}
strip <- function(text, size = 18, fill = NULL) {
  p <- ggdraw()
  if (!is.null(fill)) p <- p + draw_grob(grid::rectGrob(gp = grid::gpar(fill = fill, col = NA)))
  p + draw_label(text, size = size, fontface = "bold")
}
cols <- plot_grid(strip("IV threshold = 5e-05"), strip("IV threshold = 5e-08"),
                   strip("Power comparison"), strip("Varying N aux; N target = 6K"), nrow = 1)
qq_legend <- get_legend(make_qq("Baseline", 6000, 5e-5, TRUE))
qq_legend_row <- plot_grid(qq_legend, NULL, nrow = 1, rel_widths = c(3, 1))
bar_legend <- get_legend(make_power("Reduced", 6000, TRUE) +
    theme(legend.text = element_text(size = 18), legend.margin = margin(0, 0, 0, 0)) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)))
bar_legend_row <- plot_grid(bar_legend, NULL, nrow = 1, rel_widths = c(3, 1))
row_heading <- function(i) {
  row <- settings[i]
  plot_grid(strip(paste0("N aux = 40K | N target = ", row$n_target / 1000,
      "K | Cross-pop correlation = 0.8 | EAS background X-Y correlation = ", row$background),
      17, if (row$background == .2) "#EDF5F2" else "#F0F1F3"),
    strip(paste0("Cross-pop correlation = ", row$curve_rho), 17, "#F0F1F3"),
    nrow = 1, rel_widths = c(3, 1))
}
figure <- plot_grid(cols, row_heading(1), make_row(1), row_heading(2), make_row(2), qq_legend_row,
                    row_heading(3), make_row(3), bar_legend_row, ncol = 1,
                    rel_heights = c(.45, .4, 4.65, .4, 4.65, .65, .4, 4.65, .65))
figure <- ggdraw(figure) + draw_line(x = c(.75, .75), y = c(.04, .995),
                                   colour = "#B8B8B8", linetype = "dashed", linewidth = .5)
ggsave(file.path(out_dir, "Figure_3.pdf"), figure, width = 22, height = 17.5,
       device = grDevices::pdf, useDingbats = FALSE, bg = "white")
ggsave(file.path(out_dir, "Figure_3.png"), figure, width = 22, height = 17.5, dpi = 300, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Saved Figure_3.pdf and Figure_3.png in ", out_dir)
