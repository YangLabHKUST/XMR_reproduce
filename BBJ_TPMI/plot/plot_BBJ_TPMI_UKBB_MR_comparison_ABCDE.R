#!/usr/bin/env Rscript

# Compose panels A-E in three rows from the method result tables.
# Row 1: A discovery matrix; B BBJ XMR versus UKBB MR-APSS.
# Row 2: C four BBJ-versus-TPMI effect comparisons with one shared panel label.
# Row 3: D BBJ case studies; E significant-pair counts using each method's available fits.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else getwd()
project_dir <- normalizePath(file.path(script_dir, ".."))
output_dir <- file.path(project_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(output_dir, "BBJ_TPMI_UKBB_MR_comparison_ABCDE")

component_scripts <- file.path(script_dir, c(
  "plot_BBJ_TPMI_UKBB_MR_comparison_ABCD.R",
  "plot_BBJ_XMR_EUR_MRAPSS_effect_size_Bonferroni.R",
  "plot_BBJ_TPMI_effect_size_four_methods_Bonferroni.R"
))
supporting_scripts <- file.path(script_dir, c(
  "plot_BBJ_TPMI_method_significant_pairs.R",
  "plot_BBJ_DBP_MI_BMI_PAD_shared_columns.R"
))
required_scripts <- c(file.path(script_dir, "plot_BBJ_TPMI_UKBB_MR_comparison_ABCDE.R"),
                      component_scripts, supporting_scripts)
stopifnot(all(file.exists(required_scripts)))
component_envs <- lapply(component_scripts, function(path) {
  env <- new.env(parent = globalenv())
  env$render_only <- TRUE
  sys.source(path, envir = env)
  env
})
base <- component_envs[[1L]]
eur <- component_envs[[2L]]
cross <- component_envs[[3L]]

# Match single-cohort significance to A's cohort colours in B/C.
# Both and Neither use separate colours from E's method palette.
# EUR remains an internal data key; all figure labels use UKBB.
scatter_colours <- c(
  "Both significant" = "#7A5BA7",
  "BBJ only" = unname(base$cohort_colours["BBJ"]),
  "EUR only" = unname(base$cohort_colours["EUR"]),
  "TPMI only" = unname(base$cohort_colours["TPMI"]),
  "Neither" = "#9CAEC1"
)
recolour_scatter <- function(p, target) {
  categories <- c("Both significant", "BBJ only", paste(target, "only"), "Neither")
  p <- suppressMessages(p + scale_colour_manual(
    values = scatter_colours[categories], breaks = categories,
    labels = ifelse(categories == "EUR only", "UKBB only", categories),
    drop = FALSE, name = "Significance"
  ))
  p$layers <- lapply(p$layers, function(layer) {
    changes <- if (inherits(layer$geom, "GeomRibbon")) {
      list(fill = "#B8B8B8")
    } else if (inherits(layer$geom, "GeomLine")) {
      list(colour = "#555555")
    } else NULL
    if (!is.null(changes)) {
      return(ggproto(NULL, layer,
                     aes_params = modifyList(layer$aes_params, changes)))
    }
    layer
  })
  p
}

# Verify shared results before combining plots built by different scripts.
for (pop in c("BBJ", "EUR")) {
  old <- base$long_data[cohort == pop, .(exposure_name, outcome_name, beta, p)]
  new <- if (pop == "BBJ") {
    eur$pairs[, .(exposure_name = exposure_name_BBJ, outcome_name = outcome_name_BBJ,
                   beta_new = beta_BBJ, p_new = p_BBJ)]
  } else {
    eur$pairs[, .(exposure_name = exposure_name_EUR, outcome_name = outcome_name_EUR,
                   beta_new = beta_EUR, p_new = p_EUR)]
  }
  check <- merge(old, new, by = c("exposure_name", "outcome_name"))
  stopifnot(nrow(check) == 105L, all(abs(check$beta - check$beta_new) < 1e-10),
            all(abs(check$p - check$p_new) < 1e-12))
}
check_xmr <- merge(
  eur$pairs[, .(exposure, outcome, beta_BBJ, p_BBJ)],
  cross$pairs[method == "XMR", .(exposure, outcome, beta_new = beta_BBJ, p_new = p_BBJ)],
  by = c("exposure", "outcome")
)
stopifnot(nrow(check_xmr) == 105L,
          all(abs(check_xmr$beta_BBJ - check_xmr$beta_new) < 1e-10),
          all(abs(check_xmr$p_BBJ - check_xmr$p_new) < 1e-12))

# Panels D and E are computed directly from the input result tables.
canvas_width <- base$canvas_width
column_widths <- base$column_widths
header_height <- 0.95
top_body_height <- max(base$body_height, 7.70)
top_height <- header_height + top_body_height
middle_height <- 7.25
bottom_height <- base$bottom_height
canvas_height <- top_height + middle_height + bottom_height

panel_header <- function(label, title = NULL) {
  p <- ggdraw() + draw_label(label, x = 0.0025,
                             y = if (is.null(title)) 0.5 else 1 - 0.275 / header_height,
                             hjust = 0, fontface = "bold", size = 25, colour = "#202A33")
  if (!is.null(title)) {
    p <- p + draw_label(title, x = 0.5, y = 1 - 0.75 / header_height,
                         fontface = "bold", size = 25, colour = "#202A33")
  }
  p
}

# Panel A: association matrix and cohort legend.
panel_a_body <- ggdraw() + draw_grob(
  base$panel_a_grob, x = 0, y = 1 - base$a_content_height / top_body_height,
  width = 1, height = base$a_content_height / top_body_height
)
panel_a <- plot_grid(base$panel_a_header, panel_a_body, ncol = 1,
                     rel_heights = c(header_height, top_body_height))
panel_b_plot <- recolour_scatter(eur$figure, "EUR") +
  labs(title = NULL, y = "Effect size in UKBB (MR-APSS)") + theme(
  plot.margin = margin(8, 30, 6, 8), legend.position = "bottom",
  legend.box.margin = margin(t = 7)
)
panel_b <- plot_grid(panel_header("B"),
                     panel_b_plot, ncol = 1,
                     rel_heights = c(header_height, top_body_height))
top_row <- plot_grid(panel_a, panel_b, nrow = 1, rel_widths = column_widths)

# C has one label for all four methods, not four separately lettered subpanels.
# Match each text role to B, including annotations whose size is set in layers.
reference_text_sizes <- list()
for (layer in eur$figure$layers) {
  if (inherits(layer$geom, "GeomTextRepel")) {
    reference_text_sizes$pair <- layer$aes_params$size
  } else if (inherits(layer$geom, "GeomText")) {
    reference_text_sizes$correlation <- layer$aes_params$size
  }
}
method_plots <- lapply(cross$panels, function(p) {
  p <- recolour_scatter(p, "TPMI") + labs(tag = NULL) + theme(
    plot.title = element_text(size = eur$figure$theme$plot.title$size),
    axis.title = element_text(size = eur$figure$theme$axis.title$size),
    axis.text = element_text(size = eur$figure$theme$axis.text$size),
    legend.title = element_text(size = eur$figure$theme$legend.title$size),
    legend.text = element_text(size = eur$figure$theme$legend.text$size)
  )
  p$layers <- lapply(p$layers, function(layer) {
    text_size <- if (inherits(layer$geom, "GeomTextRepel")) {
      reference_text_sizes$pair
    } else if (inherits(layer$geom, "GeomText")) {
      reference_text_sizes$correlation
    } else NULL
    if (!is.null(text_size)) {
      return(ggproto(NULL, layer,
                     aes_params = modifyList(layer$aes_params, list(size = text_size))))
    }
    layer
  })
  p
})
panel_c_legend <- get_legend(method_plots[[1L]] + theme(legend.position = "bottom"))
method_row <- plot_grid(
  plotlist = lapply(method_plots, function(p) p + theme(legend.position = "none")),
  nrow = 1, align = "hv", axis = "tblr"
)
panel_c <- plot_grid(panel_header("C"), method_row, panel_c_legend,
                     ncol = 1, rel_heights = c(0.45, 6.15, 0.65))

# Panels D and E occupy the third row.
panel_d <- plot_grid(panel_header("D"), base$panel_c_plot,
                     ncol = 1, rel_heights = c(0.08, 0.92))
panel_e <- plot_grid(panel_header("E"), base$panel_d_plot,
                     ncol = 1, rel_heights = c(0.08, 0.92))
bottom_row <- plot_grid(panel_d, panel_e, nrow = 1, rel_widths = column_widths)
figure <- plot_grid(top_row, panel_c, bottom_row, ncol = 1,
                    rel_heights = c(top_height, middle_height, bottom_height))

ggsave(paste0(stem, ".pdf"), figure, width = canvas_width, height = canvas_height,
        device = grDevices::pdf, useDingbats = FALSE, bg = "white", limitsize = FALSE,
        title = "MR discoveries, effect comparisons, and case studies across BBJ, TPMI, and UKBB")
ggsave(paste0(stem, ".png"), figure, width = canvas_width, height = canvas_height,
        dpi = 200, bg = "white", limitsize = FALSE)

fwrite(base$panel_a_data, paste0(stem, "_panelA_data.tsv"), sep = "\t")
fwrite(eur$pairs, paste0(stem, "_panelB_data.tsv"), sep = "\t")
fwrite(eur$summary_table, paste0(stem, "_panelB_summary.tsv"), sep = "\t")
fwrite(eur$fit_data, paste0(stem, "_panelB_linear_fit.tsv"), sep = "\t")
fwrite(cross$pairs, paste0(stem, "_panelC_data.tsv"), sep = "\t")
fwrite(cross$summary_table, paste0(stem, "_panelC_summary.tsv"), sep = "\t")
fwrite(cross$fit_data, paste0(stem, "_panelC_linear_fit.tsv"), sep = "\t")
fwrite(base$panel_c_data, paste0(stem, "_panelD_data.tsv"), sep = "\t")
fwrite(base$panel_d_data, paste0(stem, "_panelE_counts.tsv"), sep = "\t")
fwrite(base$count_env$all_data, paste0(stem, "_panelE_all_data.tsv"), sep = "\t")
source_files <- unique(c(unname(base$input_files), base$count_env$input_paths,
                         base$case_env$input_paths, eur$input_paths, cross$input_paths))
fwrite(data.table(path = normalizePath(source_files), md5 = unname(tools::md5sum(source_files))),
        paste0(stem, "_inputs.tsv"), sep = "\t")
fwrite(data.table(path = normalizePath(required_scripts),
                  md5 = unname(tools::md5sum(required_scripts))),
        paste0(stem, "_component_scripts.tsv"), sep = "\t")
capture.output(sessionInfo(), file = paste0(stem, "_sessionInfo.txt"))
pdf_file <- paste0(stem, ".pdf")
stopifnot(identical(eur$input_hashes, tools::md5sum(eur$input_paths)),
          identical(cross$input_hashes, tools::md5sum(cross$input_paths)))
message(sprintf("Saved three-row figure: %.2f x %.2f inches", canvas_width, canvas_height))
message("Panel order: A/B; C (four methods); D/E.")
message("Saved: ", pdf_file)
