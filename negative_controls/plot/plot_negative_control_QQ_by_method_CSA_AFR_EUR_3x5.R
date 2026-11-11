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
default_input_dir <- normalizePath(file.path(script_dir, "..", "results"))
default_output_dir <- normalizePath(file.path(script_dir, "..", "figures"))

input_dir <- if (length(args) >= 1L) args[[1L]] else default_input_dir
output_dir <- if (length(args) >= 2L) args[[2L]] else default_output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

main_method_levels <- c(
  "TEMR",
  "CAMERA",
  "MRBEE",
  "MRBEE-TL",
  "XMR"
)
eur_method_levels <- c(
  "CAMERA (CSA + EUR FEMA IV)",
  "CAMERA (AFR + EUR FEMA IV)",
  "MRBEE (CSA + EUR FEMA IV)",
  "MRBEE (AFR + EUR FEMA IV)",
  "MR-APSS"
)
all_method_levels <- c(main_method_levels, eur_method_levels)

threshold_levels <- c("5e-08", "5e-05")
# Muted berry and dark forest green encode IV thresholds, not method identity.
threshold_colours <- c(
  "5e-08" = "#CC6677",
  "5e-05" = "#24553A"
)
threshold_shapes <- c(
  "5e-08" = 21,
  "5e-05" = 21
)
threshold_labels <- c(
  "5e-08" = expression(5 %*% 10^{-8}),
  "5e-05" = expression(5 %*% 10^{-5})
)

required_files <- c(
  # XMR: FEMA-P ranking, target LD clumping, no effect heterogeneity QC.
  XMR = file.path("crossPopIntercept0", "NC_XMR_EUR_to_AFR_CSA_FEMA_P_target_clump_noEffectHetQC_clump10000kb.tsv"),
  CAMeRa = "NC_CAMeRa_EUR_to_AFR_CSA_clump10000kb.tsv",
  CAMeRa_EUR_CSA = "NC_CAMeRa_EUR_CSA_FEMA_EUR_clump10000kb.tsv",
  CAMeRa_EUR_AFR = "NC_CAMeRa_EUR_AFR_FEMA_EUR_clump10000kb.tsv",
  MRBEE = "NC_MRBEE_TL_EUR_to_AFR_CSA_clump10000kb.tsv",
  TEMR = "NC_TEMR_EUR_to_AFR_CSA_clump10000kb.tsv",
  MRAPSS_EUR = "NC_MRAPSS_EUR_clump10000kb.tsv"
)

paths <- setNames(
  file.path(input_dir, unname(required_files)),
  names(required_files)
)
missing_files <- paths[!file.exists(paths)]
if (length(missing_files) > 0L) {
  stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"))
}

check_columns <- function(x, columns, file_label) {
  missing_columns <- setdiff(columns, names(x))
  if (length(missing_columns) > 0L) {
    stop(
      file_label,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

normalise_threshold <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  fifelse(
    is.finite(x_num) & abs(log10(x_num) - log10(5e-8)) < 1e-6,
    "5e-08",
    fifelse(
      is.finite(x_num) & abs(log10(x_num) - log10(5e-5)) < 1e-6,
      "5e-05",
      NA_character_
    )
  )
}

standardise_result <- function(
    path,
    method_label,
    target_col,
    threshold_col,
    p_col,
    method_filter = NULL) {
  x <- fread(path, na.strings = c("", "NA", "NaN", "Inf", "-Inf"))
  x[, input_row := seq_len(.N)]
  check_columns(
    x,
    c(target_col, threshold_col, p_col, "status"),
    basename(path)
  )

  if (!is.null(method_filter)) {
    check_columns(x, "method", basename(path))
    x <- x[method == method_filter]
  }

  data.table(
    target = toupper(trimws(as.character(x[[target_col]]))),
    threshold = normalise_threshold(x[[threshold_col]]),
    method = method_label,
    p_value = suppressWarnings(as.numeric(x[[p_col]])),
    status = tolower(trimws(as.character(x[["status"]]))),
    source_file = normalizePath(path),
    source_row = x$input_row,
    exposure = as.character(x$exposure),
    outcome = as.character(x$outcome)
  )
}

camera_eur_csa <- standardise_result(
  paths[["CAMeRa_EUR_CSA"]],
  eur_method_levels[[1L]],
  "target",
  "threshold",
  "p"
)
camera_eur_afr <- standardise_result(
  paths[["CAMeRa_EUR_AFR"]],
  eur_method_levels[[2L]],
  "target",
  "threshold",
  "p"
)

mrbee_eur_base <- standardise_result(
  paths[["MRBEE"]],
  "MRBEE-EUR",
  "target",
  "threshold",
  "p",
  method_filter = "MRBEE_source_specific"
)
mrbee_eur_csa <- copy(mrbee_eur_base[target == "CSA"])
mrbee_eur_csa[, `:=`(
  target = "EUR",
  method = eur_method_levels[[3L]]
)]
mrbee_eur_afr <- copy(mrbee_eur_base[target == "AFR"])
mrbee_eur_afr[, `:=`(
  target = "EUR",
  method = eur_method_levels[[4L]]
)]

mrapss_eur <- standardise_result(
  paths[["MRAPSS_EUR"]],
  eur_method_levels[[5L]],
  "target_dataset",
  "threshold",
  "beta_pvalue"
)

all_results <- rbindlist(
  list(
    standardise_result(
      paths[["XMR"]],
      "XMR",
      "target",
      "threshold",
      "p"
    ),
    standardise_result(
      paths[["CAMeRa"]],
      "CAMERA",
      "target",
      "threshold",
      "p"
    ),
    standardise_result(
      paths[["MRBEE"]],
      "MRBEE",
      "target",
      "threshold",
      "p",
      method_filter = "MRBEE_target_specific"
    ),
    standardise_result(
      paths[["MRBEE"]],
      "MRBEE-TL",
      "target",
      "threshold",
      "p",
      method_filter = "MRBEE_TL_UV"
    ),
    standardise_result(
      paths[["TEMR"]],
      "TEMR",
      "target_dataset",
      "threshold",
      "beta_pvalue"
    ),
    camera_eur_csa,
    camera_eur_afr,
    mrbee_eur_csa,
    mrbee_eur_afr,
    mrapss_eur
  ),
  use.names = TRUE
)

all_results <- all_results[
  threshold %chin% threshold_levels &
    (
      (target %chin% c("CSA", "AFR") & method %chin% main_method_levels) |
        (target == "EUR" & method %chin% eur_method_levels)
    )
]
all_results[, `:=`(
  method = factor(method, levels = all_method_levels),
  threshold = factor(threshold, levels = threshold_levels),
  valid_p = status == "ok" & is.finite(p_value) &
    p_value > 0 & p_value <= 1
)]

expected_combinations <- rbindlist(
  list(
    CJ(
      target = c("CSA", "AFR"),
      threshold = threshold_levels,
      method = main_method_levels,
      unique = TRUE
    ),
    CJ(
      target = "EUR",
      threshold = threshold_levels,
      method = eur_method_levels,
      unique = TRUE
    )
  ),
  use.names = TRUE
)
observed_combinations <- unique(
  all_results[, .(
    target,
    threshold = as.character(threshold),
    method = as.character(method)
  )]
)
missing_combinations <- expected_combinations[
  !observed_combinations,
  on = names(expected_combinations)
]
if (nrow(missing_combinations) > 0L) {
  stop(
    "Missing target/threshold/method combination(s):\n",
    paste(capture.output(print(missing_combinations)), collapse = "\n")
  )
}

summary_table <- all_results[, .(
  n_total = .N,
  n_valid_p = sum(valid_p, na.rm = TRUE),
  n_failed_or_invalid = sum(!valid_p, na.rm = TRUE),
  n_zero_p = sum(status == "ok" & p_value == 0, na.rm = TRUE),
  n_nominal_significant = sum(valid_p & p_value < 0.05, na.rm = TRUE),
  min_p = if (any(valid_p)) min(p_value[valid_p]) else NA_real_,
  max_raw_minus_log10_p = if (any(valid_p)) {
    max(-log10(pmax(p_value[valid_p], .Machine$double.xmin)))
  } else {
    NA_real_
  },
  n_capped_at_5 = if (any(valid_p)) {
    sum(-log10(pmax(p_value[valid_p], .Machine$double.xmin)) > 5)
  } else {
    0L
  }
), by = .(target, threshold, method)]

qq_data <- all_results[valid_p == TRUE]
setorder(qq_data, target, method, threshold, p_value)
qq_data[, `:=`(
  n_valid = .N,
  rank = seq_len(.N),
  expected = -log10(seq_len(.N) / (.N + 1)),
  observed_raw = -log10(pmax(p_value, .Machine$double.xmin))
), by = .(target, method, threshold)]
qq_data[, `:=`(
  observed_plot = pmin(observed_raw, 5),
  capped_at_5 = observed_raw > 5
)]

make_null_envelope <- function(panel_data) {
  reference_n <- unique(panel_data$n_valid)
  if (length(reference_n) != 1L) {
    stop(
      "Thresholds within a panel have different valid-test counts: ",
      paste(sort(reference_n), collapse = ", ")
    )
  }
  rank_ref <- seq_len(reference_n)
  envelope <- data.table(
    expected = -log10(rank_ref / (reference_n + 1)),
    lower = -log10(qbeta(0.975, rank_ref, reference_n + 1 - rank_ref)),
    upper = -log10(qbeta(0.025, rank_ref, reference_n + 1 - rank_ref))
  )
  setorder(envelope, expected)
  envelope
}

make_panel <- function(
    method_value,
    target_value,
    show_x,
    show_y,
    title_value = method_value,
    subtitle_value = NULL) {
  panel_data <- qq_data[
    method == method_value & target == target_value
  ]
  if (nrow(panel_data) == 0L) {
    stop("No valid P values for ", method_value, " / ", target_value)
  }
  null_envelope <- make_null_envelope(panel_data)

  ggplot(
    panel_data,
    aes(
      x = expected,
      y = observed_plot,
      colour = threshold,
      shape = threshold,
      fill = threshold
    )
  ) +
    geom_ribbon(
      data = null_envelope,
      aes(x = expected, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      fill = "#D9D9D9",
      alpha = 0.55
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      colour = "#555555",
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    geom_point(size = 1.8, stroke = 0.5) +
    scale_colour_manual(
      values = threshold_colours,
      guide = "none"
    ) +
    scale_fill_manual(
      values = threshold_colours,
      guide = "none"
    ) +
    scale_shape_manual(
      values = threshold_shapes,
      breaks = threshold_levels,
      labels = threshold_labels,
      drop = FALSE,
      name = "IV threshold"
    ) +
    scale_x_continuous(
      limits = c(0, 2.5),
      breaks = seq(0, 2.5, by = 0.5),
      labels = c("0", "0.5", "1", "1.5", "2", "2.5"),
      expand = expansion(mult = c(0, 0.015))
    ) +
    scale_y_continuous(
      limits = c(0, 5),
      breaks = 0:5,
      expand = expansion(mult = c(0, 0.018))
    ) +
    labs(
      title = title_value,
      subtitle = subtitle_value,
      x = if (show_x) expression(Expected~~-log[10](italic(P))) else NULL,
      y = if (show_y) expression(Observed~~-log[10](italic(P))) else NULL
    ) +
    coord_cartesian(clip = "on") +
    theme_classic(base_size = 11.5, base_family = "sans") +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12.5,
        margin = margin(b = if (is.null(subtitle_value)) 6 else 2)
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 9.5,
        colour = "#444444",
        margin = margin(b = 5)
      ),
      axis.title = element_text(size = 11.2),
      axis.text = element_text(size = 9.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 11.2, face = "bold"),
      legend.text = element_text(size = 11),
      legend.key.width = unit(1.05, "lines"),
      legend.spacing.x = unit(0.6, "lines"),
      plot.margin = margin(5, 6, 5, 6)
    ) +
    guides(
      shape = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          size = 2.6,
          colour = unname(threshold_colours[threshold_levels]),
          fill = unname(threshold_colours[threshold_levels]),
          alpha = 1
        )
      )
    )
}

top_panels <- lapply(seq_along(main_method_levels), function(i) {
  make_panel(
    main_method_levels[[i]],
    "CSA",
    show_x = FALSE,
    show_y = i == 1L,
    subtitle_value = if (main_method_levels[[i]] %in% c("CAMERA", "MRBEE")) {
      "CSA + EUR FEMA IV"
    } else NULL
  ) + theme(legend.position = "none")
})

middle_panels <- lapply(seq_along(main_method_levels), function(i) {
  make_panel(
    main_method_levels[[i]],
    "AFR",
    show_x = FALSE,
    show_y = i == 1L,
    subtitle_value = if (main_method_levels[[i]] %in% c("CAMERA", "MRBEE")) {
      "AFR + EUR FEMA IV"
    } else NULL
  ) + theme(legend.position = "none")
})

eur_panels <- list(
  make_panel(
    eur_method_levels[[1L]],
    "EUR",
    show_x = TRUE,
    show_y = TRUE,
    title_value = "CAMERA",
    subtitle_value = "CSA + EUR FEMA IV"
  ) + theme(legend.position = "none"),
  make_panel(
    eur_method_levels[[2L]],
    "EUR",
    show_x = TRUE,
    show_y = FALSE,
    title_value = "CAMERA",
    subtitle_value = "AFR + EUR FEMA IV"
  ) + theme(legend.position = "none"),
  make_panel(
    eur_method_levels[[3L]],
    "EUR",
    show_x = TRUE,
    show_y = FALSE,
    title_value = "MRBEE",
    subtitle_value = "CSA + EUR FEMA IV"
  ) + theme(legend.position = "none"),
  make_panel(
    eur_method_levels[[4L]],
    "EUR",
    show_x = TRUE,
    show_y = FALSE,
    title_value = "MRBEE",
    subtitle_value = "AFR + EUR FEMA IV"
  ) + theme(legend.position = "none"),
  make_panel(
    eur_method_levels[[5L]],
    "EUR",
    show_x = TRUE,
    show_y = FALSE,
    title_value = "MR-APSS"
  ) + theme(legend.position = "none")
)

legend_source <- make_panel("XMR", "CSA", FALSE, FALSE)
shared_legend <- get_legend(legend_source)

top_row <- plot_grid(
  plotlist = top_panels,
  nrow = 1,
  align = "hv",
  axis = "tblr"
)
middle_row <- plot_grid(
  plotlist = middle_panels,
  nrow = 1,
  align = "hv",
  axis = "tblr"
)
eur_row <- plot_grid(
  plotlist = eur_panels,
  nrow = 1,
  align = "hv",
  axis = "tblr"
)

row_header <- function(panel_label, population_label) {
  ggplot() +
    annotate(
      "text",
      x = 0.008,
      y = 0.5,
      label = panel_label,
      hjust = 0,
      vjust = 0.5,
      fontface = "bold",
      size = 5.6,
      colour = "black"
    ) +
    annotate(
      "label",
      x = 0.5,
      y = 0.5,
      label = population_label,
      hjust = 0.5,
      vjust = 0.5,
      fontface = "bold",
      size = 6.2,
      colour = "white",
      fill = "#365A73",
      linewidth = 0.45,
      label.padding = grid::unit(0.12, "lines"),
      label.r = grid::unit(0.08, "lines")
    ) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 6, 0, 6))
}

combined_plot <- plot_grid(
  row_header("A", "CSA"),
  top_row,
  row_header("B", "AFR"),
  middle_row,
  row_header("C", "EUR"),
  eur_row,
  shared_legend,
  ncol = 1,
  rel_heights = c(0.11, 1, 0.11, 1, 0.11, 1, 0.13)
)

file_stem <- "Negative_control_QQ_by_method_CSA_AFR_EUR_3x5"
png_path <- file.path(output_dir, paste0(file_stem, ".png"))
pdf_path <- file.path(output_dir, paste0(file_stem, ".pdf"))
summary_path <- file.path(output_dir, paste0(file_stem, "_summary.tsv"))
plot_data_path <- file.path(output_dir, paste0(file_stem, "_plot_data.tsv"))

ggsave(
  filename = png_path,
  plot = combined_plot,
  width = 18,
  height = 12.2,
  units = "in",
  dpi = 400,
  bg = "white"
)
ggsave(
  filename = pdf_path,
  plot = combined_plot,
  width = 18,
  height = 12.2,
  units = "in",
  device = grDevices::pdf,
  useDingbats = FALSE,
  bg = "white"
)

summary_out <- copy(summary_table)
summary_out[, `:=`(
  threshold = as.character(threshold),
  method = as.character(method)
)]
setorder(summary_out, target, method, threshold)
fwrite(summary_out, summary_path, sep = "\t", na = "NA")

plot_data_out <- qq_data[, .(
  target,
  method = as.character(method),
  threshold = as.character(threshold),
  rank,
  n_valid,
  p_value,
  expected_minus_log10_p = expected,
  observed_minus_log10_p_raw = observed_raw,
  observed_minus_log10_p_plotted = observed_plot,
  capped_at_5,
  source_file,
  source_row
)]
setorder(plot_data_out, target, method, threshold, rank)
fwrite(plot_data_out, plot_data_path, sep = "\t", na = "NA")
stopifnot(nrow(summary_table) == 30L,
          all(summary_table$n_total == 70L),
          all(summary_table[method == "XMR", n_valid_p] == 70L),
          all(summary_table[method == "CAMERA", n_valid_p] == 70L))
fwrite(all_results, file.path(output_dir, paste0(file_stem, "_all_data.tsv")), sep = "\t", na = "NA")
fwrite(summary_out, file.path(output_dir, paste0(file_stem, "_checks.tsv")), sep = "\t", na = "NA")
fwrite(data.table(
  input = names(paths), path = normalizePath(unname(paths)),
  md5 = unname(tools::md5sum(paths)),
  rows = vapply(paths, function(p) nrow(fread(p)), integer(1L))
), file.path(output_dir, paste0(file_stem, "_inputs.tsv")), sep = "\t")
capture.output(sessionInfo(), file = file.path(output_dir, paste0(file_stem, "_sessionInfo.txt")))

message("XMR source (no effect heterogeneity QC): ", paths[["XMR"]])
message("Saved PNG: ", png_path)
message("Saved PDF: ", pdf_path)
message("Saved summary: ", summary_path)
message("Saved plot data: ", plot_data_path)
