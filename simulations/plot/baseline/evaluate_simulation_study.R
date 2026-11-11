#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else {
  getwd()
}
source(file.path(SCRIPT_DIR, "..", "..", "code", "baseline", "simulation_method_utils.R"))

require_method_packages(c("data.table", "ggplot2"))
args <- commandArgs(trailingOnly = TRUE)
root <- get_arg(args, "--root", DEFAULT_WORK_ROOT)
settings_path <- get_arg(
  args, "--settings", file.path(root, "manifests/settings.tsv")
)
xmr_priority <- toupper(get_arg(args, "--xmr-priority", "FEMA_P"))
if (!xmr_priority %in% c("FEMA_P", "EUR_P", "EAS_P", "BOTH", "ALL")) {
  stop("--xmr-priority must be FEMA_P, EUR_P, EAS_P, both, or all")
}
xmr_methods <- switch(
  xmr_priority, FEMA_P = "XMR_FEMA", EUR_P = "XMR_EUR_P",
  EAS_P = "XMR_EAS_P", BOTH = c("XMR_FEMA", "XMR_EUR_P"),
  ALL = c("XMR_FEMA", "XMR_EUR_P", "XMR_EAS_P")
)
method_levels <- c(
  "XMR_FEMA", "XMR_EUR_P", "XMR_EAS_P", "CAMeRa_FEMA", "MRBEE_SOURCE", "MRBEE_TARGET",
  "MRBEE_TL", "TEMR", "MRAPSS_TARGET"
)
method_selection <- get_arg(args, "--methods", "all")
requested_methods <- if (method_selection == "all") {
  c(xmr_methods, setdiff(method_levels, c("XMR_FEMA", "XMR_EUR_P", "XMR_EAS_P")))
} else {
  unique(trimws(strsplit(method_selection, ",", fixed = TRUE)[[1L]]))
}
if (!length(requested_methods) || any(!requested_methods %in% method_levels)) {
  stop("--methods must be all or comma-separated estimator names: ",
       paste(method_levels, collapse = ", "))
}
if (has_flag(args, "--fema-qc") && has_flag(args, "--no-fema-qc")) {
  stop("Choose either --fema-qc or --no-fema-qc")
}
xmr_selection_modes <- setNames(
  paste0(c("FEMA_P", "EUR_P", "EAS_P"), "_clump",
         if (has_flag(args, "--fema-qc")) "" else "_noFEMAQC"),
  c("XMR_FEMA", "XMR_EUR_P", "XMR_EAS_P")
)
alpha <- as.numeric(get_arg(args, "--alpha", "0.05"))
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
  stop("--alpha must be in (0, 1)")
}

evaluation_dir <- file.path(root, "evaluation")
figure_dir <- file.path(root, "figures")
dir.create(evaluation_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

result_files <- Sys.glob(file.path(
  root, "method_results", "*", "*", "simulation_*_results.tsv"
))
if (!length(result_files)) {
  stop("No method result files found below ", file.path(root, "method_results"))
}
results <- data.table::rbindlist(
  lapply(result_files, data.table::fread, showProgress = FALSE),
  fill = TRUE,
  use.names = TRUE,
  idcol = "result_file_index"
)
results[, result_file := result_files[result_file_index]]
results[, result_file_index := NULL]

required <- c(
  "setting_id", "replicate", "method", "population", "threshold",
  "selection_mode", "beta_true", "rho", "n_aux", "n_target",
  "beta", "se", "p", "ci_lower", "ci_upper", "nIV", "status", "error"
)
missing <- setdiff(required, names(results))
if (length(missing)) {
  stop("Combined method results lack: ", paste(missing, collapse = ", "))
}
results <- results[method %in% requested_methods]
results <- results[!method %in% names(xmr_selection_modes) |
                     selection_mode == unname(xmr_selection_modes[method])]
if (!nrow(results)) {
  stop("No results found for requested methods: ",
       paste(requested_methods, collapse = ", "),
       "; XMR selection modes: ", paste(xmr_selection_modes, collapse = ", "))
}
message("Evaluating methods: ", paste(requested_methods, collapse = ", "))
if (any(requested_methods %in% names(xmr_selection_modes))) {
  message("XMR selection mode(s): ", paste(
    xmr_selection_modes[intersect(requested_methods, names(xmr_selection_modes))], collapse = ", "
  ))
}

key <- c(
  "setting_id", "replicate", "method", "population", "threshold",
  "selection_mode"
)
duplicates <- results[, .N, by = key][N > 1L]
if (nrow(duplicates)) {
  data.table::fwrite(
    duplicates,
    file.path(evaluation_dir, "duplicate_result_keys.tsv"),
    sep = "\t"
  )
  stop("Duplicated result keys found; see duplicate_result_keys.tsv")
}

unknown <- setdiff(unique(results$method), method_levels)
if (length(unknown)) {
  warning("Unexpected method labels: ", paste(unknown, collapse = ", "))
}
results[, method := factor(method, levels = c(method_levels, unknown))]
results[, threshold_label := threshold_label(threshold)]
results[, valid :=
  status == "ok" & is.finite(beta) & is.finite(se) & se > 0 &
  is.finite(p) & p >= 0 & p <= 1]

combined_path <- file.path(evaluation_dir, "combined_simulation_results.tsv.gz")
data.table::fwrite(results, combined_path, sep = "\t", quote = FALSE, na = "NA")

settings <- if (file.exists(settings_path)) {
  data.table::fread(settings_path)
} else {
  unique(results[, .(
    setting_id, beta_true, rho, n_aux, n_target,
    expected_replicates = max(replicate)
  )])
}
completion <- results[, .(
  observed_rows = .N,
  unique_replicates = data.table::uniqueN(replicate),
  successful_rows = sum(valid),
  failed_rows = sum(!valid)
), by = setting_id]
completion <- merge(
  settings[, .(
    setting_id,
    expected_replicates,
    n_thresholds = if ("n_thresholds" %in% names(settings)) {
      as.integer(n_thresholds)
    } else {
      2L
    }
  )],
  completion,
  by = "setting_id",
  all.x = TRUE
)
completion[, expected_methods := length(requested_methods)]
completion[, expected_rows := expected_replicates * expected_methods * n_thresholds]
completion[, complete := observed_rows == expected_rows]
data.table::fwrite(
  completion, file.path(evaluation_dir, "completion_qc.tsv"), sep = "\t"
)

wilson_interval <- function(successes, total, level = 0.95) {
  if (!total) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - level) / 2)
  proportion <- successes / total
  denominator <- 1 + z^2 / total
  centre <- (proportion + z^2 / (2 * total)) / denominator
  half_width <- z * sqrt(
    proportion * (1 - proportion) / total + z^2 / (4 * total^2)
  ) / denominator
  c(max(0, centre - half_width), min(1, centre + half_width))
}

summarize_group <- function(table, truth) {
  ok <- table[table$valid]
  rejected <- ok$p < alpha
  interval <- wilson_interval(sum(rejected), nrow(ok))
  list(
    n_total = nrow(table),
    n_ok = nrow(ok),
    success_rate = nrow(ok) / nrow(table),
    mean_beta = if (nrow(ok)) mean(ok$beta) else NA_real_,
    median_beta = if (nrow(ok)) stats::median(ok$beta) else NA_real_,
    bias = if (nrow(ok)) mean(ok$beta - truth) else NA_real_,
    empirical_sd = if (nrow(ok) > 1L) stats::sd(ok$beta) else NA_real_,
    mean_se = if (nrow(ok)) mean(ok$se) else NA_real_,
    rmse = if (nrow(ok)) {
      sqrt(mean((ok$beta - truth)^2))
    } else NA_real_,
    coverage_95 = if (nrow(ok)) {
      mean(ok$ci_lower <= truth & ok$ci_upper >= truth)
    } else NA_real_,
    rejection_rate = if (nrow(ok)) mean(rejected) else NA_real_,
    rejection_lower = interval[[1L]],
    rejection_upper = interval[[2L]],
    mean_nIV = if (nrow(ok)) mean(ok$nIV, na.rm = TRUE) else NA_real_,
    median_nIV = if (nrow(ok)) {
      as.numeric(stats::median(ok$nIV, na.rm = TRUE))
    } else NA_real_
  )
}

group_columns <- c(
  "setting_id", "beta_true", "rho", "n_aux", "n_target", "method", "population",
  "threshold", "threshold_label", "selection_mode"
)
performance <- results[, summarize_group(.SD, beta_true[[1L]]), by = group_columns]
performance[, estimand := ifelse(
  beta_true == 0, "type_I_error", "power"
)]
data.table::fwrite(
  performance,
  file.path(evaluation_dir, "method_performance.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)
data.table::fwrite(
  performance[beta_true == 0],
  file.path(evaluation_dir, "type_I_error.tsv"),
  sep = "\t"
)
data.table::fwrite(
  performance[beta_true != 0],
  file.path(evaluation_dir, "power.tsv"),
  sep = "\t"
)

plot_data <- results[valid == TRUE]
plot_root <- file.path(figure_dir, "by_setting")
dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)
setting_values <- unique(results[, .(
  setting_id, beta_true, rho, n_aux, n_target
)])

null_data <- plot_data[beta_true == 0]
qq <- NULL
if (nrow(null_data)) {
  qq <- null_data[order(p), {
    n <- .N
    rank <- seq_len(n)
    expected <- (rank - 0.5) / n
    lower_p <- stats::qbeta(0.025, rank, n - rank + 1)
    upper_p <- stats::qbeta(0.975, rank, n - rank + 1)
    list(
      expected_log10 = -log10(expected),
      observed_log10 = -log10(pmax(p, .Machine$double.xmin)),
      envelope_lower = -log10(upper_p),
      envelope_upper = -log10(pmax(lower_p, .Machine$double.xmin))
    )
  }, by = .(setting_id, beta_true, rho, n_aux, n_target, method, threshold_label)]
  data.table::fwrite(
    qq, file.path(evaluation_dir, "QQ_data.tsv.gz"), sep = "\t"
  )
}

for (i in seq_len(nrow(setting_values))) {
  setting <- setting_values[i]
  id <- as.character(setting$setting_id)
  truth <- as.numeric(setting$beta_true)
  rho_value <- as.numeric(setting$rho)
  n_aux_value <- as.integer(setting$n_aux)
  n_target_value <- as.integer(setting$n_target)
  setting_dir <- file.path(plot_root, id)
  dir.create(setting_dir, recursive = TRUE, showWarnings = FALSE)
  current <- plot_data[setting_id == id]
  title_suffix <- paste0(
    "beta=", truth, ", rho=", rho_value,
    ", N_aux=", n_aux_value, ", N_target=", n_target_value
  )

  if (nrow(current)) {
    effect_plot <- ggplot2::ggplot(
      current,
      ggplot2::aes(method, beta, fill = threshold_label)
    ) +
      ggplot2::geom_hline(yintercept = truth, linetype = "dashed") +
      ggplot2::geom_boxplot(
        position = ggplot2::position_dodge(width = 0.8),
        width = 0.68,
        outlier.shape = NA
      ) +
      ggplot2::coord_cartesian(ylim = c(-1, 1)) +
      ggplot2::labs(
        title = paste("Effect estimates:", title_suffix),
        x = NULL, y = "Estimated beta", fill = "IV P threshold"
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
        panel.grid.minor = ggplot2::element_blank()
      )
    ggplot2::ggsave(
      file.path(setting_dir, "effect_boxplot.png"),
      effect_plot, width = 12, height = 8, dpi = 300
    )
  }

  if (truth == 0) {
    current_qq <- if (is.null(qq)) {
      data.table::data.table()
    } else {
      qq[setting_id == id]
    }
    if (nrow(current_qq)) {
      qq_plot <- ggplot2::ggplot(
        current_qq,
        ggplot2::aes(expected_log10, observed_log10)
      ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = envelope_lower, ymax = envelope_upper),
          fill = "grey80"
        ) +
        ggplot2::geom_abline(linetype = "dashed") +
        ggplot2::geom_point(size = 1.2, colour = "#006D77") +
        ggplot2::facet_grid(
          rows = ggplot2::vars(method),
          cols = ggplot2::vars(threshold_label),
          scales = "free"
        ) +
        ggplot2::labs(
          title = paste("Null QQ:", title_suffix),
          x = expression(Expected~~-log[10](italic(P))),
          y = expression(Observed~~-log[10](italic(P)))
        ) +
        ggplot2::theme_bw(base_size = 9)
      ggplot2::ggsave(
        file.path(setting_dir, "QQ_plot.png"), qq_plot,
        width = 8.5, height = 15, dpi = 300, limitsize = FALSE
      )
    }

    type1 <- performance[
      setting_id == id & is.finite(rejection_rate)
    ]
    if (nrow(type1)) {
      type1_plot <- ggplot2::ggplot(
        type1,
        ggplot2::aes(method, rejection_rate, fill = threshold_label)
      ) +
        ggplot2::geom_hline(yintercept = alpha, linetype = "dashed") +
        ggplot2::geom_col(
          position = ggplot2::position_dodge(width = 0.8), width = 0.7
        ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = rejection_lower, ymax = rejection_upper),
          position = ggplot2::position_dodge(width = 0.8), width = 0.2
        ) +
        ggplot2::coord_cartesian(
          ylim = c(0, max(0.2, type1$rejection_upper, na.rm = TRUE))
        ) +
        ggplot2::labs(
          title = paste("Type I error:", title_suffix),
          x = NULL, y = "Type I error", fill = "IV P threshold"
        ) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
        )
      ggplot2::ggsave(
        file.path(setting_dir, "type_I_error.png"),
        type1_plot, width = 12, height = 8, dpi = 300
      )
    }
  } else {
    power <- performance[
      setting_id == id & is.finite(rejection_rate)
    ]
    if (nrow(power)) {
      power_plot <- ggplot2::ggplot(
        power,
        ggplot2::aes(method, rejection_rate, fill = threshold_label)
      ) +
        ggplot2::geom_col(
          position = ggplot2::position_dodge(width = 0.8), width = 0.7
        ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = rejection_lower, ymax = rejection_upper),
          position = ggplot2::position_dodge(width = 0.8), width = 0.2
        ) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(
          title = paste("Power:", title_suffix),
          x = NULL, y = "Power", fill = "IV P threshold"
        ) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
        )
      ggplot2::ggsave(
        file.path(setting_dir, "power.png"),
        power_plot, width = 12, height = 8, dpi = 300
      )
    }
  }
  message("Finished plots for setting: ", id)
}

message("Combined rows: ", nrow(results))
message("Combined results: ", combined_path)
message("Figures: ", figure_dir)
