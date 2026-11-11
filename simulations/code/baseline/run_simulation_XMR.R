#!/usr/bin/env Rscript

# XMR simulation workflow:
#   Step 1: select EUR hits; rank by FEMA, EUR, or EAS exposure P.
#   Step 2: clump all jobs together with EAS LD.
#   Step 3: save each selected IV's true class and report its composition.
#   Step 4: read Step 6 background parameters, fit, append, and report.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else {
  getwd()
}
source(file.path(SCRIPT_DIR, "simulation_method_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
xmr_priority <- toupper(get_arg(args, "--xmr-priority", "FEMA_P"))
if (!xmr_priority %in% c("FEMA_P", "EUR_P", "EAS_P", "BOTH", "ALL")) {
  stop("--xmr-priority must be FEMA_P, EUR_P, EAS_P, both, or all")
}
if (xmr_priority %in% c("BOTH", "ALL")) {
  failed <- FALSE
  priorities <- if (xmr_priority == "ALL") c("FEMA_P", "EUR_P", "EAS_P") else c("FEMA_P", "EUR_P")
  for (priority in priorities) {
    message("Running XMR clump priority: ", priority)
    child_args <- args
    inline <- grep("^--xmr-priority=", child_args)
    separate <- which(child_args == "--xmr-priority")
    if (length(inline)) child_args[inline] <- paste0("--xmr-priority=", priority)
    if (length(separate)) child_args[separate + 1L] <- priority
    status <- system2(
      file.path(R.home("bin"), "Rscript"),
      shQuote(c("--vanilla", file.path(SCRIPT_DIR, "run_simulation_XMR.R"), child_args))
    )
    failed <- failed || status != 0L
  }
  quit(save = "no", status = if (failed) 1L else 0L)
}
cfg <- parse_method_config(args, "XMR", SCRIPT_DIR)
selection_only <- has_flag(args, "--selection-only")
require_method_packages(if (selection_only) "data.table" else c("data.table", "XMR"))

het_q_cutoff <- as.numeric(get_arg(args, "--het-q-cutoff", "0.05"))
if (!is.finite(het_q_cutoff) || het_q_cutoff < 0 || het_q_cutoff > 1) {
  stop("--het-q-cutoff must be in [0, 1]")
}
if (has_flag(args, "--fema-qc") && has_flag(args, "--no-fema-qc")) {
  stop("Choose either --fema-qc or --no-fema-qc")
}
skip_fema_qc <- !has_flag(args, "--fema-qc")
result_method <- switch(xmr_priority, FEMA_P = "XMR_FEMA", EUR_P = "XMR_EUR_P", EAS_P = "XMR_EAS_P")
selection_mode <- paste0(xmr_priority, "_clump", if (skip_fema_qc) "_noFEMAQC" else "")
clump_label <- paste0(result_method, if (skip_fema_qc) "_noFEMAQC" else "")
message(
  "XMR FEMA heterogeneity QC: ",
  if (skip_fema_qc) "disabled" else paste0("enabled (BH q >= ", het_q_cutoff, ")")
)

files <- discover_summary_files(cfg)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
if (!selection_only) {
  result_file <- file.path(cfg$output_dir, "simulation_XMR_results.tsv")
  if (file.exists(result_file) && cfg$overwrite) {
    old <- read_tsv(result_file)
    replace <- old$method == result_method & old$selection_mode == selection_mode &
      old$replicate %in% files$replicate & old$threshold %in% cfg$thresholds
    old <- old[!replace]
    if (nrow(old)) {
      data.table::fwrite(old, result_file, sep = "\t", quote = FALSE, na = "NA")
    } else unlink(result_file)
  }
  completed_keys <- completed_result_keys(result_file, cfg$seed)
}
class_map <- c(
  shared = "shared", target_only = "target_only_EAS", aux_only = "aux_only_EUR",
  background = "background_only", null_effect = "null_effect"
)
class_origin <- c("foreground", "foreground", "foreground", "background_only", "null_effect")

message(
  "XMR setting=", cfg$setting_id, "; replicates=", nrow(files),
  "; thresholds=", paste(cfg$thresholds, collapse = ","),
  "; priority=", xmr_priority, "; LD=EAS; kb=", cfg$clump_kb, "; r2=", cfg$clump_r2
)

# Step 1: candidate selection and FEMA exposure effects.
analyses <- list()
jobs <- list()
candidates_by_job <- list()
summary_class_counts <- list()

for (file_index in seq_len(nrow(files))) {
  replicate <- files$replicate[file_index]
  summary_path <- files$path[file_index]
  message(sprintf(
    "Preparing replicate %03d/%03d", file_index, nrow(files)
  ))
  dat <- tryCatch(
    {
      dat <- read_simulation_summary(summary_path, cfg$expected_snps)
      if (!all(c("EFFECT_CLASS", "HAS_BACKGROUND") %in% names(dat))) {
        stop("Summary data lack EFFECT_CLASS/HAS_BACKGROUND truth annotations")
      }
      if (anyNA(dat$EFFECT_CLASS) || any(!dat$EFFECT_CLASS %in% names(class_map))) {
        stop("Unknown EFFECT_CLASS in summary data")
      }
      foreground <- dat$EFFECT_CLASS %in% c("shared", "target_only", "aux_only")
      if (anyNA(dat$HAS_BACKGROUND) ||
          any(dat$HAS_BACKGROUND != as.integer(dat$EFFECT_CLASS != "null_effect"))) {
        stop("HAS_BACKGROUND disagrees with EFFECT_CLASS")
      }
      if ("Z_FOREGROUND" %in% names(dat) &&
          (anyNA(dat$Z_FOREGROUND) || any(dat$Z_FOREGROUND != as.integer(foreground)))) {
        stop("Z_FOREGROUND disagrees with EFFECT_CLASS")
      }
      dat[, `:=`(
        IV_CLASS = unname(class_map[as.character(EFFECT_CLASS)]),
        IV_ORIGIN = class_origin[match(EFFECT_CLASS, names(class_map))]
      )]
      dat
    },
    error = function(condition) condition
  )

  family <- NULL
  if (!inherits(dat, "error")) {
    summary_class_counts[[as.character(replicate)]] <- tabulate(
      match(dat$EFFECT_CLASS, names(class_map)), length(class_map)
    )
    # Official XMR QC family: EUR exposure P < 5e-5.
    family <- dat[pval.exp.pop1 < 5e-5]
    family[, het_q := NA_real_]
    if (!skip_fema_qc && nrow(family)) {
      q_stat <- (family$b.exp.pop1 - family$b.exp.pop2)^2 /
        (family$se.exp.pop1^2 + family$se.exp.pop2^2)
      family[, het_q := stats::p.adjust(
        stats::pchisq(q_stat, df = 1L, lower.tail = FALSE),
        method = "BH"
      )]
    }
  }

  for (threshold in cfg$thresholds) {
    job_id <- make_job_id(
      clump_label, replicate, threshold, xmr_priority, "EAS"
    )
    xmr_threshold <- threshold
    candidates <- data.table::data.table()
    n_before <- NA_integer_
    status <- "failed"
    error_message <- if (inherits(dat, "error")) {
      conditionMessage(dat)
    } else {
      NA_character_
    }

    if (!inherits(dat, "error")) {
      before_qc <- family[pval.exp.pop1 < threshold]
      candidates <- if (skip_fema_qc) {
        data.table::copy(before_qc)
      } else {
        before_qc[is.finite(het_q) & het_q >= het_q_cutoff]
      }
      n_before <- nrow(before_qc)
      if (nrow(candidates)) {
        w_eur <- 1 / candidates$se.exp.pop1^2
        w_eas <- 1 / candidates$se.exp.pop2^2
        candidates[, b.exp.fema2 :=
          (w_eur * b.exp.pop1 + w_eas * b.exp.pop2) / (w_eur + w_eas)]
        candidates[, se.exp.fema2 := sqrt(1 / (w_eur + w_eas))]
        candidates[, pval.exp.fema2 := pmax(
          2 * stats::pnorm(
            abs(b.exp.fema2 / se.exp.fema2), lower.tail = FALSE
          ),
          .Machine$double.xmin
        )]
        status <- "ready"
        jobs[[length(jobs) + 1L]] <- data.table::data.table(
          job_id = job_id,
          population = "EAS",
          SNP = candidates$SNP,
          priority = switch(xmr_priority, FEMA_P = candidates$pval.exp.fema2,
                            EUR_P = candidates$pval.exp.pop1, EAS_P = candidates$pval.exp.pop2)
        )
        candidates_by_job[[job_id]] <- candidates
      } else {
        error_message <- if (skip_fema_qc) {
          "No SNP passed the EUR exposure P-value threshold"
        } else {
          "No SNP remained after XMR heterogeneity QC"
        }
      }
    }

    analyses[[length(analyses) + 1L]] <- data.table::data.table(
      replicate = replicate,
      summary_path = summary_path,
      threshold = threshold,
      selection_mode = selection_mode,
      job_id = job_id,
      model_threshold = xmr_threshold,
      n_qc_family = if (inherits(dat, "error")) NA_integer_ else nrow(family),
      n_candidates_before_qc = n_before,
      n_removed_by_qc = n_before - nrow(candidates),
      n_candidates_after_qc = if (inherits(dat, "error")) {
        NA_integer_
      } else {
        nrow(candidates)
      },
      fema_qc_applied = !skip_fema_qc,
      het_q_cutoff = if (skip_fema_qc) NA_real_ else het_q_cutoff,
      preparation_status = status,
      preparation_error = error_message
    )
  }
}

analyses <- data.table::rbindlist(analyses, fill = TRUE)
data.table::fwrite(
  analyses,
  file.path(
    cfg$output_dir,
    paste0("XMR_selection_", if (xmr_priority != "FEMA_P") paste0(xmr_priority, "_") else "",
           if (skip_fema_qc) "noFEMAQC.tsv" else "qc.tsv")
  ),
  sep = "\t", quote = FALSE, na = "NA"
)

# Step 2: one batched EAS clumping call.
selected <- data.table::data.table()
clump_error <- NA_character_
if (length(jobs)) {
  job_table <- data.table::rbindlist(jobs)
  message(
    "XMR batched EAS clumping: ",
    data.table::uniqueN(job_table$job_id), " jobs, ", nrow(job_table),
    " candidate rows"
  )
  selected_attempt <- tryCatch(
    run_vilma_clump(job_table, cfg, clump_label),
    error = function(condition) condition
  )
  if (inherits(selected_attempt, "error")) {
    clump_error <- conditionMessage(selected_attempt)
  } else {
    selected <- selected_attempt
  }
}

# Step 3: report true direct-effect classes, including for cached/completed fits.
detail_path <- file.path(cfg$output_dir, paste0(clump_label, "_IV_details.tsv"))
composition_path <- file.path(cfg$output_dir, paste0(clump_label, "_IV_composition.tsv"))
for (path in c(detail_path, composition_path)) {
  if (file.exists(path)) {
    old <- read_tsv(path)
    old <- old[!job_id %in% analyses$job_id]
    data.table::fwrite(old, path, sep = "\t", quote = FALSE, na = "NA")
  }
}
iv_counts <- list()
selection_failures <- 0L
total_runs <- nrow(analyses)
for (run_index in seq_len(total_runs)) {
  info <- analyses[run_index]
  candidates <- candidates_by_job[[info$job_id]]
  instruments <- data.table::data.table()
  counts_in_summary <- summary_class_counts[[as.character(info$replicate)]]
  selection_status <- "failed"
  selection_error <- info$preparation_error
  if (info$preparation_status == "ready" && is.na(clump_error)) {
    chosen <- selected[job_id == info$job_id][order(selection_rank)]
    matched <- match(chosen$SNP, candidates$SNP)
    if (anyNA(matched) || anyDuplicated(chosen$SNP) || anyNA(chosen$population) ||
        any(chosen$population != "EAS") ||
        !identical(as.integer(chosen$selection_rank), seq_len(nrow(chosen)))) {
      stop("Invalid selected IVs for ", info$job_id)
    }
    instruments <- data.table::copy(candidates[matched])
    instruments[, `:=`(
      setting_id = cfg$setting_id, replicate = info$replicate,
      threshold = info$threshold, selection_mode = selection_mode,
      job_id = info$job_id, selection_rank = chosen$selection_rank,
      clump_priority = chosen$priority, clump_panel = "EAS",
      clump_kb = cfg$clump_kb, clump_r2 = cfg$clump_r2
    )]
    data.table::setcolorder(instruments, c(
      "setting_id", "replicate", "threshold", "selection_mode", "job_id",
      "selection_rank", "SNP", "IV_ORIGIN", "IV_CLASS", "EFFECT_CLASS", "HAS_BACKGROUND"
    ))
    data.table::fwrite(
      instruments, detail_path, sep = "\t", quote = FALSE, na = "NA",
      append = file.exists(detail_path), col.names = !file.exists(detail_path)
    )
    selection_status <- "ok"
    selection_error <- NA_character_
  } else if (!is.null(counts_in_summary) && info$n_candidates_after_qc == 0L) {
    selection_status <- "no_candidates"
  } else if (info$preparation_status == "ready") {
    selection_error <- clump_error
  }
  selection_failed <- selection_status == "failed"
  selection_failures <- selection_failures + as.integer(selection_failed)
  selected_counts <- if (selection_failed) rep(NA_integer_, length(class_map)) else {
    tabulate(match(instruments$EFFECT_CLASS, names(class_map)), length(class_map))
  }
  n_selected <- if (selection_failed) NA_integer_ else nrow(instruments)
  iv_counts[[info$job_id]] <- data.table::data.table(
    n_foreground_IV = sum(selected_counts[1:3]),
    n_shared_IV = selected_counts[1L],
    n_target_only_IV = selected_counts[2L],
    n_aux_only_IV = selected_counts[3L],
    n_background_only_IV = selected_counts[4L],
    n_null_effect_IV = selected_counts[5L]
  )
  composition <- data.table::data.table(
    setting_id = cfg$setting_id, replicate = info$replicate,
    threshold = info$threshold, selection_mode = selection_mode, job_id = info$job_id,
    IV_ORIGIN = class_origin, IV_CLASS = unname(class_map),
    n_in_summary = counts_in_summary %||% rep(NA_integer_, length(class_map)),
    n_candidates = if (is.null(counts_in_summary)) rep(NA_integer_, length(class_map)) else {
      tabulate(match(candidates$EFFECT_CLASS, names(class_map)), length(class_map))
    },
    n_selected_IV = selected_counts, n_total_IV = n_selected,
    fraction_selected = if (!is.na(n_selected) && n_selected > 0L) {
      selected_counts / n_selected
    } else NA_real_,
    clump_panel = "EAS", clump_kb = cfg$clump_kb, clump_r2 = cfg$clump_r2,
    status = selection_status, error = selection_error, summary_path = info$summary_path
  )
  data.table::fwrite(
    composition, composition_path, sep = "\t", quote = FALSE, na = "NA",
    append = file.exists(composition_path), col.names = !file.exists(composition_path)
  )
  if (selection_only) {
    message(sprintf(
      paste0("[%d/%d] XMR IV selection rep=%03d P=%s IV=%s ",
             "foreground=%s (shared=%s target_only=%s aux_only=%s) ",
             "background_only=%s null_effect=%s status=%s"),
      run_index, total_runs, info$replicate, threshold_label(info$threshold), n_selected,
      sum(selected_counts[1:3]), selected_counts[1L], selected_counts[2L],
      selected_counts[3L], selected_counts[4L], selected_counts[5L], selection_status
    ))
    if (!is.na(selection_error)) message("  ", selection_error)
    flush.console()
  }
}
if (file.exists(detail_path)) message("Selected IV details: ", detail_path)
message("IV composition: ", composition_path)
if (selection_only) {
  message("Selection-only test: no XMR fit, no nvalid estimate, no background re-estimation.")
  quit(save = "no", status = if (selection_failures) 1L else 0L)
}

# Add origin counts to older result files without discarding completed estimates.
origin_columns <- names(iv_counts[[1L]])
existing_results <- if (file.exists(result_file)) read_tsv(result_file) else data.table::data.table()
if (nrow(existing_results)) {
  for (column in setdiff(origin_columns, names(existing_results))) {
    data.table::set(existing_results, j = column, value = NA_integer_)
  }
  for (column in origin_columns) {
    data.table::set(existing_results, j = column, value = as.integer(existing_results[[column]]))
  }
  for (run_index in seq_len(total_runs)) {
    info <- analyses[run_index]
    rows <- which(
      existing_results$setting_id == cfg$setting_id &
      existing_results$replicate == info$replicate & existing_results$threshold == info$threshold &
      existing_results$method == result_method & existing_results$population == "EAS" &
      existing_results$selection_mode == selection_mode
    )
    counts <- iv_counts[[info$job_id]]
    if (length(rows) && !is.na(counts$n_foreground_IV)) {
      current_n <- counts$n_foreground_IV + counts$n_background_only_IV + counts$n_null_effect_IV
      if (any(existing_results$status[rows] == "ok" &
              existing_results$n_clumped_IV[rows] != current_n, na.rm = TRUE)) {
        stop("Saved fit's IV count disagrees with current clump: ", info$job_id,
             "; inspect inputs before using --overwrite")
      }
      for (column in origin_columns) {
        data.table::set(existing_results, i = rows, j = column, value = counts[[column]])
      }
    }
  }
  data.table::fwrite(existing_results, result_file, sep = "\t", quote = FALSE, na = "NA")
}

# Step 4: fit, append, and report one result at a time.
for (run_index in seq_len(total_runs)) {
  info <- analyses[run_index]
  current_key <- result_key(
    info$replicate, result_method, info$threshold, selection_mode, "EAS"
  )
  if (current_key %in% completed_keys) {
    saved <- existing_results[
      setting_id == cfg$setting_id & replicate == info$replicate & threshold == info$threshold &
        method == result_method & population == "EAS" & selection_mode == info$selection_mode
    ]
    report_method_result(saved, run_index, total_runs)
    message("  Already completed; saved estimate retained.")
    next
  }

  status <- "failed"
  error_message <- info$preparation_error
  beta <- se <- p <- tau_sq <- pi0 <- NA_real_
  n_clumped <- 0L
  fit_niv <- 0L
  nvalid <- NA_integer_

  if (info$preparation_status == "ready" && is.na(clump_error)) {
    candidates <- candidates_by_job[[info$job_id]]
    iv_snps <- selected_snps(selected, info$job_id)
    matched <- match(iv_snps, candidates$SNP)
    instruments <- candidates[matched[!is.na(matched)]]
    n_clumped <- nrow(instruments)
    fit_niv <- n_clumped

    if (n_clumped <= 3L) {
      error_message <- "XMR requires more than three clumped IVs"
    } else {
      fit_attempt <- tryCatch({
        background <- load_background_parameters(cfg, info$replicate)
        tolerance <- if (n_clumped < 100L) 1e-8 else 1e-7
        fit <- XMR::fit_XMR(
          data = as.data.frame(instruments),
          C = background$C,
          Omega0 = background$Omega0,
          Threshold = info$model_threshold,
          tol1 = tolerance,
          tol2 = tolerance,
          min_thres = 1e-2,
          verbose = FALSE
        )
        beta_fit <- scalar_or_na(fit$beta)
        se_fit <- scalar_or_na(fit$beta.se)
        if (!is.finite(beta_fit) || !is.finite(se_fit) || se_fit <= 0) {
          stop("XMR returned a non-finite estimate")
        }
        list(
          beta = beta_fit,
          se = se_fit,
          p = scalar_or_na(fit$beta.pvalue),
          nIV = scalar_or_na(fit$nIV),
          nvalid = scalar_or_na(fit$nvalid),
          tau_sq = scalar_or_na(fit$tau.sq),
          pi0 = scalar_or_na(fit$pi0)
        )
      }, error = function(condition) condition)

      if (inherits(fit_attempt, "error")) {
        error_message <- conditionMessage(fit_attempt)
      } else {
        status <- "ok"
        error_message <- NA_character_
        beta <- fit_attempt$beta
        se <- fit_attempt$se
        p <- fit_attempt$p
        if (is.finite(fit_attempt$nIV)) fit_niv <- as.integer(fit_attempt$nIV)
        nvalid <- as.integer(fit_attempt$nvalid)
        tau_sq <- fit_attempt$tau_sq
        pi0 <- fit_attempt$pi0
      }
    }
  } else if (info$preparation_status == "ready") {
    error_message <- clump_error
  }

  result <- simulation_result(
    cfg, info$replicate, result_method, info$threshold,
    selection_mode, "EAS", beta = beta, se = se, p = p,
    n_iv = fit_niv, status = status, error = error_message
  )
  result[, `:=`(
    model_threshold = info$model_threshold,
    n_clumped_IV = n_clumped,
    nvalid = nvalid,
    tau_sq = tau_sq,
    pi0 = pi0,
    n_qc_family = info$n_qc_family,
    n_candidates_before_qc = info$n_candidates_before_qc,
    n_removed_by_qc = info$n_removed_by_qc,
    n_candidates_after_qc = info$n_candidates_after_qc,
    clump_kb = cfg$clump_kb,
    clump_r2 = cfg$clump_r2,
    clump_panel = "EAS"
  )]
  counts <- iv_counts[[info$job_id]]
  for (column in origin_columns) {
    data.table::set(result, j = column, value = counts[[column]])
  }
  write_checkpoint(list(result), result_file)
  report_method_result(result, run_index, total_runs)
}

final <- if (file.exists(result_file)) read_tsv(result_file) else data.table::data.table()
message(
  "XMR finished: ", sum(final$status == "ok"), "/", nrow(final),
  " successful rows; results: ", result_file
)
