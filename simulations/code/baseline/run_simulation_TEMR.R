#!/usr/bin/env Rscript

# TEMR simulation workflow:
#   Step 1: prepare target-or-auxiliary exposure candidates for every run.
#   Step 2: clump all jobs together so each ViLMA LD block is loaded once.
#   Step 3: fit each replicate/threshold, append it, and report it immediately.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else {
  getwd()
}
source(file.path(SCRIPT_DIR, "simulation_method_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
cfg <- parse_method_config(args, "TEMR", SCRIPT_DIR)
require_method_packages(c("data.table", "TEMR"))

files <- discover_summary_files(cfg)
result_file <- prepare_result_file(cfg, "simulation_TEMR_results.tsv")
completed_keys <- completed_result_keys(result_file, cfg$seed)
selection_mode <- "target_or_aux_P__combined_P__EAS_LD"

message(
  "TEMR setting=", cfg$setting_id, "; replicates=", nrow(files),
  "; thresholds=", paste(cfg$thresholds, collapse = ",")
)

# Step 1: candidate selection. Keep one compact row per planned fit.
analyses <- list()
jobs <- list()
candidates_by_job <- list()

for (file_index in seq_len(nrow(files))) {
  replicate <- files$replicate[file_index]
  summary_path <- files$path[file_index]
  message(sprintf(
    "Preparing replicate %03d/%03d", file_index, nrow(files)
  ))
  dat <- tryCatch(
    read_simulation_summary(summary_path, cfg$expected_snps),
    error = function(condition) condition
  )

  for (threshold in cfg$thresholds) {
    job_id <- make_job_id("TEMR", replicate, threshold, "combined_P", "EAS")
    candidates <- data.table::data.table()
    status <- "failed"
    error_message <- if (inherits(dat, "error")) {
      conditionMessage(dat)
    } else {
      NA_character_
    }

    if (!inherits(dat, "error")) {
      candidates <- dat[
        pval.exp.pop2 < threshold | pval.exp.pop1 < threshold
      ]
      if (nrow(candidates)) {
        combined_z <- stats::qnorm(
          pmax(candidates$pval.exp.pop2, .Machine$double.xmin) / 2,
          lower.tail = FALSE
        ) + stats::qnorm(
          pmax(candidates$pval.exp.pop1, .Machine$double.xmin) / 2,
          lower.tail = FALSE
        )
        candidates[, P.combined := pmax(
          2 * stats::pnorm(combined_z, lower.tail = FALSE),
          .Machine$double.xmin
        )]
        status <- "ready"
        jobs[[length(jobs) + 1L]] <- data.table::data.table(
          job_id = job_id,
          population = "EAS",
          SNP = candidates$SNP,
          priority = candidates$P.combined
        )
        candidates_by_job[[job_id]] <- candidates
      } else {
        error_message <- paste(
          "No target-or-auxiliary exposure SNP passed the threshold"
        )
      }
    }

    analyses[[length(analyses) + 1L]] <- data.table::data.table(
      replicate = replicate,
      threshold = threshold,
      job_id = job_id,
      n_candidates = if (inherits(dat, "error")) NA_integer_ else nrow(candidates),
      preparation_status = status,
      preparation_error = error_message
    )
  }
}

analyses <- data.table::rbindlist(analyses, fill = TRUE)
data.table::fwrite(
  analyses,
  file.path(cfg$output_dir, "TEMR_selection_qc.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

# Step 2: one batched clumping call for all replicate/threshold jobs.
selected <- data.table::data.table()
clump_error <- NA_character_
if (length(jobs)) {
  job_table <- data.table::rbindlist(jobs)
  message(
    "TEMR batched EAS clumping: ",
    data.table::uniqueN(job_table$job_id), " jobs, ", nrow(job_table),
    " candidate rows"
  )
  selected_attempt <- tryCatch(
    run_vilma_clump(job_table, cfg, "TEMR"),
    error = function(condition) condition
  )
  if (inherits(selected_attempt, "error")) {
    clump_error <- conditionMessage(selected_attempt)
  } else {
    selected <- selected_attempt
  }
}

# Step 3: fit, append, and report one result at a time.
total_runs <- nrow(analyses)
for (run_index in seq_len(total_runs)) {
  info <- analyses[run_index]
  current_key <- result_key(
    info$replicate, "TEMR", info$threshold, selection_mode, "EAS"
  )
  if (current_key %in% completed_keys) {
    message(sprintf(
      "[%d/%d] TEMR rep=%03d P=%s: already completed",
      run_index, total_runs, info$replicate,
      threshold_label(info$threshold)
    ))
    next
  }

  status <- "failed"
  error_message <- info$preparation_error
  beta <- se <- p <- statistic <- NA_real_
  n_clumped <- 0L

  if (info$preparation_status == "ready" && is.na(clump_error)) {
    candidates <- candidates_by_job[[info$job_id]]
    iv_snps <- selected_snps(selected, info$job_id)
    matched <- match(iv_snps, candidates$SNP)
    clumped <- candidates[matched[!is.na(matched)]]
    n_clumped <- nrow(clumped)

    if (n_clumped < 3L) {
      error_message <- "TEMR requires at least three clumped IVs"
    } else {
      fit_attempt <- tryCatch({
        # Official order: target EAS first, auxiliary EUR second.
        betaXG <- as.matrix(clumped[, .(b.exp.pop2, b.exp.pop1)])
        betaYG <- as.matrix(clumped[, .(b.out.pop2, b.out.pop1)])
        sebetaXG <- as.matrix(clumped[, .(se.exp.pop2, se.exp.pop1)])
        sebetaYG <- as.matrix(clumped[, .(se.out.pop2, se.out.pop1)])
        colnames(betaXG) <- colnames(betaYG) <- c("target1", "auxiliary")
        colnames(sebetaXG) <- colnames(sebetaYG) <-
          c("target1", "auxiliary")

        fit <- TEMR::TEMR(
          betaXG, betaYG, sebetaXG, sebetaYG,
          rho = NULL, meth = "CG"
        )
        target <- fit$resultall$target1
        if (is.null(target)) stop("TEMR output has no target1 result")
        beta_fit <- as.numeric(target$beta[["bbe1"]])
        statistic_fit <- as.numeric(target$stat)
        p_fit <- as.numeric(target$pvalue)
        chi_square <- -statistic_fit
        if (length(beta_fit) != 1L || length(statistic_fit) != 1L ||
            length(p_fit) != 1L || !is.finite(beta_fit) ||
            !is.finite(chi_square) || chi_square <= 0 ||
            !is.finite(p_fit) || p_fit < 0 || p_fit > 1) {
          stop("TEMR returned invalid target1 output")
        }
        # Equivalent 1-df SE: chi-square = (beta / SE)^2 and chi-square = -stat.
        se_fit <- abs(beta_fit) / sqrt(chi_square)
        if (!is.finite(se_fit) || se_fit <= 0) {
          stop("TEMR likelihood-ratio-equivalent SE is invalid")
        }
        list(
          beta = beta_fit,
          se = se_fit,
          p = max(p_fit, .Machine$double.xmin),
          statistic = statistic_fit
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
        statistic <- fit_attempt$statistic
      }
    }
  } else if (info$preparation_status == "ready") {
    error_message <- clump_error
  }

  result <- simulation_result(
    cfg, info$replicate, "TEMR", info$threshold,
    selection_mode, "EAS", beta = beta, se = se, p = p,
    n_iv = n_clumped, status = status, error = error_message
  )
  result[, `:=`(
    n_candidates = info$n_candidates,
    n_clumped_IV = n_clumped,
    auxiliary_population = "EUR",
    clump_kb = cfg$clump_kb,
    clump_r2 = cfg$clump_r2,
    clump_panel = "EAS",
    statistic = statistic
  )]
  write_checkpoint(list(result), result_file)
  report_method_result(result, run_index, total_runs)
}

final <- if (file.exists(result_file)) read_tsv(result_file) else data.table::data.table()
message(
  "TEMR finished: ", sum(final$status == "ok"), "/", nrow(final),
  " successful rows; results: ", result_file
)
