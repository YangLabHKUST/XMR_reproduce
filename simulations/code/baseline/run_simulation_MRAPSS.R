#!/usr/bin/env Rscript

# MR-APSS target workflow:
#   Step 1: prepare EAS exposure candidates for every run.
#   Step 2: clump all jobs together with EAS LD.
#   Step 3: read target background parameters, fit, append, and report.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else {
  getwd()
}
source(file.path(SCRIPT_DIR, "simulation_method_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
cfg <- parse_method_config(args, "MRAPSS", SCRIPT_DIR)
cfg$refresh_stale_clump <- TRUE
require_method_packages(c("data.table", "MRAPSS"))

files <- discover_summary_files(cfg)
result_file <- prepare_result_file(cfg, "simulation_MRAPSS_results.tsv")
selection_mode <- "target_EAS_P_and_LD"
analysis_signatures <- stats::setNames(
  rep(NA_character_, nrow(files)), as.character(files$replicate)
)
resume_signatures <- analysis_signatures
background_cache <- list()
fingerprint_parameters <- list(
  setting_id = cfg$setting_id,
  master_seed = cfg$seed,
  expected_snps = cfg$expected_snps,
  beta_true = cfg$beta_true,
  rho = cfg$rho,
  n_aux = cfg$n_aux,
  n_target = cfg$n_target,
  clump_kb = cfg$clump_kb,
  clump_r2 = cfg$clump_r2,
  clump_population = "EAS",
  python = cfg$python,
  selection_bias_correction = TRUE,
  minimum_clumped_ivs = 5L
)

message(
  "MR-APSS setting=", cfg$setting_id, "; replicates=", nrow(files),
  "; thresholds=", paste(cfg$thresholds, collapse = ",")
)

# Step 1: candidate selection.
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
  background <- if (inherits(dat, "error")) {
    dat
  } else {
    tryCatch(
      load_background_parameters(
        cfg, replicate, require_crosspop_zero = FALSE,
        summary_path = summary_path
      ),
      error = function(condition) condition
    )
  }
  replicate_name <- as.character(replicate)
  analysis_signature <- analysis_fingerprint(
    summary_path,
    background_paths = c(
      target = file.path(
        cfg$background_dir, sprintf("rep%03d", replicate),
        "background_parameters.rds"
      ),
      common_snps = cfg$common_snps,
      eas_schema = cfg$eas_schema
    ),
    code_paths = c(
      file.path(SCRIPT_DIR, "simulation_method_utils.R"),
      file.path(SCRIPT_DIR, "run_simulation_MRAPSS.R"),
      cfg$clumper
    ),
    parameters = fingerprint_parameters,
    packages = "MRAPSS"
  )
  analysis_signatures[[replicate_name]] <- analysis_signature
  if (!inherits(background, "error")) {
    background_cache[[replicate_name]] <- background
    resume_signatures[[replicate_name]] <- analysis_signature
  }

  for (threshold in cfg$thresholds) {
    job_id <- make_job_id(
      "MRAPSS_TARGET", replicate, threshold, "target_EAS", "EAS"
    )
    candidates <- data.table::data.table()
    status <- "failed"
    error_message <- if (inherits(background, "error")) {
      conditionMessage(background)
    } else {
      NA_character_
    }

    if (!inherits(dat, "error")) {
      candidates <- dat[pval.exp.pop2 < threshold]
      if (nrow(candidates) && !inherits(background, "error")) {
        status <- "ready"
        jobs[[length(jobs) + 1L]] <- data.table::data.table(
          job_id = job_id,
          population = "EAS",
          SNP = candidates$SNP,
          priority = candidates$pval.exp.pop2
        )
        candidates_by_job[[job_id]] <- candidates
      } else if (!nrow(candidates) && !inherits(background, "error")) {
        error_message <- "No EAS exposure SNP passed the threshold"
      }
    }

    analyses[[length(analyses) + 1L]] <- data.table::data.table(
      replicate = replicate,
      analysis_signature = analysis_signature,
      summary_path = summary_path,
      threshold = threshold,
      population = "EAS",
      job_id = job_id,
      n_candidates = if (inherits(dat, "error")) NA_integer_ else nrow(candidates),
      preparation_status = status,
      preparation_error = error_message
    )
  }
}

analyses <- data.table::rbindlist(analyses, fill = TRUE)
completed_keys <- completed_result_keys(
  result_file, cfg$seed, analysis_signatures = resume_signatures
)
data.table::fwrite(
  analyses,
  file.path(cfg$output_dir, "MRAPSS_selection_qc.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

# Step 2: one batched EAS clumping call.
selected <- data.table::data.table()
clump_error <- NA_character_
if (length(jobs)) {
  job_table <- data.table::rbindlist(jobs)
  message(
    "MR-APSS batched EAS clumping: ",
    data.table::uniqueN(job_table$job_id), " jobs, ", nrow(job_table),
    " candidate rows"
  )
  selected_attempt <- tryCatch(
    run_vilma_clump(job_table, cfg, "MRAPSS"),
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
    info$replicate, "MRAPSS_TARGET", info$threshold,
    selection_mode, "EAS"
  )
  if (current_key %in% completed_keys) {
    message(sprintf(
      "[%d/%d] MRAPSS_TARGET rep=%03d P=%s: already completed",
      run_index, total_runs, info$replicate,
      threshold_label(info$threshold)
    ))
    next
  }

  status <- "failed"
  error_message <- info$preparation_error
  beta <- se <- p <- pi0 <- nvalid <- NA_real_
  n_clumped <- 0L

  if (info$preparation_status == "ready" && is.na(clump_error)) {
    candidates <- candidates_by_job[[info$job_id]]
    iv_snps <- selected_snps(selected, info$job_id)
    matched <- match(iv_snps, candidates$SNP)
    clumped <- candidates[matched[!is.na(matched)]]
    n_clumped <- nrow(clumped)

    if (n_clumped < 5L) {
      error_message <- "MR-APSS requires at least five clumped IVs"
    } else {
      fit_attempt <- tryCatch({
        nuisance <- target_background(
          background_cache[[as.character(info$replicate)]]
        )
        fit_data <- data.frame(
          SNP = clumped$SNP,
          A1 = clumped$A1,
          A2 = clumped$A2,
          b.exp = clumped$b.exp.pop2,
          se.exp = clumped$se.exp.pop2,
          pval.exp = clumped$pval.exp.pop2,
          b.out = clumped$b.out.pop2,
          se.out = clumped$se.out.pop2,
          pval.out = clumped$pval.out.pop2,
          L2 = clumped$L2.pop2,
          Threshold = rep(info$threshold, n_clumped)
        )
        fit <- MRAPSS::MRAPSS(
          fit_data,
          C = nuisance$C,
          Omega = nuisance$Omega0,
          Cor.SelectionBias = TRUE
        )
        beta_fit <- scalar_or_na(fit$beta)
        se_fit <- scalar_or_na(fit$beta.se)
        if (!is.finite(beta_fit) || !is.finite(se_fit) || se_fit <= 0) {
          stop("MR-APSS returned a non-finite estimate")
        }
        list(
          beta = beta_fit,
          se = se_fit,
          p = scalar_or_na(fit$pvalue),
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
        pi0 <- fit_attempt$pi0
        nvalid <- n_clumped * pi0
      }
    }
  } else if (info$preparation_status == "ready") {
    error_message <- clump_error
  }

  result <- simulation_result(
    cfg, info$replicate, "MRAPSS_TARGET", info$threshold,
    selection_mode, "EAS", beta = beta, se = se, p = p,
    n_iv = n_clumped, status = status, error = error_message
  )
  result[, `:=`(
    analysis_signature = info$analysis_signature,
    n_candidates = info$n_candidates,
    n_clumped_IV = n_clumped,
    clump_kb = cfg$clump_kb,
    clump_r2 = cfg$clump_r2,
    clump_panel = "EAS",
    selection_bias_correction = TRUE,
    pi0 = pi0,
    nvalid = nvalid
  )]
  write_checkpoint(list(result), result_file, "analysis_signature")
  report_method_result(result, run_index, total_runs)
}

final <- if (file.exists(result_file)) read_tsv(result_file) else data.table::data.table()
message(
  "MR-APSS finished: ", sum(final$status == "ok"), "/", nrow(final),
  " successful rows; results: ", result_file
)
