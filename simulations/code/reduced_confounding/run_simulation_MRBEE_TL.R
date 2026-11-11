#!/usr/bin/env Rscript

# MRBEE simulation workflow:
#   1. select joint-FEMA hits and perform EUR LD clumping once;
#   2. apply the GxE screen and define one shared final_iv set;
#   3. fit source MRBEE, target MRBEE, and MRBEE-TL on final_iv;
#   4. source, target, and TL rows are appended and reported separately.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else {
  getwd()
}
source(file.path(SCRIPT_DIR, "simulation_method_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
cfg <- parse_method_config(args, "MRBEE_TL", SCRIPT_DIR)
cfg$refresh_stale_clump <- TRUE
require_method_packages(c("data.table", "MRBEE", "MRBEEX", "XMAP"))

# Keep the installed MRBEEX package unchanged; the simulation uses a local
# identity-LD copy with corrected half-sample inference and bounded retries.
source(file.path(SCRIPT_DIR, "MRBEE_TL_UV_Independent_simulation.R"))
environment(MRBEE_TL_UV_Independent_simulation) <- asNamespace("MRBEEX")
message("MRBEE-TL simulation uses the local half-sample/retry fixes; MRBEEX package is unchanged")

gxe_p_threshold <- as.numeric(get_arg(args, "--gxe-p-threshold", "5e-8"))
sampling_time <- as.integer(get_arg(args, "--sampling-time", "1000"))
if (!is.finite(gxe_p_threshold) || gxe_p_threshold <= 0 ||
    gxe_p_threshold >= 1) {
  stop("--gxe-p-threshold must be in (0, 1)")
}
if (is.na(sampling_time) || sampling_time < 1L) {
  stop("--sampling-time must be positive")
}

# These two short helpers only translate data into the MRBEE API.
make_mr_data <- function(dat, population) {
  if (population == "EUR") {
    data.frame(
      SNP = dat$SNP,
      bx = dat$b.exp.pop1,
      by = dat$b.out.pop1,
      bxse = dat$se.exp.pop1,
      byse = dat$se.out.pop1
    )
  } else {
    data.frame(
      SNP = dat$SNP,
      bx = dat$b.exp.pop2,
      by = dat$b.out.pop2,
      bxse = dat$se.exp.pop2,
      byse = dat$se.out.pop2
    )
  }
}

fit_mrbee <- function(dat, C) {
  if (nrow(dat) < 3L) stop("MRBEE requires at least three IVs")
  fit <- MRBEE::MRBEE.IMRP.UV(
    by = as.numeric(dat$by),
    bx = as.numeric(dat$bx),
    byse = as.numeric(dat$byse),
    bxse = as.numeric(dat$bxse),
    Rxy = C
  )
  beta <- scalar_or_na(fit$theta)
  variance <- scalar_or_na(fit$vartheta)
  if (!is.finite(beta) || !is.finite(variance) || variance <= 0) {
    stop("MRBEE returned an invalid estimate")
  }
  list(beta = beta, variance = variance, se = sqrt(variance))
}

files <- discover_summary_files(cfg)
result_file <- prepare_result_file(cfg, "simulation_MRBEE_TL_results.tsv")
analysis_signatures <- stats::setNames(
  rep(NA_character_, nrow(files)), as.character(files$replicate)
)
resume_signatures <- analysis_signatures
fingerprint_parameters <- list(
  setting_id = cfg$setting_id,
  master_seed = cfg$seed,
  expected_snps = cfg$expected_snps,
  beta_true = cfg$beta_true,
  rho = cfg$rho,
  n_aux = cfg$n_aux,
  n_target = cfg$n_target,
  aux_z2_threshold = cfg$aux_z2_threshold,
  gxe_p_threshold = gxe_p_threshold,
  gxe_adjust_method = "Sidak",
  sampling_time = sampling_time,
  clump_kb = cfg$clump_kb,
  clump_r2 = cfg$clump_r2,
  clump_population = "EUR",
  python = cfg$python,
  pip_thres = 0.5,
  ebic_delta = 0,
  ebic_gamma = 1,
  reliability_thres = 0.5
)
total_runs <- nrow(files) * length(cfg$thresholds) * 3L
run_index <- 0L
selection_qc <- list()

message(
  "MRBEE setting=", cfg$setting_id, "; replicates=", nrow(files),
  "; thresholds=", paste(cfg$thresholds, collapse = ","),
  "; estimators=SOURCE,TARGET,TL; planned fits=", total_runs
)

# Step 1: prepare one joint-FEMA EUR-clumping job per replicate/threshold.
batch_jobs <- list()
background_cache <- list()
source_background_cache <- list()
for (file_index in seq_len(nrow(files))) {
  replicate <- files$replicate[file_index]
  summary_path <- files$path[file_index]
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
  source_background <- if (inherits(dat, "error")) {
    dat
  } else {
    tryCatch(
      load_or_estimate_aux_background(cfg, replicate, summary_path),
      error = function(condition) condition
    )
  }
  replicate_name <- as.character(replicate)
  background_cache[[replicate_name]] <- background
  source_background_cache[[replicate_name]] <- source_background
  # The auxiliary cache must be checked or refreshed before its bytes enter
  # the signature. Failures still receive a signature for the current inputs.
  analysis_signature <- analysis_fingerprint(
    summary_path,
    background_paths = c(
      target = file.path(
        cfg$background_dir, sprintf("rep%03d", replicate),
        "background_parameters.rds"
      ),
      auxiliary = file.path(
        cfg$aux_background_dir, sprintf("rep%03d", replicate),
        "background_parameters.rds"
      ),
      common_snps = cfg$common_snps,
      eur_schema = cfg$eur_schema
    ),
    code_paths = c(
      file.path(SCRIPT_DIR, "simulation_method_utils.R"),
      file.path(SCRIPT_DIR, "run_simulation_MRBEE_TL.R"),
      file.path(SCRIPT_DIR, "MRBEE_TL_UV_Independent_simulation.R"),
      cfg$clumper
    ),
    parameters = fingerprint_parameters,
    packages = c("MRBEE", "MRBEEX", "XMAP")
  )
  analysis_signatures[[replicate_name]] <- analysis_signature
  if (inherits(dat, "error") || inherits(background, "error") ||
      inherits(source_background, "error")) next
  resume_signatures[[replicate_name]] <- analysis_signature

  z_eur <- dat$b.exp.pop1 / dat$se.exp.pop1
  z_eas <- dat$b.exp.pop2 / dat$se.exp.pop2
  chi_square <- z_eur^2 / source_background$C["EUR_X", "EUR_X"] +
    z_eas^2 / background$C["EAS_X", "EAS_X"]
  dat[, p_fema_joint := pmax(
    stats::pchisq(chi_square, df = 2L, lower.tail = FALSE),
    .Machine$double.xmin
  )]

  for (threshold in cfg$thresholds) {
    candidates <- dat[p_fema_joint < threshold]
    if (nrow(candidates)) {
      batch_jobs[[length(batch_jobs) + 1L]] <- data.table::data.table(
        job_id = make_job_id(
          "MRBEE_TL", replicate, threshold, "joint_FEMA", "EUR"
        ),
        population = "EUR",
        SNP = candidates$SNP,
        priority = candidates$p_fema_joint
      )
    }
  }
}
completed_keys <- completed_result_keys(
  result_file, cfg$seed, analysis_signatures = resume_signatures
)

batch_selected <- data.table::data.table()
batch_clump_error <- NA_character_
if (length(batch_jobs)) {
  batch_job_table <- data.table::rbindlist(batch_jobs)
  message(
    "MRBEE batched EUR clumping: ",
    data.table::uniqueN(batch_job_table$job_id), " jobs, ",
    nrow(batch_job_table), " candidate rows"
  )
  batch_attempt <- tryCatch(
    run_vilma_clump(batch_job_table, cfg, "MRBEE_TL"),
    error = function(condition) condition
  )
  if (inherits(batch_attempt, "error")) {
    batch_clump_error <- conditionMessage(batch_attempt)
  } else {
    batch_selected <- batch_attempt
  }
}

# Step 2: fit source, target, and TL rows one threshold at a time.
for (file_index in seq_len(nrow(files))) {
  replicate <- files$replicate[file_index]
  summary_path <- files$path[file_index]
  message(sprintf(
    "Replicate %03d/%03d: reading data and background parameters",
    file_index, nrow(files)
  ))

  dat <- tryCatch(
    read_simulation_summary(summary_path, cfg$expected_snps),
    error = function(condition) condition
  )
  background <- if (inherits(dat, "error")) {
    dat
  } else {
    background_cache[[as.character(replicate)]]
  }
  source_background <- if (inherits(dat, "error")) {
    dat
  } else {
    source_background_cache[[as.character(replicate)]]
  }

  if (!inherits(dat, "error") && !inherits(background, "error") &&
      !inherits(source_background, "error")) {
    z_eur <- dat$b.exp.pop1 / dat$se.exp.pop1
    z_eas <- dat$b.exp.pop2 / dat$se.exp.pop2
    chi_square <- z_eur^2 / source_background$C["EUR_X", "EUR_X"] +
      z_eas^2 / background$C["EAS_X", "EAS_X"]
    dat[, p_fema_joint := pmax(
      stats::pchisq(chi_square, df = 2L, lower.tail = FALSE),
      .Machine$double.xmin
    )]
  }

  for (threshold_index in seq_along(cfg$thresholds)) {
    threshold <- cfg$thresholds[threshold_index]
    joint_job <- make_job_id(
      "MRBEE_TL", replicate, threshold, "joint_FEMA", "EUR"
    )
    methods <- c("MRBEE_SOURCE", "MRBEE_TARGET", "MRBEE_TL")
    populations <- c("EUR", "EAS", "EAS")
    keys <- mapply(
      result_key,
      replicate = replicate,
      method = methods,
      threshold = threshold,
      selection_mode = "EUR_clump",
      population = populations,
      USE.NAMES = FALSE
    )

    global_error <- if (inherits(dat, "error")) {
      conditionMessage(dat)
    } else if (inherits(background, "error")) {
      conditionMessage(background)
    } else if (inherits(source_background, "error")) {
      conditionMessage(source_background)
    } else {
      NA_character_
    }
    joint_candidates <- data.table::data.table()
    if (is.na(global_error)) {
      joint_candidates <- dat[p_fema_joint < threshold]
    }

    selection_qc[[length(selection_qc) + 1L]] <- data.table::data.table(
      replicate = replicate,
      analysis_signature = analysis_signatures[[as.character(replicate)]],
      summary_path = summary_path,
      threshold = threshold,
      clump_mode = "EUR_clump",
      joint_job = joint_job,
      n_joint_candidates = if (is.na(global_error)) {
        nrow(joint_candidates)
      } else {
        NA_integer_
      },
      preparation_status = if (is.na(global_error)) "ready" else "failed",
      preparation_error = global_error
    )

    if (all(keys %in% completed_keys)) {
      for (method in methods) {
        run_index <- run_index + 1L
        message(sprintf(
          "[%d/%d] %s rep=%03d P=%s: already completed",
          run_index, total_runs, method, replicate,
          threshold_label(threshold)
        ))
      }
      next
    }

    # One joint-FEMA EUR-clumped set is screened once for GxE outliers.
    pipeline_error <- if (!is.na(global_error)) global_error else batch_clump_error
    clumped_iv <- data.table::data.table()
    final_iv <- data.table::data.table()
    gxe_removed <- character()
    gxe_completed <- FALSE

    if (is.na(pipeline_error)) {
      iv_snps <- selected_snps(batch_selected, joint_job)
      matched <- match(iv_snps, joint_candidates$SNP)
      clumped_iv <- joint_candidates[matched[!is.na(matched)]]
      if (nrow(clumped_iv) < 3L) {
        pipeline_error <- "Fewer than three joint-FEMA clumped IVs"
      }
    }

    if (is.na(pipeline_error)) {
      source_strong <- which(clumped_iv$pval.exp.pop1 < gxe_p_threshold)
      if (length(source_strong) >= 3L) {
        gxe_fit <- tryCatch(
          MRBEE::MRBEE.IMRP.UV(
            by = clumped_iv$b.exp.pop2[source_strong],
            bx = clumped_iv$b.exp.pop1[source_strong],
            byse = clumped_iv$se.exp.pop2[source_strong],
            bxse = clumped_iv$se.exp.pop1[source_strong],
            Rxy = diag(c(
              source_background$C["EUR_X", "EUR_X"],
              background$C["EAS_X", "EAS_X"]
            )),
            adjust.method = "Sidak"
          ),
          error = function(condition) condition
        )
        if (inherits(gxe_fit, "error")) {
          pipeline_error <- conditionMessage(gxe_fit)
        } else {
          delta <- as.numeric(gxe_fit$delta)
          flagged <- source_strong[which(is.finite(delta) & delta != 0)]
          gxe_removed <- clumped_iv$SNP[flagged]
          gxe_completed <- TRUE
        }
      } else {
        gxe_completed <- TRUE
      }
    }

    if (is.na(pipeline_error)) {
      final_iv <- clumped_iv[!SNP %in% gxe_removed]
      if (nrow(final_iv) < 3L) {
        pipeline_error <- "Fewer than three IVs after GxE filtering"
      }
    }

    # Source and target MRBEE use exactly the same final_iv set.
    source_dat <- make_mr_data(final_iv, "EUR")
    target_dat <- make_mr_data(final_iv, "EAS")
    source_fit <- if (!is.na(pipeline_error)) {
      simpleError(pipeline_error)
    } else {
      tryCatch(
        fit_mrbee(source_dat, source_background$C),
        error = function(condition) condition
      )
    }
    target_fit <- if (!is.na(pipeline_error)) {
      simpleError(pipeline_error)
    } else {
      tryCatch(
        fit_mrbee(target_dat, target_background(background)$C),
        error = function(condition) condition
      )
    }

    source_ok <- !inherits(source_fit, "error")
    source_result <- simulation_result(
      cfg, replicate, "MRBEE_SOURCE", threshold, "EUR_clump", "EUR",
      beta = if (source_ok) source_fit$beta else NA_real_,
      se = if (source_ok) source_fit$se else NA_real_,
      n_iv = nrow(final_iv),
      status = if (source_ok) "ok" else "failed",
      error = if (source_ok) NA_character_ else conditionMessage(source_fit)
    )

    target_ok <- !inherits(target_fit, "error")
    target_result <- simulation_result(
      cfg, replicate, "MRBEE_TARGET", threshold, "EUR_clump", "EAS",
      beta = if (target_ok) target_fit$beta else NA_real_,
      se = if (target_ok) target_fit$se else NA_real_,
      n_iv = nrow(final_iv),
      status = if (target_ok) "ok" else "failed",
      error = if (target_ok) NA_character_ else conditionMessage(target_fit)
    )

    # Transfer learning uses the source estimate as the prior information.
    tl_fit <- NULL
    if (!source_ok || !target_ok) {
      reasons <- c(
        if (!source_ok) paste0("source: ", conditionMessage(source_fit)),
        if (!target_ok) paste0("target: ", conditionMessage(target_fit))
      )
      tl_fit <- simpleError(paste(reasons, collapse = " | "))
    } else if (floor(nrow(final_iv) / 2) <= 1L) {
      tl_fit <- simpleError(paste0(
        "MRBEE-TL requires at least four final IVs: half-sample ",
        "resampling would contain only ", floor(nrow(final_iv) / 2),
        " IV"
      ))
    } else {
      # Stable in full runs, single-replicate restarts, and reordered grids.
      analysis_index <- as.integer(round(-log10(threshold) * 1e6))
      tl_seed <- derive_simulation_seed(
        cfg$seed, component_id = 5L,
        replicate = replicate, analysis_index = analysis_index
      )
      message(sprintf(
        "Starting MRBEE_TL rep=%03d P=%s final_IV=%d half_sample=%d seed=%d",
        replicate, threshold_label(threshold), nrow(final_iv),
        floor(nrow(final_iv) / 2), tl_seed
      ))
      set.seed(tl_seed)
      tl_fit <- tryCatch(
        MRBEE_TL_UV_Independent_simulation(
          by = as.numeric(target_dat$by),
          bX = matrix(as.numeric(target_dat$bx), ncol = 1L),
          byse = as.numeric(target_dat$byse),
          bXse = matrix(as.numeric(target_dat$bxse), ncol = 1L),
          Rxy = target_background(background)$C,
          theta.source = source_fit$beta,
          theta.source.cov = matrix(source_fit$variance, 1L, 1L),
          pip.thres = 0.5,
          ebic.delta = 0,
          ebic.gamma = 1,
          reliability.thres = 0.5,
          sampling.time = sampling_time
        ),
        error = function(condition) condition
      )
      if (!inherits(tl_fit, "error")) {
        tl_beta <- scalar_or_na(tl_fit$theta)
        tl_se <- scalar_or_na(tl_fit$theta.se)
        if (!is.finite(tl_beta) || !is.finite(tl_se) || tl_se <= 0) {
          tl_fit <- simpleError("MRBEE-TL returned a non-finite estimate")
        }
      }
    }

    tl_ok <- !inherits(tl_fit, "error")
    tl_result <- simulation_result(
      cfg, replicate, "MRBEE_TL", threshold, "EUR_clump", "EAS",
      beta = if (tl_ok) scalar_or_na(tl_fit$theta) else NA_real_,
      se = if (tl_ok) scalar_or_na(tl_fit$theta.se) else NA_real_,
      n_iv = nrow(final_iv),
      status = if (tl_ok) "ok" else "failed",
      error = if (tl_ok) NA_character_ else conditionMessage(tl_fit)
    )

    result_rows <- list(source_result, target_result, tl_result)
    for (result_position in seq_along(result_rows)) {
      result <- result_rows[[result_position]]
      result_method <- as.character(result$method[[1L]])
      result[, `:=`(
        analysis_signature = analysis_signatures[[as.character(replicate)]],
        source_population = "EUR",
        target_population = "EAS",
        n_source_candidates = if (is.na(global_error)) {
          nrow(joint_candidates)
        } else {
          NA_integer_
        },
        n_joint_candidates = if (is.na(global_error)) {
          nrow(joint_candidates)
        } else {
          NA_integer_
        },
        n_gxe_removed = if (gxe_completed) {
          length(gxe_removed)
        } else {
          NA_integer_
        },
        clump_kb = cfg$clump_kb,
        clump_r2 = cfg$clump_r2,
        target_clump_panel = "EUR",
        sampling_time = if (result_method == "MRBEE_TL") {
          sampling_time
        } else {
          NA_integer_
        }
      )]
      run_index <- run_index + 1L
      if (keys[result_position] %in% completed_keys) {
        message(sprintf(
          "[%d/%d] %s rep=%03d P=%s: already completed",
          run_index, total_runs, result$method, replicate,
          threshold_label(threshold)
        ))
      } else {
        write_checkpoint(list(result), result_file, "analysis_signature")
        report_method_result(result, run_index, total_runs)
      }
    }
  }
}

data.table::fwrite(
  data.table::rbindlist(selection_qc, fill = TRUE),
  file.path(cfg$output_dir, "MRBEE_TL_selection_qc.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
final <- if (file.exists(result_file)) read_tsv(result_file) else data.table::data.table()
message(
  "MRBEE finished: ", sum(final$status == "ok"), "/", nrow(final),
  " successful rows; results: ", result_file
)
