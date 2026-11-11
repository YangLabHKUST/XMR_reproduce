#!/usr/bin/env Rscript

# CAMeRa regional-FEMA workflow:
#   1. select and clump EUR and EAS exposure hits in their own LD panels;
#   2. pool the selected hits into genomic regions;
#   3. construct shared regional FEMA instruments;
#   4. estimate the EAS target effect;
#   5. append and report the result immediately.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else {
  getwd()
}
source(file.path(SCRIPT_DIR, "simulation_method_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
cfg <- parse_method_config(args, "CAMeRa", SCRIPT_DIR)
require_method_packages(c("data.table", "CAMeRa", "tibble"))

region_radius <- as.integer(get_arg(args, "--region-radius", "250000"))
if (is.na(region_radius) || region_radius < 1L) {
  stop("--region-radius must be positive")
}

# CAMeRa 0.1.0 drops matrix dimensions in one-SNP regions. Keep the official
# calculation but force a one-row matrix in that edge case.
fixed_effects_meta_analysis_fast_safe <- function(beta_mat, se_mat) {
  if (is.null(dim(beta_mat))) beta_mat <- matrix(beta_mat, nrow = 1L)
  if (is.null(dim(se_mat))) se_mat <- matrix(se_mat, nrow = 1L)
  weights <- 1 / se_mat^2
  beta <- rowSums(beta_mat * weights, na.rm = TRUE) /
    rowSums(weights, na.rm = TRUE)
  se <- sqrt(1 / rowSums(weights, na.rm = TRUE))
  z <- abs(beta / se)
  tibble::tibble(
    nstudy = apply(beta_mat, 1L, function(x) sum(is.finite(x))),
    p = stats::pnorm(z, lower.tail = FALSE),
    z = z
  )
}
assignInNamespace(
  "fixed_effects_meta_analysis_fast",
  fixed_effects_meta_analysis_fast_safe,
  ns = "CAMeRa"
)

# Convert one population/trait to the column names required by CAMeRa.
make_camera_table <- function(dat, population, trait) {
  exposure <- trait == "X"
  suffix <- if (population == "EUR") "pop1" else "pop2"
  beta_column <- paste0("b.", if (exposure) "exp." else "out.", suffix)
  se_column <- paste0("se.", if (exposure) "exp." else "out.", suffix)
  p_column <- paste0("pval.", if (exposure) "exp." else "out.", suffix)
  n_column <- paste0("N.", if (exposure) "exp." else "out.", suffix)
  data.table::data.table(
    SNP = dat$SNP,
    rsid = paste0(dat$CHR, ":", dat$BP, "_", dat$A1, "_", dat$A2),
    chr = as.character(dat$CHR),
    position = as.integer(dat$BP),
    ea = dat$A1,
    nea = dat$A2,
    beta = dat[[beta_column]],
    se = dat[[se_column]],
    p = dat[[p_column]],
    Z = dat[[beta_column]] / dat[[se_column]],
    N = dat[[n_column]],
    eaf = NA_real_,
    pop = population,
    id = paste(trait, population),
    trait = trait
  )
}

# Regional data assembly is kept as one helper because it is a direct CAMeRa
# object-construction step reused for every replicate and threshold.
build_regional_fema <- function(dat, pooled_snps) {
  populations <- c("EUR", "EAS")
  pooled_hits <- unique(dat[SNP %in% pooled_snps, .(SNP, CHR, BP)])
  if (!nrow(pooled_hits)) stop("No pooled discovery hit remained after clumping")
  regions <- merge_hit_regions(pooled_hits, region_radius)
  if (!nrow(regions)) stop("No CAMeRa region could be constructed")

  exposure <- setNames(lapply(populations, function(population) {
    make_camera_table(dat, population, "X")
  }), populations)
  outcome <- setNames(lapply(populations, function(population) {
    make_camera_table(dat, population, "Y")
  }), populations)

  exposure_regions <- list()
  outcome_regions <- list()
  raw_exposure <- list()
  raw_outcome <- list()
  retained_region <- 1L

  for (region_index in seq_len(nrow(regions))) {
    region <- regions[region_index]
    exp_region <- setNames(lapply(populations, function(population) {
      as.data.frame(exposure[[population]][
        chr == region$CHR & position >= region$start & position <= region$end
      ])
    }), paste("X", populations))
    common_rsid <- Reduce(intersect, lapply(exp_region, `[[`, "rsid"))
    if (!length(common_rsid)) next
    exp_region <- lapply(exp_region, function(table) {
      table <- table[table$rsid %in% common_rsid, , drop = FALSE]
      table[order(table$rsid), , drop = FALSE]
    })
    out_region <- setNames(lapply(populations, function(population) {
      table <- outcome[[population]][rsid %in% common_rsid]
      as.data.frame(table[order(rsid)])
    }), paste("Y", populations))

    combined_exposure <- data.table::rbindlist(exp_region, fill = TRUE)
    lead_rsid <- combined_exposure$rsid[
      which.max(abs(combined_exposure$beta / combined_exposure$se))
    ]
    region_name <- paste0("region_", retained_region)
    exposure_regions[[region_name]] <- exp_region
    outcome_regions[[region_name]] <- out_region
    raw_exposure[[region_name]] <- combined_exposure[rsid == lead_rsid]
    combined_outcome <- data.table::rbindlist(out_region, fill = TRUE)
    raw_outcome[[region_name]] <- combined_outcome[rsid == lead_rsid]
    retained_region <- retained_region + 1L
  }
  if (!length(exposure_regions)) stop("No common regional data remained")

  fit <- CAMeRa::CAMERA$new()
  fit$import_from_local(
    instrument_raw = data.table::rbindlist(raw_exposure, fill = TRUE),
    instrument_outcome = data.table::rbindlist(raw_outcome, fill = TRUE),
    instrument_regions = exposure_regions,
    instrument_outcome_regions = outcome_regions,
    exposure_ids = paste("X", populations),
    outcome_ids = paste("Y", populations),
    pops = populations
  )
  fit$fema_regional_instruments(method = "fema")
  if (is.null(fit$instrument_fema) || !nrow(fit$instrument_fema)) {
    stop("CAMeRa regional FEMA returned no instrument")
  }
  instruments <- data.table::as.data.table(fit$instrument_fema)
  missing <- setdiff(c("SNP", "rsid"), names(instruments))
  if (length(missing)) {
    stop("CAMeRa FEMA instruments lack: ", paste(missing, collapse = ", "))
  }
  list(
    fit = fit,
    instruments = instruments,
    n_regions = length(exposure_regions),
    n_pooled_hits = nrow(pooled_hits)
  )
}

files <- discover_summary_files(cfg)
result_file <- prepare_result_file(cfg, "simulation_CAMeRa_results.tsv")
completed_keys <- completed_result_keys(result_file, cfg$seed)
total_runs <- nrow(files) * length(cfg$thresholds)
run_index <- 0L
selection_qc <- list()

message(
  "CAMeRa setting=", cfg$setting_id, "; replicates=", nrow(files),
  "; thresholds=", paste(cfg$thresholds, collapse = ","),
  "; mode=regional FEMA; planned fits=", total_runs
)

# Step 1: prepare all discovery jobs, then load each ViLMA LD block once.
batch_jobs <- list()
for (file_index in seq_len(nrow(files))) {
  replicate <- files$replicate[file_index]
  dat <- tryCatch(
    read_simulation_summary(files$path[file_index], cfg$expected_snps),
    error = function(condition) condition
  )
  if (inherits(dat, "error")) next
  for (threshold in cfg$thresholds) {
    eur_candidates <- dat[pval.exp.pop1 < threshold]
    eas_candidates <- dat[pval.exp.pop2 < threshold]
    if (nrow(eur_candidates)) {
      batch_jobs[[length(batch_jobs) + 1L]] <- data.table::data.table(
        job_id = make_job_id(
          "CAMeRa", replicate, threshold, "discovery", "EUR"
        ),
        population = "EUR",
        SNP = eur_candidates$SNP,
        priority = eur_candidates$pval.exp.pop1
      )
    }
    if (nrow(eas_candidates)) {
      batch_jobs[[length(batch_jobs) + 1L]] <- data.table::data.table(
        job_id = make_job_id(
          "CAMeRa", replicate, threshold, "discovery", "EAS"
        ),
        population = "EAS",
        SNP = eas_candidates$SNP,
        priority = eas_candidates$pval.exp.pop2
      )
    }
  }
}

batch_selected <- data.table::data.table()
batch_clump_error <- NA_character_
if (length(batch_jobs)) {
  batch_job_table <- data.table::rbindlist(batch_jobs)
  message(
    "CAMeRa batched discovery clumping: ",
    data.table::uniqueN(batch_job_table$job_id), " jobs, ",
    nrow(batch_job_table), " candidate rows"
  )
  batch_attempt <- tryCatch(
    run_vilma_clump(batch_job_table, cfg, "CAMeRa_discovery"),
    error = function(condition) condition
  )
  if (inherits(batch_attempt, "error")) {
    batch_clump_error <- conditionMessage(batch_attempt)
  } else {
    batch_selected <- batch_attempt
  }
}

# Step 2: construct regional FEMA instruments and fit each planned analysis.
for (file_index in seq_len(nrow(files))) {
  replicate <- files$replicate[file_index]
  summary_path <- files$path[file_index]
  message(sprintf(
    "Replicate %03d/%03d: reading %s",
    file_index, nrow(files), basename(summary_path)
  ))
  dat <- tryCatch(
    read_simulation_summary(summary_path, cfg$expected_snps),
    error = function(condition) condition
  )

  for (threshold in cfg$thresholds) {
    run_index <- run_index + 1L
    eur_job <- make_job_id(
      "CAMeRa", replicate, threshold, "discovery", "EUR"
    )
    eas_job <- make_job_id(
      "CAMeRa", replicate, threshold, "discovery", "EAS"
    )
    current_key <- result_key(
      replicate, "CAMeRa_FEMA", threshold, "FEMA_IV", "EAS"
    )

    eur_candidates <- eas_candidates <- data.table::data.table()
    preparation_status <- "failed"
    preparation_error <- if (inherits(dat, "error")) {
      conditionMessage(dat)
    } else {
      NA_character_
    }
    if (!inherits(dat, "error")) {
      eur_candidates <- dat[pval.exp.pop1 < threshold]
      eas_candidates <- dat[pval.exp.pop2 < threshold]
      if (nrow(eur_candidates) + nrow(eas_candidates) > 0L) {
        preparation_status <- "ready"
      } else {
        preparation_error <- "No EUR or EAS exposure SNP passed the threshold"
      }
    }

    selection_qc[[length(selection_qc) + 1L]] <- data.table::data.table(
      replicate = replicate,
      summary_path = summary_path,
      threshold = threshold,
      eur_job = eur_job,
      eas_job = eas_job,
      n_candidates_EUR = if (inherits(dat, "error")) {
        NA_integer_
      } else {
        nrow(eur_candidates)
      },
      n_candidates_EAS = if (inherits(dat, "error")) {
        NA_integer_
      } else {
        nrow(eas_candidates)
      },
      preparation_status = preparation_status,
      preparation_error = preparation_error
    )

    if (current_key %in% completed_keys) {
      message(sprintf(
        "[%d/%d] CAMeRa_FEMA rep=%03d P=%s: already completed",
        run_index, total_runs, replicate, threshold_label(threshold)
      ))
      next
    }

    status <- "failed"
    error_message <- preparation_error
    beta <- se <- p <- NA_real_
    n_iv <- n_regions <- n_pooled_hits <- n_harmonised <- NA_integer_

    if (preparation_status == "ready") {
      selected <- if (is.na(batch_clump_error)) {
        batch_selected
      } else {
        simpleError(batch_clump_error)
      }

      if (inherits(selected, "error")) {
        error_message <- conditionMessage(selected)
      } else {
        pooled_snps <- unique(c(
          selected_snps(selected, eur_job),
          selected_snps(selected, eas_job)
        ))
        n_pooled_hits <- length(pooled_snps)

        regional <- tryCatch(
          build_regional_fema(dat, pooled_snps),
          error = function(condition) condition
        )
        if (inherits(regional, "error")) {
          error_message <- conditionMessage(regional)
        } else {
          n_iv <- data.table::uniqueN(regional$instruments$SNP)
          n_regions <- regional$n_regions
          n_pooled_hits <- regional$n_pooled_hits
          fit_attempt <- tryCatch({
            fit <- regional$fit
            instruments <- regional$instruments
            fit$instrument_fema <- as.data.frame(instruments)
            fit$make_outcome_local(exp = fit$instrument_fema)
            fit$harmonise(exp = fit$instrument_fema)

            harmonised <- data.table::as.data.table(fit$harmonised_dat)
            required <- c("pops", "beta.x", "se.x", "beta.y", "se.y")
            missing <- setdiff(required, names(harmonised))
            if (length(missing)) {
              stop(
                "CAMeRa harmonised data lacks: ",
                paste(missing, collapse = ", ")
              )
            }
            harmonised <- harmonised[
              !is.na(pops) & is.finite(beta.x) & is.finite(se.x) & se.x > 0 &
                is.finite(beta.y) & is.finite(se.y) & se.y > 0
            ]
            if (!nrow(harmonised)) {
              stop("No CAMeRa harmonised associations remained")
            }
            fit$harmonised_dat <- as.data.frame(harmonised)

            estimates <- data.table::as.data.table(fit$cross_estimate())
            renames <- c(
              "Estimate" = "beta",
              "Std. Error" = "se",
              "Pr(>|t|)" = "p"
            )
            for (old_name in names(renames)) {
              if (old_name %in% names(estimates)) {
                data.table::setnames(estimates, old_name, renames[[old_name]])
              }
            }
            required_estimate <- c("pops", "beta", "se")
            missing_estimate <- setdiff(required_estimate, names(estimates))
            if (length(missing_estimate)) {
              stop("CAMeRa cross_estimate returned unexpected columns")
            }
            estimate <- estimates[as.character(pops) == "EAS"]
            if (nrow(estimate) != 1L) {
              stop("CAMeRa did not return exactly one EAS target estimate")
            }
            beta_fit <- suppressWarnings(as.numeric(estimate$beta[[1L]]))
            se_fit <- suppressWarnings(as.numeric(estimate$se[[1L]]))
            p_fit <- if ("p" %in% names(estimate)) {
              suppressWarnings(as.numeric(estimate$p[[1L]]))
            } else {
              p_from_beta_se(beta_fit, se_fit)
            }
            if (!is.finite(beta_fit) || !is.finite(se_fit) || se_fit <= 0) {
              stop("CAMeRa returned a non-finite estimate")
            }
            list(
              beta = beta_fit,
              se = se_fit,
              p = p_fit,
              n_harmonised = nrow(harmonised)
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
            n_harmonised <- fit_attempt$n_harmonised
          }
        }
      }
    }

    result <- simulation_result(
      cfg, replicate, "CAMeRa_FEMA", threshold, "FEMA_IV", "EAS",
      beta = beta, se = se, p = p, n_iv = n_iv,
      status = status, error = error_message
    )
    result[, `:=`(
      n_regions = n_regions,
      n_pooled_hits = n_pooled_hits,
      n_fema_IV_before_post_clump = n_iv,
      n_harmonised_associations = n_harmonised,
      region_radius = region_radius,
      clump_kb = cfg$clump_kb,
      clump_r2 = cfg$clump_r2
    )]

    write_checkpoint(list(result), result_file)
    report_method_result(result, run_index, total_runs)
  }
}

data.table::fwrite(
  data.table::rbindlist(selection_qc, fill = TRUE),
  file.path(cfg$output_dir, "CAMeRa_selection_qc.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
final <- if (file.exists(result_file)) read_tsv(result_file) else data.table::data.table()
message(
  "CAMeRa finished: ", sum(final$status == "ok"), "/", nrow(final),
  " successful rows; results: ", result_file
)
