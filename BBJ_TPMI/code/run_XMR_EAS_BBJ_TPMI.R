#!/usr/bin/env Rscript

# Resolve repository-relative inputs independently of the launch directory.
.xmr_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(.xmr_args)) stop("Run this entry point with Rscript.")
.xmr_root <- dirname(normalizePath(sub("^--file=", "", .xmr_args[[1L]])))
while (!file.exists(file.path(.xmr_root, ".xmr-reproduce-root"))) {
  .xmr_parent <- dirname(.xmr_root)
  if (identical(.xmr_parent, .xmr_root)) stop("Repository root not found.")
  .xmr_root <- .xmr_parent
}
setwd(.xmr_root)

############################################################
# Paper analysis: XMR in BBJ and TPMI
#
# Auxiliary population: EUR
# Target populations:    BBJ and TPMI (TW)
# IV thresholds:         5e-05 and 5e-08; one master clump at 5e-05 per pair
# IV workflow:           EUR P selection -> FEMA P priority -> EAS LD clump
# No exposure-effect heterogeneity filter or additional selection-bias adjustment.
# Background parameters: re-estimated by two-stage LDSC; cross-population intercepts fixed at 0
#
# All paths and analysis settings are specified directly below.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

if (!requireNamespace("XMR", quietly = TRUE)) {
  stop("The XMR package is required")
}
if (!requireNamespace("MRAPSS", quietly = TRUE)) {
  stop("The MRAPSS package is required for its PLINK clump wrapper")
}

if (!requireNamespace("XMAP", quietly = TRUE)) {
  stop("The XMAP package is required for background LDSC estimation")
}

set.seed(1)

## Fixed paths
official_xmr_dir <- file.path(.xmr_root, "BBJ_TPMI")
formatted_dir <- file.path(official_xmr_dir, "formatted_data")
summary_exp_file <- file.path(
  official_xmr_dir,
  "raw_data/BBJ/EUR_BBJ_exposure_description.csv"
)
summary_out_file <- file.path(
  official_xmr_dir,
  "raw_data/BBJ/EUR_BBJ_outcome_description.csv"
)
ldscore_file <- file.path(.xmr_root, "negative_controls", "LDscore_eas_brit_afr.txt")

reference_dir <- file.path(dirname(official_xmr_dir), "negative_controls")
plink_bin <- file.path(reference_dir, "plink_Linux")
eas_bfile <- file.path(reference_dir, "1kg_pops/EAS_phase3_hg38_autosomes")

result_root <- file.path(.xmr_root, "BBJ_TPMI/results/masterClump5e-05")

## Fixed analysis settings
thresholds <- c(5e-5, 5e-8)
master_clump_threshold <- 5e-5
format_threshold <- function(x) sprintf("%.0e", x)
clump_kb <- 10000L
clump_r2 <- 0.001
# Resume successful analyses; failed analyses are retried. Set TRUE to rerun all.
# Background parameters are re-estimated once for every pending trait pair.
overwrite_results <- FALSE
bg_z2_threshold <- 30
bg_nblocks <- 200L

xmr_scalar_or_na <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  value <- suppressWarnings(as.numeric(x[1L]))
  if (length(value) == 0L || !is.finite(value)) NA_real_ else value
}

xmr_matrix_value_or_na <- function(x, i, j) {
  if (is.null(x) || length(dim(x)) != 2L || nrow(x) < i || ncol(x) < j) {
    return(NA_real_)
  }
  xmr_scalar_or_na(x[i, j])
}

xmr_complement_allele <- function(x) {
  chartr("ATCG", "TAGC", toupper(x))
}

run_xmr_eas <- function(target_pop) {
  target_pop <- toupper(target_pop)
  if (!target_pop %in% c("BBJ", "TW")) stop("target_pop must be BBJ or TW")

  target_label <- if (target_pop == "TW") "TPMI" else "BBJ"
  source_pop <- "EUR"
  out_dir <- file.path(result_root, target_label, "crossPopIntercept0")
  bgparas_dir <- file.path(out_dir, "bgparas")

  required_dirs <- formatted_dir
  missing_dirs <- required_dirs[!dir.exists(required_dirs)]
  if (length(missing_dirs) > 0L) {
    stop("Missing required directories: ", paste(missing_dirs, collapse = ", "))
  }

  required_files <- c(summary_exp_file, summary_out_file, ldscore_file, plink_bin)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop("Missing required files: ", paste(missing_files, collapse = ", "))
  }

  required_bfiles <- paste0(eas_bfile, c(".bed", ".bim", ".fam"))
  missing_bfiles <- required_bfiles[!file.exists(required_bfiles)]
  if (length(missing_bfiles) > 0L) {
    stop("Missing EAS PLINK reference files: ", paste(missing_bfiles, collapse = ", "))
  }

  dir.create(bgparas_dir, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(
    out_dir,
    paste0("XMR_EAS_", target_label, "_EURselected_two_thresholds_noEffectHetQC.tsv")
  )

  exposure_summary <- fread(summary_exp_file)
  outcome_summary <- fread(summary_out_file)
  if (!"UKB.code" %in% names(exposure_summary)) {
    stop("Exposure summary must contain UKB.code: ", summary_exp_file)
  }
  if (!"phecode_unique" %in% names(outcome_summary)) {
    stop("Outcome summary must contain phecode_unique: ", summary_out_file)
  }

  all_exposures <- unique(trimws(as.character(exposure_summary$UKB.code)))
  all_outcomes <- unique(trimws(as.character(outcome_summary$phecode_unique)))
  all_exposures <- all_exposures[!is.na(all_exposures) & nzchar(all_exposures)]
  all_outcomes <- all_outcomes[!is.na(all_outcomes) & nzchar(all_outcomes)]

  exposures <- all_exposures
  outcomes <- all_outcomes
  if (length(exposures) == 0L) stop("No valid exposures selected")
  if (length(outcomes) == 0L) stop("No valid outcomes selected")

  exposure_name_col <- if ("UKB.Phenotype" %in% names(exposure_summary)) {
    "UKB.Phenotype"
  } else {
    "UKB.code"
  }
  outcome_name_col <- if ("UKB.Phenotype" %in% names(outcome_summary)) {
    "UKB.Phenotype"
  } else {
    "phecode_unique"
  }

  exposure_names <- setNames(
    as.character(exposure_summary[[exposure_name_col]]),
    as.character(exposure_summary$UKB.code)
  )
  outcome_names <- setNames(
    as.character(outcome_summary[[outcome_name_col]]),
    as.character(outcome_summary$phecode_unique)
  )

  trait_pairs <- data.table(
    outcome = rep(outcomes, each = length(exposures)),
    exposure = rep(exposures, times = length(outcomes))
  )
  trait_pairs[, `:=`(
    exposure_name = unname(exposure_names[exposure]),
    outcome_name = unname(outcome_names[outcome])
  )]
  if (nrow(trait_pairs) == 0L) stop("No exposure-outcome pairs to run")

  make_analysis_id <- function(exposure, outcome, threshold) {
    paste(source_pop, target_label, exposure, outcome, format_threshold(threshold), sep = "|")
  }

  results_started <-
    !overwrite_results && file.exists(result_file) && file.info(result_file)$size > 0L
  completed <- character(0)
  if (results_started) {
    existing <- fread(result_file, na.strings = c("", "NA", "NaN", "NULL", "--"))
    if (!"analysis_id" %in% names(existing)) {
      stop("Existing result file lacks analysis_id: ", result_file)
    }
    existing <- existing[status == "ok"]
    fwrite(existing, result_file, sep = "\t", quote = FALSE, na = "NA")
    completed <- unique(as.character(existing$analysis_id))
  }

  read_gwas <- function(population, trait_code) {
    stem <- file.path(formatted_dir, paste0(population, "_", trait_code))
    candidates <- c(stem, paste0(stem, ".tsv"), paste0(stem, ".txt"), paste0(stem, ".gz"))
    filename <- candidates[file.exists(candidates)][1L]
    if (is.na(filename)) stop("Missing GWAS file: ", stem)

    dat <- fread(filename)
    required <- c("SNP", "A1", "A2", "Z", "N", "P")
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0L) {
      stop(filename, " is missing columns: ", paste(missing, collapse = ", "))
    }

    dat[, `:=`(
      SNP = as.character(SNP),
      A1 = toupper(as.character(A1)),
      A2 = toupper(as.character(A2)),
      Z = suppressWarnings(as.numeric(Z)),
      N = suppressWarnings(as.numeric(N)),
      P = suppressWarnings(as.numeric(P))
    )]
    dat <- dat[
      !is.na(SNP) & nzchar(SNP) &
        A1 %in% c("A", "C", "G", "T") & A2 %in% c("A", "C", "G", "T") & A1 != A2 &
        is.finite(Z) & is.finite(N) & N > 0 & is.finite(P) & P >= 0 & P <= 1
    ]
    dat[!duplicated(SNP)]
  }

  load_ldscore <- function() {
    ldscore <- fread(ldscore_file)
    required <- c("rsid", "allele1", "allele2", "EUR", "EAS", "EUR_EAS")
    missing <- setdiff(required, names(ldscore))
    if (length(missing) > 0L) {
      stop(ldscore_file, " is missing columns: ", paste(missing, collapse = ", "))
    }

    ldscore[, `:=`(
      rsid = as.character(rsid),
      allele1 = toupper(as.character(allele1)),
      allele2 = toupper(as.character(allele2)),
      EUR = suppressWarnings(as.numeric(EUR)),
      EAS = suppressWarnings(as.numeric(EAS)),
      EUR_EAS = suppressWarnings(as.numeric(EUR_EAS))
    )]
    palindromic <- is.na(ldscore$allele1) | is.na(ldscore$allele2) |
      ldscore$allele1 == xmr_complement_allele(ldscore$allele2)
    ldscore <- ldscore[!palindromic]

    answer <- ldscore[, .(
      SNP = rsid,
      LD_A1 = allele1,
      LD_A2 = allele2,
      L2.pop1 = EUR,
      L2.pop2 = EAS,
      L12 = EUR_EAS
    )]
    answer <- answer[
      !is.na(SNP) & nzchar(SNP) &
        is.finite(L2.pop1) & is.finite(L2.pop2) & is.finite(L12)
    ]
    answer[!duplicated(SNP)]
  }

  estimate_background <- function(exposure_code, outcome_code, gwas) {
    filename <- file.path(bgparas_dir, paste0("bgparas_", exposure_code, "_", outcome_code, ".tsv"))
    trait_names <- c(paste0(source_pop, "_", exposure_code),
                     paste0(target_pop, "_", exposure_code),
                     paste0(target_pop, "_", outcome_code))
    populations <- c(source_pop, target_pop, target_pop)
    score_columns <- c("L2.pop1", "L2.pop2", "L2.pop2")

    # Harmonize both alleles to the LD reference before taking Z products.
    align_background <- function(dat) {
      ans <- merge(dat[, .(SNP, A1, A2, Z, N)], ldscore_xmr, by = "SNP")
      same <- (ans$A1 == ans$LD_A1 & ans$A2 == ans$LD_A2) |
        (xmr_complement_allele(ans$A1) == ans$LD_A1 &
           xmr_complement_allele(ans$A2) == ans$LD_A2)
      flip <- (ans$A1 == ans$LD_A2 & ans$A2 == ans$LD_A1) |
        (xmr_complement_allele(ans$A1) == ans$LD_A2 &
           xmr_complement_allele(ans$A2) == ans$LD_A1)
      ans$Z[flip & !same] <- -ans$Z[flip & !same]
      ans[same | flip]
    }
    aligned <- lapply(gwas, align_background)
    rows <- list()
    for (i in 1:3) for (j in i:3) {
      snps <- intersect(aligned[[i]]$SNP, aligned[[j]]$SNP)
      d1 <- aligned[[i]][match(snps, SNP)]
      d2 <- aligned[[j]][match(snps, SNP)]
      l1 <- d1[[score_columns[i]]]
      l2 <- d2[[score_columns[j]]]
      cross_pop <- populations[i] != populations[j]
      lx <- if (cross_pop) d1$L12 else l1
      idx <- which(d1$Z^2 < bg_z2_threshold & d2$Z^2 < bg_z2_threshold)
      if (length(idx) <= bg_nblocks || length(snps) <= bg_nblocks) {
        stop("Too few background SNPs for LDSC: ", trait_names[i], " vs ", trait_names[j])
      }
      w1 <- 1 / pmax(l1, 1)
      w2 <- 1 / pmax(l2, 1)
      ss1 <- as.data.frame(d1[, .(Z, N)])
      ss2 <- as.data.frame(d2[, .(Z, N)])
      # For cross-population pairs, stage 1 estimates only the two marginal
      # intercepts. Their cross intercept is constrained to zero in stage 2.
      if (cross_pop) {
        marginal1 <- XMAP::estimate_gc(
          ss1[idx, ], ss1[idx, ], l1[idx], l1[idx], l1[idx],
          reg_w1 = w1[idx], reg_w2 = w1[idx], reg_wx = w1[idx],
          constrain_intercept = FALSE, nblocks = bg_nblocks)
        marginal2 <- XMAP::estimate_gc(
          ss2[idx, ], ss2[idx, ], l2[idx], l2[idx], l2[idx],
          reg_w1 = w2[idx], reg_w2 = w2[idx], reg_wx = w2[idx],
          constrain_intercept = FALSE, nblocks = bg_nblocks)
        int1 <- marginal1$tau1$coefs[1]
        int2 <- marginal2$tau1$coefs[1]
        intx <- 0
      } else {
        fit1 <- XMAP::estimate_gc(
          ss1[idx, ], ss2[idx, ], l1[idx], l2[idx], lx[idx],
          reg_w1 = w1[idx], reg_w2 = w2[idx],
          reg_wx = sqrt(w1[idx] * w2[idx]),
          constrain_intercept = FALSE, nblocks = bg_nblocks)
        int1 <- fit1$tau1$coefs[1]
        int2 <- fit1$tau2$coefs[1]
        intx <- fit1$theta$coefs[1]
      }
      fit2 <- XMAP::estimate_gc(
        ss1, ss2, l1, l2, lx,
        reg_w1 = w1, reg_w2 = w2, reg_wx = sqrt(w1 * w2),
        constrain_intercept = TRUE, int1 = int1, int2 = int2,
        intx = intx, nblocks = bg_nblocks)
      if (cross_pop && fit2$theta$coefs[1] != 0) {
        stop("Cross-population LDSC intercept was not fixed at zero")
      }
      rows[[length(rows) + 1L]] <- data.frame(
        trait1 = trait_names[i], trait2 = trait_names[j],
        pop1 = populations[i], pop2 = populations[j],
        Omega11 = fit2$tau1$coefs[2], Omega12 = fit2$theta$coefs[2],
        Omega22 = fit2$tau2$coefs[2], c1 = fit2$tau1$coefs[1],
        c12 = fit2$theta$coefs[1], c2 = fit2$tau2$coefs[1])
    }
    paras <- do.call(rbind, rows)
    fwrite(paras, filename, sep = "\t", quote = FALSE)

    expected_traits <- c(
      paste0(source_pop, "_", exposure_code),
      paste0(target_pop, "_", exposure_code),
      paste0(target_pop, "_", outcome_code)
    )

    make_matrix <- function(value_column) {
      value_matrix <- matrix(
        NA_real_,
        nrow = 3L,
        ncol = 3L,
        dimnames = list(expected_traits, expected_traits)
      )
      for (row_index in seq_len(nrow(paras))) {
        trait1 <- as.character(paras$trait1[row_index])
        trait2 <- as.character(paras$trait2[row_index])
        if (trait1 %in% expected_traits && trait2 %in% expected_traits) {
          value <- suppressWarnings(as.numeric(paras[[value_column]][row_index]))
          value_matrix[trait1, trait2] <- value
          value_matrix[trait2, trait1] <- value
        }
      }
      if (any(!is.finite(value_matrix))) {
        stop("Saved background parameters lack required entries: ", filename)
      }
      value_matrix
    }

    Omega0 <- make_matrix("Omega12")
    C <- make_matrix("c12")
    stopifnot(all(C[1, 2:3] == 0), all(C[2:3, 1] == 0))
    if (any(diag(Omega0) <= 0)) stop("Omega0 has a non-positive diagonal: ", filename)
    if (any(diag(C) <= 0)) stop("C has a non-positive diagonal: ", filename)

    for (i in seq_len(nrow(Omega0) - 1L)) {
      for (j in seq.int(i + 1L, ncol(Omega0))) {
        correlation <- Omega0[i, j] / sqrt(Omega0[i, i] * Omega0[j, j])
        if (is.finite(correlation) && abs(correlation) > 1) {
          Omega0[i, j] <- sqrt(Omega0[i, i] * Omega0[j, j]) * 0.95 * sign(correlation)
          Omega0[j, i] <- Omega0[i, j]
        }
      }
    }

    list(Omega0 = Omega0, C = C, filename = filename)
  }

  add_fema_p <- function(dat) {
    weight_source <- 1 / dat$se.exp.pop1^2
    weight_target <- 1 / dat$se.exp.pop2^2
    weight_sum <- weight_source + weight_target
    beta_fema <-
      (weight_source * dat$b.exp.pop1 + weight_target * dat$b.exp.pop2) / weight_sum
    se_fema <- sqrt(1 / weight_sum)
    dat[, pval.exp.fema := pmax(
      2 * pnorm(abs(beta_fema / se_fema), lower.tail = FALSE),
      .Machine$double.xmin
    )]
    if (any(!is.finite(dat$pval.exp.fema))) stop("Non-finite FEMA P values")
    dat
  }

  ldscore_xmr <- load_ldscore()

  prepare_pair <- function(exposure_code, outcome_code) {
    source_exposure <- read_gwas(source_pop, exposure_code)
    target_exposure <- read_gwas(target_pop, exposure_code)
    target_outcome <- read_gwas(target_pop, outcome_code)

    # Estimate on genome-wide data before IV selection; never reuse saved parameters.
    background <- estimate_background(
      exposure_code, outcome_code,
      list(source_exposure, target_exposure, target_outcome))

    target_pair <- XMR::align_data(
      dat1 = as.data.frame(target_exposure),
      dat2 = as.data.frame(target_outcome),
      trait1.name = paste0(target_pop, "_", exposure_code),
      trait2.name = paste0(target_pop, "_", outcome_code),
      h2.fix.intercept = FALSE,
      LDSC = FALSE
    )$dat
    if (ncol(target_pair) != 9L) stop("Unexpected target align_data output")
    target_pair <- as.data.table(target_pair)
    setnames(target_pair, c(
      "SNP", "A1", "A2", "b.exp.pop2", "b.out.pop2",
      "se.exp.pop2", "se.out.pop2", "pval.exp.pop2", "pval.out.pop2"
    ))

    exposure_pair <- XMR::align_data(
      dat1 = as.data.frame(target_exposure),
      dat2 = as.data.frame(source_exposure),
      trait1.name = paste0(target_pop, "_", exposure_code),
      trait2.name = paste0(source_pop, "_", exposure_code),
      h2.fix.intercept = FALSE,
      LDSC = FALSE
    )$dat
    if (ncol(exposure_pair) != 9L) stop("Unexpected exposure align_data output")
    exposure_pair <- as.data.table(exposure_pair)
    setnames(exposure_pair, c(
      "SNP", "A1", "A2", "b.exp.pop2", "b.exp.pop1",
      "se.exp.pop2", "se.exp.pop1", "pval.exp.pop2", "pval.exp.pop1"
    ))

    dat <- merge(
      target_pair,
      exposure_pair[, .(SNP, A1, A2, b.exp.pop1, se.exp.pop1, pval.exp.pop1)],
      by = c("SNP", "A1", "A2")
    )
    dat <- merge(dat, ldscore_xmr, by = "SNP")
    candidates <- dat[pval.exp.pop1 < master_clump_threshold]
    if (nrow(candidates) == 0L) stop("No EUR-exposure candidate at P < 5e-05")

    list(
      data = candidates,
      background = background,
      source_N = median(source_exposure$N, na.rm = TRUE),
      target_N = median(target_exposure$N, na.rm = TRUE)
    )
  }

  make_result <- function(
      exposure,
      exposure_name,
      outcome,
      outcome_name,
      threshold,
      status,
      error = NA_character_,
      fit = NULL,
      prepared = NULL,
      n_candidates = NA_integer_,
      n_iv = NA_integer_) {
    background <- if (is.null(prepared)) NULL else prepared$background
    data.table(
      analysis_id = make_analysis_id(exposure, outcome, threshold),
      method = "XMR",
      auxiliary = source_pop,
      target = target_label,
      target_population_code = target_pop,
      exposure = exposure,
      exposure_name = exposure_name,
      outcome = outcome,
      outcome_name = outcome_name,
      threshold = threshold,
      threshold_label = format_threshold(threshold),
      master_clump_threshold = master_clump_threshold,
      iv_selection_strategy = "nested EUR-P subsets of 5e-05 master-clumped IVs",
      beta = if (is.null(fit)) NA_real_ else xmr_scalar_or_na(fit$beta),
      se = if (is.null(fit)) NA_real_ else xmr_scalar_or_na(fit$beta.se),
      p = if (is.null(fit)) NA_real_ else xmr_scalar_or_na(fit$beta.pvalue),
      tau_sq = if (is.null(fit)) NA_real_ else xmr_scalar_or_na(fit$tau.sq),
      sigmax11 = if (is.null(fit)) NA_real_ else xmr_matrix_value_or_na(fit$SigmaX, 1, 1),
      sigmax12 = if (is.null(fit)) NA_real_ else xmr_matrix_value_or_na(fit$SigmaX, 1, 2),
      sigmax22 = if (is.null(fit)) NA_real_ else xmr_matrix_value_or_na(fit$SigmaX, 2, 2),
      n_IV = if (is.null(fit)) n_iv else xmr_scalar_or_na(fit$nIV),
      IVsignal_sum = if (is.null(fit)) NA_real_ else xmr_scalar_or_na(fit$post$IVsignal.sum),
      nvalid = if (is.null(fit)) NA_real_ else xmr_scalar_or_na(fit$nvalid),
      pi0 = if (is.null(fit)) NA_real_ else xmr_scalar_or_na(fit$pi0),
      n_candidates = n_candidates,
      source_exposure_N = if (is.null(prepared)) NA_real_ else prepared$source_N,
      target_exposure_N = if (is.null(prepared)) NA_real_ else prepared$target_N,
      selection_rule = paste("EUR exposure P <", format_threshold(threshold)),
      clump_priority = "signed inverse-variance FEMA P",
      clump_reference = "EAS",
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      dt = 0,
      modified_threshold = threshold,
      background_cross_population_intercept = 0,
      background_file = if (is.null(background)) NA_character_ else background$filename,
      status = status,
      error = as.character(error)[1L]
    )
  }

  clump_master_ivs <- function(prepared) {
    candidates <- copy(prepared$data[pval.exp.pop1 < master_clump_threshold])
    if (nrow(candidates) == 0L) stop("No EUR-exposure IV candidate")
    candidates <- add_fema_p(candidates)
    clumped <- MRAPSS::clump(
      as.data.frame(candidates),
      IV.Threshold = 1,
      SNP_col = "SNP",
      pval_col = "pval.exp.fema",
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p = 1,
      bfile = eas_bfile,
      plink_bin = plink_bin
    )
    if (is.null(clumped) || nrow(clumped) == 0L) stop("No IV remains after EAS LD clumping")
    as.data.table(clumped)
  }

  run_pair <- function(prepared, master_ivs, exposure, exposure_name, outcome, outcome_name, threshold) {
    n_candidates <- nrow(prepared$data[pval.exp.pop1 < threshold])
    clumped <- copy(master_ivs[pval.exp.pop1 < threshold])
    n_iv <- nrow(clumped)
    if (n_iv == 0L) stop("No IV remains in the nested EUR-P threshold subset")
    clumped[, Threshold := threshold]

    stopifnot(all(prepared$background$C[1, 2:3] == 0),
              all(prepared$background$C[2:3, 1] == 0))
    tolerance <- if (n_iv < 100L) 1e-8 else 1e-7
    fit <- XMR::fit_XMR(
      data = as.data.frame(clumped),
      C = prepared$background$C,
      Omega0 = prepared$background$Omega0,
      Threshold = threshold,
      tol1 = tolerance,
      tol2 = tolerance,
      min_thres = 1e-2
    )

    make_result(
      exposure = exposure,
      exposure_name = exposure_name,
      outcome = outcome,
      outcome_name = outcome_name,
      threshold = threshold,
      status = "ok",
      fit = fit,
      prepared = prepared,
      n_candidates = n_candidates,
      n_iv = n_iv
    )
  }

  append_result <- function(result) {
    fwrite(
      result,
      result_file,
      sep = "\t",
      quote = FALSE,
      na = "NA",
      append = results_started,
      col.names = !results_started
    )
    results_started <<- TRUE
  }

  message("Target: ", target_label, " (", target_pop, ")")
  message("IV thresholds: ", paste(format_threshold(thresholds), collapse = ", "))
  message("Pairs: ", nrow(trait_pairs))
  message("Re-estimated background parameters (cross-population intercept = 0): ", bgparas_dir)
  message("Results: ", result_file)
  message("Completed pairs already present: ", length(completed))

  for (pair_index in seq_len(nrow(trait_pairs))) {
    exposure <- trait_pairs$exposure[pair_index]
    exposure_name <- trait_pairs$exposure_name[pair_index]
    outcome <- trait_pairs$outcome[pair_index]
    outcome_name <- trait_pairs$outcome_name[pair_index]
    pair_ids <- make_analysis_id(exposure, outcome, thresholds)
    if (all(pair_ids %in% completed)) {
      message("[", pair_index, "/", nrow(trait_pairs), "] Skip existing: ", exposure, " -> ", outcome)
      next
    }

    message("[", pair_index, "/", nrow(trait_pairs), "] ", exposure, " -> ", outcome)
    preparation_error <- NA_character_
    prepared <- tryCatch(
      prepare_pair(exposure, outcome),
      error = function(e) {
        preparation_error <<- conditionMessage(e)
        NULL
      }
    )
    clump_error <- NA_character_
    master_ivs <- if (is.null(prepared)) NULL else tryCatch(
      clump_master_ivs(prepared),
      error = function(e) {
        clump_error <<- conditionMessage(e)
        NULL
      }
    )

    for (threshold in thresholds) {
      analysis_id <- make_analysis_id(exposure, outcome, threshold)
      if (analysis_id %in% completed) next

      failure <- function(error) make_result(
        exposure, exposure_name, outcome, outcome_name, threshold,
        status = "failed", error = error, prepared = prepared,
        n_candidates = if (is.null(prepared)) NA_integer_ else
          nrow(prepared$data[pval.exp.pop1 < threshold]),
        n_iv = if (is.null(master_ivs)) NA_integer_ else
          nrow(master_ivs[pval.exp.pop1 < threshold])
      )
      result <- if (is.null(prepared)) {
        failure(preparation_error)
      } else if (is.null(master_ivs)) {
        failure(clump_error)
      } else {
        tryCatch(
          run_pair(prepared, master_ivs, exposure, exposure_name, outcome, outcome_name, threshold),
          error = function(e) failure(conditionMessage(e))
        )
      }

      append_result(result)
      completed <- c(completed, analysis_id)
      message("  ", format_threshold(threshold), " appended [", result$status, "]")
    }
    rm(prepared, master_ivs)
    invisible(gc(verbose = FALSE))
  }

  message("XMR finished: ", result_file)
  invisible(result_file)
}

## Run the complete paper analysis for both EAS target cohorts.
run_xmr_eas("BBJ")
run_xmr_eas("TW")
