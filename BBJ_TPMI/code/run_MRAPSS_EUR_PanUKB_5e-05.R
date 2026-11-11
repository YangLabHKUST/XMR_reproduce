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
# Paper analysis: MRAPSS with official EUR exposures and
# PanUKB EUR outcomes for the XMR EAS study
#
# Exposures: 7 BBJ_TPMI EUR exposure files
# Outcomes:  15 PanUKB files named PanUKB_<abbreviation>
# Threshold: 5e-05 only
# LD:        EUR LD score and EUR 1000G reference
# Background parameters: re-estimated with two-stage EUR LDSC
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

for (package_name in c("XMR", "XMAP", "MRAPSS")) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop("The ", package_name, " package is required")
  }
}

set.seed(1)

## Fixed input paths
official_xmr_dir <- file.path(.xmr_root, "BBJ_TPMI")
exposure_dir <- file.path(official_xmr_dir, "formatted_data")
exposure_metadata_file <- file.path(
  official_xmr_dir,
  "raw_data/BBJ/EUR_BBJ_exposure_description.csv"
)
ldscore_file <- file.path(.xmr_root, "negative_controls", "LDscore_eas_brit_afr.txt")

outcome_dir <-
  file.path(.xmr_root, "auxiliary_sensitivity/data/formatted")

reference_dir <- file.path(dirname(official_xmr_dir), "negative_controls")
plink_bin <- file.path(reference_dir, "plink_Linux")
eur_bfile <- file.path(reference_dir, "1kg_pops/EUR_phase3_hg38_autosomes")

## Fixed output paths
analysis_root <- file.path(.xmr_root, "BBJ_TPMI")
background_dir <- file.path(
  analysis_root,
  "background_parameters/MRAPSS_EUR_exposure_PanUKB_outcome"
)
result_dir <- file.path(analysis_root, "results/EUR")
result_file <- file.path(
  result_dir,
  "MRAPSS_EUR_exposure_PanUKB_15outcomes_5e-05.tsv"
)

## Fixed analysis settings
threshold <- 5e-5
threshold_label <- "5e-05"
clump_kb <- 10000L
clump_r2 <- 0.001
overwrite_results <- FALSE

## The 15 outcomes and abbreviations in UKB_15_outcomes_PanUKB_manifest.xlsx.
outcome_table <- data.table(
  outcome_order = seq_len(15L),
  outcome = c(
    "BC", "PrC", "IDA", "UF", "Ep", "Cat", "Ang", "MI",
    "PAD", "UAP", "As", "COPD", "CL", "T2D", "RA"
  ),
  outcome_name = c(
    "Breast cancer",
    "Prostate cancer",
    "Iron deficiency anemia",
    "Uterine fibroid",
    "Epilepsy",
    "Cataract",
    "Angina pectoris",
    "Myocardial infarction",
    "Peripheral arterial disease",
    "Unstable angina pectoris",
    "Asthma",
    "Chronic obstructive pulmonary disease",
    "Cholelithiasis",
    "Type 2 diabetes",
    "Rheumatoid arthritis"
  ),
  panukb_code = c(
    "174.1", "185", "280.1", "218.1", "345.1", "366", "411.3",
    "411.2", "443.9", "411.1", "495", "496.21", "574.1", "250.2",
    "714.1"
  )
)
outcome_table[, outcome_file := paste0("PanUKB_", outcome)]

scalar_or_na <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  value <- suppressWarnings(as.numeric(x[1L]))
  if (length(value) == 0L || !is.finite(value)) NA_real_ else value
}

complement_allele <- function(x) {
  chartr("ATCG", "TAGC", toupper(x))
}

analysis_id_for <- function(exposure_code, outcome_abbreviation) {
  paste(
    "EUR_exposure",
    exposure_code,
    paste0("PanUKB_", outcome_abbreviation),
    threshold_label,
    sep = "|"
  )
}

resolve_file <- function(directory, stem, label) {
  base <- file.path(directory, stem)
  candidates <- c(
    base,
    paste0(base, ".tsv"),
    paste0(base, ".txt"),
    paste0(base, ".gz"),
    paste0(base, ".bgz")
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L) stop("Missing ", label, " GWAS file: ", base)
  hits[1L]
}

resolve_exposure_file <- function(exposure_code) {
  resolve_file(exposure_dir, paste0("EUR_", exposure_code), "EUR exposure")
}

resolve_outcome_file <- function(outcome_abbreviation) {
  resolve_file(outcome_dir, paste0("PanUKB_", outcome_abbreviation), "PanUKB outcome")
}

read_formatted_gwas <- function(filename, trait_label) {
  dat <- fread(filename)
  required <- c("SNP", "A1", "A2", "Z", "N", "P")
  missing_columns <- setdiff(required, names(dat))
  if (length(missing_columns) > 0L) {
    stop(
      trait_label, " file is missing columns: ",
      paste(missing_columns, collapse = ", "),
      "; file: ", filename
    )
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
      A1 %in% c("A", "C", "G", "T") &
      A2 %in% c("A", "C", "G", "T") &
      A1 != A2 &
      is.finite(Z) & is.finite(N) & N > 0 &
      is.finite(P) & P >= 0 & P <= 1
  ]
  dat[!duplicated(SNP), .(SNP, A1, A2, Z, N, P)]
}

## Validate fixed paths and build the 7 x 15 trait grid.
required_dirs <- c(exposure_dir, outcome_dir)
missing_dirs <- required_dirs[!dir.exists(required_dirs)]
if (length(missing_dirs) > 0L) {
  stop("Missing required directories: ", paste(missing_dirs, collapse = ", "))
}

required_files <- c(exposure_metadata_file, ldscore_file, plink_bin)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

required_bfiles <- paste0(eur_bfile, c(".bed", ".bim", ".fam"))
missing_bfiles <- required_bfiles[!file.exists(required_bfiles)]
if (length(missing_bfiles) > 0L) {
  stop("Missing EUR PLINK reference files: ", paste(missing_bfiles, collapse = ", "))
}

exposure_metadata <- fread(exposure_metadata_file)
required_exposure_columns <- c("UKB.code", "UKB.Phenotype")
missing_exposure_columns <- setdiff(
  required_exposure_columns,
  names(exposure_metadata)
)
if (length(missing_exposure_columns) > 0L) {
  stop(
    "Exposure metadata is missing columns: ",
    paste(missing_exposure_columns, collapse = ", ")
  )
}

exposure_table <- unique(exposure_metadata[, .(
  exposure = trimws(as.character(UKB.code)),
  exposure_name = trimws(as.character(UKB.Phenotype))
)])
exposure_table <- exposure_table[
  !is.na(exposure) & nzchar(exposure) &
    !is.na(exposure_name) & nzchar(exposure_name)
]
if (nrow(exposure_table) != 7L) {
  stop("Expected 7 official EUR exposures; found ", nrow(exposure_table))
}
if (nrow(outcome_table) != 15L || anyDuplicated(outcome_table$outcome)) {
  stop("The fixed PanUKB outcome table must contain 15 unique outcomes")
}

exposure_table[, exposure_path := vapply(
  exposure,
  resolve_exposure_file,
  character(1)
)]
outcome_table[, outcome_path := vapply(
  outcome,
  resolve_outcome_file,
  character(1)
)]

trait_pairs <- rbindlist(lapply(seq_len(nrow(exposure_table)), function(i) {
  data.table(
    exposure = exposure_table$exposure[i],
    exposure_name = exposure_table$exposure_name[i],
    exposure_path = exposure_table$exposure_path[i],
    outcome = outcome_table$outcome,
    outcome_name = outcome_table$outcome_name,
    panukb_code = outcome_table$panukb_code,
    outcome_file = outcome_table$outcome_file,
    outcome_path = outcome_table$outcome_path
  )
}))
if (nrow(trait_pairs) != 105L) stop("Expected 105 exposure-outcome pairs")

dir.create(background_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

## EUR LD score for two-stage LDSC and MRAPSS L2.
load_eur_ldscore <- function() {
  ldscore <- fread(ldscore_file)
  required <- c("rsid", "allele1", "allele2", "EUR")
  missing_columns <- setdiff(required, names(ldscore))
  if (length(missing_columns) > 0L) {
    stop(
      "LD-score file is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  ldscore[, `:=`(
    rsid = as.character(rsid),
    allele1 = toupper(as.character(allele1)),
    allele2 = toupper(as.character(allele2)),
    EUR = suppressWarnings(as.numeric(EUR))
  )]
  palindromic <-
    is.na(ldscore$allele1) | is.na(ldscore$allele2) |
    ldscore$allele1 == complement_allele(ldscore$allele2)
  ldscore <- ldscore[
    !palindromic & !duplicated(rsid) &
      !is.na(rsid) & nzchar(rsid) & is.finite(EUR)
  ]
  ldscore[, EUR_EUR := EUR]
  ldscore
}

ldscore_full <- load_eur_ldscore()
ldscore_mrapss <- ldscore_full[, .(SNP = rsid, L2 = EUR)]

align_z_to_ldscore <- function(dat, ldscore) {
  same <-
    (dat$A1 == ldscore$allele1 & dat$A2 == ldscore$allele2) |
    (
      dat$A1 == complement_allele(ldscore$allele1) &
        dat$A2 == complement_allele(ldscore$allele2)
    )
  reverse <-
    (dat$A1 == ldscore$allele2 & dat$A2 == ldscore$allele1) |
    (
      dat$A1 == complement_allele(ldscore$allele2) &
        dat$A2 == complement_allele(ldscore$allele1)
    )

  same[is.na(same)] <- FALSE
  reverse[is.na(reverse)] <- FALSE
  z <- dat$Z
  z[!(same | reverse)] <- NA_real_
  z[reverse] <- -z[reverse]
  z
}

## Two-stage EUR LDSC, following the attached MRAPSS workflow.
estimate_ldsc_entry <- function(dat1, dat2, trait1_label, trait2_label) {
  common_snps <- Reduce(intersect, list(ldscore_full$rsid, dat1$SNP, dat2$SNP))
  if (length(common_snps) < 100L) {
    stop("Too few EUR LDSC SNPs: ", trait1_label, " / ", trait2_label)
  }

  dat1_ldsc <- dat1[match(common_snps, SNP)]
  dat2_ldsc <- dat2[match(common_snps, SNP)]
  ldscore <- ldscore_full[match(common_snps, rsid)]

  z1 <- align_z_to_ldscore(dat1_ldsc, ldscore)
  z2 <- align_z_to_ldscore(dat2_ldsc, ldscore)
  usable <-
    is.finite(z1) & is.finite(z2) &
    is.finite(dat1_ldsc$N) & dat1_ldsc$N > 0 &
    is.finite(dat2_ldsc$N) & dat2_ldsc$N > 0 &
    is.finite(ldscore$EUR)
  stage1 <- usable & z1^2 < 30 & z2^2 < 30
  if (sum(stage1) < 100L || sum(usable) < 100L) {
    stop("Too few usable EUR LDSC SNPs: ", trait1_label, " / ", trait2_label)
  }

  regression_weight <- 1 / pmax(ldscore$EUR, 1)
  fit_stage1 <- XMAP::estimate_gc(
    data.frame(Z = z1[stage1], N = dat1_ldsc$N[stage1]),
    data.frame(Z = z2[stage1], N = dat2_ldsc$N[stage1]),
    ldscore$EUR[stage1],
    ldscore$EUR[stage1],
    ldscore$EUR_EUR[stage1],
    reg_w1 = regression_weight[stage1],
    reg_w2 = regression_weight[stage1],
    reg_wx = regression_weight[stage1],
    constrain_intercept = FALSE
  )

  intercepts <- c(
    as.numeric(fit_stage1$tau1$coefs[1]),
    as.numeric(fit_stage1$tau2$coefs[1]),
    as.numeric(fit_stage1$theta$coefs[1])
  )
  if (any(!is.finite(intercepts))) {
    stop("Non-finite stage-1 EUR LDSC intercept: ", trait1_label, " / ", trait2_label)
  }

  fit_stage2 <- XMAP::estimate_gc(
    data.frame(Z = z1[usable], N = dat1_ldsc$N[usable]),
    data.frame(Z = z2[usable], N = dat2_ldsc$N[usable]),
    ldscore$EUR[usable],
    ldscore$EUR[usable],
    ldscore$EUR_EUR[usable],
    reg_w1 = regression_weight[usable],
    reg_w2 = regression_weight[usable],
    reg_wx = regression_weight[usable],
    constrain_intercept = TRUE,
    intercepts[1],
    intercepts[2],
    intercepts[3]
  )

  estimates <- c(
    Omega11 = as.numeric(fit_stage2$tau1$coefs[2]),
    Omega12 = as.numeric(fit_stage2$theta$coefs[2]),
    Omega22 = as.numeric(fit_stage2$tau2$coefs[2]),
    C11 = as.numeric(fit_stage2$tau1$coefs[1]),
    C12 = as.numeric(fit_stage2$theta$coefs[1]),
    C22 = as.numeric(fit_stage2$tau2$coefs[1])
  )
  if (any(!is.finite(estimates))) {
    stop("Non-finite stage-2 EUR LDSC estimate: ", trait1_label, " / ", trait2_label)
  }
  estimates
}

self_entry_cache <- new.env(parent = emptyenv())

get_self_ldsc_entry <- function(dat, trait_label) {
  if (exists(trait_label, envir = self_entry_cache, inherits = FALSE)) {
    return(get(trait_label, envir = self_entry_cache, inherits = FALSE))
  }
  message("    LDSC diagonal: ", trait_label)
  entry <- estimate_ldsc_entry(dat, dat, trait_label, trait_label)
  assign(trait_label, entry, envir = self_entry_cache)
  entry
}

background_file_for <- function(exposure_code, outcome_abbreviation) {
  file.path(
    background_dir,
    paste0(
      "EUR_", exposure_code,
      "__PanUKB_", outcome_abbreviation,
      "_background.rds"
    )
  )
}

estimate_background <- function(
    exposure_gwas,
    outcome_gwas,
    exposure_code,
    outcome_abbreviation) {
  exposure_label <- paste0("EUR_", exposure_code)
  outcome_label <- paste0("PanUKB_", outcome_abbreviation)
  background_file <- background_file_for(exposure_code, outcome_abbreviation)

  message("    re-estimate background: ", exposure_label, " -> ", outcome_label)
  xx <- get_self_ldsc_entry(exposure_gwas, exposure_label)
  xy <- estimate_ldsc_entry(
    exposure_gwas,
    outcome_gwas,
    exposure_label,
    outcome_label
  )
  yy <- get_self_ldsc_entry(outcome_gwas, outcome_label)

  trait_labels <- c(exposure_label, outcome_label)
  Omega <- matrix(
    c(
      xx[["Omega11"]], xy[["Omega12"]],
      xy[["Omega12"]], yy[["Omega11"]]
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(trait_labels, trait_labels)
  )
  C <- matrix(
    c(
      xx[["C11"]], xy[["C12"]],
      xy[["C12"]], yy[["C11"]]
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(trait_labels, trait_labels)
  )

  if (any(!is.finite(Omega)) || any(diag(Omega) <= 0)) {
    stop("Invalid re-estimated Omega matrix: ", exposure_label, " / ", outcome_label)
  }
  if (any(!is.finite(C)) || any(diag(C) <= 0)) {
    stop("Invalid re-estimated C matrix: ", exposure_label, " / ", outcome_label)
  }

  genetic_correlation <- Omega[1, 2] / sqrt(Omega[1, 1] * Omega[2, 2])
  if (is.finite(genetic_correlation) && abs(genetic_correlation) > 1) {
    Omega[1, 2] <-
      0.95 * sign(genetic_correlation) * sqrt(Omega[1, 1] * Omega[2, 2])
    Omega[2, 1] <- Omega[1, 2]
  }
  tryCatch(
    chol(C),
    error = function(e) stop(
      "Re-estimated C matrix is not positive definite: ",
      exposure_label, " / ", outcome_label
    )
  )

  background <- list(
    C = C,
    Omega = Omega,
    file = background_file,
    exposure_label = exposure_label,
    outcome_label = outcome_label,
    xx = xx,
    xy = xy,
    yy = yy,
    reestimated = TRUE
  )
  saveRDS(background, background_file)
  background
}

prepare_pair <- function(
    exposure_gwas,
    outcome_gwas,
    exposure_code,
    outcome_abbreviation) {
  exposure_label <- paste0("EUR_", exposure_code)
  outcome_label <- paste0("PanUKB_", outcome_abbreviation)

  aligned <- XMR::align_data(
    dat1 = as.data.frame(exposure_gwas),
    dat2 = as.data.frame(outcome_gwas),
    trait1.name = exposure_label,
    trait2.name = outcome_label,
    h2.fix.intercept = FALSE,
    LDSC = FALSE
  )$dat
  if (ncol(aligned) != 9L) stop("Unexpected XMR::align_data output")

  aligned <- as.data.table(aligned)
  setnames(aligned, c(
    "SNP", "A1", "A2", "b.exp", "b.out",
    "se.exp", "se.out", "pval.exp", "pval.out"
  ))
  aligned <- merge(aligned, ldscore_mrapss, by = "SNP")
  if (nrow(aligned) == 0L) {
    stop("No aligned SNP overlaps the EUR LD-score file")
  }

  list(
    data = aligned,
    background = estimate_background(
      exposure_gwas,
      outcome_gwas,
      exposure_code,
      outcome_abbreviation
    ),
    exposure_N = median(exposure_gwas$N, na.rm = TRUE),
    outcome_N = median(outcome_gwas$N, na.rm = TRUE)
  )
}

make_result <- function(
    exposure_code,
    exposure_name,
    exposure_path,
    outcome_abbreviation,
    outcome_name,
    panukb_code,
    outcome_path,
    status,
    error = NA_character_,
    fit = NULL,
    prepared = NULL,
    n_candidates = NA_integer_,
    n_iv = NA_integer_) {
  background <- if (is.null(prepared)) NULL else prepared$background
  data.table(
    analysis_id = analysis_id_for(exposure_code, outcome_abbreviation),
    method = "MRAPSS",
    exposure_dataset = "official_EUR",
    outcome_dataset = "PanUKB_EUR",
    exposure = exposure_code,
    exposure_name = exposure_name,
    exposure_file = exposure_path,
    outcome = outcome_abbreviation,
    outcome_name = outcome_name,
    panukb_code = panukb_code,
    outcome_file = outcome_path,
    threshold = threshold,
    beta = if (is.null(fit)) NA_real_ else scalar_or_na(fit$beta),
    beta_se = if (is.null(fit)) NA_real_ else scalar_or_na(fit$beta.se),
    beta_pvalue = if (is.null(fit)) NA_real_ else scalar_or_na(fit$pvalue),
    tau_sq = if (is.null(fit)) NA_real_ else scalar_or_na(fit$tau.sq),
    pi0 = if (is.null(fit)) NA_real_ else scalar_or_na(fit$pi0),
    IVsignal_sum = if (is.null(fit)) NA_real_ else scalar_or_na(fit$post$IVsignal.sum),
    n_candidates = n_candidates,
    n_IV = n_iv,
    exposure_N = if (is.null(prepared)) NA_real_ else prepared$exposure_N,
    outcome_N = if (is.null(prepared)) NA_real_ else prepared$outcome_N,
    C11 = if (is.null(background)) NA_real_ else background$C[1, 1],
    C12 = if (is.null(background)) NA_real_ else background$C[1, 2],
    C22 = if (is.null(background)) NA_real_ else background$C[2, 2],
    Omega11 = if (is.null(background)) NA_real_ else background$Omega[1, 1],
    Omega12 = if (is.null(background)) NA_real_ else background$Omega[1, 2],
    Omega22 = if (is.null(background)) NA_real_ else background$Omega[2, 2],
    background_reestimated = if (is.null(background)) NA else background$reestimated,
    background_file = if (is.null(background)) NA_character_ else background$file,
    Cor_SelectionBias = TRUE,
    ldscore_population = "EUR",
    clump_reference = "EUR_1000G_phase3_hg38",
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    status = status,
    error = as.character(error)[1L]
  )
}

run_mrapss <- function(
    prepared,
    exposure_code,
    exposure_name,
    exposure_path,
    outcome_abbreviation,
    outcome_name,
    panukb_code,
    outcome_path) {
  candidates <- prepared$data[is.finite(pval.exp) & pval.exp < threshold]
  n_candidates <- nrow(candidates)
  if (n_candidates == 0L) {
    stop("No official-EUR exposure IV candidate at P < 5e-05")
  }

  candidates[, Threshold := threshold]
  clumped <- MRAPSS::clump(
    as.data.frame(candidates),
    IV.Threshold = threshold,
    SNP_col = "SNP",
    pval_col = "pval.exp",
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    bfile = eur_bfile,
    plink_bin = plink_bin
  )
  n_iv <- if (is.null(clumped)) 0L else nrow(clumped)
  if (n_iv < 5L) {
    stop("MRAPSS requires at least five IVs; found ", n_iv)
  }
  clumped$Threshold <- threshold

  mr_data <- as.data.frame(clumped)[, c(
    "SNP", "A1", "A2", "b.exp", "se.exp", "pval.exp",
    "b.out", "se.out", "pval.out", "L2", "Threshold"
  ), drop = FALSE]
  fit <- MRAPSS::MRAPSS(
    mr_data,
    C = prepared$background$C,
    Omega = prepared$background$Omega,
    Cor.SelectionBias = TRUE
  )
  if (!is.finite(scalar_or_na(fit$beta)) || !is.finite(scalar_or_na(fit$beta.se))) {
    stop("MRAPSS returned a non-finite estimate")
  }

  make_result(
    exposure_code = exposure_code,
    exposure_name = exposure_name,
    exposure_path = exposure_path,
    outcome_abbreviation = outcome_abbreviation,
    outcome_name = outcome_name,
    panukb_code = panukb_code,
    outcome_path = outcome_path,
    status = "ok",
    fit = fit,
    prepared = prepared,
    n_candidates = n_candidates,
    n_iv = n_iv
  )
}

results_started <-
  !overwrite_results && file.exists(result_file) && file.info(result_file)$size > 0L
completed <- character(0)
if (results_started) {
  existing_results <- fread(
    result_file,
    na.strings = c("", "NA", "NaN", "NULL", "--")
  )
  required_result_columns <- c("analysis_id", "status")
  missing_result_columns <- setdiff(
    required_result_columns,
    names(existing_results)
  )
  if (length(missing_result_columns) > 0L) {
    stop(
      "Existing result file is missing columns: ",
      paste(missing_result_columns, collapse = ", ")
    )
  }
  completed <- unique(as.character(existing_results[status == "ok", analysis_id]))
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

message("MRAPSS: official EUR exposures -> PanUKB EUR outcomes")
message("Fixed threshold: ", threshold_label)
message("Trait pairs: ", nrow(trait_pairs))
message("Outcome directory: ", outcome_dir)
message("Background mode: force re-estimation for every pending pair")
message("Background directory: ", background_dir)
message("Result file: ", result_file)
message("Successful pairs already present: ", length(completed))

cached_exposure_code <- NA_character_
cached_exposure_gwas <- NULL
cached_exposure_error <- NA_character_

for (pair_index in seq_len(nrow(trait_pairs))) {
  pair <- trait_pairs[pair_index]
  analysis_id <- analysis_id_for(pair$exposure, pair$outcome)
  if (analysis_id %in% completed) {
    message(
      "[", pair_index, "/", nrow(trait_pairs), "] Skip successful: ",
      pair$exposure, " -> PanUKB_", pair$outcome
    )
    next
  }

  message(
    "[", pair_index, "/", nrow(trait_pairs), "] ",
    pair$exposure, " -> PanUKB_", pair$outcome
  )

  if (!identical(cached_exposure_code, pair$exposure)) {
    if (!is.null(cached_exposure_gwas)) rm(cached_exposure_gwas)
    invisible(gc(verbose = FALSE))
    cached_exposure_code <- pair$exposure
    cached_exposure_error <- NA_character_
    cached_exposure_gwas <- tryCatch(
      read_formatted_gwas(
        pair$exposure_path,
        paste0("EUR_", pair$exposure)
      ),
      error = function(e) {
        cached_exposure_error <<- conditionMessage(e)
        NULL
      }
    )
  }

  preparation_error <- cached_exposure_error
  prepared <- NULL
  outcome_gwas <- NULL
  if (!is.null(cached_exposure_gwas)) {
    outcome_gwas <- tryCatch(
      read_formatted_gwas(
        pair$outcome_path,
        paste0("PanUKB_", pair$outcome)
      ),
      error = function(e) {
        preparation_error <<- conditionMessage(e)
        NULL
      }
    )
  }

  if (!is.null(cached_exposure_gwas) && !is.null(outcome_gwas)) {
    prepared <- tryCatch(
      prepare_pair(
        cached_exposure_gwas,
        outcome_gwas,
        pair$exposure,
        pair$outcome
      ),
      error = function(e) {
        preparation_error <<- conditionMessage(e)
        NULL
      }
    )
  }

  result <- if (is.null(prepared)) {
    make_result(
      exposure_code = pair$exposure,
      exposure_name = pair$exposure_name,
      exposure_path = pair$exposure_path,
      outcome_abbreviation = pair$outcome,
      outcome_name = pair$outcome_name,
      panukb_code = pair$panukb_code,
      outcome_path = pair$outcome_path,
      status = "failed",
      error = preparation_error
    )
  } else {
    tryCatch(
      run_mrapss(
        prepared = prepared,
        exposure_code = pair$exposure,
        exposure_name = pair$exposure_name,
        exposure_path = pair$exposure_path,
        outcome_abbreviation = pair$outcome,
        outcome_name = pair$outcome_name,
        panukb_code = pair$panukb_code,
        outcome_path = pair$outcome_path
      ),
      error = function(e) {
        make_result(
          exposure_code = pair$exposure,
          exposure_name = pair$exposure_name,
          exposure_path = pair$exposure_path,
          outcome_abbreviation = pair$outcome,
          outcome_name = pair$outcome_name,
          panukb_code = pair$panukb_code,
          outcome_path = pair$outcome_path,
          status = "failed",
          error = conditionMessage(e),
          prepared = prepared,
          n_candidates = sum(
            is.finite(prepared$data$pval.exp) &
              prepared$data$pval.exp < threshold
          )
        )
      }
    )
  }

  append_result(result)
  completed <- c(completed, analysis_id)
  message("  appended [", result$status, "]")

  if (!is.null(outcome_gwas)) rm(outcome_gwas)
  if (!is.null(prepared)) rm(prepared)
  invisible(gc(verbose = FALSE))
}

if (!is.null(cached_exposure_gwas)) rm(cached_exposure_gwas)
invisible(gc(verbose = FALSE))
message("MRAPSS EUR/PanUKB analysis finished: ", result_file)
