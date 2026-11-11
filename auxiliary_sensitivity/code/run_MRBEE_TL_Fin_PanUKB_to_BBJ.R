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
# MRBEE-TL with FinnGen/Pan-UKB auxiliary data and BBJ target data.
#
# Auxiliary populations: FinnGen / PanUKB
# Target population: BBJ
# Exposures: BMI / Height / Weight
# Outcomes: all non-exposure traits in the 68-trait manifest
# IV thresholds: 5e-5, 5e-6, 5e-7, 5e-8; one master clump at 5e-5 per pair
# LD clumping: 10000 kb, r2 = 0.001
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ieugwasr)
  library(XMAP)
  library(MRBEE)
  library(MRBEEX)
})

set.seed(1)

## =========================================================
## 1. Analysis parameters
## =========================================================

work_dir <- file.path(.xmr_root, "auxiliary_sensitivity")
formatted_dir <- file.path(work_dir, "data/formatted")
manifest_file <- file.path(work_dir, "data/BBJ_PanUKB_FinnGen_68_traits_data_manifest.csv")
# Background caches are shared across IV thresholds for a given trait pair.
result_dir <- file.path(work_dir, "results", "masterClump5e-05")
background_dir <- file.path(work_dir, "results", "background_MRBEE_TL_BBJ")

# Use the LD scores and PLINK references specified in DATA.md.
xmr_dir <- file.path(.xmr_root, "BBJ_TPMI")
reference_dir <- file.path(.xmr_root, "negative_controls")
ldscore_file <- file.path(xmr_dir, "LDscore_eas_brit_afr.txt")
plink_bin <- file.path(reference_dir, "plink_Linux")
eur_bfile <- file.path(reference_dir, "1kg_pops/EUR_phase3_hg38_autosomes")

combinations <- data.table(
  auxiliary = c("Fin", "PanUKB"),
  target = c("BBJ", "BBJ")
)

thresholds <- c(5e-5, 5e-6, 5e-7, 5e-8)
master_clump_threshold <- 5e-5
clump_kb <- 10000
clump_r2 <- 0.001
gxe_p_threshold <- 5e-8
sampling_time <- 1000

# Set these flags to TRUE to regenerate background parameters or result tables.
overwrite_background <- FALSE
overwrite_results <- FALSE
result_file <- file.path(result_dir, "MRBEE_TL_Fin_PanUKB_to_BBJ_68_traits.tsv")

# Read harmonized trait abbreviations from the manifest.
if (!file.exists(manifest_file)) stop("Trait manifest not found: ", manifest_file)
manifest <- fread(manifest_file)
if (!all(c("trait_name", "abbreviation") %in% names(manifest))) {
  stop("Trait manifest must contain trait_name and abbreviation")
}
exposure_names <- c("Body mass index", "Height", "Weight")
exposure_rows <- match(exposure_names, manifest$trait_name)
if (anyNA(exposure_rows)) stop("BMI, Height, or Weight is missing from the trait manifest")
exposures <- manifest$abbreviation[exposure_rows]
outcomes <- manifest[!trait_name %chin% exposure_names, abbreviation]

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(background_dir, recursive = TRUE, showWarnings = FALSE)
if (overwrite_results && file.exists(result_file)) file.remove(result_file)

## =========================================================
## 2. Read formatted GWAS data
## =========================================================

find_gwas_file <- function(prefix, trait) {
  stem <- file.path(formatted_dir, paste0(prefix, "_", trait))
  files <- c(stem, paste0(stem, ".txt"), paste0(stem, ".tsv"), paste0(stem, ".gz"))
  files <- files[file.exists(files)]
  if (length(files) == 0) stop("Missing GWAS file: ", stem)
  files[1]
}

read_gwas <- function(prefix, trait) {
  dat <- fread(find_gwas_file(prefix, trait))
  required <- c("SNP", "A1", "A2", "Z", "N", "P")
  if (!all(required %in% names(dat))) {
    stop(prefix, "_", trait, " lacks columns: ", paste(setdiff(required, names(dat)), collapse = ", "))
  }

  # Construct standardized beta and SE from the formatted Z and N columns.
  dat[, `:=`(
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    Z = as.numeric(Z),
    N = as.numeric(N),
    P = as.numeric(P)
  )]
  # Inputs are QC-filtered; construct standardized model effects.
  dat[, `:=`(beta = Z / sqrt(N), se = 1 / sqrt(N))]
  dat
}

complement_allele <- function(x) chartr("ATCG", "TAGC", x)

# Determine effect-direction flips relative to the reference alleles.
# Exclude A/T and C/G SNPs because EAF is unavailable for strand disambiguation.
allele_sign <- function(ref_a1, ref_a2, a1, a2) {
  valid <- !is.na(ref_a1) & !is.na(ref_a2) & !is.na(a1) & !is.na(a2)
  palindromic <- valid & ref_a1 == complement_allele(ref_a2)
  same <- valid & (
    (a1 == ref_a1 & a2 == ref_a2) |
    (complement_allele(a1) == ref_a1 & complement_allele(a2) == ref_a2)
  )
  reverse <- valid & (
    (a1 == ref_a2 & a2 == ref_a1) |
    (complement_allele(a1) == ref_a2 & complement_allele(a2) == ref_a1)
  )

  sign <- rep(NA_real_, length(ref_a1))
  sign[!palindromic & same] <- 1
  sign[!palindromic & !same & reverse] <- -1
  sign
}

# FEMA uses source and target exposure associations.
harmonise_exposures <- function(source, target) {
  snps <- intersect(source$SNP, target$SNP)
  source <- source[match(snps, source$SNP)]
  target <- target[match(snps, target$SNP)]
  sign <- allele_sign(source$A1, source$A2, target$A1, target$A2)
  keep <- is.finite(sign)

  data.table(
    SNP = source$SNP[keep],
    A1 = source$A1[keep],
    A2 = source$A2[keep],
    b_source = source$beta[keep],
    se_source = source$se[keep],
    z_source = source$Z[keep],
    p_source = source$P[keep],
    b_target = target$beta[keep] * sign[keep],
    se_target = target$se[keep],
    z_target = target$Z[keep] * sign[keep],
    p_target = target$P[keep]
  )
}

# After IV selection, align outcome associations to the source exposure alleles.
add_outcomes <- function(exposure_dat, source_outcome, target_outcome, snps) {
  dat <- exposure_dat[match(snps, exposure_dat$SNP)]
  source_outcome <- source_outcome[match(snps, source_outcome$SNP)]
  target_outcome <- target_outcome[match(snps, target_outcome$SNP)]

  source_sign <- allele_sign(dat$A1, dat$A2, source_outcome$A1, source_outcome$A2)
  target_sign <- allele_sign(dat$A1, dat$A2, target_outcome$A1, target_outcome$A2)
  keep <- is.finite(source_sign) & is.finite(target_sign)

  dat <- dat[keep]
  dat[, `:=`(
    b_out_source = source_outcome$beta[keep] * source_sign[keep],
    se_out_source = source_outcome$se[keep],
    b_out_target = target_outcome$beta[keep] * target_sign[keep],
    se_out_target = target_outcome$se[keep]
  )]
  dat[complete.cases(dat)]
}

## =========================================================
## 3. Estimate source and target background parameters
## =========================================================

# Read LD scores once and add within-population covariance columns.
ldscore <- fread(ldscore_file)
ldscore[, `:=`(
  rsid = as.character(rsid),
  allele1 = toupper(as.character(allele1)),
  allele2 = toupper(as.character(allele2)),
  EUR = as.numeric(EUR),
  EAS = as.numeric(EAS),
  EUR_EUR = as.numeric(EUR),
  EAS_EAS = as.numeric(EAS)
)]
ldscore <- ldscore[
  allele1 != complement_allele(allele2) & !duplicated(rsid)
]

# Align GWAS Z scores to allele1 in the LD-score table.
align_z_to_ldscore <- function(dat, ld) {
  valid <-
    !is.na(dat$A1) & !is.na(dat$A2) &
    !is.na(ld$allele1) & !is.na(ld$allele2)
  same <-
    valid & ((dat$A1 == ld$allele1 & dat$A2 == ld$allele2) |
    (complement_allele(dat$A1) == ld$allele1 & complement_allele(dat$A2) == ld$allele2)
    )
  reverse <-
    valid & ((dat$A1 == ld$allele2 & dat$A2 == ld$allele1) |
    (complement_allele(dat$A1) == ld$allele2 & complement_allele(dat$A2) == ld$allele1)
    )

  z <- dat$Z
  z[!(same | reverse)] <- NA_real_
  z[!same & reverse] <- -z[!same & reverse]
  z
}

# Fit within-population exposure-outcome LDSC and return the 2 x 2 MRBEE Rxy matrix.
estimate_Rxy <- function(exposure_dat, outcome_dat, population) {
  snps <- Reduce(intersect, list(ldscore$rsid, exposure_dat$SNP, outcome_dat$SNP))
  exposure_dat <- exposure_dat[match(snps, exposure_dat$SNP)]
  outcome_dat <- outcome_dat[match(snps, outcome_dat$SNP)]
  ld <- ldscore[match(snps, ldscore$rsid)]

  z1 <- align_z_to_ldscore(exposure_dat, ld)
  z2 <- align_z_to_ldscore(outcome_dat, ld)
  cross_column <- paste0(population, "_", population)

  valid <-
    is.finite(z1) & is.finite(z2) &
    is.finite(exposure_dat$N) & is.finite(outcome_dat$N) &
    is.finite(ld[[population]]) & is.finite(ld[[cross_column]])
  stage1 <- valid & z1^2 < 30 & z2^2 < 30
  if (sum(stage1) < 100) stop("Too few variants for background estimation")

  w1 <- 1 / pmax(ld[[population]], 1)
  w2 <- 1 / pmax(ld[[population]], 1)
  wx <- sqrt(w1 * w2)

  # LDSC stage 1: estimate intercepts freely.
  fit1 <- XMAP::estimate_gc(
    data.frame(Z = z1[stage1], N = exposure_dat$N[stage1]),
    data.frame(Z = z2[stage1], N = outcome_dat$N[stage1]),
    ld[[population]][stage1],
    ld[[population]][stage1],
    ld[[cross_column]][stage1],
    reg_w1 = w1[stage1],
    reg_w2 = w2[stage1],
    reg_wx = wx[stage1],
    constrain_intercept = FALSE
  )
  intercept <- c(fit1$tau1$coefs[1], fit1$tau2$coefs[1], fit1$theta$coefs[1])
  if (any(!is.finite(intercept))) stop("Invalid stage-1 LDSC intercept")

  # LDSC stage 2: retain stage-1 intercepts and refit using all valid SNPs.
  fit2 <- XMAP::estimate_gc(
    data.frame(Z = z1[valid], N = exposure_dat$N[valid]),
    data.frame(Z = z2[valid], N = outcome_dat$N[valid]),
    ld[[population]][valid],
    ld[[population]][valid],
    ld[[cross_column]][valid],
    reg_w1 = w1[valid],
    reg_w2 = w2[valid],
    reg_wx = wx[valid],
    constrain_intercept = TRUE,
    intercept[1], intercept[2], intercept[3]
  )

  Rxy <- matrix(c(
    fit2$tau1$coefs[1], fit2$theta$coefs[1],
    fit2$theta$coefs[1], fit2$tau2$coefs[1]
  ), nrow = 2, byrow = TRUE)

  # Cap estimated correlations outside the valid range.
  if (any(!is.finite(Rxy)) || any(diag(Rxy) <= 0)) stop("Invalid Rxy")
  rho <- Rxy[1, 2] / sqrt(Rxy[1, 1] * Rxy[2, 2])
  if (abs(rho) >= 1) {
    Rxy[1, 2] <- Rxy[2, 1] <- 0.95 * sign(rho) * sqrt(Rxy[1, 1] * Rxy[2, 2])
  }
  Rxy
}

background_file <- function(auxiliary, target, exposure, outcome) {
  file.path(background_dir, paste0(paste(auxiliary, target, exposure, outcome, sep = "_"), ".rds"))
}

get_background <- function(auxiliary, target, exposure, outcome,
                           source_exposure, source_outcome,
                           target_exposure, target_outcome) {
  filename <- background_file(auxiliary, target, exposure, outcome)
  if (file.exists(filename) && !overwrite_background) return(readRDS(filename))

  # Use EUR LD scores for FinnGen/PanUKB and EAS LD scores for BBJ.
  background <- list(
    C_source = estimate_Rxy(source_exposure, source_outcome, "EUR"),
    C_target = estimate_Rxy(target_exposure, target_outcome, "EAS"),
    filename = filename
  )
  saveRDS(background, filename)
  background
}

## =========================================================
## 4. Prepare one exposure-outcome pair
## =========================================================

prepare_pair <- function(auxiliary, target, exposure, outcome) {
  # Read the exposure and outcome GWASs for each population separately.
  source_exposure <- read_gwas(auxiliary, exposure)
  source_outcome <- read_gwas(auxiliary, outcome)
  target_exposure <- read_gwas(target, exposure)
  target_outcome <- read_gwas(target, outcome)

  background <- get_background(
    auxiliary, target, exposure, outcome,
    source_exposure, source_outcome, target_exposure, target_outcome
  )

  list(
    exposure = harmonise_exposures(source_exposure, target_exposure),
    source_outcome = source_outcome,
    target_outcome = target_outcome,
    background = background,
    source_exposure_N = median(source_exposure$N, na.rm = TRUE),
    source_outcome_N = median(source_outcome$N, na.rm = TRUE),
    target_exposure_N = median(target_exposure$N, na.rm = TRUE),
    target_outcome_N = median(target_outcome$N, na.rm = TRUE)
  )
}

## =========================================================
## 5. Run one MRBEE-TL analysis
## =========================================================

wald_p <- function(beta, se) {
  2 * pnorm(abs(beta / se), lower.tail = FALSE)
}

clump_master_ivs <- function(dat) {
  exp_dat <- dat$exposure
  C_source <- dat$background$C_source
  C_target <- dat$background$C_target

  # 1. Select candidate IVs using the joint FEMA exposure chi-square statistic.
  fema_chi2 <- exp_dat$z_source^2 / C_source[1, 1] +
    exp_dat$z_target^2 / C_target[1, 1]
  fema_p <- pchisq(fema_chi2, df = 2, lower.tail = FALSE)
  candidates <- data.table(SNP = exp_dat$SNP, P = fema_p)[P < master_clump_threshold]
  if (nrow(candidates) == 0) stop("No FEMA-significant SNP")

  # 2. Apply 10,000 kb LD clumping with the EUR reference.
  if (nrow(candidates) == 1) {
    fema_iv <- candidates$SNP
  } else {
    clumped <- ieugwasr::ld_clump(
      data.frame(rsid = candidates$SNP, pval = candidates$P),
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p = master_clump_threshold,
      bfile = eur_bfile,
      plink_bin = plink_bin
    )
    fema_iv <- clumped$rsid
  }
  if (length(fema_iv) == 0) stop("No IV remains after clumping")

  list(candidates = candidates, ivs = candidates[match(fema_iv, SNP)])
}

run_mrbee_tl <- function(dat, master_ivs, auxiliary, target, exposure, outcome, threshold) {
  exp_dat <- dat$exposure
  C_source <- dat$background$C_source
  C_target <- dat$background$C_target
  candidates <- master_ivs$candidates[P < threshold]
  if (nrow(candidates) == 0L) stop("No FEMA-significant SNP")
  fema_iv <- master_ivs$ivs[P < threshold, SNP]
  if (length(fema_iv) == 0L) stop("No IV remains in the nested FEMA threshold subset")

  # 3. Remove MR-GxE outliers among source genome-wide significant IVs.
  iv_dat <- exp_dat[match(fema_iv, exp_dat$SNP)]
  source_gws <- which(iv_dat$p_source < gxe_p_threshold)
  gxe_removed <- character(0)

  if (length(source_gws) >= 3) {
    C_gxe <- diag(c(C_source[1, 1], C_target[1, 1]))
    fit_gxe <- MRBEE::MRBEE.IMRP.UV(
      by = iv_dat$b_target[source_gws],
      bx = iv_dat$b_source[source_gws],
      byse = iv_dat$se_target[source_gws],
      bxse = iv_dat$se_source[source_gws],
      Rxy = C_gxe,
      adjust.method = "Sidak"
    )
    delta <- as.numeric(fit_gxe$delta)
    gxe_removed <- iv_dat$SNP[source_gws][which(is.finite(delta) & delta != 0)]
  }

  final_iv <- setdiff(fema_iv, gxe_removed)
  if (length(final_iv) < 3) stop("Fewer than three IVs after MR-GxE filtering")

  # 4. Attach and harmonize source and target outcome associations after IV selection.
  mr_dat <- add_outcomes(
    exp_dat, dat$source_outcome, dat$target_outcome, final_iv
  )
  if (nrow(mr_dat) < 3) stop("Fewer than three harmonised IVs")

  # 5. Estimate the source causal effect from source exposure and outcome data.
  fit_source <- MRBEE::MRBEE.IMRP.UV(
    by = mr_dat$b_out_source,
    bx = mr_dat$b_source,
    byse = mr_dat$se_out_source,
    bxse = mr_dat$se_source,
    Rxy = C_source
  )

  # 6. Estimate the target causal effect using target-population data.
  fit_target <- MRBEE::MRBEE.IMRP.UV(
    by = mr_dat$b_out_target,
    bx = mr_dat$b_target,
    byse = mr_dat$se_out_target,
    bxse = mr_dat$se_target,
    Rxy = C_target
  )

  # 7. Fit the single-exposure transfer model with the source MRBEE estimate as prior.
  fit_tl <- MRBEEX::MRBEE_TL_UV(
    by = mr_dat$b_out_target,
    bX = matrix(mr_dat$b_target, ncol = 1),
    byse = mr_dat$se_out_target,
    bXse = matrix(mr_dat$se_target, ncol = 1),
    Rxy = C_target,
    theta.source = as.numeric(fit_source$theta)[1],
    theta.source.cov = matrix(as.numeric(fit_source$vartheta)[1], nrow = 1),
    LD = "identity",
    pip.thres = 0.5,
    ebic.delta = 0,
    ebic.gamma = 1,
    reliability.thres = 0.5,
    sampling.time = sampling_time
  )

  # 8. Save source, target and transfer-learning estimates.
  result <- data.table(
    method = c("MRBEE_source_specific", "MRBEE_target_specific", "MRBEE_TL_UV"),
    population = c(auxiliary, target, target),
    beta = c(
      as.numeric(fit_source$theta)[1],
      as.numeric(fit_target$theta)[1],
      as.numeric(fit_tl$theta)[1]
    ),
    se = c(
      sqrt(as.numeric(fit_source$vartheta)[1]),
      sqrt(as.numeric(fit_target$vartheta)[1]),
      as.numeric(fit_tl$theta.se)[1]
    )
  )
  if (any(!is.finite(result$beta)) || any(!is.finite(result$se)) || any(result$se <= 0)) {
    stop("MRBEE-TL returned an invalid estimate")
  }
  result[, `:=`(
    analysis_id = paste(auxiliary, target, exposure, outcome, threshold, sep = "|"),
    auxiliary = auxiliary,
    target = target,
    exposure = exposure,
    outcome = outcome,
    threshold = threshold,
    z = beta / se,
    p = mapply(wald_p, beta, se),
    CI_low = beta - 1.96 * se,
    CI_high = beta + 1.96 * se,
    nIV = nrow(mr_dat),
    n_fema_candidates = nrow(candidates),
    n_fema_iv = length(fema_iv),
    n_gxe_removed = length(gxe_removed),
    source_exposure_N = dat$source_exposure_N,
    source_outcome_N = dat$source_outcome_N,
    target_exposure_N = dat$target_exposure_N,
    target_outcome_N = dat$target_outcome_N,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    master_clump_threshold = master_clump_threshold,
    iv_selection_strategy = "nested FEMA-P subsets; threshold-specific MR-GxE filtering",
    gxe_p_threshold = gxe_p_threshold,
    sampling_strategy = c(NA_character_, NA_character_, "half_sample_subsampling"),
    sampling_time = c(NA_integer_, NA_integer_, sampling_time),
    background_file = dat$background$filename,
    status = "ok",
    error = NA_character_
  )]

  setcolorder(result, c(
    "analysis_id", "method", "population", "auxiliary", "target",
    "exposure", "outcome", "threshold", "beta", "se", "z", "p",
    "CI_low", "CI_high", "nIV", "n_fema_candidates", "n_fema_iv",
    "n_gxe_removed", "source_exposure_N", "source_outcome_N",
    "target_exposure_N", "target_outcome_N", "clump_kb", "clump_r2",
    "master_clump_threshold", "iv_selection_strategy",
    "gxe_p_threshold", "sampling_strategy", "sampling_time",
    "background_file", "status", "error"
  ))
  result
}

## =========================================================
## 6. Run analyses and append results
## =========================================================

make_failure <- function(auxiliary, target, exposure, outcome, threshold, error) {
  data.table(
    analysis_id = paste(auxiliary, target, exposure, outcome, threshold, sep = "|"),
    method = "MRBEE_TL_UV", population = target,
    auxiliary, target, exposure, outcome, threshold,
    beta = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
    CI_low = NA_real_, CI_high = NA_real_, nIV = NA_integer_,
    n_fema_candidates = NA_integer_, n_fema_iv = NA_integer_,
    n_gxe_removed = NA_integer_, source_exposure_N = NA_real_,
    source_outcome_N = NA_real_, target_exposure_N = NA_real_,
    target_outcome_N = NA_real_, clump_kb, clump_r2,
    master_clump_threshold = master_clump_threshold,
    iv_selection_strategy = "nested FEMA-P subsets; threshold-specific MR-GxE filtering",
    gxe_p_threshold,
    sampling_strategy = "half_sample_subsampling", sampling_time,
    background_file = NA_character_, status = "failed", error = as.character(error)[1]
  )
}

append_result <- function(result) {
  has_result <- file.exists(result_file) && file.info(result_file)$size > 0
  fwrite(
    result, result_file, sep = "\t", quote = FALSE, na = "NA",
    append = has_result, col.names = !has_result
  )
}

completed <- character(0)
if (file.exists(result_file) && file.info(result_file)$size > 0) {
  completed <- unique(fread(result_file, select = "analysis_id")$analysis_id)
}

for (i in seq_len(nrow(combinations))) {
  auxiliary <- combinations$auxiliary[i]
  target <- combinations$target[i]

  for (exposure in exposures) {
    for (outcome in outcomes) {
      message("Preparing: ", auxiliary, " -> ", target, " | ", exposure, " -> ", outcome)

      pair_analysis_ids <- paste(auxiliary, target, exposure, outcome, thresholds, sep = "|")
      if (all(pair_analysis_ids %in% completed)) next

      preparation_error <- NULL
      prepared <- tryCatch(
        prepare_pair(auxiliary, target, exposure, outcome),
        error = function(e) {
          preparation_error <<- conditionMessage(e)
          NULL
        }
      )

      clump_error <- NULL
      master_ivs <- if (is.null(prepared)) NULL else tryCatch(
        clump_master_ivs(prepared),
        error = function(e) {
          clump_error <<- conditionMessage(e)
          NULL
        }
      )

      for (threshold in thresholds) {
        analysis_id <- paste(auxiliary, target, exposure, outcome, threshold, sep = "|")
        if (analysis_id %in% completed) next

        result <- if (is.null(prepared)) {
          make_failure(auxiliary, target, exposure, outcome, threshold, preparation_error)
        } else if (is.null(master_ivs)) {
          make_failure(auxiliary, target, exposure, outcome, threshold, clump_error)
        } else {
          tryCatch(
            run_mrbee_tl(prepared, master_ivs, auxiliary, target, exposure, outcome, threshold),
            error = function(e) make_failure(
              auxiliary, target, exposure, outcome, threshold, conditionMessage(e)
            )
          )
        }

        append_result(result)
        completed <- c(completed, analysis_id)
        message("  ", threshold, " appended [", result$status[1], "]")
      }

      rm(prepared, master_ivs)
      invisible(gc())
    }
  }
}

message("MRBEE-TL finished: ", result_file)
